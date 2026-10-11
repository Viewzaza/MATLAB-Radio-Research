function [T, info] = ff_signal_id(varargin)
%FF_SIGNAL_ID  Find the signals in a FieldFox spectrum and guess what they are.
%
%   T = ff_signal_id(f, P_dBm, rbw)       f in Hz, P in dBm (one trace, or sweeps x points), rbw in Hz
%   T = ff_signal_id(f, P_dBm, rbw, Name, Value, ...)
%   T = ff_signal_id("10.90.36.200")      read the trace the FieldFox shows now (SA mode) and identify it
%   T = ff_signal_id(ff)                   same, through an open FieldFox object (if a tool
%                                          holds ff.takeControl(), one sweep is triggered first)
%   T = ff_signal_id("rec.mat")            a recording (f, P, settings) from fieldfox_live or another tool
%   [T, info] = ff_signal_id(...)          also noise floor, signal regions, class probabilities
%   ff_signal_id(...)                      no output: print the table (and notes)
%
% The FieldFox gives power spectra only (no I/Q option 351), so signals are
% recognised by the width and shape of their spectrum - no demodulation:
%   1. Noise floor: robust low-order fit through the noise-only points of the
%      5-point smoothed trace, plus a running median for slow ripple. Signals
%      wider than ~60 % of the span hide the floor: widen the span.
%   2. Detection: findpeaks (height and prominence above the floor, thresholds
%      from the noise's own upper tail) on the raw and the smoothed trace. Each
%      peak grows into the region above the floor; touching regions are merged
%      (the flat top of an OFDM carrier is ONE signal) unless a dip of
%      'SplitDip' dB separates two peaks.
%   3. Measurement: 99 % occupied bandwidth obw(sxx,f,rbw,lims,99) and -26 dB
%      bandwidth powerbw(sxx,f,rbw,lims,26), with the noise power subtracted.
%   4. Identification (ff_signal_classifier, trained on synthetic FieldFox traces
%      made from the reference spectra of ff_signal_templates):
%        "template"  each reference spectrum, seen through the same RBW and
%                    detector, is fitted to the trace; smallest dB error wins
%        "ml"        bagged decision trees (fitcensemble) on bandwidth/shape features
%        "cnn"       1-D convolutional network (Deep Learning Toolbox) on the shape
%        "ensemble"  (default) weighted geometric mean of the three (weights fitted
%                    on validation traces), times a band prior: a technology not
%                    used in the band where the signal sits is made 10x less likely
%        "none"      detection and measurement only (fast)
%      ff_signal_id never trains: without ff_signal_models.mat, "ensemble" falls
%      back to the template fit (warning) - run ff_signal_classifier once (~4 min).
%   5. Band (Thailand + 3GPP band table) and channel: Wi-Fi channel
%      (wlanChannelFrequency), NR-ARFCN (TS 38.104 global raster), LTE EARFCN
%      (TS 36.101), GSM ARFCN, BLE/Bluetooth channel, UHF TV channel, FM frequency.
%
% Output table T (one row per signal, in frequency order):
%   Freq_Hz     centre of the 99 % occupied band
%   Level_dBm   highest trace point of the signal (display level, RBW-dependent)
%   OBW_Hz      99 % occupied bandwidth;  BW26_Hz  -26 dB bandwidth
%   Band        band(s) the centre falls in;  Channel  channel number (approximate)
%   Technology  best guess (e.g. "LTE 10 MHz"), "Unresolved" or "Unknown"
%   Score       probability of the best guess (0-1), from synthetic training data
%   Top3        three best guesses with probabilities
%   Method      method used (info.MethodDetail says which parts were combined)
%
% Read the result as a hypothesis, not a verdict:
%   - LTE and NR of the same bandwidth have almost the same spectrum (Wi-Fi 20 MHz
%     is close too); the band prior decides most real cases.
%   - Templates are FULLY LOADED carriers; a quiet cell or bursty Wi-Fi looks
%     weaker and more ragged. Max hold or averaging helps.
%   - "CW" means narrower than the RBW (a carrier, a harmonic, narrowband voice).
%     Lower the RBW to see its real shape.
%   - "Unresolved": fewer than 'MinPoints' trace points across its 99 % band, or
%     it runs into the edge of the span -> centre a smaller span on it.
%   - "Unknown": it fits the template of its best guess clearly worse than the
%     synthetic signals of that class did (inverter switching noise, LoRa, radar,
%     analogue video, two carriers side by side, ...). Top3 still lists the
%     closest classes. Very noisy or bursty traces are judged leniently.
%   - Classes known: NR 5/10/20/40/100 MHz, LTE 1.4-20 MHz, Wi-Fi 20/40/80 MHz,
%     BLE 1M/2M, Bluetooth BR/EDR, GSM, FM broadcast, GPS L1 C/A, DVB-T2 8 MHz,
%     CW. Anything else gets the nearest of these (or "Unknown").
%   - Wide signals need a span of ~1.5x their width so the floor is visible.
%   - The first call in a MATLAB session loads the models (~10 s); later calls
%     take ~0.1-0.5 s (more with many signals).
%   - GPS is ~-130 dBm, below the FieldFox noise floor: a visible "GPS" signal is
%     a repeater/jammer or a misclassification.
%   - Channel numbers come from the centre estimate (+- one trace point); LTE
%     lists every band whose range contains it (e.g. B5 DL / B8 UL).
%   - The accuracies printed by ff_signal_classifier are on SYNTHETIC traces:
%     an upper bound. Verify important findings by narrowing span and RBW.
%
% Name-Value options (defaults):
%   'Method'      "ensemble"  "template" | "ml" | "cnn" | "ensemble" | "none"
%   'Detector'    "auto"      trace detector: "peak" | "average" | "sample"
%                             (auto: read from the FieldFox, else "peak" = FieldFox default)
%   'Threshold'   6           dB above the noise floor to detect a signal
%   'NoiseFloor'  []          dBm (scalar or per point) instead of the estimate
%   'MinPoints'   5           fewer points across a signal -> "Unresolved"
%   'SplitDip'    10          dB dip that separates two neighbouring signals
%   'MaxSignals'  25          keep the strongest ones
%   'BandPrior'   true        use the band table in "ensemble"
%   'Combine'     "mean"      several sweeps: "mean" (power) | "max" | "last"
%   'Plot'        false       labelled spectrum
%   'Visible'     'on'        figure visibility ('off' for batch use)
%   'Prompt'      true        accepted for suite compatibility (no dialogs are used)
%   'Port'        5025        live use only
%   'SaveTo'      ""          live use: also save the trace as a recording (.mat)
%   'Models'      ""          models MAT (default: ff_signal_models.mat next to the code)
%
% Live use is PASSIVE: the instrument settings are not changed (only the data
% format, which is put back). The FieldFox must be in SA mode.
%
% See also ff_signal_classifier, ff_signal_templates, fieldfox_live.

opt = struct('Method', "ensemble", 'Detector', "auto", 'Threshold', 6, 'NoiseFloor', [], ...
    'MinPoints', 5, 'SplitDip', 10, 'MaxSignals', 25, 'BandPrior', true, 'Combine', "mean", ...
    'Plot', false, 'Visible', 'on', 'Prompt', true, 'Port', 5025, 'SaveTo', "", 'Models', "");
if nargin < 1
    error('Give a trace (f, P_dBm, rbw), an IP address, a FieldFox object or a recording file.');
end
src = varargin{1};
meta = struct('settings', [], 'idn', "", 't', NaT);
if isnumeric(src)
    if nargin < 3, error('Usage: ff_signal_id(f, P_dBm, rbw, Name, Value, ...)'); end
    opt = parseOptions(opt, varargin(4:end));
else
    opt = parseOptions(opt, varargin(2:end));
end
opt.Method = lower(string(opt.Method));
if ~ismember(opt.Method, ["template" "ml" "cnn" "ensemble" "none"])
    error('Method must be "template", "ml", "cnn", "ensemble" or "none".');
end
if isnumeric(src)
    [f, P, rbw] = varargin{1:3};
    det = opt.Detector;
else
    if ~isa(src, 'FieldFox') && endsWith(lower(string(src)), ".mat")
        [f, P, rbw, det, meta] = readRecording(src);
    else
        [f, P, rbw, det, meta] = readLive(src, opt);
    end
    if opt.Detector ~= "auto", det = opt.Detector; end
end
[f, P] = prepareTrace(f, P, opt.Combine);
det = detectorName(det);
if ~(isscalar(rbw) && rbw > 0), error('rbw must be the resolution bandwidth in Hz.'); end

%% Detect and measure
[segs, nf, noise] = detectSignals(f, P, rbw, opt);
n = numel(segs);
df = (f(end) - f(1)) / (numel(f) - 1);

%% Identify
% Never trains here (that takes minutes): without cached models "ensemble"
% falls back to the template fit; run ff_signal_classifier once to train.
classes = strings(1, 0); prob = zeros(n, 0); parts = struct(); fitBad = false(n, 1);
detail = opt.Method;
if n > 0 && opt.Method ~= "none"
    margs = {'Train', false, 'Verbose', false};
    if strlength(string(opt.Models)) > 0, margs = [margs, {'File', opt.Models}]; end
    M = ff_signal_classifier(margs{:});
    if ~M.Trained && ismember(opt.Method, ["ml" "cnn"])
        error('ff_signal_id:noModels', ['Method "%s" needs the trained models: run ' ...
            'ff_signal_classifier once (a few minutes), or use Method "template".'], opt.Method);
    end
    if ~M.Trained && opt.Method == "ensemble"
        warnOnce('ff_signal_id:noModels', ['No trained models for these templates: using the ' ...
            'template fit only. Run ff_signal_classifier once to train them.']);
    end
    classes = M.Classes;
    parts = M.predict(f, P, nf, rbw, segs, det, noise);
    [prob, detail] = combine(parts, opt.Method, M, segs, opt.BandPrior);
    % "Unknown": the trace fits the template of its best guess clearly worse
    % than genuine signals of that class did in training
    [~, kBest] = max(prob, [], 2);
    lim = M.FitLimit(min(kBest, numel(M.FitLimit)));
    fitBad = parts.FitAll(sub2ind(size(parts.FitAll), (1:n).', kBest)) > lim(:);
end

%% Table
Freq_Hz = [segs.fc].'; Level_dBm = [segs.level].'; OBW_Hz = [segs.obw].'; BW26_Hz = [segs.bw26].';
Band = strings(n, 1); Channel = strings(n, 1); Technology = strings(n, 1);
Score = nan(n, 1); Top3 = strings(n, 1); Method = repmat(opt.Method, n, 1);
if n == 0, Freq_Hz = zeros(0, 1); Level_dBm = Freq_Hz; OBW_Hz = Freq_Hz; BW26_Hz = Freq_Hz; end
for i = 1:n
    B = bandsAt(Freq_Hz(i));
    if ~isempty(B), Band(i) = strjoin(B.Name, " / "); end
    if isempty(prob)
        Technology(i) = "";
        guess = "";
    else
        [p, order] = sort(prob(i, :), 'descend');
        k = order(1:min(3, end));
        Top3(i) = strjoin(compose("%s %.2f", classes(k).', p(1:numel(k)).'), "; ");
        Score(i) = p(1);
        guess = classes(order(1));
        Technology(i) = guess;
        if segs(i).npts < opt.MinPoints || segs(i).edge
            Technology(i) = "Unresolved";
        elseif fitBad(i)
            Technology(i) = "Unknown";
        end
    end
    Channel(i) = channelLabel(Freq_Hz(i), guess, OBW_Hz(i));
end
T = table(Freq_Hz, Level_dBm, OBW_Hz, BW26_Hz, Band, Channel, Technology, Score, Top3, Method);

%% Extras
info = struct('f', f, 'P', P, 'Floor', nf, 'NoiseSigma_dB', noise.sigma, 'Noise', noise, 'RBW', rbw, ...
    'Detector', det, 'MethodDetail', detail, 'Segments', segs, 'Classes', classes, 'Probabilities', prob, ...
    'Parts', parts, 'Comb', findComb(segs, rbw, df), 'Settings', meta.settings, ...
    'IDN', meta.idn, 'Time', meta.t, 'Figure', []);
if opt.Plot, info.Figure = plotResult(f, P, nf, segs, T, rbw, opt.Visible); end
if strlength(string(opt.SaveTo)) > 0 && ~isempty(meta.settings)
    out = struct('f', f, 'P', P, 't', meta.t, 'settings', meta.settings, 'idn', meta.idn);
    save(opt.SaveTo, '-struct', 'out');
    fprintf('Trace saved to %s\n', opt.SaveTo);
end
if nargout == 0
    printResult(T, info);
    clear T
end
end

%% ===== input =================================================================

function opt = parseOptions(opt, nv)
% Case-insensitive Name, Value pairs (like the other suite tools).
if mod(numel(nv), 2), error('Options come in Name, Value pairs.'); end
names = fieldnames(opt);
for k = 1:2:numel(nv)
    j = find(strcmpi(names, nv{k}), 1);
    if isempty(j), error('Unknown option "%s".', nv{k}); end
    opt.(names{j}) = nv{k + 1};
end
end

function [f, P] = prepareTrace(f, P, how)
f = double(f(:).'); P = double(P);
if size(P, 2) ~= numel(f) && size(P, 1) == numel(f), P = P.'; end
if size(P, 2) ~= numel(f), error('P must have one column per frequency point.'); end
if numel(f) < 5, error('The trace needs at least 5 points.'); end
if any(diff(f) <= 0)
    [f, k] = sort(f);
    if any(diff(f) <= 0), error('The frequencies must be distinct.'); end
    P = P(:, k);
end
if all(~isfinite(P(:))), error('The trace has no finite values.'); end
P(isinf(P)) = NaN;
if size(P, 1) > 1
    switch lower(string(how))
        case "max",  P = max(P, [], 1);
        case "last", P = P(end, :);
        otherwise,   P = 10*log10(mean(10.^(P/10), 1, 'omitnan'));
    end
end
if any(isnan(P)), P = fillmissing(P, 'nearest'); end
end

function d = detectorName(d)
d = upper(string(d));
if startsWith(d, "AVER") || d == "RMS", d = "average";
elseif startsWith(d, "SAMP") || startsWith(d, "NEG"), d = "sample";
else, d = "peak";             % AUTO / NORMal / POSitive: peak-like
end
end

function [f, P, rbw, det, meta] = readRecording(file)
S = load(file);
if ~all(isfield(S, {'f', 'P'})), error('%s is not a recording (needs f and P).', file); end
if numel(S.f) < 5
    error(['%s is not a swept spectrum (%d frequency point): power meter and VVM logs are ' ...
        're-analysed with their own tool (''File'').'], file, numel(S.f));
end
f = S.f; P = S.P; det = "auto";
meta = struct('settings', [], 'idn', "", 't', NaT);
if isfield(S, 'settings'), meta.settings = S.settings; end
if isfield(S, 'idn'), meta.idn = string(S.idn); end
if isfield(S, 't') && ~isempty(S.t), meta.t = S.t(end); end
if isfield(S, 'settings') && isfield(S.settings, 'rbw') && isfinite(S.settings.rbw)
    rbw = S.settings.rbw;
else                                % e.g. an on-instrument recording parsed without 'RBW'
    error('%s has no settings.rbw; call ff_signal_id(f, P, rbw) instead.', file);
end
if isfield(S, 'settings') && isfield(S.settings, 'detector'), det = S.settings.detector; end
end

function [f, P, rbw, det, meta] = readLive(src, opt)
% One passive read of the trace on screen. Nothing is changed except the
% data format during the transfer, which is put back.
if isa(src, 'FieldFox')
    ff = src; own = false;
else
    fprintf('Connecting to %s:%d ...\n', string(src), opt.Port);
    ff = FieldFox(string(src), opt.Port); own = true;
end
if ff.mode() ~= "SA"
    error('The FieldFox is in %s mode. Switch it to SA (Spectrum Analyzer) mode first.', ff.mode());
end
if own
    oldFormat = ff.query("FORM:DATA?");
    restore = onCleanup(@() restoreFormat(ff, oldFormat));
end
s = ff.saSettings();
ff.sweep();                        % not in control: just waits for the current sweep
P = ff.saTrace(1);
f = linspace(s.start, s.stop, numel(P));
rbw = s.rbw; det = s.detector;
e = ff.errors();
if ~isempty(e), warning('ff_signal_id:instrument', 'Instrument reported: %s', strjoin(e, '; ')); end
meta = struct('settings', rmfield(s, 'f'), 'idn', ff.IDN, 't', datetime('now'));
fprintf('Read %d points, %s - %s, RBW %s, detector %s\n', numel(P), fmtHz(s.start), fmtHz(s.stop), fmtHz(rbw), s.detector);
end

function restoreFormat(ff, oldFormat)
try
    if startsWith(oldFormat, "ASC"), ff.write("FORM:DATA ASC,0");
    elseif contains(oldFormat, "64"), ff.write("FORM:DATA REAL,64");
    end
catch err
    warning('ff_signal_id:format', 'Could not restore FORM:DATA %s: %s', oldFormat, err.message);
end
end

%% ===== detection ===========================================================

function [segs, nf, nz] = detectSignals(f, P, rbw, opt)
N = numel(P); df = (f(end) - f(1)) / (N - 1);
w = 2*min(2, floor(N/30)) + 1;                   % 5 points (fewer on very short traces)
Ps = 10*log10(movmean(10.^(P/10), w));           % power average: much less scatter
if isempty(opt.NoiseFloor)
    [nf, isN] = noiseFloor(P, Ps);
else
    nf = double(opt.NoiseFloor(:).') .* ones(1, N);
    r = Ps - nf;
    isN = r < max(3*1.4826*mad(r, 1), 1.5);
    if nnz(isN) < 3, isN = r <= median(r); end
end
res = P - nf; rs = Ps - nf;
sig1 = max(1.4826*mad(res(isN), 1), 0.3);
nz = struct('isNoise', isN, 'sigma', sig1, 'meanRatio', mean(10.^(res(isN)/10)));
nz.mae = mean(abs(res(isN) - 10*log10(nz.meanRatio)));    % dB scatter of noise around its mean
% Thresholds from the upper tail of the noise itself: a single-sweep sample or
% average trace scatters by +-5 dB, an averaged one by a fraction of a dB.
thr1 = max(opt.Threshold + 3, prctile(res(isN), 99.5) + 3);
thr5 = max(opt.Threshold, prctile(rs(isN), 99.5) + 2);
% (padded so that a signal running into the trace edge is still found)
l1 = []; l5 = [];
if max(res) > thr1, [~, l1] = findpeaks([-1e3 res -1e3], 'MinPeakHeight', thr1, 'MinPeakProminence', thr1); end
if max(rs) > thr5,  [~, l5] = findpeaks([-1e3 rs -1e3], 'MinPeakHeight', thr5, 'MinPeakProminence', thr5); end
seeds = unique([l1(:); l5(:)]).' - 1;

% Regions above the floor; close gaps of up to 2 points (ragged plateaus)
above = rs > max(prctile(rs(isN), 90), 1.5) | res > thr1;
above = movmin(double(movmax(double(above), 5)), 5) > 0;
lab = cumsum([above(1), diff(above) == 1]) .* above;
runs = unique(lab(seeds)); runs(runs == 0) = [];

segs = struct('lo', {}, 'hi', {}, 'a', {}, 'b', {}, 'fc', {}, 'flo', {}, 'fhi', {}, ...
    'obw', {}, 'bw26', {}, 'level', {}, 'ipk', {}, 'npts', {}, 'snr', {}, 'edge', {});
sLin = max(10.^(P/10) - nz.meanRatio * 10.^(nf/10), 0);   % noise power removed
e = max(2, ceil(rbw/df));
for r = runs
    idx = find(lab == r); lo = idx(1); hi = idx(end);
    cuts = [];
    if hi - lo >= 4
        [~, lp] = findpeaks(rs(lo:hi), 'MinPeakProminence', opt.SplitDip);
        lp = lp + lo - 1;
        for j = 1:numel(lp) - 1
            [~, m] = min(rs(lp(j):lp(j + 1)));
            cuts(end + 1) = lp(j) + m - 1; %#ok<AGROW>
        end
    end
    edges = [lo, cuts + 1; cuts, hi];
    for p = 1:size(edges, 2)
        s = measure(f, P, nf, sLin, rbw, edges(1, p), edges(2, p), e, df);
        if ~isempty(s), segs(end + 1) = s; end %#ok<AGROW>
    end
end
if numel(segs) > opt.MaxSignals
    [~, k] = sort([segs.level], 'descend');
    segs = segs(sort(k(1:opt.MaxSignals)));
end
end

function s = measure(f, P, nf, sLin, rbw, lo, hi, e, df)
N = numel(P);
a = max(1, lo - e); b = min(N, hi + e);
s = [];
if b - a < 2 || ~any(sLin(lo:hi) > 0), return, end
lims = [f(a) f(b)];
% obw/powerbw refuse rbw < point spacing (wide span, small RBW); the
% bandwidths are power ratios, so the rbw scale does not change them.
r = max(rbw, max(diff(f)));
try
    [bw, flo, fhi] = obw(sLin(:), f(:), r, lims, 99);
    bw26 = powerbw(sLin(:), f(:), r, lims, 26);
catch
    bw = NaN; bw26 = NaN;
end
if ~(bw > 0), bw = (hi - lo + 1)*df; flo = f(lo) - df/2; fhi = f(hi) + df/2; end
if ~(bw26 > 0), bw26 = bw; end
[lev, ipk] = max(P(lo:hi)); ipk = ipk + lo - 1;
s = struct('lo', lo, 'hi', hi, 'a', a, 'b', b, 'fc', (flo + fhi)/2, 'flo', flo, 'fhi', fhi, ...
    'obw', bw, 'bw26', bw26, 'level', lev, 'ipk', ipk, 'npts', bw/df, 'snr', lev - nf(ipk), ...
    'edge', lo == 1 || hi == N);
end

function [nf, isN] = noiseFloor(P, Ps)
% Robust fit to the noise-only points of the smoothed trace Ps: start from
% its lowest 40 %, fit a low-order polynomial (floor tilt), keep the points
% within 3 sigma of it, repeat; then follow slow ripple with a running median
% of the noise points only (signals are bridged, not followed). Finally the
% floor is put at the median of the raw noise points. Works while signals
% cover less than ~60 % of the span.
N = numel(P); x = linspace(-1, 1, N);
isN = Ps <= prctile(Ps, 40);
for it = 1:8
    d = min(3, floor(nnz(isN)/15));
    base = polyval(polyfit(x(isN), Ps(isN), d), x);
    r = Ps - base;
    sig = max(1.4826*mad(r(isN), 1), 0.2);
    new = r < max(3*sig, 1.5);
    if isequal(new, isN) || nnz(new) < 5, break, end
    isN = new;
end
r(~isN) = NaN;
loc = movmedian(r, max(9, round(N/8)), 'omitnan');
nf = base + fillmissing(loc, 'linear', 'EndValues', 'nearest');
nf = nf + median(P(isN) - nf(isN));
end

%% ===== identification ======================================================

function [prob, used] = combine(parts, method, M, segs, useBand)
% Weighted geometric mean of the available methods (each floored at 1e-3;
% weights fitted on validation traces by ff_signal_classifier), times the band prior.
L = {parts.Template, parts.ML, parts.CNN};
names = ["template" "ml" "cnn"];
w = M.Weights;
if method ~= "ensemble", w = double(names == method); end
ok = ~cellfun(@isempty, L) & w > 0;
if ~any(ok), ok = names == "template"; w = [1 0 0]; end   % ensemble of untrained models
used = strjoin(names(ok), "+");
lp = zeros(numel(segs), numel(M.Classes));
for k = find(ok), lp = lp + w(k) * log(max(L{k}, 1e-3)); end
if method == "ensemble"
    if useBand, used = used + "+band"; end
    used = "ensemble(" + used + ")";
end
if method == "ensemble" && useBand
    for i = 1:numel(segs), lp(i, :) = lp(i, :) + log(bandPrior(segs(i).fc, M)); end
end
prob = exp(lp - max(lp, [], 2));
prob = prob ./ sum(prob, 2);
end

function pr = bandPrior(fc, M)
% A technology only exists in its own frequency ranges: outside them its prior
% is 1e-4 (a Bluetooth-shaped hump at 96 MHz is not Bluetooth). Inside a listed
% band, technologies not used there (or wider than the band) get 0.1.
hard = {   % technology, allowed ranges [low high] in Hz
    "BLE",       [2400e6 2483.5e6]
    "Bluetooth", [2400e6 2483.5e6]
    "Wi-Fi",     [2400e6 2495e6; 5150e6 5895e6; 5925e6 7125e6]
    "GNSS",      [1164e6 1300e6; 1559e6 1610e6]
    "FM",        [76e6 108e6]
    "DVB-T2",    [174e6 230e6; 470e6 862e6]
    "GSM",       [380e6 1990e6]
    "LTE",       [410e6 7125e6]
    "NR",        [410e6 7125e6] };
pr = ones(1, numel(M.Classes));
for r = 1:size(hard, 1)
    rg = hard{r, 2};
    if ~any(fc >= rg(:, 1) & fc <= rg(:, 2)), pr(M.ClassTech == hard{r, 1}) = 1e-4; end
end
B = bandsAt(fc);
if isempty(B), return, end
ok = ismember(M.ClassTech, [B.Tech{:}, "CW"]) & M.ClassBW <= 1.1*max(B.High - B.Low);
ok(M.ClassTech == "CW") = true;
pr(~ok) = min(pr(~ok), 0.1);
end

%% ===== bands and channels ==================================================

function B = bandsAt(fc)
persistent tbl
if isempty(tbl)
    c = {   % name, low (Hz), high (Hz), technologies
        "MW AM broadcast",                 526.5e3,  1606.5e3, "CW"
        "Amateur 40 m",                    7.0e6,    7.2e6,    "CW"
        "CB 27 MHz",                       26.965e6, 27.405e6, "CW"
        "FM broadcast",                    87.5e6,   108e6,    "FM"
        "Airband VHF",                     108e6,    137e6,    "CW"
        "Amateur 2 m",                     144e6,    146e6,    "CW"
        "Marine VHF",                      156e6,    162.05e6, "CW"
        "Thai CB 245 MHz",                 245e6,    246e6,    "CW"
        "Amateur 70 cm",                   430e6,    440e6,    "CW"
        "ISM/SRD 433 MHz",                 433.05e6, 434.79e6, "CW"
        "Digital TV (UHF ch 21-48)",       470e6,    694e6,    "DVB-T2"
        "n28/B28 UL",                      703e6,    748e6,    ["LTE" "NR"]
        "n28/B28 DL",                      758e6,    803e6,    ["LTE" "NR"]
        "n5/B5 UL",                        824e6,    849e6,    ["LTE" "NR"]
        "n5/B5 DL",                        869e6,    894e6,    ["LTE" "NR"]
        "n8/B8 UL (GSM 900 UL)",           880e6,    915e6,    ["LTE" "NR" "GSM"]
        "IoT/RFID 920-925 MHz",            920e6,    925e6,    "CW"
        "n8/B8 DL (GSM 900 DL)",           925e6,    960e6,    ["LTE" "NR" "GSM"]
        "GNSS L1 (GPS/Galileo/BeiDou)",    1559e6,   1610e6,   "GNSS"
        "n3/B3 UL (GSM 1800 UL)",          1710e6,   1785e6,   ["LTE" "NR" "GSM"]
        "n3/B3 DL (GSM 1800 DL)",          1805e6,   1880e6,   ["LTE" "NR" "GSM"]
        "n1/B1 UL",                        1920e6,   1980e6,   ["LTE" "NR"]
        "n1/B1 DL",                        2110e6,   2170e6,   ["LTE" "NR"]
        "n40/B40 TDD",                     2300e6,   2400e6,   ["LTE" "NR"]
        "ISM 2.4 GHz (Wi-Fi/Bluetooth)",   2400e6,   2483.5e6, ["Wi-Fi" "BLE" "Bluetooth"]
        "n41/B41 TDD",                     2496e6,   2690e6,   ["LTE" "NR"]
        "B38 TDD",                         2570e6,   2620e6,   ["LTE" "NR"]
        "n78 TDD",                         3300e6,   3800e6,   "NR"
        "Wi-Fi 5 GHz",                     5150e6,   5850e6,   "Wi-Fi"
        "Wi-Fi 6 GHz",                     5925e6,   6425e6,   "Wi-Fi" };
    tbl = struct('Name', [c{:, 1}], 'Low', [c{:, 2}], 'High', [c{:, 3}]);
    tbl.Tech = c(:, 4).';
end
k = fc >= tbl.Low & fc <= tbl.High;
if ~any(k), B = []; return, end
B = struct('Name', tbl.Name(k), 'Low', tbl.Low(k), 'High', tbl.High(k));
B.Tech = tbl.Tech(k);
end

function ch = channelLabel(fc, cls, bw)
% Channel number for the guessed technology (approximate: from the centre estimate).
ch = "";
tech = extractBefore(cls + " ", " ");
switch tech
    case "NR"
        if fc < 3e9, n = round(fc/5e3); else, n = 600000 + round((fc - 3e9)/15e3); end
        if fc > 700e6 && fc < 3.8e9, ch = sprintf("NR-ARFCN ~%d", n); end
    case "LTE"
        % band, F_low (Hz), N_offs, width (Hz)   (TS 36.101 Table 5.7.3-1)
        t = {"B1 DL", 2110e6, 0, 60e6;      "B1 UL", 1920e6, 18000, 60e6
             "B3 DL", 1805e6, 1200, 75e6;   "B3 UL", 1710e6, 19200, 75e6
             "B5 DL", 869e6, 2400, 25e6;    "B5 UL", 824e6, 20400, 25e6
             "B8 DL", 925e6, 3450, 35e6;    "B8 UL", 880e6, 21450, 35e6
             "B28 DL", 758e6, 9210, 45e6;   "B28 UL", 703e6, 27210, 45e6
             "B38", 2570e6, 37750, 50e6;    "B40", 2300e6, 38650, 100e6
             "B41", 2496e6, 39650, 194e6};
        c = strings(0);
        for i = 1:size(t, 1)
            if fc >= t{i, 2} && fc <= t{i, 2} + t{i, 4}
                c(end + 1) = sprintf("EARFCN ~%d (%s)", t{i, 3} + round((fc - t{i, 2})/100e3), t{i, 1}); %#ok<AGROW>
            end
        end
        ch = strjoin(c, " / ");
    case "Wi-Fi"
        ch = wifiChannel(fc, bw);
    case "BLE"
        k = round((fc - 2402e6)/2e6);
        if k >= 0 && k <= 39
            idx = [37, 0:10, 38, 11:36, 39];
            ch = sprintf("BLE ch %d", idx(k + 1));
            if ismember(idx(k + 1), 37:39), ch = ch + " (advertising)"; end
        end
    case "Bluetooth"
        k = round((fc - 2402e6)/1e6);
        if k >= 0 && k <= 78, ch = sprintf("BT ch %d", k); end
    case "GSM"
        if fc >= 935e6 && fc <= 960e6,           ch = sprintf("ARFCN %d (DL)", round((fc - 935e6)/200e3));
        elseif fc >= 925e6 && fc < 935e6,        ch = sprintf("ARFCN %d (DL)", 1024 + round((fc - 935e6)/200e3));
        elseif fc >= 890e6 && fc <= 915e6,       ch = sprintf("ARFCN %d (UL)", round((fc - 890e6)/200e3));
        elseif fc >= 1805e6 && fc <= 1880e6,     ch = sprintf("ARFCN %d (DL)", 512 + round((fc - 1805.2e6)/200e3));
        elseif fc >= 1710e6 && fc <= 1785e6,     ch = sprintf("ARFCN %d (UL)", 512 + round((fc - 1710.2e6)/200e3));
        end
    case "FM"
        if fc >= 87.4e6 && fc <= 108.1e6, ch = sprintf("FM %.1f MHz", round(fc/100e3)/10); end
    case "DVB-T2"
        k = 21 + round((fc - 474e6)/8e6);
        if k >= 21 && k <= 48, ch = sprintf("UHF ch %d (%g MHz)", k, 474 + 8*(k - 21)); end
    case "GPS"
        if abs(fc - 1575.42e6) < 1e6, ch = "L1 C/A"; end
end
end

function ch = wifiChannel(fc, bw)
ch = "";
if fc >= 2400e6 && fc <= 2500e6, band = 2.4; elseif fc >= 5150e6 && fc <= 5895e6, band = 5;
elseif fc >= 5925e6 && fc <= 7125e6, band = 6; else, return
end
w = 20 * 2^max(0, round(log2(max(bw, 15e6)/18e6)));    % 20 / 40 / 80 MHz channel
switch band
    case 2.4, if w == 20, list = 1:13; else, list = 3:11; end
    case 5
        if w == 20, list = [36:4:64 100:4:144 149:4:177];
        elseif w == 40, list = [38:8:62 102:8:142 151 159 167 175];
        else, list = [42 58 106 122 138 155 171];
        end
    case 6
        if w == 20, list = 1:4:233; elseif w == 40, list = 3:8:227; else, list = 7:16:215; end
end
fch = wlanChannelFrequency(list, band);
[~, k] = min(abs(fch - fc));
ch = sprintf("Wi-Fi ch %d", list(k));
if w > 20, ch = ch + sprintf(" (%d MHz)", w); end
end

function C = findComb(segs, rbw, df)
% Equally spaced narrow lines (e.g. switching-converter harmonics).
C = struct('Spacing_Hz', NaN, 'Freq_Hz', zeros(1, 0));
narrow = [segs.obw] < 3*max(rbw, df);
fn = sort([segs(narrow).fc]);
if numel(fn) < 4, return, end
d = diff(fn); sp = median(d);
tol = max(2*df, 0.03*sp);
good = abs(d - sp) < tol;
run = 0; best = 0; stop = 0;
for i = 1:numel(good)
    if good(i), run = run + 1; else, run = 0; end
    if run > best, best = run; stop = i; end
end
if best >= 3
    C.Spacing_Hz = sp; C.Freq_Hz = fn(stop - best + 1:stop + 1);
end
end

%% ===== output ==============================================================

function printResult(T, info)
if isempty(T)
    fprintf(['No signals more than the threshold above the noise floor (median %.1f dBm, RBW %s).\n' ...
        'A signal covering more than ~60 %% of the span hides the floor: widen the span, ' ...
        'or give ''NoiseFloor''.\n'], median(info.Floor), fmtHz(info.RBW));
    return
end
S = T; S.Freq_MHz = S.Freq_Hz/1e6; S.OBW_kHz = S.OBW_Hz/1e3; S.BW26_kHz = S.BW26_Hz/1e3;
disp(S(:, {'Freq_MHz', 'Level_dBm', 'OBW_kHz', 'BW26_kHz', 'Band', 'Channel', 'Technology', 'Score', 'Top3'}));
fprintf('Noise floor (median) %.1f dBm, RBW %s, detector %s, method %s.\n', median(info.Floor), ...
    fmtHz(info.RBW), info.Detector, info.MethodDetail);
if ~isnan(info.Comb.Spacing_Hz)
    fprintf(['%d equally spaced narrow lines, spacing %s: typical of a switching converter ' ...
        '(e.g. the inverter) or its harmonics.\n'], numel(info.Comb.Freq_Hz), fmtHz(info.Comb.Spacing_Hz));
end
if any(T.Technology == "Unresolved")
    fprintf(['Unresolved = too few trace points across the signal, or it runs into the edge of ' ...
        'the span: centre a smaller span on it or use more points.\n']);
end
if any(T.Technology == "Unknown")
    fprintf('Unknown = the shape fits none of the reference spectra (Top3 still lists the closest).\n');
end
end

function warnOnce(id, msg)
persistent shown
if isempty(shown), shown = strings(0); end
if ~any(shown == id), warning(id, '%s', msg); shown(end + 1) = id; end
end

function fig = plotResult(f, P, nf, segs, T, rbw, vis)
fig = figure('Name', 'ff_signal_id', 'Color', 'w', 'Visible', vis, 'Position', [60 60 1250 560]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
[sc, un] = freqUnit(max(abs(f)));
plot(ax, f/sc, P, 'Color', [0.2 0.4 0.9]);
plot(ax, f/sc, nf, '--', 'Color', [0 0.6 0]);
df = (f(end) - f(1)) / (numel(f) - 1);
for i = 1:height(T)
    lim = [segs(i).flo, max(segs(i).fhi, segs(i).flo + df)] / sc;
    xregion(ax, lim(1), lim(2), 'FaceColor', [1 0.75 0.2], 'FaceAlpha', 0.2);
    lbl = T.Technology(i);
    if ismember(lbl, ["Unresolved" "Unknown"]) && T.Top3(i) ~= ""
        lbl = lbl + ": " + extractBefore(T.Top3(i) + ";", ";") + "?";      % closest class
    elseif ~isnan(T.Score(i)) && lbl ~= ""
        lbl = sprintf('%s (%.2f)', lbl, T.Score(i));
    end
    text(ax, T.Freq_Hz(i)/sc, T.Level_dBm(i) + 1.5, lbl, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 8, 'Rotation', 90 * (height(T) > 8));
end
xlim(ax, f([1 end])/sc);
yl = ylim(ax); ylim(ax, [yl(1), yl(2) + 0.15*diff(yl)]);
xlabel(ax, "Frequency (" + un + ")"); ylabel(ax, 'Power (dBm)');
title(ax, sprintf('%d signal(s)  |  RBW %s  |  shaded = 99%% occupied bandwidth', height(T), fmtHz(rbw)));
legend(ax, {'Trace', 'Noise floor'}, 'Location', 'northeast', 'AutoUpdate', 'off');
end

function [scale, unit] = freqUnit(fmax)
if fmax >= 1e9,     scale = 1e9; unit = "GHz";
elseif fmax >= 1e6, scale = 1e6; unit = "MHz";
elseif fmax >= 1e3, scale = 1e3; unit = "kHz";
else,               scale = 1;   unit = "Hz";
end
end

function s = fmtHz(x)
if ~isfinite(x), s = '-'; return, end
[scale, unit] = freqUnit(abs(x));
s = sprintf('%.7g %s', x/scale, unit);
end
