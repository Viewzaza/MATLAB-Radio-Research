function M = ff_signal_classifier(varargin)
%FF_SIGNAL_CLASSIFIER  Train or load the signal classifiers used by ff_signal_id.
%
%   M = ff_signal_classifier                   load the cached models (train if missing, ~4 min)
%   M = ff_signal_classifier('Retrain', true)  train again
%   [f, P, rbw, truth] = M.synth("LTE 10 MHz") a synthetic FieldFox trace of one class
%   M.Accuracy                                 held-out accuracy of each method
%   M = ff_signal_classifier('Train', false)   never train (no cache: template method only)
%
% There is no I/Q on this FieldFox, so the classifiers learn from POWER
% SPECTRA. The training traces are synthetic: the reference spectra of
% ff_signal_templates passed through a model of how the analyzer measures:
%   - Gaussian RBW filter (3 dB width = RBW, 10 Hz ... 5 MHz)
%   - random point spacing (5 ... 150 points across the signal), span 1.4-3.5x the
%     visible width of the signal, random centre offset
%   - peak signal-to-noise ratio 6 ... 45 dB
%   - noise floor with ripple and tilt; signal and noise fluctuate like noise
%     (exponential power, smoothed by the VBW)
%   - detector: peak (highest value in each point's bucket), average or sample
%   - trace: clear-write, averaging (log, or power for the average detector) or max hold
%   - Wi-Fi / Bluetooth / BLE are bursty: present only part of the time
% Every synthetic trace goes through the same detection as ff_signal_id, then
%   "template"  fit error (mean |dB| difference) of each reference spectrum seen
%               through the same RBW and detector -> softmax (temperature fitted)
%   "ml"        bagged decision trees, fitcensemble (Statistics and ML Toolbox), on
%               41 features: log bandwidths (99 %, -3/-10/-20/-26 dB) and ratios,
%               flatness, symmetry, SNR, RBW, point spacing, detector and the
%               template-fit error of every class
%   "cnn"       1-D convolutional network (Deep Learning Toolbox: dlnetwork layers,
%               trainnet on CPU) on the shape resampled to 128 points, plus the
%               signal width and RBW as soft one-hot channels
%   ensemble    weighted geometric mean; weights chosen on a validation part
% Split: 70 % training, 15 % validation (temperature, weights), 15 % test (the
% accuracies and the confusion charts, ff_signal_confusion.png). "Unknown"
% limit per class: 1.25 x the 99th percentile of the training traces' own fit.
%
% Synthetic accuracy is an upper bound: real signals (partly loaded cells,
% multipath, carriers side by side) are harder, and LTE/NR of the same
% bandwidth are nearly indistinguishable from the spectrum alone.
%
% Name-Value options (defaults):
%   'Retrain'          false
%   'Train'            true    false: with no cache, return only the template method
%   'SamplesPerClass'  300     synthetic traces per class
%   'Epochs'           15      CNN training epochs
%   'File'             ff_signal_models.mat next to this file
%   'TemplateFile'     ""      templates MAT (default: ff_signal_templates' own)
%   'OutDir'           folder of 'File' (ff_signal_confusion.png)
%   'Visible'          'on'    show the confusion charts ('off': PNG only)
%   'Prompt'           true    accepted for suite compatibility (no dialogs)
%   'Seed'             7
%   'Verbose'          true
%
% Output M: Classes, ClassTech, ClassBW, Ensemble, CNN, Tau, FitLimit, Weights,
% FeatureNames, Accuracy, Confusion, Trained, Templates, and function handles
% predict (used by ff_signal_id) and synth (synthetic traces; options 'Fc',
% 'SNR', 'Detector', 'Points' across the signal, 'RBWRatio' = RBW/point spacing).
%
% See also ff_signal_id, ff_signal_templates, fitcensemble, trainnet.

opt = struct('Retrain', false, 'Train', true, 'SamplesPerClass', 300, 'Epochs', 15, 'File', "", ...
    'TemplateFile', "", 'OutDir', "", 'Visible', 'on', 'Prompt', true, 'Seed', 7, 'Verbose', true);
for k = 1:2:numel(varargin)
    if ~isfield(opt, varargin{k}), error('Unknown option "%s".', varargin{k}); end
    opt.(varargin{k}) = varargin{k + 1};
end
if strlength(string(opt.File)) == 0
    opt.File = fullfile(fileparts(mfilename('fullpath')), 'ff_signal_models.mat');
end
if strlength(string(opt.OutDir)) == 0, opt.OutDir = fileparts(opt.File); end

persistent cache cacheKey
targs = {'Verbose', true};                         % prints only when it has to build
if strlength(string(opt.TemplateFile)) > 0, targs = [targs, {'File', opt.TemplateFile}]; end
TPL = ff_signal_templates(targs{:});
key = string(opt.File) + "|" + string(TPL.Created);

if ~opt.Retrain && ~isempty(cache) && cacheKey == key && (cache.Trained || ~(opt.Train || isfile(opt.File)))
    M = cache;
    return
end
M = [];
if ~opt.Retrain && isfile(opt.File)
    S = load(opt.File, 'M');
    if isequal(S.M.Classes, TPL.Classes) && isfield(S.M, 'Weights'), M = S.M; end   % else: stale file
end
if isempty(M)
    if opt.Train || opt.Retrain
        M = trainModels(TPL, opt);
        save(opt.File, 'M', '-v7');
        if opt.Verbose, fprintf('Models saved to %s\n', opt.File); end
    else
        M = baseModel(TPL);
    end
end
M.Templates = TPL;
core = M;
M.predict = @(f, P, nf, rbw, segs, det, noise) predictSegments(core, f, P, nf, rbw, segs, det, noise);
M.synth = @(cls, varargin) synthClass(TPL, cls, varargin{:});
cache = M; cacheKey = key;
end

%% ===== prediction ==========================================================

function out = predictSegments(M, f, P, nf, rbw, segs, det, nz)
n = numel(segs); C = numel(M.Classes);
J = zeros(n, C); X = zeros(n, numel(M.FeatureNames)); Xc = zeros(128, cnnChannels, n);
for i = 1:n
    J(i, :) = templateFit(M.Templates, M.ClassIdx, f, P, nf, rbw, segs(i), det, nz);
    X(i, :) = features(f, P, nf, rbw, segs(i), det, nz, J(i, :));
    Xc(:, :, i) = cnnInput(f, P, nf, rbw, segs(i));
end
out.J = J;
out.FitAll = J / max(nz.mae, 0.5);                 % fit error in units of the noise scatter
out.Fit = min(out.FitAll, [], 2);
out.Template = softmaxRows(-(J - min(J, [], 2)) / M.Tau);
out.ML = []; out.CNN = [];
if ~isempty(M.Ensemble) && n > 0
    [~, sc] = predict(M.Ensemble, X);
    out.ML = reorder(sc, string(M.Ensemble.ClassNames), M.Classes);
end
if ~isempty(M.CNN) && n > 0
    sc = extractdata(predict(M.CNN, dlarray(Xc, "SCB"))).';  % much faster than minibatchpredict here
    out.CNN = reorder(double(sc), M.CNNClasses, M.Classes);
end
end

function p = softmaxRows(z)
p = exp(z - max(z, [], 2));
p = p ./ sum(p, 2);
end

function s = reorder(s, from, to)
[~, k] = ismember(to, from);
s = s(:, k);
end

function M = baseModel(TPL)
cls = TPL.Classes; T = TPL.Templates;
M.Classes = cls;
[~, M.ClassIdx] = ismember([T.Class], cls);
M.ClassTech = strings(size(cls)); M.ClassBW = zeros(size(cls));
for c = 1:numel(cls)
    k = [T.Class] == cls(c);
    M.ClassTech(c) = extractBefore(cls(c) + " ", " ");
    M.ClassBW(c) = max([T(k).OBW_Hz]);
end
M.ClassTech(M.ClassTech == "GPS") = "GNSS";       % band table name
M.FeatureNames = featureNames(cls);
M.Ensemble = []; M.CNN = []; M.CNNClasses = cls;
M.Tau = 0.5; M.FitLimit = 4; M.Weights = [1 1 1];
M.Accuracy = table(); M.Confusion = struct(); M.Trained = false;
M.TemplatesCreated = TPL.Created;
end

%% ===== template fit, features, CNN input ====================================

function Jc = templateFit(TPL, classIdx, f, P, nf, rbw, seg, det, nz)
% Mean |dB| error of the best fit of every template (gain, small centre shift).
N = numel(P); df = (f(end) - f(1)) / (N - 1);
w = seg.hi - seg.lo + 1;
pad = max([3, round(0.3*w), ceil(2*rbw/df)]);
idx = max(1, seg.lo - pad):min(N, seg.hi + pad);
idx = idx(1:ceil(numel(idx)/400):end);             % long traces: <= 400 points is plenty
x = f(idx) - seg.fc; y = P(idx);
nl = nz.meanRatio * 10.^(nf(idx)/10);
sl = max(10.^(y/10) - nl, 0);
[~, iv] = min(abs(log(TPL.RBWGrid) - log(rbw)));
u = subBins(df, rbw);
T = TPL.Templates;
Jt = inf(1, numel(T));
for k = 1:numel(T)
    vf = double(T(k).Views.f{iv}); vs = double(T(k).Views.s{iv});
    wk = max(T(k).OBW_Hz, rbw);
    scales = 1;
    if T(k).Technology == "FM", scales = [0.7 0.85 1 1.15]; end
    for sc = scales
        for d = [-0.06 -0.03 0 0.03 0.06] * wk
            s = viewAt(vf*sc, vs, x - d, u, det);
            G = (sl*s.') / (s*s.' + eps);
            if G <= 0, continue, end
            Jt(k) = min(Jt(k), mean(abs(y - 10*log10(G*s + nl))));
        end
    end
end
Jc = accumarray(classIdx(:), Jt(:), [], @min).';
end

function u = subBins(df, rbw)
% Offsets inside one trace point's bucket (detector acts over the bucket).
if df <= rbw/2, u = 0; return, end
Q = min(15, 2*ceil(df/rbw) + 1);
u = ((1:Q) - (Q + 1)/2) / Q * df;
end

function s = viewAt(vf, vs, x, u, det)
S = interp1(vf, vs, x(:) + u, 'linear', 0);
switch det
    case "average", s = mean(S, 2);
    case "sample",  s = S(:, (numel(u) + 1)/2);
    otherwise,      s = max(S, [], 2);
end
s = s.';
end

function names = featureNames(cls)
names = ["logOBW" "logBW26" "logBW3" "logBW10" "logRBW" "logOBWperRBW" "BW3perOBW" ...
    "BW10perOBW" "BW20perOBW" "BW26perOBW" "SNR" "Flatness" "CentreFraction" "Asymmetry" ...
    "logSpacingPerRBW" "logPoints" "Detector" "BestFit" "BestFitPerNoise", ...
    "fit_" + matlab.lang.makeValidName(cls)];
end

function x = features(f, P, nf, rbw, seg, det, nz, Jc)
N = numel(P); df = (f(end) - f(1)) / (N - 1);
lims = [f(seg.a) f(seg.b)];
sl = max(10.^(P/10) - nz.meanRatio * 10.^(nf/10), 0);
ss = movmean(sl, 3);
b = nan(1, 3); r = [3 10 20];
for i = 1:3
    try
        b(i) = powerbw(ss(:), f(:), max(rbw, max(diff(f))), lims, r(i));
    catch
    end
end
o = max(seg.obw, df); b26 = max(seg.bw26, df); b = max(b, df);
core = f >= seg.fc - 0.3*o & f <= seg.fc + 0.3*o;
flat = NaN;
if nnz(core) >= 3, flat = std(P(core)); end
in = f >= seg.flo & f <= seg.fhi;
pin = sl(in); fin = f(in);
if sum(pin) > 0
    cen = sum(pin(abs(fin - seg.fc) <= o/4)) / sum(pin);
    asym = (sum(pin.*fin)/sum(pin) - seg.fc) / o;
else
    cen = NaN; asym = NaN;
end
dcode = find(det == ["peak" "average" "sample"], 1);
Jmin = min(Jc);
x = [log10(o), log10(b26), log10(b(1)), log10(b(2)), log10(rbw), log10(o/rbw), b(1)/o, b(2)/o, ...
    b(3)/o, b26/o, seg.snr, flat, cen, asym, log10(df/rbw), log10(max(seg.npts, 0.1)), dcode, Jmin, ...
    Jmin/max(nz.mae, 0.5), Jc - Jmin];
end

function X = cnnInput(f, P, nf, rbw, seg)
% 128 x cnnChannels: the shape (dB above floor) over +-1.25x the signal
% width, normalised to its peak and also as dB/40; then the window width and
% RBW/width as soft one-hot codes (constant along the axis), so that the
% network can tell e.g. LTE 10 MHz from LTE 15 MHz (same shape, other width).
L = 128; df = (f(end) - f(1)) / (numel(f) - 1);
W = 1.25 * max([seg.bw26, seg.obw, 3*rbw]);
y = P - nf;
k = round(2*W/df/L);
if k > 1, y = 10*log10(movmean(10.^(y/10), k)); end
v = interp1(f, y, seg.fc + linspace(-W, W, L), 'linear', 0);
v = max(v(:), 0);
cw = exp(-(log10(2*W) - (3:0.125:8.5)).^2 / (2*0.1^2));    % 45 width codes
cr = exp(-(log10(rbw/(2*W)) - (-4:0.5:0)).^2 / (2*0.25^2));  % 9 RBW codes
X = [v / max(max(v), 1), min(v/40, 1.5), repmat([cw cr], L, 1)];
end

function n = cnnChannels
n = 2 + 45 + 9;
end

%% ===== synthetic FieldFox traces ===========================================

function [f, P, rbw, tr] = synthClass(TPL, cls, varargin)
% Synthetic trace of class (or template name) cls. Options: 'Fc', 'SNR',
% 'Detector', 'Points' (across the signal), 'RBWRatio' (RBW / point spacing).
o = struct('Fc', [], 'SNR', [], 'Detector', "", 'Points', [], 'RBWRatio', []);
for k = 1:2:numel(varargin), o.(varargin{k}) = varargin{k + 1}; end
T = TPL.Templates;
ks = find([T.Class] == string(cls) | [T.Name] == string(cls));
if isempty(ks), error('No template or class "%s".', cls); end
[f, P, rbw, tr] = synthTrace(TPL, ks(randi(numel(ks))), o);
end

function [f, P, rbw, tr] = synthTrace(TPL, k, o)
t = TPL.Templates(k); g = TPL.RBWGrid;
gOk = find(g <= 5e6);                                     % FieldFox max RBW 5 MHz
W = t.OBW_Hz;
if W == 0                                                 % CW: shape = RBW filter
    rbw = g(gOk(randi(numel(gOk))));
    ratio = pick(o.RBWRatio, 10^(log10(2) + rand*log10(6)));
    df = rbw / ratio;
    Wvis = 2.2 * rbw;
else
    n = pick(o.Points, 10^(log10(5) + rand*log10(30)));    % 5...150 points across
    df = min(W / n, 16e6);
    want = min(pick(o.RBWRatio, 10^(-1 + 1.6*rand)) * df, W/3);
    cand = gOk(g(gOk) <= max(W/3, g(1)));
    [~, i] = min(abs(log(g(cand)) - log(want)));
    rbw = g(cand(i));
    Wvis = W + 2.2*rbw;
end
iv = find(g == rbw, 1);
snr = pick(o.SNR, 6 + 39*rand);
vis = t.Views.f{iv}(t.Views.s{iv} > 10^(-snr/10));      % part above the noise floor
Wvis = max(Wvis, double(vis(end) - vis(1)));            % the user sets the span to see it all
Nw = min(max(round((Wvis*(1.4 + 2.1*rand) + 8*rbw) / df), 41), 1201);
span = (Nw - 1) * df;
fc = pick(o.Fc, 2e9 + 1e6*rand);
fc = max(fc, span + 1e6);
off = (2*rand - 1) * 0.25 * max(0, span - Wvis);
f = fc + off + ((0:Nw - 1) - (Nw - 1)/2) * df;

det = string(o.Detector);
if det == ""
    r = rand; det = "peak";
    if r > 0.45, det = "average"; end
    if r > 0.8, det = "sample"; end
end
r = rand;
if r < 0.5,     mode = "clear"; nSw = 1;
elseif r < 0.8, mode = "average"; nSw = randi([3 20]);
else,           mode = "maxhold"; nSw = randi([3 30]);
end
vbw = [1 1 3 10]; vbw = vbw(randi(4));                  % RBW/VBW

kap = 1;
if t.Technology == "FM", kap = 0.6 + 0.55*rand; end        % programme-dependent width
S = interp1(double(t.Views.f{iv})*kap, double(t.Views.s{iv}), f(:) - fc + subBins(df, rbw), 'linear', 0);
switch det
    case "average", sd = mean(S, 2);
    case "sample",  sd = S(:, (size(S, 2) + 1)/2);
    otherwise,      sd = max(S, [], 2);
end
A = 10^(snr/10);
x01 = linspace(0, 1, Nw).';
rip = rand*sin(2*pi*(0.2 + 1.8*rand)*x01 + 2*pi*rand) + (3*rand - 1.5)*(x01 - 0.5);
n0 = 10.^(rip/10);
cells = min(60, max(1, round(df/rbw)));                   % independent RBW cells per point
bursty = ismember(t.Technology, ["Wi-Fi" "BLE" "Bluetooth"]);
duty = 1; if bursty, duty = 0.15 + 0.85*rand; end
p = zeros(nSw, Nw);
for s = 1:nSw
    X = fluctuation(det, cells, vbw, Nw);
    on = 1;
    if bursty
        if det == "average", on = duty;
        elseif det == "sample", on = rand(Nw, 1) < duty;
        else, on = rand(Nw, 1) < 1 - (1 - duty)^cells;
        end
    end
    if W == 0, p(s, :) = A*sd + n0.*X;                     % deterministic carrier
    else,      p(s, :) = (A*sd.*on + n0) .* X;             % noise-like signal
    end
end
switch mode
    case "average"
        if det == "average", p = mean(p, 1); else, p = 10.^(mean(log10(p), 1)); end
    case "maxhold", p = max(p, [], 1);
end
L0 = -160 + 10*log10(rbw) + 15*rand;                       % DANL-like floor (dBm)
P = 10*log10(p) + L0;
tr = struct('Class', t.Class, 'Template', t.Name, 'Fc', fc, 'OBW', W, 'RBW', rbw, 'df', df, ...
    'Detector', det, 'SNR', snr, 'Mode', mode, 'Sweeps', nSw);
end

function v = pick(given, default)
if isempty(given), v = default; else, v = given; end
end

function X = fluctuation(det, cells, vbw, n)
% Power of noise-like signals relative to their mean: exponential, smoothed
% by the VBW (gamma), then max (peak) or mean (average) over the bucket.
switch det
    case "average", X = gamrnd(vbw*cells, 1/(vbw*cells), n, 1);
    case "sample",  X = gamrnd(vbw, 1/vbw, n, 1);
    otherwise,      X = max(gamrnd(vbw, 1/vbw, n, cells), [], 2);
end
end

%% ===== training ===========================================================

function M = trainModels(TPL, opt)
rng(opt.Seed);
M = baseModel(TPL);
C = numel(M.Classes); nPer = opt.SamplesPerClass; nMax = C*nPer;
nF = numel(M.FeatureNames);
X = zeros(nMax, nF); Xc = zeros(128, cnnChannels, nMax); J = zeros(nMax, C); fitOwn = zeros(nMax, 1);
y = strings(nMax, 1); n = 0;
classOf = [TPL.Templates.Class];
tAll = tic;
if opt.Verbose, fprintf('Generating %d synthetic FieldFox traces per class ...\n', nPer); end
for c = 1:C
    ks = find(classOf == M.Classes(c));
    got = 0; tries = 0;
    while got < nPer && tries < 3*nPer
        tries = tries + 1;
        [f, P, rbw, tr] = synthTrace(TPL, ks(randi(numel(ks))), struct('Fc', [], 'SNR', [], ...
            'Detector', "", 'Points', [], 'RBWRatio', []));
        [~, info] = ff_signal_id(f, P, rbw, 'Method', "none", 'Detector', tr.Detector);
        seg = pickSegment(info.Segments, tr);
        if isempty(seg) || seg.npts < 5, continue, end
        n = n + 1; got = got + 1;
        J(n, :) = templateFit(TPL, M.ClassIdx, f, P, info.Floor, rbw, seg, tr.Detector, info.Noise);
        X(n, :) = features(f, P, info.Floor, rbw, seg, tr.Detector, info.Noise, J(n, :));
        Xc(:, :, n) = cnnInput(f, P, info.Floor, rbw, seg);
        fitOwn(n) = J(n, c) / max(info.Noise.mae, 0.5);    % fit to its own class
        y(n) = M.Classes(c);
    end
    if opt.Verbose, fprintf('  %-18s %3d traces (%d tries)  %.0f s\n', M.Classes(c), got, tries, toc(tAll)); end
end
X = X(1:n, :); Xc = Xc(:, :, 1:n); J = J(1:n, :); fitOwn = fitOwn(1:n); y = categorical(y(1:n), M.Classes);
% 70 % training, 15 % validation (ensemble weights), 15 % test (reported)
cv = cvpartition(y, 'HoldOut', 0.3);
tr = training(cv); rest = find(test(cv));
cv2 = cvpartition(y(rest), 'HoldOut', 0.5);
va = false(n, 1); va(rest(training(cv2))) = true;
te = false(n, 1); te(rest(test(cv2))) = true;

% Template method: temperature with the best log-likelihood on the training part
D = J - min(J, [], 2);
yi = double(y);
taus = [0.05 0.1 0.2 0.3 0.5 0.75 1 1.5 2 3];
ll = zeros(size(taus));
for i = 1:numel(taus)
    lp = log(softmaxRows(-D(tr, :)/taus(i)) + 1e-12);
    ll(i) = mean(lp(sub2ind(size(lp), (1:nnz(tr)).', yi(tr))));
end
[~, i] = max(ll); M.Tau = taus(i);
% "Unknown" if a trace fits the template of its guessed class worse than
% 99 % of the genuine training traces of that class (+25 %)
M.FitLimit = zeros(1, C);
for c = 1:C, M.FitLimit(c) = max(2, 1.25*prctile(fitOwn(tr & y == M.Classes(c)), 99)); end
pT = softmaxRows(-D/M.Tau);

% (1) Bagged trees on the engineered features
if opt.Verbose, fprintf('Training bagged trees (fitcensemble) ...\n'); end
t0 = tic;
lrn = templateTree('MaxNumSplits', 150, 'MinLeafSize', 3, 'NumVariablesToSample', 8, 'Reproducible', true);
ens = fitcensemble(X(tr, :), y(tr), 'Method', 'Bag', 'NumLearningCycles', 60, 'Learners', lrn, ...
    'PredictorNames', cellstr(M.FeatureNames));
M.Ensemble = compact(ens);
[~, sc] = predict(M.Ensemble, X);
pML = reorder(sc, string(M.Ensemble.ClassNames), M.Classes);
if opt.Verbose, fprintf('  done in %.0f s\n', toc(t0)); end

% (2) 1-D CNN on the shape
if opt.Verbose, fprintf('Training 1-D CNN (trainnet, CPU) ...\n'); end
t0 = tic;
% Input "SCB" (128 points x channels x batch): the shape is treated as a
% 1-D spatial signal, and the fully connected layer sees every position, so
% WHERE the edges are (relative to the measured width) counts, not only
% which local patterns occur.
layers = [
    inputLayer([128 cnnChannels NaN], "SCB")
    convolution1dLayer(7, 16, 'Padding', 'same')
    batchNormalizationLayer
    reluLayer
    maxPooling1dLayer(2, 'Stride', 2)
    convolution1dLayer(5, 32, 'Padding', 'same')
    batchNormalizationLayer
    reluLayer
    maxPooling1dLayer(2, 'Stride', 2)
    convolution1dLayer(3, 32, 'Padding', 'same')
    batchNormalizationLayer
    reluLayer
    maxPooling1dLayer(2, 'Stride', 2)
    fullyConnectedLayer(64)
    reluLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(C)
    softmaxLayer];
topt = trainingOptions('adam', 'MaxEpochs', opt.Epochs, 'MiniBatchSize', 64, 'InitialLearnRate', 3e-3, ...
    'LearnRateSchedule', 'piecewise', 'LearnRateDropPeriod', max(1, round(opt.Epochs*0.7)), ...
    'LearnRateDropFactor', 0.2, 'Shuffle', 'every-epoch', 'InputDataFormats', "SCB", ...
    'Verbose', false, 'Plots', 'none', 'ExecutionEnvironment', 'cpu');
M.CNN = trainnet(Xc(:, :, tr), y(tr), layers, "crossentropy", topt);
M.CNNClasses = string(categories(y)).';
pCNN = reorder(double(minibatchpredict(M.CNN, Xc, 'InputDataFormats', "SCB")), M.CNNClasses, M.Classes);
if opt.Verbose, fprintf('  done in %.0f s\n', toc(t0)); end

% Ensemble = weighted geometric mean, weights with the best log-likelihood
% on the validation part (the models have not seen it)
L3 = cat(3, log(max(pT, 1e-3)), log(max(pML, 1e-3)), log(max(pCNN, 1e-3)));
[a, b, c] = ndgrid(0:0.25:1.5); Wg = [a(:) b(:) c(:)]; Wg(sum(Wg, 2) == 0, :) = [];
best = -inf;
for i = 1:size(Wg, 1)
    lp = softmaxRows(sum(L3(va, :, :) .* reshape(Wg(i, :), 1, 1, 3), 3));
    v = mean(log(lp(sub2ind(size(lp), (1:nnz(va)).', yi(va))) + 1e-12));
    if v > best, best = v; M.Weights = Wg(i, :); end
end
pAll = softmaxRows(sum(L3 .* reshape(M.Weights, 1, 1, 3), 3));
P = {pT, pML, pCNN, pAll};
names = ["template" "ml" "cnn" "ensemble (no band prior)"];
acc = zeros(4, 1);
for i = 1:4
    [~, k] = max(P{i}(te, :), [], 2);
    pred = categorical(M.Classes(k).', M.Classes);
    acc(i) = mean(pred == y(te));
    M.Confusion.(matlab.lang.makeValidName(extractBefore(names(i) + " ", " "))) = ...
        confusionmat(y(te), pred, 'Order', categorical(M.Classes, M.Classes));
end
M.Accuracy = table(names.', acc, repmat(nnz(te), 4, 1), 'VariableNames', {'Method', 'HeldOutAccuracy', 'NTest'});
M.Accuracy.Properties.Description = sprintf('Ensemble weights template/ml/cnn = %g/%g/%g', M.Weights);
M.Trained = true;
M.TrainedOn = datetime('now'); M.NTrain = nnz(tr);
if opt.Verbose
    fprintf(['Held-out accuracy on %d synthetic traces (template temperature %.2f, unknown if fit > %.1f...%.1f, ' ...
        'ensemble weights %g/%g/%g):\n'], nnz(te), M.Tau, min(M.FitLimit), max(M.FitLimit), M.Weights);
    disp(M.Accuracy);
    fprintf('Training took %.0f s.\n', toc(tAll));
end
saveConfusion(M, P, names, y, te, opt);
end

function seg = pickSegment(segs, tr)
seg = [];
if isempty(segs), return, end
fc = [segs.fc];
in = [segs.flo] <= tr.Fc & [segs.fhi] >= tr.Fc;
if any(in)
    k = find(in); [~, j] = max([segs(k).level]); seg = segs(k(j));
else
    [d, j] = min(abs(fc - tr.Fc));
    if d < max(tr.OBW, 3*tr.RBW)/2, seg = segs(j); end
end
end

function saveConfusion(M, P, names, y, te, opt)
fig = figure('Visible', opt.Visible, 'Color', 'w', 'Position', [20 20 1800 1000], ...
    'Name', 'ff_signal_classifier confusion charts');
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');
for i = 1:4
    [~, k] = max(P{i}(te, :), [], 2);
    pred = categorical(M.Classes(k).', M.Classes);
    set(groot, 'CurrentFigure', fig);
    nexttile(tl, i);               % the chart replaces this tile's axes
    cc = confusionchart(y(te), pred, 'RowSummary', 'row-normalized', 'FontSize', 6);
    cc.Title = sprintf('%s: %.1f %% held-out', names(i), 100*mean(pred == y(te)));
end
title(tl, 'ff\_signal\_classifier on synthetic FieldFox traces (rows = true class)');
file = fullfile(opt.OutDir, 'ff_signal_confusion.png');
exportgraphics(fig, file, 'Resolution', 90);
if opt.Verbose, fprintf('Confusion charts: %s\n', file); end
if strcmpi(opt.Visible, 'off'), close(fig); end
end
