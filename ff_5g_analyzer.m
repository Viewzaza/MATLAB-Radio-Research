function R = ff_5g_analyzer(src, varargin)
%FF_5G_ANALYZER  5G NR / LTE carrier identification and analysis with a FieldFox (no I/Q).
%
%   R = ff_5g_analyzer                                  measure (FieldFox at 10.90.36.200)
%   R = ff_5g_analyzer(ip, Name, Value, ...)
%   R = ff_5g_analyzer([], 'Demo', true)                synthetic Thai NR/LTE scene, no instrument
%   R = ff_5g_analyzer("...\nr5g_run.mat")             re-analyse a saved run (new options allowed)
%   R = ff_5g_analyzer([], 'File', "rec.mat")           a shared-format recording (fieldfox_live,
%                                                         ff_spectrum_pro, ...): spectrum analysis;
%                                                         a zero-span recording: envelope analysis
%
% Example (antenna on Port 2, no dialogs, only the 2600 and 1800 MHz bands):
%   R = ff_5g_analyzer("10.90.36.200", 'Bands', ["n41" "B3"], 'Prompt', false, 'Visible', 'off');
%   R.carriers        % one row per carrier: operator guess, technology, BW, TDD pattern, ACLR, level
%   R also has: details (per-carrier structs: zero-span fold, template ranking, ACLR, reference fit),
%   templates, coexistence, conclusions, notes, files, data (raw run, re-usable with 'File')
%
% What the FieldFox does (spectrum analyzer, Port 2 input; it never transmits):
%   1. Band scans of the Thai mobile downlink and TDD bands: average detector, 'ScanSweeps'
%      sweeps power-averaged in MATLAB (the max over the sweeps is kept for bursty carriers)
%   2. For every carrier found (strongest 'MaxCarriers'):
%      - zero span inside the carrier (RBW = VBW = 5 MHz, average detector, 'ZeroSpanTime'
%        20 ms, 'ZeroSpanSweeps' sweeps, free run or video / external trigger)
%      - a 45 ms peak-detector zero span for SSB bursts when the carrier is bursty or neither
%        continuous nor TDD (an idle NR cell sends only SSBs, every 20 ms)
%      - a wide sweep fc +- 2.5 channel bandwidths for ACLR, and the FieldFox's own channel
%        power (CHP, CALC:MEAS:DATA?) over the occupied bandwidth as a cross-check
%   Safe start: preamp as asked ('Preamp'), attenuation coupled to 'RefLevel'; on ADC overload
%   (CALC:MEAS:WAOR?) the preamp goes off, then the attenuation goes up 10 dB (max 30 dB).
%   The settings are backed up first and restored at the end, also after an error or Ctrl+C.
%
% What MATLAB computes:
%   - noise floor (DANL model and empirical), carrier detection (threshold over the floor,
%     regions split at notches / level steps only when both parts fit 3GPP bandwidths),
%     -3 dB edge fit (lsqcurvefit: RBW-convolved step + floor, nlparci intervals), OBW 99 %
%     (obw), channel power, PSD and energy per resource element (EPRE, ~RSRP at full load)
%   - channel bandwidth snapped to the TS 38.104 / 36.104 transmission bandwidths
%     (N_RB x 12 x SCS; checked against nrDLCarrierConfig.FR1BandwidthTable), EARFCN,
%     NR-ARFCN and the nearest GSCN (SSB raster)
%   - FDD / TDD from the band and from the zero-span envelope: autocorrelation period, DL duty,
%     slot string, and a match (folded-envelope correlation) against TDD templates built from
%     5G Toolbox resource grids (nrWaveformGenerator: DDDSU, DDSUU, DDDSUDDSUU, DDDDDDDSUU at
%     30 kHz, DDDSU at 15 kHz) and the LTE Toolbox (lteDuplexingInfo: UL/DL configurations 0-6
%     with special subframe 'LTESSC'); templates with identical envelopes are grouped
%   - SSB burst detection (period, burst length) for idle NR cells
%   - ACLR at +-1 and +-2 channel bandwidths vs TS 38.104 Table 6.6.3.2-1 / TS 36.104
%     6.6.2.1 (45 dB), noise-corrected, with the noise-limited bound and occupied-adjacent flags
%   - reference spectra (nrWaveformGenerator / lteTestModelTool) through the same RBW filter:
%     overlay and residual for the NR and the LTE hypothesis
%   - operator guess from the Thai band plan after the 29 Jun 2025 NBTC auction
%   - optional coexistence check with an EMI recording: UE SINR loss and throughput impact
% Results go to 'OutDir': nr5g_report.pdf, figures (PNG), nr5g_carriers.csv, nr5g_run.mat
% (raw data; re-use it with 'File') and rec\rec_*.mat (shared format: f, P, t, settings, idn).
%
% Name-Value options (defaults):
%  Scan
%   'Bands'         "thai"  ["B28" "B8" "L1500" "B3" "B1" "n40" "n41"]; more keys: "B5" "n78"
%                           "B28U" "B8U" "B3U" "B1U" (uplinks), "all"; or N x 2 [start stop] Hz
%   'Carriers'      []      skip the search: N x 1 centres or N x 2 [centre channelBW] in Hz
%                           (channelBW 0 = unknown; give it for a correct ACLR span)
%   'ScanRBW'       100e3   Hz (1-3-10 steps); 'ScanSweeps' 5; 'ScanPoints' "auto" (bin <= RBW)
%   'Threshold'     6       dB above the noise floor (PSD smoothed over 500 kHz)
%   'MinBW'         1e6     Hz, narrower detections are ignored (LTE 1.4 MHz occupies 1.08 MHz)
%   'MaxCarriers'   10      strongest carriers analysed further
%   'Preamp'        true    preamplifier (off automatically on ADC overload)
%   'RefLevel'      -30     dBm; the attenuation is coupled to it and stays constant
%   'ExtGain'       0       dB between antenna and FieldFox (POW:RF:EXTG), negative = loss;
%                           an external LNA is a positive gain
%  Zero span
%   'ZeroSpan' true; 'ZeroSpanTime' 20e-3 s; 'ZeroSpanSweeps' 6; 'ZeroSpanPoints' 1001
%   'ZeroSpanDetector' "AVER"  ("SAMP" shows +-10 dB OFDM speckle, "POS" widens the S slot)
%   'ZeroSpanOffset' 0      Hz from the carrier centre (kept 2.6 MHz inside the carrier edges)
%   'SSBCheck' true; 'SSBTime' 45e-3 s; 'SSBSweeps' 2
%   'Trigger'       "free"  "video" (level 'TriggerLevel' dBm, [] = midway DL/UL of a free
%                           sweep) or "external" (GPS 1PPS at Ref/Trig In); falls back to free
%   'LTESSC'        7       LTE special subframe configuration used for the TDD templates
%  ACLR
%   'ACLR' true; 'ACLRSweeps' 5; 'ACLRPoints' 2001; 'NoiseCorrection' true
%   'CrossCheck'    true    FieldFox channel power (CHP) over the occupied bandwidth
%  Analysis
%   'Reference'     true    toolbox reference spectra (a few s per carrier type)
%   'Templates'     "auto"  5G/LTE Toolbox templates (cached in support\nr5g); "analytic" =
%                           slot maps from TS 38.213 / 36.211 tables only (no toolbox call)
%  Coexistence
%   'EMIFile'       ""      shared-format recording of the interferer (e.g. inverter ON,
%                           recorded where a phone would be); "demo" = synthetic inverter noise
%   'EMIRefFile'    ""      the same with the interferer OFF (subtracted); "" = analyzer floor
%   'UENoiseFigure' 7       dB;  'EMICoupling' 0  dB added to the EMI level
%  Output
%   'Report' "pdf" ("docx", "html", "none"); 'Notes' ""; 'OutDir' "" (nr5g_<time>); 'Save' true
%   'Visible' 'on'; 'Prompt' true; 'Demo' false; 'File' ""; 'Port' 5025
%
% Read this before trusting the results:
%   - NO I/Q CAPTURE (option 351 absent): nothing is demodulated. No PCI, MIB/SIB, PLMN, EVM or
%     beam ID. The operator is a guess from the frequency and the 2025 NBTC band plan (some
%     block edges are inferred); technology, bandwidth and pattern are hypotheses with a score.
%   - Some TDD patterns have identical envelopes: NR 30 kHz DDDDDDDSUU (S 6:4:4), LTE config 2
%     with SSC 7 and NR 15 kHz DDDSU (S 10:2:2) all give 3.714 ms on / 1.286 ms off per 5 ms.
%     They are reported as one group; the spectrum width (LTE 18.0 / NR 18.36 or 19.08 MHz for
%     20 MHz) then decides, weakly (edge fit +-0.1-0.2 MHz over the air).
%   - The pattern is read from what is transmitted: a lightly loaded NR cell leaves D slots
%     empty (NR has no always-on CRS), so the DL duty and the slot string follow the traffic.
%   - Zero span needs the carrier about 10 dB above the floor in 5 MHz: about -76 dBm in 5 MHz
%     with the preamp (RSRP better than about -100 dBm), -60 dBm without.
%   - An idle NR cell (SSB only, ~1 ms of bursts per 20 ms) shows in a swept scan only when the
%     sweep crosses it during a burst: it is often missed with 5 sweeps. Use 'ScanSweeps' 20+ (the
%     max over the sweeps is searched) or give it in 'Carriers' for the zero-span SSB check.
%   - 'Trigger' "external" needs a 1PPS on Ref/Trig In. Without one the sweep never completes; the
%     tool then frees the trigger (tested on a mock only: if this firmware blocks commands while a
%     *OPC? is pending, the recovery may need a power-off/on of the LAN session).
%   - ACLR of a distant base station over the air is NOT a conformance test: the analyzer floor
%     limits it (max = channel power - noise in the measurement bandwidth), the adjacent channel
%     is often another operator's carrier, multipath ripples the spectrum, and the -13 dBm/MHz
%     absolute limit cannot be judged without the path loss. The 45 dB verdict is INFO only;
%     for compliance measure conducted at about -20 dBm.
%   - Swept spectra of TDD carriers are ragged (the average detector dwells tens of us per
%     point); levels are time averages including UL gaps (DL-on = average - 10 log10(duty)).
%   - The FieldFox zero-span sample rate is not specified; the time axis uses SWE:TIME? (the
%     first-to-last point time).
%   - Coexistence is a first-order SINR estimate for a UE where the EMI was recorded (same
%     antenna), throughput by the TR 36.942 truncated Shannon bound; not a link simulation.
%     The FieldFox floor (-153 dBm/Hz with preamp) is ~14 dB above a UE noise floor (-174 + NF):
%     EMI that does not raise the analyzer floor by 0.5 dB is reported as an upper bound
%     (Resolved = false). Use an LNA ('ExtGain') or a closer antenna to resolve weak EMI.

if nargin < 1, src = "10.90.36.200"; end
opt = parseOptions(varargin{:});
if ~opt.Demo && strlength(opt.File) == 0 && (ischar(src) || isstring(src)) && isscalar(string(src)) ...
        && strlength(string(src)) > 0 && (endsWith(lower(string(src)), ".mat") || isfile(string(src)))
    opt.File = string(src);
end
if opt.Demo, mode = "demo"; elseif strlength(opt.File) > 0, mode = "file"; else, mode = "instrument"; end
if mode == "instrument" && (isempty(src) || strlength(string(src)) == 0)
    error('ff_5g_analyzer:input', 'Give the FieldFox IP address, ''Demo'', true or ''File'', "...mat".');
end
if strlength(opt.OutDir) == 0
    opt.OutDir = fullfile(pwd, "nr5g_" + string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
end
if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end

switch mode
    case "instrument"
        if ~safetyGate(opt)
            fprintf('Cancelled: nothing was sent to the instrument.\n');
            R = struct('aborted', true, 'carriers', table(), 'notes', "Cancelled at the start question.");
            return
        end
        D = measureRun(string(src), opt);
    case "demo"
        D = demoRun(opt);
    otherwise
        D = loadInput(opt.File, opt);
end

fprintf('Analysing ...\n');
A = analyzeRun(D, opt);
figs = makeFigures(D, A, opt);
R = struct('carriers', A.table, 'details', A.carriers, 'templates', templateTable(A.templates), ...
    'coexistence', A.coex, 'conclusions', A.conclusions, 'notes', [D.notes; A.notes], ...
    'source', D.source, 'idn', D.idn, 'aborted', D.aborted, 'files', figs, 'outDir', opt.OutDir, ...
    'data', D, 'opt', opt);
if opt.Save
    R.files = [R.files; saveRun(D, A, opt)];
end
if opt.Report ~= "none"
    try
        R.files = [R.files; makeReport(D, A, opt, figs)];
    catch err
        warning('ff_5g_analyzer:report', 'Report not written: %s', err.message);
    end
end
printSummary(R);
end

%% ===== options ==================================================================

function opt = parseOptions(varargin)
p = inputParser;
p.FunctionName = 'ff_5g_analyzer';
defs = {'Bands', "thai"; 'Carriers', []; 'ScanRBW', 100e3; 'ScanSweeps', 5; 'ScanPoints', "auto"; ...
    'Threshold', 6; 'MinBW', 1e6; 'MaxCarriers', 10; 'Preamp', true; 'RefLevel', -30; 'ExtGain', 0; ...
    'ZeroSpan', true; 'ZeroSpanTime', 20e-3; 'ZeroSpanSweeps', 6; 'ZeroSpanPoints', 1001; ...
    'ZeroSpanDetector', "AVER"; 'ZeroSpanOffset', 0; 'SSBCheck', true; 'SSBTime', 45e-3; 'SSBSweeps', 2; ...
    'Trigger', "free"; 'TriggerLevel', []; 'LTESSC', 7; 'ACLR', true; 'ACLRSweeps', 5; 'ACLRPoints', 2001; ...
    'NoiseCorrection', true; 'CrossCheck', true; 'Reference', true; 'Templates', "auto"; ...
    'EMIFile', ""; 'EMIRefFile', ""; 'UENoiseFigure', 7; 'EMICoupling', 0; 'Report', "pdf"; 'Notes', ""; ...
    'OutDir', ""; 'Save', true; 'Visible', 'on'; 'Prompt', true; 'Demo', false; 'File', ""; 'Port', 5025};
for k = 1:size(defs, 1), addParameter(p, defs{k, 1}, defs{k, 2}); end
parse(p, varargin{:});
opt = p.Results;
opt.bandsGiven = ~any(strcmp(p.UsingDefaults, 'Bands'));
for nm = ["Report" "ZeroSpanDetector" "Trigger" "Templates" "EMIFile" "EMIRefFile" "Notes" "OutDir" "File"]
    v = opt.(nm); if isempty(v), v = ""; end
    opt.(nm) = string(v);
end
if islogical(opt.Visible) || isnumeric(opt.Visible)
    if opt.Visible, opt.Visible = 'on'; else, opt.Visible = 'off'; end
end
opt.Visible = lower(char(opt.Visible));
mustBeOneOf = @(v, set, nm) assert(any(v == set), 'ff_5g_analyzer:option', '''%s'' must be one of: %s.', nm, strjoin(set, ', '));
opt.Report = lower(opt.Report); mustBeOneOf(opt.Report, ["pdf" "docx" "html" "none"], 'Report');
opt.ZeroSpanDetector = upper(opt.ZeroSpanDetector); mustBeOneOf(opt.ZeroSpanDetector, ["AVER" "SAMP" "POS"], 'ZeroSpanDetector');
opt.Trigger = lower(opt.Trigger); if opt.Trigger == "ext", opt.Trigger = "external"; end
mustBeOneOf(opt.Trigger, ["free" "video" "external"], 'Trigger');
opt.Templates = lower(opt.Templates); mustBeOneOf(opt.Templates, ["auto" "analytic"], 'Templates');
mustBeOneOf(string(opt.Visible), ["on" "off"], 'Visible');
for nm = ["ScanSweeps" "ZeroSpanSweeps" "SSBSweeps" "ACLRSweeps" "MaxCarriers"]
    opt.(nm) = max(1, round(double(opt.(nm))));
end
opt.ZeroSpanPoints = min(10001, max(101, round(opt.ZeroSpanPoints)));
opt.ACLRPoints = min(10001, max(201, round(opt.ACLRPoints)));
opt.ZeroSpanTime = min(2, max(2e-3, opt.ZeroSpanTime));
opt.SSBTime = min(2, max(10e-3, opt.SSBTime));
assert(isscalar(opt.LTESSC) && any(opt.LTESSC == 0:9), 'ff_5g_analyzer:option', '''LTESSC'' must be 0 ... 9.');
assert(opt.ScanRBW > 0 && opt.Threshold > 0 && opt.MinBW > 0, 'ff_5g_analyzer:option', 'ScanRBW, Threshold and MinBW must be > 0.');
if isstring(opt.ScanPoints) || ischar(opt.ScanPoints)
    assert(strcmpi(string(opt.ScanPoints), "auto"), 'ff_5g_analyzer:option', '''ScanPoints'' must be "auto" or 101 ... 10001.');
    opt.ScanPoints = "auto";
else
    opt.ScanPoints = min(10001, max(101, round(opt.ScanPoints)));
end
if ~isempty(opt.Carriers)
    c = double(opt.Carriers);
    assert(ismatrix(c) && any(size(c, 2) == [1 2]) && all(c(:, 1) > 1e5 & c(:, 1) < 6.5e9), ...
        'ff_5g_analyzer:option', '''Carriers'' must be N x 1 centres or N x 2 [centre channelBW] in Hz.');
    if size(c, 2) == 1, c(:, 2) = 0; end
    opt.Carriers = c;
end
if ischar(opt.Bands) || iscellstr(opt.Bands), opt.Bands = string(opt.Bands); end
if isnumeric(opt.Bands)
    assert(size(opt.Bands, 2) == 2 && all(opt.Bands(:, 2) > opt.Bands(:, 1)) && all(opt.Bands(:) >= 1e5 & opt.Bands(:) <= 6.5e9), ...
        'ff_5g_analyzer:option', 'Numeric ''Bands'' must be N x 2 [start stop] in Hz (100 kHz ... 6.5 GHz).');
end
scanBands(opt);                             % unknown band keys fail here, before anything is measured
end

function ok = safetyGate(opt)
msg = sprintf(['ff_5g_analyzer takes control of the FieldFox (all settings are backed up and restored). ' ...
    'Connect the antenna to Port 2 (SA input, damage above +27 dBm or 50 V DC: never connect a ' ...
    'transmitter or a powered LNA output without a DC block). ExtGain = %g dB, preamp %s. ' ...
    'The FieldFox only receives.'], opt.ExtGain, string(onoff(opt.Preamp)));
fprintf('\n%s\n', msg);
ok = true;
if ~opt.Prompt, return, end
if usejava('desktop')
    ok = strcmp(questdlg(msg, 'ff_5g_analyzer', 'Measure', 'Abort', 'Abort'), 'Measure');
else
    try
        ok = strcmpi(strtrim(input('Type "yes" to measure: ', 's')), 'yes');
    catch
        ok = false;
        fprintf('No console for the start question: run with ''Prompt'', false to measure headless.\n');
    end
end
end

%% ===== measurement flow (shared by the instrument and the demo) =================

function D = newRun(source, opt)
D = struct('kind', "ff_5g_analyzer", 'version', 1, 'source', source, 'idn', "", ...
    'options', strings(1, 0), 'time', datetime('now'), 'notes', strings(0, 1), 'extGain', opt.ExtGain, ...
    'scan', struct('key', {}, 'name', {}, 'acq', {}), ...
    'zs', struct('fc', {}, 'fz', {}, 'occ', {}, 'acq', {}, 'ssb', {}), ...
    'aclr', struct('fc', {}, 'chBW', {}, 'occ', {}, 'acq', {}, 'chp', {}, 'chpNote', {}), ...
    'planned', [], 'aborted', false, 'truth', [], 'opt', opt);
end

function D = runFlow(E, D, opt)
nfail = 0;
bands = scanBands(opt);
for b = 1:numel(bands)
    fprintf('Scan %d/%d: %s %s ...\n', b, numel(bands), bands(b).name, fmtRange([bands(b).start bands(b).stop]));
    try
        a = E.sweep(planScan(bands(b), opt), opt.ScanSweeps);
        D.scan(end + 1) = struct('key', bands(b).key, 'name', bands(b).name, 'acq', a);
        nfail = 0;
    catch err
        [D, nfail, stop] = stepFailed(D, err, "Scan " + bands(b).name, nfail);
        if stop, return, end
    end
    checkpoint(D, opt);
end
[C, G] = planCarriers(D, opt);
D.planned = C;
if isempty(C)
    D.notes(end + 1, 1) = "No carrier above the threshold in the scanned bands.";
    fprintf('   no carrier found above %g dB over the floor.\n', opt.Threshold);
end
for k = 1:numel(C)
    c = C(k);
    fprintf('Carrier %d/%d: %s, about %s wide%s ...\n', k, numel(C), fmtHz6(c.fc), fmtHz(c.occ), ...
        repmat(' (bursty)', 1, c.bursty));
    ssbOnly = c.bursty;                     % an idle cell (SSB bursts only) has no ACLR to measure
    if opt.ZeroSpan
        try
            z = E.zero(planZero(c, opt, "main"), opt.ZeroSpanSweeps);
            s = [];
            if opt.SSBCheck && needSSB(z, c, D, G)
                s = E.zero(planZero(c, opt, "ssb"), opt.SSBSweeps);
                ssbOnly = ssbOnly || analyzeSSB(s, D).detected;
            end
            D.zs(end + 1) = struct('fc', c.fc, 'fz', z.settings.center, 'occ', c.occ, 'acq', z, 'ssb', s);
            nfail = 0;
        catch err
            [D, nfail, stop] = stepFailed(D, err, "Zero span at " + fmtHz6(c.fc), nfail);
            if stop, return, end
        end
    end
    if opt.ACLR && ~ssbOnly
        try
            a = E.sweep(planACLR(c, opt), opt.ACLRSweeps);
            chp = [NaN NaN]; note = "";
            if opt.CrossCheck && isfield(E, 'chp') && ~isempty(E.chp)
                [chp, note] = E.chp(c, a);
            end
            D.aclr(end + 1) = struct('fc', c.fc, 'chBW', c.chBW, 'occ', c.occ, 'acq', a, 'chp', chp, 'chpNote', note);
            nfail = 0;
        catch err
            [D, nfail, stop] = stepFailed(D, err, "ACLR sweep at " + fmtHz6(c.fc), nfail);
            if stop, return, end
        end
    end
    checkpoint(D, opt);
end
end

function [D, nfail, stop] = stepFailed(D, err, what, nfail)
nfail = nfail + 1;
D.notes(end + 1, 1) = what + " failed: " + string(err.message);
warning('ff_5g_analyzer:step', '%s failed: %s', what, err.message);
stop = nfail >= 2;                          % two failures in a row: connection lost or unit busy
if stop
    D.aborted = true;
    D.notes(end + 1, 1) = "Measurement stopped after two failed steps: analysing what was measured.";
end
end

function checkpoint(D, opt)
% Raw data after every step, so an interrupted run can be re-analysed with 'File'.
if ~opt.Save, return, end
try
    save(fullfile(opt.OutDir, "nr5g_run.mat"), "D");
catch
end
end

function bands = scanBands(opt)
bands = struct('key', {}, 'name', {}, 'start', {}, 'stop', {});
if ~isempty(opt.Carriers) && ~opt.bandsGiven, return, end      % given carriers: no search
if isnumeric(opt.Bands)
    for k = 1:size(opt.Bands, 1)
        bands(end + 1) = struct('key', "user" + k, 'name', "User " + fmtRange(opt.Bands(k, :)), ...
            'start', opt.Bands(k, 1), 'stop', opt.Bands(k, 2)); %#ok<AGROW>
    end
    return
end
B = bandTable();
keys = opt.Bands;
if any(lower(keys) == "thai"), keys = [keys(lower(keys) ~= "thai"), "B28" "B8" "L1500" "B3" "B1" "n40" "n41"]; end
if any(lower(keys) == "all"), keys = B.key.'; end
keys = unique(keys, 'stable');
for k = keys
    i = find(lower(B.key) == lower(k), 1);
    if isempty(i), error('ff_5g_analyzer:option', 'Unknown band key "%s". Known: %s.', k, strjoin(B.key, ', ')); end
    m = max(3e6, 0.03*(B.fhi(i) - B.flo(i)));
    bands(end + 1) = struct('key', B.key(i), 'name', B.name(i), 'start', max(1e5, B.flo(i) - m), ...
        'stop', min(6.5e9, B.fhi(i) + m)); %#ok<AGROW>
end
[~, o] = sort([bands.start]); bands = bands(o);
end

function cfg = planScan(band, opt)
span = band.stop - band.start;
rbw = pick13(opt.ScanRBW);
if isstring(opt.ScanPoints), pts = min(10001, max(1001, ceil(span/rbw) + 1)); else, pts = opt.ScanPoints; end
cfg = struct('label', "scan_" + band.key, 'start', band.start, 'stop', band.stop, 'rbw', rbw, ...
    'points', pts, 'det', "AVER", 'timeout', 120);
end

function cfg = planZero(c, opt, kind)
half = max(0, c.occ/2 - 2.6e6);
fz = min(max(c.fc + opt.ZeroSpanOffset, c.fc - half), c.fc + half);
cfg = struct('label', "zs", 'center', fz, 'rbw', 5e6, 'vbw', 5e6, 'points', opt.ZeroSpanPoints, ...
    'det', opt.ZeroSpanDetector, 'time', opt.ZeroSpanTime, 'trigger', opt.Trigger, ...
    'trigLevel', opt.TriggerLevel, 'timeout', opt.ZeroSpanTime + 15);
if kind == "ssb"
    cfg.label = "ssb"; cfg.center = c.fc; cfg.det = "POS"; cfg.time = opt.SSBTime;
    cfg.trigger = "free"; cfg.timeout = opt.SSBTime + 15;
end
end

function cfg = planACLR(c, opt)
bw = c.chBW*1e6;
if ~(bw > 0), bw = max(5e6, 1.15*c.occ); end
st = max(1e5, c.fc - 2.5*bw); sp = min(6.5e9, c.fc + 2.5*bw);
pts = opt.ACLRPoints;
rbw = pick13(min(max(c.occ, 1e6)/60, 2*(sp - st)/(pts - 1)));
rbw = min(3e6, max(1e3, rbw));
cfg = struct('label', "aclr", 'start', st, 'stop', sp, 'rbw', rbw, 'points', pts, 'det', "AVER", 'timeout', 180);
end

function yes = needSSB(z, c, D, G)
% SSB check only where it can show something: bursty carriers, or neither continuous nor TDD.
if c.bursty, yes = true; return, end
P = levels(z, D); N = size(P, 2); T = z.settings.sweepTime;
t = (0:N - 1).'*T/(N - 1);
r = tddAnalyze(t, P(1, :).', 0.5e-3);
fl5 = floorHzOf(z, D, G) + 10*log10(1.0645*z.settings.rbw);
isTDD = ~r.bursty && r.dyn >= 6 && r.periodicity >= 0.5 && r.onFrac >= 0.3;
isCont = r.dyn < 6 && r.hi > fl5 + 8;
yes = ~(isTDD || isCont);
end

function [C, G] = planCarriers(D, opt)
G = floorGroups(D);
if ~isempty(opt.Carriers)
    C = manualCarriers(opt.Carriers);
    return
end
C = emptyCarrier();
for i = 1:numel(D.scan)
    C = [C, detectInAcq(D.scan(i).acq, D, G, opt, "scan", i)]; %#ok<AGROW>
end
if isempty(C), return, end
[~, o] = sort([C.chP], 'descend');
C = C(o(1:min(end, opt.MaxCarriers)));
[~, o] = sort([C.fc]); C = C(o);
for k = 1:numel(C)
    s = snapCandidates(C(k).occ, C(k).band);
    s = s(s.tech ~= "NR-SSB", :);
    if ~isempty(s) && s.err(1) < 0.12, C(k).chBW = s.chBW(1); end
end
end

function C = manualCarriers(M)
C = emptyCarrier();
for k = 1:size(M, 1)
    c = newCarrier();
    c.src = "manual"; c.fc = M(k, 1); c.fcPlan = M(k, 1);
    [c.band, c.bandName, c.duplexBand] = bandOf(c.fc);
    c.chBW = M(k, 2)/1e6;
    if c.chBW > 0
        s = snapCandidates(0, c.band); s = s(s.tech ~= "NR-SSB" & abs(s.chBW - c.chBW) < 1e-6, :);
        if isempty(s), c.occ = 0.9*c.chBW*1e6; else, c.occ = min(s.occ); end
        c.occSrc = "given";
    else
        c.chBW = NaN; c.occ = 18e6; c.occSrc = "assumed";
        c.notes(end + 1) = "channel bandwidth not given: 20 MHz assumed for the measurement plan";
    end
    c.f1 = c.fc - c.occ/2; c.f2 = c.fc + c.occ/2;
    C(end + 1) = c; %#ok<AGROW>
end
end

%% ===== instrument engine =========================================================

function D = measureRun(ip, opt)
fprintf('Connecting to %s:%d ...\n', ip, opt.Port);
ff = FieldFox(ip, opt.Port);
fprintf('Connected: %s\n', ff.IDN);
if ~ff.hasOption("233"), error('ff_5g_analyzer:option', 'Needs the spectrum analyzer (option 233).'); end
ff.takeControl();                          % backs up every setting; release() restores them
cleanup = onCleanup(@() releaseInstrument(ff));
ff.setMode("SA");
D = newRun("instrument", opt);
D.idn = string(ff.IDN); D.options = ff.Options;
if opt.Preamp && ~ff.hasOption("235")
    opt.Preamp = false;
    D.notes(end + 1, 1) = "Preamp asked but option 235 is not installed: measured without it.";
end
[notes, extgInstr] = prepareSA(ff, opt);
D.notes = [D.notes; notes];
st = containers.Map('KeyType', 'char', 'ValueType', 'any');
st('preamp') = opt.Preamp; st('att') = "auto"; st('waor') = probe(ff, "CALC:MEAS:WAOR?");
if ~st('waor'), D.notes(end + 1, 1) = "CALC:MEAS:WAOR? not available: overload was not checked."; end
E.sweep = @(cfg, n) instSweep(ff, st, extgInstr, cfg, n);
E.zero = @(cfg, n) instZero(ff, st, extgInstr, cfg, n, D);
E.chp = @(c, a) instCHP(ff, c, a, extgInstr, D.extGain);
D = runFlow(E, D, opt);
try
    e = ff.errors();
    if ~isempty(e), D.notes(end + 1, 1) = "Instrument errors at the end: " + strjoin(e, "; "); end
catch
end
clear cleanup                              % restore the instrument now (prints the outcome)
end

function [notes, extgInstr] = prepareSA(ff, opt)
notes = strings(0, 1);
u = tryQuery(ff, "AMPL:UNIT?");            % display unit = TRAC:DATA? unit
if strlength(u) > 0 && ~startsWith(upper(u), "DBM")
    ff.write("AMPL:UNIT DBM");
    notes(end + 1, 1) = "Display unit was " + u + ": set to dBm.";
end
ff.write("SENS:MEAS:CHAN NONE");           % no channel measurement
ff.write("SENS:MEAS:INT NONE");            % no spectrogram
ff.write("TRIG:SOUR FREE");                % free run (zero span must not wait for a trigger)
ff.write(sprintf("SENS:POW:RF:GAIN %d", opt.Preamp));
ff.write("SENS:POW:RF:ATT:AUTO 1");        % attenuation coupled to the reference level
ff.write(sprintf("DISP:WIND:TRAC:Y:RLEV %g", opt.RefLevel));
ff.opc(30);
e = ff.errors();
if ~isempty(e), notes(end + 1, 1) = "Setup errors (ignored): " + strjoin(e, "; "); end
ff.write(sprintf("SENS:POW:RF:EXTG %.4f", opt.ExtGain));   % verified on this unit (negative = loss)
ff.opc(30);
e = ff.errors();
rb = str2double(tryQuery(ff, "SENS:POW:RF:EXTG?"));
ff.errors();
if isempty(e) && isfinite(rb) && abs(rb - opt.ExtGain) < 0.01
    extgInstr = opt.ExtGain;
else
    extgInstr = 0;
    if isfinite(rb), extgInstr = rb; end
    notes(end + 1, 1) = sprintf("POW:EXTG %g not confirmed (read back %g): ExtGain is applied in MATLAB.", opt.ExtGain, rb);
end
end

function a = instSweep(ff, st, extgInstr, cfg, n)
ff.write("SENS:SWE:TIME 0");               % PG: 0 = automatic sweep time after a zero span (UNVERIFIED)
ff.errors();
args = {'Start', cfg.start, 'Stop', cfg.stop, 'RBWAuto', false, 'RBW', cfg.rbw, 'VBW', cfg.rbw, ...
    'Points', cfg.points, 'Detector', cfg.det, 'TraceType', "CLRW", 'Average', 1, 'Preamp', st('preamp')};
att = st('att');
if isnumeric(att), args = [args, {'AttenAuto', false, 'Atten', att}]; else, args = [args, {'AttenAuto', true}]; end
ff.saConfigure(args{:});
s = ff.saSettings();
P = zeros(n, s.points); t = NaT(n, 1); over = false; k = 0;
while k < n
    k = k + 1;
    ff.sweep(cfg.timeout);
    y = ff.saTrace();
    if numel(y) ~= s.points
        error('ff_5g_analyzer:points', '%s: trace has %d points, expected %d.', cfg.label, numel(y), s.points);
    end
    if k == 1 && st('waor') && waor(ff)
        if st('preamp')                    % first remedy: preamp off, then attenuation
            fprintf('   ADC overload: preamp off\n');
            st('preamp') = false;
            ff.saConfigure('Preamp', false); s = ff.saSettings(); k = 0;
            continue
        elseif s.atten < 30
            na = min(30, 5*ceil((s.atten + 10)/5));
            fprintf('   ADC overload at %g dB attenuation: %g dB\n', s.atten, na);
            st('att') = na;
            ff.saConfigure('AttenAuto', false, 'Atten', na); s = ff.saSettings(); k = 0;
            continue
        end
        over = true;
        warning('ff_5g_analyzer:overload', '%s: ADC overload even with 30 dB attenuation: add an attenuator.', cfg.label);
    end
    P(k, :) = y; t(k) = datetime('now');
end
if st('waor') && ~over, over = waor(ff); end
a = makeAcq(cfg, s.f, P, t, s, extgInstr, over, false, NaN, "", "free");
fprintf('   %-14s %s, RBW %s, %d pts, %s, preamp %s, att %g dB%s\n', cfg.label, fmtRange([s.start s.stop]), ...
    fmtHz(s.rbw), s.points, s.detector, onoff(s.preamp), s.atten, repmat('  OVERLOAD', 1, over));
end

function a = instZero(ff, st, extgInstr, cfg, n, D)
% Zero span: FREQ:SPAN 0, SWE:TIME, RBW/VBW 5 MHz, SAMP/AVER (verified on this unit); POS and the
% trigger commands are from the programming guide (errors are checked, free run is the fallback).
notes = strings(0, 1);
ff.write("SENS:FREQ:SPAN 0");
ff.write(sprintf("SENS:FREQ:CENT %.10g", cfg.center));
ff.write(sprintf("SENS:BAND:RES %.10g", cfg.rbw));
ff.write(sprintf("SENS:BAND:VID %.10g", cfg.vbw));
ff.write(sprintf("SENS:SWE:POIN %d", cfg.points));
ff.write(sprintf("SENS:POW:RF:GAIN %d", st('preamp')));
att = st('att');
if isnumeric(att), ff.write(sprintf("SENS:POW:RF:ATT %g", att)); end
ff.opc(30);
e = ff.errors();
if ~isempty(e), notes(end + 1) = "setup: " + strjoin(e, "; "); end
det = cfg.det;
ff.write("SENS:DET:FUNC " + det);
ff.opc(30);
eD = ff.errors();
if ~isempty(eD) && det ~= "SAMP"
    notes(end + 1) = "DET " + det + " rejected in zero span (" + strjoin(eD, "; ") + "): SAMP used";
    det = "SAMP"; ff.write("SENS:DET:FUNC SAMP"); ff.opc(30); ff.errors();
end
ff.write(sprintf("SENS:SWE:TIME %.6g", cfg.time));
ff.opc(30);
eT = ff.errors();
if ~isempty(eT), notes(end + 1) = "SWE:TIME " + cfg.time + " rejected: " + strjoin(eT, "; "); end
span = ff.queryNum("SENS:FREQ:SPAN?");
if isempty(span) || span ~= 0, error('ff_5g_analyzer:zerospan', 'Zero span not accepted (span %g).', span); end
N = ff.queryNum("SENS:SWE:POIN?");
s = ff.saSettings();
T = NaN; src = "SWE:TIME?";
if isempty(eT), T = str2double(tryQuery(ff, "SENS:SWE:TIME?")); end
if ~(T > 0), T = str2double(tryQuery(ff, "SENS:SWE:MTIM?")); src = "SWE:MTIM? (may include overhead)"; end
if ~(T > 0), T = cfg.time; src = "requested (not confirmed)"; end
P = zeros(n, N); t = NaT(n, 1); k = 0;
trig = cfg.trigger;
if trig ~= "free"
    lvl = cfg.trigLevel;
    if trig == "video" && isempty(lvl)     % level midway between DL and UL from one free sweep
        ff.sweep(T + 15); y = ff.saTrace();
        if numel(y) == N
            k = 1; P(1, :) = y; t(1) = datetime('now');
            r = tddAnalyze((0:N - 1).'*T/(N - 1), y(:), 0.5e-3);
            lvl = r.thr;
        end
    end
    if trig == "video" && isempty(lvl)
        trig = "free"; notes(end + 1) = "video trigger level unknown: free run";
    else
        [ok, msg] = setTrigger(ff, trig, lvl);
        if ~ok, notes(end + 1) = msg; trig = "free"; end
    end
end
restoreTrig = onCleanup(@() setFree(ff));
while k < n
    k = k + 1;
    try
        if trig == "free", ff.sweep(T + 15); else, ff.sweep(T + 3); end
    catch err
        if trig ~= "free" && strcmp(err.identifier, 'FieldFox:timeout')
            recoverTrigger(ff, T);
            notes(end + 1) = sprintf("no %s trigger within %.1f s: free run", trig, T + 3); %#ok<AGROW>
            trig = "free"; k = k - 1;
            continue
        end
        rethrow(err)
    end
    y = ff.saTrace();
    if numel(y) ~= N, error('ff_5g_analyzer:points', 'Zero span: %d points, expected %d.', numel(y), N); end
    P(k, :) = y; t(k) = datetime('now');
end
clear restoreTrig
over = false;
if st('waor'), over = waor(ff); end
s.detector = det;
a = makeAcq(cfg, repmat(cfg.center, 1, N), P, t, s, extgInstr, over, true, T, src, trig);
a.settings.errors = strjoin(notes, "; ");
fprintf('   zero span %s: %d x %d points over %.4g s (%s), RBW %s, %s, %s%s\n', fmtHz6(cfg.center), n, N, T, ...
    src, fmtHz(s.rbw), det, trig, repmat('  OVERLOAD', 1, over));
if ~isempty(notes), warning('ff_5g_analyzer:zerospan', '%s', strjoin(notes, "; ")); end
if D.extGain ~= extgInstr, a.settings.levelNote = "ExtGain applied in MATLAB"; end
end

function [ok, msg] = setTrigger(ff, trig, lvl)
% PG: TRIGger[:SEQuence]:SOURce FREE|VIDeo|EXTernal, :VIDeo:LEVel, :SLOPe (no option needed).
msg = "";
if trig == "video"
    ff.write("TRIG:SOUR VID");
    ff.write(sprintf("TRIG:VID:LEV %.2f", lvl));
else
    ff.write("TRIG:SOUR EXT");
end
ff.write("TRIG:SLOP POS");
ff.opc(10);
e = ff.errors();
ok = isempty(e);
if ~ok
    msg = trig + " trigger not accepted (" + strjoin(e, "; ") + "): free run";
    setFree(ff);
end
end

function setFree(ff)
try
    ff.write("TRIG:SOUR FREE");
    ff.errors();
catch
end
end

function recoverTrigger(ff, T)
% The pending INIT:IMM;*OPC? completes once the trigger is free: wait, then drop its late reply.
try ff.write("TRIG:SOUR FREE"); catch, end
pause(min(3, 2*T + 0.5));
try flush(ff.Client); catch, end
try ff.query("*OPC?", 10); catch, end
pause(0.2);
try flush(ff.Client); catch, end
try ff.errors(); catch, end
end

function [chp, note] = instCHP(ff, c, a, extgInstr, extGain)
% FieldFox channel power over the occupied bandwidth on the ACLR sweep settings (verified:
% CALC:MEAS:DATA? follows FORM, CHP -> [P dBm, PSD dBm/Hz]). One single sweep.
chp = [NaN NaN]; note = "";
if ~(c.occ > 0), note = "no occupied bandwidth"; return, end
s = a.settings;
try
    ff.write("SENS:MEAS:CHAN CHP");
    ff.write(sprintf("SENS:CME:IBW %.10g", min(c.occ, s.stop - s.start)));
    ff.opc(30);
    e = ff.errors();
    ff.sweep(120);
    v = double(ff.queryBlock("CALC:MEAS:DATA?", "single", 10));
    v = v(:).'; v(v < -9e4) = NaN;
    if numel(v) >= 2, chp = v(1:2) + [extgInstr - extGain, extgInstr - extGain]; end
    if ~isempty(e), note = "CHP setup: " + strjoin(e, "; "); end
    if abs((s.start + s.stop)/2 - c.fc) > 0.02*c.occ, note = strjoin([note "span not centred on the carrier"], "; "); end
catch err
    note = "CHP failed: " + string(err.message);
end
try
    ff.write("SENS:MEAS:CHAN NONE"); ff.opc(30); ff.errors();
catch
end
end

function w = waor(ff)
try
    v = ff.queryNum("CALC:MEAS:WAOR?", 5);
    w = isscalar(v) && v == 1;
catch
    w = false;
end
end

function releaseInstrument(ff)
% FieldFox.release() only warns (FieldFox:release) when a restore step fails.
[m0, id0] = lastwarn; lastwarn('', '');
try
    ff.release();
    [~, id] = lastwarn;
    if strcmp(id, 'FieldFox:release')
        fprintf(2, ['Instrument NOT fully restored (connection lost?): its backup stays on the FieldFox and ' ...
            'is restored at the next takeControl (any suite tool).\n']);
    else
        fprintf('Instrument released (all settings restored).\n');
        lastwarn(m0, id0);
    end
catch err
    warning('ff_5g_analyzer:release', 'Could not restore the instrument: %s', err.message);
end
end

function r = tryQuery(ff, cmd)
% Reply of a query that may not exist on this firmware ("" = no reply after 5 s).
w = warning('off', 'all');
restore = onCleanup(@() warning(w));
try
    r = string(ff.query(cmd, 5));
catch
    r = "";
end
if isempty(r), r = ""; end
end

function ok = probe(ff, cmd)
r = tryQuery(ff, cmd);
e = ff.errors();                           % always drain the queue (an unknown query leaves -113)
ok = strlength(r) > 0 && isempty(e);
end

function a = makeAcq(cfg, f, P, t, s, extgInstr, over, zs, T, src, trig)
set = struct('label', cfg.label, 'start', s.start, 'stop', s.stop, 'center', (s.start + s.stop)/2, ...
    'span', s.stop - s.start, 'rbw', s.rbw, 'vbw', s.vbw, 'points', size(P, 2), 'atten', s.atten, ...
    'preamp', logical(s.preamp), 'detector', string(s.detector), 'refLevel', s.ref, 'extgInstr', extgInstr, ...
    'overload', over, 'zeroSpan', zs, 'sweepTime', T, 'timeFrom', src, 'trigger', string(trig), ...
    'errors', "", 'levelNote', "");
if zs, set.center = cfg.center; set.span = 0; end
a = struct('f', reshape(f, 1, []), 'P', P, 't', t, 'settings', set);
end

%% ===== demo engine (synthetic Thai scene + FieldFox measurement model) ============

function D = demoRun(opt)
fprintf(['DEMO: synthetic Thai NR/LTE scene through a FieldFox model (spectra analytic, zero-span ' ...
    'envelopes from 5G/LTE Toolbox waveforms).\n']);
W = demoWorld();
TT = tddTemplates(opt);
D = newRun("demo", opt);
D.idn = "SYNTHETIC DEMO (ff_5g_analyzer FieldFox model) - not a measurement";
D.options = ["233" "235"];
st = containers.Map('KeyType', 'char', 'ValueType', 'any');
st('preamp') = opt.Preamp; st('t0') = 1000*rand;
E.sweep = @(cfg, n) demoSweep(W, TT, st, cfg, n, opt);
E.zero = @(cfg, n) demoZero(W, st, cfg, n, opt);
E.chp = [];
D.notes(end + 1, 1) = "Demo: the FieldFox channel-power cross-check (CHP) is instrument-only.";
D = runFlow(E, D, opt);
D.truth = struct2table(rmfield(W, {'rip'}));
end

function W = demoWorld()
M = 1e6;
e = @(name, holder, tech, kind, fc, bw, scs, nrb, occ, tmpl, Pon, a1, a2) struct('name', string(name), ...
    'holder', string(holder), 'tech', string(tech), 'kind', string(kind), 'fc', fc, 'chBW', bw, 'scs', scs, ...
    'nrb', nrb, 'occ', occ, 'tmpl', string(tmpl), 'Pon', Pon, 'aclr1', a1, 'aclr2', a2, 'rip', rand(1, 3));
W = [ ...
    e("idle NR n28 (SSB only)", "AIS", "NR", "nrssb", 788*M, 20, 15, 106, 3.6*M, "SSB caseA L4 20ms", -64, 60, 65)
    e("LTE B3 20 MHz FDD", "AIS|True", "LTE", "ltefdd", 1845*M, 20, 15, 100, 18.0*M, "", -32, 44, 52)
    e("LTE B1 15 MHz FDD", "AIS", "LTE", "ltefdd", 2162.5*M, 15, 15, 75, 13.5*M, "", -72, 47, 55)
    e("LTE B40 20 MHz TDD", "True", "LTE", "ltetdd", 2310*M, 20, 15, 100, 18.0*M, "LTE cfg2 DSUDDDSUDD (SSC 7)", -66, 48, 52)
    e("NR n41 100 MHz TDD", "AIS", "NR", "nrtdd", 2550*M, 100, 30, 273, 273*12*30e3, "NR30 DDDSU (S 10:2:2)", -45, 46, 50)
    e("NR n41 90 MHz TDD", "True", "NR", "nrtdd", 2645*M, 90, 30, 245, 245*12*30e3, "NR30 DDDSU (S 10:2:2)", -54, 46, 50)];
W = W(:).';
end

function a = demoSweep(W, TT, st, cfg, n, opt)
N = cfg.points; f = linspace(cfg.start, cfg.stop, N); df = (cfg.stop - cfg.start)/(N - 1);
rbw = cfg.rbw; enbw = 1.0645*rbw;
pre = st('preamp');
nfl = 10^((danlOf(pre) - opt.ExtGain)/10)*enbw;        % displayed = DUT-referred (EXTG set)
Ts = max(0.02, 2.5*(cfg.stop - cfg.start)/rbw^2 + N*10e-6);
dwell = Ts/N;
dfi = min(rbw/4, max(df/4, 5e3));
fi = (cfg.start - 4*rbw - df:dfi:cfg.stop + 4*rbw + df).';
K = ceil(2.5*rbw/dfi); g = exp(-log(2)*(((-K:K)*dfi)/(rbw/2)).^2).';
nb = max(1, round(df/dfi));
spec = zeros(N, numel(W)); gated = false(1, numel(W)); tm = cell(1, numel(W));
for k = 1:numel(W)
    w = W(k);
    if w.fc + 4*w.chBW*1e6 < fi(1) || w.fc - 4*w.chBW*1e6 > fi(end), continue, end
    psd = 10^(w.Pon/10)/w.occ*demoShape(fi - w.fc, w);
    if ~any(psd), continue, end
    pr = conv(psd, g, 'same')*dfi;                      % mW in the RBW (Gaussian RBW, ENBW 1.0645 RBW)
    if nb > 1, pr = movmean(pr, nb); end                % average detector across the bucket
    spec(:, k) = interp1(fi, pr, f, 'linear', 0);
    if w.kind ~= "ltefdd", gated(k) = true; tm{k} = demoTemplate(w, TT); end
end
P = zeros(n, N); t = NaT(n, 1);
for s = 1:n
    t0 = st('t0') + 7.123*s + rand;                     % free running: random frame phase per sweep
    st('t0') = t0;
    x = repmat(nfl, 1, N);
    for k = find(any(spec, 1))
        gk = 1;
        if gated(k)
            a0 = t0 + (0:N - 1)*dwell; if W(k).kind == "nrssb", a0 = a0 + 0.37; end
            gk = (templateCum(tm{k}, a0 + dwell) - templateCum(tm{k}, a0))/dwell;
        end
        x = x + spec(:, k).'.*gk.*10.^(0.25*randn(1, N)/10);
    end
    P(s, :) = 10*log10(x.*10.^(0.3*randn(1, N)/10));
    t(s) = datetime('now');
end
sset = struct('start', cfg.start, 'stop', cfg.stop, 'rbw', rbw, 'vbw', rbw, 'atten', 0, 'preamp', pre, ...
    'detector', cfg.det, 'ref', opt.RefLevel);
a = makeAcq(cfg, f, P, t, sset, opt.ExtGain, false, false, NaN, "", "free");
a.settings.sweepTimeModel = Ts;
fprintf('   %-14s %s, RBW %s, %d pts (emulated)\n', cfg.label, fmtRange([cfg.start cfg.stop]), fmtHz(rbw), N);
end

function s = demoShape(x, w)
% Relative PSD of an emitter: flat (with a gentle multipath ripple) over the occupied band,
% leakage flat at -ACLR1 in the first and -ACLR2 in the second adjacent channel.
a = abs(x); bw = w.chBW*1e6; s = zeros(size(x));
in = a <= w.occ/2;
rip = 1.2*sin(2*pi*(x/max(w.occ, 1e6))*1.3 + 6*w.rip(1)) + 0.6*sin(2*pi*(x/max(w.occ, 1e6))*3.1 + 6*w.rip(2));
s(in) = 10.^(rip(in)/10);
if w.kind == "nrssb", return, end                     % idle cell: no leakage worth modelling
s(a > w.occ/2 & a <= 1.5*bw) = 10^(-w.aclr1/10);
s(a > 1.5*bw & a <= 2.5*bw) = 10^(-w.aclr2/10);
s(a > 2.5*bw & a <= 4*bw) = 10^(-(w.aclr2 + 12)/10);
end

function T = demoTemplate(w, TT)
if w.kind == "nrssb"
    T = ssbTemplate(15, 4, 20e-3);
else
    T = TT(find([TT.name] == w.tmpl, 1));
end
end

function a = demoZero(W, st, cfg, n, opt)
fz = cfg.center; T = cfg.time; N = cfg.points; pre = st('preamp');
danlHz = danlOf(pre) - opt.ExtGain;
k = find(arrayfun(@(w) abs(fz - w.fc) < w.occ/2, W));
if numel(k) > 1, [~, i] = max(arrayfun(@(w) w.Pon - 10*log10(w.occ), W(k))); k = k(i); end
P = zeros(n, N); t = NaT(n, 1);
if ~isempty(k)
    w = W(k);
    [x, fs, occW] = demoWaveform(w);
    pw = abs(x).^2; L = max(1, round(fs*30e-6)); env = movmean(pw, L);
    on = env > 0.2*max(env);
    x = x*sqrt(10^(w.Pon/10)/w.occ*occW/mean(pw(on)));
end
for s = 1:n
    if isempty(k)
        p = 10^(danlHz/10)*1.0645*cfg.rbw*(-log(rand(N, 1)));
        if cfg.det ~= "SAMP", p = 10^(danlHz/10)*1.0645*cfg.rbw*(1 + 0.15*randn(N, 1)); end
        p = 10*log10(abs(p));
    else
        [~, p] = zsEmulate(circshift(x, randi(numel(x))), fs, 0, cfg.rbw, T, N, cfg.det, danlHz, true);
    end
    P(s, :) = p; t(s) = datetime('now');
end
sset = struct('start', fz, 'stop', fz, 'rbw', cfg.rbw, 'vbw', cfg.vbw, 'atten', 0, 'preamp', pre, ...
    'detector', cfg.det, 'ref', opt.RefLevel);
a = makeAcq(cfg, repmat(fz, 1, N), P, t, sset, opt.ExtGain, false, true, T, "emulated", "free");
fprintf('   zero span %s: %d x %d points over %.4g s, %s (emulated)\n', fmtHz6(fz), n, N, T, cfg.det);
end

function [x, fs, occW] = demoWaveform(w)
% Toolbox waveforms for the emulated zero span (cached for the MATLAB session).
persistent cache
if isempty(cache), cache = containers.Map(); end
key = char(w.kind + "|" + w.tmpl);
if isKey(cache, key), v = cache(key); x = v{1}; fs = v{2}; occW = v{3}; return, end
switch w.kind
    case "nrtdd"
        tk = regexp(w.tmpl, 'NR(\d+) (\w+) \(S (\d+):(\d+):(\d+)\)', 'tokens', 'once');
        cfg = nrTDDConfig(20, str2double(tk{1}), tk{2}, str2double(tk(3:5)), 20, true);
        [x, info] = nrWaveformGenerator(cfg); fs = info.ResourceGrids(1).Info.SampleRate;
        occW = cfg.SCSCarriers{1}.NSizeGrid*12*cfg.SCSCarriers{1}.SubcarrierSpacing*1e3;
    case "ltetdd"
        tk = regexp(w.tmpl, 'cfg(\d).*SSC (\d)', 'tokens', 'once');
        rmc = lteRMCDL('R.9'); rmc.NDLRB = 100; rmc.DuplexMode = 'TDD';
        rmc.TDDConfig = str2double(tk{1}); rmc.SSC = str2double(tk{2}); rmc.TotSubframes = 20;
        rmc.OCNGPDSCHEnable = 'On'; rmc.OCNGPDCCHEnable = 'On';      % loaded cell (no empty subframes)
        rmc = lteRMCDL(rmc);
        [x, ~, info] = lteRMCDLTool(rmc, [1; 0; 0; 1]); fs = info.SamplingRate; occW = 100*180e3;
    case "ltefdd"
        tm = lteTestModel('1.1', '20MHz');
        [x, ~, info] = lteTestModelTool(tm); fs = info.SamplingRate; occW = 100*180e3;
    otherwise                                         % idle NR: SSB only (Case A, L = 4, 20 ms)
        cfg = nrDLCarrierConfig; cfg.FrequencyRange = 'FR1'; cfg.ChannelBandwidth = 20; cfg.NumSubframes = 40;
        cfg.SCSCarriers{1}.SubcarrierSpacing = 15; cfg.SCSCarriers{1}.NSizeGrid = 106;
        cfg.BandwidthParts{1}.SubcarrierSpacing = 15; cfg.BandwidthParts{1}.NSizeBWP = 106;
        cfg.SSBurst.BlockPattern = 'Case A'; cfg.SSBurst.TransmittedBlocks = [1 1 1 1];
        cfg.SSBurst.Period = 20; cfg.SSBurst.SubcarrierSpacingCommon = 15;
        cfg.PDSCH{1}.Enable = false; cfg.PDCCH{1}.Enable = false; cfg.CSIRS{1}.Enable = false;
        [x, info] = nrWaveformGenerator(cfg); fs = info.ResourceGrids(1).Info.SampleRate; occW = 240*15e3;
end
x = x(:, 1);
cache(key) = {x, fs, occW};
end

function E = demoEMIRecording(opt)
% Synthetic inverter-like interference for the coexistence demo: a broadband hump at 2.4-2.7 GHz
% (switching edges + a monitoring dongle) and harmonics near 800 MHz. Shared format, AVER detector.
f = linspace(700e6, 2700e6, 4001); rbw = 1e6; enbw = 1.0645*rbw;
fl = danlOf(opt.Preamp) - opt.ExtGain;
emi = -175 + 27*exp(-((f - 2.56e9)/90e6).^2) + 17*exp(-((f - 0.8e9)/40e6).^2);   % dBm/Hz, peaks -148 / -158
on = 10*log10((10.^(fl/10) + 10.^(emi/10))*enbw) + 0.3*randn(3, numel(f));
off = 10*log10(10.^(fl/10)*enbw) + 0.3*randn(3, numel(f));
mk = @(P, lbl) struct('f', f, 'P', P, 't', datetime('now') + seconds(0:2).', 'settings', ...
    struct('label', lbl, 'rbw', rbw, 'detector', "AVER", 'preamp', opt.Preamp, 'atten', 0, ...
    'extgInstr', opt.ExtGain, 'zeroSpan', false), 'idn', "SYNTHETIC DEMO EMI");
E = struct('on', mk(on, "emi_demo_on"), 'off', mk(off, "emi_demo_off"));
end

%% ===== file input ==================================================================

function D = loadInput(file, opt)
file = string(file);
if ~isfile(file), error('ff_5g_analyzer:file', 'File not found: %s', file); end
S = load(file);
if isfield(S, 'D') && isstruct(S.D) && isfield(S.D, 'kind') && string(S.D.kind) == "ff_5g_analyzer"
    D = S.D;
    D.notes(end + 1, 1) = "Re-analysis of " + file + " (measured " + string(D.time) + ").";
    fprintf('Re-analysing %s: %d scans, %d zero-span sets, %d ACLR sweeps.\n', file, numel(D.scan), numel(D.zs), numel(D.aclr));
    return
end
if ~all(isfield(S, {'f', 'P'}))
    error('ff_5g_analyzer:file', '%s is neither an ff_5g_analyzer run nor a recording with f and P.', file);
end
D = newRun("recording", opt);
if isfield(S, 'idn'), D.idn = string(S.idn); end
a = recordingAcq(S);
D.notes(end + 1, 1) = "Recording " + file + ": " + a.settings.levelNote;
[~, nm] = fileparts(file);
if a.settings.zeroSpan
    if ~(a.settings.sweepTime > 0), error('ff_5g_analyzer:file', 'Zero-span recording without a sweep time.'); end
    D.zs(1) = struct('fc', a.settings.center, 'fz', a.settings.center, 'occ', NaN, 'acq', a, 'ssb', []);
    fprintf('Zero-span recording %s: %d sweeps of %.4g s.\n', file, size(a.P, 1), a.settings.sweepTime);
else
    D.scan(1) = struct('key', "file", 'name', string(nm), 'acq', a);
    fprintf('Recording %s: %s, %d sweeps x %d points.\n', file, fmtRange(a.f([1 end])), size(a.P, 1), size(a.P, 2));
end
end

function a = recordingAcq(S)
f = double(S.f(:).'); P = double(S.P);
if size(P, 2) ~= numel(f) && size(P, 1) == numel(f), P = P.'; end
if size(P, 2) ~= numel(f), error('ff_5g_analyzer:file', 'P has %d columns but f has %d points.', size(P, 2), numel(f)); end
s = struct();
if isfield(S, 'settings') && isstruct(S.settings), s = S.settings; end
g = @(nm, d) getf(s, nm, d);
df = 0; if numel(f) > 1, df = (f(end) - f(1))/(numel(f) - 1); end
zs = logical(g('zeroSpan', df == 0)) || g('span', NaN) == 0 || df == 0;
rbw = g('rbw', NaN);
note = "";
if ~(rbw > 0)
    rbw = max(df, 1); note = "no RBW in the recording: bin spacing used";
end
T = g('sweepTime', g('swtime', g('time', NaN)));
t = NaT(size(P, 1), 1); if isfield(S, 't') && numel(S.t) == size(P, 1) && isdatetime(S.t), t = S.t(:); end
set = struct('label', string(g('label', "recording")), 'start', f(1), 'stop', f(end), 'center', g('center', mean(f([1 end]))), ...
    'span', f(end) - f(1), 'rbw', rbw, 'vbw', g('vbw', rbw), 'points', numel(f), 'atten', g('atten', 0), ...
    'preamp', logical(g('preamp', true)), 'detector', string(g('detector', "unknown")), 'refLevel', g('refLevel', NaN), ...
    'extgInstr', g('extgInstr', 0), 'overload', logical(g('overload', false)), 'zeroSpan', zs, ...
    'sweepTime', T, 'timeFrom', "recording", 'trigger', string(g('trigger', "free")), 'errors', "", ...
    'levelNote', strjoin([note "levels as recorded, ExtGain applied in MATLAB"], "; "));
if zs, set.span = 0; end
if ~isnumeric(set.atten) || ~isfinite(set.atten), set.atten = 0; end
a = struct('f', f, 'P', P, 't', t, 'settings', set);
end

function v = getf(s, nm, d)
v = d;
if isfield(s, nm) && ~isempty(s.(nm)), v = s.(nm); end
if isstring(d) || ischar(d), v = string(v); elseif islogical(d) || isnumeric(d), v = double(v); end
if isnumeric(v) && ~isscalar(v), v = v(1); end
end

%% ===== analysis ======================================================================

function A = analyzeRun(D, opt)
A = struct('notes', strings(0, 1));
A.templates = tddTemplates(opt);
A.G = floorGroups(D);
C = emptyCarrier();
for i = 1:numel(D.scan)
    C = [C, detectInAcq(D.scan(i).acq, D, A.G, opt, "scan", i)]; %#ok<AGROW>
end
C = mergeAclr(C, D, A.G, opt);
C = linkZeroSpan(C, D);
for k = 1:numel(C)
    C(k).id = k;
end
for k = 1:numel(C)
    C(k) = analyzeCarrier(C(k), C, D, A, opt);
end
A.carriers = C;
A.coex = coexistence(C, D, opt);
A.table = summaryTable(C);
A.conclusions = conclusions(C, D, A);
if D.source == "demo", A.truth = D.truth; end
end

function C = mergeAclr(C, D, G, opt)
% The ACLR sweep is centred on the carrier with a finer grid: it refines the spectral estimate
% (or adds the carrier when no scan found it, e.g. 'Carriers' given).
for k = 1:numel(D.aclr)
    X = D.aclr(k);
    cs = detectInAcq(X.acq, D, G, opt, "aclr", k);
    j = [];
    if ~isempty(cs)
        j = find([cs.f1] <= X.fc & [cs.f2] >= X.fc, 1);
        if isempty(j)
            [d, j] = min(abs([cs.fc] - X.fc));
            if d > 0.3*max(X.occ, 5e6), j = []; end
        end
    end
    i = [];
    if ~isempty(C)
        i = find([C.f1] <= X.fc & [C.f2] >= X.fc, 1);
        if isempty(i) && ~isempty(j), i = find([C.f1] <= cs(j).fc & [C.f2] >= cs(j).fc, 1); end
    end
    if ~isempty(j)
        c = cs(j);
        if ~isempty(i)                     % keep the scan's band, bursty flag and sweep spread
            c.bursty = C(i).bursty || c.bursty; c.fcPlan = C(i).fc;
            if ~isfinite(c.sweepStd), c.sweepStd = C(i).sweepStd; end
            C(i) = c;
        else
            C(end + 1) = c; %#ok<AGROW>
        end
        if ~isempty(i), C(i).aclrIdx = k; else, C(end).aclrIdx = k; end
    elseif isempty(i)                      % nothing detected: keep the planned carrier
        c = newCarrier();
        c.src = "aclr (not detected)"; c.fc = X.fc; c.fcPlan = X.fc; c.occ = X.occ;
        c.f1 = X.fc - X.occ/2; c.f2 = X.fc + X.occ/2; c.aclrIdx = k; c.chBW = X.chBW;
        [c.band, c.bandName, c.duplexBand] = bandOf(c.fc);
        c = fromPlan(c, D);
        c.notes(end + 1) = "not detected above the threshold in the ACLR sweep";
        C(end + 1) = c; %#ok<AGROW>
    else
        C(i).aclrIdx = k;
    end
end
if ~isempty(C), [~, o] = sort([C.fc]); C = C(o); end
end

function c = fromPlan(c, D)
% A carrier known only from the plan ('Carriers'): keep whether its bandwidth was given or assumed.
if ~isfield(D, 'planned') || isempty(D.planned) || ~isfield(D.planned, 'occSrc'), return, end
j = find(abs([D.planned.fc] - c.fc) < 1, 1);
if isempty(j) || strlength(D.planned(j).occSrc) == 0, return, end
c.occSrc = D.planned(j).occSrc; c.chBW = D.planned(j).chBW;
if c.occSrc == "assumed", c.notes(end + 1) = "channel bandwidth not given and not measured"; end
end

function C = linkZeroSpan(C, D)
for k = 1:numel(D.zs)
    z = D.zs(k);
    i = [];
    if ~isempty(C)
        i = find([C.f1] <= z.fc & [C.f2] >= z.fc);
        if numel(i) > 1, [~, m] = min(abs([C(i).fc] - z.fc)); i = i(m); end
    end
    if isempty(i)
        c = newCarrier();
        c.src = "zero span only"; c.fc = z.fc; c.fcPlan = z.fc; c.occ = z.occ;
        c.f1 = z.fc - z.occ/2; c.f2 = z.fc + z.occ/2;
        [c.band, c.bandName, c.duplexBand] = bandOf(c.fc);
        c = fromPlan(c, D);
        if ~(z.occ > 0), c.notes(end + 1) = "no spectrum: bandwidth unknown"; end
        c.zsIdx = k;
        C(end + 1) = c; %#ok<AGROW>
    else
        C(i).zsIdx = k;
    end
end
if ~isempty(C), [~, o] = sort([C.fc]); C = C(o); end
end

function c = analyzeCarrier(c, C, D, A, opt)
if strlength(c.band) == 0, [c.band, c.bandName, c.duplexBand] = bandOf(c.fc); end
[c.holder, c.holderConf] = operatorOf(c.f1, c.f2, c.fc);
if c.occ > 0 && c.occSrc ~= "assumed", c.snap = snapCandidates(c.occ, c.band); end
if c.occSrc == "given"                     % width not measured: every configuration of that BW fits
    S = c.snap(c.snap.tech ~= "NR-SSB" & abs(c.snap.chBW - c.chBW) < 1e-6, :);
    if isempty(S)
        S = bwTable(); S = S(S.tech ~= "NR-SSB" & abs(S.chBW - c.chBW) < 1e-6, :); S.err = zeros(height(S), 1);
        c.notes(end + 1) = sprintf("%g MHz is not a 3GPP bandwidth of band %s", c.chBW, c.band);
    end
    S.err(:) = 0; c.snap = S;
end
if c.zsIdx > 0
    c.zsr = analyzeZero(D.zs(c.zsIdx), D, A, opt);
end
c = decideTech(c);
% ACLR: the dedicated sweep, else the scan if it covers the adjacent channels
acq = []; fl = NaN;
if c.aclrIdx > 0, acq = D.aclr(c.aclrIdx).acq; end
if isempty(acq) && c.scanIdx > 0, acq = D.scan(c.scanIdx).acq; end
if ~isempty(acq), fl = floorHzOf(acq, D, A.G); end
if ~isempty(acq) && c.chBW > 0 && c.measBW > 0
    c.aclr = aclrCompute(c, acq, D, fl, C, opt);
end
if opt.Reference && ~isempty(acq) && c.chBW > 0
    try
        c.ref = refFit(c, acq, D, fl);
    catch err
        c.notes(end + 1) = "reference spectrum failed: " + string(err.message);
    end
end
if c.aclrIdx > 0 && all(isfinite(D.aclr(c.aclrIdx).chp))
    c.chp = D.aclr(c.aclrIdx).chp;
end
% raster numbers
if c.tech == "LTE", c.arfcn = "EARFCN " + earfcn(c.fc, c.band); end
if startsWith(c.tech, "NR")
    c.arfcn = "NR-ARFCN " + nrarfcn(c.fc);
    [g, fg] = gscnNearest(c.fc);
    c.gscn = sprintf("GSCN %d (%s)", g, fmtHz6(fg));
end
if c.scs > 0 && isfinite(c.psd)
    d = 1; if c.duplex == "TDD" && ~isempty(c.zsr) && c.zsr.dlDuty > 0.2, d = c.zsr.dlDuty; end
    c.epre = c.psd - 10*log10(d) + 10*log10(c.scs*1e3);
end
end

%% ---- noise floor ----------------------------------------------------------------

function d = danlOf(preamp)
% Displayed average noise level of this N9914A (dBm/Hz at 0 dB attenuation; task-context values).
if preamp, d = -153; else, d = -137; end
end

function fl = floorModel(set, D)
att = set.atten; if ~isnumeric(att) || ~isfinite(att), att = 0; end
fl = danlOf(logical(set.preamp)) + att - D.extGain;
end

function G = floorGroups(D)
% Empirical floor (5th percentile of the smoothed PSD) of every swept acquisition, by setting.
G = struct('preamp', {}, 'att', {}, 'emp', {});
acqs = [arrayfun(@(s) {s.acq}, D.scan), arrayfun(@(s) {s.acq}, D.aclr)];
for k = 1:numel(acqs)
    a = acqs{k};
    if a.settings.zeroSpan, continue, end
    G(end + 1) = struct('preamp', logical(a.settings.preamp), 'att', a.settings.atten, 'emp', empFloor(a, D)); %#ok<AGROW>
end
end

function e = empFloor(a, D)
f = a.f(:); df = (f(end) - f(1))/max(1, numel(f) - 1);
lin = mean(10.^(levels(a, D)/10), 1).';
sm = 10*log10(movmean(lin, max(1, round(500e3/df)))/(1.0645*a.settings.rbw));
e = prctile(sm, 5);
end

function fl = floorHzOf(a, D, G)
% DUT-referred analyzer floor in dBm/Hz: the lowest of the empirical floor, the DANL model + 3 dB
% and the best empirical floor of acquisitions with the same preamp / attenuation + 1.5 dB.
m = floorModel(a.settings, D);
c = m + 3;
if ~a.settings.zeroSpan, c = [c empFloor(a, D)]; end
if ~isempty(G)
    same = [G.preamp] == logical(a.settings.preamp) & [G.att] == a.settings.atten;
    if any(same), c = [c min([G(same).emp]) + 1.5]; end
end
fl = min(c);
end

function P = levels(a, D)
% DUT-referred dBm: display + external gain set on the FieldFox - total ExtGain.
P = a.P + a.settings.extgInstr - D.extGain;
end

%% ---- carrier detection ------------------------------------------------------------

function c = newCarrier()
c = struct('id', 0, 'src', "", 'band', "", 'bandName', "", 'duplexBand', "", 'holder', "", 'holderConf', "", ...
    'fc', NaN, 'f1', NaN, 'f2', NaN, 'occ', NaN, 'occCI', NaN, 'obw99', NaN, 'chP', NaN, 'psd', NaN, ...
    'epre', NaN, 'sweepStd', NaN, 'bursty', false, 'truncated', false, 'fcPlan', NaN, 'floorHz', NaN, ...
    'scanIdx', 0, 'aclrIdx', 0, 'zsIdx', 0, 'snap', table(), 'tech', "?", 'techLabel', "?", 'duplex', "?", ...
    'chBW', NaN, 'scs', NaN, 'nrb', NaN, 'measBW', NaN, 'snapErr', NaN, 'occSrc', "", ...
    'zsr', [], 'aclr', [], 'ref', [], 'chp', [NaN NaN], ...
    'pattern', "", 'patternScore', NaN, 'confidence', "low", 'evidence', strings(1, 0), 'arfcn', "", ...
    'gscn', "", 'notes', strings(1, 0));
end

function C = emptyCarrier()
C = repmat(newCarrier(), 1, 0);
end

function C = detectInAcq(a, D, G, opt, kind, idx)
% Carriers in one swept acquisition: threshold over the floor on the 500 kHz-smoothed mean PSD,
% split at notches / level steps, -3 dB edges refined by an edge fit. Regions seen only in the
% max over the sweeps are kept as bursty (idle NR: SSB only).
C = emptyCarrier();
f = a.f(:); N = numel(f);
if N < 16 || a.settings.zeroSpan, return, end
P = levels(a, D); lin = 10.^(P/10);
mlin = mean(lin, 1).'; xlin = max(lin, [], 1).';
rbw = a.settings.rbw; enbw = 1.0645*rbw; df = (f(end) - f(1))/(N - 1);
fl = floorHzOf(a, D, G);
w = max(1, round(500e3/df));
sm = 10*log10(movmean(mlin, w)/enbw); smx = 10*log10(movmean(xlin, w)/enbw);
w2 = max(1, round(max(rbw, 150e3)/df));
s2 = 10*log10(movmean(mlin, w2)/enbw); s2x = 10*log10(movmean(xlin, w2)/enbw);
pHz = 10*log10(mlin/enbw);
thr = fl + opt.Threshold;
gap = max(1, round(300e3/df));
maskM = sm > thr; maskX = smx > thr + 3;
reg = runs(maskM | maskX, gap, max(2, round(0.8*opt.MinBW/df)));
bandKey = bandOf(mean(f([1 end])));
parts = zeros(0, 3);
for r = 1:size(reg, 1)
    if mean(maskM(reg(r, 1):reg(r, 2))) < 0.5      % mostly in the max over the sweeps: bursty
        parts = [parts; reg(r, :), 1]; %#ok<AGROW>
        continue
    end
    i1 = find(maskM(reg(r, 1):reg(r, 2)), 1) + reg(r, 1) - 1;
    i2 = find(maskM(reg(r, 1):reg(r, 2)), 1, 'last') + reg(r, 1) - 1;
    pr = splitRegion(f, s2, i1, i2, opt, 0);
    parts = [parts; pr, zeros(size(pr, 1), 1)]; %#ok<AGROW>
end
for p = 1:size(parts, 1)
    i1 = parts(p, 1); i2 = parts(p, 2); isB = parts(p, 3) == 1;
    if isB, sx = s2x; py = 10*log10(xlin/enbw); ml = xlin; else, sx = s2; py = pHz; ml = mlin; end
    [fL, fR] = edges3(f, sx, i1, i2);
    if ~(fR > fL), continue, end
    win = min(max([1.5e6, 6*rbw, 8*df]), (fR - fL)/2);
    [fe1, ci1] = fitEdge(f, py, fL, +1, rbw, df, win);
    [fe2, ci2] = fitEdge(f, py, fR, -1, rbw, df, win);
    if ~(fe2 - fe1 > 0.5*(fR - fL)), fe1 = fL; fe2 = fR; ci1 = NaN; ci2 = NaN; end
    c = newCarrier();
    c.src = kind; c.f1 = fe1; c.f2 = fe2; c.fc = (fe1 + fe2)/2; c.occ = fe2 - fe1; c.occCI = hypot(ci1, ci2);
    if c.occ < 0.8*opt.MinBW, continue, end
    c.bursty = isB; c.truncated = i1 <= 2 || i2 >= N - 1; c.floorHz = fl;
    m = f >= fe1 & f <= fe2;
    pc = sum(ml(m))*df/enbw - 10^(fl/10)*c.occ;
    c.chP = 10*log10(max(pc, 1e-30)); c.psd = c.chP - 10*log10(c.occ);
    mo = f >= fe1 - 0.15*c.occ & f <= fe2 + 0.15*c.occ;
    try
        x = max(ml(mo) - 10^(fl/10)*enbw, 0);
        if nnz(x) > 4, c.obw99 = obw(x, f(mo)); end
    catch
    end
    if size(P, 1) > 1, c.sweepStd = median(std(P(:, m), 0, 1)); end
    [c.band, c.bandName, c.duplexBand] = bandOf(c.fc);
    if strlength(c.band) == 0, c.band = bandKey; end
    if kind == "scan", c.scanIdx = idx; else, c.aclrIdx = idx; end
    if c.truncated, c.notes(end + 1) = "carrier touches the sweep edge: bandwidth may be truncated"; end
    C(end + 1) = c; %#ok<AGROW>
end
end

function R = runs(mask, gap, minLen)
% [start end] index pairs of true runs, gaps <= gap merged, runs shorter than minLen dropped.
mask = mask(:); d = diff([false; mask; false]);
s = find(d == 1); e = find(d == -1) - 1;
k = 1;
while k < numel(s)
    if s(k + 1) - e(k) - 1 <= gap, e(k) = e(k + 1); s(k + 1) = []; e(k + 1) = []; else, k = k + 1; end
end
keep = (e - s + 1) >= minLen;
R = [s(keep) e(keep)];
end

function [fL, fR, plat] = edges3(f, s, i1, i2)
% Outermost -3 dB points relative to the plateau (median of the central 60 %), interpolated.
n = i2 - i1 + 1;
c1 = i1 + floor(0.2*n); c2 = i2 - floor(0.2*n);
plat = median(s(c1:max(c1, c2)));
lv = plat - 3;
iL = find(s(i1:i2) >= lv, 1) + i1 - 1; iR = find(s(i1:i2) >= lv, 1, 'last') + i1 - 1;
if isempty(iL), fL = f(i1); fR = f(i2); return, end
fL = f(iL); fR = f(iR);
if iL > 1 && s(iL) ~= s(iL - 1), fL = f(iL - 1) + (lv - s(iL - 1))*(f(iL) - f(iL - 1))/(s(iL) - s(iL - 1)); end
if iR < numel(f) && s(iR) ~= s(iR + 1), fR = f(iR) + (s(iR) - lv)*(f(iR + 1) - f(iR))/(s(iR) - s(iR + 1)); end
end

function parts = splitRegion(f, s, i1, i2, opt, depth)
% Split a region at a notch or a level step when both parts fit 3GPP bandwidths better than the
% whole (frequency-selective fading over the air also makes notches: they must not split).
parts = [i1 i2];
df = f(2) - f(1);
if depth > 4 || f(i2) - f(i1) < 2*opt.MinBW || i2 - i1 < 12, return, end
bandKey = bandOf((f(i1) + f(i2))/2);
[fL, fR] = edges3(f, s, i1, i2);
e0 = snapErr(fR - fL, bandKey);
cand = zeros(0, 3);                         % [index, prominence dB, step dB]
seg = s(i1:i2);
[~, loc, ~, prom] = findpeaks(-seg, 'MinPeakProminence', 2);
loc = loc + i1 - 1;
ok = f(loc) - f(i1) >= opt.MinBW & f(i2) - f(loc) >= opt.MinBW;
cand = [cand; loc(ok), prom(ok), zeros(nnz(ok), 1)];
wN = max(2, round(2e6/df));
if i2 - i1 > 4*wN
    mL = movmedian(s, [wN 0]); mR = movmedian(s, [0 wN]);
    j = (i1 + wN:i2 - wN).';
    dlt = mR(j + 1) - mL(j - 1);
    if max(abs(dlt)) > 6 && numel(dlt) > 2*wN + 2
        [pk, lj] = findpeaks(abs(dlt), 'MinPeakHeight', 6, 'MinPeakDistance', wN);
        cand = [cand; j(lj), zeros(numel(lj), 1), pk(:)];
    end
end
if isempty(cand), return, end
best = Inf; bj = 0; bigStep = false;
for k = 1:size(cand, 1)
    j = cand(k, 1);
    [a1, b1] = edges3(f, s, i1, j); [a2, b2] = edges3(f, s, j, i2);
    e = max(snapErr(b1 - a1, bandKey), snapErr(b2 - a2, bandKey));
    step = cand(k, 3) >= 8;
    accept = step || (e < 0.06 && e < e0 - 0.01 && (e0 > 0.03 || cand(k, 2) >= 6));
    if accept && (e < best || (step && ~bigStep))
        best = e; bj = j; bigStep = step;
    end
end
if bj > 0
    parts = [splitRegion(f, s, i1, bj, opt, depth + 1); splitRegion(f, s, bj, i2, opt, depth + 1)];
end
end

function [fe, ci] = fitEdge(f, y, fe0, side, rbw, df, win)
% One carrier edge: y = 10 log10(L * Phi(side*(f - fe)/sigma) + N), sigma from the Gaussian RBW and
% the bucket width; lsqcurvefit (Optimization Toolbox), 95 % interval from nlparci scaled for
% RBW-correlated bins. Falls back to the -3 dB point.
fe = fe0; ci = NaN;
m = abs(f - fe0) <= win;
if nnz(m) < 6, return, end
x = f(m); yy = y(m);
sig = sqrt((rbw/2.3548)^2 + df^2/12);
mdl = @(p, x) 10*log10(10.^(p(2)/10).*0.5.*erfc(-side*(x - p(1)*1e6)/(sqrt(2)*sig)) + 10.^(p(3)/10));
hi = max(yy); lo = min(yy);
p0 = [fe0/1e6, hi, lo];
lb = [(fe0 - win)/1e6, hi - 15, lo - 20]; ub = [(fe0 + win)/1e6, hi + 10, hi];
try
    o = optimoptions('lsqcurvefit', 'Display', 'off', 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-12);
    [p, ~, res, flag, ~, ~, J] = lsqcurvefit(mdl, p0, x, yy, lb, ub, o);
    if flag <= 0, return, end
    fe = p(1)*1e6;
    try
        cpar = nlparci(p, res, 'jacobian', full(J));
        ci = (cpar(1, 2) - cpar(1, 1))/2*1e6*sqrt(max(1, 1.0645*rbw/df));
    catch
    end
catch
    p = fminsearch(@(p) sum((mdl(p, x) - yy).^2), p0, optimset('Display', 'off'));
    if abs(p(1)*1e6 - fe0) <= win, fe = p(1)*1e6; end
end
end

%% ---- bandwidth tables, bands, operators -----------------------------------------------

function T = bwTable()
% Transmission bandwidth configurations, FR1 (TS 38.104 Table 5.3.2-1) and E-UTRA (TS 36.104
% Table 5.6-1). The NR rows equal nrDLCarrierConfig.FR1BandwidthTable (R2026a).
nr15 = [5 25; 10 52; 15 79; 20 106; 25 133; 30 160; 35 188; 40 216; 45 242; 50 270];
nr30 = [5 11; 10 24; 15 38; 20 51; 25 65; 30 78; 35 92; 40 106; 45 119; 50 133; 60 162; 70 189; 80 217; 90 245; 100 273];
lte = [1.4 6; 3 15; 5 25; 10 50; 15 75; 20 100];
tech = [repmat("NR", size(nr15, 1) + size(nr30, 1), 1); repmat("LTE", size(lte, 1), 1); "NR-SSB"; "NR-SSB"];
chBW = [nr15(:, 1); nr30(:, 1); lte(:, 1); NaN; NaN];
scs = [15*ones(size(nr15, 1), 1); 30*ones(size(nr30, 1), 1); 15*ones(size(lte, 1), 1); 15; 30];
nrb = [nr15(:, 2); nr30(:, 2); lte(:, 2); 20; 20];
occ = nrb.*12.*scs*1e3;
T = table(tech, chBW, scs, nrb, occ);
end

function S = snapCandidates(W, bandKey)
% 3GPP configurations allowed in the band, sorted by |W - occupied|/occupied.
T = bwTable();
B = bandTable();
i = find(B.key == bandKey, 1);
if ~isempty(i)
    okNR = T.tech == "NR" & ismember(T.chBW, B.nrbw{i});
    okLTE = T.tech == "LTE" & ismember(T.chBW, B.ltebw{i});
    okSSB = T.tech == "NR-SSB";
    keep = okNR | okLTE | okSSB;
    if B.duplex(i) == "TDD", keep = keep & ~(T.tech == "NR" & T.scs == 15 & T.chBW > 20); end
    T = T(keep, :);
end
T.err = abs(W - T.occ)./T.occ;
T = sortrows(T, 'err');
S = T;
end

function e = snapErr(W, bandKey)
S = snapCandidates(W, bandKey);
S = S(S.tech ~= "NR-SSB", :);
if isempty(S), e = Inf; else, e = S.err(1); end
end

function B = bandTable()
% Thai IMT bands (2026): DL / TDD / SDL bands and uplinks, LTE and NR channel bandwidths (MHz)
% from TS 36.101 Table 5.6.1-1 / TS 38.101-1 Table 5.3.5-1 (approximate, Rel-17).
M = 1e6;
r = {
"B28",   "700 B28/n28 DL",   758*M,  803*M,  "FDD", 28, 28, [3 5 10 15 20],     [5 10 15 20 30]
"B5",    "850 B5/n5 DL",     869*M,  884*M,  "FDD", 5,  5,  [1.4 3 5 10],       [5 10 15]
"B8",    "900 B8/n8 DL",     935*M,  960*M,  "FDD", 8,  8,  [1.4 3 5 10],       [5 10 15 20]
"L1500", "1500 B32/n75 SDL", 1452*M, 1472*M, "SDL", 32, 75, [5 10 15 20],       [5 10 15 20]
"B3",    "1800 B3/n3 DL",    1805*M, 1880*M, "FDD", 3,  3,  [1.4 3 5 10 15 20], [5 10 15 20 25 30 40]
"B1",    "2100 B1/n1 DL",    2110*M, 2170*M, "FDD", 1,  1,  [5 10 15 20],       [5 10 15 20 25 30 40 50]
"n40",   "2300 B40/n40 TDD", 2300*M, 2370*M, "TDD", 40, 40, [5 10 15 20],       [5 10 15 20 25 30 40 50 60 70]
"n41",   "2600 B41/n41 TDD", 2496*M, 2690*M, "TDD", 41, 41, [5 10 15 20],       [10 15 20 30 40 50 60 70 80 90 100]
"n78",   "3500 n78 TDD",     3300*M, 3800*M, "TDD", NaN, 78, [],                [10 15 20 30 40 50 60 70 80 90 100]
"B28U",  "700 B28/n28 UL",   703*M,  748*M,  "UL",  28, 28, [3 5 10 15 20],     [5 10 15 20 30]
"B8U",   "900 B8/n8 UL",     890*M,  915*M,  "UL",  8,  8,  [1.4 3 5 10],       [5 10 15 20]
"B3U",   "1800 B3/n3 UL",    1710*M, 1785*M, "UL",  3,  3,  [1.4 3 5 10 15 20], [5 10 15 20 25 30 40]
"B1U",   "2100 B1/n1 UL",    1920*M, 1980*M, "UL",  1,  1,  [5 10 15 20],       [5 10 15 20 25 30 40 50]
};
B = cell2table(r, 'VariableNames', {'key', 'name', 'flo', 'fhi', 'duplex', 'lteBand', 'nrBand', 'ltebw', 'nrbw'});
B.key = string(B.key); B.name = string(B.name); B.duplex = string(B.duplex);
end

function H = holderTable()
% Holdings after the 29 Jun 2025 NBTC auction (licences from 4 Aug 2025). conf: high / med /
% low ("low" = block edges inferred, not published).
M = 1e6;
r = {
758*M,  778*M,  "True", "med";   778*M,  798*M,  "AIS", "med";   798*M,  803*M,  "NT", "low"
869*M,  884*M,  "(unassigned since Aug 2025)", "med"
935*M,  960*M,  "AIS|True", "med"
1452*M, 1472*M, "True", "high"
1805*M, 1880*M, "AIS|True", "med"
2110*M, 2140*M, "True", "low";   2140*M, 2170*M, "AIS", "med"
2300*M, 2370*M, "True", "high"
2500*M, 2600*M, "AIS", "high";   2600*M, 2690*M, "True", "med"
3300*M, 3800*M, "(not assigned in 2026)", "med"
703*M,  723*M,  "True", "med";   723*M,  743*M,  "AIS", "med";   743*M,  748*M,  "NT", "low"
890*M,  915*M,  "AIS|True", "med"
1710*M, 1785*M, "AIS|True", "med"
1920*M, 1950*M, "True", "low";   1950*M, 1980*M, "AIS", "med"
};
r = reshape(r.', 4, []).';
H = cell2table(r, 'VariableNames', {'flo', 'fhi', 'holder', 'conf'});
H.holder = string(H.holder); H.conf = string(H.conf);
end

function [key, name, duplex] = bandOf(fc)
B = bandTable();
i = find(B.flo <= fc & B.fhi >= fc, 1);
if isempty(i), key = ""; name = "outside the Thai IMT bands"; duplex = "";
else, key = B.key(i); name = B.name(i); duplex = B.duplex(i); end
end

function [h, conf] = operatorOf(f1, f2, fc)
H = holderTable();
if ~(f2 > f1), f1 = fc; f2 = fc; end
ov = max(0, min(H.fhi, f2) - max(H.flo, f1));
if all(ov == 0), ov = double(H.flo <= fc & H.fhi >= fc); end
if all(ov == 0), h = "?"; conf = ""; return, end
i = find(ov >= 0.2*max(ov));
h = strjoin(unique(H.holder(i), 'stable'), " | ");
conf = H.conf(i(1));
if numel(i) > 1, conf = "low"; h = h + " (spans blocks)"; end
end

function n = earfcn(fc, bandKey)
% DL EARFCN (TS 36.101 Table 5.7.3-1): N = N_Offs + 10 (F - F_low) / 1 MHz, 100 kHz raster.
t = struct('B1', [2110 0], 'B3', [1805 1200], 'B5', [869 2400], 'B8', [925 3450], 'B28', [758 9210], ...
    'L1500', [1452 9920], 'n40', [2300 38650], 'n41', [2496 39650]);
n = "?";
k = char(bandKey);
if isempty(k) || ~isfield(t, k), return, end
v = t.(k);
n = string(round(v(2) + 10*(fc/1e6 - v(1))));
end

function n = nrarfcn(fc)
% NR-ARFCN (TS 38.104 5.4.2.1): 5 kHz global raster below 3 GHz, 15 kHz from 3 to 24.25 GHz.
if fc < 3e9, n = string(round(fc/5e3)); else, n = string(round(600000 + (fc - 3e9)/15e3)); end
end

function [g, fs] = gscnNearest(fc)
% Nearest synchronization raster point (TS 38.104 Table 5.4.3.1-1).
if fc < 3e9
    N = round((fc - [50e3 150e3 250e3])/1200e3);
    cand = N*1200e3 + [50e3 150e3 250e3]; [~, i] = min(abs(cand - fc));
    Mv = [1 3 5]; fs = cand(i); g = 3*N(i) + (Mv(i) - 3)/2;
else
    N = round((fc - 3e9)/1.44e6); fs = 3e9 + N*1.44e6; g = 7499 + N;
end
end

%% ---- technology decision ------------------------------------------------------------

function c = decideTech(c)
S = c.snap; ev = strings(1, 0);
z = c.zsr; zk = "";
if ~isempty(z), zk = z.kind; end
dup = c.duplexBand;
pick = @(S, tech) S(S.tech == tech, :);
best = table();
if ~isempty(S), best = S(S.tech ~= "NR-SSB", :); end
ssbLike = ~isempty(S) && S.tech(1) == "NR-SSB" && S.err(1) < 0.12;
isSSB = ~isempty(z) && z.ssb.detected;
if isSSB || (c.bursty && ssbLike) || (zk == "bursty" && ssbLike)
    c.tech = "NR"; c.duplex = dup; if dup == "" || dup == "SDL", c.duplex = "FDD"; end
    if ssbLike, c.scs = S.scs(1); end
    c.techLabel = "NR " + c.duplex + " (idle: SSB only)";
    c.chBW = NaN; c.nrb = NaN; c.measBW = NaN;
    if isSSB, ev(end + 1) = sprintf("SSB bursts every %.0f ms", z.ssb.period*1e3); end
    if ssbLike, ev(end + 1) = sprintf("width %s = 240 subcarriers x %d kHz", fmtHz(c.occ), c.scs); end
    if c.bursty, ev(end + 1) = "seen only in the max over the sweeps"; end
    c.confidence = conf3(numel(ev) >= 2, true);
    c.evidence = ev; return
end
if isempty(best)
    c.techLabel = "?"; c.evidence = "no bandwidth estimate";
    if zk == "tdd"                         % pattern known, width not measured
        techs = unique(z.bestTechs); c.duplex = "TDD";
        if isscalar(techs), c.tech = techs; else, c.tech = strjoin(techs, "/"); end
        c.techLabel = c.tech + " TDD (BW ?)"; c.pattern = z.best; c.patternScore = z.score;
        c.evidence = ["zero-span pattern " + z.best + sprintf(" (score %.2f)", z.score), "bandwidth not measured"];
        c.confidence = conf3(false, isscalar(techs) && z.conf ~= "low");
    elseif zk == "cont"
        c.techLabel = "? (continuous: FDD DL or loaded carrier, BW ?)"; c.duplex = dup;
    end
    return
end
if zk == "tdd"
    techs = unique(z.bestTechs);
    cand = best(ismember(best.tech, techs), :);
    if isempty(cand), cand = best; end
    tie = cand.err <= cand.err(1) + 0.005 & cand.tech == "NR" & cand.scs == 30;
    if dup == "TDD" && any(tie), cand = [cand(find(tie, 1), :); cand(setdiff(1:height(cand), find(tie, 1)), :)]; end
    if numel(techs) > 1
        ev(end + 1) = "envelope identical for " + strjoin(techs, "/") + ": technology from the spectrum width";
    else
        ev(end + 1) = "zero-span pattern matches " + techs + " (" + z.best + ", score " + sprintf("%.2f", z.score) + ")";
    end
    c = setCand(c, cand(1, :));
    c.duplex = "TDD";
    if dup ~= "TDD", ev(end + 1) = "TDD envelope in a " + dup + " band (unexpected)"; end
    ev(end + 1) = sprintf("width %s fits %s %g MHz (%.1f %% off)", fmtHz(c.occ), c.tech, c.chBW, 100*cand.err(1));
    c.confidence = conf3(dup == "TDD" && z.conf ~= "low" && cand.err(1) < 0.03, cand.err(1) < 0.06);
    c.pattern = z.best; c.patternScore = z.score;
elseif zk == "cont"
    L = pick(best, "LTE"); Nn = pick(best, "NR");
    if c.occ > 20.5e6 || isempty(L) || (~isempty(Nn) && Nn.err(1) + 0.01 < L.err(1) && L.err(1) > 0.03)
        if isempty(Nn), Nn = best; end
        c = setCand(c, Nn(1, :));
        ev(end + 1) = "continuous envelope (loaded NR or LTE); width fits NR better";
    else
        c = setCand(c, L(1, :));
        ev(end + 1) = "continuous envelope: LTE CRS is always on";
    end
    c.duplex = dup; if dup == "TDD", c.duplex = "TDD?"; ev(end + 1) = "continuous envelope in a TDD band (UL load or unsynchronised)"; end
    ev(end + 1) = sprintf("width %s fits %s %g MHz (%.1f %% off)", fmtHz(c.occ), c.tech, c.chBW, 100*c.snapErr);
    c.confidence = conf3(any(dup == ["FDD" "SDL"]) && c.snapErr < 0.03, c.snapErr < 0.06);
else
    cand = best;
    if dup == "TDD"                        % prefer NR 30 kHz in TDD bands unless LTE fits clearly better
        Nn = cand(cand.tech == "NR" & cand.scs == 30, :);
        if ~isempty(Nn) && Nn.err(1) <= cand.err(1) + 0.01, cand = Nn; end
    elseif any(dup == ["FDD" "SDL"])       % prefer LTE in FDD bands on a tie
        L = cand(cand.tech == "LTE", :);
        if ~isempty(L) && L.err(1) <= cand.err(1) + 0.005, cand = L; end
    end
    c = setCand(c, cand(1, :));
    c.duplex = dup;
    if zk == "weak", ev(end + 1) = "too weak for the zero span"; elseif zk ~= "", ev(end + 1) = "zero span: " + z.class; end
    ev(end + 1) = sprintf("width %s fits %s %g MHz (%.1f %% off); band %s", fmtHz(c.occ), c.tech, c.chBW, 100*c.snapErr, c.band);
    c.confidence = conf3(false, c.snapErr < 0.03);
end
if c.snapErr > 0.06, ev(end + 1) = sprintf("bandwidth uncertain (%.0f %% from the nearest 3GPP value)", 100*c.snapErr); c.confidence = "low"; end
if c.truncated, c.confidence = "low"; end
if c.occSrc == "given"                     % bandwidth from 'Carriers', not from the spectrum
    ev = [ev(~startsWith(ev, "width ")), sprintf("channel bandwidth %g MHz given (not measured)", c.chBW)];
    if c.confidence == "high", c.confidence = "medium"; end
end
c.techLabel = c.tech + " " + c.duplex;
c.evidence = ev;
end

function c = setCand(c, r)
c.tech = r.tech(1); c.chBW = r.chBW(1); c.scs = r.scs(1); c.nrb = r.nrb(1); c.measBW = r.occ(1);
c.snapErr = r.err(1);
end

function s = conf3(high, med)
if high, s = "high"; elseif med, s = "medium"; else, s = "low"; end
end

%% ---- TDD templates ---------------------------------------------------------------------

function TT = tddTemplates(opt)
% DL occupancy per OFDM symbol over one period. NR: resource grids of nrWaveformGenerator;
% LTE: lteDuplexingInfo. Cached in support\nr5g\tdd_templates.mat; "analytic" = tables only.
persistent mem
useTB = opt.Templates == "auto";
key = sprintf('v1_%d_%d', useTB, opt.LTESSC);
if ~isempty(mem) && isfield(mem, key), TT = mem.(key); return, end
d = fullfile(fileparts(mfilename('fullpath')), 'support', 'nr5g');
file = fullfile(d, 'tdd_templates.mat');
if useTB && isfile(file)
    try
        S = load(file);
        if isfield(S, key), TT = S.(key); mem.(key) = TT; return, end
    catch
    end
end
if useTB, fprintf('Building TDD templates with the 5G and LTE Toolboxes (once; cached in support\\nr5g) ...\n'); end
nr = {"DDDSU", [10 2 2], 30; "DDDSU", [6 4 4], 30; "DDSUU", [10 2 2], 30; "DDDSUDDSUU", [10 2 2], 30; ...
    "DDDDDDDSUU", [6 4 4], 30; "DDDSU", [10 2 2], 15};
TT = [];
for k = 1:size(nr, 1)
    TT = [TT, nrTemplate(nr{k, 1}, nr{k, 2}, nr{k, 3}, useTB)]; %#ok<AGROW>
end
for cfgNo = 0:6
    TT = [TT, lteTemplate(cfgNo, opt.LTESSC, useTB)]; %#ok<AGROW>
end
TT = groupTemplates(TT);
mem.(key) = TT;
if useTB
    try
        if ~isfolder(d), mkdir(d); end
        S = struct(); if isfile(file), S = load(file); end
        S.(key) = TT; S.info = sprintf('ff_5g_analyzer TDD templates, MATLAB %s, %s', version, string(datetime('now')));
        save(file, '-struct', 'S');
    catch err
        warning('ff_5g_analyzer:cache', 'Template cache not written: %s', err.message);
    end
end
end

function T = tmplStruct(name, tech, scs, pat, sSym, slotDur, edges, occ, src)
T = struct('name', string(name), 'tech', string(tech), 'scs', scs, 'pattern', string(pat), 'sSym', sSym, ...
    'slotDur', slotDur, 'period', edges(end), 'edges', edges, 'occ', occ, 'source', string(src), ...
    'duty', sum(diff(edges).*occ)/edges(end), 'group', 0, 'groupName', "");
end

function T = nrTemplate(pat, sSym, scs, useTB)
pat = char(pat); per = numel(pat); slotDur = 1e-3*15/scs;
name = sprintf("NR%d %s (S %d:%d:%d)", scs, pat, sSym);
occ = []; src = "analytic (TS 38.213 slot formats)";
if useTB
    try
        cfg = nrTDDConfig(20, scs, pat, sSym, 20, true);
        [~, info] = nrWaveformGenerator(cfg);
        G = info.ResourceGrids(1).ResourceGridInCarrier;
        e = sum(abs(G(:, :, 1)).^2, 1) > 0;
        nSym = per*14; nPer = floor(numel(e)/nSym);
        occ = mean(reshape(double(e(1:nPer*nSym)), nSym, nPer), 2).';
        oi = nrOFDMInfo(cfg.SCSCarriers{1}.NSizeGrid, scs);
        L = double(oi.SymbolLengths(:).')/oi.SampleRate;
        Ls = repmat(L, 1, ceil(per*slotDur/1e-3) + 1); Ls = Ls(1:nSym);
        edges = [0 cumsum(Ls)];
        src = "nrWaveformGenerator resource grid";
    catch err
        warning('ff_5g_analyzer:template', 'NR template %s from the toolbox failed (%s): analytic map used.', name, err.message);
        occ = [];
    end
end
if isempty(occ)
    occ = zeros(1, per*14);
    for s = 1:per
        switch pat(s)
            case 'D', occ((s - 1)*14 + (1:14)) = 1;
            case 'S', occ((s - 1)*14 + (1:sSym(1))) = 1;
        end
    end
    if scs == 30, L = [1112 1096*ones(1, 13)]/61.44e6*2; else, L = [2208 2192*ones(1, 6) 2208 2192*ones(1, 6)]/30.72e6; end
    edges = [0 cumsum(repmat(L, 1, per))];
end
T = tmplStruct(name, "NR", scs, pat, sSym, slotDur, edges, occ, src);
end

function T = lteTemplate(cfgNo, ssc, useTB)
pats = ["DSUUUDSUUU" "DSUUDDSUUD" "DSUDDDSUDD" "DSUUUDDDDD" "DSUUDDDDDD" "DSUDDDDDDD" "DSUUUDSUUD"];
dw = [3 9 10 11 12 3 9 10 11 6];            % DwPTS symbols, normal CP (TS 36.211 Table 4.2-1)
typ = char(pats(cfgNo + 1)); nDL = zeros(1, 10); src = "analytic (TS 36.211)";
done = false;
if useTB
    try
        enb = struct('NDLRB', 6, 'CyclicPrefix', 'Normal', 'DuplexMode', 'TDD', 'TDDConfig', cfgNo, 'SSC', ssc, 'NSubframe', 0);
        for sf = 0:9
            enb.NSubframe = sf; di = lteDuplexingInfo(enb);
            nDL(sf + 1) = di.NSymbolsDL; typ(sf + 1) = di.SubframeType(1);
        end
        done = true; src = "lteDuplexingInfo";
    catch err
        warning('ff_5g_analyzer:template', 'lteDuplexingInfo failed (%s): analytic map used.', err.message);
    end
end
if ~done
    for sf = 1:10
        switch typ(sf), case 'D', nDL(sf) = 14; case 'S', nDL(sf) = dw(ssc + 1); otherwise, nDL(sf) = 0; end
    end
end
nsf = 10;
if isequal(typ(1:5), typ(6:10)) && isequal(nDL(1:5), nDL(6:10)), nsf = 5; end
occ = zeros(1, 14*nsf);
for sf = 1:nsf, occ((sf - 1)*14 + (1:nDL(sf))) = 1; end
L = [2208 2192*ones(1, 6) 2208 2192*ones(1, 6)]/30.72e6;
edges = [0 cumsum(repmat(L, 1, nsf))];
name = sprintf("LTE cfg%d %s (SSC %d)", cfgNo, pats(cfgNo + 1), ssc);
T = tmplStruct(name, "LTE", 15, typ(1:nsf), [dw(ssc + 1) NaN NaN], 1e-3, edges, occ, src);
end

function TT = groupTemplates(TT)
% Templates with the same period and an envelope correlation >= 0.995 cannot be told apart.
n = numel(TT); g = 1:n;
for i = 1:n
    for j = i + 1:n
        if abs(TT(i).period - TT(j).period) > 2e-6 || g(j) ~= j, continue, end
        nb = round(TT(i).period/2e-6);
        a = templateBins(TT(i), nb); b = templateBins(TT(j), nb);
        a = a - mean(a); b = b - mean(b);
        r = max(real(ifft(fft(a).*conj(fft(b)))))/(norm(a)*norm(b));
        if r >= 0.995, g(j) = g(i); end
    end
end
[~, ~, gi] = unique(g, 'stable');
for k = 1:n, TT(k).group = gi(k); end
for k = 1:n, TT(k).groupName = strjoin([TT(gi == gi(k)).name], " = "); end
end

function m = templateBins(T, nb, w)
% Occupancy averaged over nb phase bins of one period (16 sub-samples per bin), optionally
% smoothed circularly over w bins (the same smoothing as the measured envelope).
K = 16;
ph = ((0:nb*K - 1) + 0.5)/(nb*K)*T.period;
idx = discretize(ph, T.edges); idx(isnan(idx)) = numel(T.occ);
v = T.occ(idx);
m = mean(reshape(v, K, nb), 1).';
if nargin > 2 && w > 1
    o = mod((0:nb - 1).' + (-floor((w - 1)/2):ceil((w - 1)/2)), nb) + 1;
    m = mean(m(o), 2);
end
end

function C = templateCum(T, t)
% Integral of the DL occupancy from 0 to t (periodic), for the demo time gating.
P = T.period; d = diff(T.edges).*T.occ; cs = [0 cumsum(d)];
k = floor(t/P); r = t - k*P;
i = discretize(r, T.edges); i(isnan(i)) = numel(T.occ);
C = k*cs(end) + cs(i) + (r - T.edges(i)).*T.occ(i);
end

function T = ssbTemplate(scs, L, period)
% SSB burst occupancy: Case A (15 kHz) symbols {2,8} + 14n, Case C (30 kHz) symbols {2,8} + 14n,
% each SSB 4 symbols, L blocks, one burst per 'period'.
sym = 1e-3*15/scs/14;
n = 0:(L/2 - 1); st = reshape([2; 8] + 14*n, 1, []);
on = false(1, round(period/sym));
for s = st, on(s + 1:s + 4) = true; end
edges = (0:numel(on))*sym;
T = tmplStruct(sprintf("SSB %d kHz L%d %g ms", scs, L, period*1e3), "NR", scs, "", [NaN NaN NaN], sym*14, ...
    edges, double(on), "TS 38.213 4.1");
end

function cfg = nrTDDConfig(bwMHz, scs, pat, sSym, nSub, loaded)
% nrDLCarrierConfig with a TDD pattern: full-band PDSCH in the D slots, the first sSym(1)
% symbols of the S slot, nothing in the U slots (R2026a has no TDD object for waveforms).
pat = char(pat);
cfg = nrDLCarrierConfig;
cfg.FrequencyRange = 'FR1'; cfg.ChannelBandwidth = bwMHz; cfg.NumSubframes = nSub;
nrb = cfg.FR1BandwidthTable{sprintf('%dkHz', scs), sprintf('%dMHz', bwMHz)};
cfg.SCSCarriers{1}.SubcarrierSpacing = scs; cfg.SCSCarriers{1}.NSizeGrid = nrb;
cfg.BandwidthParts{1}.SubcarrierSpacing = scs; cfg.BandwidthParts{1}.NSizeBWP = nrb; cfg.BandwidthParts{1}.NStartBWP = 0;
if scs == 30
    cfg.SSBurst.BlockPattern = 'Case C'; cfg.SSBurst.TransmittedBlocks = ones(1, 8);
else
    cfg.SSBurst.BlockPattern = 'Case A'; cfg.SSBurst.TransmittedBlocks = [1 1 1 1];
end
cfg.SSBurst.Period = 20; cfg.SSBurst.SubcarrierSpacingCommon = scs;
cfg.CSIRS{1}.Enable = false;
per = numel(pat); dS = find(pat == 'D') - 1; sS = find(pat == 'S') - 1;
if loaded
    pd = {nrWavegenPDSCHConfig('Label', 'PDSCH_D', 'PRBSet', 0:nrb - 1, 'SymbolAllocation', [0 14], ...
        'SlotAllocation', dS, 'Period', per, 'Modulation', '64QAM')};
    if ~isempty(sS) && sSym(1) > 0
        pd{end + 1} = nrWavegenPDSCHConfig('Label', 'PDSCH_S', 'PRBSet', 0:nrb - 1, 'SymbolAllocation', [0 sSym(1)], ...
            'SlotAllocation', sS, 'Period', per, 'Modulation', '64QAM');
    end
    cfg.PDSCH = pd;
    cfg.PDCCH{1}.SlotAllocation = sort([dS sS]); cfg.PDCCH{1}.Period = per;
    cfg.SearchSpaces{1}.NumCandidates = [8 8 4 2 1];
else
    cfg.PDSCH{1}.Enable = false; cfg.PDCCH{1}.Enable = false;
end
end

function T = templateTable(TT)
if isempty(TT), T = table(); return, end
T = table([TT.name].', [TT.tech].', [TT.period].'*1e3, 100*[TT.duty].', [TT.group].', [TT.source].', ...
    'VariableNames', {'Template', 'Tech', 'Period_ms', 'DLDuty_pct', 'Group', 'Source'});
end

%% ---- zero span ------------------------------------------------------------------------

function r = tddAnalyze(t, p, slotDur)
% Classify a zero-span envelope: TDD (periodic on/off), continuous (FDD DL), bursty (SSB-only),
% intermittent. 100 us linear-power smoothing kills the OFDM speckle.
if nargin < 3, slotDur = 0.5e-3; end
p = p(:); t = t(:);
r = struct('lo', NaN, 'hi', NaN, 'pk', NaN, 'bursty', false, 'thr', NaN, 'onFrac', NaN, 'period', NaN, ...
    'periodicity', 0, 'medianOnRun', 0, 'slotPattern', "", 'slotOnFrac', [], 'dyn', 0, 'class', "too few points");
if numel(p) < 20 || any(~isfinite(p)), return, end
dt = t(2) - t(1); T = t(end) - t(1) + dt;
w = max(1, round(0.1e-3/dt));
ps = 10*log10(movmean(10.^(p/10), w));
lo = prctile(ps, 10); hi = prctile(ps, 90); pk99 = prctile(ps, 99.5);
r.lo = lo; r.hi = hi; r.pk = pk99;
r.bursty = (hi - lo) < 6 && (pk99 - lo) > 10;       % short bursts (duty < ~10 %) above a flat floor
if r.bursty, r.thr = lo + (pk99 - lo)/2; else, r.thr = (lo + hi)/2; end
on = ps > r.thr; r.onFrac = mean(on);
if r.bursty, x = double(on) - mean(on); else, x = ps - mean(ps); end
maxLag = min(numel(x) - 1, round(min(25e-3, 0.6*T)/dt));
[ac, lags] = xcorr(x, maxLag, 'coeff'); keep = lags*dt >= 0.9e-3; ac = ac(keep); lags = lags(keep);
if numel(ac) > 3 && all(isfinite(ac))
    [pk, loc] = findpeaks(ac, 'MinPeakProminence', 0.05);
    if ~isempty(pk)
        [pmax, im] = max(pk); L = lags(loc(im));
        for k = 6:-1:2                             % prefer the smallest strong sub-multiple
            j = find(abs(lags(loc) - L/k) <= max(2, 0.03*L/k));
            if ~isempty(j) && max(pk(j)) >= 0.85*pmax, L = lags(loc(j(1))); break; end
        end
        r.period = L*dt; r.periodicity = pmax;
    end
end
d = diff([0; on; 0]); st = find(d == 1); en = find(d == -1) - 1;
if ~isempty(st), r.medianOnRun = median(en - st + 1)*dt; end
if ~isnan(r.period) && ~r.bursty && (hi - lo) > 6
    d0 = diff([1; on; 1]); so = find(d0 == -1); eo = find(d0 == 1) - 1;
    if ~isempty(so)
        [~, ii] = max(eo - so); k0 = eo(ii) + 1; nP = round(r.period/dt); k0 = mod(k0 - 1, nP) + 1;
        nS = max(1, round(r.period/slotDur)); str = blanks(nS); r.slotOnFrac = zeros(1, nS);
        for k = 1:nS                               % fold over every complete period
            idx = [];
            for m = 0:floor((numel(on) - k0 + 1)/nP) - 1
                idx = [idx, k0 + m*nP + round(((k - 1)*slotDur)/dt) + (0:max(0, round(slotDur/dt) - 1))]; %#ok<AGROW>
            end
            idx = idx(idx <= numel(on)); if isempty(idx), continue, end
            fr = mean(on(idx)); r.slotOnFrac(k) = fr;
            if fr > 0.9, str(k) = 'D'; elseif fr > 0.25, str(k) = 'S'; else, str(k) = 'U'; end
        end
        r.slotPattern = string(str);
    end
end
r.dyn = hi - lo;
if r.bursty && r.periodicity >= 0.3
    r.class = sprintf("bursty periodic, period %.1f ms, burst %.2f ms", r.period*1e3, r.medianOnRun*1e3);
elseif r.bursty
    r.class = "sporadic bursts";
elseif r.dyn < 6
    r.class = "continuous (FDD DL / always on)";
elseif r.periodicity >= 0.5 && r.onFrac >= 0.3
    r.class = sprintf("TDD periodic, period %.2f ms, DL duty %.0f %%", r.period*1e3, 100*r.onFrac);
else
    r.class = "intermittent (traffic-dependent)";
end
end

function S = matchTemplates(t, p, TT, smoothT)
% Fold the 0..1 normalised envelope on every template period and correlate (circularly, all
% shifts) with the template smoothed the same way. Score = Pearson correlation of the profiles.
t = t(:); p = p(:); dt = t(2) - t(1); Tt = t(end) - t(1) + dt;
w = max(1, round(smoothT/dt));
ps = 10*log10(movmean(10.^(p/10), w));
lo = prctile(ps, 10); hi = prctile(ps, 90);
o = min(max((ps - lo)/max(hi - lo, eps), 0), 1);
n = numel(TT);
S = struct('score', NaN(1, n), 'lag', zeros(1, n), 'nb', zeros(1, n), 'prof', {cell(1, n)}, 'tmpl', {cell(1, n)});
for c = 1:n
    Pc = TT(c).period;
    if Tt < 0.98*Pc, continue, end
    nb = max(8, round(Pc/dt));
    b = min(nb, floor(mod(t - t(1), Pc)/Pc*nb) + 1);
    prof = accumarray(b, o, [nb 1], @mean, NaN);
    ok = ~isnan(prof);
    if nnz(ok) < 0.8*nb, continue, end
    if any(~ok)
        xi = find(ok); prof(~ok) = interp1([xi - nb; xi; xi + nb], repmat(prof(ok), 3, 1), find(~ok));
    end
    m = templateBins(TT(c), nb, max(1, round(w*dt/(Pc/nb))));
    x = prof - mean(prof); y = m - mean(m);
    if norm(x) == 0 || norm(y) == 0, continue, end
    cc = real(ifft(fft(x).*conj(fft(y))))/(norm(x)*norm(y));
    [S.score(c), i] = max(cc);
    S.lag(c) = i - 1; S.nb(c) = nb; S.prof{c} = prof; S.tmpl{c} = m;
end
end

function Z = analyzeZero(zs, D, A, ~)
a = zs.acq; P = levels(a, D); [n, N] = size(P); T = a.settings.sweepTime;
TT = A.templates;
dt = T/(N - 1); t = (0:N - 1).'*dt;
fl5 = floorHzOf(a, D, A.G) + 10*log10(1.0645*a.settings.rbw);
rs = arrayfun(@(k) tddAnalyze(t, P(k, :).', 0.5e-3), 1:n);
Z = struct('kind', "", 'class', "", 'period', median([rs.period], 'omitnan'), ...
    'periodicity', median([rs.periodicity]), 'dlDuty', median([rs.onFrac]), 'dyn', median([rs.dyn]), ...
    'lo', median([rs.lo]), 'hi', median([rs.hi]), 'floor5', fl5, 'slotPattern', "", 'slotOnFrac', [], ...
    'ranking', table(), 'best', "", 'bestTechs', strings(1, 0), 'score', NaN, 'margin', NaN, 'conf', "low", ...
    'theoryDuty', NaN, 'ssb', struct('checked', false, 'detected', false, 'period', NaN, 'burst', NaN, 'level', NaN, 'note', ""), ...
    't', t, 'P', P(1, :).', 'overlay', [], 'fold', [], 'nSweeps', n, 'center', a.settings.center, ...
    'detector', a.settings.detector, 'trigger', a.settings.trigger, 'sweepTime', T);
burstyFrac = mean([rs.bursty]);
if Z.hi < fl5 + 6
    Z.kind = "weak";
    Z.class = sprintf("too weak in 5 MHz: %.1f dBm vs floor %.1f dBm", Z.hi, fl5);
elseif burstyFrac >= 0.5
    Z.kind = "bursty";
    Z.class = rs(find([rs.bursty], 1)).class;
elseif Z.dyn < 6
    Z.kind = "cont";
    Z.class = sprintf("continuous (spread %.1f dB): FDD DL / always on", Z.dyn);
elseif Z.periodicity >= 0.5 && Z.dlDuty >= 0.3
    Z.kind = "tdd";
    Z.class = sprintf("TDD periodic, period %.2f ms, DL duty %.0f %%", Z.period*1e3, 100*Z.dlDuty);
else
    Z.kind = "intermittent";
    Z.class = "intermittent (traffic-dependent)";
end
if Z.kind == "tdd" && ~isempty(TT)
    Sc = NaN(n, numel(TT)); Ss = cell(1, n);
    for k = 1:n
        Ss{k} = matchTemplates(t, P(k, :).', TT, 50e-6);
        Sc(k, :) = Ss{k}.score;
    end
    sc = mean(Sc, 1, 'omitnan');
    grp = [TT.group];
    ug = unique(grp(isfinite(sc)), 'stable');
    gs = arrayfun(@(g) max(sc(grp == g)), ug);
    [gs, o] = sort(gs, 'descend'); ug = ug(o);
    gn = arrayfun(@(g) TT(find(grp == g, 1)).groupName, ug);
    gt = arrayfun(@(g) strjoin(unique([TT(grp == g).tech]), "/"), ug);
    gp = arrayfun(@(g) TT(find(grp == g, 1)).period*1e3, ug);
    gd = arrayfun(@(g) 100*TT(find(grp == g, 1)).duty, ug);
    Z.ranking = table(gn(:), gt(:), gp(:), gd(:), gs(:), 'VariableNames', {'Templates', 'Tech', 'Period_ms', 'DLDuty_pct', 'Score'});
    if ~isempty(ug)
        Z.best = gn(1); Z.score = gs(1); Z.bestTechs = unique([TT(grp == ug(1)).tech]);
        if numel(gs) > 1, Z.margin = gs(1) - gs(2); else, Z.margin = 1; end
        if Z.score >= 0.9 && Z.margin >= 0.05, Z.conf = "high";
        elseif Z.score >= 0.75 && Z.margin >= 0.02, Z.conf = "medium"; end
        mem = find(grp == ug(1)); [~, im] = max(sc(mem)); bt = mem(im);
        Z.theoryDuty = TT(bt).duty;
        [~, kb] = max(Sc(:, bt));
        r2 = tddAnalyze(t, P(kb, :).', TT(bt).slotDur);
        Z.slotPattern = r2.slotPattern; Z.slotOnFrac = r2.slotOnFrac;
        Z.P = P(kb, :).';
        S = Ss{kb}; nb = S.nb(bt);
        mS = circshift(S.tmpl{bt}, S.lag(bt));
        b = min(nb, floor(mod(t, TT(bt).period)/TT(bt).period*nb) + 1);
        Z.overlay = struct('name', TT(bt).name, 'y', mS(b));
        alt = [];
        if numel(ug) > 1
            ma = find(grp == ug(2)); [~, ia] = max(sc(ma)); alt = ma(ia);
        end
        Z.fold = struct('phase', ((0:nb - 1).' + 0.5)/nb*TT(bt).period, 'prof', S.prof{bt}, 'tmpl', mS, ...
            'name', TT(bt).name, 'alt', [], 'altName', "");
        if ~isempty(alt) && S.nb(alt) > 0
            Z.fold.alt = struct('phase', ((0:S.nb(alt) - 1).' + 0.5)/S.nb(alt)*TT(alt).period, ...
                'prof', S.prof{alt}, 'tmpl', circshift(S.tmpl{alt}, S.lag(alt)));
            Z.fold.altName = TT(alt).name;
        end
    end
end
if ~isempty(zs.ssb)
    Z.ssb = analyzeSSB(zs.ssb, D);
elseif Z.kind == "bursty"
    Z.ssb.note = "bursts seen in the main zero span; no SSB-length sweep";
end
end

function s = analyzeSSB(a, D)
P = levels(a, D); [n, N] = size(P); T = a.settings.sweepTime; dt = T/(N - 1); t = (0:N - 1).'*dt;
s = struct('checked', true, 'detected', false, 'period', NaN, 'burst', NaN, 'level', NaN, 'note', "", 't', t, 'P', P(1, :).');
allowed = [5 10 20 40 80 160]*1e-3;
for k = 1:n
    r = tddAnalyze(t, P(k, :).', 0.5e-3);
    if r.bursty && r.periodicity >= 0.3 && any(abs(r.period - allowed) <= 0.08*allowed)
        s.detected = true; s.period = r.period; s.burst = r.medianOnRun; s.level = r.pk; s.P = P(k, :).';
        break
    end
end
if s.detected
    s.note = sprintf("SSB bursts every %.1f ms, burst %.2f ms, peak %.1f dBm (peak detector)", s.period*1e3, s.burst*1e3, s.level);
else
    s.note = "no periodic SSB bursts (loaded cell, LTE, or too weak)";
end
end

%% ---- ACLR -------------------------------------------------------------------------------

function X = aclrCompute(c, a, D, fl, C, opt)
% ACLR = main channel power / adjacent channel power, both integrated over the transmission
% bandwidth (square filter, TS 38.104 6.6.3 / TS 36.104 6.6.2); bins summed x binWidth / ENBW.
f = a.f(:); lin = mean(10.^(levels(a, D)/10), 1).'; rbw = a.settings.rbw; enbw = 1.0645*rbw;
df = (f(end) - f(1))/(numel(f) - 1);
mb = c.measBW; ch = c.chBW*1e6; offs = [-2 -1 1 2]*ch;
cover = @(x) x - mb/2 >= f(1) - df/2 && x + mb/2 <= f(end) + df/2;
cp = @(x) sum(lin(abs(f - x) <= mb/2))*df/enbw;
Nb = 10^(fl/10)*mb;
X = struct('offsets', offs, 'value', NaN(1, 4), 'lower', false(1, 4), 'flag', repmat("", 1, 4), ...
    'pMain', NaN, 'limitNoise', NaN, 'measBW', mb, 'chBW', ch, 'src', a.settings.label, 'rbw', rbw, ...
    'status', "", 'aclr1', NaN, 'note', "", 'limit', 45);
if ~cover(c.fc), X.status = "main channel outside the sweep"; return, end
Pm = cp(c.fc); PmC = Pm;
if opt.NoiseCorrection, PmC = max(Pm - Nb, 0.01*Pm); end
X.pMain = 10*log10(PmC); X.limitNoise = 10*log10(Pm/Nb);
for k = 1:4
    x = c.fc + offs(k);
    if ~cover(x), X.flag(k) = "outside span"; continue, end
    Pa = cp(x);
    if Pa < 1.26*Nb
        X.value(k) = 10*log10(PmC/(0.26*Nb)); X.lower(k) = true; X.flag(k) = "noise-limited";
    else
        PaC = Pa; if opt.NoiseCorrection, PaC = Pa - Nb; end
        X.value(k) = 10*log10(PmC/PaC);
    end
    occ = "";
    for j = 1:numel(C)
        if C(j).fc == c.fc || ~(C(j).occ > 0), continue, end
        ov = max(0, min(C(j).f2, x + mb/2) - max(C(j).f1, x - mb/2));
        if ov >= min(0.2*mb, 0.5*C(j).occ), occ = "occupied by carrier at " + fmtHz6(C(j).fc); break, end
    end
    if strlength(occ) > 0, X.flag(k) = occ;
    elseif ~X.lower(k) && X.value(k) < 25, X.flag(k) = "strong adjacent signal (another emitter?)";
    end
end
valid = isfinite(X.value) & ~startsWith(X.flag, ["occupied" "outside" "strong"]);
v1 = valid & abs(offs) == ch;
if any(v1), X.aclr1 = min(X.value(v1)); end
if ~any(valid)
    X.status = "not measurable (adjacent channels occupied / outside the sweep)";
else
    [v, i] = min(X.value(valid)); lw = X.lower(valid); lw = lw(i);
    if lw && v >= 45, X.status = sprintf(">= %.1f dB: consistent with the 45 dB limit (INFO)", v);
    elseif lw, X.status = sprintf("noise-limited: >= %.1f dB (cannot verify 45 dB)", v);
    elseif abs(v - 45) < 1, X.status = sprintf("%.1f dB: at the 45 dB limit within +-1 dB (INFO, radiated)", v);
    elseif v >= 45, X.status = sprintf("%.1f dB >= 45 dB (INFO, radiated)", v);
    else, X.status = sprintf("%.1f dB < 45 dB (INFO, radiated: check conducted)", v);
    end
end
X.note = sprintf("measurement BW %s (square), RBW %s, noise %.1f dBm in the BW, max measurable %.1f dB", ...
    fmtHz(mb), fmtHz(rbw), 10*log10(Nb), X.limitNoise);
end

%% ---- reference spectra -------------------------------------------------------------------

function F = refFit(c, a, D, fl)
% Toolbox reference PSD through the same RBW / bucket, scaled to the measured channel power;
% residual rms over the occupied band for the chosen technology and the best alternative.
hyp = c.snap(c.snap.tech ~= "NR-SSB" & c.snap.err < 0.08, :);
if isempty(hyp), F = []; return, end
pickRows = hyp(hyp.tech == c.tech & hyp.chBW == c.chBW & hyp.scs == c.scs, :);
other = hyp(hyp.tech ~= c.tech, :);
rows = [pickRows(1:min(end, 1), :); other(1:min(end, 1), :)];
if isempty(rows), rows = hyp(1, :); end
f = a.f(:); lin = mean(10.^(levels(a, D)/10), 1).'; rbw = a.settings.rbw; enbw = 1.0645*rbw;
df = (f(end) - f(1))/(numel(f) - 1);
m = abs(f - c.fc) <= c.occ/2 + 3*rbw;
meas = 10*log10(lin);
F = struct('hyp', strings(1, 0), 'rms', [], 'model', {{}}, 'f', f, 'best', "", 'note', "");
for k = 1:height(rows)
    r = rows(k, :);
    try
        [fr, psd, src] = refPSD(r.tech, r.chBW, r.scs, r.nrb, c.duplex);
    catch err
        F.note = strjoin([F.note, sprintf("%s %g MHz reference failed: %s", r.tech, r.chBW, err.message)], "; ");
        continue
    end
    dfr = fr(2) - fr(1);
    K = ceil(2.5*rbw/dfr); g = exp(-log(2)*(((-K:K)*dfr)/(rbw/2)).^2).';
    pr = conv(psd(:), g, 'same')*dfr;
    nb = max(1, round(df/dfr)); if nb > 1, pr = movmean(pr, nb); end
    mdl = interp1(c.fc + fr(:), pr, f, 'linear', 0);
    inb = abs(f - c.fc) <= r.occ/2;
    sc = (10^(c.chP/10))/max(sum(mdl(inb))*df/enbw, realmin);
    mdl = mdl*sc + 10^(fl/10)*enbw;
    md = 10*log10(mdl);
    F.hyp(end + 1) = sprintf("%s %g MHz (%d kHz) [%s]", r.tech, r.chBW, r.scs, src);
    F.rms(end + 1) = sqrt(mean((meas(m) - md(m)).^2));
    F.model{end + 1} = md;
end
if isempty(F.rms), if strlength(F.note) == 0, F = []; end, return, end
[~, i] = min(F.rms); F.best = F.hyp(i);
if numel(F.rms) > 1
    d = abs(diff(F.rms(1:2)));
    if d < 0.3, s = sprintf("reference residuals differ by %.2f dB: not decisive", d);
    else, s = sprintf("best reference fit: %s (rms %.2f dB vs %.2f dB)", F.best, min(F.rms), max(F.rms)); end
else
    s = sprintf("reference fit rms %.2f dB", F.rms);
end
F.note = strjoin([s, F.note(strlength(F.note) > 0)], "; ");
end

function [fr, psd, src] = refPSD(tech, chBW, scs, nrb, duplex)
% Fully loaded reference waveform PSD (centred, mW/Hz up to a scale). Cached per configuration.
persistent cache
if isempty(cache), cache = containers.Map(); end
key = sprintf('%s_%g_%d_%d_%s', tech, chBW, scs, nrb, duplex);
if isKey(cache, key), v = cache(key); fr = v{1}; psd = v{2}; src = v{3}; return, end
if tech == "LTE"
    bw = sprintf('%gMHz', chBW);
    if startsWith(duplex, "TDD"), tm = lteTestModel('1.1', bw, 1, 'TDD'); else, tm = lteTestModel('1.1', bw); end
    [x, ~, info] = lteTestModelTool(tm); fs = info.SamplingRate; src = "lteTestModelTool E-TM1.1";
else
    cfg = nrDLCarrierConfig; cfg.FrequencyRange = 'FR1'; cfg.ChannelBandwidth = chBW; cfg.NumSubframes = 2;
    cfg.SCSCarriers{1}.SubcarrierSpacing = scs; cfg.SCSCarriers{1}.NSizeGrid = nrb;
    cfg.BandwidthParts{1}.SubcarrierSpacing = scs; cfg.BandwidthParts{1}.NSizeBWP = nrb; cfg.BandwidthParts{1}.NStartBWP = 0;
    if scs == 30, cfg.SSBurst.BlockPattern = 'Case C'; else, cfg.SSBurst.BlockPattern = 'Case A'; end
    cfg.SSBurst.SubcarrierSpacingCommon = scs;
    cfg.PDSCH{1}.PRBSet = 0:nrb - 1; cfg.PDSCH{1}.Modulation = '256QAM'; cfg.CSIRS{1}.Enable = false;
    cfg.PDCCH{1}.Enable = false;           % the default CORESET does not fit narrow 30 kHz grids
    [x, info] = nrWaveformGenerator(cfg); fs = info.ResourceGrids(1).Info.SampleRate; src = "nrWaveformGenerator";
end
x = x(:, 1);
res = min(chBW*1e6/400, 50e3);
nfft = min(2^16, 2^nextpow2(fs/res));
nfft = min(nfft, 2^floor(log2(numel(x)/4)));
[psd, fr] = pwelch(x, hann(nfft), nfft/2, nfft, fs, 'centered');
cache(key) = {fr, psd, src};
end

%% ---- coexistence ------------------------------------------------------------------------

function T = coexistence(C, D, opt)
T = table();
if strlength(opt.EMIFile) == 0, return, end
try
    if lower(opt.EMIFile) == "demo"
        E = demoEMIRecording(opt);
        Eon = E.on; Eoff = E.off; srcNote = "synthetic demo EMI (inverter-like hump at 2.4-2.7 GHz)";
        if opt.Save
            d = fullfile(opt.OutDir, "rec"); if ~isfolder(d), mkdir(d); end
            s = Eon; save(fullfile(d, "rec_emi_demo_on.mat"), '-struct', 's');
            s = Eoff; save(fullfile(d, "rec_emi_demo_off.mat"), '-struct', 's');
        end
    else
        Eon = load(opt.EMIFile); srcNote = "EMI recording " + opt.EMIFile;
        Eoff = [];
        if strlength(opt.EMIRefFile) > 0, Eoff = load(opt.EMIRefFile); srcNote = srcNote + " minus " + opt.EMIRefFile; end
    end
catch err
    warning('ff_5g_analyzer:emi', 'EMI recording not usable: %s', err.message);
    return
end
on = recordingAcq(Eon);
off = []; if ~isempty(Eoff), off = recordingAcq(Eoff); end
Dg = struct('extGain', D.extGain);
Pon = mean(10.^(levels(on, Dg)/10), 1); fOn = on.f; enOn = 1.0645*on.settings.rbw;
if ~isempty(off), Poff = mean(10.^(levels(off, Dg)/10), 1); enOff = 1.0645*off.settings.rbw; end
rows = {};
N0 = -174 + opt.UENoiseFigure;
for k = 1:numel(C)
    c = C(k);
    if ~(c.occ > 0) || ~isfinite(c.psd), continue, end
    m = fOn >= c.f1 & fOn <= c.f2;
    if nnz(m) < 3
        rows(end + 1, :) = {c.id, c.fc/1e6, c.techLabel, NaN, NaN, N0, NaN, NaN, NaN, NaN, NaN, NaN, false, "not covered by the EMI recording"}; %#ok<AGROW>
        continue
    end
    I = mean(Pon(m))/enOn; note = srcNote;
    Iref = 10^(floorModel(on.settings, Dg)/10);
    if ~isempty(off)
        mo = off.f >= c.f1 & off.f <= c.f2;
        if nnz(mo) >= 3
            Iref = mean(Poff(mo))/enOff;        % interferer OFF (ambient + analyzer floor)
            note = note + "; OFF recording subtracted";
        else
            note = note + "; reference does not cover the carrier: DANL model subtracted";
        end
    else
        note = note + "; analyzer floor (DANL model) subtracted";
    end
    I = I - Iref;
    resolved = I >= 0.122*Iref;                 % < 0.5 dB rise: not resolved, upper bound
    if ~resolved
        I = 0.122*Iref;
        note = note + "; EMI not resolved above the floor (< 0.5 dB rise): I is an upper bound";
    end
    if contains(upper(on.settings.detector), "POS"), note = note + "; PEAK detector: EMI overestimated"; end
    if ~(c.chBW > 0), note = note + "; idle cell (SSB only): S is the SSB-burst average"; end
    Idb = 10*log10(I) + opt.EMICoupling;
    duty = 1; if c.duplex == "TDD" && ~isempty(c.zsr) && c.zsr.dlDuty > 0.2, duty = c.zsr.dlDuty; end
    S = c.psd - 10*log10(duty);
    sinr0 = S - N0; sinr1 = S - 10*log10(10^(N0/10) + 10^(Idb/10));
    th0 = tr36942(sinr0, c.tech)*c.occ*duty/1e6; th1 = tr36942(sinr1, c.tech)*c.occ*duty/1e6;
    rows(end + 1, :) = {c.id, c.fc/1e6, c.techLabel, S, Idb, N0, Idb - N0, sinr0, sinr1, sinr0 - sinr1, th0, th1, resolved, note}; %#ok<AGROW>
end
if isempty(rows), return, end
T = cell2table(rows, 'VariableNames', {'ID', 'Fc_MHz', 'Technology', 'S_DLon_dBmHz', 'I_dBmHz', 'N_UE_dBmHz', ...
    'IoverN_dB', 'SINR0_dB', 'SINR_dB', 'SINRloss_dB', 'Thr0_Mbps', 'Thr_Mbps', 'Resolved', 'Note'});
end

function e = tr36942(sinr, tech)
% Truncated Shannon bound, 3GPP TR 36.942 Annex A.2: alpha 0.6, SINR_min -10 dB, max 4.4 bit/s/Hz
% (LTE 64QAM); 7.4 bit/s/Hz used for NR 256QAM.
mx = 4.4; if startsWith(tech, "NR"), mx = 7.4; end
e = 0.6*log2(1 + 10.^(sinr/10));
e(sinr < -10) = 0; e = min(e, mx);
end

%% ---- summary ------------------------------------------------------------------------------

function T = summaryTable(C)
n = numel(C);
if n == 0
    T = table(); return
end
col = @(f) arrayfun(f, C(:));
str = @(f) arrayfun(f, C(:));
zsCls = strings(n, 1); per = NaN(n, 1); duty = NaN(n, 1); ssb = strings(n, 1);
a1 = NaN(n, 4); ast = strings(n, 1); atxt = strings(n, 1);
for k = 1:n
    z = C(k).zsr;
    if ~isempty(z)
        zsCls(k) = z.kind;
        if any(z.kind == ["tdd" "bursty" "intermittent"]), per(k) = z.period*1e3; duty(k) = 100*z.dlDuty; end
        if z.kind == "cont", duty(k) = 100; end
        if z.ssb.checked, ssb(k) = string(onoff(z.ssb.detected)); end
        if z.ssb.detected, ssb(k) = sprintf("%.0f ms", z.ssb.period*1e3); end
    end
    X = C(k).aclr;
    if ~isempty(X)
        a1(k, :) = X.value; ast(k) = X.status;
        a1(k, startsWith(X.flag, ["occupied" "outside" "strong"])) = NaN;   % not an ACLR (text column keeps it)
        lab = ["L2" "L1" "U1" "U2"]; parts = strings(1, 0);
        for j = 1:4
            if strlength(X.flag(j)) > 0 && ~X.lower(j)
                parts(end + 1) = lab(j) + " " + X.flag(j); %#ok<AGROW>
                if isfinite(X.value(j)), parts(end) = parts(end) + sprintf(" (%.1f)", X.value(j)); end
            elseif X.lower(j), parts(end + 1) = sprintf("%s >=%.1f", lab(j), X.value(j)); %#ok<AGROW>
            elseif isfinite(X.value(j)), parts(end + 1) = sprintf("%s %.1f", lab(j), X.value(j)); %#ok<AGROW>
            end
        end
        atxt(k) = strjoin(parts, " | ");
    end
end
T = table(col(@(c) c.id), str(@(c) c.band), col(@(c) c.fc/1e6), str(@(c) c.holder), str(@(c) c.techLabel), ...
    col(@(c) c.chBW), col(@(c) c.scs), col(@(c) c.nrb), col(@(c) c.occ/1e6), col(@(c) c.chP), col(@(c) c.psd), ...
    col(@(c) c.epre), zsCls, per, duty, str(@(c) c.pattern), col(@(c) c.patternScore), ssb, a1(:, 2), a1(:, 3), ...
    atxt, ast, str(@(c) c.arfcn), str(@(c) c.confidence), str(@(c) strjoin([c.evidence c.notes], "; ")), ...
    'VariableNames', {'ID', 'Band', 'Fc_MHz', 'Operator', 'Technology', 'ChBW_MHz', 'SCS_kHz', 'NRB', 'OccBW_MHz', ...
    'Level_dBm', 'PSD_dBmHz', 'EPRE_dBm', 'ZeroSpan', 'Period_ms', 'DLDuty_pct', 'TDDPattern', 'PatternScore', ...
    'SSB', 'ACLR_L1_dB', 'ACLR_U1_dB', 'ACLR', 'ACLRStatus', 'ARFCN', 'Confidence', 'Evidence'});
end

function L = conclusions(C, D, A)
L = strings(0, 1);
if D.source == "demo", L(end + 1) = "SYNTHETIC DEMO DATA - not a measurement."; end
if isempty(C)
    L(end + 1) = "No carrier found: check the antenna on Port 2, 'Preamp', 'Threshold' and the bands.";
    if D.aborted, L(end + 1) = "The measurement was stopped early: results are partial."; end
    L = L(:);
    return
end
for k = 1:numel(C)
    c = C(k);
    s = sprintf("#%d %s: %s, %s", c.id, fmtHz6(c.fc), c.techLabel, ternary(c.chBW > 0, sprintf("%g MHz", c.chBW), "BW ?"));
    if strlength(c.pattern) > 0, s = s + sprintf(", pattern %s (score %.2f)", c.pattern, c.patternScore); end
    s = s + sprintf(", %.1f dBm, operator %s, confidence %s", c.chP, c.holder, c.confidence);
    if ~isempty(c.aclr), s = s + ", ACLR " + c.aclr.status; end
    L(end + 1) = s; %#ok<AGROW>
end
if ~isempty(A.coex) && height(A.coex) > 0
    T = A.coex; ok = T.Resolved & isfinite(T.SINRloss_dB);
    if any(ok)
        loss = T.SINRloss_dB; loss(~ok) = -Inf; [mx, i] = max(loss);
        L(end + 1) = sprintf("Coexistence: worst SINR loss %.1f dB on carrier #%d (I/N %.1f dB).", mx, T.ID(i), T.IoverN_dB(i));
    end
    u = ~T.Resolved & isfinite(T.SINRloss_dB);
    if any(u)
        L(end + 1) = sprintf(['Coexistence: EMI below the analyzer floor on carrier(s) %s: SINR loss <= %.1f dB ' ...
            '(the FieldFox floor is ~%.0f dB above a UE noise floor: use an LNA or a closer antenna to resolve it).'], ...
            strjoin("#" + string(T.ID(u).'), ", "), max(T.SINRloss_dB(u)), max(T.IoverN_dB(u)) - 10*log10(0.122));
    end
end
if D.aborted, L(end + 1) = "The measurement was stopped early: results are partial."; end
L = L(:);
end

function printSummary(R)
fprintf('\n=== ff_5g_analyzer: %d carrier(s) (%s) ===\n', height(R.carriers), R.source);
if height(R.carriers) > 0
    disp(R.carriers(:, {'ID', 'Fc_MHz', 'Operator', 'Technology', 'ChBW_MHz', 'Level_dBm', 'ZeroSpan', ...
        'TDDPattern', 'ACLR_L1_dB', 'ACLR_U1_dB', 'Confidence'}));
end
for s = R.conclusions(:).', fprintf('  %s\n', s); end
if ~isempty(R.coexistence), disp(R.coexistence(:, 1:13)); end
fprintf('Output: %s\n', R.outDir);
end

%% ===== figures ========================================================================

function files = makeFigures(D, A, opt)
files = strings(0, 1);
figs0 = findall(groot, 'Type', 'figure');
try
    if ~isempty(D.scan), files(end + 1, 1) = plotScans(D, A, opt); end
    for k = 1:numel(A.carriers)
        files(end + 1, 1) = plotCarrier(A.carriers(k), D, A, opt); %#ok<AGROW>
    end
    if ~isempty(A.coex) && height(A.coex) > 0, files(end + 1, 1) = plotCoex(A, opt); end
catch err
    warning('ff_5g_analyzer:plot', 'Figure failed: %s', err.message);
    if strcmpi(opt.Visible, 'off')          % the half-drawn hidden figure would stay in the session
        f = findall(groot, 'Type', 'figure'); delete(f(~ismember(f, figs0)));
    end
end
end

function file = plotScans(D, A, opt)
n = numel(D.scan);
fig = figure('Visible', opt.Visible, 'Position', [50 50 1200 min(260*n, 1300)], 'Color', 'w');
tl = tiledlayout(fig, n, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
H = holderTable();
for i = 1:n
    a = D.scan(i).acq; ax = nexttile(tl);
    f = a.f/1e6; lin = 10.^(levels(a, D)/10); en = 1.0645*a.settings.rbw;
    mHz = 10*log10(mean(lin, 1)/en); xHz = 10*log10(max(lin, [], 1)/en);
    fl = floorHzOf(a, D, A.G);
    plot(ax, f, xHz, 'Color', [0.75 0.75 0.85]); hold(ax, 'on');
    plot(ax, f, mHz, 'Color', [0 0.3 0.7]);
    yline(ax, fl, ':', 'floor'); yline(ax, fl + opt.Threshold, '--', 'threshold');
    yl = [fl - 8, max(max(xHz), fl + 30) + 8]; ylim(ax, yl);
    for h = find(H.fhi > a.f(1) & H.flo < a.f(end)).'
        x = [max(H.flo(h), a.f(1)) min(H.fhi(h), a.f(end))]/1e6;
        patch(ax, x([1 2 2 1]), yl(1) + [0 0 3 3], opColor(H.holder(h)), 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        text(ax, mean(x), yl(1) + 4.5, H.holder(h), 'HorizontalAlignment', 'center', 'FontSize', 7);
    end
    for c = A.carriers
        if c.fc < a.f(1) || c.fc > a.f(end) || ~(c.occ > 0), continue, end
        x = [c.f1 c.f2]/1e6; y = [fl + opt.Threshold, c.psd];
        patch(ax, x([1 2 2 1]), y([1 1 2 2]), [1 0.6 0.1], 'FaceAlpha', 0.15, 'EdgeColor', [0.9 0.4 0]);
        text(ax, mean(x), c.psd + 2, sprintf('#%d %s %s', c.id, c.tech, ternary(c.chBW > 0, sprintf('%g', c.chBW), '')), ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.7 0.2 0]);
    end
    grid(ax, 'on'); ylabel(ax, 'dBm/Hz');
    title(ax, sprintf('%s  (RBW %s, %d sweeps, average detector; grey = max)', D.scan(i).name, fmtHz(a.settings.rbw), size(a.P, 1)), 'FontWeight', 'normal');
end
xlabel(ax, 'MHz');
if D.source == "demo", title(tl, 'SYNTHETIC DEMO - band scans', 'Color', 'r'); else, title(tl, 'Band scans: ' + D.idn); end
file = saveFig(fig, fullfile(opt.OutDir, "fig_scans.png"), opt);
end

function file = plotCarrier(c, D, ~, opt)
fig = figure('Visible', opt.Visible, 'Position', [60 60 1200 760], 'Color', 'w');
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);
a = []; if c.aclrIdx > 0, a = D.aclr(c.aclrIdx).acq; elseif c.scanIdx > 0, a = D.scan(c.scanIdx).acq; end
if ~isempty(a)
    f = a.f(:); lin = 10.^(levels(a, D)/10);
    m = 10*log10(mean(lin, 1)); plot(ax, f/1e6, m, 'Color', [0 0.3 0.7]); hold(ax, 'on');
    yl = [min(m) - 3, max(m) + 6];
    X = c.aclr;
    if ~isempty(X)
        for o = [0 X.offsets]
            x = (c.fc + o + [-1 1]*X.measBW/2)/1e6;
            col = [0.2 0.7 0.2]; if o ~= 0, col = [0.6 0.6 0.6]; end
            patch(ax, x([1 2 2 1]), yl([1 1 2 2]), col, 'FaceAlpha', 0.08, 'EdgeColor', col);
        end
    end
    if ~isempty(c.ref)
        sty = {'--', ':'};
        for k = 1:numel(c.ref.model)
            plot(ax, f/1e6, c.ref.model{k}, sty{min(k, 2)}, 'LineWidth', 1.1);
        end
        nP = 0; if ~isempty(X), nP = 1 + numel(X.offsets); end
        legend(ax, ["measured" repmat("", 1, nP) c.ref.hyp], 'Location', 'south', 'FontSize', 7, ...
            'Interpreter', 'none', 'AutoUpdate', 'off');
    end
    ylim(ax, yl); grid(ax, 'on'); xlabel(ax, 'MHz'); ylabel(ax, sprintf('dBm in %s', fmtHz(a.settings.rbw)));
    tt = string(sprintf('#%d %s  %s, %s', c.id, fmtHz6(c.fc), c.techLabel, ternary(c.chBW > 0, sprintf('%g MHz', c.chBW), 'BW ?')));
    if ~isempty(X), tt = tt + newline + "ACLR " + X.status; end
    title(ax, tt, 'FontWeight', 'normal', 'Interpreter', 'none', 'FontSize', 9);
else
    text(ax, 0.1, 0.5, 'no spectrum'); axis(ax, 'off');
end
ax = nexttile(tl);
z = c.zsr;
if ~isempty(z)
    plot(ax, z.t*1e3, z.P, 'Color', [0.4 0.4 0.4]); hold(ax, 'on');
    w = max(1, round(0.1e-3/(z.t(2) - z.t(1))));
    plot(ax, z.t*1e3, 10*log10(movmean(10.^(z.P/10), w)), 'b', 'LineWidth', 1);
    yline(ax, z.floor5, ':', 'floor (5 MHz)');
    if ~isempty(z.overlay)
        plot(ax, z.t*1e3, z.lo + (z.hi - z.lo)*z.overlay.y, 'r', 'LineWidth', 1.2);
        legend(ax, ["envelope" "100 us smoothed" "" z.overlay.name], 'Location', 'south', 'FontSize', 7, ...
            'Interpreter', 'none', 'AutoUpdate', 'off');
    end
    ylim(ax, [min(z.floor5, min(z.P)) - 3, max(z.P) + 8]);
    grid(ax, 'on'); xlabel(ax, 'ms'); ylabel(ax, 'dBm in 5 MHz');
    title(ax, sprintf('Zero span %s (%s, %s): %s', fmtHz6(z.center), z.detector, z.trigger, z.class), 'FontWeight', 'normal', 'Interpreter', 'none');
else
    text(ax, 0.1, 0.5, 'no zero span'); axis(ax, 'off');
end
ax = nexttile(tl);
if ~isempty(z) && ~isempty(z.fold)
    F = z.fold;
    plot(ax, F.phase*1e3, F.prof, 'k', 'LineWidth', 1.2); hold(ax, 'on');
    plot(ax, F.phase*1e3, F.tmpl, 'r', 'LineWidth', 1.1);
    leg = ["measured (folded)" F.name];
    if ~isempty(F.alt), plot(ax, F.alt.phase*1e3, F.alt.tmpl, 'Color', [0.2 0.6 0.9]); leg(end + 1) = "runner-up: " + F.altName; end
    legend(ax, leg, 'Location', 'south', 'FontSize', 7, 'Interpreter', 'none');
    ylim(ax, [-0.1 1.25]); grid(ax, 'on'); xlabel(ax, 'ms in period'); ylabel(ax, 'DL occupancy');
    title(ax, sprintf('Best: %s, score %.2f, margin %.2f (%s); slots %s', z.best, z.score, z.margin, z.conf, z.slotPattern), ...
        'FontWeight', 'normal', 'Interpreter', 'none', 'FontSize', 8);
elseif ~isempty(c.snap)
    S = c.snap(1:min(end, 6), :);
    txt = compose("%s %g MHz %d kHz %d RB: %.2f MHz (%.1f %%)", S.tech, S.chBW, S.scs, S.nrb, S.occ/1e6, 100*S.err);
    text(ax, 0.02, 0.95, ["Bandwidth candidates (occupied " + fmtHz(c.occ) + "):"; txt], 'VerticalAlignment', 'top', 'FontSize', 8, 'Interpreter', 'none');
    axis(ax, 'off');
else
    axis(ax, 'off');
end
ax = nexttile(tl);
if ~isempty(z) && z.ssb.checked
    plot(ax, z.ssb.t*1e3, z.ssb.P, 'Color', [0.5 0.2 0.6]); grid(ax, 'on');
    xlabel(ax, 'ms'); ylabel(ax, 'dBm (peak detector)');
    title(ax, "SSB check: " + z.ssb.note, 'FontWeight', 'normal', 'FontSize', 8, 'Interpreter', 'none');
else
    lines = ["Operator guess: " + c.holder + " (" + c.holderConf + ")"; "Band: " + c.bandName; ...
        "Occupied BW (edge fit): " + fmtHz(c.occ) + ternary(isfinite(c.occCI), sprintf(' +- %s', fmtHz(c.occCI)), ''); ...
        "OBW 99 %: " + fmtHz(c.obw99); sprintf("Channel power %.1f dBm, PSD %.1f dBm/Hz, EPRE %.1f dBm", c.chP, c.psd, c.epre); ...
        c.arfcn + "  " + c.gscn; "Confidence: " + c.confidence; c.evidence(:)];
    if ~isempty(c.ref), lines(end + 1) = c.ref.note; end
    text(ax, 0.02, 0.95, lines, 'VerticalAlignment', 'top', 'FontSize', 8, 'Interpreter', 'none'); axis(ax, 'off');
end
if D.source == "demo", title(tl, 'SYNTHETIC DEMO', 'Color', 'r'); end
file = saveFig(fig, fullfile(opt.OutDir, sprintf("fig_carrier_%02d.png", c.id)), opt);
end

function file = plotCoex(A, opt)
T = A.coex(isfinite(A.coex.SINR0_dB), :);
fig = figure('Visible', opt.Visible, 'Position', [70 70 1000 420], 'Color', 'w');
ax = axes(fig);
lab = compose("#%d %.0f MHz", T.ID, T.Fc_MHz); lab(~T.Resolved) = lab(~T.Resolved) + " (EMI < floor: bound)";
bar(ax, categorical(lab, lab), [T.SINR0_dB T.SINR_dB]);
grid(ax, 'on'); ylabel(ax, 'UE SINR (dB)');
legend(ax, ["without EMI" "with EMI"], 'Location', 'northeast');
title(ax, sprintf('Coexistence: UE noise figure %g dB, EMI coupling %+g dB', opt.UENoiseFigure, opt.EMICoupling), 'FontWeight', 'normal');
file = saveFig(fig, fullfile(opt.OutDir, "fig_coexistence.png"), opt);
end

function c = opColor(h)
if contains(h, "|"), c = [0.6 0.5 0.2];
elseif startsWith(h, "AIS"), c = [0.2 0.65 0.2];
elseif startsWith(h, "True"), c = [0.85 0.15 0.15];
elseif startsWith(h, "NT"), c = [0.15 0.35 0.85];
else, c = [0.6 0.6 0.6];
end
end

function file = saveFig(fig, file, opt)
exportgraphics(fig, file, 'Resolution', 120);
if strcmpi(opt.Visible, 'off'), close(fig); end
end

%% ===== saving and report ===============================================================

function files = saveRun(D, A, opt)
files = strings(0, 1);
try
    carriers = A.table;
    f = fullfile(opt.OutDir, "nr5g_run.mat");
    save(f, "D", "carriers"); files(end + 1, 1) = f;
    if height(A.table) > 0
        f = fullfile(opt.OutDir, "nr5g_carriers.csv"); writetable(A.table, f); files(end + 1, 1) = f;
    end
    d = fullfile(opt.OutDir, "rec"); if ~isfolder(d), mkdir(d); end
    recs = {};
    for i = 1:numel(D.scan), recs(end + 1, :) = {"rec_scan_" + D.scan(i).key, D.scan(i).acq}; end %#ok<AGROW>
    for i = 1:numel(D.aclr), recs(end + 1, :) = {sprintf("rec_aclr_%.0fMHz", D.aclr(i).fc/1e6), D.aclr(i).acq}; end %#ok<AGROW>
    for i = 1:numel(D.zs)
        recs(end + 1, :) = {sprintf("rec_zs_%.0fMHz", D.zs(i).fc/1e6), D.zs(i).acq}; %#ok<AGROW>
        if ~isempty(D.zs(i).ssb), recs(end + 1, :) = {sprintf("rec_ssb_%.0fMHz", D.zs(i).fc/1e6), D.zs(i).ssb}; end %#ok<AGROW>
    end
    for i = 1:size(recs, 1)
        a = recs{i, 2};
        s = struct('f', a.f, 'P', a.P, 't', a.t, 'settings', a.settings, 'idn', D.idn);
        s.settings.extGainTotal = D.extGain;
        s.settings.levelNote = "P as displayed; DUT-referred = P + extgInstr - extGainTotal";
        fn = fullfile(d, matlab.lang.makeValidName(recs{i, 1}) + ".mat");
        save(fn, '-struct', 's'); files(end + 1, 1) = fn; %#ok<AGROW>
    end
catch err
    warning('ff_5g_analyzer:save', 'Saving failed: %s', err.message);
end
end

function file = makeReport(D, A, opt, figs)
import mlreportgen.report.*
import mlreportgen.dom.*
fprintf('Writing the %s report ...\n', upper(opt.Report));
typ = char(opt.Report); if opt.Report == "html", typ = 'html-file'; end   % one self-contained file
rpt = Report(char(fullfile(opt.OutDir, "nr5g_report")), typ);
when = datetime('now', 'TimeZone', 'Asia/Bangkok', 'Format', 'yyyy-MM-dd HH:mm ''(Asia/Bangkok)''');
sub = "FieldFox measurement - " + D.idn;
if D.source == "demo", sub = "SYNTHETIC DEMO DATA - not a measurement"; end
if D.source == "recording", sub = "Re-analysis of a recording"; end
add(rpt, TitlePage('Title', '5G NR / LTE carrier analysis', 'Subtitle', char(sub), ...
    'Author', char(string(getenv('USERNAME'))), 'PubDate', char(when)));
add(rpt, TableOfContents);
ch = Chapter('Title', 'Summary');
if D.source == "demo"
    p = Paragraph('SYNTHETIC DEMO: analytic spectra and 5G/LTE Toolbox envelopes through a FieldFox model.'); p.Bold = true; p.Color = 'red'; add(ch, p);
end
add(ch, Paragraph(['Hypotheses from swept spectra and zero-span envelopes (no I/Q, nothing demodulated). ' ...
    'Operators are guesses from the Thai 2025 band plan. ACLR verdicts are INFO only (radiated).']));
if height(A.table) > 0
    add(ch, rptTable(A.table(:, {'ID', 'Fc_MHz', 'Operator', 'Technology', 'ChBW_MHz', 'Level_dBm', 'ZeroSpan', ...
        'TDDPattern', 'ACLR_L1_dB', 'ACLR_U1_dB', 'Confidence'}), 'Carriers', 40));
end
add(ch, UnorderedList(cellstr(A.conclusions)));
addImg(ch, "fig_scans.png", 'Band scans with the Thai band plan and the detected carriers');
add(rpt, ch);
for k = 1:numel(A.carriers)
    c = A.carriers(k);
    ch = Chapter('Title', char(sprintf("Carrier #%d: %s, %s", c.id, fmtHz6(c.fc), c.techLabel)));
    items = ["Operator guess: " + c.holder + " (" + c.holderConf + "), band " + c.bandName; ...
        sprintf("Occupied bandwidth (edge fit) %s, OBW 99 %% %s, channel power %.1f dBm, PSD %.1f dBm/Hz, EPRE %.1f dBm", ...
        fmtHz(c.occ), fmtHz(c.obw99), c.chP, c.psd, c.epre); ...
        "Technology " + c.techLabel + ", " + ternary(c.chBW > 0, sprintf("%g MHz, %d kHz, %d RB", c.chBW, c.scs, c.nrb), "bandwidth ?") + ", confidence " + c.confidence; ...
        c.arfcn + " " + c.gscn; "Evidence: " + strjoin(c.evidence, "; ")];
    if ~isempty(c.zsr), items(end + 1) = "Zero span: " + c.zsr.class + ternary(strlength(c.zsr.best) > 0, ...
            sprintf("; best template %s (score %.2f, margin %.2f, %s); slots %s", c.zsr.best, c.zsr.score, c.zsr.margin, c.zsr.conf, c.zsr.slotPattern), ""); end %#ok<AGROW>
    if ~isempty(c.zsr) && c.zsr.ssb.checked, items(end + 1) = "SSB: " + c.zsr.ssb.note; end %#ok<AGROW>
    if ~isempty(c.aclr), items(end + 1) = "ACLR: " + c.aclr.status + " (" + c.aclr.note + ")"; end %#ok<AGROW>
    if ~isempty(c.ref), items(end + 1) = "Reference: " + c.ref.note; end %#ok<AGROW>
    if all(isfinite(c.chp)), items(end + 1) = sprintf("FieldFox CHP cross-check: %.1f dBm (MATLAB %.1f dBm)", c.chp(1), c.chP); end %#ok<AGROW>
    if ~isempty(c.notes), items(end + 1) = "Notes: " + strjoin(c.notes, "; "); end %#ok<AGROW>
    add(ch, UnorderedList(cellstr(items)));
    if ~isempty(c.zsr) && height(c.zsr.ranking) > 0
        add(ch, rptTable(c.zsr.ranking, 'TDD template ranking', 8));
    end
    addImg(ch, sprintf("fig_carrier_%02d.png", c.id), 'Spectrum with channels and toolbox references; zero span with the best template; folded envelope; SSB check');
    add(rpt, ch);
end
if ~isempty(A.coex) && height(A.coex) > 0
    ch = Chapter('Title', 'Coexistence with the EMI recording');
    add(ch, Paragraph(sprintf(['UE at the place of the EMI recording: noise figure %g dB, EMI coupling %+g dB. ' ...
        'SINR = DL-on PSD / (kTF + EMI PSD); throughput by TR 36.942 A.2 (truncated Shannon, alpha 0.6).'], opt.UENoiseFigure, opt.EMICoupling)));
    add(ch, rptTable(A.coex(:, 1:13), 'SINR impact (Resolved = false: EMI below the analyzer floor, values are upper bounds)', 20));
    addImg(ch, "fig_coexistence.png", 'UE SINR with and without the interference');
    add(rpt, ch);
end
ch = Chapter('Title', 'Method and limits');
add(ch, UnorderedList({ ...
    'Spectra: average detector, sweeps power-averaged in MATLAB; floor = min(empirical 5th percentile, DANL model + 3 dB). Gaussian RBW, ENBW 1.0645 x RBW.', ...
    'Bandwidth: -3 dB edges refined by lsqcurvefit (RBW-convolved step); snapped to TS 38.104 Table 5.3.2-1 / TS 36.104 Table 5.6-1.', ...
    'TDD: zero span RBW 5 MHz; envelope folded on each template period and correlated with 5G Toolbox (nrWaveformGenerator grids) and LTE Toolbox (lteDuplexingInfo) templates.', ...
    'Identical envelopes (e.g. NR 30 kHz DDDDDDDSUU 6:4:4 = LTE cfg 2 SSC 7 = NR 15 kHz DDDSU) are grouped; the spectrum width decides.', ...
    'ACLR: TS 38.104 6.6.3 / TS 36.104 6.6.2 (45 dB), square measurement filter of the transmission bandwidth, noise-corrected; radiated values are INFO only.', ...
    'No I/Q: no PCI, MIB/SIB, PLMN or EVM. Operator = frequency-based guess (2025 NBTC auction; some block edges inferred).'}));
if ~isempty(D.notes), add(ch, UnorderedList(cellstr(D.notes))); end
if strlength(opt.Notes) > 0, add(ch, Paragraph("Notes: " + opt.Notes)); end
add(rpt, ch);
close(rpt);
file = string(rpt.OutputPath);
fprintf('Report: %s\n', file);

    function addImg(parent, name, cap)
        fpath = fullfile(opt.OutDir, name);
        if isfile(fpath) && any(endsWith(figs, name))
            add(parent, FormalImage('Image', char(fpath), 'Caption', cap, 'Width', '16cm'));
        end
    end
end

function bt = rptTable(T, ttl, maxRows)
import mlreportgen.report.*
import mlreportgen.dom.*
T = T(1:min(end, maxRows), :);
names = T.Properties.VariableNames;
body = cell(height(T), numel(names));
for j = 1:numel(names)
    x = T.(names{j});
    if islogical(x), c = string(x);
    elseif isnumeric(x)
        c = compose("%.4g", x);
        k = isfinite(x) & x == round(x) & abs(x) < 1e7; c(k) = compose("%d", x(k));
    else, c = string(x);
    end
    c(ismissing(c) | c == "NaN") = "-";
    body(:, j) = cellstr(c(:));
end
ft = FormalTable(replace(names, "_", " "), body);
ft.Style = {FontSize('7pt'), Border('solid'), ColSep('solid'), RowSep('solid'), Width('100%')};
ft.Header.Style = {Bold(true), BackgroundColor('#e8eef7')};
bt = BaseTable(ft);
bt.Title = ttl;
end

%% ===== helpers ==========================================================================

function [t, pdBm] = zsEmulate(x, fs, fOff, rbw, Tsweep, N, det, danlHz, addNoise)
% Swept analyser in zero span: Gaussian RBW (|H|^2 = 0.5 at +-rbw/2) at fOff, envelope power,
% N display points over Tsweep (first to last point) with a SAMP / AVER / POS detector.
x = x(:); L = numel(x);
need = ceil(Tsweep*fs*N/(N - 1)); if L < need, x = repmat(x, ceil(need/L), 1); end
x = x(1:need); n = (0:need - 1).';
if fOff ~= 0, x = x.*exp(-1j*2*pi*fOff*n/fs); end
F = ifftshift((-floor(need/2):ceil(need/2) - 1).'*fs/need);
y = ifft(fft(x).*exp(-log(2)/2*(2*F/rbw).^2));
if addNoise
    pn = 10^(danlHz/10)*rbw*sqrt(pi/(4*log(2)));
    y = y + sqrt(pn/2)*(randn(need, 1) + 1j*randn(need, 1));
end
pw = abs(y).^2;
e = round(linspace(1, need + 1, N + 1));
pd = zeros(N, 1);
for k = 1:N
    seg = pw(e(k):max(e(k), e(k + 1) - 1));
    switch det
        case "SAMP", pd(k) = seg(1);
        case "AVER", pd(k) = mean(seg);
        otherwise,   pd(k) = max(seg);
    end
end
pdBm = 10*log10(pd); t = (0:N - 1).'*Tsweep/(N - 1);
end

function v = pick13(x)
% Largest 1-3-10 step <= x (FieldFox RBW steps), at least 1 Hz.
d = 10.^floor(log10(max(x, 1)));
s = [1 3 10]*d;
v = s(find(s <= x*(1 + 1e-9), 1, 'last'));
if isempty(v), v = d; end
end

function s = onoff(b)
if b, s = 'on'; else, s = 'off'; end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function s = fmtHz(f)
if ~isfinite(f), s = '-'; return, end
a = abs(f);
if a >= 1e9, s = sprintf('%.4g GHz', f/1e9);
elseif a >= 1e6, s = sprintf('%.4g MHz', f/1e6);
elseif a >= 1e3, s = sprintf('%.4g kHz', f/1e3);
else, s = sprintf('%.4g Hz', f);
end
end

function s = fmtHz6(f)
if ~isfinite(f), s = '-'; return, end
s = sprintf('%.3f MHz', f/1e6);
end

function s = fmtRange(r)
s = sprintf('%s - %s', fmtHz(r(1)), fmtHz(r(end)));
end
