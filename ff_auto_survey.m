function R = ff_auto_survey(src, varargin)
%FF_AUTO_SURVEY  Autonomous spectrum + EMI survey and emitter inventory with a FieldFox.
%
%   R = ff_auto_survey("10.90.36.200")                 full survey 100 kHz - 6.5 GHz (~2-3 min)
%   R = ff_auto_survey(ip, 'Preset', "quick")           < 1 min;  "deep" = finer RBW, ~5-8 min
%   R = ff_auto_survey(ip, 'Range', [400e6 470e6])      part of the range only
%   R = ff_auto_survey(ip, 'Baseline', "...\survey.mat")  new / disappeared / changed emitters
%   R = ff_auto_survey(ip, 'Watch', true, 'Interval', 60, 'Duration', 3600)   survey, then monitor + alert
%   R = ff_auto_survey("demo")                          synthetic survey, no instrument
%   R = ff_auto_survey("demo", 'Watch', true, 'Duration', 5, 'Interval', 1)   synthetic watch
%   R = ff_auto_survey("...\survey.mat")                re-analyse a saved survey
%   R = ff_auto_survey("rec.mat") | (folder) | ["a.mat" "b.mat"]   any suite recording(s)
%   P = ff_auto_survey(ip, 'DryRun', true)              only the scan plan (P.plan) and time estimate
%
% What it does (live; the FieldFox only listens - no RF output):
%   1. Coarse scan of 'Range' in segments, RBW per segment (10 kHz below 30 MHz,
%      100 kHz to 1 GHz, 300 kHz to 3 GHz, 1 MHz above), 'Sweeps' sweeps with the
%      Peak (POS) and then the Average (AVER) detector, 'Repeats' times. Preamp
%      "auto": ON, but OFF again in a segment whose peak exceeds 'PreampMax' or that
%      reports ADC over-range (CALC:MEAS:WAOR?).
%   2. Noise floor per segment: 20th percentile in ~40 blocks, running minimum over
%      5 blocks (wide carriers do not lift it), then corrected to the median of the
%      noise-only points. Detection (CFAR-like): points more than max('Threshold',
%      'K' x noise spread) dB above the floor in the max-hold OR the averaged trace,
%      grown into regions, split at dips of 'SplitDip' dB. Broadband rises: floor
%      more than 'RiseDB' above the segment's own lowest floor (robust line).
%      Harmonic combs: power cepstrum (periodogram of the dB excess) + line fit.
%      >= 8 narrow detections packed together outside channelised services (IMT,
%      FM, DTV, Wi-Fi, GNSS) become one "cluster" (EMI hash / unresolved comb).
%   3. Every detection is measured (99 % OBW and -26 dB BW with the noise power
%      subtracted, -3 dB width / OBW "shoulder" (1 = flat top), spectral flatness,
%      edge slope, peak-minus-average, channel power, occupancy = % of sweeps in
%      which it is present). The 'MaxZoom' most informative ones are re-measured
%      with a span ~3x their width and RBW ~1.7 % of it (spectral lines, combs
%      inside a "hump"); the 'MaxZeroSpan' most interesting ones (EMI, unknown,
%      TDD bands) are captured in zero span: 20 ms AVER (TDD pattern: period, DDDSU
%      slots), 100 ms SAMP (bursts, duty, 50/100 Hz mains-synchronous, impulsive)
%      with a noise-calibrated contrast threshold, 100 ms POS (radar pulses).
%      A broadband rise is zoomed with a narrow span/fine RBW to look for the lines
%      (comb) of a switching converter inside it.
%   4. Classification = combined evidence, as a hypothesis with a confidence:
%        Licensed service  ff_signal_id technology (if installed) consistent with
%                          the Thailand band plan (support\survey\th_bandplan.csv),
%                          FM/DTV/IMT widths, TDD pattern, radar pulses
%        Satellite         support\survey\sat_channels.csv (frequency within Doppler
%                          + OBW/2, width), visibility and Doppler now from cached
%                          TLEs (support\survey\tle, and ff_cubesat_link's support\tle
%                          read-only; 'Online' fetches from CelesTrak, 2 h cache),
%                          ff_satsig("overhead") if present (its catalog, SGP4,
%                          Doppler-shifted f_rx), SatNOGS DB (Online, 24 h cache).
%                          A carrier in a satellite downlink band with no satellite
%                          overhead is flagged as possible interference.
%        EMI (+sub-type)   harmonic comb (spacing), broadband noise rise, mains-
%                          synchronous or impulsive bursts, any line in a GNSS band,
%                          a carrier in a satellite downlink band with no satellite
%                          overhead, possible analyzer intermodulation
%        Unknown           nothing fits well enough
%   5. Inventory table, band occupancy (per band-plan row: % of points x sweeps
%      above the threshold), figures (overview with band-plan overlay, segments,
%      waterfall, zoom thumbnails, zero span, occupancy), PDF/DOCX report.
%   6. 'Baseline': emitters matched by frequency; per-sweep levels compared with a
%      Wilcoxon rank-sum test (ranksum), Benjamini-Hochberg q -> New / Gone /
%      Changed / Same.
%   7. 'Watch': repeats the coarse scan every 'Interval' s for 'Duration' s and alerts
%      (console, watch_log.csv, optional beep and 'AlertFcn') when a pass exceeds a
%      mask built from the survey itself (max hold or floor + threshold, + 'MaskMargin'
%      dB): new or stronger emissions (many in one segment -> one summary alert), new
%      EMI (comb spacing, broadband rise), steady emitters that vanish, and anomalies:
%      per-point robust z against the earlier passes (median/MAD) and per-segment
%      features (robust z; robustcov Mahalanobis + chi2 once there is enough history,
%      or iforest with 'Anomaly' "iforest"). A change must persist 'Confirm' passes.
%
% Output R: plan, inventory (one row per emitter: Freq/Lo/Hi/OBW, levels, Excess,
% Occupancy, Class, SubType, Confidence, Band, Technology, Satellite, Evidence, Kind =
% signal | comb | rise | cluster), bandOccupancy, combs (coarse and zoom), rises,
% zoom, zeroSpan, baseline (table: Match New/Matched/In EMI region/Gone, Delta_dB,
% p, q, Status), watch (alerts, passes, known), conclusions, notes, segments
% (traces, floors), raw (all sweeps: re-analyse with ff_auto_survey(R.files.mat)),
% meta, files, outDir.
% Files in 'OutDir': survey.mat (variable "survey"), inventory.csv,
% band_occupancy.csv, fig_*.png, survey_report.pdf, rec\*.mat (fieldfox_live
% format: f, P, t, settings, idn - one file per segment and detector).
%
% Name-Value options (defaults):
%  Measurement
%   'Preset'        "full"   "quick" (< 1 min) | "full" (~2-3 min) | "deep" (~5-8 min)
%   'Range'         [100e3 6.5e9]  Hz
%   'Segments'      table()  own plan: Name, Start_Hz, Stop_Hz, RBW_Hz, Points
%   'Preamp'        "auto"   true | false | "auto" (see 1.)
%   'PreampMax'     -40      dBm peak above which "auto" switches the preamp off
%   'Detectors'     ["POS" "AVER"]
%   'Sweeps'        []       per detector and segment (quick 2, full 4, deep 8)
%   'Repeats'       []       whole coarse scan (quick 1, full 2, deep 4)
%   'RefLevel'      -20      dBm;  'Atten' "auto" or dB
%   'Zoom'          true;    'MaxZoom' [] (quick 4, full 12, deep 25); 'ZoomSweeps' 3; 'ZoomPoints' 1001
%   'ZeroSpan'      true;    'MaxZeroSpan' [] (quick 1, full 4, deep 8); 'ZeroSpanSweeps' 3
%   'Port'          5025
%  Analysis
%   'Threshold'     6        dB above the floor (minimum);  'K' 5 (x noise spread)
%   'SplitDip'      10       dB;   'RiseDB' 4 dB (broadband rise)
%   'SignalID'      "auto"   ff_signal_id method: "auto" (ensemble when its models are
%                            cached, else template) | "ensemble" | "template" | "ml" | "off"
%   'SatID'         true     satellite candidates (catalog, TLE visibility, Doppler)
%   'Online'        false    allow celestrak.org (TLE, 2 h cache) and db.satnogs.org
%                            (transmitters, 24 h cache); stops at the first error
%   'SatnogsFile'   ""       a SatNOGS DB transmitters JSON to use offline (else the 24 h cache)
%   'TleDir'        ""       TLE cache folder (default support\survey\tle)
%   'Location'      [13.7563 100.5018 10]  lat, lon (deg), alt (m) - Bangkok
%   'Time'          NaT      survey time for offline files without time stamps
%   'BandPlan'      ""       CSV like support\survey\th_bandplan.csv
%   'Baseline'      ""       another survey.mat / recording / folder
%   'Alpha'         0.05     false discovery rate for the baseline comparison
%  Watch
%   'Watch'         false;   'Interval' 60 s; 'Duration' 600 s; 'WatchSweeps' 1
%   'MaskMargin'    6        dB above the survey max hold before a pass point counts
%   'Confirm'       2        consecutive passes a change must persist before it alerts
%   'Beep'          false;   'AlertFcn' [] (called with a struct per alert)
%   'Anomaly'       "auto"   "auto" (robust z, robustcov when enough history) |
%                            "robustcov" | "iforest" | "zscore"
%  Output
%   'Report'        "pdf"    "pdf" | "docx" | "none";  'Notes' "" (text for the report)
%   'OutDir'        ""       default survey_<time> in the current folder
%   'Save'          true     false = nothing written (no files, no report)
%   'Visible'       'on'     'off' = no windows (figures are still saved)
%   'Prompt'        true     live: confirm the plan first (skipped in -batch)
%   'Demo'          false    same as src = "demo";  'Seed' 1
%   'DryRun'        false
%
% Read this before trusting the result:
%   - No I/Q (option 351 absent): nothing is demodulated. A class is a HYPOTHESIS
%     from the spectrum shape, band plan, timing and catalogs; the Evidence column
%     says why. Confidence is a heuristic score, not a calibrated probability.
%   - Levels are display levels in the RBW of each segment (not comparable between
%     segments for wide signals). Only a calibrated antenna factor gives field
%     strength; use ff_emi_study for EMI pre-compliance (CISPR limits).
%   - Most satellites are at or below the FieldFox noise (-153 dBm/Hz with preamp):
%     expect only ORBCOMM/Meteor/CubeSat beacons with a gain antenna + LNA, and only
%     during a pass. GNSS is always 17-25 dB below kTB: a visible line in a GNSS
%     band is interference.
%   - Below 30 MHz the floor is usually ambient man-made noise; a "broadband rise"
%     there needs a baseline (or the source switched off) to be called EMI.
%   - Swept peak detection sees bursty signals (Wi-Fi, TDD, radar) in only some
%     sweeps: occupancy and levels depend on sweep time; zero span shows timing.
%     Bursty EMI (mains-gated inverter combs) breaks into a different set of
%     fragments in every survey: compare EMI regions (rise / cluster / comb), not
%     single fragments; the baseline does exactly that ("In EMI region").
%   - "quick" uses 30 kHz RBW below 30 MHz: a comb with < ~45 kHz spacing is not
%     resolved there (it shows as a rise or cluster); the zoom on it resolves the lines.
%     Narrow humps (< ~5 MHz) in the 30-1000 MHz quick segment may be missed: "full".
%   - Watch compares single-sweep passes with the survey's max hold: weak changes
%     (< 'MaskMargin' above it) are not alerted; the per-point history needs 5 passes.
%   - Thailand band plan: compiled 2026-10-10 from NBTC/press, not a legal table.
%
% See also ff_signal_id, ff_emi_study, ff_cubesat_link, fieldfox_live, FieldFox.

if nargin < 1 || (isnumeric(src) && isempty(src)), src = "10.90.36.200"; end
opt = parseOptions(defaultOptions(), varargin);
opt = applyPreset(opt);
here = fileparts(mfilename('fullpath'));
opt.here = here;
opt.supp = fullfile(here, 'support', 'survey');
if strlength(string(opt.TleDir)) == 0, opt.TleDir = fullfile(opt.supp, 'tle'); end   % own cache; support\tle is read too
[kind, src] = resolveSource(src, opt);
plan = makePlan(opt);
if opt.DryRun
    R = struct('plan', plan, 'estimate_s', estimateTime(plan, opt), 'opt', opt);
    printPlan(plan, opt, R.estimate_s);
    return
end
BP = loadBandPlan(opt);
CAT = loadSatCatalog(opt);
if opt.Save, opt.OutDir = makeOutDir(opt); else, opt.OutDir = ""; end
notes = strings(0, 1);
io = [];
switch kind
    case "file"
        raw = loadRaw(src, opt);
        if opt.Watch
            notes(end + 1) = "Watch needs the instrument or the demo: ignored for a file.";
            warning('ff_auto_survey:watch', '%s', notes(end));
        end
    otherwise
        if kind == "live"
            printPlan(plan, opt, estimateTime(plan, opt));
            if opt.Prompt && isInteractive()
                input('Press Enter to start the survey (Ctrl+C to abort) ', 's');
            end
        end
        io = openIO(kind, src, opt);     % live: takeControl + guaranteed release
        raw = acquireSurvey(io, plan, opt, BP);
end
tA = tic;
R = analyzeSurvey(raw, opt, BP, CAT);
fprintf('  analysis %.1f s\n', toc(tA));
R.plan = plan;
R.notes = [R.notes; notes];
if strlength(string(opt.Baseline)) > 0
    try
        R.baseline = compareBaseline(R, opt, BP, CAT);
    catch err
        R.notes(end + 1) = "Baseline comparison failed: " + err.message;
        warning('ff_auto_survey:baseline', '%s', R.notes(end));
    end
end
if opt.Watch && kind ~= "file"
    R.watch = watchLoop(io, plan, opt, R, BP);
end
R.conclusions = conclusions(R, opt);
clear io                                  % releases the instrument now (live)
R.outDir = opt.OutDir;
if opt.Save
    R.files = saveOutputs(R, opt);
end
printSummary(R);
end

%% ===== options, plan, source ===============================================

function opt = defaultOptions()
opt = struct('Preset', "full", 'Range', [100e3 6.5e9], 'Segments', table(), 'Preamp', "auto", ...
    'PreampMax', -40, 'Detectors', ["POS" "AVER"], 'Sweeps', [], 'Repeats', [], 'RefLevel', -20, ...
    'Atten', "auto", 'Zoom', true, 'MaxZoom', [], 'ZoomSweeps', 3, 'ZoomPoints', 1001, ...
    'ZeroSpan', true, 'MaxZeroSpan', [], 'ZeroSpanSweeps', 3, 'Port', 5025, ...
    'Threshold', 6, 'K', 5, 'SplitDip', 10, 'RiseDB', 4, 'SignalID', "auto", 'SatID', true, ...
    'Online', false, 'Location', [13.7563 100.5018 10], 'Time', NaT, 'BandPlan', "", ...
    'Baseline', "", 'Alpha', 0.05, 'Watch', false, 'Interval', 60, 'Duration', 600, ...
    'WatchSweeps', 1, 'MaskMargin', 6, 'Confirm', 2, 'Beep', false, 'AlertFcn', [], 'Anomaly', "auto", 'Report', "pdf", ...
    'OutDir', "", 'Notes', "", 'Save', true, 'Visible', 'on', 'Prompt', true, 'Demo', false, ...
    'Seed', 1, 'DryRun', false, 'SatnogsFile', "", 'TleDir', "");
end

function opt = parseOptions(opt, args)
if mod(numel(args), 2), error('ff_auto_survey:args', 'Options must be Name, Value pairs.'); end
names = fieldnames(opt);
for k = 1:2:numel(args)
    i = find(strcmpi(names, string(args{k})), 1);
    if isempty(i), error('ff_auto_survey:option', 'Unknown option "%s".', string(args{k})); end
    opt.(names{i}) = args{k + 1};
end
end

function opt = applyPreset(opt)
opt.Preset = lower(string(opt.Preset));
switch opt.Preset
    case "quick", d = [2 1 4 1];
    case "full",  d = [4 2 12 4];
    case "deep",  d = [8 4 25 8];
    otherwise, error('ff_auto_survey:preset', 'Preset must be "quick", "full" or "deep".');
end
if isempty(opt.Sweeps), opt.Sweeps = d(1); end
if isempty(opt.Repeats), opt.Repeats = d(2); end
if isempty(opt.MaxZoom), opt.MaxZoom = d(3); end
if isempty(opt.MaxZeroSpan), opt.MaxZeroSpan = d(4); end
if ~opt.Zoom, opt.MaxZoom = 0; end
if ~opt.ZeroSpan, opt.MaxZeroSpan = 0; end
opt.Sweeps = max(1, round(opt.Sweeps)); opt.Repeats = max(1, round(opt.Repeats));
opt.WatchSweeps = max(1, round(opt.WatchSweeps)); opt.Confirm = max(1, round(opt.Confirm));
opt.MaskMargin = max(0, double(opt.MaskMargin));
opt.MaxZoom = max(0, round(opt.MaxZoom)); opt.MaxZeroSpan = max(0, round(opt.MaxZeroSpan));
if ~(opt.Interval > 0) || ~(opt.Duration > 0), error('ff_auto_survey:watch', 'Interval and Duration must be > 0 s.'); end
opt.Detectors = upper(string(opt.Detectors));
if any(~ismember(opt.Detectors, ["POS" "AVER" "SAMP" "NEG" "NORM"]))
    error('ff_auto_survey:det', 'Detectors must be from POS, AVER, SAMP, NEG, NORM.');
end
p = opt.Preamp;
if (isstring(p) || ischar(p)) && lower(string(p)) == "auto", opt.PreampMode = "auto";
elseif islogical(p) || isnumeric(p), opt.PreampMode = string(ternary(logical(p), "on", "off"));
else, error('ff_auto_survey:preamp', 'Preamp must be true, false or "auto".');
end
r = double(opt.Range);
if numel(r) ~= 2 || r(2) <= r(1), error('ff_auto_survey:range', 'Range must be [start stop] in Hz.'); end
if r(1) < 5e3 || r(2) > 6.5e9
    warning('ff_auto_survey:range', 'Range clipped to the SA range 5 kHz - 6.5 GHz.');
end
opt.Range = [max(r(1), 5e3), min(r(2), 6.5e9)];
if islogical(opt.Visible) || isnumeric(opt.Visible), opt.Visible = ternary(logical(opt.Visible), 'on', 'off'); end
opt.Visible = char(lower(string(opt.Visible)));
if ~ismember(opt.Visible, {'on', 'off'}), error('ff_auto_survey:visible', 'Visible must be ''on'' or ''off''.'); end
opt.Report = lower(string(opt.Report));
opt.SignalID = lower(string(opt.SignalID));
opt.Anomaly = lower(string(opt.Anomaly));
opt.Location = double(opt.Location(:).');
if numel(opt.Location) == 2, opt.Location(3) = 0; end
end

function [kind, src] = resolveSource(src, opt)
if opt.Demo, kind = "demo"; return, end
if isa(src, 'FieldFox'), kind = "live"; return, end
s = string(src);
if isscalar(s) && lower(s) == "demo", kind = "demo"; return, end
if any(endsWith(lower(s), ".mat")) || (isscalar(s) && isfolder(s))
    kind = "file"; src = s; return
end
if isscalar(s) && ~isempty(regexp(s, '^[\w\.\-:]+$', 'once')), kind = "live"; src = s; return, end
error('ff_auto_survey:source', 'Give an IP address, a FieldFox object, "demo", a .mat file or a folder.');
end

function plan = makePlan(opt)
if istable(opt.Segments) && height(opt.Segments) > 0
    plan = opt.Segments;
    need = ["Name" "Start_Hz" "Stop_Hz" "RBW_Hz"];
    if ~all(ismember(need, string(plan.Properties.VariableNames)))
        error('ff_auto_survey:segments', 'Segments needs the columns Name, Start_Hz, Stop_Hz, RBW_Hz [, Points].');
    end
    if ~ismember("Points", string(plan.Properties.VariableNames))
        plan.Points = min(10001, max(401, 2*round((plan.Stop_Hz - plan.Start_Hz)./plan.RBW_Hz) + 1));
    end
    plan.Name = string(plan.Name);
    plan = plan(:, {'Name', 'Start_Hz', 'Stop_Hz', 'RBW_Hz', 'Points'});
    return
end
switch opt.Preset
    case "quick"
        C = {"HF 0.1-30 MHz", 100e3, 30e6, 30e3, 2001; "VHF/UHF 30-1000 MHz", 30e6, 1e9, 300e3, 6001; ...
             "L/S 1-3 GHz", 1e9, 3e9, 1e6, 4001; "C 3-6.5 GHz", 3e9, 6.5e9, 3e6, 2001};
    case "full"
        C = {"HF 0.1-30 MHz", 100e3, 30e6, 10e3, 6001; "VHF 30-520 MHz", 30e6, 520e6, 100e3, 10001; ...
             "UHF 520-1000 MHz", 520e6, 1e9, 100e3, 10001; "L 1-2 GHz", 1e9, 2e9, 300e3, 10001; ...
             "S 2-3 GHz", 2e9, 3e9, 300e3, 10001; "C 3-6.5 GHz", 3e9, 6.5e9, 1e6, 7001};
    otherwise
        C = {"HF 0.1-30 MHz", 100e3, 30e6, 3e3, 10001; "VHF 30-300 MHz", 30e6, 300e6, 30e3, 10001; ...
             "UHF 300-1000 MHz", 300e6, 1e9, 100e3, 10001; "L 1-2 GHz", 1e9, 2e9, 100e3, 10001; ...
             "S 2-3 GHz", 2e9, 3e9, 100e3, 10001; "C 3-4.5 GHz", 3e9, 4.5e9, 300e3, 10001; ...
             "C 4.5-6.5 GHz", 4.5e9, 6.5e9, 300e3, 10001};
end
C{1, 2} = min(C{1, 2}, opt.Range(1));                 % allow a start below 100 kHz
rows = {};
for i = 1:size(C, 1)
    lo = max(C{i, 2}, opt.Range(1)); hi = min(C{i, 3}, opt.Range(2));
    if hi <= lo, continue, end
    frac = (hi - lo)/(C{i, 3} - C{i, 2});
    pts = min(10001, max(401, round(C{i, 5}*frac)));
    nm = C{i, 1};
    if frac < 0.999, nm = sprintf("%s-%s", fmtHz(lo), fmtHz(hi)); end
    rows(end + 1, :) = {string(nm), lo, hi, C{i, 4}, pts}; %#ok<AGROW>
end
plan = cell2table(rows, 'VariableNames', {'Name', 'Start_Hz', 'Stop_Hz', 'RBW_Hz', 'Points'});
end

function t = estimateTime(plan, opt)
tsw = 0.04 + 6e-5*plan.Points;                         % s per sweep incl. transfer (measured on this unit)
t = opt.Repeats*sum(numel(opt.Detectors)*opt.Sweeps*tsw + 0.4) + ...
    opt.MaxZoom*(opt.ZoomSweeps*0.5 + 0.6) + opt.MaxZeroSpan*1.5*(opt.ZeroSpanSweeps*0.25 + 0.6) + 3;
end

function printPlan(plan, opt, est)
fprintf('\nff_auto_survey plan (%s): %d segment(s), %d sweep(s) x %s, %d repeat(s), preamp %s\n', ...
    opt.Preset, height(plan), opt.Sweeps, strjoin(opt.Detectors, '+'), opt.Repeats, opt.PreampMode);
for i = 1:height(plan)
    fprintf('  %-22s %10s - %-10s RBW %-8s %5d points\n', plan.Name(i), fmtHz(plan.Start_Hz(i)), ...
        fmtHz(plan.Stop_Hz(i)), fmtHz(plan.RBW_Hz(i)), plan.Points(i));
end
fprintf('  zoom on up to %d, zero span on up to %d detections; estimated time %.0f s\n\n', opt.MaxZoom, opt.MaxZeroSpan, est);
end

function tf = isInteractive()
tf = usejava('desktop');
try, tf = tf && ~batchStartupOptionUsed; catch, end
end

function d = makeOutDir(opt)
d = string(opt.OutDir);
if strlength(d) == 0
    d = fullfile(pwd, "survey_" + string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
end
if ~isfolder(d), mkdir(d); end
end

%% ===== instrument / demo I/O ===============================================

function io = openIO(kind, src, opt)
io = struct('kind', kind, 'ff', [], 'idn', "", 'options', "", 'state', containers.Map(), ...
    'sa', [], 'zs', [], 'zsBegin', @() [], 'zsEnd', @() [], 'now', @() datetime('now'), 'cleanup', []);
io.state('pass') = 0;
if kind == "demo"
    rng(opt.Seed);
    [W, t0, dnote] = demoWorld(opt);
    io.idn = "SYNTHETIC DEMO (ff_auto_survey)";
    off = t0 - datetime('now');
    io.now = @() datetime('now') + off;
    st = io.state;
    io.sa = @(cfg, det, n) demoSA(W, cfg, det, n, st);
    io.zs = @(cfg, n) demoZS(W, cfg, n, st);
    io.demoNote = dnote;
    return
end
own = ~isa(src, 'FieldFox');
if own, ff = FieldFox(string(src), opt.Port); else, ff = src; end
ff.takeControl();
io.cleanup = onCleanup(@() releaseFF(ff));
io.ff = ff; io.idn = ff.IDN; io.options = strjoin(ff.Options, ',');
prepareSA(ff);
canWaor = probe(ff, "CALC:MEAS:WAOR?");
io.sa = @(cfg, det, n) liveSA(ff, cfg, det, n, canWaor, opt);
io.zs = @(cfg, n) liveZS(ff, cfg, n);
tmp = "ffsurvey_zs.sta";
io.zsBegin = @() ff.saveState(tmp);
io.zsEnd = @() zsRestore(ff, tmp);
end

function releaseFF(ff)
try
    ff.release();
catch err
    warning('ff_auto_survey:release', 'Could not restore the instrument: %s', err.message);
end
end

function prepareSA(ff)
if ff.mode() ~= "SA", ff.setMode("SA"); end
u = tryQuery(ff, "AMPL:UNIT?");
if strlength(u) > 0 && ~startsWith(upper(u), "DBM"), ff.write("AMPL:UNIT DBM"); end
for c = ["ANT" "CABL"]
    if ismember(upper(tryQuery(ff, "AMPL:CORR:" + c + "?")), ["1" "ON"]), ff.write("AMPL:CORR:" + c + " 0"); end
end
ff.write("TRIG:SOUR FREE");
ff.write("SENS:BAND:VID:AUTO 1");
ff.opc(30);
e = ff.errors();
if ~isempty(e), fprintf('  (setup: instrument reported %s - ignored)\n', strjoin(e, '; ')); end
end

function A = liveSA(ff, cfg, det, n, canWaor, opt)
args = {'Start', cfg.Start, 'Stop', cfg.Stop, 'RBWAuto', false, 'RBW', cfg.RBW, 'Points', cfg.Points, ...
    'Preamp', logical(cfg.Preamp), 'RefLevel', opt.RefLevel, 'TraceType', "CLRW", 'Average', 1, 'Detector', det};
if isnumeric(opt.Atten), args = [args, {'AttenAuto', false, 'Atten', opt.Atten}];
else, args = [args, {'AttenAuto', true}]; end
ff.saConfigure(args{:});
st = ff.saSettings();
P = zeros(n, st.points, 'single'); t = NaT(n, 1); over = false;
for k = 1:n
    ff.sweep();
    y = ff.saTrace();
    if numel(y) ~= st.points, error('ff_auto_survey:trace', 'Trace has %d points, expected %d.', numel(y), st.points); end
    P(k, :) = y; t(k) = datetime('now');
    if canWaor
        w = ff.queryNum("CALC:MEAS:WAOR?");
        over = over || (isscalar(w) && w == 1);
    end
end
A = struct('f', st.f, 'P', P, 't', t, 'rbw', st.rbw, 'det', string(det), 'preamp', logical(st.preamp), ...
    'over', over, 'settings', rmfield(st, 'f'));
end

function Z = liveZS(ff, cfg, n)
ff.write("SENS:FREQ:SPAN 0");                                % verified: zero span
ff.write(sprintf("SENS:FREQ:CENT %.10g", cfg.Center));
ff.write(sprintf("SENS:BAND:RES %.10g", cfg.RBW));
ff.write(sprintf("SENS:SWE:POIN %d", cfg.Points));
ff.write("SENS:DET:FUNC " + cfg.Det);
ff.write(sprintf("SENS:POW:RF:GAIN %d", logical(cfg.Preamp)));
ff.write(sprintf("SENS:SWE:TIME %.6g", cfg.T));              % verified in zero span
ff.opc(30);
e = ff.errors();
span = ff.queryNum("SENS:FREQ:SPAN?");
if isempty(span) || span ~= 0, error('ff_auto_survey:zerospan', 'Zero span not accepted (span %g).', span); end
N = ff.queryNum("SENS:SWE:POIN?"); rbw = ff.queryNum("SENS:BAND:RES?");
T = str2double(tryQuery(ff, "SENS:SWE:TIME?"));
if ~(T > 0), T = cfg.T; end
E = zeros(n, N, 'single');
for k = 1:n
    ff.sweep(max(30, 5*T + 10));
    y = ff.saTrace();
    if numel(y) ~= N, error('ff_auto_survey:trace', 'Zero-span trace has %d points, expected %d.', numel(y), N); end
    E(k, :) = y;
end
Z = struct('T', T, 'rbw', rbw, 'E', E, 'errors', strjoin(e, '; '));
end

function zsRestore(ff, tmp)
ff.loadState(tmp);
ff.write("MMEM:DEL """ + tmp + """");
ff.query("INIT:CONT 0;*OPC?", 30);
ff.errors();
end

function r = tryQuery(ff, cmd)
w = warning('off', 'all'); c = onCleanup(@() warning(w)); %#ok<NASGU>
try, r = string(ff.query(cmd, 5)); catch, r = ""; end
if isempty(r), r = ""; end
end

function ok = probe(ff, cmd)
r = tryQuery(ff, cmd); e = ff.errors();
ok = strlength(r) > 0 && isempty(e);
end

%% ===== acquisition ===========================================================

function raw = newRaw(io, plan, opt)
seg = struct('name', {}, 'f', {}, 'rbw', {}, 'preamp', {}, 'Pk', {}, 'Av', {}, 'tPk', {}, 'tAv', {}, ...
    'over', {}, 'settings', {});
for s = 1:height(plan)
    seg(s) = struct('name', plan.Name(s), 'f', [], 'rbw', plan.RBW_Hz(s), 'preamp', opt.PreampMode ~= "off", ...
        'Pk', [], 'Av', [], 'tPk', NaT(0, 1), 'tAv', NaT(0, 1), 'over', false, 'settings', struct());
end
raw = struct('seg', seg, 'zoom', emptyZoomRaw(), 'zs', emptyZsRaw(), 'idn', io.idn, ...
    'source', io.kind, 'options', io.options, 'tStart', io.now(), 'tEnd', NaT, 'location', opt.Location, ...
    'notes', strings(0, 1), 'synthetic', io.kind == "demo");
if isfield(io, 'demoNote'), raw.notes(end + 1) = io.demoNote; end
end

function z = emptyZoomRaw()
z = struct('target', {}, 'lo', {}, 'hi', {}, 'seg', {}, 'f', {}, 'P', {}, 't', {}, 'rbw', {}, 'preamp', {}, 'kind', {});
end

function z = emptyZsRaw()
z = struct('target', {}, 'kind', {}, 'why', {}, 'seg', {}, 'rbw', {}, 'T', {}, 'det', {}, 'E', {}, 't', {}, 'preamp', {});
end

function raw = acquireSurvey(io, plan, opt, BP)
raw = newRaw(io, plan, opt);
t0 = tic;
fprintf('Coarse scan: %d segment(s) x %d repeat(s)\n', height(plan), opt.Repeats);
raw = coarseScan(io, raw, plan, opt, opt.Repeats, opt.Sweeps, true);
if strlength(raw.abortMsg) > 0
    if ~any(arrayfun(@(s) ~isempty(s.Pk) || ~isempty(s.Av), raw.seg))
        error(char(ternary(strlength(raw.abortId) > 0, raw.abortId, "ff_auto_survey:scan")), '%s', raw.abortMsg);
    end
    raw.notes(end + 1) = "Coarse scan stopped early (" + raw.abortMsg + "): the result covers the completed segments only.";
    warning('ff_auto_survey:scan', '%s', raw.notes(end));
    raw.tEnd = io.now();
    return
end
savePartial(raw, opt);
fprintf('  coarse scan done in %.1f s\n', toc(t0));
% --- zoom on the most informative detections
C = analyzeCoarse(raw, opt);
C.emitters = clusterEmitters(C.emitters, BP, C.seg);
Ez = pickZoom(C, opt, BP);
dead = false;
for i = 1:numel(Ez)
    e = Ez(i);
    sg = raw.seg(e.seg);
    w = max(e.hi - e.lo, 2*sg.rbw);
    if ismember(e.kind, ["rise" "cluster"])   % look for spectral lines (comb) inside: narrow span, fine RBW
        span = min(3*w, max(20*sg.rbw, 0.02*e.fpk)); span = min(max(span, 50e3), 20e6); fc = e.fpk;
        rbwZ = snapRBW(max(span/500, 100));
    else
        span = min(max(3*w, 20e3), 300e6); fc = (e.lo + e.hi)/2;
        rbwZ = snapRBW(min(max(w/60, 100), 1e6));
    end
    lo = max(5e3, fc - span/2); hi = min(6.5e9, lo + span); lo = max(5e3, hi - span);
    cfg = struct('Start', lo, 'Stop', hi, 'RBW', rbwZ, 'Points', opt.ZoomPoints, 'Preamp', sg.preamp);
    try
        A = io.sa(cfg, "AVER", opt.ZoomSweeps);
        raw.zoom(end + 1) = struct('target', fc, 'lo', e.lo, 'hi', e.hi, 'seg', e.seg, 'f', A.f, 'P', A.P, ...
            't', A.t, 'rbw', A.rbw, 'preamp', A.preamp, 'kind', e.kind);
        fprintf('  zoom %2d/%d  %-12s span %-10s RBW %s\n', i, numel(Ez), fmtHz(fc), fmtHz(hi - lo), fmtHz(A.rbw));
    catch err
        raw.notes(end + 1) = sprintf("Zoom at %s failed: %s", fmtHz(fc), err.message);
        warning('ff_auto_survey:zoom', '%s', raw.notes(end));
        if isConnectionError(err), dead = true; break, end
    end
end
% --- zero span on the most interesting ones
Zt = pickZeroSpan(C, raw, opt, BP);
if dead && ~isempty(Zt)
    Zt = Zt([]); raw.notes(end + 1) = "Zero span skipped: the instrument stopped answering during the zoom.";
end
if ~isempty(Zt)
    began = false;
    try
        io.zsBegin(); began = true;
    catch err
        raw.notes(end + 1) = "Could not save the swept setup before zero span: " + err.message;
        if isConnectionError(err), dead = true; Zt = Zt([]); end
    end
    for i = 1:numel(Zt)
        z = Zt(i);
        cfg = struct('Center', z.target, 'RBW', z.rbw, 'Points', 1001, 'Det', z.det, 'T', z.T, 'Preamp', z.preamp);
        try
            Z = io.zs(cfg, opt.ZeroSpanSweeps);
            raw.zs(end + 1) = struct('target', z.target, 'kind', z.kind, 'why', z.why, 'seg', z.seg, 'rbw', Z.rbw, ...
                'T', Z.T, 'det', z.det, 'E', Z.E, 't', io.now(), 'preamp', z.preamp);
            fprintf('  zero span %-12s %s %3.0f ms RBW %-8s (%s)\n', fmtHz(z.target), z.det, 1e3*Z.T, fmtHz(Z.rbw), z.why);
        catch err
            raw.notes(end + 1) = sprintf("Zero span at %s failed: %s", fmtHz(z.target), err.message);
            warning('ff_auto_survey:zerospan', '%s', raw.notes(end));
            if isConnectionError(err), dead = true; break, end
        end
    end
    if began && ~dead
        try, io.zsEnd(); catch err, raw.notes(end + 1) = "Restore after zero span: " + err.message; end
    end
end
raw.tEnd = io.now();
fprintf('Measurement finished in %.1f s\n', toc(t0));
end

function raw = coarseScan(io, raw, plan, opt, nRep, nSw, verbose)
% Never throws: an instrument error stops the scan and is returned in raw.abortMsg/abortId
% (segments completed so far are kept). A timeout is retried once per segment.
raw.abortMsg = ""; raw.abortId = "";
for r = 1:nRep
    for s = 1:height(plan)
        for attempt = 1:2
            try
                raw = scanSegment(io, raw, plan, opt, s, r, nRep, nSw, verbose);
                break
            catch err
                if attempt == 1 && strcmp(err.identifier, 'FieldFox:timeout')
                    raw.notes(end + 1) = sprintf("%s: %s Retrying once.", plan.Name(s), err.message);
                    warning('ff_auto_survey:retry', '%s', raw.notes(end));
                    continue
                end
                raw.abortMsg = string(err.message); raw.abortId = string(err.identifier);
                return
            end
        end
    end
end
end

function raw = scanSegment(io, raw, plan, opt, s, r, nRep, nSw, verbose)
sg = raw.seg(s);
cfg = struct('Start', plan.Start_Hz(s), 'Stop', plan.Stop_Hz(s), 'RBW', plan.RBW_Hz(s), ...
    'Points', plan.Points(s), 'Preamp', sg.preamp);
tt = tic;
for d = 1:numel(opt.Detectors)
    det = opt.Detectors(d);
    A = io.sa(cfg, det, nSw);
    if d == 1 && r == 1 && opt.PreampMode == "auto" && cfg.Preamp && ...
            (A.over || max(A.P(:)) > opt.PreampMax)
        cfg.Preamp = false; sg.preamp = false;
        raw.notes(end + 1) = sprintf("%s: preamp switched off (peak %.1f dBm%s).", plan.Name(s), ...
            max(A.P(:)), ternary(A.over, ", ADC over-range", ""));
        A = io.sa(cfg, det, nSw);
    end
    if ~isempty(sg.f) && numel(A.f) ~= numel(sg.f)
        error('ff_auto_survey:axis', '%s: the trace length changed (%d -> %d).', plan.Name(s), numel(sg.f), numel(A.f));
    end
    sg.f = A.f; sg.rbw = A.rbw; sg.preamp = A.preamp; sg.over = sg.over || A.over; sg.settings = A.settings;
    if isAvgDet(det), sg.Av = [sg.Av; A.P]; sg.tAv = [sg.tAv; A.t];
    else, sg.Pk = [sg.Pk; A.P]; sg.tPk = [sg.tPk; A.t]; end
end
raw.seg(s) = sg;
if verbose
    fprintf('  [%d/%d] %-22s RBW %-8s preamp %-3s %5.1f s%s\n', r, nRep, plan.Name(s), fmtHz(sg.rbw), ...
        ternary(sg.preamp, "on", "off"), toc(tt), ternary(sg.over, "  ADC OVER-RANGE", ""));
end
end

function tf = isConnectionError(err)
tf = any(strcmp(err.identifier, {'FieldFox:timeout', 'transportlib:client:ConnectionClosed'})) || ...
    contains(lower(err.message), ["connection", "not connected", "invalid or deleted"]);
end

function savePartial(raw, opt)
if ~opt.Save || strlength(opt.OutDir) == 0, return, end
try
    survey = struct('raw', raw, 'partial', true); %#ok<NASGU>
    save(fullfile(opt.OutDir, "survey_partial.mat"), 'survey', '-v7.3');
catch err
    warning('ff_auto_survey:save', 'Partial save failed: %s', err.message);
end
end

function E = pickZoom(C, opt, BP)
E = C.emitters([]);
if opt.MaxZoom <= 0 || isempty(C.emitters), return, end
X = C.emitters;
rank = zeros(1, numel(X));
for i = 1:numel(X)
    e = X(i);
    if e.kind == "comb", rank(i) = -Inf; continue, end
    nb = e.nb;
    r = e.exc;
    if ismember(e.kind, ["rise" "cluster"]), r = e.exc + 20; end   % EMI first: resolved lines identify the source
    if nb < 30, r = r + 10; end                          % poorly resolved by the coarse RBW
    b = bandOf(BP, e.fc);
    if b.kind == "broadcast_fm", r = r - 15; end         % plenty of these; zoom only a few
    if e.hi - e.lo > 150e6 && e.kind ~= "rise", r = -Inf; end
    rank(i) = r;
end
[rs, o] = sort(rank, 'descend');
o = o(isfinite(rs));
E = X(o(1:min(end, opt.MaxZoom)));
end

function Z = pickZeroSpan(C, raw, opt, BP)
Z = struct('target', {}, 'kind', {}, 'why', {}, 'seg', {}, 'rbw', {}, 'T', {}, 'det', {}, 'preamp', {});
if opt.MaxZeroSpan <= 0 || isempty(C.emitters), return, end
X = C.emitters;
pr = -Inf(1, numel(X)); why = strings(1, numel(X)); kd = strings(1, numel(X)); rb = zeros(1, numel(X));
for i = 1:numel(X)
    e = X(i); b = bandOf(BP, e.fc); w = e.hi - e.lo;
    if e.kind == "comb"
        pr(i) = 100 + e.exc; why(i) = "comb line"; kd(i) = "B"; rb(i) = snapRBW(min(e.combSp/3, 30e3));
        X(i).fc = e.combLine;
    elseif e.kind == "rise"
        pr(i) = 80 + e.exc; why(i) = "broadband rise"; kd(i) = "B"; rb(i) = snapRBW(min(5e6, max(w/3, 10e3)));
    elseif e.kind == "cluster"
        pr(i) = 85 + e.exc; why(i) = "EMI cluster"; kd(i) = "B"; rb(i) = snapRBW(min(5e6, max(w/3, 10e3)));
        X(i).fc = e.fpk;
    elseif b.kind == "imt_tdd" && w >= 5e6
        pr(i) = 70 + e.exc; why(i) = "TDD band"; kd(i) = "A"; rb(i) = 5e6;
    elseif b.kind == "rnss" || b.kind == "ism" || b.kind == "srd" || b.kind == "other"
        pr(i) = 50 + e.exc; why(i) = "EMI-prone / unlicensed band"; kd(i) = "AB"; rb(i) = snapRBW(min(5e6, max(w/2, 1e3)));
    elseif b.kind == "radiolocation" || (b.kind == "aero_radionav" && e.par > 10)
        pr(i) = 30 + e.exc; why(i) = "radar band"; kd(i) = "C"; rb(i) = snapRBW(min(5e6, max(w, 1e3)));
    elseif ismember(b.kind, ["broadcast_fm" "broadcast_am" "broadcast_tv"])
        pr(i) = e.exc - 20; why(i) = "strong emitter"; kd(i) = "A"; rb(i) = snapRBW(min(5e6, max(w/2, 1e3)));
    else
        pr(i) = 20 + e.exc; why(i) = "unexplained / " + b.kind; kd(i) = "AB"; rb(i) = snapRBW(min(5e6, max(w/2, 1e3)));
    end
end
[~, o] = sort(pr, 'descend');
o = o(isfinite(pr(o)));
for i = o(1:min(end, opt.MaxZeroSpan))
    e = X(i);
    for k = char(kd(i))      % A: 20 ms AVER (TDD slots)  B: 100 ms SAMP (bursts, mains)  C: 100 ms POS (short pulses)
        switch k
            case 'A', T = 20e-3; det = "AVER";
            case 'C', T = 100e-3; det = "POS";
            otherwise, T = 100e-3; det = "SAMP";
        end
        Z(end + 1) = struct('target', e.fc, 'kind', string(k), 'why', why(i), 'seg', e.seg, 'rbw', rb(i), ...
            'T', T, 'det', det, 'preamp', raw.seg(e.seg).preamp); %#ok<AGROW>
    end
end
end

%% ===== offline input =========================================================

function raw = loadRaw(src, opt)
files = string(src);
if isscalar(files) && isfolder(files)
    if isfile(fullfile(files, "survey.mat"))
        files = fullfile(files, "survey.mat");
    else
        d = dir(fullfile(files, '**', '*.mat'));
        d = d(~startsWith({d.name}, 'survey_partial'));
        if isempty(d), error('ff_auto_survey:file', 'No .mat files in %s.', src); end
        files = string(fullfile({d.folder}, {d.name}));
    end
end
for i = 1:numel(files)
    if ~isfile(files(i)), error('ff_auto_survey:file', 'File not found: %s', files(i)); end
end
if isscalar(files)
    v = who('-file', files);
    if ismember('survey', v)
        L = load(files, 'survey');
        if ~isfield(L.survey, 'raw'), error('ff_auto_survey:file', '%s holds no raw survey data.', files); end
        raw = L.survey.raw;
        raw.source = "survey file " + files;
        fprintf('Loaded survey %s (%d segments, %d zooms, %d zero-span captures).\n', files, numel(raw.seg), numel(raw.zoom), numel(raw.zs));
        return
    end
end
raw = struct('seg', [], 'zoom', emptyZoomRaw(), 'zs', emptyZsRaw(), 'idn', "", 'source', "recordings", ...
    'options', "", 'tStart', NaT, 'tEnd', NaT, 'location', opt.Location, 'notes', strings(0, 1), 'synthetic', false);
seg = struct('name', {}, 'f', {}, 'rbw', {}, 'preamp', {}, 'Pk', {}, 'Av', {}, 'tPk', {}, 'tAv', {}, 'over', {}, 'settings', {});
keys = strings(0, 1); nUsed = 0;
for i = 1:numel(files)
    try
        v = who('-file', files(i));
        if ~all(ismember({'f', 'P'}, v))
            raw.notes(end + 1) = sprintf("Skipped %s: not a spectrum recording (needs f and P).", files(i)); continue
        end
        r = load(files(i), v{ismember(v, {'f', 'P', 't', 'settings', 'idn', 'rbw'})});
    catch err
        raw.notes(end + 1) = sprintf("Skipped %s: %s", files(i), err.message); continue
    end
    if numel(r.f) < 21
        raw.notes(end + 1) = sprintf("Skipped %s: not a swept spectrum recording (needs f and P).", files(i)); continue
    end
    f = double(r.f(:).'); P = double(r.P);
    if size(P, 2) ~= numel(f) && size(P, 1) == numel(f), P = P.'; end
    if size(P, 2) ~= numel(f), raw.notes(end + 1) = sprintf("Skipped %s: P does not match f.", files(i)); continue, end
    st = struct(); if isfield(r, 'settings') && isstruct(r.settings), st = r.settings; end
    rbw = getf(st, 'rbw', NaN); if isfield(r, 'rbw'), rbw = r.rbw; end
    if ~(rbw > 0), rbw = 2*(f(end) - f(1))/(numel(f) - 1); raw.notes(end + 1) = sprintf("%s: RBW unknown, assumed %s.", files(i), fmtHz(rbw)); end
    det = upper(string(getf(st, 'detector', "POS")));
    t = NaT(size(P, 1), 1);
    if isfield(r, 't') && isdatetime(r.t) && numel(r.t) == size(P, 1), t = r.t(:); end
    key = sprintf('%.9g_%.9g_%d_%.6g', f(1), f(end), numel(f), rbw);
    k = find(keys == key, 1);
    if isempty(k)
        keys(end + 1) = key; k = numel(keys);
        nm = string(getf(st, 'band', sprintf("%s-%s", fmtHz(f(1)), fmtHz(f(end)))));
        seg(k) = struct('name', nm, 'f', f, 'rbw', rbw, 'preamp', logical(getf(st, 'preamp', false)), 'Pk', [], 'Av', [], ...
            'tPk', NaT(0, 1), 'tAv', NaT(0, 1), 'over', logical(getf(st, 'overload', false)), 'settings', st);
    end
    if isAvgDet(det), seg(k).Av = [seg(k).Av; single(P)]; seg(k).tAv = [seg(k).tAv; t];
    else, seg(k).Pk = [seg(k).Pk; single(P)]; seg(k).tPk = [seg(k).tPk; t]; end
    if strlength(raw.idn) == 0 && isfield(r, 'idn'), raw.idn = string(r.idn); end
    nUsed = nUsed + 1;
end
if isempty(seg), error('ff_auto_survey:file', 'No usable spectrum recordings in the input.'); end
[~, o] = sort(arrayfun(@(s) s.f(1), seg)); raw.seg = seg(o);
tt = [vertcat(raw.seg.tPk); vertcat(raw.seg.tAv)];
if any(~isnat(tt)), raw.tStart = min(tt); raw.tEnd = max(tt); end
fprintf('Loaded %d recording(s) into %d segment(s).\n', nUsed, numel(raw.seg));
end

function v = getf(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

%% ===== analysis ===============================================================

function R = analyzeSurvey(raw, opt, BP, CAT)
fprintf('Analysing ...\n');
C = analyzeCoarse(raw, opt);
E = clusterEmitters(C.emitters, BP, C.seg);
notes = [raw.notes(:); C.notes(:)];
tSurvey = surveyTime(raw, opt);
% --- zoom and zero span, attached to the nearest detection
Zr = struct('target', {}, 'seen', {}, 'f', {}, 'x', {}, 'fl', {}, 'rbw', {}, 'obw', {}, 'lo', {}, 'hi', {}, ...
    'nLines', {}, 'comb', {}, 'tech', {}, 'techScore', {}, 'top3', {}, 'emitter', {}, 'm', {});
for i = 1:numel(raw.zoom)
    z = analyzeZoom(raw.zoom(i), opt);
    kp = ""; if isfield(raw.zoom, 'kind'), kp = string(raw.zoom(i).kind); end
    j = nearestEmitter(E, raw.zoom(i).target, raw.zoom(i).lo, raw.zoom(i).hi, kp);
    z.emitter = j;
    if j > 0, E(j) = attachZoom(E(j), z); end
    Zr(end + 1) = z; %#ok<AGROW>
end
for i = 1:numel(Zr)                  % unresolved coarse blobs that are pieces of a comb seen in a zoom
    if isempty(Zr(i).comb) || isempty(E), continue, end
    sp = Zr(i).comb.spacing;
    for j = find([E.kind] == "signal" & [E.fc] >= Zr(i).f(1) & [E.fc] <= Zr(i).f(end) & [E.obw] >= 2*sp)
        if j == Zr(i).emitter || isfinite(E(j).zoomComb), continue, end
        E(j).zoomComb = sp;
    end
end
ZS = struct('target', {}, 'kind', {}, 'why', {}, 'res', {}, 'emitter', {}, 'raw', {});
for i = 1:numel(raw.zs)
    zz = raw.zs(i);
    res = analyzeZS(zz);
    kp = ternary(zz.why == "broadband rise", "rise", ternary(zz.why == "comb line", "comb", ternary(zz.why == "EMI cluster", "cluster", "")));
    j = nearestEmitter(E, zz.target, zz.target, zz.target, kp);
    if j > 0, E(j) = attachZS(E(j), res, zz); end
    ZS(end + 1) = struct('target', zz.target, 'kind', zz.kind, 'why', zz.why, 'res', res, 'emitter', j, 'raw', i); %#ok<AGROW>
end
% --- technology labels (ff_signal_id) on the coarse traces for unzoomed emitters
[m, snote] = sidMethod(opt);
if strlength(snote) > 0, notes(end + 1) = snote; end
if m ~= "off" && ~isempty(E)
    for s = 1:numel(C.seg)
        G = C.seg(s);
        todo = find([E.seg] == s & ~[E.zoomed] & [E.kind] == "signal");
        if isempty(G.f) || isempty(todo), continue, end
        T = runSignalId(G.f, G.avMean, G.rbw, m, "average");
        for j = todo
            [tech, sc, top] = matchSid(T, E(j).lo, E(j).hi);
            if strlength(tech) > 0, E(j).tech = tech; E(j).techScore = sc; E(j).top3 = top; E(j).techFrom = "coarse"; end
        end
    end
end
% --- band plan, receiver-generated products, satellites, classification
E = intermodCheck(E);
svc = startsWith(BP.Kind, "imt") | ismember(BP.Kind, ["broadcast_tv" "broadcast_fm" "rlan"]) | contains(BP.Name, "Wi-Fi");
for j = find(arrayfun(@(e) ismember(e.kind, ["rise" "cluster"]), E))   % how much of a hump lies in carrier bands
    x = linspace(E(j).lo, E(j).hi, 200);
    in = any(BP.FLow_Hz(svc).' <= x(:) & BP.FHigh_Hz(svc).' >= x(:), 2);
    E(j).svcFrac = mean(in);
    if E(j).svcFrac > 0
        b = BP.Name(svc & BP.FLow_Hz <= E(j).hi & BP.FHigh_Hz >= E(j).lo);
        E(j).svcBands = strjoin(b(1:min(end, 3)), ", ");
    end
end
for j = find(arrayfun(@(e) e.kind == "signal", E))         % narrow signals sitting inside an EMI hump/cluster
    k = find(arrayfun(@(x) ismember(x.kind, ["rise" "cluster"]) && x.svcFrac <= 0.5 && x.seg == E(j).seg && E(j).fc >= x.lo && E(j).fc <= x.hi, E), 1);
    if ~isempty(k), E(j).inEmi = sprintf("inside the EMI %s %s-%s", E(k).kind, fmtHz(E(k).lo), fmtHz(E(k).hi)); end
end
satNote = strings(0, 1);
for j = 1:numel(E)
    E(j).band = bandOf(BP, E(j).fc);
    if opt.SatID && E(j).kind == "signal"
        [E(j).sat, sn] = satMatch(E(j), tSurvey, CAT, opt, C.seg(E(j).seg));
        satNote = [satNote; sn(:)]; %#ok<AGROW>
    end
    E(j) = classifyEmitter(E(j), opt);
end
notes = [notes; unique(satNote)];
inv = inventoryTable(E, C);
R = struct('plan', table(), 'inventory', inv, 'bandOccupancy', bandOccupancy(C, BP, inv), ...
    'combs', combTable(C, Zr), 'rises', riseTable(E), 'zoom', Zr, 'zeroSpan', zsTable(ZS, E), 'zsDetail', ZS, ...
    'baseline', struct(), 'watch', struct(), 'conclusions', strings(0, 1), 'notes', notes, ...
    'segments', C.seg, 'emitters', E, 'raw', raw, 'files', struct(), 'outDir', "", ...
    'meta', struct('source', raw.source, 'idn', raw.idn, 'time', tSurvey, 'location', opt.Location, ...
    'synthetic', raw.synthetic, 'preset', opt.Preset, 'opt', rmfield(opt, {'AlertFcn'})));
end

function C = analyzeCoarse(raw, opt)
segs = struct('name', {}, 'f', {}, 'rbw', {}, 'df', {}, 'preamp', {}, 'avMean', {}, 'pkHold', {}, 'pkMed', {}, ...
    'flA', {}, 'flP', {}, 'flS', {}, 'refA', {}, 'thrA', {}, 'thrP', {}, 'thrS', {}, 'sA', {}, 'nPk', {}, 'nAv', {}, ...
    'combs', {}, 'over', {}, 'hitFrac', {});
E = newEmitter(); E = E([]);
notes = strings(0, 1);
for s = 1:numel(raw.seg)
    sg = raw.seg(s);
    if isempty(sg.f) || (isempty(sg.Pk) && isempty(sg.Av))
        segs(s).name = sg.name; segs(s).f = []; continue
    end
    [G, Es] = analyzeSegment(sg, s, opt);
    segs(s) = G;
    E = [E, Es]; %#ok<AGROW>
end
for j = 1:numel(E), E(j).id = j; end
C = struct('seg', segs, 'emitters', E, 'notes', notes);
end

function [G, E] = analyzeSegment(sg, s, opt)
f = double(sg.f(:).'); N = numel(f); df = (f(end) - f(1))/max(N - 1, 1); rbw = double(sg.rbw);
Pk = double(sg.Pk); Av = double(sg.Av);
if isempty(Pk), Pk = Av; end
if isempty(Av), AvS = Pk; else, AvS = Av; end
avMean = 10*log10(mean(10.^(AvS/10), 1));
pkHold = max(Pk, [], 1); pkMed = median(Pk, 1);
[flA, sA, refA, fl0A] = noiseFloor(avMean, 40, f);
[flP, sP] = noiseFloor(pkHold, 40, f);
flS = noiseFloor(pkMed, 40, f);
rr = Pk - flS; rr = rr(:);                                  % spread of single peak sweeps
s0 = max(median(rr) - prctile(rr, 16), 0.2); rn = rr(rr < median(rr) + 3*s0);
sS = max(prctile(rn, 84) - median(rn), median(rn) - prctile(rn, 16));
% an averaged (RMS) trace has ~0.5-1.5 dB noise spread; more means structure the floor model
% did not follow (dense comb, steep HF floor), which must not hide steady carriers -> cap it.
% Peak traces keep their own (possibly inflated) spread: bursty hash would shatter otherwise.
thrA = max(opt.Threshold, opt.K*min(sA, 1.5)); thrP = max(opt.Threshold, opt.K*sP); thrS = max(opt.Threshold, opt.K*sS);
eA = avMean - flA; eP = pkHold - flP;
seed = eA > thrA | eP > thrP;
low = eA > max(2, thrA - 3) | eP > max(2, thrP - 3);
regs = growRegions(seed, low, max(1, round(0.5*rbw/df)));
regs = splitRegions(regs, max(eA, eP), opt.SplitDip);
G = struct('name', sg.name, 'f', f, 'rbw', rbw, 'df', df, 'preamp', sg.preamp, 'avMean', avMean, 'pkHold', pkHold, ...
    'pkMed', pkMed, 'flA', flA, 'flP', flP, 'flS', flS, 'refA', refA, 'thrA', thrA, 'thrP', thrP, 'thrS', thrS, ...
    'sA', sA, 'nPk', size(sg.Pk, 1), 'nAv', size(sg.Av, 1), 'combs', [], 'over', sg.over, 'hitFrac', []);
hitPk = (Pk - flS) > thrS;                                  % per-sweep occupancy (peak sweeps)
if ~isempty(Av), hitAv = (Av - flA) > thrA; else, hitAv = false(0, N); end
G.hitFrac = mean([hitPk; hitAv], 1);
E = newEmitter(); E = E([]);
covered = false(1, N);
for r = 1:numel(regs)
    idx = regs{r};
    if ~(numel(idx) >= 2 || max(eA(idx)) >= thrA + 3 || max(eP(idx)) >= thrP + 3), continue, end
    useAv = max(eA(idx)) >= thrA;
    if useAv, x = avMean; fl = flA; else, x = pkHold; fl = flP; end
    m = measureRegion(f, x, fl, idx, rbw, df);
    e = newEmitter();
    e.seg = s; e.kind = "signal"; e.idx = idx; e.lo = m.lo; e.hi = m.hi; e.fc = m.fc; e.obw = m.obw;
    e.bw26 = m.bw26; e.bw3 = m.bw3; e.shoulder = m.shoulder; e.flat = m.flat; e.edge = m.edge; e.fpk = m.fpk;
    e.pk = max(pkHold(idx)); e.av = max(avMean(idx)); e.chp = m.chp; e.exc = max(max(eA(idx)), max(eP(idx)));
    e.excAv = max(eA(idx)); e.par = max(pkHold(idx)) - max(avMean(idx)); e.nb = numel(idx); e.resolved = m.resolved;
    e.rbw = rbw; e.df = df; e.onlyPeak = ~useAv;
    hp = any(hitPk(:, idx), 2); ha = any(hitAv(:, idx), 2);
    e.occ = 100*mean([hp; ha]);
    e.lvl = 10*log10(mean(10.^(Pk(:, idx)/10), 2)).';
    lin = 10.^(e.lvl/10); e.tcv = std(lin)/max(mean(lin), eps);
    E(end + 1) = e; %#ok<AGROW>
    covered(idx) = true;
end
% --- fragments of one pulsed/bursty emitter (radar, TDMA) seen in few bins per sweep
% and pieces of one wide bursty emitter (Wi-Fi, TDD) split at random dips of the max hold
for pass = 1:2
    j = 1;
    while j < numel(E)
        a = E(j); b = E(j + 1);
        if pass == 1
            ok = b.lo - a.hi <= max(3*rbw, 4*df) && a.occ <= 50 && b.occ <= 50 && a.nb <= 3 && b.nb <= 3 && abs(a.pk - b.pk) <= 12;
        else
            ok = (b.idx(1) - a.idx(end))*df <= max(4*rbw, 6*df) && (a.occ < 90 || a.onlyPeak) && ...
                (b.occ < 90 || b.onlyPeak) && a.obw > 3*rbw && b.obw > 3*rbw && abs(a.pk - b.pk) <= 10;
        end
        if ok
            idx = a.idx(1):b.idx(end);
            m = measureRegion(f, pkHold, flP, idx, rbw, df);
            a.idx = idx; a.lo = m.lo; a.hi = m.hi; a.fc = m.fc; a.obw = m.obw; a.bw26 = m.bw26; a.bw3 = m.bw3;
            a.shoulder = m.shoulder; a.flat = m.flat; a.edge = m.edge; a.fpk = m.fpk; a.chp = m.chp;
            a.pk = max(pkHold(idx)); a.av = max(avMean(idx)); a.exc = max(max(eA(idx)), max(eP(idx)));
            a.excAv = max(eA(idx)); a.par = a.pk - a.av; a.onlyPeak = a.onlyPeak && b.onlyPeak;
            a.frag = a.frag + max(b.frag, 1); a.resolved = m.resolved;
            hp = any(hitPk(:, idx), 2); ha = any(hitAv(:, idx), 2); a.occ = 100*mean([hp; ha]);
            a.lvl = 10*log10(mean(10.^(Pk(:, idx)/10), 2)).';
            lin = 10.^(a.lvl/10); a.tcv = std(lin)/max(mean(lin), eps);
            E(j) = a; E(j + 1) = [];
        else
            j = j + 1;
        end
    end
end
for j = 1:numel(E), E(j).nb = numel(E(j).idx); end
% --- harmonic combs (max-hold excess, then averaged excess)
combs = combSearch(f, eP, df, rbw, thrP, 2);
if isempty(combs), combs = combSearch(f, eA, df, rbw, thrA, 2); end
G.combs = combs;
for c = 1:numel(combs)
    cb = combs(c);
    mem = false(1, numel(E));
    for j = 1:numel(E)
        if E(j).obw > cb.spacing/2 && E(j).nb > 3, continue, end
        mm = round(E(j).fpk/cb.spacing);
        mem(j) = mm >= cb.m(1) && mm <= cb.m(end) && abs(E(j).fpk - mm*cb.spacing) <= max(1.5*df, cb.spacing/6);
    end
    e = newEmitter();
    e.seg = s; e.kind = "comb"; e.lo = max(f(1), cb.fpk(1) - cb.spacing/2); e.hi = min(f(end), cb.fpk(end) + cb.spacing/2);
    e.fc = cb.fpk(1); e.fpk = cb.fpk(1); e.obw = e.hi - e.lo; e.idx = find(f >= e.lo & f <= e.hi);
    e.pk = max(pkHold(e.idx)); e.av = max(avMean(e.idx)); e.exc = max(cb.v); e.excAv = max(eA(e.idx));
    e.nb = numel(e.idx); e.rbw = rbw; e.df = df; e.lines = cb.nHarm; e.combSp = cb.spacing; e.combCI = cb.ci;
    [~, kmax] = max(cb.v); e.combLine = cb.fpk(kmax); e.members = nnz(mem);
    if any(mem), e.occ = max([E(mem).occ]); e.lvl = E(find(mem, 1)).lvl; else, e.occ = 100; end
    e.chp = 10*log10(sum(10.^(cb.lvl/10)));
    e.bw26 = NaN; e.bw3 = NaN; e.shoulder = NaN; e.flat = NaN; e.par = e.pk - e.av; e.resolved = true;
    E = [E(~mem), e];
end
% --- broadband rises of the floor
rise = movmean(fl0A, max(3, round(N/80))) - refA;
minW = max(5, min(round(0.03*N), round(5e6/df)));          % 3 % of the segment or 5 MHz
[rs, re] = runs(rise > opt.RiseDB);
for k = 1:numel(rs)
    idx = rs(k):re(k);
    if numel(idx) < minW || mean(covered(idx)) > 0.5, continue, end
    e = newEmitter();
    e.seg = s; e.kind = "rise"; e.idx = idx; e.lo = f(idx(1)); e.hi = f(idx(end)); e.fc = (e.lo + e.hi)/2;
    [~, im] = max(rise(idx)); e.fpk = f(idx(im));
    e.obw = e.hi - e.lo; e.pk = max(pkHold(idx)); e.av = mean(avMean(idx)); e.exc = max(rise(idx));
    e.excAv = e.exc; e.nb = numel(idx); e.rbw = rbw; e.df = df;
    e.chp = 10*log10(max(sum(10.^(fl0A(idx)/10) - 10.^(refA(idx)/10)), eps)*df/(1.0645*rbw));
    e.occ = 100; e.lvl = 10*log10(mean(10.^(Pk(:, idx)/10), 2)).'; e.resolved = true;
    e.bw26 = NaN; e.bw3 = NaN; e.shoulder = NaN; e.flat = NaN; e.par = NaN;
    E(end + 1) = e; %#ok<AGROW>
end
end

function [fl, sig, ref, fl0] = noiseFloor(x, nBlk, f)
% Local noise floor (dB) and noise spread (dB) of a trace. Blocks: 20th percentile,
% running minimum over 5 blocks, then corrected to the median of the noise points.
% Segments spanning more than a decade use log-spaced blocks (steep HF floors).
% fl0 = the same without the running minimum (for broadband rises).
if nargin < 2, nBlk = 40; end
x = double(x(:).'); N = numel(x);
nb = max(1, min(nBlk, floor(N/15)));
ed = round(linspace(0, N, nb + 1));
if nargin > 2 && nb >= 3 && f(1) > 0 && f(end)/f(1) > 10
    p = interp1(double(f), 0:N - 1, logspace(log10(f(1)), log10(f(end)), nb + 1), 'linear', 'extrap');
    ed = unique([0, round(p(2:end - 1)), N]);
    k = 2;
    while k <= numel(ed)                                   % at least 15 points per block
        if ed(k) - ed(k - 1) < 15 && k < numel(ed), ed(k) = []; else, k = k + 1; end
    end
    if numel(ed) > 2 && ed(end) - ed(end - 1) < 15, ed(end - 1) = []; end
    nb = numel(ed) - 1;
end
q = zeros(1, nb); c = zeros(1, nb);
for b = 1:nb
    ix = ed(b) + 1:ed(b + 1);
    q(b) = prctile(x(ix), 20); c(b) = mean(ix);
end
if nb >= 3
    qm = movmin(q, 5);
    fl = interp1(c, qm, 1:N, 'linear', 'extrap');
    fl = movmean(fl, max(3, round(min(diff(ed)))));
    fl0 = interp1(c, q, 1:N, 'linear', 'extrap');
else
    fl = repmat(min(q), 1, N); fl0 = fl;
end
% second pass: median of the noise-only points per block
r = x - fl;
s0 = max(median(r) - prctile(r, 16), 0.2);
isN = r < median(r) + 3*s0;
corr = nan(1, nb);
for b = 1:nb
    ix = ed(b) + 1:ed(b + 1); ok = isN(ix);
    if mean(ok) > 0.2, corr(b) = median(r(ix(ok))); end
end
if all(isnan(corr)), corr(:) = median(r(isN)); end
corr = fillmissing(corr, 'linear', 'EndValues', 'nearest');
if nb >= 2, cI = interp1(c, corr, 1:N, 'linear', 'extrap'); else, cI = corr; end
fl = fl + cI;
r0 = x - fl0; m0 = median(r0(r0 < median(r0) + 3*s0));
fl0 = fl0 + m0;
r = x - fl;
rn = r(r < 3*s0 + max(0, median(r)));
if numel(rn) < 10, rn = r; end
sLo = median(rn) - prctile(rn, 16); sHi = prctile(rn, 84) - median(rn);
sig = max([sLo, sHi, 0.1]);
% reference = the segment's lowest floor: line through the lowest blocks. Over a wide
% frequency ratio the line is in log10(f) with its slope limited to -28..+3 dB/decade
% (man-made noise falls ~28 dB/decade, P.372): the natural VHF/HF slope is not a "rise".
if nb >= 4
    qb = interp1(1:N, fl, c);
    if nargin > 2 && f(1) >= 30e6 && f(end)/f(1) >= 3     % VHF/UHF: antenna noise slope (HF keeps the old line)
        xb = log10(interp1(1:N, double(f), c)); xN = log10(max(double(f), f(1)));
        k = qb <= prctile(qb, 50);
        if nnz(k) >= 2, p = polyfit(xb(k), qb(k), 1); else, p = [0 min(qb)]; end
        p(1) = min(max(p(1), -28), 3);
        p(2) = prctile(qb - p(1)*xb, 25);
        ref = min(polyval(p, xN), fl);
    else
        k = qb <= prctile(qb, 30);
        if nnz(k) >= 2, p = polyfit(c(k), qb(k), 1); else, p = [0 min(qb)]; end
        ref = min(polyval(p, 1:N), fl);
    end
else
    ref = fl;
end
end

function regs = growRegions(seed, low, gap)
low = low | seed;
[s, e] = runs(low);
regs = {};
if isempty(s), return, end
S = s(1); Ee = e(1);
for k = 2:numel(s)
    if s(k) - Ee(end) - 1 <= gap, Ee(end) = e(k); else, S(end + 1) = s(k); Ee(end + 1) = e(k); end %#ok<AGROW>
end
for k = 1:numel(S)
    idx = S(k):Ee(k);
    if any(seed(idx)), regs{end + 1} = idx; end %#ok<AGROW>
end
end

function out = splitRegions(regs, e, dip)
out = {};
for r = 1:numel(regs)
    idx = regs{r};
    if numel(idx) < 7, out{end + 1} = idx; continue, end %#ok<AGROW>
    y = movmean(e(idx), 3);
    pad = min(y) - 1;
    [~, loc] = findpeaks([pad, y, pad], 'MinPeakProminence', dip);
    loc = loc - 1;
    if numel(loc) < 2, out{end + 1} = idx; continue, end %#ok<AGROW>
    cuts = zeros(1, numel(loc) - 1);
    for k = 1:numel(loc) - 1
        [~, mpos] = min(y(loc(k):loc(k + 1))); cuts(k) = loc(k) + mpos - 1;
    end
    b = [0, cuts, numel(idx)];
    for k = 1:numel(b) - 1
        sgm = idx(b(k) + 1:b(k + 1));
        if ~isempty(sgm), out{end + 1} = sgm; end %#ok<AGROW>
    end
end
end

function m = measureRegion(f, x, fl, idx, rbw, df)
N = numel(f);
i1 = max(1, idx(1) - 1); i2 = min(N, idx(end) + 1); ii = i1:i2;
lin = max(10.^(x(ii)/10) - 10.^(fl(ii)/10), 0);
if sum(lin) <= 0, lin = 10.^(x(ii)/10); end
cs = cumsum(lin)/sum(lin);
k1 = find(cs >= 0.005, 1); k2 = find(cs >= 0.995, 1);
m.lo = f(ii(k1)) - df/2; m.hi = f(ii(k2)) + df/2; m.obw = m.hi - m.lo; m.fc = (m.lo + m.hi)/2;
[pk, ip] = max(x(ii)); m.fpk = f(ii(ip)); ip = ii(ip);
m.chp = 10*log10(max(sum(lin), eps)*df/(1.0645*rbw));
a3 = find(x(ii) >= pk - 3); m.bw3 = (a3(end) - a3(1) + 1)*df;
lv = pk - 26;
jl = ip; while jl > 1 && x(jl) >= lv, jl = jl - 1; end
jr = ip; while jr < N && x(jr) >= lv, jr = jr + 1; end
m.bw26 = f(jr) - f(jl);
if max(fl(jl:jr)) > lv - 3, m.bw26 = NaN; end           % floor-limited
m.shoulder = min(1, m.bw3/max(m.obw, df));
ins = lin(k1:k2); ins = max(ins, max(ins)*1e-6);
m.flat = exp(mean(log(ins)))/mean(ins);
% edge slope (dB/MHz): from 3 dB below the peak down to 20 dB below (or the floor)
lvE = max(pk - 20, max(fl(ii)) + 3);
sl = NaN(1, 2);
if pk - 3 > lvE
    j3 = ip; while j3 > 1 && x(j3) >= pk - 3, j3 = j3 - 1; end
    j20 = j3; while j20 > 1 && x(j20) >= lvE, j20 = j20 - 1; end
    if j3 > j20, sl(1) = (pk - 3 - lvE)/((f(j3) - f(j20))/1e6); end
    j3 = ip; while j3 < N && x(j3) >= pk - 3, j3 = j3 + 1; end
    j20 = j3; while j20 < N && x(j20) >= lvE, j20 = j20 + 1; end
    if j20 > j3, sl(2) = (pk - 3 - lvE)/((f(j20) - f(j3))/1e6); end
end
m.edge = mean(sl, 'omitnan');
m.resolved = m.obw > 3*max(rbw, df);
end

function C = combSearch(f, e, df, rbw, thr, maxCombs)
% Harmonic combs: power cepstrum of the dB excess, then the lines near each multiple.
C = struct('spacing', {}, 'ci', {}, 'nHarm', {}, 'm', {}, 'fpk', {}, 'v', {}, 'lvl', {}, 'rmse', {});
x = max(e(:), 0); x(~isfinite(x)) = 0; N = numel(x);
thr = min(thr, 10);              % a dense comb inflates the noise-spread estimate; lines must also stand out locally
minSp = max(3*df, 1.5*rbw); maxSp = (f(end) - f(1))/6;     % Gaussian RBW: lines 1.5 RBW apart dip ~7 dB
if N < 64 || maxSp <= 2*minSp, return, end
for it = 1:maxCombs
    if nnz(x > thr) < 5, return, end
    [pq, qf] = periodogram(x - mean(x), hamming(N), 2^nextpow2(4*N), 1/df);
    ok = qf >= 1/maxSp & qf <= 1/minSp;
    if ~any(ok), return, end
    [pmax, i0] = max(pq.*ok); dq = qf(2) - qf(1);
    for k = 4:-1:2                                       % prefer the largest spacing (first rahmonic)
        j = find(ok & abs(qf - qf(i0)/k) <= 2*dq);
        if ~isempty(j) && max(pq(j)) >= 0.5*pmax, [~, jj] = max(pq(j)); i0 = j(jj); break, end
    end
    sp = 1/qf(i0);
    m = max(1, ceil(f(1)/sp)):floor(f(end)/sp);
    if numel(m) < 5, return, end
    hit = false(size(m)); fpk = nan(size(m)); v = nan(size(m)); lv = nan(size(m));
    for h = 1:numel(m)
        w = find(abs(f - m(h)*sp) <= max(1.5*df, sp/6));
        if isempty(w), continue, end
        [vv, iw] = max(e(w)); k = w(iw);
        nbh = abs(f - m(h)*sp) <= sp/2;
        hit(h) = vv >= thr && vv - median(e(nbh)) >= 4;
        fpk(h) = f(k); v(h) = vv; lv(h) = vv;
    end
    if nnz(hit) < 5, return, end
    span = find(hit, 1, 'last') - find(hit, 1) + 1;
    if nnz(hit)/span < 0.4, return, end
    % a real comb has lines AT the multiples and valleys half-way; dense random hash
    % (bursty EMI, many carriers) "hits" any grid -> compare with the half-shifted grid
    hitH = false(size(m));
    for h = 1:numel(m)
        w = find(abs(f - (m(h) + 0.5)*sp) <= max(1.5*df, sp/6));
        if isempty(w), continue, end
        nbh = abs(f - (m(h) + 0.5)*sp) <= sp/2;
        hitH(h) = max(e(w)) >= thr && max(e(w)) - median(e(nbh)) >= 4;
    end
    if (nnz(hit) - nnz(hitH))/numel(m) < 0.3, return, end
    mh = m(hit); fh = fpk(hit);
    spf = sum(mh.*fh)/sum(mh.^2);
    rmse = sqrt(mean((fh - mh*spf).^2));
    if rmse > max(df, rbw/2), return, end
    se = rmse/sqrt(sum(mh.^2)) + df/sqrt(12*sum(mh.^2));
    tq = tinv(0.975, max(1, nnz(hit) - 1));
    C(end + 1) = struct('spacing', spf, 'ci', spf + [-1 1]*tq*se, 'nHarm', nnz(hit), 'm', mh, 'fpk', fh, ...
        'v', v(hit), 'lvl', lv(hit), 'rmse', rmse); %#ok<AGROW>
    for h = fh, x(abs(f(:) - h) <= max(1.5*df, spf/6)) = 0; end
end
end

function z = analyzeZoom(zr, opt)
f = double(zr.f(:).'); P = double(zr.P); df = (f(end) - f(1))/max(numel(f) - 1, 1);
x = 10*log10(mean(10.^(P/10), 1));
[fl, sg] = noiseFloor(x, 12);
% an averaged trace has ~1-2 dB noise spread; more means lines/structure everywhere
% (e.g. a dense comb) - cap it so that real carriers are not hidden by the threshold
sg = min(sg, 1.5);
thr = max(opt.Threshold, opt.K*sg);
e = x - fl;
regs = growRegions(e > thr, e > max(2, thr - 3), max(1, round(0.5*zr.rbw/df)));
idx = [];
for r = 1:numel(regs)
    ir = regs{r};
    if f(ir(end)) >= zr.lo - 2*df && f(ir(1)) <= zr.hi + 2*df, idx = [idx, ir]; end %#ok<AGROW>
end
z = struct('target', zr.target, 'seen', ~isempty(idx), 'f', f, 'x', x, 'fl', fl, 'rbw', zr.rbw, 'obw', NaN, ...
    'lo', NaN, 'hi', NaN, 'nLines', 0, 'comb', [], 'tech', "", 'techScore', NaN, 'top3', "", 'emitter', 0, 'm', []);
if ~isempty(idx)
    idx = min(idx):max(idx);
    m = measureRegion(f, x, fl, idx, zr.rbw, df);
    z.m = m; z.obw = m.obw; z.lo = m.lo; z.hi = m.hi;
    [pk, ~, wd] = findpeaks(e, 'MinPeakProminence', 6, 'MinPeakHeight', thr);
    z.nLines = nnz(wd*df <= 4*zr.rbw & pk > 0);
end
cb = combSearch(f, e, df, zr.rbw, thr, 1);
if ~isempty(cb), z.comb = cb; end
[mth, ~] = sidMethod(opt);
if mth ~= "off" && z.seen
    T = runSignalId(f, P, zr.rbw, mth, "average");
    [z.tech, z.techScore, z.top3] = matchSid(T, z.lo, z.hi);
end
end

function e = attachZoom(e, z)
e.zoomed = true;
if ~isempty(z.comb), e.zoomComb = z.comb.spacing; e.lines = max(e.lines, z.comb.nHarm); end
if ~z.seen
    e.zoomNote = "not seen in the zoom (intermittent?)";
    return
end
if e.kind == "signal"
    m = z.m;
    e.obw = m.obw; e.lo = m.lo; e.hi = m.hi; e.fc = m.fc; e.bw26 = m.bw26; e.bw3 = m.bw3;
    e.shoulder = m.shoulder; e.flat = m.flat; e.edge = m.edge; e.fpk = m.fpk;
    e.resolved = m.resolved; e.rbwZoom = z.rbw;
end
e.lines = max(e.lines, z.nLines);
if strlength(z.tech) > 0 && e.kind == "signal", e.tech = z.tech; e.techScore = z.techScore; e.top3 = z.top3; e.techFrom = "zoom"; end
end

function res = analyzeZS(zz)
E = double(zz.E); [ns, N] = size(E);
t = linspace(0, zz.T, N);
res = struct('kind', zz.kind, 'label', "", 'tdd', [], 'burst', [], 'mains', false, 'impulsive', false, ...
    'continuous', false, 'period', NaN, 'pattern', "", 'duty', NaN, 'rate', NaN, 'crest', NaN, 'contrast', NaN, 'thr', NaN);
if zz.kind == "A"
    best = []; bp = -Inf;
    for k = 1:ns
        r = tddAnalyze(t, E(k, :), 0.5e-3);
        if r.periodicity > bp, bp = r.periodicity; best = r; end
    end
    res.tdd = best; res.label = replace(string(best.class), "TDD periodic", "TDD-like periodic on/off"); res.period = best.period;
    res.pattern = string(best.slotPattern); res.duty = best.onFrac;
    res.continuous = startsWith(res.label, "continuous");
    return
end
fs = (N - 1)/zz.T;
x = 10.^(E/20);
thr = noiseContrast(N, fs, zz.rbw);
[con, pp, dc, cr] = deal(nan(ns, 1)); ps = 0; fp = [];
for k = 1:ns
    lv = statelevels(x(k, :));
    con(k) = 20*log10(lv(2)/max(lv(1), eps));
    cr(k) = max(E(k, :)) - median(E(k, :));
    if con(k) >= thr
        try
            pp(k) = mean(pulseperiod(x(k, :), fs, 'StateLevels', lv));
            dc(k) = mean(dutycycle(x(k, :), fs, 'StateLevels', lv));
        catch
        end
    end
    [p1, fp] = pspectrum(x(k, :) - mean(x(k, :)), fs);
    ps = ps + p1/ns;
end
ok = fp > 2/zz.T;
[~, j] = max(ps.*ok); envPk = fp(j);
pkSNR = 10*log10(ps(j)/median(ps(ok)));
per = median(pp, 'omitnan'); cv = std(pp, 'omitnan')/mean(pp, 'omitnan');
res.contrast = median(con); res.thr = thr; res.crest = median(cr); res.duty = median(dc, 'omitnan');
res.rate = 1/per;
near = @(a, b) abs(a - b) <= 0.06*b;
isM = @(r) near(r, 50) || near(r, 100);
res.burst = struct('contrast', res.contrast, 'thr', thr, 'period', per, 'cv', cv, 'envPk', envPk, 'envSNR', pkSNR);
if res.contrast < thr
    res.continuous = true;
    res.label = sprintf("continuous / noise-like (contrast %.1f dB < %.1f dB noise limit)", res.contrast, thr);
    if pkSNR > 15 && isM(envPk)
        res.mains = true; res.label = sprintf("mains-synchronous modulation at %.0f Hz (envelope spectrum)", envPk);
    end
    return
end
if isfinite(per) && cv < 0.2 && isM(1/per)
    res.mains = true; res.label = sprintf("mains-synchronous bursts at %.1f Hz (duty %.0f%%)", 1/per, 100*res.duty);
elseif pkSNR > 15 && isM(envPk)
    res.mains = true; res.label = sprintf("mains-synchronous bursts (envelope peak %.0f Hz)", envPk);
elseif isfinite(per) && cv < 0.2
    res.label = sprintf("periodic bursts at %.4g Hz (duty %.0f%%)", 1/per, 100*res.duty);
else
    res.label = "irregular bursts";
end
if res.crest > 15 && (isnan(res.duty) || res.duty < 0.1)
    res.impulsive = true; res.label = res.label + sprintf(", impulsive (crest %.0f dB)", res.crest);
end
end

function e = attachZS(e, res, zz)
if zz.kind == "A"
    e.tdd = res.label; e.tddPeriod = res.period; e.tddPattern = res.pattern;
    if startsWith(res.label, "continuous"), e.zsCont = true; end
else
    e.zsLabel = joinS(e.zsLabel, res.label);
    e.mains = e.mains || res.mains; e.impulsive = e.impulsive || res.impulsive;
    e.zsCont = e.zsCont || res.continuous;
end
e.zsDone = true;
end

function E = intermodCheck(E)
% Harmonics (2f, 3f) and 3rd/2nd-order products of strong emitters may be made by
% the analyzer itself (or a nearby non-linearity): flag weaker lines that sit there.
if isempty(E), return, end
k = find([E.kind] == "signal" & [E.pk] > -50);
if isempty(k), return, end
[~, o] = sort([E(k).pk], 'descend'); k = k(o(1:min(end, 6)));
fs = [E(k).fc]; cand = []; lab = strings(0, 1);
for a = 1:numel(fs)
    cand = [cand, 2*fs(a), 3*fs(a)]; lab = [lab; "2x" + fmtHz(fs(a)); "3x" + fmtHz(fs(a))]; %#ok<AGROW>
    for b = 1:numel(fs)
        if a == b, continue, end
        cand = [cand, 2*fs(a) - fs(b), fs(a) + fs(b)]; %#ok<AGROW>
        lab = [lab; "2x" + fmtHz(fs(a)) + "-" + fmtHz(fs(b)); fmtHz(fs(a)) + "+" + fmtHz(fs(b))]; %#ok<AGROW>
    end
end
for j = 1:numel(E)
    if E(j).kind ~= "signal" || any(j == k) || E(j).pk > -50, continue, end
    tol = max([2*E(j).df, E(j).obw/2, 2*E(j).rbw]);
    h = find(abs(cand - E(j).fc) <= tol, 1);
    if ~isempty(h), E(j).intermod = lab(h); end
end
end

function E = clusterEmitters(E, BP, segs)
% Many narrow detections packed together outside channelised services (inverter "hash",
% a comb the RBW cannot resolve, arcing) -> one "cluster" emitter. Members more than
% 10 dB above the cluster median stay separate (a real carrier inside the hash).
if numel(E) < 8, return, end
skip = ["imt_fdd_dl" "imt_fdd_ul" "imt_tdd" "imt_sdl" "broadcast_fm" "broadcast_tv" "rlan" "rnss" ...
    "aero" "aero_radionav" "maritime" "landmobile" "broadcast_am"];   % channelised: a crowd there is traffic
drop = false(1, numel(E)); add = E([]);
for s = unique([E.seg])
    ks = find([E.seg] == s & [E.kind] == "signal");
    if numel(ks) < 8, continue, end
    [~, o] = sort([E(ks).fc]); ks = ks(o);
    ok = arrayfun(@(j) ~ismember(bandOf(BP, E(j).fc).kind, skip), ks);
    rbw = E(ks(1)).rbw; if s <= numel(segs) && isfield(segs, 'rbw') && ~isempty(segs(s).rbw), rbw = segs(s).rbw; end
    i = 1;
    while i <= numel(ks)
        if ~ok(i), i = i + 1; continue, end
        j = i;
        while j < numel(ks) && ok(j + 1)
            med = median([E(ks(i:j + 1)).obw]);
            if E(ks(j + 1)).lo - E(ks(j)).hi <= max(4*med, 10*rbw), j = j + 1; else, break, end
        end
        mem = ks(i:j);
        if numel(mem) >= 8
            pk = [E(mem).pk]; mem = mem(pk <= median(pk) + 10);
        end
        if numel(mem) >= 8
            e = newEmitter(); M = E(mem);
            e.seg = s; e.kind = "cluster"; e.idx = unique([M.idx]); e.lo = min([M.lo]); e.hi = max([M.hi]);
            e.fc = (e.lo + e.hi)/2; e.obw = e.hi - e.lo; [e.pk, im] = max([M.pk]); e.fpk = M(im).fpk;
            e.av = max([M.av]); e.exc = max([M.exc]); e.excAv = max([M.excAv]); e.par = e.pk - e.av;
            e.chp = 10*log10(sum(10.^([M.chp]/10))); e.occ = mean([M.occ]); e.lvl = M(im).lvl; e.nb = numel(e.idx);
            e.rbw = M(1).rbw; e.df = M(1).df; e.lines = numel(mem); e.resolved = true;
            e.combSp = median(diff(sort([M.fpk])));             % typical spacing (not a proven comb)
            add(end + 1) = e; %#ok<AGROW>
            drop(mem) = true;
        end
        i = j + 1;
    end
end
if isempty(add), return, end
E = [E(~drop), add];
[~, o] = sort([E.fc]); E = E(o);
for j = 1:numel(E), E(j).id = j; end
end

function j = nearestEmitter(E, f0, lo, hi, kindPref)
j = 0;
if isempty(E), return, end
lo0 = [E.lo]; hi0 = [E.hi]; tol = max([E.df], 1);
in = (f0 >= lo0 - tol & f0 <= hi0 + tol) | (hi >= lo0 & lo <= hi0);
if nargin > 4 && strlength(kindPref) > 0 && any(in & [E.kind] == kindPref), in = in & [E.kind] == kindPref; end
if any(in)
    c = find(in); [~, k] = min(abs([E(c).fc] - f0) + 1e12*([E(c).kind] == "comb" & abs([E(c).combLine] - f0) > tol(c)));
    j = c(k);
else
    [d, k] = min(abs([E.fc] - f0));
    if d <= 3*max(E(k).obw, E(k).df), j = k; end
end
end

%% ===== signal ID, band plan, satellites ======================================

function [m, note] = sidMethod(opt)
persistent cacheM cacheN
note = "";
if opt.SignalID == "off" || exist('ff_signal_id', 'file') ~= 2
    m = "off";
    if opt.SignalID ~= "off", note = "ff_signal_id not found: technology labels from the band plan only."; end
    return
end
if opt.SignalID ~= "auto", m = opt.SignalID; return, end
if isempty(cacheM)
    mf = fullfile(fileparts(which('ff_signal_id')), 'ff_signal_models.mat');
    if isfile(mf), cacheM = "ensemble"; cacheN = "";
    else, cacheM = "template"; cacheN = "ff_signal_models.mat not found: ff_signal_id uses the template fit."; end
end
m = cacheM; note = cacheN;
end

function T = runSignalId(f, P, rbw, m, det)
T = table();
try
    w = warning('off', 'all'); c = onCleanup(@() warning(w)); %#ok<NASGU>
    [~, T] = evalc('ff_signal_id(f, P, rbw, ''Method'', m, ''Detector'', det, ''Plot'', false, ''Visible'', ''off'', ''MaxSignals'', 25)');
catch
    T = table();
end
end

function [tech, sc, top] = matchSid(T, lo, hi)
tech = ""; sc = NaN; top = "";
if ~istable(T) || height(T) == 0 || ~all(ismember({'Freq_Hz', 'OBW_Hz', 'Technology'}, T.Properties.VariableNames)), return, end
a = T.Freq_Hz - max(T.OBW_Hz, 1)/2; b = T.Freq_Hz + max(T.OBW_Hz, 1)/2;
ov = max(0, min(b, hi) - max(a, lo));
in = T.Freq_Hz >= lo & T.Freq_Hz <= hi;
ov(in) = ov(in) + max(hi - lo, 1);
[o, k] = max(ov);
if o <= 0, return, end
tech = string(T.Technology(k));
if ismember('Score', T.Properties.VariableNames), sc = T.Score(k); end
if ismember('Top3', T.Properties.VariableNames), top = string(T.Top3(k)); end
end

function BP = loadBandPlan(opt)
persistent cache key
file = string(opt.BandPlan);
if strlength(file) == 0, file = fullfile(opt.supp, "th_bandplan.csv"); end
if isfile(file)
    d = dir(file); k = file + "|" + d.datenum;
    if ~isempty(key) && key == k, BP = cache; return, end
end
if ~isfile(file)
    warning('ff_auto_survey:bandplan', 'Band plan %s not found: no band labels.', file);
    BP = table(strings(0, 1), zeros(0, 1), zeros(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
        'VariableNames', {'Name', 'FLow_Hz', 'FHigh_Hz', 'Service', 'Kind', 'Holder', 'Expect', 'Note'});
    return
end
o = detectImportOptions(file, 'TextType', 'string', 'Delimiter', ',');
o = setvartype(o, intersect(o.VariableNames, {'Name', 'Service', 'Kind', 'Duplex', 'Holder', 'Expect', 'Conf', 'Note'}), 'string');
BP = readtable(file, o);
need = ["Name" "FLow_Hz" "FHigh_Hz"];
if ~all(ismember(need, string(BP.Properties.VariableNames)))
    error('ff_auto_survey:bandplan', 'Band plan %s needs the columns Name, FLow_Hz, FHigh_Hz.', file);
end
for v = ["Kind" "Expect" "Holder" "Note" "Service"]
    if ~ismember(v, string(BP.Properties.VariableNames)), BP.(v) = strings(height(BP), 1); end
    BP.(v)(ismissing(BP.(v))) = "";
end
BP.Kind(BP.Kind == "") = "other";
cache = BP; key = k;
end

function b = bandOf(BP, fc)
b = struct('name', "", 'kind', "other", 'service', "", 'expect', "", 'holder', "", 'note', "", 'all', "");
if isempty(BP) || height(BP) == 0, return, end
k = find(BP.FLow_Hz <= fc & BP.FHigh_Hz >= fc);
if isempty(k), return, end
[~, i] = min(BP.FHigh_Hz(k) - BP.FLow_Hz(k)); i = k(i);
b = struct('name', BP.Name(i), 'kind', BP.Kind(i), 'service', BP.Service(i), 'expect', BP.Expect(i), ...
    'holder', BP.Holder(i), 'note', BP.Note(i), 'all', strjoin(BP.Name(k), " | "));
end

function CAT = loadSatCatalog(opt)
file = fullfile(opt.supp, "sat_channels.csv");
if ~isfile(file), CAT = table(); return, end
o = detectImportOptions(file, 'TextType', 'string', 'Delimiter', ',');
o = setvartype(o, intersect(o.VariableNames, {'Name', 'Modulation', 'Class', 'Note'}), 'string');
CAT = readtable(file, o);
end

function [S, notes] = satMatch(e, tUTC, CAT, opt, G)
S = struct('name', "", 'score', 0, 'note', "", 'gnss', false, 'el', NaN, 'doppler', NaN, 'cands', "");
notes = strings(0, 1);
tol = max(G.df, G.rbw);
fc = e.fc; obw = e.obw;
if isempty(CAT) || height(CAT) == 0, k = [];
else
    k = find(CAT.Active == 1 & fc >= CAT.FLow_Hz - CAT.MaxDoppler_Hz - CAT.OBW_Hz/2 - tol & ...
        fc <= CAT.FHigh_Hz + CAT.MaxDoppler_Hz + CAT.OBW_Hz/2 + tol);
end
if isempty(k)
    S = satnogsHint(S, e, opt);
    if ismember(e.band.kind, ["metsat" "spaceops" "amateur_sat" "amateur" "mss"]), [S, notes] = satsigHint(S, e, tUTC, opt, tol); end
    return
end
if any(CAT.Class(k) == "gnss")
    S.gnss = true; S.name = strjoin(CAT.Name(k(CAT.Class(k) == "gnss")), ", ");
    S.note = "GNSS band (" + S.name + "): GNSS is 17-25 dB below kTB, invisible on a swept SA";
    return
end
best = -Inf; names = strings(0, 1);
for i = k(:).'
    w = CAT.OBW_Hz(i);
    res = obw > 5*max(e.rbw, e.df);                 % else the measured width is mostly the RBW
    if res && (obw > 5*w + 2*tol || obw < w/5), continue, end
    if ~res && w > 4*max(e.rbw, e.df), continue, end
    sc = 0.35; nt = CAT.Name(i);
    if CAT.NORAD(i) > 0
        [el, rr, gnote] = satGeometry(CAT.NORAD(i), tUTC, opt);
        if strlength(gnote) > 0, notes(end + 1) = gnote; end %#ok<AGROW>
        if isfinite(el)
            if el > 0
                dop = -CAT.F_Hz(i)*rr/299792458;
                sc = sc + 0.25; nt = nt + sprintf(" visible now (el %.0f deg)", el);
                if abs(fc - (CAT.F_Hz(i) + dop)) <= tol + obw/2 + 2e3
                    sc = sc + 0.2; nt = nt + sprintf(", Doppler %+.1f kHz consistent", dop/1e3);
                end
                S.el = el; S.doppler = dop;
            else
                sc = sc - 0.35; nt = nt + sprintf(" below the horizon now (el %.0f deg)", el);
            end
        else
            nt = nt + " (visibility unknown: no TLE)";
        end
    else
        sc = sc + 0.1; nt = nt + " (constellation)";
    end
    if e.pk > -70, sc = sc - 0.3; nt = nt + ", but too strong for a satellite"; end
    names(end + 1) = nt; %#ok<AGROW>
    if sc > best, best = sc; S.name = CAT.Name(i); S.score = sc; S.note = "satellite catalog: " + nt; end
end
S.cands = strjoin(names, "; ");
if best == -Inf, S.score = 0; end
S = satnogsHint(S, e, opt);
[S, n2] = satsigHint(S, e, tUTC, opt, tol);
notes = [notes(:); n2(:)];
end

function [S, notes] = satsigHint(S, e, tUTC, opt, tol)
% Cross-check with ff_satsig("overhead"): its catalog satellites above the horizon at the
% survey time, with the Doppler-shifted receive frequency. Optional (exist + try).
notes = strings(0, 1);
[O, note] = satsigOverhead(tUTC, opt);
if strlength(note) > 0, notes(end + 1) = note; end
if isempty(O) || height(O) == 0, return, end
d = abs(O.Frx_MHz*1e6 - e.fc);
k = find(d <= tol + e.obw/2 + 3e3 + 5e-6*e.fc);
if isempty(k), return, end
[~, i] = min(d(k)); i = k(i);
nt = sprintf("ff_satsig: %s (%s) above the horizon (el %.0f deg), predicted %.4f MHz incl. Doppler %+.1f kHz, %.1f kHz from the measured centre", ...
    O.Satellite(i), O.Downlink(i), O.El_deg(i), O.Frx_MHz(i), O.Doppler_Hz(i)/1e3, d(i)/1e3);
S.note = strjoin([S.note(strlength(S.note) > 0), nt], "; ");
if e.pk <= -70 || S.score > 0
    S.score = max(S.score, 0.35) + 0.25;
    if strlength(S.name) == 0 || ~contains(S.name, O.Satellite(i), 'IgnoreCase', true), S.name = string(O.Satellite(i)); end
    S.el = O.El_deg(i); S.doppler = O.Doppler_Hz(i);
end
end

function [O, note] = satsigOverhead(tUTC, opt)
% ff_satsig("overhead") once per survey time (cached); empty if absent or failing.
persistent cache key
O = []; note = "";
if ~opt.SatID || isnat(tUTC) || exist('ff_satsig', 'file') ~= 2, return, end
k = sprintf('%s_%.4f_%.4f_%d', string(tUTC, 'yyyyMMddHHmm'), opt.Location(1), opt.Location(2), opt.Online);
if ~isempty(key) && strcmp(key, k), O = cache; return, end
try
    t = tUTC; t.TimeZone = 'UTC'; t1 = tic;
    w = warning('off', 'all'); c = onCleanup(@() warning(w)); %#ok<NASGU>
    [~, O] = evalc(['ff_satsig("overhead", ''Time'', t, ''Site'', opt.Location, ''Offline'', ~opt.Online, ' ...
        '''Plot'', false, ''Quiet'', true, ''Hours'', 0.05)']);
    fprintf('  ff_satsig: satellites overhead at %s UTC (%.1f s)\n', string(tUTC, 'yyyy-MM-dd HH:mm'), toc(t1));
    if ~istable(O) || ~all(ismember({'Satellite', 'Downlink', 'El_deg', 'Frx_MHz', 'Doppler_Hz'}, O.Properties.VariableNames))
        O = []; note = "ff_satsig overhead returned an unexpected table: not used.";
    else
        O.Satellite = string(O.Satellite); O.Downlink = string(O.Downlink);
    end
catch err
    O = []; note = "ff_satsig overhead failed (" + err.message + "): catalog/TLE check only.";
end
cache = O; key = k;
end

function S = satnogsHint(S, e, opt)
% SatNOGS DB transmitters marked alive near the carrier (a weak prior: stale flags exist).
% Only for weak carriers in satellite-capable bands or with a catalog candidate.
if ~opt.Online && strlength(string(opt.SatnogsFile)) == 0, return, end
if e.pk > -70 || ~(S.score > 0 || ismember(e.band.kind, ["metsat" "spaceops" "amateur_sat" "amateur" "mss"])), return, end
Tx = satnogsTransmitters(opt);
if isempty(Tx), return, end
fc = e.fc;
k = find(Tx.alive & abs(Tx.f - fc) <= 3e-5*fc + 2e3 + min(e.obw/2, 25e3));
if isempty(k), return, end
[~, o] = sort(abs(Tx.f(k) - fc)); k = k(o);
nm = unique(Tx.name(k), 'stable');
S.note = strjoin([S.note(strlength(S.note) > 0), "SatNOGS alive transmitters near: " + strjoin(nm(1:min(end, 3)), ", ")], "; ");
if S.score <= 0, S.score = 0.1; S.name = nm(1); end
end

function [el, rr, note] = satGeometry(norad, tUTC, opt)
persistent TL
if isempty(TL), TL = containers.Map('KeyType', 'double', 'ValueType', 'any'); end
el = NaN; rr = NaN; note = "";
if isnat(tUTC), return, end
if ~isKey(TL, norad)
    [tle, note] = loadTLE(norad, opt);
    TL(norad) = tle;
end
tle = TL(norad);
if isempty(tle), return, end
try
    t = tUTC + seconds([-1 0 1]);
    r = propagateOrbit(t, tle, "PropModel", "sgp4");
    aer = eci2aer(r.', datevec(t), repmat(opt.Location, 3, 1));
    el = aer(2, 2); rr = (aer(3, 3) - aer(1, 3))/2;
    ep = tle.Epoch; if ~isempty(ep.TimeZone), ep.TimeZone = 'UTC'; ep.TimeZone = ''; end
    age = days(tUTC - ep);
    if age > 14, note = sprintf("TLE for %d is %.0f days old: visibility is approximate.", norad, age); end
catch err
    note = sprintf("Orbit of %d not computed: %s", norad, err.message);
end
end

function [tle, note] = loadTLE(norad, opt)
tle = []; note = "";
d = string(opt.TleDir); file = fullfile(d, sprintf("%d.tle", norad)); errf = file + ".err";
if opt.Online && ~(isfile(file) && fileAgeH(file) < 2) && ~(isfile(errf) && fileAgeH(errf) < 2)
    url = sprintf('https://celestrak.org/NORAD/elements/gp.php?CATNR=%d&FORMAT=TLE', norad);
    try
        txt = string(webread(url, weboptions('ContentType', 'text', 'Timeout', 20, ...
            'UserAgent', 'ff_auto_survey (MATLAB; FieldFox survey; student project)')));
        if contains(txt, newline + "1 " + norad) || startsWith(strtrim(txt), "1 ") || contains(txt, "1 " + norad)
            if ~isfolder(d), mkdir(d); end
            writelines(strtrim(splitlines(strtrim(txt))), file);
        else
            writelines(string(datetime('now')) + " " + extractBefore(txt + " ", min(200, strlength(txt) + 1)), errf);
        end
    catch err
        if ~isfolder(d), mkdir(d); end
        writelines(string(datetime('now')) + " " + err.message, errf);
        note = "CelesTrak not reachable (" + err.message + "): cached TLEs only.";
    end
end
src = ""; tmp = "";
shared = fullfile(opt.here, "support", "tle");          % ff_cubesat_link's cache (read only)
if isfile(file)
    src = file;
elseif isfile(fullfile(shared, sprintf("%d.tle", norad)))
    src = fullfile(shared, sprintf("%d.tle", norad));
else                                                   % any other cached .tle that holds this object
    dd = [dir(fullfile(d, '*.tle')); dir(fullfile(shared, '*.tle'))];
    for i = 1:numel(dd)
        L = readlines(fullfile(dd(i).folder, dd(i).name));
        j = find(startsWith(L, "1 " + string(sprintf('%5d', norad))), 1);
        if ~isempty(j) && j < numel(L)
            tmp = string(tempname) + ".tle";
            nm = "OBJ " + norad; if j > 1 && ~startsWith(L(j - 1), ["1 " "2 "]), nm = L(j - 1); end
            writelines([nm; L(j); L(j + 1)], tmp); src = tmp; break
        end
    end
end
if strlength(src) == 0, return, end
try
    tle = tleread(src); tle = tle(1);
catch err
    note = sprintf("TLE for %d unreadable: %s", norad, err.message); tle = [];
end
if strlength(tmp) > 0 && isfile(tmp), delete(tmp); end   % the temporary extract (tempname)
end

function Tx = satnogsTransmitters(opt)
persistent cache
if ~isempty(cache), Tx = cache; return, end
Tx = [];
file = string(opt.SatnogsFile);
if strlength(file) == 0
    file = fullfile(opt.supp, "satnogs_transmitters.json"); errf = file + ".err";
    stale = ~isfile(file) || fileAgeH(file) > 24;
    blocked = isfile(errf) && fileAgeH(errf) < 24;
    if stale && opt.Online && ~blocked
        fid = -1;
        try
            txt = webread('https://db.satnogs.org/api/transmitters/?format=json', weboptions('ContentType', 'text', ...
                'Timeout', 60, 'UserAgent', 'ff_auto_survey (MATLAB; FieldFox survey; student project)'));
            fid = fopen(file, 'w'); fwrite(fid, txt, 'char');
        catch err
            writelines(string(datetime('now')) + " " + err.message, errf);
        end
        if fid >= 0, fclose(fid); end                    % also after a write error
    end
end
if ~isfile(file), return, end
try
    J = jsondecode(fileread(file));
    if iscell(J), J = J(cellfun(@isstruct, J)); n = numel(J); gf = @(i, fld) fieldOr(J{i}, fld);
    else, n = numel(J); gf = @(i, fld) fieldOr(J(i), fld); end
    f = nan(n, 1); al = false(n, 1); nm = strings(n, 1);
    for i = 1:n
        v = gf(i, 'downlink_low'); if isnumeric(v) && isscalar(v), f(i) = v; end
        a = gf(i, 'alive'); al(i) = isequal(a, true);
        d = gf(i, 'description'); m = gf(i, 'mode'); s = gf(i, 'norad_cat_id');
        nm(i) = sprintf("%s %s (NORAD %s)", string(d), string(m), string(s));
    end
    Tx = table(f, al, nm, 'VariableNames', {'f', 'alive', 'name'});
    Tx = Tx(isfinite(Tx.f), :);
    cache = Tx;
catch
    Tx = [];
end
end

function v = fieldOr(s, fld)
if isfield(s, fld), v = s.(fld); else, v = []; end
if isempty(v), v = ""; if strcmp(fld, 'downlink_low'), v = []; end, end
end

function t = surveyTime(raw, opt)
t = raw.tStart;
if ~isdatetime(t) || isnat(t), t = opt.Time; end
if ~isdatetime(t) || isnat(t), t = NaT; return, end
if isempty(t.TimeZone), t.TimeZone = 'local'; end
t.TimeZone = 'UTC'; t.TimeZone = '';
end

%% ===== classification ========================================================

function e = newEmitter()
e = struct('id', 0, 'seg', 0, 'kind', "signal", 'idx', [], 'lo', NaN, 'hi', NaN, 'fc', NaN, 'fpk', NaN, 'obw', NaN, ...
    'bw26', NaN, 'bw3', NaN, 'shoulder', NaN, 'flat', NaN, 'edge', NaN, 'pk', NaN, 'av', NaN, 'chp', NaN, ...
    'exc', NaN, 'excAv', NaN, 'par', NaN, 'nb', 0, 'resolved', false, 'rbw', NaN, 'df', NaN, 'onlyPeak', false, ...
    'occ', NaN, 'lvl', [], 'tcv', NaN, 'frag', 0, 'intermod', "", 'lines', 0, 'combSp', NaN, 'combCI', [NaN NaN], 'combLine', NaN, 'members', 0, ...
    'zoomed', false, 'rbwZoom', NaN, 'zoomComb', NaN, 'zoomNote', "", 'tech', "", 'techScore', NaN, 'top3', "", ...
    'techFrom', "", 'zsDone', false, 'zsLabel', "", 'tdd', "", 'tddPeriod', NaN, 'tddPattern', "", 'mains', false, ...
    'impulsive', false, 'zsCont', false, 'band', struct('name', "", 'kind', "other", 'service', "", 'expect', "", ...
    'holder', "", 'note', "", 'all', ""), 'sat', struct('name', "", 'score', 0, 'note', "", 'gnss', false, 'el', NaN, ...
    'doppler', NaN, 'cands', ""), 'inEmi', "", 'svcFrac', 0, 'svcBands', "", 'cls', "Unknown", 'sub', "", 'conf', 0, 'evidence', "", 'scores', [0 0 0 0]);
end

function e = classifyEmitter(e, opt) %#ok<INUSD>
lic = 0; sat = 0; emi = 0; unk = 0.3; ev = strings(0, 1); sub = ""; bandFit = false;
B = e.band; k = B.kind; tech = e.tech; ts = e.techScore; if isnan(ts), ts = 0.5; end
w = e.obw;
unres = e.kind == "signal" && (~e.resolved || (tech == "CW" && w <= 5*max(e.rbw, e.df)));
if unres && ismember(tech, ["CW" "Unresolved"])            % "CW" = narrower than the RBW: no shape information
    ev(end + 1) = sprintf("narrower than ~3 x RBW (%s): shape not resolved", fmtHz(max(e.rbw, e.df)));
    tech = "";
end
known = strlength(tech) > 0 && ~ismember(tech, ["Unknown" "Unresolved"]);
% ---- EMI signatures
if e.kind == "comb"
    emi = emi + 0.7;
    sub = sprintf("harmonic comb, spacing %s (%d lines, %s-%s)", fmtHz(e.combSp), e.lines, fmtHz(e.lo), fmtHz(e.hi));
    ev(end + 1) = sprintf("cepstrum + line fit: spacing %s [%s, %s]", fmtHz(e.combSp), fmtHz(e.combCI(1)), fmtHz(e.combCI(2)));
    if e.combSp >= 2e3 && e.combSp <= 200e3, ev(end + 1) = "spacing typical of a switching converter (inverter/SMPS)"; end
end
if e.kind == "rise"
    emi = emi + 0.55;
    sub = sprintf("broadband noise rise (+%.1f dB over %s)", e.exc, fmtHz(w));
    if e.hi <= 30e6, ev(end + 1) = "below 30 MHz this may be ambient man-made/atmospheric noise: compare with a baseline"; end
end
if ismember(e.kind, ["rise" "cluster"]) && e.svcFrac > 0.5
    emi = emi - 0.3; lic = lic + 0.45;
    ev(end + 1) = sprintf("%.0f%% of it lies in carrier bands (%s): probably weak wideband carriers (LTE/NR/DVB-T2/Wi-Fi) raising the floor, not EMI", ...
        100*e.svcFrac, e.svcBands);
end
if e.kind == "cluster"
    emi = emi + 0.6;
    sub = sprintf("dense cluster of %d narrow emissions, %s-%s (EMI hash / unresolved comb)", e.lines, fmtHz(e.lo), fmtHz(e.hi));
    ev(end + 1) = sprintf("%d narrow detections packed together (median spacing %s) outside channelised services", e.lines, fmtHz(e.combSp));
    if e.hi <= 30e6, ev(end + 1) = "below 30 MHz shortwave broadcast/amateur stations also crowd together: compare with the source switched off (Baseline)"; end
end
if isfinite(e.zoomComb)
    emi = emi + 0.45; ev(end + 1) = sprintf("zoom shows a line comb, spacing %s", fmtHz(e.zoomComb));
    if strlength(sub) == 0, sub = sprintf("comb inside the signal (spacing %s)", fmtHz(e.zoomComb));
    else, sub = sub + sprintf(", line comb %s spacing inside", fmtHz(e.zoomComb)); end
    if e.zoomComb >= 2e3 && e.zoomComb <= 200e3, ev(end + 1) = "comb spacing typical of a switching converter (inverter/SMPS)"; end
end
if e.mains
    emi = emi + 0.5; ev(end + 1) = "zero span: " + e.zsLabel;
    sub = strjoin([sub(strlength(sub) > 0), "mains-synchronous"], ", ");
end
if e.impulsive && ~ismember(k, ["radiolocation" "aero_radionav"])
    emi = emi + 0.3; ev(end + 1) = "impulsive in zero span";
    if strlength(sub) == 0, sub = "impulsive interference"; end
end
dme = false;
if k == "rnss" && e.kind == "signal"
    dme = e.fc >= 1164e6 && e.fc <= 1213e6 && (e.onlyPeak || e.par > 8 || e.impulsive || e.occ < 60);
    if dme                               % DME/TACAN ground replies use 962-1213 MHz (pulse pairs)
        emi = emi + 0.15;
        ev(end + 1) = sprintf("pulsed/bursty (peak - average %.0f dB, occupancy %.0f%%) in 1164-1213 MHz, where DME/TACAN (962-1213 MHz) shares the GNSS L5/E5 band", e.par, e.occ);
    else
        emi = emi + 0.6; ev(end + 1) = "visible inside a GNSS band (GNSS itself is below the noise)";
        if e.fc <= 1213e6, ev(end + 1) = "if it turns out pulsed (zero span), DME/TACAN 962-1213 MHz is the likely source"; end
        sub = "interference in a GNSS band";
    end
elseif e.sat.gnss
    emi = emi + 0.4; ev(end + 1) = e.sat.note;
end
if e.kind == "signal" && w > 1e6 && e.shoulder < 0.55 && e.flat < 0.75 && k ~= "broadcast_fm"
    if startsWith(k, "imt")              % partly loaded or narrow (GSM/NB-IoT) carriers look like this too
        emi = emi + 0.1; ev(end + 1) = sprintf("wide (%s) with gentle edges (shoulder %.2f) - in an IMT band partly loaded or narrow carriers look like this", fmtHz(w), e.shoulder);
    else
        emi = emi + 0.35; ev(end + 1) = sprintf("wide (%s) with gentle edges (shoulder %.2f): not a flat-top carrier", fmtHz(w), e.shoulder);
        if strlength(sub) == 0, sub = "broadband emission"; end
    end
end
if strlength(e.intermod) > 0
    emi = emi + 0.2; ev(end + 1) = "sits at " + e.intermod + " of strong emitters: may be generated in the analyzer (retry with more attenuation)";
    if strlength(sub) == 0, sub = "possible harmonic/intermodulation product (" + e.intermod + ")"; end
end
if e.frag > 0, ev(end + 1) = sprintf("seen in few bins per sweep (%d fragments): pulsed or TDMA", e.frag + 1); end
if tech == "Unknown", emi = emi + 0.15; unk = unk + 0.15; ev(end + 1) = "ff_signal_id: Unknown shape"; end
if ismember(k, ["ism" "srd" "other"]) && ~known, emi = emi + 0.1; end
% ---- licensed services
licTxt = "";
if known
    hasExp = strlength(B.expect) > 0;
    ok = ~hasExp || any(startsWith(tech, split(B.expect, "|")));
    if hasExp, lic = lic + 0.55*ts*(1 - 0.5*~ok) + 0.25*ok; else, lic = lic + 0.3*ts; end
    ev(end + 1) = sprintf("ff_signal_id (%s): %s, p=%.2f", e.techFrom, tech, ts);
    if ~ok, ev(end + 1) = sprintf("%s is not expected in %s", tech, B.name); end
    licTxt = tech;
end
if ismember(e.kind, ["rise" "cluster"]) && e.svcFrac > 0.5, licTxt = "weak wideband carrier(s) raising the floor in " + e.svcBands; end
licensedKinds = ["broadcast_am" "broadcast_fm" "broadcast_tv" "imt_fdd_dl" "imt_fdd_ul" "imt_tdd" "imt_sdl" ...
    "aero" "aero_radionav" "maritime" "landmobile" "amateur" "radiolocation" "fss" "rlan" "ism" "srd"];
if e.kind == "signal" && ismember(k, licensedKinds), lic = lic + 0.15; end
if e.kind == "signal"
    if k == "broadcast_fm" && w >= 80e3 && w <= 320e3
        lic = lic + 0.4; ev(end + 1) = "FM broadcast width"; licTxt = ternary(strlength(licTxt) > 0, licTxt, "FM broadcast");
    elseif k == "broadcast_fm" && unres && w <= 4*max(e.rbw, e.df) + 320e3
        lic = lic + 0.3; ev(end + 1) = "carrier in the FM band (width not resolved at this RBW)";
        licTxt = ternary(strlength(licTxt) > 0, licTxt, "FM broadcast (unresolved)");
    elseif k == "broadcast_fm" && w <= 1.2e6
        lic = lic + 0.3; ev(end + 1) = "FM band: adjacent stations not separated at this RBW";
        licTxt = "FM broadcast (adjacent stations)";
    elseif k == "broadcast_am" && (w <= 30e3 || (unres && e.rbw <= 30e3))
        lic = lic + 0.35; ev(end + 1) = "narrow carrier in the AM band"; licTxt = ternary(strlength(licTxt) > 0, licTxt, "AM broadcast");
    elseif ismember(k, ["ism" "rlan"]) && w >= 4e6 && w <= 170e6 && (e.occ < 90 || e.tcv > 0.4 || e.onlyPeak || e.par > 6)
        lic = lic + 0.45;
        ev(end + 1) = sprintf("bursty wideband emission in a Wi-Fi band (occupancy %.0f%%, peak - average %.0f dB): RLAN-like", e.occ, e.par);
        if strlength(licTxt) == 0 || ~any(startsWith(licTxt, ["Wi-Fi" "BLE" "Bluetooth"])), licTxt = "Wi-Fi / RLAN (licence-exempt, bursty)"; end
    elseif startsWith(k, "imt") && w >= 1.2e6 && w <= 110e6 && e.shoulder >= 0.7
        lic = lic + 0.35; ev(end + 1) = sprintf("flat-top carrier of IMT width (%s, shoulder %.2f)", fmtHz(w), e.shoulder);
        licTxt = ternary(strlength(licTxt) > 0, licTxt, "IMT carrier (LTE/NR)");
    elseif k == "broadcast_tv" && w >= 5.5e6 && w <= 8.5e6 && e.shoulder >= 0.7
        lic = lic + 0.4; ev(end + 1) = "DVB-T2 channel width"; licTxt = ternary(strlength(licTxt) > 0, licTxt, "DVB-T2");
    elseif ismember(k, ["aero" "maritime" "landmobile" "amateur"]) && w <= 30e3
        lic = lic + 0.35; ev(end + 1) = "narrowband channel in a " + B.service + " band";
    elseif k == "srd" && w <= 300e3
        lic = lic + 0.35; ev(end + 1) = "narrowband emission in an SRD band (remote controls, sensors, LoRa)";
        if strlength(licTxt) == 0 || ~known, licTxt = "SRD device (licence-exempt)"; end
    elseif ismember(k, ["imt_fdd_ul" "imt_tdd"]) && (e.occ < 90 || e.onlyPeak)
        lic = lic + 0.25; ev(end + 1) = sprintf("bursty emission in an IMT uplink/TDD band (occupancy %.0f%%): handset (UE) uplink is expected", e.occ);
        licTxt = ternary(strlength(licTxt) > 0, licTxt, "IMT uplink activity (UE)");
    elseif startsWith(k, "imt")
        lic = lic + 0.2; ev(end + 1) = "emission in an IMT band (operator carriers expected)";
    end
    if dme
        lic = lic + 0.55; licTxt = "DME/TACAN (aeronautical radionavigation, pulsed)"; ev(end + 1) = "pulsed in the DME band";
    end
    bandFit = any(contains(ev, ["FM broadcast width" "narrow carrier in the AM band" "narrowband channel in a" ...
        "narrowband emission in an SRD" "RLAN-like" "DVB-T2 channel width" "flat-top carrier of IMT" ...
        "adjacent stations not separated" "carrier in the FM band" "handset (UE) uplink" "pulsed in the DME band"]));
    if k == "imt_tdd" && strlength(e.tdd) > 0 && startsWith(e.tdd, "TDD-like")
        lic = lic + 0.4; ev(end + 1) = "zero span: " + e.tdd; bandFit = true;
        licTxt = sprintf("TDD carrier %s%s", ternary(isfinite(e.tddPeriod), sprintf("(%.2f ms", 1e3*e.tddPeriod), "("), ...
            ternary(strlength(e.tddPattern) > 0, ", " + e.tddPattern + ")", ")"));
    end
    if ismember(k, ["radiolocation" "aero_radionav"]) && (e.impulsive || (e.onlyPeak && e.par > 10) || e.par > 15 || ...
            contains(e.zsLabel, "periodic bursts"))
        lic = lic + 0.4; ev(end + 1) = sprintf("pulsed (peak - average %.0f dB): radar-like", e.par);
        licTxt = "pulsed radar"; bandFit = true;
        emi = max(0, emi - 0.4);                          % pulses explain the odd shape and the impulsiveness
    end
    % ff_signal_id's "Unknown shape" (no template fits) weighs less when width and band agree
    if bandFit && tech == "Unknown", emi = max(0, emi - 0.15); unk = unk - 0.15; end
end
if strlength(e.inEmi) > 0 && ~bandFit
    emi = emi + 0.25; ev(end + 1) = e.inEmi;
    if strlength(sub) == 0, sub = "narrow emission " + e.inEmi; end
end
% ---- satellites
if e.sat.score > 0 && ~e.sat.gnss
    sat = sat + e.sat.score; ev(end + 1) = e.sat.note;
    if ismember(k, ["metsat" "spaceops" "amateur_sat" "mss" "amateur"]), sat = sat + 0.15; end
elseif strlength(e.sat.note) > 0 && ~e.sat.gnss
    ev(end + 1) = e.sat.note;
end
if e.kind == "signal" && e.sat.score <= 0 && ~e.sat.gnss && ...
        (ismember(k, ["metsat" "amateur_sat"]) || (k == "spaceops" && e.fc < 401e6))
    emi = emi + 0.4;
    ev(end + 1) = "carrier in a satellite downlink band with no catalog satellite overhead now: terrestrial interference or an uncatalogued satellite";
    if strlength(sub) == 0, sub = "unexplained carrier in a satellite downlink band"; end
end
% ---- decide
sc = [lic sat emi unk];
names = ["Licensed service" "Satellite" "EMI" "Unknown"];
[b, i] = max(sc); s2 = max(sc(setdiff(1:4, i)));
conf = min(0.95, b)*(1 - 0.5*s2/max(b, eps));
if b < 0.35, i = 4; conf = min(conf, 0.3); end
e.cls = names(i); e.conf = round(conf, 2); e.scores = sc;
switch i
    case 1
        e.sub = ternary(strlength(licTxt) > 0, licTxt, B.service + " (" + B.name + ")");
        if strlength(B.holder) > 0 && ~ismember(B.holder, ["-" "various"]), e.sub = e.sub + " - " + B.holder; end
    case 2, e.sub = e.sat.name;
    case 3, e.sub = ternary(strlength(sub) > 0, sub, "unexplained emission");
    otherwise
        if unres && strlength(tech) == 0
            e.sub = sprintf("narrow carrier (not resolved at RBW %s)", fmtHz(e.rbw));
        elseif e.kind == "signal"
            e.sub = "nearest: " + ternary(strlength(e.top3) > 0, e.top3, ternary(strlength(tech) > 0, tech, B.service));
        else
            e.sub = ternary(strlength(sub) > 0, sub, "unexplained " + e.kind);
        end
end
if e.zoomed && strlength(e.zoomNote) > 0, ev(end + 1) = e.zoomNote; end
if strlength(e.zsLabel) > 0 && ~any(contains(ev, e.zsLabel)), ev(end + 1) = "zero span: " + e.zsLabel; end
if strlength(e.tdd) > 0 && ~any(contains(ev, e.tdd)), ev(end + 1) = "zero span (20 ms): " + e.tdd; end
if strlength(B.name) > 0, ev(end + 1) = "band: " + B.all; end
e.evidence = strjoin(unique(ev, 'stable'), "; ");
end

%% ===== tables =================================================================

function T = inventoryTable(E, C)
n = numel(E);
if n == 0                                                  % typed empty table (comparisons keep working)
    z = zeros(0, 1); s = strings(0, 1);
    T = table(z, z, z, z, z, z, z, z, z, z, z, s, s, z, s, s, s, s, s, z, false(0, 1), s, z, z, s, ...
        'VariableNames', {'ID', 'Freq_Hz', 'Lo_Hz', 'Hi_Hz', 'OBW_Hz', 'BW26_Hz', 'Peak_dBm', 'Avg_dBm', 'ChanPow_dBm', ...
        'Excess_dB', 'Occupancy_pct', 'Class', 'SubType', 'Confidence', 'Band', 'Technology', 'Satellite', 'Evidence', ...
        'Segment', 'RBW_Hz', 'Zoomed', 'ZeroSpan', 'Lines', 'Spacing_Hz', 'Kind'});
    return
end
segName = strings(n, 1);
for j = 1:n, segName(j) = C.seg(E(j).seg).name; end
zsS = strings(n, 1);
for j = 1:n, zsS(j) = strjoin([E(j).tdd(strlength(E(j).tdd) > 0), E(j).zsLabel(strlength(E(j).zsLabel) > 0)], "; "); end
col = @(v) reshape(v, [], 1);
T = table((1:n).', col([E.fc]), col([E.lo]), col([E.hi]), col([E.obw]), col([E.bw26]), col([E.pk]), col([E.av]), ...
    col([E.chp]), col([E.exc]), col([E.occ]), col([E.cls]), col([E.sub]), col([E.conf]), col(sfun(@(e) e.band.name, E)), ...
    col([E.tech]), col(sfun(@(e) e.sat.name, E)), col([E.evidence]), segName, col([E.rbw]), col([E.zoomed]), zsS, ...
    col([E.lines]), col([E.combSp]), col([E.kind]), ...
    'VariableNames', {'ID', 'Freq_Hz', 'Lo_Hz', 'Hi_Hz', 'OBW_Hz', 'BW26_Hz', 'Peak_dBm', 'Avg_dBm', 'ChanPow_dBm', ...
    'Excess_dB', 'Occupancy_pct', 'Class', 'SubType', 'Confidence', 'Band', 'Technology', 'Satellite', 'Evidence', ...
    'Segment', 'RBW_Hz', 'Zoomed', 'ZeroSpan', 'Lines', 'Spacing_Hz', 'Kind'});
if n == 0, return, end
[~, o] = sort(T.Freq_Hz); T = T(o, :); T.ID = (1:n).';
end

function T = bandOccupancy(C, BP, inv)
rows = {};
for i = 1:height(BP)
    lo = BP.FLow_Hz(i); hi = BP.FHigh_Hz(i);
    nb = 0; hit = 0; mx = -Inf; flr = [];
    for s = 1:numel(C.seg)
        G = C.seg(s);
        if isempty(G.f), continue, end
        k = G.f >= lo & G.f <= hi;
        if ~any(k), continue, end
        nb = nb + nnz(k); hit = hit + sum(G.hitFrac(k)); mx = max(mx, max(G.pkHold(k))); flr = [flr, G.flA(k)]; %#ok<AGROW>
    end
    if nb < 3, continue, end
    m = inv.Freq_Hz >= lo & inv.Freq_Hz <= hi;
    cls = strjoin(unique(inv.Class(m)), ", ");
    rows(end + 1, :) = {BP.Name(i), lo, hi, BP.Kind(i), 100*hit/nb, mx, median(flr), nnz(m), cls, nb}; %#ok<AGROW>
end
if isempty(rows)
    T = table(strings(0, 1), zeros(0, 1), zeros(0, 1), strings(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), strings(0, 1), zeros(0, 1), ...
        'VariableNames', {'Band', 'Lo_Hz', 'Hi_Hz', 'Kind', 'Occupancy_pct', 'MaxLevel_dBm', 'Floor_dBm', 'Emitters', 'Classes', 'Points'});
    return
end
T = cell2table(rows, 'VariableNames', {'Band', 'Lo_Hz', 'Hi_Hz', 'Kind', 'Occupancy_pct', 'MaxLevel_dBm', 'Floor_dBm', 'Emitters', 'Classes', 'Points'});
T = sortrows(T, 'Occupancy_pct', 'descend');
end

function T = combTable(C, Zr)
rows = {};
for s = 1:numel(C.seg)
    G = C.seg(s);
    for c = 1:numel(G.combs)
        cb = G.combs(c);
        rows(end + 1, :) = {G.name, cb.spacing, cb.ci(1), cb.ci(2), cb.nHarm, cb.fpk(1), cb.fpk(end), max(cb.v)}; %#ok<AGROW>
    end
end
for i = 1:numel(Zr)                                         % combs resolved only in a zoom
    cb = Zr(i).comb;
    if isempty(cb), continue, end
    rows(end + 1, :) = {"zoom at " + fmtHz(Zr(i).target), cb.spacing, cb.ci(1), cb.ci(2), cb.nHarm, cb.fpk(1), cb.fpk(end), max(cb.v)}; %#ok<AGROW>
end
vn = {'Segment', 'Spacing_Hz', 'CIlow_Hz', 'CIhigh_Hz', 'Lines', 'First_Hz', 'Last_Hz', 'MaxExcess_dB'};
if isempty(rows)
    z = zeros(0, 1); T = table(strings(0, 1), z, z, z, z, z, z, z, 'VariableNames', vn);
else
    T = cell2table(rows, 'VariableNames', vn); T.Segment = string(T.Segment);
end
end

function T = riseTable(E)
if isempty(E), z = zeros(0, 1); T = table(z, z, z, z, 'VariableNames', {'Lo_Hz', 'Hi_Hz', 'Rise_dB', 'ExcessPower_dBm'}); return, end
k = [E.kind] == "rise";
T = table([E(k).lo].', [E(k).hi].', [E(k).exc].', [E(k).chp].', 'VariableNames', {'Lo_Hz', 'Hi_Hz', 'Rise_dB', 'ExcessPower_dBm'});
end

function T = zsTable(ZS, E)
n = numel(ZS);
lab = strings(n, 1); cls = strings(n, 1); per = nan(n, 1); duty = nan(n, 1); con = nan(n, 1); thr = nan(n, 1);
for i = 1:n
    r = ZS(i).res; lab(i) = r.label; per(i) = r.period; duty(i) = r.duty; con(i) = r.contrast; thr(i) = r.thr;
    if ZS(i).kind ~= "A" && isfinite(r.rate), per(i) = 1/r.rate; end
    if ZS(i).emitter > 0, cls(i) = E(ZS(i).emitter).cls; end
end
fq = zeros(n, 1); kd = strings(n, 1); wy = strings(n, 1);
if n > 0, fq = [ZS.target].'; kd = string([ZS.kind]).'; wy = string([ZS.why]).'; end
T = table(fq, kd, wy, lab, 1e3*per, 100*duty, con, thr, cls, ...
    'VariableNames', {'Freq_Hz', 'Capture', 'Why', 'Result', 'Period_ms', 'Duty_pct', 'Contrast_dB', 'NoiseContrast_dB', 'Class'});
end

%% ===== baseline ===============================================================

function B = compareBaseline(R, opt, BP, CAT)
src = string(opt.Baseline);
fprintf('Baseline: %s\n', src);
if isscalar(src) && endsWith(lower(src), ".mat") && ismember('survey', who('-file', src))
    L = load(src, 'survey'); Rb = L.survey;
    if ~isfield(Rb, 'emitters') || ~isfield(Rb, 'raw'), Rb = analyzeSurvey(Rb.raw, opt, BP, CAT); end
else
    o2 = opt; o2.SignalID = "off";
    Rb = analyzeSurvey(loadRaw(src, o2), o2, BP, CAT);
end
cur = R.inventory; base = Rb.inventory;
covered = @(lo, hi) any(arrayfun(@(g) ~isempty(g.f) && lo >= g.f(1) && hi <= g.f(end), R.segments));
% EMI regions (rise / cluster / comb): bursty hash shatters differently in every survey, so a
% narrow emitter inside such a region of the other survey is not "new" or "gone" by itself
regionOf = @(T, f) regionIdx(T, f);
usedB = false(height(base), 1);
rows = {};
for i = 1:height(cur)
    tol = max([2*cur.RBW_Hz(i), 0.25*cur.OBW_Hz(i)]);
    d = abs(base.Freq_Hz - cur.Freq_Hz(i));
    ov = min(base.Hi_Hz, cur.Hi_Hz(i)) - max(base.Lo_Hz, cur.Lo_Hz(i));
    c = find((ov > 0 | d <= tol) & ~usedB);
    lc = levelsAt(R.raw, R.segments, cur.Lo_Hz(i), cur.Hi_Hz(i));
    lb = levelsAt(Rb.raw, Rb.segments, cur.Lo_Hz(i), cur.Hi_Hz(i));
    [p, dl] = levelTest(lc, lb);
    if isempty(c) && regionOf(base, cur.Freq_Hz(i)) > 0
        j = regionOf(base, cur.Freq_Hz(i));
        rows(end + 1, :) = {cur.Freq_Hz(i), cur.OBW_Hz(i), "In EMI region", cur.Class(i), cur.SubType(i), base.Peak_dBm(j), cur.Peak_dBm(i), dl, p}; %#ok<AGROW>
    elseif isempty(c)
        rows(end + 1, :) = {cur.Freq_Hz(i), cur.OBW_Hz(i), "New", cur.Class(i), cur.SubType(i), NaN, cur.Peak_dBm(i), dl, p}; %#ok<AGROW>
    else
        [~, k] = min(d(c)); j = c(k); usedB(j) = true;
        rows(end + 1, :) = {cur.Freq_Hz(i), cur.OBW_Hz(i), "Matched", cur.Class(i), cur.SubType(i), base.Peak_dBm(j), cur.Peak_dBm(i), dl, p}; %#ok<AGROW>
    end
end
for j = find(~usedB).'
    if ~covered(base.Lo_Hz(j), base.Hi_Hz(j)) || regionOf(cur, base.Freq_Hz(j)) > 0, continue, end
    if ismember('Kind', base.Properties.VariableNames) && ismember(base.Kind(j), ["rise" "cluster" "comb"]) && ...
            any(ismember(cur.Kind, ["rise" "cluster" "comb"]) & cur.Lo_Hz <= base.Hi_Hz(j) & cur.Hi_Hz >= base.Lo_Hz(j))
        continue                                            % the EMI region is still there (other shape)
    end
    lc = levelsAt(R.raw, R.segments, base.Lo_Hz(j), base.Hi_Hz(j));
    lb = levelsAt(Rb.raw, Rb.segments, base.Lo_Hz(j), base.Hi_Hz(j));
    [p, dl] = levelTest(lc, lb);
    rows(end + 1, :) = {base.Freq_Hz(j), base.OBW_Hz(j), "Gone", base.Class(j), base.SubType(j), base.Peak_dBm(j), NaN, dl, p}; %#ok<AGROW>
end
vn = {'Freq_Hz', 'OBW_Hz', 'Match', 'Class', 'SubType', 'BasePeak_dBm', 'CurPeak_dBm', 'Delta_dB', 'p'};
if isempty(rows), T = cell2table(cell(0, 9), 'VariableNames', vn); else, T = cell2table(rows, 'VariableNames', vn); end
T.q = fdrBH(T.p);
st = strings(height(T), 1);
for i = 1:height(T)
    sig = T.q(i) < opt.Alpha && abs(T.Delta_dB(i)) >= 3;
    switch T.Match(i)
        case "New",  st(i) = ternary(sig || isnan(T.q(i)), "New", "New (not significant)");
        case "Gone", st(i) = ternary(sig || isnan(T.q(i)), "Gone", "Gone (not significant)");
        otherwise,   st(i) = ternary(sig, "Changed", "Same");
    end
end
T.Status = st;
T = sortrows(T, 'Freq_Hz');
B = struct('source', src, 'table', T, 'baseTime', Rb.meta.time, 'nNew', nnz(st == "New"), 'nGone', nnz(st == "Gone"), ...
    'nChanged', nnz(st == "Changed"), 'nSame', nnz(st == "Same"), 'baseInventory', base);
fprintf('  new %d, gone %d, changed %d, same %d\n', B.nNew, B.nGone, B.nChanged, B.nSame);
end

function j = regionIdx(T, f)
% Row of an EMI region (rise / cluster / comb) of inventory T that contains f, else 0.
j = 0;
if ~ismember('Kind', T.Properties.VariableNames) || height(T) == 0, return, end
k = find(ismember(T.Kind, ["rise" "cluster" "comb"]) & T.Lo_Hz <= f & T.Hi_Hz >= f, 1);
if ~isempty(k), j = k; end
end

function L = levelsAt(raw, segs, lo, hi)
L = [];
for s = 1:numel(raw.seg)
    sg = raw.seg(s);
    if isempty(sg.f) || lo < sg.f(1) || hi > sg.f(end), continue, end
    k = sg.f >= lo - 0.5*(sg.f(2) - sg.f(1)) & sg.f <= hi + 0.5*(sg.f(2) - sg.f(1));
    if ~any(k), [~, i] = min(abs(sg.f - (lo + hi)/2)); k(i) = true; end
    P = double(sg.Pk); if isempty(P), P = double(sg.Av); end
    L = [L; 10*log10(mean(10.^(P(:, k)/10), 2))]; %#ok<AGROW>
end
if isempty(L) && ~isempty(segs), L = []; end
end

function [p, d] = levelTest(a, b)
p = NaN; d = NaN;
if isempty(a) || isempty(b), return, end
d = median(a) - median(b);
if numel(a) >= 3 && numel(b) >= 3, p = ranksum(a, b); end
end

function q = fdrBH(p)
q = nan(size(p)); k = find(isfinite(p));
if isempty(k), return, end
[ps, o] = sort(p(k)); m = numel(ps);
qs = min(1, cummin(ps.*m./(1:m)', 'reverse'));
q(k(o)) = qs;
end

%% ===== watch ==================================================================

function W = watchLoop(io, plan, opt, R0, BP)
% Repeated coarse scans compared with the survey itself (R0): a mask = survey max hold
% (3-point running max) or floor + threshold, whichever is higher, + 'MaskMargin' dB.
% Mask violations -> new / stronger emitters; new combs and broadband rises from the
% pass's own analysis; known steady emitters that fall below the threshold -> gone;
% per-point robust z against the earlier passes and per-segment feature anomalies.
% An alert needs 'Confirm' consecutive passes (anomalies and errors: immediately).
nPass = max(1, floor(opt.Duration/max(opt.Interval, 0.1)));
logf = "";
if opt.Save, logf = fullfile(opt.OutDir, "watch_log.csv"); end
inv0 = R0.inventory;
if isfield(R0.baseline, 'baseInventory') && height(R0.baseline.baseInventory) > 0
    b = R0.baseline.baseInventory; v = intersect(inv0.Properties.VariableNames, b.Properties.VariableNames, 'stable');
    inv0 = [inv0(:, v); b(:, v)];
end
known = inv0(:, intersect({'Freq_Hz', 'Lo_Hz', 'Hi_Hz', 'OBW_Hz', 'RBW_Hz', 'Class', 'SubType', 'Occupancy_pct', 'Kind'}, ...
    inv0.Properties.VariableNames, 'stable'));
if ~ismember('Kind', known.Properties.VariableNames), known.Kind = repmat("signal", height(known), 1); end
known.Low = zeros(height(known), 1); known.GoneSent = false(height(known), 1);
combsKnown = R0.combs.Spacing_Hz(:);
risesKnown = [R0.rises.Lo_Hz(:), R0.rises.Hi_Hz(:)];
segs0 = R0.segments;
H = zeros(0, 6*height(plan));
traceHist = cell(1, height(plan));
maxHist = 200;                          % per-point history: last 200 passes (bounded memory in long watches)
alerts = struct('time', {}, 'pass', {}, 'type', {}, 'freq', {}, 'level', {}, 'message', {});
pend = struct('type', {}, 'lo', {}, 'hi', {}, 'n', {}, 'last', {}, 'sent', {});
passes = table();
fprintf('\nWatch: %d pass(es) every %g s, mask = survey max hold + %g dB, alert after %d consecutive pass(es) (Ctrl+C stops; the instrument is restored)\n', ...
    nPass, opt.Interval, opt.MaskMargin, opt.Confirm);
if nPass > maxHist, fprintf('  (the per-point level history keeps the last %d passes)\n', maxHist); end
for p = 1:nPass
    t0 = tic;
    io.state('pass') = p;
    raw = newRaw(io, plan, opt);
    for s = 1:numel(raw.seg), raw.seg(s).preamp = R0.raw.seg(min(s, numel(R0.raw.seg))).preamp; end
    o2 = opt; o2.PreampMode = "fixed";
    raw = coarseScan(io, raw, plan, o2, 1, opt.WatchSweeps, false);
    if strlength(raw.abortMsg) > 0
        alerts = addAlert(alerts, io.now(), p, "error", NaN, NaN, "scan failed: " + raw.abortMsg, opt, logf);
        if isConnectionError(MException(char(ternary(strlength(raw.abortId) > 0, raw.abortId, "x:y")), '%s', raw.abortMsg)), break, end
        continue
    end
    C = analyzeCoarse(raw, opt);
    cand = struct('type', {}, 'f', {}, 'lo', {}, 'hi', {}, 'lvl', {}, 'msg', {});
    for s = 1:numel(C.seg)
        G = C.seg(s);
        if isempty(G.f) || s > numel(segs0) || numel(segs0(s).f) ~= numel(G.f), continue, end
        G0 = segs0(s);
        % --- mask violations
        mask = max(movmax(G0.pkHold, 3), G0.flP + G0.thrP) + opt.MaskMargin;
        v = G.pkHold - mask;
        hit = v > 0 & (G.pkHold - G.flS) > G.thrS;
        regs = growRegions(hit, hit, max(1, round(0.5*G.rbw/G.df)));
        rows = zeros(0, 5);                                      % f, level, excess, lo, hi
        for r = 1:numel(regs)
            ix = regs{r}; [vx, im] = max(v(ix));
            rows(end + 1, :) = [G.f(ix(im)), G.pkHold(ix(im)), vx, G.f(ix(1)), G.f(ix(end))]; %#ok<AGROW>
        end
        if ~isempty(rows)
            % cluster nearby violations (a bursty or broadband source gives several per sweep)
            [~, o] = sort(rows(:, 1)); rows = rows(o, :);
            cl = cumsum([1; rows(2:end, 4) - rows(1:end - 1, 5) > max(10*G.rbw, 0.03*rows(2:end, 1))]);
            nCl = max(cl);
            if nCl > 4, cl(:) = 1; nCl = 1; end                 % many: one summary for the segment
            for q = 1:nCl
                rr = rows(cl == q, :); [~, im] = max(rr(:, 3));
                f0 = rr(im, 1); lo = min(rr(:, 4)); hi = max(rr(:, 5)); b = bandOf(BP, f0);
                tol = max(2*G.rbw, 2*G.df);
                k = find(known.Lo_Hz - tol <= hi & known.Hi_Hz + tol >= lo, 1);
                extra = ""; if size(rr, 1) > 1, extra = sprintf(" (%d violations over %s-%s)", size(rr, 1), fmtHz(lo), fmtHz(hi)); end
                if ~isempty(k) && (size(rr, 1) == 1 || known.Kind(k) ~= "signal")
                    msg = sprintf("STRONGER: %s (%s) now %.1f dBm, +%.1f dB over the survey max hold%s", fmtHz(f0), ...
                        known.SubType(k), rr(im, 2), rr(im, 3) + opt.MaskMargin, extra); typ = "change";
                else
                    msg = sprintf("NEW emission at %s, %.1f dBm (+%.1f dB over the survey max hold)%s in %s", fmtHz(f0), ...
                        rr(im, 2), rr(im, 3) + opt.MaskMargin, extra, ternary(strlength(b.name) > 0, b.name, "no band-plan entry"));
                    typ = ternary(b.kind == "rnss", "emi", "new");
                    if b.kind == "rnss", msg = msg + " - interference in a GNSS band"; end
                end
                cand(end + 1) = struct('type', typ, 'f', f0, 'lo', lo - tol, 'hi', hi + tol, 'lvl', rr(im, 2), 'msg', msg); %#ok<AGROW>
            end
        end
        % --- known steady emitters that disappeared (below floor + threshold)
        for k = find(known.Freq_Hz >= G.f(1) & known.Freq_Hz <= G.f(end) & known.Kind == "signal" & ...
                known.Occupancy_pct >= 80 & ~known.GoneSent).'
            ix = G.f >= known.Lo_Hz(k) - G.df & G.f <= known.Hi_Hz(k) + G.df;
            if ~any(ix), continue, end
            low = max(G.pkHold(ix) - G.flP(ix)) < G.thrP && max(G.avMean(ix) - G.flA(ix)) < G.thrA;   % pulsed: the average still shows it
            known.Low(k) = low*(known.Low(k) + 1);
            if known.Low(k) >= max(2, opt.Confirm)
                known.GoneSent(k) = true;
                alerts = addAlert(alerts, io.now(), p, "gone", known.Freq_Hz(k), NaN, sprintf("emitter at %s (%s) below the threshold for %d passes", ...
                    fmtHz(known.Freq_Hz(k)), known.SubType(k), known.Low(k)), opt, logf);
            end
        end
        % --- level anomaly of the (9-point smoothed) average trace vs the earlier watch
        % passes: robust z (median/MAD), runs of >= 5 points, once there are 5 passes
        y = movmean(G.avMean, 9);
        Hs = traceHist{s};
        if ~isempty(Hs) && size(Hs, 2) ~= numel(G.f), Hs = []; end
        if size(Hs, 1) >= 5
            md = median(Hs, 1); sd = max(1.4826*median(abs(Hs - md), 1), 1.5);
            z = (y - md)./sd;
            [a1, a2] = runs(abs(z) > 8);
            keep = (a2 - a1 + 1) >= 5; a1 = a1(keep); a2 = a2(keep);
            for r = 1:min(3, numel(a1))
                [~, im] = max(abs(z(a1(r):a2(r)))); fi = G.f(a1(r) + im - 1);
                if any([cand.lo] <= G.f(a2(r)) & [cand.hi] >= G.f(a1(r))), continue, end
                alerts = addAlert(alerts, io.now(), p, "anomaly", fi, y(a1(r) + im - 1), sprintf("average level %s by %.0f dB over %s-%s (robust z %.0f vs %d earlier passes)", ...
                    ternary(z(a1(r) + im - 1) > 0, "up", "down"), abs(y(a1(r) + im - 1) - md(a1(r) + im - 1)), ...
                    fmtHz(G.f(a1(r))), fmtHz(G.f(a2(r))), max(abs(z(a1(r):a2(r)))), size(Hs, 1)), opt, logf);
            end
        end
        traceHist{s} = [Hs(max(1, end - maxHist + 2):end, :); y];
    end
    % --- new combs / broadband rises from this pass's own analysis
    for e = C.emitters
        if e.kind == "comb" && ~any(abs(combsKnown - e.combSp) <= 0.02*e.combSp)
            cand(end + 1) = struct('type', "comb", 'f', e.combSp, 'lo', 0.98*e.combSp, 'hi', 1.02*e.combSp, 'lvl', e.pk, ...
                'msg', sprintf("NEW EMI: harmonic comb, spacing %s, %d lines from %s", fmtHz(e.combSp), e.lines, fmtHz(e.lo))); %#ok<AGROW>
        elseif e.kind == "rise" && ~any(risesKnown(:, 1) <= e.hi & risesKnown(:, 2) >= e.lo)
            cand(end + 1) = struct('type', "rise", 'f', e.fc, 'lo', e.lo, 'hi', e.hi, 'lvl', e.pk, ...
                'msg', sprintf("NEW EMI: broadband noise rise +%.1f dB over %s-%s", e.exc, fmtHz(e.lo), fmtHz(e.hi))); %#ok<AGROW>
        end
    end
    % --- confirmation over consecutive passes, then alert once
    % (same kind of change, overlapping frequency range, in consecutive passes)
    if ~isempty(cand)
        isR = [cand.type] == "rise"; drop = false(1, numel(cand));
        for q = find(~isR & [cand.type] ~= "comb")          % a new hump is reported once, as a rise
            drop(q) = any([cand(isR).lo] <= cand(q).hi & [cand(isR).hi] >= cand(q).lo);
        end
        cand(drop) = [];
    end
    fam = @(t) ternary(ismember(t, ["new" "emi" "change"]), "mask", t);
    for c = cand
        i = [];
        for q = 1:numel(pend)
            if fam(pend(q).type) == fam(c.type) && pend(q).lo <= c.hi && pend(q).hi >= c.lo, i = q; break, end
        end
        if isempty(i), pend(end + 1) = struct('type', c.type, 'lo', c.lo, 'hi', c.hi, 'n', 0, 'last', 0, 'sent', false); i = numel(pend); end %#ok<AGROW>
        if pend(i).last == p, continue, end                     % already counted in this pass
        if pend(i).last == p - 1, pend(i).n = pend(i).n + 1; else, pend(i).n = 1; end
        pend(i).last = p; pend(i).lo = min(pend(i).lo, c.lo); pend(i).hi = max(pend(i).hi, c.hi);
        if pend(i).n >= opt.Confirm && ~pend(i).sent
            pend(i).sent = true;
            msg = c.msg; if opt.Confirm > 1, msg = msg + sprintf(" [%d passes]", opt.Confirm); end
            alerts = addAlert(alerts, io.now(), p, ternary(ismember(c.type, ["comb" "rise"]), "emi", c.type), ...
                ternary(c.type == "comb", NaN, c.f), c.lvl, msg, opt, logf);
        end
    end
    % --- per-segment feature anomalies
    F = segFeatures(C.seg);
    [isA, why] = featureAnomaly(H, F, opt);
    for s = find(isA(:).')
        alerts = addAlert(alerts, io.now(), p, "anomaly", NaN, NaN, sprintf("segment %s: %s", C.seg(s).name, why(s)), opt, logf);
    end
    H = [H; F]; %#ok<AGROW>
    na = nnz([alerts.pass] == p);
    passes = [passes; table(p, io.now(), numel(C.emitters), numel(cand), na, ...
        'VariableNames', {'Pass', 'Time', 'Emitters', 'Candidates', 'Alerts'})]; %#ok<AGROW>
    fprintf('  pass %d/%d: %d emitters, %d candidate change(s), %d alert(s), %.1f s\n', p, nPass, numel(C.emitters), numel(cand), na, toc(t0));
    if io.kind ~= "demo" && p < nPass, pause(max(0, opt.Interval - toc(t0))); end
end
if isempty(alerts)
    A = table(NaT(0, 1), zeros(0, 1), strings(0, 1), zeros(0, 1), zeros(0, 1), strings(0, 1), ...
        'VariableNames', {'Time', 'Pass', 'Type', 'Freq_Hz', 'Level_dBm', 'Message'});
else
    A = table([alerts.time].', [alerts.pass].', [alerts.type].', [alerts.freq].', [alerts.level].', [alerts.message].', ...
        'VariableNames', {'Time', 'Pass', 'Type', 'Freq_Hz', 'Level_dBm', 'Message'});
end
unconf = sum(~[pend.sent]);
W = struct('alerts', A, 'passes', passes, 'known', known, 'features', H, 'log', logf, 'unconfirmed', unconf);
end

function tf = coveredBy(C, f)
tf = false(size(f));
for s = 1:numel(C.seg)
    if isempty(C.seg(s).f), continue, end
    tf = tf | (f >= C.seg(s).f(1) & f <= C.seg(s).f(end));
end
end

function F = segFeatures(segs)
F = nan(numel(segs), 6);
for s = 1:numel(segs)
    G = segs(s);
    if isempty(G.f), continue, end
    e = G.pkHold - G.flP;
    F(s, :) = [median(G.flA), prctile(G.flA, 90) - prctile(G.flA, 10), mean(G.hitFrac)*100, ...
        10*log10(sum(10.^(max(G.avMean - G.flA, 0)/10)) + eps), max(e), max(G.flA - G.refA)];
end
F = reshape(F.', 1, []);                                   % one row: all segments
end

function [isA, why] = featureAnomaly(H, F, opt)
nS = numel(F)/6; isA = false(1, nS); why = strings(1, nS);
names = ["floor median" "floor spread" "occupancy" "excess power" "max excess" "floor rise"];
if size(H, 1) < 5, return, end
% floors: bursty emitters (Wi-Fi, radar, hash) move occupancy and max excess a lot between
% single-sweep passes - only large, persistent-looking jumps count
md = median(H, 1, 'omitnan'); sd = max(1.4826*median(abs(H - md), 1, 'omitnan'), repmat([0.75 1 4 2 6 1], 1, nS));
z = (F - md)./sd;
for s = 1:nS
    zs = z((s - 1)*6 + (1:6));
    [m, i] = max(abs(zs));
    if m > 8, isA(s) = true; why(s) = sprintf("%s changed (robust z %.1f)", names(i), zs(i)); end
end
mode = opt.Anomaly;
ok = all(isfinite(H), 1) & all(isfinite(F)) & std(H, 0, 1) > 0;
d = nnz(ok);
try
    if (mode == "robustcov" || (mode == "auto" && size(H, 1) >= 4*d)) && d >= 2 && size(H, 1) > 2*d
        Z = (H(:, ok) - md(ok))./sd(ok); zf = (F(ok) - md(ok))./sd(ok);
        [S, mu] = robustcov(Z);
        d2 = (zf - mu)*pinv(S)*(zf - mu).';
        if d2 > chi2inv(0.999, d) && ~any(isA)
            isA(1) = true; why(1) = sprintf("joint feature anomaly (robust Mahalanobis^2 %.0f > %.0f)", d2, chi2inv(0.999, d));
        end
    elseif mode == "iforest" && size(H, 1) >= 8 && d >= 2
        mdl = iforest(H(:, ok), 'ContaminationFraction', 0.01);
        [tf, s] = isanomaly(mdl, F(ok));
        if tf && ~any(isA), isA(1) = true; why(1) = sprintf("isolation-forest anomaly (score %.2f)", s); end
    end
catch
end
end

function alerts = addAlert(alerts, t, p, typ, f, lvl, msg, opt, logf)
alerts(end + 1) = struct('time', t, 'pass', p, 'type', string(typ), 'freq', f, 'level', lvl, 'message', string(msg));
fprintf(2, '  [%s] ALERT %-7s %s\n', string(t, 'HH:mm:ss'), upper(typ), msg);
if strlength(logf) > 0
    fid = -1;
    try
        newf = ~isfile(logf);
        fid = fopen(logf, 'a');
        if newf, fprintf(fid, 'time,pass,type,freq_Hz,level_dBm,message\n'); end
        fprintf(fid, '%s,%d,%s,%.6f,%.2f,"%s"\n', string(t, 'yyyy-MM-dd HH:mm:ss'), p, typ, f, lvl, strrep(msg, '"', ''''));
    catch
    end
    if fid >= 0, fclose(fid); end                        % also after a write error
end
if opt.Beep, beep; end
if ~isempty(opt.AlertFcn)
    try
        feval(opt.AlertFcn, alerts(end));
    catch err
        warning('ff_auto_survey:alertfcn', 'AlertFcn failed: %s', err.message);
    end
end
end

%% ===== zero-span helpers (from the suite's verified research) =================

function r = tddAnalyze(t, p, slotDur)
% Classify a zero-span envelope: TDD (periodic on/off), continuous, bursty, intermittent.
if nargin < 3, slotDur = 0.5e-3; end
p = p(:); t = t(:); dt = t(2) - t(1); T = t(end) - t(1) + dt;
w = max(1, round(0.1e-3/dt));
ps = 10*log10(movmean(10.^(p/10), w));
lo = prctile(ps, 10); hi = prctile(ps, 90); pk99 = prctile(ps, 99.5);
r.lo = lo; r.hi = hi; r.pk = pk99;
r.bursty = (hi - lo) < 6 && (pk99 - lo) > 10;
if r.bursty, r.thr = lo + (pk99 - lo)/2; else, r.thr = (lo + hi)/2; end
on = ps > r.thr; r.onFrac = mean(on);
if r.bursty, x = double(on) - mean(on); else, x = ps - mean(ps); end
maxLag = min(numel(x) - 1, round(min(25e-3, 0.6*T)/dt));
r.period = NaN; r.periodicity = 0;
if maxLag > 3 && any(x ~= 0)
    [ac, lags] = xcorr(x, maxLag, 'coeff'); keep = lags*dt >= 0.9e-3; ac = ac(keep); lags = lags(keep);
    if numel(ac) > 3
        [pk, loc] = findpeaks(ac, 'MinPeakProminence', 0.05);
        if ~isempty(pk)
            [pmax, im] = max(pk); L = lags(loc(im));
            for k = 6:-1:2
                j = find(abs(lags(loc) - L/k) <= max(2, 0.03*L/k));
                if ~isempty(j) && max(pk(j)) >= 0.85*pmax, L = lags(loc(j(1))); break; end
            end
            r.period = L*dt; r.periodicity = pmax;
        end
    end
end
d = diff([0; on; 0]); st = find(d == 1); en = find(d == -1) - 1;
if isempty(st), r.medianOnRun = 0; else, r.medianOnRun = median(en - st + 1)*dt; end
r.slotPattern = ''; r.slotOnFrac = [];
if ~isnan(r.period) && ~r.bursty && (hi - lo) > 6
    d0 = diff([1; on; 1]); so = find(d0 == -1); eo = find(d0 == 1) - 1;
    if ~isempty(so)
        [~, ii] = max(eo - so); k0 = eo(ii) + 1; nP = round(r.period/dt); k0 = mod(k0 - 1, nP) + 1;
        nS = max(1, round(r.period/slotDur)); str = blanks(nS);
        r.slotOnFrac = zeros(1, nS);
        for k = 1:nS
            idx = [];
            for m = 0:floor((numel(on) - k0 + 1)/nP) - 1
                idx = [idx, k0 + m*nP + round(((k - 1)*slotDur)/dt) + (0:max(0, round(slotDur/dt) - 1))]; %#ok<AGROW>
            end
            idx = idx(idx <= numel(on)); if isempty(idx), continue, end
            fr = mean(on(idx)); r.slotOnFrac(k) = fr;
            if fr > 0.9, str(k) = 'D'; elseif fr > 0.25, str(k) = 'S'; else, str(k) = 'U'; end
        end
        r.slotPattern = strtrim(str);
    end
end
dyn = hi - lo; r.dyn = dyn;
if r.bursty && r.periodicity >= 0.3
    r.class = sprintf('bursty periodic (SSB-only idle NR / beacon), period %.1f ms, burst %.2f ms', r.period*1e3, r.medianOnRun*1e3);
elseif r.bursty
    r.class = 'sporadic bursts';
elseif dyn < 6
    r.class = 'continuous (FDD DL / always-on) - or noise only if level ~ floor';
elseif r.periodicity >= 0.5 && r.onFrac >= 0.4
    r.class = sprintf('TDD periodic, period %.2f ms, DL duty %.0f%%', r.period*1e3, 100*r.onFrac);
else
    r.class = 'intermittent (traffic-dependent)';
end
end

function thr = noiseContrast(N, fs, rbw)
% 99th percentile of the statelevels contrast of a noise-only SAMP envelope (Monte Carlo).
persistent cache
if isempty(cache), cache = containers.Map(); end
key = sprintf('%d_%.6g_%.6g', N, fs, rbw);
if isKey(cache, key), thr = cache(key); return, end
old = rng(1); c0 = onCleanup(@() rng(old)); %#ok<NASGU>
sig = 0.265*fs/rbw;
if sig > 0.5, h = exp(-0.5*((-ceil(4*sig):ceil(4*sig))/sig).^2); else, h = 1; end
c = zeros(200, 1);
for k = 1:numel(c)
    n = conv(complex(randn(1, N + numel(h)), randn(1, N + numel(h))), h, 'valid');
    lv = statelevels(abs(n(1:N)));
    c(k) = 20*log10(lv(2)/lv(1));
end
thr = prctile(c, 99);
if cache.Count >= 500, remove(cache, keys(cache)); end   % bounded (keys are continuous); recomputed identically
cache(key) = thr;
end

%% ===== demo world =============================================================

function [W, t0, note] = demoWorld(opt)
em = @(name, f0, B, P, shape, temporal, duty, from) struct('name', string(name), 'f0', f0, 'B', B, 'P', P, ...
    'shape', string(shape), 'temporal', string(temporal), 'duty', duty, 'from', from, 'sp', NaN, 'fEnd', NaN, 'slope', 0, 'rise', NaN);
W = [ ...
    em("AM 1.08 MHz", 1.08e6, 9e3, -72, "gauss", "cont", 1, 0), em("AM 1.35 MHz", 1.35e6, 9e3, -80, "gauss", "cont", 1, 0), ...
    em("FM 88.5", 88.5e6, 180e3, -42, "gauss", "cont", 1, 0), em("FM 91.0", 91.0e6, 180e3, -55, "gauss", "cont", 1, 0), ...
    em("FM 95.5", 95.5e6, 180e3, -50, "gauss", "cont", 1, 0), em("FM 99.0", 99.0e6, 180e3, -58, "gauss", "cont", 1, 0), ...
    em("FM 103.5", 103.5e6, 180e3, -52, "gauss", "cont", 1, 0), em("FM 107.0", 107.0e6, 180e3, -60, "gauss", "cont", 1, 0), ...
    em("Air band 125.1 MHz", 125.1e6, 8e3, -85, "gauss", "bursty", 0.4, 0), ...
    em("Meteor-M2-3 LRPT", 137.9e6, 110e3, -92, "rect", "cont", 1, 0), ...
    em("ORBCOMM 137.46", 137.46e6, 6e3, -96, "gauss", "bursty", 0.5, 0), ...
    em("KNACKSAT-2 TLM", 400.63e6, 12.5e3, -88, "gauss", "cont", 1, 0), ...
    em("SRD 433.92", 433.92e6, 40e3, -80, "gauss", "bursty", 0.2, 0), ...
    em("DVB-T2 530 MHz", 530e6, 7.6e6, -60, "rect", "cont", 1, 0), em("DVB-T2 578 MHz", 578e6, 7.6e6, -66, "rect", "cont", 1, 0), ...
    em("LTE B28 DL True 768", 768e6, 18e6, -65, "rect", "cont", 1, 0), em("LTE B8 DL 950", 950e6, 9e6, -62, "rect", "cont", 1, 0), ...
    em("GNSS-band interferer", 1575.42e6, 1, -85, "gauss", "cont", 1, 0), ...
    em("LTE B3 DL 1842.5", 1842.5e6, 18e6, -58, "rect", "cont", 1, 0), em("LTE B1 DL 2155", 2155e6, 13.5e6, -63, "rect", "cont", 1, 0), ...
    em("NR n41 AIS 2550", 2550e6, 98.28e6, -48, "rect", "tdd", 0.743, 0), ...
    em("Wi-Fi ch6 2437", 2437e6, 17e6, -55, "rect", "bursty", 0.3, 0), ...
    em("S-band radar 2800", 2800e6, 1.5e6, -35, "rect", "pulsed", 0.001, 0), ...
    em("Wi-Fi 5 GHz ch36", 5180e6, 17e6, -70, "rect", "bursty", 0.2, 0), ...
    em("New emitter 868 MHz", 868e6, 150e3, -72, "gauss", "bursty", 0.6, 2)];
c = em("Inverter switching comb", 20e3, 1, -62, "comb", "mains", 0.3, 0); c.sp = 20e3; c.f0 = 140e3; c.fEnd = 30e6; c.slope = 20;
h = em("Inverter broadband EMI", 70e6, 90e6, NaN, "hump", "impulsive", 1, 0); h.rise = 8;
h2 = em("New EMI hump 250 MHz", 250e6, 60e6, NaN, "hump", "impulsive", 1, 3); h2.rise = 10;
W = [W, c, h, h2];
t0 = datetime('now');
note = "SYNTHETIC DEMO DATA - not a measurement.";
% put the demo at the next good KNACKSAT-2 pass over the site, with its real Doppler
try
    [tle, ~] = loadTLE(67683, opt);
    if ~isempty(tle)
        tq = dateshift(datetime('now', 'TimeZone', 'UTC'), 'start', 'minute') + minutes(0:1440); tq.TimeZone = '';
        r = propagateOrbit(tq, tle, "PropModel", "sgp4");
        aer = eci2aer(r.', datevec(tq), repmat(opt.Location, numel(tq), 1));
        [elMax, i] = max(aer(:, 2));
        if elMax > 10
            tp = tq(i); tp.TimeZone = 'UTC'; tl = tp; tl.TimeZone = 'local'; tl.TimeZone = '';
            t0 = tl;
            r3 = propagateOrbit(tq(i) + seconds([-1 1]), tle, "PropModel", "sgp4");
            a3 = eci2aer(r3.', datevec(tq(i) + seconds([-1 1])), repmat(opt.Location, 2, 1));
            dop = -400.63e6*(a3(2, 3) - a3(1, 3))/2/299792458;
            k = find([W.name] == "KNACKSAT-2 TLM"); W(k).f0 = W(k).f0 + dop;
            note = note + sprintf(" Demo clock set to the KNACKSAT-2 pass at %s UTC (el %.0f deg, Doppler %+.1f kHz).", ...
                string(tq(i), 'yyyy-MM-dd HH:mm'), elMax, dop/1e3);
        end
    end
catch
end
end

function A = demoSA(W, cfg, det, n, st)
N = cfg.Points; f = linspace(cfg.Start, cfg.Stop, N); df = (cfg.Stop - cfg.Start)/max(N - 1, 1);
rbw = snapRBW(cfg.RBW); pass = st('pass');
noiseP = demoNoise(f, rbw, cfg.Preamp);
fe = f + df*[-0.4; -0.2; 0; 0.2; 0.4];                        % 5 sub-points per bin
isPk = ~isAvgDet(det);
tSweep = 0.04 + 6e-5*N;
P = zeros(n, N, 'single');
for k = 1:n
    if isPk, nz = max(exprnd(1, 3, N), [], 1)/1.83; else, nz = mean(exprnd(1, 16, N), 1); end
    lin = noiseP.*nz;
    for w = W
        if w.from > pass, continue, end
        [lo, hi] = reach(w, rbw);
        if hi < f(1) - df || lo > f(end) + df, continue, end
        c = demoShape(w, fe, rbw, cfg.Preamp);
        if isPk, c = max(c, [], 1); else, c = mean(c, 1); end
        switch w.temporal
            case "cont", g = 10^(0.6*randn/10);
            case {"tdd", "bursty", "mains"}, g = double(rand(1, N) < w.duty);
            case "pulsed"
                if isPk, g = double(rand(1, N) < min(1, (tSweep/N)/1e-3)); else, g = w.duty; end
            case "impulsive"
                if isPk, g = 2.5*(1 + 0.3*randn(1, N)); else, g = 1; end
            otherwise, g = 1;
        end
        lin = lin + c.*max(g, 0);
    end
    P(k, :) = 10*log10(lin);
end
lev = max(P(:));
over = (cfg.Preamp && lev > -30) || lev > -5;
st2 = struct('start', f(1), 'stop', f(end), 'points', N, 'rbw', rbw, 'vbw', rbw, 'ref', -20, 'pdiv', 10, ...
    'atten', 0, 'preamp', logical(cfg.Preamp), 'detector', string(det));
A = struct('f', f, 'P', P, 't', repmat(datetime('now'), n, 1), 'rbw', rbw, 'det', string(det), ...
    'preamp', logical(cfg.Preamp), 'over', over, 'settings', st2);
end

function noiseP = demoNoise(f, rbw, preamp)
danl = ternary(preamp, -155, -140) + 1.0*(f/1e9);           % dBm/Hz, worse at high frequency
fam = -174 + 72.5 - 27.7*log10(max(f, 1e5)/1e6) - 45;        % man-made noise (P.372, short inefficient antenna)
noiseP = (10.^(danl/10) + 10.^(fam/10))*1.0645*rbw;
end

function [lo, hi] = reach(w, rbw)
switch w.shape
    case "comb", lo = w.f0 - 5*rbw; hi = w.fEnd + 5*rbw;
    case "hump", lo = w.f0 - 1.5*w.B; hi = w.f0 + 1.5*w.B;
    otherwise, lo = w.f0 - w.B - 6*rbw; hi = w.f0 + w.B + 6*rbw;
end
end

function c = demoShape(w, fe, rbw, preamp)
sp = rbw/1.6651;
Pt = 10^(w.P/10);
switch w.shape
    case "rect"
        c = Pt/w.B*(sqrt(pi)*sp/2)*(erf((fe - w.f0 + w.B/2)/sp) - erf((fe - w.f0 - w.B/2)/sp));
    case "gauss"
        ss = max(w.B, 1)/5.15; v = sp^2 + 2*ss^2;
        c = Pt*sp/sqrt(v)*exp(-(fe - w.f0).^2/v);
    case "comb"
        c = zeros(size(fe)); M = round(fe/w.sp); m1 = ceil(w.f0/w.sp); m2 = floor(w.fEnd/w.sp);
        for d = -1:1
            mm = M + d; ok = mm >= m1 & mm <= m2;
            L = w.P - w.slope*log10(max(mm*w.sp, w.f0)/w.f0);
            c = c + ok.*10.^(L/10).*exp(-((fe - mm*w.sp)/sp).^2);
        end
    case "hump"                                             % noise rise of w.rise dB at the centre
        ss = w.B/5.15;
        c = (10^(w.rise/10) - 1)*demoNoise(fe, rbw, preamp).*exp(-(fe - w.f0).^2/(2*ss^2));
    otherwise
        c = zeros(size(fe));
end
end

function Z = demoZS(W, cfg, n, st)
N = cfg.Points; T = cfg.T; t = linspace(0, T, N); dt = T/max(N - 1, 1);
rbw = min(snapRBW(cfg.RBW), 5e6); fc = cfg.Center; pass = st('pass');
nP = demoNoise(fc, rbw, cfg.Preamp);
K = max(1, round(dt*rbw));
E = zeros(n, N, 'single');
for k = 1:n
    switch cfg.Det
        case "SAMP", nz = exprnd(1, 1, N);
        case "POS",  nz = max(exprnd(1, min(K, 50), N), [], 1);
        otherwise,   nz = mean(exprnd(1, min(K, 50), N), 1);
    end
    lin = nP*nz;
    for w = W
        if w.from > pass, continue, end
        [lo, hi] = reach(w, rbw);
        if fc < lo || fc > hi, continue, end
        c0 = demoShape(w, fc, rbw, cfg.Preamp);
        ph = rand;
        switch w.temporal
            case "tdd"                                          % DDDSU, 2.5 ms, S = 10 DL symbols
                ts = mod(t + ph*2.5e-3, 2.5e-3); sl = floor(ts/0.5e-3);
                g = double(sl <= 2 | (sl == 3 & mod(ts, 0.5e-3) < 0.5e-3*10/14));
            case "mains"
                g = double(mod(t + ph*10e-3, 10e-3) < 3e-3).*(1 + 0.3*sin(2*pi*20e3*t));
            case "bursty"
                g = double(mod(t + ph*4e-3, 4e-3) < 4e-3*w.duty);
                g = g.*double(rand(1, N) < 0.9 | g == 0);
            case "pulsed"
                if cfg.Det == "SAMP", g = double(mod(t + ph*1e-3, 1e-3) < 1e-6);
                else, g = double(floor((t + ph*1e-3)/1e-3) ~= floor((t + ph*1e-3 - dt)/1e-3))*min(1, 1e-6*rbw*3); end
            case "impulsive"
                g = 1 + 20*double(rand(1, N) < 0.01);
            otherwise
                g = ones(1, N);
        end
        if w.shape == "rect" && rbw < w.B
            if cfg.Det == "SAMP", g = g.*exprnd(1, 1, N); else, g = g.*mean(exprnd(1, min(K, 20), N), 1); end
        end
        lin = lin + c0*g;
    end
    E(k, :) = 10*log10(lin);
end
Z = struct('T', T, 'rbw', rbw, 'E', E, 'errors', "");
end

%% ===== outputs ================================================================

function files = saveOutputs(R, opt)
files = struct('mat', "", 'inventory', "", 'occupancy', "", 'figures', strings(0, 1), 'report', "", 'recordings', strings(0, 1), 'watchLog', "");
d = opt.OutDir;
survey = R; survey.files = struct(); %#ok<NASGU>
files.mat = fullfile(d, "survey.mat");
try
    save(files.mat, 'survey', '-v7.3');
    if isfile(fullfile(d, "survey_partial.mat")), delete(fullfile(d, "survey_partial.mat")); end
catch err
    warning('ff_auto_survey:save', 'survey.mat not written: %s', err.message);
end
try
    files.inventory = fullfile(d, "inventory.csv"); writetable(R.inventory, files.inventory);
    files.occupancy = fullfile(d, "band_occupancy.csv"); writetable(R.bandOccupancy, files.occupancy);
catch err
    warning('ff_auto_survey:save', 'CSV not written: %s', err.message);
end
if ~startsWith(R.raw.source, ["survey file" "recordings"])
    files.recordings = saveRecordings(R.raw, d);
end
if isfield(R.watch, 'log'), files.watchLog = R.watch.log; end
tF = tic;
try
    files.figures = makeFigures(R, opt);
catch err
    warning('ff_auto_survey:figures', 'Figures failed: %s', err.message);
end
fprintf('  figures %.1f s\n', toc(tF));
if opt.Report ~= "none"
    tF = tic;
    try
        files.report = makeReport(R, opt, files);
    catch err
        warning('ff_auto_survey:report', 'Report failed: %s', err.message);
    end
    fprintf('  report %.1f s\n', toc(tF));
end
fprintf('Results in %s\n', d);
end

function out = saveRecordings(raw, d)
out = strings(0, 1);
rd = fullfile(d, "rec"); if ~isfolder(rd), mkdir(rd); end
for s = 1:numel(raw.seg)
    sg = raw.seg(s);
    for det = ["POS" "AVER"]
        if det == "POS", P = sg.Pk; t = sg.tPk; else, P = sg.Av; t = sg.tAv; end
        if isempty(P), continue, end
        f = double(sg.f); P = double(P); idn = raw.idn; %#ok<NASGU>
        settings = struct('start', f(1), 'stop', f(end), 'points', numel(f), 'rbw', sg.rbw, 'preamp', sg.preamp, ...
            'detector', det, 'band', sg.name, 'overload', sg.over, 'tool', "ff_auto_survey"); %#ok<NASGU>
        file = fullfile(rd, sprintf("rec_seg%d_%s.mat", s, det));
        save(file, 'f', 'P', 't', 'settings', 'idn');
        out(end + 1, 1) = file; %#ok<AGROW>
    end
end
end

function files = makeFigures(R, opt)
files = strings(0, 1);
if strcmpi(opt.Visible, 'off')
    fig0 = newFig(opt, [40 40 800 600]); c0 = onCleanup(@() delete(fig0(isgraphics(fig0)))); %#ok<NASGU>
end
d = R.outDir; if strlength(d) == 0, d = opt.OutDir; end
segs = R.segments(arrayfun(@(g) ~isempty(g.f), R.segments));
if isempty(segs), return, end
inv = R.inventory;
cc = classColors();
% --- overview
fig = newFig(opt, [40 40 1500 640]);
ax = axes(fig); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on'); set(ax, 'XScale', 'log');
for s = 1:numel(segs)
    G = segs(s);
    [x, y] = decim(G.f, G.pkHold, 2500, "max"); plot(ax, x, y, 'Color', [0.72 0.8 0.95]);
    [x, y] = decim(G.f, G.avMean, 2500, "mean"); plot(ax, x, y, 'Color', [0.1 0.25 0.55]);
    [x, y] = decim(G.f, G.flA + G.thrA, 500, "mean"); plot(ax, x, y, ':', 'Color', [0.35 0.35 0.35]);
end
xlim(ax, [max(segs(1).f(1), 1e3), segs(end).f(end)]);
yl = ylim(ax); yl(2) = yl(2) + 12; ylim(ax, yl);
bandStrip(ax, R, yl);
for c = 1:4
    k = inv.Class == cc.names(c);
    if any(k), scatter(ax, inv.Freq_Hz(k), inv.Peak_dBm(k), 40, cc.rgb(c, :), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', cc.names(c)); end
end
[~, o] = sort(inv.Excess_dB, 'descend');
for i = o(1:min(end, 14)).'
    text(ax, inv.Freq_Hz(i), inv.Peak_dBm(i) + 2, shortLabel(inv(i, :)), 'FontSize', 7, 'Rotation', 35, 'Interpreter', 'none');
end
xlabel(ax, 'Frequency (Hz, log)'); ylabel(ax, 'Level (dBm in the RBW of each segment)');
title(ax, sprintf('Spectrum survey %s - %d emitters (light: max hold, dark: average, dotted: threshold)%s', ...
    string(R.meta.time, 'yyyy-MM-dd HH:mm'), height(inv), ternary(R.meta.synthetic, '  [SYNTHETIC]', '')));
h = findobj(ax, 'Type', 'Scatter'); if ~isempty(h), legend(ax, h, 'Location', 'southwest'); end
files(end + 1) = saveFig(fig, fullfile(d, "fig_overview.png"), opt);
% --- segments
fig = newFig(opt, [40 40 1400 230*numel(segs) + 60]);
tl = tiledlayout(fig, numel(segs), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:numel(segs)
    G = segs(s); ax = nexttile(tl); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
    sc = 1e-6;
    [x, y] = decim(G.f*sc, G.pkHold, 3000, "max"); plot(ax, x, y, 'Color', [0.72 0.8 0.95]);
    [x, y] = decim(G.f*sc, G.avMean, 3000, "mean"); plot(ax, x, y, 'Color', [0.1 0.25 0.55]);
    [x, y] = decim(G.f*sc, G.flA, 500, "mean"); plot(ax, x, y, '--', 'Color', [0.3 0.6 0.3]);
    plot(ax, x, y + G.thrA, ':', 'Color', [0.4 0.4 0.4]);
    xlim(ax, [G.f(1) G.f(end)]*sc); yl = ylim(ax);
    k = inv.Freq_Hz >= G.f(1) & inv.Freq_Hz <= G.f(end) & inv.Segment == G.name;
    for ci = 1:4                                           % one patch object per class (fast export)
        i = find(k & inv.Class == cc.names(ci));
        if isempty(i), continue, end
        X = [inv.Lo_Hz(i) inv.Hi_Hz(i) inv.Hi_Hz(i) inv.Lo_Hz(i)].'*sc;
        patch(ax, X, repmat(yl([1 1 2 2]).', 1, numel(i)), cc.rgb(ci, :), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end
    ylim(ax, yl);
    title(ax, sprintf('%s  RBW %s  preamp %s  floor %.1f dBm  threshold %.1f dB', G.name, fmtHz(G.rbw), ternary(G.preamp, "on", "off"), median(G.flA), G.thrA), 'Interpreter', 'none');
    ylabel(ax, 'dBm');
end
xlabel(tl, 'Frequency (MHz)');
files(end + 1) = saveFig(fig, fullfile(d, "fig_segments.png"), opt);
% --- waterfall
fig = newFig(opt, [40 40 1400 200*numel(segs) + 60]);
tl = tiledlayout(fig, numel(segs), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:numel(segs)
    sg = R.raw.seg(arrayfun(@(x) x.name == segs(s).name, R.raw.seg));
    ax = nexttile(tl);
    if isempty(sg), continue, end
    P = double(sg(1).Pk); if isempty(P), P = double(sg(1).Av); end
    imagesc(ax, segs(s).f*1e-6, 1:size(P, 1), P); axis(ax, 'xy'); colormap(ax, parula);
    clim(ax, [median(segs(s).flA) - 3, max(max(P(:)), median(segs(s).flA) + 10)]); cb = colorbar(ax); cb.Label.String = 'dBm';
    title(ax, sprintf('%s - sweeps in time order (peak detector)', segs(s).name), 'Interpreter', 'none'); ylabel(ax, 'sweep');
end
xlabel(tl, 'Frequency (MHz)');
files(end + 1) = saveFig(fig, fullfile(d, "fig_waterfall.png"), opt);
% --- zoom thumbnails
Z = R.zoom;
if ~isempty(Z)
    n = min(16, numel(Z)); nc = min(4, n); nr = ceil(n/nc);
    fig = newFig(opt, [40 40 360*nc 250*nr + 40]);
    tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:n
        z = Z(i); ax = nexttile(tl); hold(ax, 'on'); box(ax, 'on');
        sc = 1e-6; plot(ax, z.f*sc, z.x, 'Color', [0.1 0.25 0.55]); plot(ax, z.f*sc, z.fl, '--', 'Color', [0.3 0.6 0.3]);
        if isfinite(z.lo), xline(ax, [z.lo z.hi]*sc, 'r:'); end
        lab = "";
        if z.emitter > 0, e = R.emitters(z.emitter); lab = e.cls + ": " + e.sub; end
        title(ax, {sprintf('%s  RBW %s  OBW %s', fmtHz(z.target), fmtHz(z.rbw), fmtHz(z.obw)), char(extractBefore(lab + " ", min(60, strlength(lab) + 1)))}, ...
            'FontSize', 7, 'Interpreter', 'none');
        xlim(ax, [z.f(1) z.f(end)]*sc);
    end
    xlabel(tl, 'MHz'); ylabel(tl, 'dBm');
    files(end + 1) = saveFig(fig, fullfile(d, "fig_zoom.png"), opt);
end
% --- zero span
ZS = R.zsDetail;
if ~isempty(ZS)
    n = min(8, numel(ZS)); nc = min(2, n); nr = ceil(n/nc);
    fig = newFig(opt, [40 40 700*nc 220*nr + 40]);
    tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:n
        zz = R.raw.zs(ZS(i).raw); ax = nexttile(tl);
        t = linspace(0, zz.T, size(zz.E, 2))*1e3;
        plot(ax, t, zz.E(1, :), 'Color', [0.1 0.25 0.55]); grid(ax, 'on');
        title(ax, sprintf('%s %s RBW %s: %s', fmtHz(zz.target), zz.det, fmtHz(zz.rbw), ZS(i).res.label), 'FontSize', 7, 'Interpreter', 'none');
        xlim(ax, [0 t(end)]);
    end
    xlabel(tl, 'time (ms)'); ylabel(tl, 'dBm');
    files(end + 1) = saveFig(fig, fullfile(d, "fig_zerospan.png"), opt);
end
% --- occupancy
T = R.bandOccupancy;
if height(T) > 0
    T = T(1:min(end, 25), :);
    fig = newFig(opt, [40 40 1000 30*height(T) + 120]);
    ax = axes(fig); barh(ax, T.Occupancy_pct, 'FaceColor', [0.2 0.45 0.7]); grid(ax, 'on');
    set(ax, 'YTick', 1:height(T), 'YTickLabel', T.Band, 'YDir', 'reverse', 'TickLabelInterpreter', 'none', 'FontSize', 8);
    xlabel(ax, '% of points x sweeps above the detection threshold'); title(ax, 'Band occupancy (busiest first)');
    files(end + 1) = saveFig(fig, fullfile(d, "fig_occupancy.png"), opt);
end
% --- baseline
if isfield(R.baseline, 'table') && height(R.baseline.table) > 0
    B = R.baseline.table;
    fig = newFig(opt, [40 40 1300 450]); ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); set(ax, 'XScale', 'log');
    st = ["New" "Gone" "Changed" "Same"]; col = [0.85 0.2 0.1; 0.3 0.3 0.3; 0.9 0.6 0; 0.6 0.75 0.9];
    for k = 1:4
        m = startsWith(B.Status, st(k));
        if any(m), scatter(ax, B.Freq_Hz(m), B.Delta_dB(m), 45, col(k, :), 'filled', 'DisplayName', st(k)); end
    end
    yline(ax, 0, 'k-'); xlabel(ax, 'Frequency (Hz)'); ylabel(ax, 'level now - baseline (dB, median of sweeps)');
    legend(ax, 'Location', 'best'); title(ax, "Baseline comparison: " + R.baseline.source, 'Interpreter', 'none');
    files(end + 1) = saveFig(fig, fullfile(d, "fig_baseline.png"), opt);
end
% --- watch
if isfield(R.watch, 'passes') && height(R.watch.passes) > 0
    fig = newFig(opt, [40 40 1100 420]); ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
    plot(ax, R.watch.passes.Pass, R.watch.passes.Emitters, '-o', 'DisplayName', 'emitters');
    A = R.watch.alerts;
    if height(A) > 0
        [g, ty] = findgroups(A.Type);
        for k = 1:numel(ty)
            m = g == k; stem(ax, A.Pass(m) + 0.05*k, ones(nnz(m), 1)*k, 'filled', 'DisplayName', "alert: " + ty(k));
        end
    end
    xlabel(ax, 'pass'); legend(ax, 'Location', 'best'); title(ax, 'Watch: emitters per pass and alerts');
    files(end + 1) = saveFig(fig, fullfile(d, "fig_watch.png"), opt);
end
end

function bandStrip(ax, R, yl)
BP = loadBandPlan(R.meta.opt);
kinds = unique(BP.Kind);
cmap = lines(max(numel(kinds), 1));
x = xlim(ax);
lo = max(BP.FLow_Hz, x(1)); hi = min(BP.FHigh_Hz, x(2)); in = find(hi > lo);
if isempty(in), return, end
[~, ki] = ismember(BP.Kind(in), kinds);
patch(ax, [lo(in) hi(in) hi(in) lo(in)].', repmat([yl(2) - 8; yl(2) - 8; yl(2) - 2; yl(2) - 2], 1, numel(in)), ...
    reshape(cmap(ki, :), [1 numel(in) 3]), 'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
for i = in(:).'
    lo = max(BP.FLow_Hz(i), x(1)); hi = min(BP.FHigh_Hz(i), x(2));
    if hi/lo > 1.1
        text(ax, sqrt(lo*hi), yl(2) - 5, BP.Name(i), 'FontSize', 6, 'HorizontalAlignment', 'center', 'Rotation', 90, ...
            'Interpreter', 'none', 'Clipping', 'on');
    end
end
end

function s = shortLabel(row)
s = sprintf('%s %s', fmtHz(row.Freq_Hz), row.SubType);
s = char(s); if numel(s) > 38, s = [s(1:36) '..']; end
end

function cc = classColors()
cc.names = ["Licensed service" "Satellite" "EMI" "Unknown"];
cc.rgb = [0 0.45 0.74; 0.49 0.18 0.56; 0.85 0.33 0.1; 0.5 0.5 0.5];
end

function fig = newFig(opt, pos)
% Visible 'off': one hidden figure is reused (creating a canvas costs ~1-2 s each in R2026a)
persistent shared
if strcmpi(opt.Visible, 'off')
    if isempty(shared) || ~isgraphics(shared)
        shared = figure('Visible', 'off', 'Color', 'w', 'Position', pos, 'ToolBar', 'none', 'MenuBar', 'none', ...
            'HandleVisibility', 'off');
    else
        clf(shared); set(shared, 'Position', pos, 'Color', 'w');
    end
    fig = shared;
else
    fig = figure('Visible', opt.Visible, 'Color', 'w', 'Position', pos, 'ToolBar', 'none', 'MenuBar', 'none');
end
end

function [x, y] = decim(x, y, n, how)
% Fewer points for plotting (faster rendering): block max (peak traces) or block mean.
N = numel(x);
if N <= n, return, end
k = ceil(N/n); m = floor(N/k);
X = reshape(x(1:m*k), k, m); Y = reshape(y(1:m*k), k, m);
x = mean(X, 1);
if how == "max", y = max(Y, [], 1); else, y = 10*log10(mean(10.^(Y/10), 1)); end
end

function file = saveFig(fig, file, ~)
try                                    % getframe: ~10x faster than print/exportgraphics in R2026a
    drawnow; fr = getframe(fig);
    if isempty(fr.cdata) || size(fr.cdata, 1) < 50, error('ff_auto_survey:frame', 'empty frame'); end
    imwrite(fr.cdata, char(file));
catch
    try
        print(fig, char(file), '-dpng', '-r110');
    catch
        exportgraphics(fig, file, 'Resolution', 110);
    end
end
file = string(file);                   % the hidden shared figure is closed at the end of makeFigures
end

function file = makeReport(R, opt, files)
import mlreportgen.report.*
import mlreportgen.dom.*
fmt = char(opt.Report);
fprintf('Writing the %s report ...\n', upper(fmt));
rpt = Report(char(fullfile(opt.OutDir, "survey_report")), fmt);
when = datetime('now', 'TimeZone', 'Asia/Bangkok', 'Format', 'yyyy-MM-dd HH:mm ''(Asia/Bangkok)''');
sub = "FieldFox autonomous survey (" + R.meta.preset + ")";
if R.meta.synthetic, sub = "SYNTHETIC DEMO DATA - not a measurement"; end
add(rpt, TitlePage('Title', 'Spectrum survey & EMI inventory', 'Subtitle', sub, ...
    'Author', string(getenv('USERNAME')), 'PubDate', char(when)));
add(rpt, TableOfContents);
png = @(n) char(fullfile(opt.OutDir, n));
ch = Chapter('Title', 'Summary');
p = Paragraph(['Classes are HYPOTHESES from spectrum shape, band plan, timing and catalogs (no I/Q, nothing is ' ...
    'demodulated). The Evidence column says why; Confidence is a heuristic score.']);
p.Bold = true; p.Color = 'red'; add(ch, p);
add(ch, UnorderedList(cellstr(R.conclusions)));
if isfile(png("fig_overview.png")), add(ch, FormalImage('Image', png("fig_overview.png"), 'Caption', 'Survey overview with the band plan', 'Width', '16cm')); end
add(rpt, ch);
ch = Chapter('Emitter inventory');
V = R.inventory(:, {'ID', 'Freq_Hz', 'OBW_Hz', 'Peak_dBm', 'Excess_dB', 'Occupancy_pct', 'Class', 'SubType', 'Confidence', 'Band'});
add(ch, rptTable(V, 'Inventory (all emitters)', 80));
E = R.inventory(R.inventory.Class ~= "Licensed service", {'ID', 'Freq_Hz', 'Class', 'SubType', 'Evidence'});
if height(E) > 0, add(ch, rptTable(E, 'Evidence for EMI, satellite and unknown emitters', 40)); end
if isfile(png("fig_segments.png")), add(ch, FormalImage('Image', png("fig_segments.png"), 'Caption', 'Segments: floor, threshold and emitters (colour = class)', 'Width', '16cm')); end
if isfile(png("fig_zoom.png")), add(ch, FormalImage('Image', png("fig_zoom.png"), 'Caption', 'Zoomed measurements', 'Width', '16cm')); end
add(rpt, ch);
ch = Chapter('EMI, timing and occupancy');
if height(R.combs) > 0, add(ch, rptTable(R.combs, 'Harmonic combs', 10)); end
if height(R.rises) > 0, add(ch, rptTable(R.rises, 'Broadband noise rises', 10)); end
if height(R.zeroSpan) > 0, add(ch, rptTable(R.zeroSpan, 'Zero-span captures', 20)); end
if isfile(png("fig_zerospan.png")), add(ch, FormalImage('Image', png("fig_zerospan.png"), 'Caption', 'Zero span (first sweep)', 'Width', '16cm')); end
add(ch, rptTable(R.bandOccupancy(:, {'Band', 'Occupancy_pct', 'MaxLevel_dBm', 'Floor_dBm', 'Emitters', 'Classes'}), 'Band occupancy', 30));
if isfile(png("fig_waterfall.png")), add(ch, FormalImage('Image', png("fig_waterfall.png"), 'Caption', 'Waterfall per segment', 'Width', '16cm')); end
add(rpt, ch);
if isfield(R.baseline, 'table') || isfield(R.watch, 'alerts')
    ch = Chapter('Changes');
    if isfield(R.baseline, 'table')
        add(ch, Paragraph(sprintf('Baseline %s: new %d, gone %d, changed %d, same %d (Wilcoxon rank-sum on per-sweep levels, BH q < %g and |delta| >= 3 dB).', ...
            R.baseline.source, R.baseline.nNew, R.baseline.nGone, R.baseline.nChanged, R.baseline.nSame, opt.Alpha)));
        add(ch, rptTable(R.baseline.table(R.baseline.table.Status ~= "Same", :), 'Differences', 40));
        if isfile(png("fig_baseline.png")), add(ch, FormalImage('Image', png("fig_baseline.png"), 'Caption', 'Baseline comparison', 'Width', '16cm')); end
    end
    if isfield(R.watch, 'alerts')
        add(ch, rptTable(R.watch.alerts, 'Watch alerts', 60));
        if isfile(png("fig_watch.png")), add(ch, FormalImage('Image', png("fig_watch.png"), 'Caption', 'Watch', 'Width', '15cm')); end
    end
    add(rpt, ch);
end
ch = Chapter('Setup, method and caveats');
add(ch, Paragraph(sprintf('Instrument: %s. Source: %s. Location %.4f, %.4f. Survey time (UTC): %s.', R.meta.idn, R.meta.source, ...
    R.meta.location(1), R.meta.location(2), string(R.meta.time))));
if strlength(string(opt.Notes)) > 0, add(ch, Paragraph("Notes: " + string(opt.Notes))); end
if height(R.plan) > 0, add(ch, rptTable(R.plan, 'Scan plan', 20)); end
add(ch, UnorderedList({ ...
    'Floor: block 20th percentile, running minimum over 5 blocks, corrected to the median of noise-only points; detection above max(Threshold, K x noise spread).', ...
    'OBW = 99 % of the noise-subtracted power; shoulder = -3 dB width / OBW (1 = flat top). Occupancy = % of sweeps in which the emitter exceeds the threshold.', ...
    'Comb: periodogram of the dB excess (cepstrum), lines near each multiple, least-squares spacing with a 95 % interval.', ...
    'Zero span: 20 ms AVER (TDD period and slots) and 100 ms SAMP (bursts; contrast must exceed what noise alone gives, Monte Carlo).', ...
    'Satellites: sat_channels.csv (Doppler + OBW/2 window, width ratio < 5), elevation and Doppler now from SGP4 (propagateOrbit, eci2aer) with cached TLEs.', ...
    'Levels are in the RBW of each segment; no antenna factor: these are not field strengths.', ...
    'No I/Q: technologies are recognised by spectrum shape only (ff_signal_id); verify important findings with a narrower span/RBW or a demodulating receiver.'}));
if ~isempty(R.notes), add(ch, UnorderedList(cellstr(R.notes(1:min(end, 40))))); end
add(ch, UnorderedList(cellstr([files.mat; files.inventory; files.occupancy; files.figures(:)])));
add(rpt, ch);
close(rpt);
file = string(rpt.OutputPath);
fprintf('Report: %s\n', file);
end

function bt = rptTable(T, ttl, maxRows)
import mlreportgen.report.*
import mlreportgen.dom.*
T = T(1:min(end, maxRows), :);
for v = string(T.Properties.VariableNames)
    x = T.(v);
    if isdatetime(x), T.(v) = string(x, 'yyyy-MM-dd HH:mm:ss');
    elseif endsWith(v, "_Hz") && isnumeric(x), T.(v) = string(arrayfun(@fmtHz, x, 'UniformOutput', false));
    elseif isnumeric(x), T.(v) = compose("%.4g", x);
    elseif islogical(x), T.(v) = string(x);
    end
    if isstring(T.(v)), T.(v) = cellstr(T.(v)); end
end
mt = MATLABTable(T);
mt.Style = {FontSize('7pt')};
bt = BaseTable(mt);
bt.Title = ttl;
end

function txt = conclusions(R, opt)
txt = strings(0, 1);
inv = R.inventory;
if R.meta.synthetic, txt(end + 1) = "SYNTHETIC DEMO DATA - not a measurement."; end
cnt = @(c) nnz(inv.Class == c);
txt(end + 1) = sprintf("%d emitters: %d licensed service, %d satellite, %d EMI, %d unknown.", height(inv), ...
    cnt("Licensed service"), cnt("Satellite"), cnt("EMI"), cnt("Unknown"));
for i = find(inv.Kind == "comb").'
    txt(end + 1) = sprintf("EMI: harmonic comb, spacing %s, %d lines from %s to %s (switching converter / inverter signature).", ...
        fmtHz(inv.Spacing_Hz(i)), inv.Lines(i), fmtHz(inv.Lo_Hz(i)), fmtHz(inv.Hi_Hz(i))); %#ok<AGROW>
end
for i = find(startsWith(R.combs.Segment, "zoom")).'
    c = R.combs(i, :);
    txt(end + 1) = sprintf("EMI: line comb, spacing %s [%s, %s], %d lines in the %s (%s-%s)%s.", fmtHz(c.Spacing_Hz), ...
        fmtHz(c.CIlow_Hz), fmtHz(c.CIhigh_Hz), c.Lines, c.Segment, fmtHz(c.First_Hz), fmtHz(c.Last_Hz), ...
        ternary(c.Spacing_Hz >= 2e3 && c.Spacing_Hz <= 200e3, " - switching converter / inverter signature", "")); %#ok<AGROW>
end
for i = find(inv.Kind == "rise").'
    txt(end + 1) = sprintf("EMI candidate: broadband noise rise +%.1f dB over %s-%s%s.", inv.Excess_dB(i), fmtHz(inv.Lo_Hz(i)), ...
        fmtHz(inv.Hi_Hz(i)), ternary(inv.Hi_Hz(i) <= 30e6, " (below 30 MHz: may be ambient noise - check with a baseline)", "")); %#ok<AGROW>
end
for i = find(inv.Kind == "cluster").'
    txt(end + 1) = sprintf("EMI: %s, strongest %.1f dBm%s.", inv.SubType(i), inv.Peak_dBm(i), ...
        ternary(strlength(inv.ZeroSpan(i)) > 0, "; zero span: " + inv.ZeroSpan(i), "")); %#ok<AGROW>
end
k = contains(inv.ZeroSpan, "mains-synchronous");
for i = find(k).'
    txt(end + 1) = sprintf("Mains-synchronous bursts at %s: %s.", fmtHz(inv.Freq_Hz(i)), inv.ZeroSpan(i)); %#ok<AGROW>
end
for i = find(inv.Class == "EMI" & contains(inv.SubType, "GNSS")).'
    txt(end + 1) = sprintf("Interference in a GNSS band at %s (%.1f dBm): GNSS signals are below the FieldFox noise, so this line is not GNSS.", ...
        fmtHz(inv.Freq_Hz(i)), inv.Peak_dBm(i)); %#ok<AGROW>
end
for i = find(contains(inv.ZeroSpan, "TDD-like") & startsWith(inv.Band, ["2300" "2600"])).'
    txt(end + 1) = sprintf("TDD carrier at %s (%s wide, %s): %s.", fmtHz(inv.Freq_Hz(i)), fmtHz(inv.OBW_Hz(i)), inv.Band(i), inv.ZeroSpan(i)); %#ok<AGROW>
end
for i = find(inv.Class == "Satellite").'
    txt(end + 1) = sprintf("Satellite candidate %s at %s (%.1f dBm, confidence %.2f).", inv.SubType(i), fmtHz(inv.Freq_Hz(i)), inv.Peak_dBm(i), inv.Confidence(i)); %#ok<AGROW>
end
u = inv(inv.Class == "Unknown", :);
if height(u) > 0
    txt(end + 1) = sprintf("Unknown: %s - narrow the span/RBW on them (ff_signal_id) or compare with a baseline.", ...
        strjoin(sfun(@fmtHz, u.Freq_Hz(1:min(end, 8))), ", "));
end
B = R.bandOccupancy;
if height(B) > 0
    txt(end + 1) = "Busiest bands: " + strjoin(compose("%s %.0f%%", B.Band(1:min(end, 4)), B.Occupancy_pct(1:min(end, 4))), ", ") + ".";
end
if isfield(R.baseline, 'table')
    txt(end + 1) = sprintf("Against the baseline: %d new, %d gone, %d changed, %d unchanged emitters.", R.baseline.nNew, R.baseline.nGone, R.baseline.nChanged, R.baseline.nSame);
end
if isfield(R.watch, 'alerts')
    txt(end + 1) = sprintf("Watch: %d pass(es), %d alert(s) (%s).", height(R.watch.passes), height(R.watch.alerts), ...
        strjoin(unique(R.watch.alerts.Type), ", "));
end
if any(contains(R.notes, "preamp switched off")), txt(end + 1) = "Preamp was switched off in some segments (strong signals): their floor is ~15 dB higher."; end
if opt.SatID && cnt("Satellite") == 0
    txt(end + 1) = "No satellite carrier identified: most need a gain antenna + LNA and a pass overhead (see ff_cubesat_link).";
end
end

function printSummary(R)
inv = R.inventory;
fprintf('\n%-4s %-13s %-10s %8s %6s %5s  %-16s %4s  %s\n', 'ID', 'Freq', 'OBW', 'Peak', 'Exc', 'Occ%', 'Class', 'Conf', 'Sub-type');
for i = 1:height(inv)
    fprintf('%-4d %-13s %-10s %8.1f %6.1f %5.0f  %-16s %4.2f  %s\n', inv.ID(i), fmtHz(inv.Freq_Hz(i)), fmtHz(inv.OBW_Hz(i)), ...
        inv.Peak_dBm(i), inv.Excess_dB(i), inv.Occupancy_pct(i), inv.Class(i), inv.Confidence(i), ...
        extractBefore(inv.SubType(i) + " ", min(70, strlength(inv.SubType(i)) + 1)));
end
fprintf('\n'); fprintf('  - %s\n', R.conclusions); fprintf('\n');
end

%% ===== small helpers ==========================================================

function tf = isAvgDet(det)
tf = ismember(upper(string(det)), ["AVER" "RMS" "AVERAGE"]);
end

function r = snapRBW(x)
s = [10 30 100 300 1e3 3e3 10e3 30e3 100e3 300e3 1e6 3e6 5e6];
[~, i] = min(abs(log(s) - log(max(x, 10)))); r = s(i);
end

function [s, e] = runs(mask)
d = diff([false, logical(mask(:).'), false]); s = find(d == 1); e = find(d == -1) - 1;
end

function s = joinS(a, b)
s = strjoin([string(a), string(b)], "; ");
s = regexprep(s, '^(; )+|(; )+$', '');
end

function s = sfun(fn, A)
s = string(arrayfun(fn, A, 'UniformOutput', false));
end

function h = fileAgeH(p)
d = dir(p); h = (now - d.datenum)*24; %#ok<TNOW1>
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function s = fmtHz(x)
if ~isfinite(x), s = "-"; return, end
a = abs(x);
if a >= 1e9, s = sprintf("%.4g GHz", x/1e9);
elseif a >= 1e6, s = sprintf("%.4g MHz", x/1e6);
elseif a >= 1e3, s = sprintf("%.4g kHz", x/1e3);
else, s = sprintf("%.4g Hz", x); end
s = string(s);
end
