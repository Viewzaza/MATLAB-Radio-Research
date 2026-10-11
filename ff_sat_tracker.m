function R = ff_sat_tracker(varargin)
%FF_SAT_TRACKER  Doppler-tracking satellite capture with the FieldFox + measured-vs-predicted orbit check.
%
%   R = ff_sat_tracker("10.90.36.200", 'Satellite', "OrigamiSat-2")  next pass: wait for AOS, track, analyse
%   R = ff_sat_tracker(ff, 'Satellite', 67683, 'Mode', "wide")        existing FieldFox object, fixed wide span
%   R = ff_sat_tracker("demo")                                        simulated OrigamiSat-2 pass (offline)
%   R = ff_sat_tracker('Simulate', true, 'Satellite', "METEOR-M2 4")  simulated pass of another satellite
%   R = ff_sat_tracker('File', "ff_sat_tracker_20261010_155300.mat")  re-analyse a recording
%   P = ff_sat_tracker('Satellite', "KNACKSAT-2", 'PlanOnly', true)   pass list + capture plan only
%   R = ff_sat_tracker("10.90.36.200", 'Rehearse', true)              dry run of the retune loop now, on
%                                                                     the geometry of the next pass
% Example (headless):
%   R = ff_sat_tracker("demo", 'Visible', 'off', 'OutDir', "tracker_demo");
%   R = ff_sat_tracker("10.90.36.200", 'Satellite', "OrigamiSat-2", 'RxAntenna', "tracking", 'RxGain', 12, ...
%        'LNAGain', 20, 'LNANoiseFigure', 0.8, 'Prompt', false, 'Visible', 'off', 'OutDir', "pass_origami2");
%
% What it does:
%  1. Target: 'Satellite' (name, catalog Id or NORAD) -> downlink f0, occupied bandwidth, modulation class
%     and EIRP from the ff_satsig catalog (else the SatNOGS subset it ships, else 'Frequency'/'SignalBW').
%  2. Orbit: TLE from CelesTrak (gp.php?CATNR=..&FORMAT=TLE, cached 2 h in support\tracker\tle, no retry
%     for 2 h after an error); offline: any cached copy (this tool's, ff_satsig's, ff_cubesat_link's), then
%     the shipped support\tracker\tracker_sample.tle. satelliteScenario + access -> pass list.
%  3. Pass geometry at 1 s: satellite (SGP4) + states (ECEF) -> range rate -> Doppler factor and rate,
%     cross-checked against dopplershift (agreement printed). Link: transmitter with the satellite
%     antenna pattern (isotropic / Antenna Toolbox dipole / Gaussian beam, pointAt) -> txEIRP, sigstrength
%     (received isotropic power), plus the ground antenna (tracking gain or a fixed vertical whip pattern),
%     polarisation/misc losses, LNA, ITU-R P.372 man-made noise and the FieldFox DANL -> predicted C/N0.
%  4. Capture plan: "retune" = narrow span; the centre is moved to (f0 + 'FreqOffset')(1 + d(t)) every
%     'Step' s (SENS:FREQ:CENT, then INIT:IMM;*OPC? and a binary trace). "wide" = fixed centre and a span
%     covering +-(max Doppler + tolerance). Span = signal OBW + 2 x (oscillator tolerance + Doppler rate x
%     (TLE timing tolerance + Step + cycle)); RBW = 1-3-10 step >= 2 bins and >= OBW/20; preamp on.
%  5. Live: waits for AOS with a countdown (or starts at once inside a pass), takes control 'Lead' s before
%     AOS (takeControl, guaranteed release via onCleanup), records every sweep with PC timestamps
%     (before INIT, after *OPC?, after the trace read) and saves every 'SaveEvery' s, so an error or
%     Ctrl+C keeps the data recorded so far. Nothing is transmitted (SA mode only).
%  6. Analysis (live, simulated and file are identical): per sweep (wide continuous signals: per block
%     of sweeps averaged over 'Integrate' s) a sigma-clipped noise floor, the strongest peak inside the
%     predicted window (+-search width): narrow signals 3-bin peak with SNR gate 'SNRMin' and a parabolic
%     peak in dB; wide signals an OBW-wide box (matched) filter with gate 'ZMin' sigma and a half-maximum
%     centroid. Time of the carrier bin inside the sweep, channel power and C/N0. Robust track
%     cleaning, then lsqcurvefit of
%         f(t) = (f0 + df + a (t - tm)) (1 + d(t + dt; x))
%     for three models: "dt+df" (time offset + transmitter offset), "dt+df+drift" (+ oscillator drift a),
%     "dt+df+xtrack" (+ cross-track TLE position error x, km). dt is the PC-clock error PLUS the TLE
%     along-track error (they are indistinguishable; along-track ~ v dt). Parameter 1-sigma from the
%     Jacobian, residual RMS, BIC model choice, a fixed-carrier (no Doppler) comparison and a verdict.
%     Doppler-corrected waterfall, measured vs predicted Doppler, C/N0 vs elevation against the link
%     budget (median offset = EIRP/antenna error if your receive chain parameters are right), measured vs
%     predicted noise floor (excess = local EMI or a wrong LNA/antenna-noise assumption), timing.
%  7. 'Simulate': a full pass recording is synthesised (template spectra of the signal class from
%     support\satsig\satsig_templates.mat, else a CW/Gaussian model; real SGP4 geometry; known time,
%     frequency and EIRP errors; slow fading; Morse keying for CW; FieldFox RBW shape, detector noise and
%     sweep timing) and analysed by the same code, so the recovered errors can be checked against truth.
%
% Output R: source, target, tle, site, passes (table), pass, plan, rec (recording), track (table, one
% row per sweep, or per averaged block: column nSweeps), fits (table), best, link, timing, verdict
% (CONFIRMED / PROBABLE / NOT CONFIRMED / INCONCLUSIVE), sim (truth + errors), report (text),
% notes, files, figs. Recording MAT (shared suite format): f (1xN Hz), P (sweeps x N dBm), t (datetime
% UTC, mid-sweep), settings (rbw, vbw, span, points, detector, preamp, ...), idn, plus fCenter (actual
% centre per sweep), tStart/tEnd, timing, tracker (target, TLE lines, pass, plan) and analysis.
% In "retune" recordings f is the TRACKING frame: f0 + offset from the sweep centre; the absolute
% frequency of bin j of sweep k is f(j) - settings.frameCenter + fCenter(k).
%
% Name-Value options (defaults):
%  Target / orbit
%   'Satellite'      "OrigamiSat-2"  name, catalog Id ("origami2_uhf_cw") or NORAD number
%   'Downlink'       ""      catalog Id when the satellite has several downlinks (default: the most
%                            active one inside the FieldFox range)
%   'Frequency'      []      Hz, nominal downlink (overrides the catalog)
%   'SignalBW'       []      Hz, occupied bandwidth (default: catalog OBW, else 1 kHz)
%   'TLE'            ""      TLE file instead of CelesTrak (entry with the NORAD, else the first)
%   'Offline'        false   never go online (cached TLE of any age, else the shipped sample)
%  Site / time / pass
%   'Site'           [13.7563 100.5018 10]  [lat lon alt_m]  SET YOUR LAB LOCATION
%   'SiteName'       "Bangkok";  'TimeZone' "Asia/Bangkok" (printed times)
%   'Pass'           "next"  "next" = in progress or next, "now" = in progress only, k = k-th pass
%   'Hours'          24      pass search window;  'MinElevation' 5 deg (mask = AOS/LOS)
%   'MinPeakElevation' 15    deg; "next" skips lower passes
%   'Start'          []      reference time for Simulate/PlanOnly (datetime; default now)
%   'Wait'           true    block until AOS with a countdown; false = return the plan if not in a pass
%   'MaxWait'        3       h; a later pass returns the plan instead of waiting
%   'Lead'           30      s recorded before AOS and after LOS (TLE timing errors)
%   'Duration'       []      s, stop after this (default: until LOS + Lead)
%   'ClockShift'     0       s added to the PC clock for all predictions (known clock error)
%   'Rehearse'       false   dry run NOW on the geometry of the chosen pass (ClockShift set so the run
%                            is centred on TCA; 'Duration' default 60 s): tests retune timing; no real
%                            satellite signal is expected (use a signal generator if you want one)
%  Capture
%   'Port'           5025    SCPI socket port
%   'Mode'           "retune"  or "wide"
%   'Step'           2       s between centre updates in "retune" (0 = every sweep)
%   'Span'           []      Hz (auto);  'RBW' [] Hz (auto, 1-3-10);  'Points' 401
%   'Preamp'         true;   'Detector' "AVER" (RMS; "POS" catches short bursts but biases the floor)
%   'RefLevel'       -40 dBm; 'Atten' [] (auto)
%   'FreqOffset'     0       Hz, known transmitter offset (e.g. df from an earlier pass) applied to tuning
%   'FreqTolPPM'     5;  'FreqTolHz' 500   satellite oscillator tolerance for the span/search window
%   'TimeTol'        []      s TLE timing tolerance ([] = 5 s + 3 s per day of TLE age)
%   'SaveTo'         ""      recording file (default ff_sat_tracker_<time>.mat in OutDir or pwd)
%   'SaveEvery'      30      s between intermediate saves;  'LivePlot' true (waterfall while recording)
%  Analysis
%   'SNRMin'         7       dB, narrow signals (OBW <= 3 RBW): 3-bin smoothed peak over the mean floor
%   'ZMin'           6       wide signals: OBW-wide matched-filter (box) peak, in noise standard deviations
%   'Integrate'      []      s, sweeps averaged (same centre, linear power) before detection. [] = auto:
%                            wide continuous signals 5 s (limited to 5 % OBW Doppler smear), else 0
%                            (single sweeps; bursty/keyed signals are never averaged automatically)
%   'MaxTimeOffset'  60      s (fit bound for dt)
%  Link budget (prediction of C/N0 and of the noise floor)
%   'EIRP'           []      dBm (catalog value, else 30 dBm ASSUMED)
%   'TxAntenna'      "auto"  "isotropic" (tumbling/unknown attitude; default below 2 GHz) | "dipole"
%                            (Antenna Toolbox, nadir-stabilised, axis across the orbit plane) | "beam"
%                            (Gaussian, 'TxGain' dBi; default above 2 GHz)
%   'TxGain'         18      dBi (beam);  'TxPointing' "station" (pointAt the ground station) | "nadir"
%   'RxAntenna'      "tracking"  gain 'RxGain' toward the satellite (rotator) | "whip" (fixed vertical
%                            dipole pattern from the Antenna Toolbox; null at zenith)
%   'RxGain'         0       dBi;  'PolLoss' 3 dB;  'MiscLoss' 1 dB (feeder before the LNA, pointing)
%   'LNAGain'        0 dB;  'LNANoiseFigure' 1 dB;  'CableLoss' 0 dB (LNA to FieldFox)
%   'Environment'    "residential"  P.372 man-made noise ("city","residential","rural","quiet")
%   'AntennaTemp'    []      K (overrides P.372);  'DANL' [] dBm/Hz ([] = -153 preamp on, -137 off)
%  Simulation ('Simulate', true or first argument "demo")
%   'Seed' 1; 'SimTimeError' 3 s; 'SimFreqOffset' [] Hz (+2.5 ppm); 'SimCrossTrack' 0 km;
%   'SimEIRPError' -2 dB; 'SimFading' 1.5 dB (slow log-normal); 'SimDuty' [] ([on period] s; CW
%   default = Morse-like keying); 'SimInterferer' [] ([offset_Hz level_dBm] fixed carrier);
%   'SimSweepTime' [] s (model); 'SimNoiseExcess' 0 dB (local EMI on top of the predicted floor)
%  Output
%   'File' ""; 'PlanOnly' false; 'Plot' true; 'Visible' 'on'|'off' ('off': closed once saved); 'OutDir' "" (no files unless live);
%   'Prompt' true (confirm before taking control; dialog only in an interactive desktop session);
%   'Quiet' false; 'PDF' false (Report Generator PDF in OutDir)
%
% Honest limits:
%  - No I/Q on this FieldFox: the carrier frequency comes from the swept trace (resolution ~ bin width /
%    10 for a clean CW with parabolic interpolation, ~ RBW/3 for wide signals). Nothing is demodulated.
%  - The fitted df includes the FieldFox time-base error (use an external 10 MHz GPSDO reference for
%    sub-ppm claims); dt includes the PC clock error (keep Windows time NTP-synchronised).
%  - Morse/beacon signals are only seen when the sweep crosses the carrier while it is keyed; bursty
%    packets (9k6 FSK every 60 s) give few points: use "wide" + more sweeps, or ff_satsig MAXH windows.
%  - Short or low passes make dt and df degenerate (a straight Doppler line): read the 1-sigma values.
%  - The link budget uses ASSUMED EIRPs where none is published (catalog column Source); polarisation
%    and tumbling are a fixed loss. The C/N0 offset is a calibration check, not a measurement of EIRP
%    unless the receive chain (antenna gain, LNA, cables) is known.
%  - The AVER (RMS) detector is assumed for the noise floor; other detectors are flagged.
%  - Simulated recordings are idealised (template spectra, no adjacent signals unless 'SimInterferer').
%
% See also ff_satsig, ff_cubesat_link, ff_ccsds_linksim, FieldFox.

if nargin == 0
    help('ff_sat_tracker'); if nargout > 0, R = []; end, return
end
args = varargin; src = [];
if isa(args{1}, 'FieldFox') || mod(numel(args), 2) == 1
    src = args{1}; args = args(2:end);
end
opt = parseOpts(args);
if (ischar(src) || isstring(src)) && any(strcmpi(string(src), ["demo", "simulate", "sim"]))
    opt.Simulate = true; src = [];
end
notes = strings(0, 1);
R = struct('source', "", 'target', [], 'tle', [], 'site', struct('lla', opt.Site, 'name', opt.SiteName), ...
    'passes', table(), 'pass', [], 'plan', [], 'rec', [], 'track', table(), 'fits', table(), 'best', [], ...
    'link', [], 'timing', [], 'verdict', "", 'sim', [], 'report', strings(0, 1), 'notes', strings(0, 1), ...
    'files', strings(0, 1), 'figs', gobjects(0), 'status', "");

if strlength(opt.File) > 0
    %% ---- re-analysis of a recording ------------------------------------------------------
    [rec, tgt, tle, notes] = loadRecording(opt, notes);
    R.source = "file " + opt.File;
    tA = min(rec.tStart); tB = max(rec.tEnd);
    eph = ephemeris(tle, tA - seconds(300), tB + seconds(300), tgt, opt, rec.settings.preamp);
    pass = passFromEph(eph, tA, tB, opt);
    plan = planFromRec(rec, eph, tgt, tle, pass, opt);
else
    %% ---- target, orbit, pass, plan ---------------------------------------------------------
    [tgt, notes] = resolveTarget(opt, notes);
    [tle, notes] = getTLE(tgt, opt, notes);
    if opt.Simulate || opt.PlanOnly
        if isempty(opt.Start), tNow = datetime('now', 'TimeZone', 'UTC'); else, tNow = toUTC(opt.Start, opt.TimeZone); end
    else
        tNow = nowUTC(opt);
    end
    passes = findPasses(tle, tNow - minutes(20), opt);
    R.passes = passes;
    [pass, notes] = choosePass(passes, tNow, opt, notes);
    if isempty(pass)
        R.target = tgt; R.tle = tle; R.notes = notes(:); R.status = "no pass";
        say(opt, '%s\n', notes); return
    end
    if opt.Rehearse && ~opt.Simulate
        if isempty(opt.Duration), opt.Duration = 60; end
        opt.ClockShift = opt.ClockShift + seconds(pass.TCA - nowUTC(opt)) - opt.Duration/2;
        tNow = nowUTC(opt);
        notes(end+1) = sprintf(['REHEARSAL: the PC clock is shifted by %+.0f s so this run replays the pass of %s ' ...
            'around TCA. No real satellite signal is expected.'], opt.ClockShift, localStr(pass.TCA, opt));
    end
    if opt.Simulate, tNow = pass.AOS - seconds(opt.Lead); end
    plan = makePlan(tgt, tle, pass, tNow, opt);
    eph = plan.eph;
    notes = [notes(:); plan.notes(:)];
    R.target = tgt; R.tle = tle; R.pass = pass; R.plan = rmfield(plan, 'eph');
    if opt.PlanOnly
        R.status = "plan"; R.notes = notes(:);
        if ~opt.Quiet, printPlan(R, opt); end
        return
    end
    if opt.Simulate
        R.source = "simulated";
        [rec, R.sim] = simulateRec(plan, tgt, tle, opt);
    else
        if isempty(src), error('ff_sat_tracker:usage', 'Give an IP address or a FieldFox object, or use ''Simulate'' / ''File''.'); end
        R.source = "live " + string(srcName(src));
        wait = seconds(plan.tStart - nowUTC(opt));
        if wait > 0 && (~opt.Wait || wait > 3600*opt.MaxWait)
            R.status = "waiting not allowed"; R.notes = [notes(:); sprintf(['The pass starts in %.1f min (AOS %s): not ' ...
                'waiting (''Wait'' false or more than ''MaxWait'' %.1f h). Plan returned.'], wait/60, localStr(pass.AOS, opt), opt.MaxWait)];
            if ~opt.Quiet, printPlan(R, opt); fprintf('%s\n', R.notes(end)); end
            return
        end
        if ~opt.Quiet, printPlan(R, opt); end
        [rec, notes, stopped] = acquireLive(src, plan, tgt, tle, opt, notes);
        if stopped == "declined"
            R.status = "declined"; R.notes = notes(:); return
        end
    end
end

%% ---- analysis, figures, report ----------------------------------------------------------
if isempty(rec.P) || size(rec.P, 1) < 1
    R.notes = [notes(:); "No sweeps were recorded."]; R.status = "no data"; R.rec = rec;
    say(opt, 'No sweeps were recorded.\n'); return
end
A = analyzeRec(rec, eph, tgt, plan, opt);
R.target = tgt; R.tle = tle; R.pass = pass; R.plan = rmfield(plan, 'eph');
R.rec = rec; R.track = A.track; R.fits = A.fits; R.best = A.best; R.link = A.link;
R.timing = A.timing; R.verdict = A.verdict; R.status = "analysed";
if ~isempty(R.sim), R.sim = simCompare(R.sim, A); end
R.notes = [notes(:); A.notes(:)];
R.report = buildReport(R, opt);
if ~opt.Quiet, fprintf('%s\n', R.report); end
if opt.Plot
    R.figs = makeFigures(R, A, eph, opt);
end
R.files = saveOutputs(R, A, opt);
if strcmpi(opt.Visible, 'off') && (any(endsWith(R.files, ".png")) || nargout == 0)
    delete(R.figs(isgraphics(R.figs)));                         % hidden figures: none left once saved
end
end

%% ======================================================================================
%% options
%% ======================================================================================
function opt = parseOpts(args)
opt = struct('Satellite', "OrigamiSat-2", 'Downlink', "", 'Frequency', [], 'SignalBW', [], 'TLE', "", ...
    'Offline', false, 'Site', [13.7563 100.5018 10], 'SiteName', "Bangkok", 'TimeZone', "Asia/Bangkok", ...
    'Pass', "next", 'Hours', 24, 'MinElevation', 5, 'MinPeakElevation', 15, 'Start', [], 'Wait', true, ...
    'MaxWait', 3, 'Lead', 30, 'Duration', [], 'ClockShift', 0, 'Rehearse', false, 'Port', 5025, ...
    'Mode', "retune", 'Step', 2, 'Span', [], 'RBW', [], 'Points', 401, 'Preamp', true, 'Detector', "AVER", ...
    'RefLevel', -40, 'Atten', [], 'FreqOffset', 0, 'FreqTolPPM', 5, 'FreqTolHz', 500, 'TimeTol', [], ...
    'SaveTo', "", 'SaveEvery', 30, 'LivePlot', true, 'SNRMin', 7, 'ZMin', 6, 'Integrate', [], ...
    'MaxTimeOffset', 60, 'EIRP', [], ...
    'TxAntenna', "auto", 'TxGain', 18, 'TxPointing', "station", 'RxAntenna', "tracking", 'RxGain', 0, ...
    'PolLoss', 3, 'MiscLoss', 1, 'LNAGain', 0, 'LNANoiseFigure', 1, 'CableLoss', 0, ...
    'Environment', "residential", 'AntennaTemp', [], 'DANL', [], 'Simulate', false, 'Seed', 1, ...
    'SimTimeError', 3, 'SimFreqOffset', [], 'SimCrossTrack', 0, 'SimEIRPError', -2, 'SimFading', 1.5, ...
    'SimDuty', [], 'SimInterferer', [], 'SimSweepTime', [], 'SimNoiseExcess', 0, 'File', "", ...
    'PlanOnly', false, 'Plot', true, 'Visible', 'on', 'OutDir', "", 'Prompt', true, 'Quiet', false, 'PDF', false);
if mod(numel(args), 2) ~= 0, error('ff_sat_tracker:option', 'Options must be Name, Value pairs.'); end
fn = fieldnames(opt);
for k = 1:2:numel(args)
    j = find(strcmpi(fn, string(args{k})), 1);
    if isempty(j)
        error('ff_sat_tracker:option', 'Unknown option "%s". Valid: %s', string(args{k}), strjoin(fn, ', '));
    end
    opt.(fn{j}) = args{k + 1};
end
for s = ["Downlink", "TLE", "SiteName", "TimeZone", "Mode", "Detector", "TxAntenna", "TxPointing", ...
        "RxAntenna", "Environment", "SaveTo", "File", "OutDir"]
    opt.(s) = string(opt.(s));
end
opt.Mode = lower(opt.Mode); opt.TxAntenna = lower(opt.TxAntenna); opt.RxAntenna = lower(opt.RxAntenna);
opt.TxPointing = lower(opt.TxPointing); opt.Environment = lower(opt.Environment); opt.Detector = upper(opt.Detector);
opt.Visible = char(string(opt.Visible));
if ~ismember(opt.Mode, ["retune", "wide"]), error('ff_sat_tracker:option', '''Mode'' must be "retune" or "wide".'); end
if ~ismember(opt.TxAntenna, ["auto", "isotropic", "dipole", "beam"]), error('ff_sat_tracker:option', '''TxAntenna'': auto | isotropic | dipole | beam.'); end
if ~ismember(opt.RxAntenna, ["tracking", "whip"]), error('ff_sat_tracker:option', '''RxAntenna'': tracking | whip.'); end
if ~ismember(opt.Detector, ["AVER", "SAMP", "POS", "NEG", "NORM", "AUTO"]), error('ff_sat_tracker:option', '''Detector'': AVER | SAMP | POS | NEG | NORM | AUTO.'); end
if opt.Points < 3 || opt.Points ~= round(opt.Points), error('ff_sat_tracker:option', '''Points'' must be an integer >= 3.'); end
for s = ["Wait", "Rehearse", "Preamp", "Offline", "Simulate", "PlanOnly", "Plot", "Prompt", "Quiet", "PDF", "LivePlot"]
    opt.(s) = logical(opt.(s));
end
if isduration(opt.ClockShift), opt.ClockShift = seconds(opt.ClockShift); end
if isduration(opt.Duration), opt.Duration = seconds(opt.Duration); end
if strlength(opt.OutDir) > 0 && ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
opt.Given = strings(1, 0);
for k = 1:2:numel(args), opt.Given(end+1) = fn{strcmpi(fn, string(args{k}))}; end
end

function say(opt, varargin)
if ~opt.Quiet, fprintf(varargin{:}); end
end

function d = supportDir(sub)
d = fullfile(fileparts(mfilename('fullpath')), 'support', 'tracker');
if nargin > 0, d = fullfile(d, sub); end
if ~isfolder(d), mkdir(d); end
end

function t = nowUTC(opt)
t = datetime('now', 'TimeZone', 'UTC') + seconds(opt.ClockShift);
end

function t = toUTC(t, zone)
if ~isdatetime(t), t = datetime(t); end
if isempty(t.TimeZone), t.TimeZone = char(zone); end
t.TimeZone = 'UTC';
end

function s = localStr(t, opt, fmt)
if nargin < 3, fmt = 'yyyy-MM-dd HH:mm:ss'; end
t.TimeZone = char(opt.TimeZone); s = string(t, fmt);
end

function s = srcName(src)
if isa(src, 'FieldFox'), s = src.IP; else, s = string(src); end
end

%% ======================================================================================
%% target and TLE
%% ======================================================================================
function [tgt, notes] = resolveTarget(opt, notes)
q = opt.Satellite; C = table(); S = table();
if exist('ff_satsig', 'file') == 2
    try
        [C, info] = ff_satsig("catalog"); S = info.satnogs;
    catch err
        notes(end+1) = "ff_satsig catalog not available: " + string(err.message);
    end
end
tgt = struct('name', "", 'norad', NaN, 'id', "", 'f0', NaN, 'obw', NaN, 'modulation', "", 'class', "", ...
    'template', "", 'eirp_dBm', NaN, 'eirpSource', "", 'burst', "", 'status', "", 'source', "");
isNum = isnumeric(q) || ~isnan(str2double(string(q)));
rows = [];
if ~isempty(C)
    if isNum
        rows = find(C.NORAD == double(string(q)));
    else
        qs = string(q);
        rows = find(strcmpi(C.Id, qs));
        if isempty(rows)
            qn = lower(regexprep(qs, '[\s_\-\(\)]', ''));
            hay = lower(regexprep(C.Satellite + " " + C.Id, '[\s_\-\(\)]', ''));
            rows = find(contains(hay, qn) & isfinite(C.NORAD));
            if ~isempty(rows), rows = find(C.NORAD == C.NORAD(rows(1))); end
        end
    end
    if strlength(opt.Downlink) > 0
        rows = rows(strcmpi(C.Id(rows), opt.Downlink));
        if isempty(rows), error('ff_sat_tracker:target', 'Downlink "%s" not found for "%s" (see ff_satsig("catalog")).', opt.Downlink, string(q)); end
    end
end
if ~isempty(rows)
    Ck = C(rows, :);
    inRange = Ck.F0_Hz >= 100e3 & Ck.F0_Hz <= 6.5e9;
    score = Ck.Active + 10*inRange;
    if ~isempty(opt.Frequency)
        score = score - abs(Ck.F0_Hz - opt.Frequency)/1e6;    % borrow BW/class of the nearest downlink
    end
    [~, i] = max(score); r = Ck(i, :);
    tgt.name = r.Satellite; tgt.norad = r.NORAD; tgt.id = r.Id; tgt.f0 = r.F0_Hz; tgt.obw = r.OBW_Hz;
    tgt.modulation = r.Modulation; tgt.class = r.Class; tgt.template = r.Template; tgt.eirp_dBm = r.EIRP_dBm;
    tgt.burst = r.Burst; tgt.status = r.Status; tgt.source = "ff_satsig catalog (" + r.Source + ")";
    if r.Active == 0, notes(end+1) = "Catalog status of " + r.Id + ": " + r.Status + " (no signal expected)."; end
    if height(Ck) > 1 && strlength(opt.Downlink) == 0
        notes(end+1) = "Downlink " + r.Id + " chosen; others: " + strjoin(Ck.Id(Ck.Id ~= r.Id), ", ") + " (use 'Downlink').";
    end
elseif ~isempty(S)
    if isNum, k = find(S.norad == double(string(q)));
    else, k = find(contains(lower(S.sat_name), lower(string(q)))); end
    if ~isempty(k)
        k = k(1); tgt.name = S.sat_name(k); tgt.norad = S.norad(k); tgt.f0 = S.downlink_low_Hz(k);
        tgt.modulation = S.mode(k); tgt.id = "satnogs:" + S.uuid(k);
        m = upper(S.mode(k)); b = S.baud(k);
        if contains(m, "CW")
            tgt.obw = 150; tgt.class = "CW"; tgt.template = "cw16";
        elseif contains(m, "AFSK")
            tgt.obw = 12e3; tgt.class = "AFSK-FM"; tgt.template = "afsk1k2";
        elseif any(m == ["FM", "FMN"])
            tgt.obw = 12e3; tgt.class = "FM"; tgt.template = "";      % Gaussian PSD model in 'Simulate'
        elseif isfinite(b) && b > 0
            tgt.obw = max(2*b, 2e3); tgt.class = "FSK/GMSK"; tgt.template = "";
        else
            tgt.obw = 10e3;
        end
        tgt.source = "SatNOGS DB subset (CC BY-SA 4.0); first alive transmitter";
        notes(end+1) = sprintf("%s: frequency from SatNOGS (%s, %s): check it.", tgt.name, fmtHz(tgt.f0), tgt.modulation);
    end
end
if isNum && ~isfinite(tgt.norad), tgt.norad = double(string(q)); tgt.name = "NORAD " + string(q); end
if ~isNum && ~isfinite(tgt.norad)
    tgt.name = string(q);
    tgt.norad = noradFromCaches(string(q));
    if ~isfinite(tgt.norad) && strlength(opt.TLE) == 0
        error('ff_sat_tracker:target', '"%s" is not in the catalog or the TLE caches: give the NORAD number or a ''TLE'' file.', string(q));
    end
end
if ~isempty(opt.Frequency), tgt.f0 = opt.Frequency; tgt.source = tgt.source + " (frequency from 'Frequency')"; end
if ~isempty(opt.SignalBW), tgt.obw = opt.SignalBW; end
if ~isfinite(tgt.f0), error('ff_sat_tracker:target', 'No downlink frequency known for "%s": give ''Frequency''.', string(q)); end
if ~isfinite(tgt.obw), tgt.obw = 1e3; notes(end+1) = "Signal bandwidth unknown: 1 kHz assumed ('SignalBW')."; end
if ~isempty(opt.EIRP)
    tgt.eirp_dBm = opt.EIRP; tgt.eirpSource = "'EIRP' option";
elseif isfinite(tgt.eirp_dBm)
    tgt.eirpSource = "catalog (often ASSUMED, see Status/Source)";
else
    tgt.eirp_dBm = 30; tgt.eirpSource = "30 dBm ASSUMED (no published EIRP)";
end
if tgt.f0 < 100e3 || tgt.f0 > 6.5e9
    notes(end+1) = sprintf("f0 = %s is outside the N9914A SA range (100 kHz - 6.5 GHz): use a down-converter / LNB and 'Frequency' = IF.", fmtHz(tgt.f0));
end
end

function n = noradFromCaches(name)
n = NaN; base = fileparts(mfilename('fullpath'));
dirs = [string(supportDir('tle')), fullfile(base, "support", "satsig", "tle"), fullfile(base, "support", "tle"), ...
    fullfile(base, "support", "satsig", "tle_sample"), string(supportDir())];
qn = regexprep(name, '[\s_\-\(\)]', '');
for d = dirs
    for f = dir(fullfile(d, '*.tle'))'
        E = parseTLE(fileread(fullfile(f.folder, f.name)));
        for e = E
            if strcmpi(regexprep(e.name, '[\s_\-\(\)]', ''), qn), n = e.norad; return, end
        end
    end
end
end

function E = parseTLE(txt)
L = strtrim(splitlines(string(txt))); L = L(strlength(L) > 0);
E = struct('name', {}, 'norad', {}, 'l1', {}, 'l2', {}, 'epoch', {});
for k = 1:numel(L) - 1
    a = char(L(k)); b = char(L(k + 1));
    if numel(a) >= 69 && numel(b) >= 69 && a(1) == '1' && b(1) == '2' && strcmp(a(3:7), b(3:7))
        nm = "SAT " + string(strtrim(a(3:7)));
        if k > 1 && ~startsWith(L(k - 1), ["1 ", "2 "]), nm = regexprep(L(k - 1), '^0 ', ''); end
        yy = str2double(a(19:20)); doy = str2double(a(21:32)); yr = 2000 + yy - 100*(yy >= 57);
        E(end+1) = struct('name', strtrim(string(nm)), 'norad', str2double(a(3:7)), 'l1', string(a(1:69)), ...
            'l2', string(b(1:69)), 'epoch', datetime(yr, 1, 1, 'TimeZone', 'UTC') + days(doy - 1)); %#ok<AGROW>
    end
end
end

function [tle, notes] = getTLE(tgt, opt, notes)
% CelesTrak CATNR (2 h cache, no retry for 2 h after an error) -> any cached copy -> shipped sample.
n = tgt.norad; base = fileparts(mfilename('fullpath'));
if strlength(opt.TLE) > 0
    if ~isfile(opt.TLE), error('ff_sat_tracker:tle', 'TLE file not found: %s', opt.TLE); end
    E = parseTLE(fileread(opt.TLE));
    if isempty(E), error('ff_sat_tracker:tle', 'No TLE in %s.', opt.TLE); end
    k = find([E.norad] == n, 1); if isempty(k), k = 1; end
    tle = E(k); tle.source = "TLE file " + opt.TLE; return
end
if ~isfinite(n) || n > 99999
    error('ff_sat_tracker:tle', 'NORAD %g: no TLE possible (CelesTrak gives no TLE above 99999; give an OMM via ''TLE'' after conversion).', n);
end
cache = fullfile(supportDir('tle'), sprintf('%d.tle', n)); errF = cache + ".err";
tle = [];
if isfile(cache) && fileAgeH(cache) < 2
    E = parseTLE(fileread(cache)); E = E([E.norad] == n);
    if ~isempty(E), tle = E(1); tle.source = sprintf("CelesTrak, cached %.0f min ago", 60*fileAgeH(cache)); end
end
if isempty(tle) && ~opt.Offline
    if isfile(errF) && fileAgeH(errF) < 2
        notes(end+1) = "CelesTrak failed less than 2 h ago (" + strtrim(string(fileread(errF))) + "): not retried.";
    else
        url = sprintf('https://celestrak.org/NORAD/elements/gp.php?CATNR=%d&FORMAT=TLE', n);
        try
            say(opt, 'Fetching %s\n', url);
            txt = webread(url, weboptions('ContentType', 'text', 'Timeout', 20, ...
                'UserAgent', 'ff_sat_tracker/1.0 (MATLAB R2026a; student CubeSat research, Bangkok)'));
            E = parseTLE(txt); E = E([E.norad] == n);
            if isempty(E), error('no TLE for %d in the reply: "%s"', n, extractBefore(string(txt) + string(blanks(60)), 61)); end
            writelines([E(1).name; E(1).l1; E(1).l2], cache);
            if isfile(errF), delete(errF); end
            tle = E(1); tle.source = "CelesTrak (downloaded now)";
        catch err
            writelines(string(datetime('now')) + " " + extractBefore(string(err.message) + newline, newline), errF);
            notes(end+1) = "CelesTrak download failed: " + extractBefore(string(err.message) + newline, newline);
        end
    end
end
if isempty(tle)
    % any cached copy (newest epoch): this tool, ff_satsig, ff_cubesat_link, then the shipped samples
    files = [dir(fullfile(supportDir('tle'), '*.tle')); dir(fullfile(base, 'support', 'satsig', 'tle', '*.tle')); ...
        dir(fullfile(base, 'support', 'tle', '*.tle')); dir(fullfile(supportDir(), '*.tle')); ...
        dir(fullfile(base, 'support', 'satsig', 'tle_sample', '*.tle'))];
    best = [];
    for f = files'
        E = parseTLE(fileread(fullfile(f.folder, f.name))); E = E([E.norad] == n);
        if ~isempty(E) && (isempty(best) || E(1).epoch > best.epoch)
            best = E(1); best.source = "cached/shipped file " + string(f.name) + " (" + string(f.folder) + ")";
        end
    end
    if isempty(best)
        error('ff_sat_tracker:tle', 'No TLE for NORAD %d (CelesTrak unavailable and no cached copy). Give ''TLE''.', n);
    end
    tle = best;
end
age = days(datetime('now', 'TimeZone', 'UTC') - tle.epoch);
if abs(age) > 7
    notes(end+1) = sprintf("TLE epoch %s is %.0f days old: pass times can be off by minutes, Doppler fits by tens of seconds.", string(tle.epoch, 'yyyy-MM-dd HH:mm'), age);
end
end

function h = fileAgeH(file)
h = Inf; d = dir(file);
if ~isempty(d), h = hours(datetime('now') - datetime(d.datenum, 'ConvertFrom', 'datenum')); end
end

function file = writeTLE(tle)
file = string(tempname) + ".tle";
writelines([tle.name; tle.l1; tle.l2], file);
end

%% ======================================================================================
%% passes and geometry
%% ======================================================================================
function T = findPasses(tle, tA, opt)
file = writeTLE(tle); c = onCleanup(@() delete(file));
sc = satelliteScenario(tA, tA + hours(opt.Hours), 10);
sat = satellite(sc, char(file), 'OrbitPropagator', 'sgp4');
gs = groundStation(sc, opt.Site(1), opt.Site(2), 'Altitude', opt.Site(3), 'MinElevationAngle', opt.MinElevation);
I = accessIntervals(access(sat(1), gs));
[~, el, ~, tt] = aer(gs, sat(1));
tt = tt(:); el = el(:); if isempty(tt.TimeZone), tt.TimeZone = 'UTC'; end
n = height(I);
AOS = NaT(n, 1, 'TimeZone', 'UTC'); LOS = AOS; TCA = AOS; MaxEl = NaN(n, 1);
for k = 1:n
    a = I.StartTime(k); b = I.EndTime(k);
    if isempty(a.TimeZone), a.TimeZone = 'UTC'; b.TimeZone = 'UTC'; end
    AOS(k) = a; LOS(k) = b;
    in = find(tt >= a & tt <= b);
    if isempty(in), continue, end
    [MaxEl(k), i] = max(el(in)); TCA(k) = tt(in(i));
end
Duration_min = minutes(LOS - AOS);
T = table(AOS, TCA, LOS, MaxEl, Duration_min);
end

function [pass, notes] = choosePass(T, tNow, opt, notes)
pass = [];
if isempty(T), notes(end+1) = sprintf("No pass above %g deg in the next %g h.", opt.MinElevation, opt.Hours); return, end
future = find(T.LOS > tNow + seconds(30));
p = opt.Pass;
if isnumeric(p)
    if p < 1 || p > numel(future), notes(end+1) = sprintf("Pass index %d: only %d passes in the window.", p, numel(future)); return, end
    k = future(p);
else
    switch lower(string(p))
        case "now"
            k = find(T.AOS <= tNow & T.LOS > tNow, 1);
            if isempty(k), notes(end+1) = "No pass in progress ('Pass', ""now"")."; return, end
        case "next"
            k = future(find(T.MaxEl(future) >= opt.MinPeakElevation, 1));
            if isempty(k), notes(end+1) = sprintf("No pass with max elevation >= %g deg in the next %g h.", opt.MinPeakElevation, opt.Hours); return, end
            low = future(future < k);
            if ~isempty(low), notes(end+1) = sprintf("%d lower pass(es) skipped ('MinPeakElevation' %g deg).", numel(low), opt.MinPeakElevation); end
        otherwise
            error('ff_sat_tracker:option', '''Pass'' must be "next", "now" or an index.');
    end
end
pass = table2struct(T(k, :));
pass.index = k;
end

function eph = ephemeris(tle, tA, tB, tgt, opt, preamp)
% 1 s geometry, Doppler (states/ECEF, checked with dopplershift) and the link prediction.
f0 = tgt.f0; c0 = physconst('LightSpeed');
file = writeTLE(tle); cl = onCleanup(@() delete(file));
tA = dateshift(tA, 'start', 'second'); tB = tA + seconds(ceil(seconds(tB - tA)));   % whole-second grid
sc = satelliteScenario(tA, tB, 1);
sat = satellite(sc, char(file), 'OrbitPropagator', 'sgp4'); sat = sat(1);
gs = groundStation(sc, opt.Site(1), opt.Site(2), 'Altitude', opt.Site(3), 'Name', char(opt.SiteName));
[p, v, tt] = states(sat, 'CoordinateFrame', 'ecef');
tt = tt(:); if isempty(tt.TimeZone), tt.TimeZone = 'UTC'; end
[az, el, rng] = aer(gs, sat);
eph.t0 = tt(1); eph.t = tt; eph.tau = seconds(tt - tt(1));
eph.p = p.'; eph.v = v.'; eph.g = lla2ecef(opt.Site(:).');
eph.az = az(:); eph.el = el(:); eph.rng = rng(:);
rho = eph.p - eph.g; eph.dfac = -sum(rho.*eph.v, 2)./vecnorm(rho, 2, 2)/c0;
eph.fD = f0*eph.dfac; eph.rate = gradient(eph.fD, eph.tau); eph.f0 = f0;
eph.speed = vecnorm(eph.v + cross(repmat([0 0 7.2921159e-5], numel(tt), 1), eph.p, 2), 2, 2);
eph.dsCheck_Hz = NaN;
try
    fDt = dopplershift(sat, gs, Frequency = f0); fDt = fDt(:);
    ok = isfinite(fDt) & eph.el > 0;
    if any(ok), eph.dsCheck_Hz = max(abs(fDt(ok) - eph.fD(ok))); end
catch
end
% satellite antenna and EIRP pattern (Satellite Communications + Antenna Toolbox)
kind = opt.TxAntenna;
if kind == "auto", if f0 > 2e9, kind = "beam"; else, kind = "isotropic"; end, end
lam = c0/f0;
switch kind
    case "isotropic"
        tx = transmitter(sat, 'Frequency', f0, 'Power', tgt.eirp_dBm - 30, 'SystemLoss', 0, ...
            'Antenna', phased.IsotropicAntennaElement);
        txt = "isotropic (tumbling / unknown attitude)";
    case "dipole"
        d = design(dipole, f0); pk = max(pattern(d, f0, 0:5:360, -90:5:90), [], 'all');
        tx = transmitter(sat, 'Frequency', f0, 'Power', tgt.eirp_dBm - pk - 30, 'SystemLoss', 0, ...
            'Antenna', d, 'MountingAngles', [0 0 90]);
        txt = sprintf("Antenna Toolbox dipole (%.1f dBi peak), nadir-stabilised, axis across the orbit plane", pk);
    case "beam"
        G = opt.TxGain; D = lam/pi*sqrt(10^(G/10)/0.6);
        if opt.TxPointing == "station", pointAt(sat, gs); end
        tx = transmitter(sat, 'Frequency', f0, 'Power', tgt.eirp_dBm - G - 30, 'SystemLoss', 0);
        gaussianAntenna(tx, 'DishDiameter', D, 'ApertureEfficiency', 0.6);
        txt = sprintf("Gaussian beam %.1f dBi (%.2f m), pointing %s", G, D, opt.TxPointing);
end
lk = link(tx, receiver(gs));
eph.rip_dBm = reshape(sigstrength(lk), [], 1) + 30;     % received isotropic power (no Rx antenna)
eph.eirp_dBm = reshape(txEIRP(lk), [], 1) + 30;
eph.txText = txt;
% receive side
[eph.grx, eph.rxText] = rxGain(eph.el, f0, opt);
eph.noise = noiseModel(f0, preamp, opt);
eph.C_dBm = eph.rip_dBm + eph.grx - opt.PolLoss - opt.MiscLoss + opt.LNAGain - opt.CableLoss;
eph.C_dBm(eph.el < 0) = NaN;
eph.CN0 = eph.C_dBm - eph.noise.N0in;
u = [true; diff(eph.tau) > 0];                                 % interp1 needs strictly increasing times
if ~all(u)
    for fn = ["t", "tau", "p", "v", "az", "el", "rng", "dfac", "fD", "rate", "speed", "rip_dBm", "eirp_dBm", "grx", "C_dBm", "CN0"]
        x = eph.(fn); eph.(fn) = x(u, :);
    end
end
end

function [g, txt] = rxGain(el, f0, opt)
switch opt.RxAntenna
    case "tracking"
        g = opt.RxGain*ones(size(el));
        txt = sprintf("tracking antenna, %.1f dBi toward the satellite", opt.RxGain);
    case "whip"
        d = design(dipole, f0); e = -90:1:90;
        pv = pattern(d, f0, 0, e); pv = double(pv(:));
        g = interp1(e(:), pv, max(min(el, 90), -90));
        txt = "fixed vertical whip (Antenna Toolbox dipole pattern, null at zenith)";
end
end

function N = noiseModel(f0, preamp, opt)
% Noise density at the FieldFox input: (antenna noise + LNA noise) x LNA gain / cable loss + DANL.
if isempty(opt.DANL), if preamp, danl = -153; else, danl = -137; end, else, danl = opt.DANL; end
fM = f0/1e6;
if ~isempty(opt.AntennaTemp)
    Ta = opt.AntennaTemp; txt = sprintf("antenna %g K (option)", Ta);
elseif fM <= 1000
    tab = struct('city', [76.8 27.7], 'residential', [72.5 27.7], 'rural', [67.2 27.7], 'quiet', [53.6 28.6]);
    if ~isfield(tab, opt.Environment), error('ff_sat_tracker:option', '''Environment'': city | residential | rural | quiet.'); end
    c = tab.(opt.Environment);
    Fa = max(c(1) - c(2)*log10(fM), 52 - 23*log10(fM));         % man-made vs galactic (P.372)
    Ta = max(290*10^(Fa/10), 30);
    txt = sprintf("ITU-R P.372 %s man-made noise Fa = %.1f dB (Ta %.0f K%s)", opt.Environment, Fa, Ta, ...
        ternary(fM > 250, ", extrapolated above 250 MHz", ""));
else
    Ta = 150; txt = "antenna 150 K assumed above 1 GHz (sky + ground spill-over)";
end
k = 1.380649e-23;
Tl = 0; if opt.LNAGain > 0, Tl = 290*(10^(opt.LNANoiseFigure/10) - 1); end
Next = 10*log10(k*(Ta + Tl)*1e3) + opt.LNAGain - opt.CableLoss;
N = struct('danl', danl, 'Ta', Ta, 'Tlna', Tl, 'Next', Next, 'N0in', 10*log10(10^(Next/10) + 10^(danl/10)), 'text', txt);
end

function v = interpEph(eph, tau, field)
v = interp1(eph.tau, eph.(field), tau, 'linear', NaN);
end

function d = dfacX(eph, tq, xtKm)
% Doppler factor -v_r/c at times tq (s from eph.t0); optional cross-track shift of the satellite (km).
tq = tq(:); c0 = 299792458;
P = interp1(eph.tau, eph.p, tq, 'spline'); V = interp1(eph.tau, eph.v, tq, 'spline');
if nargin > 2 && xtKm ~= 0
    Vi = V + cross(repmat([0 0 7.2921159e-5], numel(tq), 1), P, 2);
    nh = cross(P, Vi, 2); nh = nh./vecnorm(nh, 2, 2);
    P = P + 1e3*xtKm*nh;
end
rho = P - eph.g;
d = -sum(rho.*V, 2)./vecnorm(rho, 2, 2)/c0;
end

function pass = passFromEph(eph, tA, tB, opt)
in = eph.t >= tA & eph.t <= tB;
if ~any(in), in = true(size(eph.t)); end
el = eph.el; el(~in) = -Inf; [mx, i] = max(el);
vis = find(eph.el >= opt.MinElevation);
if isempty(vis), aos = tA; los = tB; else, aos = eph.t(vis(1)); los = eph.t(vis(end)); end
pass = struct('AOS', aos, 'TCA', eph.t(i), 'LOS', los, 'MaxEl', mx, 'Duration_min', minutes(los - aos), 'index', NaN);
end

%% ======================================================================================
%% capture plan
%% ======================================================================================
function plan = makePlan(tgt, tle, pass, tNow, opt)
notes = strings(0, 1);
f0 = tgt.f0; obw = tgt.obw;
tStart = max(tNow, pass.AOS - seconds(opt.Lead)); tEnd = pass.LOS + seconds(opt.Lead);
if ~isempty(opt.Duration), tEnd = min(tEnd, tStart + seconds(opt.Duration)); end
eph = ephemeris(tle, tStart - seconds(300), tEnd + seconds(300), tgt, opt, opt.Preamp);
% refine AOS/LOS/TCA from the 1 s geometry where the window covers them
vis = find(eph.el >= opt.MinElevation);
if ~isempty(vis) && vis(1) > 1, pass.AOS = eph.t(vis(1)); end
if ~isempty(vis) && vis(end) < numel(eph.t), pass.LOS = eph.t(vis(end)); end
[mx, i] = max(eph.el); if i > 1 && i < numel(eph.el), pass.TCA = eph.t(i); pass.MaxEl = mx; end
in = eph.t >= tStart & eph.t <= tEnd;
if ~any(in), in = true(size(eph.t)); end
maxDop = max(abs(eph.fD(in))); maxRate = max(abs(eph.rate(in)));
ageD = days(tStart - tle.epoch);
timeTol = opt.TimeTol; if isempty(timeTol), timeTol = 5 + 3*abs(ageD); end
tolF = f0*opt.FreqTolPPM*1e-6 + opt.FreqTolHz;
N = opt.Points; cyc = 0.3;
if opt.Mode == "retune"
    halfW = tolF + maxRate*(timeTol + opt.Step + cyc) + obw/2;
    span = niceUp(max(2*halfW + obw/2, 3*obw));
    center0 = (f0 + opt.FreqOffset)*(1 + interpEph(eph, seconds(tStart - eph.t0), 'dfac'));
else
    halfW = maxDop + tolF + obw/2;
    span = niceUp(2*halfW + obw);
    center0 = f0 + opt.FreqOffset;
end
if ~isempty(opt.Span), span = opt.Span; end
halfW = min(halfW, span/2);
rbw = rbwStep(max([2*span/(N - 1), obw/20, 10]));
if ~isempty(opt.RBW), rbw = opt.RBW; end
if isempty(opt.SimSweepTime), swt = min(max(0.03 + 0.5*span/rbw^2, 0.03), 3); else, swt = opt.SimSweepTime; end
smear = maxRate*swt;
if opt.Mode == "retune" && smear > rbw
    notes(end+1) = sprintf("Doppler moves up to %.0f Hz during one sweep (%.2f s est.) > RBW %s: lines smear near TCA.", smear, swt, fmtHz(rbw));
end
if opt.Mode == "wide" && span/rbw > 2e4
    notes(end+1) = "Wide span / RBW ratio is large: sweeps will be slow; consider 'Mode' ""retune"".";
end
if center0 - span/2 < 100e3 || center0 + span/2 > 6.5e9
    notes(end+1) = "Span reaches outside the 100 kHz - 6.5 GHz SA range.";
end
plan = struct('mode', opt.Mode, 'f0', f0, 'freqOffset', opt.FreqOffset, 'center0', center0, 'span', span, ...
    'rbw', rbw, 'points', N, 'step', opt.Step, 'detector', opt.Detector, 'preamp', opt.Preamp, ...
    'tStart', tStart, 'tEnd', tEnd, 'lead', opt.Lead, 'maxDoppler', maxDop, 'maxRate', maxRate, ...
    'halfW', halfW, 'tolF', tolF, 'timeTol', timeTol, 'tleAgeDays', ageD, 'sweepTimeEst', swt, ...
    'clockShift', opt.ClockShift, 'pass', pass, 'notes', notes, 'dsCheck_Hz', eph.dsCheck_Hz, 'eph', eph);
end

function plan = planFromRec(rec, eph, tgt, tle, pass, opt)
s = rec.settings;
in = eph.t >= min(rec.tStart) & eph.t <= max(rec.tEnd); if ~any(in), in = true(size(eph.t)); end
tolF = tgt.f0*opt.FreqTolPPM*1e-6 + opt.FreqTolHz;
halfW = s.span/2;
if isfield(s, 'searchHalf') && isfinite(s.searchHalf), halfW = s.searchHalf; end
plan = struct('mode', s.mode, 'f0', tgt.f0, 'freqOffset', s.freqOffset, 'center0', rec.fCenter(1), ...
    'span', s.span, 'rbw', s.rbw, 'points', numel(rec.f), 'step', s.step, 'detector', s.detector, ...
    'preamp', s.preamp, 'tStart', min(rec.tStart), 'tEnd', max(rec.tEnd), 'lead', NaN, ...
    'maxDoppler', max(abs(eph.fD(in))), 'maxRate', max(abs(eph.rate(in))), 'halfW', halfW, 'tolF', tolF, ...
    'timeTol', NaN, 'tleAgeDays', days(min(rec.tStart) - tle.epoch), 'sweepTimeEst', s.sweepTime, ...
    'clockShift', 0, 'pass', pass, 'notes', strings(0, 1), 'dsCheck_Hz', eph.dsCheck_Hz, 'eph', eph);
end

function printPlan(R, opt)
p = R.plan; t = R.target; s = R.pass;
fprintf('\nFF_SAT_TRACKER plan: %s (NORAD %d), downlink %s, OBW %s (%s)\n', t.name, t.norad, fmtHz(t.f0), fmtHz(t.obw), t.modulation);
fprintf('  TLE epoch %s UTC (%.1f d old at start), %s\n', string(R.tle.epoch, 'yyyy-MM-dd HH:mm'), p.tleAgeDays, R.tle.source);
fprintf('  Pass: AOS %s  TCA %s  LOS %s (%s), max el %.1f deg\n', localStr(s.AOS, opt, 'HH:mm:ss'), ...
    localStr(s.TCA, opt, 'HH:mm:ss'), localStr(s.LOS, opt, 'HH:mm:ss'), localStr(s.TCA, opt, 'yyyy-MM-dd z'), s.MaxEl);
fprintf('  Record %s .. %s; max Doppler +-%s, max rate %.1f Hz/s\n', localStr(p.tStart, opt, 'HH:mm:ss'), ...
    localStr(p.tEnd, opt, 'HH:mm:ss'), fmtHz(p.maxDoppler), p.maxRate);
fprintf('  Mode %s: span %s, RBW %s, %d points, step %g s, search +-%s, detector %s, preamp %d\n', p.mode, ...
    fmtHz(p.span), fmtHz(p.rbw), p.points, p.step, fmtHz(p.halfW), p.detector, p.preamp);
if ~isempty(R.passes)
    fprintf('  Passes in the window (local time):\n');
    for k = 1:height(R.passes)
        fprintf('   %2d  AOS %s  max el %5.1f deg  %4.1f min\n', k, localStr(R.passes.AOS(k), opt, 'MM-dd HH:mm:ss'), ...
            R.passes.MaxEl(k), R.passes.Duration_min(k));
    end
end
end

%% ======================================================================================
%% live acquisition
%% ======================================================================================
function [rec, notes, stopped] = acquireLive(src, plan, tgt, tle, opt, notes)
stopped = "";
if isa(src, 'FieldFox'), ff = src; else, ff = FieldFox(string(src), opt.Port); end
say(opt, 'Connected: %s\n', ff.IDN);
if opt.Prompt && usejava('desktop') && feature('ShowFigureWindows')
    a = questdlg(sprintf(['Take control of the FieldFox (%s) for the %s pass (AOS %s)? ' ...
        'The instrument state is backed up and restored afterwards. Nothing is transmitted.'], ff.IP, tgt.name, ...
        localStr(plan.pass.AOS, opt, 'HH:mm:ss')), 'ff_sat_tracker', 'Yes', 'No', 'Yes');
    if ~strcmp(a, 'Yes'), stopped = "declined"; rec = struct('P', []); notes(end+1) = "User declined."; return, end
end
% countdown until 'setup' time (take control ~20 s before the recording starts)
setup = 20; lastPrint = [];
while true
    remS = seconds(plan.tStart - nowUTC(opt)) - setup;
    if remS <= 0, break, end
    if isempty(lastPrint) || toc(lastPrint) >= ternary(remS > 120, 60, 10)
        say(opt, 'AOS %s local: recording starts in %s (max el %.0f deg)\n', localStr(plan.pass.AOS, opt, 'HH:mm:ss'), ...
            string(seconds(remS + setup), 'hh:mm:ss'), plan.pass.MaxEl);
        lastPrint = tic;
    end
    pause(min(remS, 1));
end
wasInControl = ff.InControl;
if ~wasInControl, ff.takeControl(); end
guard = onCleanup(@() releaseIf(ff, wasInControl));
N = plan.points;
ff.saConfigure(Center = plan.center0, Span = plan.span, RBW = plan.rbw, VBW = plan.rbw, Points = N, ...
    Preamp = plan.preamp, Detector = plan.detector, TraceType = "CLRW", RefLevel = opt.RefLevel);
if isempty(opt.Atten), ff.write("SENS:POW:RF:ATT:AUTO 1"); else, ff.write(sprintf("SENS:POW:RF:ATT %g", opt.Atten)); end
ff.write("SENS:POW:RF:EXTG 0");                                    % displayed levels = FieldFox input
e = ff.errors(); if ~isempty(e), notes(end+1) = "Instrument: " + strjoin(e, "; "); end
st = ff.saSettings();
span = st.stop - st.start; rbw = st.rbw; N = st.points;
swt = NaN;
try
    swt = ff.queryNum("SENS:SWE:TIME?");      % reported sweep time (in-sweep timing of the carrier bin)
catch
end
if isempty(swt), swt = NaN; end
if abs(span - plan.span) > 1 || abs(rbw - plan.rbw) > 0.5
    notes(end+1) = sprintf("FieldFox set span %s / RBW %s (asked %s / %s).", fmtHz(span), fmtHz(rbw), fmtHz(plan.span), fmtHz(plan.rbw));
end
frel = linspace(-span/2, span/2, N);
% recording buffers
nMax = max(100, ceil(seconds(plan.tEnd - plan.tStart)/0.05) + 10); nMax = min(nMax, 50000);
P = zeros(min(nMax, 2000), N, 'single'); fC = zeros(size(P, 1), 1); tS = NaT(size(P, 1), 1, 'TimeZone', 'UTC'); tE = tS;
cyc = zeros(size(P, 1), 1); swp = cyc; rd = cyc; ret = false(size(P, 1), 1);
if strlength(opt.SaveTo) > 0, saveFile = opt.SaveTo;
else
    d = opt.OutDir; if strlength(d) == 0, d = string(pwd); end
    saveFile = fullfile(d, "ff_sat_tracker_" + string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) + ".mat");
end
settings = struct('rbw', rbw, 'vbw', rbw, 'span', span, 'points', N, 'detector', plan.detector, 'preamp', plan.preamp, ...
    'atten', st.atten, 'ref', st.ref, 'sweepTime', swt, 'mode', plan.mode, 'step', plan.step, 'f0', tgt.f0, ...
    'frameCenter', ternary(plan.mode == "retune", tgt.f0, plan.center0), 'freqOffset', plan.freqOffset, ...
    'searchHalf', plan.halfW, 'extGain', 0, 'clockShift', opt.ClockShift, 'source', "live");
meta = trackerMeta(tgt, tle, plan, opt, "live");
% live waterfall
lp = [];
if opt.LivePlot && strcmp(opt.Visible, 'on') && feature('ShowFigureWindows')
    lp = livePlotInit(frel, span, tgt);
end
say(opt, 'Recording %s .. %s (local) in %s mode, centre %s, span %s, RBW %s\n', localStr(plan.tStart, opt, 'HH:mm:ss'), ...
    localStr(plan.tEnd, opt, 'HH:mm:ss'), plan.mode, fmtHz(plan.center0), fmtHz(span), fmtHz(rbw));
while nowUTC(opt) < plan.tStart, pause(0.05); end
k = 0; lastRet = -Inf; fc = st.start + span/2; lastSave = tic; lastMsg = tic; cycEst = 0.3;
t0 = plan.eph.t0; tmo = 15; if isfinite(swt), tmo = max(15, 5*swt + 10); end
try
    while nowUTC(opt) < plan.tEnd
        k = k + 1;
        if k > size(P, 1)                                       % grow buffers
            m = size(P, 1); P(2*m, N) = 0; fC(2*m) = 0; tS(2*m) = NaT; tE(2*m) = NaT; cyc(2*m) = 0; swp(2*m) = 0; rd(2*m) = 0; ret(2*m) = false;
        end
        c1 = tic; tNow = nowUTC(opt);
        if plan.mode == "retune" && seconds(tNow - t0) - lastRet >= plan.step
            tq = seconds(tNow - t0) + cycEst/2;
            fc = round((tgt.f0 + plan.freqOffset)*(1 + interpEph(plan.eph, tq, 'dfac')));
            ff.write(sprintf("SENS:FREQ:CENT %.0f", fc));
            lastRet = seconds(tNow - t0); ret(k) = true;
        end
        tS(k) = nowUTC(opt);
        ff.query("INIT:IMM;*OPC?", tmo);
        tE(k) = nowUTC(opt); swp(k) = seconds(tE(k) - tS(k));
        y = ff.saTrace();
        rd(k) = seconds(nowUTC(opt) - tE(k));
        if k == 1 || mod(k, 200) == 0
            e = ff.errors(); if ~isempty(e), notes(end+1) = sprintf("Instrument at sweep %d: %s", k, strjoin(e, "; ")); end %#ok<AGROW>
        end
        if numel(y) ~= N, error('ff_sat_tracker:points', 'Trace has %d points, expected %d.', numel(y), N); end
        P(k, :) = single(y); fC(k) = fc;
        cyc(k) = toc(c1); cycEst = 0.8*cycEst + 0.2*cyc(k);
        if ~isempty(lp), livePlotUpdate(lp, P, k, tS, plan, tgt, fC); end
        if toc(lastMsg) > 10
            tq = seconds(tE(k) - t0);
            [mx, j] = max(y);
            say(opt, '%s  el %5.1f deg  Doppler %+8.0f Hz  sweep %4d  cycle %.2f s  peak %.1f dBm at %+.0f Hz (floor %.1f)\n', ...
                localStr(tE(k), opt, 'HH:mm:ss'), interpEph(plan.eph, tq, 'el'), interpEph(plan.eph, tq, 'fD'), k, cycEst, mx, ...
                fc + frel(j) - (tgt.f0 + plan.freqOffset)*(1 + interpEph(plan.eph, tq, 'dfac')), median(y));
            lastMsg = tic;
        end
        if toc(lastSave) > opt.SaveEvery
            saveRec(saveFile, packRec(P, fC, tS, tE, cyc, swp, rd, ret, k, frel, settings, ff.IDN, meta));
            lastSave = tic;
        end
    end
catch err
    notes(end+1) = sprintf("Acquisition stopped at sweep %d: %s", k, err.message);
    warning('ff_sat_tracker:acquire', 'Acquisition stopped at sweep %d: %s', k, err.message);
    k = k - 1;                                              % the failed sweep has no data
end
clear guard                                                 % restore the instrument now
rec = packRec(P, fC, tS, tE, cyc, swp, rd, ret, max(k, 0), frel, settings, ff.IDN, meta);
if k > 0
    saveRec(saveFile, rec);
    notes(end+1) = "Recording saved: " + saveFile;
end
end

function releaseIf(ff, wasInControl)
if wasInControl, return, end
try
    ff.release();
catch err
    warning('ff_sat_tracker:release', 'Release failed: %s', err.message);
end
end

function meta = trackerMeta(tgt, tle, plan, opt, source)
p = rmfield(plan, 'eph');
meta = struct('version', "ff_sat_tracker 1.0", 'source', source, 'target', tgt, 'tle', tle, ...
    'site', struct('lla', opt.Site, 'name', opt.SiteName), 'pass', plan.pass, 'plan', p, ...
    'link', struct('EIRP_dBm', tgt.eirp_dBm, 'TxAntenna', opt.TxAntenna, 'TxGain', opt.TxGain, ...
    'TxPointing', opt.TxPointing, 'RxAntenna', opt.RxAntenna, 'RxGain', opt.RxGain, 'PolLoss', opt.PolLoss, ...
    'MiscLoss', opt.MiscLoss, 'LNAGain', opt.LNAGain, 'LNANoiseFigure', opt.LNANoiseFigure, ...
    'CableLoss', opt.CableLoss, 'Environment', opt.Environment));
end

function rec = packRec(P, fC, tS, tE, cyc, swp, rd, ret, k, frel, settings, idn, meta)
i = 1:k;
rec = struct();
rec.f = settings.frameCenter + frel;
rec.P = P(i, :);
rec.t = tS(i) + (tE(i) - tS(i))/2;
rec.settings = settings;
rec.idn = string(idn);
rec.fCenter = fC(i);
rec.tStart = tS(i); rec.tEnd = tE(i);
rec.timing = struct('cycle_s', cyc(i), 'sweep_s', swp(i), 'read_s', rd(i), 'retuned', ret(i));
rec.tracker = meta;
end

function saveRec(file, rec)
try
    S = rec;
    save(file, '-struct', 'S', '-v7');
catch err
    warning('ff_sat_tracker:save', 'Could not save %s: %s', file, err.message);
end
end

function lp = livePlotInit(frel, span, tgt)
lp.fig = figure('Name', 'ff_sat_tracker live', 'NumberTitle', 'off', 'Color', 'w');
lp.ax = axes(lp.fig);
lp.im = imagesc(lp.ax, frel/1e3, 1:300, NaN(300, numel(frel)));
set(lp.ax, 'YDir', 'normal'); colormap(lp.ax, turbo); colorbar(lp.ax);
xlabel(lp.ax, 'Offset from the sweep centre (kHz)'); ylabel(lp.ax, 'Sweep (last 300)');
title(lp.ax, sprintf('%s %s - tracking frame', tgt.name, fmtHz(tgt.f0)));
lp.span = span;
end

function livePlotUpdate(lp, P, k, tS, plan, tgt, fC) %#ok<INUSD>
if ~isvalid(lp.fig), return, end
a = max(1, k - 299); Z = NaN(300, size(P, 2)); Z(end - (k - a):end, :) = P(a:k, :);
lp.im.CData = Z;
lo = prctile(Z(:), 5); hi = max(Z(:)); if isfinite(lo) && hi > lo, clim(lp.ax, [lo, hi]); end
title(lp.ax, sprintf('%s  %s  centre %s', tgt.name, string(tS(k), 'HH:mm:ss'), fmtHz(fC(k))));
drawnow limitrate
end

%% ======================================================================================
%% simulation
%% ======================================================================================
function [rec, sim] = simulateRec(plan, tgt, tle, opt)
rng(opt.Seed);
eph = plan.eph; f0 = tgt.f0; N = plan.points; span = plan.span; rbw = plan.rbw; enbw = 1.0645*rbw;
frel = linspace(-span/2, span/2, N);
dfo = opt.SimFreqOffset; if isempty(dfo), dfo = round(2.5e-6*f0); end
truth = struct('dt_s', opt.SimTimeError, 'df_Hz', dfo, 'ppm', dfo/f0*1e6, 'xtrack_km', opt.SimCrossTrack, ...
    'eirpError_dB', opt.SimEIRPError, 'fading_dB', opt.SimFading, 'noiseExcess_dB', opt.SimNoiseExcess);
[u, D1, shapeTxt] = signalShape(tgt, rbw);
floorLin = 10^((eph.noise.N0in + 10*log10(enbw) + opt.SimNoiseExcess)/10);
Tsw = plan.sweepTimeEst;
M = min(max(round(Tsw/N*rbw), 1), 64);                      % RBW-independent samples per bucket (AVER)
% slow fading (1 s AR(1), 5 s correlation)
nF = numel(eph.tau); a = exp(-1/5); w = randn(nF, 1); fade = zeros(nF, 1); fade(1) = w(1);
for i = 2:nF, fade(i) = a*fade(i - 1) + sqrt(1 - a^2)*w(i); end
fade = opt.SimFading*fade;
% keying: CW = Morse-like 75 ms elements, 55 % key-down; 'SimDuty' [on period] overrides
duty = opt.SimDuty; burstTxt = "";
if isempty(duty) && tgt.class ~= "CW" && isBursty(tgt.burst)
    bt = string(tgt.burst);
    t1 = regexp(bt, 'every\s+(\d+(?:\.\d+)?)\s*s', 'tokens', 'once');
    t2 = regexp(bt, '(\d+(?:\.\d+)?)\s*ms bursts in (?:a )?(\d+(?:\.\d+)?)\s*ms', 'tokens', 'once');
    if ~isempty(t2)
        duty = [str2double(t2(1)), str2double(t2(2))]/1e3;
    elseif ~isempty(t1)
        duty = [0.5, str2double(t1(1))];
    else
        duty = [1, 30];
    end
    burstTxt = sprintf("; bursts %.3g s every %g s (from the catalog Burst text '%s'; length ASSUMED)", duty(1), duty(2), bt);
end
keyed = isempty(duty) && tgt.class == "CW";
if keyed, keyEl = 0.075; nK = ceil((eph.tau(end) + 10)/keyEl); keyState = rand(nK, 1) < 0.55; end
dutyPhase = rand*10;
tA = seconds(plan.tStart - eph.t0); tB = seconds(plan.tEnd - eph.t0);
nMax = ceil((tB - tA)/(Tsw + 0.05)) + 5;
P = zeros(nMax, N, 'single'); fC = zeros(nMax, 1); tS = zeros(nMax, 1); tE = tS; cyc = tS; swp = tS; rd = tS; ret = false(nMax, 1);
t = tA; k = 0; lastRet = -Inf; fc = plan.center0;
while t < tB && k < nMax
    k = k + 1;
    pre = 0.004 + 0.002*rand;                               % trigger latency
    if plan.mode == "retune" && t - lastRet >= plan.step
        fc = round((f0 + plan.freqOffset)*(1 + interpEph(eph, t + pre + Tsw/2, 'dfac')));
        lastRet = t; ret(k) = true; pre = pre + 0.003;
    end
    tp = t + pre + (0:N-1)/(N-1)*Tsw;
    ftrue = (f0 + dfo)*(1 + dfacX(eph, tp + truth.dt_s, truth.xtrack_km)).';
    C = interp1(eph.tau, eph.C_dBm, tp, 'linear', NaN) + truth.eirpError_dB + interp1(eph.tau, fade, tp, 'linear', 0);
    C(~isfinite(C)) = -300;
    on = ones(1, N);
    if keyed, on = double(keyState(min(floor(tp/keyEl) + 1, nK)).');
    elseif ~isempty(duty), on = double(mod(tp + dutyPhase, duty(2)) < duty(1)); end
    fa = fc + frel;
    sig = 10.^(C/10).*interp1(u, D1, fa - ftrue, 'linear', 0).*on;
    if ~isempty(opt.SimInterferer)
        sig = sig + 10^(opt.SimInterferer(2)/10)*exp(-4*log(2)*((fa - f0 - opt.SimInterferer(1))/rbw).^2);
    end
    noise = floorLin*mean(-log(rand(M, N)), 1);
    if plan.detector == "POS", noise = floorLin*max(-log(rand(M, N)), [], 1); end
    if plan.detector == "SAMP", noise = floorLin*(-log(rand(1, N))); end
    P(k, :) = single(10*log10(sig + noise));
    fC(k) = fc; tS(k) = t + pre; tE(k) = t + pre + Tsw;
    rd(k) = 0.02 + 0.01*rand; swp(k) = Tsw + pre; cyc(k) = swp(k) + rd(k) + 0.002;
    t = t + cyc(k);
end
i = 1:k;
settings = struct('rbw', rbw, 'vbw', rbw, 'span', span, 'points', N, 'detector', plan.detector, 'preamp', plan.preamp, ...
    'atten', 0, 'ref', opt.RefLevel, 'sweepTime', Tsw, 'mode', plan.mode, 'step', plan.step, 'f0', f0, ...
    'frameCenter', ternary(plan.mode == "retune", f0, plan.center0), 'freqOffset', plan.freqOffset, ...
    'searchHalf', plan.halfW, 'extGain', 0, 'clockShift', 0, 'source', "simulated");
rec = struct();
rec.f = settings.frameCenter + frel; rec.P = P(i, :);
rec.tStart = eph.t0 + seconds(tS(i)); rec.tEnd = eph.t0 + seconds(tE(i));
rec.t = rec.tStart + (rec.tEnd - rec.tStart)/2;
rec.settings = settings; rec.idn = "SIMULATED (ff_sat_tracker)"; rec.fCenter = fC(i);
rec.timing = struct('cycle_s', cyc(i), 'sweep_s', swp(i), 'read_s', rd(i), 'retuned', ret(i));
rec.tracker = trackerMeta(tgt, tle, plan, opt, "simulated");
rec.truth = truth;
sim = struct('truth', truth, 'shape', shapeTxt + burstTxt, 'keyed', keyed, 'duty', duty, 'samplesPerBucket', M);
end

function [u, D1, txt] = signalShape(tgt, rbw)
% Displayed power vs offset u for 1 mW of signal through a Gaussian RBW (|H|^2 = 0.5 at +-RBW/2).
file = fullfile(fileparts(mfilename('fullpath')), 'support', 'satsig', 'satsig_templates.mat');
key = tgt.template; T = [];
if strlength(key) > 0 && isfile(file)
    S = load(file, 'tpl');
    if isfield(S.tpl, key), T = S.tpl.(key); end
end
if ~isempty(T)
    f = double(T.f(:)); p = max(double(T.psd(:)), 0);
    keep = find(p >= 1e-5*max(p));
    lo = f(keep(1)) - 4*rbw; hi = f(keep(end)) + 4*rbw;
    h = min(T.df, rbw/5); n = ceil((hi - lo)/h) + 1;
    if n > 4e5, h = (hi - lo)/4e5; n = 400001; end
    u = lo + (0:n-1)'*h;
    pu = interp1(f, p, u, 'linear', 0); pu = pu/sum(pu);
    if rbw > 2*h
        m = ceil(3*rbw/h); g = exp(-4*log(2)*((-m:m)'*h/rbw).^2);
        D1 = conv(pu, g, 'same');
    else
        D1 = pu/h*1.0645*rbw;
    end
    txt = "ff_satsig template " + key + " (" + string(T.note) + ")";
else
    if tgt.class == "CW" || tgt.obw <= rbw
        u = (-4*rbw:rbw/40:4*rbw)'; D1 = exp(-4*log(2)*(u/rbw).^2); txt = "pure tone through the RBW";
    else
        s = tgt.obw/4; h = min(rbw, s)/10; u = (-(3*tgt.obw + 4*rbw):h:(3*tgt.obw + 4*rbw))';
        pu = exp(-u.^2/(2*s^2)); pu = pu/sum(pu);
        m = ceil(3*rbw/h); g = exp(-4*log(2)*((-m:m)'*h/rbw).^2); D1 = conv(pu, g, 'same');
        txt = sprintf("Gaussian PSD model (sigma = OBW/4 = %s)", fmtHz(s));
    end
end
end

function sim = simCompare(sim, A)
t = sim.truth; b = A.best;
sim.recovered = struct('dt_s', NaN, 'df_Hz', NaN, 'xtrack_km', NaN, 'eirp_dB', NaN);
if ~isempty(b)
    sim.recovered.dt_s = b.dt_s; sim.recovered.df_Hz = b.df_Hz; sim.recovered.xtrack_km = b.xtrack_km;
    sim.err = struct('dt_s', b.dt_s - t.dt_s, 'df_Hz', b.df_Hz - t.df_Hz, ...
        'dt_sigmas', (b.dt_s - t.dt_s)/b.dt_sigma_s, 'df_sigmas', (b.df_Hz - t.df_Hz)/b.df_sigma_Hz);
end
if isfield(A.link, 'offset_dB'), sim.recovered.eirp_dB = A.link.offset_dB; end
end

%% ======================================================================================
%% recordings from file
%% ======================================================================================
function [rec, tgt, tle, notes] = loadRecording(opt, notes)
if ~isfile(opt.File), error('ff_sat_tracker:file', 'File not found: %s', opt.File); end
S = load(opt.File);
for v = ["f", "P", "t"]
    if ~isfield(S, v), error('ff_sat_tracker:file', '%s has no "%s" (suite format f, P, t, settings).', opt.File, v); end
end
rec = struct(); rec.f = double(S.f(:)).'; rec.P = S.P;
if size(rec.P, 2) ~= numel(rec.f) && size(rec.P, 1) == numel(rec.f), rec.P = rec.P.'; end
K = size(rec.P, 1);
t = S.t(:); if ~isdatetime(t), error('ff_sat_tracker:file', 't must be a datetime vector.'); end
t = toUTC(t, opt.TimeZone); rec.t = t;
s = struct(); if isfield(S, 'settings'), s = S.settings; end
if ~isfield(s, 'rbw') || ~isfinite(s.rbw)
    if ~isempty(opt.RBW), s.rbw = opt.RBW; else, s.rbw = 2*(rec.f(end) - rec.f(1))/(numel(rec.f) - 1); end
    notes(end+1) = sprintf("RBW not in the file: %s assumed ('RBW').", fmtHz(s.rbw));
end
def = struct('vbw', s.rbw, 'span', rec.f(end) - rec.f(1), 'points', numel(rec.f), 'detector', "AVER", 'preamp', true, ...
    'atten', NaN, 'ref', NaN, 'sweepTime', NaN, 'mode', "wide", 'step', NaN, 'f0', NaN, ...
    'frameCenter', mean(rec.f([1 end])), 'freqOffset', 0, 'searchHalf', NaN, 'extGain', 0, 'clockShift', 0, 'source', "file");
for fn = string(fieldnames(def))'
    if ~isfield(s, fn), s.(fn) = def.(fn); end
end
s.detector = upper(string(s.detector)); s.mode = string(s.mode); s.preamp = logical(s.preamp);
if isfield(S, 'fCenter'), rec.fCenter = double(S.fCenter(:)); else, rec.fCenter = repmat(s.frameCenter, K, 1); end
if isfield(S, 'tStart') && isfield(S, 'tEnd')
    rec.tStart = toUTC(S.tStart(:), opt.TimeZone); rec.tEnd = toUTC(S.tEnd(:), opt.TimeZone);
else
    rec.tStart = t; rec.tEnd = t;
    notes(end+1) = "No sweep start/end times in the file: the in-sweep time of the carrier is not corrected.";
end
if isfield(S, 'timing'), rec.timing = S.timing; else, rec.timing = struct('cycle_s', [diff(seconds(t - t(1))); NaN], 'sweep_s', NaN(K, 1), 'read_s', NaN(K, 1), 'retuned', false(K, 1)); end
rec.idn = ""; if isfield(S, 'idn'), rec.idn = string(S.idn); end
if isfield(S, 'truth'), rec.truth = S.truth; end
useMeta = isfield(S, 'tracker') && isfield(S.tracker, 'target') && ~ismember("Satellite", opt.Given);
if useMeta
    tgt = S.tracker.target; tle = S.tracker.tle;
    if ~isempty(opt.EIRP), tgt.eirp_dBm = opt.EIRP; tgt.eirpSource = "'EIRP' option"; end
    if ~isempty(opt.SignalBW), tgt.obw = opt.SignalBW; end
    notes(end+1) = "Target and TLE from the recording (" + tle.source + ").";
    if strlength(opt.TLE) > 0, [tle, notes] = getTLE(tgt, opt, notes); end
else
    [tgt, notes] = resolveTarget(opt, notes);
    [tle, notes] = getTLE(tgt, opt, notes);
end
if ~isfinite(s.f0), s.f0 = tgt.f0; end
rec.settings = s;
rec.tracker = struct('target', tgt, 'tle', tle);
if isfield(S, 'tracker'), rec.tracker = S.tracker; rec.tracker.target = tgt; rec.tracker.tle = tle; end
end

%% ======================================================================================
%% analysis
%% ======================================================================================
function A = analyzeRec(rec0, eph, tgt, plan, opt)
notes = strings(0, 1);
s = rec0.settings; N = numel(rec0.f);
f0 = tgt.f0; rbw = s.rbw; enbw = 1.0645*rbw; obw = tgt.obw;
frel = rec0.f - s.frameCenter; df = (frel(end) - frel(1))/(N - 1);
wide = obw > 3*rbw;
[rec, Tint, nSw, dur] = integrateSweeps(rec0, tgt, plan, opt, wide);
K = size(rec.P, 1);
if Tint > 0
    notes(end+1) = sprintf("Sweeps averaged over %.1f s blocks before detection (%d blocks, median %d sweeps each; 'Integrate').", ...
        Tint, K, round(median(nSw)));
end
Tsw = dur; if isfinite(s.sweepTime) && s.sweepTime > 0, Tsw = min(dur, s.sweepTime); end
tauMid = seconds(rec.t - eph.t0);
fExpMid = (f0 + s.freqOffset)*(1 + interpEph(eph, tauMid, 'dfac'));
halfW = plan.halfW; if ~isfinite(halfW), halfW = s.span/2; end
detCorr = 0; detOK = ismember(s.detector, ["AVER", "SAMP"]);    % floor = mean noise power (unbiased for both)
if ~detOK, notes(end+1) = "Detector " + s.detector + ": the noise floor (and C/N0) is biased; use AVER for calibrated C/N0."; end
% ---- per-sweep (per-block) detection ----
% fixed (terrestrial/EMI) carriers inside the window are found (constant absolute frequency while the
% predicted Doppler moves) and masked, then the detection is repeated (up to 3 lines)
lines = zeros(0, 1); lineW = max(1.5*rbw, 2*df); tpl = [];
if wide, [tu, tD] = signalShape(tgt, rbw); tpl = struct('u', tu, 'D', tD); end
for it = 1:4
    [fMeas, snr, flo, chP, wid, tdet, ok, zs, sgr] = detectAll(rec, frel, fExpMid, halfW, rbw, obw, opt, wide, ...
        tauMid, Tsw, lines, lineW, tpl);
    fl = fixedLine(fMeas, ok, tdet, eph, f0, rbw, df, lines);
    if isnan(fl) || it == 4, break, end
    lines(end+1, 1) = fl; %#ok<AGROW>
end
if ~isempty(lines)
    notes(end+1) = sprintf("Fixed carrier(s) at %s (constant frequency while the predicted Doppler moved): masked +-%s as terrestrial/EMI.", ...
        strjoin(arrayfun(@fmtHz, lines), ", "), fmtHz(lineW));
end
fExp = (f0 + s.freqOffset)*(1 + interpEph(eph, tdet, 'dfac'));
el = interpEph(eph, tdet, 'el'); az = interpEph(eph, tdet, 'az'); rngKm = interpEph(eph, tdet, 'rng')/1e3;
% gross outliers vs a running median of the residual to the prediction
res0 = fMeas - fExp; good = ok;
if nnz(ok) >= 5
    r = res0(ok); rm = movmedian(r, 9); dv = abs(r - rm);
    sc = max(1.4826*median(dv), max(df, rbw/3));
    g = dv <= 6*sc; idx = find(ok); good(idx(~g)) = false;
end
% ---- fits ----
F = fitDoppler(tdet(good), fMeas(good), eph, f0, s.freqOffset, opt, max(df, rbw/4));
keepAll = false(K, 1); gi = find(good); if ~isempty(F.keep), keepAll(gi(F.keep)) = true; end
fFit = NaN(K, 1); if ~isempty(F.best), fFit = F.best.model(tdet); end
% ---- link: measured vs predicted ----
N0meas = flo - 10*log10(enbw) + detCorr;
CN0meas = chP - N0meas;
CN0pred = interpEph(eph, tdet, 'CN0');
floorPred = eph.noise.N0in + 10*log10(enbw) - detCorr;
Lk = struct('N0in_pred', eph.noise.N0in, 'floorPred_dBm', floorPred, 'floorMeas_dBm', median(flo, 'omitnan'), ...
    'noiseExcess_dB', median(flo, 'omitnan') - floorPred, 'noise', eph.noise, 'tx', eph.txText, 'rx', eph.rxText, ...
    'EIRP_dBm', tgt.eirp_dBm, 'EIRPsource', tgt.eirpSource, 'offset_dB', NaN, 'offsetIQR_dB', NaN, 'n', 0);
use = keepAll & isfinite(CN0meas) & isfinite(CN0pred);
% near the threshold only the up-faded sweeps are detected (selection bias): use sweeps predicted
% >= 6 dB above the detection threshold for the offset when there are enough of them
if wide        % box sum of C over N0*OBW must exceed ZMin x its relative noise std
    cn0Thr = 10*log10(obw*opt.ZMin*median(sgr, 'omitnan'));
else
    cn0Thr = opt.SNRMin + 10*log10(max(enbw, obw));
end
Lk.cn0Threshold = cn0Thr; strong = use & CN0pred >= cn0Thr + 6;
Lk.offsetBasis = "sweeps predicted >= 6 dB above the detection threshold";
if nnz(strong) < 10, strong = use; Lk.offsetBasis = "all used sweeps (few strong ones: biased high near the threshold)"; end
if nnz(strong) >= 3
    dlt = CN0meas(strong) - CN0pred(strong);
    Lk.offset_dB = median(dlt); Lk.offsetIQR_dB = iqr(dlt); Lk.n = nnz(strong);
end
% predicted detectability (SNR of the 3-bin smoothed peak, approx.) vs elevation
Lk.elBins = (0:10:80)'; Lk.cn0MeasBin = NaN(9, 1); Lk.cn0PredBin = NaN(9, 1); Lk.detFrac = NaN(9, 1);
elAll = interpEph(eph, tauMid, 'el');
for b = 1:9
    inb = elAll >= Lk.elBins(b) & elAll < Lk.elBins(b) + 10;
    if any(inb), Lk.detFrac(b) = mean(keepAll(inb)); end
    u2 = use & el >= Lk.elBins(b) & el < Lk.elBins(b) + 10;
    if any(u2), Lk.cn0MeasBin(b) = median(CN0meas(u2)); Lk.cn0PredBin(b) = median(CN0pred(u2)); end
end
% ---- timing ----
tm = rec0.timing;
Tm = struct('sweeps', size(rec0.P, 1), 'units', K, 'integrate_s', Tint, 'cycle_mean_s', mean(tm.cycle_s, 'omitnan'), 'cycle_p95_s', prctile(tm.cycle_s, 95), ...
    'sweep_mean_s', mean(tm.sweep_s, 'omitnan'), 'read_mean_s', mean(tm.read_s, 'omitnan'), ...
    'overhead_mean_s', mean(tm.cycle_s - tm.sweep_s, 'omitnan'), 'retunes', nnz(tm.retuned), ...
    'maxGap_s', max([0; diff(seconds(rec0.t - rec0.t(1)))]), 'reportedSweepTime_s', s.sweepTime, ...
    'smear_Hz', plan.maxRate*mean(Tsw, 'omitnan'));
% ---- verdict ----
nG = nnz(keepAll); swingPred = NaN; swingMeas = NaN;
if nG >= 2, swingMeas = max(fMeas(keepAll)) - min(fMeas(keepAll)); swingPred = max(fExp(keepAll)) - min(fExp(keepAll)); end
gate = ternary(wide, sprintf("%g sigma (box filter)", opt.ZMin), sprintf("%g dB SNR", opt.SNRMin));
if nG < 5 || isempty(F.best)
    verdict = sprintf("INCONCLUSIVE: %d usable detections (need >= 5). No signal above %s in the search window, or too short a capture.", nG, gate);
else
    rf = F.best.rms_Hz; rx = F.rmsFixed_Hz;
    tk = tdet(keepAll); fe = fExp(keepAll);                    % curvature of the predicted Doppler that a
    cl = polyfit(tk - mean(tk), fe, 1);                         % fixed carrier + drift cannot follow
    curv = sqrt(mean((fe - polyval(cl, tk - mean(tk))).^2));
    stats = sprintf("residual %s RMS vs %s for a fixed carrier; predicted Doppler curvature over the detections %s RMS", ...
        fmtHz(rf), fmtHz(rx), fmtHz(curv));
    est = sprintf("dt = %+.2f +- %.2f s, df = %+.0f +- %.0f Hz (%+.2f ppm)", F.best.dt_s, F.best.dt_sigma_s, F.best.df_Hz, F.best.df_sigma_Hz, F.best.ppm);
    if curv < 2*min(rf, rx)
        verdict = sprintf("INCONCLUSIVE: the Doppler curvature over the detections is below the measurement scatter (%s): longer/higher pass, more 'Integrate', or a narrower signal.", stats);
    elseif rf < rx/3
        verdict = sprintf("DOPPLER SIGNATURE CONFIRMED: %d detections follow the TLE Doppler (%s); %s.", nG, stats, est);
    elseif rx < rf/2
        verdict = sprintf("NOT CONFIRMED: a fixed carrier fits much better (%s): terrestrial emitter / EMI / another satellite.", stats);
    elseif rf < rx/1.5 && nG >= 10
        verdict = sprintf("PROBABLE: the TLE Doppler model fits better than a fixed carrier but the scatter is large (%s); %s.", stats, est);
    else
        verdict = sprintf("NOT CONFIRMED: %s: the detections may be another emitter (terrestrial/EMI) or a different satellite.", stats);
    end
end
if nG >= 5 && ~isempty(F.best) && (abs(F.best.dt_s) > 0.9*opt.MaxTimeOffset)
    notes(end+1) = "The fitted time offset sits at the bound 'MaxTimeOffset': TLE much older/wrong, or PC clock error.";
end
if ~isempty(F.best)
    F.best.raw = F.raw; F.best.rmsFixed_Hz = F.rmsFixed_Hz; F.best.swingMeas_Hz = swingMeas; F.best.swingPred_Hz = swingPred;
end
track = table(rec.t, tdet, el, az, rngKm, fExp, fMeas, fFit, snr, zs, flo, chP, CN0meas, CN0pred, wid, nSw, ok, keepAll, ...
    'VariableNames', ["t", "tau_s", "el_deg", "az_deg", "range_km", "fPred_Hz", "fMeas_Hz", "fFit_Hz", "SNR_dB", "Z", ...
    "floor_dBm", "chanPow_dBm", "CN0meas_dBHz", "CN0pred_dBHz", "width_Hz", "nSweeps", "detected", "used"]);
A = struct('track', track, 'fits', F.table, 'best', [], 'fit', F, 'link', Lk, 'timing', Tm, 'verdict', verdict, ...
    'notes', notes, 'swingPred_Hz', swingPred, 'swingMeas_Hz', swingMeas, 'raw', F.raw, 'binHz', df, 'Tsw', Tsw, ...
    'wide', wide, 'gate', gate);
A.best = F.best;
end

function [fMeas, snr, flo, chP, wid, tdet, ok, zs, sgr] = detectAll(rec, frel, fExpMid, halfW, rbw, obw, opt, wide, tauMid, Tsw, lines, lineW, tpl)
K = size(rec.P, 1); N = numel(frel);
fMeas = NaN(K, 1); snr = NaN(K, 1); flo = NaN(K, 1); chP = NaN(K, 1); wid = NaN(K, 1); tdet = tauMid; ok = false(K, 1);
zs = NaN(K, 1); sgr = NaN(K, 1);
for k = 1:K
    fa = rec.fCenter(k) + frel;
    mask = false(1, N);
    for L = lines(:).', mask = mask | abs(fa - L) <= lineW; end
    if wide
        d = detectSweepWide(double(rec.P(k, :)), fa, fExpMid(k), halfW, rbw, obw, opt.ZMin, mask, tpl);
    else
        d = detectSweep(double(rec.P(k, :)), fa, fExpMid(k), halfW, rbw, obw, opt.SNRMin, mask);
    end
    flo(k) = d.floor; snr(k) = d.snr; zs(k) = d.z; sgr(k) = d.sigRel;
    if d.ok
        fMeas(k) = d.f; chP(k) = d.chP; wid(k) = d.width; ok(k) = true;
        tdet(k) = tauMid(k) + ((d.j - 1)/(N - 1) - 0.5)*Tsw(k);
    end
end
end

function fl = fixedLine(fMeas, ok, tdet, eph, f0, rbw, df, lines)
% A detection cluster at one absolute frequency (+-1 RBW) lasting >= 30 s while the predicted
% Doppler moved by > 5 RBW cannot be the satellite.
fl = NaN; i = find(ok); if numel(i) < 5, return, end
w = max(rbw, 2*df); fm = fMeas(i);
edges = (min(fm) - w):w:(max(fm) + 2*w);
if numel(edges) < 3 || numel(edges) > 1e6, return, end
[c, ~, bin] = histcounts(fm, edges);
c2 = c + [c(2:end) 0];
[cm, b] = max(c2);
if cm < max(5, 0.15*numel(i)), return, end
sel = bin == b | bin == b + 1; tt = tdet(i(sel));
if max(tt) - min(tt) < 30, return, end
dd = interpEph(eph, tt, 'dfac'); swing = f0*(max(dd) - min(dd));
if ~(swing > 5*w), return, end
fl = median(fm(sel));
if any(abs(lines - fl) < w), fl = NaN; end
end

function d = detectSweep(y, fa, fExp, halfW, rbw, obw, snrMin, mask)
% Floor = MEAN noise power (median-clipped, so it is unbiased for power-averaging and sample
% detectors), strongest 3-bin-smoothed peak in the window, region, centroid / parabolic peak.
N = numel(y); p = 10.^(y/10);
if nargin > 7 && any(mask) && ~all(mask), p(mask) = median(p(~mask)); end
m0 = median(p);
for it = 1:4
    k = p < 4*m0; if nnz(k) < N/4, break, end
    m0 = median(p(k));
end
m = mean(p(p < 10*m0));
d = struct('ok', false, 'f', NaN, 'snr', NaN, 'floor', 10*log10(m), 'chP', NaN, 'width', NaN, 'j', NaN, 'z', NaN, 'sigRel', NaN);
ps = movmean(p, 3);
win = find(abs(fa - fExp) <= halfW);
if isempty(win), return, end
[pk, i] = max(ps(win)); j = win(i);
d.snr = 10*log10(pk/m); d.j = j;
if d.snr < snrMin, return, end
thr = m*10^(max(1.5, min(3, d.snr - 3))/10);
a = j; while a > 1 && ps(a - 1) > thr, a = a - 1; end
b = j; while b < N && ps(b + 1) > thr, b = b + 1; end
df = (fa(end) - fa(1))/(N - 1);
w = max(p(a:b) - m, 0);
if sum(w) <= 0, return, end
fc = sum(fa(a:b).*w)/sum(w);
if obw <= 2*rbw                                                % narrow line: parabola in dB (Gaussian RBW)
    [~, jj] = max(p(a:b)); jj = jj + a - 1;
    if jj > 1 && jj < N
        yl = 10*log10(p(jj - 1:jj + 1)); den = yl(1) - 2*yl(2) + yl(3);
        if den < 0
            del = 0.5*(yl(1) - yl(3))/den;
            if abs(del) <= 1, fc = fa(jj) + del*df; d.j = jj; end
        end
    end
end
inC = abs(fa - fc) <= max(obw/2 + rbw, 1.5*rbw); inC(a:b) = true;
noi = ~inC & p < 10*m0; if nnz(noi) >= N/4, m = mean(p(noi)); d.floor = 10*log10(m); d.snr = 10*log10(pk/m); end
d.chP = 10*log10(max(sum(max(p(inC) - m, 0))*df/(1.0645*rbw), 1e-30));
d.f = fc; d.width = (b - a + 1)*df;
d.ok = d.width <= 5*obw + 8*rbw + 3*df;                       % broadband junk is not the carrier
end

function d = detectSweepWide(y, fa, fExp, halfW, rbw, obw, zMin, mask, tpl)
% Wide signal (OBW > 3 RBW): box (matched) filter of OBW width on linear power. Noise statistics
% from the bins outside the signal: floor m (mean) and relative std r (robust, 1.4826 MAD/median);
% the box mean of B bins with a correlation length Lc = ENBW/bin has std m r sqrt(Lc/B).
N = numel(y); p = 10.^(y/10); dfb = (fa(end) - fa(1))/(N - 1);
if nargin > 7 && any(mask) && ~all(mask), p(mask) = median(p(~mask)); end
B = min(max(1, round(obw/dfb)), floor(N/3)); Lc = max(1, 1.0645*rbw/dfb);
d = struct('ok', false, 'f', NaN, 'snr', NaN, 'floor', NaN, 'chP', NaN, 'width', NaN, 'j', NaN, 'z', NaN, 'sigRel', NaN);
noi = p < 3*median(p);
for it = 1:2
    m = mean(p(noi)); r = 1.4826*mad(p(noi), 1)/median(p(noi));
    sg = m*r*sqrt(min(1, Lc/B));
    pw = movmean(p, B);
    win = find(abs(fa - fExp) <= halfW);
    if isempty(win) || ~isfinite(sg) || sg <= 0, d.floor = 10*log10(m); return, end
    [pk, i] = max(pw(win)); j = win(i);
    ex = abs(fa - fa(j)) <= obw/2 + 2*rbw;                      % exclude the candidate, re-estimate
    nn = noi & ~ex;
    if nnz(nn) < max(20, N/5), break, end
    noi = nn;
end
d.floor = 10*log10(m); d.z = (pk - m)/sg; d.snr = 10*log10(max(pk, m*1.0001)/m); d.sigRel = sg/m; d.j = j;
if d.z < zMin, return, end
h = m + 0.5*(pk - m);                                          % half-maximum region of the box output
a = j; while a > 1 && pw(a - 1) > h, a = a - 1; end
b = j; while b < N && pw(b + 1) > h, b = b + 1; end
w = pw(a:b) - m; fc = sum(fa(a:b).*w)/sum(w);
if nargin > 8 && ~isempty(tpl), fc = tplFit(p - m, fa, fc, tpl, obw, dfb); end   % least-squares template shift
[~, jc] = min(abs(fa - fc)); d.j = jc;
inC = abs(fa - fc) <= obw/2 + rbw;
d.chP = 10*log10(max(sum(p(inC) - m)*dfb/(1.0645*rbw), 1e-30));
d.f = fc; d.width = (b - a + 1)*dfb;
d.ok = d.width <= 2.5*obw + 4*rbw && d.width >= 0.3*obw;       % half-max width of box*signal ~ OBW
end

function fc = tplFit(e, fa, fc0, tpl, obw, dfb)
% Shift of the displayed signal template (ff_satsig spectrum through the RBW) that maximises the
% least-squares fit to the noise-subtracted trace (matched filter), parabolic refinement.
lagMax = obw/2; reg = abs(fa - fc0) <= obw + lagMax;
fc = fc0; if nnz(reg) < 5, return, end
x = fa(reg); x = x(:); ev = e(reg); ev = ev(:);
lags = -lagMax:dfb/2:lagMax;
D = interp1(tpl.u, tpl.D, x - fc0 - lags, 'linear', 0);
c = sum(ev.*D, 1); sd = sum(D.^2, 1); q = sign(c).*c.^2./max(sd, realmin);
[~, i] = max(q); del = 0;
if i > 1 && i < numel(q)
    den = q(i - 1) - 2*q(i) + q(i + 1); if den < 0, del = max(-1, min(1, 0.5*(q(i - 1) - q(i + 1))/den)); end
end
fc = fc0 + lags(i) + del*dfb/2;
end

function [rec, Tint, nSw, sd] = integrateSweeps(rec, tgt, plan, opt, wide)
% Average consecutive sweeps (same centre) in linear power over Tint seconds.
K = size(rec.P, 1); nSw = ones(K, 1); sd = seconds(rec.tEnd - rec.tStart);
Tint = opt.Integrate;
if isempty(Tint)
    if wide && ~isBursty(tgt.burst) && tgt.class ~= "CW"
        Tint = min(5, 0.05*tgt.obw/max(plan.maxRate, 1e-3));
    else
        Tint = 0;
    end
end
cyc = median(seconds(diff(rec.t)), 'omitnan'); if ~isfinite(cyc), cyc = Inf; end
if Tint <= 0 || K < 4 || Tint < 1.5*cyc, Tint = 0; return, end
ts = seconds(rec.t - rec.t(1)); g = zeros(K, 1); n = 1; a = 1;
for k = 1:K
    if k > a && (ts(k) - ts(a) >= Tint || rec.fCenter(k) ~= rec.fCenter(a)), n = n + 1; a = k; end
    g(k) = n;
end
L = 10.^(double(rec.P)/10);
P = zeros(n, size(L, 2));
for c = 1:size(L, 2), P(:, c) = accumarray(g, L(:, c), [n 1], @mean); end
nSw = accumarray(g, 1);
first = accumarray(g, (1:K)', [n 1], @min); last = accumarray(g, (1:K)', [n 1], @max);
tm = rec.t(1) + seconds(accumarray(g, ts, [n 1], @mean));
r2 = rec;
r2.P = single(10*log10(P)); r2.t = tm; r2.tStart = rec.tStart(first); r2.tEnd = rec.tEnd(last);
r2.fCenter = rec.fCenter(first);
sd = accumarray(g, sd, [n 1], @mean);
rec = r2;
end

function F = fitDoppler(tq, fm, eph, f0, fOff, opt, binHz)
% lsqcurvefit of f(t) = (f0 + df + a (t - tm)) (1 + d(t + dt; x)), three models, robust re-fit, BIC.
tq = tq(:); fm = fm(:); n = numel(tq);
F = struct('table', table(), 'best', [], 'keep', [], 'rmsFixed_Hz', NaN, 'raw', struct('mean_Hz', NaN, 'rms_Hz', NaN), 'models', []);
if n >= 1
    r0 = fm - (f0 + fOff)*(1 + dfacX(eph, tq, 0));
    F.raw = struct('mean_Hz', mean(r0), 'rms_Hz', sqrt(mean(r0.^2)), 'std_Hz', std(r0));
end
if n < 5, return, end
tm = mean(tq);
names = ["dt+df", "dt+df+drift", "dt+df+xtrack"]; masks = logical([1 1 0 0; 1 1 1 0; 1 1 0 1]);
dtB = opt.MaxTimeOffset; dfB = max(20e3, 50e-6*f0);
lbA = [-dtB, -dfB/1e3, -20, -150]; ubA = -lbA;
% coarse start: grid over dt, df = median residual
best = Inf; x0 = [0 0 0 0];
for dt = -dtB:0.5:dtB
    r = fm - f0*(1 + dfacX(eph, tq + dt, 0)); dfv = median(r); s = sqrt(mean((r - dfv).^2));
    if s < best, best = s; x0 = [dt, dfv/1e3, 0, 0]; end
end
% fitted as f - f0 (values ~kHz) so finite-difference Jacobians keep full precision
model = @(q, t) 1e3*q(2) + q(3)*(t - tm) + (f0 + 1e3*q(2) + q(3)*(t - tm)).*dfacX(eph, t + q(1), q(4));
fm0 = fm - f0;
o = optimoptions('lsqcurvefit', 'Display', 'off', 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-9, ...
    'OptimalityTolerance', 1e-10, 'MaxIterations', 200);
keep = true(n, 1); rows = {};
models = struct('name', {}, 'q', {}, 'model', {}, 'res', {});
for m = 1:3
    fr = masks(m, :);
    if nnz(keep) < nnz(fr) + 3, continue, end
    expand = @(x) subsasgnFree(x0.*fr, fr, x);
    fun = @(x, t) model(expand(x), t);
    xs = x0(fr);
    for it = 1:(1 + 2*(m == 1))                                % outlier rejection on the base model only
        [x, resn, ~, ~, ~, ~, J] = lsqcurvefit(fun, xs, tq(keep), fm0(keep), lbA(fr), ubA(fr), o);
        res = fm0 - fun(x, tq);
        if m > 1, break, end
        sc = 1.4826*median(abs(res(keep) - median(res(keep))));
        nk = abs(res) <= max(4*sc, binHz);
        if isequal(nk, keep) || nnz(nk) < 5, break, end
        keep = nk; xs = x;
    end
    nk = nnz(keep); mf = nnz(fr);
    sig = NaN(1, mf);
    if nk > mf
        Jf = full(J); cv = (Jf.'*Jf)\eye(mf)*resn/(nk - mf); sig = sqrt(abs(diag(cv))).';
    end
    q = expand(x); qs = zeros(1, 4); qs(fr) = sig;
    rms = sqrt(resn/nk); bic = nk*log(max(resn/nk, 1e-12)) + mf*log(nk);
    rows(end+1, :) = {names(m), nk, q(1), qs(1), 1e3*q(2), 1e3*qs(2), 1e3*q(2)/f0*1e6, q(3), qs(3), q(4), qs(4), rms, bic}; %#ok<AGROW>
    models(end+1) = struct('name', names(m), 'q', q, 'model', @(t) f0 + model(q, t(:)), 'res', res); %#ok<AGROW>
    if m == 1, x0(fr) = x; end
end
if isempty(rows), return, end
T = cell2table(rows, 'VariableNames', ["Model", "N", "dt_s", "dt_sigma_s", "df_Hz", "df_sigma_Hz", "ppm", ...
    "drift_Hz_s", "drift_sigma", "xtrack_km", "xtrack_sigma_km", "rms_Hz", "BIC"]);
b = 1;                                                         % a richer model must win clearly (dBIC > 6)
for m = 2:height(T), if T.BIC(m) < T.BIC(b) - 6, b = m; end, end
T.Chosen = false(height(T), 1); T.Chosen(b) = true;
F.table = T; F.keep = keep; F.models = models;
F.best = table2struct(T(b, :)); F.best.model = models(b).model; F.best.q = models(b).q;
F.best.alongTrack_km = F.best.dt_s*mean(eph.speed)/1e3;
c = polyfit(tq(keep) - tm, fm(keep), 1); F.rmsFixed_Hz = sqrt(mean((fm(keep) - polyval(c, tq(keep) - tm)).^2));
end

function q = subsasgnFree(q, fr, x)
q(fr) = x;
end

%% ======================================================================================
%% report, figures, files
%% ======================================================================================
function L = buildReport(R, opt)
t = R.target; p = R.plan; s = R.pass; T = R.timing; Lk = R.link; tr = R.track;
L = strings(0, 1);
L(end+1) = sprintf("FF_SAT_TRACKER (%s): %s, NORAD %d, %s downlink %s, OBW %s", R.source, t.name, t.norad, t.id, fmtHz(t.f0), fmtHz(t.obw));
L(end+1) = sprintf("  modulation: %s", t.modulation);
L(end+1) = sprintf("TLE: epoch %s UTC, %.2f days before the capture; %s", string(R.tle.epoch, 'yyyy-MM-dd HH:mm:ss'), p.tleAgeDays, R.tle.source);
L(end+1) = sprintf("Pass: AOS %s, TCA %s (max el %.1f deg), LOS %s (local %s)", localStr(s.AOS, opt, 'HH:mm:ss'), ...
    localStr(s.TCA, opt, 'HH:mm:ss'), s.MaxEl, localStr(s.LOS, opt, 'HH:mm:ss'), localStr(s.TCA, opt, 'yyyy-MM-dd'));
L(end+1) = sprintf("Capture: %s mode, span %s, RBW %s, %d points, %s detector, preamp %d, step %g s, search +-%s", ...
    p.mode, fmtHz(p.span), fmtHz(p.rbw), p.points, p.detector, p.preamp, p.step, fmtHz(p.halfW));
L(end+1) = sprintf("  %d sweeps; cycle %.3f s mean (p95 %.3f), sweep+OPC %.3f s, trace read %.3f s, overhead %.3f s; %d retunes; max gap %.2f s", ...
    T.sweeps, T.cycle_mean_s, T.cycle_p95_s, T.sweep_mean_s, T.read_mean_s, T.overhead_mean_s, T.retunes, T.maxGap_s);
L(end+1) = sprintf("  Doppler max +-%s, max rate %.1f Hz/s (moves %.0f Hz per sweep); dopplershift vs states check %.3f Hz", ...
    fmtHz(p.maxDoppler), p.maxRate, T.smear_Hz, p.dsCheck_Hz);
nd = nnz(tr.detected); nu = nnz(tr.used);
gate = sprintf("%g dB SNR", opt.SNRMin); if any(isfinite(tr.Z)), gate = sprintf("%g sigma (OBW box filter)", opt.ZMin); end
if T.integrate_s > 0
    L(end+1) = sprintf("Detections: %d of %d blocks (%.1f s, %d sweeps) above %s, %d used after outlier rejection", nd, height(tr), T.integrate_s, T.sweeps, gate, nu);
else
    L(end+1) = sprintf("Detections: %d of %d sweeps above %s, %d used after outlier rejection", nd, height(tr), gate, nu);
end
raw = R.fits;
if ~isempty(R.best)
    L(end+1) = sprintf("TLE prediction alone (with 'FreqOffset' %g Hz): residual mean %+.0f Hz, RMS %.0f Hz", p.freqOffset, R.best.raw.mean_Hz, R.best.raw.rms_Hz);
    L(end+1) = "Fits  f(t) = (f0 + df + a(t - tm)) (1 + d(t + dt; x)):";
    for k = 1:height(raw)
        L(end+1) = sprintf("  %-13s N %4d  dt %+7.2f +- %5.2f s  df %+8.1f +- %6.1f Hz (%+6.2f ppm)  a %+6.2f Hz/s  x %+6.2f km  RMS %7.1f Hz  BIC %8.1f%s", ...
            raw.Model(k), raw.N(k), raw.dt_s(k), raw.dt_sigma_s(k), raw.df_Hz(k), raw.df_sigma_Hz(k), raw.ppm(k), ...
            raw.drift_Hz_s(k), raw.xtrack_km(k), raw.rms_Hz(k), raw.BIC(k), ternary(raw.Chosen(k), "  <- chosen", "")); %#ok<AGROW>
    end
    b = R.best;
    L(end+1) = sprintf("  -> transmitter frequency %s (%+.0f Hz = %+.2f ppm vs %s; includes the FieldFox time-base error)", ...
        fmtHz(t.f0 + b.df_Hz), b.df_Hz, b.ppm, fmtHz(t.f0));
    L(end+1) = sprintf("  -> time offset %+.2f s = PC clock error + TLE along-track error (%+.1f km along track if the clock is right)", b.dt_s, b.alongTrack_km);
    L(end+1) = sprintf("  -> fixed-carrier (no Doppler) fit RMS %s; measured Doppler swing %s vs predicted %s", ...
        fmtHz(b.rmsFixed_Hz), fmtHz(R.best.swingMeas_Hz), fmtHz(R.best.swingPred_Hz));
end
L(end+1) = "Verdict: " + R.verdict;
L(end+1) = sprintf("Link: EIRP %.1f dBm (%s); satellite antenna %s; ground: %s", Lk.EIRP_dBm, Lk.EIRPsource, Lk.tx, Lk.rx);
L(end+1) = sprintf("  noise at the FieldFox input %.1f dBm/Hz (DANL %.0f, external %.1f; %s)", Lk.N0in_pred, Lk.noise.danl, Lk.noise.Next, Lk.noise.text);
L(end+1) = sprintf("  noise floor in RBW: measured %.1f dBm vs predicted %.1f dBm (excess %+.1f dB%s)", Lk.floorMeas_dBm, ...
    Lk.floorPred_dBm, Lk.noiseExcess_dB, ternary(Lk.noiseExcess_dB > 3, ": local noise/EMI or LNA/antenna-noise not modelled", ""));
if Lk.n >= 3
    L(end+1) = sprintf("  C/N0 measured - predicted: median %+.1f dB (IQR %.1f dB, %d %s) -> EIRP ~ %.1f dBm if the receive chain is right", ...
        Lk.offset_dB, Lk.offsetIQR_dB, Lk.n, Lk.offsetBasis, Lk.EIRP_dBm + Lk.offset_dB);
end
v = isfinite(Lk.detFrac);
if any(v)
    L(end+1) = "  elevation  detect%  C/N0 meas  C/N0 pred (dB-Hz)";
    for b2 = find(v)'
        L(end+1) = sprintf("   %2d-%2d deg  %5.0f   %8.1f   %8.1f", Lk.elBins(b2), Lk.elBins(b2) + 10, 100*Lk.detFrac(b2), Lk.cn0MeasBin(b2), Lk.cn0PredBin(b2)); %#ok<AGROW>
    end
end
if ~isempty(R.sim)
    tr0 = R.sim.truth;
    L(end+1) = sprintf("Simulation truth: dt %+.2f s, df %+.0f Hz (%+.2f ppm), x-track %.1f km, EIRP error %+.1f dB, noise excess %+.1f dB; signal %s", ...
        tr0.dt_s, tr0.df_Hz, tr0.ppm, tr0.xtrack_km, tr0.eirpError_dB, tr0.noiseExcess_dB, R.sim.shape);
    if isfield(R.sim, 'err')
        L(end+1) = sprintf("  recovered: dt error %+.2f s (%.1f sigma), df error %+.1f Hz (%.1f sigma), EIRP offset %+.1f dB", ...
            R.sim.err.dt_s, abs(R.sim.err.dt_sigmas), R.sim.err.df_Hz, abs(R.sim.err.df_sigmas), R.sim.recovered.eirp_dB);
    end
end
for n = R.notes(:)', L(end+1) = "Note: " + n; end %#ok<AGROW>
L = L(:);
end

function figs = makeFigures(R, ~, eph, opt)
figs = gobjects(0);
tr = R.track; rec = R.rec; s = rec.settings; f0 = R.target.f0;
tc = R.pass.TCA; tmin = minutes(rec.t - tc);
N = numel(rec.f); frel = rec.f - s.frameCenter; K = size(rec.P, 1);
% raw waterfall: offset from f0 (absolute frame, every sweep at its own centre)
fa0 = rec.fCenter - f0;                                        % K x 1 centre offsets
lo = min(fa0) + frel(1); hi = max(fa0) + frel(end);
step = max((frel(end) - frel(1))/(N - 1), (hi - lo)/1500);
gx = lo:step:hi;
gx2 = frel(1):max((frel(end) - frel(1))/(N - 1), (frel(end) - frel(1))/1500):frel(end);
fmod = NaN(K, 1);
if ~isempty(R.best), fmod = R.best.model(seconds(rec.t - eph.t0)); end
fexpMid = (f0 + s.freqOffset)*(1 + interpEph(eph, seconds(rec.t - eph.t0), 'dfac'));
gD = max(1, ceil(K/1500)); nr = ceil(K/gD);                   % display rows (linear-power average)
W1 = zeros(nr, numel(gx), 'single'); C1 = W1; W2 = zeros(nr, numel(gx2), 'single'); C2 = W2;
for k = 1:K
    r = ceil(k/gD); y = 10.^(double(rec.P(k, :))/10);
    v = interp1(fa0(k) + frel, y, gx, 'linear', NaN); m = isfinite(v);
    W1(r, m) = W1(r, m) + v(m); C1(r, m) = C1(r, m) + 1;
    ref = fmod(k); if ~isfinite(ref), ref = fexpMid(k); end
    v = interp1(rec.fCenter(k) + frel - ref, y, gx2, 'linear', NaN); m = isfinite(v);
    W2(r, m) = W2(r, m) + v(m); C2(r, m) = C2(r, m) + 1;
end
W1 = 10*log10(W1./C1); W1(C1 == 0) = NaN; W2 = 10*log10(W2./C2); W2(C2 == 0) = NaN;
tmin = accumarray(ceil((1:K)'/gD), tmin(:), [nr 1], @mean);
cl = [prctile(W1(:), 5), max(W1(:))]; if ~all(isfinite(cl)) || cl(2) <= cl(1), cl = [-160 -60]; end
fig = figure('Visible', opt.Visible, 'Name', 'ff_sat_tracker Doppler', 'Color', 'w', 'Position', [60 60 1250 820]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s  %s  pass %s (max el %.0f deg)  [%s]', R.target.name, fmtHz(f0), localStr(tc, opt, 'yyyy-MM-dd HH:mm z'), R.pass.MaxEl, R.source), 'Interpreter', 'none');
[fsc, fun] = fscale(max(abs([lo hi])));
ax = nexttile(tl);
im = imagesc(ax, gx/fsc, tmin, W1); set(im, 'AlphaData', ~isnan(W1)); set(ax, 'YDir', 'normal'); clim(ax, cl); colormap(ax, turbo);
hold(ax, 'on');
te = minutes(eph.t - tc); inr = te >= min(tmin) - 0.5 & te <= max(tmin) + 0.5;
plot(ax, ((f0 + s.freqOffset)*(1 + eph.dfac(inr)) - f0)/fsc, te(inr), 'w--', 'LineWidth', 1);
u = tr.used;
plot(ax, (tr.fMeas_Hz(u) - f0)/fsc, minutes(tr.t(u) - tc), 'm.', 'MarkerSize', 6);
if ~isempty(R.best), plot(ax, (R.best.model(eph.tau(inr)) - f0)/fsc, te(inr), 'k-', 'LineWidth', 0.8); end
xlabel(ax, "Offset from f0 (" + fun + ")"); ylabel(ax, 'Minutes from TCA'); title(ax, 'Waterfall: TLE Doppler (dashed), detections (dots), fit (line)');
cb = colorbar(ax); cb.Label.String = 'dBm';
ax = nexttile(tl);
[fs2, fu2] = fscale(max(abs(gx2)));
im = imagesc(ax, gx2/fs2, tmin, W2); set(im, 'AlphaData', ~isnan(W2)); set(ax, 'YDir', 'normal'); clim(ax, cl); colormap(ax, turbo);
xlabel(ax, "Offset from the fitted Doppler curve (" + fu2 + ")"); ylabel(ax, 'Minutes from TCA');
title(ax, ternary(isempty(R.best), 'Doppler-corrected waterfall (TLE prediction)', 'Doppler-corrected waterfall (fitted model)'));
colorbar(ax);
ax = nexttile(tl);
plot(ax, nanIfEmpty(minutes(tr.t(u) - tc)), nanIfEmpty(tr.fMeas_Hz(u) - tr.fPred_Hz(u)), '.', 'Color', [0.85 0.33 0.1]); hold(ax, 'on');
if ~isempty(R.best), plot(ax, minutes(tr.t(u) - tc), tr.fMeas_Hz(u) - tr.fFit_Hz(u), '.', 'Color', [0 0.45 0.74]); end
yline(ax, 0, 'k:'); grid(ax, 'on');
xlabel(ax, 'Minutes from TCA'); ylabel(ax, 'Residual (Hz)');
if ~isempty(R.best)
    legend(ax, sprintf('measured - TLE (RMS %.0f Hz)', R.best.raw.rms_Hz), sprintf('measured - fit %s (RMS %.1f Hz)', R.best.Model, R.best.rms_Hz), 'Location', 'best');
else
    legend(ax, 'measured - TLE', 'Location', 'best');
end
title(ax, 'Doppler residuals');
ax = nexttile(tl);
yyaxis(ax, 'left');
h1 = plot(ax, te(inr), eph.fD(inr)/1e3, 'k-'); hold(ax, 'on');
h2 = plot(ax, nanIfEmpty(minutes(tr.t(u) - tc)), nanIfEmpty((tr.fMeas_Hz(u) - f0)/1e3), '.', 'Color', [0.85 0.33 0.1]);
ylabel(ax, 'Doppler (kHz)');
yyaxis(ax, 'right'); h3 = plot(ax, te(inr), eph.el(inr), '-', 'Color', [0.5 0.5 0.5]); ylabel(ax, 'Elevation (deg)');
xlabel(ax, 'Minutes from TCA'); grid(ax, 'on'); legend(ax, [h1 h2 h3], 'TLE prediction', 'measured', 'elevation', 'Location', 'best');
title(ax, 'Measured vs predicted Doppler');
figs(end+1) = fig;
% link / noise / timing
fig = figure('Visible', opt.Visible, 'Name', 'ff_sat_tracker link', 'Color', 'w', 'Position', [80 80 1250 420]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s: link budget check and timing', R.target.name));
ax = nexttile(tl);
ev = eph.el >= 0 & inr;
[es, ord] = sort(eph.el(ev)); cp = eph.CN0(ev); plot(ax, es, cp(ord), 'k-'); hold(ax, 'on');
plot(ax, nanIfEmpty(tr.el_deg(u)), nanIfEmpty(tr.CN0meas_dBHz(u)), '.', 'Color', [0 0.45 0.74]);
yline(ax, R.link.cn0Threshold, 'r:', 'detection threshold (approx.)');
xlabel(ax, 'Elevation (deg)'); ylabel(ax, 'C/N0 (dB-Hz)'); grid(ax, 'on');
legend(ax, 'predicted (link budget)', 'measured', 'Location', 'best'); title(ax, 'C/N0 vs elevation');
ax = nexttile(tl);
ttr = minutes(tr.t - tc);
plot(ax, ttr, tr.floor_dBm, '.', 'Color', [0.5 0.5 0.5]); hold(ax, 'on');
yline(ax, R.link.floorPred_dBm, 'k--', 'predicted floor');
yyaxis(ax, 'right');
if any(isfinite(tr.Z))
    plot(ax, ttr, tr.Z, '.', 'Color', [0.85 0.33 0.1]); ylabel(ax, 'Box-filter detection statistic (sigma)');
else
    plot(ax, ttr, tr.SNR_dB, '.', 'Color', [0.85 0.33 0.1]); ylabel(ax, 'Peak SNR in window (dB)');
end
yyaxis(ax, 'left'); ylabel(ax, 'Noise floor in RBW (dBm)'); xlabel(ax, 'Minutes from TCA'); grid(ax, 'on');
title(ax, 'Noise floor (EMI check) and SNR');
ax = nexttile(tl);
c = rec.timing.cycle_s; c = c(isfinite(c));
if ~isempty(c)
    histogram(ax, 1e3*c, 30); hold(ax, 'on');
    xline(ax, 1e3*mean(rec.timing.sweep_s, 'omitnan'), 'r--', 'sweep + *OPC?');
end
xlabel(ax, 'Cycle time per sweep (ms)'); ylabel(ax, 'Sweeps'); grid(ax, 'on');
title(ax, sprintf('Timing: %d sweeps, %d retunes', R.timing.sweeps, R.timing.retunes));
figs(end+1) = fig;
end

function files = saveOutputs(R, ~, opt)
files = strings(0, 1);
live = startsWith(R.source, "live");
if strlength(opt.OutDir) == 0 && ~live, return, end
d = opt.OutDir; if strlength(d) == 0, d = string(pwd); end
stamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
f = fullfile(d, "tracker_report.txt"); writelines(R.report, f); files(end+1) = f;
for k = 1:numel(R.figs)
    if isvalid(R.figs(k))
        nm = ["tracker_doppler.png", "tracker_link.png"];
        f = fullfile(d, nm(min(k, 2))); exportgraphics(R.figs(k), f, 'Resolution', 130); files(end+1) = f; %#ok<AGROW>
    end
end
rec = R.rec; rec.analysis = struct('track', R.track, 'fits', R.fits, 'verdict', R.verdict, 'link', R.link, 'timing', R.timing);
if live && isfield(R, 'notes') && any(startsWith(R.notes, "Recording saved: "))
    f = extractAfter(R.notes(find(startsWith(R.notes, "Recording saved: "), 1)), "Recording saved: ");
else
    f = fullfile(d, "ff_sat_tracker_" + lower(extractBefore(R.source + " ", " ")) + "_" + stamp + ".mat");
end
saveRec(f, rec); files(end+1) = f;
writetable(R.track, fullfile(d, "tracker_track.csv")); files(end+1) = fullfile(d, "tracker_track.csv");
if opt.PDF
    try
        files(end+1) = makePDF(R, files, d);
    catch err
        warning('ff_sat_tracker:pdf', 'PDF report failed: %s', err.message);
    end
end
files = files(:);
end

function f = makePDF(R, files, d)
import mlreportgen.report.*
import mlreportgen.dom.*
f = fullfile(d, "tracker_report.pdf");
rpt = Report(char(f), 'pdf');
add(rpt, TitlePage('Title', 'FieldFox satellite Doppler tracking', 'Subtitle', sprintf('%s - %s', R.target.name, R.source), ...
    'Author', 'ff_sat_tracker'));
ch = Chapter('Summary');
add(ch, Preformatted(char(strjoin(R.report, newline))));
add(rpt, ch);
ch = Chapter('Figures');
for p = reshape(files(endsWith(files, ".png")), 1, [])
    im = Image(char(p)); im.Width = '6.5in'; im.Height = []; add(ch, im);
end
add(rpt, ch);
if ~isempty(R.fits)
    ch = Chapter('Doppler fits');
    T = R.fits(:, ["Model", "N", "dt_s", "dt_sigma_s", "df_Hz", "df_sigma_Hz", "ppm", "rms_Hz", "BIC", "Chosen"]);
    add(ch, MATLABTable(T)); add(rpt, ch);
end
close(rpt);
end

%% ======================================================================================
%% small helpers
%% ======================================================================================
function v = nanIfEmpty(v)
if isempty(v), v = NaN; end
end

function b = isBursty(txt)
% Catalog Burst text -> bursty (packets/TDMA/periodic beacons) or continuous.
t = lower(string(txt));
b = strlength(t) > 0 && ~contains(t, ["continuous", "constant"]) && ...
    contains(t, ["burst", "tdma", "packet", "every", "periodic", "beacon", "event", "test"]);
end

function v = niceUp(x)
e = 10^floor(log10(x)); m = [1 2 5 10]; v = e*m(find(m*e >= x, 1));
end

function r = rbwStep(x)
steps = [10 30 100 300 1e3 3e3 10e3 30e3 100e3 300e3 1e6 3e6 5e6];
r = steps(find(steps >= x, 1)); if isempty(r), r = 5e6; end
end

function [sc, unit] = fscale(x)
if x >= 2e6, sc = 1e6; unit = "MHz"; elseif x >= 2e3, sc = 1e3; unit = "kHz"; else, sc = 1; unit = "Hz"; end
end

function s = fmtHz(f)
if ~isfinite(f), s = "n/a"; return, end
a = abs(f);
if a >= 1e9
    s = sprintf('%.6f GHz', f/1e9);
elseif a >= 1e6
    s = sprintf('%.4f MHz', f/1e6);
elseif a >= 1e3
    s = sprintf('%.2f kHz', f/1e3);
else
    s = sprintf('%.1f Hz', f);
end
s = string(s);
end

function o = ternary(c, a, b)
if c, o = a; else, o = b; end
end
