function R = ff_emi_study(ip, varargin)
%FF_EMI_STUDY  Solar-inverter EMI study with a FieldFox: OFF-ON-OFF scans, statistics, report.
%
%   R = ff_emi_study                               measure (FieldFox at 10.90.36.200)
%   R = ff_emi_study(ip, Name, Value, ...)
%   R = ff_emi_study([], 'Data', "demo")           synthetic demo, no instrument needed
%   R = ff_emi_study([], 'Data', "...\emi_study.mat")             re-analyse a saved study
%   R = ff_emi_study([], 'Data', ["off.mat" "on.mat" "off2.mat"]) recordings (fieldfox_live)
%   R = ff_emi_study([], 'Data', "...\emi_study_20261009_101500")  a study folder
%   R = ff_emi_study([], 'DryRun', true)           only the scan plan (R.bands), to edit
%
% Example (edit the plan, measure, no dialogs):
%   P = ff_emi_study([], 'DryRun', true);  B = P.bands([1 2 5], :);  B.Points(1) = 4001;
%   R = ff_emi_study("10.90.36.200", 'Bands', B, 'Sweeps', 10, 'Notes', "roof, 3 m, biconical");
%
% Output R: bands, protocol, summary (one row per band and detector), emissions
% (one row per flagged emission), comb, zeroSpan, limits, range, conclusions,
% files, outDir and study (all raw sweeps).
%
% Protocol (A-B-A): every band is scanned with the inverter OFF, then ON, then OFF
% again. The second OFF shows whether the ambient (radio stations, phones, other
% equipment) changed during the test. Each band is swept 'Sweeps' times with the
% Peak detector (DET POS) and then with the Average detector (DET AVER): firmware
% A.12.45 has one detector for all traces, so they are separate sweep sets.
%
% Analysis, for every frequency point of every band and detector:
%   - Wilcoxon rank-sum test (ranksum) OFF vs ON, Benjamini-Hochberg false
%     discovery rate over the band (mafdr) -> q value
%   - effect = ON - mean(OFF1, OFF2) in dB, sweep by sweep, with a bootstrap
%     confidence interval (bootci)
%   - flagged ("new emission") when q < Alpha and the CI lower bound > MinEffect dB
%   - harmonic comb of the switching frequency: cepstrum (periodogram of the ON-OFF
%     dB difference), then a least-squares fit of the harmonic lines (fitlm)
%   - flagged peaks are labelled with ff_signal_id (if installed) and a band table,
%     to tell licensed services (FM, LTE, ...) from inverter emissions
%   - zero span (level vs time) on the 'ZeroSpanK' strongest new emissions, inverter
%     ON and OFF: statelevels, pulsewidth, pulseperiod, dutycycle and the envelope
%     spectrum (pspectrum) -> burst rate and duty cycle
%   - with a transducer: levels in dBuV or dBuV/m and an INDICATIVE comparison with
%     CISPR 11/32 class B limits
%   - interference range: distance at which each emission falls below the noise of a
%     victim receiver, propagationModel("freespace") and ("close-in") + pathloss
% Results go to 'OutDir': emi_report.pdf, figures (PNG), emissions.csv, emi_study.mat
% (everything; re-use it with 'Data') and rec\*.mat recordings (f, P, t, settings,
% idn - the fieldfox_live format, one file per state, band and detector).
%
% Default scan plan (R.bands; edit the table and pass it back with 'Bands'):
%   Conducted  150 kHz - 30 MHz   RBW 9 kHz    switching noise of the inverter
%   VHF        30 - 300 MHz       RBW 120 kHz
%   UHF        300 MHz - 1 GHz    RBW 120 kHz
%   L/S        1 - 3 GHz          RBW 1 MHz
%   GNSS L1    1565 - 1585 MHz    RBW 100 kHz
%   Points = 2 per RBW (401 to 10001).
%
% Name-Value options (defaults):
%  Measurement
%   'Bands'          table   Name, Start_Hz, Stop_Hz, RBW_Hz [, Points]
%   'Sweeps'         20      sweeps per band, detector and inverter state (at least 3;
%                            fewer than ~20 cannot flag a single narrowband line)
%   'Detectors'      ["POS" "AVER"]  Peak and Average detector sweep sets
%   'Preamp'         false   preamplifier (option 235): lower noise floor, overloads sooner
%   'RefLevel'       -20     dBm (the auto attenuation follows it)
%   'Atten'          "auto"  or a value in dB
%   'Protocol'       "ABA"   "ABA" = OFF, ON, OFF;  "AB" = OFF, ON
%   'Prompt'         true    ask before each step; false = run the steps back to back
%   'SwitchFcn'      []      function called with "OFF" or "ON" before each step
%                            (e.g. to drive a relay or smart plug)
%   'SettleTime'     0       s to wait after each switch
%   'ZeroSpanK'      3       strongest new emissions examined in zero span (0 = none),
%                            the strongest of each band first
%   'ZeroSpanFreqs'  []      Hz, extra frequencies always examined in zero span
%                            (e.g. a harmonic you already know); must lie in a band
%   'ZeroSpanTime'   0.1     s per zero-span sweep (5 mains cycles at 50 Hz)
%   'ZeroSpanPoints' 4001
%   'ZeroSpanSweeps' 5
%   'ZeroSpanRBW'    []      Hz; [] = the RBW of the band
%   'Port'           5025
%  Analysis
%   'Test'           "ranksum"  or "ttest2" (Welch t-test on the dB values)
%   'Alpha'          0.05    false discovery rate; intervals are 1-Alpha (95 %)
%   'MinEffect'      3       dB, smallest rise that counts as an emission
%   'Boot'           1000    bootstrap resamples
%   'SignalID'       "auto"  ff_signal_id method for the labels: "auto" (= "ensemble" when
%                            its trained models are cached, else "template"), "template",
%                            "ml", "cnn", "ensemble", "none" (detection only) or "off"
%   'Transducer'     ""      CSV file (or N x 2 matrix): Freq_Hz, Factor_dB. Antenna
%                            factor incl. cable loss (dB/m), or LISN + limiter + cable (dB).
%                            Level = dBm + 107 + factor
%   'TransducerType' "antenna"  "antenna" -> dBuV/m, "lisn" -> dBuV
%   'LimitStandard'  "CISPR32"  or "CISPR11" (class B, indicative only)
%   'MeasDistance'   3       m from the antenna to the inverter
%   'RxGain'         0       dBi, measuring antenna (used when there is no transducer)
%   'VictimNF'       5       dB, noise figure of a victim receiver
%   'INR'            -6      dB, allowed interference-to-noise ratio
%   'VictimGain'     0       dBi
%  Output
%   'Report'         "pdf"   "pdf", "docx" or "none"
%   'Notes'          ""      text for the report (site, inverter model, antenna, ...)
%   'OutDir'         ""      default emi_study_<time> in the current folder
%   'Visible'        'on'    'off' = no windows (figures are still saved)
%   'Data'           ""      offline input, see above
%   'DryRun'         false
%
% Safe to stop: the instrument settings are saved first and restored at the end,
% also after an error or Ctrl+C. Every sweep set is saved at once (rec\*.mat) and
% the raw study goes to emi_study.mat before the analysis starts. If the
% connection fails half-way, what was measured (at least OFF1 and ON) is still
% analysed. The FieldFox only listens (no RF output).
%
% Read this before trusting the numbers:
%   - INDICATIVE PRE-SCAN ONLY. There is no quasi-peak or CISPR-average detector
%     (option 361 is not installed), the RBW is a 3 dB filter (CISPR uses 6 dB) and a
%     roof or garden is not a test site. Peak vs a QP limit and RMS average vs an AV
%     limit read high, so a "pass" is meaningful but a "fail" needs a real EMC lab.
%   - Below ~30 MHz an antenna near the inverter is in the near field (closer than
%     wavelength/2pi: 1.6 m at 30 MHz, 320 m at 150 kHz) and the noise travels mostly
%     on the cables: use a LISN there, and do not trust the range estimate.
%   - The range estimate treats the emission as a small source in free space or on
%     a typical urban path (close-in model, n = 2.9); real cables radiate differently.
%     It is only made when the antenna was in the far field (MeasDistance >
%     wavelength/2pi): with the default 3 m, from about 16 MHz up.
%   - Zero span: the envelope follows changes up to about RBW/2 (4.5 kHz at 9 kHz
%     RBW), so it shows mains-related bursts (50/100 Hz), not the switching itself;
%     the switching frequency comes from the harmonic comb. The zero-span sweep-time
%     command (SWE:TIME) is not yet verified on firmware A.12.45: the table column
%     TimeFrom says where the time axis came from.
%   - The comb can only be resolved when the switching frequency is above about
%     2.5 frequency steps and the RBW (11 kHz with the default plan); for a slower
%     inverter use a narrower conducted band with a smaller RBW.
%   - A significant change is not proof that the inverter caused it: look at the
%     OFF1/OFF2 drift, the signal labels, and repeat at another time of day.
%   - Safety: switch the inverter with its own AC breaker and DC isolator, never by
%     unplugging DC connectors. Protect the FieldFox input (max +27 dBm, 50 V DC)
%     with a transient limiter when you use a LISN.

if nargin < 1 || isempty(ip) || all(strlength(string(ip)) == 0), ip = "10.90.36.200"; end
opt = struct('Bands', defaultBands(), 'Sweeps', 20, 'Detectors', ["POS" "AVER"], 'Preamp', false, ...
    'RefLevel', -20, 'Atten', "auto", 'Protocol', "ABA", 'Prompt', true, 'SwitchFcn', [], ...
    'SettleTime', 0, 'ZeroSpanK', 3, 'ZeroSpanFreqs', [], 'ZeroSpanTime', 0.1, 'ZeroSpanPoints', 4001, ...
    'ZeroSpanSweeps', 5, 'ZeroSpanRBW', [], 'Port', 5025, 'Test', "ranksum", 'Alpha', 0.05, ...
    'MinEffect', 3, 'Boot', 1000, 'SignalID', "auto", 'Transducer', "", 'TransducerType', "antenna", ...
    'LimitStandard', "CISPR32", 'MeasDistance', 3, 'RxGain', 0, 'VictimNF', 5, 'INR', -6, ...
    'VictimGain', 0, 'Report', "pdf", 'Notes', "", 'OutDir', "", 'Visible', 'on', 'Data', "", ...
    'DryRun', false);
if mod(numel(varargin), 2), error('Options must be Name, Value pairs.'); end
names = fieldnames(opt);
for k = 1:2:numel(varargin)
    j = find(strcmpi(names, varargin{k}), 1);
    if isempty(j), error('Unknown option "%s".', string(varargin{k})); end
    opt.(names{j}) = varargin{k + 1};
end
opt.Visible = char(matlab.lang.OnOffSwitchState(opt.Visible));
opt.Prompt = logical(opt.Prompt);
opt.Detectors = upper(string(opt.Detectors));
opt.Protocol = upper(string(opt.Protocol));
opt.Test = lower(string(opt.Test));
opt.Report = lower(string(opt.Report));
opt.TransducerType = lower(string(opt.TransducerType));
opt.LimitStandard = upper(string(opt.LimitStandard));
opt.SignalID = lower(string(opt.SignalID));
if ~ismember(opt.Protocol, ["ABA" "AB"]), error('Protocol must be "ABA" or "AB".'); end
if ~ismember(opt.Test, ["ranksum" "ttest2"]), error('Test must be "ranksum" or "ttest2".'); end
if ~ismember(opt.Report, ["pdf" "docx" "none"]), error('Report must be "pdf", "docx" or "none".'); end
if opt.Sweeps < 3, error('Sweeps must be at least 3 (5 or more for useful statistics).'); end
bands = checkBands(opt.Bands);

if opt.DryRun
    R.bands = bands;
    R.states = protocolStates(opt.Protocol);
    R.sweepsTotal = height(bands) * numel(opt.Detectors) * opt.Sweeps * numel(R.states);
    disp(bands);
    fprintf('%d sweeps in total (%s), plus zero span on up to %d emissions.\n', ...
        R.sweepsTotal, strjoin(R.states, ' -> '), opt.ZeroSpanK);
    return
end

if strlength(string(opt.OutDir)) == 0, opt.OutDir = fullfile(pwd, "emi_study_" + stamp()); end
opt.OutDir = string(opt.OutDir);
if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end

%% 1. Data: measure, or load / synthesize
matFile = fullfile(opt.OutDir, "emi_study.mat");
if any(strlength(string(opt.Data)) > 0)
    S = loadData(string(opt.Data), opt, bands);
    if S.source == "demo", S = saveRecordings(S, opt.OutDir); end
else
    S = measure(ip, opt, bands);
    study = S;
    save(matFile, 'study');               % raw data first: re-analyse with 'Data' if anything below fails
end
if isempty(S.seg) || ~any([S.seg.state] == "ON") || ~any([S.seg.state] == "OFF1")
    warning('ff_emi_study:incomplete', 'No complete OFF/ON pair was measured: nothing to analyse.');
    R.study = S; R.outDir = opt.OutDir; R.files = string(matFile);
    return
end

%% 2. Analysis
fprintf('Analysing (%s test, FDR %.2g, bootstrap %d) ...\n', opt.Test, opt.Alpha, opt.Boot);
tr = loadTransducer(opt.Transducer);
res = analyzeStudy(S, opt, tr);
[res, E] = findEmissions(res, opt);
[res, E, limT] = applyLimits(res, E, opt, tr);
[zsT, Z] = zeroSpanMetrics(S, res);
[rangeT, curves] = rangeEstimate(E, opt);
summary = summaryTable(res);
combT = combTable(res);
conc = conclusions(S, res, E, combT, zsT, limT, rangeT, opt, tr);

%% 3. Files, figures and report
study = S; results = res; emissions = E; comb = combT; zeroSpan = zsT; limits = limT;
range = rangeT; conclusionsText = conc;
save(matFile, 'study', 'results', 'summary', 'emissions', 'comb', 'zeroSpan', 'limits', 'range', 'conclusionsText');
csvFile = fullfile(opt.OutDir, "emissions.csv");
writetable(E, csvFile);
files = makeFigures(S, res, Z, curves, opt);
files = [matFile; csvFile; files; [S.seg.file].'];
files = files(strlength(files) > 0);

R.study = S; R.results = res; R.bands = S.bands; R.protocol = S.protocol; R.summary = summary;
R.emissions = E; R.comb = combT; R.zeroSpan = zsT; R.limits = limT; R.range = rangeT;
R.conclusions = conc; R.outDir = opt.OutDir; R.files = files;
if opt.Report ~= "none"
    try
        rep = makeReport(R, opt, files);
        R.files = [rep; R.files];
    catch err
        warning('ff_emi_study:report', 'Report not written: %s', err.message);
    end
end

fprintf('\n===== Conclusions (automatic - check them against the figures) =====\n');
fprintf('- %s\n', conc);
fprintf('\nResults in %s\n', opt.OutDir);
if nargout == 0, clear R, end
end

%% ===== measurement ===========================================================

function S = measure(ip, opt, bands)
fprintf('Connecting to %s:%d ...\n', ip, opt.Port);
ff = FieldFox(ip, opt.Port);
fprintf('Connected: %s\n', ff.IDN);
if ~ff.hasOption("233"), error('ff_emi_study needs the spectrum analyzer (option 233).'); end
if opt.Preamp && ~ff.hasOption("235")
    warning('No preamplifier (option 235): measuring without it.');
    opt.Preamp = false;
end
ff.takeControl();                        % backs up every setting; release() restores them
cleanup = onCleanup(@() releaseInstrument(ff));
ff.setMode("SA");
S = newStudy("instrument", bands, opt);
S.idn = string(ff.IDN); S.options = ff.Options;
S.notes = [S.notes; prepareSA(ff)];
canWaor = probe(ff, "CALC:MEAS:WAOR?");   % ADC over-range flag of the last sweep
if ~canWaor, S.notes(end + 1, 1) = "CALC:MEAS:WAOR? not available: overload was not checked."; end

states = protocolStates(opt.Protocol);
segs = {}; zs = {};
targets = table();
try
    for s = 1:numel(states)
        if ~switchInverter(states(s), s, opt)
            S.aborted = true;
            fprintf('Stopped by the user before step %s.\n', states(s));
            break
        end
        t0 = datetime('now');
        for b = 1:height(bands)
            first = true;
            for det = opt.Detectors
                seg = acquire(ff, bands(b, :), b, det, states(s), opt, canWaor, first);
                seg.file = saveRecording(seg, S.idn, opt.OutDir);
                segs{end + 1} = seg; %#ok<AGROW>
                first = false;
            end
        end
        S.protocol = [S.protocol; {s, states(s), inverterOf(states(s)), t0, datetime('now')}];
        S.seg = [segs{:}];
        if states(s) == "ON", targets = pickTargets(S, opt); end
        if ~isempty(targets) && any(states(s) == ["ON" "OFF2"])
            zs = [zs, zeroSpan(ff, targets, states(s), opt)]; %#ok<AGROW>
        end
    end
    e = ff.errors();
    if ~isempty(e), S.notes(end + 1, 1) = "Instrument errors at the end: " + strjoin(e, "; "); end
catch err
    % Keep what was measured (an OFF1/ON pair can still be analysed); the
    % onCleanup below restores the instrument in any case.
    if isempty(segs), rethrow(err), end
    S.aborted = true;
    S.notes(end + 1, 1) = "Measurement stopped by an error: " + err.message;
    warning('ff_emi_study:measure', 'Measurement stopped by an error (%s): analysing what was measured.', err.message);
end
S.seg = [segs{:}];
S.zs = [zs{:}];
S.targets = targets;
clear cleanup                            % restore the instrument now
fprintf('Instrument released (all settings restored).\n');
end

function notes = prepareSA(ff)
% Plain dBm traces of a swept SA: the transducer is applied in MATLAB.
notes = strings(0, 1);
u = tryQuery(ff, "AMPL:UNIT?");          % PG: display unit, also the TRAC:DATA? unit (not yet run on this unit)
if strlength(u) > 0 && ~startsWith(upper(u), "DBM")
    ff.write("AMPL:UNIT DBM");
    notes(end + 1, 1) = "Display unit was " + u + ": set to dBm for the study.";
    if ~startsWith(upper(tryQuery(ff, "AMPL:UNIT?")), "DBM")
        notes(end + 1, 1) = "WARNING: the unit could not be set to dBm: levels are in " + u + ".";
        warning('ff_emi_study:unit', '%s', notes(end));
    end
end
for c = ["ANT" "CABL"]                   % PG: on-instrument antenna/cable corrections (not yet run on this unit)
    if ismember(upper(tryQuery(ff, "AMPL:CORR:" + c + "?")), ["1" "ON"])
        ff.write("AMPL:CORR:" + c + " 0");
        notes(end + 1, 1) = "Instrument " + c + " correction was ON: switched off (MATLAB applies the transducer).";
    end
end
ff.write("SENS:MEAS:CHAN NONE");         % no channel measurement (research: SENS:MEAS:CHAN)
ff.write("SENS:MEAS:INT NONE");          % no spectrogram display (research: SENS:MEAS:INT)
ff.write("TRIG:SOUR FREE");              % PG: a video/external trigger would stall FFT and zero-span sweeps (not yet run on this unit)
ff.opc(30);
e = ff.errors();
ff.write("SENS:BAND:VID:AUTO 1");        % UNVERIFIED on A.12.45: VBW coupled to RBW for every band
ff.opc(30);
eV = ff.errors();
if ~isempty(eV), notes(end + 1, 1) = "SENS:BAND:VID:AUTO not accepted (" + strjoin(eV, "; ") + "): VBW left as it was."; end
if ~isempty(e)
    notes(end + 1, 1) = "Instrument errors during setup (ignored): " + strjoin(e, "; ");
    warning('ff_emi_study:setup', '%s', notes(end));
end
end

function seg = acquire(ff, band, b, det, state, opt, canWaor, configure)
if configure
    cfg = {'Start', band.Start_Hz, 'Stop', band.Stop_Hz, 'RBWAuto', false, 'RBW', band.RBW_Hz, ...
        'Points', band.Points, 'Preamp', opt.Preamp, 'RefLevel', opt.RefLevel, ...
        'TraceType', "CLRW", 'Average', 1};
    if isnumeric(opt.Atten), cfg = [cfg, {'AttenAuto', false, 'Atten', opt.Atten}];
    else, cfg = [cfg, {'AttenAuto', true}]; end
    ff.saConfigure(cfg{:});
end
ff.saConfigure('Detector', det);
st = ff.saSettings();
n = opt.Sweeps;
P = zeros(n, st.points); t = NaT(n, 1); over = false;
tic;
for k = 1:n
    ff.sweep();
    y = ff.saTrace();
    if numel(y) ~= st.points, error('Trace has %d points, expected %d.', numel(y), st.points); end
    P(k, :) = y; t(k) = datetime('now');
    if canWaor
        w = ff.queryNum("CALC:MEAS:WAOR?");
        over = over || (isscalar(w) && w == 1);
    end
end
msg = '';
if over, msg = '  ADC OVER-RANGE: use less preamp / more attenuation'; end
fprintf('  [%-4s] %-22s %-5s %3d sweeps  %.2f s/sweep  RBW %s%s\n', state, band.Name, det, n, ...
    toc/n, fmtHz(st.rbw), msg);
settings = rmfield(st, 'f');
settings.binStep = (st.stop - st.start) / max(st.points - 1, 1);
settings.state = state; settings.band = band.Name; settings.bandIndex = b;
settings.detectorLabel = detLabel(det); settings.overload = over;
seg = makeSeg(state, b, det, st.f, P, t, settings, over);
end

function zs = zeroSpan(ff, targets, state, opt)
% Level vs time at each target frequency (zero span). The swept setup is saved
% to a .sta file first and loaded afterwards, so the zero-span sweep time cannot
% leak into the next swept band (only verified MMEM commands are used for that).
zs = {};
tmp = "emi_zs.sta";
saved = false;
try
    ff.saveState(tmp); saved = true;
catch err
    warning('ff_emi_study:zerospan', 'Could not save the swept setup (%s).', err.message);
end
for i = 1:height(targets)
    fc = targets.Freq_Hz(i);
    try
        ff.write("SENS:FREQ:SPAN 0");                              % research: zero span = FREQ:SPAN 0
        ff.write(sprintf("SENS:FREQ:CENT %.10g", fc));
        ff.write(sprintf("SENS:BAND:RES %.10g", targets.RBW_Hz(i)));
        ff.write(sprintf("SENS:SWE:POIN %d", opt.ZeroSpanPoints));
        ff.write("SENS:DET:FUNC SAMP");                            % envelope samples
        ff.opc(30);
        e = ff.errors();
        ff.write(sprintf("SENS:SWE:TIME %.6g", opt.ZeroSpanTime)); % UNVERIFIED on A.12.45: zero-span sweep time
        ff.opc(30);
        eT = ff.errors();
        if ~isempty(eT), e = [e; "SWE:TIME " + opt.ZeroSpanTime + " rejected: " + eT]; end %#ok<AGROW>
        span = ff.queryNum("SENS:FREQ:SPAN?");
        if isempty(span) || span ~= 0
            warning('ff_emi_study:zerospan', 'Zero span not accepted at %s (span %g): skipped.', fmtHz(fc), span);
            continue
        end
        N = ff.queryNum("SENS:SWE:POIN?");
        env = zeros(opt.ZeroSpanSweeps, N);
        for k = 1:opt.ZeroSpanSweeps
            ff.sweep();
            y = ff.saTrace();
            if numel(y) ~= N, error('Trace has %d points, expected %d.', numel(y), N); end
            env(k, :) = y;
        end
        mt = str2double(tryQuery(ff, "SENS:SWE:MTIM?"));          % measured sweep time (research: after INIT:IMM;*OPC?)
        T = NaN;
        if isempty(eT), T = str2double(tryQuery(ff, "SENS:SWE:TIME?")); end   % UNVERIFIED query
        src = "SWE:TIME?";
        if ~(T > 0) && mt > 0
            T = mt; src = "SWE:MTIM? (may include overhead)";
        elseif ~(T > 0)
            T = opt.ZeroSpanTime; src = "requested (not confirmed)";
        end
        z = struct('freq', fc, 'band', targets.Band(i), 'state', state, 'rbw', ff.queryNum("SENS:BAND:RES?"), ...
            'T', T, 'mtime', mt, 'timeSource', src, 'env', env, 't', datetime('now'), ...
            'errors', strjoin([e; ff.errors()], "; "));
    catch err
        warning('ff_emi_study:zerospan', 'Zero span at %s failed: %s', fmtHz(fc), err.message);
        continue
    end
    zs{end + 1} = z; %#ok<AGROW>
    fprintf('  [%-4s] zero span %s: %d x %d points over %.4g s (%s, RBW %s)\n', state, fmtHz(fc), ...
        opt.ZeroSpanSweeps, N, T, src, fmtHz(z.rbw));
    if strlength(z.errors) > 0, warning('ff_emi_study:zerospan', 'Instrument: %s', z.errors); end
end
if saved
    ff.loadState(tmp);
    ff.write("MMEM:DEL """ + tmp + """");
    ff.query("INIT:CONT 0;*OPC?", 30);                             % single sweep again
    ff.write("FORM:DATA REAL,32");       % what FieldFox.m last set (saTrace), in case the .sta changed it
    ff.errors();
end
end

function targets = pickTargets(S, opt)
% Zero-span frequencies: 'ZeroSpanFreqs' plus the 'ZeroSpanK' strongest new
% emissions from a quick OFF1 vs ON screening (Welch t-test), one per band first.
rows = zeros(0, 4);
det = primaryDet(S);
for b = 1:height(S.bands)
    off = getSeg(S, "OFF1", b, det); on = getSeg(S, "ON", b, det);
    if isempty(off) || isempty(on), continue, end
    [~, p] = ttest2(on.P, matchAxis(off, on.f), 'Vartype', 'unequal');
    p(isnan(p)) = 1;
    eff = mean(on.P, 1) - mean(matchAxis(off, on.f), 1);
    [s, e] = runs(fdrBH(p) < opt.Alpha & eff > opt.MinEffect, 1);
    for r = 1:numel(s)
        [~, i] = max(eff(s(r):e(r)));
        k = s(r) + i - 1;
        rows(end + 1, :) = [on.f(k), b, on.settings.rbw, eff(k)]; %#ok<AGROW>
    end
end
T = array2table(rows, 'VariableNames', {'Freq_Hz', 'Band', 'RBW_Hz', 'Effect_dB'});
T = sortrows(T, 'Effect_dB', 'descend');
keep = false(height(T), 1);
for pass = 1:2                            % pass 1: the strongest per band; pass 2: fill up
    for i = 1:height(T)
        if keep(i) || nnz(keep) >= opt.ZeroSpanK, continue, end
        if pass == 1 && any(T.Band(keep) == T.Band(i)), continue, end
        if ~any(abs(T.Freq_Hz(keep) - T.Freq_Hz(i)) < 3*T.RBW_Hz(i)), keep(i) = true; end
    end
end
T = T(keep, :);
fu = opt.ZeroSpanFreqs(:);
for i = 1:numel(fu)                       % user frequencies: RBW of the band they fall in
    b = find(S.bands.Start_Hz <= fu(i) & S.bands.Stop_Hz >= fu(i), 1);
    if isempty(b), warning('ff_emi_study:zerospan', '%s is outside the scan plan: skipped.', fmtHz(fu(i))); continue, end
    if any(abs(T.Freq_Hz - fu(i)) < 3*S.bands.RBW_Hz(b)), continue, end   % already a target
    T = [T; {fu(i), b, S.bands.RBW_Hz(b), NaN}]; %#ok<AGROW>
end
if ~isempty(opt.ZeroSpanRBW), T.RBW_Hz(:) = opt.ZeroSpanRBW; end
if isempty(T), fprintf('  No new emission found for zero span.\n'); end
targets = T;
end

function ok = switchInverter(state, step, opt)
want = inverterOf(state);
switch state
    case "OFF1", msg = "switch the inverter OFF (AC breaker and DC isolator).";
    case "ON",   msg = "switch the inverter ON and wait until it exports power.";
    otherwise,   msg = "switch the inverter OFF again.";
end
msg = sprintf('Step %d (%s): %s Keep the antenna/LISN, cables and people exactly as they are.', step, state, msg);
fprintf('\n%s\n', msg);
if ~isempty(opt.SwitchFcn), opt.SwitchFcn(want); end
ok = true;
if opt.Prompt
    if usejava('desktop')
        ok = strcmp(questdlg(msg, 'ff_emi_study', 'Continue', 'Abort', 'Continue'), 'Continue');
    else
        ok = ~strcmpi(strtrim(input('Press Enter when done (q = abort): ', 's')), 'q');
    end
end
if ok && opt.SettleTime > 0, pause(opt.SettleTime); end
end

function releaseInstrument(ff)
try
    ff.release();
catch err
    warning('ff_emi_study:release', 'Could not restore the instrument: %s', err.message);
end
end

function r = tryQuery(ff, cmd)
% Reply of a query that may not exist on this firmware ("" = no reply after 5 s;
% the tcpclient time-out warning is expected then and is not shown).
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
e = ff.errors();                         % always drain the queue (an unknown query leaves -113)
ok = strlength(r) > 0 && isempty(e);
end

%% ===== data: study, recordings, demo =========================================

function S = newStudy(source, bands, opt)
S = struct('source', source, 'bands', bands, 'seg', [], 'zs', [], 'targets', table(), ...
    'protocol', table('Size', [0 5], 'VariableTypes', {'double', 'string', 'string', 'datetime', 'datetime'}, ...
    'VariableNames', {'Step', 'State', 'Inverter', 'Start', 'End'}), 'idn', "", 'options', strings(1, 0), ...
    'notes', strings(0, 1), 'aborted', false, 'created', datetime('now', 'TimeZone', 'Asia/Bangkok'), ...
    'sweeps', opt.Sweeps, 'detectors', opt.Detectors, 'preamp', opt.Preamp, 'refLevel', opt.RefLevel, ...
    'protocolName', opt.Protocol, 'userNotes', string(opt.Notes));
S.atten = opt.Atten;
end

function seg = makeSeg(state, b, det, f, P, t, settings, over)
seg = struct('state', string(state), 'band', b, 'det', string(det), 'f', f(:).', 'P', P, ...
    't', t(:), 'settings', settings, 'overload', logical(over), 'file', "");
end

function S = loadData(data, opt, bands)
if isscalar(data) && lower(data) == "demo"
    S = demoStudy(opt, bands);
    return
end
if isscalar(data) && isfolder(data)    % a study folder: its emi_study.mat, else its rec\*.mat
    if isfile(fullfile(data, "emi_study.mat"))
        data = fullfile(data, "emi_study.mat");
    else
        d = [dir(fullfile(data, "rec_*.mat")); dir(fullfile(data, "rec", "rec_*.mat"))];
        if isempty(d), error('No emi_study.mat or rec_*.mat recordings in %s.', data); end
        data = string(fullfile({d.folder}, {d.name}));
    end
end
for i = 1:numel(data)
    if ~isfile(data(i)), error('File not found: %s', data(i)); end
end
if isscalar(data) && ismember('study', who('-file', data))
    L = load(data, 'study');
    S = L.study;
    if isempty(S.seg), error('%s holds no measured sweeps.', data); end
    S.source = "study file " + data;
    fprintf('Loaded study %s (%d sweep sets).\n', data, numel(S.seg));
    return
end
recs = cell(1, numel(data));
for i = 1:numel(data)
    recs{i} = load(data(i));
    if ~all(isfield(recs{i}, {'f', 'P'})), error('%s is not a recording (needs f and P).', data(i)); end
    if numel(recs{i}.f) < 5
        error(['%s is not a swept spectrum (%d frequency point): power meter and VVM logs are ' ...
            're-analysed with their own tool (''File'').'], data(i), numel(recs{i}.f));
    end
end
S = newStudy("recordings", table(), opt);
S.idn = recField(recs{1}, 'idn', "");
own = all(cellfun(@(r) isfield(r, 'settings') && isfield(r.settings, 'state'), recs));
if own                                   % recordings written by ff_emi_study
    names = strings(0, 1); segs = cell(1, numel(recs)); B = zeros(0, 4);
    for i = 1:numel(recs)
        st = recs{i}.settings;
        b = find(names == string(st.band), 1);
        if isempty(b)
            names(end + 1, 1) = string(st.band); b = numel(names); %#ok<AGROW>
            B(b, :) = [recs{i}.f(1), recs{i}.f(end), st.rbw, numel(recs{i}.f)];
        end
        segs{i} = makeSeg(st.state, b, detCode(st.detector), recs{i}.f, recs{i}.P, recField(recs{i}, 't', NaT(size(recs{i}.P, 1), 1)), st, st.overload);
        segs{i}.file = data(i);
    end
    S.bands = table(names, B(:, 1), B(:, 2), B(:, 3), B(:, 4), 'VariableNames', {'Name', 'Start_Hz', 'Stop_Hz', 'RBW_Hz', 'Points'});
    S.seg = [segs{:}];
    S = protocolFromSegs(S);
    return
end
% Generic recordings (e.g. fieldfox_live): 1 file = split at the largest time gap,
% 2 files = OFF, ON; 3 files = OFF, ON, OFF again.
if numel(recs) > 3, error('Give 1, 2 or 3 recordings (OFF, ON, OFF again).'); end
f = recs{1}.f(:).';
st = recField(recs{1}, 'settings', struct());
rbw = recField(st, 'rbw', f(2) - f(1));
det = "REC";
if isfield(st, 'detector'), det = detCode(st.detector); end
if isscalar(recs)
    t = recField(recs{1}, 't', NaT(0, 1));
    if numel(t) ~= size(recs{1}.P, 1) || numel(t) < 4, error('One recording needs time stamps to split it into OFF and ON.'); end
    [gap, i] = max(seconds(diff(t)));
    parts = {recs{1}.P(1:i, :), recs{1}.P(i + 1:end, :)}; times = {t(1:i), t(i + 1:end)};
    S.notes(end + 1, 1) = sprintf("One recording split at its largest time gap (%.0f s): first part = OFF, second = ON.", gap);
else
    parts = cellfun(@(r) r.P, recs, 'UniformOutput', false);
    times = cellfun(@(r) recField(r, 't', NaT(size(r.P, 1), 1)), recs, 'UniformOutput', false);
end
states = ["OFF1" "ON" "OFF2"];
segs = cell(1, numel(parts));
for i = 1:numel(parts)
    fi = f;
    if numel(recs) > 1, fi = recs{i}.f(:).'; end
    if numel(fi) ~= numel(f) || max(abs(fi - f)) > (f(2) - f(1))/2
        error('The recordings have different frequency axes: use the same span and points.');
    end
    segs{i} = makeSeg(states(i), 1, det, f, parts{i}, times{i}, st, false);
    if numel(recs) > 1, segs{i}.file = data(i); end
end
S.seg = [segs{:}];
S.bands = table("Recorded " + fmtHz(f(1)) + " - " + fmtHz(f(end)), f(1), f(end), rbw, numel(f), ...
    'VariableNames', {'Name', 'Start_Hz', 'Stop_Hz', 'RBW_Hz', 'Points'});
S.protocolName = "recordings";
S.sweeps = size(parts{1}, 1); S.detectors = det;
S = protocolFromSegs(S);
end

function S = protocolFromSegs(S)
states = unique([S.seg.state], 'stable');
for s = 1:numel(states)
    t = vertcat(S.seg([S.seg.state] == states(s)).t);
    t = t(~isnat(t));
    if isempty(t), t = NaT; end
    S.protocol = [S.protocol; {s, states(s), inverterOf(states(s)), min(t), max(t)}];
end
end

function v = recField(s, name, default)
if isfield(s, name), v = s.(name); else, v = default; end
if isstring(default) && ~isstring(v), v = string(v); end
end

function S = demoStudy(opt, bands)
% Synthetic study: ambient FM/AM/LTE carriers, inverter = 20 kHz harmonic comb
% below 30 MHz, broadband rise 30-120 MHz and a narrowband line at 148 MHz.
rs = RandStream('mt19937ar', 'Seed', 7);
S = newStudy("demo", bands, opt);
S.idn = "SYNTHETIC DEMO DATA (no instrument)";
S.notes(end + 1, 1) = "SYNTHETIC DEMO DATA: the spectra are simulated, not measured.";
states = protocolStates(opt.Protocol);
t0 = datetime('now') - minutes(30);
segs = {}; zs = {}; targets = table();
for s = 1:numel(states)
    tStart = t0;
    for b = 1:height(bands)
        f = linspace(bands.Start_Hz(b), bands.Stop_Hz(b), bands.Points(b));
        for det = opt.Detectors
            P = zeros(opt.Sweeps, numel(f)); t = NaT(opt.Sweeps, 1);
            for k = 1:opt.Sweeps
                P(k, :) = demoTrace(f, bands.RBW_Hz(b), det, states(s) == "ON", opt.Preamp, rs);
                t0 = t0 + seconds(0.3); t(k) = t0;
            end
            st = struct('start', f(1), 'stop', f(end), 'points', numel(f), 'rbw', bands.RBW_Hz(b), ...
                'vbw', bands.RBW_Hz(b), 'ref', opt.RefLevel, 'pdiv', 10, 'atten', 0, 'preamp', opt.Preamp, ...
                'detector', det, 'binStep', f(2) - f(1), 'state', states(s), 'band', bands.Name(b), ...
                'bandIndex', b, 'detectorLabel', detLabel(det), 'overload', false);
            segs{end + 1} = makeSeg(states(s), b, det, f, P, t, st, false); %#ok<AGROW>
        end
    end
    S.protocol = [S.protocol; {s, states(s), inverterOf(states(s)), tStart, t0}];
    S.seg = [segs{:}];
    if states(s) == "ON", targets = pickTargets(S, opt); end
    if ~isempty(targets) && any(states(s) == ["ON" "OFF2"])
        for i = 1:height(targets)
            env = zeros(opt.ZeroSpanSweeps, opt.ZeroSpanPoints);
            for k = 1:opt.ZeroSpanSweeps
                env(k, :) = demoEnvelope(targets.Freq_Hz(i), targets.RBW_Hz(i), opt.ZeroSpanPoints, ...
                    opt.ZeroSpanTime, states(s) == "ON", opt.Preamp, rs);
            end
            zs{end + 1} = struct('freq', targets.Freq_Hz(i), 'band', targets.Band(i), 'state', states(s), ...
                'rbw', targets.RBW_Hz(i), 'T', opt.ZeroSpanTime, 'mtime', opt.ZeroSpanTime, ...
                'timeSource', "synthetic", 'env', env, 't', t0, 'errors', ""); %#ok<AGROW>
        end
    end
    t0 = t0 + minutes(2);
end
S.zs = [zs{:}];
S.targets = targets;
end

function y = demoTrace(f, rbw, det, on, preamp, rs)
fl = -145 - 15*preamp + 10*log10(rbw);             % noise floor in the RBW (dBm)
g = @(f0, bw) exp(-2.7726*((f - f0)/bw).^2);       % 3 dB Gaussian RBW shape
p = 10^(fl/10) * ones(size(f));
for fm = [88.5 91 95.5 99 103.5 107]*1e6, p = p + 10^(-4.5) * g(fm, max(rbw, 100e3)); end
p = p + 1e-7 * g(1.2e6, rbw);
lte = f >= 925e6 & f <= 935e6;
p(lte) = p(lte) + 10^((-95 + 10*log10(rbw/9e3) + randn(rs))/10);
if on
    lvl = -55 - 20*log10(max(f, 150e3)/150e3);
    dx = f - round(f/20e3)*20e3;
    p = p + (f < 30e6) .* 10.^(lvl/10) .* exp(-2.7726*(dx/rbw).^2);
    h = f >= 30e6 & f <= 120e6;
    p(h) = p(h) + 10.^((fl + 1 + 6*(1 - (f(h) - 30e6)/90e6))/10);
    p = p + 1e-7 * g(148e6, rbw);
end
if det == "AVER"
    y = 10*log10(p) + 0.4*randn(rs, size(f));
else
    y = 10*log10(p) + 5*10^(fl/10)./p + 1.2*randn(rs, size(f));
end
end

function y = demoEnvelope(fc, rbw, N, T, on, preamp, rs)
fl = -145 - 15*preamp + 10*log10(rbw);
t = (0:N - 1)/(N - 1)*T;
p = 10.^((fl + 1.5*randn(rs, 1, N))/10);
if on && (fc < 120e6 || abs(fc - 148e6) < 2*rbw)
    lvl = -70;
    if fc < 30e6, lvl = -55 - 20*log10(max(fc, 150e3)/150e3); end
    burst = mod(t + 0.01*rand(rs), 0.01) < 0.003;    % 100 Hz bursts, 30 % duty
    p = p + burst .* 10^(lvl/10) .* (1 + 0.3*sin(2*pi*20e3*t));
end
y = 10*log10(p);
end

function S = saveRecordings(S, outDir)
for i = 1:numel(S.seg)
    S.seg(i).file = saveRecording(S.seg(i), S.idn, outDir);
end
end

function file = saveRecording(seg, idn, outDir)
% Same format as fieldfox_live: f, P, t, settings, idn.
d = fullfile(outDir, "rec");
if ~isfolder(d), mkdir(d); end
file = fullfile(d, sprintf("rec_%s_b%d_%s.mat", seg.state, seg.band, seg.det));
out = struct('f', seg.f, 'P', seg.P, 't', seg.t, 'settings', seg.settings, 'idn', idn);
save(file, '-struct', 'out');
end

%% ===== analysis ==============================================================

function res = analyzeStudy(S, opt, tr)
dets = unique([S.seg.det], 'stable');
res = {};
for b = 1:height(S.bands)
    for det = dets
        off1 = getSeg(S, "OFF1", b, det); on = getSeg(S, "ON", b, det); off2 = getSeg(S, "OFF2", b, det);
        if isempty(off1) || isempty(on), continue, end
        f = on.f;
        P2 = [];
        if ~isempty(off2), P2 = matchAxis(off2, f); end
        tic;
        r = binStats(f, matchAxis(off1, f), on.P, P2, opt);
        r.band = b; r.bandName = S.bands.Name(b); r.det = det; r.detLabel = detLabel(det);
        r.f = f; r.df = f(2) - f(1); r.rbw = on.settings.rbw;
        r.nSweeps = [size(off1.P, 1), size(on.P, 1), size(P2, 1)];
        [r.corr, r.unit] = transducerCorr(tr, f, opt);
        r.overload = off1.overload || on.overload || (~isempty(off2) && off2.overload);
        fprintf('  %-22s %-8s %5d points: %4d flagged, floor %+.1f dB  (%.1f s)\n', r.bandName, r.detLabel, ...
            numel(f), nnz(r.flag), r.floorChange, toc);
        res{end + 1} = r; %#ok<AGROW>
    end
end
res = [res{:}];
end

function r = binStats(f, off1, on, off2, opt)
off = [off1; off2];
nb = numel(f);
if opt.Test == "ttest2"
    [~, p] = ttest2(on, off, 'Vartype', 'unequal');
else
    p = ones(1, nb);
    for k = 1:nb, p(k) = ranksum(off(:, k), on(:, k)); end
end
p(isnan(p)) = 1;
r.p = p;
r.q = fdrBH(p);
% Smallest p the rank-sum test can give with these sweep counts: if it is above
% Alpha/points, an isolated narrowband emission can never pass the FDR.
r.pmin = 0;
if opt.Test == "ranksum", r.pmin = ranksum(1:size(off, 1), size(off, 1) + (1:size(on, 1))); end
% A-B-A effect, sweep by sweep: ON minus the mean of the OFF sweeps before and after
n = min(size(on, 1), size(off1, 1));
if ~isempty(off2), n = min(n, size(off2, 1)); end
ref = off1(1:n, :);
if ~isempty(off2), ref = (ref + off2(1:n, :))/2; end
d = on(1:n, :) - ref;
r.effect = mean(d, 1);
m = median(d, 2);                     % per-sweep shift of the whole band (ignores single carriers)
r.floorChange = mean(m);
if n >= 3
    ci = bootci(opt.Boot, {@mean, d}, 'Alpha', opt.Alpha, 'Type', 'per');
    r.floorCI = bootci(opt.Boot, {@mean, m}, 'Alpha', opt.Alpha, 'Type', 'per').';
else
    ci = nan(2, nb); r.floorCI = [NaN NaN];
end
r.ciLo = ci(1, :); r.ciHi = ci(2, :);
lo = r.ciLo; lo(isnan(lo)) = r.effect(isnan(lo));   % fewer than 3 sweep pairs: no CI, use the effect
r.flag = r.q < opt.Alpha & lo > opt.MinEffect;
r.offMean = powMean(off); r.onMean = powMean(on);
r.offMax = max(off, [], 1); r.onMax = max(on, [], 1);
r.floorOff = median(r.offMean); r.floorOn = median(r.onMean);
r.drift = false(1, nb); r.driftFrac = NaN;
if ~isempty(off2)                     % OFF1 vs OFF2: did the ambient change during the test?
    [~, p2] = ttest2(off1, off2, 'Vartype', 'unequal');
    p2(isnan(p2)) = 1;
    r.drift = fdrBH(p2) < opt.Alpha & abs(mean(off2, 1) - mean(off1, 1)) > opt.MinEffect;
    r.driftFrac = mean(r.drift);
end
end

function [res, E] = findEmissions(res, opt)
% Groups flagged points into emissions, finds the harmonic comb, labels sources.
E = cell(numel(res), 1);
[method, note] = signalIdMethod(opt);
if strlength(note) > 0, fprintf('  %s\n', note); end
warned = false;
for i = 1:numel(res)
    r = res(i);
    res(i).comb = combDetect(r, opt);
    sp = res(i).comb.spacing;
    [s, e] = runs(r.flag, 1);
    if isempty(s), continue, end
    ids = table();
    if method ~= "off"
        try                                % optional suite module: never let it stop the study
            t0 = tic;
            ids = ff_signal_id(r.f, r.onMean, r.rbw, 'Method', method, 'Detector', idDetector(r.det), ...
                'Visible', 'off', 'Prompt', false);
            fprintf('  ff_signal_id %-22s %-8s %3d signals (%s, %.1f s)\n', r.bandName, r.detLabel, height(ids), method, toc(t0));
        catch err
            if ~warned, warning('ff_emi_study:signalid', 'ff_signal_id failed: %s', err.message); end
            warned = true;
        end
    end
    n = numel(s);
    T = emptyEmissions(n);
    T.Band(:) = r.bandName; T.Detector(:) = r.detLabel; T.Unit(:) = r.unit; T.RBW_Hz(:) = r.rbw;
    for j = 1:n
        idx = s(j):e(j);
        [~, i0] = max(r.onMean(idx)); k = idx(i0);
        T.Freq_Hz(j) = r.f(k);
        T.Width_Hz(j) = numel(idx)*r.df;
        T.Shape(j) = "broadband";
        if T.Width_Hz(j) <= 3*max(r.df, r.rbw)
            T.Shape(j) = "narrowband";
            nb = abs(r.f - r.f(k)) <= 25*max(r.df, r.rbw);      % a piece of a wide hump, not a line?
            if median(r.effect(nb)) > opt.MinEffect/2 && isnan(sp), T.Shape(j) = "narrowband on a broadband rise"; end
        end
        T.Kind(j) = "raised";
        if r.offMean(k) < r.floorOff + 6, T.Kind(j) = "new"; end
        T.Level_ON(j) = r.onMean(k) + r.corr(k); T.Level_OFF(j) = r.offMean(k) + r.corr(k);
        T.Effect_dB(j) = r.effect(k); T.CI_low_dB(j) = r.ciLo(k); T.CI_high_dB(j) = r.ciHi(k);
        T.q(j) = min(r.q(idx)); T.Drift(j) = any(r.drift(idx));
        if ~isnan(sp) && abs(r.f(k) - round(r.f(k)/sp)*sp) <= max(r.df, r.rbw)
            T.Harmonic(j) = round(r.f(k)/sp);
        end
        T.Service(j) = serviceBand(r.f(k), max(r.df, r.rbw)/2);
        [T.SignalID(j), T.ID_Score(j)] = matchId(ids, r.f(k), max(3*r.df, r.rbw));
        T.LikelySource(j) = likelySource(T(j, :));
    end
    E{i} = sortrows(T, 'Effect_dB', 'descend');
end
E = vertcat(E{:});
if isempty(E), E = emptyEmissions(0); end
end

function T = emptyEmissions(n)
T = table(strings(n, 1), strings(n, 1), zeros(n, 1), zeros(n, 1), strings(n, 1), strings(n, 1), ...
    zeros(n, 1), zeros(n, 1), strings(n, 1), zeros(n, 1), zeros(n, 1), zeros(n, 1), zeros(n, 1), ...
    false(n, 1), nan(n, 1), strings(n, 1), strings(n, 1), nan(n, 1), strings(n, 1), zeros(n, 1), ...
    'VariableNames', {'Band', 'Detector', 'Freq_Hz', 'Width_Hz', 'Shape', 'Kind', 'Level_ON', 'Level_OFF', ...
    'Unit', 'Effect_dB', 'CI_low_dB', 'CI_high_dB', 'q', 'Drift', 'Harmonic', 'Service', 'SignalID', ...
    'ID_Score', 'LikelySource', 'RBW_Hz'});
end

function [m, note] = signalIdMethod(opt)
% "auto": the full ensemble only when its trained models are cached (training
% takes ~5 min); otherwise the template fit (no training).
m = opt.SignalID; note = "";
if m == "off", return, end
if ~exist('ff_signal_id', 'file')
    m = "off"; note = "ff_signal_id not found: emissions are labelled with the band table only.";
    return
end
if m == "auto"
    cache = fullfile(fileparts(which('ff_signal_classifier')), 'ff_signal_models.mat');
    if exist('ff_signal_classifier', 'file') && isfile(cache)
        m = "ensemble";
    else
        m = "template"; note = "ff_signal_id: no trained models cached, using the template fit (run ff_signal_classifier once for the full ensemble).";
    end
end
end

function d = idDetector(det)
switch det
    case "AVER", d = "average";
    case "SAMP", d = "sample";
    otherwise, d = "peak";
end
end

function [tech, score] = matchId(ids, f0, tol)
tech = ""; score = NaN;
if isempty(ids) || height(ids) == 0, return, end
try
    [d, j] = min(abs(ids.Freq_Hz - f0));
    if ismember('OBW_Hz', ids.Properties.VariableNames), tol = max(tol, ids.OBW_Hz(j)/2); end
    if d <= tol, tech = string(ids.Technology(j)); score = double(ids.Score(j)); end
catch
end
end

function src = likelySource(row)
% "CW" only means narrower than the RBW (an inverter line looks like that too); a
% service is claimed only for a confident ff_signal_id label (score >= 0.7).
service = strlength(row.SignalID) > 0 && ~ismember(row.SignalID, ["Unknown" "Unresolved" "CW"]) && row.ID_Score >= 0.7;
if ~isnan(row.Harmonic)
    src = "inverter: switching harmonic " + row.Harmonic;
elseif row.Drift
    src = "ambient drift (OFF1 differs from OFF2): check";
elseif service && row.Kind == "raised"
    src = "licensed service (" + row.SignalID + "): level change, probably traffic";
elseif service
    src = "looks like " + row.SignalID + " but new when ON: check";
elseif strlength(row.Service) > 0 && row.Kind == "raised"
    src = "existing signal in " + row.Service + " band: check";
else
    src = "inverter (likely)";
end
end

function c = combDetect(r, opt)
% Harmonic comb in the ON-OFF difference: cepstrum, then a fit of the lines.
c = struct('found', false, 'spacing', NaN, 'ci', [NaN NaN], 'nHarm', 0, 'fMax', NaN, 'm', [], ...
    'fpk', [], 'q', [], 'pq', [], 'q0', NaN, 'cepSNR', NaN, 'rmse', NaN);
f = r.f; df = r.df;
minSp = max(2.5*df, r.rbw); maxSp = (f(end) - f(1))/8;
if nnz(r.flag) < 5 || maxSp <= 2*minSp, return, end
x = max(r.effect, 0); x(~isfinite(x)) = 0; x = x - mean(x);
% Power cepstrum = spectrum of the (dB) spectrum. With the spectrum sampled every
% df Hz, the "quefrency" axis is in seconds: a comb of spacing D shows at 1/D.
[pq, qf] = periodogram(x(:), hamming(numel(x)), 2^nextpow2(8*numel(x)), 1/df);
ok = qf >= 1/maxSp & qf <= 1/minSp;
if ~any(ok), return, end
pk = pq .* ok;
[pmax, i0] = max(pk);
dq = qf(2) - qf(1);
for k = 4:-1:2                        % prefer the first rahmonic (= the largest spacing)
    j = find(ok & abs(qf - qf(i0)/k) <= 2*dq);
    if ~isempty(j) && max(pq(j)) >= 0.5*pmax
        [~, jj] = max(pq(j)); i0 = j(jj);
        break
    end
end
sp = 1/qf(i0);
c.q = qf; c.pq = pq; c.q0 = qf(i0); c.cepSNR = pq(i0)/median(pq(ok));
% A harmonic counts when a significant line near m*sp stands out of its neighbourhood
m = max(1, ceil(f(1)/sp)):floor(f(end)/sp);
hit = false(size(m)); fpk = nan(size(m));
for h = 1:numel(m)
    w = find(abs(f - m(h)*sp) <= max(df, sp/4));
    if isempty(w), continue, end
    [v, iw] = max(r.effect(w)); k = w(iw);
    nbh = abs(f - m(h)*sp) <= sp/2 + df;
    hit(h) = v >= opt.MinEffect && r.q(k) < opt.Alpha && v - min(r.effect(nbh)) >= opt.MinEffect;
    fpk(h) = f(k);
end
if nnz(hit) < 5 || nnz(hit)/(find(hit, 1, 'last') - find(hit, 1) + 1) < 0.5, return, end
mdl = fitlm(m(hit).', fpk(hit).', 'Intercept', false);
% A real comb has its lines on the multiples within one frequency cell; peaks of
% noise on a broadband rise scatter over the whole search window.
if mdl.RMSE > max(df, r.rbw), return, end
c.found = true;
c.spacing = mdl.Coefficients.Estimate(1);
c.ci = coefCI(mdl, opt.Alpha);
c.nHarm = nnz(hit); c.m = m(hit); c.fpk = fpk(hit); c.fMax = max(fpk(hit)); c.rmse = mdl.RMSE;
end

function [res, E, T] = applyLimits(res, E, opt, tr)
% INDICATIVE comparison with CISPR class B limits (only with a transducer).
T = table();
if ~isempty(E) && height(E) > 0, E.Margin_dB = nan(height(E), 1); end
rows = {};
for i = 1:numel(res)
    r = res(i);
    L = nan(size(r.f)); name = "";
    if ~isempty(tr) && r.unit ~= "dBm", [L, name] = cisprLimit(r.f, r.unit, r.det, opt); end
    res(i).limit = L; res(i).limitName = name;
    if ~any(isfinite(L)), continue, end
    lvl = r.onMean + r.corr; amb = r.offMean + r.corr;
    if r.det ~= "AVER", lvl = r.onMax + r.corr; amb = r.offMax + r.corr; end   % peak: max hold
    mg = lvl - L;
    [worst, k] = max(mg);
    fw = NaN; ff = NaN;                    % worst margin among the flagged (inverter) points
    m2 = mg; m2(~r.flag) = NaN;
    if any(isfinite(m2)), [fw, kf] = max(m2); ff = r.f(kf); end
    rows(end + 1, :) = {r.bandName, r.detLabel, name, r.unit, fw, ff, worst, r.f(k), nnz(mg > 0), max(amb - L)}; %#ok<AGROW>
    if ~isempty(E) && height(E) > 0
        sel = find(E.Band == r.bandName & E.Detector == r.detLabel);
        for j = sel.'
            [~, kk] = min(abs(r.f - E.Freq_Hz(j)));
            E.Margin_dB(j) = lvl(kk) - L(kk);
        end
    end
end
if ~isempty(rows)
    T = cell2table(rows, 'VariableNames', {'Band', 'Detector', 'Limit', 'Unit', 'InverterMargin_dB', ...
        'InverterFreq_Hz', 'WorstMargin_dB', 'WorstFreq_Hz', 'PointsOverLimit', 'AmbientWorstMargin_dB'});
end
end

function [L, name] = cisprLimit(f, unit, det, opt)
% Class B limits. Conducted = AC mains port (CISPR 11 group 1 = CISPR 32), dBuV.
% Radiated = dBuV/m at 10 m (30 MHz-1 GHz) and 3 m (CISPR 32 above 1 GHz),
% moved to MeasDistance with 20*log10(d_limit/d) (far field only).
L = nan(size(f));
qp = det ~= "AVER";                       % peak (or unknown) data vs the QP/PK limit
d = opt.MeasDistance;
if qp, dn = "QP"; else, dn = "AV"; end
if unit == "dBuV"
    a = f >= 150e3 & f < 500e3; x = log10(f(a)/150e3)/log10(500/150);
    if qp, L(a) = 66 - 10*x; L(f >= 500e3 & f < 5e6) = 56; L(f >= 5e6 & f <= 30e6) = 60;
    else,  L(a) = 56 - 10*x; L(f >= 500e3 & f < 5e6) = 46; L(f >= 5e6 & f <= 30e6) = 50; end
    name = opt.LimitStandard + " class B mains port, " + dn;
else
    if qp
        L(f >= 30e6 & f <= 230e6) = 30; L(f > 230e6 & f <= 1e9) = 37;
        L = L + 20*log10(10/d);
        dn = "QP (PK above 1 GHz)";
    end
    if opt.LimitStandard == "CISPR32"
        hi = [50 54] + 20*qp;             % AV 50/54, PK 70/74 dBuV/m at 3 m
        L(f > 1e9 & f <= 3e9) = hi(1) + 20*log10(3/d);
        L(f > 3e9 & f <= 6e9) = hi(2) + 20*log10(3/d);
    end
    name = sprintf("%s class B radiated at %g m, %s", opt.LimitStandard, d, dn);
end
end

function [T, Z] = zeroSpanMetrics(S, res)
% Burst timing from the zero-span envelopes (Signal Processing Toolbox).
T = table(); Z = struct('freq', {}, 'state', {}, 't', {}, 'env', {}, 'lev', {}, 'fp', {}, 'ps', {});
if isempty(S.zs), return, end
sp = NaN;
for i = 1:numel(res)
    if isfield(res(i), 'comb') && res(i).comb.found, sp = res(i).comb.spacing; break, end
end
rows = {};
for i = 1:numel(S.zs)
    z = S.zs(i);
    [ns, N] = size(z.env);
    fs = (N - 1)/z.T;
    x = 10.^(z.env/20);                   % linear envelope (sqrt of mW)
    thr = noiseContrast(N, fs, z.rbw);    % contrast that noise alone can reach
    [con, pw, pp, dc] = deal(nan(ns, 1));
    ps = 0; lev = [NaN NaN];
    for k = 1:ns
        lv = statelevels(x(k, :));
        con(k) = 20*log10(lv(2)/lv(1));
        if k == 1, lev = 20*log10(lv); end
        if con(k) >= thr                  % only a real on/off pattern has pulses
            try
                pw(k) = mean(pulsewidth(x(k, :), fs, 'StateLevels', lv));
                pp(k) = mean(pulseperiod(x(k, :), fs, 'StateLevels', lv));
                dc(k) = mean(dutycycle(x(k, :), fs, 'StateLevels', lv));
            catch
            end
        end
        [p1, fp] = pspectrum(x(k, :) - mean(x(k, :)), fs);
        ps = ps + p1/ns;
    end
    ok = fp > 2/z.T;
    [~, j] = max(ps .* ok);
    envPk = fp(j);
    period = median(pp, 'omitnan');
    cv = std(pp, 'omitnan')/mean(pp, 'omitnan');   % sweep-to-sweep period spread
    src = "requested";
    if isfield(z, 'timeSource'), src = z.timeSource; end
    rows(end + 1, :) = {z.freq, S.bands.Name(z.band), z.state, z.rbw, z.T, src, ns, median(con), thr, ...
        1e3*median(pw, 'omitnan'), 1e3*period, 1/period, 100*median(dc, 'omitnan'), envPk, ...
        interpretRate(1/period, envPk, median(con), sp, thr, cv)}; %#ok<AGROW>
    Z(end + 1) = struct('freq', z.freq, 'state', z.state, 't', (0:N - 1)/fs, 'env', z.env(1, :), ...
        'lev', lev, 'fp', fp, 'ps', ps); %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'Freq_Hz', 'Band', 'State', 'RBW_Hz', 'SweepTime_s', 'TimeFrom', 'Sweeps', ...
    'Contrast_dB', 'NoiseContrast_dB', 'PulseWidth_ms', 'Period_ms', 'Rate_Hz', 'DutyCycle_pct', 'EnvPeak_Hz', 'Interpretation'});
end

function thr = noiseContrast(N, fs, rbw)
% 99th percentile of the statelevels contrast of a noise-only envelope with the
% same length, sample rate and RBW (Monte Carlo). Noise alone gives 7-9 dB, so a
% fixed threshold would report bursts in pure noise.
persistent cache
if isempty(cache), cache = containers.Map(); end
key = sprintf('%d_%.6g_%.6g', N, fs, rbw);
if isKey(cache, key), thr = cache(key); return, end
old = rng(1);
restoreRng = onCleanup(@() rng(old));
sig = 0.265*fs/rbw;                       % Gaussian RBW filter (3 dB = RBW), in samples
if sig > 0.5
    h = exp(-0.5*((-ceil(4*sig):ceil(4*sig))/sig).^2);
else
    h = 1;                                % samples further apart than 1/RBW: independent
end
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

function s = interpretRate(rate, envPk, contrast, sp, thr, cv)
if contrast < thr || isnan(rate)
    if contrast < thr
        s = sprintf("noise-like: contrast %.1f dB is within what noise alone gives (< %.1f dB)", contrast, thr);
    else
        s = "no regular bursts";
        if ~isnan(envPk), s = s + sprintf(", envelope peak %.4g Hz", envPk); end
    end
    return
end
if cv > 0.2
    s = sprintf("irregular bursts (period varies %.0f%% between sweeps)", 100*cv);
    return
end
near = @(a, b) abs(a - b) <= 0.05*b;
if near(rate, 100)
    s = "bursts at 100 Hz = 2 x 50 Hz mains: follows the rectified mains voltage";
elseif near(rate, 50)
    s = "bursts at 50 Hz: mains synchronous";
elseif ~isnan(sp) && near(rate, sp)
    s = "repeats at the switching frequency (comb spacing)";
else
    s = sprintf("bursts repeat at %.4g Hz", rate);
end
end

function [T, curves] = rangeEstimate(E, opt)
% Distance at which each emission falls below the noise of a victim receiver.
% Far-field formulas only: the measuring antenna and the victim must be farther
% than lambda/2pi (1.6 m at 30 MHz, 320 m at 150 kHz); closer is left out (NaN).
T = table(); curves = struct('label', {}, 'd', {}, 'mF', {}, 'mC', {});
if isempty(E) || height(E) == 0, return, end
P = E(E.Detector == "Peak", :);
if isempty(P), P = E; end
c0 = physconst('LightSpeed');
P.Valid = P.Unit ~= "dBuV" & opt.MeasDistance > c0./P.Freq_Hz/(2*pi);
P = sortrows(P, {'Valid', 'Effect_dB'}, {'descend', 'descend'});
[~, iu] = unique(round(P.Freq_Hz/1e3), 'stable');
P = P(sort(iu(1:min(5, end))), :);
d = logspace(0, 5, 301);                  % 1 m to 100 km
rx = rxsite("cartesian", "AntennaPosition", [d; zeros(size(d)); 2*ones(size(d))]);
rx0 = rxsite("cartesian", "AntennaPosition", [opt.MeasDistance; 0; 2]);
pmF = propagationModel("freespace");
% close-in (Sun et al. 2016), n = 2.9 from 1 m. The toolbox returns the mean of
% NumDataPoints shadow-fading draws: with the default Sigma this is the median
% loss +-0.2 dB (keep Sigma at its default: a tiny Sigma biases it by -5.5 dB).
pmC = propagationModel("close-in");
rs = rng; rng(1);                         % repeatable +-0.2 dB jitter
restore = onCleanup(@() rng(rs));
rows = {};
for i = 1:height(P)
    fe = P.Freq_Hz(i);
    lam = c0/fe;
    exc = 10*log10(max(10^(P.Level_ON(i)/10) - 10^(P.Level_OFF(i)/10), eps));   % inverter part only
    tx = txsite("cartesian", "AntennaPosition", [0; 0; 2], "TransmitterFrequency", fe);
    switch P.Unit(i)
        case "dBm",    eirp = exc - opt.RxGain + pathloss(pmF, rx0, tx);
        case "dBuV/m", eirp = exc + 20*log10(opt.MeasDistance) - 104.77;    % E = sqrt(30 P)/d
        otherwise,     eirp = NaN;        % conducted (LISN): no radiated estimate
    end
    thr = -174 + 10*log10(P.RBW_Hz(i)) + opt.VictimNF + opt.INR;
    mF = eirp + opt.VictimGain - pathloss(pmF, rx, tx) - thr;
    mC = eirp + opt.VictimGain - pathloss(pmC, rx, tx) - thr;
    near = d < lam/(2*pi);
    mF(near) = NaN; mC(near) = NaN;
    note = "far field beyond " + sprintf('%.3g', lam/(2*pi)) + " m";
    if P.Unit(i) == "dBuV"
        note = "conducted (LISN) measurement: no radiated estimate";
    elseif ~P.Valid(i)
        note = sprintf("NOT estimated: measured in the near field (%g m < lambda/2pi = %.3g m); the noise couples through the cables", ...
            opt.MeasDistance, lam/(2*pi));
        eirp = NaN; mF(:) = NaN; mC(:) = NaN;
    end
    rows(end + 1, :) = {fe, P.Band(i), P.Level_ON(i), P.Unit(i), eirp, thr, reach(mF, d), reach(mC, d), note}; %#ok<AGROW>
    if P.Valid(i), curves(end + 1) = struct('label', string(fmtHz(fe)), 'd', d, 'mF', mF, 'mC', mC); end %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'Freq_Hz', 'Band', 'Level_ON', 'Unit', 'EIRP_dBm', ...
    'Threshold_dBm', 'Range_FreeSpace_m', 'Range_CloseIn_m', 'Note'});
end

function r = reach(margin, d)
% Largest distance with margin >= 0 (0 = below the threshold even at 1 m, Inf = beyond 100 km).
i = find(margin >= 0, 1, 'last');
if isempty(i) || all(isnan(margin)), r = 0;
elseif i == numel(d), r = Inf;
else, r = d(i);
end
if all(isnan(margin)), r = NaN; end
end

function T = summaryTable(res)
n = numel(res);
T = table(strings(n, 1), strings(n, 1), zeros(n, 1), zeros(n, 1), strings(n, 1), zeros(n, 1), zeros(n, 1), ...
    strings(n, 1), zeros(n, 1), zeros(n, 1), zeros(n, 1), zeros(n, 1), nan(n, 1), false(n, 1), ...
    'VariableNames', {'Band', 'Detector', 'Points', 'RBW_Hz', 'Sweeps_OFF1_ON_OFF2', 'FloorOFF', 'FloorON', ...
    'Unit', 'FloorChange_dB', 'CI_low_dB', 'CI_high_dB', 'FlaggedPoints', 'DriftPct', 'Overload'});
for i = 1:n
    r = res(i);
    T.Band(i) = r.bandName; T.Detector(i) = r.detLabel; T.Points(i) = numel(r.f); T.RBW_Hz(i) = r.rbw;
    T.Sweeps_OFF1_ON_OFF2(i) = strjoin(string(r.nSweeps), "/");
    T.FloorOFF(i) = median(r.offMean + r.corr, 'omitnan'); T.FloorON(i) = median(r.onMean + r.corr, 'omitnan');
    T.Unit(i) = r.unit; T.FloorChange_dB(i) = r.floorChange;
    T.CI_low_dB(i) = r.floorCI(1); T.CI_high_dB(i) = r.floorCI(2);
    T.FlaggedPoints(i) = nnz(r.flag); T.DriftPct(i) = 100*r.driftFrac; T.Overload(i) = r.overload;
end
end

function T = combTable(res)
rows = cell(0, 9);
for i = 1:numel(res)
    c = res(i).comb;
    if c.found
        rows(end + 1, :) = {res(i).bandName, res(i).detLabel, c.spacing, c.ci(1), c.ci(2), c.nHarm, ...
            c.fMax, c.cepSNR, c.rmse}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'Band', 'Detector', 'Spacing_Hz', 'CI_low_Hz', 'CI_high_Hz', ...
    'Harmonics', 'HighestHarmonic_Hz', 'CepstrumSNR', 'LineRMSE_Hz'});
end

function txt = conclusions(S, res, E, combT, zsT, limT, rangeT, opt, tr)
txt = strings(0, 1);
if isSynthetic(S)
    txt(end + 1) = "SYNTHETIC DEMO DATA - these are not measurements.";
end
txt(end + 1) = sprintf("Data: %s; protocol %s; %s.", S.idn, strjoin(unique(S.protocol.State, 'stable'), "-"), S.source);
if S.aborted, txt(end + 1) = "The protocol was stopped early: only the steps measured are analysed."; end
weak = arrayfun(@(r) r.pmin > opt.Alpha/numel(r.f), res);
if any(weak)
    txt(end + 1) = sprintf("Too few sweeps for the rank-sum test in %s: a single narrowband emission cannot reach the FDR threshold there (only wide changes can). Use 20 or more sweeps, or 'Test' ""ttest2"".", ...
        strjoin(unique([res(weak).bandName]), ", "));
end
pd = primaryDet(S);
for i = find([res.det] == pd)
    r = res(i);
    ci = r.floorCI;
    if ci(1) > 1
        txt(end + 1) = sprintf("%s: inverter ON raised the noise floor by %.1f dB (%.0f%% CI %.1f to %.1f dB, %s detector).", ...
            r.bandName, r.floorChange, 100*(1 - opt.Alpha), ci(1), ci(2), r.detLabel); %#ok<AGROW>
    elseif nnz(r.flag) == 0
        txt(end + 1) = sprintf("%s: no significant change (floor %+.1f dB, CI %.1f to %.1f dB; floor %.1f %s in %s RBW).", ...
            r.bandName, r.floorChange, ci(1), ci(2), r.floorOff + median(r.corr, 'omitnan'), r.unit, fmtHz(r.rbw)); %#ok<AGROW>
    end
    if ~isempty(E) && height(E) > 0
        Eb = E(E.Band == r.bandName & E.Detector == r.detLabel, :);
        if height(Eb) > 0
            [~, j] = max(Eb.Effect_dB);
            nh = nnz(~isnan(Eb.Harmonic));
            extra = "";
            if nh > 0, extra = sprintf(" (%d on the harmonic comb)", nh); end
            txt(end + 1) = sprintf("%s: %d emission(s) flagged%s, %d new; strongest %+.1f dB (CI %.1f to %.1f) at %s, %s.", ...
                r.bandName, height(Eb), extra, nnz(Eb.Kind == "new"), Eb.Effect_dB(j), Eb.CI_low_dB(j), ...
                Eb.CI_high_dB(j), fmtHz(Eb.Freq_Hz(j)), Eb.LikelySource(j)); %#ok<AGROW>
        end
    end
    if contains(r.bandName, "GNSS", 'IgnoreCase', true) && nnz(r.flag) == 0
        txt(end + 1) = sprintf("GNSS L1: no rise above the analyzer floor (%.1f dBm in %s). GNSS receivers work far below this, so this does not prove there is no effect on GPS.", ...
            r.floorOff, fmtHz(r.rbw)); %#ok<AGROW>
    end
    if r.driftFrac > 0.05
        txt(end + 1) = sprintf("%s: OFF1 and OFF2 differ in %.0f%% of the points: the ambient changed during the test; treat those points with care.", ...
            r.bandName, 100*r.driftFrac); %#ok<AGROW>
    end
end
for i = 1:height(combT)
    txt(end + 1) = sprintf("Harmonic comb in %s (%s): spacing %.4g kHz (CI %.4g to %.4g kHz), %d lines up to %s. This is the inverter switching frequency (or a multiple of it).", ...
        combT.Band(i), combT.Detector(i), combT.Spacing_Hz(i)/1e3, combT.CI_low_Hz(i)/1e3, combT.CI_high_Hz(i)/1e3, ...
        combT.Harmonics(i), fmtHz(combT.HighestHarmonic_Hz(i))); %#ok<AGROW>
end
zsOn = [];
if height(zsT) > 0, zsOn = find(zsT.State == "ON").'; end
for i = zsOn
    s = sprintf("Zero span at %s (inverter ON): %s", fmtHz(zsT.Freq_Hz(i)), zsT.Interpretation(i));
    if ~isnan(zsT.DutyCycle_pct(i)), s = s + sprintf("; duty cycle %.0f%%, pulse width %.3g ms", zsT.DutyCycle_pct(i), zsT.PulseWidth_ms(i)); end
    j = find(zsT.State ~= "ON" & abs(zsT.Freq_Hz - zsT.Freq_Hz(i)) < 1, 1);
    if ~isempty(j), s = s + sprintf(" (inverter OFF: %s)", zsT.Interpretation(j)); end
    txt(end + 1) = s + "."; %#ok<AGROW>
end
for i = find([res.overload])
    txt(end + 1) = sprintf("ADC over-range in %s (%s): strong signals - repeat with the preamp off or more attenuation.", ...
        res(i).bandName, res(i).detLabel); %#ok<AGROW>
end
if ~isempty(tr) && any([res.unit] == "dBm")
    txt(end + 1) = "The transducer table does not cover " + strjoin(unique([res([res.unit] == "dBm").bandName]), ", ") + ...
        ": left in dBm (no limit comparison there).";
end
if isempty(tr)
    txt(end + 1) = "No transducer given: levels are dBm at the analyzer input, so no limit comparison was made.";
elseif height(limT) > 0
    for i = 1:height(limT)
        s = sprintf("INDICATIVE (no quasi-peak detector), %s %s vs %s: ", limT.Band(i), limT.Detector(i), limT.Limit(i));
        if isnan(limT.InverterMargin_dB(i))
            s = s + "no inverter emission flagged";
        else
            s = s + sprintf("the flagged inverter emissions reach %s the limit (at %s)", aboveBelow(limT.InverterMargin_dB(i)), ...
                fmtHz(limT.InverterFreq_Hz(i)));
        end
        s = s + sprintf("; whole ON spectrum incl. ambient: %s at %s, %d points over", ...
            aboveBelow(limT.WorstMargin_dB(i)), fmtHz(limT.WorstFreq_Hz(i)), limT.PointsOverLimit(i));
        if limT.AmbientWorstMargin_dB(i) > 0
            s = s + sprintf(". With the inverter OFF the level (ambient signals or the analyzer noise) is itself up to %.1f dB over the limit: this setup cannot show compliance there", ...
                limT.AmbientWorstMargin_dB(i));
        end
        txt(end + 1) = s + "."; %#ok<AGROW>
    end
end
if height(rangeT) > 0
    v = rangeT(~isnan(rangeT.Range_FreeSpace_m), :);
    if height(v) == 0
        txt(end + 1) = "Interference range: not estimated (the emissions were measured in the near field or with a LISN).";
    else
        [~, j] = max(v.Range_FreeSpace_m);
        txt(end + 1) = sprintf("Interference range (rough): the emission at %s stays above the noise of a nearby receiver (I/N %g dB) up to ~%s in free space or ~%s with the close-in model.", ...
            fmtHz(v.Freq_Hz(j)), opt.INR, fmtDist(v.Range_FreeSpace_m(j)), fmtDist(v.Range_CloseIn_m(j)));
    end
end
txt(end + 1) = "Caveats: indicative pre-scan only (no QP/CISPR-average detector, 3 dB RBW, not a test site); below 30 MHz the near field and cable coupling dominate.";
txt = txt(:);
end

%% ===== figures ===============================================================

function files = makeFigures(S, res, Z, curves, opt)
% Each figure on its own: one that fails (or is closed while drawing) is skipped.
files = strings(0, 1);
jobs = {@() figOverview(res, opt)};
for b = 1:height(S.bands), jobs{end + 1} = @() figBand(S, res, b, opt); end %#ok<AGROW>
for i = find(arrayfun(@(r) r.comb.found, res)), jobs{end + 1} = @() figComb(res(i), opt); end %#ok<AGROW>
fz = unique([Z.freq], 'stable');
for i = 1:numel(fz), jobs{end + 1} = @() figZeroSpan(Z([Z.freq] == fz(i)), i, opt); end %#ok<AGROW>
if ~isempty(curves), jobs{end + 1} = @() figRange(curves, opt); end
for k = 1:numel(jobs)
    before = findall(groot, 'Type', 'figure');
    try
        f = jobs{k}();
        files = [files; f(:)]; %#ok<AGROW>
    catch err
        warning('ff_emi_study:figure', 'A figure was skipped: %s', err.message);
        if strcmpi(opt.Visible, 'off'), delete(setdiff(findall(groot, 'Type', 'figure'), before)); end
    end
end
end

function file = figOverview(res, opt)
% Floor change per band with CI
fig = newFig(opt, [80 80 900 420]);
labels = unique([res.bandName], 'stable'); dets = unique([res.det], 'stable');
Y = nan(numel(labels), numel(dets)); L = Y; H = Y;
for i = 1:numel(res)
    a = find(labels == res(i).bandName); b = find(dets == res(i).det);
    Y(a, b) = res(i).floorChange; L(a, b) = res(i).floorChange - res(i).floorCI(1); H(a, b) = res(i).floorCI(2) - res(i).floorChange;
end
ax = axes(fig);
hb = bar(ax, [Y; nan(size(Y))]);          % extra NaN row keeps the grouping with one band
hold(ax, 'on'); grid(ax, 'on');
for b = 1:numel(dets)
    errorbar(ax, hb(b).XEndPoints(1:numel(labels)), Y(:, b), L(:, b), H(:, b), 'k.', 'LineWidth', 1, ...
        'HandleVisibility', 'off');
end
set(ax, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'XLim', [0.5 numel(labels) + 0.5]); xtickangle(ax, 20);
ylabel(ax, 'ON - OFF noise floor (dB)');
legend(ax, detLabels(dets), 'Location', 'best');
title(ax, sprintf('Noise-floor change per band (%.0f%% bootstrap CI)', 100*(1 - opt.Alpha)));
file = saveFig(fig, fullfile(opt.OutDir, "fig_overview.png"), opt);
end

function file = figBand(S, res, b, opt)
% Spectra, difference with CI, waterfall
pd = primaryDet(S);
rb = res([res.band] == b);
if isempty(rb), file = strings(0, 1); return, end
r = rb(1); if any([rb.det] == pd), r = rb([rb.det] == pd); end
[sc, un] = freqUnit(r.f(end)); x = r.f/sc;
fig = newFig(opt, [60 40 1200 950]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact');
ax1 = nexttile(tl); hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
sty = {'-', '--'}; hl = gobjects(0); names = strings(0);
for j = 1:numel(rb)
    hl(end + 1) = plot(ax1, x, rb(j).offMean + rb(j).corr, sty{min(j, 2)}, 'Color', [0.2 0.4 0.9]); %#ok<AGROW>
    hl(end + 1) = plot(ax1, x, rb(j).onMean + rb(j).corr, sty{min(j, 2)}, 'Color', [0.85 0.2 0.1]); %#ok<AGROW>
    names = [names, "OFF " + rb(j).detLabel, "ON " + rb(j).detLabel]; %#ok<AGROW>
    if any(isfinite(rb(j).limit))
        hl(end + 1) = plot(ax1, x, rb(j).limit, sty{min(j, 2)}, 'Color', 'k', 'LineWidth', 1.3); %#ok<AGROW>
        names(end + 1) = "Limit " + rb(j).detLabel + " (indicative)"; %#ok<AGROW>
    end
end
if any(r.flag)
    hl(end + 1) = plot(ax1, x(r.flag), r.onMean(r.flag) + r.corr(r.flag), 'kx', 'MarkerSize', 6);
    names(end + 1) = "Flagged (" + r.detLabel + ")";
end
legend(ax1, hl, names, 'Location', 'best', 'FontSize', 7);
ylabel(ax1, "Level (" + r.unit + ")");
title(ax1, sprintf('%s   RBW %s   %s sweeps (OFF1/ON/OFF2)', r.bandName, fmtHz(r.rbw), strjoin(string(r.nSweeps), '/')));
ax2 = nexttile(tl); hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
ok = isfinite(r.ciLo);
fill(ax2, [x(ok) fliplr(x(ok))], [r.ciLo(ok) fliplr(r.ciHi(ok))], [0.75 0.8 1], 'EdgeColor', 'none');
plot(ax2, x, r.effect, 'Color', [0.1 0.2 0.6]);
plot(ax2, x(r.flag), r.effect(r.flag), 'r.', 'MarkerSize', 8);
if any(r.drift), plot(ax2, x(r.drift), r.effect(r.drift), 'o', 'Color', [0.9 0.6 0], 'MarkerSize', 4); end
yline(ax2, 0, 'k'); yline(ax2, opt.MinEffect, 'k--');
ylabel(ax2, 'ON - OFF (dB)');
title(ax2, sprintf('%s detector: effect with %.0f%% CI; red = flagged (q < %g, CI > %g dB), orange = OFF1/OFF2 drift', ...
    r.detLabel, 100*(1 - opt.Alpha), opt.Alpha, opt.MinEffect), 'FontWeight', 'normal', 'FontSize', 9);
ax3 = nexttile(tl);
W = []; edges = 0; lab = strings(0, 1);
for st = unique([S.seg.state], 'stable')
    g = getSeg(S, st, b, r.det);
    if isempty(g), continue, end
    W = [W; matchAxis(g, r.f) + r.corr]; %#ok<AGROW>
    edges(end + 1) = size(W, 1); lab(end + 1) = st; %#ok<AGROW>
end
imagesc(ax3, x([1 end]), [1 size(W, 1)], W);
set(ax3, 'YDir', 'normal'); colormap(ax3, turbo); cb = colorbar(ax3); cb.Label.String = r.unit;
for e = edges(2:end - 1), yline(ax3, e + 0.5, 'w-', 'LineWidth', 1.5); end
set(ax3, 'YTick', (edges(1:end - 1) + edges(2:end))/2 + 0.5, 'YTickLabel', lab);
xlabel(ax3, "Frequency (" + un + ")"); title(ax3, 'All sweeps in time order', 'FontWeight', 'normal');
linkaxes([ax1 ax2 ax3], 'x'); xlim(ax1, x([1 end]));
file = saveFig(fig, fullfile(opt.OutDir, sprintf("fig_band%d.png", b)), opt);
end

function file = figComb(r, opt)
% Harmonic comb: difference, cepstrum, line fit
c = r.comb;
fig = newFig(opt, [100 60 1100 800]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
plot(ax, r.f/1e6, r.effect, 'Color', [0.1 0.2 0.6]);
[~, loc] = ismember(c.fpk, r.f);
plot(ax, c.fpk/1e6, r.effect(loc), 'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
xlabel(ax, 'Frequency (MHz)'); ylabel(ax, 'ON - OFF (dB)');
title(ax, sprintf('%s (%s): %d harmonic lines of %.4g kHz', r.bandName, r.detLabel, c.nHarm, c.spacing/1e3));
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
ok = c.q > 0 & c.q <= 4*c.q0;
plot(ax, 1e6*c.q(ok), 10*log10(c.pq(ok)), 'Color', [0.3 0.3 0.3]);
xline(ax, 1e6*c.q0, 'r--', sprintf('%.4g us = 1/(%.4g kHz)', 1e6*c.q0, 1e-3/c.q0));
xlabel(ax, 'Quefrency (\mus)'); ylabel(ax, 'Cepstrum (dB)');
title(ax, 'Cepstrum of the ON-OFF difference (periodogram of the dB spectrum)', 'FontWeight', 'normal');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
plot(ax, c.m, c.fpk/1e6, 'k.');
plot(ax, c.m, c.m*c.spacing/1e6, 'r-');
xlabel(ax, 'Harmonic number'); ylabel(ax, 'Line frequency (MHz)');
title(ax, sprintf('fitlm: spacing %.5g kHz, %.0f%% CI %.5g to %.5g kHz, RMSE %.3g kHz', c.spacing/1e3, ...
    100*(1 - opt.Alpha), c.ci(1)/1e3, c.ci(2)/1e3, c.rmse/1e3), 'FontWeight', 'normal');
file = saveFig(fig, fullfile(opt.OutDir, sprintf("fig_comb_band%d_%s.png", r.band, r.det)), opt);
end

function file = figZeroSpan(zi, i, opt)
% Zero span: envelope and envelope spectrum, ON vs OFF
fig = newFig(opt, [120 80 1100 650]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');
ax1 = nexttile(tl); hold(ax1, 'on'); grid(ax1, 'on');
ax2 = nexttile(tl); hold(ax2, 'on'); grid(ax2, 'on');
col = containers.Map({'ON', 'OFF2', 'OFF1'}, {[0.85 0.2 0.1], [0.2 0.4 0.9], [0.2 0.4 0.9]});
for k = 1:numel(zi)
    c = col(char(zi(k).state));
    plot(ax1, 1e3*zi(k).t, zi(k).env, 'Color', c, 'DisplayName', "Inverter " + inverterOf(zi(k).state));
    plot(ax2, zi(k).fp, 10*log10(zi(k).ps), 'Color', c, 'DisplayName', "Inverter " + inverterOf(zi(k).state));
    if zi(k).state == "ON"
        yline(ax1, zi(k).lev, ':', 'Color', c, 'HandleVisibility', 'off');
    end
end
xlabel(ax1, 'Time (ms)'); ylabel(ax1, 'Level (dBm)'); legend(ax1, 'Location', 'best');
title(ax1, sprintf('Zero span at %s: first sweep (dotted = statelevels)', fmtHz(zi(1).freq)));
set(ax2, 'XScale', 'log'); xlabel(ax2, 'Envelope frequency (Hz)'); ylabel(ax2, 'Envelope spectrum (dB)');
xline(ax2, [50 100], ':', {'50 Hz', '100 Hz'}, 'HandleVisibility', 'off'); legend(ax2, 'Location', 'best');
title(ax2, 'Envelope spectrum (pspectrum, mean of the sweeps)', 'FontWeight', 'normal');
file = saveFig(fig, fullfile(opt.OutDir, sprintf("fig_zerospan_%d.png", i)), opt);
end

function file = figRange(curves, opt)
% Interference range
fig = newFig(opt, [140 100 900 500]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
co = lines(numel(curves));
for i = 1:numel(curves)
    semilogx(ax, curves(i).d, curves(i).mF, '-', 'Color', co(i, :), 'DisplayName', curves(i).label + " free space");
    semilogx(ax, curves(i).d, curves(i).mC, '--', 'Color', co(i, :), 'DisplayName', curves(i).label + " close-in");
end
set(ax, 'XScale', 'log'); yline(ax, 0, 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel(ax, 'Distance from the inverter (m)'); ylabel(ax, 'Level above victim threshold (dB)');
title(ax, sprintf('Interference range: threshold = kTB + NF %g dB + I/N %g dB', opt.VictimNF, opt.INR));
legend(ax, 'Location', 'eastoutside', 'FontSize', 7);
file = saveFig(fig, fullfile(opt.OutDir, "fig_range.png"), opt);
end

function fig = newFig(opt, pos)
fig = figure('Visible', opt.Visible, 'Color', 'w', 'Position', pos);
end

function file = saveFig(fig, file, opt)
exportgraphics(fig, file, 'Resolution', 130);
if strcmpi(opt.Visible, 'off'), close(fig); end
end

%% ===== report ================================================================

function file = makeReport(R, opt, files)
import mlreportgen.report.*
import mlreportgen.dom.*
S = R.study;
fprintf('Writing the %s report ...\n', upper(opt.Report));
rpt = Report(fullfile(opt.OutDir, "emi_report"), char(opt.Report));
when = datetime('now', 'TimeZone', 'Asia/Bangkok', 'Format', 'yyyy-MM-dd HH:mm ''(Asia/Bangkok)''');
sub = "FieldFox OFF-ON-OFF pre-scan";
if isSynthetic(S), sub = "SYNTHETIC DEMO DATA - not a measurement"; end
add(rpt, TitlePage('Title', 'Solar-inverter EMI study', 'Subtitle', sub, ...
    'Author', string(getenv('USERNAME')), 'PubDate', char(when)));
add(rpt, TableOfContents);
png = @(name) char(fullfile(opt.OutDir, name));

ch = Chapter('Title', 'Conclusions');
warn = Paragraph(['INDICATIVE PRE-SCAN ONLY: no quasi-peak or CISPR-average detector (option 361 is ' ...
    'not installed), 3 dB RBW instead of the CISPR 6 dB bandwidth, and not a calibrated test site.']);
warn.Bold = true; warn.Color = 'red';
add(ch, warn);
add(ch, Paragraph('These sentences are generated automatically from the statistics. Check them against the figures.'));
add(ch, UnorderedList(cellstr(R.conclusions)));
if isfile(png("fig_overview.png"))
    add(ch, FormalImage('Image', png("fig_overview.png"), 'Caption', 'Noise-floor change per band and detector', 'Width', '15cm'));
end
add(ch, rptTable(R.summary, 'Summary per band and detector', 20));
add(rpt, ch);

ch = Chapter('Setup and protocol');
add(ch, Paragraph(sprintf('Instrument: %s. Source: %s. Options: %s.', S.idn, S.source, strjoin(S.options, ','))));
if strlength(S.userNotes) > 0, add(ch, Paragraph("Notes: " + S.userNotes)); end
add(ch, Paragraph(sprintf(['Protocol: inverter OFF, ON%s. Every band was swept %d times per detector ' ...
    '(%s); the detectors are separate sweep sets because firmware A.12.45 has one detector for all traces. ' ...
    'Preamp %s, reference level %g dBm, attenuation %s.'], ...
    repmat(', OFF again', 1, S.protocolName == "ABA"), S.sweeps, strjoin(detLabels(S.detectors), ' and '), ...
    string(S.preamp), S.refLevel, string(S.atten))));
add(ch, rptTable(S.bands, 'Scan plan', 20));
if height(S.protocol) > 0, add(ch, rptTable(S.protocol, 'Protocol steps (local time)', 10)); end
if ~isempty(S.notes), add(ch, UnorderedList(cellstr(S.notes))); end
tr = "none (levels in dBm at the analyzer input)";
if isnumeric(opt.Transducer) && ~isempty(opt.Transducer)
    tr = opt.TransducerType + " factor table (" + size(opt.Transducer, 1) + " rows)";
elseif istable(opt.Transducer)
    tr = opt.TransducerType + " factor table (" + height(opt.Transducer) + " rows)";
elseif any(strlength(string(opt.Transducer)) > 0)
    tr = opt.TransducerType + " factor from " + string(opt.Transducer);
end
add(ch, Paragraph("Transducer: " + tr + "."));
add(rpt, ch);

ch = Chapter('Results per band');
for b = 1:height(S.bands)
    sec = Section('Title', char(S.bands.Name(b)));
    f = png(sprintf("fig_band%d.png", b));
    if isfile(f), add(sec, FormalImage('Image', f, 'Caption', ['Spectra OFF/ON, ON-OFF difference with ' ...
            'bootstrap CI, and all sweeps as a waterfall'], 'Width', '16cm')); end
    Eb = R.emissions([], :);
    if height(R.emissions) > 0, Eb = R.emissions(R.emissions.Band == S.bands.Name(b), :); end
    if height(Eb) > 0
        add(sec, Paragraph(sprintf('%d flagged emissions; the strongest are listed.', height(Eb))));
        add(sec, rptTable(emissionView(sortrows(Eb, 'Effect_dB', 'descend')), 'Flagged emissions', 15));
    else
        add(sec, Paragraph('No flagged emissions in this band.'));
    end
    add(ch, sec);
end
add(rpt, ch);

ch = Chapter('Switching-frequency harmonic comb');
add(ch, Paragraph(['Method: the ON-OFF difference (dB) is treated as a signal; its periodogram (the power ' ...
    'cepstrum) peaks at a quefrency of 1/spacing. Lines near each multiple of the spacing that are significant ' ...
    'and stand out of their neighbourhood are fitted with fitlm (no intercept): the slope is the spacing.']));
if height(R.comb) > 0
    add(ch, rptTable(R.comb, 'Detected combs', 10));
    for f = files(contains(files, "fig_comb_")).'
        add(ch, FormalImage('Image', f, 'Caption', 'Harmonic comb', 'Width', '15cm'));
    end
else
    add(ch, Paragraph('No harmonic comb was detected.'));
end
add(rpt, ch);

ch = Chapter('Zero-span time analysis');
if height(R.zeroSpan) > 0
    add(ch, Paragraph(['Zero span shows the level at one frequency against time. statelevels, pulsewidth, ' ...
        'pulseperiod and dutycycle are only used when the on/off contrast is above the 99th percentile that noise alone reaches (NoiseContrast_dB, simulated for the same RBW, points and sweep time) and the period repeats within 20 % between sweeps. The envelope can only follow ' ...
        'changes slower than about RBW/2.']));
    add(ch, rptTable(R.zeroSpan, 'Zero-span metrics', 20));
    for f = files(contains(files, "fig_zerospan_")).'
        add(ch, FormalImage('Image', f, 'Caption', 'Zero-span envelope and envelope spectrum', 'Width', '15cm'));
    end
else
    add(ch, Paragraph('No zero-span data (offline recordings, or no new emission found).'));
end
add(rpt, ch);

ch = Chapter('Indicative limit comparison');
p = Paragraph(['Indicative only: peak data compared with quasi-peak limits and RMS-average data with average ' ...
    'limits both read high; distance scaling 20 log10(d) is valid only in the far field.']);
p.Bold = true; add(ch, p);
if height(R.limits) > 0
    add(ch, rptTable(R.limits, ['Margins in dB, positive = above the limit. Inverter = worst flagged point; Worst = whole ON ' 'spectrum incl. ambient; Ambient = inverter OFF'], 20));
else
    add(ch, Paragraph('Not done: give a transducer (antenna factor or LISN CSV) to get dBuV/m or dBuV.'));
end
add(rpt, ch);

ch = Chapter('Interference range estimate');
add(ch, Paragraph(sprintf(['The inverter part of each emission (ON minus OFF power) is turned into an EIRP and ' ...
    'propagated with propagationModel("freespace") and ("close-in", n = 2.9). The range is where it falls below ' ...
    'kTB + %g dB noise figure + %g dB I/N in the measurement RBW. 0 = below the threshold already at 1 m; ' ...
    'Inf = beyond 100 km. Below 30 MHz the near field and cable coupling make this invalid.'], opt.VictimNF, opt.INR)));
if height(R.range) > 0
    add(ch, rptTable(R.range, 'Range estimate', 10));
    if isfile(png("fig_range.png")), add(ch, FormalImage('Image', png("fig_range.png"), 'Caption', 'Level above the victim threshold vs distance', 'Width', '15cm')); end
end
add(rpt, ch);

ch = Chapter('Methods and files');
add(ch, UnorderedList({ ...
    sprintf('Per-point test: %s (OFF1+OFF2 vs ON), Benjamini-Hochberg FDR at %g.', opt.Test, opt.Alpha), ...
    sprintf('Effect: mean over sweeps of ON - (OFF1+OFF2)/2 in dB; %d-resample percentile bootstrap CI (bootci).', opt.Boot), ...
    sprintf('Flag: q < %g and CI lower bound > %g dB. Band floor change: mean of the per-sweep median difference.', opt.Alpha, opt.MinEffect), ...
    'Spectra are power averages of the sweeps (mW, not dB).', ...
    'Source labels: ff_signal_id when available, otherwise only the band table; a match only says the frequency fits.'}));
add(ch, UnorderedList(cellstr(files(1:min(end, 40)))));
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
    if isstring(T.(v)), T.(v) = cellstr(T.(v)); end    % strings would be shown in quotes
end
mt = MATLABTable(T);
mt.Style = {FontSize('7pt')};
bt = BaseTable(mt);
bt.Title = ttl;
end

function V = emissionView(E)
V = table(E.Freq_Hz, E.Detector, E.Shape + " / " + E.Kind, compose("%.1f", E.Level_ON) + " " + E.Unit, ...
    compose("%+.1f [%.1f, %.1f]", E.Effect_dB, E.CI_low_dB, E.CI_high_dB), compose("%.2g", E.q), ...
    E.SignalID, E.LikelySource, 'VariableNames', {'Freq_Hz', 'Detector', 'Type', 'Level_ON', ...
    'Effect_dB_CI', 'q', 'SignalID', 'LikelySource'});
if ismember('Margin_dB', E.Properties.VariableNames) && any(isfinite(E.Margin_dB))
    V.Margin_dB = compose("%+.1f", E.Margin_dB);
end
end

%% ===== small helpers =========================================================

function B = defaultBands()
Name = ["Conducted 150k-30M"; "VHF 30-300M"; "UHF 300M-1G"; "L/S 1-3G"; "GNSS L1"];
Start_Hz = [150e3; 30e6; 300e6; 1e9; 1565e6];
Stop_Hz = [30e6; 300e6; 1e9; 3e9; 1585e6];
RBW_Hz = [9e3; 120e3; 120e3; 1e6; 100e3];
B = table(Name, Start_Hz, Stop_Hz, RBW_Hz);
B.Points = min(10001, max(401, ceil(2*(Stop_Hz - Start_Hz)./RBW_Hz) + 1));
end

function B = checkBands(B)
if ~istable(B) || ~all(ismember({'Name', 'Start_Hz', 'Stop_Hz', 'RBW_Hz'}, B.Properties.VariableNames))
    error('Bands must be a table with Name, Start_Hz, Stop_Hz, RBW_Hz (and optionally Points).');
end
B.Name = string(B.Name);
auto = min(10001, max(401, ceil(2*(B.Stop_Hz - B.Start_Hz)./B.RBW_Hz) + 1));
if ~ismember('Points', B.Properties.VariableNames), B.Points = auto; end
bad = isnan(B.Points) | B.Points < 2;
B.Points(bad) = auto(bad);
B.Points = min(10001, round(B.Points));
if any(B.Stop_Hz <= B.Start_Hz), error('Each band needs Stop_Hz > Start_Hz.'); end
if any(B.Start_Hz < 5e3 | B.Stop_Hz > 6.5e9), warning('The N9914A covers about 5 kHz (usable) to 6.5 GHz.'); end
B = B(:, {'Name', 'Start_Hz', 'Stop_Hz', 'RBW_Hz', 'Points'});
end

function s = protocolStates(p)
if p == "AB", s = ["OFF1" "ON"]; else, s = ["OFF1" "ON" "OFF2"]; end
end

function s = inverterOf(state)
s = "OFF"; if state == "ON", s = "ON"; end
end

function tf = isSynthetic(S)
tf = S.source == "demo" || contains(S.idn, "SYNTHETIC");
end

function d = primaryDet(S)
d = S.seg(1).det;
if any([S.seg.det] == "POS"), d = "POS"; end
end

function s = detLabel(det)
switch upper(string(det))
    case "POS", s = "Peak";
    case "AVER", s = "Average";
    case "SAMP", s = "Sample";
    case "NORM", s = "Normal";
    case "NEG", s = "Negative";
    otherwise, s = "Recorded";
end
end

function s = detLabels(dets)
s = string(arrayfun(@detLabel, dets, 'UniformOutput', false));
end

function c = detCode(det)
d = upper(string(det));
c = "REC";
for k = ["POS" "AVER" "SAMP" "NORM" "NEG"]
    if startsWith(d, k), c = k; end
end
end

function g = getSeg(S, state, b, det)
g = [];
for i = 1:numel(S.seg)
    if S.seg(i).state == state && S.seg(i).band == b && S.seg(i).det == det, g = S.seg(i); return, end
end
end

function P = matchAxis(seg, f)
P = seg.P;
if numel(seg.f) == numel(f) && max(abs(seg.f - f)) <= (f(2) - f(1))/2, return, end
warning('ff_emi_study:axis', 'Frequency axes differ (%s): interpolating.', seg.state);
P = interp1(seg.f, P.', f, 'linear', 'extrap').';
end

function y = powMean(P)
y = 10*log10(mean(10.^(P/10), 1));            % average in mW, not in dB
end

function q = fdrBH(p)
% Benjamini-Hochberg adjusted p-values (mafdr, Bioinformatics Toolbox, if present).
if exist('mafdr', 'file')
    q = reshape(mafdr(p(:), 'BHFDR', true), size(p));
else
    [ps, o] = sort(p(:)); m = numel(ps);
    qs = min(1, cummin(ps .* m ./ (1:m).', 'reverse'));
    q = zeros(size(p)); q(o) = qs;
end
end

function [s, e] = runs(mask, gap)
% Start/end indices of true runs, joining runs separated by <= gap false points.
mask = mask(:).';
d = diff([false mask false]);
s = find(d == 1); e = find(d == -1) - 1;
if isempty(s), return, end
keep = [true, s(2:end) - e(1:end - 1) - 1 > gap];
e = e([keep(2:end) true]); s = s(keep);
end

function tr = loadTransducer(x)
tr = [];
if isempty(x) || (~isnumeric(x) && all(strlength(string(x)) == 0)), return, end
if isnumeric(x), M = x; elseif istable(x), M = x{:, 1:2}; else, M = readmatrix(x); end
M = M(all(isfinite(M(:, 1:2)), 2), 1:2);
if size(M, 1) < 2, error('The transducer needs at least 2 rows of Freq_Hz, Factor_dB.'); end
if max(M(:, 1)) < 1e5
    M(:, 1) = M(:, 1)*1e6;
    fprintf('Transducer frequencies look like MHz: converted to Hz.\n');
end
M = sortrows(M, 1);
tr = struct('f', M(:, 1).', 'v', M(:, 2).');
end

function [c, unit] = transducerCorr(tr, f, opt)
c = zeros(size(f)); unit = "dBm";
if isempty(tr), return, end
v = interp1(log10(tr.f), tr.v, log10(f), 'linear');   % NaN outside the table
if mean(isnan(v)) > 0.5, return, end                  % table does not cover this band: stay in dBm
c = 107 + v;                                          % dBm -> dBuV (50 ohm), + factor
unit = "dBuV/m";
if opt.TransducerType == "lisn", unit = "dBuV"; end
end

function name = serviceBand(fHz, tol)
% Common services in Thailand. A match only says the frequency fits a band.
persistent bands
if isempty(bands)
    bands = {
        'AM broadcast (MW)',         526.5e3, 1606.5e3
        'FM radio',                   87.5e6,   108e6
        'Airband VHF',                 118e6,   137e6
        'Amateur 2 m',                 144e6,   148e6
        'Marine VHF',                  156e6,   162.1e6
        'UHF land mobile',             380e6,   470e6
        'ISM 433 MHz',              433.05e6, 434.79e6
        'Digital TV (UHF)',            470e6,   694e6
        'LTE/5G 700 (B28/n28)',        703e6,   803e6
        'GSM/LTE 850 (B5)',            824e6,   894e6
        'GSM/LTE 900 (B8)',            880e6,   960e6
        'IoT/RFID 920-925 MHz',        920e6,   925e6
        'GNSS L1 (GPS/Galileo/BDS)',  1559e6,  1610e6
        'LTE 1800 (B3)',              1710e6,  1880e6
        '3G/LTE 2100 (B1)',           1920e6,  2170e6
        'LTE TDD 2300 (B40)',         2300e6,  2400e6
        'Wi-Fi/Bluetooth 2.4 GHz',    2400e6, 2483.5e6
        '5G NR 2600 (n41)',           2496e6,  2690e6
        '5G NR 3500 (n78)',           3300e6,  3800e6
        'Wi-Fi 5 GHz',                5150e6,  5850e6 };
end
hit = fHz >= [bands{:, 2}] - tol & fHz <= [bands{:, 3}] + tol;
name = string(strjoin(bands(hit, 1), ' / '));
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

function s = aboveBelow(m)
if m > 0, s = sprintf("%.1f dB above", m); else, s = sprintf("%.1f dB below", -m); end
end

function s = fmtDist(d)
if isnan(d), s = '-';
elseif d == 0, s = '< 1 m';
elseif isinf(d), s = '> 100 km';
elseif d >= 1000, s = sprintf('%.3g km', d/1000);
else, s = sprintf('%.3g m', d);
end
end

function s = stamp()
s = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
end
