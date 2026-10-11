function R = ff_srs4_analyzer(ip, varargin)
%FF_SRS4_ANALYZER  Exam-grade analysis of an SRS-4 (RRC-PSK CubeSat) S-band downlink with a FieldFox.
%
%   R = ff_srs4_analyzer                                  measure (FieldFox at 10.90.36.200)
%   R = ff_srs4_analyzer(ip, Name, Value, ...)
%   R = ff_srs4_analyzer([], 'Demo', true)                synthetic SRS-4, no instrument needed
%   R = ff_srs4_analyzer([], 'File', "...\srs4_run.mat")  re-analyse a saved run (new options allowed)
%   R = ff_srs4_analyzer([], 'File', "rec.mat")           a fieldfox_live / shared-format recording
%                                                         of the downlink (spectrum analysis only)
%
% Example (over the air, antennas 3 m apart, no dialogs):
%   R = ff_srs4_analyzer("10.90.36.200", 'Distance', 3, 'TxAntennaGain', 6.1, ...
%         'RxAntennaGain', "patch", 'Prompt', false, 'Visible', 'off');
%   R.limits                      % PASS / FAIL / INFO table against the SRS-4 datasheet
%
% The radio: Satlab SRS-4, TX 2200-2290 MHz, BPSK/QPSK/8PSK, root-raised-cosine (RRC)
% roll-off 0.2/0.25/0.35, 10 kBd-5 MBd, 20-33 dBm, CCSDS framing (ASM, RS(255,223),
% randomizer, convolutional code r=1/2 K=7). Default here: 1 MBd (sl-radio --tx-rate=1000000),
% one QPSK RS+CC frame = 4144 coded bits = 2072 symbols = 2.072 ms.
%
% What the FieldFox does (spectrum analyzer, Port 2 input; it never transmits):
%   1. Band search: the strongest wideband emission in 'Band' (coarse sweep)
%   2. Zoom: span 3 x Rs, RBW about 1 % of Rs (1-3-10 steps), average detector,
%      'Averages' single sweeps power-averaged in MATLAB (they are kept for the bootstrap)
%   3. Mask sweep: wide enough for the SFCG/ECSS mask (+-3 Rs or +-8 Rs)
%   4. Its own channel power (CHP) and 99 % / 99.9 % OBW (CALC:MEAS:DATA?) as a cross-check
%   5. 2nd-harmonic window 2 fc +- 3 Rs (attenuation set to balance the analyzer's own 2nd
%      harmonic against its noise; if something shows, it is re-measured with up to +10 dB:
%      an analyzer-made harmonic drops 1 dB per dB, a real one stays) and a coarse
%      30 MHz - 6.5 GHz spurious scan (peak detector, max hold)
%   6. 'LineScan' (optional): RBW 10-100 Hz near fc for idle-frame spectral lines
%   7. Zero span at fc: envelope samples (RBW 5 MHz, sample detector, 'EnvelopeTime')
%      and a long sweep for burst on/off timing ('BurstTime'); the TX-off level between
%      bursts measures the zero-span noise that sets the envelope SNR
% Safe start: preamp OFF, attenuation AUTO, overload flag (CALC:MEAS:WAOR?) checked on the
% first sweep of every set; on overload the attenuation goes up 10 dB (to 30 dB max).
% The settings are backed up first and restored at the end, also after an error or Ctrl+C.
%
% What MATLAB computes:
%   - Rs = -3 dB bandwidth. The RRC transmit filter makes the power spectrum a raised
%     cosine (RC), which is exactly 3 dB down at fc +- Rs/2 for ANY roll-off.
%   - OBW 99 % / 99.9 % (obw), -26 dB emission bandwidth, roll-off from OBW99/Rs (RC table)
%   - least-squares fit (lsqcurvefit) of RC (x) Gaussian RBW + noise: fc, Rs, alpha, level,
%     noise, with 95 % confidence intervals from nlparci (scaled for RBW-correlated bins)
%     and from a bootstrap over the sweeps (the more honest one; t-scaled spread)
%   - reference: PSD of ccsdsTMWaveformGenerator (QPSK, RS+CC, roll-off 0.2/0.25/0.35)
%     through the same RBW filter; best roll-off = smallest residual RMS
%   - frequency error vs 'ExpectedFreq' in Hz and ppm, with the FieldFox reference error
%   - channel power (bandpower and the FieldFox CHP), C/N0, Es/N0, Eb/N0 for BPSK/QPSK/8PSK
%   - ECSS-E-ST-50-05C Table 5-5 OBW limit and the SFCG/ECSS mask (Fig. 6-5): margin
%   - 2nd harmonic and spurious in dBc against the 60 dBc SFDR spec, with detection limits
%   - envelope statistics (time below -10 dB, std of the dB envelope, PAPR) compared with
%     BPSK/QPSK/8PSK simulated by ccsdsTMWaveformGenerator at the measured SNR and RBW
%   - burst length (statelevels, pulsewidth) vs frame time for every modulation and FEC
%   - Over the air: EIRP and TX power from Friis (fspl), antenna gains given in dBi or
%     designed with Antenna Toolbox ("patch" = patchMicrostrip, "dipole")
% Results go to 'OutDir': srs4_report.pdf, figures (PNG), srs4_limits.csv, srs4_run.mat
% (raw data + results; re-use it with 'File') and rec\rec_<step>.mat (shared format:
% f, P, t, settings, idn - as displayed by the FieldFox).
%
% Name-Value options (defaults):
%  Signal
%   'Band'            [2200e6 2290e6]  Hz, search range (SRS-4 TX band)
%   'ExpectedRate'    1e6     Bd (sl-radio --tx-rate); [] = estimate from the band search
%   'ExpectedFreq'    []      Hz, nominal tx.freq; [] = unknown (frequency error not judged)
%   'ExpectedRolloff' []      0.2, 0.25 or 0.35 if you know it
%  Levels and over-the-air geometry
%   'ExtGain'         0       dB between DUT and FieldFox, negative = loss (pad + cable),
%                             set on the instrument (POW:EXTG); e.g. -40.6 for a 40 dB pad
%   'Distance'        []      m, TX antenna to FieldFox antenna (over the air); [] = conducted
%   'TxAntennaGain'   0       dBi, or "patch" / "dipole" (Antenna Toolbox design at fc; the
%                             designed air-substrate patch has ~10 dBi, a CubeSat patch such
%                             as the Satlab SAS-2 has 6.1 dBi: give the number if you know it)
%   'RxAntennaGain'   0       dBi, or "patch" / "dipole"
%   'PolLoss'         0       dB polarisation mismatch (e.g. 3 for RHCP patch -> linear antenna)
%  Measurement
%   'Averages'        20      zoom sweeps (power averaged; >= 3 for the bootstrap)
%   'Points'          "auto"  zoom trace points (101-10001; auto = 2 per RBW, 1001-10001)
%   'Mask' 'CrossCheck' 'Harmonics' 'WideScan' 'ZeroSpan'   true    steps on/off
%   'LineScan'        false   idle-frame line scan (slow: RBW 10-100 Hz)
%   'MaskSweeps' 10, 'HarmonicSweeps' 10, 'WideSweeps' 3, 'WideRBW' 1e6
%   'EnvelopeSweeps' 10, 'EnvelopeTime' 2e-3 s (nudged so the 10001 samples do not strobe
%   one symbol phase), 'BurstSweeps' 3, 'BurstTime' 0.5 s, 'LineSweeps' 3
%   'Port'            5025
%  Analysis
%   'FitRange'        15      dB below the top used by the RC fit, symmetric about fc (the
%                             skirts are excluded: truncated RRC filters and PA regrowth bias
%                             a deeper fit - 25 dB gave Rs -2 % and alpha +0.05 in simulation)
%   'Boot'            200     bootstrap resamples (0 = off)
%   'RefPpm'          []      FieldFox reference accuracy; [] = 0.7 ppm (internal, spec)
%                             or 0.05 ppm when ROSC:SOUR? reports EXT (GPSDO assumed)
%  Output
%   'Report'          "pdf"   "pdf", "docx" or "none"
%   'Notes'           ""      text for the report (setup, antennas, radio settings ...)
%   'OutDir'          ""      default srs4_<time> in the current folder
%   'Visible'         'on'    'off' = no windows (figures are still saved)
%   'Prompt'          true    false = no safety question before measuring
%   'Demo'            false   synthetic SRS-4 (ccsdsTMWaveformGenerator + FieldFox model)
%   'DemoModulation'  "QPSK"  "BPSK" / "QPSK" / "8PSK" for the demo
%   'File'            ""      re-analyse a saved run or recording
%
% Read this before trusting the numbers:
%   - NO I/Q CAPTURE (option 351 is not installed): constellation, EVM and the modulation
%     order cannot be measured. A power spectrum is identical for BPSK, QPSK and 8PSK with
%     the same RRC filter; only the zero-span envelope statistics hint at BPSK vs QPSK/8PSK,
%     and QPSK vs 8PSK is barely separable. The FieldFox zero-span sample rate is not
%     specified, so even that inference is indirect.
%   - Frequency: the PSK carrier is suppressed, so fc comes from the spectral fit, not a
%     counter. The FieldFox internal reference is +-0.7 ppm (spec, +-0.4 typ, +1 ppm/year
%     ageing): no better than the SRS-4 spec of 0.5 ppm. Lock it to an external 10 MHz
%     GPSDO (ROSC:SOUR EXT) before you claim anything at the ppm level. Statistical
%     precision of fc with 20 sweeps was about +-0.4 ppm (95 %) at 1 MBd in simulation and
%     grows with Rs: raise 'Averages' for tighter intervals.
%   - Over the air: antenna gains, polarisation and multipath make EIRP good to +-3 dB at
%     best; the analyzer noise floor limits the mask/spurious dynamic range. Compliance with
%     the 60 dBc SFDR and the -60 dB mask points needs a conducted test (30-40 dB pad, first
%     pad rated >= 10 W) or a filter in front of the FieldFox.
%   - The FieldFox makes its own 2nd harmonic (-60 dBc at -30 dBm mixer level, SHI +30 dBm);
%     the 3rd harmonic (6.6-6.87 GHz) is above the 6.5 GHz range of the N9914A.
%   - The RC fit assumes an ideal raised cosine; with a 10-symbol truncated RRC and a PA at
%     5 dB back-off the 15 dB fit was biased by alpha +0.006, Rs -0.5 %, fc +70 Hz in
%     simulation (noise-free). Treat roll-off claims as +-0.03.
%   - SFCG mask breakpoints were read from ECSS-E-ST-50-05C Rev.2 Fig. 6-5 (2011): check
%     them against SFCG Rec. 21-2R5 before quoting.
%   - SAFETY: the FieldFox input is damaged above +27 dBm (0.5 W) or 50 V DC; the SRS-4
%     delivers up to +33 dBm (2 W). Never connect its TX port without a power pad.

if nargin < 1 || isempty(ip) || all(strlength(string(ip)) == 0), ip = "10.90.36.200"; end
opt = struct('Band', [2200e6 2290e6], 'ExpectedRate', 1e6, 'ExpectedFreq', [], 'ExpectedRolloff', [], ...
    'ExtGain', 0, 'Distance', [], 'TxAntennaGain', 0, 'RxAntennaGain', 0, 'PolLoss', 0, ...
    'Averages', 20, 'Points', "auto", 'Mask', true, 'CrossCheck', true, 'Harmonics', true, ...
    'WideScan', true, 'ZeroSpan', true, 'LineScan', false, 'MaskSweeps', 10, 'HarmonicSweeps', 10, ...
    'WideSweeps', 3, 'WideRBW', 1e6, 'EnvelopeSweeps', 10, 'EnvelopeTime', 2e-3, 'BurstSweeps', 3, ...
    'BurstTime', 0.5, 'LineSweeps', 3, 'Port', 5025, 'FitRange', 15, 'Boot', 200, 'RefPpm', [], ...
    'Report', "pdf", 'Notes', "", 'OutDir', "", 'Visible', 'on', 'Prompt', true, 'Demo', false, ...
    'DemoModulation', "QPSK", 'File', "");
if mod(numel(varargin), 2), error('ff_srs4_analyzer:args', 'Options must be Name, Value pairs.'); end
names = fieldnames(opt);
for k = 1:2:numel(varargin)
    j = find(strcmpi(names, varargin{k}), 1);
    if isempty(j), error('ff_srs4_analyzer:args', 'Unknown option "%s".', string(varargin{k})); end
    opt.(names{j}) = varargin{k + 1};
end
opt = checkOptions(opt);

if strlength(opt.OutDir) == 0, opt.OutDir = fullfile(pwd, "srs4_" + stamp()); end
if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
rawFile = fullfile(opt.OutDir, "srs4_run.mat");

%% 1. Data: measure, synthesize or load
if opt.Demo
    D = demoRun(opt);
elseif strlength(opt.File) > 0
    D = loadRun(opt.File, opt);
else
    if ~safetyGate(opt)
        fprintf('Stopped at the safety question: nothing was measured.\n');
        R = struct('outDir', opt.OutDir, 'files', strings(0, 1));
        if nargout == 0, clear R, end
        return
    end
    D = measureRun(ip, opt);
end
run_data = D;
save(rawFile, 'run_data');                % raw data first: re-analyse with 'File' if anything below fails
recFiles = saveRecordings(D, opt.OutDir);

%% 2. Analysis
fprintf('Analysing ...\n');
A = analyzeRun(D, opt);
results = A; limits = A.limits; conclusionsText = A.conclusions;
save(rawFile, 'run_data', 'results', 'limits', 'conclusionsText');
csvFile = fullfile(opt.OutDir, "srs4_limits.csv");
writetable(A.limits, csvFile);

%% 3. Figures and report
figs = makeFigures(D, A, opt);
R.data = D; R.results = A; R.limits = A.limits; R.ebno = A.power.ebno; R.frames = A.frames;
R.conclusions = A.conclusions; R.outDir = opt.OutDir;
R.files = [rawFile; csvFile; figs; recFiles];
if opt.Report ~= "none"
    try
        rep = makeReport(R, opt, figs);
        R.files = [rep; R.files];
    catch err
        warning('ff_srs4_analyzer:report', 'Report not written: %s', err.message);
    end
end
fprintf('\n===== SRS-4 downlink vs datasheet (automatic - check the figures) =====\n');
disp(A.limits(:, {'Parameter', 'Measured', 'Limit', 'Status'}));
fprintf('- %s\n', A.conclusions);
fprintf('\nResults in %s\n', opt.OutDir);
if nargout == 0, clear R, end
end

%% ===== options and safety ===================================================

function opt = checkOptions(opt)
opt.Visible = char(matlab.lang.OnOffSwitchState(opt.Visible));
opt.Prompt = logical(opt.Prompt); opt.Demo = logical(opt.Demo);
for nm = ["Mask" "CrossCheck" "Harmonics" "WideScan" "ZeroSpan" "LineScan"]
    opt.(nm) = logical(opt.(nm));
end
opt.Report = lower(string(opt.Report)); opt.File = string(opt.File); opt.OutDir = string(opt.OutDir);
opt.Notes = string(opt.Notes); opt.DemoModulation = upper(string(opt.DemoModulation));
if ~ismember(opt.Report, ["pdf" "docx" "none"]), error('ff_srs4_analyzer:args', 'Report must be "pdf", "docx" or "none".'); end
if ~ismember(opt.DemoModulation, ["BPSK" "QPSK" "8PSK"]), error('ff_srs4_analyzer:args', 'DemoModulation must be BPSK, QPSK or 8PSK.'); end
if numel(opt.Band) ~= 2 || opt.Band(2) <= opt.Band(1) || opt.Band(1) < 2e6 || opt.Band(2) > 6.5e9
    error('ff_srs4_analyzer:args', 'Band must be [start stop] in Hz inside 2 MHz - 6.5 GHz.');
end
if ~isempty(opt.ExpectedRate) && ~(isscalar(opt.ExpectedRate) && opt.ExpectedRate >= 1e3 && opt.ExpectedRate <= 20e6)
    error('ff_srs4_analyzer:args', 'ExpectedRate must be a symbol rate in Bd (1 kBd - 20 MBd) or [].');
end
if ~isempty(opt.ExpectedRolloff) && ~(isscalar(opt.ExpectedRolloff) && opt.ExpectedRolloff > 0 && opt.ExpectedRolloff <= 1)
    error('ff_srs4_analyzer:args', 'ExpectedRolloff must be in (0, 1] or [].');
end
if ~isempty(opt.ExpectedFreq) && (opt.ExpectedFreq < opt.Band(1) || opt.ExpectedFreq > opt.Band(2))
    warning('ff_srs4_analyzer:args', 'ExpectedFreq %.6g MHz is outside Band.', opt.ExpectedFreq/1e6);
end
if ~isempty(opt.Distance) && ~(isscalar(opt.Distance) && opt.Distance > 0)
    error('ff_srs4_analyzer:args', 'Distance must be a positive number of metres or [].');
end
for nm = ["TxAntennaGain" "RxAntennaGain"]
    g = opt.(nm);
    if ~(isnumeric(g) && isscalar(g)) && ~any(strcmpi(string(g), ["patch" "dipole"]))
        error('ff_srs4_analyzer:args', '%s must be a gain in dBi, "patch" or "dipole".', nm);
    end
end
if ~(isnumeric(opt.ExtGain) && isscalar(opt.ExtGain) && abs(opt.ExtGain) <= 100)
    error('ff_srs4_analyzer:args', 'ExtGain must be in dB (-100 to 100; negative = loss).');
end
if isstring(opt.Points) || ischar(opt.Points)
    if ~strcmpi(opt.Points, "auto"), error('ff_srs4_analyzer:args', 'Points must be "auto" or 101-10001.'); end
    opt.Points = "auto";
elseif ~(isscalar(opt.Points) && opt.Points >= 101 && opt.Points <= 10001)
    error('ff_srs4_analyzer:args', 'Points must be "auto" or 101-10001.');
else
    opt.Points = round(opt.Points);
end
opt.Averages = max(1, round(opt.Averages));
opt.Boot = max(0, round(opt.Boot));
for nm = ["MaskSweeps" "HarmonicSweeps" "WideSweeps" "EnvelopeSweeps" "BurstSweeps" "LineSweeps"]
    opt.(nm) = max(1, round(opt.(nm)));
end
end

function ok = safetyGate(opt)
msg = sprintf(['SAFETY: the FieldFox SA input (Port 2) is damaged above +27 dBm (0.5 W) or 50 V DC, ' ...
    'and the SRS-4 transmits up to +33 dBm (2 W). Conducted: use 30-40 dB of attenuation (first pad ' ...
    'rated >= 10 W) and enter it as ExtGain (now %g dB). Over the air: keep the antennas at least 1 m ' ...
    'apart. Never connect the SRS-4 TX port directly.'], opt.ExtGain);
fprintf('\n%s\n', msg);
if isempty(opt.Distance) && opt.ExtGain > -20
    fprintf(['NOTE: no Distance and ExtGain > -20 dB: if this is a conducted test, the pad is missing ' ...
        'from ExtGain (levels will be wrong) or missing altogether (DANGER).\n']);
end
ok = true;
if opt.Prompt
    if usejava('desktop')
        ok = strcmp(questdlg(msg, 'ff_srs4_analyzer', 'Measure', 'Abort', 'Abort'), 'Measure');
    else
        try
            ok = strcmpi(strtrim(input('Type "yes" to measure: ', 's')), 'yes');
        catch
            ok = false;
            fprintf('No console for the safety question: run with ''Prompt'', false to measure headless.\n');
        end
    end
end
end

%% ===== measurement flow (shared by the instrument and the demo) =============

function D = newRun(source, opt)
D = struct('source', source, 'idn', "", 'options', strings(1, 0), 'time', datetime('now'), ...
    'notes', strings(0, 1), 'extGain', opt.ExtGain, 'extgInstr', 0, 'refSource', "unknown", ...
    'acq', struct(), 'inst', struct('done', false), 'search', [], 'plan', struct(), ...
    'aborted', false, 'truth', [], 'opt', opt);
end

function D = runFlow(E, D, opt)
% E.sweep(cfg, n), E.zero(cfg, n), E.cross(settings, Rs, n), E.setAtt(att), E.getAtt()
try
    D = runFlowSteps(E, D, opt);
catch err
    if ~isfield(D.acq, 'coarse') || strcmp(err.identifier, 'ff_srs4_analyzer:level'), rethrow(err), end
    D.aborted = true;
    D.notes(end + 1, 1) = "Measurement stopped: " + err.message;
    warning('ff_srs4_analyzer:measure', 'Measurement stopped (%s): analysing what was measured.', err.message);
end
end

function D = runFlowSteps(E, D, opt)
fprintf('1/7 Band search %s ...\n', fmtRange(opt.Band));
D.acq.coarse = E.sweep(planCoarse(opt, D), 3);
Sx = searchEmission(D.acq.coarse, D, opt);
D.search = Sx;
D = saveStep(D, opt);
if Sx.PinInput > 20
    error('ff_srs4_analyzer:level', ['%.1f dBm at the FieldFox input: too close to the +27 dBm damage level. ' ...
        'Add attenuation before measuring.'], Sx.PinInput);
end
if ~Sx.found
    D.notes(end + 1, 1) = "No emission above the noise floor in the band: is the SRS-4 transmitting?";
    warning('ff_srs4_analyzer:nosignal', '%s', D.notes(end));
    return
end
fprintf('   strongest emission: %s, -3 dB width %s, %.1f dBm%s\n', fmtHz(Sx.fc), fmtHz(Sx.bw3), Sx.Pch, ...
    repmat(' (NARROWBAND: CW-like)', 1, ~Sx.wide));
Rs0 = opt.ExpectedRate;
if isempty(Rs0)
    Rs0 = max(Sx.bw3, 10e3);
    D.notes(end + 1, 1) = "ExpectedRate not given: zoom based on the band-search width " + fmtHz(Rs0) + ".";
end
fc0 = Sx.fc;
ref = 5*ceil((toDisplay(Sx.Pch, D) + 5)/5);         % reference level ~5 dB above the channel power

fprintf('2/7 Zoom: %d sweeps ...\n', opt.Averages);
z = E.sweep(planZoom(fc0, Rs0, ref, opt), opt.Averages);
D.acq.zoom = z;
D.plan = struct('Rs0', Rs0, 'fc0', fc0, 'ref', ref);
q = quickEdges(z, D, Rs0);
if isfinite(q.Rs3) && (isempty(opt.ExpectedRate) || abs(q.Rs3/Rs0 - 1) > 0.2)
    if ~isempty(opt.ExpectedRate)
        D.notes(end + 1, 1) = sprintf("-3 dB width %s differs from ExpectedRate %s by more than 20 %%: zoomed again.", ...
            fmtHz(q.Rs3), fmtHz(opt.ExpectedRate));
        warning('ff_srs4_analyzer:rate', '%s', D.notes(end));
    end
    Rs0 = q.Rs3; fc0 = q.fc3;
    fprintf('   re-zoom with Rs ~ %s\n', fmtHz(Rs0));
    z = E.sweep(planZoom(fc0, Rs0, ref, opt), opt.Averages);
    q = quickEdges(z, D, Rs0);
end
if isfinite(q.fc3), fc0 = q.fc3; end
D.acq.zoom = z;
D.plan = struct('Rs0', Rs0, 'fc0', fc0, 'ref', ref);
D = saveStep(D, opt);

steps = {
    @() crossStep(), opt.CrossCheck, '3/7 FieldFox channel power and OBW'
    @() maskStep(), opt.Mask, '4/7 Mask sweep'
    @() harmStep(), opt.Harmonics, '5/7 2nd harmonic'
    @() wideStep(), opt.WideScan, '5/7 Spurious scan 30 MHz - 6.5 GHz'
    @() lineStep(), opt.LineScan, '6/7 Idle-frame line scan'
    @() zeroStep(), opt.ZeroSpan, '7/7 Zero span (envelope and bursts)'};
for k = 1:size(steps, 1)
    if ~steps{k, 2}, continue, end
    fprintf('%s ...\n', steps{k, 3});
    try
        steps{k, 1}();
    catch err
        if strcmp(err.identifier, 'FieldFox:timeout') || contains(err.message, 'connection', 'IgnoreCase', true)
            D.notes(end + 1, 1) = string(steps{k, 3}) + " failed: " + err.message + " (later steps skipped)";
            warning('ff_srs4_analyzer:step', '%s', D.notes(end));
            D.aborted = true;
            break
        end
        D.notes(end + 1, 1) = string(steps{k, 3}) + " failed: " + err.message;
        warning('ff_srs4_analyzer:step', '%s', D.notes(end));
    end
    D = saveStep(D, opt);
end

    function crossStep()
        D.inst = E.cross(z.settings, Rs0, min(opt.Averages, 10));
        if strlength(D.inst.note) > 0, D.notes(end + 1, 1) = "FieldFox CHP/OBW: " + D.inst.note; end
    end
    function maskStep()
        D.acq.mask = E.sweep(planMask(fc0, Rs0, ref, opt), opt.MaskSweeps);
    end
    function harmStep()
        if 2*fc0 + 3*Rs0 > 6.5e9
            D.notes(end + 1, 1) = "2nd harmonic above 6.5 GHz: not measured."; return
        end
        att0 = E.getAtt();
        restore = onCleanup(@() E.setAtt(att0));
        % The analyzer makes its own 2nd harmonic (mixer level - SHI, SHI +30 dBm) and attenuation
        % raises its noise: about Pin + 27 dB of attenuation balances the two (DANL -135 dBm/Hz).
        attH = min(30, max([0, z.settings.atten, 5*round((Sx.PinInput + 27)/5)]));
        E.setAtt(attH);
        D.acq.harm = E.sweep(planHarm(fc0, Rs0, ref, opt), opt.HarmonicSweeps);
        if quickExcess(D.acq.harm, 2*fc0, Rs0) > 1 && attH < 30
            fprintf('   something at 2 fc: repeating with more attenuation (analyzer distortion test)\n');
            E.setAtt(min(30, attH + 10));
            D.acq.harm2 = E.sweep(planHarm(fc0, Rs0, ref, opt), opt.HarmonicSweeps);
        end
    end
    function wideStep()
        D.acq.wide = E.sweep(planWide(ref, opt), opt.WideSweeps);
    end
    function lineStep()
        D.acq.line = E.sweep(planLine(fc0, Rs0, ref, opt), opt.LineSweeps);
    end
    function zeroStep()
        D.acq.env = E.zero(planZero(fc0, opt.EnvelopeTime, "SAMP", ref, Rs0), opt.EnvelopeSweeps);
        D.acq.burst = E.zero(planZero(fc0, opt.BurstTime, "AVER", ref, Rs0), opt.BurstSweeps);
    end
end

function D = saveStep(D, opt)
% Raw data after every step: an interrupted run can still be analysed with 'File'.
run_data = D;
try
    save(fullfile(opt.OutDir, "srs4_run.mat"), 'run_data');
catch
end
end

function c = planCoarse(opt, D)
span = diff(opt.Band);
rbw = 100e3;
if ~isempty(opt.ExpectedRate), rbw = snapRBW(min(100e3, max(opt.ExpectedRate/10, 1e3))); end
c = struct('label', "coarse", 'start', opt.Band(1), 'stop', opt.Band(2), 'rbw', rbw, ...
    'points', autoPoints(span, rbw, 1001, 10001), 'det', "AVER", 'ref', 10 - D.extgInstr, 'timeout', 120);
end

function c = planZoom(fc, Rs, ref, opt)
span = 3*Rs;
rbw = snapRBW(0.01*Rs);
n = opt.Points;
if ~isnumeric(n), n = autoPoints(span, rbw, 1001, 10001); end
c = struct('label', "zoom", 'center', fc, 'span', span, 'rbw', rbw, 'vbw', min(3*rbw, 5e6), ...
    'points', n, 'det', "AVER", 'ref', ref, 'timeout', 120);
end

function c = planMask(fc, Rs, ref, opt) %#ok<INUSD>
xEnd = 8; if Rs >= 2e6, xEnd = 3; end            % SFCG/ECSS curve ends (offset / Rs)
span = min(2*1.1*xEnd*Rs, 400e6);
rbw = snapRBW(max(0.01*Rs, span/4500));
c = struct('label', "mask", 'center', fc, 'span', span, 'rbw', rbw, 'vbw', min(3*rbw, 5e6), ...
    'points', autoPoints(span, rbw, 1001, 10001), 'det', "AVER", 'ref', ref, 'timeout', 120);
end

function c = planHarm(fc, Rs, ref, opt) %#ok<INUSD>
span = 6*Rs;
rbw = snapRBW(0.02*Rs);
c = struct('label', "harm", 'center', 2*fc, 'span', span, 'rbw', rbw, 'vbw', min(3*rbw, 5e6), ...
    'points', autoPoints(span, rbw, 1001, 4001), 'det', "AVER", 'ref', ref - 20, 'timeout', 120);
end

function c = planWide(ref, opt)
c = struct('label', "wide", 'start', 30e6, 'stop', 6.5e9, 'rbw', snapRBW(opt.WideRBW), ...
    'points', 10001, 'det', "POS", 'ref', ref, 'timeout', 180);
end

function c = planLine(fc, Rs, ref, opt) %#ok<INUSD>
% Idle frames repeat every T_frame, so their spectrum is a comb with spacing 1/T_frame:
% 241 Hz (BPSK RS+CC) ... 1652 Hz (8PSK uncoded) at 1 MBd.
bits = [1816 2072 3632 4144];
spMin = Rs/max(bits); spMax = 3*Rs/min(bits);
rbw = min(100, max(10, snapRBW(spMin/5)));
span = min(200e3, max(5e3, 25*spMax));
c = struct('label', "line", 'center', fc + 0.2*Rs, 'span', span, 'rbw', rbw, 'vbw', 3*rbw, ...
    'points', autoPoints(span, rbw, 1001, 10001), 'det', "AVER", 'ref', ref - 10, 'timeout', 600);
end

function c = planZero(fc, T, det, ref, Rs)
% Envelope samples must not strobe at a fixed symbol phase: with a sample spacing of exactly
% k symbols (0.2 us at 5 MBd!) every sample sees the same point of the pulse. The sweep time
% is nudged so that the spacing in symbols has a golden-ratio fraction (phases fill evenly).
N = 10001;
if det == "SAMP" && nargin > 4 && isfinite(Rs)
    u = T/(N - 1)*Rs;                     % sample spacing in symbols
    g = 0.618034;
    if u >= 1, u = floor(u) + g; else, u = 1/(floor(1/u) + g); end
    T = round(u*(N - 1)/Rs, 7);           % zero-span sweep-time resolution 100 ns
end
c = struct('label', "env", 'center', fc, 'rbw', 5e6, 'vbw', 5e6, 'points', N, 'time', T, ...
    'det', det, 'ref', ref, 'timeout', max(60, 5*T + 30));
if det ~= "SAMP", c.label = "burst"; end
end

function q = quickEdges(a, D, Rs0)
% -3 dB edges of the averaged trace (used to steer the measurement only).
f = a.f; lin = mean(10.^(levels(a, D)/10), 1);
n = numel(f); edge = [1:max(2, round(0.1*n)), min(n - 1, round(0.9*n)):n];
N = median(lin(edge));
s = max(lin - N, 0);
q = struct('Rs3', NaN, 'fc3', NaN);
if ~any(s > 0), return, end
c0 = sum(f.*s)/sum(s);
top = median(s(abs(f - c0) < 0.25*Rs0));
if ~(top > 0), return, end
[lo, hi] = edgePair(f, s, top, 3, c0, a.settings.rbw, Rs0);
q.Rs3 = hi - lo; q.fc3 = (lo + hi)/2;
end

function x = quickExcess(a, f0, Rs)
% dB rise of the average level within +-Rs of f0 over the outer part of the window.
lin = mean(10.^(a.P/10), 1); d = abs(a.f - f0);
inner = d <= Rs; outer = d >= 2*Rs;
x = 10*log10(mean(lin(inner))/median(lin(outer)));
if ~isfinite(x), x = 0; end
end

%% ===== instrument engine =====================================================

function D = measureRun(ip, opt)
fprintf('Connecting to %s:%d ...\n', ip, opt.Port);
ff = FieldFox(ip, opt.Port);
fprintf('Connected: %s\n', ff.IDN);
if ~ff.hasOption("233"), error('ff_srs4_analyzer:option', 'Needs the spectrum analyzer (option 233).'); end
ff.takeControl();                        % backs up every setting; release() restores them
cleanup = onCleanup(@() releaseInstrument(ff));
ff.setMode("SA");
D = newRun("instrument", opt);
D.idn = string(ff.IDN); D.options = ff.Options;
[notes, D.extgInstr, D.refSource] = prepareSA(ff, opt);
D.notes = [D.notes; notes];
st = containers.Map();
st('att') = "auto"; st('waor') = probe(ff, "CALC:MEAS:WAOR?");
if ~st('waor'), D.notes(end + 1, 1) = "CALC:MEAS:WAOR? not available: overload was not checked."; end
E.sweep = @(cfg, n) instSweep(ff, st, D.extgInstr, cfg, n);
E.zero = @(cfg, n) instZero(ff, st, D.extgInstr, cfg, n);
E.cross = @(s, Rs, n) instCross(ff, s, Rs, n);
E.setAtt = @(att) setAtt(st, att);
E.getAtt = @() st('att');
D = runFlow(E, D, opt);                 % keeps what was measured if a step fails
try
    e = ff.errors();
    if ~isempty(e), D.notes(end + 1, 1) = "Instrument errors at the end: " + strjoin(e, "; "); end
catch
end
clear cleanup                            % restore the instrument now
fprintf('Instrument released (all settings restored).\n');
end

function setAtt(st, att)
st('att') = att; %#ok<NASGU> (containers.Map is a handle)
end

function [notes, extgInstr, refSrc] = prepareSA(ff, opt)
notes = strings(0, 1);
u = tryQuery(ff, "AMPL:UNIT?");          % PG: display unit = TRAC:DATA? unit
if strlength(u) > 0 && ~startsWith(upper(u), "DBM")
    ff.write("AMPL:UNIT DBM");
    notes(end + 1, 1) = "Display unit was " + u + ": set to dBm.";
end
ff.write("SENS:MEAS:CHAN NONE");         % no channel measurement
ff.write("SENS:MEAS:INT NONE");          % no spectrogram
ff.write("TRIG:SOUR FREE");              % free run (zero span must not wait for a trigger)
ff.write("SENS:POW:RF:GAIN 0");          % safe start: preamp OFF
ff.write("SENS:POW:RF:ATT:AUTO 1");      % safe start: attenuation coupled to the reference level
ff.opc(30);
e = ff.errors();
if ~isempty(e), notes(end + 1, 1) = "Setup errors (ignored): " + strjoin(e, "; "); end
ff.write("SENS:AVER:TYPE POW");          % PG: power (RMS) averaging for the average detector; UNVERIFIED on A.12.45
ff.opc(30);
e = ff.errors();
if ~isempty(e), notes(end + 1, 1) = "AVER:TYPE POW not accepted (" + strjoin(e, "; ") + "): MATLAB averages power anyway."; end
% External gain (pad/cable): PG [:SENSe]:POWer[:RF]:EXTGain, negative = loss. UNVERIFIED on A.12.45:
% read back; if it is not taken, MATLAB applies ExtGain to the data instead.
ff.write(sprintf("SENS:POW:RF:EXTG %.4f", opt.ExtGain));
ff.opc(30);
e = ff.errors();
rb = str2double(tryQuery(ff, "SENS:POW:RF:EXTG?"));
ff.errors();
if isempty(e) && isfinite(rb) && abs(rb - opt.ExtGain) < 0.01
    extgInstr = opt.ExtGain;
else
    extgInstr = 0;
    if isfinite(rb), extgInstr = rb; end
    notes(end + 1, 1) = sprintf("POW:EXTG %g not confirmed (read back %g, errors: %s): ExtGain is applied in MATLAB.", ...
        opt.ExtGain, rb, strjoin(e, "; "));
end
refSrc = upper(erase(tryQuery(ff, "SENS:ROSC:SOUR?"), '"'));   % PG: INT / EXT; UNVERIFIED query
ff.errors();
if strlength(refSrc) == 0, refSrc = "unknown"; end
end

function a = instSweep(ff, st, extgInstr, cfg, n)
args = {'RBWAuto', false, 'RBW', cfg.rbw, 'Points', cfg.points, 'Detector', cfg.det, ...
    'TraceType', "CLRW", 'Average', 1, 'Preamp', false};
if isfield(cfg, 'start'), args = [args, {'Start', cfg.start, 'Stop', cfg.stop}];
else, args = [args, {'Center', cfg.center, 'Span', cfg.span}]; end
if isfield(cfg, 'vbw'), args = [args, {'VBW', cfg.vbw}]; end
if isfield(cfg, 'ref') && isfinite(cfg.ref), args = [args, {'RefLevel', cfg.ref}]; end
att = st('att');
if isnumeric(att), args = [args, {'AttenAuto', false, 'Atten', att}]; else, args = [args, {'AttenAuto', true}]; end
ff.saConfigure(args{:});
s = ff.saSettings();
P = zeros(n, s.points); t = NaT(n, 1); over = false; k = 0;
while k < n
    k = k + 1;
    ff.sweep(cfg.timeout);
    y = ff.saTrace();
    if numel(y) ~= s.points, error('ff_srs4_analyzer:points', '%s: trace has %d points, expected %d.', cfg.label, numel(y), s.points); end
    if k == 1 && st('waor') && waor(ff)
        if s.atten < 30                  % automatic attenuation increase, then start the set again
            na = min(30, 5*ceil((s.atten + 10)/5));
            fprintf('   ADC overload at %g dB attenuation: %g dB\n', s.atten, na);
            st('att') = na;
            ff.saConfigure('AttenAuto', false, 'Atten', na);
            s = ff.saSettings(); k = 0;
            continue
        end
        over = true;
        warning('ff_srs4_analyzer:overload', '%s: ADC overload even with 30 dB attenuation: add an external attenuator.', cfg.label);
    end
    P(k, :) = y; t(k) = datetime('now');
end
if st('waor') && ~over, over = waor(ff); end
a = makeAcq(cfg, s.f, P, t, s, extgInstr, over, false, NaN, "");
fprintf('   %-6s %s, RBW %s, %d pts, %s, att %g dB%s\n', cfg.label, fmtRange([s.start s.stop]), fmtHz(s.rbw), ...
    s.points, s.detector, s.atten, repmat('  OVERLOAD', 1, over));
end

function a = instZero(ff, st, extgInstr, cfg, n)
% Zero span: FREQ:SPAN 0, SWE:TIME (accepted on this unit), sample or average detector.
ff.write("SENS:FREQ:SPAN 0");
ff.write(sprintf("SENS:FREQ:CENT %.10g", cfg.center));
ff.write(sprintf("SENS:BAND:RES %.10g", cfg.rbw));
ff.write(sprintf("SENS:BAND:VID %.10g", cfg.vbw));
ff.write(sprintf("SENS:SWE:POIN %d", cfg.points));
if isfinite(cfg.ref), ff.write(sprintf("DISP:WIND:TRAC:Y:RLEV %g", cfg.ref)); end
att = st('att');
if isnumeric(att), ff.write(sprintf("SENS:POW:RF:ATT %g", att)); end
ff.opc(30);
e = ff.errors();
det = cfg.det;
ff.write("SENS:DET:FUNC " + det);
ff.opc(30);
eD = ff.errors();
if ~isempty(eD) && det ~= "SAMP"
    e = [e; "DET " + det + " rejected in zero span (" + strjoin(eD, "; ") + "): SAMP used"];
    det = "SAMP"; ff.write("SENS:DET:FUNC SAMP"); ff.opc(30); ff.errors();
end
ff.write(sprintf("SENS:SWE:TIME %.6g", cfg.time));
ff.opc(30);
eT = ff.errors();
if ~isempty(eT), e = [e; "SWE:TIME " + cfg.time + " rejected: " + strjoin(eT, "; ")]; end
span = ff.queryNum("SENS:FREQ:SPAN?");
if isempty(span) || span ~= 0, error('ff_srs4_analyzer:zerospan', 'Zero span not accepted (span %g).', span); end
N = ff.queryNum("SENS:SWE:POIN?");
s = ff.saSettings();
T = NaN; src = "SWE:TIME?";
if isempty(eT), T = str2double(tryQuery(ff, "SENS:SWE:TIME?")); end
if ~(T > 0)
    T = str2double(tryQuery(ff, "SENS:SWE:MTIM?")); src = "SWE:MTIM? (may include overhead)";
end
if ~(T > 0), T = cfg.time; src = "requested (not confirmed)"; end
P = zeros(n, N); t = NaT(n, 1);
for k = 1:n
    ff.sweep(cfg.timeout);
    y = ff.saTrace();
    if numel(y) ~= N, error('ff_srs4_analyzer:points', 'Zero span: %d points, expected %d.', numel(y), N); end
    P(k, :) = y; t(k) = datetime('now');
end
over = false;
if st('waor'), over = waor(ff); end
s.detector = det;
a = makeAcq(cfg, repmat(cfg.center, 1, N), P, t, s, extgInstr, over, true, T, src);
a.settings.errors = strjoin(e, "; ");
fprintf('   zero span %s: %d x %d points over %.4g s (%s), RBW %s, %s\n', fmtHz(cfg.center), n, N, T, src, ...
    fmtHz(s.rbw), det);
if ~isempty(e), warning('ff_srs4_analyzer:zerospan', '%s', strjoin(e, "; ")); end
end

function X = instCross(ff, s, Rs, n)
% FieldFox's own channel power and OBW on the zoom settings. CALC:MEAS:DATA? follows FORM
% (binary when REAL,32): CHP -> [P dBm, PSD dBm/Hz], OBW -> [OBW Hz, power dBm] (verified on
% this unit). Averaging on the instrument (TRAC1:TYPE AVG, INIT:REST) for this step only.
X = struct('done', false, 'ibw', NaN, 'chp', [NaN NaN], 'obw99', [NaN NaN], 'obw999', [NaN NaN], 'note', "");
ibw = min(1.35*Rs, s.stop - s.start);
notes = strings(0, 1);
try
    ff.write("SENS:MEAS:CHAN CHP");
    ff.write(sprintf("SENS:CME:IBW %.10g", ibw));
    ff.write(sprintf("SENS:AVER:COUN %d", n));
    ff.write("TRAC1:TYPE AVG");
    ff.opc(30); e = ff.errors();
    if ~isempty(e), notes(end + 1) = "CHP setup: " + strjoin(e, "; "); end
    ff.restartAverage(max(120, 10*n));
    X.chp = firstTwo(ff.queryBlock("CALC:MEAS:DATA?", "single", 10));
    X.ibw = ibw;
    for p = [99 99.9]
        ff.write("SENS:MEAS:CHAN OBW");
        ff.write(sprintf("SENS:OBW:PPOW %g", p));       % PG: always set it (default 10 %)
        ff.opc(30); e = ff.errors();
        ff.restartAverage(max(120, 10*n));
        v = firstTwo(ff.queryBlock("CALC:MEAS:DATA?", "single", 10));
        if ~isempty(e)                    % the percentage may not be what was asked: discard
            notes(end + 1) = sprintf("OBW %g %% setup rejected (%s): FieldFox OBW discarded", p, strjoin(e, "; ")); %#ok<AGROW>
            v = [NaN NaN];
        end
        if p == 99, X.obw99 = v; else, X.obw999 = v; end
    end
    X.done = true;
catch err
    notes(end + 1) = "stopped: " + err.message;
end
try
    ff.write("SENS:MEAS:CHAN NONE"); ff.write("TRAC1:TYPE CLRW"); ff.write("SENS:AVER:COUN 1");
    ff.opc(30); ff.errors();
    s2 = ff.saSettings();
    if abs(s2.start - s.start) > 1 || abs(s2.stop - s.stop) > 1 || abs(s2.rbw - s.rbw) > 1
        notes(end + 1) = "the channel measurement changed span/RBW (restored)";
        ff.saConfigure('Start', s.start, 'Stop', s.stop, 'RBW', s.rbw);
    end
catch err
    notes(end + 1) = "cleanup: " + err.message;
end
X.note = strjoin(notes, "; ");
fprintf('   CHP %.2f dBm (IBW %s), OBW99 %s, OBW99.9 %s\n', X.chp(1), fmtHz(X.ibw), fmtHz(X.obw99(1)), fmtHz(X.obw999(1)));
end

function v = firstTwo(v)
v = double(v(:).'); v(v < -9e4) = NaN;
if numel(v) < 2, v = [v NaN NaN]; end
v = v(1:2);
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
try
    ff.release();
catch err
    warning('ff_srs4_analyzer:release', 'Could not restore the instrument: %s', err.message);
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
e = ff.errors();                         % always drain the queue (an unknown query leaves -113)
ok = strlength(r) > 0 && isempty(e);
end

function a = makeAcq(cfg, f, P, t, s, extgInstr, over, zs, T, src)
set = struct('label', cfg.label, 'start', s.start, 'stop', s.stop, 'center', (s.start + s.stop)/2, ...
    'span', s.stop - s.start, 'rbw', s.rbw, 'vbw', s.vbw, 'points', size(P, 2), 'atten', s.atten, ...
    'preamp', s.preamp, 'detector', string(s.detector), 'refLevel', s.ref, 'extgInstr', extgInstr, ...
    'overload', over, 'zeroSpan', zs, 'sweepTime', T, 'timeFrom', src, 'errors', "");
if zs, set.center = cfg.center; set.span = 0; end
a = struct('f', reshape(f, 1, []), 'P', P, 't', t, 'settings', set);
end

%% ===== demo engine (synthetic SRS-4 + FieldFox measurement model) =========

function D = demoRun(opt)
fprintf('DEMO: synthetic SRS-4 downlink (ccsdsTMWaveformGenerator, %s RS+CC) through a FieldFox model.\n', ...
    opt.DemoModulation);
T = demoTruth(opt);
D = newRun("demo", opt);
D.idn = "SYNTHETIC DEMO (ff_srs4_analyzer FieldFox model) - not a measurement";
D.options = ["233" "235"]; D.extgInstr = opt.ExtGain; D.refSource = "INT (emulated)";
st = containers.Map(); st('att') = "auto";
E.sweep = @(cfg, n) demoSweep(T, st, cfg, n);
E.zero = @(cfg, n) demoZero(T, st, cfg, n);
E.cross = @(s, Rs, n) demoCross(T, st, s, Rs, n);
E.setAtt = @(att) setAtt(st, att);
E.getAtt = @() st('att');
D = runFlow(E, D, opt);
D.truth = rmfield(T, {'y', 'S', 'fn', 'S2', 'fn2', 'yr'});
end

function T = demoTruth(opt)
T.Rs = 1e6; if ~isempty(opt.ExpectedRate), T.Rs = opt.ExpectedRate; end
T.alpha = 0.2; if ~isempty(opt.ExpectedRolloff), T.alpha = opt.ExpectedRolloff; end
T.fnom = 2245e6; if ~isempty(opt.ExpectedFreq), T.fnom = opt.ExpectedFreq; end
T.fc = T.fnom*(1 + 0.4e-6);              % +0.40 ppm (about +0.9 kHz)
T.mod = opt.DemoModulation; T.M = find(["BPSK" "QPSK" "8PSK"] == T.mod);
T.Ptx = 27;                               % dBm at the SRS-4 TX connector
[T.Gtx, ~] = antGain(opt.TxAntennaGain, T.fc);
[T.Grx, ~] = antGain(opt.RxAntennaGain, T.fc);
T.ota = ~isempty(opt.Distance) || opt.ExtGain == 0;
if ~isempty(opt.Distance)
    T.Pin = T.Ptx + T.Gtx + T.Grx - fspl(opt.Distance, physconst('LightSpeed')/T.fc) - opt.PolLoss + opt.ExtGain;
elseif opt.ExtGain ~= 0
    T.Pin = T.Ptx + opt.ExtGain;          % conducted through the pad
else
    T.Pin = -12;                          % over the air, unknown geometry
end
T.extg = opt.ExtGain;
T.danl = -135;                            % dBm/Hz at 0 dB attenuation (-139 typ + 4 dB at 2.1-2.8 GHz)
T.K = 10;                                 % independent samples per bin of the average detector
T.h2 = -63; T.spurs = [T.fc + 26e6, -66]; % dBc
T.ambient = zeros(0, 3);
if T.ota                                  % centre, width, total dBm at the input
    T.ambient = [98e6 0.2e6 -45; 945e6 10e6 -58; 1842.5e6 15e6 -60; 2140e6 10e6 -62; 2437e6 20e6 -66; 2655e6 20e6 -62];
end
T.bits = 4144;                            % RS+CC frame (518 bytes)
T.Tframe = T.bits/(T.M*T.Rs); T.lineSpacing = 1/T.Tframe;
T.burstOn = 24*T.Tframe; T.burstOff = 20e-3;   % 1 idle + 20 data + 3 idle frames, then 20 ms off
% waveform: ccsdsTMWaveformGenerator, PA (Rapp, about 5 dB back-off), Bussgang split
sps = 8; T.fs = sps*T.Rs;
tm = ccsdsTMWaveformGenerator("WaveformSource", "synchronization and channel coding", ...
    "HasRandomizer", true, "HasASM", true, "ChannelCoding", "concatenated", "ConvolutionalCodeRate", "1/2", ...
    "RSMessageLength", 223, "RSInterleavingDepth", 1, "Modulation", T.mod, ...
    "PulseShapingFilter", "root raised cosine", "RolloffFactor", T.alpha, "FilterSpanInSymbols", 10, ...
    "SamplesPerSymbol", sps);
rs = RandStream('mt19937ar', 'Seed', 11);
x = zeros(0, 1);
while numel(x) < 60*sps*1000
    x = [x; tm(randi(rs, [0 1], tm.NumInputBits, 1))]; %#ok<AGROW>
end
x = x/rms(x);
nl = comm.MemorylessNonlinearity('Method', 'Rapp model', 'Smoothness', 2, 'OutputSaturationLevel', 10^(5/20));
y = nl(x);
G = (x'*y)/(x'*x); d = y - G*x;           % Bussgang: y = G x + uncorrelated distortion d
h = rcosdesign(T.alpha, 10, sps, 'sqrt');
fn = (-4096:4095).'/8192*sps;             % normalised frequency f/Rs
Hl = abs(freqz(h, 1, fn, sps)).^2;        % the generator's own RRC filter: smooth linear part
[Sd, fd] = pwelch(d, hann(4096), 2048, 4096, sps, 'centered', 'psd');
S = abs(G)^2*Hl/(sum(Hl)*(fn(2) - fn(1))) + interp1(fd, movmean(Sd, 25), fn, 'linear', 0);
S = S./(1 + (abs(fn)/3.6).^24);         % sps = 8 ends at +-4 Rs: taper instead of a step
T.S = S/(sum(S)*(fn(2) - fn(1))); T.fn = fn;     % unit power, per unit of f/Rs
[S2, f2] = pwelch(x.^2, hann(4096), 2048, 4096, sps, 'centered', 'psd');   % 2nd harmonic (squared)
T.S2 = movmean(S2, 9)/(sum(movmean(S2, 9))*(f2(2) - f2(1))); T.fn2 = f2;
T.y = y/rms(y); T.yr = [];
end

function att = demoAtt(T, st, cfg)
att = st('att');
if ~isnumeric(att)                        % FieldFox-like coupling to the reference level
    refIn = 0; if isfield(cfg, 'ref') && isfinite(cfg.ref), refIn = cfg.ref + T.extg; end
    att = min(30, max(0, 5*ceil((refIn + 10)/5)));
end
end

function a = demoSweep(T, st, cfg, n)
if isfield(cfg, 'start'), f = linspace(cfg.start, cfg.stop, cfg.points);
else, f = linspace(cfg.center - cfg.span/2, cfg.center + cfg.span/2, cfg.points); end
att = demoAtt(T, st, cfg);
over = demoOverload(T, att);
if over && ~isnumeric(st('att'))          % same rule as the instrument engine
    while over && att < 30
        att = min(30, att + 10); st('att') = att; over = demoOverload(T, att);
    end
end
e = demoExpected(T, f, cfg.rbw, cfg.det, att);
P = demoDetector(e.sig, e.noise, cfg.det, T.K, cfg.rbw, f, n) - T.extg;
s = struct('start', f(1), 'stop', f(end), 'rbw', cfg.rbw, 'vbw', getf(cfg, 'vbw', cfg.rbw), 'atten', att, ...
    'preamp', false, 'detector', cfg.det, 'ref', getf(cfg, 'ref', NaN));
t = datetime('now') + seconds((0:n - 1).'*0.2);
a = makeAcq(cfg, f, P, t, s, T.extg, over, false, NaN, "");
end

function o = demoOverload(T, att)
Ptot = 10*log10(10^(T.Pin/10) + sum(10.^(T.ambient(:, 3)/10)));
o = Ptot - att > -5;                      % emulated ADC over-range at -5 dBm mixer level
end

function e = demoExpected(T, f, rbw, det, att)
% Expected trace (mW in the RBW, input-referred) and noise: PSDs smoothed by the Gaussian
% RBW power response G = exp(-4 ln2 (df/RBW)^2) (area = ENBW = 1.0645 RBW) and, for the
% average/sample detectors, the bin width; peak detector: maximum over the bin.
df = (f(end) - f(1))/max(numel(f) - 1, 1);
sig = zeros(size(f));
Pin = 10^(T.Pin/10);
lineMode = rbw < T.lineSpacing/3;
% main signal (or idle-frame lines)
if lineMode
    k = ceil((f(1) - 3*rbw - T.fc)/T.lineSpacing):floor((f(end) + 3*rbw - T.fc)/T.lineSpacing);
    fk = T.fc + k*T.lineSpacing;
    Pk = Pin*interp1(T.fn, T.S, (fk - T.fc)/T.Rs, 'linear', 0)/T.Rs*T.lineSpacing;
    for i = 1:numel(fk)
        m = abs(f - fk(i)) < 4*rbw;
        sig(m) = sig(m) + Pk(i)*exp(-4*log(2)*((f(m) - fk(i))/rbw).^2);
    end
else
    sig = sig + smoothPsd(T.fc + T.fn*T.Rs, Pin*T.S/T.Rs, f, rbw, det);
end
sig = sig + smoothPsd(2*T.fc + T.fn2*T.Rs, Pin*10^(T.h2/10)*T.S2/T.Rs, f, rbw, det);
for i = 1:size(T.spurs, 1)
    sig = sig + cwTrace(f, T.spurs(i, 1), Pin*10^(T.spurs(i, 2)/10), rbw, df, det);
end
for i = 1:size(T.ambient, 1)
    fa = T.ambient(i, 1) + linspace(-0.6, 0.6, 241)*T.ambient(i, 2);
    pa = 10^(T.ambient(i, 3)/10)/T.ambient(i, 2)*(abs(fa - T.ambient(i, 1)) <= T.ambient(i, 2)/2);
    sig = sig + smoothPsd(fa, pa, f, rbw, det);
end
e.sig = sig;
e.noise = 10^((T.danl + att)/10)*1.0645*rbw*ones(size(f));
end

function y = smoothPsd(fx, px, f, rbw, det)
% Trace of a PSD px(fx) (mW/Hz) seen through the RBW filter at the points f.
y = zeros(size(f));
lo = min(fx) - 4*rbw; hi = max(fx) + 4*rbw;
m = f >= lo & f <= hi;
if ~any(m), return, end
df = (f(end) - f(1))/max(numel(f) - 1, 1);
h = min([rbw/4, max(df/2, eps), median(diff(fx))]);
h = max(h, (hi - lo)/2e5);
g = lo:h:hi;
pg = interp1(fx, px, g, 'linear', 0);
u = (-ceil(3*rbw/h):ceil(3*rbw/h))*h;   % odd length: no half-step shift
G = exp(-4*log(2)*(u/rbw).^2)*h;          % area = ENBW
tr = conv(pg, G, 'same');
w = max(1, round(df/h));
if w > 1
    if det == "POS", tr = movmax(tr, w); else, tr = movmean(tr, w); end
end
y(m) = interp1(g, tr, f(m), 'linear', 0);
end

function y = cwTrace(f, f0, P, rbw, df, det)
d = abs(f - f0);
if det == "POS" || df > rbw, d = max(0, d - df/2); end
y = P*exp(-4*log(2)*(d/rbw).^2);
end

function P = demoDetector(sig, noise, det, K, rbw, f, n)
% Detector statistics: average -> Gamma(K) power fluctuation correlated over the RBW,
% sample -> exponential, peak -> noise read about 5 dB high.
N = numel(f); df = (f(end) - f(1))/max(N - 1, 1);
w = max(1, rbw/df);
P = zeros(n, N);
for k = 1:n
    z = randn(1, N + 8*ceil(w));
    if w > 1.5
        g = exp(-0.5*((-ceil(3*w):ceil(3*w))/(w/2.35)).^2);
        z = conv(z, g/sqrt(sum(g.^2)), 'same');
    end
    z = z(1:N);
    switch det
        case "SAMP", fl = -log(max(normcdf(z), 1e-12));
        case "POS",  fl = 10^(0.5)*exp(0.25*z - 0.25^2/2);
        otherwise,   sg = sqrt(log(1 + 1/K)); fl = exp(sg*z - sg^2/2);
    end
    if det == "POS"
        P(k, :) = 10*log10(sig.*exp(0.05*z) + noise.*fl);
    else
        P(k, :) = 10*log10((sig + noise).*fl);
    end
end
end

function a = demoZero(T, st, cfg, n)
N = cfg.points; Tsw = cfg.time; dt = Tsw/(N - 1); rbw = min(cfg.rbw, 5e6);
att = demoAtt(T, st, cfg);
yr = demoFiltered(T, rbw);
L = numel(yr);
Pin = 10^(T.Pin/10);
Pn = 10^((T.danl + att)/10)*1.0645*rbw;
Pamb = 0;                                 % ambient inside a 5 MHz RBW around fc: none in the demo
P = zeros(n, N);
period = T.burstOn + T.burstOff;
for k = 1:n
    t0 = 10*rand;
    tk = t0 + (0:N - 1)*dt;
    on = mod(tk, period) < T.burstOn;
    i0 = mod(round(tk*T.fs), L) + 1;
    if cfg.det == "SAMP"
        if dt*rbw > 2
            nz = sqrt(Pn/2)*(randn(1, N) + 1i*randn(1, N));
        else                              % correlated noise: filter white noise with the same RBW
            m = round((N - 1)*dt*T.fs) + 1;
            nz = gaussFilter(sqrt(Pn/2)*(randn(1, m) + 1i*randn(1, m)), rbw, T.fs);
            nz = nz(min(m, round((0:N - 1)*dt*T.fs) + 1));
        end
        p = abs(sqrt(Pin)*yr(i0).'.*on + nz).^2 + Pamb;
    else                                  % average detector: mean power over each bin
        mw = min(64, max(1, round(dt*T.fs)));
        ps = zeros(1, N);
        for j = 0:mw - 1
            ps = ps + abs(yr(mod(i0 - 1 + j*max(1, floor(dt*T.fs/mw)), L) + 1).').^2/mw;
        end
        Kn = max(1, dt*rbw);
        p = Pin*ps.*on + Pn*gamrnd(Kn, 1/Kn, 1, N) + Pamb;
    end
    P(k, :) = 10*log10(p) - T.extg;
end
s = struct('start', cfg.center, 'stop', cfg.center, 'rbw', rbw, 'vbw', cfg.vbw, 'atten', att, ...
    'preamp', false, 'detector', cfg.det, 'ref', cfg.ref);
t = datetime('now') + seconds((0:n - 1).'*Tsw);
a = makeAcq(cfg, repmat(cfg.center, 1, N), P, t, s, T.extg, false, true, Tsw, "SWE:TIME? (emulated)");
end

function yr = demoFiltered(T, rbw)
persistent cache key
k = sprintf('%.6g_%.6g_%.6g', rbw, T.Rs, T.alpha) + T.mod;
if ~isempty(key) && key == k, yr = cache; return, end
yr = gaussFilter(T.y.', rbw, T.fs).';
cache = yr; key = k;
end

function y = gaussFilter(x, rbw, fs)
% Gaussian filter with 3 dB power bandwidth rbw (two-sided +-rbw/2): sigma_t = sqrt(ln2)/(pi rbw).
sg = sqrt(log(2))/(pi*rbw)*fs;
if sg < 0.3, y = x; return, end
u = -ceil(4*sg):ceil(4*sg);
h = exp(-0.5*(u/sg).^2); h = h/sum(h);
y = conv(x, h, 'same');
end

function X = demoCross(T, st, s, Rs, n)
cfg = struct('label', "cross", 'start', s.start, 'stop', s.stop, 'points', s.points, 'rbw', s.rbw, ...
    'det', "AVER", 'ref', s.refLevel);
a = demoSweep(T, st, cfg, n);
lin = mean(10.^(a.P/10), 1); f = a.f; df = f(2) - f(1); enbw = 1.0645*s.rbw;
ibw = min(1.35*Rs, s.stop - s.start);
m = abs(f - (s.start + s.stop)/2) <= ibw/2;
Pch = 10*log10(sum(lin(m))*df/enbw);
c = cumsum(lin)/sum(lin);
o = @(p) interp1(c + (1:numel(c))*1e-12, f, 0.5 + p/200) - interp1(c + (1:numel(c))*1e-12, f, 0.5 - p/200);
Ptot = 10*log10(sum(lin)*df/enbw);
X = struct('done', true, 'ibw', ibw, 'chp', [Pch, Pch - 10*log10(ibw)], 'obw99', [o(99), Ptot], ...
    'obw999', [o(99.9), Ptot], 'note', "emulated (demo)");
end

%% ===== loading and saving =====================================================

function D = loadRun(file, opt)
if ~isfile(file), error('ff_srs4_analyzer:file', 'File not found: %s', file); end
S = load(file);
if isfield(S, 'run_data')
    D = S.run_data;
    D.source = "file: " + D.source;
    if opt.ExtGain ~= 0 && opt.ExtGain ~= D.extGain
        D.notes(end + 1, 1) = sprintf("ExtGain changed from %g dB (measurement) to %g dB (re-analysis).", D.extGain, opt.ExtGain);
        D.extGain = opt.ExtGain;
    end
    fprintf('Loaded run (%s) from %s\n', D.source, file);
    return
end
if ~isfield(S, 'f') || ~isfield(S, 'P'), error('ff_srs4_analyzer:file', '%s has neither run_data nor f/P.', file); end
D = newRun("recording: " + string(file), opt);
if isfield(S, 'idn'), D.idn = string(S.idn); end
set = struct('rbw', NaN, 'start', S.f(1), 'stop', S.f(end));
if isfield(S, 'settings'), set = S.settings; end
rbw = getf(set, 'rbw', NaN);
if ~isfinite(rbw)
    rbw = 3*(S.f(end) - S.f(1))/(numel(S.f) - 1);
    D.notes(end + 1, 1) = "RBW not in the recording: assumed " + fmtHz(rbw) + ".";
end
s = struct('start', S.f(1), 'stop', S.f(end), 'rbw', rbw, 'vbw', getf(set, 'vbw', NaN), ...
    'atten', getf(set, 'atten', NaN), 'preamp', getf(set, 'preamp', false), ...
    'detector', string(getf(set, 'detector', "unknown")), 'ref', getf(set, 'ref', NaN));
t = NaT(size(S.P, 1), 1); if isfield(S, 't'), t = S.t(:); end
cfg = struct('label', "zoom");
D.acq.zoom = makeAcq(cfg, S.f, S.P, t, s, 0, false, false, NaN, "");
Rs0 = opt.ExpectedRate;
if isempty(Rs0)
    q = quickEdges(D.acq.zoom, D, (S.f(end) - S.f(1))/3);
    Rs0 = q.Rs3;
end
D.plan = struct('Rs0', Rs0, 'fc0', NaN, 'ref', NaN);
D.notes(end + 1, 1) = "Recording: only the spectrum analysis is possible (no mask, harmonic or zero-span data).";
if ~contains(lower(s.detector), "aver")
    D.notes(end + 1, 1) = "Detector " + s.detector + ": the average (RMS) detector is needed for unbiased power and shape.";
end
end

function files = saveRecordings(D, outDir)
% Shared recording format (fieldfox_live): f, P (dBm as displayed), t, settings, idn.
files = strings(0, 1);
if ~isfield(D, 'acq') || isempty(fieldnames(D.acq)), return, end
d = fullfile(outDir, "rec");
if ~isfolder(d), mkdir(d); end
for nm = string(fieldnames(D.acq)).'
    a = D.acq.(nm);
    settings = a.settings; settings.extGainTotal = D.extGain;
    settings.levelNote = "P as displayed; DUT-referred = P + extgInstr - extGainTotal";
    f = a.f; P = a.P; t = a.t; idn = D.idn;
    files(end + 1, 1) = fullfile(d, "rec_" + nm + ".mat"); %#ok<AGROW>
    save(files(end), 'f', 'P', 't', 'settings', 'idn');
end
end

function a = getAcq(D, name)
a = [];
if isfield(D, 'acq') && isfield(D.acq, name), a = D.acq.(name); end
end

function P = levels(a, D)
% DUT-referred dBm: display + external gain already applied by the instrument - total ExtGain.
P = a.P + a.settings.extgInstr - D.extGain;
end

function p = toDisplay(pCorr, D)
p = pCorr + D.extGain - D.extgInstr;
end

%% ===== analysis =================================================================

function A = analyzeRun(D, opt)
A = struct('ok', false, 'notes', D.notes, 'search', D.search);
z = getAcq(D, "zoom");
A.power = struct('ebno', table()); A.frames = table();
if isempty(z) || ~isfield(D.plan, 'Rs0') || ~isfinite(D.plan.Rs0)
    A.limits = limitRow("Downlink spectrum", "not measured", "-", "NOT MEASURED", "no emission found / no zoom data");
    A.conclusions = "No zoomed spectrum of the downlink: check that the SRS-4 transmits in the band.";
    return
end
Rs0 = D.plan.Rs0;
mk = getAcq(D, "mask");
A.floor = noiseFloor(z, mk, D, Rs0);
A.spec = spectrumAnalysis(z, D, opt, Rs0, A.floor);
A.ref = referenceOverlay(z, D, A.spec, A.floor, opt);
A.power = powerAnalysis(z, D, A.spec, A.floor, opt);
A.inst = instCompare(D, A);
A.freq = freqAnalysis(A.spec, D, opt);
A.ecss = ecssObw(A.spec);
if isempty(mk), mk = z; end
A.mask = maskAnalysis(mk, D, A.spec, A.floor);
A.harm = harmAnalysis(getAcq(D, "harm"), getAcq(D, "harm2"), D, A.spec, A.power);
A.wide = wideAnalysis(getAcq(D, "wide"), D, A.spec, A.power, opt);
A.burst = burstAnalysis(getAcq(D, "burst"), D, A.spec, A.power, A.floor);
A.env = envelopeAnalysis(getAcq(D, "env"), D, A.spec, A.power, A.floor, A.burst);
A.lines = lineAnalysis(getAcq(D, "line"), D, A.spec);
[A.frameRs, A.frameRsSource, relU] = frameRate(A, opt);
A.frames = frameTable(A.frameRs, relU, A.burst);
A.limits = limitsTable(A, D, opt);
A.conclusions = conclusions(A, D, opt);
A.ok = true;
end

function Sx = searchEmission(a, D, opt)
% Strongest wideband emission of the coarse sweep: runs of bins > floor + 6 dB.
f = a.f; lin = mean(10.^(levels(a, D)/10), 1);
rbw = a.settings.rbw; enbw = 1.0645*rbw; df = (f(end) - f(1))/(numel(f) - 1);
fl = prctile(lin, 25);
sm = movmean(lin, max(1, round(rbw/df)));
[s, e] = runs(10*log10(sm) > 10*log10(fl) + 6, 2);
rows = zeros(0, 6);
for r = 1:numel(s)
    i = s(r):e(r);
    ex = max(lin(i) - fl, 0);
    pw = sum(ex)*df/enbw;
    [pk, j] = max(sm(i));
    above = i(sm(i) >= pk/2);
    w3 = (f(above(end)) - f(above(1))) + df;
    rows(end + 1, :) = [sum(f(i).*ex)/max(sum(ex), eps), w3, 10*log10(pw), f(i(j)), f(i(1)), f(i(end))]; %#ok<AGROW>
end
T = array2table(rows, 'VariableNames', {'Centroid_Hz', 'Width3dB_Hz', 'Power_dBm', 'Peak_Hz', 'Start_Hz', 'Stop_Hz'});
T.Wideband = T.Width3dB_Hz > 3*rbw;
T = sortrows(T, 'Power_dBm', 'descend');
Sx = struct('found', false, 'fc', NaN, 'bw3', NaN, 'Pch', NaN, 'wide', false, 'table', T, ...
    'floor_dBmHz', 10*log10(fl/enbw), 'PinInput', 10*log10(sum(lin)*df/enbw) + D.extGain);
if isempty(T), return, end
cand = 1:height(T);
if ~isempty(opt.ExpectedFreq)
    tol = 2*max([getf(opt, 'ExpectedRate', 0), 3*rbw, 100e3]);
    near = find(T.Start_Hz - tol <= opt.ExpectedFreq & T.Stop_Hz + tol >= opt.ExpectedFreq);
    if ~isempty(near), cand = near; end
end
k = cand(find(T.Wideband(cand), 1));
if isempty(k), k = cand(1); end
Sx.found = true; Sx.fc = T.Centroid_Hz(k); Sx.bw3 = T.Width3dB_Hz(k); Sx.Pch = T.Power_dBm(k);
Sx.wide = T.Wideband(k);
end

function F = noiseFloor(z, mk, D, Rs0)
% Analyzer noise density (DUT-referred, mW/Hz) from the outer part of the mask sweep when it
% has the zoom's attenuation (it reaches far from the carrier), else from the zoom span.
% Spectral regrowth there would make it an upper bound.
a = z; src = "zoom span edges (PA regrowth may dominate there: N0 is an upper bound)";
F.upperBound = true;
if ~isempty(mk) && mk.settings.atten == z.settings.atten, a = mk; src = "mask sweep edges"; F.upperBound = false; end
f = a.f; lin = mean(10.^(levels(a, D)/10), 1); enbw = 1.0645*a.settings.rbw;
c = (f(1) + f(end))/2; d = abs(f - c);
edge = d > 0.8*max(d) & d > 0.75*Rs0;
if nnz(edge) < 5, edge = d > 0.8*max(d); end
F.N0 = median(lin(edge))/enbw;
F.N0dB = 10*log10(F.N0);
F.source = src + sprintf(" (|f - fc| > %.2f Rs)", min(d(edge))/Rs0);
F.atten = a.settings.atten;
end

function S = spectrumAnalysis(z, D, opt, Rs0, F)
f = z.f; P = levels(z, D); n = size(P, 1);
rbw = z.settings.rbw; enbw = 1.0645*rbw; df = (f(end) - f(1))/(numel(f) - 1);
lin = mean(10.^(P/10), 1); pdB = 10*log10(lin);
Nbin = F.N0*enbw; NdB = 10*log10(Nbin);
sig = lin - Nbin;
sp = max(sig, 0);
c0 = sum(f.*sp)/sum(sp);
top = median(sig(abs(f - c0) < 0.25*Rs0));
if ~(top > 0), error('ff_srs4_analyzer:nosignal', 'No signal above the noise floor in the zoom.'); end
S = struct('f', f, 'pdB', pdB, 'rbw', rbw, 'enbw', enbw, 'df', df, 'nSweeps', n, 'Nbin_dB', NdB, ...
    'top_dB', 10*log10(top), 'snrTop', 10*log10(top) - NdB);
% -3 dB and -26 dB bandwidths (noise subtracted, local regression across each edge)
[S.f3lo, S.f3hi] = edgePair(f, sig, top, 3, c0, rbw, Rs0);
S.Rs3 = S.f3hi - S.f3lo; S.fc3 = (S.f3lo + S.f3hi)/2;
[S.f26lo, S.f26hi] = edgePair(f, sig, top, 26, c0, rbw, Rs0);
S.ebw26 = S.f26hi - S.f26lo;
S.ebw26Limited = S.snrTop < 32;
% occupied bandwidth on the noise-subtracted trace, inside the region where the signal shows
sm = movmean(lin, max(1, round(rbw/df)));
on = find(sm > Nbin*10^(0.3));
use = false(size(f));
if ~isempty(on), use(f >= f(on(1)) - 2*rbw & f <= f(on(end)) + 2*rbw) = true; end
s2 = sp; s2(~use) = 0;
[S.obw99, S.obw99lo, S.obw99hi] = obw(s2(:), f(:), [], 99);
[S.obw999, S.obw999lo, S.obw999hi] = obw(s2(:), f(:), [], 99.9);
S.obw99raw = obw(lin(:), f(:), [], 99);  % noise included, as the FieldFox computes it
S.fcObw = (S.obw99lo + S.obw99hi)/2;
S.obw99Limited = S.snrTop < 20; S.obw999Limited = S.snrTop < 30;
% raised-cosine fit
S.fit = fitSpectrum(f, P, rbw, df, Rs0, c0, S, NdB, opt);
S.fc = S.fit.fc; S.Rs = S.fit.Rs; S.alpha = S.fit.alpha;
S.alphaObw = alphaFromObw(S.obw99/S.Rs);
S.alphaObw3 = alphaFromObw(S.obw99/S.Rs3);
end

function [lo, hi] = edgePair(f, s, top, xdb, c0, rbw, Rs0)
lev = top*10^(-xdb/10);
df = (f(end) - f(1))/(numel(f) - 1);
sm = movmean(s, max(1, round(rbw/df)));
L = 10*log10(lev);
sdb = 10*log10(max(s, top*1e-6));
m = max(2, round(0.03*Rs0/df));
il = find(f < c0 & sm < lev, 1, 'last');
ih = find(f > c0 & sm < lev, 1, 'first');
lo = cross(il, il + 1); hi = cross(ih - 1, ih);
    function x = cross(i1, i2)
        x = NaN;
        if isempty(i1) || isempty(i2) || i1 < 1 || i2 > numel(f), return, end
        j = max(1, i1 - m + 1):min(numel(f), i2 + m - 1);
        if numel(j) < 3, return, end
        u = (f(j) - f(i1))/Rs0;
        if numel(j) >= 6
            p = polyfit(u, sdb(j), 2);    % the RC edge is curved in dB: a line would bias the width
            r = roots(p - [0 0 L]);
            r = real(r(abs(imag(r)) < 1e-12));
            r = r(r >= u(1) & r <= u(end));
            if ~isempty(r), [~, k] = min(abs(r - (f(i2) - f(i1))/2/Rs0)); x = f(i1) + r(k)*Rs0; return, end
        end
        p = polyfit(u, sdb(j), 1);
        if p(1) == 0, return, end
        x = f(i1) + (L - p(2))/p(1)*Rs0;
        if x < f(max(1, i1 - m)) || x > f(min(numel(f), i2 + m)), x = NaN; end
    end
end

function a = alphaFromObw(r)
% Roll-off whose ideal raised-cosine PSD has OBW99/Rs = r.
persistent Atab Rtab
if isempty(Atab)
    Atab = 0.02:0.005:1; Rtab = zeros(size(Atab));
    fn = linspace(-1.1, 1.1, 22001);
    for k = 1:numel(Atab)
        c = cumsum(rcShape(fn, 0, 1, Atab(k))); c = c/c(end) + (1:numel(c))*1e-13;
        Rtab(k) = interp1(c, fn, 0.995) - interp1(c, fn, 0.005);
    end
end
a = interp1(Rtab, Atab, r, 'linear', NaN);
end

function s = rcShape(f, fc, Rs, a)
% Raised-cosine power spectrum (= |RRC|^2), peak 1, -3 dB at fc +- Rs/2.
d = abs(f - fc); f1 = (1 - a)*Rs/2; f2 = (1 + a)*Rs/2;
s = double(d <= f1);
m = d > f1 & d <= f2;
s(m) = 0.5*(1 + cos(pi/(a*Rs)*(d(m) - f1)));
end

function y = rcTrace(f, fc, Rs, a, rbw, binw)
% Raised cosine seen through the Gaussian RBW (and the bin of an average detector), peak ~1.
h = min(rbw/4, (1 + a)*Rs/400);
if binw > h, h = min(h, binw/2); end
lo = fc - (1 + a)*Rs/2 - 4*rbw - binw; hi = fc + (1 + a)*Rs/2 + 4*rbw + binw;
g = lo:h:hi;
u = (-ceil(3*rbw/h):ceil(3*rbw/h))*h;   % odd length: no half-step shift
G = exp(-4*log(2)*(u/rbw).^2); G = G/sum(G);
t = conv(rcShape(g, fc, Rs, a), G, 'same');
w = round(binw/h);
if w > 1, t = movmean(t, w); end
y = interp1(g, t, f, 'linear', 0);
end

function Ft = fitSpectrum(f, P, rbw, df, Rs0, c0, S, NdB, opt)
% lsqcurvefit of 10 log10(P0 RC(f - fc; Rs, alpha) (x) RBW + N0) to the averaged trace,
% only where it is within FitRange dB of the top (the skirts are where truncated RRC
% filters and PA regrowth differ from the ideal RC).
lin = mean(10.^(P/10), 1); pdB = 10*log10(lin);
topDb = S.top_dB;
% fit region: symmetric about the centre of the -FitRange crossings, so that the boundary
% (where the model mismatch is largest) cannot pull fc to one side
sig = lin - 10^(NdB/10);
[lo, hi] = edgePair(f, sig, 10^(topDb/10), opt.FitRange, c0, rbw, Rs0);
if isfinite(lo) && isfinite(hi)
    use = abs(f - (lo + hi)/2) <= (hi - lo)/2;
else
    use = pdB > 10*log10(10^(topDb/10) + 10^(NdB/10)) - opt.FitRange;
end
fixN = (topDb - NdB) > opt.FitRange + 3;          % noise not inside the fit region
binw = df*(df > rbw/2);
fm = c0;
r0 = S.Rs3/Rs0; if ~isfinite(r0) || r0 < 0.6 || r0 > 1.4, r0 = 1; end
a0 = alphaFromObw(S.obw99/(r0*Rs0)); if ~isfinite(a0), a0 = 0.25; end
a0 = min(max(a0, 0.05), 0.9);
x0 = [0, r0, a0, topDb, NdB];
lb = [-0.3, 0.6*r0, 0.02, topDb - 10, NdB - 15];
ub = [0.3, 1.4*r0, 1, topDb + 10, topDb - 3];
if fixN
    mdl = @(x, ff) rcModel(ff, fm, Rs0, [x NdB], rbw, binw);
    x0 = x0(1:4); lb = lb(1:4); ub = ub(1:4);
else
    mdl = @(x, ff) rcModel(ff, fm, Rs0, x, rbw, binw);
end
x0 = min(max(x0, lb), ub);
if nnz(use) < numel(x0) + 3
    error('ff_srs4_analyzer:fit', ['Only %d trace points lie within %g dB of the top: too few for the raised-cosine ' ...
        'fit. Is this an RRC-PSK downlink zoomed with RBW << Rs? (A CW line or a coarse recording of another ' ...
        'signal cannot be fitted; check ''ExpectedRate''.)'], nnz(use), opt.FitRange);
end
o = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxFunctionEvaluations', 4000);
[x, ~, res, ef, ~, ~, J] = lsqcurvefit(mdl, x0, f(use), pdB(use), lb, ub, o);
wn = warning('off', 'stats:nlparci:IllConditionedJacobian');
try
    ci = nlparci(x, res, 'jacobian', full(J));
catch
    ci = NaN(numel(x), 2);                         % degenerate Jacobian: no model-based interval
end
warning(wn);
infl = sqrt(max(1, 1.0645*rbw/df));                 % neighbouring bins share the RBW: fewer independent residuals
ci = x(:) + (ci - x(:))*infl;
Ft = struct('exitflag', ef, 'fixedN0', fixN, 'use', use, 'nUsed', nnz(use), 'inflate', infl);
Ft.fc = fm + x(1)*Rs0; Ft.Rs = x(2)*Rs0; Ft.alpha = x(3); Ft.P0dB = x(4);
Ft.N0dB = NdB; if ~fixN, Ft.N0dB = x(5); end
Ft.ciFc = fm + ci(1, :)*Rs0; Ft.ciRs = ci(2, :)*Rs0; Ft.ciAlpha = ci(3, :);
xm = x; if fixN, xm = [x NdB]; end
Ft.model = rcModel(f, fm, Rs0, xm, rbw, binw);
Ft.resid = pdB - Ft.model;
Ft.rms = rms(Ft.resid(use));
% bootstrap over the sweeps (resample whole sweeps, refit)
n = size(P, 1);
Ft.boot = zeros(0, 3); Ft.bciFc = [NaN NaN]; Ft.bciRs = [NaN NaN]; Ft.bciAlpha = [NaN NaN];
if opt.Boot > 0 && n >= 3
    rs = RandStream('mt19937ar', 'Seed', 7);
    ob = optimoptions(o, 'MaxFunctionEvaluations', 600);
    B = zeros(opt.Boot, 3);
    Pl = 10.^(P/10);
    for b = 1:opt.Boot
        k = randi(rs, n, n, 1);
        pb = 10*log10(mean(Pl(k, :), 1));
        xb = lsqcurvefit(mdl, x, f(use), pb(use), lb, ub, ob);
        B(b, :) = [fm + xb(1)*Rs0, xb(2)*Rs0, xb(3)];
    end
    Ft.boot = B;
    % normal interval from the bootstrap spread with a small-sample t factor: percentile
    % intervals from only 10-30 sweeps are too narrow (simulation)
    k = tinv(0.975, n - 1)*sqrt(n/(n - 1));
    est = [Ft.fc Ft.Rs Ft.alpha]; sd = std(B, 0, 1);
    Ft.bciFc = est(1) + [-1 1]*k*sd(1); Ft.bciRs = est(2) + [-1 1]*k*sd(2);
    Ft.bciAlpha = est(3) + [-1 1]*k*sd(3);
end
end

function y = rcModel(f, fm, Rs0, x, rbw, binw)
y = 10*log10(10^(x(4)/10)*rcTrace(f, fm + x(1)*Rs0, x(2)*Rs0, x(3), rbw, binw) + 10^(x(5)/10));
end

function F = referenceOverlay(~, ~, S, Fl, opt)
% ccsdsTMWaveformGenerator (QPSK, RS+CC) PSDs at the legal roll-offs, scaled to the fitted
% fc/Rs, through the same RBW filter, level fitted; smallest residual RMS = best roll-off.
alphas = unique([0.2 0.25 0.35, opt.ExpectedRolloff]);
f = S.f; pdB = S.pdB; rbw = S.rbw; df = S.df;
Nbin = Fl.N0*S.enbw;
% level from the flat top; shape compared over the RC-fit region (the skirts below it are
% dominated by the PA and the real filter length, which the reference does not know)
reg = S.fit.use;
topR = abs(f - S.fc) < 0.3*S.Rs;
F = struct('alpha', alphas, 'rms', NaN(size(alphas)), 'best', NaN, 'curve', [], 'resid', [], 'region', reg, ...
    'note', "not computed (too few points in the fit region)");
if nnz(reg) < 10 || nnz(topR) < 3, return, end
curves = cell(size(alphas));
for k = 1:numel(alphas)
    [Sr, fn] = refPsd(alphas(k));
    tr = smoothPsd(S.fc + fn*S.Rs, Sr/S.Rs, f, rbw, "AVER");
    tr = max(tr, realmin);
    c = 10*log10(mean(10.^(pdB(topR)/10) - Nbin)/mean(tr(topR)));
    curves{k} = 10*log10(10^(c/10)*tr + Nbin);
    F.rms(k) = rms(pdB(reg) - curves{k}(reg));
end
[~, i] = min(F.rms);
F.best = alphas(i); F.curve = curves{i}; F.resid = pdB - F.curve;
F.curves = curves;
F.note = sprintf("ccsdsTMWaveformGenerator QPSK RS+CC, 10-symbol RRC, %d frames, pwelch", 12);
if df > rbw/2, F.note = F.note + "; bins wider than RBW/2"; end
end

function [S, fn] = refPsd(alpha)
% Normalised PSD (unit power, per unit of f/Rs) of the SRS-4-like CCSDS waveform.
persistent cache
if isempty(cache), cache = containers.Map(); end
key = sprintf('%.4f', alpha);
if isKey(cache, key), v = cache(key); S = v{1}; fn = v{2}; return, end
tm = ccsdsTMWaveformGenerator("WaveformSource", "synchronization and channel coding", ...
    "HasRandomizer", true, "HasASM", true, "ChannelCoding", "concatenated", "ConvolutionalCodeRate", "1/2", ...
    "RSMessageLength", 223, "RSInterleavingDepth", 1, "Modulation", "QPSK", ...
    "PulseShapingFilter", "root raised cosine", "RolloffFactor", alpha, "FilterSpanInSymbols", 10, ...
    "SamplesPerSymbol", 8);
rs = RandStream('mt19937ar', 'Seed', 3);
x = zeros(0, 1);
for k = 1:12, x = [x; tm(randi(rs, [0 1], tm.NumInputBits, 1))]; end %#ok<AGROW>
[S, fn] = pwelch(x, hann(4096), 2048, 4096, 8, 'centered', 'psd');
S = movmean(S, 3); S = S/(sum(S)*(fn(2) - fn(1)));
cache(key) = {S, fn};
end

function W = powerAnalysis(~, D, S, Fl, opt)
f = S.f; lin = 10.^(S.pdB/10); psd = lin/S.enbw;
ibw = min((1 + S.alpha)*S.Rs*1.05, f(end) - f(1));
band = S.fc + [-1 1]*ibw/2;
band = [max(band(1), f(1)) min(band(2), f(end))];
Ptot = bandpower(psd(:), f(:), band, 'psd');
Psig = Ptot - Fl.N0*diff(band);
W.ibw = diff(band); W.ChP_dBm = 10*log10(Ptot);
W.C_dBm = 10*log10(max(Psig, realmin));
W.Cfit_dBm = S.fit.P0dB - 10*log10(S.enbw) + 10*log10(S.Rs);   % flat-top PSD x Rs
W.N0_dBmHz = Fl.N0dB;
W.CN0 = W.C_dBm - Fl.N0dB;
W.EsN0 = W.CN0 - 10*log10(S.Rs);
W.nfAnalyzer = Fl.N0dB + 174;
mods = ["BPSK"; "QPSK"; "8PSK"]; m = [1; 2; 3];
rc = 223/255/2;                           % CCSDS RS(255,223) + CC r=1/2 = 0.4373
W.ebno = table(mods, m, W.EsN0 - 10*log10(m), W.EsN0 - 10*log10(m*rc), S.Rs*m/1e3, S.Rs*m*0.419/1e3, ...
    'VariableNames', {'Modulation', 'BitsPerSymbol', 'EbN0_channel_dB', 'EbN0_info_RSCC_dB', ...
    'ChannelRate_kbps', 'Payload_RSCC_kbps'});
W.table7 = table7(S.Rs);
% transmitter power: conducted (ExtGain) or over the air (Friis)
W.Ptx = NaN; W.EIRP = NaN; W.txSource = "not available (give Distance for over the air or ExtGain for conducted)";
W.Gtx = NaN; W.Grx = NaN; W.fspl = NaN; W.gainNote = "";
if ~isempty(opt.Distance)
    lam = physconst('LightSpeed')/S.fc;
    [W.Gtx, gt] = antGain(opt.TxAntennaGain, S.fc);
    [W.Grx, gr] = antGain(opt.RxAntennaGain, S.fc);
    W.fspl = fspl(opt.Distance, lam);
    W.EIRP = W.C_dBm - W.Grx + W.fspl + opt.PolLoss;
    W.Ptx = W.EIRP - W.Gtx;
    W.txSource = sprintf("Friis: d = %g m, FSPL %.1f dB, Gtx %.1f dBi (%s), Grx %.1f dBi (%s), pol. loss %g dB", ...
        opt.Distance, W.fspl, W.Gtx, gt, W.Grx, gr, opt.PolLoss);
    if opt.Distance < 10*lam
        W.gainNote = sprintf("Distance %.2g m < 10 wavelengths (%.2g m): near-field, Friis is not reliable.", opt.Distance, 10*lam);
    end
elseif D.extGain ~= 0
    W.Ptx = W.C_dBm; W.txSource = sprintf("conducted, ExtGain %g dB", D.extGain);
end
end

function T = table7(Rs)
% SRS-4 datasheet Table 7: input power (dBm) for < 1 % PER, columns CC+RS / CC / RS / none.
rates = [100e3 1e6 4e6 5e6];
tab = {[-122 -121 -114 -112; -119 -118 -111 -109], [-112.5 -111 -104 -102; -108.5 -108 -101 -99], ...
    [-106 -105 -98 -96; -102 -101.5 NaN -93], [-105.5 -104 NaN -95; NaN -100.5 NaN -92]};
[~, i] = min(abs(log(rates/Rs)));
Rt = rates(i); p = tab{i};
rcode = [223/255/2, 1/2, 223/255, 1];
cn0 = p + 174 - 2.05;                     % NF 2.05 dB
mods = ["BPSK"; "QPSK"];
eb = cn0 - 10*log10([1; 2]*Rt*rcode);
T = table(mods, repmat(Rt/1e3, 2, 1), p(:, 1), p(:, 2), p(:, 3), p(:, 4), eb(:, 1), eb(:, 4), ...
    'VariableNames', {'Modulation', 'Rate_kBd', 'P_CCRS_dBm', 'P_CC_dBm', 'P_RS_dBm', 'P_none_dBm', ...
    'EbN0_CCRS_dB', 'EbN0_none_dB'});
end

function [g, desc] = antGain(spec, fc)
% Antenna gain in dBi: a number, or an Antenna Toolbox design at fc (boresight gain).
persistent cache
if isempty(cache), cache = containers.Map(); end
if isnumeric(spec), g = spec; desc = "given"; return, end
key = lower(string(spec)) + sprintf('_%.4g', fc);
if isKey(cache, key), v = cache(key); g = v{1}; desc = v{2}; return, end
switch lower(string(spec))
    case "patch"
        ant = design(patchMicrostrip, fc);            % air substrate, finite ground plane
        g = pattern(ant, fc, 0, 90, 'Type', 'realizedgain');
        desc = "patchMicrostrip designed at fc (air substrate; a CubeSat FR4/ceramic patch gives 5-7 dBi)";
    case "dipole"
        ant = design(dipole, fc);
        g = pattern(ant, fc, 0, 0, 'Type', 'realizedgain');
        desc = "half-wave dipole designed at fc";
    otherwise
        error('ff_srs4_analyzer:args', 'Antenna gain must be a number (dBi), "patch" or "dipole".');
end
cache(key) = {g, desc};
end

function C = instCompare(D, A)
C = struct('done', false, 'chp', NaN, 'chpML', NaN, 'diff', NaN, 'obw99', NaN, 'obw999', NaN, 'obw99ML', NaN);
if ~isfield(D, 'inst') || ~D.inst.done, return, end
x = D.inst; S = A.spec;
corr = D.extgInstr - D.extGain;           % instrument results are in display units
C.done = true;
C.chp = x.chp(1) + corr;
m = abs(S.f - (S.f(1) + S.f(end))/2) <= x.ibw/2;
C.chpML = 10*log10(sum(10.^(S.pdB(m)/10))*S.df/S.enbw);
C.diff = C.chp - C.chpML;
C.obw99 = x.obw99(1); C.obw999 = x.obw999(1); C.obw99ML = S.obw99raw;
end

function Q = freqAnalysis(S, D, opt)
Q.fc = S.fc;
ci = S.fit.bciFc; Q.ciSource = "bootstrap over sweeps";
if any(~isfinite(ci)), ci = S.fit.ciFc; Q.ciSource = "nlparci (RBW-scaled)"; end
Q.ci = ci; Q.Ustat = diff(ci)/2;
Q.refSource = D.refSource;
if ~isempty(opt.RefPpm), Q.refPpm = opt.RefPpm; Q.refNote = "given (RefPpm)";
elseif startsWith(D.refSource, "EXT"), Q.refPpm = 0.05; Q.refNote = "external 10 MHz (GPSDO assumed; give RefPpm)";
else, Q.refPpm = 0.7; Q.refNote = "FieldFox internal reference, spec +-0.7 ppm (+1 ppm/year ageing)";
end
Q.expected = opt.ExpectedFreq;
Q.nominalGuess = round(S.fc/10e3)*10e3;
Q.err = NaN; Q.ppm = NaN; Q.Uppm = NaN; Q.status05 = "INFO"; Q.status25 = "INFO";
f0 = opt.ExpectedFreq;
if isempty(f0), return, end
Q.err = S.fc - f0; Q.ppm = Q.err/f0*1e6;
Q.Uppm = Q.refPpm + Q.Ustat/f0*1e6;       % worst case: reference + statistical (95 %)
Q.status05 = judge(abs(Q.ppm), Q.Uppm, 0.5);
Q.status25 = judge(abs(Q.ppm), Q.Uppm, 2.5);
end

function s = judge(x, U, lim)
if x + U <= lim, s = "PASS";
elseif x - U > lim, s = "FAIL";
else, s = "INCONCLUSIVE";
end
end

function E = ecssObw(S)
% ECSS-E-ST-50-05C Rev.2 Table 5-5, 2200-2290 MHz, category A telemetry.
Rs = S.Rs;
if Rs < 10e3, lim = 300e3; rule = "Rs < 10 ksps: 300 kHz";
elseif Rs <= 60e3, lim = min(1200e3, 30*Rs); rule = "10-60 ksps: 1200 kHz or 30 Rs, whichever is smaller";
elseif Rs < 2e6, lim = min(6e6, max(1200e3, 12*Rs)); rule = "60 ksps-2 Msps: 1200 kHz or 12 Rs, whichever is larger, up to 6 MHz";
else, lim = min(6e6, 1.1*Rs); rule = "Rs >= 2 Msps: 1.1 Rs, up to 6 MHz";
end
E = struct('limit', lim, 'rule', rule, 'obw', S.obw99, 'margin', lim - S.obw99, 'status', "PASS");
if S.obw99 > lim, E.status = "FAIL"; end
end

function M = maskAnalysis(a, D, S, Fl)
% SFCG/ECSS spectral mask (ECSS-E-ST-50-05C Rev.2 Fig. 6-5), relative to the spectrum peak,
% linear in dB vs linear frequency offset / Rs.
f = a.f; lin = mean(10.^(levels(a, D)/10), 1); rbw = a.settings.rbw;
df = (f(end) - f(1))/(numel(f) - 1);
pdB = 10*log10(movmean(lin, max(1, round(rbw/df))));   % smoothed over one RBW: no noise spikes
pk = max(pdB);
x = abs(f - S.fc)/S.Rs;
if S.Rs >= 2e6
    xs = [0 0.5 1.4 3 1e9]; name = "suppressed carrier, Rs >= 2 Msps";
else
    xs = [0 0.5 3 8 1e9]; name = "Rs < 2 Msps (and PCM/PM/SP-L)";
end
ms = [0 0 -30 -60 -60];
mask = interp1(xs, ms, x);
noiseRel = Fl.N0dB + 10*log10(1.0645*rbw) - pk;
rel = pdB - pk;
testable = mask > noiseRel + 3;
margin = mask - rel;
viol = testable & margin < 0;
M = struct('name', name, 'x', x, 'rel', rel, 'mask', mask, 'noiseRel', noiseRel, 'testable', testable, ...
    'span', max(x), 'xEnd', xs(4), 'minMargin', NaN, 'worstX', NaN, 'status', "PASS", 'source', a.settings.label);
slope = testable & x > xs(2);   % inside +-0.5 Rs the mask is 0 dB and the peak itself defines it
if any(slope)
    [M.minMargin, i] = min(margin(slope)); xt = x(slope); M.worstX = xt(i);
end
if any(viol)
    M.status = "FAIL";
elseif ~all(testable(x <= xs(4)))
    xl = min(x(~testable));
    M.status = sprintf("PASS (to %.0f dB; noise-limited beyond +-%.2g Rs)", -noiseRel - 3, xl);
end
if max(x) < xs(4) && M.status ~= "FAIL"   % the mask was not checked where it matters most
    M.status = sprintf("INCONCLUSIVE (passes within +-%.2g Rs only; the mask runs to +-%g Rs: needs the mask sweep)", max(x), xs(4));
end
end

function H = harmAnalysis(a, a2, D, S, W)
H = struct('done', false, 'dBc', NaN, 'limit_dBc', NaN, 'detected', false, 'line', false, ...
    'internal_dBc', NaN, 'attTest', "", 'status', "NOT MEASURED", 'note', "");
if isempty(a), return, end
r = h2one(a, D, S, W);
H = mergeStruct(H, r);
H.done = true;
if H.detected
    H.status = "FAIL"; if H.dBc <= -60, H.status = "PASS"; end
else
    H.status = "PASS"; if H.limit_dBc > -60, H.status = "INCONCLUSIVE"; end
end
if H.detected && H.dBc < H.internal_dBc + 6
    H.status = "INCONCLUSIVE";            % within 6 dB of what the analyzer itself produces
    H.attTest = "level close to the analyzer's own 2nd harmonic. ";
end
if ~isempty(a2)
    % an analyzer-made 2nd harmonic falls 1 dB (relative) per dB of extra attenuation, a real one stays
    r2 = h2one(a2, D, S, W);
    da = a2.settings.atten - a.settings.atten;
    if H.detected && r2.detected
        ch = r2.dBc - H.dBc;
        if ch < -max(2, 0.5*da)
            H.attTest = H.attTest + sprintf("dropped %.1f dB with +%g dB attenuation: generated INSIDE the analyzer", -ch, da);
            H.status = "INCONCLUSIVE";
        else
            H.attTest = sprintf("changed %+.1f dB with +%g dB attenuation: a real emission of the DUT", ch, da);
            H.status = "FAIL"; if H.dBc <= -60, H.status = "PASS"; end
        end
    elseif H.detected
        H.attTest = H.attTest + sprintf("vanished with +%g dB attenuation: analyzer-made (or below the higher noise)", da);
        H.status = "INCONCLUSIVE";
    end
end
end

function r = h2one(a, D, S, W)
f = a.f; lin = mean(10.^(levels(a, D)/10), 1); rbw = a.settings.rbw; enbw = 1.0645*rbw;
df = (f(end) - f(1))/(numel(f) - 1);
x = f - 2*S.fc; bw = (1 + S.alpha)*S.Rs;   % squaring doubles the occupied bandwidth
inner = abs(x) <= bw; outer = abs(x) >= 1.6*bw;
if nnz(outer) < 10, outer = abs(x) >= 0.8*max(abs(x)); end
Nb = median(lin(outer));
ex = sum(lin(inner) - Nb)*df/enbw;
nIndep = max(1, nnz(inner)*df/enbw);
sig = std(lin(outer))*sqrt(nIndep);          % RBW-correlated bins: nIndep independent groups
sig = max(sig, Nb*nnz(inner)*df/enbw*0.02);        % 2 % floor: noise-level estimate error
r.detected = ex > 3*sig;
r.dBc = 10*log10(max(ex, realmin)) - W.C_dBm;
r.limit_dBc = 10*log10(3*sig) - W.C_dBm;
pd = 10*log10(lin);
r.line = max(pd(abs(x) < 3*rbw)) - median(pd(inner)) > 6 && r.detected;
mixer = W.C_dBm + D.extGain - a.settings.atten;     % input-referred level minus attenuation
r.internal_dBc = mixer - 30;                         % SHI +30 dBm (spec <= 4 GHz)
r.atten = a.settings.atten;
r.note = sprintf("mixer %.0f dBm -> analyzer's own 2nd harmonic about %.0f dBc", mixer, r.internal_dBc);
end

function W2 = wideAnalysis(a, D, S, W, opt)
W2 = struct('done', false, 'spurs', table(), 'worst_dBc', NaN, 'limit_dBc', NaN, 'status', "NOT MEASURED", 'note', "");
if isempty(a), return, end
f = a.f; pk = max(levels(a, D), [], 1);
df = (f(end) - f(1))/(numel(f) - 1);
fl = movmedian(pk, max(5, round(100e6/df)));
ex = abs(f - S.fc) < 2.5*max(S.obw99, S.Rs) + 2*a.settings.rbw;   % fundamental + out-of-band domain
d = pk - fl; d(ex) = 0;
[~, loc] = findpeaks(d, 'MinPeakHeight', 10, 'MinPeakProminence', 8, 'MinPeakDistance', 3);
rows = {};
for i = loc
    rows(end + 1, :) = {f(i), pk(i), pk(i) - W.C_dBm, classify(f(i), S, opt)}; %#ok<AGROW>
end
if isempty(rows)
    T = table('Size', [0 4], 'VariableTypes', {'double', 'double', 'double', 'string'}, ...
        'VariableNames', {'Freq_Hz', 'Level_dBm', 'Rel_dBc', 'Likely'});
else
    T = cell2table(rows, 'VariableNames', {'Freq_Hz', 'Level_dBm', 'Rel_dBc', 'Likely'});
    T = sortrows(T, 'Rel_dBc', 'descend');
end
W2.spurs = T; W2.done = true;
W2.limit_dBc = median(fl(~ex)) + 10 - W.C_dBm;
own = ~startsWith(T.Likely, "ambient");
if any(own), W2.worst_dBc = max(T.Rel_dBc(own)); end
if isfinite(W2.worst_dBc) && W2.worst_dBc > -60
    W2.status = "FAIL";
elseif W2.limit_dBc > -60
    W2.status = "INCONCLUSIVE";
else
    W2.status = "PASS";
end
W2.note = sprintf("RBW %s peak detector, max hold of %d sweeps; detection limit %.0f dBc", fmtHz(a.settings.rbw), ...
    size(a.P, 1), W2.limit_dBc);
end

function s = classify(f, S, opt)
if abs(f - 2*S.fc) < 2*S.Rs, s = "2nd harmonic"; return, end
if f >= opt.Band(1) && f <= opt.Band(2), s = "in band (DUT spur or another emitter)"; return, end
B = [87.5 108 "ambient: FM broadcast"; 470 790 "ambient: DVB-T2 TV"; 703 960 "ambient: LTE/GSM/NR 700-900";
    1427 1518 "ambient: LTE L-band"; 1710 1880 "ambient: LTE/GSM 1800"; 1920 2170 "ambient: UMTS/LTE 2100";
    2300 2400 "ambient: LTE 2300"; 2400 2483.5 "ambient: Wi-Fi/ISM 2.4 GHz"; 2496 2690 "ambient: LTE/NR 2600";
    3300 3800 "ambient: NR 3.5 GHz"; 5150 5850 "ambient: Wi-Fi 5 GHz"];
for i = 1:size(B, 1)
    if f >= str2double(B(i, 1))*1e6 && f <= str2double(B(i, 2))*1e6, s = B(i, 3); return, end
end
s = "unknown (check with the DUT off; change attenuation to test for analyzer spurs)";
end

function V = envelopeAnalysis(a, D, S, W, Fl, B)
% Zero-span envelope statistics vs BPSK/QPSK/8PSK hypotheses simulated with the generator.
V = struct('done', false, 'status', "NOT MEASURED", 'decision', "not measured", 'table', table(), 'post', [NaN NaN NaN], 'note', "");
if isempty(a), return, end
P = levels(a, D); [ns, N] = size(P);
rbw = a.settings.rbw; T = a.settings.sweepTime; dt = T/(N - 1);
frac = rbwFraction(S, rbw);
Pon = 10*log10(10^(W.C_dBm/10)*frac + Fl.N0*1.0645*rbw);
Pn = Fl.N0dB + 10*log10(1.0645*rbw);
V.noiseFrom = "spectrum noise floor";
if isfield(B, 'noiseMeas') && isfinite(B.noiseMeas) && B.rbw == rbw
    Pn = B.noiseMeas; V.noiseFrom = "zero-span TX-off level";   % measured, not predicted
end
feats = NaN(ns, 3); keep = false(ns, 1); samples = {};
for k = 1:ns
    y = P(k, :);
    good = movmedian(y, 51) > Pon - 6;    % drop the parts where the transmitter was off
    if nnz(good) < 0.5*N, continue, end
    p = 10.^(y(good)/10); p = p/mean(p);
    feats(k, :) = [mean(p < 0.1), std(10*log10(p)), 10*log10(prctile(p, 99.99))];
    keep(k) = true; samples{end + 1} = p; %#ok<AGROW>
end
V.nUsed = nnz(keep); V.nSweeps = ns;
snr = 10^((W.C_dBm + 10*log10(frac) - Pn)/10);
V.snr_dB = 10*log10(snr); V.rbw = rbw; V.dt = dt; V.Pon = Pon; V.Pn = Pn;
if V.nUsed == 0
    V.status = "INCONCLUSIVE"; V.note = "No zero-span sweep with the transmitter on for most of the sweep.";
    return
end
m = mean(feats(keep, :), 1); se = std(feats(keep, :), 0, 1)/sqrt(V.nUsed);
if V.nUsed < 2, se = [0.01 0.2 0.3]; end
V.meas = m; V.measSE = se; V.samples = [samples{:}];
mods = ["BPSK" "QPSK" "8PSK"];
H = zeros(3, 3); Hs = zeros(3, 3); V.simSamples = cell(1, 3);
for i = 1:3
    [H(i, :), Hs(i, :), V.simSamples{i}] = simEnvelope(mods(i), S.Rs, S.alpha, rbw, snr, dt, N);
end
floorErr = [0.01 0.3 0.4];               % model error: FieldFox detector/IF details unknown
ll = zeros(1, 3);
for i = 1:3
    sg = sqrt(se.^2 + Hs(i, :).^2 + floorErr.^2);
    ll(i) = -0.5*sum(((m(1:2) - H(i, 1:2))./sg(1:2)).^2) - sum(log(sg(1:2)));
end
post = exp(ll - max(ll)); post = post/sum(post);
V.post = post;
V.table = table(["measured"; mods.'], [m(1); H(:, 1)]*100, [m(2); H(:, 2)], [m(3); H(:, 3)], [NaN; post.'], ...
    'VariableNames', {'Case', 'Below10dB_pct', 'StdDbEnv_dB', 'PAPR_1e4_dB', 'Posterior'});
if m(2) < 1 && m(1) < 0.002
    V.decision = "constant envelope (GMSK/FSK-like): not an RRC-PSK signal";
elseif post(1) > 0.9
    V.decision = sprintf("BPSK (posterior %.2f)", post(1));
elseif post(2) + post(3) > 0.9
    V.decision = sprintf("QPSK or 8PSK (posterior %.2f; QPSK %.2f vs 8PSK %.2f is a weak distinction)", ...
        post(2) + post(3), post(2), post(3));
else
    V.decision = sprintf("inconclusive (BPSK %.2f, QPSK %.2f, 8PSK %.2f)", post);
end
[~, ib] = max(post);
if m(3) < H(ib, 3) - 0.7
    V.decision = V.decision + sprintf("; PAPR %.1f dB is below the ideal %.1f dB: peaks compressed (PA near saturation)", m(3), H(ib, 3));
end
V.status = "INFO"; V.done = true;
V.note = sprintf("%d of %d sweeps used, SNR in RBW %.1f dB (noise: %s), RBW %s, sample spacing %.3g us = %.3f symbols", ...
    V.nUsed, ns, V.snr_dB, V.noiseFrom, fmtHz(rbw), dt*1e6, dt*S.Rs);
if abs(dt*S.Rs - round(dt*S.Rs)) < 0.05 && dt*S.Rs >= 0.95
    V.note = V.note + " - WARNING: whole-symbol spacing strobes one symbol phase";
end
if rbw < (1 + S.alpha)*S.Rs
    V.note = V.note + "; RBW narrower than the signal: the envelope is smoothed (simulated the same way)";
end
end

function fr = rbwFraction(S, rbw)
% Fraction of the signal power inside a Gaussian RBW centred on fc.
fx = linspace(-1, 1, 4001)*(1 + S.alpha)*S.Rs/2;
s = rcShape(fx, 0, S.Rs, S.alpha);
fr = sum(s.*exp(-4*log(2)*(fx/rbw).^2))/sum(s);
end

function [mu, sd, smp] = simEnvelope(mod, Rs, alpha, rbw, snr, dt, N)
% Envelope features of a CCSDS waveform (RS coding, RRC) + noise through a Gaussian RBW,
% sampled every dt like the zero-span trace; blocks of N samples = one "sweep".
persistent cache
if isempty(cache) || cache.Count >= 24, cache = containers.Map(); end   % keys are measured values: bounded
key = sprintf('%s_%.4g_%.3f_%.4g_%.3g_%.3g_%d', mod, Rs, alpha, rbw, 10*log10(snr), dt, N);
if isKey(cache, key), v = cache(key); [mu, sd, smp] = v{:}; return, end
sps = min(32, max(4, 2*ceil(1.5*max(rbw, (1 + alpha)*Rs)/Rs)));
fs = sps*Rs;
tm = ccsdsTMWaveformGenerator("WaveformSource", "synchronization and channel coding", ...
    "HasRandomizer", true, "HasASM", true, "ChannelCoding", "RS", "RSMessageLength", 223, ...
    "RSInterleavingDepth", 1, "Modulation", mod, "PulseShapingFilter", "root raised cosine", ...
    "RolloffFactor", min(max(alpha, 0.05), 1), "FilterSpanInSymbols", 10, "SamplesPerSymbol", sps);
rs = RandStream('mt19937ar', 'Seed', 5);
x = zeros(0, 1);
while numel(x) < 3e5, x = [x; tm(randi(rs, [0 1], tm.NumInputBits, 1))]; end %#ok<AGROW>
x = x.'/rms(x);
enbw = 1.0645*rbw;
nz = sqrt(fs/(snr*enbw)/2)*(randn(rs, size(x)) + 1i*randn(rs, size(x)));
y = gaussFilter(x, rbw, fs) + gaussFilter(nz, rbw, fs);
% sample at the measured spacing with a fractional start: no artificial symbol strobing
tq = 1 + rand(rs)*sps + (0:floor((numel(y) - 2 - sps)/(dt*fs)))*dt*fs;
p = abs(interp1(y, tq, 'linear')).^2;
nb = max(3, floor(numel(p)/min(N, numel(p))));
L = floor(numel(p)/nb);
F = zeros(nb, 3);
for b = 1:nb
    q = p((b - 1)*L + (1:L)); q = q/mean(q);
    F(b, :) = [mean(q < 0.1), std(10*log10(q)), 10*log10(prctile(q, 99.99))];
end
mu = mean(F, 1); sd = std(F, 0, 1)/sqrt(nb);
smp = p(1:min(end, 2e4))/mean(p);
cache(key) = {mu, sd, smp};
end

function B = burstAnalysis(a, D, S, W, Fl)
B = struct('done', false, 'state', "not measured", 'burstLen', NaN, 'gap', NaN, 'period', NaN, ...
    'duty', NaN, 'nBursts', 0, 'status', "NOT MEASURED", 'note', "", 'lev', [NaN NaN], 'thr', NaN, ...
    'noiseMeas', NaN, 'rbw', NaN);
if isempty(a), return, end
B.rbw = a.settings.rbw;
P = levels(a, D); [ns, N] = size(P);
T = a.settings.sweepTime; fs = (N - 1)/T;
frac = rbwFraction(S, a.settings.rbw);
Pon = 10*log10(10^(W.C_dBm/10)*frac + Fl.N0*1.0645*a.settings.rbw);
Pn = Fl.N0dB + 10*log10(1.0645*a.settings.rbw);
B.Pon = Pon; B.Pn = Pn; B.T = T; B.done = true;
if Pon - Pn < 10
    B.state = "contrast too low"; B.status = "INCONCLUSIVE";
    B.note = sprintf("expected on/off contrast only %.1f dB in the zero-span RBW", Pon - Pn); return
end
w = []; sep = []; per = []; dc = []; on = 0; offLev = [];
for k = 1:ns
    x = P(k, :);                          % dB: a linear histogram cannot resolve the noise state
    lv = statelevels(x);
    if k == 1, B.lev = lv; end
    % noise-calibrated: the low state must sit on the analyzer noise and the high state at
    % the expected channel level, otherwise "two levels" are just envelope fluctuations
    bursty = lv(2) - lv(1) > 10 && lv(1) < Pn + 6 && lv(2) > Pon - 6;
    if bursty
        offLev(end + 1) = 10*log10(mean(10.^(x(x < mean(lv))/10))); %#ok<AGROW>   % TX-off noise (mean power)
        try
            row = @(v) reshape(v, 1, []);
            w = [w, row(pulsewidth(x, fs, 'StateLevels', lv))]; %#ok<AGROW>
            sep = [sep, row(pulsesep(x, fs, 'StateLevels', lv))]; %#ok<AGROW>
            per = [per, row(pulseperiod(x, fs, 'StateLevels', lv))]; %#ok<AGROW>
            dc = [dc, row(dutycycle(x, fs, 'StateLevels', lv))]; %#ok<AGROW>
        catch
        end
    elseif median(P(k, :)) > Pon - 6
        on = on + 1;
    end
end
B.thr = (Pn + Pon)/2;
if ~isempty(w)
    B.state = "bursts"; B.burstLen = median(w); B.nBursts = numel(w);
    B.gap = median(sep); B.period = median(per); B.duty = 100*median(dc);
    B.note = sprintf("%d complete bursts in %d sweeps of %.3g s; burst %.3f ms (spread %.3f ms), gap %.3f ms", ...
        numel(w), ns, T, 1e3*B.burstLen, 1e3*(max(w) - min(w)), 1e3*B.gap);
    if ~isempty(offLev)
        B.noiseMeas = median(offLev);
        B.note = B.note + sprintf("; TX-off noise %.1f dBm vs %.1f dBm predicted from the spectrum", B.noiseMeas, Pn);
    end
elseif on == ns
    B.state = "continuous"; B.note = sprintf("TX on for all %d sweeps of %.3g s (idle frames between data?)", ns, T);
else
    B.state = "no complete burst"; B.note = "bursts longer than the sweep, or TX off";
end
B.status = "INFO";
B.resolution = T/(N - 1);
end

function [Rs, src, relU] = frameRate(A, opt)
% The most exact symbol rate for frame-time arithmetic: idle-frame lines, else the sl-radio
% setting (exact) when the fit agrees within 2 %, else the fit with its uncertainty.
S = A.spec;
if A.lines.found && abs(A.lines.RsLines/S.Rs - 1) < 0.05
    Rs = A.lines.RsLines; src = "idle-frame lines"; relU = 0.001;
elseif ~isempty(opt.ExpectedRate) && abs(S.Rs/opt.ExpectedRate - 1) < 0.02
    Rs = opt.ExpectedRate; src = "ExpectedRate (fit agrees within 2 %)"; relU = 0;
else
    Rs = S.Rs; src = "RC fit"; relU = max(abs(S.fit.bciRs - S.Rs))/S.Rs;
    if ~isfinite(relU), relU = 0.01; end
end
end

function T = frameTable(Rs, relU, B)
% SRS-4 frame: 2-byte header + 217-byte payload + CRC32C = 223 bytes; RS adds 32 parity;
% the ASM (4 bytes) is added after the randomizer; CC r=1/2 doubles everything incl. the ASM.
mods = ["BPSK" "QPSK" "8PSK"]; M = [1 2 3];
cod = ["none (CRC only)" "RS(255,223)" "CC r=1/2" "RS + CC"]; bits = [227 259 454 518]*8;
rows = {};
for i = 1:3
    for j = 1:4
        Tf = bits(j)/(M(i)*Rs);
        nF = NaN; ok = "";
        if isfield(B, 'burstLen') && isfinite(B.burstLen)
            nF = B.burstLen/Tf;
            tol = max(0.1, 3*B.resolution/Tf) + nF*relU;     % timing resolution + Rs uncertainty
            ok = string(abs(nF - round(nF)) <= tol && round(nF) >= 2);
        end
        rows(end + 1, :) = {mods(i), cod(j), bits(j), 1e3*Tf, 1/Tf, 1736/Tf/1e3, nF, ok}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'Modulation', 'Coding', 'BitsPerFrame', 'Tframe_ms', ...
    'LineSpacing_Hz', 'Payload_kbps', 'FramesInBurst', 'WholeFrames'});
end

function L = lineAnalysis(a, D, S)
L = struct('done', false, 'found', false, 'spacing', NaN, 'Tframe', NaN, 'match', "", 'RsLines', NaN, ...
    'status', "NOT MEASURED", 'note', "", 'loc', []);
if isempty(a), return, end
f = a.f; pdB = 10*log10(mean(10.^(levels(a, D)/10), 1));
df = (f(end) - f(1))/(numel(f) - 1); rbw = a.settings.rbw;
bits = [1816 2072 3632 4144]; M = [1 2 3];
spMin = S.Rs/max(bits); spMax = 3*S.Rs/min(bits);
base = movmedian(pdB, max(5, round(3*spMax/df)));
d = pdB - base;
[~, loc] = findpeaks(d, 'MinPeakProminence', 4, 'MinPeakDistance', max(1, round(0.6*spMin/df)));
L.done = true; L.status = "INFO"; L.d = d; L.loc = loc;
if numel(loc) < 5
    L.note = "no regular lines: data traffic, encryption or the RBW too wide (lines appear only for repeated idle frames)";
    return
end
x = f(loc);
s0 = median(diff(x));
k = round((x - x(1))/s0);
p = polyfit(k, x, 1);
r = x - polyval(p, k);
good = abs(r) < 0.25*p(1);
if nnz(good) >= 5, p = polyfit(k(good), x(good), 1); end
L.found = true; L.spacing = p(1); L.Tframe = 1/p(1);
mods = ["BPSK" "QPSK" "8PSK"]; cods = ["no FEC" "RS" "CC" "RS+CC"];
[ii, jj] = ndgrid(1:3, 1:4);
e = abs(log(L.spacing./(M(ii).*S.Rs./bits(jj))));
[L.err, b] = min(e(:));
L.RsLines = L.spacing*bits(jj(b))/M(ii(b));
same = find(e(:) <= L.err + 0.005);       % equal bits/symbol per frame cannot be told apart
L.match = strjoin(compose("%s %s", mods(ii(same)).', cods(jj(same)).'), " = ");
L.note = sprintf("%d lines, spacing %.2f Hz -> T_frame %.4f ms; matches %s (%.1f %% off) -> Rs %s", ...
    numel(loc), L.spacing, 1e3*L.Tframe, L.match, 100*(exp(L.err) - 1), fmtHz(L.RsLines));
if rbw > L.spacing/3, L.note = L.note + " (RBW marginal)"; end
end

function T = limitsTable(A, D, opt)
S = A.spec; W = A.power; Q = A.freq;
r = limitRow("Centre frequency (fit)", sprintf("%s [%s ... %s]", fmtHz6(S.fc), fmtHz6(Q.ci(1)), fmtHz6(Q.ci(2))), ...
    "2200-2290 MHz", "PASS", "95 % CI: " + Q.ciSource);
if S.fc - S.obw99/2 < opt.Band(1) || S.fc + S.obw99/2 > opt.Band(2), r.Status = "FAIL"; end
r.Limit = fmtRange(opt.Band);
T = r;
% symbol rate
ci = S.fit.bciRs; if any(~isfinite(ci)), ci = S.fit.ciRs; end
st = "INFO"; lim = "10 kBd - 5 MBd"; note = sprintf("-3 dB width %s", fmtHz(S.Rs3));
if ~isempty(opt.ExpectedRate)
    tol = max(0.02*opt.ExpectedRate, max(abs(ci - S.Rs)));
    st = "PASS"; if abs(S.Rs - opt.ExpectedRate) > tol, st = "FAIL"; end
    lim = sprintf("%s (expected) +- %.1f %%", fmtHz(opt.ExpectedRate), 100*tol/opt.ExpectedRate);
end
T = [T; limitRow("Symbol rate Rs (fit)", sprintf("%s Bd [%s ... %s]", fmtNum(S.Rs), fmtNum(ci(1)), fmtNum(ci(2))), lim, st, note)];
% roll-off
ci = S.fit.bciAlpha; if any(~isfinite(ci)), ci = S.fit.ciAlpha; end
st = "INFO"; lim = "0.20 / 0.25 / 0.35";
if ~isempty(opt.ExpectedRolloff)
    st = "PASS"; if abs(S.alpha - opt.ExpectedRolloff) > max(0.03, max(abs(ci - S.alpha))), st = "FAIL"; end
    lim = sprintf("%.2f (expected) +- 0.03", opt.ExpectedRolloff);
end
T = [T; limitRow("Roll-off alpha", sprintf("fit %.3f [%.3f ... %.3f], OBW lookup %.3f, template %.2f", S.alpha, ci(1), ci(2), ...
    S.alphaObw, A.ref.best), lim, st, "template = best ccsdsTMWaveformGenerator reference")];
% occupied bandwidth (datasheet values at alpha 0.2)
aRef = A.ref.best; if ~isfinite(aRef), aRef = S.alpha; end
is02 = aRef == 0.2;
st = "INFO"; if is02, st = pf(S.obw99/S.Rs <= 1.07*1.02); end
if S.obw99Limited, st = "INCONCLUSIVE"; end
T = [T; limitRow("OBW 99 %", sprintf("%s = %.3f Rs", fmtHz(S.obw99), S.obw99/S.Rs), "1.07 Rs (alpha 0.2, +2 % meas.)", st, ...
    sprintf("ideal RC: %.3f Rs at alpha %.2f", idealObw(aRef, 99), aRef))];
st = "INFO"; if is02, st = pf(S.obw999/S.Rs <= 1.18*1.02); end
if S.obw999Limited, st = "INCONCLUSIVE"; end
T = [T; limitRow("OBW 99.9 %", sprintf("%s = %.3f Rs", fmtHz(S.obw999), S.obw999/S.Rs), "1.18 Rs (alpha 0.2, +2 % meas.)", st, ...
    sprintf("ideal RC: %.3f Rs; more = PA regrowth", idealObw(aRef, 99.9)))];
st = "INFO"; if S.ebw26Limited, st = "INCONCLUSIVE"; end
T = [T; limitRow("-26 dB emission bandwidth", sprintf("%s = %.3f Rs", fmtHz(S.ebw26), S.ebw26/S.Rs), "-", st, "")];
% frequency error
if isempty(Q.expected)
    T = [T; limitRow("Frequency error", "ExpectedFreq not given", "0.5 ppm (20 C) / 2.5 ppm", "INFO", ...
        sprintf("if tx.freq = %s: %+.0f Hz = %+.2f ppm (guess)", fmtHz6(Q.nominalGuess), S.fc - Q.nominalGuess, ...
        (S.fc - Q.nominalGuess)/Q.nominalGuess*1e6))];
else
    m = sprintf("%+.0f Hz = %+.3f ppm (+- %.2f ppm)", Q.err, Q.ppm, Q.Uppm);
    T = [T; limitRow("Frequency error, initial", m, "<= 0.5 ppm", Q.status05, Q.refNote)];
    T = [T; limitRow("Frequency error, over temperature", m, "<= 2.5 ppm", Q.status25, "")];
end
% power
if isfinite(W.Ptx)
    u = 0.8; if ~isempty(opt.Distance), u = 3; end
    st = "PASS"; if W.Ptx + u < 20 || W.Ptx - u > 33, st = "FAIL"; elseif W.Ptx - u < 20 || W.Ptx + u > 33, st = "INCONCLUSIVE"; end
    T = [T; limitRow("TX output power", sprintf("%.1f dBm (+- %.1f dB)", W.Ptx, u), "20-33 dBm", st, W.txSource)];
else
    T = [T; limitRow("TX output power", sprintf("%.1f dBm at the analyzer", W.C_dBm), "20-33 dBm", "INFO", W.txSource)];
end
% spectral limits
T = [T; limitRow("ECSS OBW limit (Table 5-5)", fmtHz(A.ecss.obw), fmtHz(A.ecss.limit), A.ecss.status, A.ecss.rule)];
T = [T; limitRow("SFCG/ECSS spectral mask", sprintf("min margin %.1f dB at %.2f Rs", A.mask.minMargin, A.mask.worstX), ...
    A.mask.name, firstWord(A.mask.status), A.mask.status)];
H = A.harm;
if H.done
    m = sprintf("< %.1f dBc (not detected)", H.limit_dBc); if H.detected, m = sprintf("%.1f dBc", H.dBc); end
    T = [T; limitRow("2nd harmonic", m, "<= -60 dBc", H.status, H.note + ". " + H.attTest)];
else
    T = [T; limitRow("2nd harmonic", "not measured", "<= -60 dBc", "NOT MEASURED", "")];
end
X = A.wide;
if X.done
    m = sprintf("worst %.1f dBc", X.worst_dBc); if ~isfinite(X.worst_dBc), m = "none above the detection limit"; end
    T = [T; limitRow("Spurious (SFDR)", m, ">= 60 dBc", X.status, X.note)];
else
    T = [T; limitRow("Spurious (SFDR)", "not measured", ">= 60 dBc", "NOT MEASURED", "")];
end
% other information
V = A.env;
T = [T; limitRow("Modulation (envelope)", V.decision, "BPSK/QPSK/8PSK", V.status, V.note)];
B = A.burst; m = B.state; if isfinite(B.burstLen), m = sprintf("burst %.3f ms", 1e3*B.burstLen); end
T = [T; limitRow("Burst / frame timing", m, "n x T_frame", B.status, B.note)];
L = A.lines;
if L.done, T = [T; limitRow("Idle-frame lines", L.note, "1/T_frame", L.status, "")]; end
C = A.inst;
if C.done
    T = [T; limitRow("FieldFox CHP vs MATLAB", sprintf("%.2f vs %.2f dBm (%+.2f dB)", C.chp, C.chpML, C.diff), ...
        "agree within 0.5 dB", pf(abs(C.diff) <= 0.5), sprintf("FieldFox OBW99 %s vs MATLAB %s (noise incl.)", ...
        fmtHz(C.obw99), fmtHz(C.obw99ML)))];
end
ov = false;
for nm = string(fieldnames(D.acq)).'
    ov = ov || D.acq.(nm).settings.overload;
end
T = [T; limitRow("Analyzer overload (WAOR)", string(ov), "none", pf(~ov), "")];
end

function r = limitRow(p, m, l, s, n)
r = table(string(p), string(m), string(l), string(s), string(n), ...
    'VariableNames', {'Parameter', 'Measured', 'Limit', 'Status', 'Note'});
end

function s = pf(tf)
if tf, s = "PASS"; else, s = "FAIL"; end
end

function w = firstWord(s)
w = extractBefore(s + " ", " ");
end

function r = idealObw(a, p)
fn = linspace(-1.1, 1.1, 22001);
c = cumsum(rcShape(fn, 0, 1, a)); c = c/c(end) + (1:numel(c))*1e-13;
r = interp1(c, fn, 0.5 + p/200) - interp1(c, fn, 0.5 - p/200);
end

function txt = conclusions(A, D, ~)
S = A.spec; W = A.power; Q = A.freq;
txt = strings(0, 1);
if D.source == "demo"
    txt(end + 1) = sprintf("SYNTHETIC DEMO: truth fc %s, Rs %s, alpha %.2f, %s - compare with the estimates below.", ...
        fmtHz6(D.truth.fc), fmtHz(D.truth.Rs), D.truth.alpha, D.truth.mod);
end
txt(end + 1) = sprintf("Downlink at %s, Rs = %s (-3 dB width %s), roll-off %.3f (fit), %.2f (best reference).", ...
    fmtHz6(S.fc), fmtHz(S.Rs), fmtHz(S.Rs3), S.alpha, A.ref.best);
txt(end + 1) = sprintf("OBW99 %s = %.3f Rs, OBW99.9 %s = %.3f Rs (datasheet at alpha 0.2: 1.07 / 1.18 Rs).", ...
    fmtHz(S.obw99), S.obw99/S.Rs, fmtHz(S.obw999), S.obw999/S.Rs);
txt(end + 1) = sprintf("C = %.1f dBm, N0 = %.1f dBm/Hz (analyzer, NF about %.0f dB) -> C/N0 %.1f dB-Hz, Es/N0 %.1f dB.", ...
    W.C_dBm, W.N0_dBmHz, W.nfAnalyzer, W.CN0, W.EsN0);
if A.floor.upperBound
    txt(end) = txt(end) + " N0 came from the zoom edges (no mask sweep): an upper bound, so C/N0 is a lower bound.";
end
if isfinite(W.Ptx), txt(end + 1) = sprintf("Transmitter power %.1f dBm (%s).", W.Ptx, W.txSource); end
if isempty(Q.expected)
    txt(end + 1) = sprintf("Frequency error not judged (no ExpectedFreq); the FieldFox reference alone is +-%.2g ppm = +-%.0f Hz.", ...
        Q.refPpm, Q.refPpm*S.fc*1e-6);
else
    txt(end + 1) = sprintf("Frequency error %+.0f Hz (%+.3f ppm), uncertainty +-%.2f ppm: %s vs 0.5 ppm, %s vs 2.5 ppm.", ...
        Q.err, Q.ppm, Q.Uppm, Q.status05, Q.status25);
end
txt(end + 1) = "Mask: " + A.mask.status + "; ECSS OBW " + A.ecss.status + ".";
if A.harm.done, txt(end + 1) = "2nd harmonic: " + A.harm.status + " (" + A.harm.note + ")."; end
if A.wide.done, txt(end + 1) = "Spurious: " + A.wide.status + " (" + A.wide.note + ")."; end
if A.env.done, txt(end + 1) = "Envelope statistics suggest " + A.env.decision + "."; end
if A.burst.done, txt(end + 1) = "Zero span: " + A.burst.state + ". " + A.burst.note; end
if isfield(A.burst, 'burstLen') && isfinite(A.burst.burstLen)
    k = A.frames.WholeFrames == "true";
    if any(k)
        txt(end + 1) = "Burst length is a whole number of frames for: " + strjoin(A.frames.Modulation(k) + " " + A.frames.Coding(k), ", ") + ".";
    end
end
if A.lines.found, txt(end + 1) = "Idle-frame lines: " + A.lines.note + "."; end
if isfinite(A.mask.noiseRel) && A.mask.noiseRel > -60
    txt(end + 1) = sprintf("Dynamic range: the floor (analyzer noise or PA regrowth) is only %.0f dB below the peak, so -60 dB points cannot be verified here (conducted test with a pad, or a notch filter, needed).", -A.mask.noiseRel);
end
txt(end + 1) = "No I/Q: EVM and constellation are not measurable with this FieldFox.";
end

%% ===== figures ===================================================================

function files = makeFigures(D, A, opt)
files = strings(0, 1);
if ~A.ok
    if isfield(D.acq, 'coarse'), files(end + 1, 1) = figSearch(D, A, opt); end
    return
end
files(end + 1, 1) = figSpectrum(D, A, opt);
files(end + 1, 1) = figFit(A, opt);
if isfield(D.acq, 'coarse'), files(end + 1, 1) = figSearch(D, A, opt); end
if A.harm.done || A.wide.done, files(end + 1, 1) = figSpurs(D, A, opt); end
if A.env.done || A.burst.done, files(end + 1, 1) = figEnvelope(D, A, opt); end
if A.lines.done, files(end + 1, 1) = figLines(D, A, opt); end
end

function file = figSpectrum(D, A, opt)
S = A.spec; M = A.mask;
fig = newFig(opt, [80 60 1100 820]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);
x = (S.f - S.fc)/1e3;
plot(ax, x, S.pdB, 'Color', [0 0.3 0.8], 'DisplayName', sprintf('measured (%d sweeps, RBW %s)', S.nSweeps, fmtHz(S.rbw)));
hold(ax, 'on');
plot(ax, x, S.fit.model, 'k--', 'LineWidth', 1.2, 'DisplayName', sprintf('RC fit: Rs %s, \\alpha %.3f', fmtHz(S.Rs), S.alpha));
if ~isempty(A.ref.curve)
    plot(ax, x, A.ref.curve, ':', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5, ...
        'DisplayName', sprintf('ccsdsTMWaveformGenerator QPSK RS+CC, \\alpha %.2f', A.ref.best));
end
yline(ax, S.Nbin_dB, ':', 'Color', [0.5 0.5 0.5], 'DisplayName', 'noise floor');
v = [S.f3lo S.f3hi]; v = v(isfinite(v)); if ~isempty(v), xline(ax, (v - S.fc)/1e3, 'g-', 'HandleVisibility', 'off'); end
xline(ax, ([S.obw99lo S.obw99hi] - S.fc)/1e3, 'm--', 'HandleVisibility', 'off');
plot(ax, NaN, NaN, 'g-', 'DisplayName', sprintf('-3 dB: %s', fmtHz(S.Rs3)));
plot(ax, NaN, NaN, 'm--', 'DisplayName', sprintf('OBW99: %s', fmtHz(S.obw99)));
grid(ax, 'on'); xlabel(ax, sprintf('Offset from fc = %s (kHz)', fmtHz6(S.fc))); ylabel(ax, 'dBm in RBW (DUT-referred)');
ttl = 'SRS-4 downlink spectrum'; if D.source == "demo", ttl = [ttl ' (SYNTHETIC DEMO)']; end
title(ax, ttl); legend(ax, 'Location', 'southoutside', 'NumColumns', 3, 'FontSize', 8);
yl = [S.Nbin_dB - 5, S.top_dB + 8]; if all(isfinite(yl)) && yl(2) > yl(1), ylim(ax, yl); end
ax = nexttile(tl);
[xs, i] = sort(M.x .* sign(getMaskF(D, A) - S.fc));
plot(ax, xs, M.rel(i), 'Color', [0 0.3 0.8], 'DisplayName', 'measured (rel. to peak)');
hold(ax, 'on');
plot(ax, xs, M.mask(i), 'r', 'LineWidth', 1.5, 'DisplayName', "mask: " + M.name);
yline(ax, M.noiseRel, ':', 'Color', [0.5 0.5 0.5], 'DisplayName', 'noise floor');
grid(ax, 'on'); xlabel(ax, 'Offset / Rs'); ylabel(ax, 'dB relative to peak');
title(ax, "SFCG/ECSS mask: " + M.status, 'FontSize', 9);
legend(ax, 'Location', 'southoutside', 'NumColumns', 3, 'FontSize', 8);
ylim(ax, [max(-90, M.noiseRel - 10), 5]);
xe = min(max(M.x), M.xEnd*1.1); if isfinite(xe) && xe > 0, xlim(ax, [-1 1]*xe); end
file = saveFig(fig, fullfile(opt.OutDir, "fig_spectrum.png"), opt);
end

function f = getMaskF(D, A)
a = getAcq(D, "mask"); if isempty(a), a = getAcq(D, "zoom"); end
f = a.f; if numel(f) ~= numel(A.mask.x), f = A.spec.f; end
end

function file = figFit(A, opt)
S = A.spec;
fig = newFig(opt, [100 80 1100 700]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl, [1 3]);
x = (S.f - S.fc)/S.Rs;
plot(ax, x, S.fit.resid, 'k', 'DisplayName', sprintf('measured - RC fit (rms %.2f dB in fit region)', S.fit.rms));
hold(ax, 'on');
if ~isempty(A.ref.resid)
    plot(ax, x, A.ref.resid, 'Color', [0.85 0.33 0.1], 'DisplayName', ...
        sprintf('measured - reference \\alpha %.2f (rms %.2f dB)', A.ref.best, min(A.ref.rms)));
end
u = S.fit.use; plot(ax, x(u), zeros(1, nnz(u)) - 3, '.', 'Color', [0.3 0.7 0.3], 'MarkerSize', 3, 'DisplayName', 'fit region');
grid(ax, 'on'); ylim(ax, [-4 4]); xlabel(ax, 'Offset / Rs'); ylabel(ax, 'Residual (dB)');
legend(ax, 'Location', 'northeast', 'FontSize', 8); title(ax, 'Residuals: model mismatch shows where the real filter/PA differs');
B = S.fit.boot;
lab = {'fc - fit (Hz)', 'Rs (kBd)', 'roll-off \alpha'};
for k = 1:3
    ax = nexttile(tl);
    if size(B, 1) >= 10
        v = B(:, k); if k == 1, v = v - S.fc; elseif k == 2, v = v/1e3; end
        histogram(ax, v, 25, 'FaceColor', [0.3 0.5 0.8]);
    else
        text(ax, 0.1, 0.5, 'no bootstrap (fewer than 3 sweeps or Boot = 0)', 'Units', 'normalized');
    end
    xlabel(ax, lab{k}); grid(ax, 'on');
end
title(tl, sprintf('Bootstrap over sweeps (%d resamples)', size(B, 1)));
file = saveFig(fig, fullfile(opt.OutDir, "fig_fit.png"), opt);
end

function file = figSearch(D, A, opt)
a = D.acq.coarse;
fig = newFig(opt, [120 100 1000 380]);
ax = axes(fig);
plot(ax, a.f/1e6, 10*log10(mean(10.^(levels(a, D)/10), 1)), 'Color', [0 0.3 0.8]);
hold(ax, 'on'); grid(ax, 'on');
if ~isempty(A.search) && A.search.found
    xline(ax, A.search.fc/1e6, 'r--', sprintf('%s (%.1f dBm)', fmtHz6(A.search.fc), A.search.Pch));
end
xlabel(ax, 'Frequency (MHz)'); ylabel(ax, 'dBm in RBW'); title(ax, sprintf('Band search, RBW %s', fmtHz(a.settings.rbw)));
file = saveFig(fig, fullfile(opt.OutDir, "fig_search.png"), opt);
end

function file = figSpurs(D, A, opt)
fig = newFig(opt, [140 120 1100 760]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
a = getAcq(D, "wide");
ax = nexttile(tl);
if ~isempty(a)
    pk = max(levels(a, D), [], 1);
    plot(ax, a.f/1e9, pk, 'Color', [0 0.3 0.8], 'DisplayName', 'max hold, peak detector'); hold(ax, 'on');
    yline(ax, A.power.C_dBm - 60, 'r--', 'DisplayName', 'carrier - 60 dBc');
    T = A.wide.spurs;
    if height(T) > 0
        plot(ax, T.Freq_Hz/1e9, T.Level_dBm, 'rv', 'MarkerFaceColor', 'r', 'DisplayName', 'peaks');
        for i = 1:min(8, height(T))
            text(ax, T.Freq_Hz(i)/1e9, T.Level_dBm(i) + 3, sprintf('%.0f dBc', T.Rel_dBc(i)), 'FontSize', 7, 'HorizontalAlignment', 'center');
        end
    end
    legend(ax, 'Location', 'northeast', 'FontSize', 8);
    title(ax, "Spurious scan: " + A.wide.status + " (" + A.wide.note + ")", 'FontSize', 9);
end
grid(ax, 'on'); xlabel(ax, 'Frequency (GHz)'); ylabel(ax, 'dBm in RBW');
ax = nexttile(tl);
a = getAcq(D, "harm");
if ~isempty(a)
    plot(ax, (a.f - 2*A.spec.fc)/1e3, 10*log10(mean(10.^(levels(a, D)/10), 1)), 'DisplayName', ...
        sprintf('att %g dB', a.settings.atten)); hold(ax, 'on');
    a2 = getAcq(D, "harm2");
    if ~isempty(a2)
        plot(ax, (a2.f - 2*A.spec.fc)/1e3, 10*log10(mean(10.^(levels(a2, D)/10), 1)), 'DisplayName', ...
            sprintf('att %g dB', a2.settings.atten));
    end
    xline(ax, [-1 1]*(1 + A.spec.alpha)*A.spec.Rs/1e3, 'k:', 'HandleVisibility', 'off');
    legend(ax, 'FontSize', 8);
    m = sprintf('< %.1f dBc', A.harm.limit_dBc); if A.harm.detected, m = sprintf('%.1f dBc', A.harm.dBc); end
    title(ax, sprintf('2nd harmonic window (2 fc = %s): %s, %s', fmtHz6(2*A.spec.fc), m, A.harm.status), 'FontSize', 9);
end
grid(ax, 'on'); xlabel(ax, 'Offset from 2 fc (kHz)'); ylabel(ax, 'dBm in RBW');
file = saveFig(fig, fullfile(opt.OutDir, "fig_spurs.png"), opt);
end

function file = figEnvelope(D, A, opt)
fig = newFig(opt, [160 140 1200 420]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
V = A.env; co = lines(3); mods = ["BPSK" "QPSK" "8PSK"];
ax1 = nexttile(tl); ax2 = nexttile(tl);
if V.done
    e = -30:0.5:12;
    histogram(ax1, 10*log10(V.samples), e, 'Normalization', 'pdf', 'FaceColor', [0.6 0.6 0.6], 'DisplayName', 'measured');
    hold(ax1, 'on'); hold(ax2, 'on');
    [c, x] = ccdfCurve(V.samples); semilogy(ax2, x, c, 'k', 'LineWidth', 1.5, 'DisplayName', 'measured');
    for i = 1:3
        h = histcounts(10*log10(V.simSamples{i}), e, 'Normalization', 'pdf');
        plot(ax1, e(1:end - 1) + 0.25, h, 'Color', co(i, :), 'LineWidth', 1.2, 'DisplayName', ...
            sprintf('%s (p = %.2f)', mods(i), V.post(i)));
        [c, x] = ccdfCurve(V.simSamples{i}); semilogy(ax2, x, c, 'Color', co(i, :), 'DisplayName', mods(i));
    end
    set(ax2, 'YScale', 'log'); ylim(ax2, [1e-4 1]); xlim(ax2, [0 10]);
    legend(ax1, 'Location', 'northwest', 'FontSize', 7); legend(ax2, 'FontSize', 7);
    title(ax1, "Envelope: " + V.decision, 'FontSize', 8);
end
xlabel(ax1, 'Instantaneous power / mean (dB)'); ylabel(ax1, 'pdf'); grid(ax1, 'on');
xlabel(ax2, 'dB above mean'); ylabel(ax2, 'CCDF'); title(ax2, 'Peak-to-average (CCDF)', 'FontSize', 8); grid(ax2, 'on');
ax = nexttile(tl);
a = getAcq(D, "burst");
if ~isempty(a)
    P = levels(a, D); t = (0:size(P, 2) - 1)*a.settings.sweepTime/(size(P, 2) - 1);
    plot(ax, 1e3*t, P(1, :), 'Color', [0 0.3 0.8]); hold(ax, 'on');
    yline(ax, A.burst.Pn, ':', 'noise'); yline(ax, A.burst.Pon, '--', 'expected on');
    title(ax, "Zero span: " + A.burst.state, 'FontSize', 8);
end
xlabel(ax, 'Time (ms)'); ylabel(ax, 'dBm (RBW)'); grid(ax, 'on');
file = saveFig(fig, fullfile(opt.OutDir, "fig_envelope.png"), opt);
end

function [c, x] = ccdfCurve(p)
x = 0:0.1:12; d = 10*log10(p(:));
c = arrayfun(@(v) mean(d > v), x);
end

function file = figLines(D, A, opt)
a = getAcq(D, "line"); L = A.lines;
fig = newFig(opt, [180 160 1000 400]);
ax = axes(fig);
plot(ax, (a.f - A.spec.fc)/1e3, 10*log10(mean(10.^(levels(a, D)/10), 1))); hold(ax, 'on'); grid(ax, 'on');
pl = 10*log10(mean(10.^(levels(a, D)/10), 1));
if ~isempty(L.loc), plot(ax, (a.f(L.loc) - A.spec.fc)/1e3, pl(L.loc) + 1, 'rv', 'MarkerSize', 4); end
xlabel(ax, 'Offset from fc (kHz)'); ylabel(ax, 'dBm in RBW');
title(ax, sprintf('Idle-frame line scan, RBW %s: %s', fmtHz(a.settings.rbw), L.note), 'FontSize', 8);
file = saveFig(fig, fullfile(opt.OutDir, "fig_lines.png"), opt);
end

function fig = newFig(opt, pos)
fig = figure('Visible', opt.Visible, 'Color', 'w', 'Position', pos);
end

function file = saveFig(fig, file, opt)
exportgraphics(fig, file, 'Resolution', 130);
if strcmpi(opt.Visible, 'off'), close(fig); end
end

%% ===== report ====================================================================

function file = makeReport(R, opt, figs)
import mlreportgen.report.*
import mlreportgen.dom.*
D = R.data; A = R.results;
fprintf('Writing the %s report ...\n', upper(opt.Report));
rpt = Report(char(fullfile(opt.OutDir, "srs4_report")), char(opt.Report));
when = datetime('now', 'TimeZone', 'Asia/Bangkok', 'Format', 'yyyy-MM-dd HH:mm ''(Asia/Bangkok)''');
sub = "FieldFox N9914A measurement - " + D.idn;
if D.source == "demo", sub = "SYNTHETIC DEMO DATA - not a measurement"; end
add(rpt, TitlePage('Title', 'SRS-4 S-band downlink analysis', 'Subtitle', char(sub), ...
    'Author', char(string(getenv('USERNAME'))), 'PubDate', char(when)));
add(rpt, TableOfContents);

ch = Chapter('Title', 'Summary');
if D.source == "demo"
    p = Paragraph('SYNTHETIC DEMO: ccsdsTMWaveformGenerator waveform through a FieldFox measurement model.'); p.Bold = true; p.Color = 'red'; add(ch, p);
end
add(ch, Paragraph('Automatic verdicts against the Satlab SRS-4 datasheet (SLDS-SRS4-1.2) and ECSS/SFCG limits. Check them against the figures.'));
add(ch, rptTable(R.limits(:, {'Parameter', 'Measured', 'Limit', 'Status'}), 'Results vs limits', 40, {'23%', '39%', '19%', '19%'}));
add(ch, UnorderedList(cellstr(R.conclusions)));
if ~A.ok                                  % no downlink found: band search and notes only
    if ~isempty(D.notes), add(ch, UnorderedList(cellstr(D.notes))); end
    addImg(ch, "fig_search.png", 'Band search');
    add(rpt, ch);
    close(rpt);
    file = string(rpt.OutputPath);
    fprintf('Report: %s\n', file);
    return
end
addImg(ch, "fig_spectrum.png", 'Spectrum with RC fit and CCSDS reference (top); SFCG/ECSS mask (bottom)');
add(rpt, ch);

ch = Chapter('Setup and method');
add(ch, Paragraph(sprintf('Instrument: %s. Source: %s. ExtGain %g dB (%s). Frequency reference: %s.', D.idn, D.source, ...
    D.extGain, ternary(D.extgInstr == D.extGain, "set on the instrument", "applied in MATLAB"), D.refSource)));
if strlength(opt.Notes) > 0, add(ch, Paragraph("Notes: " + opt.Notes)); end
add(ch, rptTable(acqTable(D), 'Acquisitions (as configured on the analyzer)', 20));
if ~isempty(D.notes), add(ch, UnorderedList(cellstr(D.notes))); end
add(ch, UnorderedList({ ...
    'Safe start: preamp OFF, attenuation AUTO; ADC overload (CALC:MEAS:WAOR?) checked on the first sweep of each set, attenuation +10 dB per overload (max 30 dB).', ...
    'Average detector; sweeps power-averaged in MATLAB (mW, not dB). Levels are DUT-referred: display + instrument external gain - ExtGain.', ...
    'RBW filter modelled as Gaussian, noise bandwidth 1.0645 x RBW (assumption; cross-checked by the FieldFox channel power).'}));
add(rpt, ch);

ch = Chapter('Spectrum: symbol rate, roll-off, bandwidths');
S = A.spec;
add(ch, Paragraph(['Theory: the SRS-4 transmit filter is root-raised-cosine (RRC), so the transmitted power spectrum is ' ...
    '|H_RRC(f)|^2 = raised cosine: flat to (1-alpha)Rs/2, a cosine roll-off to (1+alpha)Rs/2. It is exactly 3 dB down at ' ...
    'fc +- Rs/2 for any alpha, so the -3 dB bandwidth equals the symbol rate. The occupied bandwidth (99 %) is 1.073 Rs at alpha 0.2.']));
T = table(["-3 dB edges"; "OBW midpoint"; "RC fit (nlparci, RBW-scaled)"; "RC fit (bootstrap)"], ...
    [S.fc3; S.fcObw; S.fc; mean(S.fit.bciFc)], [string(fmtHz(S.Rs3)); "-"; ci2(S.fit.ciRs, @fmtHz); ci2(S.fit.bciRs, @fmtHz)], ...
    ["-"; sprintf("%.3f (OBW lookup)", S.alphaObw); ci2(S.fit.ciAlpha, @(v) sprintf('%.3f', v)); ci2(S.fit.bciAlpha, @(v) sprintf('%.3f', v))], ...
    'VariableNames', {'Method', 'fc_Hz', 'Rs', 'Alpha'});
add(ch, rptTable(T, 'Estimates (95 % intervals)', 10));
add(ch, Paragraph(sprintf(['Fit: lsqcurvefit on %d points within %g dB of the top, noise %s; residual rms %.2f dB. Reference: %s, ' ...
    'residual rms %s dB for alpha %s. OBW 99 %%: %s, 99.9 %%: %s, -26 dB: %s.'], S.fit.nUsed, opt.FitRange, ...
    ternary(S.fit.fixedN0, "fixed from the floor", "fitted"), S.fit.rms, A.ref.note, mat2str(round(A.ref.rms, 2)), ...
    mat2str(A.ref.alpha), fmtHz(S.obw99), fmtHz(S.obw999), fmtHz(S.ebw26))));
addImg(ch, "fig_fit.png", 'Fit residuals and bootstrap distributions');
add(rpt, ch);

ch = Chapter('Frequency error');
Q = A.freq;
add(ch, Paragraph(sprintf(['fc = %s, 95 %% interval %s ... %s (%s). Reference: %s (+-%.2g ppm = +-%.0f Hz at fc). ' ...
    'The PSK carrier is suppressed: there is no line to count, fc comes from the spectral fit.'], fmtHz6(Q.fc), ...
    fmtHz6(Q.ci(1)), fmtHz6(Q.ci(2)), Q.ciSource, Q.refNote, Q.refPpm, Q.refPpm*Q.fc*1e-6)));
if ~isempty(Q.expected)
    add(ch, Paragraph(sprintf('Expected %s: error %+.0f Hz = %+.3f ppm, uncertainty +-%.2f ppm -> %s (0.5 ppm), %s (2.5 ppm).', ...
        fmtHz6(Q.expected), Q.err, Q.ppm, Q.Uppm, Q.status05, Q.status25)));
end
add(ch, Paragraph(['For ppm-level claims lock the FieldFox to a GPS-disciplined 10 MHz reference (Ref In, ROSC:SOUR EXT). ' ...
    'From orbit the Doppler shift is up to +-53 kHz (24 ppm) at 2245 MHz for a 500 km LEO, far more than the 0.5 ppm spec.']));
add(rpt, ch);

ch = Chapter('Power, C/N0 and Eb/N0');
W = A.power;
add(ch, Paragraph(sprintf(['Channel power %.2f dBm in %s (bandpower), signal C = %.2f dBm after removing the noise ' ...
    '(fit: %.2f dBm). Analyzer noise density N0 = %.1f dBm/Hz (%s) -> effective noise figure about %.0f dB. ' ...
    'C/N0 = %.1f dB-Hz, Es/N0 = C/N0 - 10 log10(Rs) = %.1f dB.'], W.ChP_dBm, fmtHz(W.ibw), W.C_dBm, W.Cfit_dBm, ...
    W.N0_dBmHz, A.floor.source, W.nfAnalyzer, W.CN0, W.EsN0)));
add(ch, rptTable(W.ebno, 'Eb/N0 at the analyzer: Eb/N0 = Es/N0 - 10 log10(bits/symbol x code rate)', 5));
add(ch, Paragraph(['Exam context - SRS-4 receiver Table 7 (sensitivity for < 1 % PER, NF 2.05 dB): required C/N0 = P + 174 - NF; ' ...
    'required Eb/N0 = C/N0 - 10 log10(information rate).']));
add(ch, rptTable(W.table7, 'SRS-4 Table 7 at the nearest rate', 5));
add(ch, Paragraph("Transmitter: " + W.txSource + ternary(isfinite(W.Ptx), sprintf(" -> %.1f dBm", W.Ptx), "") + ". " + W.gainNote));
if A.inst.done
    add(ch, Paragraph(sprintf('FieldFox CHP %.2f dBm vs MATLAB %.2f dBm over the same %s (%+.2f dB). FieldFox OBW 99 %%: %s, 99.9 %%: %s.', ...
        A.inst.chp, A.inst.chpML, fmtHz(D.inst.ibw), A.inst.diff, fmtHz(A.inst.obw99), fmtHz(A.inst.obw999))));
end
add(rpt, ch);

ch = Chapter('Spectral mask and bandwidth limits');
add(ch, Paragraph(sprintf(['ECSS-E-ST-50-05C Rev.2 Table 5-5 (category A, 2200-2290 MHz): %s -> limit %s; measured OBW99 %s: %s. ' ...
    'SFCG/ECSS mask (Fig. 6-5, %s): %s; minimum margin %.1f dB at %.2f Rs; the noise floor is %.0f dB below the peak.'], ...
    A.ecss.rule, fmtHz(A.ecss.limit), fmtHz(A.ecss.obw), A.ecss.status, A.mask.name, A.mask.status, A.mask.minMargin, ...
    A.mask.worstX, -A.mask.noiseRel)));
add(ch, Paragraph('Mask breakpoints read from the ECSS figure (0 dB to 0.5 Rs, -30 dB at 1.4 or 3 Rs, -60 dB at 3 or 8 Rs): check against SFCG Rec. 21-2R5.'));
add(rpt, ch);

ch = Chapter('Spurious emissions and 2nd harmonic');
if A.harm.done
    add(ch, Paragraph(sprintf('2nd harmonic window: %s. %s. %s', ternary(A.harm.detected, sprintf("%.1f dBc", A.harm.dBc), ...
        sprintf("not detected (< %.1f dBc)", A.harm.limit_dBc)), A.harm.note, A.harm.attTest)));
    if A.harm.line, add(ch, Paragraph('A spectral line at exactly 2 fc: squaring removes BPSK modulation (+-1 -> 1), so this points to BPSK.')); end
end
if A.wide.done
    add(ch, Paragraph(A.wide.note + ". Ambient labels are band guesses: repeat with the SRS-4 off to be sure."));
    if height(A.wide.spurs) > 0, add(ch, rptTable(A.wide.spurs, 'Peaks outside the out-of-band domain (2.5 x OBW)', 15)); end
end
addImg(ch, "fig_spurs.png", 'Spurious scan and 2nd harmonic window');
add(ch, Paragraph(['Limits: the 3rd harmonic (6.6-6.87 GHz) is beyond 6.5 GHz. The FieldFox makes its own 2nd harmonic ' ...
    '(-60 dBc at -30 dBm mixer level): use a high-pass filter (> 3 GHz) for a conducted harmonic test.']));
add(rpt, ch);

ch = Chapter('Zero span: envelope and bursts');
if A.env.done
    add(ch, Paragraph(char("Without I/Q the constellation is invisible, but the envelope is not: RRC-BPSK passes through zero " + ...
        "often (a symbol change +1 -> -1), QPSK/8PSK much less. Hypotheses are simulated with ccsdsTMWaveformGenerator " + ...
        "at the measured SNR, RBW and sample spacing; posterior with a model-error floor. " + A.env.note)));
    add(ch, rptTable(A.env.table, 'Envelope statistics', 5));
    add(ch, Paragraph("Decision: " + A.env.decision));
end
if A.burst.done, add(ch, Paragraph("Bursts: " + A.burst.state + ". " + A.burst.note)); end
add(ch, rptTable(A.frames, char("SRS-4 frame time by modulation and coding (223-byte block + ASM), Rs from " + A.frameRsSource), 12));
addImg(ch, "fig_envelope.png", 'Envelope statistics and zero-span timeline');
add(rpt, ch);

if A.lines.done
    ch = Chapter('Idle-frame spectral lines');
    add(ch, Paragraph(char("Repeated identical idle frames make the signal periodic with T_frame, so its spectrum is a comb " + ...
        "with spacing 1/T_frame; Rs = bits_per_frame / (bits_per_symbol x T_frame). " + A.lines.note)));
    addImg(ch, "fig_lines.png", 'Line scan');
    add(rpt, ch);
end

ch = Chapter('Exam notes and limitations');
add(ch, UnorderedList({ ...
    'RRC: the PSD is a raised cosine; -3 dB width = Rs for any alpha; null-to-null = (1+alpha) Rs.', ...
    'OBW99/Rs = 1.073 (0.20), 1.103 (0.25), 1.167 (0.35) for the ideal RC; PA regrowth widens OBW99.9 first.', ...
    'Frame (default 217-byte payload, RS+CC): 518 bytes = 4144 coded bits; T_frame = 4144 / (Rs log2 M) = 2.072 ms for QPSK at 1 MBd.', ...
    'Payload rate with RS+CC = 0.419 x Rs log2 M (0.84 Mbit/s for QPSK at 1 MBd).', ...
    'C/N0 = C - N0; Es/N0 = C/N0 - 10 log10 Rs; Eb/N0 = Es/N0 - 10 log10(log2 M x code rate).', ...
    'Friis: P_rx = P_tx + G_tx + G_rx - FSPL - L_pol; FSPL = 20 log10(4 pi d / lambda).', ...
    'No I/Q capture: EVM, constellation and modulation order are not measurable; the envelope test is indirect.', ...
    'FieldFox reference +-0.7 ppm: lock to a GPSDO for frequency-error claims at the 0.5 ppm level.', ...
    'Over the air the noise floor limits mask/spurious verification; compliance needs a conducted test with a >= 10 W pad.'}));
add(ch, UnorderedList(cellstr(R.files(1:min(end, 30)))));
add(rpt, ch);
close(rpt);
file = string(rpt.OutputPath);
fprintf('Report: %s\n', file);

    function addImg(parent, name, cap, w)
        if nargin < 4, w = '16cm'; end
        fpath = fullfile(opt.OutDir, name);
        if isfile(fpath) && any(endsWith(figs, name))
            add(parent, FormalImage('Image', char(fpath), 'Caption', cap, 'Width', w));
        end
    end
end

function T = acqTable(D)
rows = {};
for nm = string(fieldnames(D.acq)).'
    s = D.acq.(nm).settings;
    span = string(fmtHz(s.span)); if s.zeroSpan, span = sprintf("0 (%.3g s)", s.sweepTime); end
    rows(end + 1, :) = {nm, fmtHz6(s.center), span, fmtHz(s.rbw), s.points, s.detector, s.atten, size(D.acq.(nm).P, 1), s.overload}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'Step', 'Center', 'Span', 'RBW', 'Points', 'Detector', 'Atten_dB', 'Sweeps', 'Overload'});
end

function bt = rptTable(T, ttl, maxRows, widths)
% Report table from plain text cells (a MATLABTable would show strings in quotes).
import mlreportgen.report.*
import mlreportgen.dom.*
T = T(1:min(end, maxRows), :);
names = T.Properties.VariableNames;
body = cell(height(T), numel(names));
for j = 1:numel(names)
    x = T.(names{j});
    if isdatetime(x), c = string(x, 'yyyy-MM-dd HH:mm:ss');
    elseif endsWith(names{j}, "_Hz") && isnumeric(x), c = string(arrayfun(@fmtHz6, x, 'UniformOutput', false));
    elseif islogical(x), c = string(x);
    elseif isnumeric(x)
        c = compose("%.4g", x);
        k = isfinite(x) & x == round(x) & abs(x) < 1e7; c(k) = compose("%d", x(k));
    else, c = string(x);
    end
    c(ismissing(c)) = "-";
    body(:, j) = cellstr(c(:));
end
ft = FormalTable(replace(names, "_", " "), body);
ft.Style = {FontSize('7pt'), Border('solid'), ColSep('solid'), RowSep('solid'), Width('100%')};
ft.Header.Style = {Bold(true), BackgroundColor('#e8eef7')};
if nargin > 3 && numel(widths) == numel(names)
    grp = TableColSpecGroup;
    specs = TableColSpec.empty(0, 1);
    for k = 1:numel(widths)
        sp = TableColSpec; sp.Style = {Width(widths{k})}; specs(k) = sp;
    end
    grp.ColSpecs = specs;
    ft.ColSpecGroups = grp;
end
bt = BaseTable(ft);
bt.Title = ttl;
end

function s = ci2(ci, fmt)
if any(~isfinite(ci)), s = "-"; return, end
s = string(fmt(ci(1))) + " ... " + string(fmt(ci(2)));
end

%% ===== small helpers =========================================================

function v = snapRBW(x)
% FieldFox RBW in 1-3-10 steps (10 Hz - 5 MHz), the largest step not above x.
st = unique([reshape([1; 3]*10.^(1:6), 1, []), 5e6]);
st = st(st <= 5e6);
v = st(find(st <= x*1.0001, 1, 'last'));
if isempty(v), v = st(1); end
end

function n = autoPoints(span, rbw, nmin, nmax)
n = min(nmax, max(nmin, ceil(2*span/rbw) + 1));
if mod(n, 2) == 0, n = n + 1; end
n = min(n, nmax);
end

function [s, e] = runs(mask, gap)
% Start/end indices of true-runs; gaps up to 'gap' samples are bridged.
mask = mask(:).';
if gap > 0
    d = diff([1 mask 1]);
    gs = find(d == -1); ge = find(d == 1) - 1;
    for i = 1:numel(gs)
        if ge(i) - gs(i) + 1 <= gap && gs(i) > 1 && ge(i) < numel(mask), mask(gs(i):ge(i)) = true; end
    end
end
d = diff([0 mask 0]);
s = find(d == 1); e = find(d == -1) - 1;
end

function v = getf(s, name, default)
v = default;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); end
end

function s = mergeStruct(s, r)
for k = string(fieldnames(r)).', s.(k) = r.(k); end
end

function s = ternary(c, a, b)
if c, s = string(a); else, s = string(b); end
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
s = sprintf('%.5g %s', x/scale, unit);
end

function s = fmtHz6(x)
if ~isfinite(x), s = '-'; return, end
s = sprintf('%.6f MHz', x/1e6);
end

function s = fmtNum(x)
s = sprintf('%.6g', x);
end

function s = fmtRange(b)
if diff(b) < 1e-3*abs(mean(b)), s = sprintf('%s +- %s', fmtHz6(mean(b)), fmtHz(diff(b)/2)); return, end
s = sprintf('%s - %s', fmtHz(b(1)), fmtHz(b(2)));
end

function s = stamp()
s = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
end
