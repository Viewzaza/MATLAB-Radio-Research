function R = ff_ccsds_linksim(varargin)
%FF_CCSDS_LINKSIM  Monte-Carlo link simulator of the Satlab SRS-4 CCSDS telemetry chain (BER / PER vs Eb/N0).
%
%   R = ff_ccsds_linksim                       QPSK + RS + CC (SRS-4 default), ideal receiver, 'Quick' (~1 min)
%   R = ff_ccsds_linksim(Name, Value, ...)
%   R = ff_ccsds_linksim('Receiver', "realistic", 'Impairments', "srs4")  AFC + timing + carrier loops + ASM sync
%   R = ff_ccsds_linksim('Preset', "full")     exam-quality curves (8 samples/symbol, 0.25 dB steps, more frames)
%   R = ff_ccsds_linksim('Pass', "...\cubesat_link.mat")  PER and throughput over the ff_cubesat_link passes
%   R = ff_ccsds_linksim('Measured', "...\srs4_run.mat") BER/PER at the C/N0 measured by ff_srs4_analyzer
%   R = ff_ccsds_linksim('Pass', "demo", 'Measured', "demo")   both modes with synthetic inputs
%
% Example (headless; the realistic receiver against an uncompensated LEO Doppler):
%   R = ff_ccsds_linksim('Receiver', "realistic", 'Impairments', "leo", 'AFCRange', 100e3, ...
%         'Visible', 'off', 'Prompt', false);
%   R.summary, R.curves, R.validation.fig13, R.validation.table7
%
% No instrument and no network are used: this is the computational core behind the FieldFox
% measurements (ff_srs4_analyzer) and the pass planner (ff_cubesat_link).
%
% Transmitter = the SRS-4 downlink frame, built with ccsdsTMWaveformGenerator
% (Satellite Communications Toolbox), CCSDS 131.0-B-3:
%   223-byte transfer frame (SRS-4: 2-byte header + 217-byte payload + CRC32C)
%   -> RS(255,223), interleaving depth 1, dual basis (= Satlab)       ['RS', 'concatenated']
%   -> CCSDS randomizer -> ASM 1ACFFC1D
%   -> convolutional code r = 1/2, K = 7 (171,133 octal), G2 inverted ['CC', 'concatenated']
%   -> BPSK / QPSK / 8PSK -> root-raised-cosine (RRC) filter, roll-off 'Rolloff', span 10 symbols.
%   QPSK + RS + CC: 4144 coded bits = 2072 symbols = 2.072 ms per frame at 1 MBd (verified).
%   8PSK: the generator maps the continuous code-bit stream (frames are not symbol aligned:
%   4144/3 is not an integer - an ASSUMPTION for the SRS-4). Its 8PSK path is ~25x slower, so
%   the same code bits + pskmod + comm.RaisedCosineTransmitFilter are used; identity with the
%   generator is checked on every run (R.config.tx8check ~1e-16).
%   The convolutional encoder runs continuously across frames, so the ASM is encoded too:
%   code bits 13-64 of every frame are a fixed 52-bit pattern (the first 12 depend on the
%   previous frame) - that pattern is what the frame synchroniser searches for.
% Channel: AWGN (awgn) at the requested Eb/N0, plus optional impairments (below).
% Receiver ('Receiver'):
%   "ideal"      genie: perfect carrier phase, frequency, symbol timing and frame alignment.
%                Matched RRC filter (comm.RaisedCosineReceiveFilter), symbol-spaced samples.
%                Doppler / DopplerRate / TimingOffset are NOT applied (a genie removes them);
%                phase noise and IQ imbalance ARE applied (random, cannot be removed).
%   "realistic"  random carrier phase and symbol timing in every batch, plus the impairments:
%                AFC = comm.CoarseFrequencyCompensator (FFT of the M-th power of the matched-
%                filtered acquisition frames; rejected when the two halves disagree - noise;
%                clamped to +-'AFCRange') -> matched RRC filter (2 sps)
%                -> comm.SymbolSynchronizer (Gardner) -> comm.CarrierSynchronizer (2nd-order
%                PLL) -> ASM frame sync on the bit stream (all M phase rotations; resolves the
%                PSK phase ambiguity; acquisition over 4 idle frames, then flywheel +- 2
%                symbols) -> as below. Frames the synchroniser loses count as packet errors.
%   Both: soft demodulation (pskdemod, exact LLR; max-log for 8PSK) -> Viterbi (vitdec, soft
%   'unquant', traceback 42, per frame with 48 bits of context each side; the G2 inversion
%   is undone by negating every second LLR) -> derandomizer -> ccsdsRSDecode.
% Measured: BER of the 1784 transfer-frame bits, packet error rate (PER = frame error rate:
% one SRS-4 packet per frame), raw demodulator BER (hard decisions on the code bits), RS
% decoder failures, frames lost / symbol slips / rejected AFC estimates; 95 % confidence
% intervals (berconfint).
% Theory: berawgn (uncoded PSK), bercoding (CC soft union bound with distspec; RS hard),
% binocdf (RS t = 16 packet error), Shannon limit for the same spectral efficiency.
%
% Eb/N0 definitions (the axis 'EbNo' and all tables use the FIRST one):
%   Eb/N0 (info)    per transfer-frame bit (1784 per frame), ASM + coding overhead included:
%                   Es/N0 = Eb/N0 + 10 log10(1784 x bits/symbol / code bits per frame)
%   Eb/N0 (payload) per SRS-4 payload bit (217 bytes), = info + 10 log10(1784/1736) = +0.12 dB.
%                   ff_cubesat_link and the Table 7 conversion use this one.
%   Eb/N0 (channel) per modulated code bit = Es/N0 - 10 log10(bits/symbol): the x axis of the
%                   SRS-4 datasheet Figure 13 (uncoded BER of the demodulator).
% In AWGN with the ideal receiver, curves in Eb/N0 do not depend on the symbol rate.
% 'SymbolRate' only matters for the realistic receiver (Doppler, phase noise in Hz) and for
% the conversions to dBm (Table 7), bit rate and throughput.
%
% Validation against the SRS-4 datasheet (rev 1.2):
%   Figure 13: demodulator implementation loss < 0.3 dB for Eb/N0 1-10 dB (uncoded BER,
%              on hardware). Reproduced with an uncoded BPSK/QPSK sweep through the
%              REALISTIC receiver: loss = Eb/N0 - berawgn^-1(raw BER), per point, with its
%              95 % interval. The ideal receiver must give ~0 dB (simulator calibration).
%   Table 7:   input power for PER < 1 % (217-byte payload). Converted with NF 'NoiseFigure':
%              Eb/N0 = P + 174 - NF - 10 log10(payload bit rate), compared with the simulated
%              Eb/N0 at PER = 1 %. The difference = what the hardware loses beyond this model
%              (implementation loss, NF tolerance, test set-up). 'Table7Sweep' fills all columns.
%   Example result (quick preset, ideal receiver, 'Table7Sweep', 10 Oct 2026): Figure 13 loss
%   0.08-0.18 dB (BPSK, realistic receiver: PASS); Table 7 minus model, mean per FEC: none
%   +0.05 dB, CC +0.21 dB, RS +2.64 dB, RS+CC +1.07 dB. Demodulator, Viterbi and NF 2.05 dB
%   explain the uncoded and CC columns; the RS columns are 1-2.6 dB worse than hard-decision
%   RS(255,223) with t = 16 predicts - an open question for the SRS-4 (not in this model).
%
% Extra modes:
%   'Pass'      ff_cubesat_link result (struct or cubesat_link.mat): its per-second Eb/N0 for
%               your configuration -> PER(t) from the simulated curve -> throughput(t), data
%               volume per pass; compared with the margin-threshold model of ff_cubesat_link.
%               The configuration (modulation, FEC, symbol rate) is taken from that result
%               unless you set 'Modulation'/'Coding'/'SymbolRate' yourself. The curves are
%               AWGN; add 'Receiver' "realistic", 'Impairments' "leo" to include Doppler.
%   'Measured'  ff_srs4_analyzer result (struct or srs4_run.mat): its C/N0 (at the FieldFox, whose
%               noise floor is ~20 dB above an SRS-4 receiver) and channel power -> Eb/N0 for
%               your configuration -> BER/PER, margin to PER 1 %, the extra path loss the link
%               could take and the equivalent distance. Also for an SRS-4-class receiver
%               (NF 2.05 dB) at the same input power.
%
% Output R: settings, config (frame structure, rates, Eb/N0 conversions), impairments,
% curves (table per Eb/N0 point: BER, PER, raw BER, CIs, frames, frames lost, slips, rejected
% AFC estimates, source sim/cache/skipped), theory,
% summary (Eb/N0 at PER 1 % and BER 1e-5 with CI, coding gain, Shannon gap, losses),
% ideal (comparison sweep, realistic receiver only), validation.fig13 / .table7,
% pass, measured, notes, files, outDir. Files in 'OutDir': linksim.mat, curves.csv,
% [table7.csv, fig13.csv, pass.csv, measured.csv], fig_*.png.
% Curves are cached in support\linksim\linksim_cache.mat (per configuration and Eb/N0 point):
% a second run re-uses them and only adds frames where the stop rule is not met yet.
%
% Name-Value options (defaults):
%  Link
%   'Modulation'     "QPSK"   "BPSK" | "QPSK" | "8PSK" (8PSK: TX only on the SRS-4, no Table 7)
%   'Coding'         "concatenated"  "none" | "CC" | "RS" | "concatenated" (= "CC+RS")
%   'Rolloff'        0.2      RRC roll-off (SRS-4: 0.2, 0.25, 0.35)
%   'SymbolRate'     1e6      Bd (sl-radio --tx-rate=1000000)
%   'EbNo'           []       dB (info), [] = preset grid for the coding (sorted ascending)
%  Monte Carlo
%   'Preset'         "quick"  "quick": 4 samples/symbol, 0.5 dB steps, <= 600 frames and 30
%                             packet errors per point, 60 s per sweep;
%                             "full": 8 samples/symbol, 0.25 dB, <= 20000 frames, 100 errors
%   'Frames'         []       max frames per Eb/N0 point (preset)
%   'MaxErrors'      []       stop a point after this many packet (frame) errors (preset)
%   'MinFrames'      []       at least this many frames per point (preset: 50 / 200)
%   'TimeLimit'      []       s for the main sweep (preset: 60 / Inf), shared out over the open
%                             points; the ideal-comparison and Figure 13 sweeps get half each
%   'SamplesPerSymbol' []     even, >= 2 (preset 4 / 8)
%   'UseParallel'    false    parfor over Eb/N0 points (identical results). Measured: 'full', 7
%                             points x 1500 frames: 82 s serial, 32 s on 6 workers (+30 s to
%                             start the pool); 'quick' is slower in parallel - keep false
%   'Seed'           1        random seed (results are reproducible)
%   'Cache'          true     true | false | "refresh" (ignore and overwrite cached points)
%   'CacheFile'      ""       default support\linksim\linksim_cache.mat next to this file
%  Receiver and impairments
%   'Receiver'       "ideal"  "ideal" | "realistic"
%   'CompareIdeal'   true     realistic: also run the ideal receiver -> synchronisation loss
%   'Impairments'    "none"   preset: "none"; "srs4" = 0.5 ppm frequency error (1.12 kHz at
%                             2245 MHz), 0.5 ppm symbol clock, ASSUMED phase noise (below) and
%                             IQ imbalance 0.2 dB / 1 deg; "leo" = "srs4" + 45 kHz Doppler and
%                             -800 Hz/s (S band, 500 km, uncompensated). The options below override.
%   'Doppler'        []       Hz, carrier frequency offset at the start of every batch
%   'DopplerRate'    []       Hz/s, frequency ramp within a batch
%   'PhaseNoise'     []       [] / "none" | "typical" = [-80 -88 -95 -115] dBc/Hz at [1k 10k 100k 1M] Hz
%                             (an ASSUMPTION: the SRS-4 datasheet gives no phase noise) | a 2 x N
%                             matrix [dBc/Hz; offset Hz] for comm.PhaseNoise (offsets < fs/2)
%   'TimingOffset'   []       [tau ppm]: fixed delay tau (symbols, NaN = random every batch) and
%                             symbol-clock offset in ppm; [] = [NaN 0] (random timing phase)
%   'IQImbalance'    []       [amplitude_dB phase_deg] (iqimbal), transmitter side, not corrected
%   'AFCRange'       15e3     Hz, AFC clamp (SRS-4 datasheet 15 kHz; firmware 2.0: 100 kHz)
%   'LoopBandwidth'  [0.005 0.005]  normalised (x Rs) bandwidth of the timing and carrier loops
%  Validation
%   'Validate'       true     Figure 13 sweep (uncoded, realistic receiver, BPSK/QPSK) + Table 7
%   'Table7Sweep'    false    true = simulate all 4 FEC columns for BPSK and QPSK (8 sweeps)
%   'NoiseFigure'    2.05     dB, SRS-4 receiver NF used for Table 7 and 'Measured'
%  Modes
%   'Pass'           []       ff_cubesat_link struct | MAT file | "demo" | true (choose a file)
%   'Measured'       []       ff_srs4_analyzer struct | MAT file | "demo" | true (choose a file)
%   'ImplLoss'       0        dB subtracted from the Eb/N0 of 'Pass' and 'Measured' (hardware
%                             loss beyond the simulated receiver, e.g. the Table 7 difference)
%   'MarginTarget'   3        dB, ff_cubesat_link's threshold model (for the comparison)
%  Output
%   'PlotEbNo'       []       dB (info) for constellation / eye / spectrum; [] = PER nearest 10 %
%   'Scopes'         false    true = also comm.ConstellationDiagram + spectrumAnalyzer (desktop only)
%   'Plot'           true     figures (PNG files are written even with 'Visible' 'off')
%   'OutDir'         ""       default linksim_<yyyymmdd_HHMMSS> in the current folder
%   'Save'           true     write MAT/CSV/PNG files
%   'Visible'        'on'     'off' = no windows
%   'Prompt'         true     ask before a long run / for a file (never in -batch)
%   'Verbose'        true     progress in the Command Window
%
% Read this before trusting the numbers:
%   - Monte Carlo resolution: a point with n frames and 0 errors only says PER < ~3/n
%     (95 %); BER below ~10/(bits simulated) is not measured. CIs are given for every point.
%   - The realistic receiver is a generic coherent demodulator with SRS-4-like parameters,
%     not Satlab's firmware: loop bandwidths, acquisition and AGC (here: ideal gain) differ.
%   - The decoder works frame by frame (with context), so a cycle slip or symbol slip costs
%     the frame it happens in, like a real receiver with ASM re-acquisition.
%   - RS interleaving depth 1 (SRS-4 default frame): Viterbi error bursts are NOT spread over
%     several RS codewords - the concatenated "theory" (independent byte errors) is optimistic.
%   - 8PSK + CC has no closed-form theory here; Table 7 has no 8PSK (the SRS-4 RX is BPSK/QPSK).
%     The 8th-power AFC fails at the low Es/N0 of coded 8PSK (estimates rejected, see AFCRejects):
%     a real 8PSK ground receiver needs Doppler pre-compensation or data-aided acquisition.
%   - The FieldFox has no I/Q capture: the modulation order, EVM and the BER of a real SRS-4
%     cannot be measured with it. 'Measured' PREDICTS BER/PER from the measured C/N0 only.
%     For ppm-level frequency claims use an external 10 MHz GPSDO (FieldFox +-0.7 ppm vs the
%     SRS-4 0.5 ppm spec).
%
% See also ff_srs4_analyzer, ff_cubesat_link, ccsdsTMWaveformGenerator, ccsdsRSDecode, vitdec,
% berawgn, bercoding, berconfint, comm.SymbolSynchronizer, comm.CarrierSynchronizer.

[opt, usingDefaults] = parseOptions(varargin{:});
interactive = usejava('desktop') && ~batchStartupOptionUsed;
notes = strings(0, 1);
oldVis = get(groot, 'DefaultFigureVisible');
if strcmpi(opt.Visible, 'off'), set(groot, 'DefaultFigureVisible', 'off'); end
restoreVis = onCleanup(@() set(groot, 'DefaultFigureVisible', oldVis));

%% 0. Inputs of the extra modes (they may define the configuration)
passIn = []; measIn = [];
if ~isempty(opt.Pass) && ~isequal(opt.Pass, false)
    passIn = loadPass(opt.Pass, opt, interactive);
    if ~isempty(passIn.user)
        u = passIn.user;
        if ismember('Modulation', usingDefaults), opt.Modulation = u.Modulation; end
        if ismember('Coding', usingDefaults), opt.Coding = u.Coding; end
        if ismember('SymbolRate', usingDefaults), opt.SymbolRate = u.SymbolRate; end
        notes(end+1, 1) = sprintf("Pass mode: configuration %s, %s, %s from %s.", opt.Modulation, ...
            opt.Coding, rateText(opt.SymbolRate), passIn.source);
    end
end
if ~isempty(opt.Measured) && ~isequal(opt.Measured, false)
    measIn = loadMeasured(opt.Measured, opt, interactive);
end

%% 1. Configuration
run = presetRun(opt);
imp = impairmentSet(opt);
cfg = linkConfig(opt.Modulation, opt.Coding, opt.Rolloff, run.sps, opt.SymbolRate, opt.Receiver, imp);
grid = opt.EbNo;
if isempty(grid), grid = defaultGrid(cfg, run); end
grid = unique(double(grid(:)).');
cache = cacheOpen(opt);
if opt.Verbose
    fprintf('ff_ccsds_linksim: %s %s, roll-off %.2f, %s, %s receiver, preset %s (%d samples/symbol)\n', ...
        cfg.mod, cfg.codingName, cfg.rolloff, rateText(cfg.Rs), cfg.receiver, run.preset, run.sps);
    fprintf('  frame: %d info bits -> %d code bits = %.1f symbols = %.3f ms; payload %.1f kbit/s\n', ...
        cfg.Ni, cfg.Lc, cfg.symPerFrame, cfg.Tframe*1e3, cfg.RbPay/1e3);
end
if run.preset == "full" && opt.Prompt && interactive && isinf(run.TimeLimit)
    est = estimateMinutes(cfg, grid, run);
    if est > 5
        a = input(sprintf('Preset "full" may take ~%.0f min (cached points are re-used). Continue? [y/N] ', est), 's');
        if ~strncmpi(a, 'y', 1), error('ff_ccsds_linksim:cancelled', 'Cancelled by the user.'); end
    end
end

%% 2. Main sweep (+ ideal comparison)
[S, cache] = sweep(cfg, imp, grid, run, opt, cache, "main");
SI = [];
if cfg.realistic && opt.CompareIdeal
    cfgI = linkConfig(opt.Modulation, opt.Coding, opt.Rolloff, run.sps, opt.SymbolRate, "ideal", imp);
    runI = presetRun(opt, "ideal"); runI.TimeLimit = runI.TimeLimit/2;   % secondary sweep: half budget
    [SI, cache] = sweep(cfgI, imp, grid, runI, opt, cache, "ideal comparison");
end
th = theory(cfg, S.T.EbNo_dB);
summ = summarize(cfg, S, SI, opt.NoiseFigure);

%% 3. Validation against the datasheet
V = struct('fig13', table(), 'fig13Note', "", 'table7', table(), 'table7Note', "", 'calibration', table());
if opt.Validate
    [V.fig13, V.fig13Note, cache, V.fig13Sweep] = validateFig13(cfg, imp, run, opt, cache);
end
V.calibration = calibrationTable(S, SI, cfg);
T7 = table7Compare(cfg, summ.EbNo_PER1, summ.EbNo_PER1_CI, opt.NoiseFigure);
if opt.Table7Sweep
    [T7, cache] = table7Sweep(opt, imp, run, cache);
end
V.table7 = T7;
V.table7Note = table7Note(T7, opt.NoiseFigure);

%% 4. Pass and measurement predictions
PS = []; MS = [];
if ~isempty(passIn), PS = passPrediction(passIn, S, cfg, summ, opt); end
if ~isempty(measIn), MS = measuredPrediction(measIn, S, cfg, summ, opt); end
cacheSave(cache, opt);

%% 5. Output
R = struct();
R.settings = rmfield(opt, {'Pass', 'Measured'});
R.settings.PassSource = ""; R.settings.MeasuredSource = "";
if ~isempty(passIn), R.settings.PassSource = passIn.source; end
if ~isempty(measIn), R.settings.MeasuredSource = measIn.source; end
R.config = rmfield(cfg, {'pat', 'PN', 'trellis', 'txPSD'});
R.impairments = imp;
R.curves = S.T;
R.theory = th;
R.summary = summ;
R.ideal = [];
if ~isempty(SI), R.ideal = SI.T; end
R.validation = rmfieldIf(V, 'fig13Sweep');
R.pass = PS; R.measured = MS;
R.notes = [notes; S.notes; summ.notes; V.fig13Note; V.table7Note];
if ~isempty(PS), R.notes = [R.notes; PS.notes]; end
if ~isempty(MS), R.notes = [R.notes; MS.notes]; end
R.notes = R.notes(strlength(R.notes) > 0);
R.files = strings(0, 1); R.outDir = "";

if opt.Save
    if strlength(opt.OutDir) == 0
        opt.OutDir = fullfile(pwd, "linksim_" + string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
    end
    if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
    R.outDir = string(opt.OutDir);
end
if opt.Plot
    figs = makeFigures(cfg, imp, S, SI, summ, V, PS, MS, opt, interactive);
    R.files = [R.files; figs];
end
if opt.Save
    f = fullfile(opt.OutDir, "curves.csv"); writetable(S.T, f); R.files(end+1, 1) = f;
    if ~isempty(V.fig13), f = fullfile(opt.OutDir, "fig13.csv"); writetable(V.fig13, f); R.files(end+1, 1) = f; end
    if ~isempty(V.table7), f = fullfile(opt.OutDir, "table7.csv"); writetable(V.table7, f); R.files(end+1, 1) = f; end
    if ~isempty(PS), f = fullfile(opt.OutDir, "pass.csv"); writetable(PS.table, f); R.files(end+1, 1) = f; end
    if ~isempty(MS), f = fullfile(opt.OutDir, "measured.csv"); writetable(MS.table, f); R.files(end+1, 1) = f; end
    f = fullfile(opt.OutDir, "linksim.mat");
    R.files = [string(f); R.files(strlength(R.files) > 0)];
    save(f, 'R');
end
if opt.Verbose, printSummary(R, cfg); end
end

%% ===== options ===================================================================

function [opt, usingDefaults] = parseOptions(varargin)
p = inputParser; p.FunctionName = 'ff_ccsds_linksim';
add = @(n, v) addParameter(p, n, v);
add('Modulation', "QPSK"); add('Coding', "concatenated"); add('Rolloff', 0.2); add('SymbolRate', 1e6);
add('EbNo', []); add('Preset', "quick"); add('Frames', []); add('MaxErrors', []); add('MinFrames', []);
add('TimeLimit', []); add('SamplesPerSymbol', []); add('UseParallel', false); add('Seed', 1);
add('Cache', true); add('CacheFile', "");
add('Receiver', "ideal"); add('CompareIdeal', true); add('Impairments', "none");
add('Doppler', []); add('DopplerRate', []); add('PhaseNoise', []); add('TimingOffset', []);
add('IQImbalance', []); add('AFCRange', 15e3); add('LoopBandwidth', [0.005 0.005]);
add('Validate', true); add('Table7Sweep', false); add('NoiseFigure', 2.05);
add('Pass', []); add('Measured', []); add('ImplLoss', 0); add('MarginTarget', 3);
add('PlotEbNo', []); add('Scopes', false); add('Plot', true); add('OutDir', ""); add('Save', true);
add('Visible', 'on'); add('Prompt', true); add('Verbose', true);
parse(p, varargin{:});
opt = p.Results; usingDefaults = p.UsingDefaults;
opt.Modulation = upper(string(opt.Modulation));
if ~ismember(opt.Modulation, ["BPSK" "QPSK" "8PSK"])
    error('ff_ccsds_linksim:modulation', 'Modulation must be "BPSK", "QPSK" or "8PSK".');
end
opt.Coding = canonicalCoding(opt.Coding);
opt.Receiver = lower(string(opt.Receiver));
if ~ismember(opt.Receiver, ["ideal" "realistic"])
    error('ff_ccsds_linksim:receiver', 'Receiver must be "ideal" or "realistic".');
end
opt.Preset = lower(string(opt.Preset));
if ~ismember(opt.Preset, ["quick" "full"]), error('ff_ccsds_linksim:preset', 'Preset must be "quick" or "full".'); end
opt.Impairments = lower(string(opt.Impairments));
if ~ismember(opt.Impairments, ["none" "srs4" "leo"])
    error('ff_ccsds_linksim:impairments', 'Impairments must be "none", "srs4" or "leo".');
end
opt.OutDir = string(opt.OutDir); opt.CacheFile = string(opt.CacheFile);
opt.Visible = char(opt.Visible);
if ~ismember(lower(opt.Visible), {'on', 'off'}), error('ff_ccsds_linksim:visible', 'Visible must be ''on'' or ''off''.'); end
validateattributes(opt.Rolloff, {'numeric'}, {'scalar', '>', 0, '<=', 1}, 'ff_ccsds_linksim', 'Rolloff');
validateattributes(opt.SymbolRate, {'numeric'}, {'scalar', 'positive', 'finite'}, 'ff_ccsds_linksim', 'SymbolRate');
validateattributes(opt.LoopBandwidth, {'numeric'}, {'numel', 2, 'positive', '<', 0.2}, 'ff_ccsds_linksim', 'LoopBandwidth');
end

function c = canonicalCoding(c)
c = lower(strrep(string(c), ' ', ''));
switch c
    case {"none", "uncoded", "nofec"}, c = "none";
    case {"cc", "conv", "convolutional", "viterbi"}, c = "CC";
    case {"rs", "reedsolomon"}, c = "RS";
    case {"concatenated", "cc+rs", "rs+cc", "rscc", "ccrs", "all"}, c = "concatenated";
    otherwise
        error('ff_ccsds_linksim:coding', 'Coding must be "none", "CC", "RS" or "concatenated" (got "%s").', c);
end
end

function run = presetRun(opt, receiver)
if nargin < 2, receiver = opt.Receiver; end
if opt.Preset == "full"
    run = struct('preset', "full", 'sps', 8, 'Frames', 20000, 'MaxErrors', 100, 'MinFrames', 200, ...
        'TimeLimit', Inf, 'step', 0.25, 'fig13', 1:10, 'fig13Frames', 4000);
else
    run = struct('preset', "quick", 'sps', 4, 'Frames', 600, 'MaxErrors', 30, 'MinFrames', 50, ...
        'TimeLimit', 60, 'step', 0.5, 'fig13', [1 3 5 7 9], 'fig13Frames', 400);
end
if ~isempty(opt.SamplesPerSymbol), run.sps = opt.SamplesPerSymbol; end
if run.sps < 2 || mod(run.sps, 2) ~= 0
    error('ff_ccsds_linksim:sps', 'SamplesPerSymbol must be even and >= 2.');
end
if ~isempty(opt.Frames), run.Frames = opt.Frames; run.fig13Frames = min(run.fig13Frames, opt.Frames); end
if ~isempty(opt.MaxErrors), run.MaxErrors = opt.MaxErrors; end
if ~isempty(opt.MinFrames), run.MinFrames = opt.MinFrames; end
run.MinFrames = min(run.MinFrames, run.Frames);
if ~isempty(opt.TimeLimit), run.TimeLimit = opt.TimeLimit; end
run.receiver = string(receiver);
if run.receiver == "realistic"
    run.F = 30; run.nAcq = 4;            % data frames per batch, acquisition (idle) frames
else
    run.F = 25; run.nAcq = 0;
end
run.seed = opt.Seed; run.useParallel = logical(opt.UseParallel); run.verbose = logical(opt.Verbose);
end

function imp = impairmentSet(opt)
typicalPN = [-80 -88 -95 -115; 1e3 1e4 1e5 1e6];
imp = struct('preset', opt.Impairments, 'f0', 0, 'fdot', 0, 'pn', zeros(2, 0), 'tau', NaN, 'ppm', 0, ...
    'iq', zeros(1, 0), 'afc', opt.AFCRange, 'loop', double(opt.LoopBandwidth(:).'));
switch opt.Impairments
    case "srs4"
        imp.f0 = 0.5e-6*2245e6; imp.pn = typicalPN; imp.ppm = 0.5; imp.iq = [0.2 1];
    case "leo"
        imp.f0 = 45e3; imp.fdot = -800; imp.pn = typicalPN; imp.ppm = 0.5; imp.iq = [0.2 1];
end
if ~isempty(opt.Doppler), imp.f0 = opt.Doppler; end
if ~isempty(opt.DopplerRate), imp.fdot = opt.DopplerRate; end
if ~isempty(opt.PhaseNoise)
    pn = opt.PhaseNoise;
    if (isstring(pn) || ischar(pn)) && lower(string(pn)) == "typical", imp.pn = typicalPN;
    elseif (isstring(pn) || ischar(pn)) || isequal(pn, false), imp.pn = zeros(2, 0);
    else
        if size(pn, 1) ~= 2, error('ff_ccsds_linksim:pn', 'PhaseNoise must be [dBc/Hz; offset Hz] (2 x N).'); end
        imp.pn = double(pn);
    end
end
if ~isempty(opt.TimingOffset)
    to = double(opt.TimingOffset(:).');
    imp.tau = to(1); if numel(to) > 1, imp.ppm = to(2); end
end
if ~isempty(opt.IQImbalance)
    if isequal(opt.IQImbalance, false), imp.iq = zeros(1, 0); else, imp.iq = double(opt.IQImbalance(1:2)); end
end
if ~isnan(imp.tau), imp.tau = mod(imp.tau, 1); end
end

%% ===== link configuration ===========================================================

function cfg = linkConfig(modName, coding, rolloff, sps, Rs, receiver, imp)
cfg.mod = string(modName); cfg.coding = string(coding); cfg.rolloff = rolloff; cfg.sps = sps; cfg.Rs = Rs;
cfg.receiver = string(receiver); cfg.realistic = cfg.receiver == "realistic";
cfg.fs = Rs*sps; cfg.span = 10;
switch cfg.mod
    case "BPSK", cfg.M = 2; cfg.phOff = pi;   cfg.order = 'gray';            cfg.carOff = 0;    cfg.nRot = 2;
    case "QPSK", cfg.M = 4; cfg.phOff = pi/4; cfg.order = [0 2 3 1];         cfg.carOff = pi/4; cfg.nRot = 4;
    case "8PSK", cfg.M = 8; cfg.phOff = pi/4; cfg.order = [0 4 6 2 3 7 5 1]; cfg.carOff = 0;    cfg.nRot = 8;
end
cfg.b = log2(cfg.M);
cfg.isCC = any(cfg.coding == ["CC" "concatenated"]);
cfg.isRS = any(cfg.coding == ["RS" "concatenated"]);
names = struct('none', "none (uncoded)", 'CC', "CC r1/2 K7", 'RS', "RS(255,223)", 'concatenated', "RS(255,223) + CC r1/2 K7");
cfg.codingName = names.(char(cfg.coding));
cfg.genCoding = "none";
if cfg.isCC && cfg.isRS, cfg.genCoding = "concatenated";
elseif cfg.isCC, cfg.genCoding = "convolutional";
elseif cfg.isRS, cfg.genCoding = "RS";
end
cfg.Ni = 1784; cfg.Np = 217*8;                    % transfer-frame bits, SRS-4 payload bits
cfg.Lt = cfg.Ni; if cfg.isRS, cfg.Lt = 2040; end  % bits after the ASM (before CC)
cfg.Lf = 32 + cfg.Lt;
cfg.Lc = cfg.Lf*(1 + cfg.isCC);                   % code bits per frame
cfg.symPerFrame = cfg.Lc/cfg.b;
cfg.Tframe = cfg.symPerFrame/Rs;
cfg.RbInfo = cfg.Ni/cfg.Tframe; cfg.RbPay = cfg.Np/cfg.Tframe;
cfg.esFactor = 10*log10(cfg.Ni*cfg.b/cfg.Lc);     % Es/N0 = Eb/N0 + esFactor
cfg.chFactor = 10*log10(cfg.Ni/cfg.Lc);           % Eb/N0(channel) = Eb/N0 + chFactor
cfg.payFactor = 10*log10(cfg.Ni/cfg.Np);          % Eb/N0(payload) = Eb/N0 + payFactor
cfg.eta = cfg.Ni*cfg.b/cfg.Lc;                    % information bits per symbol
cfg.shannon = 10*log10((2^cfg.eta - 1)/cfg.eta);  % unconstrained complex AWGN limit
cfg.trellis = poly2trellis(7, [171 133]); cfg.tblen = 42; cfg.Lpre = 96; cfg.Lpost = 96;
asm = [0 0 0 1 1 0 1 0 1 1 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 0 1].';
cfg.asm = asm;
if cfg.isCC
    c = convenc(asm, cfg.trellis); c(2:2:end) = ~c(2:2:end);
    cfg.pat = logical(c(13:64)); cfg.patOff = 12;  % code bits 13..64 do not depend on the previous frame
else
    cfg.pat = logical(asm); cfg.patOff = 0;
end
% randomizer sequence: the generator's output for all-zero data (RS codeword of zeros = zeros)
g = makeGenerator(cfg, false, true);
[~, e] = g(zeros(cfg.Ni, 1));
cfg.PN = logical(e(33:end));
% reference power and matched-filter gain from a noiseless burst
nb = 6; bits = randi(RandStream('mlfg6331_64', 'Seed', 4242), [0 1], cfg.Ni*nb, 1);   % deterministic
O1 = txObjects(cfg);
[wf, enc] = O1.tm(bits);
cfg.fast8 = false; cfg.txGain = 1; cfg.tx8check = NaN;
if cfg.M == 8
    % The generator's 8PSK path is ~25x slower than BPSK/QPSK. Same code bits (BPSK generator),
    % same CCSDS 8PSK mapping (pskmod) and the same RRC (comm.RaisedCosineTransmitFilter) give
    % the identical waveform: checked here on every run, the generator is used if they differ.
    c8 = cfg; c8.fast8 = true;
    [w2, e2] = txWaveform(c8, txObjects(c8), bits);
    n = min(numel(wf), numel(w2)); k = (w2(1:n)'*wf(1:n))/(w2(1:n)'*w2(1:n));
    cfg.tx8check = max(abs(k*w2(1:n) - wf(1:n)))/rms(wf);
    if cfg.tx8check < 1e-6 && isequal(e2(1:numel(enc)), enc), cfg.fast8 = true; cfg.txGain = k; end
end
cfg.Ps = mean(abs(wf).^2);
mf = comm.RaisedCosineReceiveFilter('RolloffFactor', rolloff, 'FilterSpanInSymbols', cfg.span, ...
    'InputSamplesPerSymbol', sps, 'DecimationFactor', 1);
z = mf(wf); ys = z(cfg.span*sps + 1:sps:end);
nb3 = floor(numel(enc)/cfg.b)*cfg.b;
ref = pskmod(double(enc(1:nb3)), cfg.M, cfg.phOff, cfg.order, 'InputType', 'bit');
n = min(numel(ys), numel(ref)) - 1;
cfg.A = (ref(1:n)'*ys(1:n))/(ref(1:n)'*ref(1:n));
cfg.mfEVM = rms(ys(1:n)/cfg.A - ref(1:n));        % ISI of the truncated RRC pair (span 10)
nw = min(numel(wf), 2^15);
[p, f] = pwelch(wf(1:nw), hann(1024), 512, 1024, cfg.fs, 'centered');
cfg.txPSD = struct('f', f/Rs, 'dB', 10*log10(p/max(p)));
cfg.obw99 = obw(wf, cfg.fs)/Rs; cfg.obw999 = obw(wf, cfg.fs, [], 99.9)/Rs;
cfg.key = cacheKeyStr(cfg, imp);
cfg.keyNum = mod(sum(double(char(cfg.key)).*(1:strlength(cfg.key))), 1e6);
end

function O = txObjects(cfg)
% transmitter: ccsdsTMWaveformGenerator (BPSK/QPSK, or 8PSK reference), or the verified 8PSK fast path
if isfield(cfg, 'fast8') && cfg.fast8
    O.tm = makeGenerator(cfg, false, false, "BPSK");          % code bits only
    O.txf = comm.RaisedCosineTransmitFilter('RolloffFactor', cfg.rolloff, 'FilterSpanInSymbols', cfg.span, ...
        'OutputSamplesPerSymbol', cfg.sps);
else
    O.tm = makeGenerator(cfg, true, false);
end
end

function [wf, enc] = txWaveform(cfg, O, bits)
% one burst from the encoder reset state (frames are continuous inside a burst)
reset(O.tm);
if cfg.fast8
    [~, enc] = O.tm(bits);
    nb = floor(numel(enc)/3)*3;
    reset(O.txf);
    wf = cfg.txGain*O.txf(pskmod(double(enc(1:nb)), 8, cfg.phOff, cfg.order, 'InputType', 'bit'));
else
    [wf, enc] = O.tm(bits);
end
end

function g = makeGenerator(cfg, shaped, pnOnly, modName)
coding = cfg.genCoding;
if nargin < 4, modName = cfg.mod; end
if pnOnly
    coding = "none"; if cfg.isRS, coding = "RS"; end
    modName = "BPSK";
end
args = {'WaveformSource', 'synchronization and channel coding', 'HasRandomizer', true, 'HasASM', true, ...
    'ChannelCoding', char(coding), 'Modulation', char(modName)};
if shaped
    args = [args, {'PulseShapingFilter', 'root raised cosine', 'RolloffFactor', cfg.rolloff, ...
        'FilterSpanInSymbols', cfg.span, 'SamplesPerSymbol', cfg.sps}];
else
    args = [args, {'PulseShapingFilter', 'none'}];
end
if any(coding == ["RS" "concatenated"])
    args = [args, {'RSMessageLength', 223, 'RSInterleavingDepth', 1, 'IsRSMessageShortened', false}];
else
    args = [args, {'NumBytesInTransferFrame', 223}];
end
if any(coding == ["convolutional" "concatenated"])
    args = [args, {'ConvolutionalCodeRate', '1/2', 'InvertCCPath2', true}];
end
g = ccsdsTMWaveformGenerator(args{:});
end

function k = cacheKeyStr(cfg, imp)
k = sprintf('v2|%s|%s|a%.3f|sps%d|span%d|%s', cfg.mod, cfg.coding, cfg.rolloff, cfg.sps, cfg.span, cfg.receiver);
if ~isempty(imp.pn), k = [k sprintf('|pn%s@%s|Rs%.6g', mat2str(imp.pn(1, :)), mat2str(imp.pn(2, :)), cfg.Rs)]; end
if ~isempty(imp.iq), k = [k sprintf('|iq%.3g,%.3g', imp.iq)]; end
if cfg.realistic
    k = [k sprintf('|Rs%.6g|f%.6g|fd%.6g|tau%.4g|ppm%.4g|afc%.6g|bw%.4g,%.4g', cfg.Rs, imp.f0, imp.fdot, ...
        imp.tau, imp.ppm, imp.afc, imp.loop)];
end
k = string(k);
end

function g = defaultGrid(cfg, run)
s = run.step;
switch cfg.coding
    case "none", r = [0 11];
    case "CC", r = [0 6];
    case "RS", r = [3 8.5];
    otherwise, r = [0 4];
end
if cfg.M == 8, r = r + [2 3.5]; end
if cfg.realistic, r(2) = r(2) + 1; end
g = r(1):s:r(2);
end

function m = estimateMinutes(cfg, grid, run)
perFrame = 3e-3*cfg.Lc/4144*(1 + cfg.realistic)*run.sps/4 + 1e-3;
m = numel(grid)*run.Frames*perFrame/60*0.4;
end

%% ===== Monte Carlo ===============================================================

function names = countNames()
names = ["frames" "frameErr" "bits" "bitErr" "rawBits" "rawErr" "rsFail" "lost" "slips" "seconds" "afcRejects"];
end

function [S, cache] = sweep(cfg, imp, grid, run, opt, cache, label)
NC = numel(countNames()); nP = numel(grid);
[entry, ie] = cacheGet(cache, cfg.key, opt.Cache);
C = zeros(nP, NC); src = repmat("sim", nP, 1); need = true(nP, 1);
for i = 1:nP
    j = find(abs(entry.ebno - grid(i)) < 1e-9, 1);
    if ~isempty(j)
        C(i, :) = entry.C(j, :);
        if enough(C(i, :), run), need(i) = false; src(i) = "cache"; end
    end
end
snaps = cell(nP, 1);
S.notes = strings(0, 1);
t0 = tic;
if run.verbose
    fprintf('  sweep "%s": %d Eb/N0 points (%d cached), stop at %d packet errors or %d frames\n', ...
        label, nP, nnz(~need), run.MaxErrors, run.Frames);
end
usePar = run.useParallel && nnz(need) > 1 && poolReady();
if usePar
    pool = gcp; idx = find(need); nI = numel(idx);
    budget = run.TimeLimit*pool.NumWorkers/nI;
    gI = grid(idx); C0 = C(idx, :); Cn = zeros(nI, NC); Sn = cell(nI, 1);
    parfor q = 1:nI
        [Cn(q, :), Sn{q}] = simPoint(cfg, imp, gI(q), run, C0(q, :), budget, true);
    end
    C(idx, :) = C(idx, :) + Cn; snaps(idx) = Sn;
    S.notes(end+1, 1) = sprintf("Sweep %s: parfor on %d workers.", label, pool.NumWorkers);
else
    stopHigher = false;
    for i = 1:nP
        if stopHigher && C(i, 1) == 0, src(i) = "skipped"; continue, end
        if need(i)
            rest = nnz(need(i:end));
            budget = (run.TimeLimit - toc(t0))/max(rest, 1);
            if budget <= 0 && C(i, 1) > 0, src(i) = "time limit"; continue, end
            [cn, snaps{i}] = simPoint(cfg, imp, grid(i), run, C(i, :), max(budget, 0), true);
            C(i, :) = C(i, :) + cn;
            if run.verbose, printPoint(grid(i), C(i, :), cfg); end
        end
        if C(i, 2) == 0 && C(i, 1) >= run.MinFrames, stopHigher = true; end
    end
    if any(src == "skipped")
        S.notes(end+1, 1) = sprintf("Sweep %s: points above %.2f dB skipped (no packet errors at a lower Eb/N0).", ...
            label, grid(find(src == "skipped", 1)) - run.step);
    end
    if any(src == "time limit")
        S.notes(end+1, 1) = sprintf("Sweep %s: time limit reached; %d point(s) kept their cached counts only.", ...
            label, nnz(src == "time limit"));
    end
end
% cache update (accumulated counts)
for i = 1:nP
    if C(i, 1) == 0, continue, end
    j = find(abs(entry.ebno - grid(i)) < 1e-9, 1);
    if isempty(j), entry.ebno(end+1) = grid(i); entry.C(end+1, :) = C(i, :);
    else, entry.C(j, :) = C(i, :);
    end
end
entry.key = cfg.key;
if isempty(ie), cache(end+1) = entry; else, cache(ie) = entry; end
S.T = countsTable(cfg, grid, C, src);
cacheSave(cache, opt);                          % an interrupted run keeps the finished sweeps
S.snaps = snaps; S.cfg = cfg; S.run = run; S.label = string(label); S.elapsed = toc(t0);
end

function ok = enough(c, run)
ok = c(1) >= run.Frames || (c(2) >= run.MaxErrors && c(1) >= run.MinFrames);
end

function ok = poolReady()
ok = false;
try
    ok = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    if ok, p = gcp('nocreate'); if isempty(p), parpool; end, end
catch e
    warning('ff_ccsds_linksim:pool', 'No parallel pool (%s): running serially.', e.message); ok = false;
end
end

function printPoint(x, c, cfg)
fprintf('    Eb/N0 %5.2f dB: %6d frames, %5d packet errors (PER %8.2e), BER %8.2e, raw BER %8.2e', ...
    x, c(1), c(2), c(2)/max(c(1), 1), c(4)/max(c(3), 1), c(6)/max(c(5), 1));
if cfg.realistic, fprintf(', lost %d, slips %d, AFC rejected %d', c(8), c(9), c(11)); end
fprintf(', %.1f s\n', c(10));
end

function [cn, snap] = simPoint(cfg, imp, ebno, run, c0, budget, wantSnap)
% Simulates batches at one Eb/N0 until the stop rule (counted together with the cached counts c0).
seed = mod(run.seed*1000003 + round((ebno + 100)*1000)*7 + c0(1)*13 + cfg.keyNum, 2^31 - 2) + 1;
st = RandStream('mlfg6331_64', 'Seed', seed);
O = rxObjects(cfg, imp, run, seed);
esn0 = ebno + cfg.esFactor;
cn = zeros(1, numel(countNames())); snap = [];
t0 = tic;
while true
    [cb, sb] = runBatch(cfg, imp, esn0, run, st, O, wantSnap && isempty(snap));
    cn = cn + cb;
    if isempty(snap) && ~isempty(sb), snap = sb; end
    tot = c0 + cn;
    if enough(tot, run) || toc(t0) >= budget, break, end
end
cn(10) = cn(10) + toc(t0);
end

function O = rxObjects(cfg, imp, run, seed)
O.tx = txObjects(cfg);
O.pn = [];
if ~isempty(imp.pn)
    keep = imp.pn(2, :) < cfg.fs/2;
    if any(keep)
        O.pn = comm.PhaseNoise('Level', imp.pn(1, keep), 'FrequencyOffset', imp.pn(2, keep), ...
            'SampleRate', cfg.fs, 'RandomStream', 'mt19937ar with seed', 'Seed', mod(seed, 2^31));
    end
end
O.mfEye = comm.RaisedCosineReceiveFilter('RolloffFactor', cfg.rolloff, 'FilterSpanInSymbols', cfg.span, ...
    'InputSamplesPerSymbol', cfg.sps, 'DecimationFactor', 1);
if cfg.realistic
    names = struct('BPSK', 'BPSK', 'QPSK', 'QPSK', 'x8PSK', '8PSK');
    m = names.(matlab.lang.makeValidName(char(cfg.mod), 'Prefix', 'x'));
    nA = round(max(run.nAcq, 1)*cfg.symPerFrame*cfg.sps);
    O.nA = nA; O.afc = imp.afc;
    nA2 = 2*nA/cfg.sps;                                      % acquisition block at 2 sps
    O.cfc = comm.CoarseFrequencyCompensator('Modulation', m, 'Algorithm', 'FFT-based', ...
        'SampleRate', 2*cfg.Rs, 'FrequencyResolution', 2*cfg.Rs/2^nextpow2(nA2));
    resH = 2*cfg.Rs/2^nextpow2(nA2/2);
    O.cfcHalf = comm.CoarseFrequencyCompensator('Modulation', m, 'Algorithm', 'FFT-based', ...
        'SampleRate', 2*cfg.Rs, 'FrequencyResolution', resH);
    O.afcTol = 3*resH + 1e-3*cfg.Rs;
    O.mf = comm.RaisedCosineReceiveFilter('RolloffFactor', cfg.rolloff, 'FilterSpanInSymbols', cfg.span, ...
        'InputSamplesPerSymbol', cfg.sps, 'DecimationFactor', cfg.sps/2);
    O.ss = comm.SymbolSynchronizer('Modulation', 'PAM/PSK/QAM', 'TimingErrorDetector', 'Gardner (non-data-aided)', ...
        'SamplesPerSymbol', 2, 'DampingFactor', 1, 'NormalizedLoopBandwidth', imp.loop(1));
    O.cs = comm.CarrierSynchronizer('Modulation', m, 'ModulationPhaseOffset', 'Custom', ...
        'CustomPhaseOffset', cfg.carOff, 'SamplesPerSymbol', 1, 'DampingFactor', 0.707, ...
        'NormalizedLoopBandwidth', imp.loop(2));
    O.vfd = [];
else
    O.mf = O.mfEye;
end
end

function [c, s] = runBatch(cfg, imp, esn0, run, st, O, wantSnap)
c = zeros(1, numel(countNames())); s = [];
nF = run.nAcq + run.F + 1;                          % + 1 guard frame (filter / decoder tail)
bits = randi(st, [0 1], cfg.Ni*nF, 1);
[wf, enc] = txWaveform(cfg, O.tx, bits);
N = numel(wf); fs = cfg.fs;
x = wf;
if ~isempty(imp.iq), x = iqimbal(x, imp.iq(1), imp.iq(2)); end
if ~isempty(O.pn), x = O.pn(x); end
xDisp = x;                                          % before the carrier / timing offsets (display)
if cfg.realistic
    n = (0:N-1).'/fs;
    x = x.*exp(1j*(2*pi*rand(st) + 2*pi*(imp.f0*n + 0.5*imp.fdot*n.^2)));
    tau = imp.tau; if isnan(tau), tau = rand(st); end
    ramp = imp.ppm*1e-6*(0:N-1).';
    d = 2 + tau*cfg.sps + ramp - min(0, min(ramp));
    vfd = dsp.VariableFractionalDelay('InterpolationMethod', 'Farrow', 'MaximumDelay', ceil(max(d)) + 4);
    x = vfd(x, d);
end
y = awgn(x, esn0 - 10*log10(cfg.sps), 10*log10(cfg.Ps), st);
nv = 10^(-esn0/10);
if cfg.realistic
    [ys, fk, rk, info] = rxRealistic(cfg, O, y, nF, run.nAcq);
    c(9) = info.slips; c(11) = info.afcRejected;
else
    reset(O.mf);
    z = O.mf(y);
    ys = z(cfg.span*cfg.sps + 1:cfg.sps:end)/cfg.A;
    fk = (0:nF-1).'*cfg.Lc + 1; rk = zeros(nF, 1);
end
% decode the data frames
L = cell(1, cfg.nRot);
data = run.nAcq + (1:run.F);
for k = data
    c(1) = c(1) + 1;
    if isnan(fk(k)), c(2) = c(2) + 1; c(8) = c(8) + 1; continue, end
    r = rk(k) + 1;
    if isempty(L{r}), L{r} = demodLLR(ys*exp(-1j*2*pi*rk(k)/cfg.nRot), cfg, nv); end
    [info, fail, hard, ok] = decodeFrame(cfg, L{r}, fk(k));
    if ~ok, c(2) = c(2) + 1; c(8) = c(8) + 1; continue, end
    tx = bits((k-1)*cfg.Ni + (1:cfg.Ni));
    be = nnz(info(:) ~= tx);
    c(3) = c(3) + cfg.Ni; c(4) = c(4) + be; c(2) = c(2) + (be > 0); c(7) = c(7) + fail;
    c(5) = c(5) + cfg.Lc; c(6) = c(6) + nnz(hard(:) ~= logical(enc((k-1)*cfg.Lc + (1:cfg.Lc))));
end
if wantSnap
    k1 = data(1); r = rk(k1);
    i0 = max(1, floor((fk(k1) - 1)/cfg.b) + 1); if isnan(i0), i0 = 1; end
    s.sym = ys(i0:min(end, i0 + 2999))*exp(-1j*2*pi*r/cfg.nRot);
    reset(O.mfEye);
    ze = O.mfEye(xDisp + (y - x))/cfg.A;            % same noise, perfect carrier and timing
    e0 = cfg.span*cfg.sps + 1;
    s.eye = ze(e0:min(end, e0 + 400*cfg.sps - 1));
    nw = min(N, 2^15);
    [p, f] = pwelch(y(1:nw), hann(1024), 512, 1024, fs, 'centered');
    s.psd = struct('f', f/cfg.Rs, 'dB', 10*log10(p));
    s.rx = y(1:min(N, 2^14));
end
end

function L = demodLLR(ys, cfg, nv)
if cfg.M == 8
    L = pskdemod(ys, cfg.M, cfg.phOff, cfg.order, 'OutputType', 'approxllr', 'NoiseVariance', nv);
else
    L = pskdemod(ys, cfg.M, cfg.phOff, cfg.order, 'OutputType', 'llr', 'NoiseVariance', nv);
end
L = L(:);
end

function [info, fail, hard, ok] = decodeFrame(cfg, L, f)
NB = numel(L); info = []; fail = false; hard = [];
ok = f >= 1 && f + cfg.Lc - 1 <= NB;
if ~ok, return, end
if cfg.isCC
    a = f - cfg.Lpre; z = f + cfg.Lc + cfg.Lpost - 1;
    seg = zeros(z - a + 1, 1);
    src = max(a, 1):min(z, NB);
    seg(src - a + 1) = L(src);
    seg(2:2:end) = -seg(2:2:end);                   % undo the G2 inversion (LLR sign)
    dec = vitdec(seg, cfg.trellis, cfg.tblen, 'trunc', 'unquant');
    fb = logical(dec(cfg.Lpre/2 + (1:cfg.Lf)));
else
    fb = L(f + (0:cfg.Lf-1)) < 0;
end
tf = xor(fb(33:end), cfg.PN);
if cfg.isRS
    [info, ne] = ccsdsRSDecode(tf, 223);
    fail = ne < 0;
else
    info = tf;
end
hard = L(f + (0:cfg.Lc-1)) < 0;
end

function [ys, fk, rk, info] = rxRealistic(cfg, O, y, nF, nAcq)
N = numel(y); fs = cfg.fs;
info = struct('slips', 0, 'fEst', NaN, 'afcRejected', 0);
% AFC: FFT of the M-th power over the acquisition frames, clamped to the AFC range. Estimated on
% the matched-filtered 2-sps signal (on the raw 4-sps samples the out-of-band noise makes the QPSK
% 4th-power line vanish below ~3 dB Eb/N0 - verified), corrected at the input rate, then filtered again.
nA = min(N, O.nA);
reset(O.mf); za = O.mf(y(1:nA));
reset(O.cfc);
[~, fe] = O.cfc(za);
% consistency check (a real AFC also needs a stable estimate): the two halves of the block must
% agree, otherwise the M-th power line is buried in noise (8PSK, low SNR) and the estimate is a
% random bin up to +-fs/(2M) away - then no correction is applied and the loops must pull in.
h = floor(numel(za)/2);
reset(O.cfcHalf); [~, f1] = O.cfcHalf(za(1:h));
reset(O.cfcHalf); [~, f2] = O.cfcHalf(za(h+1:2*h));
if abs(f1 - f2) > O.afcTol, fe = 0; info.afcRejected = 1; end
fe = max(min(fe, O.afc), -O.afc); info.fEst = fe;
y = y.*exp(-1j*2*pi*fe*(0:N-1).'/fs);
reset(O.mf); reset(O.ss); reset(O.cs);
z = O.mf(y)/cfg.A;
z = O.ss(z);
ys = O.cs(z);
[fk, rk, info.slips] = frameSync(cfg, ys, nF, nAcq);
end

function [fk, rk, slips] = frameSync(cfg, ys, nF, nAcq)
% ASM search on hard bits for every phase rotation (resolves the PSK ambiguity).
fk = NaN(nF, 1); rk = zeros(nF, 1); slips = 0;
q = 1 - 2*double(cfg.pat(:)); Lq = numel(q);
NB = numel(ys)*cfg.b;
if NB < cfg.Lc*max(nAcq, 1) + Lq, return, end
Cr = zeros(NB - Lq + 1, cfg.nRot);
for r = 0:cfg.nRot - 1
    hb = pskdemod(ys*exp(-1j*2*pi*r/cfg.nRot), cfg.M, cfg.phOff, cfg.order, 'OutputType', 'bit');
    Cr(:, r+1) = conv(1 - 2*double(hb(:)), flipud(q), 'valid');
end
nC = size(Cr, 1);
P = min(cfg.Lc, nC);
met = zeros(P, cfg.nRot);
for j = 0:max(nAcq, 1) - 1
    idx = (1:P).' + j*cfg.Lc; ok = idx <= nC;
    met(ok, :) = met(ok, :) + Cr(idx(ok), :);
end
[mx, li] = max(met(:));
if mx < 0.4*max(nAcq, 1)*Lq, return, end           % acquisition failed: all frames lost
[p1, r1] = ind2sub(size(met), li);
f1 = p1 - cfg.patOff;
k1 = 1 + round((f1 - (cfg.span*cfg.b + 1))/cfg.Lc);
if k1 < 1 || k1 > nF, return, end
fk(k1) = f1; rk(k1) = r1 - 1;
for k = k1 + 1:nF
    e = fk(k-1) + cfg.Lc; pe = e + cfg.patOff;
    if pe + 2*cfg.b > nC || pe - 2*cfg.b < 1, break, end
    ce = Cr(pe, rk(k-1) + 1);
    best = -Inf; fb = e; rb = rk(k-1);
    for jj = -2:2
        pj = pe + jj*cfg.b;
        [v, r] = max(Cr(pj, :));
        if v > best, best = v; fb = e + jj*cfg.b; rb = r - 1; end
    end
    if best - ce >= Lq/4
        fk(k) = fb; rk(k) = rb; slips = slips + (fb ~= e);
    else
        fk(k) = e; rk(k) = rk(k-1);
    end
end
end

function T = countsTable(cfg, x, C, src)
x = x(:);
[per, perCI] = berconfint(C(:, 2), max(C(:, 1), 1), 0.95);
[ber, berCI] = berconfint(C(:, 4), max(C(:, 3), 1), 0.95);
none = C(:, 1) == 0;
per(none) = NaN; perCI(none, :) = NaN; ber(none) = NaN; berCI(none, :) = NaN;
raw = C(:, 6)./C(:, 5); raw(C(:, 5) == 0) = NaN;
T = table(x, x + cfg.payFactor, x + cfg.esFactor, x + cfg.chFactor, C(:, 1), C(:, 2), per, perCI(:, 1), perCI(:, 2), ...
    C(:, 3), C(:, 4), ber, berCI(:, 1), berCI(:, 2), C(:, 5), C(:, 6), raw, C(:, 7), C(:, 8), C(:, 9), C(:, 11), C(:, 10), src(:), ...
    'VariableNames', {'EbNo_dB', 'EbNo_payload_dB', 'EsN0_dB', 'EbNo_channel_dB', 'Frames', 'PacketErrors', ...
    'PER', 'PER_lo', 'PER_hi', 'Bits', 'BitErrors', 'BER', 'BER_lo', 'BER_hi', 'RawBits', 'RawErrors', 'RawBER', ...
    'RSFailures', 'FramesLost', 'Slips', 'AFCRejects', 'Seconds', 'Source'});
end

%% ===== cache =====================================================================

function f = cacheFileName(opt)
f = opt.CacheFile;
if strlength(f) == 0
    f = fullfile(fileparts(mfilename('fullpath')), 'support', 'linksim', 'linksim_cache.mat');
end
end

function cache = cacheOpen(opt)
cache = struct('key', {}, 'ebno', {}, 'C', {});
if isequal(opt.Cache, false), return, end
f = cacheFileName(opt);
if isfile(f)
    try
        s = load(f, 'cache');
        if isfield(s, 'cache') && isstruct(s.cache) && all(isfield(s.cache, {'key', 'ebno', 'C'}))
            cache = s.cache;
            cache = cache(arrayfun(@(e) size(e.C, 2) == numel(countNames()), cache));  % older formats
        end
    catch e
        warning('ff_ccsds_linksim:cache', 'Cache %s not readable (%s): starting empty.', f, e.message);
    end
end
end

function [entry, ie] = cacheGet(cache, key, mode)
% ie = index of the key in the cache (also for "refresh", so the entry is overwritten, not duplicated)
entry = struct('key', key, 'ebno', zeros(1, 0), 'C', zeros(0, numel(countNames())));
ie = [];
if isempty(cache), return, end
ie = find([cache.key] == key, 1);
if isempty(ie), return, end
if (isstring(mode) || ischar(mode)) && lower(string(mode)) == "refresh", return, end
if isequal(mode, false), return, end
entry = cache(ie);
end

function cacheSave(cache, opt)
if isequal(opt.Cache, false), return, end
f = cacheFileName(opt);
try
    d = fileparts(f); if ~isfolder(d), mkdir(d); end
    save(f, 'cache');
catch e
    warning('ff_ccsds_linksim:cache', 'Could not write the cache %s (%s).', f, e.message);
end
end

%% ===== theory and summaries ===========================================================

function th = theory(cfg, x)
x = x(:); M = cfg.M; n = numel(x);
th = table(x, berPSK(x, M), NaN(n, 1), NaN(n, 1), ...
    'VariableNames', {'EbNo_dB', 'BER_uncodedPSK', 'BER_system', 'PER_system'});
spect = distspec(cfg.trellis, 7);
xrs = x - 10*log10(cfg.Lf/cfg.Lt);                 % per RS / CC info bit, ASM overhead removed
pch = berPSK(x + cfg.chFactor, M);
rsPER = @(p) binocdf(16, 255, -expm1(8*log1p(-min(p, 0.5))), 'upper');
switch cfg.coding
    case "none"
        th.BER_system = pch; th.PER_system = -expm1(cfg.Ni*log1p(-pch));
        lab = "exact: berawgn at the channel Eb/N0 (ASM overhead included)";
    case "CC"
        if M <= 4
            th.BER_system = min(bercoding(xrs, 'conv', 'soft', 1/2, spect), 0.5);
            lab = "bercoding soft-decision union bound (pessimistic above BER ~1e-4)";
        else
            lab = "8PSK + CC: no closed form, simulation only";
        end
    case "RS"
        if M <= 4, th.BER_system = bercoding(xrs, 'RS', 'hard', 255, 223); end
        th.PER_system = rsPER(pch);
        lab = "bercoding RS hard decision; PER = P(> 16 byte errors), binocdf";
    otherwise
        if M <= 4
            pcc = min(bercoding(x + 10*log10(cfg.Ni/cfg.Lf), 'conv', 'soft', 1/2, spect), 0.5);
            th.PER_system = rsPER(pcc);
            lab = "approx: CC soft union bound (pessimistic near threshold) -> independent byte errors (optimistic for I = 1: Viterbi bursts) -> RS t = 16";
        else
            lab = "8PSK + CC: no closed form, simulation only";
        end
end
th.Properties.Description = lab;
end

function s = summarize(cfg, S, SI, NF)
T = S.T; s = struct(); s.notes = strings(0, 1);
[s.EbNo_PER1, s.EbNo_PER1_CI] = crossing(T, 'PER', 1e-2);
[s.EbNo_BER5, s.EbNo_BER5_CI] = crossing(T, 'BER', 1e-5);
s.EbNo_PER1_payload = s.EbNo_PER1 + cfg.payFactor;
s.uncodedEbNo_BER5 = invTheory(1e-5, cfg.M);
s.codingGain_BER5 = s.uncodedEbNo_BER5 - s.EbNo_BER5;
s.shannonEbNo = cfg.shannon;
s.shannonGap_BER5 = s.EbNo_BER5 - cfg.shannon;
s.reqCN0_PER1_dBHz = s.EbNo_PER1 + 10*log10(cfg.RbInfo);
s.sensitivity_PER1_dBm = s.reqCN0_PER1_dBHz - 174 + NF;   % SRS-4-class receiver
s.NF = NF;
% demodulator loss from the raw BER (theory: berawgn at the channel Eb/N0)
[s.rawLoss_dB, s.rawLoss_points, s.rawLoss_range] = rawLoss(T, cfg.M);
s.syncLoss_PER1 = NaN;
if ~isempty(SI)
    s.idealEbNo_PER1 = crossing(SI.T, 'PER', 1e-2);
    s.syncLoss_PER1 = s.EbNo_PER1 - s.idealEbNo_PER1;
end
if isnan(s.EbNo_PER1)
    s.notes(end+1, 1) = "PER = 1 % is not crossed inside the simulated Eb/N0 range: widen 'EbNo'.";
end
s.obw99_Rs = cfg.obw99; s.obw999_Rs = cfg.obw999;
end

function [x0, ci] = crossing(T, var, target)
% Eb/N0 where an error-rate curve crosses target (log-linear interpolation) and its 95 % interval
% from the berconfint bounds. Zero-error points use 0.5/N (a pseudo-count) to stay on a log scale;
% when the crossing lies next to a zero-error point the lower CI end is widened to the last point
% that still had errors (the curve between them is not measured).
ok = T.Frames > 0;
x = T.EbNo_dB(ok); n = T.Frames(ok);
if var == "BER", n = T.Bits(ok); end
fl = 0.5./max(n, 1);
v = max(T.(var)(ok), fl); lo = max(T.(var + "_lo")(ok), fl); hi = max(T.(var + "_hi")(ok), fl);
[x0, i] = cross1(x, v, target); ci = [cross1(x, lo, target) cross1(x, hi, target)];
raw = T.(var)(ok);
if ~isempty(i) && ~(raw(i+1) > 0)                         % upper bracketing point has no errors
    ci(1) = min(ci(1), x(i));
end
end

function [x0, i] = cross1(x, v, target)
x0 = NaN; i = [];
if numel(x) < 2, return, end
i = find(v(1:end-1) >= target & v(2:end) < target, 1);
if isempty(i), return, end
a = log10(v(i)); b = log10(v(i+1));
x0 = x(i) + (log10(target) - a)/(b - a)*(x(i+1) - x(i));
end

function xi = invTheory(p, M)
% Eb/N0 (channel) at which uncoded M-PSK theory gives BER p
[xs, lt] = berTable(M);
ok = isfinite(lt); [lt, iu] = unique(lt(ok)); xs = xs(ok); xs = xs(iu);
xi = interp1(lt, xs, log10(p), 'linear', NaN);
end

function b = berPSK(x, M)
% berawgn(x, 'psk', M, 'nondiff'). BPSK/QPSK: closed form, direct. 8PSK: berawgn integrates
% numerically (~2 ms per point, profiled: 95 s per run before), so a 0.05 dB table is computed
% once per MATLAB session and interpolated in log10 (error < 1e-3 relative).
if M <= 4, b = berawgn(x, 'psk', M, 'nondiff'); return, end
[xs, lt] = berTable(M);
ok = isfinite(lt);
b = 10.^interp1(xs(ok), lt(ok), x, 'pchip', NaN);
b(x < xs(1)) = berawgn(x(x < xs(1)), 'psk', M, 'nondiff');
b(x > xs(end)) = 0;
end

function [xs, lt] = berTable(M)
persistent tab
if isempty(tab), tab = containers.Map('KeyType', 'double', 'ValueType', 'any'); end
if ~isKey(tab, M)
    xs = (-8:0.05:24).';
    t = berawgn(xs, 'psk', M, 'nondiff');
    lt = log10(t); lt(t <= 1e-300) = -Inf;
    tab(M) = {xs, lt};
end
v = tab(M); xs = v{1}; lt = v{2};
end

function [m, n, rng] = rawLoss(T, M)
ok = T.RawErrors >= 20 & T.RawBER < 0.3;
m = NaN; n = nnz(ok); rng = [NaN NaN];
if n == 0, return, end
rng = [min(T.EbNo_channel_dB(ok)) max(T.EbNo_channel_dB(ok))];
xth = invTheory(T.RawBER(ok), M);
m = mean(T.EbNo_channel_dB(ok) - xth, 'omitnan');
end

function C = calibrationTable(S, SI, cfg)
C = table();
add = @(name, T) table(string(name), mean(T.EbNo_channel_dB(T.RawErrors >= 20 & T.RawBER < 0.3)), ...
    rawLoss(T, cfg.M), nnz(T.RawErrors >= 20 & T.RawBER < 0.3), ...
    'VariableNames', {'Sweep', 'MeanChannelEbNo_dB', 'RawLoss_dB', 'Points'});
C = [C; add(S.label + " (" + S.cfg.receiver + ")", S.T)];
if ~isempty(SI), C = [C; add(SI.label + " (ideal)", SI.T)]; end
end

%% ===== datasheet validation ========================================================

function [T, note, cache, Sv] = validateFig13(cfg, imp, run, opt, cache)
T = table(); note = ""; Sv = [];
if cfg.M == 8
    note = "Figure 13 validation skipped: the SRS-4 receiver (and Figure 13) is BPSK/QPSK only.";
    return
end
cfgV = linkConfig(cfg.mod, "none", cfg.rolloff, run.sps, cfg.Rs, "realistic", imp);
runV = presetRun(opt, "realistic");
runV.Frames = run.fig13Frames; runV.MinFrames = min(runV.MinFrames, runV.Frames);
runV.TimeLimit = runV.TimeLimit/2;                       % secondary sweep: half budget
runV.MaxErrors = max(runV.MaxErrors, 50);
x = run.fig13 - cfgV.chFactor;                          % info Eb/N0 for the requested channel Eb/N0
[Sv, cache] = sweep(cfgV, imp, x, runV, opt, cache, "Figure 13 (uncoded, realistic)");
V = Sv.T;
[raw, rawCI] = berconfint(V.RawErrors, max(V.RawBits, 1), 0.95);
th = berPSK(V.EbNo_channel_dB, cfgV.M);
loss = V.EbNo_channel_dB - invTheory(raw, cfgV.M);
lossLo = V.EbNo_channel_dB - invTheory(rawCI(:, 1), cfgV.M);  % lower BER bound -> smallest loss
lossHi = V.EbNo_channel_dB - invTheory(rawCI(:, 2), cfgV.M);  % upper BER bound -> largest loss
lossLo(rawCI(:, 1) == 0) = -Inf;
verdict = repmat("INFO", height(V), 1);
ok = V.RawErrors >= 20;
verdict(ok & lossHi <= 0.3) = "PASS";                          % whole 95 % interval below 0.3 dB
verdict(ok & lossLo > 0.3) = "FAIL";                           % whole interval above 0.3 dB
verdict(ok & lossLo <= 0.3 & lossHi > 0.3) = "UNCERTAIN";      % interval straddles 0.3 dB: more frames
verdict(~ok) = "TOO FEW ERRORS";
T = table(V.EbNo_channel_dB, V.RawBits, V.RawErrors, raw, rawCI(:, 1), rawCI(:, 2), th, loss, lossLo, lossHi, verdict, ...
    'VariableNames', {'EbNo_channel_dB', 'Bits', 'Errors', 'BER', 'BER_lo', 'BER_hi', 'BER_theory', ...
    'Loss_dB', 'Loss_lo_dB', 'Loss_hi_dB', 'Verdict_vs_0p3dB'});
good = ok & isfinite(loss);
note = sprintf("Figure 13 check (%s uncoded, realistic receiver, impairments ""%s""): demodulator loss %.2f dB mean (%.2f to %.2f dB over %d points with >= 20 errors) vs datasheet < 0.3 dB.", ...
    cfg.mod, imp.preset, mean(loss(good)), min(loss(good)), max(loss(good)), nnz(good));
if ~any(good), note = "Figure 13 check: too few raw bit errors for a loss estimate (raise 'Frames')."; end
end

function T = table7Data()
rates = [100e3 1e6 4e6 5e6];
P = cat(3, [-122 -121 -114 -112; -119 -118 -111 -109], [-112.5 -111 -104 -102; -108.5 -108 -101 -99], ...
    [-106 -105 -98 -96; -102 -101.5 NaN -93], [-105.5 -104 NaN -95; NaN -100.5 NaN -92]);
T = struct('rates', rates, 'P', P, 'mods', ["BPSK" "QPSK"], 'cols', ["concatenated" "CC" "RS" "none"]);
end

function T = table7Compare(cfg, x1, ci, NF)
% SRS-4 Table 7 (PER < 1 %, 217-byte payload) converted to Eb/N0 per payload bit and compared.
T = table();
D = table7Data();
im = find(D.mods == cfg.mod); ic = find(D.cols == cfg.coding);
if isempty(im) || isempty(ic), return, end
for ir = 1:numel(D.rates)
    P = D.P(im, ic, ir);
    if isnan(P), continue, end
    Rs = D.rates(ir); rbPay = cfg.Np*Rs*cfg.b/cfg.Lc;
    ebDs = P + 174 - NF - 10*log10(rbPay);
    ebSim = x1 + cfg.payFactor;
    pSim = ebSim + 10*log10(rbPay) - 174 + NF;
    T = [T; table(cfg.mod, cfg.coding, Rs/1e3, P, ebDs, ebSim, ci(1) + cfg.payFactor, ci(2) + cfg.payFactor, ...
        ebDs - ebSim, pSim, string(cfg.receiver), ...
        'VariableNames', {'Modulation', 'Coding', 'Rate_kBd', 'P_table_dBm', 'EbNo_table_dB', 'EbNo_sim_dB', ...
        'EbNo_sim_lo_dB', 'EbNo_sim_hi_dB', 'Difference_dB', 'P_sim_dBm', 'Receiver'})]; %#ok<AGROW>
end
end

function [T, cache] = table7Sweep(opt, imp, run, cache)
T = table();
for m = ["BPSK" "QPSK"]
    for c = ["none" "CC" "RS" "concatenated"]
        cfgk = linkConfig(m, c, opt.Rolloff, run.sps, opt.SymbolRate, opt.Receiver, imp);
        [Sk, cache] = sweep(cfgk, imp, defaultGrid(cfgk, run), presetRun(opt), opt, cache, "Table 7: " + m + " " + c);
        [x1, ci] = crossing(Sk.T, 'PER', 1e-2);
        T = [T; table7Compare(cfgk, x1, ci, opt.NoiseFigure)]; %#ok<AGROW>
        cacheSave(cache, opt);
    end
end
end

function note = table7Note(T, NF)
note = "";
if isempty(T), note = "Table 7 comparison: no Table 7 column for this modulation/coding (8PSK is TX only)."; return, end
d = T.Difference_dB(isfinite(T.Difference_dB));
if isempty(d), note = "Table 7 comparison: the simulated PER did not cross 1 %."; return, end
note = sprintf("Table 7 (NF %.2f dB, payload Eb/N0): the SRS-4 needs %.2f to %.2f dB (mean %.2f dB) more Eb/N0 than this %s-receiver model at PER 1 %% - hardware implementation loss + NF tolerance + test set-up.", ...
    NF, min(d), max(d), mean(d), T.Receiver(1));
cods = unique(T.Coding, 'stable');
if numel(cods) > 1                                   % 'Table7Sweep': where is the loss?
    m = arrayfun(@(c) mean(T.Difference_dB(T.Coding == c & isfinite(T.Difference_dB))), cods);
    note = note + " Mean difference per FEC: " + strjoin(compose("%s %+.2f dB", cods, m), ", ") + ".";
    noRS = ismember(cods, ["none" "CC"]);
    if all(abs(m(noRS)) < 0.5) && all(m(~noRS) > 0.5)
        note = note + " Uncoded and CC agree within 0.5 dB (demodulator + Viterbi + NF consistent with the datasheet); the RS modes need more: the extra loss sits in the SRS-4 RS path, which this model does not explain.";
    end
end
end

%% ===== pass mode =================================================================

function P = loadPass(src, opt, interactive)
P = struct('passes', [], 'user', [], 'source', "");
if isequal(src, true)
    if ~(opt.Prompt && interactive)
        error('ff_ccsds_linksim:pass', '''Pass'', true needs a dialog: give a struct or a MAT file instead.');
    end
    [f, d] = uigetfile('cubesat_link.mat', 'ff_cubesat_link result');
    if isequal(f, 0), error('ff_ccsds_linksim:cancelled', 'No file chosen.'); end
    src = string(fullfile(d, f));
end
if (isstring(src) || ischar(src)) && lower(string(src)) == "demo"
    P = demoPass(opt); return
end
if isstring(src) || ischar(src)
    s = load(src);
    if isfield(s, 'R'), Rc = s.R; else, error('ff_ccsds_linksim:pass', '%s has no variable R (ff_cubesat_link).', src); end
    P.source = string(src);
else
    Rc = src; P.source = "ff_cubesat_link result (struct)";
end
if ~isfield(Rc, 'pass') || isempty(Rc.pass)
    error('ff_ccsds_linksim:pass', 'The ff_cubesat_link result has no passes (R.pass).');
end
if ~isfield(Rc.pass, 'ebnoUser')
    error('ff_ccsds_linksim:pass', 'R.pass has no ebnoUser (Eb/N0 of your configuration): re-run ff_cubesat_link.');
end
for k = 1:numel(Rc.pass)
    p = Rc.pass(k);
    P.passes(k).t = p.t(:); P.passes(k).el = p.el(:); P.passes(k).ebno = p.ebnoUser(:);
    P.passes(k).maxEl = max(p.el);
end
if isfield(Rc, 'user') && ~isempty(Rc.user)
    u = Rc.user;
    if ~ismember(string(u.Modulation), ["BPSK" "QPSK" "8PSK"])
        error('ff_ccsds_linksim:pass', ['The ff_cubesat_link configuration is %s: this simulator models the ' ...
            'S-band SRS-4 (BPSK/QPSK/8PSK) only.'], string(u.Modulation));
    end
    P.user = struct('Modulation', string(u.Modulation), 'Coding', canonicalCoding(u.FEC), ...
        'SymbolRate', double(u.SymbolRate_Bd), 'Rb', double(u.Rb_bps), 'EbNoReq', double(u.EbN0_req_dB));
end
end

function P = demoPass(opt)
% Three synthetic passes (500 km circular orbit, no Earth rotation) with a small ground antenna
% (14 dBic helix, Tsys 300 K: G/T -10.8 dB/K) so that the waterfall shows during the pass.
mu = 3.986004418e14; Re = 6371e3; h = 500e3; r = Re + h; w = sqrt(mu/r^3);
f = 2245e6; lam = physconst('LightSpeed')/f;
eirp_dBW = 33 + 6 - 1 - 30; GT = 14 - 10*log10(300); losses = 3 + 1;
modName = opt.Modulation; cod = opt.Coding; Rs = opt.SymbolRate;
bps = find(["BPSK" "QPSK" "8PSK"] == modName);                % 1, 2, 3 bits/symbol
Lt = 1784; if any(cod == ["RS" "concatenated"]), Lt = 2040; end
Lc = (32 + Lt)*(1 + any(cod == ["CC" "concatenated"]));
Rb = 1736*Rs*bps/Lc;                                     % payload bit rate (217-byte payload)
t0 = datetime(2026, 10, 10, 2, 0, 0, 'TimeZone', 'UTC');
elMax = [25 50 85];
for k = 1:3
    g0 = acos(Re*cosd(elMax(k))/r) - deg2rad(elMax(k));
    tt = (-700:700).';
    th = w*tt; cg = cos(th).*cos(g0);
    g = acos(cg);
    el = atan2d(cg - Re/r, sin(g));
    rho = sqrt(Re^2 + r^2 - 2*Re*r*cg);
    keep = el >= 10;
    cn0 = eirp_dBW - fspl(rho(keep), lam) + GT + 228.6 - losses;
    P.passes(k).t = t0 + hours(1.6*(k-1)) + seconds(tt(keep));
    P.passes(k).el = el(keep);
    P.passes(k).ebno = cn0 - 10*log10(Rb);
    P.passes(k).maxEl = max(el(keep));
end
P.user = struct('Modulation', modName, 'Coding', cod, 'SymbolRate', Rs, 'Rb', Rb, 'EbNoReq', NaN);
P.source = "demo (500 km orbit, 14 dBic helix, G/T -10.8 dB/K, EIRP 38 dBm)";
end

function PS = passPrediction(P, S, cfg, summ, opt)
PS = struct('notes', strings(0, 1), 'table', table(), 'series', []);
n = numel(P.passes);
rows = cell(n, 1);
thr = summ.EbNo_PER1_payload + opt.MarginTarget;
for k = 1:n
    p = P.passes(k);
    x = p.ebno - cfg.payFactor - opt.ImplLoss;            % info Eb/N0 seen by the simulated receiver
    [per, flag, perHi] = curveAt(S.T, 'PER', x);
    dt = seconds(median(diff(p.t))); if isempty(dt) || ~isfinite(dt), dt = 1; end
    tput = cfg.RbPay*(1 - per);
    vol = sum(tput)*dt/8e6;
    volThr = nnz(p.ebno >= thr)*dt*cfg.RbPay/8e6;
    s = struct('t', p.t, 'el', p.el, 'ebnoPayload', p.ebno, 'ebnoInfo', x, 'PER', per, 'PER_hi', perHi, ...
        'flag', flag, 'throughput', tput);
    PS.series = [PS.series; s];
    aos = p.t(1); if isdatetime(aos) && ~isempty(aos.TimeZone), aos.TimeZone = 'UTC'; end   % ff_cubesat_link: UTC
    rows{k} = table(k, aos, p.maxEl, numel(p.t)*dt/60, min(p.ebno), max(p.ebno), mean(per), ...
        nnz(per < 0.01)*dt/60, vol, volThr, nnz(flag == 1)*dt/60, ...
        'VariableNames', {'Pass', 'AOS_UTC', 'MaxEl_deg', 'Duration_min', 'MinEbNo_dB', 'MaxEbNo_dB', 'MeanPER', ...
        'TimePER1_min', 'Volume_MB', 'VolumeThreshold_MB', 'AboveSimRange_min'});
end
PS.table = vertcat(rows{:});
PS.notes(end+1, 1) = sprintf("Pass mode (%s): %s %s at %s, payload %.1f kbit/s; Eb/N0 (payload) from the pass minus %.2f dB (payload->info) minus ImplLoss %.2f dB.", ...
    P.source, cfg.mod, cfg.coding, rateText(cfg.Rs), cfg.RbPay/1e3, cfg.payFactor, opt.ImplLoss);
PS.notes(end+1, 1) = sprintf("Threshold model for comparison: full rate when Eb/N0(payload) >= %.2f dB (PER 1 %% point) + %.1f dB margin.", ...
    summ.EbNo_PER1_payload, opt.MarginTarget);
if any(arrayfun(@(s) any(s.flag == 1), PS.series))
    PS.notes(end+1, 1) = "Above the simulated Eb/N0 range the PER is held at the last point (0.5/frames when no errors were simulated; the 95 % upper bound is in series.PER_hi): the throughput there is ~full rate.";
end
end

function [v, flag, vHi] = curveAt(T, var, x)
% Error rate at Eb/N0 x from the simulated curve (log-linear interpolation). Zero-error points
% use the pseudo-count 0.5/N; vHi = the same interpolation of the 95 % upper bound (berconfint),
% the honest number to quote where no errors were simulated. flag: -1 below the simulated range
% (held at the first point), +1 above it (held at the last point), 0 inside.
ok = T.Frames > 0;
xs = T.EbNo_dB(ok); n = T.Frames(ok); if var == "BER", n = T.Bits(ok); end
vs = max(T.(var)(ok), 0.5./max(n, 1)); hs = max(T.(var + "_hi")(ok), 0.5./max(n, 1));
flag = zeros(size(x)); v = NaN(size(x)); vHi = v;
if isempty(xs), return, end
flag(x < xs(1)) = -1; flag(x > xs(end)) = 1;
if isscalar(xs)
    v(:) = vs; vHi(:) = hs; return
end
xc = min(max(x, xs(1)), xs(end));
v = 10.^interp1(xs, log10(vs), xc, 'linear');
vHi = 10.^interp1(xs, log10(hs), xc, 'linear');
end

%% ===== measured mode ==============================================================

function M = loadMeasured(src, opt, interactive)
if isequal(src, true)
    if ~(opt.Prompt && interactive)
        error('ff_ccsds_linksim:measured', '''Measured'', true needs a dialog: give a struct or a MAT file instead.');
    end
    [f, d] = uigetfile('srs4_run.mat', 'ff_srs4_analyzer result');
    if isequal(f, 0), error('ff_ccsds_linksim:cancelled', 'No file chosen.'); end
    src = string(fullfile(d, f));
end
if (isstring(src) || ischar(src)) && lower(string(src)) == "demo"
    lam = physconst('LightSpeed')/2245.01e6;
    A.power = struct('C_dBm', -62, 'N0_dBmHz', -161, 'CN0', 99, 'EsN0', 39, 'fspl', fspl(3, lam));
    A.spec = struct('Rs', 1.0002e6, 'fc', 2245.01e6);
    M = struct('A', A, 'source', "demo (C = -62 dBm, N0 = -161 dBm/Hz at the FieldFox, 3 m over the air)");
    return
end
if isstring(src) || ischar(src)
    s = load(src);
    if isfield(s, 'results'), A = s.results;
    elseif isfield(s, 'R') && isfield(s.R, 'results'), A = s.R.results;
    else, error('ff_ccsds_linksim:measured', '%s has no ''results'' (ff_srs4_analyzer srs4_run.mat).', src);
    end
    M.source = string(src);
else
    if isfield(src, 'results'), A = src.results; else, A = src; end
    M.source = "ff_srs4_analyzer result (struct)";
end
if ~isfield(A, 'power') || ~isfield(A.power, 'CN0') || ~isfinite(A.power.CN0)
    error('ff_ccsds_linksim:measured', 'The analyzer result has no C/N0 (results.power.CN0): was the downlink found?');
end
M.A = A;
end

function MS = measuredPrediction(Min, S, cfg, summ, opt)
A = Min.A; MS = struct('notes', strings(0, 1), 'table', table());
Rs = cfg.Rs; fc = 2245e6;
if isfield(A, 'spec') && isfield(A.spec, 'Rs') && isfinite(A.spec.Rs), Rs = A.spec.Rs; end
if isfield(A, 'spec') && isfield(A.spec, 'fc') && isfinite(A.spec.fc), fc = A.spec.fc; end
rbInfo = cfg.Ni*Rs*cfg.b/cfg.Lc;
C = A.power.C_dBm;
cn0 = [A.power.CN0; C + 174 - opt.NoiseFigure];
who = ["FieldFox (as measured)"; sprintf("SRS-4-class receiver, NF %.2f dB", opt.NoiseFigure)];
x = cn0 - 10*log10(rbInfo) - opt.ImplLoss;
[per, fp, perHi] = curveAt(S.T, 'PER', x); [ber, fb, berHi] = curveAt(S.T, 'BER', x);
margin = x - summ.EbNo_PER1;
d0 = NaN;
if isfield(A.power, 'fspl') && isfinite(A.power.fspl)
    d0 = 10^(A.power.fspl/20)*physconst('LightSpeed')/fc/(4*pi);
end
dmax = d0*10.^(margin/20);
note = repmat("", 2, 1);
note(fp == 1) = "above the simulated range: quote PER/BER_upper95 (no errors simulated there)";
note(fp == -1) = "below the simulated range: PER held at the lowest point";
MS.table = table(who, cn0, x, ber, berHi, per, perHi, margin, dmax, note, ...
    'VariableNames', {'Receiver', 'CN0_dBHz', 'EbNo_info_dB', 'BER', 'BER_upper95', 'PER', 'PER_upper95', ...
    'MarginToPER1_dB', 'EquivalentRange_m', 'Note'});
MS.Rs = Rs; MS.fc = fc; MS.C_dBm = C; MS.distance_m = d0; MS.fb = fb;
MS.notes(end+1, 1) = sprintf("Measured (%s): C = %.1f dBm, Rs = %s -> %s %s info rate %.1f kbit/s. The FieldFox C/N0 is limited by its own noise floor (N0 %.1f dBm/Hz): the SRS-4-class row is what an SRS-4 receiver would see at the same input power.", ...
    Min.source, C, rateText(Rs), cfg.mod, cfg.coding, rbInfo/1e3, A.power.N0_dBmHz);
if isfinite(d0)
    MS.notes(end+1, 1) = sprintf("Equivalent range = measured distance %.3g m x 10^(margin/20) (free space, same antennas): how far the link would still give PER < 1 %%.", d0);
end
if isfield(A, 'floor') && isfield(A.floor, 'upperBound') && A.floor.upperBound
    MS.notes(end+1, 1) = "ff_srs4_analyzer took N0 from the zoom-span edges (no mask sweep; PA regrowth may dominate there): N0 is an upper bound, so the FieldFox C/N0, Eb/N0 and margin are lower bounds.";
end
end

%% ===== figures ===================================================================

function files = makeFigures(cfg, imp, S, SI, summ, V, PS, MS, opt, interactive)
files = strings(0, 1);
T = S.T;
ttl = sprintf('SRS-4 CCSDS TM: %s, %s, \\alpha = %.2f, %s receiver', cfg.mod, cfg.codingName, cfg.rolloff, cfg.receiver);
% 1. BER / PER with theory
f1 = figure('Name', 'Link BER/PER', 'Color', 'w', 'Position', [60 60 1150 520]);
tl = tiledlayout(f1, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); title(tl, ttl);
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); set(ax, 'YScale', 'log');
xt = linspace(min(T.EbNo_dB) - 1, max(max(T.EbNo_dB) + 1, 11), 300).';
plot(ax, xt, berPSK(xt, cfg.M), 'k--', 'DisplayName', sprintf('uncoded %s theory (berawgn)', cfg.mod));
tht = theory(cfg, xt);
if any(isfinite(tht.BER_system)), plot(ax, xt, tht.BER_system, '-', 'Color', [0.5 0.5 0.5], 'DisplayName', 'system BER theory'); end
if any(isfinite(tht.PER_system)), plot(ax, xt, tht.PER_system, ':', 'Color', [0.6 0.2 0.2], 'DisplayName', 'system PER theory'); end
ok = T.Frames > 0;
errBars(ax, T.EbNo_dB(ok), T.BER(ok), T.BER_lo(ok), T.BER_hi(ok), 'o-', [0 0.35 0.75], 'simulated BER (95 % CI)');
errBars(ax, T.EbNo_dB(ok), T.PER(ok), T.PER_lo(ok), T.PER_hi(ok), 's-', [0.85 0.3 0.1], 'simulated PER (95 % CI)');
if ~isempty(SI)
    oi = SI.T.Frames > 0;
    plot(ax, SI.T.EbNo_dB(oi), SI.T.PER(oi), 's--', 'Color', [0.85 0.3 0.1 0.5], 'DisplayName', 'PER, ideal receiver');
end
yline(ax, 1e-2, ':', 'PER 1 %', 'HandleVisibility', 'off'); yline(ax, 1e-5, ':', 'BER 10^{-5}', 'HandleVisibility', 'off');
xline(ax, cfg.shannon, '-.', sprintf('Shannon %.2f dB', cfg.shannon), 'Color', [0.3 0.6 0.3], 'HandleVisibility', 'off');
ylim(ax, [1e-7 1]); xlim(ax, [xt(1) xt(end)]);
xlabel(ax, 'E_b/N_0 per information bit [dB] (ASM + FEC overhead included)'); ylabel(ax, 'error rate');
legend(ax, 'Location', 'southwest'); title(ax, 'Coded BER and packet error rate');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); set(ax, 'YScale', 'log');
xc = linspace(-4, 11, 300).';
plot(ax, xc, berPSK(xc, cfg.M), 'k--', 'DisplayName', 'theory (berawgn)');
plot(ax, xc + 0.3, berPSK(xc, cfg.M), 'k:', 'DisplayName', 'theory + 0.3 dB (datasheet limit)');
o = T.RawBits > 0;
plot(ax, T.EbNo_channel_dB(o), T.RawBER(o), 'o', 'Color', [0 0.35 0.75], 'DisplayName', "raw BER, this sweep (" + cfg.receiver + ")");
if ~isempty(V.fig13)
    F = V.fig13;
    errBars(ax, F.EbNo_channel_dB, F.BER, F.BER_lo, F.BER_hi, 'd-', [0.85 0.3 0.1], 'uncoded, realistic receiver');
end
ylim(ax, [1e-6 0.5]); xlim(ax, [-4 11]);
xlabel(ax, 'E_b/N_0 per channel bit [dB]'); ylabel(ax, 'demodulator (uncoded) BER');
legend(ax, 'Location', 'southwest'); title(ax, 'Demodulator loss (datasheet Figure 13: < 0.3 dB)');
files(end+1, 1) = saveFig(f1, opt, "fig_ber_per.png");
% 2. PER vs Table 7
if ~isempty(V.table7)
    f2 = figure('Name', 'Table 7', 'Color', 'w', 'Position', [80 80 900 480]);
    ax = axes(f2); hold(ax, 'on'); grid(ax, 'on'); set(ax, 'YScale', 'log');
    errBars(ax, T.EbNo_dB(ok) + cfg.payFactor, T.PER(ok), T.PER_lo(ok), T.PER_hi(ok), 's-', [0.85 0.3 0.1], 'simulated PER');
    T7 = V.table7(V.table7.Modulation == cfg.mod & V.table7.Coding == cfg.coding, :);
    [g, xg] = findgroups(round(T7.EbNo_table_dB, 1));       % rates with the same Eb/N0 share a line
    cl = lines(numel(xg) + 1);
    for i = 1:numel(xg)
        r = T7(g == i, :);
        xline(ax, mean(r.EbNo_table_dB), '--', 'LineWidth', 1.5, 'Color', cl(i + 1, :), 'DisplayName', ...
            sprintf('Table 7: %s kBd (%s dBm)', strjoin(compose('%g', r.Rate_kBd), '/'), strjoin(compose('%.1f', r.P_table_dBm), '/')));
    end
    xlim(ax, [min([T.EbNo_dB(ok) + cfg.payFactor; T7.EbNo_table_dB]) - 0.5, max([T.EbNo_dB(ok) + cfg.payFactor; T7.EbNo_table_dB]) + 0.5]);
    yline(ax, 1e-2, ':', 'PER 1 %', 'HandleVisibility', 'off');
    xlabel(ax, sprintf('E_b/N_0 per payload bit [dB] (Table 7 converted with NF %.2f dB)', opt.NoiseFigure));
    ylabel(ax, 'packet error rate'); ylim(ax, [1e-4 1]);
    title(ax, sprintf('%s %s: simulation vs SRS-4 Table 7', cfg.mod, cfg.codingName));
    legend(ax, 'Location', 'southwest');
    files(end+1, 1) = saveFig(f2, opt, "fig_table7.png");
end
% 3. constellation, eye, spectrum at the plot point
ip = plotIndex(S, opt.PlotEbNo);
sn = [];
if ~isempty(ip), sn = S.snaps{ip}; end
if isempty(sn) && ~isempty(ip)                       % all points came from the cache: one fresh batch
    O = rxObjects(cfg, imp, S.run, 7);
    [~, sn] = runBatch(cfg, imp, T.EbNo_dB(ip) + cfg.esFactor, S.run, RandStream('mlfg6331_64', 'Seed', 7), O, true);
end
% clean reference: one batch at Es/N0 20 dB through the same receiver (impairments still applied)
esHigh = 20;
O = rxObjects(cfg, imp, S.run, 11);
[~, sh] = runBatch(cfg, imp, esHigh, S.run, RandStream('mlfg6331_64', 'Seed', 11), O, true);
if ~isempty(sn)
    xp = T.EbNo_dB(ip);
    ref = pskmod((0:cfg.M-1).', cfg.M, cfg.phOff, cfg.order);
    f3 = figure('Name', 'Constellation', 'Color', 'w', 'Position', [80 80 1000 500]);
    tl = tiledlayout(f3, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s symbols after the %s receiver (impairments "%s")', cfg.mod, cfg.receiver, imp.preset));
    pts = {sn.sym, sh.sym}; ttl2 = {sprintf('operating point: E_b/N_0 %.2f dB, E_s/N_0 %.2f dB', xp, xp + cfg.esFactor), ...
        sprintf('reference: E_s/N_0 %d dB (shows ISI, phase noise, IQ imbalance, loop jitter)', esHigh)};
    for i = 1:2
        ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
        plot(ax, real(pts{i}), imag(pts{i}), '.', 'MarkerSize', 4, 'Color', [0 0.35 0.75]);
        plot(ax, real(ref), imag(ref), 'r+', 'MarkerSize', 14, 'LineWidth', 2);
        axis(ax, 'equal'); xlim(ax, [-1.8 1.8]); ylim(ax, [-1.8 1.8]);
        xlabel(ax, 'in-phase'); ylabel(ax, 'quadrature'); title(ax, ttl2{i}, 'FontWeight', 'normal');
    end
    files(end+1, 1) = saveFig(f3, opt, "fig_constellation.png");
    up = 4;                                              % smoother traces for the eye (display only)
    ze = resample(sh.eye, up, 1); ze = ze(8*up*cfg.sps + 1:end - 8*up*cfg.sps);
    he = eyediagram(ze, 2*cfg.sps*up, 2, 0);   % sample 1 = symbol instant -> plotted at t = 0 (verified)
    set(he, 'Name', 'Eye diagram');
    axs = findobj(he, 'Type', 'axes');
    if ~isempty(axs)
        title(axs(end), sprintf('Eye after the matched RRC filter, E_s/N_0 %d dB (in-phase)', esHigh));
        if numel(axs) > 1, title(axs(1), 'quadrature'); end
    end
    files(end+1, 1) = saveFig(he, opt, "fig_eye.png");
    f4 = figure('Name', 'Spectrum', 'Color', 'w', 'Position', [100 100 900 450]);
    ax = axes(f4); hold(ax, 'on'); grid(ax, 'on');
    plot(ax, cfg.txPSD.f, cfg.txPSD.dB, 'LineWidth', 1.5, 'DisplayName', 'transmitted (pwelch)');
    plot(ax, sn.psd.f, sn.psd.dB - max(sn.psd.dB), 'DisplayName', sprintf('received, E_b/N_0 %.2f dB', xp));
    xline(ax, [-1 1]*cfg.obw99/2, '--', 'HandleVisibility', 'off');
    xline(ax, [-1 1]*(1 + cfg.rolloff)/2, ':', 'HandleVisibility', 'off');
    xlim(ax, [-1.5 1.5]); ylim(ax, [-70 5]);
    xlabel(ax, 'frequency offset / R_s'); ylabel(ax, 'PSD [dB, peak = 0]');
    title(ax, sprintf('OBW99 = %.3f R_s (datasheet 1.07 R_s at \\alpha 0.2), OBW99.9 = %.3f R_s (1.18)', cfg.obw99, cfg.obw999));
    legend(ax, 'Location', 'south');
    files(end+1, 1) = saveFig(f4, opt, "fig_spectrum.png");
    if opt.Scopes && interactive && strcmpi(opt.Visible, 'on')
        try
            cd = comm.ConstellationDiagram('ReferenceConstellation', ref, 'Title', 'SRS-4 link simulator'); cd(sn.sym);
            sa = spectrumAnalyzer('SampleRate', cfg.fs, 'Title', 'Received (simulated)'); sa(sn.rx);
        catch e
            warning('ff_ccsds_linksim:scopes', 'Scopes not available: %s', e.message);
        end
    end
end
% 4. pass
if ~isempty(PS)
    f5 = figure('Name', 'Pass', 'Color', 'w', 'Position', [120 60 1000 750]);
    tl = tiledlayout(f5, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Pass prediction: %s %s at %s (payload %.0f kbit/s)', cfg.mod, cfg.coding, rateText(cfg.Rs), cfg.RbPay/1e3));
    a = gobjects(4, 1); for i = 1:4, a(i) = nexttile(tl); hold(a(i), 'on'); grid(a(i), 'on'); end
    for k = 1:numel(PS.series)
        s = PS.series(k); m = minutes(s.t - s.t(1)); nm = sprintf('pass %d', k);
        plot(a(1), m, s.el, 'DisplayName', nm); plot(a(2), m, s.ebnoPayload, 'DisplayName', nm);
        plot(a(3), m, s.PER, 'DisplayName', nm); plot(a(4), m, s.throughput/1e3, 'DisplayName', nm);
    end
    yline(a(2), summ.EbNo_PER1_payload, '--', 'PER 1 %', 'HandleVisibility', 'off');
    set(a(3), 'YScale', 'log'); ylim(a(3), [1e-4 1]);
    fl = 0.5/max(T.Frames);
    yline(a(3), fl, ':', 'simulation floor 0.5/frames (no errors seen)', 'HandleVisibility', 'off', 'FontSize', 8);
    ylabel(a(1), 'elevation [deg]'); ylabel(a(2), 'E_b/N_0 payload [dB]'); ylabel(a(3), 'PER'); ylabel(a(4), 'throughput [kbit/s]');
    xlabel(a(4), 'minutes after AOS'); legend(a(1), 'Location', 'best');
    files(end+1, 1) = saveFig(f5, opt, "fig_pass.png");
end
% 5. measured
if ~isempty(MS)
    f6 = figure('Name', 'Measured', 'Color', 'w', 'Position', [140 80 900 480]);
    ax = axes(f6); hold(ax, 'on'); grid(ax, 'on'); set(ax, 'YScale', 'log');
    errBars(ax, T.EbNo_dB(ok), T.PER(ok), T.PER_lo(ok), T.PER_hi(ok), 's-', [0.85 0.3 0.1], 'simulated PER');
    errBars(ax, T.EbNo_dB(ok), T.BER(ok), T.BER_lo(ok), T.BER_hi(ok), 'o-', [0 0.35 0.75], 'simulated BER');
    xl = [min(T.EbNo_dB) - 1, max(T.EbNo_dB) + 1];
    cl = [0.2 0.5 0.2; 0.5 0.2 0.6];
    for i = 1:height(MS.table)
        xv = MS.table.EbNo_info_dB(i);
        lab = sprintf('%s: E_b/N_0 %.1f dB, margin %.1f dB to PER 1 %%', MS.table.Receiver(i), xv, MS.table.MarginToPER1_dB(i));
        if xv >= xl(1) && xv <= xl(2)
            xline(ax, xv, '-', 'LineWidth', 1.5, 'Color', cl(i, :), 'DisplayName', lab);
        else                                             % off scale: say so in the legend
            plot(ax, NaN, NaN, '-', 'LineWidth', 1.5, 'Color', cl(i, :), 'DisplayName', lab + " (off scale)");
        end
    end
    xlim(ax, xl); ylim(ax, [1e-7 1]);
    xlabel(ax, 'E_b/N_0 per information bit [dB]'); ylabel(ax, 'error rate'); legend(ax, 'Location', 'southwest');
    title(ax, sprintf('Measured C/N0 -> predicted BER/PER (%s %s, R_s %s)', cfg.mod, cfg.coding, rateText(MS.Rs)));
    files(end+1, 1) = saveFig(f6, opt, "fig_measured.png");
end
end

function errBars(ax, x, y, lo, hi, style, color, name)
x = x(:); y = y(:); lo = lo(:); hi = hi(:);
ok = isfinite(y) & y > 0;
h = plot(ax, x(ok), y(ok), style, 'Color', color, 'MarkerFaceColor', color, 'DisplayName', name);
z = isfinite(hi) & ~(y > 0);                          % zero-error points: show the upper bound
if any(z), plot(ax, x(z), hi(z), 'v', 'Color', color, 'HandleVisibility', 'off'); end
okb = ok & isfinite(lo) & isfinite(hi) & lo > 0;
if any(okb)
    errorbar(ax, x(okb), y(okb), y(okb) - lo(okb), hi(okb) - y(okb), 'LineStyle', 'none', 'Color', color, ...
        'HandleVisibility', 'off');
end
if isempty(h), plot(ax, NaN, NaN, style, 'Color', color, 'DisplayName', name); end
end

function ip = plotIndex(S, want)
T = S.T; ok = find(T.Frames > 0 & cellfun(@(s) ~isempty(s), S.snaps));
if isempty(ok), ok = find(T.Frames > 0); end
ip = [];
if isempty(ok), return, end
if ~isempty(want)
    [~, j] = min(abs(T.EbNo_dB(ok) - want)); ip = ok(j); return
end
[~, j] = min(abs(log10(max(T.PER(ok), 1e-6)) + 1)); ip = ok(j);
end

function f = saveFig(h, opt, name)
% PNG via exportgraphics (no drawnow: profiled at ~4 s per figure in -batch); hidden figures are closed
f = "";
if opt.Save
    f = string(fullfile(opt.OutDir, name));
    try
        exportgraphics(h, f, 'Resolution', 130);
    catch e
        warning('ff_ccsds_linksim:figure', 'Could not save %s: %s', name, e.message); f = "";
    end
end
if strcmpi(opt.Visible, 'off'), close(h); end
end

%% ===== text ======================================================================

function printSummary(R, cfg)
s = R.summary;
fprintf('\n=== ff_ccsds_linksim summary: %s %s, %s receiver ===\n', cfg.mod, cfg.codingName, cfg.receiver);
fprintf('  Eb/N0 at PER 1 %%   : %s dB (info), %s dB (payload)\n', ciText(s.EbNo_PER1, s.EbNo_PER1_CI), ...
    ciText(s.EbNo_PER1_payload, s.EbNo_PER1_CI + cfg.payFactor));
fprintf('  Eb/N0 at BER 1e-5  : %s dB; coding gain vs uncoded %s %.2f dB; Shannon limit %.2f dB (gap %.2f dB)\n', ...
    ciText(s.EbNo_BER5, s.EbNo_BER5_CI), cfg.mod, s.codingGain_BER5, s.shannonEbNo, s.shannonGap_BER5);
fprintf('  required C/N0      : %.2f dB-Hz at %s (receiver NF %.2f dB -> sensitivity %.1f dBm)\n', s.reqCN0_PER1_dBHz, ...
    rateText(cfg.Rs), s.NF, s.sensitivity_PER1_dBm);
fprintf('  demodulator loss   : %.2f dB from the raw BER (%d points, channel Eb/N0 %.1f to %.1f dB)\n', ...
    s.rawLoss_dB, s.rawLoss_points, s.rawLoss_range);
if isfinite(s.syncLoss_PER1), fprintf('  synchronisation loss at PER 1 %%: %.2f dB (realistic - ideal)\n', s.syncLoss_PER1); end
fprintf('  OBW99 %.3f Rs, OBW99.9 %.3f Rs (datasheet 1.07 / 1.18 Rs at alpha 0.2)\n', s.obw99_Rs, s.obw999_Rs);
if ~isempty(R.validation.fig13)
    disp(R.validation.fig13(:, {'EbNo_channel_dB', 'Errors', 'BER', 'BER_theory', 'Loss_dB', 'Loss_lo_dB', 'Loss_hi_dB', 'Verdict_vs_0p3dB'}));
end
if ~isempty(R.validation.table7)
    disp(R.validation.table7(:, {'Modulation', 'Coding', 'Rate_kBd', 'P_table_dBm', 'EbNo_table_dB', 'EbNo_sim_dB', 'Difference_dB'}));
end
if ~isempty(R.pass), disp(R.pass.table); end
if ~isempty(R.measured), disp(R.measured.table); end
for i = 1:numel(R.notes), fprintf('  - %s\n', R.notes(i)); end
if strlength(R.outDir) > 0, fprintf('  files in %s\n', R.outDir); end
end

function s = ciText(x, ci)
if ~isfinite(x), s = "n/a"; return, end
s = sprintf('%.2f [%.2f, %.2f]', x, ci(1), ci(2));
end

function s = rateText(r)
if r >= 1e6, s = sprintf('%g MBd', r/1e6); else, s = sprintf('%g kBd', r/1e3); end
end

function s = rmfieldIf(s, f)
if isfield(s, f), s = rmfield(s, f); end
end
