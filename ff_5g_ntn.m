function R = ff_5g_ntn(varargin)
%FF_5G_NTN  5G NR Non-Terrestrial Network (NTN) downlink from a LEO satellite: pass geometry,
%           link-level BLER, link budget, throughput over a pass, Doppler/delay vs NR numerology,
%           and coexistence with a measured SRS-4 (FieldFox) spectrum.
%
%   R = ff_5g_ntn                                   600 km circular LEO, overhead pass over Bangkok,
%                                                   band n256 (S-band), Quick link-level preset
%                                                   (about 1 min of simulation, 2-3 min in total)
%   R = ff_5g_ntn(Name, Value, ...)
%   R = ff_5g_ntn('TLE', 67683)                     real orbit: KNACKSAT-2 from CelesTrak (cached 2 h)
%   R = ff_5g_ntn('TLE', "my.tle", 'Pass', "first") a TLE file; first pass in 'SearchHours'
%   R = ff_5g_ntn('Interference', "demo")           + coexistence with a synthetic SRS-4 spectrum
%   R = ff_5g_ntn('Interference', "srs4_run.mat")   + coexistence with a MEASURED spectrum: an
%                                                   ff_srs4_analyzer run (file or its R) or any
%                                                   recording in the suite format (f, P, settings)
%   R = ff_5g_ntn('File', "...\ntn5g_run.mat", ...) reuse the BLER curves of an earlier run (no
%                                                   link-level simulation): fast what-if budgets
%
% Example (headless, exam figures and PDF, no dialogs):
%   R = ff_5g_ntn('MaxElevation', 60, 'Interference', "demo", 'Visible', 'off', 'Prompt', false);
%   R.challenges                   % Doppler / delay vs NR numerology, and the Rel-17 NTN answers
%   R.throughput                   % expected throughput and data volume over the pass
%   R2 = ff_5g_ntn('File', R, 'EIRPDensity', 28, 'Visible', 'off', 'Report', "none");  % what-if
%
% No instrument is used. The only network access is the optional TLE download from celestrak.org
% (gp.php?CATNR=..&FORMAT=TLE, cached 2 h in support\tle, not retried for 2 h after an error).
%
% What it computes
%  1. Geometry (Satellite Communications Toolbox): satelliteScenario + satellite + groundStation,
%     aer, dopplershift (Doppler and Doppler rate), latency, states; a circular orbit is placed
%     so that it passes over the site with the requested 'MaxElevation' (lla2eci + reckon), a TLE
%     uses SGP4 and access/accessIntervals to find the pass. Per time step: elevation, slant
%     range, one-way delay (timing advance = 2 x delay), Doppler, Doppler rate, FSPL (fspl),
%     off-nadir (scan) angle, P.618 atmospheric loss (p618PropagationLosses, if the maps are in
%     support\p618; otherwise a fixed 0.5 dB with a note).
%  2. Link budget per time step for a HANDHELD UE (TR 38.821 assumptions, see below):
%        SNR = EIRP density - scan loss - FSPL - atmosphere - polarisation - shadowing
%              - scintillation + G/T + 228.6        (dB, per resource element = per Hz)
%     cross-checked against satelliteCNR(satelliteCNRConfig) at three elevations.
%  3. Link level (5G Toolbox): slot-by-slot PDSCH, 1 layer: nrDLSCH (LDPC) -> nrPDSCH + DM-RS
%     (3 symbols, for Doppler) -> nrOFDMModulate -> nrTDLChannel with an NTN delay profile
%     (TR 38.811 NTN-TDL-A..D, NTN-TDLA100, NTN-TDLC5) and SatelliteDopplerShift = the CFO the UE
%     sees (or nrCDLChannel NTN-CDL-A..D with SatelliteElevationAngle) -> AWGN ->
%     CP-correlation fine CFO estimate (+-SCS/2 range, averaged over slots) ->
%     nrTimingEstimate -> nrOFDMDemodulate -> nrChannelEstimate -> nrEqualizeMMSE ->
%     nrPDSCHDecode -> nrDLSCHDecoder. BLER and throughput vs SNR per MCS (nrPDSCHMCSTables,
%     TS 38.214 Table 5.1.3.1-1, TBS from nrTBS), SCS 15 and 30 kHz, two cases:
%        "UE pre-comp": the UE removes the predicted Doppler (GNSS position + SIB19 ephemeris);
%                       residual = |Doppler rate| x 'PrecompInterval'/2 + fc x 'PrecompVelError'/c
%        "none":        the full Doppler at 'LinkElevation' reaches the receiver
%     CFO sweep: BLER vs carrier frequency offset with and without the fine CFO estimator at
%     'CFOSweepSNR' -> the CFO tolerance of each numerology in Hz. Expect about 0.05 SCS without
%     the estimator (the DM-RS channel estimate follows the phase ramp only while it turns less
%     than ~90 deg between DM-RS symbols up to 5 symbols apart: f < 1/(4 x 5 x Tsym)) and just
%     under SCS/2 with it (the CP correlation is ambiguous beyond +-SCS/2: an integer-SCS offset
%     needs the SSB search). Both scale with the SCS: the numerology choice for LEO Doppler.
%  4. Throughput over the pass: ideal link adaptation (best MCS for the SNR of each time step)
%     using the simulated BLER curves, scaled to 'ChannelBW' with nrTBS; Shannon bound for
%     reference; data volume per pass. Without pre-compensation: outage while |Doppler| exceeds
%     the simulated CFO tolerance.
%  5. Doppler/delay challenges vs NR numerology (R.numerology, R.challenges): Doppler / SCS,
%     CFO tolerance, Doppler drift per slot, RTT vs the NR timing-advance range, differential
%     delay in the beam vs CP, delay drift, HARQ processes, and the 3GPP Rel-17 NTN solutions.
%  6. Coexistence ('Interference'): the interferer's power that falls into the n256 victim
%     channel (TX leakage, integrated from the measured spectrum minus the analyzer noise floor,
%     which is the median PSD more than 3 Rs from the interferer's carrier) plus the power in the
%     adjacent channels divided by the UE adjacent-channel selectivity ('ACS'), scaled by the
%     geometry/coupling options -> I/N, SINR and throughput loss over the pass, I/N and volume
%     loss vs separation, and the free-space separation that brings I/N down to 'INR'. Parts
%     of the channel outside the measured span are filled with the SFCG/ECSS spectral mask (an
%     upper bound for a mask-compliant transmitter) and flagged.
% Results go to 'OutDir': ntn5g_report.pdf, figures (PNG), ntn5g_pass.csv, ntn5g_bler.csv and
% ntn5g_run.mat (variable results = R; re-use it with 'File').
%
% Name-Value options (defaults):
%  Orbit and site
%   'Altitude'        600e3   m, circular orbit (TR 38.821 LEO-600); ignored with 'TLE'
%   'Inclination'     53      deg (circular orbit); must exceed |site latitude| + offset
%   'MaxElevation'    90      deg, peak elevation of the constructed pass (circular orbit)
%   'TLE'             ""      NORAD number (CelesTrak, cached 2 h) or a TLE file; e.g. 67683
%                             (KNACKSAT-2, a real Thai CubeSat orbit, used only for geometry)
%   'Start'           now     UTC datetime: pass construction time / start of the TLE search
%   'SearchHours'     24      h searched for passes (TLE)
%   'Pass'            "best"  "best" (highest elevation), "first", or the pass number (TLE)
%   'Site'            [13.7563 100.5018 10]   lat (deg), lon (deg), altitude (m): Bangkok
%   'SiteName'        "Bangkok"
%   'MinElevation'    10      deg, elevation mask (TR 38.821 LEO minimum)
%   'TimeStep'        5       s
%  Band and link budget
%   'Band'            "n256"  "n256" (S-band, DL 2170-2200 MHz) or "n255" (L-band, DL 1525-1559)
%   'CarrierFreq'     []      Hz; [] = top channel of the band (2197.5 MHz for n256 / 5 MHz):
%                             the channel next to the SRS-4 band 2200-2290 MHz (worst case)
%   'ChannelBW'       5e6     Hz: 5, 10, 15 or 20 MHz (TS 38.101-5 NTN bands)
%   'SCS'             [15 30] kHz subcarrier spacings to compare (FR1 NTN: 15 and 30)
%   'EIRPDensity'     []      dBW/MHz; [] = TR 38.821 Table 6.1.1.1-1 Set-1 S-band:
%                             34 (LEO-600, altitude <= 900 km), 40 (LEO-1200)
%   'SatGain'         30      dBi, satellite max antenna gain (TR 38.821, S-band LEO)
%   'BeamDiameter'    []      km; [] = 50 (LEO-600) / 90 (LEO-1200), for the differential delay
%   'ScanLossExp'     1       scan loss = -10 n log10(cos(off-nadir)): earth-fixed beam steered
%                             by a planar array (n = 1: projected aperture); 0 = none
%   'UEGain'          0       dBi   } handheld UE, TR 38.821 Table 6.1.1.1-3
%   'UENF'            7       dB    }
%   'UETemp'          290     K antenna temperature
%   'PolLoss'         3       dB (circular satellite -> linear handheld)
%   'ShadowMargin'    3       dB
%   'ScintLoss'       2.2     dB (ionospheric; Bangkok is near the equatorial anomaly, where
%                             post-sunset scintillation at L/S-band can be much worse)
%   'AtmLoss'         "p618"  tropospheric loss: "p618" (ITU-R P.618 maps) or a fixed dB value
%   'P618Exceedance'  1       % of an average year
%  Link level
%   'Preset'          "quick" "quick": MCS [0 4 9 13], SNR -8:2:14 dB, 5 MHz, <= 30 slots/point
%                             "full":  MCS 0:2:16, SNR -8:1:16 dB, up to 10 MHz, <= 100 slots
%                             (Full takes 10-30 min; use 'UseParallel', true)
%   'MCS' 'SNR' 'SimBW' 'MaxSlots' 'MaxErrors'   override the preset ('SNR' in dB, ascending)
%   'Channel'         "NTN-TDL-C"  NTN-TDL-A/B (NLOS), NTN-TDL-C/D (LOS), NTN-TDLA100, NTN-TDLC5
%                             (nrTDLChannel), or NTN-CDL-A..D (nrCDLChannel, SatelliteElevationAngle =
%                             the LinkElevation; the CFO is then applied to the waveform because
%                             nrCDLChannel has no SatelliteDopplerShift; co-polarised isotropic
%                             antennas (RX: a half-wave ULA), so the elevation changes only the
%                             angles, which matter with RxAntennas > 1; the 2 RX elements are
%                             correlated by the LOS, so CDL gains less diversity than TDL)
%   'DelaySpread'     100e-9  s (NTN-TDL/CDL-A..D only; TR 38.811 values depend on scenario/elevation)
%   'UESpeed'         3       km/h (fading Doppler; the satellite Doppler is separate)
%   'RxAntennas'      1       UE receive antennas
%   'LinkElevation'   30      deg: the Doppler and Doppler rate of the simulated case
%   'PrecompInterval' 0.1     s between UE Doppler pre-compensation updates
%   'PrecompVelError' 1       m/s, satellite velocity error of the broadcast ephemeris
%   'FineCFO'         true    CP-correlation fine CFO estimator in the receiver
%   'CFOSweep'        true    BLER vs CFO (SCS 15/30, fine estimator off/on)
%   'CFOSweepSNR'     10      dB    'CFOSweepMCS'  4
%   'Seed'            1
%   'UseParallel'     "auto"  parfor over the simulation jobs: "auto" = only if a pool is
%                             already open, true = open one (up to 8 workers), N = open one with
%                             N workers, false = serial. Results are identical (seeded per job).
%                             Measured, Quick: 69 s serial, 39 s on 4 workers (+20 s pool start):
%                             worth it for Full, hardly for Quick
%  Coexistence
%   'Interference'    []      "demo", a MAT file, an ff_srs4_analyzer result, or struct(f, P)
%                             (from an ff_srs4_analyzer run: the sweep that covers the victim
%                             channel, preferring the average detector and RBW <= ChannelBW/50,
%                             i.e. its mask sweep when fc is near 2200 MHz, else its wide scan)
%   'InterfCoupling'  0       dB added to the measured level (extra path, body loss, ...)
%   'InterfMeasDistance'   [] m, interferer -> measuring antenna (from the run if it has it)
%   'InterfVictimDistance' [] m, interferer -> 5G UE; with both: 20 log10(meas/victim) dB
%   'InterfMeasGain'  0       dBi, measuring antenna (the UE has 'UEGain')
%   'ACS'             33      dB, UE adjacent-channel selectivity (TS 38.101-1 Table 7.5-1 value
%                             for <= 20 MHz channels; check TS 38.101-5 for NTN)
%   'INR'             -6      dB, protection criterion I/N (1 dB desensitisation)
%  Output
%   'Report'          "pdf"   "pdf", "docx" or "none"
%   'Notes'           ""      text for the report
%   'OutDir'          ""      default ntn5g_<time> in the current folder
%   'Visible'         'on'    'off' = no windows (figures are still saved)
%   'Prompt'          true    ask before downloading a TLE and before a long Full run
%                             (only in an interactive MATLAB desktop)
%   'File'            ""      earlier ntn5g_run.mat or R: reuse its link-level curves
%
% Output R: options, orbit, pass (summary), timeline (table per time step), budget, linkLevel,
% cfoSweep, simInfo, throughput, numerology, challenges, interference (summary, throughput,
% vsDistance, ...), conclusions, assumptions, notes, files, outDir, runtime_s.
%
% Read this before trusting the numbers:
%   - A MODEL, not a measurement. The satellite parameters are 3GPP study values (TR 38.821
%     Set-1, S-band); a real NTN satellite (and certainly a CubeSat TLE used for geometry) has
%     its own EIRP, beam layout and payload type. Losses are planning margins.
%   - Earth-fixed beam centred on the UE (no beam-edge loss) and a regenerative payload (gNB on
%     board): the RTT is the service link only. A transparent payload adds the feeder link
%     (TR 38.821 quotes about 26 ms worst-case RTT for transparent LEO-600).
%   - Link level: 1 layer, no HARQ retransmissions (NTN may disable HARQ feedback), PDSCH in
%     symbols 2-13 (2-symbol CORESET reserved), no SSB/CSI-RS overhead, practical channel and
%     timing estimation, SNR per resource element and receive antenna. The BLER curves come from
%     'SimBW' and are applied to 'ChannelBW' (LDPC block size effects ignored). A few tens of
%     slots per point (Quick) give BLER to about +-0.1: use Full for smooth curves.
%   - The satellite Doppler over the simulated slots is held constant (the drift over 30 slots
%     is < 20 Hz); the fading drops are 1 s apart (independent for UESpeed >= 1 km/h).
%   - "none" models a UE without any Doppler pre-compensation or SSB-based integer-CFO search:
%     a real UE acquires the DL Doppler during cell search and then tracks it; the hard NTN
%     problems are UL frequency/timing (many UEs, different Doppler) and initial access. The
%     simulation shows WHY a +-SCS/2 fine estimator alone cannot cope with LEO Doppler.
%   - Coexistence: the FieldFox has no I/Q; only the power spectrum is used. Leakage below the
%     analyzer noise floor cannot be measured: R.interference gives the noise-subtracted
%     estimate AND the upper bound (with the floor). Over-the-air levels depend on antenna
%     gains, polarisation and multipath (+-3 dB at best). A measurement taken with a pad and
%     ExtGain is DUT-referred (transmitter output): then 'InterfCoupling' must hold the TX
%     antenna gain - path loss (+ UE gain - measuring gain). The UE noise density of a handheld
%     (about -167 dBm/Hz) is usually far below the analyzer floor: the run says when the upper
%     bound cannot prove the I/N criterion. Only the victim and the first adjacent channels are
%     counted (no blocking, intermodulation or second-adjacent terms).
%   - Free-space distance scaling is valid in the far field and without obstructions.
%
% 5G Toolbox NTN support in R2026a (checked in MATLAB): nrTDLChannel DelayProfile NTN-TDL-A..D,
% NTN-TDLA100, NTN-TDLC5 and property SatelliteDopplerShift; nrCDLChannel NTN-CDL-A..D and
% SatelliteElevationAngle (no SatelliteDopplerShift). There is no nrNTN* function or NTN
% satellite/timing-advance object, so geometry comes from the Satellite Communications Toolbox
% and the NTN timing/HARQ figures are computed here from TS 38.211/38.213 formulas.
%
% Sources: 3GPP TR 38.821 (NTN solutions: parameter sets, link budgets), TR 38.811 (NTN channel
% models), TS 38.101-5 (NTN UE: bands n255/n256), TS 38.211/38.213/38.214 (numerology, TA,
% MCS/TBS), ITU-R P.618-14, ECSS-E-ST-50-05C / SFCG Rec. 21-2 (spectral mask). Values marked
% "TR 38.821" are typical study values, quoted from the TR; check them against the current
% release before using them in an exam answer.
%
% See also ff_srs4_analyzer, ff_cubesat_link, fieldfox_live, nrTDLChannel, satelliteScenario.

tStart = tic;
opt = parseOptions(varargin);
here = fileparts(mfilename('fullpath'));
if strlength(string(opt.OutDir)) == 0, opt.OutDir = fullfile(pwd, "ntn5g_" + stamp()); end
opt.OutDir = string(opt.OutDir);
if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
notes = strings(0, 1);
c0 = 299792458;

%% 1. Geometry
fprintf('1/6 Pass geometry (%s) ...\n', orbitLabel(opt));
[T, orb, notes] = passGeometry(opt, here, notes);
orb = satelliteParameters(orb, opt);
[T, atm, notes] = atmosphere(T, opt, here, notes);

%% 2. Link budget
fprintf('2/6 Link budget per time step (%s) ...\n', orb.eirpSource);
[T, B] = linkBudget(T, opt, orb, atm);
vis = T.Visible;
fprintf('   pass %s - %s (UTC+7), max elevation %.1f deg, SNR %.1f .. %.1f dB\n', ...
    string(T.Time(find(vis, 1)), 'HH:mm:ss', 'en_US'), string(T.Time(find(vis, 1, 'last')), 'HH:mm:ss'), ...
    max(T.Elevation_deg), min(T.SNR_dB(vis)), max(T.SNR_dB(vis)));

%% 3. Link level
cond = linkConditions(T, opt, c0);
notes = [notes; cond.notes];
[LL, CF, simInfo, notes, opt] = linkLevelOrReuse(opt, cond, notes);

%% 4. Throughput over the pass
fprintf('4/6 Throughput over the pass, 5/6 Doppler/delay vs numerology ...\n');
[T, TP, tol] = passThroughput(T, LL, CF, opt);
lowSnr = T.Visible & T.SNR_dB < min(LL.SNR_dB);
if any(lowSnr)
    notes(end + 1, 1) = sprintf(['%.0f %% of the pass has an SNR below the simulated range (%g dB): counted as ' ...
        'zero throughput (conservative). Extend ''SNR'' downwards for low passes or small EIRP.'], ...
        100*nnz(lowSnr)/nnz(T.Visible), min(LL.SNR_dB));
end

%% 5. Doppler / delay vs numerology
[NUM, CH] = challengeTables(T, opt, orb, cond, tol, c0);

%% 6. Coexistence
IF = struct('done', false);
if ~isempty(opt.Interference) && ~(isstring(opt.Interference) && all(strlength(opt.Interference) == 0))
    fprintf('6/6 Coexistence with the measured spectrum ...\n');
    try
        [IF, T] = coexistence(opt, T, LL, B, orb);
    catch err
        warning('ff_5g_ntn:interference', 'Coexistence analysis failed: %s (%s)', err.message, where(err));
        IF = struct('done', false, 'error', string(err.message));
        notes(end + 1, 1) = "Coexistence analysis failed: " + err.message;
    end
end

%% Results
R = struct();
R.options = opt; R.orbit = orb; R.pass = passSummary(T, orb, opt);
R.timeline = T; R.budget = B.table; R.linkLevel = LL; R.cfoSweep = CF; R.simInfo = simInfo;
R.throughput = TP; R.numerology = NUM; R.challenges = CH; R.interference = IF;
R.assumptions = assumptionList(opt, orb, atm, cond);
R.notes = [notes; atm.notes];
R.conclusions = conclusions(R, opt, cond, tol);
R.outDir = opt.OutDir;
fprintf('\n'); fprintf('  - %s\n', R.conclusions); fprintf('\n');

files = strings(0, 1);
figs0 = findall(groot, 'Type', 'figure');
try
    files = makeFigures(R, opt, cond, tol);
catch err
    warning('ff_5g_ntn:figures', 'Figures failed: %s (%s)', err.message, where(err));
    if strcmpi(opt.Visible, 'off')          % the half-drawn hidden figure would stay in the session
        f = findall(groot, 'Type', 'figure'); delete(f(~ismember(f, figs0)));
    end
end
csv1 = fullfile(opt.OutDir, "ntn5g_pass.csv"); writetable(T(T.Visible, :), csv1);
csv2 = fullfile(opt.OutDir, "ntn5g_bler.csv"); writetable(LL, csv2);
R.files = [files; csv1; csv2];
R.runtime_s = toc(tStart);
matFile = fullfile(opt.OutDir, "ntn5g_run.mat");
results = R; results.options.Interference = summarizeInput(opt.Interference);
if isstruct(results.options.File), results.options.File = "struct"; end
save(matFile, 'results');
R.files = [R.files; matFile];
if opt.Report ~= "none"
    try
        R.files = [makeReport(R, opt); R.files];
    catch err
        warning('ff_5g_ntn:report', 'Report failed (%s). Figures and data are in %s.', err.message, opt.OutDir);
    end
end
R.runtime_s = toc(tStart);
fprintf('Done in %.0f s. Results in %s\n', R.runtime_s, opt.OutDir);
end

%% ===== options ==================================================================

function opt = parseOptions(args)
opt = struct('Altitude', 600e3, 'Inclination', 53, 'MaxElevation', 90, 'TLE', "", 'Start', [], ...
    'SearchHours', 24, 'Pass', "best", 'Site', [13.7563 100.5018 10], 'SiteName', "Bangkok", ...
    'MinElevation', 10, 'TimeStep', 5, ...
    'Band', "n256", 'CarrierFreq', [], 'ChannelBW', 5e6, 'SCS', [15 30], 'EIRPDensity', [], ...
    'SatGain', 30, 'BeamDiameter', [], 'ScanLossExp', 1, 'UEGain', 0, 'UENF', 7, 'UETemp', 290, ...
    'PolLoss', 3, 'ShadowMargin', 3, 'ScintLoss', 2.2, 'AtmLoss', "p618", 'P618Exceedance', 1, ...
    'Preset', "quick", 'MCS', [], 'SNR', [], 'SimBW', [], 'MaxSlots', [], 'MaxErrors', [], ...
    'Channel', "NTN-TDL-C", 'DelaySpread', 100e-9, 'UESpeed', 3, 'RxAntennas', 1, ...
    'LinkElevation', 30, 'PrecompInterval', 0.1, 'PrecompVelError', 1, 'FineCFO', true, ...
    'CFOSweep', true, 'CFOSweepSNR', 10, 'CFOSweepMCS', 4, 'Seed', 1, 'UseParallel', "auto", ...
    'Interference', [], 'InterfCoupling', 0, 'InterfMeasDistance', [], 'InterfVictimDistance', [], ...
    'InterfMeasGain', 0, 'ACS', 33, 'INR', -6, ...
    'File', "", 'Report', "pdf", 'Notes', "", 'OutDir', "", 'Visible', 'on', 'Prompt', true);
if mod(numel(args), 2), error('ff_5g_ntn:args', 'Options must be Name, Value pairs.'); end
names = fieldnames(opt);
for k = 1:2:numel(args)
    j = find(strcmpi(names, args{k}), 1);
    if isempty(j), error('ff_5g_ntn:args', 'Unknown option "%s".', string(args{k})); end
    opt.(names{j}) = args{k + 1};
end
opt.Visible = char(matlab.lang.OnOffSwitchState(opt.Visible));
opt.Prompt = logical(opt.Prompt);
opt.Band = lower(string(opt.Band)); opt.Preset = lower(string(opt.Preset)); opt.Report = lower(string(opt.Report));
opt.Channel = upper(string(opt.Channel)); opt.Pass = string(opt.Pass);
opt.SiteName = string(opt.SiteName); opt.Notes = string(opt.Notes);
if ~ismember(opt.Band, ["n256" "n255"]), error('ff_5g_ntn:args', 'Band must be "n256" or "n255".'); end
if ~ismember(opt.Preset, ["quick" "full"]), error('ff_5g_ntn:args', 'Preset must be "quick" or "full".'); end
if ~ismember(opt.Report, ["pdf" "docx" "none"]), error('ff_5g_ntn:args', 'Report must be "pdf", "docx" or "none".'); end
chans = ["NTN-TDL-A" "NTN-TDL-B" "NTN-TDL-C" "NTN-TDL-D" "NTN-TDLA100" "NTN-TDLC5" ...
    "NTN-CDL-A" "NTN-CDL-B" "NTN-CDL-C" "NTN-CDL-D"];
if ~ismember(opt.Channel, chans), error('ff_5g_ntn:args', 'Channel must be one of %s.', strjoin(chans, ', ')); end
if opt.ChannelBW < 1000, opt.ChannelBW = opt.ChannelBW*1e6; end           % MHz given
if ~ismember(opt.ChannelBW, [5 10 15 20]*1e6), error('ff_5g_ntn:args', 'ChannelBW must be 5, 10, 15 or 20 MHz.'); end
opt.SCS = unique(opt.SCS(:).');
if isempty(opt.SCS) || any(~ismember(opt.SCS, [15 30])), error('ff_5g_ntn:args', 'SCS must be 15 and/or 30 (kHz).'); end
if numel(opt.Site) < 2, error('ff_5g_ntn:args', 'Site must be [lat lon] or [lat lon altitude_m].'); end
if numel(opt.Site) < 3, opt.Site(3) = 0; end
if opt.MinElevation < 0 || opt.MinElevation >= 80, error('ff_5g_ntn:args', 'MinElevation must be 0-80 deg.'); end
if opt.MaxElevation <= opt.MinElevation || opt.MaxElevation > 90
    error('ff_5g_ntn:args', 'MaxElevation must lie between MinElevation and 90 deg.');
end
if opt.RxAntennas < 1 || opt.RxAntennas ~= round(opt.RxAntennas), error('ff_5g_ntn:args', 'RxAntennas must be a positive integer.'); end
DL = bandEdges(opt.Band);
if isempty(opt.CarrierFreq), opt.CarrierFreq = DL(2) - opt.ChannelBW/2; end
if opt.CarrierFreq - opt.ChannelBW/2 < DL(1) - 1 || opt.CarrierFreq + opt.ChannelBW/2 > DL(2) + 1
    warning('ff_5g_ntn:band', 'The %g MHz channel at %.2f MHz is not inside the %s DL band %g-%g MHz.', ...
        opt.ChannelBW/1e6, opt.CarrierFreq/1e6, opt.Band, DL/1e6);
end
% TLE: NORAD number or file
opt.NORAD = []; opt.TLEFile = "";
tl = opt.TLE;
if isnumeric(tl) && ~isempty(tl), opt.NORAD = tl(1);
elseif strlength(string(tl)) > 0
    if all(isstrprop(char(strtrim(string(tl))), 'digit')), opt.NORAD = str2double(tl);
    else, opt.TLEFile = string(tl); end
end
if isempty(opt.Start), opt.Start = dateshift(datetime('now', 'TimeZone', 'UTC'), 'start', 'minute'); end
if ~isdatetime(opt.Start), opt.Start = datetime(opt.Start, 'TimeZone', 'UTC'); end
if isempty(opt.Start.TimeZone), opt.Start.TimeZone = 'UTC'; end
opt.Start.TimeZone = 'UTC';
% presets (fields the user did not set come from the preset)
opt = applyPreset(opt, opt.Preset, ["MCS" "SNR" "SimBW" "MaxSlots" "MaxErrors"]);
opt = checkLinkLevel(opt);
end

function opt = applyPreset(opt, preset, fields)
if preset == "quick"
    P = struct('MCS', [0 4 9 13], 'SNR', -8:2:14, 'SimBW', 5e6, 'MaxSlots', 30, 'MaxErrors', 6, ...
        'CFOPoints', [0 0.02 0.05 0.1 0.15 0.2 0.3 0.45 0.55]);
else
    P = struct('MCS', 0:2:16, 'SNR', -8:1:16, 'SimBW', min(opt.ChannelBW, 10e6), 'MaxSlots', 100, ...
        'MaxErrors', 20, 'CFOPoints', [0:0.025:0.3 0.35:0.05:0.6 0.8 1]);
end
if ~isfield(opt, 'PresetFields'), opt.PresetFields = strings(1, 0); end
for f = fields
    if isempty(opt.(f)) || ismember(f, opt.PresetFields)
        opt.(f) = P.(f); opt.PresetFields = union(opt.PresetFields, f);
    end
end
opt.CFOPoints = P.CFOPoints; opt.Preset = preset;
end

function opt = checkLinkLevel(opt)
if opt.SimBW < 1000, opt.SimBW = opt.SimBW*1e6; end
if ~ismember(opt.SimBW, [5 10 15 20]*1e6), error('ff_5g_ntn:args', 'SimBW must be 5, 10, 15 or 20 MHz.'); end
opt.SNR = sort(unique(opt.SNR(:).'));
opt.MCS = unique(opt.MCS(:).');
if any(opt.MCS < 0 | opt.MCS > 28), error('ff_5g_ntn:args', 'MCS indices must be 0-28 (TS 38.214 Table 5.1.3.1-1).'); end
if numel(opt.SNR) < 2, error('ff_5g_ntn:args', 'SNR needs at least two values.'); end
for f = ["MaxSlots" "MaxErrors"]
    v = opt.(f);
    if ~isscalar(v) || v < 1 || v ~= round(v), error('ff_5g_ntn:args', '%s must be a positive integer.', f); end
end
if ~isscalar(opt.CFOSweepMCS) || ~ismember(opt.CFOSweepMCS, 0:28), error('ff_5g_ntn:args', 'CFOSweepMCS must be 0-28.'); end
end

function DL = bandEdges(band)
% TS 38.101-5 Table 5.2.2-1 (FDD): n256 DL 2170-2200 / UL 1980-2010; n255 DL 1525-1559 / UL 1626.5-1660.5 MHz
if band == "n256", DL = [2170 2200]*1e6; else, DL = [1525 1559]*1e6; end
end

function n = nrbFor(bw, scs)
% TS 38.101-1 Table 5.3.2-1 transmission bandwidth configuration (NRB), FR1
B = [5 10 15 20]*1e6; N15 = [25 52 79 106]; N30 = [11 24 38 51];
i = find(B == bw, 1);
if scs == 15, n = N15(i); else, n = N30(i); end
end

function s = orbitLabel(opt)
if ~isempty(opt.NORAD), s = sprintf('TLE NORAD %d', opt.NORAD);
elseif strlength(opt.TLEFile) > 0, s = "TLE " + opt.TLEFile;
else, s = sprintf('circular %.0f km, i = %g deg, peak elevation %g deg', opt.Altitude/1e3, opt.Inclination, opt.MaxElevation);
end
end

%% ===== geometry ==================================================================

function [T, orb, notes] = passGeometry(opt, here, notes)
Re = 6378137; mu = 3.986004418e14; c0 = 299792458;
site = opt.Site(:).';
orb = struct('name', "", 'source', "", 'circular', false, 'altitude_km', NaN, 'inclination', NaN, ...
    'tle', strings(0, 1), 'passes', table());
if ~isempty(opt.NORAD) || strlength(opt.TLEFile) > 0
    [tleFile, orb.source, notes] = getTLE(opt, here, notes);
    orb.tle = readlines(tleFile); orb.tle = orb.tle(strlength(strtrim(orb.tle)) > 0);
    t0 = opt.Start;
    sc = satelliteScenario(t0, t0 + hours(opt.SearchHours), 30);
    sat = satellite(sc, tleFile, 'OrbitPropagator', 'sgp4');
    if numel(sat) > 1
        notes(end + 1, 1) = sprintf('TLE file has %d satellites: the first (%s) is used.', numel(sat), sat(1).Name);
        sat = sat(1);
    end
    gs = groundStation(sc, site(1), site(2), 'Altitude', site(3), 'MinElevationAngle', opt.MinElevation, 'Name', opt.SiteName);
    AI = accessIntervals(access(sat, gs));
    if height(AI) == 0
        error('ff_5g_ntn:nopass', 'No pass above %g deg in %g h from %s UTC. Increase SearchHours or lower MinElevation.', ...
            opt.MinElevation, opt.SearchHours, string(t0));
    end
    [~, elc, ~, tc] = aer(gs, sat);
    AI.MaxElevation = zeros(height(AI), 1);
    for k = 1:height(AI)
        m = tc >= AI.StartTime(k) & tc <= AI.EndTime(k);
        if any(m), AI.MaxElevation(k) = max(elc(m)); else, AI.MaxElevation(k) = opt.MinElevation; end
    end
    if opt.Pass == "best", [~, k] = max(AI.MaxElevation);
    elseif opt.Pass == "first", k = 1;
    else
        k = str2double(opt.Pass);
        if ~isfinite(k) || k < 1 || k > height(AI), error('ff_5g_ntn:args', 'Pass must be "best", "first" or 1-%d.', height(AI)); end
    end
    orb.passes = AI(:, intersect(["StartTime" "EndTime" "Duration" "MaxElevation"], AI.Properties.VariableNames, 'stable'));
    for v = intersect(["StartTime" "EndTime"], orb.passes.Properties.VariableNames)    % local time like the timeline
        x = orb.passes.(v); if isempty(x.TimeZone), x.TimeZone = 'UTC'; end
        x.TimeZone = 'Asia/Bangkok'; orb.passes.(v) = x;
    end
    orb.passIndex = k;
    sc2 = satelliteScenario(AI.StartTime(k) - seconds(60), AI.EndTime(k) + seconds(60), opt.TimeStep);
    sat2 = satellite(sc2, tleFile, 'OrbitPropagator', 'sgp4'); sat2 = sat2(1);
    orb.name = string(sat2.Name);
    oe = orbitalElements(sat2);
    orb.inclination = oe.Inclination; orb.period_min = oe.Period/60;
    try                                                                   % TLE age at the pass
        ep = oe.Epoch; if isempty(ep.TimeZone), ep.TimeZone = 'UTC'; end
        ts = AI.StartTime(k); if isempty(ts.TimeZone), ts.TimeZone = 'UTC'; end
        orb.tleEpoch = ep; ageD = days(ts - ep);
        if abs(ageD) > 3
            notes(end + 1, 1) = sprintf(['TLE epoch %s UTC is %.1f days from the pass: SGP4 errors grow by roughly ' ...
                'km per day in LEO, so pass times can be off by seconds to minutes (geometry statistics stay valid).'], ...
                string(ep, 'yyyy-MM-dd HH:mm'), ageD);
        end
    catch
    end
else
    h = opt.Altitude; a = Re + h; inc = opt.Inclination; n = sqrt(mu/a^3);
    psiMask = acosd(Re/a*cosd(opt.MinElevation)) - opt.MinElevation;     % Earth central angle at the mask
    half = ceil(deg2rad(psiMask)/n) + 90;                                 % s, mask crossing -> closest approach
    tmid = opt.Start + seconds(half);
    psi = acosd(Re/a*cosd(opt.MaxElevation)) - opt.MaxElevation;          % central angle at closest approach
    if abs(site(1)) + psi >= inc
        error('ff_5g_ntn:args', 'Inclination %g deg is too low for a pass over latitude %g deg (needs > %.1f).', ...
            inc, site(1), abs(site(1)) + psi);
    end
    beta = asind(min(1, cosd(inc)/cosd(site(1))));                       % ground-track heading (ascending)
    [plat, plon] = reckon(site(1), site(2), psi, beta + 90);              % sub-satellite point at closest approach
    u = asind(sind(plat)/sind(inc));                                      % argument of latitude there
    r = lla2eci([plat plon h], datevec(tmid));
    raan = mod(atan2d(r(2), r(1)) - atan2d(sind(u)*cosd(inc), cosd(u)), 360);
    nu0 = mod(u - rad2deg(n*half), 360);
    sc2 = satelliteScenario(opt.Start, opt.Start + seconds(2*half), opt.TimeStep);
    sat2 = satellite(sc2, a, 0, inc, raan, 0, nu0, 'OrbitPropagator', 'two-body-keplerian', ...
        'Name', sprintf('LEO-%.0f', h/1e3));
    orb.name = string(sat2.Name); orb.circular = true; orb.inclination = inc; orb.period_min = 2*pi/n/60;
    orb.source = sprintf(['Circular orbit %.0f km, i = %g deg, RAAN %.2f deg, placed to pass over %s ' ...
        'with ~%g deg peak elevation (two-body Keplerian)'], h/1e3, inc, raan, opt.SiteName, opt.MaxElevation);
end
gs2 = groundStation(sc2, site(1), site(2), 'Altitude', site(3), 'MinElevationAngle', opt.MinElevation, 'Name', opt.SiteName);
[az, el, rng, t] = aer(gs2, sat2);
[fD, ~, dI] = dopplershift(sat2, gs2, Frequency = opt.CarrierFreq);
lat = latency(sat2, gs2);
geo = states(sat2, 'CoordinateFrame', 'geographic');
t = t(:); t.TimeZone = 'Asia/Bangkok';
rate = dI.DopplerRate(:); N = numel(t);
if numel(rate) == N - 1                    % R2026a: finite difference between samples (N-1 values)
    rate = [rate(1); (rate(1:end - 1) + rate(2:end))/2; rate(end)];       % -> at the sample times
end
T = table(t, az(:), el(:), rng(:)/1e3, geo(3, :).'/1e3, geo(1, :).', geo(2, :).', rng(:)/c0*1e3, ...
    fD(:), rate, dI.RelativeVelocity(:), ...
    'VariableNames', ["Time" "Azimuth_deg" "Elevation_deg" "Range_km" "SatAlt_km" "SatLat_deg" ...
    "SatLon_deg" "Delay_ms" "Doppler_Hz" "DopplerRate_Hzps" "RelVelocity_mps"]);
T.Visible = T.Elevation_deg >= opt.MinElevation & isfinite(T.Doppler_Hz);
if nnz(T.Visible) < 3
    error('ff_5g_ntn:nopass', 'The pass has fewer than 3 time steps above %g deg: reduce TimeStep.', opt.MinElevation);
end
d = abs(lat(:) - rng(:)/c0); d = d(T.Visible & isfinite(d));
if ~isempty(d) && max(d) > 1e-6
    notes(end + 1, 1) = sprintf('latency() differs from range/c by up to %.1f us.', max(d)*1e6);
end
orb.altitude_km = mean(T.SatAlt_km(T.Visible));
end

function [file, src, notes] = getTLE(opt, here, notes)
if strlength(opt.TLEFile) > 0
    if ~isfile(opt.TLEFile), error('ff_5g_ntn:tle', 'TLE file not found: %s', opt.TLEFile); end
    file = opt.TLEFile; src = "TLE file " + file; return
end
norad = opt.NORAD;
if norad > 99999, error('ff_5g_ntn:tle', 'NORAD > 99999 cannot be a TLE (use an OMM file with satellite()).'); end
d = fullfile(here, "support", "tle");
if ~isfolder(d), mkdir(d); end
file = fullfile(d, sprintf('%d.tle', norad));
errFile = file + ".err";
age = fileAgeHours(file);
if isfile(file) && age < 2
    src = sprintf('CelesTrak NORAD %d, cached %.0f min ago (refreshed every 2 h)', norad, age*60); return
end
if isfile(errFile) && fileAgeHours(errFile) < 2
    notes(end + 1, 1) = "CelesTrak failed less than 2 h ago (" + strtrim(string(fileread(errFile))) + "): not retried.";
elseif ~askUser(opt, sprintf('Download the TLE of NORAD %d from celestrak.org? [y/n] ', norad))
    notes(end + 1, 1) = "TLE download declined.";
else
    url = sprintf('https://celestrak.org/NORAD/elements/gp.php?CATNR=%d&FORMAT=TLE', norad);
    try
        fprintf('   fetching %s\n', url);
        txt = webread(url, weboptions('ContentType', 'text', 'Timeout', 20));
        L = strtrim(splitlines(string(txt))); L = L(strlength(L) > 0);
        i1 = find(startsWith(L, "1 ") & strlength(L) >= 69, 1);
        if isempty(i1) || i1 + 1 > numel(L) || ~startsWith(L(i1 + 1), "2 ")
            error('ff_5g_ntn:tle', 'CelesTrak returned no TLE for NORAD %d (%s).', norad, extractBefore(string(txt) + newline, newline));
        end
        nm = "SAT " + norad; if i1 > 1, nm = L(i1 - 1); end
        writelines([nm; L(i1); L(i1 + 1)], file);
        if isfile(errFile), delete(errFile); end
        src = sprintf('CelesTrak NORAD %d (downloaded now)', norad); return
    catch err
        msg = string(err.message); msg = extractBefore(msg + newline, newline);
        writelines(string(datetime('now')) + " " + msg, errFile);     % CelesTrak etiquette: no retry for 2 h
        notes(end + 1, 1) = "CelesTrak download failed: " + msg;
    end
end
if isfile(file)
    src = sprintf('CelesTrak NORAD %d, cached %.1f h ago (no fresh download)', norad, fileAgeHours(file));
    return
end
error('ff_5g_ntn:tle', 'No TLE for NORAD %d (no cache, no download). Give a TLE file or use the circular orbit.', norad);
end

function h = fileAgeHours(file)
h = Inf; d = dir(file);
if ~isempty(d), h = hours(datetime('now') - datetime(d.datenum, 'ConvertFrom', 'datenum')); end
end

function ok = askUser(opt, question)
ok = true;
if opt.Prompt && usejava('desktop')
    a = input(question, 's');
    ok = isempty(a) || any(lower(a(1)) == 'yj');
end
end

function orb = satelliteParameters(orb, opt)
% TR 38.821 Table 6.1.1.1-1 (Set-1, S-band): LEO-600 34 dBW/MHz, 50 km beam; LEO-1200 40 dBW/MHz, 90 km
if orb.altitude_km <= 900, cls = "LEO-600"; e = 34; bd = 50; else, cls = "LEO-1200"; e = 40; bd = 90; end
orb.class = cls;
orb.eirpDensity = e; orb.eirpSource = sprintf('TR 38.821 Set-1 S-band %s value (%g dBW/MHz)', cls, e);
if ~isempty(opt.EIRPDensity), orb.eirpDensity = opt.EIRPDensity; orb.eirpSource = sprintf('user value %g dBW/MHz', opt.EIRPDensity); end
if opt.Band == "n255", orb.eirpSource = orb.eirpSource + " (S-band value reused for L-band: assumption)"; end
orb.beamDiameter_km = bd;
if ~isempty(opt.BeamDiameter), orb.beamDiameter_km = opt.BeamDiameter; end
end

function [T, atm, notes] = atmosphere(T, opt, here, notes)
atm = struct('source', "", 'grid', table(), 'notes', strings(0, 1));
el = T.Elevation_deg; vis = T.Visible;
if isnumeric(opt.AtmLoss)
    T.Atm_dB = repmat(opt.AtmLoss, height(T), 1);
    atm.source = sprintf('fixed %.2f dB (user)', opt.AtmLoss);
    return
end
p = fullfile(here, 'support', 'p618');
if isfolder(p), addpath(genpath(p)); end
try
    ws = warning('off', 'satcom:p618PropagationLosses:InvalidXPDFrequency');   % XPD is not used here
    restore = onCleanup(@() warning(ws));
    eg = unique([max(5, floor(min(el(vis)))) 10 15 20 30 45 60 90]);
    eg = eg(eg >= max(5, floor(min(el(vis)))));
    cfg = p618Config;
    cfg.Frequency = max(opt.CarrierFreq, 1e9); cfg.Latitude = opt.Site(1); cfg.Longitude = opt.Site(2);
    pe = opt.P618Exceedance;
    cfg.GasAnnualExceedance = pe; cfg.CloudAnnualExceedance = pe; cfg.RainAnnualExceedance = pe;
    cfg.ScintillationAnnualExceedance = pe; cfg.TotalAnnualExceedance = pe;
    cfg.AntennaDiameter = 0.1; cfg.AntennaEfficiency = 0.5;            % handheld: small aperture (scintillation)
    G = zeros(numel(eg), 5);
    for k = 1:numel(eg)
        cfg.ElevationAngle = eg(k);
        pl = p618PropagationLosses(cfg, 'StationHeight', max(opt.Site(3), 0)/1e3);
        G(k, :) = [pl.At pl.Ag pl.Ac pl.Ar pl.As];
    end
    atm.grid = array2table([eg(:) G], 'VariableNames', ["Elevation_deg" "Total_dB" "Gas_dB" "Cloud_dB" "Rain_dB" "Scint_dB"]);
    e = min(max(el, eg(1)), 90);
    T.Atm_dB = interp1(eg(:), G(:, 1), e, 'pchip');
    atm.source = sprintf('ITU-R P.618-14 (p618PropagationLosses, %g %% of the year): gas + cloud + rain + tropospheric scintillation', pe);
catch err
    T.Atm_dB = repmat(0.5, height(T), 1);
    atm.source = "fixed 0.5 dB (P.618 not available)";
    atm.notes(end + 1, 1) = "P.618 not available (" + string(err.message) + "): fixed 0.5 dB tropospheric loss used. " + ...
        "Install the ITU-R P.618 maps in support\p618 for elevation-dependent values.";
end
end

%% ===== link budget ================================================================

function [T, B] = linkBudget(T, opt, orb, atm)
Re = 6371e3; c0 = 299792458;
fc = opt.CarrierFreq;
Tsys = opt.UETemp + 290*(10^(opt.UENF/10) - 1);
GT = opt.UEGain - 10*log10(Tsys);
Lfix = opt.PolLoss + opt.ShadowMargin + opt.ScintLoss;
sinTh = (Re + opt.Site(3))./(Re + T.SatAlt_km*1e3).*cosd(T.Elevation_deg);
T.OffNadir_deg = asind(min(max(sinTh, 0), 1));
T.ScanLoss_dB = -10*opt.ScanLossExp*log10(max(cosd(T.OffNadir_deg), 1e-3));
T.FSPL_dB = fspl(T.Range_km*1e3, c0/fc);
T.SNR_dB = (orb.eirpDensity - 60) - T.ScanLoss_dB - T.FSPL_dB - T.Atm_dB - Lfix + GT + 228.6;
T.SNR_dB(~T.Visible) = NaN;
scs = opt.SCS(1); nrb = nrbFor(opt.ChannelBW, scs); Bocc = nrb*12*scs*1e3;
B = struct('Tsys', Tsys, 'GT', GT, 'Lfix', Lfix, 'Bocc', Bocc, 'nrb', nrb);
% budget table at the peak, at LinkElevation and at the mask (TR 38.821 style), with satelliteCNR
vis = find(T.Visible);
[~, iMax] = max(T.Elevation_deg);
iLink = vis(find(T.Elevation_deg(vis) >= opt.LinkElevation, 1));
if isempty(iLink), iLink = iMax; end
iMin = vis(1);
idx = [iMax iLink iMin];
terms = ["Elevation" "Slant range" "Satellite EIRP density" "EIRP in occupied BW" "Scan loss" ...
    "Free-space path loss" "Atmospheric loss (P.618)" "Polarisation loss" "Shadowing margin" ...
    "Scintillation (ionosphere)" "UE G/T" "Boltzmann constant" "C/N0" "Occupied bandwidth" ...
    "SNR (this tool)" "SNR (satelliteCNR)"];
units = ["deg" "km" "dBW/MHz" "dBW" "dB" "dB" "dB" "dB" "dB" "dB" "dB/K" "dBW/K/Hz" "dB-Hz" "MHz" "dB" "dB"];
V = zeros(numel(terms), 3);
for j = 1:3
    i = idx(j);
    eirp = orb.eirpDensity + 10*log10(Bocc/1e6);
    cn0 = eirp - T.ScanLoss_dB(i) - T.FSPL_dB(i) - T.Atm_dB(i) - Lfix + GT + 228.6;
    cfg = satelliteCNRConfig;
    cfg.TransmitterPower = eirp - opt.SatGain; cfg.TransmitterSystemLoss = 0;
    cfg.TransmitterAntennaGain = opt.SatGain - T.ScanLoss_dB(i);
    cfg.Distance = T.Range_km(i); cfg.Frequency = fc/1e9;
    cfg.MiscellaneousLoss = T.Atm_dB(i) + Lfix; cfg.GainToNoiseTemperatureRatio = GT;
    cfg.ReceiverSystemLoss = 0; cfg.Bandwidth = Bocc/1e6; cfg.BitRate = Bocc/1e6; cfg.SymbolRate = Bocc/1e6;
    cnr = satelliteCNR(cfg);
    V(:, j) = [T.Elevation_deg(i); T.Range_km(i); orb.eirpDensity; eirp; -T.ScanLoss_dB(i); -T.FSPL_dB(i); ...
        -T.Atm_dB(i); -opt.PolLoss; -opt.ShadowMargin; -opt.ScintLoss; GT; 228.6; cn0; Bocc/1e6; ...
        T.SNR_dB(i); cnr];
end
B.table = table(terms(:), units(:), V(:, 1), V(:, 2), V(:, 3), 'VariableNames', ...
    ["Term" "Unit" "PeakElevation" "LinkElevation" "MinElevation"]);
B.cnrCheck = max(abs(V(end, :) - V(end - 1, :)));
B.atmSource = atm.source;
end

%% ===== link level =================================================================

function cond = linkConditions(T, opt, c0)
vis = find(T.Visible);
i = vis(find(T.Elevation_deg(vis) >= opt.LinkElevation, 1));      % rising side: approaching, Doppler > 0
cond.notes = strings(0, 1);
low = isempty(i);
if low                                                              % never reaches LinkElevation:
    [~, j] = max(abs(T.Doppler_Hz(vis))); i = vis(j);               % worst Doppler (near the mask)
    cond.notes(end + 1, 1) = sprintf(['Peak elevation %.1f deg < LinkElevation %g deg: the link-level Doppler is the ' ...
        'largest of the pass (%.1f deg elevation) and the pre-compensation residual uses the largest Doppler rate.'], ...
        max(T.Elevation_deg), opt.LinkElevation, T.Elevation_deg(i));
end
cond.index = i; cond.el = T.Elevation_deg(i);
cond.fD = abs(T.Doppler_Hz(i)); cond.rate = abs(T.DopplerRate_Hzps(i));
if low, cond.rate = max(abs(T.DopplerRate_Hzps(vis))); end
cond.resid = cond.rate*opt.PrecompInterval/2 + opt.CarrierFreq*opt.PrecompVelError/c0;
cond.fdUE = opt.UESpeed/3.6*opt.CarrierFreq/c0;
cond.maxDoppler = max(abs(T.Doppler_Hz(vis)));
cond.maxRate = max(abs(T.DopplerRate_Hzps(vis)));
cond.residMax = cond.maxRate*opt.PrecompInterval/2 + opt.CarrierFreq*opt.PrecompVelError/c0;
end

function [LL, CF, info, notes, opt] = linkLevelOrReuse(opt, cond, notes)
reuse = isstruct(opt.File) || strlength(string(opt.File)) > 0;
if reuse
    S = opt.File;
    if ~isstruct(S)
        if ~isfile(S), error('ff_5g_ntn:file', 'File not found: %s', S); end
        S = load(S);
        if isfield(S, 'results'), S = S.results; end
    end
    if ~isfield(S, 'linkLevel') || ~istable(S.linkLevel)
        error('ff_5g_ntn:file', '''File'' must be an ff_5g_ntn result (ntn5g_run.mat or R).');
    end
    LL = S.linkLevel; CF = table(); if isfield(S, 'cfoSweep'), CF = S.cfoSweep; end
    info = struct('reused', true, 'seconds', 0, 'parallel', false, 'jobs', 0);
    if isfield(S, 'simInfo'), info = S.simInfo; info.reused = true; end
    have = unique(LL.SCS_kHz).';
    if ~all(ismember(opt.SCS, have))
        notes(end + 1, 1) = "Reused curves only have SCS " + strjoin(string(have), "/") + " kHz: SCS set to that.";
    end
    opt.SCS = have;
    if isfield(S, 'options') && isfield(S.options, 'SimBW'), opt.SimBW = S.options.SimBW; end
    notes(end + 1, 1) = "Link-level curves reused from an earlier run (no new simulation).";
    fprintf('3/6 Link level: reusing %d BLER points from the earlier run.\n', height(LL));
    return
end
if opt.Preset == "full" && ~askUser(opt, 'The Full link-level preset takes 10-30 min. Run it? (n = Quick) [y/n] ')
    opt = checkLinkLevel(applyPreset(opt, "quick", ["MCS" "SNR" "SimBW" "MaxSlots" "MaxErrors"]));
    notes(end + 1, 1) = "Full preset declined at the prompt: Quick preset used.";
end
[LL, CF, info] = linkLevel(opt, cond);
end

function [LL, CF, info] = linkLevel(opt, cond)
mt = nrPDSCHMCSTables; mt = mt.QAM64Table;
jobs = [];
for s = 1:numel(opt.SCS)
    scs = opt.SCS(s); nrb = nrbFor(opt.SimBW, scs); seed = opt.Seed*100 + s;
    base = struct('kind', "snr", 'scs', scs, 'nrb', nrb, 'mcs', 0, 'mod', "QPSK", 'tcr', 0.1, ...
        'points', opt.SNR, 'cfo', 0, 'snr', NaN, 'fine', logical(opt.FineCFO), 'precomp', "UE pre-comp", ...
        'seed', seed, 'maxSlots', opt.MaxSlots, 'maxErrors', opt.MaxErrors, 'slotsPerDrop', 5, ...
        'stopAtZero', true, 'stopAtOne', false, 'channel', opt.Channel, 'ds', opt.DelaySpread, ...
        'fdUE', cond.fdUE, 'nRx', opt.RxAntennas, 'fc', opt.CarrierFreq, 'el', cond.el);
    for m = opt.MCS
        J = base; r = mt(mt.MCSIndex == m, :);
        J.mcs = m; J.mod = string(r.Modulation); J.tcr = r.TargetCodeRate; J.cfo = cond.resid;
        jobs = [jobs J]; %#ok<AGROW>
    end
    J = base; r = mt(mt.MCSIndex == min(opt.MCS), :);                     % no pre-compensation: most robust MCS only
    J.mcs = min(opt.MCS); J.mod = string(r.Modulation); J.tcr = r.TargetCodeRate;
    J.cfo = cond.fD; J.precomp = "none"; J.stopAtZero = true;
    jobs = [jobs J]; %#ok<AGROW>
    if opt.CFOSweep
        r = mt(mt.MCSIndex == opt.CFOSweepMCS, :);
        for fine = [false true]
            J = base; J.kind = "cfo"; J.mcs = opt.CFOSweepMCS; J.mod = string(r.Modulation); J.tcr = r.TargetCodeRate;
            J.points = opt.CFOPoints*scs*1e3; J.snr = opt.CFOSweepSNR; J.fine = fine; J.precomp = "sweep";
            J.stopAtZero = false; J.stopAtOne = ~fine;
            jobs = [jobs J]; %#ok<AGROW>
        end
    end
end
nJ = numel(jobs);
M = 0;
up = opt.UseParallel; nW = min(8, feature('numcores'));
if isnumeric(up) && isscalar(up) && up > 1, nW = min(round(up), feature('numcores')); up = true; end
if (islogical(up) || isnumeric(up)) && up, up = "true"; elseif (islogical(up) || isnumeric(up)), up = "false"; end
up = lower(string(up));
if up ~= "false" && nJ >= 4 && license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
    try
        pool = gcp('nocreate');
        if isempty(pool) && up == "true", pool = parpool('Processes', nW); end
        if ~isempty(pool), M = pool.NumWorkers; end
    catch err
        warning('ff_5g_ntn:parallel', 'No parallel pool (%s): running serially.', err.message);
    end
end
fprintf('3/6 Link level: %d jobs (PDSCH, %s, %g MHz), %s ...\n', nJ, opt.Channel, opt.SimBW/1e6, ...
    ternary(M > 0, sprintf('parfor on %d workers', M), 'serial'));
out = cell(1, nJ);
t0 = tic;
if M > 0
    parfor (k = 1:nJ, M)
        out{k} = runJob(jobs(k));
    end
else
    for k = 1:nJ                      % plain loop: in order (parfor with 0 workers runs backwards)
        out{k} = runJob(jobs(k));
    end
end
info = struct('reused', false, 'seconds', toc(t0), 'parallel', M > 0, 'workers', M, 'jobs', nJ);
% assemble tables
LL = table(); CF = table();
for k = 1:nJ
    J = jobs(k); o = out{k}; n = numel(J.points);
    if J.kind == "snr"
        slotDur = 1e-3*15/J.scs;
        t = table(repmat(J.scs, n, 1), repmat(J.precomp, n, 1), repmat(J.mcs, n, 1), repmat(J.mod, n, 1), ...
            repmat(J.tcr, n, 1), J.points(:), o.bler(:), o.slots(:), o.errors(:), o.simulated(:), ...
            repmat(o.tbs, n, 1), (1 - o.bler(:))*o.tbs/slotDur/1e6, repmat(J.cfo, n, 1), ...
            'VariableNames', ["SCS_kHz" "Precomp" "MCS" "Modulation" "CodeRate" "SNR_dB" "BLER" "Slots" ...
            "Errors" "Simulated" "TBS_bits" "Thr_Mbps" "CFO_Hz"]);
        LL = [LL; t]; %#ok<AGROW>
    else
        t = table(repmat(J.scs, n, 1), repmat(J.fine, n, 1), J.points(:), J.points(:)/(J.scs*1e3), o.bler(:), ...
            o.slots(:), o.errors(:), o.simulated(:), repmat(J.snr, n, 1), repmat(J.mcs, n, 1), ...
            'VariableNames', ["SCS_kHz" "FineCFO" "CFO_Hz" "CFO_norm" "BLER" "Slots" "Errors" "Simulated" "SNR_dB" "MCS"]);
        CF = [CF; t]; %#ok<AGROW>
    end
end
fprintf('   link level done in %.0f s\n', info.seconds);
end

function o = runJob(J)
% One PDSCH link-level job: an SNR sweep (fixed CFO) or a CFO sweep (fixed SNR).
tj = tic;
[car, pd] = pdschSetup(J.scs, J.nrb, J.mod);
[~, pinfo] = nrPDSCHIndices(car, pd);
tbs = nrTBS(char(J.mod), 1, J.nrb, pinfo.NREPerPRB, J.tcr);
enc = nrDLSCH('MultipleHARQProcesses', false, 'TargetCodeRate', J.tcr);
dec = nrDLSCHDecoder('MultipleHARQProcesses', false, 'TargetCodeRate', J.tcr, 'TransportBlockLength', tbs, ...
    'LDPCDecodingAlgorithm', 'Normalized min-sum', 'MaximumLDPCIterationCount', 8);
oi = nrOFDMInfo(car);
isCDL = startsWith(J.channel, "NTN-CDL");
if isCDL                                                   % NTN-CDL: 1 TX, nRx isotropic RX (ULA), LOS at J.el
    tdl = nrCDLChannel('DelayProfile', char(J.channel), 'DelaySpread', J.ds, 'CarrierFrequency', J.fc, ...
        'SampleRate', oi.SampleRate, 'MaximumDopplerShift', J.fdUE, 'Seed', J.seed, ...
        'SatelliteElevationAngle', min(max(J.el, 0), 90));
    % co-polarised isotropic elements: the default arrays (TX 45 deg slant, RX 0 deg) lose 3 dB of
    % polarisation mismatch, which the link budget already holds as PolLoss (measured: 0.0 dB mean gain)
    tdl.TransmitAntennaArray.Size = [1 1 1 1 1]; tdl.TransmitAntennaArray.Element = 'isotropic';
    tdl.TransmitAntennaArray.PolarizationAngles = [0 90];
    tdl.ReceiveAntennaArray.Size = [1 J.nRx 1 1 1]; tdl.ReceiveAntennaArray.Element = 'isotropic';
    tdl.ReceiveAntennaArray.PolarizationAngles = [0 90];
else
    tdl = nrTDLChannel('DelayProfile', char(J.channel), 'SampleRate', oi.SampleRate, 'NumTransmitAntennas', 1, ...
        'NumReceiveAntennas', J.nRx, 'MaximumDopplerShift', J.fdUE, 'Seed', J.seed);
    if ismember(J.channel, ["NTN-TDL-A" "NTN-TDL-B" "NTN-TDL-C" "NTN-TDL-D"]), tdl.DelaySpread = J.ds; end
end
ci = info(tdl); dly = ci.ChannelFilterDelay; pad = ceil(ci.MaximumChannelDelay) + 16;
nSmp = 0;                                                  % sample counter for the CDL frequency offset
rs = RandStream('mt19937ar', 'Seed', J.seed);
np = numel(J.points);
o = struct('bler', NaN(1, np), 'errors', zeros(1, np), 'slots', zeros(1, np), 'simulated', false(1, np), 'tbs', tbs);
ones1 = 0;
for p = 1:np
    if J.kind == "snr", snrDb = J.points(p); cfo = J.cfo; else, snrDb = J.snr; cfo = J.points(p); end
    reset(rs);                                                             % common random numbers per point
    N0 = 1/sqrt(2*J.nRx*double(oi.Nfft)*10^(snrDb/10));                    % SNR per RE and receive antenna
    acc = 0; nErr = 0; s = 0;
    while s < J.maxSlots
        if mod(s, J.slotsPerDrop) == 0                                     % new fading drop, 1 s later
            release(tdl); tdl.InitialTime = s/J.slotsPerDrop;
            if ~isCDL, tdl.SatelliteDopplerShift = cfo; end
        end
        car.NSlot = mod(s, car.SlotsPerFrame);
        [ind, pinfo] = nrPDSCHIndices(car, pd);
        dmi = nrPDSCHDMRSIndices(car, pd); dms = nrPDSCHDMRS(car, pd);
        setTransportBlock(enc, randi(rs, [0 1], tbs, 1));
        cw = enc(char(pd.Modulation), 1, pinfo.G, 0);
        grid = nrResourceGrid(car, 1);
        grid(ind) = nrPDSCH(car, pd, cw); grid(dmi) = dms;
        tx = nrOFDMModulate(car, grid);
        rx = tdl([tx; zeros(pad, 1)]);
        if isCDL                                           % CFO as a continuous phase ramp
            rx = rx.*exp(1j*2*pi*cfo*(nSmp + (0:size(rx, 1) - 1).')/oi.SampleRate);
            nSmp = nSmp + size(tx, 1);
        end
        rx = rx + N0*complex(randn(rs, size(rx)), randn(rs, size(rx)));
        if J.fine
            [rx, acc] = cpFineCFO(rx, oi, car, dly, J.scs*1e3, acc);
        end
        off = nrTimingEstimate(car, rx, dmi, dms);
        if off > pad, off = dly; end                                       % implausible peak: keep the nominal delay
        rxg = nrOFDMDemodulate(car, rx(1 + off:end, :));
        [H, nVar] = nrChannelEstimate(car, rxg, dmi, dms, 'CDMLengths', pd.DMRS.CDMLengths);
        [r, h] = nrExtractResources(ind, rxg, H);
        eq = nrEqualizeMMSE(r, h, nVar);
        llr = nrPDSCHDecode(car, pd, eq, nVar);
        reset(dec);
        [~, err] = dec(llr, char(pd.Modulation), 1, 0);
        nErr = nErr + double(err); s = s + 1;
        if nErr >= J.maxErrors, break, end
    end
    o.errors(p) = nErr; o.slots(p) = s; o.bler(p) = nErr/s; o.simulated(p) = true;
    if J.stopAtZero && nErr == 0                                            % monotonic: higher SNR -> BLER 0
        o.bler(p + 1:end) = 0; break
    end
    if J.stopAtOne
        if o.bler(p) == 1, ones1 = ones1 + 1; else, ones1 = 0; end
        if ones1 >= 2, o.bler(p + 1:end) = 1; break, end                   % ICI only grows with the CFO
    end
end
if J.kind == "snr"
    fprintf('   SCS %2d kHz, MCS %2d (%s %.2f), %-11s CFO %7.0f Hz: BLER %s  (%.0f s)\n', J.scs, J.mcs, J.mod, ...
        J.tcr, J.precomp, J.cfo, strjoin(compose('%.2f', o.bler(o.simulated)), ' '), toc(tj));
else
    fprintf('   SCS %2d kHz, CFO sweep at %g dB, fine CFO %-3s: BLER %s  (%.0f s)\n', J.scs, J.snr, ...
        string(matlab.lang.OnOffSwitchState(J.fine)), strjoin(compose('%.2f', o.bler(o.simulated)), ' '), toc(tj));
end
end

function [rx, acc] = cpFineCFO(rx, oi, car, dly, scsHz, acc)
% CP correlation (van de Beek): r(n) r*(n+Nfft) over the second half of every CP (free of ISI and
% of the transmit windowing). Range +-SCS/2; the correlation is accumulated over the slots of the
% point (an AFC that averages), so the estimate improves with time like a real UE tracking loop.
N = double(oi.Nfft);
spsf = car.SlotsPerSubframe;
k0 = mod(car.NSlot, spsf)*14;
L = double(oi.SymbolLengths(k0 + (1:14))); cp = double(oi.CyclicPrefixLengths(k0 + (1:14)));
st = [0 cumsum(L(1:end - 1))];
a = 0;
for k = 1:14
    q = ceil(cp(k)/2);
    i0 = st(k) + dly + (cp(k) - q + 1:cp(k));
    i0 = i0(i0 >= 1 & i0 + N <= size(rx, 1));
    a = a + sum(rx(i0, :).*conj(rx(i0 + N, :)), 'all');
end
acc = acc + a;
f = -angle(acc)/(2*pi)*scsHz;
n = (0:size(rx, 1) - 1).';
rx = rx.*exp(-1j*2*pi*f*n/oi.SampleRate);
end

%% ===== throughput over the pass ======================================================

function [T, TP, tol] = passThroughput(T, LL, CF, opt)
vis = T.Visible; dt = seconds(median(diff(T.Time)));
TP = table(); tol = struct('scs', [], 'fineOff', [], 'fineOn', [], 'source', strings(0, 1));
for scs = opt.SCS
    nrbF = nrbFor(opt.ChannelBW, scs);
    curves = curvesFor(LL, scs, opt);
    [thr, mcs] = adaptThroughput(T.SNR_dB, curves);
    thr(~vis) = 0; mcs(~vis) = NaN;
    % tolerance from the CFO sweep (fine estimator on) for the "no pre-compensation" case
    [tOff, tOn, src] = cfoTolerance(CF, scs, opt);
    tol.scs(end + 1) = scs; tol.fineOff(end + 1) = tOff; tol.fineOn(end + 1) = tOn; tol.source(end + 1) = src;
    tUse = tOn; if ~opt.FineCFO, tUse = tOff; end
    ok = abs(T.Doppler_Hz) <= tUse;
    thrNo = thr.*ok;
    Bocc = nrbF*12*scs*1e3;
    shannon = Bocc*log2(1 + 10.^(T.SNR_dB/10))/1e6; shannon(~vis) = 0;
    tag = sprintf('SCS%d', scs);
    T.("MCS_" + tag) = mcs; T.("Thr_" + tag + "_Mbps") = thr; T.("ThrNoPrecomp_" + tag + "_Mbps") = thrNo;
    if scs == opt.SCS(1), T.Shannon_Mbps = shannon; end
    row = table(scs, nrbF, max([curves.rate]), mean(thr(vis)), max(thr(vis)), sum(thr(vis))*dt/8, ...
        mean(thrNo(vis)), sum(thrNo(vis))*dt/8, 100*mean(~ok(vis)), mean(shannon(vis)), tUse, ...
        'VariableNames', ["SCS_kHz" "NRB" "PeakRate_Mbps" "MeanThr_Mbps" "MaxThr_Mbps" "Volume_MB" ...
        "MeanThrNoPrecomp_Mbps" "VolumeNoPrecomp_MB" "OutageNoPrecomp_pct" "MeanShannon_Mbps" "CFOTolerance_Hz"]);
    TP = [TP; row]; %#ok<AGROW>
end
end

function [thr, mcs] = adaptThroughput(snr, curves)
% Ideal link adaptation: per time step the MCS with the highest (1 - BLER) x rate.
snr = snr(:); Tm = zeros(numel(snr), numel(curves));
for m = 1:numel(curves)
    Tm(:, m) = (1 - blerAt(curves(m).x, curves(m).y, snr))*curves(m).rate;
end
[thr, j] = max(Tm, [], 2);
mcs = reshape([curves(j).mcs], [], 1); mcs = double(mcs); mcs(thr <= 0 | isnan(snr)) = NaN;
thr(isnan(snr)) = 0;
end

function b = blerAt(x, y, s)
% log-linear interpolation of a BLER curve; below the simulated range = failure (conservative)
x = x(:); y = min(max(y(:), 1e-4), 1);
[x, i] = unique(x); y = y(i);
b = ones(size(s));
in = s >= x(1);
if isscalar(x), b(in) = y; else, b(in) = 10.^interp1(x, log10(y), min(s(in), x(end)), 'linear'); end
b(b <= 1.0001e-4) = 0;
end

function [tOff, tOn, src] = cfoTolerance(CF, scs, opt)
% Largest CFO with BLER <= 0.1 at CFOSweepSNR (linear interpolation between sweep points).
tOff = NaN; tOn = NaN; src = "assumed";
if ~isempty(CF) && any(CF.SCS_kHz == scs)
    for fine = [false true]
        C = sortrows(CF(CF.SCS_kHz == scs & CF.FineCFO == fine, :), 'CFO_Hz');
        if isempty(C), continue, end
        x = C.CFO_Hz; y = C.BLER; i = find(y > 0.1, 1);
        if isempty(i), t = x(end);
        elseif i == 1, t = 0;
        else, t = x(i - 1) + (0.1 - y(i - 1))/(y(i) - y(i - 1))*(x(i) - x(i - 1));
        end
        if fine, tOn = t; else, tOff = t; end
    end
    src = sprintf("simulated (BLER <= 0.1, MCS %d, %g dB)", opt.CFOSweepMCS, opt.CFOSweepSNR);
end
if isnan(tOn), tOn = 0.4*scs*1e3; src = "assumed 0.4 SCS (no CFO sweep)"; end
if isnan(tOff), tOff = 0.1*scs*1e3; end
end

%% ===== numerology / challenges =========================================================

function [NUM, CH] = challengeTables(T, opt, orb, cond, tol, c0)
vis = T.Visible;
fc = opt.CarrierFreq;
dMax = max(T.Delay_ms(vis)); dMin = min(T.Delay_ms(vis)); rtt = 2*dMax;
drift = max(abs(T.RelVelocity_mps(vis)))/c0;                            % s/s
elMin = min(T.Elevation_deg(vis));
diffDelay = orb.beamDiameter_km*1e3*cosd(elMin)/c0;                     % s, across the beam at the mask
NUM = table();
for scs = opt.SCS
    mu = log2(scs/15);
    cpUs = 144*2^-mu/(15e3*2048)*1e6;                                    % normal CP (not the long first one)
    k = find(tol.scs == scs, 1);
    row = table(scs, 1e3/scs, cpUs, 1/2^mu, cond.maxDoppler/1e3, cond.maxDoppler/(scs*1e3), ...
        0.1*scs, tol.fineOff(k)/1e3, tol.fineOn(k)/1e3, tol.fineOn(k)/fc*1e6, cond.residMax, ...
        cond.maxRate*1e-3/2^mu, 2.003/2^mu, rtt, diffDelay*1e6, cpUs*1e-6/drift, ceil(rtt/(1/2^mu)) + 2, ...
        'VariableNames', ["SCS_kHz" "Symbol_us" "CP_us" "Slot_ms" "MaxDoppler_kHz" "Doppler_over_SCS" ...
        "CFOtol_rule_kHz" "CFOtol_noFine_kHz" "CFOtol_fine_kHz" "CFOtol_fine_ppm" "ResidualPrecomp_Hz" ...
        "DriftPerSlot_Hz" "MaxTA_ms" "RTT_ms" "DiffDelayBeam_us" "CPconsumed_s" "HARQneeded"]);
    NUM = [NUM; row]; %#ok<AGROW>
end
ppm = cond.maxDoppler/fc*1e6;
s1 = opt.SCS(1);
ch = { ...
 "Doppler (common to the beam)", sprintf('+-%.1f kHz = %.1f ppm = %.2f x SCS%d', cond.maxDoppler/1e3, ppm, cond.maxDoppler/(s1*1e3), s1), ...
    "UE oscillator +-0.1 ppm after DL sync; CP-based fine CFO covers only +-SCS/2; ~0.1 SCS keeps ICI below ~-15 dB", ...
    "UE pre-compensates with its GNSS position + satellite ephemeris (SIB19), mainly for UL; DL: network may pre-compensate at beam centre, UE tracks with SSB/TRS";
 "Doppler rate", sprintf('up to %.0f Hz/s = %.2f Hz per 1 ms slot, %.0f Hz per 100 ms', cond.maxRate, cond.maxRate*1e-3, cond.maxRate*0.1), ...
    "fixed base station: ~0 (UE motion only)", ...
    "continuous pre-compensation update from the ephemeris; residual after updates every " + opt.PrecompInterval + " s: " + sprintf('%.0f Hz', cond.residMax);
 "Propagation delay / RTT", sprintf('%.2f-%.2f ms one way, RTT up to %.1f ms (service link only)', dMin, dMax, rtt), ...
    sprintf('max timing advance %.2f ms (SCS15) / %.2f ms (SCS30)', 2.003, 1.0015), ...
    "UE-specific TA from GNSS + ephemeris, common TA broadcast (SIB19), K_offset scheduling offset, extended RAR window";
 "Differential delay in the beam", sprintf('%.0f us across a %.0f km beam at %.0f deg elevation', diffDelay*1e6, orb.beamDiameter_km, elMin), ...
    "CP 4.7 us (SCS15) / 2.3 us (SCS30); PRACH formats sized for the cell radius", ...
    "UE applies its own TA before PRACH, so the gNB sees only the residual";
 "Delay drift", sprintf('up to %.1f us/s: the SCS15 CP is used up in %.2f s', drift*1e6, 4.69e-6/drift), ...
    "~0 (static UE-BS distance)", "UE updates its TA autonomously (common TA drift and drift rate in SIB19)";
 "HARQ round trip", sprintf('RTT %.1f ms -> about %d HARQ processes at 1 ms slots', rtt, ceil(rtt) + 2), ...
    "16 HARQ processes", "up to 32 HARQ processes; HARQ feedback can be disabled per process (then RLC ARQ / blind repetition)";
 "Path loss and SNR", sprintf('FSPL %.1f-%.1f dB, SNR %.1f to %.1f dB for a 0 dBi handheld', min(T.FSPL_dB(vis)), max(T.FSPL_dB(vis)), ...
    min(T.SNR_dB(vis)), max(T.SNR_dB(vis))), "cellular links often 10-30 dB SNR", ...
    "low MCS (QPSK), slot aggregation/repetition, larger satellite antennas (Rel-18 coverage enhancements)";
 "Multipath (delay spread)", sprintf('%s, DS %.0f ns: much shorter than the CP', opt.Channel, opt.DelaySpread*1e9), ...
    "TDL-A..E 30-300 ns+", "LOS-dominated: SCS 30 kHz halves the CP but doubles the CFO tolerance (the better trade for LEO)"};
CH = cell2table(ch, 'VariableNames', ["Challenge" "ThisPass" "TerrestrialNR" "Rel17NTN"]);
CH.Challenge = string(CH.Challenge); CH.ThisPass = string(CH.ThisPass);
CH.TerrestrialNR = string(CH.TerrestrialNR); CH.Rel17NTN = string(CH.Rel17NTN);
end

%% ===== coexistence ======================================================================

function [IF, T] = coexistence(opt, T, LL, B, orb) %#ok<INUSD>
S = loadInterference(opt);
fv = opt.CarrierFreq; bw = opt.ChannelBW; Bocc = B.Bocc;
vic = fv + [-1 1]*Bocc/2;
adjL = fv - bw + [-1 1]*Bocc/2; adjU = fv + bw + [-1 1]*Bocc/2;
f = S.f(:); P = S.P(:);
enbw = 1.065*S.rbw;                                     % Gaussian RBW filter: ENBW / RBW = 1.065
lin = 10.^(P/10)/enbw;                                  % mW/Hz
% emission parameters (centre, symbol rate = -3 dB width, flat-top PSD)
[fcI, RsI, pkPsd] = emissionShape(f, P, lin, S);
% analyzer floor: median PSD away from the emission (main lobe + 3rd/5th-order regrowth < 3 Rs)
off = abs(f - fcI) > 3*RsI + S.rbw;
if nnz(off) >= max(20, 0.2*numel(f))
    flPsd = median(lin(off)); flSrc = sprintf('median of the %d points > 3 Rs from the carrier', nnz(off));
else
    q = sort(lin); flPsd = median(q(1:max(3, ceil(0.1*numel(q)))));
    flSrc = 'lowest 10 % of the points (span mostly occupied: floor biased low, estimate biased high)';
    S.notes(end + 1, 1) = "Analyzer floor taken from the " + flSrc + ".";
end
fl = 10*log10(flPsd*enbw);
maskF = @(x) maskPsd(x, fcI, RsI, pkPsd);
[rawV, exV, mkV, covV, detV] = bandPower(f, lin, flPsd, maskF, vic);
[rawL, exL, mkL, covL] = bandPower(f, lin, flPsd, maskF, adjL);
[rawU, exU, mkU, covU] = bandPower(f, lin, flPsd, maskF, adjU);
if covV < 0.5
    S.notes(end + 1, 1) = sprintf(['Only %.0f %% of the victim channel lies inside the measured span: the leakage is ' ...
        'mostly the SFCG/ECSS mask fill (a worst case for a mask-compliant transmitter), not a measurement. ' ...
        'Use a spectrum that covers %.1f-%.1f MHz (average detector), e.g. ff_srs4_analyzer with ''WideScan'' or a ' ...
        'fieldfox_live / ff_spectrum_pro recording.'], 100*covV, (fv - 1.5*bw)/1e6, (fv + 1.5*bw)/1e6);
end
acs = 10^(opt.ACS/10);
leakUB = rawV + mkV; leakEst = exV + mkV;
adjUB = (rawL + mkL + rawU + mkU)/acs; adjEst = (exL + mkL + exU + mkU)/acs;
% coupling: measured level -> level at the 5G UE
cpl = opt.InterfCoupling + opt.UEGain - opt.InterfMeasGain; cplNote = string(sprintf('%+.1f dB (InterfCoupling + UEGain - InterfMeasGain)', cpl));
dMeas = opt.InterfMeasDistance; if isempty(dMeas), dMeas = S.measDistance; end
dVic = opt.InterfVictimDistance;
cplD = 0;
if ~isempty(dMeas) && ~isempty(dVic) && dMeas > 0 && dVic > 0
    cplD = 20*log10(dMeas/dVic); cplNote = cplNote + sprintf(', %+.1f dB free space %g m -> %g m', cplD, dMeas, dVic);
end
k = 10^((cpl + cplD)/10);
N = 10^((-228.6 + 10*log10(B.Tsys) + 10*log10(Bocc) + 30)/10);         % mW, UE noise in Bocc
I_est = (leakEst + adjEst)*k; I_ub = (leakUB + adjUB)*k;
IN_est = 10*log10(I_est/N); IN_ub = 10*log10(max(I_ub, 1e-30)/N);
IN0_ub = 10*log10(max((leakUB + adjUB)*10^(cpl/10), 1e-30)/N);         % at the measuring distance
dReq = NaN;
if ~isempty(dMeas) && dMeas > 0
    dReq = dMeas*10^(max(IN0_ub - opt.INR, 0)/20);
end
% SINR and throughput over the pass
vis = T.Visible;
T.SINR_est_dB = T.SNR_dB - 10*log10(1 + I_est/N);
T.SINR_ub_dB = T.SNR_dB - 10*log10(1 + I_ub/N);
TL = table();
dt = seconds(median(diff(T.Time)));
for scs = opt.SCS
    curves = curvesFor(LL, scs, opt);
    thr0 = T.("Thr_SCS" + scs + "_Mbps");
    [te, ~] = adaptThroughput(T.SINR_est_dB, curves); te(~vis) = 0;
    [tu, ~] = adaptThroughput(T.SINR_ub_dB, curves); tu(~vis) = 0;
    T.("ThrInterfEst_SCS" + scs + "_Mbps") = te; T.("ThrInterfUB_SCS" + scs + "_Mbps") = tu;
    TL = [TL; table(scs, sum(thr0(vis))*dt/8, sum(te(vis))*dt/8, sum(tu(vis))*dt/8, ...
        100*(1 - sum(te(vis))/max(sum(thr0(vis)), eps)), 100*(1 - sum(tu(vis))/max(sum(thr0(vis)), eps)), ...
        'VariableNames', ["SCS_kHz" "Volume_MB" "VolumeInterfEst_MB" "VolumeInterfUB_MB" "LossEst_pct" "LossUB_pct"])]; %#ok<AGROW>
end
% I/N and pass volume loss vs interferer-UE separation (free space from the measuring distance)
VD = table();
if ~isempty(dMeas) && dMeas > 0
    d = logspace(log10(dMeas), log10(max([100*dMeas, 3*dReq, 1e3])), 40).';
    k0 = 10^(cpl/10)*(dMeas./d).^2;
    curves = curvesFor(LL, opt.SCS(1), opt);
    thr0 = T.("Thr_SCS" + opt.SCS(1) + "_Mbps"); v0 = max(sum(thr0(vis)), eps);
    lossE = zeros(size(d)); lossU = lossE;
    for j = 1:numel(d)
        te = adaptThroughput(T.SNR_dB - 10*log10(1 + (leakEst + adjEst)*k0(j)/N), curves); te(~vis) = 0;
        tu = adaptThroughput(T.SNR_dB - 10*log10(1 + (leakUB + adjUB)*k0(j)/N), curves); tu(~vis) = 0;
        lossE(j) = 100*(1 - sum(te(vis))/v0); lossU(j) = 100*(1 - sum(tu(vis))/v0);
    end
    VD = table(d, 10*log10((leakEst + adjEst)*k0/N), 10*log10(max((leakUB + adjUB)*k0, 1e-30)/N), ...
        lossE, lossU, 'VariableNames', ["Distance_m" "IN_est_dB" "IN_ub_dB" "LossEst_pct" "LossUB_pct"]);
end
dBm = @(x) 10*log10(max(x, 1e-30));
% can this measurement show compliance at all? analyzer floor referred to the UE vs the UE noise density
n0 = -228.6 + 10*log10(B.Tsys) + 30;                                    % dBm/Hz
flUE = dBm(flPsd) + cpl + cplD;
if flUE - n0 > opt.INR
    dd = flUE - n0; rel = ternary(dd >= 0, sprintf('%.1f dB above', dd), sprintf('only %.1f dB below', -dd));
    S.notes(end + 1, 1) = sprintf(['The analyzer floor referred to the UE (%.1f dBm/Hz) is %s the UE noise ' ...
        'density (%.1f dBm/Hz): the upper bound cannot prove I/N <= %g dB. Lower the floor (preamp, smaller RBW, ' ...
        'directional antenna, shorter distance) or measure conducted and add the path loss.'], flUE, rel, n0, opt.INR);
end
fmtE = @(x) ternary(x > 0, sprintf("%.1f dBm", dBm(x)), "not detected above the analyzer floor");
IF = struct('done', true, 'source', S.source, 'notes', S.notes, 'f', f.', 'P', P.', 'rbw', S.rbw, ...
    'detector', S.detector, 'floor_dBm', fl, 'floorPsd_dBmHz', dBm(flPsd), 'floorSource', string(flSrc), ...
    'fc', fcI, 'Rs', RsI, 'leakDetected', detV, ...
    'peakPsd_dBmHz', dBm(pkPsd), 'victim', vic, 'adjLower', adjL, 'adjUpper', adjU, ...
    'coverage', [covV covL covU], 'leakEst_dBm', dBm(leakEst), 'leakUB_dBm', dBm(leakUB), ...
    'leakMask_dBm', dBm(mkV), 'adjEst_dBm', dBm(adjEst), 'adjUB_dBm', dBm(adjUB), 'ACS', opt.ACS, ...
    'coupling_dB', cpl + cplD, 'couplingNote', cplNote, 'noise_dBm', dBm(N), 'I_est_dBm', dBm(I_est), ...
    'I_ub_dBm', dBm(I_ub), 'IN_est_dB', IN_est, 'IN_ub_dB', IN_ub, 'INR', opt.INR, ...
    'desense_est_dB', 10*log10(1 + I_est/N), 'desense_ub_dB', 10*log10(1 + I_ub/N), ...
    'measDistance', dMeas, 'requiredDistance_m', dReq, 'throughput', TL, 'vsDistance', VD);
IF.summary = table(["Victim channel"; "Interferer"; "Analyzer floor"; "TX leakage in channel (estimate)"; ...
    "TX leakage in channel (upper bound)"; "  of which SFCG-mask fill (outside span)"; ...
    "Adjacent power / ACS (estimate)"; "Adjacent power / ACS (upper bound)"; "Coupling to the UE"; ...
    "UE noise (Bocc)"; "I/N estimate"; "I/N upper bound"; "SINR loss (estimate / upper bound)"; ...
    "Separation for I/N <= INR (upper bound)"], ...
    [sprintf("%.2f-%.2f MHz (%.2f MHz occupied)", vic/1e6, Bocc/1e6); ...
     sprintf("fc %.4f MHz, Rs %.0f kBd (-3 dB width), %s", fcI/1e6, RsI/1e3, S.source); ...
     sprintf("%.1f dBm in RBW %s (%.1f dBm/Hz, %s)", fl, fmtHz(S.rbw), dBm(flPsd), flSrc); ...
     fmtE(leakEst) + sprintf(" (measured coverage %.0f %%)", 100*covV); sprintf("%.1f dBm", dBm(leakUB)); ...
     ternary(mkV > 0, sprintf("%.1f dBm", dBm(mkV)), "none (channel inside the measured span)"); ...
     fmtE(adjEst) + sprintf(" (ACS %g dB)", opt.ACS); sprintf("%.1f dBm", dBm(adjUB)); cplNote; ...
     sprintf("%.1f dBm (T_sys %.0f K)", dBm(N), B.Tsys); ...
     ternary(I_est > 0, sprintf("%.1f dB", IN_est), "no interference measurable above the analyzer floor"); ...
     sprintf("%.1f dB", IN_ub); ...
     sprintf("%.2f / %.2f dB", 10*log10(1 + I_est/N), 10*log10(1 + I_ub/N)); ...
     ternary(isnan(dReq), "unknown (give InterfMeasDistance)", sprintf("%.0f m (free space, INR %g dB)", dReq, opt.INR))], ...
    'VariableNames', ["Item" "Value"]);
end

function curves = curvesFor(LL, scs, opt)
% BLER curves (UE pre-comp case) per MCS, with the peak rate of the FULL channel (nrTBS at its NRB)
mt = nrPDSCHMCSTables; mt = mt.QAM64Table;
L = LL(LL.SCS_kHz == scs & LL.Precomp == "UE pre-comp", :);
nrbF = nrbFor(opt.ChannelBW, scs); slotDur = 1e-3*15/scs;
curves = struct('mcs', {}, 'x', {}, 'y', {}, 'rate', {});
for m = unique(L.MCS).'
    Lm = sortrows(L(L.MCS == m, :), 'SNR_dB'); r = mt(mt.MCSIndex == m, :);
    [car, pd] = pdschSetup(scs, nrbF, string(r.Modulation));
    [~, pinfo] = nrPDSCHIndices(car, pd);
    tbsF = nrTBS(char(pd.Modulation), 1, nrbF, pinfo.NREPerPRB, r.TargetCodeRate);
    curves(end + 1) = struct('mcs', m, 'x', Lm.SNR_dB, 'y', Lm.BLER, 'rate', tbsF/slotDur/1e6); %#ok<AGROW>
end
end

function [car, pd] = pdschSetup(scs, nrb, modulation)
% Carrier + PDSCH used everywhere: 1 layer, all PRBs, symbols 2-13 (2-symbol CORESET reserved),
% DM-RS type A at symbol 2 with 2 additional positions (3 DM-RS symbols, for Doppler)
car = nrCarrierConfig('SubcarrierSpacing', scs, 'NSizeGrid', nrb);
pd = nrPDSCHConfig;
pd.Modulation = char(modulation); pd.NumLayers = 1; pd.PRBSet = 0:nrb - 1; pd.SymbolAllocation = [2 12];
pd.DMRS.DMRSTypeAPosition = 2; pd.DMRS.DMRSAdditionalPosition = 2;
end

function [fc, Rs, pk] = emissionShape(f, P, lin, S)
[~, i] = max(P);
top = P(i);
lo = i; while lo > 1 && P(lo - 1) > top - 3, lo = lo - 1; end
hi = i; while hi < numel(P) && P(hi + 1) > top - 3, hi = hi + 1; end
fc = (f(lo) + f(hi))/2; Rs = max(f(hi) - f(lo), S.rbw);
if isfinite(S.fcHint), fc = S.fcHint; end
if isfinite(S.RsHint), Rs = S.RsHint; end
m = abs(f - fc) <= max(0.25*Rs, S.rbw);
pk = median(lin(m));
end

function p = maskPsd(f, fc, Rs, pk)
% SFCG Rec. 21-2 / ECSS-E-ST-50-05C Fig. 6-5 mask relative to the peak PSD (suppressed carrier):
% < 2 Msps: 0 dB to 0.5 Rs, -30 dB at 3 Rs, -60 dB at 8 Rs; >= 2 Msps: -30 dB at 1.4 Rs, -60 at 3 Rs
if Rs < 2e6, x = [0 0.5 3 8 1e6]; else, x = [0 0.5 1.4 3 1e6]; end
m = interp1(x, [0 0 -30 -60 -60], abs(f - fc)/Rs, 'linear', -60);
p = pk*10.^(m/10);
end

function [raw, ex, mk, cov, det] = bandPower(f, lin, flPsd, maskF, band)
% power (mW) in [band(1) band(2)]: measured part (raw, and with the analyzer floor subtracted) +
% mask-filled rest. The floor is subtracted from the band total (noise fluctuations cancel); an
% excess below 20 % of the floor power in the band is treated as "not detected" (ex = 0, det = false).
fg = linspace(band(1), band(2), 4001).'; dfg = fg(2) - fg(1);
df = median(diff(f));
covd = fg >= f(1) - df/2 & fg <= f(end) + df/2;
cov = mean(covd);
raw = 0; ex = 0; mk = 0; det = true;
if any(covd)
    raw = sum(interp1(f, lin, fg(covd), 'linear', 'extrap'))*dfg;
    flP = flPsd*nnz(covd)*dfg;
    ex = raw - flP;
    if ex < 0.2*flP, ex = 0; det = false; end
end
if any(~covd), mk = sum(maskF(fg(~covd)))*dfg; end
end

function S = loadInterference(opt)
x = opt.Interference;
S = struct('f', [], 'P', [], 'rbw', NaN, 'detector', "unknown", 'source', "", 'notes', strings(0, 1), ...
    'fcHint', NaN, 'RsHint', NaN, 'measDistance', []);
if (isstring(x) || ischar(x)) && strcmpi(string(x), "demo")
    S = demoInterference(opt, S); return
end
if isstring(x) || ischar(x)
    file = string(x);
    if ~isfile(file), error('ff_5g_ntn:interference', 'Interference file not found: %s', file); end
    x = load(file); src = "file " + file;
else
    src = "workspace struct";
end
if isfield(x, 'run_data'), D = x.run_data; A = []; if isfield(x, 'results'), A = x.results; end
elseif isfield(x, 'data') && isstruct(x.data) && isfield(x.data, 'acq'), D = x.data; A = []; if isfield(x, 'results'), A = x.results; end
else, D = []; A = [];
end
if ~isempty(D)
    S = fromRun(D, A, opt, S); S.source = "ff_srs4_analyzer run (" + src + ", " + S.source + ")"; return
end
if ~isfield(x, 'f') || ~isfield(x, 'P'), error('ff_5g_ntn:interference', 'The interference data needs f and P (or an ff_srs4_analyzer run).'); end
S.f = double(x.f(:)).'; P = double(x.P);
if size(P, 2) ~= numel(S.f), P = P.'; end
if size(P, 2) ~= numel(S.f), error('ff_5g_ntn:interference', 'P must have one column per frequency in f.'); end
if numel(unique(S.f)) < 10
    error('ff_5g_ntn:interference', 'A zero-span (or < 10 point) recording has no spectrum: use a swept recording.');
end
[S.f, i] = sort(S.f); P = P(:, i);
S.P = 10*log10(mean(10.^(P/10), 1, 'omitnan'));                         % power average over the sweeps
if any(~isfinite(S.P)), error('ff_5g_ntn:interference', 'The recording has frequency points without data.'); end
if isfield(x, 'settings') && isstruct(x.settings)
    st = x.settings;
    if isfield(st, 'rbw'), S.rbw = double(st.rbw); end
    if isfield(st, 'detector'), S.detector = string(st.detector); end
    if isfield(st, 'extgInstr') && isfield(st, 'extGainTotal')                 % ff_srs4_analyzer rec_*.mat
        S.P = S.P + st.extgInstr - st.extGainTotal;
        if st.extGainTotal ~= 0
            S.notes(end + 1, 1) = sprintf(['Recording made with ExtGain %g dB: converted to DUT-referred levels ' ...
                '(transmitter output). Put TX antenna gain - path loss (+ UE gain - measuring gain) into InterfCoupling.'], ...
                st.extGainTotal);
        end
    end
end
if ~isfinite(S.rbw)
    S.rbw = 3*(S.f(end) - S.f(1))/(numel(S.f) - 1);
    S.notes(end + 1, 1) = "RBW not in the recording: assumed 3 bins (" + fmtHz(S.rbw) + ").";
end
S.source = "recording (" + src + ")";
S = detectorNote(S);
end

function S = fromRun(D, A, opt, S)
% a run in which the analyzer found no downlink (SRS-4 off, wrong band) has no interferer:
% its coarse sweep would otherwise be read as an "emission" at the highest noise peak
noDL = ~isempty(A) && isfield(A, 'ok') && ~A.ok;
if isfield(D, 'search') && isstruct(D.search) && isfield(D.search, 'found') && ~D.search.found, noDL = true; end
if noDL
    error('ff_5g_ntn:interference', ['The ff_srs4_analyzer run found no downlink (SRS-4 off or outside its ' ...
        '''Band''): there is no interferer to assess. Use a run with the SRS-4 transmitting.']);
end
names = string(fieldnames(D.acq)).';
names = names(~startsWith(names, "harm"));
fv = opt.CarrierFreq; want = fv + [-1.5 1.5]*opt.ChannelBW; chan = fv + [-0.5 0.5]*opt.ChannelBW;
best = ""; bestScore = -Inf;
% Rank the run's spectra: victim-channel coverage first, then the average detector and an RBW
% that resolves the channel edge (the 30 MHz-6.5 GHz max-hold scan at RBW 1 MHz smears the
% interferer's main lobe into the channel: +20-30 dB of false leakage), then the adjacent channels.
for nm = names
    a = D.acq.(nm);
    if ~isfield(a, 'settings') || (isfield(a.settings, 'zeroSpan') && a.settings.zeroSpan), continue, end
    cv = @(b) max(0, min(a.f(end), b(2)) - max(a.f(1), b(1)))/diff(b);
    isAvg = isfield(a.settings, 'detector') && contains(lower(string(a.settings.detector)), ["aver" "rms"]);
    score = 100*cv(chan) + 20*cv(want) + 15*isAvg + 10*(a.settings.rbw <= opt.ChannelBW/50) ...
        - log10(max(a.settings.rbw, 1));
    if score > bestScore, best = nm; bestScore = score; end
end
if best == "", error('ff_5g_ntn:interference', 'The run has no swept spectrum.'); end
a = D.acq.(best);
eg = 0; if isfield(D, 'extGain'), eg = D.extGain; end
ei = 0; if isfield(a.settings, 'extgInstr'), ei = a.settings.extgInstr; end
P = a.P + ei - eg;                                                      % DUT-referred (= port level when ExtGain = 0)
S.f = double(a.f(:)).'; S.P = 10*log10(mean(10.^(P/10), 1)); S.rbw = a.settings.rbw;
S.detector = string(a.settings.detector); S.source = "acquisition '" + best + "'";
if eg ~= 0
    S.notes(end + 1, 1) = sprintf(['The run used ExtGain %g dB: levels are DUT-referred (transmitter output), not ' ...
        'received. Put TX antenna gain - path loss (+ UE gain - measuring gain) into InterfCoupling.'], eg);
end
if isfield(D, 'opt') && isfield(D.opt, 'Distance') && ~isempty(D.opt.Distance), S.measDistance = D.opt.Distance; end
if ~isempty(A) && isfield(A, 'spec') && isstruct(A.spec)
    if isfield(A.spec, 'fc') && isfinite(A.spec.fc), S.fcHint = A.spec.fc; end
    if isfield(A.spec, 'Rs') && isfinite(A.spec.Rs), S.RsHint = A.spec.Rs; end
end
S = detectorNote(S);
end

function S = detectorNote(S)
d = lower(S.detector);
if contains(d, ["pos" "peak" "max"])
    S.notes(end + 1, 1) = "Peak / max-hold detector: levels read high for noise-like signals (conservative).";
elseif ~contains(d, ["aver" "rms"])
    S.notes(end + 1, 1) = "Detector " + S.detector + ": the average (RMS) detector gives unbiased power.";
end
end

function S = demoInterference(opt, S)
% Synthetic SRS-4 (CCSDS TM, QPSK, RS+CC, RRC 0.2, 1 MBd) at 2201.0 MHz through a mild PA
% (modified Rapp), seen by a FieldFox model: 10 kHz RBW, average detector, -135 dBm/Hz floor.
Rs = 1e6; sps = 24; fs = Rs*sps; fc = 2201.0e6; Prx = -25;
tm = ccsdsTMWaveformGenerator('ChannelCoding', 'concatenated', 'Modulation', 'QPSK', 'RolloffFactor', 0.2, ...
    'SamplesPerSymbol', sps, 'HasRandomizer', true, 'HasASM', true);
rs = RandStream('mt19937ar', 'Seed', opt.Seed);
x = zeros(0, 1);
for k = 1:6, x = [x; tm(randi(rs, [0 1], tm.NumInputBits, 1))]; end %#ok<AGROW>
x = x/rms(x);
pa = comm.MemorylessNonlinearity('Method', 'Modified Rapp model', 'Smoothness', 1.5, 'OutputSaturationLevel', 1.1, 'LinearGain', 0);
y = pa(x);
[pxx, fb] = pwelch(y, hann(8192), 4096, 8192, fs, 'centered', 'psd');
pxx = pxx/(sum(pxx)*(fb(2) - fb(1)))*10^(Prx/10);                      % mW/Hz, total Prx dBm
rbw = 10e3; f = (2185e6:10e3:2215e6);
psd = interp1(fb + fc, pxx, f, 'linear', 0);
floorPsd = 10^(-135/10);
nSw = 10;
P = zeros(nSw, numel(f));
for k = 1:nSw
    P(k, :) = 10*log10((psd + floorPsd.*max(0.2, 1 + 0.3*randn(rs, 1, numel(f))))*1.065*rbw);
end
S.f = f; S.P = 10*log10(mean(10.^(P/10), 1)); S.rbw = rbw; S.detector = "AVER (demo)";
S.source = sprintf("SYNTHETIC demo: SRS-4 at %.1f MHz, 1 MBd QPSK RRC 0.2, %g dBm received, 3 m away", fc/1e6, Prx);
S.measDistance = 3;
S.notes(end + 1, 1) = "Synthetic SRS-4 spectrum (ccsdsTMWaveformGenerator + comm.MemorylessNonlinearity): not a measurement.";
end

function s = summarizeInput(x)
if isempty(x), s = "";
elseif isstring(x) || ischar(x), s = string(x);
else, s = "struct";
end
end

%% ===== summaries ===========================================================================

function P = passSummary(T, orb, opt)
vis = T.Visible; tv = T.Time(vis);
[mx, i] = max(T.Elevation_deg);
P = struct('satellite', orb.name, 'start', tv(1), 'stop', tv(end), 'duration_s', seconds(tv(end) - tv(1)), ...
    'maxElevation', mx, 'timeMaxElevation', T.Time(i), 'rangeMin_km', min(T.Range_km(vis)), ...
    'rangeMax_km', max(T.Range_km(vis)), 'delayMin_ms', min(T.Delay_ms(vis)), 'delayMax_ms', max(T.Delay_ms(vis)), ...
    'dopplerMax_Hz', max(abs(T.Doppler_Hz(vis))), 'dopplerRateMax_Hzps', max(abs(T.DopplerRate_Hzps(vis))), ...
    'snrMin_dB', min(T.SNR_dB(vis)), 'snrMax_dB', max(T.SNR_dB(vis)), 'minElevation', opt.MinElevation, ...
    'boundDoppler_Hz', NaN, 'boundRange_km', NaN);
if orb.circular
    P.boundDoppler_Hz = max(dopplerShiftCircularOrbit(opt.MinElevation, opt.Altitude, opt.Site(3), opt.CarrierFreq));
    P.boundRange_km = slantRangeCircularOrbit(opt.MinElevation, opt.Altitude, opt.Site(3))/1e3;
end
end

function A = assumptionList(opt, orb, atm, cond)
A = [ ...
    "Orbit: " + orb.source;
    sprintf("Satellite EIRP density %g dBW/MHz (%s), max gain %g dBi, earth-fixed beam centred on the UE, scan loss exponent %g", ...
        orb.eirpDensity, orb.eirpSource, opt.SatGain, opt.ScanLossExp);
    sprintf("Handheld UE: %g dBi, NF %g dB, antenna %g K (TR 38.821 Table 6.1.1.1-3); polarisation %g dB, shadowing %g dB, scintillation %g dB", ...
        opt.UEGain, opt.UENF, opt.UETemp, opt.PolLoss, opt.ShadowMargin, opt.ScintLoss);
    "Troposphere: " + atm.source;
    sprintf("Band %s, carrier %.2f MHz, channel %g MHz, SCS %s kHz", opt.Band, opt.CarrierFreq/1e6, opt.ChannelBW/1e6, strjoin(string(opt.SCS), '/'));
    sprintf("Link level: %s (DS %g ns), UE speed %g km/h (fading Doppler %.1f Hz), %d Rx, PDSCH symbols 2-13, 3 DM-RS symbols, no HARQ retransmission, practical estimation", ...
        opt.Channel, opt.DelaySpread*1e9, opt.UESpeed, cond.fdUE, opt.RxAntennas);
    sprintf("Doppler at %.1f deg: %.1f kHz, rate %.0f Hz/s; residual after UE pre-compensation %.0f Hz (update %g s, ephemeris velocity error %g m/s)", ...
        cond.el, cond.fD/1e3, cond.rate, cond.resid, opt.PrecompInterval, opt.PrecompVelError);
    "Throughput: ideal link adaptation over the simulated BLER curves, full-buffer FDD downlink, no control/SSB overhead beyond the 2-symbol CORESET"];
end

function txt = conclusions(R, opt, cond, tol)
p = R.pass; txt = strings(0, 1);
txt(end + 1) = sprintf(['Pass of %s: %s-%s (UTC+7), %.0f s above %g deg, peak %.1f deg; slant range %.0f-%.0f km, ' ...
    'one-way delay %.2f-%.2f ms (RTT up to %.1f ms).'], p.satellite, string(p.start, 'HH:mm:ss'), string(p.stop, 'HH:mm:ss'), ...
    p.duration_s, opt.MinElevation, p.maxElevation, p.rangeMin_km, p.rangeMax_km, p.delayMin_ms, p.delayMax_ms, 2*p.delayMax_ms);
txt(end + 1) = sprintf('Doppler at %.2f MHz: up to +-%.1f kHz (%.1f ppm), rate up to %.0f Hz/s', opt.CarrierFreq/1e6, ...
    p.dopplerMax_Hz/1e3, p.dopplerMax_Hz/opt.CarrierFreq*1e6, p.dopplerRateMax_Hzps);
if isfinite(p.boundDoppler_Hz)
    txt(end) = txt(end) + sprintf(' (dopplerShiftCircularOrbit bound at %g deg: %.1f kHz)', opt.MinElevation, p.boundDoppler_Hz/1e3);
end
txt(end) = txt(end) + ".";
txt(end + 1) = sprintf('Handheld SNR (EIRP density %g dBW/MHz): %.1f dB at the mask to %.1f dB at the peak.', ...
    R.orbit.eirpDensity, p.snrMin_dB, p.snrMax_dB);
for i = 1:height(R.throughput)
    r = R.throughput(i, :);
    txt(end + 1) = sprintf(['SCS %d kHz, %g MHz: mean %.2f Mbit/s (peak %.2f), %.1f MB per pass with UE Doppler ' ...
        'pre-compensation; without it %.1f MB (outage %.0f %% of the pass, CFO tolerance %.1f kHz).'], ...
        r.SCS_kHz, opt.ChannelBW/1e6, r.MeanThr_Mbps, r.MaxThr_Mbps, r.Volume_MB, r.VolumeNoPrecomp_MB, ...
        r.OutageNoPrecomp_pct, r.CFOTolerance_Hz/1e3); %#ok<AGROW>
end
for i = 1:numel(tol.scs)
    txt(end + 1) = sprintf(['CFO tolerance SCS %d kHz (%s): +-%.2f kHz without and +-%.2f kHz with the CP fine ' ...
        'estimator, against +-%.1f kHz LEO Doppler: pre-compensation (residual ~%.0f Hz) is essential.'], ...
        tol.scs(i), tol.source(i), tol.fineOff(i)/1e3, tol.fineOn(i)/1e3, cond.maxDoppler/1e3, cond.residMax); %#ok<AGROW>
end
L = R.linkLevel; n = L(L.Precomp == "none" & L.Simulated, :);
if ~isempty(n)
    txt(end + 1) = sprintf('Link level without pre-compensation (simulated CFO %.1f kHz): BLER %.2f at %g dB SNR (MCS %d).', ...
        n.CFO_Hz(1)/1e3, n.BLER(end), n.SNR_dB(end), n.MCS(1));
end
txt(end + 1) = sprintf(['Timing: RTT %.1f ms exceeds the NR timing-advance range (2.0 ms at SCS 15, 1.0 ms at 30): ' ...
    'UE-specific TA from GNSS + ephemeris, K_offset and up to 32 HARQ processes (Rel-17).'], 2*p.delayMax_ms);
if isfield(R.interference, 'done') && R.interference.done
    I = R.interference;
    est = ternary(isfinite(I.IN_est_dB), sprintf('%.1f dB (estimate)', I.IN_est_dB), 'not measurable above the analyzer floor (estimate)');
    txt(end + 1) = sprintf(['Coexistence: I/N %s / %.1f dB (upper bound) -> SINR loss %.2f / %.2f dB, ' ...
        'pass volume loss %.0f / %.0f %% (SCS %d).'], est, I.IN_ub_dB, I.desense_est_dB, I.desense_ub_dB, ...
        I.throughput.LossEst_pct(1), I.throughput.LossUB_pct(1), I.throughput.SCS_kHz(1));
    if I.coverage(1) < 0.5
        txt(end) = txt(end) + sprintf(' Mostly SFCG-mask fill: only %.0f %% of the victim channel was measured.', 100*I.coverage(1));
    end
    if isfinite(I.requiredDistance_m)
        txt(end + 1) = sprintf('Separation for I/N <= %g dB: about %.0f m in free space (upper bound).', opt.INR, I.requiredDistance_m);
    end
end
txt = txt(:);
end

%% ===== figures ==============================================================================

function files = makeFigures(R, opt, cond, ~)
files = strings(0, 1);
T = R.timeline; vis = T.Visible; up = T.Elevation_deg > 0;
% 1. pass geometry
fig = newFig(opt, [60 60 1500 860]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s over %s: %s', R.pass.satellite, opt.SiteName, string(R.pass.start, 'yyyy-MM-dd')));
pax = polaraxes(tl); pax.Layout.Tile = 1;
polarplot(pax, deg2rad(T.Azimuth_deg(up)), 90 - T.Elevation_deg(up), 'Color', [0.6 0.6 0.6]); hold(pax, 'on');
polarplot(pax, deg2rad(T.Azimuth_deg(vis)), 90 - T.Elevation_deg(vis), 'b', 'LineWidth', 2);
pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise'; pax.RLim = [0 90];
pax.RTick = [0 30 60 80]; pax.RTickLabel = {'90', '60', '30', '10'};
title(pax, 'Sky track (az, el)');
ax = nexttile(tl, 2);
yyaxis(ax, 'left'); plot(ax, T.Time, T.Elevation_deg, 'LineWidth', 1.5); ylabel(ax, 'Elevation (deg)');
yline(ax, opt.MinElevation, ':');
yyaxis(ax, 'right'); plot(ax, T.Time, T.Range_km, 'LineWidth', 1.5); ylabel(ax, 'Slant range (km)');
grid(ax, 'on'); title(ax, 'Elevation and range'); xlabel(ax, 'Time (UTC+7)');
ax = nexttile(tl, 3);
yyaxis(ax, 'left'); plot(ax, T.Time, T.Doppler_Hz/1e3, 'LineWidth', 1.5); ylabel(ax, 'Doppler (kHz)');
yline(ax, [-1 1]*opt.SCS(1)/2, ':');
yyaxis(ax, 'right'); plot(ax, T.Time, T.DopplerRate_Hzps, 'LineWidth', 1.5); ylabel(ax, 'Doppler rate (Hz/s)');
grid(ax, 'on'); title(ax, sprintf('Doppler at %.1f MHz (dotted: +-SCS%d/2)', opt.CarrierFreq/1e6, opt.SCS(1)));
ax = nexttile(tl, 4);
plot(ax, T.Time, T.Delay_ms, T.Time, 2*T.Delay_ms, 'LineWidth', 1.5); hold(ax, 'on');
yline(ax, 2.003, '--', 'max TA SCS15'); yline(ax, 1.0015, '--', 'max TA SCS30');
grid(ax, 'on'); ylabel(ax, 'ms'); legend(ax, 'one-way delay', 'RTT (service link)', 'Location', 'best');
title(ax, 'Delay vs NR timing-advance range');
gx = geoaxes(tl); gx.Layout.Tile = 5;
geoplot(gx, T.SatLat_deg(up), T.SatLon_deg(up), 'Color', [0.6 0.6 0.6]); hold(gx, 'on');
geoplot(gx, T.SatLat_deg(vis), T.SatLon_deg(vis), 'b', 'LineWidth', 2);
geoplot(gx, opt.Site(1), opt.Site(2), 'r^', 'MarkerFaceColor', 'r');
geobasemap(gx, 'grayland'); title(gx, 'Ground track (visible part in blue)');
ax = nexttile(tl, 6);
plot(ax, T.Time, T.SNR_dB, 'k', 'LineWidth', 1.5); hold(ax, 'on');
plot(ax, T.Time, -T.ScanLoss_dB, T.Time, -(T.FSPL_dB - min(T.FSPL_dB(vis))), T.Time, -T.Atm_dB);
grid(ax, 'on'); ylabel(ax, 'dB'); legend(ax, 'SNR (handheld)', 'scan loss', 'extra FSPL vs zenith', 'atmosphere', 'Location', 'best');
title(ax, 'Link budget per time step');
files(end + 1) = saveFig(fig, fullfile(opt.OutDir, "fig_pass.png"), opt);

% 2. BLER / throughput vs SNR
L = R.linkLevel; nS = numel(opt.SCS);
fig = newFig(opt, [80 80 1400 420*nS]);
tl = tiledlayout(fig, nS, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('PDSCH link level, %s, %g MHz simulated', opt.Channel, opt.SimBW/1e6));
for s = 1:nS
    scs = opt.SCS(s);
    ax = nexttile(tl); hold(ax, 'on');
    Lu = L(L.SCS_kHz == scs & L.Precomp == "UE pre-comp", :);
    for m = unique(Lu.MCS).'
        q = sortrows(Lu(Lu.MCS == m, :), 'SNR_dB');
        semilogy(ax, q.SNR_dB, max(q.BLER, 1e-3), '-o', 'DisplayName', sprintf('MCS %d %s %.2f', m, q.Modulation(1), q.CodeRate(1)));
    end
    Ln = sortrows(L(L.SCS_kHz == scs & L.Precomp == "none", :), 'SNR_dB');
    if ~isempty(Ln)
        semilogy(ax, Ln.SNR_dB, max(Ln.BLER, 1e-3), 'k--x', 'DisplayName', sprintf('MCS %d, no pre-comp (%.1f kHz)', Ln.MCS(1), Ln.CFO_Hz(1)/1e3));
    end
    set(ax, 'YScale', 'log'); ylim(ax, [1e-3 1.05]); grid(ax, 'on'); yline(ax, 0.1, ':', 'HandleVisibility', 'off');
    xlabel(ax, 'SNR per RE (dB)'); ylabel(ax, 'BLER (0 plotted as 1e-3)'); legend(ax, 'Location', 'southwest', 'FontSize', 7);
    title(ax, sprintf('BLER, SCS %d kHz (UE pre-comp residual %.0f Hz)', scs, cond.resid));
    ax = nexttile(tl); hold(ax, 'on');
    for m = unique(Lu.MCS).'
        q = sortrows(Lu(Lu.MCS == m, :), 'SNR_dB');
        plot(ax, q.SNR_dB, q.Thr_Mbps, '-o', 'DisplayName', sprintf('MCS %d', m));
    end
    x = linspace(min(L.SNR_dB), max(L.SNR_dB), 100);
    nrb = nrbFor(opt.SimBW, scs);
    plot(ax, x, nrb*12*scs*1e3*log2(1 + 10.^(x/10))/1e6, 'k:', 'DisplayName', 'Shannon (occupied BW)');
    grid(ax, 'on'); xlabel(ax, 'SNR per RE (dB)'); ylabel(ax, 'Throughput (Mbit/s)');
    legend(ax, 'Location', 'northwest', 'FontSize', 7); title(ax, sprintf('Throughput, SCS %d kHz, %g MHz', scs, opt.SimBW/1e6));
end
files(end + 1) = saveFig(fig, fullfile(opt.OutDir, "fig_bler.png"), opt);

% 3. throughput over the pass
fig = newFig(opt, [100 100 1300 700]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Expected DL throughput over the pass (%g MHz, ideal link adaptation)', opt.ChannelBW/1e6));
ax = nexttile(tl); hold(ax, 'on');
for scs = opt.SCS
    plot(ax, T.Time, T.("Thr_SCS" + scs + "_Mbps"), 'LineWidth', 1.5, 'DisplayName', sprintf('SCS %d, UE pre-comp', scs));
    plot(ax, T.Time, T.("ThrNoPrecomp_SCS" + scs + "_Mbps"), '--', 'DisplayName', sprintf('SCS %d, no pre-comp', scs));
    if isfield(R.interference, 'done') && R.interference.done
        plot(ax, T.Time, T.("ThrInterfUB_SCS" + scs + "_Mbps"), ':', 'LineWidth', 1.5, 'DisplayName', sprintf('SCS %d, + interference (UB)', scs));
    end
end
plot(ax, T.Time, T.Shannon_Mbps, 'k:', 'DisplayName', 'Shannon bound');
grid(ax, 'on'); ylabel(ax, 'Mbit/s'); legend(ax, 'Location', 'northwest', 'FontSize', 8);
ax = nexttile(tl); hold(ax, 'on');
dt = seconds(median(diff(T.Time)));
for scs = opt.SCS
    plot(ax, T.Time, cumsum(T.("Thr_SCS" + scs + "_Mbps"))*dt/8, 'LineWidth', 1.5, 'DisplayName', sprintf('SCS %d', scs));
end
yyaxis(ax, 'right'); plot(ax, T.Time, T.("MCS_SCS" + opt.SCS(1)), '.', 'DisplayName', sprintf('MCS (SCS %d)', opt.SCS(1)));
ylabel(ax, 'selected MCS'); yyaxis(ax, 'left'); ylabel(ax, 'cumulative data (MB)');
grid(ax, 'on'); xlabel(ax, 'Time (UTC+7)'); legend(ax, 'Location', 'northwest');
files(end + 1) = saveFig(fig, fullfile(opt.OutDir, "fig_throughput.png"), opt);

% 4. CFO tolerance
CF = R.cfoSweep;
if ~isempty(CF)
    fig = newFig(opt, [120 120 1100 520]);
    ax = axes(fig); hold(ax, 'on');
    for scs = unique(CF.SCS_kHz).'
        for fine = [false true]
            q = sortrows(CF(CF.SCS_kHz == scs & CF.FineCFO == fine, :), 'CFO_Hz');
            if isempty(q), continue, end
            st = '-o'; if fine, st = '-s'; end
            semilogy(ax, q.CFO_Hz/1e3, max(q.BLER, 1e-3), st, 'LineWidth', 1.2, ...
                'DisplayName', sprintf('SCS %d kHz, fine CFO %s', scs, string(matlab.lang.OnOffSwitchState(fine))));
        end
    end
    set(ax, 'YScale', 'log'); ylim(ax, [1e-3 1.05]); grid(ax, 'on');
    xline(ax, cond.resid/1e3, 'g-', sprintf('residual after pre-comp %.0f Hz', cond.resid), 'HandleVisibility', 'off');
    xline(ax, min(cond.fD, max(CF.CFO_Hz))/1e3, 'r--', sprintf('LEO Doppler at %.0f deg: %.1f kHz%s', cond.el, cond.fD/1e3, ...
        ternary(cond.fD > max(CF.CFO_Hz), ' (off scale)', '')), 'HandleVisibility', 'off');
    yline(ax, 0.1, ':', 'HandleVisibility', 'off');
    xlabel(ax, 'Carrier frequency offset at the UE (kHz)'); ylabel(ax, 'BLER (0 plotted as 1e-3)');
    legend(ax, 'Location', 'southeast');
    title(ax, sprintf('CFO tolerance vs numerology: MCS %d at %g dB SNR (%s)', opt.CFOSweepMCS, opt.CFOSweepSNR, opt.Channel));
    files(end + 1) = saveFig(fig, fullfile(opt.OutDir, "fig_cfo.png"), opt);
end

% 5. coexistence
I = R.interference;
if isfield(I, 'done') && I.done
    fig = newFig(opt, [140 140 1300 760]);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Coexistence: interferer spectrum vs the 5G NTN victim channel');
    ax = nexttile(tl, [1 2]); hold(ax, 'on');
    psd = I.P - 10*log10(1.065*I.rbw);
    plot(ax, I.f/1e6, psd, 'b', 'DisplayName', 'measured PSD (dBm/Hz)');
    fm = linspace(min([I.f I.adjLower(1)]), max([I.f I.adjUpper(2)]), 2000);
    plot(ax, fm/1e6, 10*log10(maskPsd(fm, I.fc, I.Rs, 10^(I.peakPsd_dBmHz/10))), 'm--', 'DisplayName', 'SFCG/ECSS mask (fill)');
    yline(ax, I.floorPsd_dBmHz, 'k:', 'analyzer floor', 'HandleVisibility', 'off');
    yl = ylim(ax);
    patch(ax, I.victim([1 2 2 1])/1e6, yl([1 1 2 2]), [1 0 0], 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName', '5G victim channel');
    patch(ax, I.adjLower([1 2 2 1])/1e6, yl([1 1 2 2]), [1 0.6 0], 'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
        'DisplayName', 'adjacent channels (ACS)');
    patch(ax, I.adjUpper([1 2 2 1])/1e6, yl([1 1 2 2]), [1 0.6 0], 'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    ylim(ax, yl);
    if 2200e6 > fm(1) && 2200e6 < fm(end), xline(ax, 2200, 'k-', 'SRS-4 TX band edge', 'HandleVisibility', 'off'); end
    grid(ax, 'on'); xlabel(ax, 'Frequency (MHz)'); ylabel(ax, 'dBm/Hz'); legend(ax, 'Location', 'northwest', 'FontSize', 8);
    title(ax, sprintf('I/N %.1f dB (estimate) / %.1f dB (upper bound), coupling %+.1f dB', I.IN_est_dB, I.IN_ub_dB, I.coupling_dB));
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, T.Time, T.SNR_dB, 'k', 'LineWidth', 1.5, 'DisplayName', 'SNR');
    plot(ax, T.Time, T.SINR_est_dB, 'b--', 'DisplayName', 'SINR (estimate)');
    plot(ax, T.Time, T.SINR_ub_dB, 'r:', 'LineWidth', 1.5, 'DisplayName', 'SINR (upper-bound interference)');
    grid(ax, 'on'); ylabel(ax, 'dB'); xlabel(ax, 'Time (UTC+7)'); legend(ax, 'Location', 'best');
    title(ax, sprintf('Pass volume loss SCS %d: %.0f %% (estimate) / %.0f %% (upper bound)', I.throughput.SCS_kHz(1), ...
        I.throughput.LossEst_pct(1), I.throughput.LossUB_pct(1)));
    ax = nexttile(tl); hold(ax, 'on');
    V = I.vsDistance;
    if ~isempty(V)
        yyaxis(ax, 'left');
        h1 = plot(ax, V.Distance_m, V.IN_est_dB, '-', 'LineWidth', 1.5);
        h2 = plot(ax, V.Distance_m, V.IN_ub_dB, '--', 'LineWidth', 1.5);
        yline(ax, opt.INR, ':', sprintf('I/N = %g dB', opt.INR)); ylabel(ax, 'I/N (dB)');
        yyaxis(ax, 'right');
        h3 = plot(ax, V.Distance_m, V.LossEst_pct, '-', 'LineWidth', 1);
        h4 = plot(ax, V.Distance_m, V.LossUB_pct, '-.', 'LineWidth', 1);
        ylabel(ax, 'pass volume loss (%)'); ylim(ax, [0 105]);
        set(ax, 'XScale', 'log'); grid(ax, 'on'); xlabel(ax, 'Interferer - UE distance (m), free space');
        xlim(ax, V.Distance_m([1 end]));
        legend(ax, [h1 h2 h3 h4], 'I/N estimate', 'I/N upper bound', 'loss estimate', 'loss upper bound', ...
            'Location', 'northeast', 'FontSize', 7);
        title(ax, sprintf('Separation (measured at %g m)', I.measDistance));
    else
        axis(ax, 'off');
        text(ax, 0.05, 0.5, {'I/N vs distance needs the measuring distance:', 'give ''InterfMeasDistance'' (m).'});
    end
    files(end + 1) = saveFig(fig, fullfile(opt.OutDir, "fig_coexistence.png"), opt);
end
files = files(:);
end

function fig = newFig(opt, pos)
fig = figure('Visible', opt.Visible, 'Color', 'w', 'Position', pos);
end

function file = saveFig(fig, file, opt)
exportgraphics(fig, file, 'Resolution', 130);
if strcmpi(opt.Visible, 'off'), close(fig); end
file = string(file);
end

%% ===== report =================================================================================

function file = makeReport(R, opt)
import mlreportgen.report.*
import mlreportgen.dom.*
fprintf('Writing the %s report ...\n', upper(opt.Report));
rpt = Report(char(fullfile(opt.OutDir, "ntn5g_report")), char(opt.Report));
when = datetime('now', 'TimeZone', 'Asia/Bangkok', 'Format', 'yyyy-MM-dd HH:mm ''(Asia/Bangkok)''');
add(rpt, TitlePage('Title', '5G NR NTN downlink from a LEO satellite', ...
    'Subtitle', sprintf('%s over %s, band %s: geometry, link level, budget, throughput%s', R.pass.satellite, ...
    opt.SiteName, opt.Band, ternary(isfield(R.interference, 'done') && R.interference.done, ', coexistence', '')), ...
    'Author', string(getenv('USERNAME')), 'PubDate', char(when)));
add(rpt, TableOfContents);
png = @(name) char(fullfile(opt.OutDir, name));
img = @(ch, name, cap) addImg(ch, png(name), cap);

ch = Chapter('Title', 'Conclusions');
w = Paragraph(['MODEL RESULTS: 3GPP study parameters (TR 38.821), not a real NTN satellite. The interference part ' ...
    'uses a power spectrum only (no I/Q).']); w.Bold = true; w.Color = 'red'; add(ch, w);
add(ch, UnorderedList(cellstr(R.conclusions)));
if strlength(opt.Notes) > 0, add(ch, Paragraph("Notes: " + opt.Notes)); end
add(ch, rptTable(R.throughput, 'Throughput and data volume per pass', 10));
add(rpt, ch);

ch = Chapter('Title', 'Scenario and assumptions');
add(ch, UnorderedList(cellstr(R.assumptions)));
if ~isempty(R.notes), add(ch, Paragraph('Notes from the run:')); add(ch, UnorderedList(cellstr(R.notes))); end
if istable(R.orbit.passes) && height(R.orbit.passes) > 0, add(ch, rptTable(R.orbit.passes, 'Passes found (TLE search)', 15)); end
add(rpt, ch);

ch = Chapter('Title', 'Pass geometry and link budget');
img(ch, "fig_pass.png", 'Sky track, elevation/range, Doppler, delay vs timing advance, ground track, SNR');
add(ch, rptTable(R.budget, 'Link budget (TR 38.821 style) at three elevations; last row: satelliteCNR cross-check', 20));
add(ch, img2(png("fig_throughput.png"), 'Throughput over the pass'));
add(rpt, ch);

ch = Chapter('Title', 'Link-level simulation (5G Toolbox)');
add(ch, Paragraph(sprintf(['PDSCH slots through nrTDLChannel (%s) with SatelliteDopplerShift = the CFO at the UE, AWGN, ' ...
    'CP-based fine CFO estimation, practical timing/channel estimation, MMSE equalisation and LDPC decoding. %d jobs, %.0f s%s.'], ...
    opt.Channel, R.simInfo.jobs, R.simInfo.seconds, ternary(isfield(R.simInfo, 'parallel') && R.simInfo.parallel, ' (parfor)', ''))));
img(ch, "fig_bler.png", 'BLER and throughput vs SNR per MCS');
img(ch, "fig_cfo.png", 'BLER vs carrier frequency offset: numerology and fine CFO estimation');
add(rpt, ch);

ch = Chapter('Title', 'Doppler and delay vs NR numerology');
add(ch, rptTable(R.numerology, 'Numerology table', 5));
add(ch, rptTable(R.challenges, 'NTN challenges, terrestrial NR limits and the Rel-17 answers', 10));
add(rpt, ch);

if isfield(R.interference, 'done') && R.interference.done
    ch = Chapter('Title', 'Coexistence with the measured spectrum');
    add(ch, Paragraph("Source: " + R.interference.source));
    if ~isempty(R.interference.notes), add(ch, UnorderedList(cellstr(R.interference.notes))); end
    img(ch, "fig_coexistence.png", 'Interferer PSD, victim and adjacent channels; SINR over the pass');
    add(ch, rptTable(R.interference.summary, 'Interference budget', 20));
    add(ch, rptTable(R.interference.throughput, 'Data volume with interference', 5));
    add(rpt, ch);
end

ch = Chapter('Title', 'Files');
add(ch, UnorderedList(cellstr(R.files)));
add(rpt, ch);
close(rpt);
file = string(rpt.OutputPath);
fprintf('Report: %s\n', file);
end

function addImg(ch, file, cap)
if isfile(file), add(ch, img2(file, cap)); end
end

function im = img2(file, cap)
import mlreportgen.report.*
if isfile(file)
    im = FormalImage('Image', file, 'Caption', cap, 'Width', '16cm');
else
    im = mlreportgen.dom.Paragraph(['(missing figure: ' char(cap) ')']);
end
end

function bt = rptTable(T, ttl, maxRows)
import mlreportgen.report.*
import mlreportgen.dom.*
T = T(1:min(end, maxRows), :);
for v = string(T.Properties.VariableNames)
    x = T.(v);
    if isdatetime(x), T.(v) = string(x, 'yyyy-MM-dd HH:mm:ss');
    elseif isduration(x), T.(v) = string(x);
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

%% ===== small helpers ============================================================================

function s = fmtHz(x)
if ~isfinite(x), s = "NaN"; return, end
a = abs(x);
if a >= 1e9, s = sprintf('%.4g GHz', x/1e9);
elseif a >= 1e6, s = sprintf('%.4g MHz', x/1e6);
elseif a >= 1e3, s = sprintf('%.4g kHz', x/1e3);
else, s = sprintf('%.4g Hz', x);
end
s = string(s);
end

function s = where(err)
% 'function line N' of the first stack frame inside this file (for warnings)
s = "";
for k = 1:numel(err.stack)
    if contains(err.stack(k).file, mfilename), s = sprintf('%s line %d', err.stack(k).name, err.stack(k).line); return, end
end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function s = stamp()
s = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
end
