function R = ff_cubesat_link(varargin)
%FF_CUBESAT_LINK  CubeSat pass plan, time-varying link budget and Doppler table for a Satlab SRS-4 downlink.
%
%   R = ff_cubesat_link                        KNACKSAT-2 (NORAD 67683) over Bangkok, next 24 h
%   R = ff_cubesat_link(Name, Value, ...)
%   R = ff_cubesat_link('Offline', true)       no internet: cached TLE, else support\tle\thai_cubesats.tle
%   R = ff_cubesat_link('Altitude', 500)       circular sun-synchronous orbit (no TLE needed)
%   R = ff_cubesat_link('Band', "UHF")         437 MHz amateur-band case (9k6 GMSK, estimate)
%
% Example (headless, your lab site, report and a FieldFox Doppler-tracking table):
%   R = ff_cubesat_link('Site', [13.7290 100.7750 5], 'SiteName', "Lab roof", ...
%        'FieldFoxTrack', true, 'Visible', 'off', 'Prompt', false);
%   R.passes, R.best, R.doppler.summary, R.budget
%
% SET 'Site' TO YOUR LAB: the default is central Bangkok (13.7563 N, 100.5018 E, 10 m).
% A few km moves pass times by seconds only, but the elevation mask, the horizon and the
% P.618 rain climate are site-specific.
%
% What it does (no instrument is used; only celestrak.org is contacted, for the TLE):
%  1. Orbit: TLE from CelesTrak (gp.php?CATNR=..&FORMAT=TLE, cached 2 h in support\tle,
%     no retry for 2 h after an error), a TLE file, or Keplerian elements written as a TLE
%     (SGP4, epoch = 'Start', no drag). TLE age is reported: old TLE = wrong pass times.
%  2. Survey: satelliteScenario + satellite (SGP4) + groundStation -> access, accessIntervals,
%     accessPercentage, states (ground track).
%  3. Every pass in its own 1 s scenario: aer (az, el, slant range), dopplershift (+ Doppler
%     rate), the CubeSat antenna (Antenna Toolbox patchMicrostrip / Gaussian proxy /
%     isotropic) with pointAt(sat, gs), a Gaussian dish on a gimbal that pointAt's the
%     satellite, transmitter/receiver/link -> ebno, sigstrength, txEIRP, linkIntervals,
%     linkPercentage. Atmosphere: ITU-R P.618 (p618Config + p618PropagationLosses: gas,
%     cloud, rain, scintillation) vs elevation if the maps in support\p618 exist, else a
%     fixed loss.
%  4. Required Eb/N0 for every SRS-4 configuration (BPSK/QPSK/8PSK x none/CC/RS/CC+RS x
%     symbol rate) from the SRS-4 datasheet Table 7 (PER < 1 %):
%         Eb/N0_req = Psens - (-174 + NF + 10 log10(Rb)),  NF = 2.05 dB,
%     Rb = net payload rate (datasheet throughput: 217-byte payload; RS(255,223); CC r1/2).
%     Between table rates: linear in log10(Rs). 8PSK (TX only, no table) = QPSK value + the
%     8PSK-vs-QPSK theory penalty: an ESTIMATE. Theory columns (berawgn, bercoding with
%     distspec(poly2trellis(7,[171 133])), RS via binocdf) give the implementation loss.
%  5. Margin vs time for every configuration; availability; data volume per pass and per
%     day for: your configuration, the best fixed one (day), the best per pass, the
%     fastest one that closes >= 95 % of the time, and an adaptive (per-second best feasible
%     configuration) upper bound. The SRS-4 has no ACM: adaptive needs the ground to
%     reconfigure the radio, so treat it as a bound. Strategies only use configurations that
%     fit 'MaxOBW', the ECSS OBW limit and 'BusRate' (column Allowed).
%  6. Static check at the highest elevation with satelliteCNRConfig / satelliteCNR (the
%     classic link-budget table) - it must agree with the time-varying result.
%  7. Doppler pre-compensation CSV (UTC, local, az, el, range, range rate, Doppler, rate,
%     receiver tuning, uplink pre-compensated TX frequency) and checks against the SRS-4
%     receiver: demodulator pull-in ~1 % Rs, AFC 15 kHz (datasheet) / 100 kHz (fw 2.0 for
%     Rs > 30 kBd), including the SRS-4 frequency error (2.5 ppm worst case over temperature).
%     Tolerance = how many seconds of timing error (or table step) the link survives.
%  8. Optional: FieldFox Doppler-tracking table (centre frequency vs time, span, RBW),
%     eclipse (aeroDataPackage ephemeris, in a SEPARATE scenario), satelliteScenarioViewer.
%  9. Figures (ground track geoplot, skyplot, link vs time, Doppler, data volume, required
%     Eb/N0), CSV/MAT files and a PDF report (MATLAB Report Generator).
%
% Output R: settings, orbit (TLE, age), site, survey, passes (one row per pass), pass
% (time series per pass), requirements (Eb/N0 per configuration), theory, user (your
% configuration), strategies (volume per pass), best, budget (static link budget table),
% static (satelliteCNR), doppler (table, summary, limits), fieldfox, atmosphere, antenna,
% eclipse, notes, files, outDir. Files in 'OutDir': passes.csv, doppler_table.csv,
% strategies.csv, link_budget.csv, requirements.csv, [fieldfox_track.csv], tle_used.tle,
% fig_*.png, cubesat_link.mat, cubesat_link_report.pdf. A 24 h run takes ~1-2 min (each
% pass is ~1-10 s, so parfor would not pay for starting a pool).
%
% Name-Value options (defaults; S band unless noted):
%  Orbit
%   'NORAD'          67683   KNACKSAT-2 (NAPA-1 = 46320). Fetched from CelesTrak.
%   'TLE'            ""      TLE file instead (the entry matching 'NORAD', else the first)
%   'Offline'        false   true = never go online: cached TLE (any age) or the sample file
%   'Altitude'       []      km: circular orbit from elements instead of a TLE
%   'Inclination'    []      deg; [] = sun-synchronous for 'Altitude'
%   'RAAN'           []      deg; [] = from 'LTAN'
%   'LTAN'           10.5    h, local time of the ascending node (SSO)
%   'TrueAnomaly'    0       deg at 'Start'
%  Site and time
%   'Site'           [13.7563 100.5018 10]  [lat lon alt_m]  SET YOUR LAB LOCATION
%   'SiteName'       "Bangkok"
%   'TimeZone'       "Asia/Bangkok"   local time in tables (UTC+7)
%   'Start'          now     datetime or text; without a zone it is read in 'TimeZone'
%   'Duration'       24      h (or a duration)
%   'MinElevation'   10      deg, elevation mask
%   'SampleTime'     10      s, survey step (pass search)
%   'PassSampleTime' 1       s, step inside each pass
%   'MaxPasses'      20      passes analysed in detail (chronological)
%  SRS-4 transmitter (CubeSat)
%   'Band'           "S"     or "UHF" (437 MHz, 1 W, 9600 Bd GMSK-like, isotropic, Yagi)
%   'Frequency'      2245e6  Hz (SRS-4 TX 2200-2290 MHz; exact channel unknown)
%   'TxPower_dBm'    33      dBm (SRS-4: 20-33 dBm). sl-radio --tx-gain=50 is a gain setting, not
%                            dBm: measure the output (FieldFox channel power through an attenuator,
%                            max +27 dBm at Port 2) and enter that value
%   'TxLoss'         1       dB, radio-to-antenna cable/connectors
%   'SymbolRate'     1e6     Bd, your setting (sl-radio --tx-rate=1000000)
%   'Modulation'     "QPSK"  "BPSK" | "QPSK" | "8PSK" (your setting)
%   'FEC'            "CC+RS" "none" | "CC" | "RS" | "CC+RS" (your setting)
%   'Rolloff'        0.2     RRC roll-off (0.2, 0.25, 0.35) -> occupied bandwidth
%   'SymbolRates'    [0.1 0.25 0.5 1 2 2.5 4 5]*1e6  candidates for the best configuration
%   'CubeSatAntenna' "patch" "patch" = design(patchMicrostrip, f) (about 10 dBi, linear),
%                            "gaussian" = wide-beam proxy of 'CubeSatGain' dBi, "isotropic"
%   'CubeSatGain'    6       dBi, peak gain of the "gaussian" proxy
%   'CubeSatPointing' "target"  pointAt(sat, gs) during passes; "nadir" = body-fixed nadir
%  Ground station receiver
%   'DishDiameter'   2.4     m, Gaussian dish on a gimbal (UHF: 1.48 m = 14 dBi Yagi proxy)
%   'DishEfficiency' 0.55
%   'Tsys'           200     K, system noise temperature (UHF: 700 K, urban noise)
%   'GT'             []      dB/K; [] = dish gain - 10 log10(Tsys)
%   'RxSystemLoss'   1       dB, pointing/radome/feed (demodulator loss is in Table 7)
%   'PolarizationLoss' 3     dB (linear patch -> circular dish; 0.3 for matched RHCP)
%   'Exceedance'     1       % of an average year for P.618 (1 % = 99 % availability)
%   'AtmLoss'        0.5     dB, fixed fallback when P.618 is not available (and at UHF)
%   'MarginTarget'   3       dB required above Eb/N0_req
%   'MaxOBW'         Inf     Hz, your licensed bandwidth (NBTC/ITU filing): wider configurations
%                            are excluded from the "best" strategies
%   'ECSSMask'       true    also apply the ECSS-E-ST-50-05C Table 5-5 OBW limit (2200-2290 MHz,
%                            e.g. 1.1 Rs above 2 MBd: roll-off 0.35 fails there)
%   'BusRate'        100e6   bit/s the satellite can feed the SRS-4 (Ethernet; RS-422 12.5e6,
%                            CAN 1e6): faster configurations are excluded
%   'NoiseFigure'    2.05    dB, SRS-4 RX NF used to turn Table 7 into Eb/N0
%  Doppler
%   'UplinkFrequency' 2067.5e6  Hz, SRS-4 receive channel (Table 7 frequency)
%   'UplinkSymbolRate' []    Bd; [] = 'SymbolRate' (sl-radio --rx-rate)
%   'FreqErrorPPM'   2.5     SRS-4 frequency error (0.5 initial, 2.5 max over temperature)
%   'UHFImplLoss'    3       dB added to coherent-BPSK theory for the UHF estimate
%   'FieldFoxTrack'  false   true = write fieldfox_track.csv (centre frequency vs time)
%   'TrackStep'      1       s, step of the FieldFox track
%  Output
%   'Eclipse'        false   true = eclipse/sunlight per pass (needs aeroDataPackage data)
%   'Viewer'         false   true = open satelliteScenarioViewer (needs a display/WebGL)
%   'Report'         "pdf"   "pdf", "docx" or "none"
%   'Notes'          ""      text for the report
%   'OutDir'         ""      default cubesat_link_<time> in the current folder
%   'Visible'        'on'    'off' = no windows (figures are still saved)
%   'Prompt'         true    ask to confirm the site when it is the default (ignored in -batch)
%
% Use the FieldFox track (with 'FieldFoxTrack', true) in a capture loop, e.g.
%   T = readtable(fullfile(R.outDir, "fieldfox_track.csv"));
%   ff.saConfigure('Center', T.Center_Hz(1), 'Span', R.fieldfox.trackSpan_Hz, 'RBW', R.fieldfox.rbw_Hz);
%   % then at each T.UTC: ff.write("SENS:FREQ:CENT " + T.Center_Hz(k)); P = ff.saTrace; ...
%
% Read this before trusting the numbers:
%   - It is a PREDICTION. The TLE drag term, the attitude (pointAt is perfect), antenna
%     patterns (the Antenna Toolbox patch has a 1-wavelength ground plane, larger than a 1U
%     face) and the ground G/T dominate the result: change them and look at the margin.
%   - Table 7 is the SRS-4 RECEIVER (uplink band, 2067.5 MHz). Using it for the downlink
%     assumes the ground demodulator is as good as the SRS-4 (implementation loss < 0.3 dB).
%     The margin itself does not depend on how Rb is defined: Eb/N0 - Eb/N0_req =
%     C/N0 - (Psens + 174 - NF).
%   - P.618 gives losses exceeded 'Exceedance' % of an average year, not today's weather;
%     the sky-noise increase (tsky) is reported, not added to Tsys. Below 1 GHz (UHF) P.618
%     does not apply: fixed loss, no ionospheric scintillation / Faraday rotation model.
%   - Doppler sign: positive while the satellite approaches (received f = f + Doppler).
%     A FieldFox has +-0.7 ppm reference error (+-1.6 kHz at 2245 MHz) and the SRS-4 up to
%     2.5 ppm: absolute frequency (ppm) claims need an external 10 MHz GPSDO; the Doppler
%     RATE and the time of closest approach do not.
%   - No I/Q on the FieldFox N9914A: over the air it can confirm power, occupied bandwidth,
%     symbol rate (-3 dB width) and Doppler, not the modulation order, EVM or BER.
%   - The UHF case is a generic 9k6 GMSK-like estimate, not the SRS-4.
%
% Toolboxes: Satellite Communications, Aerospace (tleread), Antenna, Phased Array,
% Communications (berawgn, bercoding), Statistics (binocdf), Navigation (skyplot), Mapping
% (geoplot basemap 'darkwater' - installed, no download; scircle1), Report Generator.

here = fileparts(mfilename('fullpath'));
opt = parseOptions(varargin);
if strlength(opt.OutDir) == 0, opt.OutDir = fullfile(pwd, "cubesat_link_" + stamp()); end
opt.OutDir = string(opt.OutDir);
if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
opt = askSite(opt);
notes = strings(0, 1);
if isequal(opt.Site, [13.7563 100.5018 10])
    notes(end+1) = "Site = default central Bangkok; set 'Site' to the lab [lat lon alt_m].";
end
t0run = tic;

%% 1. Orbit (TLE) and SRS-4 requirements
[orb, notes] = getOrbit(opt, here, notes);
fprintf('Orbit: %s (NORAD %d), %s; TLE epoch %s UTC (age %.1f days at Start)\n', orb.name, orb.norad, ...
    orb.source, string(orb.epoch, 'yyyy-MM-dd HH:mm'), orb.ageDays);
[REQ, theo] = srs4Requirements(opt);
u = find(REQ.Modulation == opt.Modulation & REQ.FEC == opt.FEC & abs(REQ.SymbolRate_Bd - opt.SymbolRate) < 1, 1);
user = table2struct(REQ(u, :));
fprintf('Your configuration: %s -> payload %.3g kbit/s, Eb/N0_req %.2f dB (%s)\n', REQ.Config(u), ...
    REQ.Rb_bps(u)/1e3, REQ.EbN0_req_dB(u), REQ.Source(u));

%% 2. Survey: passes over 'Duration'
S = survey(opt, orb);
fprintf('Survey %s .. %s UTC: %d passes above %g deg (in view %.2f %% of the time)\n', ...
    string(opt.Start, 'yyyy-MM-dd HH:mm'), string(opt.Start + opt.Duration, 'yyyy-MM-dd HH:mm'), ...
    height(S.access), opt.MinElevation, S.accessPct);

%% 3. Atmosphere and antennas
[atm, notes] = atmosphere(opt, here, notes);
[ant, notes] = cubesatAntenna(opt, notes);

%% 4. Passes in detail
nP = min(height(S.access), opt.MaxPasses);
if height(S.access) > opt.MaxPasses
    notes(end+1) = sprintf('%d passes found; only the first %d (MaxPasses) are analysed, so per-day volumes are low.', ...
        height(S.access), opt.MaxPasses);
end
P = [];
for k = 1:nP
    fprintf('Pass %d/%d: %s UTC ... ', k, nP, string(S.access.StartTime(k), 'MM-dd HH:mm:ss'));
    tp = tic;
    p = analysePass(numel(P) + 1, S.access.StartTime(k), S.access.EndTime(k), opt, orb, ant, atm, user);
    if isempty(p)
        fprintf('fewer than 2 samples above the mask, skipped\n');
    else
        P = [P, p]; %#ok<AGROW>
        fprintf('max el %.1f deg, %.0f s (%.1f s)\n', p.maxEl, seconds(p.los - p.aos), toc(tp));
    end
end

%% 5. Link evaluation, static check, Doppler, FieldFox plan, eclipse
[P, E] = evaluate(P, REQ, u, opt);
try
    [budget, stat] = staticCheck(P, REQ(u, :), opt, ant);
catch err
    budget = table(); stat = struct('ok', false);
    notes(end+1) = "Static satelliteCNR check failed: " + err.message;
end
D = dopplerAnalysis(P, opt);
FF = fieldfoxPlan(P, opt);
EC = struct('ok', false, 'note', "not requested ('Eclipse', true)", 'intervals', table(), 'fraction', NaN, ...
    'passInEclipse', []);
if opt.Eclipse, EC = eclipseAnalysis(opt, orb, P); end
if opt.Eclipse && ~EC.ok, notes(end+1) = EC.note; end
passes = passTable(P, E, D, EC, opt);

%% 6. Results, files, figures, report
R.settings = opt; R.orbit = orb; R.site = struct('Name', opt.SiteName, 'LatLonAlt', opt.Site);
R.survey = rmfield(S, {'sc', 'sat', 'gs'}); R.passes = passes; R.pass = P;
R.requirements = REQ; R.theory = theo; R.user = user; R.strategies = E.perPass; R.best = E.best;
R.budget = budget; R.static = stat; R.doppler = D; R.fieldfox = FF; R.atmosphere = atm;
R.eclipse = EC; R.antenna = rmfield(ant, 'obj');
R.notes = [notes(:); atm.notes(:); E.notes(:); D.notes(:)];
R.outDir = opt.OutDir;
files = writeFiles(R, opt);
files = [files; makeFigures(R, S, opt)];
R.files = files;
matFile = fullfile(opt.OutDir, "cubesat_link.mat");
save(matFile, 'R');
R.files = [matFile; R.files];
if opt.Report ~= "none"
    try
        rep = makeReport(R, opt);
        R.files = [rep; R.files];
    catch err
        warning('ff_cubesat_link:report', 'Report not written: %s', err.message);
    end
end
if opt.Viewer
    try
        R.viewer = satelliteScenarioViewer(S.sc, 'ShowDetails', true);
    catch err
        warning('ff_cubesat_link:viewer', 'Viewer not opened: %s', err.message);
    end
end
printSummary(R, opt);
fprintf('Done in %.0f s. Results in %s\n', toc(t0run), opt.OutDir);
if nargout == 0, clear R, end
end

%% ===== options ===============================================================

function opt = parseOptions(args)
opt = struct('NORAD', 67683, 'TLE', "", 'Offline', false, 'Altitude', [], 'Inclination', [], 'RAAN', [], ...
    'LTAN', 10.5, 'TrueAnomaly', 0, 'Site', [13.7563 100.5018 10], 'SiteName', "Bangkok", ...
    'TimeZone', "Asia/Bangkok", 'Start', [], 'Duration', 24, 'MinElevation', 10, 'SampleTime', 10, ...
    'PassSampleTime', 1, 'MaxPasses', 20, 'Band', "S", 'Frequency', [], 'TxPower_dBm', [], 'TxLoss', 1, ...
    'SymbolRate', [], 'Modulation', [], 'FEC', [], 'Rolloff', 0.2, 'SymbolRates', [], 'CubeSatAntenna', [], ...
    'CubeSatGain', 6, 'CubeSatPointing', "target", 'DishDiameter', [], 'DishEfficiency', 0.55, 'Tsys', [], ...
    'GT', [], 'RxSystemLoss', 1, 'PolarizationLoss', 3, 'Exceedance', 1, 'AtmLoss', 0.5, 'MarginTarget', 3, ...
    'NoiseFigure', 2.05, 'UplinkFrequency', 2067.5e6, 'UplinkSymbolRate', [], 'FreqErrorPPM', 2.5, ...
    'UHFImplLoss', 3, 'FieldFoxTrack', false, 'TrackStep', 1, 'Eclipse', false, 'Viewer', false, ...
    'MaxOBW', Inf, 'ECSSMask', true, 'BusRate', 100e6, ...
    'Report', "pdf", 'Notes', "", 'OutDir', "", 'Visible', 'on', 'Prompt', true);
if mod(numel(args), 2), error('ff_cubesat_link:args', 'Options must be Name, Value pairs.'); end
names = fieldnames(opt);
for k = 1:2:numel(args)
    j = find(strcmpi(names, args{k}), 1);
    if isempty(j), error('ff_cubesat_link:args', 'Unknown option "%s".', string(args{k})); end
    opt.(names{j}) = args{k + 1};
end
opt.Band = upper(string(opt.Band));
if ~ismember(opt.Band, ["S" "UHF"]), error('ff_cubesat_link:args', 'Band must be "S" or "UHF".'); end
if opt.Band == "S"
    dflt = struct('Frequency', 2245e6, 'TxPower_dBm', 33, 'SymbolRate', 1e6, 'Modulation', "QPSK", ...
        'FEC', "CC+RS", 'SymbolRates', [0.1 0.25 0.5 1 2 2.5 4 5]*1e6, 'CubeSatAntenna', "patch", ...
        'DishDiameter', 2.4, 'Tsys', 200);
else
    dflt = struct('Frequency', 437e6, 'TxPower_dBm', 30, 'SymbolRate', 9600, 'Modulation', "GMSK", ...
        'FEC', "none", 'SymbolRates', [1200 2400 4800 9600 19200], 'CubeSatAntenna', "isotropic", ...
        'DishDiameter', 1.48, 'Tsys', 700);
end
for f = string(fieldnames(dflt)).'
    if isempty(opt.(f)), opt.(f) = dflt.(f); end
end
if isempty(opt.UplinkSymbolRate), opt.UplinkSymbolRate = opt.SymbolRate; end
opt.Modulation = upper(string(opt.Modulation));
opt.FEC = upper(string(opt.FEC)); if opt.FEC == "NONE", opt.FEC = "none"; end
opt.CubeSatAntenna = lower(string(opt.CubeSatAntenna));
opt.CubeSatPointing = lower(string(opt.CubeSatPointing));
opt.Report = lower(string(opt.Report));
opt.TLE = string(opt.TLE); opt.SiteName = string(opt.SiteName); opt.TimeZone = string(opt.TimeZone);
opt.Notes = string(opt.Notes);
opt.Visible = char(matlab.lang.OnOffSwitchState(opt.Visible));
opt.Prompt = logical(opt.Prompt); opt.Offline = logical(opt.Offline);
opt.FieldFoxTrack = logical(opt.FieldFoxTrack); opt.Eclipse = logical(opt.Eclipse); opt.Viewer = logical(opt.Viewer);
if ~ismember(opt.CubeSatAntenna, ["patch" "gaussian" "isotropic"])
    error('ff_cubesat_link:args', 'CubeSatAntenna must be "patch", "gaussian" or "isotropic".');
end
if ~ismember(opt.CubeSatPointing, ["target" "nadir"]), error('ff_cubesat_link:args', 'CubeSatPointing must be "target" or "nadir".'); end
if ~ismember(opt.Report, ["pdf" "docx" "none"]), error('ff_cubesat_link:args', 'Report must be "pdf", "docx" or "none".'); end
if opt.Band == "S"
    if ~ismember(opt.Modulation, ["BPSK" "QPSK" "8PSK"]), error('ff_cubesat_link:args', 'Modulation must be BPSK, QPSK or 8PSK.'); end
    if ~ismember(opt.FEC, ["none" "CC" "RS" "CC+RS"]), error('ff_cubesat_link:args', 'FEC must be "none", "CC", "RS" or "CC+RS".'); end
    rs = unique([opt.SymbolRates(:); opt.SymbolRate]);
    if any(rs < 10e3 | rs > 5e6), error('ff_cubesat_link:args', 'SRS-4 symbol rates are 10 kBd (fw 1.7+) to 5 MBd.'); end
    if ~ismember(opt.Rolloff, [0.2 0.25 0.35]), error('ff_cubesat_link:args', 'Rolloff must be 0.2, 0.25 or 0.35.'); end
    if opt.Frequency < 2200e6 || opt.Frequency > 2290e6
        warning('ff_cubesat_link:band', 'Frequency %.1f MHz is outside the SRS-4 TX band 2200-2290 MHz.', opt.Frequency/1e6);
    end
    if opt.TxPower_dBm < 20 || opt.TxPower_dBm > 33
        warning('ff_cubesat_link:power', 'SRS-4 TX power is 20-33 dBm; %.1f dBm used anyway.', opt.TxPower_dBm);
    end
else
    opt.Modulation = "GMSK"; opt.FEC = "none";
end
site = double(opt.Site(:)).';
if numel(site) == 2, site(3) = 0; end
if numel(site) ~= 3 || abs(site(1)) > 90, error('ff_cubesat_link:args', 'Site must be [lat lon alt_m].'); end
opt.Site = site;
if isempty(opt.Start)
    opt.Start = dateshift(datetime('now', 'TimeZone', 'UTC'), 'start', 'minute');
else
    st = opt.Start;
    if ~isdatetime(st), st = datetime(string(st)); end
    if isempty(st.TimeZone), st.TimeZone = opt.TimeZone; end
    st.TimeZone = 'UTC'; opt.Start = st;
end
if isnumeric(opt.Duration), opt.Duration = hours(opt.Duration); end
if opt.Duration < seconds(2*opt.SampleTime), error('ff_cubesat_link:args', 'Duration must be at least 2 x SampleTime.'); end
if opt.PassSampleTime > opt.SampleTime, error('ff_cubesat_link:args', 'PassSampleTime must not exceed SampleTime.'); end
if opt.MinElevation < 0 || opt.MinElevation >= 90, error('ff_cubesat_link:args', 'MinElevation must be 0-90 deg.'); end
if ~isempty(opt.Altitude) && (opt.Altitude < 150 || opt.Altitude > 5000)
    error('ff_cubesat_link:args', 'Altitude must be 150-5000 km (SGP4 near-Earth model).');
end
if isempty(opt.GT)
    lam = physconst('LightSpeed')/opt.Frequency;
    g = 10*log10(opt.DishEfficiency*(pi*opt.DishDiameter/lam)^2);
    opt.GT = g - 10*log10(opt.Tsys);
    opt.GTSource = sprintf('%.2f m dish, %.0f %% efficiency (%.1f dBi), Tsys %.0f K', opt.DishDiameter, ...
        100*opt.DishEfficiency, g, opt.Tsys);
else
    opt.GTSource = "given ('GT')";
end
end

function opt = askSite(opt)
if ~opt.Prompt || batchStartupOptionUsed || ~usejava('desktop'), return; end
if ~isequal(opt.Site, [13.7563 100.5018 10]), return; end
fprintf('Ground station: %s (%.4f, %.4f, %g m) - the DEFAULT, not your lab.\n', opt.SiteName, opt.Site);
s = input('Press Enter to keep it, or type [lat lon alt_m]: ', 's');
v = sscanf(regexprep(s, '[\[\],;]', ' '), '%f').';
if numel(v) >= 2
    if numel(v) == 2, v(3) = 0; end
    opt.Site = v(1:3); opt.SiteName = "Site";
end
end

%% ===== orbit =================================================================

function [orb, notes] = getOrbit(opt, here, notes)
tleDir = fullfile(here, "support", "tle");
sample = fullfile(tleDir, "thai_cubesats.tle");
if ~isempty(opt.Altitude)
    [name, l1, l2, orb] = keplerTLE(opt);
    orb.source = sprintf('Keplerian elements (%.0f km, i = %.2f deg, RAAN %.2f deg), SGP4 without drag', ...
        opt.Altitude, orb.inc, orb.raan);
elseif strlength(opt.TLE) > 0
    if ~isfile(opt.TLE), error('ff_cubesat_link:tle', 'TLE file not found: %s', opt.TLE); end
    [name, l1, l2] = pickTLE(readlines(opt.TLE), opt.NORAD, true);
    orb.source = "TLE file " + opt.TLE;
else
    if opt.NORAD > 99999, error('ff_cubesat_link:tle', 'NORAD > 99999 cannot be a TLE: give an OMM file via satellite() yourself.'); end
    if ~isfolder(tleDir), mkdir(tleDir); end
    cache = fullfile(tleDir, sprintf('%d.tle', opt.NORAD));
    [name, l1, l2, orb.source, notes] = celestrak(opt, cache, sample, notes);
end
orb.name = string(strtrim(name)); orb.norad = str2double(l1(3:7));
orb.file = fullfile(opt.OutDir, "tle_used.tle");
writelines([orb.name; string(l1); string(l2)], orb.file);
tl = tleread(orb.file);
orb.epoch = tl(1).Epoch; if isempty(orb.epoch.TimeZone), orb.epoch.TimeZone = 'UTC'; end
orb.ageDays = days(opt.Start - orb.epoch);
orb.inclination = tl(1).Inclination; orb.meanMotion_revday = tl(1).MeanMotion*86400/360;   % tleread: deg/s
orb.bstar = tl(1).BStar; orb.lines = [orb.name; string(l1); string(l2)];
if abs(orb.ageDays) > 14
    notes(end+1) = sprintf(['TLE is %.0f days from Start: SGP4 errors grow to tens of km (minutes in pass time). ' ...
        'Use a fresh TLE or a Start near the epoch.'], orb.ageDays);
elseif abs(orb.ageDays) > 3
    notes(end+1) = sprintf('TLE is %.1f days old at Start: expect several seconds to a minute of pass-time error for a low, draggy CubeSat.', orb.ageDays);
end
end

function [name, l1, l2, src, notes] = celestrak(opt, cache, sample, notes)
age = fileAgeHours(cache);
errFile = cache + ".err";
if opt.Offline
    if isfile(cache)
        [name, l1, l2] = pickTLE(readlines(cache), opt.NORAD, false);
        src = sprintf('cached CelesTrak TLE (%.1f h old, offline)', age); return
    end
    [name, l1, l2] = pickTLE(readlines(sample), opt.NORAD, false);
    src = "sample TLE " + sample + " (offline)"; return
end
if isfile(cache) && age < 2
    [name, l1, l2] = pickTLE(readlines(cache), opt.NORAD, false);
    src = sprintf('CelesTrak, cached %.0f min ago (refresh every 2 h)', age*60); return
end
if isfile(errFile) && fileAgeHours(errFile) < 2
    notes(end+1) = "CelesTrak failed less than 2 h ago (" + strtrim(fileread(errFile)) + "); not retried.";
else
    url = sprintf('https://celestrak.org/NORAD/elements/gp.php?CATNR=%d&FORMAT=TLE', opt.NORAD);
    try
        fprintf('Fetching %s\n', url);
        txt = webread(url, weboptions('ContentType', 'text', 'Timeout', 20));
        [name, l1, l2] = pickTLE(splitlines(string(txt)), opt.NORAD, false);
        writelines([string(strtrim(name)); string(l1); string(l2)], cache);
        if isfile(errFile), delete(errFile); end
        src = "CelesTrak (downloaded now)"; return
    catch err
        msg = string(err.message); msg = extractBefore(msg + newline, newline);
        writelines(string(datetime('now')) + " " + msg, errFile);   % no retry for 2 h (CelesTrak etiquette)
        notes(end+1) = "CelesTrak download failed: " + msg;
    end
end
if isfile(cache)
    [name, l1, l2] = pickTLE(readlines(cache), opt.NORAD, false);
    src = sprintf('cached CelesTrak TLE, %.1f h old (download failed)', fileAgeHours(cache));
else
    try
        [name, l1, l2] = pickTLE(readlines(sample), opt.NORAD, false);
    catch
        error('ff_cubesat_link:tle', ['No TLE for NORAD %d: CelesTrak gave none (%s) and the sample file has ' ...
            'none either. Check the number, or give ''TLE'' or ''Altitude''.'], opt.NORAD, notes(end));
    end
    src = "sample TLE " + sample + " (download failed)";
end
end

function h = fileAgeHours(file)
h = Inf;
d = dir(file);
if ~isempty(d), h = hours(datetime('now') - datetime(d.datenum, 'ConvertFrom', 'datenum')); end
end

function [name, l1, l2] = pickTLE(lines, norad, firstIfMissing)
% Entries: optional name line, line 1, line 2 (checksums verified).
lines = strtrim(string(lines)); lines = lines(strlength(lines) > 0);
E = struct('name', {}, 'l1', {}, 'l2', {});
for k = 1:numel(lines) - 1
    a = char(lines(k)); b = char(lines(k + 1));
    if numel(a) >= 69 && numel(b) >= 69 && a(1) == '1' && b(1) == '2' && strcmp(a(3:7), b(3:7))
        nm = "SAT " + string(a(3:7));
        if k > 1 && ~startsWith(lines(k - 1), ["1 " "2 "]), nm = lines(k - 1); end
        if tleChecksum(a) ~= a(69) || tleChecksum(b) ~= b(69)
            warning('ff_cubesat_link:checksum', 'TLE checksum mismatch for %s.', nm);
        end
        E(end+1) = struct('name', regexprep(nm, '^0 ', ''), 'l1', a(1:69), 'l2', b(1:69)); %#ok<AGROW>
    end
end
if isempty(E), error('ff_cubesat_link:tle', 'No TLE found (reply: "%s").', extractBefore(strjoin(lines, " ") + "  ", 60)); end
j = find(arrayfun(@(e) str2double(e.l1(3:7)) == norad, E), 1);
if isempty(j)
    if ~firstIfMissing, error('ff_cubesat_link:tle', 'NORAD %d is not in the TLE data; give ''TLE'' or another ''NORAD''.', norad); end
    j = 1;
end
name = E(j).name; l1 = E(j).l1; l2 = E(j).l2;
end

function c = tleChecksum(line)
s = line(1:68);
c = char('0' + mod(sum(s(isstrprop(s, 'digit')) - '0') + sum(s == '-'), 10));
end

function [name, l1, l2, orb] = keplerTLE(opt)
% Circular orbit as a TLE (epoch = Start) so that every scenario propagates the same SGP4 orbit.
mu = 3.986004418e14; Re = 6378137; J2 = 1.08262668e-3;
a = Re + opt.Altitude*1e3;
n = sqrt(mu/a^3);
inc = opt.Inclination;
if isempty(inc)                                   % sun-synchronous: RAAN drift = 360 deg / year
    wSun = 2*pi/(365.2422*86400);
    inc = acosd(-2*wSun*a^3.5/(3*J2*Re^2*sqrt(mu)));
end
raan = opt.RAAN;
if isempty(raan), raan = mod(sunRA(opt.Start) + 15*(opt.LTAN - 12), 360); end
t = opt.Start;
doy = day(t, 'dayofyear') + hours(timeofday(t))/24;
l1 = sprintf('1 99999U %-8s %02d%012.8f %10s %8s %8s 0 %4d', '26999A', mod(year(t), 100), doy, ...
    ' .00000000', ' 00000-0', ' 00000-0', 999);
l2 = sprintf('2 99999 %8.4f %8.4f %07d %8.4f %8.4f %11.8f%5d', inc, raan, 0, 0, mod(opt.TrueAnomaly, 360), ...
    n*86400/(2*pi), 1);
l1 = [l1 tleChecksum(l1)]; l2 = [l2 tleChecksum(l2)];
name = sprintf('KEPLER-%.0fKM', opt.Altitude);
orb = struct('inc', inc, 'raan', raan);
end

function ra = sunRA(t)
% Low-precision solar right ascension (deg), good to ~0.01 deg (Astronomical Almanac).
n = juliandate(t) - 2451545.0;
L = 280.460 + 0.9856474*n; g = 357.528 + 0.9856003*n;
lam = L + 1.915*sind(g) + 0.020*sind(2*g); ep = 23.439 - 4e-7*n;
ra = mod(atan2d(cosd(ep)*sind(lam), cosd(lam)), 360);
end

%% ===== SRS-4 requirements =====================================================

function [REQ, theo] = srs4Requirements(opt)
NF = opt.NoiseFigure;
fecs = ["none" "CC" "RS" "CC+RS"];
eff = [217/227, 217/227/2, 217/259, 217/259/2];        % payload / channel bits (datasheet throughput)
frameBits = 8*223;                                      % header + payload + CRC under one RS block
theo = theoryTable(fecs, eff, frameBits, opt);
obwFac = interp1([0.2 0.25 0.35], [1.072 1.098 1.163], opt.Rolloff);
if opt.Band == "UHF"
    rs = unique([opt.SymbolRates(:); opt.SymbolRate]).';
    e = theo.EbN0_dB(theo.Modulation == "BPSK" & theo.FEC == "none") + opt.UHFImplLoss;
    n = numel(rs);
    REQ = table(compose("GMSK %g Bd", rs(:)), repmat("GMSK", n, 1), repmat("none", n, 1), rs(:), ones(n, 1), ...
        ones(n, 1), rs(:), repmat(e, n, 1), e - 174 + NF + 10*log10(rs(:)), repmat(e - opt.UHFImplLoss, n, 1), ...
        repmat(opt.UHFImplLoss, n, 1), 1.0*rs(:), repmat(opt.MaxOBW, n, 1), rs(:) <= opt.BusRate & rs(:) <= opt.MaxOBW, ...
        repmat(sprintf("estimate: coherent BPSK theory + %g dB", opt.UHFImplLoss), n, 1), 'VariableNames', reqNames());
    return
end
RsTab = [100e3 1e6 4e6 5e6];
P7.BPSK = [-112 -102 -96 -95; -121 -111 -105 -104; -114 -104 -98 NaN; -122 -112.5 -106 -105.5];   % none; CC; RS; CC+RS
P7.QPSK = [-109 -99 -93 -92; -118 -108 -101.5 -100.5; -111 -101 NaN NaN; -119 -108.5 -102 NaN];
d8 = theo.EbN0_dB(theo.Modulation == "8PSK" & theo.FEC == "none") - theo.EbN0_dB(theo.Modulation == "QPSK" & theo.FEC == "none");
rs = unique([opt.SymbolRates(:); opt.SymbolRate]).';
mods = ["BPSK" "QPSK" "8PSK"]; bps = [1 2 3];
C = {};
for m = 1:3
    for f = 1:4
        for r = rs
            Rb = r*bps(m)*eff(f);
            if mods(m) == "8PSK"                    % SRS-4 RX has no 8PSK: no Table 7 data
                e = fromTable(P7.QPSK(f, :) + 174 - NF - 10*log10(RsTab*2*eff(f)), RsTab, r) + d8;
                s = sprintf("estimate: QPSK + %.1f dB 8PSK theory penalty (no Table 7 data)", d8);
                th = theo.EbN0_dB(theo.Modulation == "8PSK" & theo.FEC == fecs(f));
                if isempty(th), th = NaN; end
            else
                [e, s] = fromTable(P7.(mods(m))(f, :) + 174 - NF - 10*log10(RsTab*bps(m)*eff(f)), RsTab, r);
                th = theo.EbN0_dB(theo.Modulation == mods(m) & theo.FEC == fecs(f));
            end
            lbl = sprintf('%s %s %s', mods(m), fecs(f), rateLabel(r));
            obwMax = opt.MaxOBW;
            if opt.ECSSMask, obwMax = min(obwMax, ecssObwLimit(r)); end
            C(end+1, :) = {string(lbl), mods(m), fecs(f), r, bps(m), eff(f), Rb, e, e - 174 + NF + 10*log10(Rb), ...
                th, e - th, obwFac*r, obwMax, obwFac*r <= obwMax && Rb <= opt.BusRate, string(s)}; %#ok<AGROW>
        end
    end
end
REQ = cell2table(C, 'VariableNames', reqNames());
end

function n = reqNames()
n = {'Config', 'Modulation', 'FEC', 'SymbolRate_Bd', 'BitsPerSymbol', 'Efficiency', 'Rb_bps', 'EbN0_req_dB', ...
    'Psens_equiv_dBm', 'Theory_EbN0_dB', 'ImplLoss_dB', 'OBW99_Hz', 'OBWmax_Hz', 'Allowed', 'Source'};
end

function L = ecssObwLimit(rs)
% ECSS-E-ST-50-05C Rev.2 Table 5-5: max OBW, category-A telemetry, 2200-2290 MHz (from SFCG Rec 21-2).
if rs < 10e3, L = 300e3;
elseif rs < 60e3, L = min(1200e3, 30*rs);
elseif rs < 2e6, L = min(max(1200e3, 12*rs), 6e6);
else, L = min(1.1*rs, 6e6);
end
end

function [v, src] = fromTable(E, RsTab, r)
ok = isfinite(E); x = log10(RsTab(ok)); y = E(ok); xq = log10(r);
j = find(abs(RsTab(ok) - r) < 1, 1);
if ~isempty(j)
    v = y(j); src = "Table 7";
elseif xq > x(1) && xq < x(end)
    v = interp1(x, y, xq); src = "Table 7, interpolated in log(Rs)";
elseif xq < x(1)
    v = y(1); src = "held from 100 kBd (below Table 7)";
else
    v = y(end); hi = RsTab(ok); hi = hi(end);
    if hi < 5e6
        src = sprintf("held from %s (Table 7 N/A: SRS-4 RS decoder limit)", rateLabel(hi));
    else
        src = "held from 5 MBd";
    end
end
end

function s = rateLabel(r)
if r >= 1e6, s = sprintf('%g MBd', r/1e6); else, s = sprintf('%g kBd', r/1e3); end
end

function T = theoryTable(fecs, eff, N, opt)
% Ideal coherent AWGN theory, Eb/N0 per payload bit for PER = 1 % (frame of N bits).
spect = distspec(poly2trellis(7, [171 133]), 7);
ccBer = @(x) min(bercoding(x, 'conv', 'soft', 1/2, spect), 0.5);
perU = @(pb) -expm1(N*log1p(-min(pb, 0.5)));
perRS = @(pb) binocdf(16, 255, -expm1(8*log1p(-min(pb, 0.5))), 'upper');
method = ["berawgn, independent bit errors over the frame"; ...
    "bercoding soft-decision union bound (an upper bound on BER: ~0.3-0.5 dB pessimistic)"; ...
    "berawgn -> RS(255,223) t = 16 symbol errors, binocdf (no frame-sync loss)"; ...
    "bercoding bound -> RS t = 16, ideal interleaving"];
rows = {};
for M = [2 4 8]
    mname = ["BPSK" "QPSK" "8PSK"]; mname = mname(log2(M));
    pb = @(x) berawgn(x, 'psk', M, 'nondiff');
    for f = 1:4
        if M == 8 && f > 1, continue, end          % coded 8PSK: no closed form used here
        k = 10*log10(eff(f));                      % Eb_channel(or CC input)/Eb_payload
        switch fecs(f)
            case "none",  per = @(x) perU(pb(x + k));
            case "RS",    per = @(x) perRS(pb(x + k));
            case "CC",    per = @(x) perU(ccBer(x + 10*log10(217/227)));
            case "CC+RS", per = @(x) perRS(ccBer(x + 10*log10(217/259)));
        end
        g = @(x) log10(max(per(x), 1e-300)) + 2;
        rows(end+1, :) = {mname, fecs(f), fzero(g, [-4 25]), method(f)}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'Modulation', 'FEC', 'EbN0_dB', 'Method'});
if opt.Band == "UHF", T = T(T.FEC == "none", :); end
end

%% ===== scenarios ===============================================================

function S = survey(opt, orb)
sc = satelliteScenario(opt.Start, opt.Start + opt.Duration, opt.SampleTime);
sat = satellite(sc, orb.file, "OrbitPropagator", "sgp4");
gs = groundStation(sc, opt.Site(1), opt.Site(2), "Altitude", opt.Site(3), "Name", opt.SiteName, ...
    "MinElevationAngle", opt.MinElevation);
ac = access(sat, gs);
T = accessIntervals(ac);
if height(T) > 0
    if isempty(T.StartTime.TimeZone), T.StartTime.TimeZone = 'UTC'; T.EndTime.TimeZone = 'UTC'; end
end
[pos, ~, tt] = states(sat, "CoordinateFrame", "geographic");
if isempty(tt.TimeZone), tt.TimeZone = 'UTC'; end
S = struct('access', T, 'accessPct', accessPercentage(ac), 't', tt, 'lat', pos(1, :), 'lon', pos(2, :), ...
    'alt', pos(3, :), 'elements', orbitalElements(sat), 'sc', sc, 'sat', sat, 'gs', gs);
end

function [atm, notes] = atmosphere(opt, here, notes)
atm = struct('method', sprintf('fixed %.2f dB (AtmLoss)', opt.AtmLoss), 'el', [0 90], 'At', opt.AtmLoss*[1 1], ...
    'Ag', NaN, 'Ac', NaN, 'Ar', NaN, 'As', NaN, 'tsky', NaN, 'notes', strings(0, 1));
mapDir = fullfile(here, "support", "p618");
if opt.Frequency < 1e9
    atm.notes(end+1) = "UHF: ITU-R P.618 covers 1-55 GHz, fixed AtmLoss used (ionospheric effects not modelled).";
    return
end
if isfolder(mapDir), addpath(mapDir); end
if isempty(which('iturp618v14Maps.mat'))
    atm.notes(end+1) = "ITU-R P.618 maps not installed (support\p618): fixed AtmLoss used.";
    return
end
cfg = p618Config;
cfg.Frequency = opt.Frequency; cfg.Latitude = opt.Site(1); cfg.Longitude = opt.Site(2);
cfg.TotalAnnualExceedance = opt.Exceedance; cfg.GasAnnualExceedance = opt.Exceedance;
cfg.CloudAnnualExceedance = opt.Exceedance; cfg.RainAnnualExceedance = opt.Exceedance;
cfg.ScintillationAnnualExceedance = opt.Exceedance;
cfg.AntennaDiameter = opt.DishDiameter; cfg.AntennaEfficiency = opt.DishEfficiency;
elg = unique([5 6 7 8 10 12.5 15 20 25 30 40 50 60 75 90 max(5, opt.MinElevation)]);
w = warning('off', 'satcom:p618PropagationLosses:InvalidXPDFrequency');   % XPD is not used (< 4 GHz)
c = onCleanup(@() warning(w));
try
    L = zeros(numel(elg), 6);
    for k = 1:numel(elg)
        cfg.ElevationAngle = elg(k);
        [pl, ~, tsky] = p618PropagationLosses(cfg);
        L(k, :) = [pl.At pl.Ag pl.Ac pl.Ar pl.As tsky];
    end
    atm.method = sprintf('ITU-R P.618-14 (p618PropagationLosses), exceeded %g %% of an average year', opt.Exceedance);
    atm.el = elg; atm.At = L(:, 1).'; atm.Ag = L(:, 2).'; atm.Ac = L(:, 3).'; atm.Ar = L(:, 4).';
    atm.As = L(:, 5).'; atm.tsky = L(:, 6).';
    if opt.MinElevation < 5, atm.notes(end+1) = "P.618 values below 5 deg elevation are held at the 5 deg value."; end
catch err
    atm.notes(end+1) = "P.618 failed (" + err.message + "): fixed AtmLoss used.";
end
end

function L = atmLoss(atm, el)
L = interp1(atm.el, atm.At, min(max(el, atm.el(1)), atm.el(end)), 'pchip');
end

function [ant, notes] = cubesatAntenna(opt, notes)
lam = physconst('LightSpeed')/opt.Frequency;
ant = struct('kind', opt.CubeSatAntenna, 'obj', [], 'D', NaN, 'eff', 0.6, 'peakGain_dBi', NaN, 'text', "");
if ant.kind == "patch" && opt.Band == "UHF"
    notes(end+1) = "A 437 MHz patch (34 cm) does not fit a CubeSat: isotropic antenna used instead.";
    ant.kind = "isotropic";
end
switch ant.kind
    case "patch"
        ant.obj = design(patchMicrostrip, opt.Frequency);
        ant.peakGain_dBi = pattern(ant.obj, opt.Frequency, 0, 90);
        ant.text = sprintf('Antenna Toolbox design(patchMicrostrip, %.1f MHz): %.1f dBi peak, linear polarisation', ...
            opt.Frequency/1e6, ant.peakGain_dBi);
    case "gaussian"
        ant.D = lam/pi*sqrt(10^(opt.CubeSatGain/10)/ant.eff);
        ant.peakGain_dBi = opt.CubeSatGain;
        ant.text = sprintf('Gaussian proxy: %.1f dBi peak (%.3f m aperture, 60 %%)', opt.CubeSatGain, ant.D);
    case "isotropic"
        ant.obj = phased.IsotropicAntennaElement;
        ant.peakGain_dBi = 0;
        ant.text = "isotropic (0 dBi; phased.IsotropicAntennaElement)";
end
end

function p = analysePass(k, aos, los, opt, orb, ant, atm, user)
pad = seconds(2*opt.SampleTime);
sc = satelliteScenario(aos - pad, los + pad, opt.PassSampleTime);
sat = satellite(sc, orb.file, "OrbitPropagator", "sgp4");
gs = groundStation(sc, opt.Site(1), opt.Site(2), "Altitude", opt.Site(3), "Name", opt.SiteName, ...
    "MinElevationAngle", opt.MinElevation);
if opt.CubeSatPointing == "target", pointAt(sat, gs); end
txArgs = {"Name", "SRS-4", "Frequency", opt.Frequency, "Power", opt.TxPower_dBm - 30, ...
    "SystemLoss", opt.TxLoss, "BitRate", user.Rb_bps/1e6};
if ant.kind == "gaussian"
    tx = transmitter(sat, txArgs{:});
    gaussianAntenna(tx, "DishDiameter", ant.D, "ApertureEfficiency", ant.eff);
else
    tx = transmitter(sat, txArgs{:}, "Antenna", ant.obj);
end
g = gimbal(gs);
pointAt(g, sat);
fixedLoss = opt.PolarizationLoss + atmLoss(atm, opt.MinElevation);
rx = receiver(g, "Name", "GS", "SystemLoss", opt.RxSystemLoss, "PreReceiverLoss", 0, ...
    "GainToNoiseTemperatureRatio", opt.GT, "RequiredEbNo", user.EbN0_req_dB + opt.MarginTarget + fixedLoss);
gaussianAntenna(rx, "DishDiameter", opt.DishDiameter, "ApertureEfficiency", opt.DishEfficiency);
lk = link(tx, rx);
[eb, t] = ebno(lk);
rip = sigstrength(lk);
eirp = txEIRP(lk);
[az, el, r] = aer(gs, sat);
[fD, ~, info] = dopplershift(sat, gs, Frequency = opt.Frequency);
dr = info.DopplerRate;                           % N-1 values between samples
rate = NaN(size(fD));
if numel(dr) >= 1
    rate(1) = dr(1); rate(end) = dr(end);
    if numel(fD) > 2, rate(2:end-1) = (dr(1:end-1) + dr(2:end))/2; end
end
if isempty(t.TimeZone), t.TimeZone = 'UTC'; end
v = el >= opt.MinElevation & isfinite(fD);
p = [];
if nnz(v) < 2, return, end                       % grazing pass: fewer than 2 samples above the mask
i = find(v, 1):find(v, 1, 'last');
lam = physconst('LightSpeed')/opt.Frequency;
grx = 10*log10(opt.DishEfficiency*(pi*opt.DishDiameter/lam)^2);
p.k = k;
p.t = t(i); p.az = az(i); p.el = el(i); p.rangeKm = r(i)/1e3;
p.rangeRate = -fD(i)*physconst('LightSpeed')/opt.Frequency;
p.fD = fD(i); p.rate = fillmissing(rate(i), 'nearest');   % the AOS/LOS samples border NaN (no LOS)
p.cn0Free = eb(i) + 10*log10(user.Rb_bps);    % C/N0 from the link object (free space, antennas, losses)
p.atmLoss = atmLoss(atm, p.el);
p.cn0 = p.cn0Free - p.atmLoss - opt.PolarizationLoss;
p.eirp = eirp(i); p.rip = rip(i);
p.prx_dBm = p.rip + 30 + grx - p.atmLoss - opt.PolarizationLoss;
p.aos = p.t(1); p.los = p.t(end);
[p.maxEl, im] = max(p.el); p.tca = p.t(im); p.iMax = im;
p.partial = v(1) || v(end);
p.linkPct = linkPercentage(lk);                % toolbox view: free space + fixed losses at the mask
p.linkIntervals = linkIntervals(lk);
end

%% ===== evaluation ================================================================

function [P, E] = evaluate(P, REQ, u, opt)
M = height(REQ); nP = numel(P); dt = opt.PassSampleTime;
Rb = REQ.Rb_bps(:).'; req = REQ.EbN0_req_dB(:).'; allowed = REQ.Allowed(:).';
okU = zeros(nP, 1); vol = zeros(nP, M); okT = zeros(nP, M); tot = zeros(nP, 1); volA = zeros(nP, 1);
for k = 1:nP
    C = P(k).cn0(:) - 10*log10(Rb);
    mg = C - req;
    ok = mg >= opt.MarginTarget;
    okU(k) = sum(ok(:, u))*dt;                  % your configuration, allowed or not
    ok = ok & allowed;                          % strategies only use allowed configurations
    vol(k, :) = sum(ok, 1).*Rb*dt/8e6;
    okT(k, :) = sum(ok, 1)*dt; tot(k) = numel(P(k).t)*dt;
    [best, jb] = max(ok.*Rb, [], 2);
    jb(best == 0) = 0;
    volA(k) = sum(best)*dt/8e6;
    P(k).ebnoUser = C(:, u).'; P(k).marginUser = mg(:, u).';
    P(k).adaptRb = best.'; P(k).adaptCfg = jb.';
end
E.notes = strings(0, 1);
scale = min(1, 86400/seconds(opt.Duration));     % never extrapolate a short window to a day
unit = "MB/day";
if hours(opt.Duration) < 24, unit = sprintf('MB in %.3g h', hours(opt.Duration)); end
E.best = struct('user', REQ.Config(u), 'userMBperDay', NaN, 'userAvail_pct', NaN, 'dayFixed', "", ...
    'dayFixedMBperDay', NaN, 'dayFixedAvail_pct', NaN, 'robust', "", 'robustMBperDay', NaN, ...
    'perPassMBperDay', NaN, 'adaptiveMBperDay', NaN, 'linkTime_s', sum(tot), 'scaleToDay', scale, 'unit', unit, ...
    'dayFixedSource', "", 'volumeMatrix_MB', []);
E.perPass = table();
if nP == 0, E.notes(end+1) = "No pass above the elevation mask in the analysed window."; return, end
avail = 100*sum(okT, 1)/sum(tot);
volU = okU*Rb(u)/8e6;
E.best.userMBperDay = sum(volU)*scale; E.best.userAvail_pct = 100*sum(okU)/sum(tot);
if ~allowed(u)
    E.notes(end+1) = sprintf(['Your configuration %s is outside the allowed set (OBW99 %.3g MHz > %.3g MHz limit, ' ...
        'or %.3g Mbit/s > BusRate): it is shown, but not used by the strategies.'], REQ.Config(u), ...
        REQ.OBW99_Hz(u)/1e6, REQ.OBWmax_Hz(u)/1e6, Rb(u)/1e6);
end
pp = (1:nP).';
if ~any(vol(:) > 0)
    E.notes(end+1) = "No allowed configuration closes the link with the target margin.";
    E.best.dayFixed = "none"; E.best.robust = "none"; E.best.dayFixedMBperDay = 0; E.best.perPassMBperDay = 0;
    E.best.adaptiveMBperDay = 0; E.best.dayFixedSource = "";
    E.perPass = table(pp, [P.maxEl].', tot, volU, 100*okU./tot, zeros(nP, 1), repmat("none", nP, 1), zeros(nP, 1), ...
        volA, 'VariableNames', {'Pass', 'MaxEl_deg', 'Time_s', 'User_MB', 'UserAvail_pct', 'DayBest_MB', ...
        'PassBestConfig', 'PassBest_MB', 'Adaptive_MB'});
    return
end
[~, jDay] = max(sum(vol, 1));
[vp, jPass] = max(vol, [], 2);
cfgPass = REQ.Config(jPass); cfgPass(vp == 0) = "none";
cand = find(avail >= 95 & allowed);
if isempty(cand), jRob = []; else, [~, i] = max(Rb(cand)); jRob = cand(i); end
E.best.dayFixed = REQ.Config(jDay); E.best.dayFixedMBperDay = sum(vol(:, jDay))*scale;
E.best.dayFixedAvail_pct = avail(jDay); E.best.dayFixedSource = REQ.Source(jDay);
if ~isempty(jRob)
    E.best.robust = REQ.Config(jRob); E.best.robustMBperDay = sum(vol(:, jRob))*scale;
else
    E.best.robust = "none";
end
E.best.perPassMBperDay = sum(vp)*scale;
E.best.adaptiveMBperDay = sum(volA)*scale;
E.best.volumeMatrix_MB = vol;
E.perPass = table(pp, [P.maxEl].', tot, volU, 100*okU./tot, vol(:, jDay), cfgPass, vp, volA, ...
    'VariableNames', {'Pass', 'MaxEl_deg', 'Time_s', 'User_MB', 'UserAvail_pct', 'DayBest_MB', ...
    'PassBestConfig', 'PassBest_MB', 'Adaptive_MB'});
if contains(REQ.Source(jDay), "estimate")
    E.notes(end+1) = "The best fixed configuration uses an ESTIMATED Eb/N0 (" + REQ.Source(jDay) + ").";
end
if scale < 1, E.notes(end+1) = sprintf('Per-day values are averages over %.1f h.', hours(opt.Duration)); end
if hours(opt.Duration) < 24
    E.notes(end+1) = sprintf('Window is %.3g h (< 1 day): the "per day" fields hold totals for the window.', hours(opt.Duration));
end
end

function [B, st] = staticCheck(P, U, opt, ant)
B = table(); st = struct('ok', false);
if isempty(P), return, end
[~, k] = max([P.maxEl]); i = P(k).iMax;
pdbw = opt.TxPower_dBm - 30;
gtx = P(k).eirp(i) - pdbw + opt.TxLoss;           % gain toward the station at that instant
c = satelliteCNRConfig;                           % (gains must be >= 0: a negative gain goes into the losses)
c.TransmitterPower = pdbw; c.TransmitterSystemLoss = opt.TxLoss; c.TransmitterAntennaGain = max(gtx, 0);
c.Distance = P(k).rangeKm(i); c.Frequency = opt.Frequency/1e9;
c.MiscellaneousLoss = P(k).atmLoss(i) + opt.PolarizationLoss + max(-gtx, 0);
c.GainToNoiseTemperatureRatio = opt.GT; c.ReceiverSystemLoss = opt.RxSystemLoss;
c.BitRate = U.Rb_bps/1e6; c.SymbolRate = U.SymbolRate_Bd/1e6; c.Bandwidth = U.OBW99_Hz/1e6;
[cnr, info] = satelliteCNR(c);
tv = P(k).ebnoUser(i);
st = struct('ok', true, 'pass', k, 'time', P(k).t(i), 'elevation', P(k).el(i), 'config', c, 'CNR_dB', cnr, ...
    'info', info, 'EbN0_static', info.ReceivedEbNo, 'EbN0_timeVarying', tv, 'difference', info.ReceivedEbNo - tv);
items = ["TX power"; "TX line loss"; "CubeSat antenna gain toward station"; "EIRP"; "Slant range"; ...
    "Free-space path loss"; "Atmospheric loss"; "Polarisation loss"; ...
    "Received isotropic power"; "Ground G/T"; "RX system loss"; "Boltzmann constant"; "C/N0"; ...
    "Payload bit rate"; "Eb/N0 available"; "Eb/N0 required (" + U.Source + ")"; "Margin"; "Target margin"; ...
    "CNR in OBW99 (" + sprintf('%.3g MHz', U.OBW99_Hz/1e6) + ")"; "Carrier power at dish output"];
vals = [pdbw; -opt.TxLoss; gtx; info.TransmitterEIRP; P(k).rangeKm(i); -info.FSPL; -P(k).atmLoss(i); ...
    -opt.PolarizationLoss; info.ReceivedIsotropicPower; opt.GT; -opt.RxSystemLoss; 228.6; ...
    info.CarrierToNoiseDensityRatio; 10*log10(U.Rb_bps); info.ReceivedEbNo; U.EbN0_req_dB; ...
    info.ReceivedEbNo - U.EbN0_req_dB; opt.MarginTarget; cnr; P(k).prx_dBm(i)];
units = ["dBW"; "dB"; "dBi"; "dBW"; "km"; "dB"; "dB"; "dB"; "dBW"; "dB/K"; "dB"; "dB (-k, dBW/K/Hz)"; "dB-Hz"; ...
    "dB-bit/s"; "dB"; "dB"; "dB"; "dB"; "dB"; "dBm"];
B = table(items, vals, units, 'VariableNames', {'Item', 'Value', 'Unit'});
st.antenna = ant.text;
end

function D = dopplerAnalysis(P, opt)
D = struct('table', table(), 'summary', table(), 'notes', strings(0, 1));
if isempty(P), return, end
fUp = opt.UplinkFrequency; ratio = fUp/opt.Frequency;
rows = cell(numel(P), 1); S = cell(numel(P), 1);
errUp = opt.FreqErrorPPM*1e-6*fUp; errDn = opt.FreqErrorPPM*1e-6*opt.Frequency;
pullUp = 0.01*opt.UplinkSymbolRate; pullDn = 0.01*opt.SymbolRate;
afc15 = 15e3; afc100 = 100e3*(opt.UplinkSymbolRate > 30e3) + min(100e3, opt.UplinkSymbolRate)*(opt.UplinkSymbolRate <= 30e3);
tol = @(range, err, rmax) max((range - err)/rmax, 0);
for k = 1:numel(P)
    p = P(k); n = numel(p.t);
    loc = p.t; loc.TimeZone = opt.TimeZone;
    up = p.fD*ratio;
    rows{k} = table(repmat(k, n, 1), p.t(:), loc(:), p.az(:), p.el(:), p.rangeKm(:), p.rangeRate(:), p.fD(:), ...
        p.rate(:), opt.Frequency + p.fD(:), up(:), fUp - up(:), 'VariableNames', {'Pass', 'UTC', 'Local', 'Az_deg', ...
        'El_deg', 'Range_km', 'RangeRate_m_s', 'Doppler_Hz', 'DopplerRate_Hz_s', 'RxTune_Hz', 'UplinkDoppler_Hz', ...
        'UplinkTx_Hz'});
    rUp = max(abs(p.rate))*ratio; rDn = max(abs(p.rate));
    S{k} = table(k, p.maxEl, max(abs(p.fD))/1e3, rDn, max(abs(up))/1e3, ...
        100*mean(abs(up) + errUp <= pullUp), 100*mean(abs(up) + errUp <= afc15), 100*mean(abs(up) + errUp <= afc100), ...
        tol(pullUp, errUp, rUp), tol(afc15, errUp, rUp), tol(afc100, errUp, rUp), tol(pullDn, errDn, rDn), ...
        'VariableNames', {'Pass', 'MaxEl_deg', 'MaxDoppler_kHz', 'MaxRate_Hz_s', 'MaxUplinkDoppler_kHz', ...
        'NoComp_PullIn_pct', 'NoComp_AFC15_pct', 'NoComp_AFC100_pct', 'TimeTol_PullIn_s', 'TimeTol_AFC15_s', ...
        'TimeTol_AFC100_s', 'DownTimeTol_PullIn_s'});
end
D.table = vertcat(rows{:});
D.table.UTC.Format = 'yyyy-MM-dd HH:mm:ss'; D.table.Local.Format = 'yyyy-MM-dd HH:mm:ss';
D.summary = vertcat(S{:});
D.limits = struct('pullInUp_Hz', pullUp, 'pullInDown_Hz', pullDn, 'afcDatasheet_Hz', afc15, 'afcFw2_Hz', afc100, ...
    'freqErrUp_Hz', errUp, 'freqErrDown_Hz', errDn, 'uplinkFrequency', fUp);
if opt.Band == "UHF"
    D.notes(end+1) = "UHF: the SRS-4 pull-in/AFC columns do not apply (different radio); use the Doppler columns only.";
    return
end
worst = min(D.summary.NoComp_AFC100_pct);
D.notes(end+1) = sprintf(['Uplink to the SRS-4 (%.1f MHz, %s): pull-in %.1f kHz, AFC 15 kHz (datasheet) / %.0f kHz (fw 2.0); ' ...
    'SRS-4 frequency error %.1f kHz (%.1f ppm). Max uplink Doppler %.1f kHz -> ground pre-compensation is required ' ...
    'unless AFC 100 kHz (in range %.0f %% of the worst pass).'], fUp/1e6, rateLabel(opt.UplinkSymbolRate), pullUp/1e3, ...
    afc100/1e3, errUp/1e3, opt.FreqErrorPPM, max(D.summary.MaxUplinkDoppler_kHz), worst);
D.notes(end+1) = sprintf(['With pre-compensation from this table, the timing error (or table step) the SRS-4 tolerates is ' ...
    '%.1f s (pull-in only) or %.0f s (AFC 15 kHz) at the worst Doppler rate.'], min(D.summary.TimeTol_PullIn_s), ...
    min(D.summary.TimeTol_AFC15_s));
if errUp >= pullUp
    D.notes(end+1) = "The SRS-4 frequency error alone exceeds the demodulator pull-in at this symbol rate: keep the AFC enabled.";
end
end

function F = fieldfoxPlan(P, opt)
f = opt.Frequency; Rs = opt.SymbolRate;
obw999 = Rs*(1.18*(opt.Rolloff == 0.2) + (1 + opt.Rolloff)*(opt.Rolloff ~= 0.2));
if opt.Band == "UHF", obw999 = 1.5*Rs; end
err = opt.FreqErrorPPM*1e-6*f + 0.7e-6*f;           % SRS-4 + FieldFox reference
rbws = [10 30 100 300 1e3 3e3 10e3 30e3 100e3 300e3 1e6 3e6];
F = struct('obw999_Hz', obw999, 'freqUncertainty_Hz', err, 'rbw_Hz', rbws(find(rbws <= max(Rs/100, 10), 1, 'last')), ...
    'trackSpan_Hz', NaN, 'fixedCenter_Hz', f, 'fixedSpan_Hz', NaN, 'table', table(), 'file', "");
if isempty(P), return, end
rmax = max(arrayfun(@(p) max(abs(p.rate)), P)); dmax = max(arrayfun(@(p) max(abs(p.fD)), P));
F.trackSpan_Hz = niceUp(obw999 + 2*(rmax*opt.TrackStep + err));
F.fixedSpan_Hz = niceUp(obw999 + 2*(dmax + err));
F.maxDoppler_Hz = dmax; F.maxRate_Hz_s = rmax;
rows = cell(numel(P), 1);
for k = 1:numel(P)
    p = P(k);
    tq = p.t(1):seconds(opt.TrackStep):p.t(end);
    x = seconds(p.t - p.t(1)); xq = seconds(tq - p.t(1));
    d = interp1(x, p.fD, xq, 'pchip'); e = interp1(x, p.el, xq, 'pchip');
    loc = tq; loc.TimeZone = opt.TimeZone;
    rows{k} = table(repmat(k, numel(tq), 1), tq(:), loc(:), e(:), f + d(:), d(:), ...
        repmat(F.trackSpan_Hz, numel(tq), 1), repmat(F.rbw_Hz, numel(tq), 1), 'VariableNames', ...
        {'Pass', 'UTC', 'Local', 'El_deg', 'Center_Hz', 'Doppler_Hz', 'Span_Hz', 'RBW_Hz'});
end
F.table = vertcat(rows{:});
F.table.UTC.Format = 'yyyy-MM-dd HH:mm:ss'; F.table.Local.Format = 'yyyy-MM-dd HH:mm:ss';
end

function v = niceUp(x)
d = 10^(floor(log10(x)) - 1); v = ceil(x/d)*d;      % round up to 2 significant digits
end

function EC = eclipseAnalysis(opt, orb, P)
% Separate scenario: eclipse() without ephemeris data breaks every later call on its scenario.
EC = struct('ok', false, 'note', "", 'intervals', table(), 'fraction', NaN, 'passInEclipse', []);
try
    sc = satelliteScenario(opt.Start, opt.Start + opt.Duration, 60);
    s = satellite(sc, orb.file, "OrbitPropagator", "sgp4");
    ec = eclipse(s);
    T = eclipseIntervals(ec);
    EC.intervals = T; EC.ok = true;
    if height(T) > 0
        if isempty(T.StartTime.TimeZone), T.StartTime.TimeZone = 'UTC'; T.EndTime.TimeZone = 'UTC'; end
        EC.fraction = sum(seconds(T.EndTime - T.StartTime))/seconds(opt.Duration);
        EC.passInEclipse = arrayfun(@(p) any(T.StartTime <= p.los & T.EndTime >= p.aos), P);
    else
        EC.fraction = 0; EC.passInEclipse = false(1, numel(P));
    end
    EC.note = sprintf('In eclipse %.1f %% of the time (eclipse/eclipseIntervals, separate scenario).', 100*EC.fraction);
catch err
    EC.note = "Eclipse not computed: " + regexprep(string(err.message), '<[^>]*>', '') + ...
        " Run aeroDataPackage once to install the ephemeris.";
end
end

function T = passTable(P, E, D, EC, opt)
T = table();
if isempty(P), return, end
n = numel(P);
aos = [P.aos].'; los = [P.los].'; tca = [P.tca].';
la = aos; la.TimeZone = opt.TimeZone; ll = los; ll.TimeZone = opt.TimeZone;
fmt = 'yyyy-MM-dd HH:mm:ss';
aos.Format = fmt; los.Format = fmt; tca.Format = fmt; la.Format = fmt; ll.Format = fmt;
T = table((1:n).', aos, los, la, ll, seconds(los - aos), [P.maxEl].', tca, arrayfun(@(p) p.az(1), P).', ...
    arrayfun(@(p) p.az(end), P).', arrayfun(@(p) min(p.rangeKm), P).', D.summary.MaxDoppler_kHz, ...
    D.summary.MaxRate_Hz_s, arrayfun(@(p) max(p.ebnoUser), P).', arrayfun(@(p) min(p.marginUser), P).', ...
    E.perPass.UserAvail_pct, E.perPass.User_MB, E.perPass.Adaptive_MB, [P.partial].', ...
    'VariableNames', {'Pass', 'AOS_UTC', 'LOS_UTC', 'AOS_Local', 'LOS_Local', 'Duration_s', 'MaxEl_deg', 'TCA_UTC', ...
    'AzAOS_deg', 'AzLOS_deg', 'MinRange_km', 'MaxDoppler_kHz', 'MaxRate_Hz_s', 'MaxEbN0_dB', 'MinMargin_dB', ...
    'UserAvail_pct', 'User_MB', 'Adaptive_MB', 'Partial'});
if EC.ok && numel(EC.passInEclipse) == n, T.InEclipse = EC.passInEclipse(:); end
end

%% ===== files, figures, report ======================================================

function files = writeFiles(R, opt)
d = opt.OutDir; files = strings(0, 1);
if height(R.passes) > 0
    files(end+1) = fullfile(d, "passes.csv"); writetable(R.passes, files(end));
    files(end+1) = fullfile(d, "doppler_table.csv"); writetable(R.doppler.table, files(end));
    files(end+1) = fullfile(d, "strategies.csv"); writetable(R.strategies, files(end));
    files(end+1) = fullfile(d, "link_budget.csv"); writetable(R.budget, files(end));
end
files(end+1) = fullfile(d, "requirements.csv"); writetable(R.requirements, files(end));
if opt.FieldFoxTrack && height(R.fieldfox.table) > 0
    files(end+1) = fullfile(d, "fieldfox_track.csv"); writetable(R.fieldfox.table, files(end));
end
files = files(:);
end

function files = makeFigures(R, S, opt)
files = strings(0, 1);
P = R.pass;
figs = {@() figGroundTrack(R, S, opt), @() figRequirements(R, opt)};
if ~isempty(P)
    figs = [figs, {@() figSky(P, opt), @() figLink(R, opt), @() figDoppler(R, opt), @() figVolume(R, opt)}];
end
for k = 1:numel(figs)
    old = findall(groot, 'Type', 'figure');
    try
        files(end+1) = figs{k}(); %#ok<AGROW>
    catch err
        warning('ff_cubesat_link:fig', 'Figure %s not written: %s', func2str(figs{k}), err.message);
        if strcmpi(opt.Visible, 'off')                  % no half-drawn hidden figure left behind
            new = findall(groot, 'Type', 'figure'); delete(new(~ismember(new, old)));
        end
    end
end
files = files(:);
end

function file = figGroundTrack(R, S, opt)
fig = newFig(opt, [80 80 1000 600]);
gx = geoaxes(fig, 'Basemap', 'darkwater');          % installed with MATLAB: no map download
hold(gx, 'on');
[la, lo] = breakWrap(S.lat, S.lon);
geoplot(gx, la, lo, '-', 'Color', [0.55 0.55 0.6], 'LineWidth', 0.6, 'DisplayName', 'Ground track');
c = lines(max(1, numel(R.pass)));
for k = 1:numel(R.pass)
    m = S.t >= R.pass(k).aos & S.t <= R.pass(k).los;
    [a, b] = breakWrap(S.lat(m), S.lon(m));
    geoplot(gx, a, b, '-', 'Color', c(k, :), 'LineWidth', 3, 'DisplayName', sprintf('Pass %d (%.0f deg)', k, R.pass(k).maxEl));
end
h = mean(S.alt); Re = 6371e3;
lambda = acosd(Re/(Re + h)*cosd(opt.MinElevation)) - opt.MinElevation;
[cl, cn] = scircle1(opt.Site(1), opt.Site(2), lambda);
geoplot(gx, cl, cn, 'r--', 'LineWidth', 1.2, 'DisplayName', sprintf('Visible above %g deg', opt.MinElevation));
geoplot(gx, opt.Site(1), opt.Site(2), 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 9, 'DisplayName', opt.SiteName);
geolimits(gx, opt.Site(1) + [-1 1]*min(80, 2.2*lambda), opt.Site(2) + [-1 1]*min(170, 3*lambda));
legend(gx, 'Location', 'eastoutside', 'FontSize', 8);
title(gx, sprintf('%s ground track, %s .. +%.0f h UTC (h = %.0f km)', R.orbit.name, ...
    string(opt.Start, 'yyyy-MM-dd HH:mm'), hours(opt.Duration), h/1e3), 'Interpreter', 'none');
file = saveFig(fig, fullfile(opt.OutDir, "fig_groundtrack.png"), opt);
end

function [la, lo] = breakWrap(la, lo)
la = la(:).'; lo = lo(:).';
j = find(abs(diff(lo)) > 180);
for k = numel(j):-1:1
    la = [la(1:j(k)) NaN la(j(k)+1:end)]; lo = [lo(1:j(k)) NaN lo(j(k)+1:end)];
end
end

function file = figSky(P, opt)
fig = newFig(opt, [80 80 640 600]);
n = max(arrayfun(@(p) numel(p.t), P));
az = NaN(n, numel(P)); el = az;                  % rows = time, columns = passes; right-aligned so that
for k = 1:numel(P)                                 % the last row (marker + label) is each pass's LOS
    m = numel(P(k).t); az(n-m+1:n, k) = P(k).az; el(n-m+1:n, k) = P(k).el;
end
lbl = "P" + (1:numel(P));
try
    skyplot(az, el, lbl, 'MaskElevation', opt.MinElevation, 'GroupData', categorical(lbl));  % Navigation Toolbox
catch
    pax = polaraxes(fig); hold(pax, 'on');
    pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise'; pax.RLim = [0 90];
    pax.RTick = 0:30:90; pax.RTickLabel = {'90', '60', '30', '0'};
    for k = 1:numel(P), polarplot(pax, deg2rad(P(k).az), 90 - P(k).el, 'LineWidth', 1.5); end
    legend(pax, lbl);
end
title(sprintf('Sky plot from %s (mask %g deg)', opt.SiteName, opt.MinElevation));
file = saveFig(fig, fullfile(opt.OutDir, "fig_skyplot.png"), opt);
end

function file = figLink(R, opt)
P = R.pass; U = R.user;
fig = newFig(opt, [80 80 1100 720]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
c = lines(numel(P));
ax = gobjects(4, 1);
for i = 1:4, ax(i) = nexttile(tl); hold(ax(i), 'on'); grid(ax(i), 'on'); end
for k = 1:numel(P)
    x = minutes(P(k).t - P(k).t(1)); nm = sprintf('Pass %d', k);
    plot(ax(1), x, P(k).el, 'Color', c(k, :), 'DisplayName', nm);
    plot(ax(2), x, P(k).ebnoUser, 'Color', c(k, :), 'DisplayName', nm);
    plot(ax(3), x, P(k).marginUser, 'Color', c(k, :), 'DisplayName', nm);
    stairs(ax(4), x, P(k).adaptRb/1e6, 'Color', c(k, :), 'DisplayName', nm);
end
yline(ax(2), U.EbN0_req_dB, 'k--', 'Required', 'HandleVisibility', 'off');
yline(ax(2), U.EbN0_req_dB + opt.MarginTarget, 'r--', 'Required + target', 'HandleVisibility', 'off');
yline(ax(3), opt.MarginTarget, 'r--', 'Target', 'HandleVisibility', 'off');
yline(ax(3), 0, 'k-', 'HandleVisibility', 'off');
yline(ax(4), U.Rb_bps/1e6, 'k--', 'Your configuration', 'HandleVisibility', 'off');
ylabel(ax(1), 'Elevation (deg)'); ylabel(ax(2), 'Eb/N0 available (dB)'); ylabel(ax(3), 'Margin (dB)');
ylabel(ax(4), 'Adaptive payload rate (Mbit/s)');
for i = 1:4, xlabel(ax(i), 'Minutes since AOS'); end
title(ax(1), 'Elevation'); title(ax(2), "Eb/N0, " + U.Config);
title(ax(3), "Margin, " + U.Config); title(ax(4), 'Best feasible configuration (upper bound)');
legend(ax(1), 'Location', 'best', 'FontSize', 7);
title(tl, sprintf('%s -> %s: %.1f dBm, %s, G/T %.1f dB/K, %s', R.orbit.name, opt.SiteName, opt.TxPower_dBm, ...
    R.antenna.kind, opt.GT, R.atmosphere.method), 'Interpreter', 'none', 'FontSize', 9);
file = saveFig(fig, fullfile(opt.OutDir, "fig_link.png"), opt);
end

function file = figDoppler(R, opt)
P = R.pass; L = R.doppler.limits;
fig = newFig(opt, [80 80 1000 700]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');
a1 = nexttile(tl); hold(a1, 'on'); grid(a1, 'on');
a2 = nexttile(tl); hold(a2, 'on'); grid(a2, 'on');
c = lines(numel(P));
ratio = opt.UplinkFrequency/opt.Frequency;
for k = 1:numel(P)
    x = minutes(P(k).t - P(k).tca);
    plot(a1, x, P(k).fD/1e3, '-', 'Color', c(k, :), 'DisplayName', sprintf('Pass %d downlink', k));
    if opt.Band == "S"
        plot(a1, x, P(k).fD*ratio/1e3, ':', 'Color', c(k, :), 'HandleVisibility', 'off');
    end
    plot(a2, x, P(k).rate, 'Color', c(k, :), 'DisplayName', sprintf('Pass %d', k));
end
if opt.Band == "S"
    yline(a1, [-1 1]*L.afcDatasheet_Hz/1e3, 'k--', 'AFC 15 kHz', 'HandleVisibility', 'off', 'LabelVerticalAlignment', 'top');
    yline(a1, [-1 1]*L.pullInUp_Hz/1e3, 'r:', 'pull-in 1 % Rs', 'HandleVisibility', 'off', ...
        'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
    if L.afcFw2_Hz <= 1.2*max(abs([P.fD]))*ratio
        yline(a1, [-1 1]*L.afcFw2_Hz/1e3, 'b--', 'AFC fw 2.0', 'HandleVisibility', 'off');
    end
end
ylabel(a1, 'Doppler (kHz)'); xlabel(a1, 'Minutes from closest approach');
ylabel(a2, 'Doppler rate (Hz/s)'); xlabel(a2, 'Minutes from closest approach');
title(a1, sprintf('Doppler at %.1f MHz (solid); uplink %.1f MHz dotted', opt.Frequency/1e6, opt.UplinkFrequency/1e6));
legend(a1, 'Location', 'eastoutside', 'FontSize', 7);
title(a2, 'Doppler rate (dopplershift info.DopplerRate)');
file = saveFig(fig, fullfile(opt.OutDir, "fig_doppler.png"), opt);
end

function file = figVolume(R, opt)
T = R.strategies;
fig = newFig(opt, [80 80 900 520]);
ax = axes(fig);
b = bar(ax, categorical("P" + T.Pass, "P" + T.Pass), [T.User_MB T.DayBest_MB T.PassBest_MB T.Adaptive_MB]);
legend(ax, ["Yours: " + R.best.user, "Best fixed: " + R.best.dayFixed, "Best per pass", "Adaptive (bound)"], ...
    'Location', 'northoutside', 'FontSize', 8, 'NumColumns', 2);
ylabel(ax, 'Payload per pass (MB)'); grid(ax, 'on');
title(ax, sprintf('Yours %.1f, best fixed %.1f, adaptive %.1f %s (margin >= %g dB)', R.best.userMBperDay, ...
    R.best.dayFixedMBperDay, R.best.adaptiveMBperDay, R.best.unit, opt.MarginTarget), 'FontSize', 9);
b(1).FaceColor = [0.3 0.3 0.3];
file = saveFig(fig, fullfile(opt.OutDir, "fig_volume.png"), opt);
end

function file = figRequirements(R, opt)
Q = R.requirements;
fig = newFig(opt, [80 80 900 560]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
mods = unique(Q.Modulation, 'stable'); fecs = unique(Q.FEC, 'stable');
c = lines(numel(fecs)); st = ["-o" "-s" "--d"];
for m = 1:numel(mods)
    for f = 1:numel(fecs)
        q = Q(Q.Modulation == mods(m) & Q.FEC == fecs(f), :);
        if isempty(q), continue, end
        plot(ax, q.SymbolRate_Bd, q.EbN0_req_dB, st(min(m, 3)), 'Color', c(f, :), 'MarkerSize', 4, ...
            'DisplayName', mods(m) + " " + fecs(f));
    end
end
set(ax, 'XScale', 'log'); xlabel(ax, 'Symbol rate (Bd)'); ylabel(ax, 'Required Eb/N0 per payload bit (dB)');
legend(ax, 'Location', 'eastoutside', 'FontSize', 7);
if opt.Band == "S"
    title(ax, sprintf('SRS-4 Table 7 (PER < 1 %%, NF %.2f dB) -> Eb/N0; 8PSK dashed = estimate', opt.NoiseFigure));
else
    title(ax, 'UHF estimate (theory + implementation loss)');
end
file = saveFig(fig, fullfile(opt.OutDir, "fig_requirements.png"), opt);
end

function fig = newFig(opt, pos)
fig = figure('Visible', opt.Visible, 'Color', 'w', 'Position', pos);
end

function file = saveFig(fig, file, opt)
exportgraphics(fig, file, 'Resolution', 130);
if strcmpi(opt.Visible, 'off'), close(fig); end
end

function printSummary(R, opt)
fprintf('\n===== %s from %s (%.4f, %.4f) =====\n', R.orbit.name, opt.SiteName, opt.Site(1), opt.Site(2));
if height(R.passes) == 0
    fprintf('No pass above %g deg in the window.\n', opt.MinElevation);
else
    V = R.passes(:, {'Pass', 'AOS_Local', 'Duration_s', 'MaxEl_deg', 'MaxDoppler_kHz', 'MaxEbN0_dB', 'MinMargin_dB', 'User_MB'});
    disp(V);
    B = R.best;
    fprintf('Your configuration %s: %.1f %s, margin >= %g dB for %.0f %% of the pass time\n', B.user, ...
        B.userMBperDay, B.unit, opt.MarginTarget, B.userAvail_pct);
    fprintf('Best fixed: %s (%.1f %s); robust (>= 95 %% closed): %s; adaptive bound %.1f %s\n', ...
        B.dayFixed, B.dayFixedMBperDay, B.unit, B.robust, B.adaptiveMBperDay, B.unit);
    if R.static.ok
        fprintf('Static check (satelliteCNR) at %.1f deg: Eb/N0 %.2f dB vs time-varying %.2f dB\n', ...
            R.static.elevation, R.static.EbN0_static, R.static.EbN0_timeVarying);
    end
end
if ~isempty(R.notes), fprintf('Notes:\n'); fprintf('  - %s\n', R.notes); end
end

function file = makeReport(R, opt)
import mlreportgen.report.*
import mlreportgen.dom.*
fprintf('Writing the %s report ...\n', upper(opt.Report));
rpt = Report(char(fullfile(opt.OutDir, "cubesat_link_report")), char(opt.Report));
when = datetime('now', 'TimeZone', opt.TimeZone, 'Format', 'yyyy-MM-dd HH:mm z');
add(rpt, TitlePage('Title', sprintf('%s link plan: SRS-4 %s band', R.orbit.name, opt.Band), ...
    'Subtitle', sprintf('%s (%.4f, %.4f), %s + %.0f h', opt.SiteName, opt.Site(1), opt.Site(2), ...
    string(opt.Start, 'yyyy-MM-dd HH:mm ''UTC'''), hours(opt.Duration)), ...
    'Author', string(getenv('USERNAME')), 'PubDate', char(when)));
add(rpt, TableOfContents);
png = @(n) char(fullfile(opt.OutDir, n));
ch = Chapter('Title', 'Summary');
p = Paragraph('PREDICTION from a TLE and datasheet values, not a measurement. Check the assumptions in the last chapter.');
p.Bold = true; p.Color = 'red'; add(ch, p);
B = R.best;
L = {sprintf('Orbit: %s (NORAD %d), %s; TLE epoch %s UTC, %.1f days before Start.', R.orbit.name, R.orbit.norad, ...
    R.orbit.source, string(R.orbit.epoch, 'yyyy-MM-dd HH:mm'), R.orbit.ageDays), ...
    sprintf('%d passes above %g deg; in view %.2f %% of the time.', height(R.passes), opt.MinElevation, R.survey.accessPct), ...
    sprintf('Your configuration %s: %.1f %s (margin >= %g dB %.0f %% of pass time).', B.user, B.userMBperDay, B.unit, ...
    opt.MarginTarget, B.userAvail_pct), ...
    sprintf('Best fixed configuration: %s, %.1f %s. Fastest that closes >= 95 %%: %s. Best per pass: %.1f, adaptive bound: %.1f (same unit).', ...
    B.dayFixed, B.dayFixedMBperDay, B.unit, B.robust, B.perPassMBperDay, B.adaptiveMBperDay)};
add(ch, UnorderedList(L));
if R.eclipse.ok, add(ch, Paragraph(R.eclipse.note)); end
if ~isempty(R.notes), add(ch, Paragraph('Notes:')); add(ch, UnorderedList(cellstr(R.notes))); end
if strlength(opt.Notes) > 0, add(ch, Paragraph("User notes: " + opt.Notes)); end
add(rpt, ch);
ch = Chapter('Orbit and passes');
add(ch, Paragraph(sprintf('TLE used (%s):', R.orbit.source)));
add(ch, Preformatted(char(strjoin(R.orbit.lines, newline))));
for n = ["fig_groundtrack.png" "fig_skyplot.png"]
    if isfile(png(n)), add(ch, FormalImage('Image', png(n), 'Width', '15cm')); end
end
if height(R.passes) > 0
    add(ch, rptTable(R.passes(:, {'Pass', 'AOS_Local', 'LOS_Local', 'Duration_s', 'MaxEl_deg', 'AzAOS_deg', ...
        'AzLOS_deg', 'MinRange_km', 'MaxDoppler_kHz', 'MinMargin_dB', 'User_MB'}), ...
        sprintf('Passes (local time %s)', opt.TimeZone), 30));
end
add(rpt, ch);
ch = Chapter('Link budget');
add(ch, Paragraph(sprintf(['Transmitter: SRS-4 %.1f dBm, %.1f dB line loss, %s, pointing %s. Receiver: G/T %.1f dB/K (%s) ' ...
    'on a gimbal tracking the satellite, %.1f dB system loss, %.1f dB polarisation loss. Atmosphere: %s.'], ...
    opt.TxPower_dBm, opt.TxLoss, R.antenna.text, opt.CubeSatPointing, opt.GT, opt.GTSource, opt.RxSystemLoss, ...
    opt.PolarizationLoss, R.atmosphere.method)));
if R.static.ok
    add(ch, rptTable(R.budget, sprintf('Static link budget at the highest elevation (%.1f deg), %s (satelliteCNR)', ...
        R.static.elevation, R.user.Config), 30));
    add(ch, Paragraph(sprintf(['Cross-check: satelliteCNR gives Eb/N0 %.2f dB, the time-varying link (ebno) %.2f dB ' ...
        'at the same instant (difference %.3f dB).'], R.static.EbN0_static, R.static.EbN0_timeVarying, R.static.difference)));
end
if isfile(png("fig_link.png")), add(ch, FormalImage('Image', png("fig_link.png"), 'Width', '16cm')); end
add(rpt, ch);
ch = Chapter('SRS-4 configurations and data volume');
add(ch, Paragraph(['Required Eb/N0 per payload bit from the SRS-4 datasheet Table 7 (PER < 1 %, NF 2.05 dB): ' ...
    'Eb/N0 = Psens + 174 - NF - 10 log10(Rb). Theory = ideal coherent AWGN for the same PER; ImplLoss = Table 7 - theory.']));
if isfile(png("fig_requirements.png")), add(ch, FormalImage('Image', png("fig_requirements.png"), 'Width', '15cm')); end
add(ch, rptTable(R.theory, 'Theory used for the comparison (PER = 1 %, 223-byte frame)', 10));
add(ch, Paragraph(['The CC theory is a union bound (pessimistic), so a slightly negative implementation loss is ' ...
    'expected. The ~2.6 dB gap of RS-only is plausibly frame synchronisation: at the required point the channel ' ...
    'BER is ~4e-3, where an exact 32-bit ASM match fails for ~12 % of frames (interpretation, not a datasheet statement).']));
Q = R.requirements(:, {'Config', 'Rb_bps', 'EbN0_req_dB', 'Theory_EbN0_dB', 'ImplLoss_dB', 'OBW99_Hz', 'Allowed', 'Source'});
add(ch, rptTable(Q, 'Required Eb/N0 per configuration (Allowed = OBW and bus-rate limits)', 120));
if height(R.strategies) > 0
    add(ch, rptTable(R.strategies, 'Payload per pass (MB) for each strategy', 30));
    if isfile(png("fig_volume.png")), add(ch, FormalImage('Image', png("fig_volume.png"), 'Width', '15cm')); end
end
add(rpt, ch);
ch = Chapter('Doppler and frequency plan');
if height(R.doppler.summary) > 0
    add(ch, UnorderedList(cellstr(R.doppler.notes)));
    add(ch, rptTable(R.doppler.summary, 'Doppler per pass: % of time inside the SRS-4 ranges without compensation; time tolerance with compensation', 30));
    if isfile(png("fig_doppler.png")), add(ch, FormalImage('Image', png("fig_doppler.png"), 'Width', '15cm')); end
    F = R.fieldfox;
    add(ch, Paragraph(sprintf(['FieldFox capture: fixed span %.3g MHz at %.4f MHz, or Doppler tracking with span %.3g kHz ' ...
        '(step %g s), RBW %g kHz. Frequency uncertainty (SRS-4 + FieldFox reference) +-%.1f kHz.'], F.fixedSpan_Hz/1e6, ...
        F.fixedCenter_Hz/1e6, F.trackSpan_Hz/1e3, opt.TrackStep, F.rbw_Hz/1e3, F.freqUncertainty_Hz/1e3)));
end
add(rpt, ch);
ch = Chapter('Assumptions and files');
add(ch, UnorderedList({ ...
    'SGP4 propagation from the TLE; old TLEs shift pass times (drag). Keplerian input has no drag term.', ...
    'pointAt is ideal attitude; the Antenna Toolbox patch is linearly polarised and on a 1-wavelength ground plane.', ...
    'Table 7 is the SRS-4 receiver; the downlink requirement assumes an equally good ground demodulator.', ...
    '8PSK and UHF requirements are estimates (no measured data).', ...
    'P.618 losses are statistics (exceedance of an average year), not a weather forecast.', ...
    'Adaptive = per-second best feasible configuration; the SRS-4 has no ACM, so it is an upper bound.', ...
    'Absolute frequency claims need a GPS-disciplined 10 MHz reference on the FieldFox.'}));
add(ch, UnorderedList(cellstr(R.files)));
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
    if isdatetime(x), T.(v) = string(x, 'MM-dd HH:mm:ss');
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

function s = stamp()
s = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
end
