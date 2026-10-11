function R = ff_gnss_margin(ip, varargin)
%FF_GNSS_MARGIN  GPS L1 interference margin: satellite link budget + FieldFox L1 noise.
%
%   R = ff_gnss_margin                       measure with the FieldFox at 10.90.36.200
%   R = ff_gnss_margin(ip, Name, Value, ...)
%   R = ff_gnss_margin("", ...)              no instrument: satellites only, or a saved
%                                            spectrum given with 'Spectrum'
%   R = ff_gnss_margin('Spectrum', s, ...)   same as ff_gnss_margin("", 'Spectrum', s, ...)
%
% Examples
%   R = ff_gnss_margin("10.90.36.200", 'Lat', 13.7275, 'Lon', 100.7782, 'Alt', 5);
%   R = ff_gnss_margin("", 'Spectrum', "gnss_l1_..._off.mat", 'SpectrumOn', "gnss_l1_..._on.mat");
%   R = ff_gnss_margin("", 'Antenna', [10 -2; 30 1; 60 3; 90 4], 'PolLoss', 0);  % datasheet dBic
%
% How much can a solar inverter's noise at 1575.42 MHz weaken GPS reception?
%   1. Satellites (Satellite Communications Toolbox): the GPS constellation from an
%      almanac (satellite), your site (groundStation), an L1 transmitter on every
%      satellite (EIRP ~26.8 dBW) and a receiver whose antenna is an Antenna
%      Toolbox patch designed for L1 (design(patchMicrostrip,1575.42e6)) pointing
%      at the zenith. link + sigstrength/ebno give the received power and C/N0 of
%      every satellite above the mask, aer its azimuth/elevation/range and
%      accessIntervals a 24 h visibility timeline. For a SEM almanac the Navigation
%      Toolbox (semread + gnssconstellation + lookangles) recomputes az/el as a
%      cross-check of PRNs, time zone and site (same almanac orbit model).
%   2. Receiver chain (RF Toolbox rfbudget): antenna -> [cable, e.g. a .s2p from
%      ff_vna] -> LNA -> receiver. Its noise figure F gives the system noise
%      temperature Tsys = Tant + 290*(F-1) and N0 = k*Tsys at the antenna.
%   3. FieldFox (SA mode, tracking generator off: no RF output; all settings are
%      saved first and restored on every exit): 1575.42 MHz +-10 MHz, RBW 100 kHz,
%      preamp on, 0 dB attenuation, average (RMS) detector, 'Sweeps' sweeps averaged
%      in mW, ADC over-range checked after each state (if the firmware supports it).
%      Noise density = trace / (1.06*RBW). Mode "offon" measures with the inverter
%      OFF, then ON. The extra noise I0 = ON - OFF (weighted with the C/A sinc^2
%      spectrum), divided by the chain gain to refer it to the antenna, lowers the
%      C/N0 of every satellite by 10*log10(1 + I0/N0).
%   4. Table, console report and one figure: sky plot coloured by C/N0, L1
%      spectrum OFF/ON with the C/A spectrum, C/N0 bars OFF vs ON, visibility.
%
% Set-up: connect the FieldFox (port 2) where the GNSS receiver would be:
% GNSS antenna -> [cable] -> LNA -> FieldFox ('MeasuredAt' "receiver"). An active
% GNSS antenna needs DC power through a bias tee (DC block towards the FieldFox).
%
% Name-Value options (defaults):
%  Site and satellites
%   'Lat','Lon','Alt' 13.7563, 100.5018, 10  site (deg, deg, m above the ellipsoid).
%                                     The default is central Bangkok: SET YOUR LAB.
%   'Time'          now               datetime; no time zone = Asia/Bangkok. Default
%                                     with a measurement/recording: its start time
%   'MaskAngle'     10                deg; satellites below it are not used
%   'Almanac'       shipped gpsAlmanac.txt   SEM almanac (.al3/.txt), TLE file or
%                                     RINEX GPS navigation file (.rnx/.YYn)
%   'EIRP'          26.8              dBW per satellite (L1 C/A)
%   'Antenna'       "patch"           "patch", "isotropic" (0 dBi), your own Antenna
%                                     Toolbox element (boresight +z = zenith, x =
%                                     north) or a datasheet table [el_deg gain_dBi]
%   'PolLoss'       3                 dB; linear patch receiving the RHCP signal
%                                     (0 for a datasheet table in dBic)
%   'AntennaTemp'   100               K; sky + ground noise the antenna sees (50-150)
%   'TimelineHours' 24                length of the visibility timeline
%  Receiver chain
%   'Cable'         []                .s2p file, sparameters object, ff_vna result
%                                     or a loss in dB; [] = no cable
%   'CableAfterLNA' false             true: antenna -> LNA -> cable (active antenna)
%   'LNAGain'       28                dB (LNAGain 0 and LNANF 0 = no LNA)
%   'LNANF'         1.5               dB
%   'ReceiverNF'    3                 dB, GNSS receiver front end
%   'MeasuredAt'    "receiver"        where the FieldFox is: "receiver" (after the
%                                     cable/LNA, instead of the GNSS receiver) or
%                                     "antenna" (at the antenna, 0 dB gain)
%  Measurement
%   'Mode'          "offon"           "offon" (inverter OFF, then ON) or "single"
%   'Sweeps'        50                sweeps averaged per state
%   'RBW'           100e3             Hz
%   'Span'          20e6              Hz
%   'Points'        401
%   'MinRise'       0.1               dB; smallest noise rise trusted as real
%                                     (FieldFox repeatability between OFF and ON)
%   'Spectrum'      []                offline data instead of the FieldFox: struct
%                                     with f (Hz), P (dBm, sweeps x points),
%                                     optional P_on and rbw; or a .mat recording
%                                     (f, P, settings) of this tool / fieldfox_live
%   'SpectrumOn'    []                second recording (inverter ON) for 'Spectrum'
%  Thresholds and output
%   'TrackThreshold' 25               dB-Hz, typical tracking threshold
%   'AcqThreshold'   35               dB-Hz, typical acquisition threshold
%   'Visible'       'on'              figure visibility ('off': deleted once saved)
%   'Prompt'        true             false = no input() prompts ("offon" measures
%                                     twice back to back)
%   'Save'          true              recordings + results (.mat/.csv) + figure .png
%   'OutDir'        pwd
%   'Port'          5025
%
% Output R: sats (one row per satellite above the mask: PRN, Az/El/Range, Pr,
% C/N0 nominal and with interference, margins), noise (measurement), margin
% (allowed noise rise), chain (rfbudget, NF, Tsys, N0), timeline, scenario
% (satelliteScenario: satelliteScenarioViewer(R.scenario) shows it in 3-D),
% navCrossCheck_deg (satcom vs Navigation Toolbox az/el), files, figure.
% Recordings hold f, P, t, settings, idn like fieldfox_live.
%
% Read before trusting the numbers:
%  * Detection limit: the FieldFox adds its own noise and cannot tell it apart
%    from external noise. Only the OFF->ON change is blamed on the inverter; the
%    OFF measurement is the reference. "single" counts everything above the
%    expected antenna+chain noise as interference: a pessimistic upper bound.
%  * Without an LNA in front the FieldFox (preamp on) is far less sensitive than
%    the antenna's own noise, so only strong interference shows. Measure after
%    an LNA / active antenna. With nothing connected you see only the FieldFox's
%    own floor, well below the expected noise (the report says so).
%  * The shipped almanac is from 2020: positions and PRNs are only indicative.
%    Use a current SEM almanac (navcen.uscg.gov, GPS almanacs), a GPS TLE file
%    (celestrak.org) or a RINEX navigation file for the real sky.
%  * design(patchMicrostrip) is an ideal, linearly polarised air-substrate patch
%    (~10 dBi peak, 3 dB polarisation loss) whose low-elevation gain depends on
%    azimuth (deep dips near the horizon). Real ceramic RHCP GNSS patches give
%    ~3-5 dBic at the zenith, nearly independent of azimuth: for realistic C/N0
%    give the datasheet gain table.
%  * EIRP, antenna temperature and the thresholds are typical values;
%    interference is assumed to enter through the antenna. Narrow (CW) spurs can
%    hurt the C/A code more than this spectral estimate: they are listed
%    separately. FieldFox absolute accuracy is about +-1.5 dB.

if nargin < 1, ip = "10.90.36.200"; end
if isempty(ip), ip = ""; end
opt = struct('Lat', 13.7563, 'Lon', 100.5018, 'Alt', 10, ...
    'Time', datetime('now', 'TimeZone', 'Asia/Bangkok'), 'MaskAngle', 10, ...
    'Almanac', fullfile(matlabroot, 'toolbox', 'shared', 'orbit', 'orbitdata', 'gpsAlmanac.txt'), ...
    'EIRP', 26.8, 'Antenna', "patch", 'PolLoss', 3, 'AntennaTemp', 100, 'TimelineHours', 24, ...
    'Cable', [], 'CableAfterLNA', false, 'LNAGain', 28, 'LNANF', 1.5, 'ReceiverNF', 3, ...
    'MeasuredAt', "receiver", 'Mode', "offon", 'Sweeps', 50, 'RBW', 100e3, 'Span', 20e6, ...
    'Points', 401, 'MinRise', 0.1, 'Spectrum', [], 'SpectrumOn', [], ...
    'TrackThreshold', 25, 'AcqThreshold', 35, 'Visible', 'on', 'Prompt', true, ...
    'Save', true, 'OutDir', pwd, 'Port', 5025);
names = fieldnames(opt);
if (ischar(ip) || isstring(ip)) && any(strcmpi(names, ip)) && mod(numel(varargin), 2) == 1
    varargin = [{ip} varargin]; ip = "";       % called as ff_gnss_margin('Spectrum', s, ...)
end
ip = string(ip);
if mod(numel(varargin), 2), error('Options must be Name, Value pairs.'); end
timeGiven = false;
for k = 1:2:numel(varargin)
    j = find(strcmpi(names, varargin{k}), 1);
    if isempty(j), error('Unknown option "%s".', varargin{k}); end
    opt.(names{j}) = varargin{k + 1};
    timeGiven = timeGiven || names{j} == "Time";
end
opt.Mode = lower(string(opt.Mode));
opt.MeasuredAt = lower(string(opt.MeasuredAt));
if ~ismember(opt.Mode, ["offon" "single"]), error('''Mode'' must be "offon" or "single".'); end
if ~ismember(opt.MeasuredAt, ["receiver" "antenna"]), error('''MeasuredAt'' must be "receiver" or "antenna".'); end
if ~(isscalar(opt.Sweeps) && opt.Sweeps >= 1), error('''Sweeps'' must be at least 1.'); end
if ~(isscalar(opt.TimelineHours) && opt.TimelineHours > 0), error('''TimelineHours'' must be > 0.'); end

L1 = 1575.42e6;                 % GPS L1 carrier (Hz)
Rc = 1.023e6;                   % C/A chip rate (Hz): main lobe = L1 +- Rc
Rb = 50;                        % navigation data rate (bit/s): C/N0 = Eb/N0 + 10*log10(Rb)
kB = physconst('Boltzmann');
stamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
files = strings(0, 1);
if opt.Save && ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
if abs(opt.Lat - 13.7563) < 1e-9 && abs(opt.Lon - 100.5018) < 1e-9
    fprintf('Site = default (central Bangkok). Pass ''Lat'', ''Lon'', ''Alt'' of your lab.\n');
end
if ~isfile(opt.Almanac), error('Almanac file not found: %s', opt.Almanac); end
C = rxChain(opt, L1, Rc, kB);                   % RF Toolbox chain (checks 'Cable' early)
A = rxAntenna(opt.Antenna, L1);                 % Antenna Toolbox (checks 'Antenna' early)

%% 1. L1 noise spectrum: FieldFox (no RF output) or a saved spectrum
m = [];
if ~isempty(opt.Spectrum)
    m = loadSpectrum(opt.Spectrum, opt.SpectrumOn, opt.RBW);
elseif strlength(ip) > 0
    m = measureL1(ip, opt, L1);
    if opt.Save, files = [files; saveRecordings(m, opt.OutDir, stamp)]; end
end
if ~timeGiven && ~isempty(m) && ~isempty(m.toff) && ~isnat(m.toff(1))
    opt.Time = m.toff(1);                       % satellites as during the measurement
end
T = opt.Time;
if isempty(T.TimeZone), T.TimeZone = 'Asia/Bangkok'; end

%% 2. Satellites (Satellite Communications + Antenna + Navigation Toolbox)
fprintf('Satellite scenario (%s, %g h)... ', string(T, 'yyyy-MM-dd HH:mm'), opt.TimelineHours);
tic;
G = gnssLinks(opt, T, A, C.Tsys_K, L1, Rb, kB);
fprintf('%.0f s\n', toc);
sats = G.sats;
if ~isempty(sats), C.budget.AvailableInputPower = max(sats.Pr_dBW) + 30; end
N = l1Noise(m, opt, C, L1, Rc, kB);

%% 3. Per-satellite C/N0 and margins
sats.CN0_Interf_dBHz = sats.CN0_Nominal_dBHz - N.degradation_dB;
sats.MarginTrack_dB  = sats.CN0_Interf_dBHz - opt.TrackThreshold;
sats.MarginAcq_dB    = sats.CN0_Interf_dBHz - opt.AcqThreshold;
M = margins(sats, opt, C, N, kB);

%% 4. Report, figure, files
report(opt, T, G, C, N, M, sats);
fig = plotResults(opt, T, G, C, N, M, sats, L1, Rc);
if opt.Save
    base = fullfile(opt.OutDir, "gnss_margin_" + stamp);
    out = struct('sats', sats, 'noise', rmfield(N, 'spurs'), 'spurs', N.spurs, 'margin', M, ...
        'chain', rmfield(C, 'budget'), 'timeline', G.timeline, 'site', [opt.Lat opt.Lon opt.Alt], ...
        'time', T, 'almanac', string(opt.Almanac), 'almanacEpoch', G.epoch);
    save(base + ".mat", '-struct', 'out');
    writetable(sats, base + ".csv");
    exportgraphics(fig, base + ".png", 'Resolution', 130);
    files = [files; base + ".mat"; base + ".csv"; base + ".png"];
    fprintf('Saved: %s\n', strjoin(files, ', '));
end
if strcmpi(char(matlab.lang.OnOffSwitchState(opt.Visible)), 'off') && (opt.Save || nargout == 0)
    delete(fig);                                % hidden figure: do not leave it in memory
end

R = struct('sats', sats, 'noise', N, 'margin', M, 'chain', C, 'timeline', G.timeline, ...
    'scenario', G.scenario, 'antenna', G.antenna, 'time', T, 'site', [opt.Lat opt.Lon opt.Alt], ...
    'almanac', string(opt.Almanac), 'almanacEpoch', G.epoch, 'almanacAge_days', G.age, ...
    'cn0CrossCheck_dB', G.crossCheck, 'navCrossCheck_deg', G.navCheck, 'measurement', m, ...
    'files', files, 'figure', fig);
end

%% ===== receiver chain ===================================================

function C = rxChain(opt, f0, Rc, kB)
% Antenna -> [cable] -> LNA -> receiver as an rfbudget (Friis): cascaded gain and
% noise figure referred to the antenna terminal.
cab = cableElement(opt.Cable, f0);
lna = amplifier('Name', 'LNA', 'Gain', opt.LNAGain, 'NF', opt.LNANF);
rcv = amplifier('Name', 'Receiver', 'Gain', 0, 'NF', opt.ReceiverNF);
if isempty(cab),         els = {lna, rcv};
elseif opt.CableAfterLNA, els = {lna, cab, rcv};
else,                    els = {cab, lna, rcv};
end
C.budget = rfbudget([els{:}], f0, -130, 2*Rc);
C.stages = string(cellfun(@(e) e.Name, els, 'UniformOutput', false));
C.NF_dB = C.budget.NF(end);
C.Tant_K = opt.AntennaTemp;
C.Tsys_K = opt.AntennaTemp + 290*(10^(C.NF_dB/10) - 1);
C.N0_dBWHz = 10*log10(kB*C.Tsys_K);
if opt.MeasuredAt == "receiver"          % FieldFox replaces the receiver: all stages before it
    C.Gmeas_dB = C.budget.TransducerGain(end - 1);
    nfMeas = C.budget.NF(end - 1);
else
    C.Gmeas_dB = 0; nfMeas = 0;
end
% Noise the antenna and chain alone should give at the FieldFox input
C.expected_dBmHz = 10*log10(kB*(opt.AntennaTemp + 290*(10^(nfMeas/10) - 1))) + 30 + C.Gmeas_dB;
end

function e = cableElement(c, f0)
if isempty(c), e = []; return, end
if isnumeric(c)
    e = attenuator('Name', 'Cable', 'Attenuation', c);
    return
end
if isstruct(c) && isfield(c, 'S'), c = c.S; end            % ff_vna result
if ischar(c) || isstring(c), c = sparameters(char(c)); end
if ~isa(c, 'sparameters') || c.NumPorts ~= 2
    error('''Cable'' must be a 2-port .s2p file, an sparameters object, an ff_vna result or a loss in dB.');
end
fr = c.Frequencies;
if f0 < fr(1) || f0 > fr(end)
    error('The cable S-parameters (%.4g-%.4g MHz) do not cover 1575.42 MHz.', fr(1)/1e6, fr(end)/1e6);
end
% Value at L1 with |S| and phase interpolated separately: nport's linear complex
% interpolation shrinks |S21| of a long cable between measured points.
S0 = zeros(2);
for k = 1:4
    v = squeeze(c.Parameters(mod(k - 1, 2) + 1, ceil(k/2), :));
    S0(k) = interp1(fr, abs(v), f0)*exp(1j*interp1(fr, unwrap(angle(v)), f0));
end
e = nport(sparameters(repmat(S0, 1, 1, 3), f0 + [-5e6 0 5e6], c.Impedance));
e.Name = 'Cable';
end

%% ===== satellites =======================================================

function G = gnssLinks(opt, T, A, Tsys, f0, Rb, kB)
step = max(60, opt.TimelineHours*3600/720);   % <= 720 time steps keeps accessIntervals quick
sc = satelliteScenario(T, T + hours(opt.TimelineHours), step);
alm = string(opt.Almanac);
if ~isempty(regexpi(alm, '\.(rnx|\d\dn|nav)$', 'once'))
    sat = satellite(sc, rinexread(alm));       % RINEX navigation file (Navigation Toolbox)
else
    sat = satellite(sc, alm);                  % SEM almanac or TLE: satcom detects the format
end
names = sat.Name(:);
prn = prnNumber(names);
[G.epoch, healthy] = almanacInfo(sat);
G.age = days(T - G.epoch);
if isnan(G.age)
    fprintf('\n'); warning('ff_gnss_margin:almanac', 'Almanac epoch unknown: check that it is recent.');
elseif abs(G.age) > 30
    fprintf('\n');
    warning('ff_gnss_margin:almanac', ['The almanac is %.0f days old (%s): satellite positions and PRNs ' ...
        'are only indicative.\nUse a current SEM almanac (navcen.uscg.gov, GPS almanacs), a GPS TLE file ' ...
        '(celestrak.org) or a RINEX navigation file with ''Almanac''.'], G.age, char(G.epoch, 'yyyy-MM-dd'));
end

gs = groundStation(sc, opt.Lat, opt.Lon, 'Altitude', opt.Alt, 'MinElevationAngle', opt.MaskAngle, 'Name', "Site");
[az, el, rng] = aer(gs, sat, T);
vis = el(:) >= opt.MaskAngle & healthy;
ant = A.ant; Gmax = A.Gmax; G.antennaName = A.name;
iso = phased.IsotropicAntennaElement;
sats = table(prn(vis), names(vis), az(vis), el(vis), rng(vis)/1e3, ...
    'VariableNames', {'PRN', 'Satellite', 'Az_deg', 'El_deg', 'Range_km'});
G.crossCheck = NaN;
if any(vis)
    % EIRP from an isotropic satellite antenna; receiver G/T at the antenna's peak
    % gain. ebno uses the pattern gain towards each satellite relative to the peak.
    tx = transmitter(sat(vis), 'Frequency', f0, 'Power', opt.EIRP, 'BitRate', Rb/1e6, ...
        'Antenna', iso, 'SystemLoss', 0);
    Gdir = 0;                                   % gain added outside satcom (fallback / table)
    rxArgs = {'SystemLoss', opt.PolLoss, 'PreReceiverLoss', opt.PolLoss, 'Name', "GNSS receiver"};
    if isempty(A.gainFcn)
        try
            % Roll 180 deg: the antenna's +z (boresight) points down in the station's
            % NED frame by default; this turns it to the zenith.
            rx = receiver(gs, 'Antenna', ant, 'MountingAngles', [0; 0; 180], ...
                'GainToNoiseTemperatureRatio', Gmax - 10*log10(Tsys), rxArgs{:});
        catch err
            warning('ff_gnss_margin:antenna', ['The receiver does not accept this antenna (%s): its ' ...
                'pattern gain towards each satellite is added to an isotropic receiver instead.'], err.message);
            A.gainFcn = @(a, e) arrayfun(@(x, y) pattern(ant, f0, mod(x + 180, 360) - 180, y), a, e);
            G.antennaName = G.antennaName + " (pattern added)";
        end
    end
    if ~isempty(A.gainFcn)
        rx = receiver(gs, 'Antenna', iso, 'GainToNoiseTemperatureRatio', -10*log10(Tsys), rxArgs{:});
        Gdir = reshape(A.gainFcn(sats.Az_deg, sats.El_deg), [], 1);
    end
    lnk = link(tx, rx);
    cn0 = reshape(ebno(lnk, T), [], 1) + 10*log10(Rb) + Gdir;
    [~, pri] = sigstrength(lnk, T);            % power at the receiver input (incl. antenna gain)
    pri = reshape(pri, [], 1) + Gdir;
    G.crossCheck = max(abs(cn0 - (pri - 10*log10(kB*Tsys))));
    sats.Pr_dBW = pri;
    sats.CN0_Nominal_dBHz = cn0;
else
    sats.Pr_dBW = zeros(0, 1);
    sats.CN0_Nominal_dBHz = zeros(0, 1);
end
G.sats = sortrows(sats, 'El_deg', 'descend');

tl = accessIntervals(access(gs, sat(healthy)));
tl.StartTime.TimeZone = T.TimeZone;
tl.EndTime.TimeZone = T.TimeZone;
G.timeline = tl;
[~, elAll] = aer(gs, sat(healthy));
nv = sum(elAll >= opt.MaskAngle, 1);
G.nVisible = [min(nv) max(nv)];
G.scenario = sc;
G.antenna = ant;
G.Gmax = Gmax;
G.navCheck = navCrossCheck(alm, T, opt, prn, az, el);
end

function d = navCrossCheck(alm, T, opt, prn, az, el)
% Largest az/el difference (deg) between satcom aer and the Navigation Toolbox
% (semread + gnssconstellation + lookangles) for the satellites above the mask.
% Both use the almanac orbit model, so this checks PRNs, time zone and site.
d = NaN;
if ~isempty(regexpi(alm, '\.(rnx|\d\dn|nav|tle)$', 'once')), return, end   % SEM only (TLE: semread fails -> NaN)
try
    nav = semread(alm);
    [p, ~, id] = gnssconstellation(T, nav, 'GNSSFileType', "SEM");
    [az2, el2] = lookangles([opt.Lat opt.Lon opt.Alt], p, opt.MaskAngle);
    [~, ia, ib] = intersect(prn(:), id(:));
    up = el(ia) >= opt.MaskAngle;
    if any(up)
        dAz = mod(az(ia(up)) - az2(ib(up)) + 180, 360) - 180;
        d = max(abs([el(ia(up)) - el2(ib(up)); dAz.*cosd(el(ia(up)))]));
    end
catch
end
end

function p = prnNumber(names)
p = nan(numel(names), 1);
for k = 1:numel(names)
    tok = regexp(names(k), 'PRN\D*(\d+)', 'tokens', 'once');
    if ~isempty(tok), p(k) = str2double(tok(1)); end
end
end

function [epoch, healthy] = almanacInfo(sat)
% Almanac reference time (GPS week + time of applicability, or TLE epoch) and health.
healthy = true(numel(sat), 1);
epoch = NaT('TimeZone', 'UTC');
for k = 1:numel(sat)
    try
        oe = orbitalElements(sat(k));
    catch
        continue
    end
    if isfield(oe, 'SatelliteHealth'), healthy(k) = oe.SatelliteHealth == 0; end
    if isnat(epoch)
        if isfield(oe, 'GPSWeekNumber')
            epoch = datetime(1980, 1, 6, 'TimeZone', 'UTC') + days(7*oe.GPSWeekNumber) + seconds(oe.GPSTimeOfApplicability);
        elseif isfield(oe, 'Epoch')
            epoch = oe.Epoch;
            if isempty(epoch.TimeZone), epoch.TimeZone = 'UTC'; end
        end
    end
end
end

function A = rxAntenna(a, f0)
% Receive antenna, its peak gain (dBi) and, for a gain table, the gain towards
% (az, el). The L1 patch is designed and solved once per MATLAB session (~15 s).
persistent patch patchGmax
A = struct('ant', phased.IsotropicAntennaElement, 'Gmax', 0, 'name', "isotropic 0 dBi", 'gainFcn', []);
if isnumeric(a)
    if size(a, 2) ~= 2 || size(a, 1) < 2 || any(diff(a(:, 1)) <= 0)
        error('A gain table ''Antenna'' is [el_deg gain_dBi] with increasing elevations (>= 2 rows).');
    end
    A.Gmax = max(a(:, 2));
    A.gainFcn = @(~, e) interp1(a(:, 1), a(:, 2), e, 'linear', 'extrap');
    A.name = sprintf("datasheet gain table (%g dBi at %g deg ... %g dBi at %g deg)", a(1, 2), a(1, 1), a(end, 2), a(end, 1));
elseif ischar(a) || isstring(a)
    switch lower(string(a))
        case "patch"
            if isempty(patch)
                patch = design(patchMicrostrip, f0);
                g = pattern(patch, f0, -180:5:180, 0:2:90);
                patchGmax = max(g(:));
            end
            A.ant = patch; A.Gmax = patchGmax; A.name = "patchMicrostrip designed for L1";
        case "isotropic"
        otherwise
            error('''Antenna'' must be "patch", "isotropic", an antenna object or a gain table.');
    end
else
    A.ant = a; A.name = string(class(a));
    try
        g = pattern(a, f0, -180:5:180, 0:2:90);
        A.Gmax = max(g(:));
    catch err
        warning('ff_gnss_margin:antenna', 'Could not compute the pattern of the %s antenna (%s): peak gain 0 dBi assumed.', A.name, err.message);
    end
end
end

%% ===== measurement ======================================================

function m = measureL1(ip, opt, f0)
fprintf('Connecting to %s:%d ...\n', ip, opt.Port);
ff = FieldFox(ip, opt.Port);
fprintf('Connected: %s\n', ff.IDN);
if ~ff.hasOption("233"), error('ff_gnss_margin:noSA', 'This FieldFox has no spectrum analyzer (option 233).'); end
ff.takeControl();                           % saves the full state; restored on every exit
restore = onCleanup(@() releaseFF(ff));
args = {'Center', f0, 'Span', opt.Span, 'Points', opt.Points, 'RBW', opt.RBW, 'VBW', opt.RBW, ...
    'AttenAuto', false, 'Atten', 0, 'RefLevel', -40, 'Detector', "AVER", 'TraceType', "CLRW", 'Average', 1};
if ff.hasOption("235")
    args = [args, {'Preamp', true}];
else
    warning('ff_gnss_margin:preamp', 'No preamplifier (option 235): the FieldFox noise floor is much higher.');
end
ff.saConfigure(args{:});
% Tracking-generator source off (no RF output; the saved state restores it).
% SOUR:ENAB is in the programming guide, not yet verified on firmware A.12.45.
ff.write("SOUR:ENAB OFF");
e = ff.errors();
if ~isempty(e), warning('ff_gnss_margin:source', 'SOUR:ENAB OFF not accepted (%s).', strjoin(e, '; ')); end
s = ff.saSettings();
m = struct('f', s.f, 'rbw', s.rbw, 'settings', rmfield(s, 'f'), 'idn', ff.IDN, 'source', "FieldFox " + ip);
fprintf('SA: %.3f-%.3f MHz, %d points, RBW %g kHz, preamp %d, atten %g dB, detector %s\n', ...
    s.start/1e6, s.stop/1e6, s.points, s.rbw/1e3, s.preamp, s.atten, s.detector);
checkSpectrum(m.f, m.rbw, f0, 1.023e6);
if opt.Mode == "offon"
    ask(opt, 'Switch the solar inverter OFF.');
    [m.Poff, m.toff, o1] = takeSweeps(ff, opt.Sweeps, numel(m.f), "OFF");
    ask(opt, 'Switch the solar inverter ON and wait until it delivers power.');
    [m.Pon, m.ton, o2] = takeSweeps(ff, opt.Sweeps, numel(m.f), "ON");
else
    [m.Poff, m.toff, o1] = takeSweeps(ff, opt.Sweeps, numel(m.f), "L1");
    m.Pon = []; m.ton = NaT(0, 1); o2 = NaN;
end
m.overRange = [o1 o2];
e = ff.errors();
if ~isempty(e), warning('ff_gnss_margin:instrument', 'FieldFox reported: %s', strjoin(e, '; ')); end
clear restore                               % release now (also runs on errors)
end

function [P, t, ovr] = takeSweeps(ff, n, nPts, label)
P = zeros(n, nPts); t = NaT(n, 1);
fprintf('Measuring %s: %d sweeps ', label, n);
for k = 1:n
    ff.sweep();
    y = ff.saTrace();
    if numel(y) ~= nPts
        error('ff_gnss_margin:trace', 'Trace has %d points, expected %d (timeout or wrong reply).', numel(y), nPts);
    end
    P(k, :) = y;
    t(k) = datetime('now');
    if mod(k, 10) == 0, fprintf('.'); end
end
fprintf(' done\n');
ovr = overRange(ff);
if ovr == 1
    warning('ff_gnss_margin:overrange', ['%s: ADC over-range. A strong signal somewhere overloads the ' ...
        'preamp, so these levels are not reliable. Add a GNSS band-pass filter.'], label);
end
end

function ovr = overRange(ff)
% ADC over-range flag of the last sweep (1 = overload). CALC:MEAS:WAOR? is in the
% 2022 programming guide (p.234) but not verified on firmware A.12.45: NaN if
% unsupported (no reply -> 3 s timeout, error -113 in the queue).
ovr = NaN;
try
    v = ff.queryNum("CALC:MEAS:WAOR?", 3);
    if isscalar(v), ovr = v; end
catch
end
e = ff.errors();
if ~isempty(e)
    ovr = NaN;
    e = e(~startsWith(e, "-113"));             % "Undefined header" = WAOR unsupported
    if ~isempty(e), warning('ff_gnss_margin:instrument', 'FieldFox reported: %s', strjoin(e, '; ')); end
end
end

function ask(opt, msg)
if opt.Prompt
    input(sprintf('%s Press Enter to measure... ', msg), 's');
else
    fprintf('%s (Prompt=false: measuring now)\n', msg);
end
end

function releaseFF(ff)
try
    ff.release();
    fprintf('FieldFox settings restored.\n');
catch err
    warning('ff_gnss_margin:release', 'Could not restore the FieldFox: %s', err.message);
end
end

function files = saveRecordings(m, outDir, stamp)
% Same layout as fieldfox_live recordings: f, P, t, settings, idn.
if isempty(m.Pon), sets = {"", m.Poff, m.toff};
else, sets = {"_off", m.Poff, m.toff; "_on", m.Pon, m.ton};
end
files = strings(0, 1);
for k = 1:size(sets, 1)
    out = struct('f', m.f, 'P', sets{k, 2}, 't', sets{k, 3}, 'settings', m.settings, 'idn', m.idn);
    files(end + 1, 1) = fullfile(outDir, "gnss_l1_" + stamp + sets{k, 1} + ".mat"); %#ok<AGROW>
    save(files(end), '-struct', 'out');
end
end

function m = loadSpectrum(spec, specOn, rbw)
a = readSpectrum(spec, rbw);
m = struct('f', a.f, 'rbw', a.rbw, 'settings', a.settings, 'idn', a.idn, 'source', a.source, ...
    'Poff', a.P, 'toff', a.t, 'Pon', a.Pon, 'ton', NaT(0, 1), 'overRange', [NaN NaN]);
if ~isempty(specOn)
    b = readSpectrum(specOn, rbw);
    if numel(b.f) ~= numel(a.f) || max(abs(b.f - a.f)) > 1
        error('''Spectrum'' and ''SpectrumOn'' have different frequency axes.');
    end
    m.Pon = b.P; m.ton = b.t;
end
fprintf('Spectrum from %s: %.3f-%.3f MHz, %d points, RBW %g kHz\n', m.source, m.f(1)/1e6, m.f(end)/1e6, numel(m.f), m.rbw/1e3);
checkSpectrum(m.f, m.rbw, 1575.42e6, 1.023e6);
end

function s = readSpectrum(x, rbw)
s.source = "struct";
if ischar(x) || isstring(x), s.source = string(x); x = load(x); end
if ~isfield(x, 'f') || ~isfield(x, 'P'), error('''Spectrum'' needs f (Hz) and P (dBm, sweeps x points).'); end
s.f = x.f(:).';
s.P = rowIfVector(x.P);
if size(s.P, 2) ~= numel(s.f), error('''Spectrum'': P must have numel(f) columns.'); end
s.Pon = [];
if isfield(x, 'P_on') && ~isempty(x.P_on), s.Pon = rowIfVector(x.P_on); end
s.rbw = NaN;
if isfield(x, 'rbw'), s.rbw = x.rbw;
elseif isfield(x, 'settings') && isfield(x.settings, 'rbw'), s.rbw = x.settings.rbw;
end
if isnan(s.rbw)
    warning('ff_gnss_margin:rbw', 'The spectrum has no RBW: assuming %g Hz (''RBW'' option).', rbw);
    s.rbw = rbw;
end
s.t = NaT(0, 1); if isfield(x, 't'), s.t = x.t; end
s.settings = struct(); if isfield(x, 'settings'), s.settings = x.settings; end
s.idn = ""; if isfield(x, 'idn'), s.idn = string(x.idn); end
end

function P = rowIfVector(P)
if isvector(P), P = P(:).'; end
end

%% ===== interference =====================================================

function N = l1Noise(m, opt, C, f0, Rc, kB)
% Noise density at the FieldFox input (averaged in mW), the extra noise I0
% referred to the antenna and the C/N0 degradation 10*log10(1 + I0/N0).
N = struct('measured', false, 'mode', "none", 'source', "", 'f', [], 'psdOff_dBmHz', [], 'psdOn_dBmHz', [], ...
    'main_dBmHz', [NaN NaN], 'wide_dBmHz', [NaN NaN], 'mainPower_dBm', [NaN NaN], 'widePower_dBm', [NaN NaN], ...
    'rise_dB', NaN, 'expected_dBmHz', C.expected_dBmHz, 'vsExpected_dB', NaN, 'I0_dBmHz', -Inf, ...
    'I0overN0_dB', -Inf, 'degradation_dB', 0, 'detectLimit_dB', NaN, 'significant', false, ...
    'upperBound', false, 'sweeps', [0 0], 'overRange', [NaN NaN], 'spurs', table());
if isempty(m), return, end
f = m.f(:).';                              % checked by checkSpectrum when loaded/measured
main = abs(f - f0) <= Rc;
wide = abs(f - f0) <= 10e6;
nbw = 1.06*m.rbw;                         % noise bandwidth of the Gaussian RBW filter
w = sinc((f - f0)/Rc).^2; w = w/sum(w);    % C/A power spectrum, normalised over the span
bwWide = min(f(end), f0 + 10e6) - max(f(1), f0 - 10e6);
hasOn = ~isempty(m.Pon);
S1 = 10.^(m.Poff/10)/nbw;                  % mW/Hz, sweeps x points
x1 = S1*w.';                               % C/A-weighted density of each sweep
N.measured = true;
N.source = m.source;
N.mode = "single"; if hasOn, N.mode = "offon"; end
N.f = f;
N.psdOff_dBmHz = 10*log10(mean(S1, 1));
N.sweeps(1) = size(S1, 1);
N.overRange = m.overRange;
lin1 = mean(S1, 1);
N.main_dBmHz(1) = 10*log10(mean(lin1(main)));
N.wide_dBmHz(1) = 10*log10(mean(lin1(wide)));
Sexp = 10^(C.expected_dBmHz/10);
N.vsExpected_dB = 10*log10(mean(x1)/Sexp);
% No analyzer can show less noise than kT (-174 dBm/Hz) at its input. A lower
% reading means an amplitude offset was active or another program (Data Link,
% front panel) changed the settings during the measurement.
if N.main_dBmHz(1) < -176
    warning('ff_gnss_margin:belowThermal', ['Measured %.1f dBm/Hz is below thermal noise (-174 dBm/Hz): ' ...
        'the FieldFox settings were probably changed during the measurement (another program or the ' ...
        'front panel). Repeat the measurement with nothing else connected to the FieldFox.'], N.main_dBmHz(1));
end
if hasOn
    S2 = 10.^(m.Pon/10)/nbw;
    x2 = S2*w.';
    lin2 = mean(S2, 1);
    N.psdOn_dBmHz = 10*log10(lin2);
    N.sweeps(2) = size(S2, 1);
    N.main_dBmHz(2) = 10*log10(mean(lin2(main)));
    N.wide_dBmHz(2) = 10*log10(mean(lin2(wide)));
    dS = mean(x2) - mean(x1);
    se = sqrt(var(x1)/numel(x1) + var(x2)/numel(x2));
    N.rise_dB = 10*log10(mean(x2)/mean(x1));
    ylast = N.psdOn_dBmHz;
else
    % No OFF reference: everything above the expected antenna+chain noise counts
    % (FieldFox noise included), so this is an upper bound.
    dS = mean(x1) - Sexp;
    se = std(x1)/sqrt(numel(x1));
    N.upperBound = true;
    ylast = N.psdOff_dBmHz;
end
N.mainPower_dBm = N.main_dBmHz + 10*log10(2*Rc);
N.widePower_dBm = N.wide_dBmHz + 10*log10(bwWide);
dSmin = max(2*se, (10^(opt.MinRise/10) - 1)*mean(x1));
G = 10^(C.Gmeas_dB/10);
N0 = kB*C.Tsys_K*1e3;                      % mW/Hz at the antenna
I0 = max(dS, 0)/G;
N.I0_dBmHz = 10*log10(I0);
N.I0overN0_dB = 10*log10(I0/N0);
N.degradation_dB = 10*log10(1 + I0/N0);
N.detectLimit_dB = 10*log10(1 + dSmin/G/N0);
N.significant = dS > dSmin;

% Narrow spurs (CW-like) in the last measured state: trace level (dBm in the
% RBW at the FieldFox) >= 6 dB above the median floor
ytr = ylast + 10*log10(nbw);
pk = []; loc = [];
if max(ytr) > median(ytr) + 6
    [pk, loc] = findpeaks(ytr, f, 'MinPeakProminence', 6, 'MinPeakHeight', median(ytr) + 6, ...
        'SortStr', 'descend', 'NPeaks', 10);
end
rise = nan(numel(pk), 1);
if hasOn, rise = pk(:) - interp1(f, N.psdOff_dBmHz + 10*log10(nbw), loc(:)); end
N.spurs = table(loc(:)/1e6, pk(:), rise, abs(loc(:) - f0) <= Rc, strings(numel(pk), 1), ...
    'VariableNames', {'Freq_MHz', 'Level_dBm', 'Rise_dB', 'InMainLobe', 'Technology'});
if ~isempty(pk) && exist('ff_signal_id', 'file')
    try
        id = ff_signal_id(f, ytr, m.rbw, 'Detector', "average", 'Visible', 'off', 'Prompt', false);
        for k = 1:numel(pk)
            [d, j] = min(abs(id.Freq_Hz - loc(k)));
            if d <= 2*m.rbw, N.spurs.Technology(k) = string(id.Technology(j)); end
        end
    catch
    end
end
end

function checkSpectrum(f, rbw, f0, Rc)
if f(1) > f0 - Rc || f(end) < f0 + Rc
    error('The spectrum (%.3f-%.3f MHz) must cover the C/A main lobe %.3f-%.3f MHz.', ...
        f(1)/1e6, f(end)/1e6, (f0 - Rc)/1e6, (f0 + Rc)/1e6);
end
n = nnz(abs(f - f0) <= Rc);
if n < 5
    error(['Only %d spectrum points fall in the C/A main lobe (+-1.023 MHz): measure a span of ' ...
        '~20 MHz around 1575.42 MHz.'], n);
end
if rbw > Rc/3
    warning('ff_gnss_margin:rbw', 'RBW %g kHz smears the C/A main lobe: use 100 kHz or less.', rbw/1e3);
end
end

function M = margins(sats, opt, C, N, kB)
% How much the noise may rise before satellites are lost. A rise of r dB in the
% total noise density lowers every C/N0 by r dB.
cn0 = sort(sats.CN0_Nominal_dBHz, 'descend');
on = sats.CN0_Interf_dBHz;
N0 = kB*C.Tsys_K*1e3;
toI0 = @(r) 10*log10(N0*(10.^(r/10) - 1)) + 0./(r > 0);   % dBm/Hz at the antenna, NaN if r <= 0
M.nVisible = numel(cn0);
M.nTrack = [nnz(cn0 >= opt.TrackThreshold) nnz(on >= opt.TrackThreshold)];
M.nAcq = [nnz(cn0 >= opt.AcqThreshold) nnz(on >= opt.AcqThreshold)];
M.riseTrackAll_dB = NaN; M.riseFixAcq_dB = NaN; M.riseFixTrack_dB = NaN;
if ~isempty(cn0), M.riseTrackAll_dB = cn0(end) - opt.TrackThreshold; end
if numel(cn0) >= 4                         % a position fix needs 4 satellites
    M.riseFixAcq_dB = cn0(4) - opt.AcqThreshold;
    M.riseFixTrack_dB = cn0(4) - opt.TrackThreshold;
end
M.I0maxTrackAll_dBmHz = toI0(M.riseTrackAll_dB);
M.I0maxFixAcq_dBmHz = toI0(M.riseFixAcq_dB);
M.I0maxFixTrack_dBmHz = toI0(M.riseFixTrack_dB);
M.measuredRise_dB = N.degradation_dB;
end

%% ===== report and figure ================================================

function report(opt, T, G, C, N, M, sats)
fprintf('\n===== GPS L1 interference margin =====\n');
fprintf('Site %.5f, %.5f, %g m   Time %s\n', opt.Lat, opt.Lon, opt.Alt, string(T, 'yyyy-MM-dd HH:mm:ss z'));
fprintf('Almanac %s (epoch %s, %.0f days before the analysis time)\n', string(opt.Almanac), ...
    char(G.epoch, 'yyyy-MM-dd'), G.age);
fprintf('Antenna: %s, peak %.1f dBi, pointing at the zenith; polarisation loss %g dB; EIRP %g dBW\n', ...
    G.antennaName, G.Gmax, opt.PolLoss, opt.EIRP);
if ~isnan(G.navCheck)
    fprintf('Navigation Toolbox (gnssconstellation + lookangles) agrees with satcom aer within %.3f deg\n', G.navCheck);
end
b = C.budget;
fprintf('Receiver chain (rfbudget, cumulative from the antenna; power = strongest satellite):\n');
for k = 1:numel(C.stages)
    fprintf('  %-9s gain %6.1f dB   NF %5.2f dB   signal %7.1f dBm\n', C.stages(k), ...
        b.TransducerGain(k), b.NF(k), b.OutputPower(k));
end
fprintf('System NF %.2f dB, Tsys = %g K antenna + %.0f K = %.0f K  ->  N0 = %.1f dBW/Hz (%.1f dBm/Hz)\n', ...
    C.NF_dB, C.Tant_K, C.Tsys_K - C.Tant_K, C.Tsys_K, C.N0_dBWHz, C.N0_dBWHz + 30);

if N.measured
    fprintf('\nL1 measurement (%s, %s): FieldFox at the %s, gain antenna->FieldFox %.1f dB\n', ...
        N.mode, N.source, opt.MeasuredAt, C.Gmeas_dB);
    lab = ["OFF" "ON"]; if N.mode == "single", lab = ["measured" "-"]; end
    fprintf('  Noise density, C/A main lobe +-1.023 MHz: %s %.1f dBm/Hz', lab(1), N.main_dBmHz(1));
    if N.mode == "offon", fprintf(', %s %.1f dBm/Hz (%+.2f dB)', lab(2), N.main_dBmHz(2), diff(N.main_dBmHz)); end
    fprintf('\n  Noise density, +-10 MHz:                  %s %.1f dBm/Hz', lab(1), N.wide_dBmHz(1));
    if N.mode == "offon", fprintf(', %s %.1f dBm/Hz (%+.2f dB)', lab(2), N.wide_dBmHz(2), diff(N.wide_dBmHz)); end
    fprintf('\n  Power in the main lobe %.1f dBm, in +-10 MHz %.1f dBm (%s)\n', N.mainPower_dBm(1), N.widePower_dBm(1), lab(1));
    fprintf('  Expected noise of antenna + chain at the FieldFox: %.1f dBm/Hz; measured %s is %+.1f dB from it\n', ...
        N.expected_dBmHz, lab(1), N.vsExpected_dB);
    if N.vsExpected_dB < -3
        fprintf(['  -> Less noise than the antenna+chain must produce: is the antenna/LNA connected and powered\n' ...
            '     (bias tee)? Are ''LNAGain''/''Cable'' right? (Nothing connected = FieldFox floor only.)\n']);
    elseif N.vsExpected_dB > 6
        fprintf(['  -> More noise than expected: the FieldFox''s own noise floor dominates (more LNA gain\n' ...
            '     helps) or interference is already present. Weak external noise is hidden.\n']);
    end
    if any(N.overRange == 1), fprintf('  !! ADC over-range during the measurement: levels unreliable.\n'); end
    what = "Extra noise ON-OFF"; if N.upperBound, what = "Noise above expected (UPPER BOUND)"; end
    fprintf('  %s, C/A-weighted, at the antenna: %.1f dBm/Hz = %.1f dB re N0\n', what, N.I0_dBmHz, N.I0overN0_dB);
    fprintf('  -> C/N0 degradation %.2f dB. Detection limit of this measurement: %.2f dB', N.degradation_dB, N.detectLimit_dB);
    if N.significant, fprintf(' (significant)\n'); else, fprintf(' (NOT significant: no measurable effect)\n'); end
    if ~isempty(N.spurs)
        fprintf('  Narrow spurs (>= 6 dB above the floor):\n');
        disp(N.spurs);
    end
else
    fprintf('\nNo L1 measurement (give an IP address or ''Spectrum''): nominal C/N0 only.\n');
end

fprintf('\nSatellites above %g deg at %s (%d; %d-%d over the next %g h):\n', opt.MaskAngle, ...
    string(T, 'HH:mm'), height(sats), G.nVisible(1), G.nVisible(2), opt.TimelineHours);
disp(sats(:, {'PRN', 'Az_deg', 'El_deg', 'Range_km', 'Pr_dBW', 'CN0_Nominal_dBHz', ...
    'CN0_Interf_dBHz', 'MarginTrack_dB', 'MarginAcq_dB'}));
if G.crossCheck > 0.05
    warning('ff_gnss_margin:check', 'ebno and sigstrength C/N0 differ by %.2f dB.', G.crossCheck);
end
fprintf('Trackable (>= %g dB-Hz): %d -> %d satellites; acquirable (>= %g dB-Hz): %d -> %d\n', ...
    opt.TrackThreshold, M.nTrack(1), M.nTrack(2), opt.AcqThreshold, M.nAcq(1), M.nAcq(2));
fprintf('Allowed noise rise: %.1f dB keeps all tracked, %.1f dB keeps 4 acquirable (fix from cold start)\n', ...
    M.riseTrackAll_dB, M.riseFixAcq_dB);
fprintf('  = interference up to %.1f / %.1f dBm/Hz at the antenna (%.1f / %.1f dBm/Hz at the FieldFox).\n', ...
    M.I0maxTrackAll_dBmHz, M.I0maxFixAcq_dBmHz, M.I0maxTrackAll_dBmHz + C.Gmeas_dB, M.I0maxFixAcq_dBmHz + C.Gmeas_dB);
if N.measured
    fprintf('Measured degradation %.2f dB -> remaining margin %.1f dB (all tracked), %.1f dB (4 acquirable)\n', ...
        N.degradation_dB, M.riseTrackAll_dB - N.degradation_dB, M.riseFixAcq_dB - N.degradation_dB);
end
end

function fig = plotResults(opt, T, G, C, N, M, sats, f0, Rc)
fig = figure('Name', 'GNSS L1 interference margin', 'NumberTitle', 'off', 'Color', 'w', ...
    'Position', [60 60 1400 860], 'Visible', opt.Visible);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('GPS L1 interference margin  |  %.4f, %.4f  |  %s', opt.Lat, opt.Lon, ...
    string(T, 'yyyy-MM-dd HH:mm z')), 'FontWeight', 'bold');
lab = "G" + compose('%02d', sats.PRN);
lab(isnan(sats.PRN)) = sats.Satellite(isnan(sats.PRN));
lab = matlab.lang.makeUniqueStrings(lab);
onName = "Inverter ON"; offName = "Inverter OFF";
if N.mode == "single", offName = "Measured"; onName = "With measured noise (upper bound)";
elseif ~N.measured, offName = "Nominal"; onName = "No measurement";
end

% Sky plot coloured by C/N0 with interference
edges = [-Inf opt.TrackThreshold opt.AcqThreshold max(45, opt.AcqThreshold + 5) Inf];
names = {sprintf('< %g dB-Hz (lost)', edges(2)), sprintf('%g-%g (track only)', edges(2), edges(3)), ...
    sprintf('%g-%g dB-Hz', edges(3), edges(4)), sprintf('>= %g dB-Hz', edges(4))};
cols = [0.85 0.1 0.1; 1 0.55 0; 0.85 0.75 0; 0.1 0.65 0.2];
if ~isempty(sats)
    cat = discretize(sats.CN0_Interf_dBHz, edges, 'categorical', names);
    used = ismember(names, string(cat));
    sp = skyplot(tl, sats.Az_deg, sats.El_deg, lab, 'GroupData', removecats(cat), ...
        'MaskElevation', opt.MaskAngle);
    sp.ColorOrder = cols(used, :);
    sp.Layout.Tile = 1;
    legend(sp, names(used), 'Location', 'southoutside');
else
    ax = nexttile(tl, 1); axis(ax, 'off');
    text(ax, 0.5, 0.5, 'No satellite above the mask', 'HorizontalAlignment', 'center');
end

% L1 spectrum as noise density at the FieldFox input
ax = nexttile(tl, 2); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
if N.measured
    fM = N.f/1e6;
    g = 10^(C.Gmeas_dB/10);
    xregion(ax, (f0 - Rc)/1e6, (f0 + Rc)/1e6, 'FaceColor', [1 0.85 0.3], 'FaceAlpha', 0.25);
    h = plot(ax, fM, N.psdOff_dBmHz, 'Color', [0.2 0.45 0.9], 'LineWidth', 1.2); leg = offName;
    if ~isempty(N.psdOn_dBmHz)
        h(end + 1) = plot(ax, fM, N.psdOn_dBmHz, 'Color', [0.9 0.2 0.2], 'LineWidth', 1.2); leg(end + 1) = onName;
    end
    h(end + 1) = yline(ax, C.expected_dBmHz, '--', 'Color', [0 0.6 0], 'LineWidth', 1.2);
    leg(end + 1) = "Expected antenna + chain noise";
    if ~isempty(sats)
        ca = 10*log10(10^(max(sats.Pr_dBW)/10)*1e3*g/Rc*sinc((N.f - f0)/Rc).^2);
        h(end + 1) = plot(ax, fM, ca, 'k', 'LineWidth', 1); leg(end + 1) = "C/A signal of the strongest satellite";
    end
    base = 10.^(N.psdOff_dBmHz/10);
    if isfinite(M.I0maxFixAcq_dBmHz)
        h(end + 1) = plot(ax, fM, 10*log10(base + g*10^(M.I0maxFixAcq_dBmHz/10)), 'm:', 'LineWidth', 1.4);
        leg(end + 1) = "Limit: 4 satellites still acquirable";
    end
    if isfinite(M.I0maxTrackAll_dBmHz)
        h(end + 1) = plot(ax, fM, 10*log10(base + g*10^(M.I0maxTrackAll_dBmHz/10)), 'm--', 'LineWidth', 1);
        leg(end + 1) = "Limit: all satellites still tracked";
    end
    legend(ax, h, leg, 'Location', 'southoutside', 'NumColumns', 2);
    xlim(ax, fM([1 end]));
    yl = [min([N.psdOff_dBmHz(:); C.expected_dBmHz]) - 20, ...      % keep the sinc^2 nulls off-scale
        max([N.psdOff_dBmHz(:); N.psdOn_dBmHz(:); C.expected_dBmHz]) + 10];
    if isfinite(M.I0maxTrackAll_dBmHz), yl(2) = max(yl(2), M.I0maxTrackAll_dBmHz + C.Gmeas_dB + 5); end
    if ~isempty(sats), yl(1) = min(yl(1), 10*log10(10^(max(sats.Pr_dBW)/10)*1e3*g/Rc) - 10); end
    ylim(ax, yl);
    xlabel(ax, 'Frequency (MHz)'); ylabel(ax, 'Noise density at the FieldFox (dBm/Hz)');
    title(ax, sprintf('L1 spectrum (RBW noise BW 1.06xRBW), C/A main lobe shaded; degradation %.2f dB', N.degradation_dB));
else
    axis(ax, 'off');
    text(ax, 0.5, 0.5, {'No L1 measurement', 'give an IP address or ''Spectrum'''}, 'HorizontalAlignment', 'center');
end

% C/N0 bars
ax = nexttile(tl, 3); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
if ~isempty(sats)
    x = categorical(lab, lab);
    if N.measured
        bb = bar(ax, x, [sats.CN0_Nominal_dBHz sats.CN0_Interf_dBHz], 'grouped');
        bb(2).FaceColor = [0.9 0.2 0.2];
    else
        bb = bar(ax, x, sats.CN0_Nominal_dBHz); onName = strings(1, 0);
    end
    bb(1).FaceColor = [0.2 0.45 0.9];
    yline(ax, opt.TrackThreshold, 'r--', sprintf('tracking %g', opt.TrackThreshold), 'LabelHorizontalAlignment', 'left');
    yline(ax, opt.AcqThreshold, '--', sprintf('acquisition %g', opt.AcqThreshold), 'Color', [0.9 0.5 0], ...
        'LabelHorizontalAlignment', 'left');
    legend(ax, bb, [offName, onName], 'Location', 'southoutside', 'NumColumns', 2);
    ylim(ax, [min([sats.CN0_Interf_dBHz; opt.TrackThreshold]) - 5, max(sats.CN0_Nominal_dBHz) + 5]);
end
ylabel(ax, 'C/N0 (dB-Hz)');
title(ax, 'C/N0 per satellite (highest elevation first)');

% Visibility timeline
ax = nexttile(tl, 4); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
tlb = G.timeline;
tg = unique(string(tlb.Target));
[~, order] = sort(prnNumber(tg)); tg = tg(order);
plot(ax, [T T + hours(opt.TimelineHours)], [NaN NaN]);    % datetime x axis even if empty
for k = 1:numel(tg)
    r = tlb(string(tlb.Target) == tg(k), :);
    for j = 1:height(r)
        plot(ax, [r.StartTime(j) r.EndTime(j)], [k k], '-', 'LineWidth', 5, 'Color', [0.2 0.45 0.9]);
    end
end
yticks(ax, 1:numel(tg)); yticklabels(ax, erase(tg, "PRN:"));
ylim(ax, [0 numel(tg) + 1]);
xlim(ax, [T, T + hours(opt.TimelineHours)]);
ylabel(ax, 'Satellite (PRN)');
title(ax, sprintf('Above %g deg in the %g h from the analysis time (%d-%d at a time)', opt.MaskAngle, ...
    opt.TimelineHours, G.nVisible(1), G.nVisible(2)));
end
