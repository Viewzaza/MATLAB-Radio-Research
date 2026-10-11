function varargout = ff_satsig(cmd, varargin)
%FF_SATSIG  Satellite signal library and Doppler-aware AUTOMATIC satellite signal identification.
%
%   C = ff_satsig("catalog")                    satellite downlinks: NORAD, f0, BW, modulation, framing,
%                                                burst pattern, polarisation, status, source
%   C = ff_satsig("catalog", 'Refresh', true)   also refresh the SatNOGS DB subset (bulk JSON, cached 24 h)
%   D = ff_satsig("define", "OrigamiSat-1")      definition sheet + synthetic waveform (toolbox generators)
%                                                + spectrum, spectrogram, zero-span envelope, constellation
%   T = ff_satsig("overhead")                    catalog satellites above the horizon now (Bangkok), with
%                                                az/el/range and the Doppler-shifted frequency + rate
%   [T, info] = ff_satsig("identify", "rec.mat")                     identify a FieldFox recording
%   [T, info] = ff_satsig("identify", "10.90.36.200", 'Bands', [437.48e6 437.53e6], 'Duration', 600)
%   [T, info] = ff_satsig("identify", "10.90.36.200")                plan bands from what is overhead now
%   [T, info] = ff_satsig("demo")                simulated OrigamiSat-2 pass -> automatic identification
%   [T, info] = ff_satsig("demo", 'Demo', "both")   OrigamiSat-2 + KNACKSAT-2 (max-hold) passes
%   [T, info] = ff_satsig("demo", 'Demo', "weak")   OrigamiSat-2 17 dB weaker: found by Doppler stacking
%   ff_satsig("templates")                       rebuild support\satsig\satsig_templates.mat (~1 min)
%
% Example (headless, for a report):
%   [T, info] = ff_satsig("demo", 'Visible', 'off', 'OutDir', "satsig_demo");
%   D = ff_satsig("define", "OrigamiSat-1", 'Visible', 'off', 'OutDir', "satsig_define");
%
% IDENTIFY - how it works (no I/Q on this FieldFox, so nothing is demodulated):
%  1. Detection in every sweep: sigma-clipped noise floor, findpeaks on the raw and 3-point smoothed
%     trace, regions grown down to the floor, noise-subtracted power-weighted centroid, width at
%     X = min(10, SNR-3) dB below the peak, channel power.
%  2. Tracking: detections are chained over time (gate = resolution + max Doppler rate x gap; slope
%     prediction once a track has 3 points). CW and AFSK of one satellite share one track. After
%     scoring, tracks whose best hypothesis is the same satellite are merged and re-scored (beacons
%     every 60 s near TCA break the linear gate).
%  3. Candidates: every catalog downlink (and SatNOGS transmitter with 'IncludeSatNOGS') whose
%     frequency +- (max Doppler + OBW/2 + oscillator tolerance) covers the track. Constellation rows
%     (ORBCOMM, Iridium, GNSS) expand to every member satellite. Re-entered/decommissioned rows are
%     excluded (listed in the notes) unless 'IncludeInactive'.
%  4. Prediction: TLEs from CelesTrak (GROUP/CATNR, FORMAT=TLE, cached 2 h, no retry for 2 h after an
%     error; offline: cache or the shipped sample), satelliteScenario + satellite (SGP4) + states
%     (ECEF) -> elevation, range and Doppler factor -v_r/c (identical to dopplershift within 0.03 Hz).
%  5. Score (log-likelihood, normalised to a posterior over all hypotheses incl. "fixed frequency /
%     terrestrial" and "unknown satellite of class X"):
%       frequency : multi-sweep track -> lsqcurvefit of f(t) = (f0+df)(1 + d(t+dt)) (dt = TLE timing
%                   error, df = oscillator offset) with priors on dt and df; residual RMS vs a fixed
%                   (linear drift) carrier and a polynomial "unknown satellite"; single sweeps ->
%                   distance to the predicted Doppler-shifted carrier
%       visibility: fraction of detections with the satellite above the mask
%       bandwidth : measured width vs the RBW-smoothed template width at the same level (log-normal)
%       shape     : correlation of the Doppler-aligned mean spectrum with the template
%       level     : implausibly strong for the satellite (e.g. a "GNSS" signal above the floor) -> penalty
%       prior     : status (active / scheduled / event-based), source (catalog > SatNOGS > generic)
%       The fixed-carrier hypothesis allows a linear drift with a 2 Hz/s prior (terrestrial oscillators).
%  6. Weak-signal search (Doppler stacking, 'Stack'): for every catalog downlink that can fall in the
%     band and whose satellite is up, each sweep is shifted by the predicted Doppler f0*d(t+dt) and the
%     noise-normalised power is averaged (tracks attributed to other emitters masked). A signal that follows
%     that satellite's Doppler adds up, noise and fixed carriers do not: ~5*log10(N) dB gain for N
%     sweeps with the signal on. Reported with a robust z (threshold 'StackZ' + look-elsewhere term);
%     "inconclusive" when the Doppler changes by < 5 RBW over the stacked sweeps.
%  Every result is a HYPOTHESIS with a confidence, not a decode. A fixed-frequency signal inside a
%  satellite band is flagged as possible interference (EMI); any visible line in a GNSS band too.
%
% Name-Value options (defaults; which subcommands use them):
%   'Site'           [13.7563 100.5018 10]  [lat lon alt_m] (all). SET YOUR LAB LOCATION.
%   'SiteName'       "Bangkok"
%   'TimeZone'       "Asia/Bangkok"  zone of recordings whose t has no zone, and of printed times
%   'Time'           now (UTC)       overhead: time to evaluate (datetime; without zone = 'TimeZone')
%   'Offline'        false   true = never go online (TLE/SatNOGS caches, else the shipped samples)
%   'Plot'           true    figures;  'Visible' 'on'/'off' (hidden figures are closed once saved);  'OutDir' "" = do not save
%   'Prompt'         true    live: confirm before taking control of the instrument (dialog only
%                            in an interactive desktop session)
%   'IncludeInactive' false  also consider re-entered / decommissioned downlinks
%   'IncludeSatNOGS' false   catalog/overhead/identify: add SatNOGS DB transmitters (stale 'alive'
%                            flags exist: lower prior). Demo: true if the subset file exists.
%   'ElevationMask'  0       deg
%   'Quiet'          false   no console printing
%  catalog:  'Refresh' false (SatNOGS bulk download, max once per 24 h); 'Search' ""; 'Export' "" (CSV)
%  define:   'RxGain' 0 dBi; 'Seed' 1
%  overhead: 'Constellations' true (ORBCOMM, Iridium, GNSS members); 'Hours' 24 (pass list); 'Search' ""
%  identify: 'Threshold' 6 dB; 'FreqTolPPM' 5 + 'FreqTolHz' 500 (satellite oscillator 1-sigma);
%            'TimeTol' [] s (TLE timing 1-sigma; [] = 5 s + 3 s per day of TLE age); 'RxGain' 0 dB
%            (antenna + LNA - cable, for the level plausibility test); 'PriorFixed' 0.5;
%            'TrackGap' 180 s; 'MinDetections' [] (2 if >= 10 sweeps, else 1);
%            'Demo' "" | "origami" | "knacksat" | "both" | "weak"; 'Seed' 1;
%            'Stack' true (Doppler-stacking weak-signal search); 'StackZ' 6
%   live:    'Bands' [] (k x 2 [start stop] Hz); 'Targets' [] (catalog ids/names -> bands);
%            'Duration' 300 s; 'Interval' 0 s; 'RBW' [] (auto 1-3-10 <= span/(points-1));
%            'Points' 401; 'Preamp' true; 'Detector' "POS"; 'TraceType' "CLRW" | "MAXH";
%            'HoldTime' 10 s (MAXH window); 'Passive' false (read what the screen shows, change
%            nothing); 'Port' 5025; 'SaveTo' "" (auto ff_satsig_<time>_b<k>.mat in pwd)
%
% Live use: ff.takeControl() + guaranteed release (onCleanup) restores the instrument; nothing is
% transmitted (SA mode only). Antenna on Port 2 (SA input). Recordings use the suite format
% f (1xN Hz), P (sweeps x N dBm), t (datetime), settings (rbw, detector, ...), idn, so
% ff_satsig("identify", file) works on fieldfox_live / ff_spectrum_pro recordings too.
%
% Honest limits:
%  - No I/Q: classes are inferred from spectrum width/shape and timing only. CW vs a terrestrial
%    carrier is separated by the Doppler track, not by the spectrum.
%  - Swept SA + bursty beacons: a packet is seen only if the sweep crosses its band while it is on.
%    Use 'TraceType',"MAXH" windows for sparse beacons (KNACKSAT-2: every 60 s).
%  - GNSS is 17-25 dB below kT: a GNSS-band "signal" on the FieldFox is interference/repeater.
%  - EIRPs are mostly unpublished (ASSUMED values are marked); the level test is deliberately loose.
%  - Demo traces are synthetic (template spectra + noise + real SGP4 geometry): the accuracy shown
%    there is an upper bound. Verify on real passes.
%  - Confidence = posterior over the hypotheses listed (incl. "fixed" and "unknown satellite of class
%    X"); a satellite missing from the catalog/SatNOGS/TLE set cannot be named. "CHECK: poor Doppler
%    fit" marks a winner that still fits badly (weak or fragmented track).
%  - Doppler stacking assumes the signal is on in a reasonable fraction of the stacked sweeps; its dt
%    search is coarse (9 steps over +-2 sigma of the TLE timing error).
%
% Data: CelesTrak (TLE), SatNOGS DB (CC BY-SA 4.0, attribution required), IARU coordination pages,
% team documents (see the Source column). Cached in support\satsig.
%
% See also ff_cubesat_link, ff_signal_id, ff_srs4_analyzer, ff_ccsds_linksim, ff_gnss_margin.

if nargin < 1 || strlength(string(cmd)) == 0
    help('ff_satsig'); return
end
cmd = lower(string(cmd));
switch cmd
    case "catalog"
        [C, info] = cmdCatalog(varargin{:});
        if nargout == 0, printCatalog(C, info); else, out = {C, info}; varargout = out(1:nargout); end
    case "define"
        D = cmdDefine(varargin{:});
        if nargout > 0, varargout{1} = D; elseif isfield(D, 'fig'), dropHidden([D.fig]); end
    case "overhead"
        [T, P, info] = cmdOverhead(varargin{:});
        if nargout == 0, printOverhead(T, P, info); else, out = {T, P, info}; varargout = out(1:nargout); end
        if nargout < 3 && isfield(info, 'fig'), dropHidden(info.fig); end
    case {"identify", "id"}
        if isempty(varargin), error('ff_satsig:usage', 'ff_satsig("identify", data, ...): data = recording file/struct, IP address or FieldFox object.'); end
        [T, info] = cmdIdentify(varargin{:});
        if nargout == 0, printIdentify(T, info); else, out = {T, info}; varargout = out(1:nargout); end
        if nargout < 2, dropHidden(info.figs); end
    case "demo"
        args = varargin;
        if ~any(strcmpi(string(args(1:2:end)), "Demo")), args = [args, {'Demo', "origami"}]; end
        [T, info] = cmdIdentify("demo", args{:});
        if nargout == 0, printIdentify(T, info); else, out = {T, info}; varargout = out(1:nargout); end
        if nargout < 2, dropHidden(info.figs); end
    case "templates"
        tpl = buildTemplates(fullfile(supportDir(), 'satsig_templates.mat'));
        if nargout > 0, varargout{1} = tpl; end
    otherwise
        error('ff_satsig:cmd', 'Unknown subcommand "%s". Use catalog | define | overhead | identify | demo | templates.', cmd);
end
end

%% ======================================================================================
%% options and paths
%% ======================================================================================
function opt = baseOptions()
opt = struct('Site', [13.7563 100.5018 10], 'SiteName', "Bangkok", 'TimeZone', "Asia/Bangkok", ...
    'Time', [], 'Offline', false, 'Plot', true, 'Visible', 'on', 'Prompt', true, 'OutDir', "", ...
    'IncludeInactive', false, 'IncludeSatNOGS', false, 'ElevationMask', 0, 'Quiet', false);
end

function opt = parseOpts(opt, args)
if mod(numel(args), 2) ~= 0, error('ff_satsig:option', 'Options must be Name, Value pairs.'); end
fn = fieldnames(opt);
for k = 1:2:numel(args)
    j = find(strcmpi(fn, string(args{k})), 1);
    if isempty(j)
        error('ff_satsig:option', 'Unknown option "%s". Valid here: %s', string(args{k}), strjoin(fn, ', '));
    end
    opt.(fn{j}) = args{k + 1};
end
opt.Visible = char(string(opt.Visible));
opt.OutDir = string(opt.OutDir);
if strlength(opt.OutDir) > 0 && ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
end

function d = supportDir()
d = fullfile(fileparts(mfilename('fullpath')), 'support', 'satsig');
if ~isfolder(d), mkdir(d); end
end

function say(opt, varargin)
if ~opt.Quiet, fprintf(varargin{:}); end
end

function t = toUTC(t, zone)
if isempty(t.TimeZone), t.TimeZone = char(zone); end
t.TimeZone = 'UTC';
end

%% ======================================================================================
%% CATALOG
%% ======================================================================================
function [C, info] = cmdCatalog(varargin)
opt = baseOptions(); opt.Refresh = false; opt.Search = ""; opt.Export = "";
opt = parseOpts(opt, varargin);
info = struct('notes', strings(0, 1), 'satnogs', []);
C = builtinCatalog(opt);
[S, note] = satnogsTable(opt, logical(opt.Refresh));
info.notes(end+1) = note; info.satnogs = S;
if ~isempty(S)
    C.SatNOGS = strings(height(C), 1);
    for k = 1:height(C)
        if isfinite(C.NORAD(k))
            n = nnz(S.norad == C.NORAD(k) & abs(S.downlink_low_Hz - C.F0_Hz(k)) < 5e3);
            C.SatNOGS(k) = sprintf('%d alive tx within 5 kHz', n);
        end
    end
end
if opt.IncludeSatNOGS && ~isempty(S), C = [C; satnogsRows(S, C)]; end
if strlength(string(opt.Search)) > 0
    q = string(opt.Search);
    C = C(contains(C.Id, q, 'IgnoreCase', true) | contains(C.Satellite, q, 'IgnoreCase', true) | ...
        contains(C.Modulation, q, 'IgnoreCase', true), :);
end
if strlength(string(opt.Export)) > 0
    writetable(C, opt.Export); info.notes(end+1) = "Catalog written to " + string(opt.Export);
end
end

function C = builtinCatalog(opt)
% One row per downlink. Active: prior that it is on the air now (0 = re-entered/decommissioned).
% Template: synthesis key (synthWave) used for the definition sheet and the bandwidth/shape test.
n = ["Id","Satellite","NORAD","Group","TleMatch","GeoLon","F0_Hz","FAlt_Hz","BandLo_Hz","BandHi_Hz", ...
    "OBW_Hz","Modulation","Rate","Framing","Burst","Polarization","EIRP_dBm","PrxMin_dBm","Alt_km", ...
    "Status","Active","Class","Template","Emission","Confidence","Source"];
o2 = "OrigamiSat-2 (FO-126) JS1YRU"; o1 = "OrigamiSat-1 (FO-98) JS1YAX"; k2 = "KNACKSAT-2";
iss = "ISS (ZARYA)"; amU = "435-438 MHz class prior (SatNOGS 2026-10-10: 615 alive UHF transmitters)";
R = {
"origami2_uhf_cw", o2, 68795, "", "", NaN, 437.505e6, NaN, NaN, NaN, 150, "CW Morse 16 wpm (80 cpm), 28-byte hex telemetry", 16, "Morse: callsign + ORIGAMI2 + hex TLM", "periodic beacon", "linear (assumed)", NaN, NaN, 535, "in orbit (2026-088D, 523x544 km, 97.5 deg); SatNOGS good obs to 2026-10-09 incl. INSTED/KMUTNB Bangkok", 1, "CW", "cw16", "", "high", "IARU 932; SatNOGS DB 68795; CelesTrak SATCAT; ANS-140"
"origami2_uhf_afsk", o2, 68795, "", "", NaN, 437.505e6, NaN, NaN, NaN, 10e3, "AFSK 1k2 Bell-202 (1200/2200 Hz, 1700 Hz centre) on FM", 1200, "AX.25 UI (HDLC, NRZI)", "periodic beacon", "linear (assumed)", NaN, NaN, 535, "in orbit; SatNOGS good obs to 2026-10-09", 1, "AFSK-FM", "afsk1k2", "", "high", "IARU 932; SatNOGS DB 68795"
"origami2_c_cw", "OrigamiSat-2 C-band beacon", 68795, "", "", NaN, 5840e6, NaN, NaN, NaN, 200, "CW", NaN, "none", "scheduled passes", "CP, 18.0 dBic deployable reflectarray (deployed 2026-05-23)", 48, NaN, 535, "EIRP 48 dBm ASSUMED (1 W + 18 dBic, beam on station); received by JA1OGZ (IC-905) per ANS-144", 0.5, "CW", "cw16", "", "medium", "IARU 932; SatNOGS DB; AMSAT ANS-144; IEEE TAP 2026 (antenna)"
"origami2_c_bpsk", "OrigamiSat-2 C-band data (BPSK)", 68795, "", "", NaN, 5840e6, NaN, NaN, NaN, 116e3, "BPSK 100 kbps (100k/500k/1M/5M/10M/20M selectable)", 100e3, "unknown (CCSDS-like assumed)", "scheduled passes", "CP 18.0 dBic", 48, NaN, 535, "EIRP 48 dBm ASSUMED (1 W + 18 dBic, beam on station); listed active; no public reception report found", 0.3, "PSK-RRC", "ccsds_bpsk100k", "", "medium", "IARU 932"
"origami2_c_qpsk", "OrigamiSat-2 C-band data (QPSK 20 Mbps)", 68795, "", "", NaN, 5840e6, NaN, NaN, NaN, 11.6e6, "QPSK 20 Mbps (10 Msym/s), RRC 0.35 assumed", 10e6, "unknown", "scheduled passes", "CP 18.0 dBic", 48, NaN, 535, "EIRP 48 dBm ASSUMED (1 W + 18 dBic, beam on station); listed active; no public reception report found", 0.3, "PSK-RRC", "ccsds_qpsk20M", "", "medium", "IARU 932; SSC23 abstract"
"origami1_uhf_cw", o1, 43933, "", "", NaN, 437.505e6, 437.499e6, NaN, NaN, 500, "CW beacon (A1A)", 20, "Morse (content not published)", "continuous beacon", "unknown", NaN, NaN, 500, "RE-ENTERED 2022-04-30", 0, "CW", "cw20", "500HA1A", "high", "SatNOGS DB 43933; CelesTrak SATCAT; IARU 591; amsat-bb 2019-01"
"origami1_uhf_afsk", o1, 43933, "", "", NaN, 437.505e6, 437.499e6, NaN, NaN, 20e3, "AFSK 1k2 on FM (F2D)", 1200, "AX.25 (assumed)", "constant use", "unknown", NaN, NaN, 500, "RE-ENTERED 2022-04-30", 0, "AFSK-FM", "afsk1k2", "20K0F2D", "high", "SatNOGS DB 43933; IARU 591; amsat-bb 2019-01"
"origami1_c_fsk", "OrigamiSat-1 5.84 GHz downlink", 43933, "", "", NaN, 5840e6, NaN, NaN, NaN, 420e3, "FSK 115 kbps (F1D)", 115e3, "unknown", "scheduled (10 s on / 5 s off reported, unverified)", "unknown", NaN, NaN, 500, "RE-ENTERED 2022-04-30", 0, "FSK/GMSK", "fsk115k", "420KF1D", "medium", "SatNOGS DB 43933; amsat-bb 2019-01"
"knacksat2_tlm", k2, 67683, "amateur", "", NaN, 400.630e6, NaN, NaN, NaN, 12.5e3, "FSK 9k6 (G3RUH scrambler)", 9600, "AX.25", "beacon every 60 s", "linear, rotating (RHCP recommended at the GS)", 30, NaN, 353, "in orbit (ISS deploy 2026-02-03); SatNOGS good obs to 2026-10-09; not an amateur frequency", 1, "FSK/GMSK", "g3ruh9k6", "", "high", "KNACKSAT-2 reception resource (team PDF); SatNOGS DB 67683"
"knacksat2_aprs", "KNACKSAT-2 APRS digipeater HS0K", 67683, "amateur", "", NaN, 145.825e6, NaN, NaN, NaN, 12.5e3, "9k6 GMSK AX.25 APRS digipeater", 9600, "AX.25 APRS (no path)", "on schedule only", "linear", 27, NaN, 353, "IARU coordinated; activation windows only (verify)", 0.3, "FSK/GMSK", "gmsk9k6", "", "medium", "IARU 961; team PDF"
"napa1", "NAPA-1 (RTAF-SAT 1)", 46320, "", "", NaN, 401.125e6, NaN, NaN, NaN, 1.5e3, "BPSK 1200 Bd", 1200, "unknown", "periodic", "unknown", NaN, NaN, 291, "in orbit, 289x293 km and decaying fast; SatNOGS good obs 2026-10-08", 1, "PSK-RRC", "bpsk1k2", "", "medium", "SatNOGS DB 46320; CelesTrak SATCAT"
"iss_rpt", iss + " NA1SS voice repeater", 25544, "stations", "", NaN, 437.800e6, NaN, NaN, NaN, 16e3, "NBFM voice (uplink 145.990 MHz, PL 67 Hz)", NaN, "none", "when in use", "linear", 40, NaN, 420, "active per ARISS status 2026-10-07 (EIRP ASSUMED ~10 W)", 0.6, "FM-voice", "nbfm", "16K0F3E", "high", "ARISS current status page"
"iss_aprs", iss + " RS0ISS APRS", 25544, "stations", "", NaN, 437.825e6, NaN, NaN, NaN, 12e3, "AFSK 1k2 on FM (APRS)", 1200, "AX.25 APRS", "bursty packets", "linear", NaN, NaN, 420, "APRS on 437.825 MHz (NOT 145.825) per ARISS 2026-10-07", 0.6, "AFSK-FM", "afsk1k2", "", "high", "ARISS current status page"
"iss_sstv", iss + " SSTV / crew voice", 25544, "stations", "", NaN, 145.800e6, NaN, NaN, NaN, 12e3, "FM SSTV (PD120/Robot36) or voice", NaN, "none", "event-based", "linear", NaN, NaN, 420, "event-based only", 0.1, "FM-voice", "sstv", "", "medium", "SatNOGS DB ISS; ARISS"
"iss_hamtv", iss + " HamTV", 25544, "stations", "", NaN, 2395e6, NaN, NaN, NaN, 2.7e6, "DVB-S QPSK 2 Msym/s FEC 1/2 RRC 0.35", 2e6, "MPEG-TS", "test transmissions", "RHCP", 40, NaN, 420, "test signal per ARISS 2026-10-07", 0.5, "PSK-RRC", "dvbs2_2M", "", "medium", "ARISS; BATC wiki; ARRL"
"meteor_m23_lrpt", "METEOR-M2 3 LRPT", 57166, "weather", "", NaN, 137.900e6, 137.100e6, NaN, NaN, 96e3, "OQPSK 72 ksym/s RRC 0.6 (80k mode exists)", 72e3, "CCSDS CADU: conv k7 r1/2 + RS(255,223) I=4, ASM 1ACFFC1D", "continuous (over land)", "RHCP", 37, NaN, 815, "operational; SatNOGS good obs 2026-10-09", 1, "OQPSK-RRC", "lrpt_oqpsk", "", "high (freq) / medium (alpha 0.6)", "SatNOGS DB + network; CelesTrak SATCAT"
"meteor_m24_lrpt", "METEOR-M2 4 LRPT", 59051, "weather", "", NaN, 137.900e6, 137.100e6, NaN, NaN, 96e3, "OQPSK 72 ksym/s RRC 0.6 (80k mode exists)", 72e3, "CCSDS CADU: conv k7 r1/2 + RS(255,223) I=4, ASM 1ACFFC1D", "continuous (over land)", "RHCP", 37, NaN, 820, "operational; SatNOGS good obs 2026-10-10", 1, "OQPSK-RRC", "lrpt_oqpsk", "", "high (freq) / medium (alpha 0.6)", "SatNOGS DB + network; CelesTrak SATCAT"
"noaa15_apt", "NOAA 15 APT", 25338, "", "", NaN, 137.620e6, NaN, NaN, NaN, 34e3, "APT: AM 2.4 kHz subcarrier on FM (dev 17 kHz)", 4160, "APT, 2 lines/s", "OFF", "RHCP", NaN, NaN, 807, "DECOMMISSIONED 2025-08-19 20:37 UTC", 0, "APT-FM", "apt", "", "high", "NOAA OSPO POES status"
"noaa18_apt", "NOAA 18 APT", 28654, "", "", NaN, 137.9125e6, NaN, NaN, NaN, 34e3, "APT: AM 2.4 kHz subcarrier on FM (dev 17 kHz)", 4160, "APT, 2 lines/s", "OFF", "RHCP", NaN, NaN, 854, "DECOMMISSIONED 2025-06-06 17:40 UTC", 0, "APT-FM", "apt", "", "high", "NOAA OSPO POES status"
"noaa19_apt", "NOAA 19 APT", 33591, "", "", NaN, 137.100e6, NaN, NaN, NaN, 34e3, "APT: AM 2.4 kHz subcarrier on FM (dev 17 kHz)", 4160, "APT, 2 lines/s", "OFF", "RHCP", NaN, NaN, 870, "DECOMMISSIONED 2025-08-13 18:55 UTC", 0, "APT-FM", "apt", "", "high", "NOAA OSPO POES status"
"orbcomm", "ORBCOMM OG2 subscriber downlinks", NaN, "orbcomm", "^ORBCOMM", NaN, 137.500e6, NaN, 137.200e6, 137.800e6, 6e3, "SDPSK (+/-90 deg) 4800 bps RRC ~0.4", 4800, "ORBCOMM", "continuous", "RHCP", 35, NaN, 750, "operational constellation (channel per satellite not fixed here)", 1, "PSK-RRC", "orbcomm", "", "medium", "ESA 'The Orbcomm Experience'; NTIA; SatNOGS"
"iridium_ra", "Iridium ring alert (ch 7)", NaN, "iridium-NEXT", "^IRIDIUM", NaN, 1626.270833e6, NaN, NaN, NaN, 31.5e3, "DE-QPSK 25 ksym/s RRC 0.4", 25e3, "Iridium simplex", "8.28 ms bursts in 90 ms frames", "RHCP", 48, NaN, 780, "operational (EIRP 48 dBm = ROUGH guess)", 1, "PSK-RRC", "iridium", "", "high (plan) / low (EIRP)", "sigidwiki Iridium"
"iridium_duplex", "Iridium duplex/simplex channels", NaN, "iridium-NEXT", "^IRIDIUM", NaN, 1621.25e6, NaN, 1616.0e6, 1626.5e6, 31.5e3, "DE-QPSK 25 ksym/s RRC 0.4 (TDMA)", 25e3, "Iridium TDMA", "8.28 ms bursts in 90 ms frames", "RHCP", 48, NaN, 780, "operational (EIRP 48 dBm = ROUGH guess)", 1, "PSK-RRC", "iridium", "", "high", "1616 + 0.020833(2n-1) MHz, n = 1..240"
"gps_l1", "GPS L1 C/A", NaN, "gnss", "^GPS", NaN, 1575.42e6, NaN, NaN, NaN, 2.046e6, "BPSK(1) CDMA 1.023 Mcps", 1.023e6, "LNAV 50 bps", "continuous", "RHCP", NaN, -128.5, 20200, "operational; ~18 dB BELOW kT: never visible on a swept SA", 1, "GNSS", "gps_l1", "", "high", "IS-GPS-200"
"gps_l5", "GPS L5", NaN, "gnss", "^GPS", NaN, 1176.45e6, NaN, NaN, NaN, 20.46e6, "QPSK(10) 10.23 Mcps", 10.23e6, "CNAV", "continuous", "RHCP", NaN, -124.9, 20200, "operational (block IIF/III); below the noise floor", 1, "GNSS", "gps_l5", "", "high", "IS-GPS-705"
"gal_e1", "Galileo E1 OS", NaN, "gnss", "^GSAT0", NaN, 1575.42e6, NaN, NaN, NaN, 4.092e6, "CBOC(6,1,1/11)", 1.023e6, "I/NAV", "continuous", "RHCP", NaN, -127, 23222, "operational; below the noise floor", 1, "GNSS", "gal_e1", "", "high", "Galileo OS SIS ICD"
"gal_e5", "Galileo E5 AltBOC", NaN, "gnss", "^GSAT0", NaN, 1191.795e6, NaN, NaN, NaN, 51.15e6, "AltBOC(15,10)", 10.23e6, "F/NAV + I/NAV", "continuous", "RHCP", NaN, -125, 23222, "operational; below the noise floor", 1, "GNSS", "gal_e5", "", "high", "Galileo OS SIS ICD"
"bds_b1i", "BeiDou B1I", NaN, "gnss", "^BEIDOU", NaN, 1561.098e6, NaN, NaN, NaN, 4.092e6, "BPSK(2) 2.046 Mcps", 2.046e6, "D1/D2", "continuous", "RHCP", NaN, -133, 21528, "operational (GEO/IGSO/MEO over Asia); below the noise floor", 1, "GNSS", "bds_b1i", "", "high", "BDS ICD B1I"
"bds_b3i", "BeiDou B3I", NaN, "gnss", "^BEIDOU", NaN, 1268.52e6, NaN, NaN, NaN, 20.46e6, "BPSK(10) 10.23 Mcps", 10.23e6, "D1/D2", "continuous", "RHCP", NaN, -133, 21528, "operational; below the noise floor", 1, "GNSS", "bds_b3i", "", "high", "BDS ICD B3I"
"glo_l1", "GLONASS L1OF (FDMA k=-7..+6)", NaN, "gnss", "^COSMOS", NaN, 1602e6, NaN, 1598.0625e6, 1605.375e6, 1.022e6, "BPSK 0.511 Mcps, 1602 + 0.5625k MHz", 0.511e6, "50 bps + 100 Hz meander", "continuous", "RHCP", NaN, -131, 19100, "operational; below the noise floor", 1, "GNSS", "glo_l1", "", "high", "GLONASS ICD"
"glo_l2", "GLONASS L2OF (FDMA)", NaN, "gnss", "^COSMOS", NaN, 1246e6, NaN, 1242.9375e6, 1248.625e6, 1.022e6, "BPSK 0.511 Mcps, 1246 + 0.4375k MHz", 0.511e6, "", "continuous", "RHCP", NaN, -131, 19100, "operational; below the noise floor", 1, "GNSS", "glo_l1", "", "high", "GLONASS ICD"
"qzss_l6", "QZSS L6 (CLAS/MADOCA)", NaN, "gnss", "^QZS", NaN, 1278.75e6, NaN, NaN, NaN, 5.115e6, "CASM 2.5575 Mcps (BPSK-like main lobe 5.1 MHz)", 2.5575e6, "L6D/L6E 2 kbps", "continuous", "RHCP", NaN, -125, 35786, "operational; visible from Bangkok; below the noise floor", 1, "GNSS", "qzss_l6", "", "high", "IS-QZSS-L6"
"gk2a_lrit", "GEO-KOMPSAT-2A LRIT (128.2E)", 43823, "weather", "", NaN, 1692.14e6, NaN, NaN, NaN, 170e3, "BPSK 128 ksym/s RRC 0.5", 128e3, "CCSDS conv + RS, ASM 1ACFFC1D, 64 kbps", "continuous", "linear?", NaN, NaN, 35786, "operational; el 54 deg from Bangkok; needs a dish/grid + LNA", 1, "PSK-RRC", "gk2a_lrit", "", "medium", "sigidwiki GK-2A LRIT"
"inmarsat_l", "Inmarsat-4 F1 L-band (143.5E)", 28628, "", "", NaN, 1542e6, NaN, 1525e6, 1559e6, 2.4e3, "BPSK 1200 sym/s (STD-C / EGC) and others", 1200, "Inmarsat-C", "continuous", "RHCP", NaN, NaN, 35786, "SATCAT status P (partial); channel plan not verified", 0.8, "PSK-RRC", "bpsk1k2", "", "low", "CelesTrak SATCAT"
"thuraya_gmr1", "Thuraya GMR-1 (98.5E)", NaN, "", "", 98.5, 1542e6, NaN, 1525e6, 1559e6, 31.25e3, "pi/4-DQPSK 23.4 ksym/s", 23.4e3, "GMR-1 TDMA", "TDMA", "LHCP", NaN, NaN, 35786, "operational (Osmocom GMR-1)", 0.8, "PSK-RRC", "thuraya", "", "medium", "Osmocom GMR wiki"
"amsat_uhf_cw", "Unknown amateur CubeSat (CW beacon)", NaN, "", "", NaN, 436.5e6, NaN, 435e6, 438e6, 150, "CW 12-25 wpm", 20, "Morse", "periodic beacon", "linear", 30, NaN, 500, amU + ": CW 80", 0.5, "CW", "cw16", "", "medium", "SatNOGS DB statistics"
"amsat_uhf_afsk", "Unknown amateur CubeSat (AFSK 1k2)", NaN, "", "", NaN, 436.5e6, NaN, 435e6, 438e6, 10e3, "AFSK 1k2 on FM", 1200, "AX.25", "periodic beacon", "linear", 30, NaN, 500, amU + ": AFSK 44", 0.5, "AFSK-FM", "afsk1k2", "", "medium", "SatNOGS DB statistics"
"amsat_uhf_fsk9k6", "Unknown amateur CubeSat (FSK/GFSK 9k6)", NaN, "", "", NaN, 436.5e6, NaN, 435e6, 438e6, 12e3, "FSK / GFSK 9k6 (G3RUH)", 9600, "AX.25 / AX.100", "periodic beacon", "linear", 30, NaN, 500, amU + ": FSK 97 + GFSK 70", 0.5, "FSK/GMSK", "g3ruh9k6", "", "medium", "SatNOGS DB statistics"
"amsat_uhf_gmsk4k8", "Unknown amateur CubeSat (GMSK 4k8)", NaN, "", "", NaN, 436.5e6, NaN, 435e6, 438e6, 5e3, "GMSK 4k8", 4800, "AX.25 / AX.100", "periodic beacon", "linear", 30, NaN, 500, amU + ": GMSK 63", 0.5, "FSK/GMSK", "gmsk4k8", "", "medium", "SatNOGS DB statistics"
"amsat_uhf_bpsk1k2", "Unknown amateur CubeSat (BPSK 1k2)", NaN, "", "", NaN, 436.5e6, NaN, 435e6, 438e6, 1.5e3, "BPSK 1k2", 1200, "AX.25 / FUNcube", "periodic beacon", "linear", 30, NaN, 500, amU + ": BPSK 41", 0.5, "PSK-RRC", "bpsk1k2", "", "medium", "SatNOGS DB statistics"
"amsat_uhf_lora", "Unknown amateur CubeSat (LoRa 125 kHz)", NaN, "", "", NaN, 436.5e6, NaN, 435e6, 438e6, 125e3, "LoRa CSS SF7-12, BW 125 kHz", NaN, "LoRa", "periodic beacon", "linear", 30, NaN, 500, amU + ": LoRa 39", 0.5, "LoRa", "lora125", "", "medium", "SatNOGS DB statistics"
"amsat_uhf_fm", "Unknown amateur satellite (FM voice)", NaN, "", "", NaN, 436.5e6, NaN, 435e6, 438e6, 16e3, "NBFM voice transponder", NaN, "none", "when in use", "linear", 30, NaN, 500, amU + ": FM 61", 0.5, "FM-voice", "nbfm", "", "medium", "SatNOGS DB statistics"
"amsat_vhf_afsk", "Unknown amateur satellite (VHF AFSK/FSK)", NaN, "", "", NaN, 145.9e6, NaN, 145.8e6, 146.0e6, 10e3, "AFSK 1k2 / FSK 9k6", 1200, "AX.25", "periodic beacon", "linear", 30, NaN, 500, "145.8-146 MHz amateur-satellite class prior", 0.5, "AFSK-FM", "afsk1k2", "", "medium", "SatNOGS DB statistics"
"sband_ccsds", "S-band CCSDS downlink (Satlab SRS-4 class)", NaN, "", "", NaN, 2245e6, NaN, 2200e6, 2290e6, 1.35e6, "QPSK/OQPSK/BPSK RRC 0.35, 0.1-2 Msym/s (here 1 Msym/s)", 1e6, "CCSDS TM, ASM 1ACFFC1D, CC/RS", "scheduled passes", "RHCP typ.", 33, NaN, 500, "class (see ff_srs4_analyzer)", 0.5, "PSK-RRC", "srs4_qpsk1M", "", "medium", "Satlab SRS-4 datasheet; CCSDS 131.0-B"
"dvbs2_lband_if", "DVB-S2 carrier at an LNB L-band IF (conducted)", NaN, "", "", NaN, 1550e6, NaN, 950e6, 2150e6, 36e6, "DVB-S2 QPSK/8PSK 30 Msym/s RRC 0.2 (e.g. Thaicom Ku via LNB)", 30e6, "DVB-S2 BBFRAME / MPEG-TS", "continuous", "n/a (cable)", NaN, NaN, 35786, "conducted measurement of an LNB output, not over the air", 0.2, "PSK-RRC", "dvbs2_30M", "", "medium", "ETSI EN 302 307"
};
C = cell2table(R, 'VariableNames', cellstr(n));
for v = ["Id","Satellite","Group","TleMatch","Modulation","Framing","Burst","Polarization","Status","Class","Template","Emission","Confidence","Source"]
    C.(v) = string(C.(v));
end
hg = opt.Site(3);
C.MaxDoppler_Hz = arrayfun(@(f, a) maxDoppler(f, a, hg), C.F0_Hz, C.Alt_km);
C.Kind = repmat("sat", height(C), 1);
C.Kind(C.TleMatch ~= "") = "constellation";
C.Kind(isfinite(C.GeoLon)) = "geo";
C.Kind(~isfinite(C.NORAD) & C.TleMatch == "" & ~isfinite(C.GeoLon)) = "generic";
C.SrcPrior = ones(height(C), 1);
C.SrcPrior(C.Kind == "generic") = 0.3;
end

function d = maxDoppler(f0, altkm, hg)
if ~isfinite(altkm) || altkm > 30000, d = 50; return, end
d = max(abs(dopplerShiftCircularOrbit(0:0.5:90, altkm*1e3, hg, f0)));
end

function r = maxDopplerRate(f0, altkm)
% Doppler rate at the closest approach of an overhead pass (circular orbit, no Earth rotation).
if ~isfinite(altkm) || altkm > 30000, r = 0; return, end
Re = 6371e3; h = altkm*1e3; w = sqrt(3.986004418e14/(Re + h)^3);
t = linspace(-30, 30, 3001); th = w*t;
Rg = sqrt(Re^2 + (Re + h)^2 - 2*Re*(Re + h)*cos(th));
r = max(abs(gradient(gradient(Rg, t), t)))*f0/physconst('LightSpeed');
end

function printCatalog(C, info)
fprintf('\n%-18s %6s %12s %9s  %-38s %s\n', 'Id', 'NORAD', 'f0 (MHz)', 'OBW', 'Modulation', 'Status');
for k = 1:height(C)
    fprintf('%-18s %6s %12.6f %9s  %-38s %s\n', C.Id(k), numstr(C.NORAD(k)), C.F0_Hz(k)/1e6, fmtHz(C.OBW_Hz(k)), ...
        extractBefore(C.Modulation(k) + string(blanks(38)), 39), extractBefore(C.Status(k) + string(blanks(60)), 61));
end
for n = info.notes(:)', fprintf('Note: %s\n', n); end
end

%% ----- SatNOGS DB (bulk, cached 24 h; offline: shipped subset) ------------------------
function [S, note] = satnogsTable(opt, refresh)
file = fullfile(supportDir(), 'satnogs_subset.csv');
note = "";
if refresh && ~opt.Offline
    d = fullfile(supportDir(), 'satnogs'); if ~isfolder(d), mkdir(d); end
    fT = fullfile(d, 'transmitters.json'); fS = fullfile(d, 'satellites.json');
    [okT, nT] = fetchFile("https://db.satnogs.org/api/transmitters/?format=json", fT, 24);
    [okS, nS] = fetchFile("https://db.satnogs.org/api/satellites/?format=json", fS, 24);
    note = "SatNOGS: " + nT + "; " + nS;
    if okT && okS
        try
            S = buildSatnogsSubset(fT, fS);
            writetable(S, file);
            note = note + sprintf("; subset rebuilt (%d alive transmitters on in-orbit satellites)", height(S));
        catch err
            note = note + "; subset rebuild failed: " + string(err.message);
        end
    end
end
if ~isfile(file)
    S = table(); note = strtrim(note + " No SatNOGS subset file."); return
end
io = detectImportOptions(file, 'Delimiter', ',', 'TextType', 'string');
txt = intersect(["uuid", "sat_name", "description", "mode", "service", "type", "iaru", "updated"], string(io.VariableNames));
num = intersect(["norad", "downlink_low_Hz", "downlink_high_Hz", "baud"], string(io.VariableNames));
io = setvartype(io, cellstr(txt), 'string'); io = setvartype(io, cellstr(num), 'double');
S = readtable(file, io);
need = ["uuid", "norad", "sat_name", "description", "downlink_low_Hz", "downlink_high_Hz", "mode", "baud", "service", "type", "iaru", "updated"];
if ~all(ismember(need, string(S.Properties.VariableNames)))
    S = table(); note = strtrim(note + " SatNOGS subset file has unexpected columns (rebuild with 'Refresh', true)."); return
end
for v = txt, S.(v)(ismissing(S.(v))) = ""; end
d = dir(file);
if strlength(note) == 0
    note = sprintf("SatNOGS DB subset: %d transmitters (file of %s; CC BY-SA 4.0, db.satnogs.org)", height(S), ...
        string(datetime(d.datenum, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd')));
end
end

function S = buildSatnogsSubset(fT, fS)
T = jsondecode(fileread(fT)); X = jsondecode(fileread(fS));
if ~iscell(T), T = num2cell(T); end
if ~iscell(X), X = num2cell(X); end
sid = strings(numel(X), 1); sst = sid; snm = sid; snor = nan(numel(X), 1);
for k = 1:numel(X)
    sid(k) = string(gf(X{k}, 'sat_id', "")); sst(k) = string(gf(X{k}, 'status', ""));
    snm(k) = replace(string(gf(X{k}, 'name', "")), ",", " "); snor(k) = double(gf(X{k}, 'norad_cat_id', NaN));
end
rows = cell(0, 12);
for k = 1:numel(T)
    t = T{k};
    j = find(sid == string(gf(t, 'sat_id', "")), 1);
    if isempty(j) || sst(j) ~= "in orbit", continue, end
    nor = double(gf(t, 'norad_cat_id', NaN)); if ~isfinite(nor), nor = snor(j); end
    f = double(gf(t, 'downlink_low', NaN));
    if ~isfinite(nor) || ~isfinite(f) || f < 100e6 || f > 6.5e9, continue, end
    if ~isequal(gf(t, 'alive', false), true) || string(gf(t, 'status', "")) ~= "active", continue, end
    ty = string(gf(t, 'type', ""));
    if ~ismember(ty, ["Transmitter", "Transceiver", "Transponder"]), continue, end
    rows(end+1, :) = {string(gf(t, 'uuid', "")), nor, snm(j), ...
        extractBefore(replace(string(gf(t, 'description', "")), [",", newline], [" ", " "]) + string(blanks(60)), 61), ...
        f, double(gf(t, 'downlink_high', NaN)), string(gf(t, 'mode', "")), double(gf(t, 'baud', NaN)), ...
        string(gf(t, 'service', "")), ty, string(isequal(gf(t, 'iaru_coordination', ""), 'IARU Coordinated')), ...
        extractBefore(string(gf(t, 'updated', "")) + "          ", 11)}; %#ok<AGROW>
end
S = cell2table(rows, 'VariableNames', {'uuid','norad','sat_name','description','downlink_low_Hz','downlink_high_Hz', ...
    'mode','baud','service','type','iaru','updated'});
S.description = strtrim(S.description);
S = sortrows(S, 'downlink_low_Hz');
end

function v = gf(s, name, d)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = d; end
if isstring(d) && ~(ischar(v) || isstring(v)), v = d; end
end

function R = satnogsRows(S, C)
% SatNOGS transmitters as catalog rows (lower prior; duplicates of catalog rows dropped).
n = height(S); keep = true(n, 1);
for k = 1:n
    keep(k) = ~any(C.NORAD == S.norad(k) & abs(C.F0_Hz - S.downlink_low_Hz(k)) < 3e3);
end
S = S(keep, :); n = height(S);
R = C([], :);
if n == 0, return, end
cls = strings(n, 1); key = cls; bw = nan(n, 1);
for k = 1:n, [cls(k), key(k), bw(k)] = modeMap(S.mode(k), S.baud(k)); end
mod = strtrim(S.mode + " " + string(arrayfun(@(b) numstr(b), S.baud)) + " Bd");
R = table("satnogs_" + extractBefore(S.uuid + "        ", 9), S.sat_name, S.norad, repmat("satnogs", n, 1), ...
    strings(n, 1), nan(n, 1), S.downlink_low_Hz, nan(n, 1), nan(n, 1), nan(n, 1), bw, mod, S.baud, ...
    repmat("", n, 1), repmat("unknown", n, 1), repmat("unknown", n, 1), nan(n, 1), nan(n, 1), repmat(600, n, 1), ...
    "SatNOGS alive flag (" + S.updated + "): " + S.description, repmat(0.6, n, 1), cls, key, strings(n, 1), ...
    repmat("low", n, 1), repmat("SatNOGS DB (CC BY-SA 4.0)", n, 1), ...
    'VariableNames', C.Properties.VariableNames(1:26));
R.MaxDoppler_Hz = maxDoppler(1, 600, 10)*R.F0_Hz;          % scales with f0
R.Kind = repmat("sat", n, 1); R.SrcPrior = repmat(0.5, n, 1);
if ismember("SatNOGS", C.Properties.VariableNames), R.SatNOGS = repmat("SatNOGS row", n, 1); end
end

function [cls, key, bw] = modeMap(mode, baud)
m = upper(string(mode)); b = baud; if ~isfinite(b), b = 0; end
cls = "unknown"; key = ""; bw = max(2*b, 10e3);
if m == "CW", cls = "CW"; key = "cw20"; bw = 150;
elseif startsWith(m, "AFSK"), cls = "AFSK-FM"; key = "afsk1k2"; bw = 10e3;
elseif ismember(m, ["FM", "FMN"]), cls = "FM-voice"; key = "nbfm"; bw = 16e3;
elseif m == "SSTV", cls = "FM-voice"; key = "sstv"; bw = 12e3;
elseif m == "LORA", cls = "LoRa"; key = "lora125"; bw = 125e3;
elseif m == "APT", cls = "APT-FM"; key = "apt"; bw = 34e3;
elseif m == "LRPT", cls = "OQPSK-RRC"; key = "lrpt_oqpsk"; bw = 96e3;
elseif contains(m, ["FSK", "GMSK", "MSK", "AX.100", "FFSK"])
    cls = "FSK/GMSK"; bw = max(1.3*b + 2e3, 3e3);
    if b <= 1800, key = "gmsk1k2"; elseif b <= 6000, key = "gmsk4k8"; elseif b <= 14400, key = "g3ruh9k6";
    elseif b <= 28800, key = "gmsk19k2"; else, key = "gmsk38k4"; end
elseif contains(m, "BPSK") || m == "PSK"
    cls = "PSK-RRC"; bw = max(1.35*b, 1e3);
    if b <= 3000, key = "bpsk1k2"; else, key = "bpsk9k6"; end
end
end

%% ----- downloads with caching (CelesTrak 2 h, SatNOGS 24 h) ------------------------------
function [ok, note] = fetchFile(url, file, maxAgeH)
ok = false; err = string(file) + ".err";
if isfile(file) && fileAgeH(file) < maxAgeH
    ok = true; note = sprintf("%s cached %.1f h ago", fname(file), fileAgeH(file)); return
end
if isfile(err) && fileAgeH(err) < min(6, maxAgeH)
    note = sprintf("%s: last download failed %.1f h ago, not retried", fname(file), fileAgeH(err));
    ok = isfile(file); return
end
try
    fprintf('Downloading %s\n', url);
    tmp = string(file) + ".part";
    websave(tmp, url, weboptions('Timeout', 120, 'UserAgent', userAgent()));
    jsondecode(fileread(tmp));                       % must parse
    movefile(tmp, file, 'f');
    if isfile(err), delete(err); end
    ok = true; note = sprintf("%s downloaded now", fname(file));
catch e
    writelines(string(datetime('now')) + " " + string(e.message), err);
    note = sprintf("%s download failed (%s)", fname(file), e.message);
    ok = isfile(file);
end
end

function [txt, src] = fetchText(url, cache, maxAgeH, offline)
txt = ""; err = string(cache) + ".err";
if isfile(cache) && fileAgeH(cache) < maxAgeH
    txt = string(fileread(cache)); src = sprintf("cache (%.0f min old)", 60*fileAgeH(cache)); return
end
if offline
    if isfile(cache), txt = string(fileread(cache)); src = sprintf("cache (%.1f h old, offline)", fileAgeH(cache));
    else, src = "none (offline)"; end
    return
end
if isfile(err) && fileAgeH(err) < maxAgeH
    src = "CelesTrak failed less than 2 h ago: not retried";
    if isfile(cache), txt = string(fileread(cache)); src = src + sprintf("; cache %.1f h old", fileAgeH(cache)); end
    return
end
try
    fprintf('Fetching %s\n', url);
    r = string(webread(url, weboptions('ContentType', 'text', 'Timeout', 30, 'UserAgent', userAgent())));
    if isempty(parseTLE(r)), error('no TLE in reply: "%s"', extractBefore(r + string(blanks(80)), 81)); end
    writelines(r, cache); txt = r; src = "CelesTrak (downloaded now)";
    if isfile(err), delete(err); end
catch e
    writelines(string(datetime('now')) + " " + string(e.message), err);
    src = "CelesTrak failed (" + string(e.message) + ")";
    if isfile(cache), txt = string(fileread(cache)); src = src + sprintf("; cache %.1f h old", fileAgeH(cache)); end
end
end

function u = userAgent()
u = "ff_satsig/1.0 (MATLAB R2026a; student CubeSat research, Bangkok)";
end

function h = fileAgeH(file)
h = Inf; d = dir(file);
if ~isempty(d), h = hours(datetime('now') - datetime(d.datenum, 'ConvertFrom', 'datenum')); end
end

function s = fname(file)
[~, a, b] = fileparts(file); s = string(a) + string(b);
end

%% ----- TLEs ------------------------------------------------------------------------------
function E = parseTLE(txt, group)
if nargin < 2, group = ""; end
L = strtrim(splitlines(string(txt))); L = L(strlength(L) > 0);
E = struct('name', {}, 'norad', {}, 'l1', {}, 'l2', {}, 'epoch', {}, 'group', {});
for k = 1:numel(L) - 1
    a = char(L(k)); b = char(L(k + 1));
    if numel(a) >= 69 && numel(b) >= 69 && a(1) == '1' && b(1) == '2' && strcmp(a(3:7), b(3:7))
        nm = "SAT " + string(strtrim(a(3:7)));
        if k > 1 && ~startsWith(L(k - 1), ["1 ", "2 "]), nm = regexprep(L(k - 1), '^0 ', ''); end
        E(end+1) = struct('name', strtrim(nm), 'norad', str2double(a(3:7)), 'l1', string(a(1:69)), ...
            'l2', string(b(1:69)), 'epoch', tleEpoch(a), 'group', string(group)); %#ok<AGROW>
    end
end
end

function t = tleEpoch(l1)
yy = str2double(l1(19:20)); doy = str2double(l1(21:32));
yr = 2000 + yy; if yy >= 57, yr = 1900 + yy; end
t = datetime(yr, 1, 1, 'TimeZone', 'UTC') + days(doy - 1);
end

function [E, notes] = getTLEs(opt, groups, norads)
% Group files (cached 2 h), CATNR for NORADs not in a group, shipped sample as the last resort.
notes = strings(0, 1);
tdir = fullfile(supportDir(), 'tle'); if ~isfolder(tdir), mkdir(tdir); end
sdir = fullfile(supportDir(), 'tle_sample');
E = parseTLE("");
groups = unique(string(groups(strlength(string(groups)) > 0)));
for g = groups(:)'
    url = sprintf('https://celestrak.org/NORAD/elements/gp.php?GROUP=%s&FORMAT=TLE', g);
    [txt, src] = fetchText(url, fullfile(tdir, g + ".tle"), 2, logical(opt.Offline));
    Eg = parseTLE(txt, g);
    if isempty(Eg) && isfile(fullfile(sdir, g + ".tle"))
        Eg = parseTLE(fileread(fullfile(sdir, g + ".tle")), g); src = src + " -> shipped sample (" + g + ")";
    end
    notes(end+1) = sprintf("TLE group %s: %d entries, %s", g, numel(Eg), src); %#ok<AGROW>
    E = [E, Eg]; %#ok<AGROW>
end
norads = unique(norads(isfinite(norads)));
have = [E.norad]; todo = setdiff(norads(:)', have); missing = [];
maxQ = 15;                                       % politeness: at most 15 CATNR queries per call
if numel(todo) > maxQ
    notes(end+1) = sprintf("%d NORAD IDs are not in the TLE groups: only the first %d are queried by CATNR.", numel(todo), maxQ);
    todo = todo(1:maxQ);
end
for n = todo
    url = sprintf('https://celestrak.org/NORAD/elements/gp.php?CATNR=%d&FORMAT=TLE', n);
    [txt, src] = fetchText(url, fullfile(tdir, sprintf('catnr_%d.tle', n)), 2, logical(opt.Offline));
    En = parseTLE(txt, "catnr"); En = En([En.norad] == n);
    if isempty(En)
        for f = dir(fullfile(sdir, '*.tle'))'
            Es = parseTLE(fileread(fullfile(f.folder, f.name)), erase(string(f.name), ".tle"));
            En = Es([Es.norad] == n);
            if ~isempty(En), src = src + " -> shipped sample"; En = En(1); break, end
        end
    end
    if isempty(En), missing(end+1) = n; %#ok<AGROW>
    else, notes(end+1) = sprintf("TLE NORAD %d: %s", n, src); E = [E, En]; end %#ok<AGROW>
end
if ~isempty(missing)
    notes(end+1) = "No TLE for NORAD " + strjoin(string(missing), ", ") + ternary(logical(opt.Offline), " (offline, not in the cache or shipped sample)", "");
end
% de-duplicate (newest epoch wins; keep the group tag of the first occurrence)
if ~isempty(E)
    [~, ord] = sort([E.epoch], 'descend'); E = E(ord);
    [~, ia] = unique([E.norad], 'stable'); E = E(ia);
end
end

function file = writeTLEfile(E)
file = string(tempname) + ".tle";
L = strings(3*numel(E), 1);
for k = 1:numel(E), L(3*k-2:3*k) = [E(k).name; E(k).l1; E(k).l2]; end
writelines(L, file);
end

%% ----- orbit prediction: satelliteScenario + SGP4 + states (ECEF) -------------------------
function P = predictSats(E, tA, tB, opt, stepS)
% Geometry of each TLE entry seen from opt.Site between tA and tB (UTC). Doppler factor
% d = -v_r/c (received f = f0 (1 + d)); identical to satelliteScenario dopplershift.
if nargin < 5, stepS = max(1, ceil(seconds(tB - tA)/2500)); end
P = struct('norad', {}, 'name', {}, 'group', {}, 't', {}, 'tau', {}, 'el', {}, 'az', {}, 'rng', {}, ...
    'dfac', {}, 'rate', {}, 'epoch', {});
if isempty(E), return, end
sc = satelliteScenario(tA, tB, stepS);
file = writeTLEfile(E); c = onCleanup(@() delete(file));
sat = satellite(sc, char(file), 'OrbitPropagator', 'sgp4');
[p, v, tt] = states(sat, 'CoordinateFrame', 'ecef');
lat = opt.Site(1); lon = opt.Site(2); h = opt.Site(3);
g = lla2ecef([lat lon h])'; c0 = physconst('LightSpeed'); wgs = wgs84Ellipsoid;
tau = seconds(tt - tA); tau = tau(:);
for k = 1:numel(sat)
    pk = p(:, :, k); vk = v(:, :, k);
    rho = pk - g; rng = vecnorm(rho); rr = sum(rho.*vk, 1)./rng;
    [az, el] = ecef2aer(pk(1, :), pk(2, :), pk(3, :), lat, lon, h, wgs);
    d = -rr(:)/c0;
    P(k) = struct('norad', E(k).norad, 'name', E(k).name, 'group', E(k).group, 't', tt(:), 'tau', tau, ...
        'el', el(:), 'az', az(:), 'rng', rng(:), 'dfac', d, 'rate', gradient(d, tau), 'epoch', E(k).epoch);
end
end

function P = predictGeo(lonDeg, tA, tB, opt, name)
tau = [0; seconds(tB - tA)]; tt = tA + seconds(tau);
[az, el, r] = geodetic2aer(0, lonDeg, 35786e3, opt.Site(1), opt.Site(2), opt.Site(3), wgs84Ellipsoid);
P = struct('norad', NaN, 'name', name, 'group', "geo", 't', tt, 'tau', tau, 'el', [el; el], 'az', [az; az], ...
    'rng', [r; r], 'dfac', [0; 0], 'rate', [0; 0], 'epoch', tA);
end

%% ======================================================================================
%% OVERHEAD
%% ======================================================================================
function [T, Pass, info] = cmdOverhead(varargin)
opt = baseOptions(); opt.Constellations = true; opt.Hours = 24; opt.Search = "";
opt = parseOpts(opt, varargin);
t = opt.Time; if isempty(t), t = datetime('now', 'TimeZone', 'UTC'); end
if ~isdatetime(t), t = datetime(t); end
t = toUTC(t, opt.TimeZone);
C = builtinCatalog(opt); info.notes = strings(0, 1);
if opt.IncludeSatNOGS
    [S, note] = satnogsTable(opt, false); info.notes(end+1) = note;
    if ~isempty(S), C = [C; satnogsRows(S, C)]; end
end
if ~opt.IncludeInactive, C = C(C.Active > 0, :); end
if strlength(string(opt.Search)) > 0
    q = string(opt.Search);
    C = C(contains(C.Id, q, 'IgnoreCase', true) | contains(C.Satellite, q, 'IgnoreCase', true), :);
end
if ~opt.Constellations, C = C(C.Kind ~= "constellation", :); end
grp = unique(C.Group(C.Kind == "constellation" | (C.Kind == "sat" & C.Group ~= "satnogs")));
if opt.IncludeSatNOGS, grp = unique([grp; "satnogs"]); end
[E, n1] = getTLEs(opt, grp, C.NORAD(C.Kind == "sat" & C.Group ~= "satnogs"));
info.notes = [info.notes(:); n1(:)];
% satellites needed: specific NORADs + constellation members
keepE = ismember([E.norad], C.NORAD(C.Kind == "sat"));
for k = find(C.Kind == "constellation")'
    keepE = keepE | ([E.group] == C.Group(k) & ~cellfun(@isempty, regexp(cellstr([E.name]), C.TleMatch(k), 'once')));
end
E = E(keepE);
say(opt, 'Propagating %d satellites (SGP4) for %s UTC ...\n', numel(E), string(t, 'yyyy-MM-dd HH:mm:ss'));
P = predictSats(E, t - seconds(30), t + seconds(30), opt, 10);
rows = {};
for k = 1:numel(P)
    [~, i0] = min(abs(P(k).tau - 30));
    el = P(k).el(i0);
    if el < opt.ElevationMask, continue, end
    ri = find((C.Kind == "sat" & C.NORAD == P(k).norad) | ...
        (C.Kind == "constellation" & C.Group == P(k).group & ~cellfun(@isempty, regexp(cellstr(repmat(P(k).name, height(C), 1)), C.TleMatch, 'once'))));
    for r = ri(:)'
        f0 = C.F0_Hz(r); fd = f0*P(k).dfac(i0);
        rows(end+1, :) = {P(k).name, P(k).norad, C.Id(r), f0/1e6, P(k).az(i0), el, P(k).rng(i0)/1e3, ...
            -P(k).dfac(i0)*physconst('LightSpeed')/1e3, fd, f0*P(k).rate(i0), (f0 + fd)/1e6, C.Status(r), ...
            days(t - P(k).epoch)}; %#ok<AGROW>
    end
end
for r = find(C.Kind == "geo")'
    Pg = predictGeo(C.GeoLon(r), t, t + seconds(1), opt, C.Satellite(r));
    if Pg.el(1) >= opt.ElevationMask
        rows(end+1, :) = {C.Satellite(r), NaN, C.Id(r), C.F0_Hz(r)/1e6, Pg.az(1), Pg.el(1), Pg.rng(1)/1e3, 0, 0, 0, ...
            C.F0_Hz(r)/1e6, C.Status(r), NaN}; %#ok<AGROW>
    end
end
vn = {'Satellite','NORAD','Downlink','F0_MHz','Az_deg','El_deg','Range_km','RangeRate_km_s','Doppler_Hz', ...
    'DopplerRate_Hz_s','Frx_MHz','Status','TLEage_days'};
if isempty(rows), T = cell2table(cell(0, numel(vn)), 'VariableNames', vn);
else, T = sortrows(cell2table(rows, 'VariableNames', vn), {'F0_MHz', 'El_deg'}, {'ascend', 'descend'}); end
% next passes of the specific (non-constellation) satellites
Es = E(ismember([E.norad], C.NORAD(C.Kind == "sat" & C.Group ~= "satnogs" & C.Alt_km < 30000)));
Pass = passTable(Es, t, hours(opt.Hours), opt);
info.time = t; info.site = opt.Site; info.siteName = opt.SiteName; info.catalog = C; info.timeZone = opt.TimeZone;
if opt.Plot && height(T) > 0
    fig = figure('Visible', opt.Visible, 'Name', 'ff_satsig overhead', 'Position', [60 60 760 700]);
    u = unique(T.Satellite, 'stable'); [~, iu] = unique(T.Satellite, 'stable');
    lab = u; if numel(u) > 40, lab = repmat("", size(u)); end
    skyplot(fig, T.Az_deg(iu), T.El_deg(iu), lab);
    title(sprintf('Catalog satellites above %g deg, %s, %s UTC', opt.ElevationMask, opt.SiteName, string(t, 'yyyy-MM-dd HH:mm')));
    info.fig = fig;
    if strlength(opt.OutDir) > 0, exportgraphics(fig, fullfile(opt.OutDir, "overhead_skyplot.png"), 'Resolution', 120); end
    closeSaved(fig, opt);
end
if strlength(opt.OutDir) > 0
    writetable(T, fullfile(opt.OutDir, "overhead.csv")); writetable(Pass, fullfile(opt.OutDir, "passes.csv"));
end
end

function Pass = passTable(E, t, dur, opt)
vn = {'Satellite','NORAD','AOS_UTC','LOS_UTC','MaxEl_deg','TMax_UTC','Duration_min'};
Pass = cell2table(cell(0, numel(vn)), 'VariableNames', vn);
if isempty(E), return, end
sc = satelliteScenario(t, t + dur, 20);
file = writeTLEfile(E); c = onCleanup(@() delete(file));
sat = satellite(sc, char(file), 'OrbitPropagator', 'sgp4');
gs = groundStation(sc, opt.Site(1), opt.Site(2), 'Altitude', opt.Site(3), 'MinElevationAngle', opt.ElevationMask, 'Name', opt.SiteName);
[p, ~, tt] = states(sat, 'CoordinateFrame', 'ecef');
rows = {};
for k = 1:numel(sat)
    I = accessIntervals(access(sat(k), gs));
    [~, el] = ecef2aer(p(1, :, k), p(2, :, k), p(3, :, k), opt.Site(1), opt.Site(2), opt.Site(3), wgs84Ellipsoid);
    for j = 1:height(I)
        in = tt >= I.StartTime(j) & tt <= I.EndTime(j);
        if ~any(in), continue, end
        [m, im] = max(el(in)); ti = tt(in);
        rows(end+1, :) = {E(k).name, E(k).norad, I.StartTime(j), I.EndTime(j), m, ti(im), minutes(I.EndTime(j) - I.StartTime(j))}; %#ok<AGROW>
    end
end
if ~isempty(rows), Pass = sortrows(cell2table(rows, 'VariableNames', vn), 'AOS_UTC'); end
end

function printOverhead(T, P, info)
fprintf('\nAbove the horizon at %s UTC (%s local), %s [%.4f %.4f]:\n', string(info.time, 'yyyy-MM-dd HH:mm:ss'), ...
    string(datetime(info.time, 'TimeZone', info.timeZone), 'HH:mm'), info.siteName, info.site(1), info.site(2));
isG = ismember(T.Downlink, info.catalog.Id(info.catalog.Class == "GNSS"));
if height(T) == 0, fprintf('  (none of the catalog satellites)\n'); end
fprintf('%-26s %-17s %11s %6s %5s %7s %9s %8s %13s\n', 'Satellite', 'Downlink', 'f0 MHz', 'Az', 'El', 'Range', 'Doppler', 'Hz/s', 'f_rx MHz');
for k = find(~isG)'
    fprintf('%-26s %-17s %11.6f %6.1f %5.1f %7.0f %+9.0f %+8.1f %13.6f\n', extractBefore(T.Satellite(k) + string(blanks(26)), 27), ...
        T.Downlink(k), T.F0_MHz(k), T.Az_deg(k), T.El_deg(k), T.Range_km(k), T.Doppler_Hz(k), T.DopplerRate_Hz_s(k), T.Frx_MHz(k));
end
if any(isG)
    fam = ["gps", "GPS"; "gal", "Galileo"; "glo", "GLONASS"; "bds", "BeiDou"; "qzss", "QZSS"];
    parts = strings(0, 1);
    for j = 1:size(fam, 1)
        n = numel(unique(T.Satellite(isG & startsWith(T.Downlink, fam(j, 1)))));
        if n > 0, parts(end+1) = sprintf('%s %d', fam(j, 2), n); end %#ok<AGROW>
    end
    fprintf(['GNSS satellites above the mask: %s (rows in T; their signals are 17-25 dB below kT, so a ' ...
        'visible line in a GNSS band is interference)\n'], strjoin(parts, ', '));
end
if height(P) > 0
    fprintf('\nNext passes (UTC; local = UTC+7 for Bangkok):\n');
    for k = 1:min(height(P), 25)
        ip = ""; if P.AOS_UTC(k) <= info.time + seconds(1), ip = "  (in progress)"; end
        fprintf('  %-26s %s - %s  max el %4.1f deg%s\n', extractBefore(string(P.Satellite(k)) + string(blanks(26)), 27), ...
            string(P.AOS_UTC(k), 'MM-dd HH:mm:ss'), string(P.LOS_UTC(k), 'HH:mm:ss'), P.MaxEl_deg(k), ip);
    end
end
for n = info.notes(:)', fprintf('Note: %s\n', n); end
end

%% ======================================================================================
%% DEFINE
%% ======================================================================================
function D = cmdDefine(name, varargin)
if nargin < 1 || strlength(string(name)) == 0, name = ""; end
opt = baseOptions(); opt.RxGain = 0; opt.Seed = 1;
opt = parseOpts(opt, varargin);
C = builtinCatalog(opt);
q = string(name);
if q == "", printCatalog(C, struct('notes', "Give a name or Id, e.g. ff_satsig(""define"", ""OrigamiSat-1"")")); D = []; return, end
hit = strcmpi(C.Id, q);
if ~any(hit)
    qn = lower(regexprep(q, '[\s_\-]', ''));
    hay = lower(regexprep(C.Id + " " + C.Satellite + " " + C.Modulation, '[\s_\-]', ''));
    hit = contains(hay, qn);
end
if ~any(hit), error('ff_satsig:define', 'No catalog row matches "%s". See ff_satsig("catalog").', q); end
rows = find(hit)';
D = struct([]);
for r = rows
    row = C(r, :);
    W = synthWave(row.Template, opt.Seed, synthText(row));
    d = defineSheet(row, W, opt);
    if opt.Plot, d.fig = plotDefinition(row, W, d, opt); end
    if isempty(D), D = d; else, D(end+1) = d; end %#ok<AGROW>
end
end

function txt = synthText(row)
% Callsign and illustrative payload for the CW / AX.25 synthesis of this row.
c = regexp(row.Satellite, '[A-Z]{1,2}\d[A-Z]{1,4}', 'match', 'once');
if strlength(c) == 0, c = "CUBESAT"; end
if contains(row.Id, "origami2"), name = "ORIGAMI2"; elseif contains(row.Id, "origami1"), name = "ORIGAMI1";
else, name = upper(regexprep(extractBefore(row.Satellite + " ", " "), '[^A-Za-z0-9]', '')); end
txt = struct('call', string(c), 'name', string(name));
end

function d = defineSheet(row, W, opt)
x = W.x(:); fs = W.fs;
on = find(abs(x) > 1e-3*max(abs(x)));                     % trim the silent padding, keep keying gaps
xa = x(on(1):on(end));
[o, flo, fhi] = obw(complex(x), fs);
try
    bw3 = powerbw(complex(x), fs);
catch
    bw3 = NaN;
end
a = abs(xa); papr = 10*log10(max(a.^2)/mean(a.^2)); cv = std(a)/mean(a);
f0 = row.F0_Hz; alt = row.Alt_km; maxD = row.MaxDoppler_Hz; rateD = maxDopplerRate(f0, alt);
obwNom = row.OBW_Hz;
span = niceUp(max(1.3*(obwNom + 2*maxD), 3*obwNom));
rbw = rbwStep(min(span/400, max(obwNom/3, 10)));
det = "AVER"; if ismember(row.Class, ["CW", "AFSK-FM", "FSK/GMSK"]) || contains(lower(row.Burst), ["beacon", "burst", "tdma", "packet"]), det = "POS"; end
ff = struct('Center_Hz', f0, 'Span_Hz', span, 'RBW_Hz', rbw, 'Points', 401, 'Detector', det, 'Preamp', true, ...
    'TraceType', "CLRW (record every sweep for the Doppler track); MAXH 10 s windows for sparse beacons", ...
    'ZeroSpan', sprintf('centre on the Doppler-corrected frequency, RBW %s, SAMP/POS, sweep 2-10 s', fmtHz(rbwStep(min(5e6, max(obwNom, 10)), true))));
% link estimate (FieldFox DANL -153 dBm/Hz preamp; with a 20 dB / 1 dB NF LNA -170 dBm/Hz)
el = [90; 45; 20; 10; 5]; eirp = row.EIRP_dBm; basis = "catalog value; see Status for ASSUMED/ROUGH";
if ~isfinite(eirp) && alt < 3000, eirp = 30; basis = "ASSUMED 30 dBm (1 W, 0 dBi)"; end
if isfinite(row.PrxMin_dBm)
    Prx = repmat(row.PrxMin_dBm, size(el)) + opt.RxGain; Rg = nan(size(el)); basis = "ICD minimum received power";
elseif isfinite(eirp)
    Rg = slantRangeCircularOrbit(el, alt*1e3, opt.Site(3)); Rg = Rg(:);
    Prx = eirp - fspl(Rg, physconst('LightSpeed')/f0) - 3 + opt.RxGain;   % 3 dB polarisation/pointing
else
    Rg = nan(size(el)); Prx = nan(size(el));
end
inR = Prx - max(0, 10*log10(max(obwNom, 1)/rbw));
SNRff = inR - (-153 + 10*log10(1.0645*rbw)); SNRlna = inR - (-170 + 10*log10(1.0645*rbw));
link = table(el, Rg/1e3, Prx, SNRff, SNRlna, 'VariableNames', {'El_deg', 'Range_km', 'Prx_dBm', 'SNR_FieldFox_dB', 'SNR_withLNA_dB'});
st = row.Status; if row.Active <= 0, st = st + "  ->  NOT receivable now (link table below is hypothetical)"; end
sheet = struct('Id', row.Id, 'Satellite', row.Satellite, 'NORAD', row.NORAD, 'Status', st, ...
    'Frequency', fmtHz(f0) + ternary(isfinite(row.FAlt_Hz), " (alt " + fmtHz(row.FAlt_Hz) + ")", ""), ...
    'Band', ternary(isfinite(row.BandLo_Hz), fmtHz(row.BandLo_Hz) + " - " + fmtHz(row.BandHi_Hz), ""), ...
    'Emission', decodeEmission(row.Emission), 'Modulation', row.Modulation, 'Rate', row.Rate, 'Framing', row.Framing, ...
    'Burst', row.Burst, 'Polarization', row.Polarization, 'EIRP', ternary(isfinite(eirp), sprintf('%.0f dBm (%s)', eirp, basis), "unknown"), ...
    'Orbit', sprintf('%g km', alt), 'MaxDoppler_Hz', maxD, 'MaxDopplerRate_Hz_s', rateD, ...
    'SynthOBW99_Hz', o, 'SynthBW3dB_Hz', bw3, 'SynthPAPR_dB', papr, 'SynthEnvelopeCV', cv, 'Synthesis', W.note, 'Message', W.msg, ...
    'Hint', classHint(row.Class), 'Confidence', row.Confidence, 'Source', row.Source);
d = struct('id', row.Id, 'sheet', sheet, 'wave', W, 'obw', [o flo fhi], 'fieldfox', ff, 'link', link, 'fig', []);
if ~opt.Quiet
    fprintf('\n==== %s  [%s] ====\n', row.Satellite, row.Id);
    fn = fieldnames(sheet);
    for k = 1:numel(fn)
        v = sheet.(fn{k}); if isnumeric(v), v = numstr(v); end
        if strlength(string(v)) > 0, fprintf('  %-20s %s\n', fn{k}, string(v)); end
    end
    fprintf('  FieldFox setup       centre %s, span %s, RBW %s, %d pts, %s detector, preamp on\n', fmtHz(f0), fmtHz(span), fmtHz(rbw), 401, det);
    fprintf('                       zero span: %s\n', ff.ZeroSpan);
    fprintf('  Link (0 dBi + RxGain %g dB, -3 dB pol.):\n', opt.RxGain); disp(link);
end
end

function fig = plotDefinition(row, W, d, opt)
x = W.x(:); fs = W.fs;
fig = figure('Visible', opt.Visible, 'Name', "ff_satsig define " + row.Id, 'Position', [40 40 1250 820]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% 1) spectrum
ax = nexttile(tl);
L = max(64, min(2^nextpow2(fs/max(d.obw(1)/150, fs/2^16)), 2^floor(log2(numel(x)/4))));
[p, f] = pwelch(x, hann(L), L/2, L, fs, 'centered');
[sc, unit] = fscale(max(abs(d.obw(2:3))) * 3 + 1);
pd = 10*log10(p/max(p));
plot(ax, f*sc, pd, 'Color', [0.6 0.6 0.6]); hold(ax, 'on');
rbw = d.fieldfox.RBW_Hz;
if rbw > 2*(f(2) - f(1))
    g = exp(-4*log(2)*((-ceil(3*rbw/(f(2)-f(1))):ceil(3*rbw/(f(2)-f(1))))'*(f(2)-f(1))/rbw).^2);
    ps = conv(p, g, 'same'); plot(ax, f*sc, 10*log10(ps/max(ps)), 'b', 'LineWidth', 1.4);
    legend(ax, 'Welch PSD', sprintf('FieldFox view (RBW %s)', fmtHz(rbw)), 'Location', 'south');
end
xline(ax, d.obw(2:3)*sc, 'r--', 'HandleVisibility', 'off');
lim = min(max(4*max(abs(d.obw(2:3))), 5*(f(2) - f(1))), fs/2);
xlim(ax, [-lim lim]*sc); ylim(ax, [-60 3]); grid(ax, 'on');
xlabel(ax, "Offset from carrier (" + unit + ")"); ylabel(ax, 'dB (normalised)');
title(ax, sprintf('Spectrum: OBW99 %s (obw on complex baseband)', fmtHz(d.obw(1))));
% 2) spectrogram (centred, works for real-valued keying too)
ax = nexttile(tl);
fr = max(d.obw(1)/40, fs/4096); if W.kind == "cw", fr = 2*W.Rs; end
L = 2^nextpow2(max(16, fs/fr)); L = min(L, 2^floor(log2(numel(x)/8)));
[S, F, Tt] = spectrogram(x, hann(L), round(0.75*L), max(L, 256), fs, 'centered');
S = abs(S).^2; [tsc, tun] = tscale(numel(x)/fs);
imagesc(ax, Tt*tsc, F*sc, 10*log10(S/max(S(:)) + 1e-12)); axis(ax, 'xy'); clim(ax, [-50 0]); colormap(ax, turbo);
ylim(ax, [-lim lim]*sc);
xlabel(ax, "Time (" + tun + ")"); ylabel(ax, "Frequency (" + unit + ")"); title(ax, sprintf('Spectrogram (%d-point Hann, %s resolution)', L, fmtHz(fs/L)));
% 3) zero-span envelope as the FieldFox would show it (Gaussian RBW >= OBW, sample detector)
ax = nexttile(tl);
rz = min(5e6, max(d.obw(1), 10)); rz = min(rz, fs/2.5);
if W.kind == "cw", rz = rbwStep(max(3*W.Rs, 10), true); end    % resolves the dots
X = fft(x); ff_ = (-floor(numel(x)/2):ceil(numel(x)/2)-1)'*fs/numel(x);
Hf = ifftshift(exp(-2*log(2)*(ff_/rz).^2)); y = ifft(X.*Hf);
env = 20*log10(abs(y)/max(abs(y)) + 1e-6);
tt = (0:numel(x)-1)'/fs;
plot(ax, tt*tsc, env, 'k'); ylim(ax, [-40 3]); grid(ax, 'on');
xlabel(ax, "Time (" + tun + ")"); ylabel(ax, 'dB rel. peak');
title(ax, sprintf('Zero-span envelope (RBW %s): PAPR %.1f dB, CV %.2f', fmtHz(rz), d.sheet.SynthPAPR_dB, d.sheet.SynthEnvelopeCV));
% 4) constellation / discriminator / code autocorrelation / keying
ax = nexttile(tl);
plotKind(ax, W);
sgtitle(tl, sprintf('%s  [%s]  %s  -  %s', row.Satellite, fmtHz(row.F0_Hz), row.Modulation, row.Status), 'Interpreter', 'none', 'FontSize', 10);
if strlength(opt.OutDir) > 0, exportgraphics(fig, fullfile(opt.OutDir, "define_" + row.Id + ".png"), 'Resolution', 110); end
closeSaved(fig, opt);
end

function plotKind(ax, W)
x = W.x(:); fs = W.fs;
switch W.kind
    case {"psk", "burst-psk"}
        y = mfSymbols(W);
        plot(ax, real(y), imag(y), '.', 'MarkerSize', 6); axis(ax, 'equal'); grid(ax, 'on');
        m = max(abs(y)) * 1.2; if ~isfinite(m) || m == 0, m = 1; end
        xlim(ax, [-m m]); ylim(ax, [-m m]);
        title(ax, sprintf('Constellation after RRC matched filter (%s)', ternary(W.oqpsk, "OQPSK: Q delayed Ts/2", "symbol-spaced")));
        xlabel(ax, 'I'); ylabel(ax, 'Q');
    case "fm"
        n = min(numel(x) - 1, round(fs*min(0.02, 40/W.Rs)));
        a = find(abs(x) > 0.1*max(abs(x)), 1); if isempty(a), a = 1; end
        seg = x(a:min(end, a + n));
        fi = diff(unwrap(angle(seg)))*fs/(2*pi);
        plot(ax, (0:numel(fi)-1)/fs*1e3, fi/1e3, 'b'); grid(ax, 'on');
        xlabel(ax, 'Time (ms)'); ylabel(ax, 'Instantaneous frequency (kHz)');
        title(ax, 'FM discriminator output (constellation undefined for constant-envelope FM)');
    case "chirp"
        n = min(numel(x) - 1, round(fs*3*2^10/125e3));
        fi = diff(unwrap(angle(x(1:n+1))))*fs/(2*pi);
        plot(ax, (0:n-1)/fs*1e3, fi/1e3); grid(ax, 'on'); xlabel(ax, 'Time (ms)'); ylabel(ax, 'kHz');
        title(ax, 'LoRa-like chirps: instantaneous frequency');
    case "gnss"
        n = round(fs*1e-3); x1 = x(1:n);
        r = abs(ifft(fft(x1).*conj(fft(x1)))); r = r/max(r);
        lag = (0:n-1)'/fs*W.Rs; lag(lag > lag(end)/2) = lag(lag > lag(end)/2) - n/fs*W.Rs;
        [lag, o] = sort(lag); plot(ax, lag, r(o)); grid(ax, 'on'); xlim(ax, [-20 20]);
        xlabel(ax, 'Lag (chips)'); ylabel(ax, '|R| normalised');
        title(ax, 'Code autocorrelation over 1 ms (acquisition principle; below noise on an SA)');
    otherwise   % cw
        n = min(numel(x), round(fs*6));
        plot(ax, (0:n-1)/fs, abs(x(1:n)), 'k'); grid(ax, 'on'); ylim(ax, [-0.1 1.2]);
        xlabel(ax, 'Time (s)'); ylabel(ax, 'Key (envelope)');
        title(ax, "Morse keying: " + W.msg, 'Interpreter', 'none');
end
end

function y = mfSymbols(W)
x = W.x(:); sps = W.sps; beta = W.beta;
h = rcosdesign(beta, 10, sps, 'sqrt');
z = conv(x, h); z = z(10*sps + 1:end - 10*sps);           % tx + rx filter delay
if numel(z) < 20*sps, y = complex(0); return, end
e = zeros(sps, 1);
for k = 1:sps, e(k) = mean(abs(z(k:sps:end)).^2); end
[~, k0] = max(e);
if W.oqpsk
    hs = round(sps/2); i = k0:sps:numel(z) - hs;
    y = real(z(i)) + 1j*imag(z(i + hs));
    if mean(abs(real(z(i)))) < mean(abs(imag(z(i + hs))))/2
        i = k0 + hs:sps:numel(z) - hs; y = real(z(i)) + 1j*imag(z(i + hs));
    end
else
    y = z(k0:sps:end);
end
y = y(abs(y) > 0.2*max(abs(y)));
y = y(1:min(end, 3000));
end

function h = classHint(cls)
switch cls
    case "CW", h = "Narrower than the RBW; on/off keying in zero span (dot = 1.2/wpm s). Separate from a terrestrial carrier by its Doppler track.";
    case "AFSK-FM", h = "~10 kHz flat-topped FM; constant envelope in zero span; packets ~0.3-1 s.";
    case "FSK/GMSK", h = "~1.0-1.3 x bit rate wide, rounded top; constant envelope; bursts.";
    case "PSK-RRC", h = "Flat top with RRC skirts, OBW = (1+alpha) Rs; envelope CV 0.2-0.4 (QPSK).";
    case "OQPSK-RRC", h = "Same spectrum as QPSK but envelope CV ~0.16 (OQPSK avoids zero crossings).";
    case "GNSS", h = "Below the noise floor on a swept SA. A visible line in a GNSS band is interference.";
    case "LoRa", h = "Flat 125 kHz block; chirps visible in a spectrogram/zero span.";
    case "APT-FM", h = "~34 kHz FM with discrete lines (2.4 kHz subcarrier).";
    case "FM-voice", h = "8-16 kHz FM, varying with audio; constant envelope.";
    otherwise, h = "";
end
end

function s = decodeEmission(e)
e = upper(strtrim(string(e))); s = "";
if strlength(e) < 7, return, end
c = char(e); bwp = c(1:4); L = regexp(bwp, '[HKMG]', 'match', 'once');
if isempty(L), return, end
mult = struct('H', 1, 'K', 1e3, 'M', 1e6, 'G', 1e9);
bw = str2double(strrep(bwp, L, '.'))*mult.(L);
m1 = containers.Map({'N','A','H','R','J','B','C','F','G','D','P','K','L','M','Q','V','W','X'}, ...
    {'unmodulated carrier','AM double sideband','SSB full carrier','SSB reduced carrier','SSB suppressed carrier', ...
    'independent sidebands','vestigial sideband','frequency modulation','phase modulation','AM + angle modulation', ...
    'pulse unmodulated','pulse AM','pulse width','pulse position','pulse angle','pulse combination','combination','other'});
m2 = containers.Map({'0','1','2','3','7','8','9','X'}, {'no modulating signal','digital, no subcarrier', ...
    'digital with modulating subcarrier','analogue, single channel','two or more digital channels', ...
    'two or more analogue channels','composite digital + analogue','other'});
m3 = containers.Map({'N','A','B','C','D','E','F','W','X'}, {'no information','telegraphy for aural reception (Morse)', ...
    'automatic telegraphy','facsimile','data / telemetry / telecommand','telephony (voice)','television', ...
    'combination','other'});
g = @(m, k) ternary(isKey(m, k), string(m(k)), "?");
s = sprintf('%s = necessary bandwidth %s; %s; %s; %s', e, fmtHz(bw), g(m1, c(5)), g(m2, c(6)), g(m3, c(7)));
end

%% ======================================================================================
%% SYNTHESIS (toolbox generators) and PSD TEMPLATES
%% ======================================================================================
function keys = templateKeys()
keys = ["cw16","cw20","afsk1k2","g3ruh9k6","gmsk1k2","gmsk4k8","gmsk9k6","gmsk19k2","gmsk38k4","bpsk1k2", ...
    "bpsk9k6","lora125","apt","nbfm","sstv","lrpt_oqpsk","ccsds_bpsk100k","ccsds_qpsk20M","srs4_qpsk1M", ...
    "gk2a_lrit","fsk115k","dvbs2_2M","dvbs2_30M","iridium","orbcomm","thuraya","gps_l1","gps_l5","gal_e1", ...
    "gal_e5","bds_b1i","bds_b3i","glo_l1","qzss_l6"];
end

function W = synthWave(key, seed, txt)
% Synthetic complex-baseband waveform for a template key (toolbox generators where they exist).
% txt.call / txt.name personalise the CW and AX.25 content (illustrative, not the real telemetry).
if nargin < 2, seed = 1; end
if nargin < 3, txt = struct('call', "JS1YRU", 'name', "ORIGAMI2"); end
key = string(key); s0 = rng; rng(seed); restore = onCleanup(@() rng(s0));
W = struct('key', key, 'x', [], 'fs', NaN, 'Rs', NaN, 'kind', "psk", 'note', "", 'sps', NaN, 'beta', NaN, ...
    'oqpsk', false, 'msg', "");
switch key
    case {"cw16", "cw20"}
        wpm = str2double(extractAfter(key, 2)); fs = 2000;
        hx = reshape(dec2hex(randi([0 255], 1, 6), 2)', 1, []);
        msg = txt.call + " " + txt.name + " " + string(hx);
        [env, tdot] = morseEnv(msg, wpm, fs);
        W.x = complex([zeros(round(0.3*fs), 1); env; zeros(round(0.3*fs), 1)]);
        W.fs = fs; W.Rs = 1/tdot; W.kind = "cw";
        W.msg = msg + sprintf(" (%g wpm, dot %.0f ms; illustrative content)", wpm, 1e3*tdot);
        W.note = "custom Morse keyer (PARIS timing, 5 ms raised-cosine key shaping)";
    case "afsk1k2"
        fs = 48e3; spb = fs/1200;
        b = ax25Bits("CQ", txt.call, txt.name + " BAT=7.94V ADCS=2 T=+21C SEQ=0042");
        nz = nrziEnc(b);
        tone = 1200 + 1000*(1 - repelem(nz, spb));
        aud = cos(2*pi*cumsum(tone)/fs);
        g = zeros(round(0.05*fs), 1);
        W.x = [g; fmBB(aud, fs, 3000); g]; W.fs = fs; W.Rs = 1200; W.kind = "fm";
        W.msg = sprintf('AX.25 UI frame %s>CQ, %d bits incl. flags/stuffing (illustrative payload)', txt.call, numel(b));
        W.note = "AX.25 UI + FCS (comm.CRCGenerator, CRC-16/X.25 check 0x906E), HDLC stuffing, NRZI, Bell-202 AFSK, FM dev 3 kHz (complex baseband)";
    case "g3ruh9k6"
        sps = 10; R = 9600; fs = R*sps;
        b = [repmat([0 1 1 1 1 1 1 0]', 40, 1); ax25Bits("CQ", "HS0K", "KNACKSAT-2 TLM " + string(reshape(dec2hex(randi([0 255], 1, 90), 2)', 1, [])))];
        scr = comm.Scrambler('CalculationBase', 2, 'Polynomial', '1 + z^-12 + z^-17', 'InitialConditions', zeros(1, 17));
        bs = scr(nrziEnc(b));
        h = rcosdesign(0.5, 6, sps, 'normal'); h = h/max(h);
        m = conv(upsample(2*double(bs) - 1, sps), h, 'same'); m = m/max(abs(m));
        g = zeros(round(0.02*fs), 1);
        W.x = [g; fmBB(m, fs, 3000); g]; W.fs = fs; W.Rs = R; W.kind = "fm";
        W.msg = sprintf('%d bits (40 flags TX delay + AX.25 frame)', numel(b));
        W.note = "AX.25 + NRZI + G3RUH scrambler (comm.Scrambler 1+z^-12+z^-17) + raised-cosine shaping + FM dev 3 kHz";
    case {"gmsk1k2", "gmsk4k8", "gmsk9k6", "gmsk19k2", "gmsk38k4"}
        R = str2double(replace(extractAfter(key, 4), "k", ".")) * 1e3; sps = 10; fs = R*sps;
        scr = comm.Scrambler('CalculationBase', 2, 'Polynomial', '1 + z^-12 + z^-17', 'InitialConditions', zeros(1, 17));
        b = scr(nrziEnc([repmat([0 1 1 1 1 1 1 0]', 20, 1); ax25Bits("CQ", "CUBESAT", string(blanks(120)) + "TLM")]));
        gm = comm.GMSKModulator('BitInput', true, 'BandwidthTimeProduct', 0.5, 'SamplesPerSymbol', sps);
        W.x = gm(double(b)); W.fs = fs; W.Rs = R; W.kind = "fm";
        W.note = "comm.GMSKModulator BT 0.5 on scrambled NRZI AX.25 bits";
    case {"bpsk1k2", "bpsk9k6"}
        R = 1200; sps = 40; if key == "bpsk9k6", R = 9600; sps = 10; end
        sym = pskmod(randi([0 1], 800, 1), 2);
        tx = comm.RaisedCosineTransmitFilter('RolloffFactor', 0.35, 'FilterSpanInSymbols', 10, 'OutputSamplesPerSymbol', sps);
        W.x = tx([sym; zeros(10, 1)]); W.fs = R*sps; W.Rs = R; W.sps = sps; W.beta = 0.35;
        W.note = "pskmod + comm.RaisedCosineTransmitFilter (RRC 0.35)";
    case "lora125"
        W.x = loraCSS(10, 125e3, 4, [zeros(1, 8) randi([0 1023], 1, 24)]); W.fs = 500e3; W.Rs = 125e3/2^10; W.kind = "chirp";
        W.note = "custom LoRa-like CSS chirps SF10 BW 125 kHz (no LoRa PHY in MATLAB)";
    case "apt"
        fsW = 20800; img = [];
        for k = 1:2
            syncA = repmat([0 0 1 1], 1, 14);
            img = [img syncA(1:39) 0.5 + 0.4*sin(2*pi*(1:2041)/300 + k)]; %#ok<AGROW>
        end
        aud = ammod(0.13 + 0.87*repelem(img, 5)', 2400, fsW);
        aud = resample(aud, 10, 1); fs = 208e3;
        W.x = fmBB(aud/max(abs(aud)), fs, 17e3); W.fs = fs; W.Rs = 4160; W.kind = "fm";
        W.note = "ammod 2.4 kHz AM subcarrier (APT lines) + FM dev 17 kHz";
    case "nbfm"
        fs = 48e3; v = bandpass(randn(fs, 1), [300 3000], fs); v = v/max(abs(v));
        W.x = fmBB(v, fs, 5e3); W.fs = fs; W.Rs = 3000; W.kind = "fm";
        W.note = "band-limited noise 'voice' (bandpass 300-3000 Hz) + FM dev 5 kHz";
    case "sstv"
        fs = 48e3; n = 0.5e-3*fs; tones = [];
        for line = 1:4
            tones = [tones; repmat(1200, round(5e-3*fs), 1); repelem(1500 + 800*rand(300, 1), n)]; %#ok<AGROW>
        end
        W.x = fmBB(cos(2*pi*cumsum(tones)/fs), fs, 3.5e3); W.fs = fs; W.Rs = 2000; W.kind = "fm";
        W.note = "SSTV-like audio (1200 Hz sync, 1500-2300 Hz pixels) + FM dev 3.5 kHz";
    case "lrpt_oqpsk"
        tm = ccsdsTMWaveformGenerator("ChannelCoding", "concatenated", "RSInterleavingDepth", 4, ...
            "Modulation", "OQPSK", "RolloffFactor", 0.6, "SamplesPerSymbol", 8);
        W.x = runTM(tm, 2); W.fs = 72e3*8; W.Rs = 72e3; W.sps = 8; W.beta = 0.6; W.oqpsk = true;
        W.note = "ccsdsTMWaveformGenerator: concatenated (conv k7 r1/2 + RS(255,223) I=4), OQPSK, RRC 0.6";
    case "ccsds_bpsk100k"
        tm = ccsdsTMWaveformGenerator("ChannelCoding", "none", "Modulation", "BPSK", "RolloffFactor", 0.35, "SamplesPerSymbol", 8);
        W.x = runTM(tm, 4); W.fs = 800e3; W.Rs = 100e3; W.sps = 8; W.beta = 0.35;
        W.note = "ccsdsTMWaveformGenerator BPSK uncoded RRC 0.35 (OrigamiSat-2 C-band 100 kbps mode; framing assumed)";
    case "ccsds_qpsk20M"
        tm = ccsdsTMWaveformGenerator("ChannelCoding", "none", "Modulation", "QPSK", "RolloffFactor", 0.35, "SamplesPerSymbol", 4);
        W.x = runTM(tm, 4); W.fs = 40e6; W.Rs = 10e6; W.sps = 4; W.beta = 0.35;
        W.note = "ccsdsTMWaveformGenerator QPSK uncoded RRC 0.35, 10 Msym/s (20 Mbps)";
    case "srs4_qpsk1M"
        tm = ccsdsTMWaveformGenerator("ChannelCoding", "convolutional", "Modulation", "QPSK", "RolloffFactor", 0.35, "SamplesPerSymbol", 4);
        W.x = runTM(tm, 3); W.fs = 4e6; W.Rs = 1e6; W.sps = 4; W.beta = 0.35;
        W.note = "ccsdsTMWaveformGenerator QPSK + conv r1/2, RRC 0.35, 1 Msym/s (SRS-4 class)";
    case "gk2a_lrit"
        tm = ccsdsTMWaveformGenerator("ChannelCoding", "convolutional", "Modulation", "BPSK", "RolloffFactor", 0.5, "SamplesPerSymbol", 8);
        W.x = runTM(tm, 3); W.fs = 128e3*8; W.Rs = 128e3; W.sps = 8; W.beta = 0.5;
        W.note = "ccsdsTMWaveformGenerator BPSK + conv r1/2, RRC 0.5, 128 ksym/s (LRIT 64 kbps)";
    case "fsk115k"
        m = comm.FSKModulator('ModulationOrder', 2, 'FrequencySeparation', 254e3, 'SymbolRate', 115e3, ...
            'SamplesPerSymbol', 16, 'ContinuousPhase', true, 'BitInput', true);
        W.x = m(randi([0 1], 1500, 1)); W.fs = 115e3*16; W.Rs = 115e3; W.kind = "fm";
        W.note = "comm.FSKModulator 2-FSK 115 kbps, +/-127 kHz (from 420KF1D: B = 2M + 2DK, K = 1.2)";
    case {"dvbs2_2M", "dvbs2_30M"}
        if key == "dvbs2_2M", Rs = 2e6; b = 0.35; sps = 4; else, Rs = 30e6; b = 0.2; sps = 2; end
        s2 = dvbs2WaveformGenerator("MODCOD", 4, "RolloffFactor", b, "SamplesPerSymbol", sps);
        np = s2.MinNumPackets; pk = [repmat(int2bit(71, 8), 1, np); randi([0 1], 187*8, np)];
        W.x = s2(pk(:)); W.fs = Rs*sps; W.Rs = Rs; W.sps = sps; W.beta = b;
        W.note = sprintf('dvbs2WaveformGenerator MODCOD 4 (QPSK 1/2), RRC %.2f, TS packets with 0x47 sync', b);
    case "iridium"
        Rs = 25e3; sps = 10; fs = Rs*sps;
        sym = pskmod(mod(cumsum(randi([0 3], 207, 1)), 4), 4, pi/4);
        tx = comm.RaisedCosineTransmitFilter('RolloffFactor', 0.4, 'FilterSpanInSymbols', 10, 'OutputSamplesPerSymbol', sps);
        bur = tx([sym; zeros(10, 1)]);
        fr = zeros(round(0.09*fs), 1); x = [fr; fr];
        x(1:numel(bur)) = bur; x(numel(fr) + (1:numel(bur))) = bur;
        W.x = x; W.fs = fs; W.Rs = Rs; W.sps = sps; W.beta = 0.4; W.kind = "burst-psk";
        W.note = "DE-QPSK (pskmod on cumulative symbols) + RRC 0.4; 8.28 ms burst per 90 ms frame";
    case "orbcomm"
        ph = cumsum((2*randi([0 1], 2000, 1) - 1)*pi/2);
        tx = comm.RaisedCosineTransmitFilter('RolloffFactor', 0.4, 'FilterSpanInSymbols', 10, 'OutputSamplesPerSymbol', 10);
        W.x = tx([exp(1j*ph); zeros(10, 1)]); W.fs = 48e3; W.Rs = 4800; W.sps = 10; W.beta = 0.4;
        W.note = "SDPSK (+/-90 deg per bit) + RRC 0.4, 4800 bps";
    case "thuraya"
        tx = comm.RaisedCosineTransmitFilter('RolloffFactor', 0.35, 'FilterSpanInSymbols', 10, 'OutputSamplesPerSymbol', 8);
        W.x = tx([dpskmod(randi([0 3], 2000, 1), 4, pi/4); zeros(10, 1)]); W.fs = 23.4e3*8; W.Rs = 23.4e3; W.sps = 8; W.beta = 0.35;
        W.note = "dpskmod pi/4-DQPSK + RRC 0.35, 23.4 ksym/s (GMR-1)";
    case "gps_l1"
        g = gpsWaveformGenerator("SignalType", "legacy", "PRNID", 5, "SampleRate", 10.23e6);
        x = g(randi([0 1], 1, 1)); W.x = x(1:round(10e-3*10.23e6)); W.fs = 10.23e6; W.Rs = 1.023e6; W.kind = "gnss";
        W.note = "gpsWaveformGenerator legacy (L1 C/A, PRN 5)";
    case "gps_l5"
        g = gpsWaveformGenerator("SignalType", "l5", "PRNID", 5, "SampleRate", 30.69e6);
        x = g(randi([0 1], 2, 1)); W.x = x(1:round(4e-3*30.69e6)); W.fs = 30.69e6; W.Rs = 10.23e6; W.kind = "gnss";
        W.note = "gpsWaveformGenerator l5 (PRN 5)";
    case "gal_e1"
        g = galileoWaveformGenerator("SignalType", "e1", "SVID", 11, "SampleRate", 20.46e6);
        W.x = g(randi([0 1], 2, 1)); W.fs = 20.46e6; W.Rs = 1.023e6; W.kind = "gnss";
        W.note = "galileoWaveformGenerator e1 (CBOC, SVID 11)";
    case "gal_e5"
        g = galileoWaveformGenerator("SignalType", "e5", "SVID", 11, "SampleRate", 122.76e6);
        x = g({randi([0 1], 1, 1), randi([0 1], 5, 1)}); W.x = x(1:round(2e-3*122.76e6)); W.fs = 122.76e6; W.Rs = 10.23e6; W.kind = "gnss";
        W.note = "galileoWaveformGenerator e5 (AltBOC, SVID 11)";
    case "bds_b1i"
        gs = comm.GoldSequence('FirstPolynomial', 'z^11+z^10+z^9+z^8+z^7+z+1', 'SecondPolynomial', 'z^11+z^9+z^8+z^5+z^4+z^3+z^2+z+1', ...
            'FirstInitialConditions', [0 1 0 1 0 1 0 1 0 1 0], 'SecondInitialConditions', [0 1 0 1 0 1 0 1 0 1 0], 'Index', 0, 'SamplesPerFrame', 2046);
        c = 1 - 2*double(gs());
        W.x = complex(repelem(repmat(c, 5, 1), 10)); W.fs = 20.46e6; W.Rs = 2.046e6; W.kind = "gnss";
        W.note = "comm.GoldSequence 11-stage, 2046 chips at 2.046 Mcps (B1I-like; not the ICD code phase)";
    case "bds_b3i"
        pn = comm.PNSequence('Polynomial', 'z^13 + z^4 + z^3 + z + 1', 'InitialConditions', ones(1, 13), 'SamplesPerFrame', 10230);
        c = 1 - 2*double(pn());
        W.x = complex(repelem(repmat(c, 2, 1), 4)); W.fs = 40.92e6; W.Rs = 10.23e6; W.kind = "gnss";
        W.note = "comm.PNSequence 13-stage, 10230 chips at 10.23 Mcps (B3I-like)";
    case "glo_l1"
        pn = comm.PNSequence('Polynomial', 'z^9 + z^5 + 1', 'InitialConditions', ones(1, 9), 'SamplesPerFrame', 511);
        c = 1 - 2*double(pn());
        W.x = complex(repelem(repmat(c, 10, 1), 10)); W.fs = 5.11e6; W.Rs = 0.511e6; W.kind = "gnss";
        W.note = "comm.PNSequence z^9+z^5+1 (GLONASS ST code, 511 chips, 0.511 Mcps); no GLONASS generator in MATLAB";
    case "qzss_l6"
        pn = comm.PNSequence('Polynomial', 'z^10 + z^3 + 1', 'InitialConditions', ones(1, 10), 'SamplesPerFrame', 5115);
        c = 1 - 2*double(pn());
        W.x = complex(repelem(c, 8)); W.fs = 20.46e6; W.Rs = 2.5575e6; W.kind = "gnss";
        W.note = "BPSK chips 2.5575 Mcps (spectrum of CASM L6; casmmod builds the full code-shift keying)";
    otherwise
        error('ff_satsig:template', 'Unknown template "%s".', key);
end
W.x = W.x(:);
end

function x = runTM(tm, n)
x = [];
for k = 1:n, x = [x; tm(randi([0 1], tm.NumInputBits, 1))]; end %#ok<AGROW>
end

function y = fmBB(m, fs, dev)
y = exp(1j*2*pi*dev*cumsum(m(:))/fs);
end

function bits = ax25Bits(dst, src, info)
% AX.25 UI frame -> FCS (comm.CRCGenerator, CRC-16/X.25) -> HDLC bit stuffing + flags (LSB first).
payload = [encCall(dst, false) encCall(src, true) 3 240 double(char(info))];
crc = comm.CRCGenerator('Polynomial', 'z^16 + z^12 + z^5 + 1', 'InitialConditions', 1, 'DirectMethod', true, ...
    'ReflectInputBytes', true, 'ReflectChecksums', true, 'FinalXOR', 1);
cw = crc(reshape(int2bit(payload, 8), [], 1));
fcs = bit2int(cw(end-15:end), 16);
bytes = [payload bitand(fcs, 255) bitshift(fcs, -8)];
b = reshape(int2bit(bytes, 8, false), 1, []);
o = zeros(1, ceil(numel(b)*1.25)); n = 0; run = 0;
for k = 1:numel(b)
    n = n + 1; o(n) = b(k);
    if b(k) == 1, run = run + 1; if run == 5, n = n + 1; o(n) = 0; run = 0; end, else, run = 0; end
end
flag = [0 1 1 1 1 1 1 0];
bits = [repmat(flag, 1, 8) o(1:n) repmat(flag, 1, 2)]';
end

function out = encCall(call, last)
c = double(upper(char(pad(string(call), 6)))); c = c(1:6);
out = [bitshift(c, 1), bitor(96, double(last))];
end

function y = nrziEnc(b)
y = zeros(size(b)); s = 0;
for k = 1:numel(b), if b(k) == 0, s = 1 - s; end, y(k) = s; end
end

function [env, tdot] = morseEnv(msg, wpm, fs)
M = containers.Map( ...
    {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z', ...
    '0','1','2','3','4','5','6','7','8','9','=','/'}, ...
    {'.-','-...','-.-.','-..','.','..-.','--.','....','..','.---','-.-','.-..','--','-.','---','.--.','--.-', ...
    '.-.','...','-','..-','...-','.--','-..-','-.--','--..','-----','.----','..---','...--','....-','.....', ...
    '-....','--...','---..','----.','-...-','-..-.'});
tdot = 1.2/wpm; u = [];
for ch = upper(char(msg))
    if ch == ' ', u = [u zeros(1, 4)]; continue, end %#ok<AGROW>
    if ~isKey(M, ch), continue, end
    for el = M(ch), u = [u ones(1, 1 + 2*(el == '-')) 0]; end %#ok<AGROW>
    u = [u 0 0]; %#ok<AGROW>
end
env = repelem(u', round(tdot*fs));
w = round(0.005*fs); r = 0.5 - 0.5*cos(pi*(0:w-1)'/w);
env = conv(env, r/sum(r), 'same');
end

function y = loraCSS(SF, BW, os, syms)
N = 2^SF; fs = BW*os; n = (0:N*os-1)'; y = zeros(N*os*numel(syms), 1);
for k = 1:numel(syms)
    f = mod(syms(k) + n/os, N)*BW/N - BW/2;
    y((k-1)*N*os + (1:N*os)) = exp(1j*2*pi*cumsum(f)/fs);
end
end

function tpl = buildTemplates(file)
keys = templateKeys(); tpl = struct();
for k = keys
    tic; W = synthWave(k, 1); x = W.x;
    o = obw(complex(x), W.fs);
    rb = max(o/300, W.fs/2^18);
    L = 2^nextpow2(W.fs/rb); L = min(L, 2^floor(log2(numel(x)/4))); L = max(L, 64);
    [p, f] = pwelch(x, hann(L), L/2, L, W.fs, 'centered');
    tpl.(k) = struct('f', single(f), 'psd', single(p/max(p)), 'df', W.fs/L, 'obw', o, 'fs', W.fs, 'Rs', W.Rs, ...
        'kind', W.kind, 'note', W.note);
    fprintf('template %-15s OBW99 %10s  df %8s  (%.1f s)\n', k, fmtHz(o), fmtHz(W.fs/L), toc);
end
save(file, 'tpl', '-v7');
fprintf('Saved %s\n', file);
end

function tpl = loadTemplates()
persistent T
if isempty(T)
    f = fullfile(supportDir(), 'satsig_templates.mat');
    if isfile(f)
        S = load(f, 'tpl'); T = S.tpl;
    else
        warning('ff_satsig:templates', 'Template file missing: building %s now (~1 min).', f);
        T = buildTemplates(f);
    end
end
tpl = T;
end

function [u, D1] = tplDisplay(T, rbw)
% Displayed power (mW per 1 mW of signal) vs offset u for a Gaussian RBW (|H|^2 = 0.5 at +-RBW/2).
f = double(T.f(:)); p = double(T.psd(:)); p = max(p, 0);
keep = find(p >= 1e-5*max(p));
lo = f(keep(1)) - 4*rbw; hi = f(keep(end)) + 4*rbw;
h = min(T.df, rbw/5); n = ceil((hi - lo)/h) + 1;
if n > 4e5, h = (hi - lo)/4e5; n = 400001; end
u = lo + (0:n-1)'*h;
pu = interp1(f, p, u, 'linear', 0); pu = pu/sum(pu);          % power fraction per grid cell
if rbw > 2*h
    m = ceil(3*rbw/h); g = exp(-4*log(2)*((-m:m)'*h/rbw).^2);
    D1 = conv(pu, g, 'same');
else
    D1 = pu/h*1.0645*rbw;
end
end

function W = tplWidths(tpl, key, rbw)
% -X dB width (X = 3:0.5:20) of the RBW-smoothed template (cached).
persistent M
if isempty(M), M = containers.Map(); end
k = sprintf('%s|%.6g', key, rbw);
if isKey(M, k), W = M(k); return, end
[u, D1] = tplDisplay(tpl.(key), rbw); h = u(2) - u(1); s = D1/max(D1);
X = 3:0.5:20; W = zeros(size(X));
for i = 1:numel(X)
    a = find(s >= 10^(-X(i)/10)); W(i) = (a(end) - a(1))*h + h;
end
M(k) = W;
end

%% ======================================================================================
%% IDENTIFY
%% ======================================================================================
function [T, info] = cmdIdentify(data, varargin)
opt = baseOptions();
opt.Threshold = 6; opt.FreqTolPPM = 5; opt.FreqTolHz = 500; opt.TimeTol = []; opt.RxGain = 0;
opt.PriorFixed = 0.5; opt.TrackGap = 180; opt.MinDetections = []; opt.Demo = ""; opt.Seed = 1;
opt.Bands = []; opt.Targets = []; opt.Duration = 300; opt.Interval = 0; opt.RBW = []; opt.Points = 401;
opt.Preamp = true; opt.Detector = "POS"; opt.TraceType = "CLRW"; opt.HoldTime = 10; opt.Passive = false;
opt.Port = 5025; opt.SaveTo = ""; opt.Stack = true; opt.StackZ = 6;
userSN = any(strcmpi(string(varargin(1:2:end)), "IncludeSatNOGS"));
opt = parseOpts(opt, varargin);
info = struct('notes', strings(0, 1), 'bands', [], 'detections', [], 'tracks', [], 'hypotheses', {{}}, ...
    'stack', [], 'stackProfiles', [], 'truth', [], 'figs', gobjects(0), 'files', strings(0, 1), 'settings', opt);
tpl = loadTemplates();
C = builtinCatalog(opt);
isDemo = (isstring(data) || ischar(data)) && strcmpi(string(data), "demo");
if isDemo
    if strlength(string(opt.Demo)) == 0, opt.Demo = "origami"; end
    if ~userSN, opt.IncludeSatNOGS = isfile(fullfile(supportDir(), 'satnogs_subset.csv')); end
    [R, truth, dn] = simulateDemo(opt.Demo, opt, tpl);
    info.truth = truth; info.notes = [info.notes(:); dn(:)];
elseif isstruct(data)
    R = normalizeRec(data, opt);
elseif isa(data, 'FieldFox') || isLiveAddress(data)
    [R, ln] = acquireLive(data, opt, C); info.notes = [info.notes(:); ln(:)];
else
    files = string(data); R = [];
    for f = files(:)'
        if ~isfile(f), error('ff_satsig:file', 'Recording not found: %s', f); end
        S = load(f); r = normalizeRec(S, opt); r.source = f; R = [R, r]; %#ok<AGROW>
    end
end
if isempty(R), error('ff_satsig:nodata', 'No sweeps to analyse.'); end
% hypothesis rows: catalog (+ SatNOGS)
if opt.IncludeSatNOGS
    [S, note] = satnogsTable(opt, false); info.notes(end+1) = note;
    if ~isempty(S), C = [C; satnogsRows(S, C)]; end
end
excluded = C(C.Active <= 0, :);
if ~opt.IncludeInactive, C = C(C.Active > 0, :); else, C.Active(C.Active <= 0) = 0.02; end
rows = {}; allTr = [];
for b = 1:numel(R)
    band = bandInfo(R(b));
    [Dt, bn] = detectBand(R(b), band, opt); info.notes = [info.notes(:); bn(:)];
    tr = trackDetections(Dt, band, opt);
    minDet = opt.MinDetections; if isempty(minDet), minDet = 1 + (band.nSweeps >= 10); end
    nSingle = nnz(arrayfun(@(t) numel(t.idx) < minDet, tr));
    tr = tr(arrayfun(@(t) numel(t.idx) >= minDet, tr));
    if nSingle > 0, info.notes(end+1) = sprintf("Band %d: %d track(s) with fewer than %d detections ignored (noise/glitches).", b, nSingle, minDet); end
    sig = arrayfun(@(t) trackSignificant(t, Dt, band), tr);
    if any(~sig)
        info.notes(end+1) = sprintf(['Band %d: %d sparse track(s) at the detection threshold dropped (chance alignments of ' ...
            'noise peaks are likely in long recordings); the Doppler stacking still searches every candidate satellite.'], b, nnz(~sig));
        tr = tr(sig);
    end
    say(opt, 'Band %d: %s - %s, %d sweeps, %d detections, %d tracks\n', b, fmtHz(band.f(1)), fmtHz(band.f(end)), ...
        band.nSweeps, numel(Dt.tau), numel(tr));
    % excluded (inactive) rows that would have matched
    for k = 1:numel(tr)
        fm = median(Dt.fc(tr(k).idx));
        ex = excluded(abs(excluded.F0_Hz - fm) <= excluded.MaxDoppler_Hz + excluded.OBW_Hz/2 + 20e3 | ...
            abs(excluded.FAlt_Hz - fm) <= excluded.MaxDoppler_Hz + 20e3, :);
        for e = 1:height(ex)
            note = sprintf("%s (%s) matches %s by frequency but is excluded: %s", ex.Satellite(e), ex.Id(e), fmtHz(fm), ex.Status(e));
            if ~any(info.notes == note), info.notes(end+1) = note; end
        end
    end
    [H, pn, P] = scoreBand(tr, Dt, band, C, tpl, opt); info.notes = [info.notes(:); pn(:)];
    [tr, H, mn] = mergeTracks(tr, H, Dt, band, C, tpl, opt, P); info.notes = [info.notes(:); mn(:)];
    for k = 1:numel(tr)
        tr(k).band = b; tr(k).hyp = H{k};
        rows(end+1, :) = trackRow(tr(k), Dt, band, H{k}, b, size(rows, 1) + 1); %#ok<AGROW>
    end
    if opt.Stack
        [Sb, prof, sn] = stackSearch(band, Dt, C, opt, b, tr); info.notes = [info.notes(:); sn(:)];
        if isempty(info.stack), info.stack = Sb; else, info.stack = [info.stack; Sb]; end
        info.stackProfiles = [info.stackProfiles, prof];
    end
    info.bands = [info.bands, band]; info.detections = [info.detections, Dt]; info.hypotheses = [info.hypotheses, H];
    if isempty(allTr), allTr = tr; else, allTr = [allTr, tr]; end %#ok<AGROW>
end
vn = {'Track','Band','Freq_MHz','Drift_kHz','Rate_Hz_s','NDet','OnFrac','Level_dBm','SNR_dB','BW_Hz','First_UTC', ...
    'Last_UTC','Best','NORAD','Confidence','Second','FitRMS_Hz','TimeOffset_s','FreqOffset_Hz','El_deg','Evidence','Flag'};
if isempty(rows), T = cell2table(cell(0, numel(vn)), 'VariableNames', vn);
else, T = cell2table(rows, 'VariableNames', vn); end
info.tracks = allTr; info.recordings = R;
if isDemo, info.notes = [info.notes(:); demoCheck(T, info.truth, info.stack)]; end
if opt.Plot, info.figs = plotIdentify(T, info, opt); end
if strlength(opt.OutDir) > 0
    writetable(T, fullfile(opt.OutDir, "satsig_identify.csv"));
    if ~isempty(info.stack) && height(info.stack) > 0, writetable(info.stack, fullfile(opt.OutDir, "satsig_stack.csv")); end
    Tsave = T; infoSave = rmfield(info, 'figs'); %#ok<NASGU>
    save(fullfile(opt.OutDir, "satsig_identify.mat"), 'Tsave', 'infoSave');
    info.files = [info.files; fullfile(opt.OutDir, "satsig_identify.csv"); fullfile(opt.OutDir, "satsig_identify.mat")];
    if isDemo
        for b = 1:numel(R)
            rec = R(b); f = rec.f; P = rec.P; t = rec.t; settings = rec.settings; idn = rec.idn; %#ok<NASGU>
            fn = fullfile(opt.OutDir, sprintf('satsig_demo_rec_b%d.mat', b));
            save(fn, 'f', 'P', 't', 'settings', 'idn'); info.files(end+1) = fn;
        end
    end
end
end

function tf = isLiveAddress(x)
tf = (isstring(x) || ischar(x)) && isscalar(string(x)) && ...
    (~isempty(regexp(char(x), '^\d{1,3}(\.\d{1,3}){3}$', 'once')) || strcmpi(string(x), "localhost"));
end

function R = normalizeRec(S, opt)
R = [];
for k = 1:numel(S)
    s = S(k);
    if ~all(isfield(s, {'f', 'P'})), error('ff_satsig:rec', 'A recording needs fields f (Hz) and P (dBm, sweeps x points).'); end
    f = double(s.f(:))'; P = double(s.P);
    if size(P, 2) ~= numel(f) && size(P, 1) == numel(f), P = P.'; end
    if size(P, 2) ~= numel(f), error('ff_satsig:rec', 'P must be sweeps x numel(f).'); end
    if isfield(s, 't') && ~isempty(s.t)
        t = s.t(:);
        if isnumeric(t), t = datetime(t, 'ConvertFrom', 'datenum'); end
        t = toUTC(t, opt.TimeZone);
    else
        t0 = opt.Time; if isempty(t0), error('ff_satsig:rec', 'The recording has no time stamps: give ''Time'' (start) for the prediction.'); end
        t = toUTC(datetime(t0), opt.TimeZone) + seconds((0:size(P, 1) - 1)');
    end
    if numel(t) ~= size(P, 1), t = t(1) + seconds((0:size(P, 1) - 1)'*mean(diff(seconds(t - t(1))))); end
    st = struct(); if isfield(s, 'settings'), st = s.settings; end
    rbw = NaN; if isfield(st, 'rbw'), rbw = double(st.rbw); end
    if ~isfinite(rbw), rbw = rbwStep((f(end) - f(1))/max(numel(f) - 1, 1)); end
    det = "POS"; if isfield(st, 'detector'), det = upper(string(st.detector)); end
    idn = ""; if isfield(s, 'idn'), idn = string(s.idn); end
    ok = all(isfinite(P), 2);
    r = struct('f', f, 'P', P(ok, :), 't', t(ok), 'settings', st, 'rbw', rbw, 'detector', det, 'idn', idn, 'source', "struct");
    R = [R, r]; %#ok<AGROW>
end
end

function band = bandInfo(r)
band = struct('f', r.f, 'df', (r.f(end) - r.f(1))/max(numel(r.f) - 1, 1), 'rbw', r.rbw, 'span', r.f(end) - r.f(1), ...
    't0', r.t(1), 'tau', seconds(r.t - r.t(1)), 'nSweeps', size(r.P, 1), 'P', r.P, 'nf', [], 'detector', r.detector, 'source', r.source);
end

%% ----- detection ---------------------------------------------------------------------------
function [D, notes] = detectBand(r, band, opt)
notes = strings(0, 1);
D = struct('s', [], 'tau', [], 'fc', [], 'fpk', [], 'lvl', [], 'snr', [], 'w', [], 'X', [], 'pw', [], 'edge', [], 'nf', [], 'thr', []);
ws = warning('off', 'signal:findpeaks:largeMinPeakHeight'); restoreW = onCleanup(@() warning(ws)); %#ok<NASGU>
f = band.f(:); df = band.df; rbw = band.rbw; N = numel(f);
if N < 11, notes(end+1) = "Too few trace points for detection."; return, end
enbw = 1.0645*rbw; gapBins = max(2, ceil(rbw/df));
for s = 1:size(r.P, 1)
    p = r.P(s, :)';
    m = median(p); sd = 1.4826*median(abs(p - m));
    for it = 1:6
        k = p < m + 2.5*max(sd, 0.2); m = median(p(k)); sd = 1.4826*median(abs(p(k) - m));
    end
    sd = max(sd, 0.25); nf = m;
    D.nf(end+1, 1) = nf;
    plin = 10.^(p/10); nlin = 10^(nf/10);
    ps = 10*log10(movmean(plin, 3));
    pk = []; h1 = nf + max(opt.Threshold, 3.5*sd); h2 = nf + max(opt.Threshold - 1, 2.5*sd);   % 3-pt mean: lower noise
    if max(p) > h1, [~, a] = findpeaks(p, 'MinPeakHeight', h1, 'MinPeakProminence', 3); pk = [pk; a(:)]; end
    if max(ps) > h2, [~, a] = findpeaks(ps, 'MinPeakHeight', h2, 'MinPeakProminence', 2); pk = [pk; a(:)]; end
    if isempty(pk), continue, end
    pk = unique(pk);
    lev = nf + max(1.5, 1.2*sd);
    reg = zeros(0, 2);
    for i = pk'
        a = i; while a > 1 && ps(a - 1) > lev, a = a - 1; end
        b = i; while b < N && ps(b + 1) > lev, b = b + 1; end
        reg(end+1, :) = [max(1, min(a, i - 1)), min(N, max(b, i + 1))]; %#ok<AGROW>
    end
    reg = sortrows(reg); mr = reg(1, :);
    for i = 2:size(reg, 1)
        if reg(i, 1) <= mr(end, 2) + gapBins, mr(end, 2) = max(mr(end, 2), reg(i, 2)); else, mr(end+1, :) = reg(i, :); end %#ok<AGROW>
    end
    for i = 1:size(mr, 1)
        ix = (mr(i, 1):mr(i, 2))';
        lin = max(plin(ix) - nlin, 0);
        if sum(lin) <= 0, continue, end
        [lv, ipk] = max(p(ix)); snr = lv - nf;
        if snr < opt.Threshold - 1, continue, end
        X = min(10, max(3, snr - 3));
        if numel(ix) > 6, pw = ps(ix); else, pw = p(ix); end
        ab = find(pw >= max(pw) - X);
        w = (ab(end) - ab(1) + 1)*df;
        if numel(ix) <= 2, P = 10*log10(max(lin)); else, P = 10*log10(sum(lin)*df/enbw); end
        D.s(end+1, 1) = s; D.tau(end+1, 1) = band.tau(s); D.fc(end+1, 1) = sum(f(ix).*lin)/sum(lin);
        D.fpk(end+1, 1) = f(ix(ipk)); D.lvl(end+1, 1) = lv; D.snr(end+1, 1) = snr; D.w(end+1, 1) = w;
        D.X(end+1, 1) = X; D.pw(end+1, 1) = P; D.edge(end+1, 1) = ix(1) == 1 || ix(end) == N;
        D.thr(end+1, 1) = h1 - nf;
    end
end
end

function ok = trackSignificant(t, D, band)
% A track is kept when it is dense (>= 20 % of the sweeps it spans), long (>= 10 detections) or
% clearly above the per-sweep detection threshold (median excess >= 3 dB).
i = t.idx(:); tau = D.tau(i);
spanned = nnz(band.tau >= tau(1) & band.tau <= tau(end));
ok = numel(i) >= 10 || numel(i)/max(spanned, 1) >= 0.2 || median(D.snr(i) - D.thr(i)) >= 3;
end

%% ----- tracking ----------------------------------------------------------------------------
function tr = trackDetections(D, band, opt)
tr = struct('idx', {}, 'lastTau', {}, 'lastF', {});
if isempty(D.tau), return, end
fmid = mean(band.f); maxRate = 4e-7*fmid + 20; base = max(3*band.df, 2*band.rbw);
for s = unique(D.s)'
    ii = find(D.s == s); t = D.tau(ii(1));
    act = find(arrayfun(@(T) t - T.lastTau <= opt.TrackGap && t > T.lastTau, tr));
    cost = inf(numel(ii), numel(act));
    for a = 1:numel(act)
        [fp, g] = predictTrack(tr(act(a)), D, t, base, maxRate);
        cost(:, a) = abs(D.fc(ii) - fp)/g;
    end
    cost(cost > 1) = inf; used = false(numel(ii), 1);
    while any(isfinite(cost(:)))
        [~, j] = min(cost(:)); [r, c] = ind2sub(size(cost), j);
        k = act(c); tr(k).idx(end+1) = ii(r); tr(k).lastTau = t; tr(k).lastF = D.fc(ii(r));
        used(r) = true; cost(r, :) = inf; cost(:, c) = inf;
    end
    for r = find(~used)'
        tr(end+1) = struct('idx', ii(r), 'lastTau', t, 'lastF', D.fc(ii(r))); %#ok<AGROW>
    end
end
end

function [fp, g] = predictTrack(T, D, t, base, maxRate)
dt = t - T.lastTau;
k = T.idx(D.tau(T.idx) >= T.lastTau - 120);
if numel(k) >= 3 && range(D.tau(k)) >= 10
    c = polyfit(D.tau(k) - T.lastTau, D.fc(k), 1);
    r = D.fc(k) - polyval(c, D.tau(k) - T.lastTau);
    fp = polyval(c, dt); g = base + 0.3*maxRate*dt + 2*std(r);
else
    fp = T.lastF; g = base + maxRate*dt;
end
end

%% ----- scoring -----------------------------------------------------------------------------
function [H, notes, P] = scoreBand(tr, D, band, C, tpl, opt, P)
% Score every track of one band. P (predictions) is computed here unless given (re-scoring).
notes = strings(0, 1); H = cell(1, numel(tr));
tA = band.t0 - seconds(400); tB = band.t0 + seconds(max(band.tau) + 400);
if nargin < 7, P = predictSats([], tA, tB, opt); havP = false; else, havP = true; end
if isempty(tr), return, end
% frequency-matching rows per track
match = cell(1, numel(tr));
for k = 1:numel(tr)
    fm = median(D.fc(tr(k).idx));
    osc = opt.FreqTolPPM*1e-6*C.F0_Hz + opt.FreqTolHz;
    win = C.MaxDoppler_Hz + C.OBW_Hz/2 + 4*osc;
    spec = abs(C.F0_Hz - fm) <= win | abs(C.FAlt_Hz - fm) <= win;
    inBand = isfinite(C.BandLo_Hz) & fm >= C.BandLo_Hz - C.MaxDoppler_Hz - C.OBW_Hz/2 & fm <= C.BandHi_Hz + C.MaxDoppler_Hz + C.OBW_Hz/2;
    match{k} = find(spec | inBand);
end
if ~havP
    Cm = C(unique(vertcat(match{:})), :);
    grp = unique([Cm.Group(Cm.Kind == "constellation"); Cm.Group(Cm.Kind == "sat" & Cm.Group ~= "")]);
    needAll = Cm.NORAD(Cm.Kind == "sat");
    needCat = Cm.NORAD(Cm.Kind == "sat" & Cm.Group ~= "satnogs");     % CATNR queries only for catalog rows
    if ~isempty(needAll) || ~isempty(grp)
        [E, n1] = getTLEs(opt, grp, needCat); notes = [notes; n1(:)];
        keepE = ismember([E.norad], needAll);
        for r = find(Cm.Kind == "constellation")'
            keepE = keepE | ([E.group] == Cm.Group(r) & ~cellfun(@isempty, regexp(cellstr([E.name]), Cm.TleMatch(r), 'once')));
        end
        E = E(keepE);
        nMiss = numel(setdiff(unique(needAll(:)'), [E.norad]));
        if nMiss > 0, notes(end+1) = sprintf("%d candidate satellite(s) without a TLE were skipped (SatNOGS rows use GROUP=satnogs only).", nMiss); end
        if ~isempty(E)
            say(opt, '  predicting %d satellites (satelliteScenario, SGP4) ...\n', numel(E));
            P = predictSats(E, tA, tB, opt);
            for k = 1:numel(P), P(k).tau = seconds(P(k).t - band.t0); end     % same time base as the tracks
        end
    end
end
for k = 1:numel(tr)
    H{k} = scoreTrack(tr(k), D, band, C(match{k}, :), P, tpl, opt);
end
end

function [tr, H, notes] = mergeTracks(tr, H, D, band, C, tpl, opt, P)
% Association by hypothesis: two tracks whose best hypothesis is the same satellite (e.g. a beacon
% seen every 60 s near TCA, where the Doppler rate breaks the linear track gate) are merged when the
% merged track is still best explained by that satellite.
notes = strings(0, 1);
again = true;
while again
    again = false; n = numel(tr);
    lab = strings(n, 1);
    for k = 1:n
        h = H{k}(1);
        if ismember(h.kind, ["sat", "constellation"]) && h.post >= 0.5, lab(k) = h.label; end
    end
    for a = 1:n - 1
        for b = a + 1:n
            if lab(a) == "" || lab(a) ~= lab(b), continue, end
            m = tr(a); m.idx = sort([tr(a).idx(:); tr(b).idx(:)])';
            if numel(unique(D.s(m.idx))) < numel(m.idx), continue, end     % two detections in one sweep
            m.lastTau = max(tr(a).lastTau, tr(b).lastTau);
            Hm = scoreBand(m, D, band, C, tpl, opt, P);
            okFit = Hm{1}(1).rms <= 1.5*max(H{a}(1).rms, H{b}(1).rms) + 2*max(band.df, band.rbw);
            if Hm{1}(1).label == lab(a) && Hm{1}(1).post >= 0.5 && okFit
                notes(end+1) = sprintf("Band: tracks at %s and %s merged (both best explained by %s; merged fit rms %.0f Hz).", ...
                    fmtHz(median(D.fc(tr(a).idx))), fmtHz(median(D.fc(tr(b).idx))), lab(a), Hm{1}(1).rms); %#ok<AGROW>
                tr(a) = m; H{a} = Hm{1}; tr(b) = []; H(b) = [];
                again = true; break
            end
        end
        if again, break, end
    end
end
end

function H = scoreTrack(trk, D, band, Ck, P, tpl, opt)
i = trk.idx(:); tau = D.tau(i); y = D.fc(i); N = numel(i);
Neff = min(N, 30); NeffW = min(N, 6); sig0 = max(band.df/sqrt(12), 5);
dur = tau(end) - tau(1); doFit = N >= 4 && dur >= 30;
trkS = struct('X', D.X(i), 'w', D.w(i), 'snr', D.snr(i), 'edge', D.edge(i), 'pw', D.pw(i));
spec = alignedSpectrum(trk, D, band);
H = struct('label', {}, 'kind', {}, 'norad', {}, 'rows', {}, 'keys', {}, 'logPrior', {}, 'logLoc', {}, ...
    'logRes', {}, 'logPen', {}, 'logVis', {}, 'logWidth', {}, 'logShape', {}, 'logLevel', {}, 'total', {}, ...
    'post', {}, 'dt', {}, 'dfo', {}, 'rms', {}, 'el', {}, 'curveTau', {}, 'curveF', {}, 'modes', {}, 'f0', {});
resTerm = @(rss) -(Neff/2)*log(rss/N + sig0^2);
pen = @(kp) -(kp/2)*log(max(Neff, 1));
% --- fixed frequency (terrestrial / not a satellite)
ratePen = 0;
if N >= 3 && dur >= 10
    c = polyfit(tau - mean(tau), y, 1); res = y - polyval(c, tau - mean(tau)); kp = 2;
    se2 = (sum(res.^2)/max(N - 2, 1) + sig0^2)/sum((tau - mean(tau)).^2);    % slope variance
    ratePen = -0.5*c(1)^2/(2^2 + se2);              % terrestrial oscillators drift << 2 Hz/s
else
    res = y - mean(y); kp = 1;
end
llw = widthLL("cw16", trkS, band, tpl); llu = widthLL("", trkS, band, tpl);
lw = (NeffW/N)*sum(logsumexp2(llw, llu) - log(2));
rc = shapeCorr(spec, "cw16", band, tpl);
H(end+1) = mkHyp("Fixed-frequency carrier (terrestrial / not a predicted satellite)", "fixed", NaN, [], "", ...
    log(opt.PriorFixed), -log(max(band.span, 1)), resTerm(sum(res.^2)), pen(kp) + ratePen, 0, lw, ...
    log(0.5*exp(1.5*(rc - 1)) + 0.5*exp(-0.75)), 0, NaN, NaN, sqrt(mean(res.^2)), NaN, [], [], "", mean(y));
% --- rows grouped into hypotheses
gen = find(Ck.Kind == "generic"); nGen = numel(gen);
for r = gen'
    ord = 0; if N >= 8 && dur >= 60, ord = 3; elseif N >= 5 && dur >= 30, ord = 2; end
    if ord > 0, [c, ~, mu] = polyfit(tau, y, ord); res = y - polyval(c, tau, [], mu); else, res = y - mean(y); end
    lw = (NeffW/N)*sum(widthLL(Ck.Template(r), trkS, band, tpl));
    B = Ck.BandHi_Hz(r) - Ck.BandLo_Hz(r) + 2*Ck.MaxDoppler_Hz(r);
    H(end+1) = mkHyp(Ck.Satellite(r) + " (class only, no orbit)", "generic", NaN, r, Ck.Template(r), ...
        log(Ck.Active(r)*Ck.SrcPrior(r)/max(nGen, 1)), -log(B), resTerm(sum(res.^2)), pen(ord + 1), 0, lw, ...
        1.5*(shapeCorr(spec, Ck.Template(r), band, tpl) - 1), 0, NaN, NaN, sqrt(mean(res.^2)), NaN, [], [], "", mean(y)); %#ok<AGROW>
end
% specific satellites: group rows of one NORAD with the same nominal frequency
sp = find(Ck.Kind == "sat");
done = false(size(sp));
for a = 1:numel(sp)
    if done(a), continue, end
    r = sp(a);
    same = sp(Ck.NORAD(sp) == Ck.NORAD(r) & abs(Ck.F0_Hz(sp) - Ck.F0_Hz(r)) < 2e3);
    done(ismember(sp, same)) = true;
    j = find([P.norad] == Ck.NORAD(r), 1);
    if isempty(j)
        continue                                   % no TLE (noted by getTLEs)
    end
    H(end+1) = satHyp(Ck, same, P(j), tau, y, trkS, spec, band, tpl, opt, doFit, resTerm, pen, 1); %#ok<AGROW>
end
% constellation members
for r = find(Ck.Kind == "constellation")'
    mem = find(arrayfun(@(p) p.group == Ck.Group(r) && ~isempty(regexp(p.name, Ck.TleMatch(r), 'once')), P));
    if isempty(mem), continue, end
    vis = arrayfun(@(p) any(interp1(p.tau, p.el, tau, 'linear', -90) >= opt.ElevationMask - 2), P(mem));
    mem = mem(vis);
    for j = mem
        H(end+1) = satHyp(Ck, r, P(j), tau, y, trkS, spec, band, tpl, opt, doFit, resTerm, pen, numel(mem)); %#ok<AGROW>
    end
end
% GEO by longitude
for r = find(Ck.Kind == "geo")'
    Pg = predictGeo(Ck.GeoLon(r), band.t0, band.t0 + seconds(max(tau) + 1), opt, Ck.Satellite(r));
    res = y - mean(y);
    lw = (NeffW/N)*sum(widthLL(Ck.Template(r), trkS, band, tpl));
    H(end+1) = mkHyp(Ck.Satellite(r), "geo", NaN, r, Ck.Template(r), log(Ck.Active(r)), ...
        -log(Ck.BandHi_Hz(r) - Ck.BandLo_Hz(r)), resTerm(sum(res.^2)), pen(1), log(0.002 + 0.998*(Pg.el(1) >= opt.ElevationMask)), ...
        lw, 1.5*(shapeCorr(spec, Ck.Template(r), band, tpl) - 1), 0, NaN, NaN, sqrt(mean(res.^2)), Pg.el(1), [], [], "", mean(y)); %#ok<AGROW>
end
tot = [H.logPrior] + [H.logLoc] + [H.logRes] + [H.logPen] + [H.logVis] + [H.logWidth] + [H.logShape] + [H.logLevel];
pp = exp(tot - max(tot)); pp = pp/sum(pp);
for k = 1:numel(H), H(k).total = tot(k); H(k).post = pp(k); end
[~, o] = sort(pp, 'descend'); H = H(o);
end

function h = satHyp(Ck, rows, Pj, tau, y, trkS, spec, band, tpl, opt, doFit, resTerm, pen, nMem)
r = rows(1); N = numel(tau); NeffW = min(N, 6);
f0 = Ck.F0_Hz(r); isBandRow = isfinite(Ck.BandLo_Hz(r));
osc = opt.FreqTolPPM*1e-6*f0 + opt.FreqTolHz;
age = abs(days(band.t0 - Pj.epoch));
sigT = opt.TimeTol; if isempty(sigT), sigT = min(300, 5 + 3*age); end
Fd = griddedInterpolant(Pj.tau, Pj.dfac, 'pchip', 'nearest');
Fe = griddedInterpolant(Pj.tau, Pj.el, 'linear', 'nearest');
Fr = griddedInterpolant(Pj.tau, Pj.rng, 'linear', 'nearest');
Fq = griddedInterpolant(Pj.tau, Pj.rate, 'linear', 'nearest');
if isBandRow
    lo = (Ck.BandLo_Hz(r) - f0)/1e3; hi = (Ck.BandHi_Hz(r) - f0)/1e3;
else
    lo = -4*osc/1e3; hi = 4*osc/1e3;
end
model = @(q, t) (f0 + q(2)*1e3).*(1 + Fd(t + q(1)*10));
dtMax = min(360, max(120, 4*sigT));
if doFit
    best = []; bestR = Inf;
    o = optimoptions('lsqcurvefit', 'Display', 'off', 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-9, 'MaxFunctionEvaluations', 600);
    for dt0 = unique([-min(sigT, dtMax/2), 0, min(sigT, dtMax/2)])
        fch = mean(y./(1 + Fd(tau + dt0)));
        q0 = [dt0/10, min(max((fch - f0)/1e3, lo), hi)];
        try
            [q, rn] = lsqcurvefit(model, q0, tau, y, [-dtMax/10, lo], [dtMax/10, hi], o);
        catch
            q = q0; rn = sum((y - model(q0, tau)).^2);
        end
        rn = rn + sum(y*0) ; pr = rn/max(sig0Of(band), 1)^2*0;
        if rn + pr < bestR, bestR = rn + pr; best = q; end
    end
    q = best; dt = q(1)*10; dfo = q(2)*1e3; kp = 2;
    sigLoc = osc; logPt = -0.5*(dt/sigT)^2;
else
    dt = 0; dfo = mean(y - f0*(1 + Fd(tau))); dfo = min(max(dfo, lo*1e3), hi*1e3); kp = 1;
    sigLoc = sqrt(osc^2 + (f0*max(abs(Fq(tau)))*sigT)^2); logPt = 0;
end
q = [dt/10, dfo/1e3];
res = y - model(q, tau);
if isBandRow, logLoc = -log(Ck.BandHi_Hz(r) - Ck.BandLo_Hz(r));
else, logLoc = -0.5*(dfo/sigLoc)^2 - log(sigLoc*sqrt(2*pi)); end
el = Fe(tau + dt); frac = mean(el >= opt.ElevationMask - 1);
% level plausibility (loose): measured channel power vs a generous prediction
eirp = Ck.EIRP_dBm(r); prx = NaN(N, 1);
if isfinite(Ck.PrxMin_dBm(r)), prx = repmat(Ck.PrxMin_dBm(r) + 6, N, 1);
else
    if ~isfinite(eirp) && Ck.Alt_km(r) < 3000, eirp = 33; end
    if isfinite(eirp), prx = eirp - fspl(max(Fr(tau + dt), 1e3), physconst('LightSpeed')/f0) + 3; end
end
prx = prx + opt.RxGain;
ex = trkS.pw - prx - 20; ex(~isfinite(ex)) = 0;
logLevel = -(NeffW/N)*sum(0.5*max(0, ex/5).^2);
keys = unique(Ck.Template(rows)); keys = keys(keys ~= "");
if isempty(keys), keys = ""; end
ll = zeros(N, numel(keys));
for m = 1:numel(keys), ll(:, m) = widthLL(keys(m), trkS, band, tpl); end
lw = (NeffW/N)*sum(logsumexpRows(ll) - log(numel(keys)));
rho = max(arrayfun(@(k) shapeCorr(spec, k, band, tpl), keys));
[~, best] = max(ll, [], 2);
modes = strjoin(arrayfun(@(m) string(sprintf('%s x%d', modeName(Ck, rows, keys(m)), nnz(best == m))), 1:numel(keys)), ', ');
ct = linspace(tau(1) - 60, tau(end) + 60, 200)'; cf = model(q, ct);
lab = Ck.Satellite(r);
if isBandRow || nMem > 1, lab = Pj.name + " - " + Ck.Satellite(r); end
lab = lab + " " + fmtHz(f0) + " (" + strjoin(unique(Ck.Class(rows)), "/") + ")";
h = mkHyp(lab, ternary(nMem > 1 || isBandRow, "constellation", "sat"), Pj.norad, rows, strjoin(keys, ","), ...
    log(max(Ck.Active(r)*Ck.SrcPrior(r)/nMem, 1e-6)), logLoc, resTerm(sum(res.^2)), pen(kp) + logPt, ...
    log(0.002 + 0.998*frac), lw, 1.5*(rho - 1), logLevel, dt, dfo, sqrt(mean(res.^2)), median(el), ct, cf, modes, f0);
end

function s = sig0Of(band)
s = max(band.df/sqrt(12), 5);
end

function n = modeName(Ck, rows, key)
k = rows(find(Ck.Template(rows) == key, 1));
if isempty(k), n = "?"; else, n = Ck.Class(k); end
end

function h = mkHyp(label, kind, norad, rows, keys, lp, lloc, lres, lpen, lvis, lw, lsh, llev, dt, dfo, rms, el, ct, cf, modes, f0)
h = struct('label', string(label), 'kind', string(kind), 'norad', norad, 'rows', rows, 'keys', string(keys), ...
    'logPrior', lp, 'logLoc', lloc, 'logRes', lres, 'logPen', lpen, 'logVis', lvis, 'logWidth', lw, 'logShape', lsh, ...
    'logLevel', llev, 'total', NaN, 'post', NaN, 'dt', dt, 'dfo', dfo, 'rms', rms, 'el', el, 'curveTau', ct, ...
    'curveF', cf, 'modes', string(modes), 'f0', f0);
end

function ll = widthLL(key, trkS, band, tpl)
N = numel(trkS.w);
if key == "" || ~isfield(tpl, key)
    ll = repmat(-log(log(max(band.span/band.df, 3))), N, 1); return
end
W = tplWidths(tpl, key, band.rbw);
We = interp1(3:0.5:20, W, trkS.X, 'linear', 'extrap') + band.df;
s = log(1.4) + 0.5./max(trkS.snr - 3, 1) + 0.7*trkS.edge;
ll = -0.5*((log(trkS.w) - log(We))./s).^2 - log(s*sqrt(2*pi));
end

function spec = alignedSpectrum(trk, D, band)
i = trk.idx(:); df = band.df;
H = max(1.5*max(D.w(i)), 6*df); H = min(H, band.span/2);
v = (-ceil(H/df):ceil(H/df))'*df; acc = zeros(size(v)); n = 0;
for k = i'
    p = 10.^(band.P(D.s(k), :)'/10) - 10^(D.nf(D.s(k))/10);
    a = interp1(band.f(:), p, D.fc(k) + v, 'linear', NaN);
    if all(isfinite(a)), acc = acc + max(a, 0); n = n + 1; end
end
if n == 0, spec = struct('v', v, 's', []); else, spec = struct('v', v, 's', acc/n); end
end

function rho = shapeCorr(spec, key, band, tpl)
rho = 0.5;
if isempty(spec.s) || key == "" || ~isfield(tpl, key) || std(spec.s) == 0, return, end
[u, D1] = tplDisplay(tpl.(key), band.rbw);
st = interp1(u, D1, spec.v, 'linear', 0);
if band.df > band.rbw                       % bin-max of a peak detector broadens by ~1 bin
    st = movmax(st, 3);
end
if std(st) == 0, return, end
c = corrcoef(spec.s, st); rho = c(1, 2);
if ~isfinite(rho), rho = 0; end
end

function s = logsumexp2(a, b)
m = max(a, b); s = m + log(exp(a - m) + exp(b - m));
end

function s = logsumexpRows(L)
m = max(L, [], 2); s = m + log(sum(exp(L - m), 2));
end

function row = trackRow(trk, D, band, H, b, kk)
i = trk.idx(:); tau = D.tau(i); y = D.fc(i); N = numel(i);
rate = NaN; if N >= 3 && range(tau) > 0, c = polyfit(tau, y, 1); rate = c(1); end
nsw = nnz(band.tau >= tau(1) & band.tau <= tau(end));
h1 = H(1); sec = ""; if numel(H) > 1, sec = string(sprintf('%s (%.2f)', H(2).label, H(2).post)); end
ev = sprintf('loc %.1f res %.1f vis %.1f bw %.1f shape %.1f lvl %.1f prior %.1f', h1.logLoc, h1.logRes, h1.logVis, ...
    h1.logWidth, h1.logShape, h1.logLevel, h1.logPrior);
if h1.modes ~= "", ev = h1.modes + "; " + ev; end
if N < 4 || range(tau) < 30, ev = "no Doppler fit (" + N + " det.); " + ev; end
flag = emiFlag(median(y), h1);
if h1.kind == "fixed" && h1.rms > 3*max(band.df, band.rbw)
    chk = "CHECK: no hypothesis fits well (fixed-carrier rms " + round(h1.rms) + " Hz): weak/fragmented track";
    if startsWith(flag, "EMI"), flag = flag + "; " + chk; else, flag = chk; end
end
if ismember(h1.kind, ["sat", "constellation"]) && h1.rms > 5*max(band.df, band.rbw)
    flag = strtrim(flag + " CHECK: poor Doppler fit (rms " + round(h1.rms) + " Hz): weak/fragmented track or wrong candidate");
end
row = {kk, b, median(y)/1e6, (max(y) - min(y))/1e3, rate, N, N/max(nsw, 1), max(D.lvl(i)), median(D.snr(i)), ...
    median(D.w(i)), band.t0 + seconds(tau(1)), band.t0 + seconds(tau(end)), h1.label, h1.norad, h1.post, sec, ...
    h1.rms, h1.dt, h1.dfo, h1.el, string(ev), flag};
end

function flag = emiFlag(f, h)
B = [137e6 138e6 "137 MHz MetSat / ORBCOMM downlink"; 145.8e6 146e6 "2 m amateur-satellite"; ...
    400.15e6 401e6 "400 MHz space operation / MetSat"; 401e6 403e6 "401-403 MHz MetAids/EESS"; ...
    435e6 438e6 "70 cm amateur-satellite"; 1164e6 1300e6 "GNSS L5/E5/L2/G2/E6/B3"; 1525e6 1559e6 "MSS L-band downlink"; ...
    1559e6 1610e6 "GNSS L1/E1/B1/G1"; 1610e6 1626.5e6 "Iridium / MSS"; 1670e6 1710e6 "MetSat L-band (GK-2A)"; ...
    2200e6 2290e6 "S-band space operation downlink"; 2390e6 2450e6 "13 cm amateur-satellite (HamTV 2395)"; ...
    5830e6 5850e6 "5.8 GHz amateur-satellite (OrigamiSat-2)"];
lo = double(B(:, 1)); hi = double(B(:, 2)); k = find(f >= lo & f <= hi, 1);
flag = "";
if isempty(k), return, end
if any(double(B(k, 1)) == [1164e6 1559e6])
    flag = "EMI: visible signal in " + B(k, 3) + " (GNSS is below the noise floor) - interference/repeater";
elseif h.kind == "fixed"
    flag = "Possible interference in " + B(k, 3) + " (fixed frequency, no satellite Doppler)";
end
end

function n = demoCheck(T, truth, St)
n = strings(0, 1);
for k = 1:numel(truth)
    j = find(T.NORAD == truth(k).norad);
    if isempty(j)
        n(end+1) = sprintf("DEMO CHECK %s: no track attributed (single-sweep detection failed).", truth(k).name); %#ok<AGROW>
    else
        [~, m] = max(T.Confidence(j)); j = j(m);
        n(end+1) = sprintf("DEMO CHECK %s: track %d, confidence %.3f; fitted dt %.1f s (true %.1f), df %.0f Hz (true %.0f), rms %.0f Hz.", ...
            truth(k).name, T.Track(j), T.Confidence(j), T.TimeOffset_s(j), truth(k).dt, T.FreqOffset_Hz(j), truth(k).dfo, T.FitRMS_Hz(j)); %#ok<AGROW>
    end
    if ~isempty(St) && height(St) > 0
        i = find(St.NORAD == truth(k).norad & St.Detected);
        if ~isempty(i)
            [~, m] = max(St.PeakZ(i)); i = i(m);
            n(end+1) = sprintf("DEMO CHECK %s (Doppler stacking): z = %.1f, df %.0f Hz (true %.0f), dt %.0f s (true %.1f).", ...
                truth(k).name, St.PeakZ(i), St.FreqOffset_Hz(i), truth(k).dfo, St.TimeOffset_s(i), truth(k).dt); %#ok<AGROW>
        end
    end
end
n = n(:);
end

function printIdentify(T, info)
fprintf('\n%3s %4s %12s %7s %4s %6s %7s  %-58s %6s %7s %6s %7s\n', '#', 'band', 'f (MHz)', 'drift', 'N', 'SNR', 'BW', ...
    'Best hypothesis', 'conf*', 'rms Hz', 'dt s', 'df Hz');
for k = 1:height(T)
    fprintf('%3d %4d %12.6f %6.1fk %4d %6.1f %7s  %-58s %6.3f %7.0f %6.1f %7.0f\n', k, T.Band(k), T.Freq_MHz(k), T.Drift_kHz(k), ...
        T.NDet(k), T.SNR_dB(k), fmtHz(T.BW_Hz(k)), extractBefore(T.Best(k) + string(blanks(58)), 59), min(T.Confidence(k), 0.999), ...
        T.FitRMS_Hz(k), T.TimeOffset_s(k), T.FreqOffset_Hz(k));
    if strlength(T.Second(k)) > 0, fprintf('      2nd: %s\n', T.Second(k)); end
    if strlength(T.Flag(k)) > 0, fprintf('      FLAG: %s\n', T.Flag(k)); end
end
if height(T) == 0, fprintf('No emissions detected above the threshold.\n'); end
if height(T) > 0, fprintf(['* model posterior over the listed hypotheses (shown capped at 0.999); it is only as good as the ' ...
        'candidate list, the TLEs and the EIRP/width templates.\n']); end
St = info.stack;
if ~isempty(St) && height(St) > 0
    St = sortrows(St, 'PeakZ', 'descend');
    fprintf('\nWeak-signal search (Doppler stacking along each predicted pass; other emitters masked):\n');
    fprintf('%4s %-44s %12s %6s %7s %8s %6s  %s\n', 'band', 'Satellite / downlink', 'f0 (MHz)', 'sweeps', 'z', 'df Hz', 'dt s', 'result');
    for k = 1:min(height(St), 12)
        fprintf('%4d %-44s %12.6f %6d %7.1f %8.0f %6.0f  %s\n', St.Band(k), extractBefore(St.Satellite(k) + string(blanks(44)), 45), ...
            St.F0_MHz(k), St.NSweeps(k), St.PeakZ(k), St.FreqOffset_Hz(k), St.TimeOffset_s(k), St.Result(k));
    end
end
for n = unique(info.notes(:), 'stable')', fprintf('Note: %s\n', n); end
end

%% ----- weak-signal search: Doppler stacking ----------------------------------------------------
function [S, prof, notes] = stackSearch(band, D, C, opt, b, tr)
% For every catalog downlink that can fall into the band and whose satellite is up during the
% recording: shift each sweep by the predicted Doppler f0*d(t + dt) and average the noise-normalised
% excess power (bins of already-tracked emissions are masked). A satellite signal adds up along its
% own Doppler curve, noise and other emissions do not: ~5*log10(N) dB of detection gain for N sweeps
% with the signal on. dt (TLE timing) is scanned; the peak offset estimates the oscillator offset.
notes = strings(0, 1);
prof = struct('label', {}, 'u', {}, 'z', {}, 'band', {});
vn = {'Band', 'Satellite', 'NORAD', 'Downlink', 'F0_MHz', 'NSweeps', 'DopplerSpan_Hz', 'PeakZ', 'FreqOffset_Hz', ...
    'TimeOffset_s', 'Width_Hz', 'Detected', 'Result'};
S = cell2table(cell(0, numel(vn)), 'VariableNames', vn);
if band.nSweeps < 10 || isempty(D.nf) || numel(D.nf) ~= band.nSweeps, return, end
f = band.f(:)'; df = band.df; osc = opt.FreqTolPPM*1e-6*C.F0_Hz + opt.FreqTolHz;
reach = C.MaxDoppler_Hz + C.OBW_Hz/2 + 4*osc;
rows = find(ismember(C.Kind, ["sat", "constellation"]) & C.Class ~= "GNSS" & ~isfinite(C.BandLo_Hz) & ...
    C.F0_Hz + reach >= f(1) & C.F0_Hz - reach <= f(end) & C.OBW_Hz < 0.3*band.span);
if isempty(rows), return, end
Cr = C(rows, :);
grp = unique([Cr.Group(Cr.Kind == "constellation"); Cr.Group(Cr.Kind == "sat" & Cr.Group ~= "")]);
[E, ~] = getTLEs(opt, grp, Cr.NORAD(Cr.Kind == "sat" & Cr.Group ~= "satnogs"));
keepE = ismember([E.norad], Cr.NORAD(Cr.Kind == "sat"));
for r = find(Cr.Kind == "constellation")'
    keepE = keepE | ([E.group] == Cr.Group(r) & ~cellfun(@isempty, regexp(cellstr([E.name]), Cr.TleMatch(r), 'once')));
end
E = E(keepE);
if isempty(E), return, end
P = predictSats(E, band.t0 - seconds(120), band.t0 + seconds(max(band.tau) + 120), opt);
for k = 1:numel(P), P(k).tau = seconds(P(k).t - band.t0); end
% noise-normalised excess power; detections of tracks attributed to OTHER emitters are masked per
% candidate (fixed carriers / unknown classes for everybody), the candidate's own track is kept
X0 = 10.^(band.P/10)./10.^(D.nf(:)/10) - 1;
owner = zeros(numel(D.s), 1);                    % 0 = not in a track (never masked), -1 = fixed/unknown
for k = 1:numel(tr)
    h = tr(k).hyp(1);
    if ismember(h.kind, ["sat", "constellation"]) && h.post >= 0.5, owner(tr(k).idx) = h.norad; else, owner(tr(k).idx) = -1; end
end
% (satellite, downlink) pairs that are up during the recording
pairs = zeros(0, 2);
for r = 1:height(Cr)
    if Cr.Kind(r) == "sat", js = find([P.norad] == Cr.NORAD(r));
    else, js = find(arrayfun(@(p) p.group == Cr.Group(r) && ~isempty(regexp(p.name, Cr.TleMatch(r), 'once')), P)); end
    for j = js(:)'
        if nnz(interp1(P(j).tau, P(j).el, band.tau, 'linear', -90) >= opt.ElevationMask) >= 10, pairs(end+1, :) = [r j]; end %#ok<AGROW>
    end
end
% satellites already attributed to a track (annotated, still searched: the track may be a fragment)
known = nan(1, numel(tr));
for k = 1:numel(tr)
    h = tr(k).hyp(1);
    if ismember(h.kind, ["sat", "constellation"]) && h.post >= 0.5, known(k) = h.norad; end
end
if isempty(pairs), return, end
key = string(arrayfun(@(j) P(j).norad, pairs(:, 2))) + "|" + string(round(Cr.F0_Hz(pairs(:, 1))/1e3)) + "|" + Cr.Template(pairs(:, 1));
[~, iu] = unique(key, 'stable'); pairs = pairs(iu, :);
if size(pairs, 1) > 40
    notes(end+1) = sprintf("Band %d: %d stacking candidates, only the 40 with the most sweeps above the mask are searched.", b, size(pairs, 1));
    nUp = arrayfun(@(j) nnz(interp1(P(j).tau, P(j).el, band.tau, 'linear', -90) >= opt.ElevationMask), pairs(:, 2));
    [~, o] = sort(nUp, 'descend'); pairs = pairs(o(1:40), :);
end
rowsOut = {};
zThr = opt.StackZ + max(0, log10(9*size(pairs, 1)));        % look-elsewhere: candidates x dt trials
for q = 1:size(pairs, 1)
    r = pairs(q, 1); Pj = P(pairs(q, 2)); f0 = Cr.F0_Hz(r);
    X = X0;
    for k = find(owner ~= 0 & owner ~= Pj.norad)'
        X(D.s(k), abs(f - D.fc(k)) <= 0.75*D.w(k) + 3*band.rbw) = NaN;      % incl. skirts below -X dB
    end
    Fd = griddedInterpolant(Pj.tau, Pj.dfac, 'pchip', 'nearest');
    use = find(interp1(Pj.tau, Pj.el, band.tau, 'linear', -90) >= opt.ElevationMask);
    sigT = opt.TimeTol; if isempty(sigT), sigT = min(60, 5 + 3*abs(days(band.t0 - Pj.epoch))); end
    dts = unique(round(linspace(-2*sigT, 2*sigT, 9)));
    U = min(4*osc(rows(r)) + Cr.OBW_Hz(r)/2 + 3*band.rbw, band.span/2);
    u = (-ceil(U/df):ceil(U/df))*df;
    wsm = max(1, round(Cr.OBW_Hz(r)/df));                        % matched (box) smoothing over the OBW
    best = struct('z', -Inf, 'u', NaN, 'dt', NaN, 'w', NaN, 'zz', []);
    for dt = dts
        fs_ = f0*(1 + Fd(band.tau(use) + dt));
        acc = zeros(size(u)); cnt = acc;
        for i = 1:numel(use)
            v = interp1(f, X(use(i), :), u + fs_(i), 'linear', NaN);
            ok = isfinite(v); acc(ok) = acc(ok) + v(ok); cnt(ok) = cnt(ok) + 1;
        end
        pr = acc./max(cnt, 1); pr(cnt < 0.5*numel(use)) = NaN;
        if nnz(isfinite(pr)) < max(10, 1.5*wsm), continue, end
        m = median(pr, 'omitnan'); sd = 1.4826*median(abs(pr - m), 'omitnan');     % per-bin noise
        nIndep = min(16, max(1, wsm*df/max(df, band.rbw)));                       % independent bins averaged (capped)
        pr = movmean(pr, wsm, 'omitnan', 'Endpoints', 'fill');
        if ~(sd > 0), continue, end
        zz = (pr - m)/(sd/sqrt(nIndep)); [zmax, im] = max(zz);
        if zmax > best.z
            a = find(zz >= zmax/2); w = (numel(a))*df;
            best = struct('z', zmax, 'u', u(im), 'dt', dt, 'w', w, 'zz', zz);
        end
    end
    if ~isfinite(best.z), continue, end
    span = f0*range(Fd(band.tau(use) + best.dt));
    det = best.z >= zThr && span >= 5*band.rbw;
    if best.z < zThr, res = sprintf("not found (z < %.1f)", zThr);
    elseif ~det, res = "inconclusive: Doppler change < 5 RBW (a fixed carrier would stack too)";
    else
        res = "DETECTED along the predicted Doppler curve";
        % tracks that lie on this curve are probably fragments of the same weak signal
        fc = @(t) (f0 + best.u)*(1 + Fd(t + best.dt));
        on = [];
        for k = 1:numel(tr)
            i = tr(k).idx;
            if known(k) ~= Pj.norad && mean(abs(D.fc(i) - fc(D.tau(i))) <= 3*band.rbw + 3*df + 300) >= 0.6, on(end+1) = k; end %#ok<AGROW>
        end
        kk = find(known == Pj.norad);
        if ~isempty(kk), res = res + " (also attributed to track " + strjoin(string(kk), ", ") + ")"; end
        if ~isempty(on)
            notes(end+1) = sprintf("Band %d: Doppler stacking detects %s (z %.1f); track(s) %s of this band lie on its curve and are probably fragments of it.", ...
                b, Pj.name + " - " + Cr.Id(r), best.z, strjoin(string(on), ", ")); %#ok<AGROW>
        end
    end
    lab = Pj.name + " - " + Cr.Id(r);
    rowsOut(end+1, :) = {b, lab, Pj.norad, Cr.Id(r), f0/1e6, numel(use), span, best.z, best.u, best.dt, best.w, det, res}; %#ok<AGROW>
    prof(end+1) = struct('label', lab, 'u', u, 'z', best.zz, 'band', b); %#ok<AGROW>
end
if ~isempty(rowsOut), S = cell2table(rowsOut, 'VariableNames', vn); end
end

%% ----- identify plots ------------------------------------------------------------------------
function figs = plotIdentify(T, info, opt)
figs = gobjects(0);
for b = 1:numel(info.bands)
    band = info.bands(b); D = info.detections(b);
    fig = figure('Visible', opt.Visible, 'Name', sprintf('ff_satsig identify band %d', b), 'Position', [30 30 1300 760]);
    tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact');
    ax = nexttile(tl, [1 2]);
    fc0 = mean(band.f); [sc, unit] = fscale(band.span/2);
    imagesc(ax, (band.f - fc0)*sc, band.tau/60, band.P); axis(ax, 'xy'); colormap(ax, turbo);
    cb = colorbar(ax); cb.Label.String = 'dBm';
    lo = prctile(band.P(:), 5); clim(ax, [lo, max(lo + 10, prctile(band.P(:), 99.9))]);
    hold(ax, 'on');
    plot(ax, (D.fc - fc0)*sc, D.tau/60, 'w.', 'MarkerSize', 6);
    tk = find(T.Band == b);
    for k = tk(:)'
        tr = info.tracks(k); h1 = tr.hyp(1);
        if ~isempty(h1.curveTau)
            plot(ax, (h1.curveF - fc0)*sc, h1.curveTau/60, 'm--', 'LineWidth', 1.3);
        end
        ii = tr.idx; [~, im] = max(D.snr(ii));
        text(ax, (D.fc(ii(im)) - fc0)*sc, D.tau(ii(im))/60, sprintf('  #%d %s (%.0f%%)', k, shortLabel(h1.label), 100*h1.post), ...
            'Color', 'w', 'FontSize', 8, 'Interpreter', 'none', 'FontWeight', 'bold');
    end
    xlabel(ax, sprintf('Frequency - %s (%s)', fmtHz(fc0), unit)); ylabel(ax, sprintf('Minutes after %s UTC', string(band.t0, 'yyyy-MM-dd HH:mm:ss')));
    title(ax, sprintf('Band %d: waterfall, detections (white), fitted Doppler curve of the best hypothesis (magenta)', b));
    ax2 = nexttile(tl);
    k1 = tk(1:min(end, 6)); labs = strings(0, 1); vals = [];
    for k = k1(:)'
        Hk = info.tracks(k).hyp; n = min(3, numel(Hk));
        for m = 1:n, labs(end+1) = sprintf('#%d %s', k, shortLabel(Hk(m).label)); vals(end+1) = Hk(m).post; end %#ok<AGROW>
    end
    if ~isempty(vals)
        barh(ax2, flipud(vals(:))); yticks(ax2, 1:numel(vals)); yticklabels(ax2, flipud(labs(:)));
        ax2.TickLabelInterpreter = 'none'; ax2.FontSize = 7; xlim(ax2, [0 1]); grid(ax2, 'on');
    end
    title(ax2, 'Posterior of the top 3 hypotheses per track'); xlabel(ax2, 'Confidence');
    figs(end+1) = fig; %#ok<AGROW>
    if strlength(opt.OutDir) > 0, exportgraphics(fig, fullfile(opt.OutDir, sprintf('identify_band%d.png', b)), 'Resolution', 110); end
    closeSaved(fig, opt);
end
% Doppler fit detail of the best satellite tracks
sel = find(arrayfun(@(t) ~isempty(t.hyp(1).curveTau) && numel(t.idx) >= 4, info.tracks));
if ~isempty(sel)
    sel = sel(1:min(end, 3));
    fig = figure('Visible', opt.Visible, 'Name', 'ff_satsig Doppler fit', 'Position', [60 60 1100 300*numel(sel) + 80]);
    tl = tiledlayout(fig, numel(sel), 2, 'TileSpacing', 'compact');
    for k = sel(:)'
        tr = info.tracks(k); D = info.detections(tr.band); h1 = tr.hyp(1); i = tr.idx;
        ax = nexttile(tl);
        plot(ax, D.tau(i)/60, (D.fc(i) - h1.f0)/1e3, 'ko', 'MarkerSize', 4); hold(ax, 'on');
        plot(ax, h1.curveTau/60, (h1.curveF - h1.f0)/1e3, 'm-', 'LineWidth', 1.3);
        Pfd = (h1.curveF./(1 + h1.dfo/h1.f0)); %#ok<NASGU>
        grid(ax, 'on'); xlabel(ax, 'min'); ylabel(ax, sprintf('f - %s (kHz)', fmtHz(h1.f0)));
        title(ax, sprintf('#%d %s: dt %.1f s, df %.0f Hz', k, shortLabel(h1.label), h1.dt, h1.dfo), 'Interpreter', 'none', 'FontSize', 9);
        legend(ax, 'measured centroid', 'fitted prediction', 'Location', 'best');
        ax = nexttile(tl);
        fit = interp1(h1.curveTau, h1.curveF, D.tau(i));
        stem(ax, D.tau(i)/60, D.fc(i) - fit, 'filled', 'MarkerSize', 3); grid(ax, 'on');
        xlabel(ax, 'min'); ylabel(ax, 'residual (Hz)'); title(ax, sprintf('Residual RMS %.0f Hz', h1.rms), 'FontSize', 9);
    end
    figs(end+1) = fig;
    if strlength(opt.OutDir) > 0, exportgraphics(fig, fullfile(opt.OutDir, 'identify_dopplerfit.png'), 'Resolution', 110); end
    closeSaved(fig, opt);
end
% Doppler-stacked profiles (weak-signal search)
St = info.stack;
if ~isempty(St) && height(St) > 0
    [~, o] = sort(St.PeakZ, 'descend'); o = o(1:min(end, 6));
    fig = figure('Visible', opt.Visible, 'Name', 'ff_satsig Doppler stacking', 'Position', [80 80 1000 520]);
    ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
    for k = o(:)'
        pr = info.stackProfiles(k);
        plot(ax, pr.u/1e3, pr.z, 'LineWidth', 1.1, 'DisplayName', sprintf('%s (z %.1f)', shortLabel(pr.label), St.PeakZ(k)));
    end
    yline(ax, opt.StackZ, 'r--', 'HandleVisibility', 'off');
    xlabel(ax, 'Offset from the predicted Doppler-shifted frequency (kHz)'); ylabel(ax, 'z (robust)');
    title(ax, 'Weak-signal search: spectra averaged along each predicted pass (other emitters masked)');
    legend(ax, 'Location', 'northeast', 'Interpreter', 'none');
    figs(end+1) = fig;
    if strlength(opt.OutDir) > 0, exportgraphics(fig, fullfile(opt.OutDir, 'identify_stacking.png'), 'Resolution', 110); end
    closeSaved(fig, opt);
end
end

function closeSaved(fig, opt)
% A hidden figure is closed once saved: hidden figures must not pile up in a long session.
if strcmpi(opt.Visible, 'off') && strlength(opt.OutDir) > 0, close(fig); end
end

function dropHidden(h)
% Hidden figures the caller does not receive (no output requested) are deleted.
h = h(isgraphics(h));
if ~isempty(h), delete(h(strcmp(get(h, {'Visible'}), 'off'))); end
end

function s = shortLabel(l)
s = extractBefore(string(l) + string(blanks(46)), 47); s = strtrim(s);
end

%% ======================================================================================
%% LIVE ACQUISITION
%% ======================================================================================
function [R, notes] = acquireLive(src, opt, C)
notes = strings(0, 1);
if isa(src, 'FieldFox'), ff = src;
else
    fprintf('Connecting to %s:%d ...\n', string(src), opt.Port);
    try
        ff = FieldFox(string(src), opt.Port);
    catch err
        error('ff_satsig:connect', 'Cannot connect to a FieldFox at %s:%d (%s).', string(src), opt.Port, err.message);
    end
end
stamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
saveTo = string(opt.SaveTo); if strlength(saveTo) == 0, saveTo = fullfile(pwd, "ff_satsig_" + stamp); end
saveTo = erase(saveTo, ".mat");
if opt.Passive
    if ff.mode() ~= "SA", error('ff_satsig:mode', 'The FieldFox is in %s mode. Switch it to SA mode first (Passive changes nothing).', ff.mode()); end
    oldFormat = ff.query("FORM:DATA?"); restore = onCleanup(@() restoreFormat(ff, oldFormat)); %#ok<NASGU>
    s = ff.saSettings();
    cfg = struct('start', s.start, 'stop', s.stop, 'rbw', s.rbw, 'points', s.points, 'det', string(s.detector));
    notes(end+1) = "Passive: the instrument settings were not changed (only FORM:DATA during transfers, restored).";
else
    bands = opt.Bands;
    if isempty(bands), [bands, pn] = planBands(opt, C); notes = [notes; pn(:)]; end
    if size(bands, 2) ~= 2, error('ff_satsig:bands', '''Bands'' must be k x 2 [start stop] in Hz.'); end
    cfg = struct('start', {}, 'stop', {}, 'rbw', {}, 'points', {}, 'det', {});
    for k = 1:size(bands, 1)
        sp = bands(k, 2) - bands(k, 1);
        if sp <= 0, error('ff_satsig:bands', 'Band %d: stop must be above start.', k); end
        rb = opt.RBW; if isempty(rb), rb = rbwStep(max(sp/(opt.Points - 1), 10)); end
        cfg(k) = struct('start', bands(k, 1), 'stop', bands(k, 2), 'rbw', rb, 'points', opt.Points, 'det', upper(string(opt.Detector)));
    end
    if opt.Prompt && usejava('desktop') && feature('ShowFigureWindows')
        a = questdlg(sprintf(['ff_satsig will take control of the FieldFox (state is backed up and restored), ' ...
            'sweep %d band(s) for %g s. Antenna on Port 2. Continue?'], numel(cfg), opt.Duration), 'ff_satsig', 'Continue', 'Cancel', 'Continue');
        if ~strcmp(a, 'Continue'), error('ff_satsig:cancel', 'Cancelled by the user.'); end
    end
    ff.takeControl();
    rel = onCleanup(@() safeRelease(ff)); %#ok<NASGU>
end
nb = numel(cfg);
R = repmat(struct('f', [], 'P', [], 't', datetime.empty(0, 1), 'settings', struct(), 'rbw', NaN, 'detector', "", ...
    'idn', ff.IDN, 'source', "live " + string(ff.IP)), 1, nb);
for k = 1:nb, R(k).t.TimeZone = 'UTC'; end
t0 = tic; cyc = 0; stopReason = "duration reached"; lastY = []; nDup = 0;
try
    while toc(t0) < opt.Duration
        cyc = cyc + 1;
        for k = 1:nb
            if ~opt.Passive && (nb > 1 || cyc == 1)
                ff.saConfigure(Start = cfg(k).start, Stop = cfg(k).stop, RBW = cfg(k).rbw, Points = cfg(k).points, ...
                    Detector = cfg(k).det, Preamp = logical(opt.Preamp), AttenAuto = true, TraceType = "CLRW", Average = 1);
                if cyc == 1
                    cfg(k).rbw = ff.queryNum("SENS:BAND:RES?"); cfg(k).points = ff.queryNum("SENS:SWE:POIN?");
                    cfg(k).start = ff.queryNum("SENS:FREQ:STAR?"); cfg(k).stop = ff.queryNum("SENS:FREQ:STOP?");
                end
            end
            ta = datetime('now', 'TimeZone', 'UTC');
            if ~opt.Passive && upper(string(opt.TraceType)) == "MAXH"
                ff.write("TRAC1:TYPE CLRW"); ff.write("TRAC1:TYPE MAXH");
                th = tic; while toc(th) < opt.HoldTime, ff.sweep(); end
                tm = ta + (datetime('now', 'TimeZone', 'UTC') - ta)/2;
            else
                el = ff.sweep(); tm = ta + seconds(el/2);
            end
            y = ff.saTrace(1);
            if opt.Passive
                st = ff.queryNum("SENS:FREQ:STAR?"); so = ff.queryNum("SENS:FREQ:STOP?");
                if abs(st - cfg(k).start) > 1 || abs(so - cfg(k).stop) > 1 || numel(y) ~= cfg(k).points
                    stopReason = "the span was changed on the instrument (passive mode stops)"; break
                end
            end
            if numel(y) ~= cfg(k).points, error('ff_satsig:trace', 'Trace has %d points, expected %d.', numel(y), cfg(k).points); end
            if opt.Passive && isequal(y, lastY), nDup = nDup + 1; continue, end    % free run: same sweep read twice
            lastY = y;
            R(k).P(end+1, :) = y(:)'; R(k).t(end+1, 1) = tm;
        end
        if startsWith(stopReason, "the span"), break, end
        if mod(cyc, 5) == 0 || ~opt.Quiet && mod(cyc, 10) == 0
            saveLive(R, cfg, saveTo); fprintf('  %d cycles, %.0f s\n', cyc, toc(t0));
        end
        if opt.Passive, pause(max(opt.Interval, 0.5)); elseif opt.Interval > 0, pause(opt.Interval); end
    end
catch err
    if ~strcmp(err.identifier, 'FieldFox:timeout') && ~strcmp(err.identifier, 'ff_satsig:trace'), rethrow(err); end
    stopReason = "instrument error: " + string(err.message);
end
for k = 1:nb
    R(k).f = linspace(cfg(k).start, cfg(k).stop, cfg(k).points); R(k).rbw = cfg(k).rbw; R(k).detector = cfg(k).det;
    R(k).settings = struct('start', cfg(k).start, 'stop', cfg(k).stop, 'points', cfg(k).points, 'rbw', cfg(k).rbw, ...
        'detector', cfg(k).det, 'preamp', logical(opt.Preamp), 'center', (cfg(k).start + cfg(k).stop)/2, 'tool', "ff_satsig");
end
files = saveLive(R, cfg, saveTo);
notes(end+1) = sprintf("Live: %d cycles in %.0f s (%s). Saved %s", cyc, toc(t0), stopReason, strjoin(files, ", "));
if nDup > 0, notes(end+1) = sprintf("Passive: %d repeated reads of an unchanged trace were dropped.", nDup); end
keep = arrayfun(@(r) size(r.P, 1) > 0, R); R = R(keep);
if isempty(R), error('ff_satsig:nodata', 'No sweep was recorded (%s).', stopReason); end
end

function files = saveLive(R, cfg, base)
files = strings(0, 1);
for k = 1:numel(R)
    if isempty(R(k).P), continue, end
    f = linspace(cfg(k).start, cfg(k).stop, cfg(k).points); P = R(k).P; t = R(k).t; idn = R(k).idn; %#ok<NASGU>
    settings = struct('start', cfg(k).start, 'stop', cfg(k).stop, 'points', cfg(k).points, 'rbw', cfg(k).rbw, ...
        'detector', cfg(k).det, 'center', (cfg(k).start + cfg(k).stop)/2, 'tool', "ff_satsig"); %#ok<NASGU>
    fn = sprintf('%s_b%d.mat', base, k);
    save(fn, 'f', 'P', 't', 'settings', 'idn'); files(end+1) = fn; %#ok<AGROW>
end
end

function safeRelease(ff)
try
    ff.release();
catch err
    warning('ff_satsig:release', 'Could not restore the FieldFox: %s', err.message);
end
end

function restoreFormat(ff, oldFormat)
try
    if startsWith(oldFormat, "ASC")
        ff.write("FORM:DATA ASC,0");
    elseif contains(oldFormat, "64")
        ff.write("FORM:DATA REAL,64");
    else
        ff.write("FORM:DATA REAL,32");
    end
    ff.syncFormat();
catch err
    warning('ff_satsig:format', 'Could not restore FORM:DATA %s: %s', oldFormat, err.message);
end
end

function [bands, notes] = planBands(opt, C)
% Bands for the targets (or for the active catalog satellites visible now / during 'Duration').
notes = strings(0, 1);
if ~isempty(opt.Targets)
    q = string(opt.Targets); hit = false(height(C), 1);
    for s = q(:)'
        hit = hit | strcmpi(C.Id, s) | contains(C.Satellite, s, 'IgnoreCase', true);
    end
    rows = find(hit & C.Kind == "sat" & C.Active > 0);
    if isempty(rows), error('ff_satsig:targets', 'No active satellite downlink matches ''Targets''.'); end
else
    o = opt; o.Plot = false; o.Quiet = true;
    [T, Pass] = cmdOverhead('Site', o.Site, 'Offline', o.Offline, 'Constellations', false, 'Plot', false, ...
        'Quiet', true, 'Hours', max(opt.Duration/3600, 0.05), 'ElevationMask', opt.ElevationMask);
    vis = unique([T.NORAD; Pass.NORAD]);
    rows = find(ismember(C.NORAD, vis) & C.Kind == "sat" & C.Active >= 0.5 & C.F0_Hz <= 6.5e9);
    if isempty(rows)
        [~, Pn] = cmdOverhead('Site', o.Site, 'Offline', o.Offline, 'Constellations', false, 'Plot', false, 'Quiet', true, 'Hours', 12);
        msg = "No active catalog satellite is up now. Next passes:";
        for k = 1:min(height(Pn), 6), msg = msg + newline + sprintf('  %s  %s UTC  max el %.0f', Pn.Satellite(k), string(Pn.AOS_UTC(k), 'MM-dd HH:mm'), Pn.MaxEl_deg(k)); end
        error('ff_satsig:nothingUp', '%s\nGive ''Bands'' or ''Targets'' to record anyway.', msg);
    end
end
lo = zeros(0, 1); hi = lo;
for r = rows'
    osc = opt.FreqTolPPM*1e-6*C.F0_Hz(r) + opt.FreqTolHz;
    hw = 1.15*(C.MaxDoppler_Hz(r) + C.OBW_Hz(r)/2) + 2*osc;
    lo(end+1, 1) = C.F0_Hz(r) - hw; hi(end+1, 1) = C.F0_Hz(r) + hw; %#ok<AGROW>
    notes(end+1) = "Planned band for " + C.Id(r) + ": " + fmtHz(lo(end)) + " - " + fmtHz(hi(end)); %#ok<AGROW>
end
[lo, o] = sort(lo); hi = hi(o); bands = [lo(1) hi(1)];
for k = 2:numel(lo)
    if lo(k) <= bands(end, 2), bands(end, 2) = max(bands(end, 2), hi(k)); else, bands(end+1, :) = [lo(k) hi(k)]; end %#ok<AGROW>
end
if size(bands, 1) > 4, notes(end+1) = "More than 4 bands: only the first 4 are swept."; bands = bands(1:4, :); end
end

%% ======================================================================================
%% DEMO: synthetic FieldFox recordings of real (SGP4) passes
%% ======================================================================================
function [R, truth, notes] = simulateDemo(which, opt, tpl)
rng(opt.Seed);
which = lower(string(which)); if which == "true", which = "origami"; end
list = strings(0, 1);
if any(which == ["origami", "both"]), list(end+1) = "origami"; end
if which == "weak", list(end+1) = "weak"; end
if any(which == ["knacksat", "both"]), list(end+1) = "knacksat"; end
if isempty(list), error('ff_satsig:demo', '''Demo'' must be "origami", "knacksat", "both" or "weak".'); end
notes = "DEMO: synthetic FieldFox traces (template spectra + Gaussian RBW + peak-detected noise) on real SGP4 geometry from the shipped/cached TLE. Upper-bound accuracy.";
R = []; truth = struct('name', {}, 'norad', {}, 'dt', {}, 'dfo', {});
for w = list(:)'
    if w == "origami" || w == "weak"
        c = struct('name', "OrigamiSat-2", 'norad', 68795, 'f0', 437.505e6, 'dt', 20, 'dfo', 1350, 'drift', 0.15, ...
            'span', 50e3, 'rbw', 100, 'interval', 4, 'eirp', 30, 'rxg', 6, 'N0', -168, 'fix', [437.4885e6 -128 1.0]);
    else
        c = struct('name', "KNACKSAT-2", 'norad', 67683, 'f0', 400.630e6, 'dt', -12, 'dfo', -900, 'drift', -0.1, ...
            'span', 50e3, 'rbw', 100, 'interval', 10, 'eirp', 30, 'rxg', 6, 'N0', -168, 'fix', [400.648e6 -133 0.4]);
    end
    if w == "weak", c.eirp = c.eirp - 17; c.fix(2) = c.fix(2) - 12; end      % below the single-sweep threshold
    o = opt; o.Offline = true;
    [E, ~] = getTLEs(o, strings(0, 1), c.norad);
    if isempty(E), error('ff_satsig:demo', 'No TLE for NORAD %d in the shipped sample.', c.norad); end
    [aos, los] = demoPass(E(1), opt);
    t = (aos - seconds(60):seconds(c.interval):los + seconds(60))';
    Pp = predictSats(E(1), aos - seconds(500), los + seconds(500), opt, 1);
    tq = seconds(t - (aos - seconds(500))) + c.dt;
    d = interp1(Pp.tau, Pp.dfac, tq, 'pchip'); el = interp1(Pp.tau, Pp.el, tq); rg = interp1(Pp.tau, Pp.rng, tq);
    tau = seconds(t - t(1));
    fT = (c.f0 + c.dfo)*(1 + d) + c.drift*tau;
    f = linspace(c.f0 - c.span/2, c.f0 + c.span/2, 401); df = f(2) - f(1);
    Prx = c.eirp - fspl(rg, physconst('LightSpeed')/c.f0) - 3 + c.rxg + 2*randn(size(rg));
    P = zeros(numel(t), numel(f)); phase = 30*rand; ph60 = 60*rand;
    [uC, DC] = tplDisplay(tpl.cw16, c.rbw);
    for s = 1:numel(t)
        sig = zeros(size(f));
        on = el(s) > 0;
        if w ~= "knacksat"
            ph = mod(tau(s) + phase, 30);
            if on && ph < 18 && rand < 0.6, sig = sig + 10^(Prx(s)/10)*binMax(uC, DC, f - fT(s), df); end
            if on && ph >= 24 && ph < 25.5
                [u, D1] = tplDisplay(tpl.afsk1k2, c.rbw); sig = sig + 10^(Prx(s)/10)*binMax(u, D1, f - fT(s), df);
            end
        else
            tb = ph60 + 60*floor((tau(s) + 5 - ph60)/60);       % a beacon inside the 10 s max-hold window?
            if abs(tb - tau(s)) <= 5
                tqb = seconds(t(1) - (aos - seconds(500))) + tb + c.dt;
                if interp1(Pp.tau, Pp.el, tqb) > 0
                    fb = (c.f0 + c.dfo)*(1 + interp1(Pp.tau, Pp.dfac, tqb, 'pchip')) + c.drift*tb;
                    pb = c.eirp - fspl(interp1(Pp.tau, Pp.rng, tqb), physconst('LightSpeed')/c.f0) - 3 + c.rxg + 2*randn;
                    [u, D1] = tplDisplay(tpl.g3ruh9k6, c.rbw); sig = sig + 10^(pb/10)*binMax(u, D1, f - fb, df);
                end
            end
        end
        if rand < c.fix(3)
            sig = sig + 10^((c.fix(2) + 0.5*randn)/10)*binMax(uC, DC, f - (c.fix(1) + 0.05*tau(s)), df);
        end
        nl = 10^((c.N0 + 10*log10(1.0645*c.rbw))/10);
        K = 1 + (w == "knacksat")*9;                            % max hold over ~10 sweeps
        noise = nl*max(-log(rand(numel(f), 2*K)), [], 2)';
        P(s, :) = 10*log10(sig + noise);
    end
    settings = struct('start', f(1), 'stop', f(end), 'points', numel(f), 'rbw', c.rbw, 'detector', "POS", ...
        'preamp', true, 'center', c.f0, 'tool', "ff_satsig demo", 'traceType', ternary(w == "knacksat", "MAXH 10 s", "CLRW"));
    r = struct('f', f, 'P', P, 't', t, 'settings', settings, 'rbw', c.rbw, 'detector', "POS", ...
        'idn', "SIMULATED (ff_satsig demo)", 'source', "demo " + c.name);
    R = [R, r]; %#ok<AGROW>
    truth(end+1) = struct('name', c.name, 'norad', c.norad, 'dt', c.dt, 'dfo', c.dfo); %#ok<AGROW>
    notes(end+1) = sprintf("DEMO %s pass %s - %s UTC (TLE epoch %s); truth: TLE timing error %+g s, oscillator %+g Hz, drift %g Hz/s; fixed carrier at %s.", ...
        c.name, string(aos, 'yyyy-MM-dd HH:mm:ss'), string(los, 'HH:mm:ss'), string(E(1).epoch, 'yyyy-MM-dd HH:mm'), c.dt, c.dfo, c.drift, fmtHz(c.fix(1))); %#ok<AGROW>
    if w == "knacksat", notes(end+1) = "DEMO KNACKSAT-2: 60 s beacon seen through 10 s MAXH windows (idealised: the sweep is assumed to catch the packet)."; end %#ok<AGROW>
end
end

function v = binMax(u, D1, x, df)
% Peak-detector value per trace bin: max of the displayed response over the bin.
v = interp1(u, D1, x, 'linear', 0);
for o = [-0.4 -0.2 0.2 0.4]*df
    v = max(v, interp1(u, D1, x + o, 'linear', 0));
end
end

function [aos, los] = demoPass(E, opt)
% First pass after the TLE epoch with max elevation >= 30 deg (else the highest in 2 days).
t0 = dateshift(E.epoch, 'start', 'minute');
sc = satelliteScenario(t0, t0 + days(2), 30);
file = writeTLEfile(E); c = onCleanup(@() delete(file));
sat = satellite(sc, char(file), 'OrbitPropagator', 'sgp4');
gs = groundStation(sc, opt.Site(1), opt.Site(2), 'Altitude', opt.Site(3));
I = accessIntervals(access(sat, gs));
[p, ~, tt] = states(sat, 'CoordinateFrame', 'ecef');
[~, el] = ecef2aer(p(1, :), p(2, :), p(3, :), opt.Site(1), opt.Site(2), opt.Site(3), wgs84Ellipsoid);
if isempty(I) || height(I) == 0, error('ff_satsig:demo', 'No pass of %s over the site within 2 days of the TLE epoch.', E.name); end
mx = zeros(height(I), 1);
for j = 1:height(I), mx(j) = max([el(tt >= I.StartTime(j) & tt <= I.EndTime(j)), -90]); end
j = find(mx >= 30, 1); if isempty(j), [~, j] = max(mx); end
aos = I.StartTime(j); los = I.EndTime(j);
if isempty(aos.TimeZone), aos.TimeZone = 'UTC'; los.TimeZone = 'UTC'; end
end

%% ======================================================================================
%% small helpers
%% ======================================================================================
function r = rbwStep(x, up)
% FieldFox RBW steps 1-3-10: largest step <= x (or smallest >= x with up = true), 10 Hz - 5 MHz.
if nargin < 2, up = false; end
s = [1 3] .* 10.^(0:7)'; s = sort(s(:)); s = s(s >= 10 & s <= 5e6);
if up, k = find(s >= x, 1); if isempty(k), k = numel(s); end
else, k = find(s <= x, 1, 'last'); if isempty(k), k = 1; end, end
r = s(k);
end

function v = niceUp(x)
m = [1 2 3 5 10]; e = 10^floor(log10(x)); v = e*m(find(m*e >= x, 1));
end

function [sc, unit] = tscale(T)
if T >= 2, sc = 1; unit = "s"; else, sc = 1e3; unit = "ms"; end
end

function [sc, unit] = fscale(x)
if x >= 1e6, sc = 1e-6; unit = "MHz"; elseif x >= 1e3, sc = 1e-3; unit = "kHz"; else, sc = 1; unit = "Hz"; end
end

function s = fmtHz(f)
if ~isfinite(f), s = "-"; return, end
a = abs(f);
if a >= 1e9, s = sprintf('%.6g GHz', f/1e9); elseif a >= 1e6, s = sprintf('%.7g MHz', f/1e6);
elseif a >= 1e3, s = sprintf('%.4g kHz', f/1e3); else, s = sprintf('%.3g Hz', f); end
s = string(s);
end

function s = numstr(v)
if isempty(v) || ~isfinite(v), s = "-"; else, s = string(sprintf('%.6g', v)); end
end

function o = ternary(c, a, b)
if c, o = a; else, o = b; end
end
