function R = ff_vna(ip, varargin)
%FF_VNA  S-parameters with a Keysight FieldFox, analysed with RF / SI / Antenna Toolbox.
%
%   R = ff_vna                              S11 of port 1, FieldFox at 10.90.36.200
%   R = ff_vna(ip, Name, Value, ...)        e.g. ff_vna("10.90.36.200", 'Ports', 2)
%   R = ff_vna(ff, ...)                     use an open FieldFox object
%   R = ff_vna([], 'File', "dut.s2p")       analyse a saved file (no instrument)
%   R = ff_vna([], 'Demo', "antenna")       synthetic data: "antenna" | "filter"
%
% Measurement (instrument):
%   'Method' "na"        vector S-parameters in NA mode (options 210/211):
%                        port 1 only (S11), or both ports (S11 S21 S12 S22).
%   'Method' "tracking"  scalar |S21| with the spectrum analyzer + tracking
%                        generator (source = port 1, SA input = port 2). High
%                        dynamic range for filters, but no phase and no S11.
%   QuickCal (option 112, no cal kit needed): 'Cal' "quick1" (port 1) or
%   "quick2" (full 2-port, needs a thru). You are told what to connect.
%   The instrument state is saved first and restored at the end, also after an
%   error, so the FieldFox returns to the spectrum analyzer as you left it.
%   That restore also brings back the old calibration, so a QuickCal made here
%   is used for this run only. To reuse one cal for several DUTs, hold control:
%     ff = FieldFox("10.90.36.200"); ff.takeControl();
%     A = ff_vna(ff, 'Cal', "quick1");  B = ff_vna(ff);   % same Start/Stop/Points/IFBW/Power
%     ff.release();
%
% Analysis (figures + R.metrics):
%   rfplot of all parameters; Smith chart (smithplot) of S11/S22; return loss
%   and VSWR (vswr) with the best match and the -10 dB bandwidth; 2-port:
%   insertion loss, 3 dB bandwidth (powerbw), phase, group delay (groupdelay);
%   S-parameter quality: passivity/causality/reciprocity (tdr.sParamQuality);
%   TDR impedance vs distance (Signal Integrity tdr) and distance-to-fault
%   (windowed IFFT of S11); rational model (rational, TendsToZero=false) and
%   its fit error; L/Pi/Tee matching network (matchingnetwork) at 'MatchFrequency'
%   with the predicted matched return loss; fixture removal (deembedsparams);
%   Touchstone export (rfwrite) and an rfbudget-ready nport element (R.nport);
%   comparison with a simulated antenna (Antenna Toolbox design, sparameters,
%   impedance, resonantFrequency) overlaid on rfplot and the Smith chart.
%
% !! RF OUTPUT: NA mode and the tracking generator send RF out of port 1 (NA
% with 2 ports: also port 2), default -10 dBm. Disconnect anything that must
% not be driven (LISN, limiter, receiver). Port damage level +27 dBm, 50 VDC.
%
% Name-Value options (defaults):
%   'Method'          "na"     "na" | "tracking"
%   'Ports'           1        1 = S11, 2 = S11 S21 S12 S22 (NA, option 211)
%   'Start', 'Stop'   1e6, 1e9 Hz (0.5-1.5 x DesignFrequency if 'Antenna' is set)
%   'Points'          401      2 ... 10001
%   'IFBW'            1e3      NA IF bandwidth, Hz (10 ... 100e3; lower = less noise, slower)
%   'RBW'             []       tracking: SA resolution bandwidth, Hz ([] = auto)
%   'Power'           -10      source power, dBm (NA -45...+3, tracking -45...0)
%   'Average'         1        sweeps averaged
%   'Cal'             "none"   "none" | "quick1" | "quick2" (QuickCal)
%   'Connector'       "Type N -M-,50"  DUT connector(s) for QuickCal (1 or 2 strings)
%   'Normalize'       true     tracking: measure a thru first (needs Prompt)
%   'File'            ""       .s1p/.s2p/.sNp, or a .mat saved by ff_vna (offline)
%   'Demo'            ""       "antenna" | "filter" synthetic data (offline)
%   'Fixture'         ""       .s2p of a cable/fixture on port 1 to remove (file
%                              port 1 at the instrument, port 2 at the DUT)
%   'Fixture2'        ""       same for port 2 (default: mirror of 'Fixture')
%   'VelocityFactor'  0.66     of the cable, for TDR/DTF distances (solid PE coax)
%   'MaxDistance'     []       m shown in TDR/DTF ([] = alias-free range)
%   'MatchFrequency'  []       Hz: design a matching network for S11 there
%   'MatchTopology'   "L"      "L" | "Pi" | "Tee"
%   'Antenna'         ""       "dipole" | "monopole" | "patchMicrostrip" | antenna object
%   'DesignFrequency' []       Hz, antenna design frequency ([] = measured best match)
%   'AntennaPoints'   25       simulated frequencies (dipole ~20 s, patch ~90 s)
%   'Analyze'         true     false = measure/load only, no analysis or figures
%   'SaveTouchstone'  []       true (default when measuring) | false | file name
%                              (a fixture-corrected measurement also keeps *_raw.sNp)
%   'SaveFigures'     false    PNG of every figure in OutDir
%   'OutDir'          pwd
%   'Visible'         'on'     'off' = figures hidden (scripts); they are in R.figures:
%                              delete(R.figures) when done (no output: deleted at the end)
%   'Prompt'          true     false = never wait for the keyboard (QuickCal and the
%                              tracking thru step are skipped); or a handle
%                              @(msg) ... called with each instruction, whose text
%                              reply is used: @(msg) "s" = unattended quick1 with
%                              the port open and the optional load skipped
%   'TcpPort'         5025
%
% Result R: S (sparameters), f (Hz), metrics (table Metric/Value/Unit), files
% (saved Touchstone files), cal, settings, idn, fit, quality, tdr, dtf, match,
% antenna, nport, budget, figures, otherFiles, notes. "tracking" also saves a
% .mat recording (f, P, t, settings, idn, S21_dB) like fieldfox_live.
% Read R.notes: they say when a result means nothing (e.g. open port, no
% transmission). Measured data are saved before the analysis, and an analysis
% step that fails only adds a note.
%
% Accuracy notes:
%  - Without a user cal the FieldFox applies its factory CalReady correction at
%    the test ports: a cable you add is part of the result (QuickCal at the
%    cable end, or remove it with 'Fixture').
%  - QuickCal is best with Type-N or 7/16 connectors, not recommended for female
%    SMA/3.5 mm. A cal is valid only for its frequencies, points, IFBW and power.
%  - "tracking" is magnitude only. Without a thru (Prompt=false) S21 is the
%    trace minus the source power, so cable losses are included.
%  - TDR comes from a rational fit of band-limited data: it needs a low start
%    frequency (<= Stop/100). Resolution ~ c*vf/(2*span). DTF ignores cable loss.
%  - The antenna simulation is an ideal antenna in free space fed with 50 ohm (no
%    balun, cable, mast or ground): a few % resonance difference is normal.
%  - Matching networks use ideal lossless L and C.

if nargin < 1, ip = []; end
opt = struct('Method', "na", 'Ports', 1, 'Start', [], 'Stop', [], 'Points', 401, ...
    'IFBW', 1e3, 'RBW', [], 'Power', -10, 'Average', 1, 'Cal', "none", ...
    'Connector', "Type N -M-,50", 'Normalize', true, 'File', "", 'Demo', "", ...
    'Fixture', "", 'Fixture2', "", 'VelocityFactor', 0.66, 'MaxDistance', [], ...
    'MatchFrequency', [], 'MatchTopology', "L", 'Antenna', "", 'DesignFrequency', [], ...
    'AntennaPoints', 25, 'Analyze', true, 'SaveTouchstone', [], 'SaveFigures', false, ...
    'OutDir', pwd, 'Visible', 'on', 'Prompt', true, 'TcpPort', 5025);
names = fieldnames(opt);
if (ischar(ip) || isStringScalar(ip)) && any(strcmpi(ip, names)) && mod(numel(varargin), 2) == 1
    varargin = [{ip} varargin]; ip = [];          % ff_vna('File', "x.s2p") without the ip
end
if mod(numel(varargin), 2) ~= 0, error('ff_vna:option', 'Options must be Name, Value pairs.'); end
for k = 1:2:numel(varargin)
    i = find(strcmpi(names, varargin{k}), 1);
    if isempty(i)
        error('ff_vna:option', 'Unknown option "%s". Valid options: %s.', varargin{k}, strjoin(names, ', '));
    end
    opt.(names{i}) = varargin{k + 1};
end
opt.Method = lower(string(opt.Method)); opt.Cal = lower(string(opt.Cal));
opt.Demo = lower(string(opt.Demo)); opt.File = string(opt.File);
if ~any(opt.Method == ["na" "tracking"]), error('ff_vna:option', 'Method must be "na" or "tracking".'); end
if ~any(opt.Ports == [1 2]), error('ff_vna:option', 'Ports must be 1 or 2.'); end
if ~any(opt.Cal == ["none" "quick1" "quick2"]), error('ff_vna:option', 'Cal must be "none", "quick1" or "quick2".'); end
if opt.Cal == "quick2", opt.Ports = 2; end
if opt.Points < 2 || opt.Points > 10001, error('ff_vna:option', 'Points must be 2 ... 10001.'); end
if ~exist(opt.OutDir, 'dir'), mkdir(opt.OutDir); end

R = struct('S', [], 'f', [], 'metrics', newMetrics(), 'files', strings(0, 1), ...
    'method', opt.Method, 'source', "", 'idn', "", 'settings', struct(), ...
    'cal', struct('text', "unknown (offline data)"), 'scalar', false, 'fit', [], ...
    'quality', [], 'tdr', [], 'dtf', [], 'match', [], 'antenna', [], 'nport', [], ...
    'budget', [], 'figures', gobjects(0), 'otherFiles', strings(0, 1), 'notes', strings(0, 1));

%% Get the data
live = strlength(opt.File) == 0 && strlength(opt.Demo) == 0;
if strlength(opt.File) > 0
    [S, R] = loadFile(opt.File, R);
elseif strlength(opt.Demo) > 0
    [S, R] = demoData(opt, R);
else
    [S, R] = measure(ip, opt, R);
end

%% Fixture removal, export (before the analysis, so a measurement is never lost), analysis
Sraw = S;
if strlength(string(opt.Fixture)) > 0
    if R.scalar
        warning('ff_vna:fixture', 'Fixture removal needs vector (NA) data: ignored for scalar S21.');
    else
        try
            [S, R] = removeFixture(S, opt, R);
        catch err
            warning('ff_vna:fixture', 'Fixture not removed: %s', err.message);
            R.notes(end + 1, 1) = "Fixture NOT removed: " + err.message;
        end
    end
end
R.S = S; R.f = S.Frequencies;
saveTs = opt.SaveTouchstone;
if isempty(saveTs), saveTs = live; end
named = isstring(saveTs) || ischar(saveTs);
if ~R.scalar && (named || logical(saveTs))
    ext = sprintf('.s%dp', S.NumPorts);
    if named
        [fd, fn, fe] = fileparts(string(saveTs));
        if strlength(fd) == 0, fd = string(opt.OutDir); end
        if ~strcmpi(fe, ext), fe = ext; end
        file = fullfile(fd, fn + fe);
    else
        file = fullfile(opt.OutDir, "ff_vna_" + stamp() + ext);
    end
    if exist(file, 'file'), fprintf('Replacing %s\n', file); end
    rfwrite(S, file, 'ForceOverwrite', true);       % (no "replace?" dialog)
    R.files(end + 1, 1) = file;
    fprintf('Touchstone saved: %s\n', file);
    if live && ~isequal(Sraw.Parameters, S.Parameters)       % keep the raw measurement too
        file = replace(file, ext, "_raw" + ext);
        rfwrite(Sraw, file, 'ForceOverwrite', true);
        R.files(end + 1, 1) = file;
    end
end
hidden = any(strcmpi(string(opt.Visible), ["off" "false" "0"]));
if opt.Analyze
    figs0 = findall(groot, 'Type', 'figure');
    try
        R = analyse(S, R, opt);
    catch err
        where = "";
        if ~isempty(err.stack), where = sprintf(" (in %s, line %d)", err.stack(1).name, err.stack(1).line); end
        warning('ff_vna:analysis', 'Analysis stopped: %s%s. The data in R.S are complete.', err.message, where);
        R.notes(end + 1, 1) = "Analysis stopped: " + err.message + where;
        if hidden                                   % its hidden figures are not in R: delete them
            figs = findall(groot, 'Type', 'figure');
            delete(figs(~ismember(figs, figs0)));
        end
    end
end
if opt.SaveFigures
    st = stamp();
    for k = 1:numel(R.figures)
        if ~isvalid(R.figures(k)), continue, end
        file = fullfile(opt.OutDir, "ff_vna_" + st + "_" + k + ".png");
        exportgraphics(R.figures(k), file, 'Resolution', 150);
        R.otherFiles(end + 1, 1) = file;
    end
end
printSummary(R);
if hidden && nargout == 0                           % hidden figures nobody holds (except ans): delete
    delete(R.figures(isvalid(R.figures)));
end
end

%% ===== measurement ========================================================

function [S, R] = measure(ip, opt, R)
if isa(ip, 'FieldFox')
    ff = ip;
else
    if isempty(ip) || strlength(string(ip)) == 0, ip = "10.90.36.200"; end
    fprintf('Connecting to %s:%d ...\n', ip, opt.TcpPort);
    ff = FieldFox(string(ip), opt.TcpPort);
end
fprintf('Connected: %s\n', ff.IDN);
need = "210";
if opt.Method == "tracking", need = ["210" "233"]; elseif opt.Ports == 2, need = ["210" "211"]; end
if opt.Cal ~= "none", need(end + 1) = "112"; end
missing = need(~arrayfun(@(o) ff.hasOption(o), need));
if ~isempty(missing)
    error('ff_vna:option', 'This measurement needs FieldFox option(s) %s.', strjoin(missing, ', '));
end
R.idn = ff.IDN; R.source = ff.IP;
mine = ~ff.InControl;                  % a caller that holds control keeps it (and a QuickCal)
if mine, ff.takeControl(); end
restore = onCleanup(@() releaseSafe(ff, mine));
if opt.Method == "tracking"
    [S, R] = measureTracking(ff, opt, R);
else
    [S, R] = measureNA(ff, opt, R);
end
clear restore            % instrument back to how it was, now
end

function [S, R] = measureNA(ff, opt, R)
pw = min(max(opt.Power, -45), 3);
if pw ~= opt.Power, warning('ff_vna:power', 'NA source power limited to %g dBm (range -45...+3).', pw); end
ports = "port 1"; if opt.Ports == 2, ports = "ports 1 AND 2"; end
rfWarning(opt, sprintf('NA mode sends %.1f dBm RF out of %s for this measurement.', pw, ports));

ff.setMode("NA");
[f1, f2] = sweepRange(opt);
w = @(varargin) ff.write(sprintf(varargin{:}));
w('SENS:FREQ:STAR %.10g', f1); w('SENS:FREQ:STOP %.10g', f2);
w('SENS:SWE:POIN %d', opt.Points); w('SENS:BWID %.10g', opt.IFBW);
w('SOUR:POW %.1f', pw); w('SENS:AVER:COUN %d', max(1, round(opt.Average)));
params = ["S11" "S21" "S12" "S22"];
if opt.Ports == 1, params = "S11"; end
w('CALC:PAR:COUN %d', numel(params));
for k = 1:numel(params), w('CALC:PAR%d:DEF %s', k, params(k)); end
ff.opc(30);
warnErrors(ff, 'NA setup');

st.start = readNum(ff, "SENS:FREQ:STAR?", f1); st.stop = readNum(ff, "SENS:FREQ:STOP?", f2);
st.points = readNum(ff, "SENS:SWE:POIN?", opt.Points); st.ifbw = readNum(ff, "SENS:BWID?", opt.IFBW);
st.power = readNum(ff, "SOUR:POW?", pw); st.average = readNum(ff, "SENS:AVER:COUN?", max(1, round(opt.Average)));
st.params = params; st.ports = opt.Ports; st.method = "na"; st.time = datetime('now');
fprintf('NA: %s - %s, %d points, IFBW %s, %.1f dBm, %s\n', fmtHz(st.start), fmtHz(st.stop), ...
    st.points, fmtHz(st.ifbw), st.power, strjoin(params, ' '));

calDone = "";
if opt.Cal ~= "none", calDone = quickCal(ff, opt); end
R.cal = calState(ff, calDone);
fprintf('Correction: %s\n', R.cal.text);

% Single sweeps; with averaging the FieldFox needs one INIT per average in hold mode
tmo = 30 + 3*st.points/st.ifbw*min(2, numel(params));
if st.average > 1
    ff.write("SENS:AVER:CLE");                 % UNVERIFIED on this unit (PG: CAT/NA/VVM)
    warnErrors(ff, 'SENS:AVER:CLE (restart averaging)');
end
tic;
for k = 1:max(1, st.average), opcQuery(ff, "INIT:IMM", tmo); end
fprintf('Swept %d x in %.1f s\n', max(1, st.average), toc);

f = ff.naFrequencies(); f = f(:);
if numel(f) ~= st.points
    error('ff_vna:data', 'FREQ:DATA? returned %d values, %d points expected (timeout or instrument busy).', numel(f), st.points);
end
D = complex(zeros(opt.Ports, opt.Ports, numel(f)));
ij = [1 1; 2 1; 1 2; 2 2];
for k = 1:numel(params)
    s = ff.naSData(k);
    if numel(s) ~= numel(f)
        error('ff_vna:data', '%s: %d values received, %d expected (timeout or instrument busy).', params(k), numel(s), numel(f));
    end
    D(ij(k, 1), ij(k, 2), :) = s;
end
warnErrors(ff, 'NA data read');
S = sparameters(D, f, 50);
R.settings = st; R.method = "na";
end

function [S, R] = measureTracking(ff, opt, R)
pw = min(max(opt.Power, -45), 0);
if pw ~= opt.Power, warning('ff_vna:power', 'Tracking generator power limited to %g dBm (range -45...0).', pw); end
rfWarning(opt, sprintf(['The tracking generator sends %.1f dBm RF out of port 1 ' ...
    '(port 2 = spectrum analyzer input).'], pw));
[f1, f2] = sweepRange(opt);
args = {'Start', f1, 'Stop', f2, 'Points', opt.Points, 'Preamp', false, 'AttenAuto', true, ...
    'RefLevel', pw + 5, 'TraceType', "CLRW", 'Average', 1};
if isempty(opt.RBW), args = [args {'RBWAuto', true}]; else, args = [args {'RBW', opt.RBW}]; end
ff.saConfigure(args{:});
rfOff = onCleanup(@() sourceOff(ff));
ff.write("SOUR:ENAB 1"); ff.write("SOUR:TRAC TSWP"); ff.write(sprintf("SOUR:POW %.1f", pw));
ff.opc(30);
e = ff.errors();
if ~isempty(e), error('ff_vna:tracking', 'Tracking generator setup failed: %s', strjoin(e, '; ')); end

s = ff.saSettings();
f = s.f;
P = zeros(0, numel(f)); t = NaT(0, 1); rows = strings(0, 1);
yNorm = []; thru = false;
if opt.Normalize && prompting(opt)
    ask(opt, 'TRACKING step 1: connect a THRU (cable/adapter, no DUT) from port 1 to port 2, then press Enter.');
    yThru = acquire(ff, opt.Average, numel(f));
    P(end + 1, :) = yThru; t(end + 1, 1) = datetime('now'); rows(end + 1, 1) = "thru";
    thru = true;
    % The FieldFox's own normalisation is only a cross-check (what TRAC:DATA? returns
    % with SOUR:NORM ON is not verified on this unit); S21 comes from DUT - thru.
    ff.write("SOUR:NORM 1"); ff.opc(30);
    e = ff.errors(); normOK = isempty(e);
    if ~normOK, R.notes(end + 1, 1) = "SOUR:NORM not accepted (" + strjoin(e, '; ') + "): MATLAB normalisation only."; end
    ask(opt, 'TRACKING step 2: connect the DUT between port 1 (source) and port 2 (input), then press Enter.');
    if normOK
        yNorm = acquire(ff, opt.Average, numel(f));         % instrument-normalised trace (dB)
        ff.write("SOUR:NORM 0"); ff.opc(30);
    end
else
    yThru = pw*ones(1, numel(f));
    R.notes(end + 1, 1) = "Tracking without a thru: S21 = trace - source power (test cables included, about +-1 dB).";
end
yDut = acquire(ff, opt.Average, numel(f));
P(end + 1, :) = yDut; t(end + 1, 1) = datetime('now'); rows(end + 1, 1) = "dut";
clear rfOff
S21_dB = yDut - yThru;
if ~isempty(yNorm)
    ok = S21_dB > min(S21_dB) + 20;              % well above the noise floor
    if any(ok)
        R.notes(end + 1, 1) = sprintf("FieldFox SOUR:NORM trace vs MATLAB normalisation: median difference %.2f dB.", ...
            median(abs(yNorm(ok) - S21_dB(ok))));
    else
        R.notes(end + 1, 1) = "SOUR:NORM accepted, but no transmission above the noise to compare with.";
    end
end

st = rmfield(s, 'f'); st.method = "tracking"; st.sourcePower = pw; st.rows = rows;
st.thru = thru; st.time = datetime('now');
R.settings = st; R.method = "tracking"; R.scalar = true;
R.cal = struct('text', "tracking generator, " + string(ifElse(thru, "normalised to a thru", "not normalised")));
rec = struct('f', f, 'P', P, 't', t, 'settings', st, 'idn', ff.IDN, 'S21_dB', S21_dB);
file = fullfile(opt.OutDir, "ff_vna_track_" + stamp() + ".mat");
save(file, '-struct', 'rec');
R.otherFiles(end + 1, 1) = file;
fprintf('Tracking recording saved: %s\n', file);
S = scalarS21(f, S21_dB);
end

function y = acquire(ff, n, npts)
% n single sweeps averaged in power (dBm)
Y = zeros(max(1, n), npts);
for k = 1:max(1, n)
    opcQuery(ff, "INIT:IMM", 300);
    yk = ff.saTrace();
    if numel(yk) ~= npts
        error('ff_vna:data', 'SA trace: %d values received, %d expected (timeout or instrument busy).', numel(yk), npts);
    end
    Y(k, :) = yk;
end
y = 10*log10(mean(10.^(Y/10), 1));
end

function how = quickCal(ff, opt)
% Guided QuickCal (option 112). The cal set is used for this measurement;
% release() reloads the saved instrument state afterwards.
two = opt.Cal == "quick2";
ports = 1; if two, ports = [1 2]; end
if ~prompting(opt)
    % Unattended, nobody can confirm that the ports are open (a DUT left on the
    % port during the internal-standard step gives a wrong cal) or connect a thru.
    warning('ff_vna:cal', ['QuickCal "%s" needs someone to connect things and cannot run with ' ...
        'Prompt=false: skipped. For a scripted cal pass a function handle as Prompt.'], opt.Cal);
    how = "skipped (needs Prompt)"; return
end
conn = string(opt.Connector);
if isscalar(conn), conn = [conn conn]; end
for p = ports, ff.write(sprintf('SENS:CORR:COLL:CONN %d,"%s"', p, conn(p))); end
% Insertable pair (one male, one female): internal standards on port 1 only (QSOLT),
% otherwise on both ports (SOLR, any thru adapter).
insertable = two && xor(contains(conn(1), "-M-"), contains(conn(2), "-M-"));
intPorts = ports; if insertable, intPorts = 1; end
ff.write("SENS:CORR:COLL:METH:QCAL:CAL " + strjoin(string(ports), ","));
e = ff.errors();
if ~isempty(e), error('ff_vna:cal', 'QuickCal setup failed: %s', strjoin(e, '; ')); end

fprintf('\nQuickCal (%s), DUT connector: %s\n', opt.Cal, strjoin(conn(ports), ' / '));
ask(opt, sprintf('CAL step 1: leave port %s OPEN (nothing connected), then press Enter.', strjoin(string(ports), ' and ')));
for p = intPorts
    % The verified PG example measures INT on port 1 only (insertable thru); INT 2
    % for a non-insertable pair is UNVERIFIED on this unit, so check it at once.
    opcQuery(ff, sprintf("SENS:CORR:COLL:INT %d", p), 120);
    e = ff.errors();
    if ~isempty(e), warning('ff_vna:cal', 'QuickCal internal standards, port %d: %s', p, strjoin(e, '; ')); end
end
loads = 0;
if prompting(opt)
    for p = ports
        r = ask(opt, sprintf('CAL (optional): connect a 50 ohm LOAD to port %d and press Enter, or type s + Enter to skip.', p));
        if ~startsWith(lower(strtrim(r)), "s")
            opcQuery(ff, sprintf("SENS:CORR:COLL:LOAD %d", p), 120); loads = loads + 1;
        end
    end
end
if two
    ask(opt, 'CAL step 3: connect port 1 to port 2 with a THRU (cable or adapter), then press Enter.');
    opcQuery(ff, "SENS:CORR:COLL:THRU 1,2", 120);
end
opcQuery(ff, "SENS:CORR:COLL:SAVE 0", 60);
e = ff.errors();
if ~isempty(e), warning('ff_vna:cal', 'QuickCal reported: %s', strjoin(e, '; ')); end
how = sprintf("QuickCal %s (%d load%s)", opt.Cal, loads, ifElse(loads == 1, "", "s"));
fprintf('%s done.\n', how);
ask(opt, 'Now connect the DUT, then press Enter to measure.');
end

function c = calState(ff, done)
c.corr = boolQuery(ff, "SENS:CORR?");          % all correction incl. CalReady
c.user = boolQuery(ff, "SENS:CORR:USER?");     % user calibration
c.performed = string(done);
if c.user && strlength(c.performed) > 0
    c.text = "user cal ON: " + c.performed + ", reference plane where the standards were connected";
elseif c.user
    c.text = "user cal ON (made earlier, details cannot be read over SCPI)";
elseif c.corr
    c.text = "factory CalReady only: reference plane at the FieldFox test ports";
else
    c.text = "correction OFF (raw data)";
end
if startsWith(c.performed, "skipped"), c.text = c.text + " (QuickCal skipped)"; end
end

function tf = boolQuery(ff, cmd)
r = upper(strtrim(erase(string(ff.query(cmd)), '"')));
tf = any(r == ["1" "+1" "ON"]);
end

function opcQuery(ff, cmd, tmo)
% <cmd>;*OPC? that fails loudly when the instrument does not finish in time.
r = string(ff.query(cmd + ";*OPC?", tmo));
if isempty(r) || strtrim(r) ~= "1"
    error('ff_vna:timeout', 'No completion reply to "%s" within %.0f s.', cmd, tmo);
end
end

function v = readNum(ff, cmd, default)
v = ff.queryNum(cmd);
if isempty(v) || ~isfinite(v(1)), v = default; else, v = v(1); end
end

function warnErrors(ff, where)
e = ff.errors();
if ~isempty(e), warning('ff_vna:instrument', '%s: instrument reported %s', where, strjoin(e, '; ')); end
end

function sourceOff(ff)
try
    ff.write("SOUR:ENAB 0");
catch err
    warning('ff_vna:source', 'Could not switch the tracking generator off: %s', err.message);
end
end

function releaseSafe(ff, mine)
if ~isvalid(ff), return, end
resync(ff);
if ~mine, return, end
try
    ff.release();
    fprintf('Instrument restored (mode %s).\n', ff.mode());
catch err
    warning('ff_vna:release', 'Could not restore the instrument: %s', err.message);
end
end

function resync(ff)
% After a timeout the late reply is still on its way: skip everything up to the
% answer of a fresh *IDN?, so that release() reads its own replies. Gives up
% after two silent read timeouts (dead connection).
ws = warning('off', 'transportlib:client:ReadlineWarning');
restoreWarn = onCleanup(@() warning(ws));
try
    flush(ff.Client, "input");
    ff.write("*IDN?");
    r = ""; silent = 0;
    while ~startsWith(r, ff.IDN) && silent < 2
        r = string(readline(ff.Client));
        if isempty(r), r = ""; silent = silent + 1; end
        r = strtrim(r);
    end
catch err
    warning('ff_vna:resync', 'Connection problem: %s', err.message);
end
end

function rfWarning(opt, msg)
fprintf(['\n*** RF OUTPUT ***  %s\n    Disconnect LISNs, limiters and receivers that must not be driven.' ...
    '\n    Port damage level: +27 dBm, 50 VDC.\n'], msg);
if prompting(opt)
    r = ask(opt, 'Press Enter to continue, or type n + Enter to cancel.');
    if startsWith(lower(strtrim(r)), "n"), error('ff_vna:cancel', 'Cancelled.'); end
end
end

function tf = prompting(opt)
tf = isa(opt.Prompt, 'function_handle') || (islogical(opt.Prompt) || isnumeric(opt.Prompt)) && opt.Prompt;
end

function r = ask(opt, msg)
r = "";
if isa(opt.Prompt, 'function_handle')
    fprintf('%s\n', msg);
    r = string(opt.Prompt(msg));
elseif prompting(opt)
    r = string(input(char(string(msg) + " "), 's'));
else
    fprintf('%s  [Prompt=false: not waiting]\n', msg);
end
end

function [f1, f2] = sweepRange(opt)
def = [1e6 1e9];
hasAnt = ~(isstring(opt.Antenna) || ischar(opt.Antenna)) || strlength(string(opt.Antenna)) > 0;
if hasAnt && ~isempty(opt.DesignFrequency), def = [0.5 1.5]*opt.DesignFrequency; end
f1 = opt.Start; f2 = opt.Stop;
if isempty(f1), f1 = def(1); end
if isempty(f2), f2 = def(2); end
f1 = max(f1, 30e3); f2 = min(f2, 6.5e9);
if f2 <= f1, error('ff_vna:option', 'Stop must be above Start (30 kHz ... 6.5 GHz).'); end
end

%% ===== offline data =======================================================

function [S, R] = loadFile(file, R)
[~, ~, ext] = fileparts(file);
if strcmpi(ext, ".mat")
    L = load(file);
    if isfield(L, 'S') && isa(L.S, 'sparameters')
        S = L.S;
    elseif isfield(L, 'R') && isfield(L.R, 'S') && isa(L.R.S, 'sparameters')
        S = L.R.S;
    elseif isfield(L, 'S21_dB') && isfield(L, 'f')
        S = scalarS21(L.f, L.S21_dB); R.scalar = true;
    else
        error('ff_vna:file', '%s has no sparameters object S and no tracking data (f, S21_dB).', file);
    end
    if isfield(L, 'idn'), R.idn = string(L.idn); end
    if isfield(L, 'settings'), R.settings = L.settings; end
else
    S = sparameters(char(file));
end
if S.NumPorts > 2
    S = snp2smp(S, [1 2]);
    R.notes(end + 1, 1) = "File has more than 2 ports: ports 1-2 analysed, the others terminated in 50 ohm.";
end
R.method = "file"; R.source = file;
fprintf('Loaded %s: %d-port, %d points, %s - %s\n', file, S.NumPorts, numel(S.Frequencies), ...
    fmtHz(S.Frequencies(1)), fmtHz(S.Frequencies(end)));
end

function [S, R] = demoData(opt, R)
% Synthetic "measurements" built from RF Toolbox models, seen through a 50 cm cable.
[f1, f2] = sweepRange(opt);
f = linspace(f1, f2, opt.Points).';
f0 = opt.DesignFrequency; if isempty(f0), f0 = 433e6; end
rs = RandStream('mt19937ar', 'Seed', 1);
noise = @(n) 2e-3*complex(randn(rs, n), randn(rs, n));      % VNA trace noise
cable = sparameters(txlineCoaxial('LineLength', 0.5, 'LossTangent', 2e-4, 'SigmaCond', 5.8e7), f);
switch opt.Demo
    case "antenna"          % series RLC resonance (60 ohm, Q = 8) at f0
        L = 8*60/(2*pi*f0); C = 1/((2*pi*f0)^2*L);
        ZL = 60 + 1j*2*pi*f*L + 1./(1j*2*pi*f*C);
        g = gammain(cable, ZL) + noise([numel(f) 1]);
        S = sparameters(reshape(g, 1, 1, []), f, 50);
    case "filter"           % 5th-order Chebyshev band-pass, +-3% around f0
        flt = rffilter('ResponseType', 'Bandpass', 'FilterType', 'Chebyshev', 'FilterOrder', 5, ...
            'PassbandFrequency', f0*[0.97 1.03], 'PassbandAttenuation', 0.5);
        S = cascadesparams(cable, sparameters(flt, f), cable);
        S = sparameters(S.Parameters + reshape(noise([4*numel(f) 1]), 2, 2, []), f, 50);
    otherwise
        error('ff_vna:option', 'Demo must be "antenna" or "filter".');
end
R.method = "demo"; R.source = "Demo " + opt.Demo;
R.cal = struct('text', "synthetic data");
fprintf('Demo data: %s at %s, %s - %s\n', opt.Demo, fmtHz(f0), fmtHz(f1), fmtHz(f2));
end

function S = scalarS21(f, s21dB)
% Scalar (magnitude-only) S21 as a 2-port; entries that were not measured are NaN.
D = nan(2, 2, numel(f));
D(2, 1, :) = 10.^(s21dB(:)/20);
S = sparameters(D, f(:), 50);
end

function [S, R] = removeFixture(S, opt, R)
f = S.Frequencies;
F1 = sparameters(char(opt.Fixture));
if f(1) < F1.Frequencies(1) || f(end) > F1.Frequencies(end)
    warning('ff_vna:fixture', 'Fixture file covers %s - %s only: outside points become NaN.', ...
        fmtHz(F1.Frequencies(1)), fmtHz(F1.Frequencies(end)));
end
F1 = rfinterp1(F1, f);
if S.NumPorts == 1
    % 1-port: invert the fixture 2-port (port 2 at the DUT)
    p = F1.Parameters; gm = S.Parameters;
    gd = (gm - p(1,1,:)) ./ (p(1,2,:).*p(2,1,:) + p(2,2,:).*(gm - p(1,1,:)));
    S = sparameters(gd, f, 50);
else
    if strlength(string(opt.Fixture2)) > 0, F2 = rfinterp1(sparameters(char(opt.Fixture2)), f); else, F2 = F1; end
    S = deembedsparams(S, F1, snp2smp(F2, [2 1]));   % port 2 side is mirrored
end
R.notes(end + 1, 1) = "Fixture removed: " + string(opt.Fixture);
end

%% ===== analysis ===========================================================

function R = analyse(S, R, opt)
f = S.Frequencies; n = S.NumPorts;
[fs, fu] = freqUnit(f(end));
M = R.metrics;
M = addMetric(M, "Start frequency", f(1), "Hz");
M = addMetric(M, "Stop frequency", f(end), "Hz");
M = addMetric(M, "Points", numel(f), "");
vis = opt.Visible;
oldVis = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', vis);          % toolbox plots that open their own figure
restoreVis = onCleanup(@() set(groot, 'DefaultFigureVisible', oldVis));
newFig = @(name) figure('Name', "ff_vna - " + name, 'NumberTitle', 'off', 'Color', 'w', 'Visible', vis);
label = R.source;
rlAx = [];

%% Antenna simulation first (slow), so it can be overlaid
if hasAntenna(opt)
    if R.scalar
        R.notes(end + 1, 1) = "Antenna comparison needs S11 (NA data): skipped.";
    else
        try
            R.antenna = simulateAntenna(S, opt);
        catch err
            warning('ff_vna:antenna', 'Antenna simulation failed: %s', err.message);
            R.notes(end + 1, 1) = "Antenna simulation failed: " + err.message;
        end
    end
end

%% All parameters
fig = newFig("S-parameters"); R.figures(end + 1) = fig;
rfplot(S); grid on; hold on
if ~isempty(R.antenna)
    h = rfplot(R.antenna.S, 1, 1, 'db', '--');
    h.DisplayName = "Simulated " + R.antenna.name + " S11";
end
title(sprintf('%s  |  %s', label, R.cal.text), 'Interpreter', 'none', 'FontSize', 9);

%% Reflection: return loss, VSWR, Smith chart, rational fit
if ~R.scalar
    for p = 1:n
        g = rfparam(S, p, p); nm = sprintf('S%d%d', p, p);
        rl = -20*log10(abs(g));
        [rlBest, ib] = max(rl);
        sw = vswr(g);
        [lo, hi] = bandAbove(rl, f, ib, 10);
        M = addMetric(M, nm + " best match frequency", f(ib), "Hz");
        M = addMetric(M, nm + " best return loss", rlBest, "dB");
        M = addMetric(M, nm + " min VSWR", sw(ib), "");
        M = addMetric(M, nm + " -10 dB bandwidth", hi - lo, "Hz");
        M = addMetric(M, nm + " -10 dB low edge", lo, "Hz");
        M = addMetric(M, nm + " -10 dB high edge", hi, "Hz");
        if rlBest < 3
            R.notes(end + 1, 1) = sprintf(['%s: |%s| is close to 1 at all frequencies (return loss <= %.1f dB): ' ...
                'port %d is open, shorted or has nothing connected, so match, bandwidth and ' ...
                'matching-network results mean nothing here.'], nm, nm, rlBest, p);
        end
    end
    % Rational model. S11 does not tend to zero at high frequency; the tolerance
    % and pole limit stop it from fitting measurement noise (default: 100s of poles).
    try
        R.fit = rational(S, 'TendsToZero', false, 'Tolerance', -30, 'MaxPoles', 40);
        M = addMetric(M, "Rational fit poles", numel(R.fit.Poles), "");
        M = addMetric(M, "Rational fit error", R.fit.ErrDB, "dB");
        if R.fit.ErrDB > -20
            R.notes(end + 1, 1) = sprintf(['Rational fit error %.1f dB is poor: the data are mostly noise ' ...
                '(e.g. |S21| at the noise floor), so do not use the model.'], R.fit.ErrDB);
        end
    catch err
        R.notes(end + 1, 1) = "Rational fit failed: " + err.message;
    end
    try
        q = tdr.sParamQuality(S);
        R.quality = q;
        M = addMetric(M, "Passivity", q.PassivityPercentage, "%");
        M = addMetric(M, "Causality", q.CausalityPercentage, "%");
        if n == 2, M = addMetric(M, "Reciprocity", q.ReciprocityPercentage, "%"); end
    catch err
        R.notes(end + 1, 1) = "sParamQuality failed: " + err.message;
    end

    fig = newFig("Return loss and VSWR"); R.figures(end + 1) = fig;
    ax1 = subplot(2, 1, 1, 'Parent', fig); hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
    ax2 = subplot(2, 1, 2, 'Parent', fig); hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
    for p = 1:n
        g = rfparam(S, p, p);
        plot(ax1, f/fs, -20*log10(abs(g)), 'LineWidth', 1.3, 'DisplayName', sprintf('S%d%d measured', p, p));
        plot(ax2, f/fs, min(vswr(g), 20), 'LineWidth', 1.3, 'DisplayName', sprintf('VSWR port %d', p));
    end
    if ~isempty(R.fit)
        r = freqresp(R.fit, f);                  % K x 1 for 1 port, 2 x 2 x K for 2 ports
        if ndims(r) == 3, r = squeeze(r(1, 1, :)); end
        plot(ax1, f/fs, -20*log10(abs(r(:))), 'k:', 'LineWidth', 1.2, 'DisplayName', 'S11 rational fit');
    end
    rl1 = -20*log10(abs(rfparam(S, 1, 1)));
    [~, ib] = max(rl1);
    [lo, hi] = bandAbove(rl1, f, ib, 10);
    if hi > lo, xregion(ax1, lo/fs, hi/fs, 'FaceColor', [0.3 0.8 0.3], 'FaceAlpha', 0.15, 'DisplayName', 'S11 >= 10 dB'); end
    plot(ax1, f(ib)/fs, rl1(ib), 'kv', 'MarkerFaceColor', 'y', 'DisplayName', 'Best match');
    yline(ax1, 10, '--', '10 dB (VSWR 1.9)', 'HandleVisibility', 'off');
    if ~isempty(R.antenna)
        plot(ax1, R.antenna.S.Frequencies/fs, -20*log10(abs(rfparam(R.antenna.S, 1, 1))), '--', ...
            'LineWidth', 1.2, 'DisplayName', "Simulated " + R.antenna.name);
    end
    ylabel(ax1, 'Return loss (dB)'); legend(ax1, 'Location', 'best');
    title(ax1, 'Return loss (higher = better match)');
    ylabel(ax2, 'VSWR'); xlabel(ax2, "Frequency (" + fu + ")"); ylim(ax2, [1 10]);
    yline(ax2, 2, '--', 'HandleVisibility', 'off'); legend(ax2, 'Location', 'best');
    rlAx = ax1;

    fig = newFig("Smith chart"); R.figures(end + 1) = fig;
    hs = smithplot(S, 1, 1); labels = "S11";
    if n == 2, add(hs, f, rfparam(S, 2, 2)); labels(end + 1) = "S22"; end
    if ~isempty(R.antenna)
        add(hs, R.antenna.S.Frequencies, rfparam(R.antenna.S, 1, 1));
        labels(end + 1) = "Simulated " + R.antenna.name;
    end
    hs.LegendLabels = cellstr(labels);
    R.nport = nport(S);                       % rfbudget-ready element of the measured network
end

%% Transmission
if n == 2
    s21 = rfparam(S, 2, 1); g21 = 20*log10(abs(s21));
    [pk, ipk] = max(g21);
    [bw, lo, hi] = powerbw(abs(s21).^2, f, [], 3);
    M = addMetric(M, "S21 peak frequency", f(ipk), "Hz");
    M = addMetric(M, "Insertion loss at peak", -pk, "dB");
    M = addMetric(M, "S21 3 dB bandwidth", bw, "Hz");
    M = addMetric(M, "S21 3 dB low edge", lo, "Hz");
    M = addMetric(M, "S21 3 dB high edge", hi, "Hz");
    if pk < -40
        R.notes(end + 1, 1) = sprintf(['|S21| stays below %.0f dB: nothing (or a blocking DUT) connected between ' ...
            'port 1 and port 2; bandwidth and group delay are noise.'], ceil(pk));
    elseif lo <= f(1) || hi >= f(end)
        R.notes(end + 1, 1) = "S21 3 dB band reaches the sweep edge: widen Start/Stop for the true bandwidth.";
    end
    fig = newFig("Transmission"); R.figures(end + 1) = fig;
    nRows = 1 + 2*~R.scalar;
    ax = subplot(nRows, 1, 1, 'Parent', fig); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    plot(ax, f/fs, g21, 'LineWidth', 1.3, 'DisplayName', '|S21|');
    if ~R.scalar, plot(ax, f/fs, 20*log10(abs(rfparam(S, 1, 2))), '--', 'DisplayName', '|S12|'); end
    xregion(ax, lo/fs, hi/fs, 'FaceColor', [0.3 0.6 1], 'FaceAlpha', 0.15, 'DisplayName', '3 dB band');
    ylabel(ax, 'dB'); legend(ax, 'Location', 'best');
    title(ax, sprintf('Insertion loss %.2f dB at %s, 3 dB bandwidth %s', -pk, fmtHz(f(ipk)), fmtHz(bw)));
    if ~R.scalar
        gd = groupdelay(S, 2, 1, 'Aperture', 2*median(diff(f)));
        inBand = f >= lo & f <= hi;
        M = addMetric(M, "S21 group delay (median in 3 dB band)", median(gd(inBand)), "s");
        try
            R.budget = rfbudget(R.nport, f(ipk), -30, 1e6);    % e.g. -30 dBm in, 1 MHz signal
            M = addMetric(M, "rfbudget transducer gain at peak", R.budget.TransducerGain(end), "dB");
        catch err
            R.notes(end + 1, 1) = "rfbudget failed: " + err.message;
        end
        ax = subplot(nRows, 1, 2, 'Parent', fig);
        set(fig, 'CurrentAxes', ax);           % rfplot draws into the current axes
        rfplot(S, 2, 1, 'angle'); grid(ax, 'on'); title(ax, 'S21 phase');
        ax = subplot(nRows, 1, 3, 'Parent', fig);
        plot(ax, f/fs, gd*1e9, 'LineWidth', 1.2); grid(ax, 'on');
        ylabel(ax, 'Group delay (ns)'); xlabel(ax, "Frequency (" + fu + ")");
        if any(inBand), ylim(ax, [min(0, min(gd(inBand)*1e9)) max(1, 3*max(gd(inBand)*1e9))]); end
    else
        xlabel(ax, "Frequency (" + fu + ")");
    end
end

%% Time domain on S11: TDR (Signal Integrity tdr) and windowed-IFFT distance to fault
if ~R.scalar
    try
        [R, M] = timeDomain(S, R, M, opt, newFig);
    catch err
        warning('ff_vna:timedomain', 'TDR/DTF failed: %s', err.message);
        R.notes(end + 1, 1) = "TDR/DTF failed: " + err.message;
    end
end

%% Matching network for the measured S11
if ~R.scalar && ~isempty(opt.MatchFrequency)
    try
        [R, M] = designMatch(S, R, M, opt);
        if ~isempty(rlAx) && isvalid(rlAx)
            plot(rlAx, f/fs, -20*log10(abs(R.match.gammaIn)), 'm-.', 'LineWidth', 1.3, ...
                'DisplayName', "With " + R.match.topology + " match (predicted)");
            xline(rlAx, opt.MatchFrequency/fs, 'm:', 'HandleVisibility', 'off');
        end
        smithplot(R.match.network, 'CircuitIndex', 1);    % opens its own figure
        R.figures(end + 1) = gcf;
        set(gcf, 'Name', 'ff_vna - Matching network (impedance path)', 'NumberTitle', 'off');
    catch err
        warning('ff_vna:match', 'Matching network not designed: %s', err.message);
        R.notes(end + 1, 1) = "Matching network not designed: " + err.message;
    end
end

%% Antenna comparison numbers
if ~isempty(R.antenna)
    A = R.antenna;
    M = addMetric(M, "Antenna sim best match", A.fSimBest, "Hz");
    M = addMetric(M, "Antenna resonantFrequency (X=0)", A.fRes, "Hz");
    M = addMetric(M, "Measured best match (sim band)", A.fMeasBest, "Hz");
    M = addMetric(M, "Resonance error meas vs sim", A.errorPct, "%");
    M = addMetric(M, "Sim |Z| at design frequency", abs(A.Z0), "ohm");
    M = addMetric(M, "Sim vs meas |S11| RMS difference", A.rmsDb, "dB");
    if isnan(A.fRes)
        R.notes(end + 1, 1) = "resonantFrequency found no reactance zero inside the simulated band " + ...
            "(only the band edge): use the simulated best match instead.";
    end
end
R.metrics = M;
end

function [R, M] = timeDomain(S, R, M, opt, newFig)
f = S.Frequencies; K = numel(f);
c0 = 299792458; vf = opt.VelocityFactor;
if max(abs(diff(f, 2))) > 1e-6*(f(2) - f(1))        % DTF needs equal steps
    f = linspace(f(1), f(end), K).';
    S = rfinterp1(S, f);
end
S1 = S; if S.NumPorts == 2, S1 = snp2smp(S, 1); end     % S11 with port 2 in 50 ohm
s11 = rfparam(S1, 1, 1);
df = f(2) - f(1); span = f(end) - f(1);
dRes = c0*vf/(2*span); dAlias = c0*vf/(2*df);
dMax = dAlias;
if ~isempty(opt.MaxDistance), dMax = min(opt.MaxDistance, dAlias); end

% Distance to fault: Kaiser window, zero-padded IFFT (band-pass mode). The IFFT is
% periodic: the last bins are small negative distances (a reflection right at the
% reference plane, e.g. after fixture removal, spreads to both sides of 0).
win = kaiser(K, 6);
nfft = 2^nextpow2(16*K);
h = ifft(s11.*win/mean(win), nfft)*nfft/K;        % |h| = reflection coefficient
d = (0:nfft - 1).'*c0*vf/(2*nfft*df);
guard = min(3*dRes, dAlias/4);
d(d >= dAlias - guard) = d(d >= dAlias - guard) - dAlias;
dMax = min(dMax, dAlias - guard);
[d, o] = sort(d); h = h(o);
keep = d <= dMax;
d = d(keep); gDb = 20*log10(abs(h(keep)) + eps);
[pk, ip] = max(gDb);
R.dtf = struct('distance', d, 'gammaDb', gDb, 'resolution', dRes, 'aliasFree', dAlias, 'velocityFactor', vf);
M = addMetric(M, "DTF strongest reflection at", d(ip), "m");
M = addMetric(M, "DTF strongest reflection level", pk, "dB");
M = addMetric(M, "DTF resolution", dRes, "m");
M = addMetric(M, "DTF alias-free range", dAlias, "m");

% TDR: step response of a rational fit (SI Toolbox), rise time from the bandwidth
tr = 0.5/f(end); te = 2*dMax/(c0*vf);
T = tdr(S1, RiseTime=tr, SampleTime=max(tr/4, te/20000), EndTime=te);
dT = T.Time*c0*vf/2; Z = T.ImpedanceTDR(:, 1);
rho = min(max((Z - 50)./(Z + 50), -0.99), 0.99);   % fit ripple can give |rho| > 1 near an open/short
Z = 50*(1 + rho)./(1 - rho);                        % so Z is shown within 0.25 ... 9950 ohm
R.tdr = T;
M = addMetric(M, "TDR fit error", T.Fit(1).ErrDB, "dB");
iNear = find(dT >= 2*dRes, 1); if isempty(iNear), iNear = 1; end
M = addMetric(M, "TDR impedance just after the port", Z(iNear), "ohm");
if f(1) > f(end)/100
    R.notes(end + 1, 1) = "TDR needs a low start frequency (<= Stop/100); here it is only indicative - trust the DTF plot.";
end

fig = newFig("Time domain"); R.figures(end + 1) = fig;
ax = subplot(2, 1, 1, 'Parent', fig);
plot(ax, dT, Z, 'LineWidth', 1.2); grid(ax, 'on');
ylim(ax, [0 200]); yline(ax, 50, '--', '50 \Omega');
ylabel(ax, 'Impedance (\Omega)');
title(ax, sprintf('TDR (Signal Integrity tdr, rise time %.3g ns, vf %.2f)', tr*1e9, vf));
ax = subplot(2, 1, 2, 'Parent', fig);
plot(ax, d, gDb, 'LineWidth', 1.2); grid(ax, 'on'); hold(ax, 'on');
plot(ax, d(ip), pk, 'kv', 'MarkerFaceColor', 'y');
ylabel(ax, '|\Gamma| (dB)'); xlabel(ax, 'Distance from the reference plane (m)');
title(ax, sprintf('Distance to fault (Kaiser-windowed IFFT): resolution %.2f m, range %.1f m', dRes, dAlias));
end

function [R, M] = designMatch(S, R, M, opt)
f = S.Frequencies; f0 = opt.MatchFrequency;
if f0 < f(1) || f0 > f(end), error('ff_vna:match', 'MatchFrequency is outside the measured range.'); end
S1 = S; if S.NumPorts == 2, S1 = snp2smp(S, 1); end
g = rfparam(S1, 1, 1);
topo = string(opt.MatchTopology);
mn = matchingnetwork('SourceImpedance', 50, 'LoadImpedance', S1, 'CenterFrequency', f0, ...
    'Bandwidth', min(0.1*f0, (f(end) - f(1))/2), 'Components', topo);
tbl = circuitDescriptions(mn);
Smn = sparameters(mn, f, 50, 1);           % best design; port 2 goes to the load
gc = min(abs(g), 0.9999).*exp(1j*angle(g));   % avoid infinite impedance at |S11| = 1
gin = gammain(Smn, 50*(1 + gc)./(1 - gc));
rlIn = -20*log10(abs(gin)); rl0 = -20*log10(abs(g));
[~, i0] = min(abs(f - f0));
[lo, hi] = bandAbove(rlIn, f, i0, 10);
R.match = struct('network', mn, 'circuits', tbl, 'topology', topo, 'frequency', f0, ...
    'S', Smn, 'gammaIn', gin, 'circuit', exportCircuits(mn, 1));
M = addMetric(M, "Match: return loss before at f0", rl0(i0), "dB");
M = addMetric(M, "Match: predicted return loss after at f0", rlIn(i0), "dB");
M = addMetric(M, "Match: predicted -10 dB bandwidth", hi - lo, "Hz");
fprintf('\nMatching network (%s) for %s, best design first (component 1 is at the source side):\n', topo, fmtHz(f0));
disp(tbl);
end

function A = simulateAntenna(S, opt)
f = S.Frequencies; s11 = rfparam(S, 1, 1);
[~, ib] = min(abs(s11));
f0 = opt.DesignFrequency;
if isempty(f0)
    f0 = f(ib);
    fprintf('No DesignFrequency: designing the antenna for the measured best match %s.\n', fmtHz(f0));
end
ant = opt.Antenna;
if isstring(ant) || ischar(ant)
    name = string(ant);
    fprintf('Designing %s for %s (Antenna Toolbox) ...\n', name, fmtHz(f0));
    ant = design(feval(char(name)), f0);
else
    name = string(class(ant));
end
lo = 0.7*f0; hi = 1.3*f0;
if f(1) < f0 && f(end) > f0, lo = max(lo, f(1)); hi = min(hi, f(end)); end
fa = linspace(lo, hi, opt.AntennaPoints).';
fprintf('Simulating %s at %d frequencies (can take a minute) ...\n', name, numel(fa));
tic;
Ssim = sparameters(ant, fa, 50);
Z0 = impedance(ant, f0);
try
    fRes = resonantFrequency(ant, [lo hi]);
    fRes = fRes(fRes > lo*1.005 & fRes < hi*0.995);     % a band edge is no resonance
    [~, k] = min(abs(fRes - f0)); fRes = fRes(k);       % X = 0 closest to the design frequency
    if isempty(fRes), fRes = NaN; end
catch
    fRes = NaN;
end
fprintf('Simulation took %.0f s.\n', toc);
fFine = linspace(lo, hi, 2001);                  % best match between the sim points
[~, is] = min(interp1(fa, 20*log10(abs(rfparam(Ssim, 1, 1))), fFine, 'spline'));
inB = find(f >= lo & f <= hi);
fMeas = NaN;
if ~isempty(inB), [~, im] = min(abs(s11(inB))); fMeas = f(inB(im)); end
Sm = rfinterp1(S, fa);
d = 20*log10(abs(rfparam(Sm, 1, 1))) - 20*log10(abs(rfparam(Ssim, 1, 1)));
A = struct('name', name, 'object', ant, 'designFrequency', f0, 'S', Ssim, 'Z0', Z0, ...
    'fRes', fRes, 'fSimBest', fFine(is), 'fMeasBest', fMeas, ...
    'errorPct', 100*(fMeas - fFine(is))/fFine(is), 'rmsDb', sqrt(mean(d(~isnan(d)).^2)));
end

function tf = hasAntenna(opt)
a = opt.Antenna;
tf = ~(isstring(a) || ischar(a)) || strlength(string(a)) > 0;
end

%% ===== helpers ============================================================

function [lo, hi] = bandAbove(y, f, i0, thr)
% Edges (linear interpolation) of the contiguous region around i0 where y >= thr.
lo = NaN; hi = NaN;
if y(i0) < thr, return, end
a = i0; while a > 1 && y(a - 1) >= thr, a = a - 1; end
b = i0; while b < numel(y) && y(b + 1) >= thr, b = b + 1; end
lo = f(a); hi = f(b);
if a > 1, lo = interp1(y([a - 1 a]), f([a - 1 a]), thr); end
if b < numel(y), hi = interp1(y([b + 1 b]), f([b + 1 b]), thr); end
end

function M = newMetrics()
M = table(strings(0, 1), zeros(0, 1), strings(0, 1), 'VariableNames', {'Metric', 'Value', 'Unit'});
end

function M = addMetric(M, name, value, unit)
if isempty(value), value = NaN; end
M = [M; {string(name), double(value(1)), string(unit)}];
end

function printSummary(R)
fprintf('\n=== ff_vna: %s (%s) ===\n', R.source, R.method);
fprintf('Calibration: %s\n', R.cal.text);
if height(R.metrics) > 0
    T = R.metrics;
    txt = strings(height(T), 1);
    for k = 1:height(T)
        if T.Unit(k) == "Hz", txt(k) = fmtHz(T.Value(k));
        elseif T.Unit(k) == "s", txt(k) = sprintf('%.4g ns', T.Value(k)*1e9);
        else, txt(k) = strtrim(sprintf('%.4g %s', T.Value(k), T.Unit(k)));
        end
        fprintf('  %-42s %s\n', T.Metric(k), txt(k));
    end
end
for k = 1:numel(R.notes), fprintf('Note: %s\n', R.notes(k)); end
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
s = sprintf('%.6g %s', x/scale, unit);
end

function v = ifElse(c, a, b)
if c, v = a; else, v = b; end
end

function s = stamp()
s = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
end
