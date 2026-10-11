function R = ff_vvm_logger(ip, varargin)
%FF_VVM_LOGGER  Magnitude and phase vs time at one frequency (FieldFox vector voltmeter, option 308).
%
%   ff_vvm_logger                                    S11 at the FieldFox's current VVM frequency
%   ff_vvm_logger(ip, 'Parameter', "S21", 'Frequency', 100e6, 'Zero', true, 'Duration', 600)
%   R = ff_vvm_logger('File', "vvm_20261008_101500.mat")    re-analyse a saved log
%   R = ff_vvm_logger('Demo', true)                          synthetic data, no instrument
%
% Use it to test the phase and amplitude stability of a cable, connector,
% adapter or antenna feed: leave it, warm it, bend it, and see the drift.
%
% !! RF OUTPUT: VVM mode transmits a CW signal from Port 1 (S11 and S21, up to
% +3 dBm). Connect only the cable or device under test - not an antenna you may
% not transmit with, nor a receiver that cannot take it. S21 = Port 1 -> device
% -> Port 2. The RF stops when the tool ends (the previous mode is restored).
%
% Calibration: VVM uses the calibration made in CAT or NA mode (calibrate there
% first). Without one, the factory "CalReady" correction at the test ports is
% used. 'Zero', true makes the readings relative: the first reading becomes
% 0 dB / 0 deg (or 180 deg with 'ZeroRef', "SHORT").
%
% Results:
%   magnitude (dB) and unwrapped phase (deg) vs time
%   drift   straight-line fit with fitlm (Statistics and Machine Learning
%           Toolbox): deg/hour and dB/hour with 95% confidence intervals; the
%           phase drift is also given as a delay change (ps/hour, = -dphase /
%           (360 f); for S11 it is the round trip)
%   Allan deviation of the phase with allanvar (Navigation Toolbox), the
%           standard tool for sensor stability: at averaging time tau it is the
%           typical difference between two consecutive tau-long averages.
%           Falling as tau^-1/2 = white noise (averaging helps); flat = flicker
%           floor; rising = drift (temperature, movement). The readings are put
%           on a uniform time grid first (interp1); phase steps > 180 deg between
%           two readings cannot be unwrapped - use a shorter Interval then.
%
% Name-Value options (defaults):
%   'Parameter'  "S11"    "S11" | "S21" ("AB" | "BA" need an external source)
%   'Frequency'  []       CW frequency (Hz); [] = keep the FieldFox VVM setting
%   'Zero'       false    relative readings (zero at the start)
%   'ZeroRef'    "OPEN"   "OPEN" (0 deg) | "SHORT" (180 deg)
%   'Power'      -15      source power, -45 to +3 dBm (low by default; [] = keep the instrument's)
%   'IFBW'       []       IF bandwidth (Hz), 10 Hz - 100 kHz: lower = less noise; [] = keep
%   'Average'    []       averages; [] = keep
%   'Interval'   1        s between readings (0 = as fast as possible)
%   'Duration'   60       s to log (Inf = until the window is closed)
%   'Save'       true     MAT, CSV and PNG in OutDir (vvm_<time>.*)
%   'OutDir'     pwd
%   'File'       ""       analyse a saved vvm_*.mat instead of measuring
%   'Demo'       false    synthetic readings instead of the instrument
%   'Visible'    'on'     'off' for headless runs (the PNG is still saved)
%   'Prompt'     true     ask before the RF output is switched on; false = never ask
%   'Port'       5025
%
% Output R: t, mag_dB, phase_deg (as read, -180..180), phaseUnwrapped_deg,
% f (Hz), settings, drift (struct), adev (table Tau_s, ADEV_deg, N), idn, files.
% The MAT file uses the shared layout f, P, t, settings, idn - here P is the
% magnitude in dB (not dBm) and the phase is stored next to it.

if nargin < 1, ip = ""; end
opt = struct('Parameter', "S11", 'Frequency', [], 'Zero', false, 'ZeroRef', "OPEN", 'Power', -15, ...
    'IFBW', [], 'Average', [], 'Interval', 1, 'Duration', 60, 'Save', true, 'OutDir', pwd, ...
    'File', "", 'Demo', false, 'Visible', 'on', 'Prompt', true, 'Port', 5025);
[ip, opt] = parseArgs(ip, varargin, opt);
opt.Parameter = upper(string(opt.Parameter)); opt.ZeroRef = upper(string(opt.ZeroRef));
live = strlength(string(opt.File)) == 0 && ~opt.Demo;
if ~ismember(opt.Parameter, ["S11" "S21" "AB" "BA"])
    error('ff_vvm_logger:parameter', 'Parameter must be S11, S21, AB or BA.');
end
if ~ismember(opt.ZeroRef, ["OPEN" "SHORT"]), error('ff_vvm_logger:zeroref', 'ZeroRef must be OPEN or SHORT.'); end
if ~isempty(opt.Power) && (opt.Power < -45 || opt.Power > 3)
    error('ff_vvm_logger:power', 'Power must be between -45 and +3 dBm.');
end
if live && ~isfinite(opt.Duration) && strcmp(opt.Visible, 'off')
    error('ff_vvm_logger:duration', 'Duration = Inf needs a visible window (close it to stop).');
end

t = NaT(0, 1); mag = zeros(0, 1); ph = zeros(0, 1); settings = struct(); idn = "";
% A hidden figure is only made after the measurement: MATLAB renders pending
% figures during pause(), which can stall the readings by seconds.
fig = gobjects(0); h = struct(); figGuard = []; %#ok<NASGU>
if strcmp(opt.Visible, 'on'), [fig, h] = makeFigure(); end

if live
    fprintf('Connecting to %s:%d ...\n', ip, opt.Port);
    ff = FieldFox(ip, opt.Port);
    guard = instrumentGuard(ff);                        %#ok<NASGU> release + disconnect on any exit
    idn = ff.IDN;
    modes = ff.query("INST:CAT?");
    if ~contains(upper(modes), "VVM")
        error('ff_vvm_logger:noVVM', ['This FieldFox lists no VVM mode (INST:CAT? -> %s).\n' ...
            'The vector voltmeter needs option 308.'], modes);
    end
    fprintf(['\n*** RF OUTPUT ***  VVM mode transmits a CW signal from Port 1 (%s).\n' ...
        'Connect only the cable/device under test. The RF stops when this tool ends.\n\n'], opt.Parameter);
    warning('ff_vvm_logger:rf', 'VVM mode: RF output on Port 1.');
    if opt.Prompt && ~confirm('Switch the FieldFox to VVM mode and turn the RF output on?', false)
        R = []; closeFig(); return
    end
    ff.takeControl();
    ff.setMode("VVM");
    configureVVM(ff, opt);
    settings = readVVMSettings(ff);
    settings.interval = opt.Interval; settings.duration = opt.Duration;
    settings.sync = startSync(ff);
    settings.unitsCheck = unitsCheck(ff);
    if opt.Zero, zeroVVM(ff, opt.ZeroRef); end
    settings.zero = opt.Zero; settings.zeroRef = opt.ZeroRef;
    fprintf('VVM: %s at %s, power %g dBm, IFBW %s, zero %s\n', settings.parameter, ...
        fmtHz(settings.frequency), settings.power_dBm, fmtHz(settings.ifbw), string(opt.Zero));
    if hasFig(), setappdata(fig, 'busy', true); figBusy = busyGuard(fig); end %#ok<NASGU> closable after an error
    logLoop(ff);
    if hasFig(), setappdata(fig, 'busy', false); end
    if opt.Zero, ff.write("SENS:CORR:ZERO:STAT OFF"); end
    clear guard                                         % restores the FieldFox now (RF off)
elseif strlength(string(opt.File)) > 0
    [t, mag, ph, settings, idn] = loadLog(opt.File);
else
    [t, mag, ph, settings, idn] = demoLog(opt);
end

%% Analysis
if numel(ph) < 2, warning('ff_vvm_logger:empty', 'Not enough readings.'); R = []; closeFig(); return, end
phU = rad2deg(unwrap(deg2rad(ph)));
[drift, adev] = stability(t, mag, phU, settings.frequency);
if isempty(fig)
    [fig, h] = makeFigure();
    figGuard = hiddenFigureGuard(fig, opt.Visible);     %#ok<NASGU> no hidden figure left on errors
end
showResults();

R = struct('t', t, 'mag_dB', mag, 'phase_deg', ph, 'phaseUnwrapped_deg', phU, ...
    'f', settings.frequency, 'settings', settings, 'drift', drift, 'adev', adev, ...
    'idn', idn, 'files', strings(0, 1));
if opt.Save && (live || opt.Demo)
    R.files = saveLog();
elseif opt.Save && hasFig()
    [~, n] = fileparts(string(opt.File));
    png = fullfile(opt.OutDir, n + "_analysis.png");
    exportgraphics(fig, png, 'Resolution', 150); R.files = png;
    fprintf('Saved %s\n', png);
end
fprintf(['Phase drift %+.3g deg/h [%.3g, %.3g] (= %+.3g ps/h delay), magnitude drift %+.3g dB/h, ' ...
    'phase std %.3g deg\n'], drift.phase_deg_per_h, drift.phase_CI95, drift.delay_ps_per_h, ...
    drift.mag_dB_per_h, drift.phase_std_deg);
if median(mag, 'omitnan') < -60
    warning('ff_vvm_logger:weak', ['Median magnitude %.1f dB: the signal is at the noise floor, so the phase ' ...
        'is mostly noise and the drift / Allan numbers mean little. Is the device connected?'], median(mag, 'omitnan'));
end
if stopRequested(fig) || strcmp(opt.Visible, 'off'), closeFig(); end

%% ===== nested helpers (share the variables above) ==========================

    function logLoop(ff)
        tStart = tic; tPlot = tic; k = 0; tNext = 0;
        fprintf('Logging every %g s for %g s (close the window to stop) ...\n', opt.Interval, opt.Duration);
        while ~stopRequested(fig) && toc(tStart) < opt.Duration
            while toc(tStart) < tNext && ~stopRequested(fig)          % wait, keep the window alive
                pause(min(0.05, tNext - toc(tStart)));
            end
            if stopRequested(fig) || toc(tStart) >= opt.Duration, break, end
            if settings.sync, ff.sweep(60); end
            v = readValues(ff, "TRAC:DATA? FDATA", 2);  % magnitude (dB), phase (deg)
            k = k + 1;
            t(k, 1) = datetime('now'); mag(k, 1) = v(1); ph(k, 1) = v(2);
            tNext = nextSlot(tNext, opt.Interval, toc(tStart));
            if hasFig() && toc(tPlot) > 0.5             % plotting must not slow the readings
                updateTraces(); drawnow limitrate; tPlot = tic;
            end
        end
        fprintf('%d readings.\n', numel(ph));
    end

    function ok = hasFig()
        ok = ~isempty(fig) && isvalid(fig);
    end

    function [fig, h] = makeFigure()
        fig = figure('Name', 'FieldFox VVM logger', 'NumberTitle', 'off', 'Color', 'w', ...
            'Position', [80 80 1200 760], 'Visible', opt.Visible, 'CloseRequestFcn', @closeRequest);
        tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        h.mag = nexttile(tl); hold(h.mag, 'on'); grid(h.mag, 'on'); box(h.mag, 'on');
        h.magLine = plot(h.mag, NaT, NaN, '.-', 'Color', [0.3 0.5 1]);
        ylabel(h.mag, 'Magnitude (dB)');
        h.ph = nexttile(tl); hold(h.ph, 'on'); grid(h.ph, 'on'); box(h.ph, 'on');
        h.phLine = plot(h.ph, NaT, NaN, '.-', 'Color', [0.85 0.33 0.1]);
        ylabel(h.ph, 'Phase, unwrapped (deg)');
        h.adev = nexttile(tl); h.txt = nexttile(tl); axis(h.txt, 'off');
        title(h.adev, 'Allan deviation of the phase');
    end

    function updateTraces()
        if ~hasFig(), return, end
        set(h.magLine, 'XData', t, 'YData', mag);
        set(h.phLine, 'XData', t, 'YData', rad2deg(unwrap(deg2rad(ph))));
        title(h.mag, sprintf('%s magnitude, last %.3f dB', settings.parameter, mag(end)));
        title(h.ph, sprintf('%s phase at %s, last %.2f deg (%d readings)', settings.parameter, ...
            fmtHz(settings.frequency), ph(end), numel(ph)));
    end

    function showResults()
        if ~hasFig(), return, end
        updateTraces();
        th = hours(t - t(1));
        plot(h.ph, t, drift.phaseFit(th), 'k--', 'LineWidth', 1.2);
        plot(h.mag, t, drift.magFit(th), 'k--', 'LineWidth', 1.2);
        legend(h.ph, {'Phase', sprintf('Fit %+.3g deg/h', drift.phase_deg_per_h)}, 'Location', 'best');
        legend(h.mag, {'Magnitude', sprintf('Fit %+.3g dB/h', drift.mag_dB_per_h)}, 'Location', 'best');
        cla(h.adev);
        if ~isempty(adev)
            loglog(h.adev, adev.Tau_s, adev.ADEV_deg, 'o-', 'LineWidth', 1.3); hold(h.adev, 'on');
            loglog(h.adev, adev.Tau_s, adev.ADEV_deg(1)*sqrt(adev.Tau_s(1)./adev.Tau_s), 'k:');
            legend(h.adev, {'Phase ADEV', 'White noise (\tau^{-1/2})'}, 'Location', 'best');
            grid(h.adev, 'on');
        end
        xlabel(h.adev, '\tau (s)'); ylabel(h.adev, 'Allan deviation (deg)');
        title(h.adev, 'Allan deviation of the phase');
        txt = {sprintf('%s at %s', settings.parameter, fmtHz(settings.frequency))
            sprintf('Readings          %d over %.1f min', numel(ph), minutes(t(end) - t(1)))
            sprintf('Phase drift       %+.4g deg/h', drift.phase_deg_per_h)
            sprintf('   95%% CI         [%.4g, %.4g]', drift.phase_CI95)
            sprintf('   as delay       %+.4g ps/h', drift.delay_ps_per_h)
            sprintf('Phase change      %+.4g deg (fit, whole run)', drift.phase_change_deg)
            sprintf('Phase std (raw)   %.4g deg', drift.phase_std_deg)
            sprintf('Phase std (detr.) %.4g deg', drift.phase_resid_std_deg)
            sprintf('Mag. drift        %+.4g dB/h', drift.mag_dB_per_h)
            sprintf('   95%% CI         [%.4g, %.4g]', drift.mag_CI95)
            sprintf('Mag. std          %.4g dB', drift.mag_std_dB)
            sprintf('Mag. mean         %.4g dB', mean(mag, 'omitnan'))};
        if ~isempty(adev)
            [mn, i] = min(adev.ADEV_deg);
            txt{end + 1} = sprintf('ADEV minimum      %.3g deg at tau = %.3g s', mn, adev.Tau_s(i));
        end
        if isfield(settings, 'zero') && settings.zero, txt{end + 1} = 'Relative (zeroed) readings'; end
        cla(h.txt); axis(h.txt, 'off');
        text(h.txt, 0, 1, txt, 'VerticalAlignment', 'top', 'FontName', 'Consolas', 'FontSize', 9, ...
            'Interpreter', 'none');
        drawnow;
    end

    function files = saveLog()
        if ~exist(opt.OutDir, 'dir'), mkdir(opt.OutDir); end
        base = fullfile(opt.OutDir, "vvm_" + stamp());
        out = struct('f', settings.frequency, 'P', mag, 't', t, 'settings', settings, 'idn', idn, ...
            'phase_deg', ph, 'phaseUnwrapped_deg', phU, 'drift', rmfield(drift, {'phaseFit', 'magFit'}), 'adev', adev);
        save(base + ".mat", '-struct', 'out');
        writetable(table(t, seconds(t - t(1)), mag, ph, phU, 'VariableNames', ...
            {'Time', 'Elapsed_s', 'Mag_dB', 'Phase_deg', 'PhaseUnwrapped_deg'}), base + ".csv");
        files = [base + ".mat"; base + ".csv"; base + ".png"];
        if hasFig(), exportgraphics(fig, files(3), 'Resolution', 150); else, files(3) = []; end
        fprintf('Saved %s.mat / .csv / .png\n', base);
    end

    function closeFig()
        if hasFig(), delete(fig); end
    end
end

%% ===== local functions =====================================================
% Graphics callbacks and cleanup guards are local (not nested) functions on
% purpose: a nested-function handle held by a figure keeps the whole workspace,
% including the instrument connection, alive after an error.

function configureVVM(ff, opt)
ff.write("CALC:PAR:DEF " + opt.Parameter);
if ~isempty(opt.Frequency), ff.write(sprintf("SENS:FREQ:CENT %.10g", opt.Frequency)); end
if ~isempty(opt.Power),     ff.write(sprintf("SOUR:POW %g", opt.Power)); end
if ~isempty(opt.IFBW),      ff.write(sprintf("SENS:BWID %g", opt.IFBW)); end   % listed for VVM, unverified on this unit
if ~isempty(opt.Average),   ff.write(sprintf("SENS:AVER:COUN %d", opt.Average)); end
ff.opc(30);
e = ff.errors();
if ~isempty(e), warning('ff_vvm_logger:config', 'FieldFox reported: %s', strjoin(e, '; ')); end
end

function s = readVVMSettings(ff)
s.mode      = "VVM";
s.parameter = erase(ff.query("CALC:PAR:DEF?"), '"');
s.frequency = ff.queryNum("SENS:FREQ:CENT?");
s.power_dBm = queryNumSafe(ff, "SOUR:POW?");
s.ifbw      = queryNumSafe(ff, "SENS:BWID?");
s.average   = queryNumSafe(ff, "SENS:AVER:COUN?");
s.units     = "P = |" + s.parameter + "| in dB (not dBm), phase in degrees";
s.start     = datetime('now');
ff.errors();                                            % a refused query must not linger in the queue
end

function sync = startSync(ff)
% One fresh measurement per reading (INIT:IMM;*OPC?) if VVM accepts it,
% otherwise let it free-run and read its latest value.
ff.sweep(60);
e = ff.errors();
sync = isempty(e);
if ~sync
    fprintf('VVM did not accept a single-measurement trigger (%s): free-running instead.\n', strjoin(e, '; '));
    ff.write("INIT:CONT 1");
end
end

function c = unitsCheck(ff)
% FDATA should be magnitude in dB and phase in deg; SDATA is the raw complex
% value. Printing both once lets you confirm the units on the instrument.
fd = readValues(ff, "TRAC:DATA? FDATA", 2);
sd = readValues(ff, "TRAC:DATA? SDATA", 2);
z = complex(sd(1), sd(2));
c = [fd(1) fd(2) 20*log10(abs(z)) rad2deg(angle(z))];
fprintf('Units check: FDATA %.3f / %.2f   SDATA -> %.3f dB / %.2f deg\n', c);
e = ff.errors();
if ~isempty(e), fprintf('  (SDATA not available: %s)\n', strjoin(e, '; ')); c(3:4) = NaN; end
end

function zeroVVM(ff, ref)
% Relative mode: the present reading becomes the reference (0 dB, 0 or 180 deg).
if ref == "SHORT", ref = "SHOR"; end
ff.write("SENS:CORR:ZERO:REF " + ref);
ff.write("SENS:CORR:ZERO:STAT ON");
ff.opc(30);
e = ff.errors();
if ~isempty(e), warning('ff_vvm_logger:zero', 'Zeroing: %s', strjoin(e, '; ')); end
end

function [d, adev] = stability(t, mag, phU, f)
% Linear drift (fitlm) and Allan deviation (allanvar) of the phase.
th = hours(t - t(1));
mp = fitlm(th, phU); cp = coefCI(mp);
mm = fitlm(th, mag); cm = coefCI(mm);
d.phase_deg_per_h   = mp.Coefficients.Estimate(2);
d.phase_CI95        = cp(2, :);
d.delay_ps_per_h    = -d.phase_deg_per_h/(360*f)*1e12;     % phase = -360 f tau
d.phase_change_deg  = d.phase_deg_per_h*th(end);
d.phase_std_deg     = std(phU, 'omitnan');
d.phase_resid_std_deg = std(mp.Residuals.Raw, 'omitnan');
d.mag_dB_per_h      = mm.Coefficients.Estimate(2);
d.mag_CI95          = cm(2, :);
d.mag_std_dB        = std(mag, 'omitnan');
d.phaseFit          = @(x) predict(mp, x(:));
d.magFit            = @(x) predict(mm, x(:));

adev = [];
ts = seconds(t - t(1));
dt = median(diff(ts));
if numel(phU) < 8 || ~(dt > 0), return, end
tu = (0:floor(ts(end)/dt)).'*dt;                         % uniform time grid
pu = interp1(ts, phU, tu, 'linear');
pu = pu(isfinite(pu));
try
    [av, tau] = allanvar(pu, 'octave', 1/dt);            % Navigation Toolbox
catch
    [av, tau] = allanOverlapping(pu, 1/dt);              % same formula, if the toolbox is missing
end
adev = table(tau(:), sqrt(av(:)), floor(numel(pu)./(tau(:)/dt)), 'VariableNames', {'Tau_s', 'ADEV_deg', 'N'});
end

function [avar, tau] = allanOverlapping(x, fs)
% Overlapping Allan variance of samples x at rate fs, octave-spaced tau
% (the formula used by allanvar).
x = x(:);
m = 2.^(0:floor(log2((numel(x) - 1)/2))).';
theta = [0; cumsum(x)]/fs;                              % integral of x, starting at 0
n = numel(theta);
avar = zeros(size(m));
for i = 1:numel(m)
    k = m(i);
    d = theta(1 + 2*k:n) - 2*theta(1 + k:n - k) + theta(1:n - 2*k);
    avar(i) = sum(d.^2)/(2*(k/fs)^2*(n - 2*k));
end
tau = m/fs;
end

function [t, mag, ph, settings, idn] = loadLog(file)
S = load(file);
if ~all(isfield(S, {'f', 'P', 't', 'phase_deg'}))
    error('ff_vvm_logger:file', '%s is not a VVM log (needs f, P, t, phase_deg).', file);
end
t = S.t(:); mag = double(S.P(:)); ph = double(S.phase_deg(:));
settings = getfieldOr(S, 'settings', struct());
settings.frequency = S.f(1);
if ~isfield(settings, 'parameter'), settings.parameter = "S?"; end
idn = string(getfieldOr(S, 'idn', ""));
fprintf('Loaded %d readings from %s\n', numel(ph), file);
end

function [t, mag, ph, settings, idn] = demoLog(opt)
% Synthetic cable test: slow thermal drift + random walk + white noise; the
% phase is wrapped to +-180 deg like the instrument reports it. No instrument.
n = max(120, round(min(opt.Duration, 3600)/max(opt.Interval, 0.1)));
dt = max(opt.Interval, 0.1);
t = datetime('now') + seconds((0:n - 1).'*dt);
rng(11);
th = (0:n - 1).'*dt/3600;
phU = 177 + 30*th + cumsum(0.01*randn(n, 1)) + 0.05*randn(n, 1);
ph = mod(phU + 180, 360) - 180;
mag = -0.25 - 0.4*th + 0.003*randn(n, 1);
settings = struct('mode', "VVM (synthetic demo)", 'parameter', opt.Parameter, ...
    'frequency', defaultTo(opt.Frequency, 100e6), 'power_dBm', -10, 'ifbw', 1e3, 'average', 1, ...
    'units', "P = |S| in dB (not dBm), phase in degrees", 'start', t(1), 'interval', dt, ...
    'duration', n*dt, 'sync', true, 'zero', false, 'zeroRef', "OPEN");
idn = "Synthetic demo (no instrument)";
end

function v = readValues(ff, cmd, n)
% Numbers from a reply that may be ASCII ("a,b") or a binary block ("#..."):
% TRACe:DATA? follows FORM:DATA, which another program may have left on REAL.
c = ff.Client;
ff.write(cmd);
b = read(c, 1, "uint8");
if b == '#'
    nd = double(read(c, 1, "uint8")) - 48;
    len = NaN;
    if nd >= 1 && nd <= 9, len = str2double(char(read(c, nd, "uint8"))); end
    if ~(len >= 0), error('ff_vvm_logger:block', 'Malformed binary block from the FieldFox.'); end
    bytes = uint8(read(c, len, "uint8"));
    read(c, 1, "uint8");                                % trailing LF
    if len == 4*n, v = double(typecast(bytes(:).', 'single'));
    else, v = double(typecast(bytes(:).', 'double')); end
elseif b == 10                                          % empty reply
    v = [];
else
    v = sscanf(strtrim(erase([char(b) char(readline(c))], '"')), '%f,').';
end
v(abs(v) > 1e30) = NaN;                                 % 9.91E37 = "no data"
if numel(v) < n, v(end + 1:n) = NaN; end
end

function v = queryNumSafe(ff, cmd)
try
    v = ff.queryNum(cmd);
    if isempty(v), v = NaN; end
catch
    v = NaN;
end
end

function closeRequest(src, ~)
% Closing the window while measuring = "stop"; otherwise just close it.
if isappdata(src, 'busy') && getappdata(src, 'busy')
    setappdata(src, 'stop', true);
else
    delete(src);
end
end

function s = stopRequested(fig)
s = ~isempty(fig) && (~isvalid(fig) || (isappdata(fig, 'stop') && getappdata(fig, 'stop')));
end

function tNext = nextSlot(tNext, interval, now)
% Readings on a fixed time grid; after a stall skip to the next free slot (no burst).
tNext = tNext + interval;
if interval > 0 && tNext < now, tNext = ceil(now/interval)*interval; end
end

function g = instrumentGuard(ff)
g = onCleanup(@() closeInstrument(ff));
end

function closeInstrument(ff)
% Restore the FieldFox (backup .sta, mode, format, sweep) and drop the connection.
if ~isvalid(ff), return, end
try
    flush(ff.Client, "input");                          % drop half-read replies after an error
    ff.release();
catch err
    warning('FieldFox:release', 'Could not restore the FieldFox: %s', err.message);
end
delete(ff);
end

function g = hiddenFigureGuard(fig, visible)
g = onCleanup(@() deleteIfHidden(fig, visible));
end

function deleteIfHidden(fig, visible)
if strcmp(visible, 'off') && isvalid(fig), delete(fig); end
end

function g = busyGuard(fig)
% After an error or Ctrl+C while measuring, the window must close normally again.
g = onCleanup(@() notBusy(fig));
end

function notBusy(fig)
if isvalid(fig), setappdata(fig, 'busy', false); end
end

function ok = confirm(question, default)
r = strtrim(input(question + " [y/n] ", 's'));
if isempty(r), ok = default; else, ok = lower(r(1)) == 'y'; end
end

function [ip, opt] = parseArgs(ip, args, opt)
names = fieldnames(opt);
if (ischar(ip) || isstring(ip)) && any(strcmpi(names, string(ip)))
    args = [{ip} args]; ip = "";                        % called as ff_vvm_logger('File', x)
end
if mod(numel(args), 2), error('Options must be Name, Value pairs.'); end
for k = 1:2:numel(args)
    i = find(strcmpi(names, string(args{k})));
    if isempty(i), error('Unknown option "%s". Valid: %s', args{k}, strjoin(names, ', ')); end
    opt.(names{i}) = args{k + 1};
end
if isempty(ip) || strlength(string(ip)) == 0, ip = "10.90.36.200"; end
ip = string(ip);
end

function v = getfieldOr(S, name, default)
if isfield(S, name), v = S.(name); else, v = default; end
end

function v = defaultTo(v, d)
if isempty(v), v = d; end
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

function s = stamp()
s = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
end
