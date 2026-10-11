function R = ff_power_meter(ip, varargin)
%FF_POWER_METER  Channel-power logger using the FieldFox built-in power meter (option 310).
%
%   ff_power_meter                                   10.90.36.200, keep the CPM settings
%   ff_power_meter(ip, 'Center', 100e6, 'ChannelBW', 1e6, 'Duration', 600)
%   R = ff_power_meter('File', "cpm_20261008_101500.mat")   re-analyse a saved log
%   R = ff_power_meter('Demo', true)                         synthetic data, no instrument
%
% The FieldFox "Channel Power Meter" (CPM mode) integrates the power in one
% channel, Center +/- ChannelBW/2, with its own receiver - no power sensor
% needed. The input is Port 2 (the SA RF input; max +27 dBm, 50 V DC: use an
% attenuator / DC block with a LISN or a strong source).
% This tool reads that value every Interval seconds and shows
%   Strip chart   power vs time, plus a moving power average
%   Histogram     distribution of the readings with Normal and kernel densities
%                 fitted by fitdist (Statistics and Machine Learning Toolbox)
%   Statistics    power mean, median, spread, percentiles, fitted mu/sigma with
%                 95% confidence intervals (paramci)
%   SA check      at the end, one SA-mode measurement of the same channel
%                 (RMS detector, power-averaged sweeps, the CPM attenuation)
%                 integrated with bandpower
%                 (Signal Processing Toolbox) - a cross-check of the CPM value
% The FieldFox is restored exactly as it was (backup .sta) on every exit.
% Close the window to stop early.
%
% Name-Value options (defaults):
%   'Center'     []      channel centre (Hz); [] = keep the FieldFox CPM setting
%   'ChannelBW'  []      integration bandwidth (Hz), 10 Hz - 100 MHz; [] = keep
%   'Average'    []      CPM averages (1 = off); [] = keep
%   'Atten'      []      RF attenuation: CPM accepts only 10 or 30 dB; [] = keep
%   'Interval'   1       s between readings
%   'Duration'   60      s to log (Inf = until the window is closed)
%   'CompareSA'  true    SA-mode cross-check at the end
%   'SASweeps'   10      SA sweeps averaged for the cross-check
%   'Save'       true    save MAT, CSV and PNG in OutDir (cpm_<time>.*)
%   'OutDir'     pwd
%   'File'       ""      analyse a saved cpm_*.mat instead of measuring
%   'Demo'       false   synthetic readings instead of the instrument
%   'Visible'    'on'    'off' for headless runs (the PNG is still saved)
%   'Prompt'     true    ask before switching the instrument to CPM; false = never ask
%   'Port'       5025
%
% Output R: t (datetime), P_dBm, f (centre, Hz), settings, stats (table),
% dist (fitted NormalDistribution), sa (SA cross-check, [] if not done), idn, files.
% The MAT file uses the shared recording layout - f (1x1 centre), P (readings x 1,
% dBm), t, settings, idn - so the other tools can replay it.
%
% What the number means: the total power that reaches Port 2 inside the channel.
% With nothing connected you see the analyzer's own noise in that bandwidth
% (+10 dB per decade of ChannelBW). There is no CISPR quasi-peak detector on this
% unit, so use it for trends and comparisons (e.g. solar inverter OFF vs ON),
% not as a compliance result. The SA cross-check assumes a Gaussian RBW filter
% (noise bandwidth 1.065 x RBW), like fieldfox_live.

if nargin < 1, ip = ""; end
opt = struct('Center', [], 'ChannelBW', [], 'Average', [], 'Atten', [], 'Interval', 1, ...
    'Duration', 60, 'CompareSA', true, 'SASweeps', 10, 'Save', true, 'OutDir', pwd, ...
    'File', "", 'Demo', false, 'Visible', 'on', 'Prompt', true, 'Port', 5025);
[ip, opt] = parseArgs(ip, varargin, opt);
live = strlength(string(opt.File)) == 0 && ~opt.Demo;
if live && ~isfinite(opt.Duration) && strcmp(opt.Visible, 'off')
    error('ff_power_meter:duration', 'Duration = Inf needs a visible window (close it to stop).');
end
if ~isempty(opt.ChannelBW) && (opt.ChannelBW < 10 || opt.ChannelBW > 100e6)
    error('ff_power_meter:bw', 'ChannelBW must be between 10 Hz and 100 MHz in CPM mode.');
end
if ~isempty(opt.Atten) && ~ismember(opt.Atten, [10 30])
    error('ff_power_meter:atten', 'CPM mode accepts only 10 or 30 dB attenuation.');
end

t = NaT(0, 1); P = zeros(0, 1); settings = struct(); idn = ""; sa = [];
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
    if ~contains(upper(modes), "CPM")
        error('ff_power_meter:noCPM', ['This FieldFox lists no CPM mode (INST:CAT? -> %s).\n' ...
            'The built-in channel power meter needs option 310.'], modes);
    end
    m0 = ff.mode();
    if opt.Prompt && m0 ~= "CPM" && ~confirm(sprintf(['The FieldFox is in %s mode. Switch it to CPM ' ...
            'for this measurement (everything is restored afterwards)?'], m0), true)
        R = []; closeFig(); return
    end
    ff.takeControl();
    ff.setMode("CPM");
    configureCPM(ff, opt);
    settings = readCPMSettings(ff);
    settings.interval = opt.Interval; settings.duration = opt.Duration;
    settings.sync = startSync(ff);
    fprintf('CPM: %s, channel %s wide, averaging %s, atten %g dB\n', fmtHz(settings.center), ...
        fmtHz(settings.channelBW), onOff(settings.averageOn, settings.averageCount), settings.atten_dB);
    if hasFig(), setappdata(fig, 'busy', true); figBusy = busyGuard(fig); end %#ok<NASGU> closable after an error
    logLoop(ff);
    if opt.CompareSA && ~isempty(P) && ~stopRequested(fig)
        sa = compareSA(ff, settings.center, settings.channelBW, opt.SASweeps, settings.atten_dB);
    end
    if hasFig(), setappdata(fig, 'busy', false); end
    clear guard                                         % restores the FieldFox now
elseif strlength(string(opt.File)) > 0
    [t, P, settings, idn, sa] = loadLog(opt.File);
else
    [t, P, settings, idn] = demoLog(opt);
end

%% Analysis
if isempty(P), warning('ff_power_meter:empty', 'No readings.'); R = []; closeFig(); return, end
[stats, pd] = powerStats(P);
if ~isempty(sa)
    sa.cpmMean_dBm = stats.Value(stats.Quantity == "Mean (power)");
    fprintf('SA cross-check: SA %.2f dBm, CPM mean %.2f dBm (difference %+.2f dB)\n', ...
        sa.channelPower_dBm, sa.cpmMean_dBm, sa.cpmMean_dBm - sa.channelPower_dBm);
end
if isempty(fig)
    [fig, h] = makeFigure();
    figGuard = hiddenFigureGuard(fig, opt.Visible);     %#ok<NASGU> no hidden figure left on errors
end
showResults();

R = struct('t', t, 'P_dBm', P, 'f', settings.center, 'settings', settings, 'stats', stats, ...
    'dist', pd, 'sa', sa, 'idn', idn, 'files', strings(0, 1));
if opt.Save && (live || opt.Demo)
    R.files = saveLog();
elseif opt.Save && hasFig()
    png = fullfile(opt.OutDir, erase(string(getFileName(opt.File)), ".mat") + "_analysis.png");
    exportgraphics(fig, png, 'Resolution', 150); R.files = png;
    fprintf('Saved %s\n', png);
end
disp(stats);
if stopRequested(fig) || strcmp(opt.Visible, 'off'), closeFig(); end

%% ===== nested helpers (share the variables above) ==========================

    function logLoop(ff)
        tStart = tic; tPlot = tic; tHist = tic; k = 0; tNext = 0;
        fprintf('Logging every %g s for %g s (close the window to stop) ...\n', opt.Interval, opt.Duration);
        while ~stopRequested(fig) && toc(tStart) < opt.Duration
            while toc(tStart) < tNext && ~stopRequested(fig)          % wait, keep the window alive
                pause(min(0.05, tNext - toc(tStart)));
            end
            if stopRequested(fig) || toc(tStart) >= opt.Duration, break, end
            if settings.sync, ff.sweep(60); end
            v = readValues(ff, "TRAC?", 1);
            k = k + 1;
            t(k, 1) = datetime('now'); P(k, 1) = v(1);
            tNext = nextSlot(tNext, opt.Interval, toc(tStart));
            if hasFig() && toc(tPlot) > 0.5             % plotting must not slow the readings
                updateStrip();
                if toc(tHist) > 3, updateHist(); tHist = tic; end
                drawnow limitrate; tPlot = tic;
            end
        end
        fprintf('%d readings.\n', numel(P));
    end

    function ok = hasFig()
        ok = ~isempty(fig) && isvalid(fig);
    end

    function [fig, h] = makeFigure()
        fig = figure('Name', 'FieldFox channel power meter', 'NumberTitle', 'off', 'Color', 'w', ...
            'Position', [80 80 1200 760], 'Visible', opt.Visible, 'CloseRequestFcn', @closeRequest);
        tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        h.strip = nexttile(tl, [1 3]); hold(h.strip, 'on'); grid(h.strip, 'on'); box(h.strip, 'on');
        h.pts = plot(h.strip, NaT, NaN, '.-', 'Color', [0.3 0.5 1]);
        h.avg = plot(h.strip, NaT, NaN, 'k', 'LineWidth', 1.5);
        ylabel(h.strip, 'Channel power (dBm)');
        legend(h.strip, {'Reading', 'Moving power average (10)'}, 'Location', 'best', 'AutoUpdate', 'off');
        h.hist = nexttile(tl); h.sa = nexttile(tl); h.txt = nexttile(tl); axis(h.txt, 'off');
        title(h.hist, 'Distribution'); title(h.sa, 'SA-mode cross-check');
        axis(h.sa, 'off');
    end

    function updateStrip()
        if ~hasFig(), return, end
        set(h.pts, 'XData', t, 'YData', P);
        set(h.avg, 'XData', t, 'YData', 10*log10(movmean(10.^(P/10), 10, 'omitnan')));
        title(h.strip, sprintf('%s  |  channel %s  |  last %.2f dBm  (%d readings)', ...
            fmtHz(settings.center), fmtHz(settings.channelBW), P(end), numel(P)));
    end

    function updateHist()
        if ~hasFig(), return, end
        x = P(isfinite(P));
        cla(h.hist); hold(h.hist, 'on'); grid(h.hist, 'on'); box(h.hist, 'on');
        histogram(h.hist, x, 'Normalization', 'pdf', 'FaceColor', [0.6 0.75 1]);
        if numel(x) >= 3 && std(x) > 0
            pn = fitdist(x, 'Normal'); pk = fitdist(x, 'Kernel');
            xx = linspace(min(x) - 2*std(x), max(x) + 2*std(x), 200);
            plot(h.hist, xx, pdf(pn, xx), 'r', 'LineWidth', 1.5);
            plot(h.hist, xx, pdf(pk, xx), 'k--', 'LineWidth', 1.2);
            legend(h.hist, {'Readings', sprintf('Normal \\mu=%.2f \\sigma=%.2f', pn.mu, pn.sigma), 'Kernel'}, ...
                'Location', 'best');
        end
        xlabel(h.hist, 'Channel power (dBm)'); ylabel(h.hist, 'Probability density');
        title(h.hist, 'Distribution');
    end

    function showResults()
        if ~hasFig(), return, end
        updateStrip(); updateHist();
        rows = compose('%-22s %10.6g %s', stats.Quantity, stats.Value, stats.Unit);
        txt = [{sprintf('Channel %s, BW %s', fmtHz(settings.center), fmtHz(settings.channelBW))}; rows];
        if ~isempty(sa)
            axis(h.sa, 'on'); cla(h.sa); hold(h.sa, 'on'); grid(h.sa, 'on'); box(h.sa, 'on');
            [sc, su] = freqUnit(max(abs(sa.f)));
            plot(h.sa, sa.f/sc, sa.trace_dBm, 'b');
            xregion(h.sa, sa.channel(1)/sc, sa.channel(2)/sc, 'FaceColor', [1 0.85 0.3], 'FaceAlpha', 0.3);
            xlabel(h.sa, "Frequency (" + su + ")"); ylabel(h.sa, 'dBm (RMS detector)');
            title(h.sa, sprintf('SA: %.2f dBm in channel (RBW %s)', sa.channelPower_dBm, fmtHz(sa.rbw)));
            txt = [txt; {''; sprintf('SA bandpower %10.2f dBm', sa.channelPower_dBm); ...
                sprintf('CPM mean     %10.2f dBm', sa.cpmMean_dBm); ...
                sprintf('CPM - SA     %+10.2f dB', sa.cpmMean_dBm - sa.channelPower_dBm); ...
                sprintf('SA 99%% OBW   %s', fmtHz(sa.obw_Hz))}];
        end
        cla(h.txt); axis(h.txt, 'off');
        text(h.txt, 0, 1, txt, 'VerticalAlignment', 'top', 'FontName', 'Consolas', 'FontSize', 9, ...
            'Interpreter', 'none');
        drawnow;
    end

    function files = saveLog()
        if ~exist(opt.OutDir, 'dir'), mkdir(opt.OutDir); end
        base = fullfile(opt.OutDir, "cpm_" + stamp());
        out = struct('f', settings.center, 'P', P, 't', t, 'settings', settings, 'idn', idn, ...
            'stats', stats, 'sa', sa);
        save(base + ".mat", '-struct', 'out');
        writetable(table(t, seconds(t - t(1)), P, 'VariableNames', {'Time', 'Elapsed_s', 'Power_dBm'}), base + ".csv");
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

function configureCPM(ff, opt)
if ~isempty(opt.Center),    ff.write(sprintf("SENS:FREQ %.10g", opt.Center)); end
if ~isempty(opt.ChannelBW), ff.write(sprintf("SENS:FREQ:SPAN %.10g", opt.ChannelBW)); end
if ~isempty(opt.Average)
    ff.write(sprintf("SENS:AVER:ENAB %d", opt.Average > 1));
    if opt.Average > 1, ff.write(sprintf("SENS:AVER:COUN %d", opt.Average)); end
end
if ~isempty(opt.Atten), ff.write(sprintf("SENS:POW:RF:ATT %g", opt.Atten)); end
ff.write("SENS:AMPL:UNIT DBM");                    % readings in dBm whatever the screen shows
ff.opc(30);
e = ff.errors();
if ~isempty(e), warning('ff_power_meter:config', 'FieldFox reported: %s', strjoin(e, '; ')); end
end

function s = readCPMSettings(ff)
s.mode         = "CPM";
s.center       = ff.queryNum("SENS:FREQ?");
s.channelBW    = ff.queryNum("SENS:FREQ:SPAN?");
s.averageOn    = queryNumSafe(ff, "SENS:AVER:ENAB?") == 1;
s.averageCount = queryNumSafe(ff, "SENS:AVER:COUN?");
s.atten_dB     = queryNumSafe(ff, "SENS:POW:RF:ATT?");
s.unit         = "dBm";
s.start        = datetime('now');
end

function sync = startSync(ff)
% One fresh measurement per reading (INIT:IMM;*OPC?) if CPM accepts it,
% otherwise let the meter free-run and read its latest value.
ff.sweep(60);
e = ff.errors();
sync = isempty(e);
if ~sync
    fprintf('CPM did not accept a single-measurement trigger (%s): free-running instead.\n', strjoin(e, '; '));
    ff.write("INIT:CONT 1");
end
end

function sa = compareSA(ff, fc, bw, nSweeps, atten)
% Same channel in SA mode: RMS detector, power-averaged sweeps, bandpower over
% the channel. Same attenuation as CPM, so that even a noise-only channel
% compares like with like (the attenuation sets the analyzer noise floor).
if bw < 1e3
    fprintf('SA cross-check skipped: channel narrower than 1 kHz.\n'); sa = []; return
end
rbw = min(max(10^floor(log10(bw/20)), 10), 1e6);       % >= 20 RBWs across the channel
span = 2*bw;
pts = min(10001, max(401, ceil(span/rbw) + 1));         % points no further apart than the RBW
fprintf('SA cross-check: span %s, RBW %s, %d points, %d sweeps ...\n', fmtHz(span), fmtHz(rbw), pts, nSweeps);
if isfinite(atten), att = {'Atten', atten}; else, att = {'AttenAuto', true}; end
ff.saConfigure(att{:}, Center=fc, Span=span, RBW=rbw, Points=pts, Detector="AVER", Preamp=false, ...
    Average=1, TraceType="CLRW");
s = ff.saSettings();
Y = zeros(nSweeps, s.points);
for k = 1:nSweeps
    ff.sweep(120);
    Y(k, :) = ff.saTrace();
end
sxx = mean(10.^(Y/10), 1);                              % mW per trace point, power average
lims = [fc - bw/2, fc + bw/2];
pmw = bandpower(sxx/(1.065*s.rbw), s.f, lims, 'psd');   % PSD = point power / noise bandwidth
[obwHz, flo, fhi] = obw(sxx(:), s.f(:), s.rbw, lims, 99);
sa = struct('f', s.f, 'trace_dBm', 10*log10(sxx), 'rbw', s.rbw, 'channel', lims, ...
    'channelPower_dBm', 10*log10(pmw), 'obw_Hz', obwHz, 'obwEdges_Hz', [flo fhi], ...
    'sweeps', nSweeps, 'detector', s.detector, 'atten_dB', s.atten, 'cpmMean_dBm', NaN);
end

function [T, pd] = powerStats(P)
x = P(isfinite(P));
pm = 10*log10(mean(10.^(x/10)));                        % average in mW, not in dB
q = prctile(x, [10 50 90]);
if numel(x) >= 3 && std(x) > 0
    pd = fitdist(x, 'Normal'); ci = paramci(pd);
else
    pd = []; ci = nan(2);
end
name = ["Readings"; "Mean (power)"; "Mean of dB values"; "Median"; "Std deviation"; ...
    "Minimum"; "Maximum"; "10th percentile"; "90th percentile"; "Max - mean (power)"; ...
    "Fitted mu"; "mu 95% CI low"; "mu 95% CI high"; "Fitted sigma"];
val = [numel(x); pm; mean(x); q(2); std(x); min(x); max(x); q(1); q(3); max(x) - pm; ...
    getMu(pd); ci(1, 1); ci(2, 1); getSigma(pd)];
unit = ["-"; "dBm"; "dBm"; "dBm"; "dB"; "dBm"; "dBm"; "dBm"; "dBm"; "dB"; "dBm"; "dBm"; "dBm"; "dB"];
T = table(name, val, unit, 'VariableNames', {'Quantity', 'Value', 'Unit'});
end

function v = getMu(pd), if isempty(pd), v = NaN; else, v = pd.mu; end, end
function v = getSigma(pd), if isempty(pd), v = NaN; else, v = pd.sigma; end, end

function [t, P, settings, idn, sa] = loadLog(file)
S = load(file);
need = {'f', 'P', 't'};
if ~all(isfield(S, need)) || size(S.P, 2) ~= 1
    error('ff_power_meter:file', '%s is not a power-meter log (needs f, P with one column, t).', file);
end
t = S.t(:); P = double(S.P(:));
settings = getfieldOr(S, 'settings', struct());
if ~isfield(settings, 'center'), settings.center = S.f(1); end
if ~isfield(settings, 'channelBW'), settings.channelBW = NaN; end
idn = string(getfieldOr(S, 'idn', ""));
sa = getfieldOr(S, 'sa', []);
fprintf('Loaded %d readings from %s\n', numel(P), file);
end

function [t, P, settings, idn] = demoLog(opt)
% Synthetic log: a -40 dBm channel with slow fading and short bursts (e.g. an
% inverter switching on and off). No instrument involved.
n = max(30, round(min(opt.Duration, 3600)/opt.Interval));
t = datetime('now') + seconds((0:n - 1).'*opt.Interval);
rng(7);
P = -40 + 0.3*randn(n, 1) + 0.8*sin(2*pi*(0:n - 1).'/max(n/3, 1));
burst = rand(n, 1) < 0.05;
P(burst) = P(burst) + 6 + 2*rand(nnz(burst), 1);
settings = struct('mode', "CPM (synthetic demo)", 'center', defaultTo(opt.Center, 100e6), ...
    'channelBW', defaultTo(opt.ChannelBW, 1e6), 'averageOn', false, 'averageCount', 1, ...
    'atten_dB', NaN, 'unit', "dBm", 'start', t(1), 'interval', opt.Interval, ...
    'duration', opt.Duration, 'sync', true);
idn = "Synthetic demo (no instrument)";
end

function v = readValues(ff, cmd, n)
% Numbers from a reply that may be ASCII ("a,b") or a binary block ("#..."):
% the FieldFox applies FORM:DATA to some queries but not to others.
c = ff.Client;
ff.write(cmd);
b = read(c, 1, "uint8");
if b == '#'
    nd = double(read(c, 1, "uint8")) - 48;
    len = NaN;
    if nd >= 1 && nd <= 9, len = str2double(char(read(c, nd, "uint8"))); end
    if ~(len >= 0), error('ff_power_meter:block', 'Malformed binary block from the FieldFox.'); end
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

function ok = confirm(question, default)
r = strtrim(input(question + " [y/n] ", 's'));
if isempty(r), ok = default; else, ok = lower(r(1)) == 'y'; end
end

function [ip, opt] = parseArgs(ip, args, opt)
names = fieldnames(opt);
if (ischar(ip) || isstring(ip)) && any(strcmpi(names, string(ip)))
    args = [{ip} args]; ip = "";                        % called as ff_power_meter('File', x)
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

function n = getFileName(file)
[~, n, e] = fileparts(string(file)); n = n + e;
end

function s = onOff(on, n)
if on, s = sprintf('ON (%g)', n); else, s = 'OFF'; end
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
