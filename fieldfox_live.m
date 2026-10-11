function fieldfox_live(ip, varargin)
%FIELDFOX_LIVE  Real-time spectrum monitor for a Keysight FieldFox, in MATLAB.
%
%   fieldfox_live                          connect to 10.90.36.200
%   fieldfox_live("10.90.36.200")
%   fieldfox_live(ip, Name, Value, ...)
%
% Reads the trace the FieldFox is measuring in SA mode over the LAN (SCPI
% socket, port 5025). It does NOT change the analyzer's settings: set
% frequency, span, RBW and reference level on the FieldFox and the app follows.
%
% Window:
%   Spectrum       live (blue), average (black), max hold (red), baseline (green)
%   Waterfall      last sweeps, newest at the top
%   Signals table  peaks of the average trace, labelled with the band they fall in
%   Measurements   noise floor, channel power and 99% occupied bandwidth
%                  (channel = "Ch center / Ch width" boxes; empty = whole span)
%
% Baseline (EMI check, e.g. solar inverter OFF vs ON):
%   1. Inverter OFF -> press "Set baseline" (averages the next BaselineSweeps sweeps)
%   2. Inverter ON  -> bins where the average rises more than Threshold dB above
%      the baseline AND above the highest level seen during the baseline are
%      marked with red x and listed as new emissions.
%
% Buttons: Pause, Max hold, Reset, Set baseline, Clear baseline, Record, Snapshot.
% Close the window to stop.
%
% Name-Value options (defaults):
%   'Record'         false   start recording at once (fieldfox_rec_<time>.mat)
%   'AutoBaseline'   false   take the baseline from the first sweeps automatically
%   'BaselineFile'   ""      use the average of a saved recording as baseline
%   'BaselineSweeps' 20      sweeps averaged into the baseline
%   'AvgCount'       10      sweeps in the running average
%   'Threshold'      6       dB above baseline to call an emission "new"
%   'History'        200     sweeps kept for the waterfall
%   'MaxSweeps'      Inf     stop by itself after this many sweeps
%   'SnapshotAtEnd'  false   save a snapshot when the app stops
%   'OutDir'         pwd     where recordings and snapshots go
%   'Port'           5025
%
% A recording (.mat) contains f (Hz), P (dBm, sweeps x points), t (datetime),
% settings and idn; a long one continues in a new file every 256 MB of data.
% Needs only base MATLAB: no Keysight software required.

if nargin < 1 || strlength(string(ip)) == 0, ip = "10.90.36.200"; end
opt = struct('Port', 5025, 'Record', false, 'AutoBaseline', false, 'BaselineFile', "", ...
    'BaselineSweeps', 20, 'AvgCount', 10, 'Threshold', 6, 'History', 200, ...
    'MaxSweeps', Inf, 'SnapshotAtEnd', false, 'OutDir', pwd, 'Visible', 'on');
for k = 1:2:numel(varargin)
    if ~isfield(opt, varargin{k}), error('Unknown option "%s".', varargin{k}); end
    opt.(varargin{k}) = varargin{k + 1};
end

%% Connect
fprintf('Connecting to %s:%d ...\n', ip, opt.Port);
ff = tcpclient(ip, opt.Port, 'Timeout', 10, 'ConnectTimeout', 5);
configureTerminator(ff, "LF");
idn = query(ff, '*IDN?');
idnParts = split(idn, ",");
model = idnParts(min(2, numel(idnParts)));
fprintf('Connected: %s\n', idn);

instMode = erase(query(ff, 'INST:SEL?'), '"');
if instMode ~= "SA"
    error('The FieldFox is in %s mode. Switch it to SA (Spectrum Analyzer) mode first.', instMode);
end

% The trace is read as ASCII. If another program left the format on binary,
% switch to ASCII and put it back when the app closes.
oldFormat = query(ff, 'FORM:DATA?');
if ~startsWith(oldFormat, "ASC")
    writeline(ff, 'FORM:DATA ASC,0');
end
restore = onCleanup(@() restoreFormat(ff, oldFormat));

%% State
s        = readSettings(ff);
[fScale, fUnit] = freqUnit(s.stop);
H        = nan(opt.History, s.points);   % waterfall history, newest in row 1
nH       = 0;
maxHold  = -inf(1, s.points);
avgTr    = nan(1, s.points);
yPrev    = [];
base     = emptyBaseline();
rec      = [];
recChunk = 2^25;                         % values per recording file (256 MB), then a new file
nSweeps  = 0;
rate     = NaN;
running  = true;
loopActive = false;
forceSettings = false;
newBins  = false(1, s.points);

%% Window
fig = figure('Name', "FieldFox live - " + model + " @ " + ip, 'NumberTitle', 'off', ...
    'Color', 'w', 'Position', [60 60 1400 820], 'Visible', opt.Visible, ...
    'CloseRequestFcn', @onClose);

axS = axes(fig, 'Position', [0.05 0.50 0.62 0.44]);
hold(axS, 'on'); grid(axS, 'on'); box(axS, 'on');
hBase = plot(axS, NaN, NaN, '--', 'Color', [0 0.6 0], 'LineWidth', 1.2);
hMax  = plot(axS, NaN, NaN, 'r', 'LineWidth', 0.8);
hLive = plot(axS, NaN, NaN, 'Color', [0.3 0.5 1], 'LineWidth', 0.8);
hAvg  = plot(axS, NaN, NaN, 'k', 'LineWidth', 1.4);
hNew  = plot(axS, NaN, NaN, 'rx', 'MarkerSize', 9, 'LineWidth', 2);
hPk   = plot(axS, NaN, NaN, 'kv', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
hCh   = xregion(axS, s.start/fScale, s.stop/fScale, 'FaceColor', [1 0.85 0.3], 'FaceAlpha', 0.12);
ylabel(axS, 'Power (dBm)');
legend(axS, [hLive hAvg hMax hBase hNew], {'Live', 'Average', 'Max hold', 'Baseline', 'New emission'}, ...
    'Location', 'northeast', 'AutoUpdate', 'off');

axW = axes(fig, 'Position', [0.05 0.16 0.62 0.26]);
hW = imagesc(axW, [0 1], [1 opt.History], H);
colormap(axW, turbo); cb = colorbar(axW); cb.Label.String = 'dBm';
ylabel(axW, 'Sweeps ago (newest at top)'); set(axW, 'YDir', 'reverse');
linkaxes([axS axW], 'x');

tbl = uitable(fig, 'Units', 'normalized', 'Position', [0.70 0.42 0.285 0.52], ...
    'ColumnName', {'Freq', 'dBm', 'vs base', 'Possible source'}, ...
    'ColumnWidth', {85, 50, 70, 175}, 'RowName', [], 'Data', cell(0, 4));
measTxt = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.70 0.16 0.285 0.24], 'BackgroundColor', [0.96 0.96 0.96], ...
    'HorizontalAlignment', 'left', 'FontName', 'Consolas', 'FontSize', 9, 'String', '');

uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.70 0.105 0.09 0.03], ...
    'String', "Ch center (" + fUnit + ")", 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
edCenter = uicontrol(fig, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.79 0.105 0.06 0.035], 'String', '');
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.855 0.105 0.07 0.03], ...
    'String', "Ch width (" + fUnit + ")", 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
edWidth = uicontrol(fig, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.925 0.105 0.06 0.035], 'String', '');

btn = @(style, str, x, cb) uicontrol(fig, 'Style', style, 'String', str, 'Units', 'normalized', ...
    'Position', [x 0.03 0.08 0.05], 'FontSize', 10, 'Callback', cb);
bPause = btn('togglebutton', 'Pause',          0.05, []);
bMax   = btn('togglebutton', 'Max hold',       0.14, @(~, ~) set(hMax, 'Visible', onOff(bMax.Value)));
btn('pushbutton',   'Reset',          0.23, @(~, ~) resetTraces());
btn('pushbutton',   'Set baseline',   0.32, @(~, ~) startBaseline());
btn('pushbutton',   'Clear baseline', 0.41, @(~, ~) clearBaseline());
bRec   = btn('togglebutton', 'Record',         0.50, @(~, ~) toggleRecord());
btn('pushbutton',   'Snapshot',       0.59, @(~, ~) snapshot());
status = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.70 0.03 0.285 0.05], ...
    'BackgroundColor', 'w', 'HorizontalAlignment', 'left', 'FontSize', 9, 'String', '');
bMax.Value = 1;

try
    applySettingsToPlots();
    if strlength(string(opt.BaselineFile)) > 0, loadBaseline(opt.BaselineFile); end
    if opt.AutoBaseline, startBaseline(); end
    if opt.Record, bRec.Value = 1; toggleRecord(); end
catch err
    shutdown();                       % no window (holding the connection) left behind
    rethrow(err);
end

%% Live loop
tLast = tic; tSettings = tic; tSlow = tic;
loopActive = true;
try
    liveLoop();
catch err
    loopActive = false;
    shutdown();
    rethrow(err);
end
loopActive = false;
if isvalid(fig)
    updateMeasurements(); drawnow;
    if opt.SnapshotAtEnd, snapshot(); end
end
shutdown();
fprintf('Stopped after %d sweeps.\n', nSweeps);

%% ===== nested helpers (share the variables above) ==========================

    function liveLoop()
        % Ctrl+C is not caught by try/catch, and the window callbacks keep this
        % workspace (and the connection) alive: this guard shuts down instead.
        interrupted = onCleanup(@() shutdownIfInterrupted());
        while running
            if bPause.Value
                pause(0.1);
                continue
            end

            % Re-read settings once a second so changes made on the FieldFox show up
            if forceSettings || toc(tSettings) > 1
                checkSettings();
                tSettings = tic; forceSettings = false;
            end

            y = sscanf(char(query(ff, 'TRAC:DATA?')), '%f,').';
            if numel(y) ~= s.points
                forceSettings = true;
                continue
            end
            if isequal(y, yPrev)          % same sweep as last time: wait for a new one
                pause(0.01);
                continue
            end
            yPrev = y;
            processSweep(y);

            if toc(tSlow) > 0.5           % table and measurements twice a second
                updateMeasurements();
                drawnow;                  % also lets button clicks run
                tSlow = tic;
            end
            dt = toc(tLast); tLast = tic;
            if isnan(rate), rate = 1/dt; else, rate = 0.9*rate + 0.1/dt; end
            updateStatus();
            drawnow limitrate;

            if nSweeps >= opt.MaxSweeps, running = false; end
        end
        loopActive = false;
    end

    function shutdownIfInterrupted()
        if loopActive, loopActive = false; shutdown(); end
    end

    function processSweep(y)
        nSweeps = nSweeps + 1;
        H = [y; H(1:end - 1, :)];
        nH = min(nH + 1, opt.History);
        maxHold = max(maxHold, y);
        nA = min(opt.AvgCount, nH);
        avgTr = 10*log10(mean(10.^(H(1:nA, :)/10), 1));    % average in mW, not dB

        if base.collecting
            base.acc(end + 1, :) = y;
            if size(base.acc, 1) >= opt.BaselineSweeps, finishBaseline(); end
        end
        newBins = false(1, s.points);
        if base.ready
            newBins = avgTr > base.avg + opt.Threshold & avgTr > base.max;
        end
        if ~isempty(rec)
            rec.P(end + 1, :) = y;
            rec.t(end + 1, 1) = datetime('now');
            if numel(rec.P) >= recChunk   % memory stays bounded: save, continue in a new file
                saveRecording(); rec = newRecording(s);
                fprintf('Recording continues in %s\n', rec.file);
            end
        end

        set(hLive, 'YData', y);
        set(hAvg,  'YData', avgTr);
        set(hMax,  'YData', maxHold);
        set(hNew,  'XData', s.f(newBins)/fScale, 'YData', avgTr(newBins));
        set(hW,    'CData', H, 'AlphaData', ~isnan(H));
    end

    function updateMeasurements()
        if nH == 0, return, end
        lin = 10.^(avgTr/10);                              % mW per trace point
        floorDb = median(avgTr, 'omitnan');

        % Peaks of the average trace, 10 dB above the noise floor
        isPk = islocalmax(avgTr, 'MinProminence', 6) & avgTr > floorDb + 10;
        pk = find(isPk);
        [~, order] = sort(avgTr(pk), 'descend');
        pk = pk(order(1:min(15, end)));
        set(hPk, 'XData', s.f(pk)/fScale, 'YData', avgTr(pk));
        rows = cell(numel(pk), 4);
        for i = 1:numel(pk)
            rows{i, 1} = sprintf('%.4f', s.f(pk(i))/fScale);
            rows{i, 2} = sprintf('%.1f', avgTr(pk(i)));
            if base.ready
                d = avgTr(pk(i)) - base.avg(pk(i));
                rows{i, 3} = sprintf('%+.1f%s', d, repmat(' NEW', 1, newBins(pk(i))));
            else
                rows{i, 3} = '';
            end
            rows{i, 4} = bandName(s.f(pk(i)), s.binStep/2);
        end
        tbl.Data = rows;
        tbl.ColumnName{1} = char("Freq (" + fUnit + ")");

        % Channel power and 99% occupied bandwidth inside the channel.
        % Channel power = sum of trace points * (point spacing / noise bandwidth
        % of the RBW filter, ~1.065 x RBW for a Gaussian filter). It is an
        % estimate: exact only with an average/RMS detector and points <= RBW apart.
        [c1, c2] = channelEdges();
        inCh = s.f >= c1 & s.f <= c2;
        chMw = sum(lin(inCh)) * s.binStep / (1.065*s.rbw);
        cumP = cumsum(lin(inCh)) / sum(lin(inCh));
        fCh = s.f(inCh);
        obwHz = NaN;
        if nnz(inCh) > 2
            lo = fCh(find(cumP >= 0.005, 1));
            hi = fCh(find(cumP >= 0.995, 1));
            obwHz = hi - lo;
        end
        set(hCh, 'Value', [c1 c2]/fScale, 'Visible', onOff(c1 > s.start || c2 < s.stop));

        txt = {
            sprintf('Noise floor (median)  %8.1f dBm', floorDb)
            sprintf('Peak                  %8.1f dBm @ %.4f %s', max(avgTr), s.f(find(avgTr == max(avgTr), 1))/fScale, fUnit)
            sprintf('Channel %s - %s', fmtHz(c1), fmtHz(c2))
            sprintf('  Channel power  ~ %8.1f dBm', 10*log10(chMw))
            sprintf('  99%% occ. BW      %s', fmtHz(obwHz))
            sprintf('RBW %s  VBW %s  %d pts', fmtHz(s.rbw), fmtHz(s.vbw), s.points)
            sprintf('Point spacing %s', fmtHz(s.binStep))
            };
        if base.ready
            txt{end + 1} = sprintf('New emissions: %d points (> +%g dB)', nnz(newBins), opt.Threshold);
            if any(newBins)
                [dMax, iMax] = max((avgTr - base.avg) .* newBins);
                txt{end + 1} = sprintf('  strongest +%.1f dB @ %.4f %s', dMax, s.f(iMax)/fScale, fUnit);
            end
        end
        measTxt.String = txt;
    end

    function [c1, c2] = channelEdges()
        c  = str2double(edCenter.String) * fScale;
        bw = str2double(edWidth.String) * fScale;
        if isfinite(c) && isfinite(bw) && bw > 0
            c1 = max(c - bw/2, s.start); c2 = min(c + bw/2, s.stop);
        else
            c1 = s.start; c2 = s.stop;
        end
    end

    function updateStatus()
        txt = sprintf('%.1f sweeps/s   #%d', rate, nSweeps);
        if base.collecting
            txt = sprintf('%s   BASELINE %d/%d', txt, size(base.acc, 1), opt.BaselineSweeps);
        end
        if ~isempty(rec), txt = sprintf('%s   REC %d', txt, size(rec.P, 1)); end
        status.String = txt;
        title(axS, sprintf('%s live  |  %s - %s  |  RBW %s', model, fmtHz(s.start), fmtHz(s.stop), fmtHz(s.rbw)));
    end

    function checkSettings()
        sNew = readSettings(ff);
        if sNew.start ~= s.start || sNew.stop ~= s.stop || sNew.points ~= s.points
            if ~isempty(rec), saveRecording(); rec = newRecording(sNew); end
            if base.ready || base.collecting
                fprintf('Frequency settings changed: baseline cleared.\n');
                base = emptyBaseline();
            end
            s = sNew;
            H = nan(opt.History, s.points); nH = 0; yPrev = [];
            resetTraces();
            applySettingsToPlots();
        else
            s.rbw = sNew.rbw; s.vbw = sNew.vbw; s.ref = sNew.ref; s.pdiv = sNew.pdiv;
            applyLevels();
        end
    end

    function applySettingsToPlots()
        [fScale, fUnit] = freqUnit(s.stop);
        x = s.f/fScale;
        n = nan(1, s.points);
        set([hLive hAvg hMax], 'XData', x, 'YData', n);
        if base.ready, set(hBase, 'XData', x, 'YData', base.avg); else, set(hBase, 'XData', NaN, 'YData', NaN); end
        set(hW, 'XData', x([1 end]), 'YData', [1 opt.History], 'CData', H, 'AlphaData', ~isnan(H));
        xlabel(axW, "Frequency (" + fUnit + ")");
        xlim(axS, x([1 end]));
        applyLevels();
    end

    function applyLevels()
        if isfinite(s.ref) && isfinite(s.pdiv)
            lim = [s.ref - 10*s.pdiv, s.ref];       % same scale as the FieldFox screen
            ylim(axS, lim); clim(axW, lim);
        end
    end

    function resetTraces()
        maxHold = -inf(1, s.points);
        set(hMax, 'YData', nan(1, s.points));
    end

    %% Baseline
    function b = emptyBaseline()
        b = struct('collecting', false, 'ready', false, 'acc', [], 'avg', [], 'max', [], 'f', []);
    end

    function startBaseline()
        base = emptyBaseline();
        base.collecting = true;
        base.acc = zeros(0, s.points);
        set(hBase, 'XData', NaN, 'YData', NaN);
        fprintf('Baseline: averaging the next %d sweeps...\n', opt.BaselineSweeps);
    end

    function finishBaseline()
        base.avg = 10*log10(mean(10.^(base.acc/10), 1));
        base.max = max(base.acc, [], 1);
        base.f = s.f;
        base.collecting = false; base.ready = true;
        set(hBase, 'XData', s.f/fScale, 'YData', base.avg);
        fprintf('Baseline ready (%d sweeps).\n', size(base.acc, 1));
    end

    function clearBaseline()
        base = emptyBaseline();
        newBins = false(1, s.points);
        set(hBase, 'XData', NaN, 'YData', NaN);
        set(hNew, 'XData', NaN, 'YData', NaN);
    end

    function loadBaseline(file)
        S = load(file, 'f', 'P');
        if numel(S.f) ~= s.points || abs(S.f(1) - s.start) > s.binStep/2 || abs(S.f(end) - s.stop) > s.binStep/2
            warning('Baseline file %s has a different frequency axis than the FieldFox now; not used.', file);
            return
        end
        base = emptyBaseline();
        base.acc = S.P;
        finishBaseline();
        fprintf('Baseline loaded from %s.\n', file);
    end

    %% Recording and snapshots
    function toggleRecord()
        if bRec.Value
            rec = newRecording(s);
            fprintf('Recording to %s\n', rec.file);
        elseif ~isempty(rec)
            saveRecording();
            rec = [];
        end
    end

    function r = newRecording(st)
        r.f = st.f; r.P = zeros(0, st.points); r.t = NaT(0, 1);
        r.settings = rmfield(st, 'f'); r.idn = idn;
        r.file = fullfile(opt.OutDir, "fieldfox_rec_" + stamp() + ".mat");
    end

    function saveRecording()
        if isempty(rec) || isempty(rec.P), return, end
        out = struct('f', rec.f, 'P', rec.P, 't', rec.t, 'settings', rec.settings, 'idn', rec.idn);
        save(rec.file, '-struct', 'out');
        fprintf('Saved %d sweeps to %s\n', size(rec.P, 1), rec.file);
    end

    function snapshot()
        b = fullfile(opt.OutDir, "fieldfox_snap_" + stamp());
        try
            imwrite(frame2im(getframe(fig)), b + ".png");
        catch
            exportgraphics(axS, b + ".png", 'Resolution', 150);
        end
        T = table(s.f.', get(hLive, 'YData').', avgTr.', maxHold.', ...
            'VariableNames', {'Freq_Hz', 'Live_dBm', 'Average_dBm', 'MaxHold_dBm'});
        if base.ready
            T.Baseline_dBm = base.avg.';
            T.NewEmission = newBins.';
        end
        writetable(T, b + ".csv");
        fprintf('Snapshot saved: %s.png / .csv\n', b);
    end

    %% Closing
    function onClose(~, ~)
        running = false;
        if ~loopActive, shutdown(); end
    end

    function shutdown()
        if ~isempty(rec), saveRecording(); rec = []; end
        restore = [];                 % runs the onCleanup: trace format back
        if isvalid(fig), delete(fig); end
    end
end

%% ===== local functions =====================================================

function r = query(ff, cmd)
writeline(ff, cmd);
r = strtrim(readline(ff));
end

function s = readSettings(ff)
s.start  = str2double(query(ff, 'SENS:FREQ:STAR?'));
s.stop   = str2double(query(ff, 'SENS:FREQ:STOP?'));
s.points = str2double(query(ff, 'SENS:SWE:POIN?'));
s.rbw    = str2double(query(ff, 'SENS:BAND:RES?'));
s.vbw    = str2double(query(ff, 'SENS:BAND:VID?'));
s.ref    = str2double(query(ff, 'DISP:WIND:TRAC:Y:RLEV?'));
s.pdiv   = str2double(query(ff, 'DISP:WIND:TRAC:Y:PDIV?'));
s.f      = linspace(s.start, s.stop, s.points);
s.binStep = (s.stop - s.start) / max(s.points - 1, 1);
end

function restoreFormat(ff, oldFormat)
if isvalid(ff) && ~startsWith(oldFormat, "ASC")
    writeline(ff, "FORM:DATA " + oldFormat);
end
end

function name = bandName(fHz, tol)
% Common services in Thailand. A match only says the frequency fits a band.
persistent bands
if isempty(bands)
    bands = {
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
        'Wi-Fi 5 GHz',                5150e6,  5850e6
        'Wi-Fi 6E 6 GHz',             5925e6,  6425e6 };
end
hit = fHz >= [bands{:, 2}] - tol & fHz <= [bands{:, 3}] + tol;
if any(hit), name = strjoin(bands(hit, 1), ' / '); else, name = ''; end
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

function s = onOff(v)
if v, s = 'on'; else, s = 'off'; end
end

function s = stamp()
s = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
end
