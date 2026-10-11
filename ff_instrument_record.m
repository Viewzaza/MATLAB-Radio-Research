function R = ff_instrument_record(ip, varargin)
%FF_INSTRUMENT_RECORD  Long unattended spectrum recording on the FieldFox itself (SA Record/Playback, option 236).
%
%   R = ff_instrument_record                                 10.90.36.200: record 10 min, current SA settings
%   R = ff_instrument_record(ip, 'Start', 150e3, 'Stop', 30e6, 'RBW', 10e3, 'Duration', 8*3600, 'Interval', 10)
%   R = ff_instrument_record(ip, 'Session', "latest")        download + parse a session already on the FieldFox
%   R = ff_instrument_record('parse', folder)                offline: a session folder copied from the FieldFox
%
% The FieldFox writes the sweeps to its own drive, so the recording carries on
% if the PC sleeps or the LAN drops: start it, leave it, and fetch it later
% with 'Session'. What the tool does:
%   1. backs up the instrument state (.sta) and sets the SA - only the
%      settings you give, the others stay as they are on the FieldFox
%   2. creates a new session (RECP:SESS:NEW -> "AutoSessionNN") and records the
%      raw SA sweep every Interval s for Duration s. The FieldFox also pauses by
%      itself after Duration (RECP:CONF:RTLS), whatever happens to the PC.
%   3. shows the progress and a live trace (close the window to stop early),
%      then stops and closes the session
%   4. downloads the session files (MMEM:CAT? / MMEM:DATA?, or FTP as a
%      fallback) into OutDir\RecordPlayback\GPSA\<session>, same folder layout
%   5. parses them into the shared recording layout (f, P, t, settings, idn)
%      and draws a waterfall, the average / max-hold / min spectra and the
%      power in the span vs time (bandpower, Signal Processing Toolbox).
%      With a known RBW the average spectrum is also passed to ff_signal_id.
% The FieldFox is restored (backup .sta, mode, sweep, data format) on every
% exit, including errors and Ctrl+C. The session stays on the FieldFox: delete
% it there (File > Manage Files) once you have it on the PC.
%
% Name-Value options (defaults):
%   'Duration'   600     s to record (max 1e6 s; Inf = until the window is closed)
%   'Interval'   0       s between recorded sweeps (0 = every sweep), max 10000
%   'Start' 'Stop' 'Center' 'Span' 'RBW' 'Points' 'RefLevel' 'Atten' 'Preamp' 'Detector'
%                []      SA settings; [] = keep the FieldFox setting
%   'Storage'    "INT"   "INT" | "SD" | "USB": the drive the FieldFox records to
%   'Session'    ""      no recording: download + parse this session ("latest" = newest)
%   'Download'   true    download and parse after recording
%   'Transfer'   "auto"  "scpi" (MMEM:DATA?) | "ftp" (anonymous FTP) | "auto" (SCPI, then FTP)
%   'Preview'    5       s between live-trace previews while recording (0 = none)
%   'Identify'   true    run ff_signal_id on the average spectrum (needs the RBW)
%   'Save'       true    MAT + PNG in OutDir (ffrec_<session>.*)
%   'OutDir'     pwd
%   'Visible'    'on'    'off' for headless runs (progress goes to the console)
%   'Prompt'     true    ask before recording; false = never ask
%   'Port'       5025
% With 'parse', also 'RBW' (Hz, if you know it) and 'Session' (which session in
% the folder; default the newest) are used.
%
% Output R: f (1xN Hz), P (sweeps x N, dBm), t (datetime), recNo, settings, idn,
% session, remote (path on the FieldFox), folder (local copy), level (table of
% power vs time), signals (ff_signal_id table or []), files.
%
% Caveats:
%  - The sweeps are not gap-free (no real-time option 350): each record is one
%    whole sweep, then the analyzer sweeps on. Short bursts can fall in between.
%  - "RAWM" records the raw detector trace in dBm; trace averaging and max hold
%    on the screen are not recorded. There is no quasi-peak detector (no option
%    361): for inverter EMI this is an indicative pre-scan, not a CISPR result.
%  - Times come from the FieldFox clock, which may differ from the PC clock;
%    after a recording the offset is printed and saved (settings.clockOffset_s).
%  - The RBW is not stored in the session files: offline, give 'RBW' to get
%    absolute power vs time; otherwise the mean level per point is shown.
%  - The session file layout was decoded from real files (A.12.45). The
%    download location on the FieldFox and the FTP paths are not documented by
%    Keysight; the tool searches the usual places and reports what it used.

if nargin < 1, ip = ""; end
opt = struct('Duration', 600, 'Interval', 0, 'Start', [], 'Stop', [], 'Center', [], 'Span', [], ...
    'RBW', [], 'Points', [], 'RefLevel', [], 'Atten', [], 'Preamp', [], 'Detector', [], ...
    'Storage', "INT", 'Session', "", 'Download', true, 'Transfer', "auto", 'Preview', 5, ...
    'Identify', true, 'Save', true, 'OutDir', pwd, 'Visible', 'on', 'Prompt', true, 'Port', 5025);
task = "record"; folder = "";
if (ischar(ip) || isstring(ip)) && strcmpi(string(ip), "parse")
    if isempty(varargin), error('ff_instrument_record:parse', 'Usage: ff_instrument_record(''parse'', folder)'); end
    folder = string(varargin{1}); varargin(1) = []; ip = ""; task = "parse";
end
[ip, opt] = parseArgs(ip, varargin, opt);
opt.Storage = upper(string(opt.Storage)); opt.Transfer = lower(string(opt.Transfer));
opt.Session = string(opt.Session);
if task == "record" && strlength(opt.Session) > 0, task = "session"; end
if ~ismember(opt.Storage, ["INT" "SD" "USB"]), error('ff_instrument_record:storage', 'Storage must be INT, SD or USB.'); end
if ~ismember(opt.Transfer, ["auto" "scpi" "ftp"]), error('ff_instrument_record:transfer', 'Transfer must be auto, scpi or ftp.'); end
if task == "record"
    if ~(opt.Duration > 0) || (isfinite(opt.Duration) && opt.Duration > 1e6)
        error('ff_instrument_record:duration', 'Duration must be > 0 and at most 1e6 s (or Inf).');
    end
    if ~isfinite(opt.Duration) && strcmp(opt.Visible, 'off')
        error('ff_instrument_record:duration', 'Duration = Inf needs a visible window (close it to stop).');
    end
    if opt.Interval < 0 || opt.Interval > 1e4, error('ff_instrument_record:interval', 'Interval must be 0 - 10000 s.'); end
end

idn = ""; session = ""; remote = ""; sessDir = ""; files = strings(0, 1);
live = struct();                                        % SA settings read while recording
fig = gobjects(0); h = struct(); figGuard = []; %#ok<SETNU> set by showSession()

if task == "parse"
    [sessDir, session] = findLocalSession(folder, opt.Session);
else
    fprintf('Connecting to %s:%d ...\n', ip, opt.Port);
    ff = FieldFox(ip, opt.Port);
    job = containers.Map({'recording', 'cdir', 'get'}, {false, "", ""});
    guard = instrumentGuard(ff, job);                   %#ok<NASGU> stop the session + restore on any exit
    idn = ff.IDN;
    if ~ff.hasOption("236")
        error('ff_instrument_record:no236', ['This FieldFox has no option 236 (*OPT? -> %s): ' ...
            'SA Record/Playback is not available.'], strjoin(ff.Options, ','));
    end
    if task == "record"
        ok = recordSession(ff, job);
        if ~ok, R = []; clear guard; return, end
    else
        [remote, session] = findRemoteSession(ff, opt.Storage, opt.Session);
        if sessionOpen(ff)
            warning('ff_instrument_record:open', ['A Record/Playback session is open on the FieldFox ' ...
                '(maybe still recording): the files may be incomplete. Stop it on the FieldFox and run again ' ...
                'for the complete data.']);
        end
    end
    if strlength(session) == 0
        warning('ff_instrument_record:nosession', ['Could not find the session folder on the FieldFox. ' ...
            'Copy RecordPlayback\\GPSA\\<session> with the FieldFox File Manager and use ''parse''.']);
    elseif task == "session" || opt.Download
        [sessDir, files] = downloadSession(ff, job, ip, opt, remote, session);
    end
    clear guard                                         % restores the FieldFox now
end

R = struct('f', [], 'P', [], 't', NaT(0, 1), 'recNo', [], 'settings', struct(), 'idn', idn, ...
    'session', session, 'remote', remote, 'folder', sessDir, 'level', [], 'signals', [], 'files', files);
if strlength(sessDir) == 0, return, end

%% Parse and show
[f, P, t, recNo, hdr] = readSession(sessDir);
rbw = getfieldOr(live, 'rbw', NaN);                    % read back from the FieldFox when recording
if ~isfinite(rbw), rbw = defaultTo(opt.RBW, NaN); end
settings = struct('start', f(1), 'stop', f(end), 'center', hdr.center_Hz, 'span', hdr.span_Hz, ...
    'points', hdr.points, 'rbw', rbw, 'source', "FieldFox Record/Playback (RAWM)", 'session', session, ...
    'sessionStart', hdr.time, 'remote', remote, 'unit', "dBm");
if ~isempty(fieldnames(live)), settings = mergeStruct(rmfield(live, 'f'), settings); end
if isfield(live, 'recordStart') && ~isempty(t)
    settings.clockOffset_s = seconds(t(1) - live.recordStart);   % + up to one sweep of latency
    fprintf('FieldFox clock minus PC clock: about %+.0f s (the times below are FieldFox times)\n', settings.clockOffset_s);
end
R.f = f; R.P = P; R.t = t; R.recNo = recNo; R.settings = settings;
fprintf('%s: %d sweeps, %d points, %s - %s', session, size(P, 1), hdr.points, fmtHz(f(1)), fmtHz(f(end)));
if ~isempty(t), fprintf(', %s to %s', string(t(1)), string(t(end))); end
fprintf('\n');
if isempty(P), warning('ff_instrument_record:empty', 'The session holds no sweeps.'); return, end

R.level = levelVsTime(f, P, t, rbw);
if opt.Identify && isfinite(rbw) && exist('ff_signal_id', 'file')
    try
        R.signals = ff_signal_id(f, 10*log10(mean(10.^(P/10), 1)), rbw, ...
            'Detector', string(getfieldOr(live, 'detector', "auto")));   % detector read back while recording
    catch err
        fprintf('ff_signal_id skipped: %s\n', err.message);
    end
end
showSession();
if opt.Save
    if ~exist(opt.OutDir, 'dir'), mkdir(opt.OutDir); end
    base = fullfile(opt.OutDir, "ffrec_" + matlab.lang.makeValidName(session));
    out = struct('f', f, 'P', P, 't', t, 'settings', settings, 'idn', idn, 'recNo', recNo);
    save(base + ".mat", '-struct', 'out');
    exportgraphics(fig, base + ".png", 'Resolution', 150);
    R.files = [R.files; base + ".mat"; base + ".png"];
    fprintf('Saved %s.mat / .png\n', base);
end
if ~isempty(R.signals), disp(R.signals); end
if strcmp(opt.Visible, 'off'), delete(fig); end

%% ===== nested helpers (share the variables above) ==========================

    function ok = recordSession(ff, job)
        ok = false;
        if opt.Prompt && ~confirm(sprintf(['Record %s on the FieldFox (%s storage), one sweep every %s? ' ...
                'The instrument is restored afterwards.'], fmtDur(opt.Duration), opt.Storage, ...
                intervalText(opt.Interval)), true)
            return
        end
        ff.takeControl();                               % backup .sta + single sweep
        ff.setMode("SA");
        sa = {'Start', 'Stop', 'Center', 'Span', 'RBW', 'Points', 'RefLevel', 'Atten', 'Preamp', 'Detector'};
        args = {};
        for i = 1:numel(sa)
            if ~isempty(opt.(sa{i})), args = [args, sa(i), {opt.(sa{i})}]; end %#ok<AGROW>
        end
        if ~isempty(args), ff.saConfigure(args{:}); end
        live = ff.saSettings();
        live.sweepTime = ff.sweep(120);
        step = max(opt.Interval, live.sweepTime);
        fprintf('SA %s - %s, %d points, RBW %s, detector %s; one sweep takes %.2f s\n', fmtHz(live.start), ...
            fmtHz(live.stop), live.points, fmtHz(live.rbw), live.detector, live.sweepTime);
        if sessionOpen(ff)
            error('ff_instrument_record:open', ['A Record/Playback session is already open on the FieldFox. ' ...
                'Stop and close it there (Record/Playback > Stop, Close Session) and run again.']);
        end
        before = listSessions(ff, opt.Storage);
        ff.write("RECP:SESS:SDEV " + opt.Storage);
        ff.write("RECP:SESS:NEW");
        ff.opc(60);
        job('recording') = true;                        % from here the guard stops/closes the session
        throwOnErrors(ff, 'Creating the session');
        [remote, session] = newSession(listSessions(ff, opt.Storage), before);
        trl = queryNumSafe(ff, "RECP:SESS:TRL?");       % trace record limit of the session
        ff.errors();
        ff.write("RECP:CONF:RSO RAWM");
        ff.write(sprintf("RECP:CONF:RTIN %g", opt.Interval));
        ff.write(sprintf("RECP:CONF:RTLS %g", finiteOr(opt.Duration, 0)));
        ff.opc(30);
        throwOnErrors(ff, 'Configuring the recording');
        nEst = finiteOr(floor(opt.Duration/step) + 1, Inf);
        if isfinite(trl) && trl > 0 && nEst > trl
            warning('ff_instrument_record:limit', ['About %d sweeps expected but the session holds at most %d: ' ...
                'the FieldFox stops recording after about %s. Use a longer Interval.'], nEst, trl, fmtDur(trl*step));
        end
        live.interval = opt.Interval; live.duration = opt.Duration; live.storage = opt.Storage;
        live.traceLimit = trl; live.recordStart = datetime('now');
        ff.write("INIT:CONT 1");                        % the recorder takes the sweeps as they run
        ff.write("RECP:ACT:REC");
        ff.opc(30);
        throwOnErrors(ff, 'Starting the recording');
        if ~sessionOpen(ff), warning('ff_instrument_record:state', 'The FieldFox reports no open session after REC.'); end
        fprintf('Recording to %s (%s) for %s ...\n', session, remote, fmtDur(opt.Duration));
        progressLoop(ff, step);
        ff.write("RECP:ACT:STOP"); ff.opc(60);
        ff.write("RECP:SESS:CLOS"); ff.opc(60);
        job('recording') = false;                       %#ok<NASGU>
        e = ff.errors();
        if ~isempty(e), warning('ff_instrument_record:stop', 'FieldFox reported: %s', strjoin(e, '; ')); end
        live.recordStop = datetime('now');
        fprintf('Recording stopped after %s; about %d sweeps.\n', fmtDur(seconds(live.recordStop - live.recordStart)), ...
            floor(seconds(live.recordStop - live.recordStart)/step) + 1);
        ok = true;
    end

    function progressLoop(ff, step)
        if strcmp(opt.Visible, 'on'), makeProgressFigure(); setappdata(fig, 'busy', true); end
        progressGuard = figureGuard(fig); %#ok<NASGU> the window also goes after an error or Ctrl+C
        tStart = tic; tPrev = -Inf; tLog = 0; tCheck = 0;
        logEvery = min(max(10, finiteOr(opt.Duration/10, 60)), 600);
        while toc(tStart) < opt.Duration && ~stopRequested(fig)
            pause(0.2);
            el = toc(tStart);
            if opt.Preview > 0 && el - tPrev >= opt.Preview && hasFig()
                y = ff.saTrace(); tPrev = el;
                if numel(y) == numel(live.f), set(h.live, 'YData', y); end
            end
            if el >= tCheck + 30                        % is the session still open?
                tCheck = el;
                if ~sessionOpen(ff), warning('ff_instrument_record:state', 'The session was closed on the FieldFox.'); break, end
            end
            if hasFig()
                frac = min(el/opt.Duration, 1);
                set(h.bar, 'XData', [0 frac frac 0]);
                h.txt.String = sprintf('%s   %s of %s   about %d sweeps', session, fmtDur(el), ...
                    fmtDur(opt.Duration), floor(el/step) + 1);
                drawnow limitrate
            end
            if el >= tLog + logEvery
                tLog = el;
                fprintf('  %s / %s elapsed, about %d sweeps\n', fmtDur(el), fmtDur(opt.Duration), floor(el/step) + 1);
            end
        end
        if hasFig(), delete(fig); end
        fig = gobjects(0);
    end

    function makeProgressFigure()
        fig = figure('Name', 'FieldFox recording', 'NumberTitle', 'off', 'Color', 'w', ...
            'Position', [80 80 1000 560], 'Visible', opt.Visible, 'CloseRequestFcn', @closeRequest);
        tl = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        ax = nexttile(tl); axis(ax, [0 1 0 1]); axis(ax, 'off');
        h.bar = patch(ax, [0 0 0 0], [0.35 0.35 0.65 0.65], [0.3 0.6 1], 'EdgeColor', 'none');
        rectangle(ax, 'Position', [0 0.35 1 0.3]);
        h.txt = text(ax, 0, 0.9, '', 'FontSize', 11, 'Interpreter', 'none');
        text(ax, 0, 0.05, 'Close this window to stop the recording early.', 'Color', [0.4 0.4 0.4]);
        ax = nexttile(tl, [3 1]); grid(ax, 'on'); box(ax, 'on'); hold(ax, 'on');
        [h.sc, su] = freqUnit(max(abs(live.f)));
        h.live = plot(ax, live.f/h.sc, nan(size(live.f)), 'Color', [0 0.3 0.8]);
        xlabel(ax, "Frequency (" + su + ")"); ylabel(ax, 'dBm');
        title(ax, 'Live trace (what is being recorded)');
    end

    function showSession()
        fig = figure('Name', "FieldFox recording " + session, 'NumberTitle', 'off', 'Color', 'w', ...
            'Position', [60 60 1250 820], 'Visible', opt.Visible);
        figGuard = hiddenFigureGuard(fig, opt.Visible);   % no hidden figure left on errors
        tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        [sc, su] = freqUnit(max(abs(f)));
        el = seconds(t - t(1));
        [tsc, tu] = timeUnit(max(el));
        ax = nexttile(tl, [1 2]);
        imagesc(ax, f/sc, el/tsc, P);
        colormap(ax, turbo); cb = colorbar(ax); cb.Label.String = 'dBm';
        lim = prctile(P(isfinite(P)), [1 99.9]);
        if numel(lim) == 2 && lim(2) > lim(1), clim(ax, lim); end
        xlabel(ax, "Frequency (" + su + ")");
        ylabel(ax, "Time since " + string(t(1)) + " (" + tu + ")");
        title(ax, sprintf('%s: %d sweeps (waterfall, first sweep at the top)', session, size(P, 1)), 'Interpreter', 'none');
        ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
        plot(ax, f/sc, 10*log10(mean(10.^(P/10), 1)), 'k', 'LineWidth', 1.2);
        plot(ax, f/sc, max(P, [], 1), 'r');
        plot(ax, f/sc, min(P, [], 1), 'Color', [0.3 0.6 1]);
        legend(ax, {'Average (power)', 'Max hold', 'Min hold'}, 'Location', 'best');
        xlabel(ax, "Frequency (" + su + ")"); ylabel(ax, 'dBm'); title(ax, 'Spectra over the whole recording');
        ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
        plot(ax, R.level.Time, R.level.(2), '.-', 'Color', [0 0.45 0.74]);
        plot(ax, R.level.Time, R.level.Peak_dBm, '.', 'Color', [0.85 0.33 0.1]);
        ylabel(ax, 'dBm');
        legend(ax, {strrep(R.level.Properties.VariableNames{2}, '_', ' '), 'Highest point'}, 'Location', 'best');
        if isfinite(rbw)
            title(ax, sprintf('Power in the span vs time (bandpower, RBW %s)', fmtHz(rbw)));
        else
            title(ax, 'Mean level per point vs time (RBW unknown)');
        end
        drawnow;
    end

    function ok = hasFig()
        ok = ~isempty(fig) && isvalid(fig);
    end
end

%% ===== local functions =====================================================
% Graphics callbacks and cleanup guards are local (not nested) functions on
% purpose: a nested-function handle held by a figure keeps the whole workspace,
% including the instrument connection, alive after an error.

function g = instrumentGuard(ff, job)
g = onCleanup(@() closeInstrument(ff, job));
end

function closeInstrument(ff, job)
% Stop/close our session if it is still open, go back to the start-up folder,
% restore the FieldFox (backup .sta, mode, format, sweep) and disconnect.
if ~isvalid(ff), return, end
try
    flush(ff.Client, "input");                          % drop half-read replies after an error
    if job('recording')
        ff.write("RECP:ACT:STOP"); ff.opc(60);
        ff.write("RECP:SESS:CLOS"); ff.opc(60);
        fprintf('Recording stopped; what was recorded stays on the FieldFox.\n');
    end
    if strlength(job('cdir')) > 0, ff.write("MMEM:CDIR """ + job('cdir') + """"); end
    ff.release();
catch err
    warning('FieldFox:release', 'Could not restore the FieldFox: %s', err.message);
end
delete(ff);
end

function g = figureGuard(fig)
g = onCleanup(@() deleteIfValid(fig));
end

function deleteIfValid(fig)
if ~isempty(fig) && isvalid(fig), delete(fig); end
end

function g = hiddenFigureGuard(fig, visible)
g = onCleanup(@() deleteIfHidden(fig, visible));
end

function deleteIfHidden(fig, visible)
if strcmp(visible, 'off') && isvalid(fig), delete(fig); end
end

function closeRequest(src, ~)
% Closing the window while recording = "stop"; otherwise just close it.
if isappdata(src, 'busy') && getappdata(src, 'busy')
    setappdata(src, 'stop', true);
else
    delete(src);
end
end

function s = stopRequested(fig)
s = ~isempty(fig) && (~isvalid(fig) || (isappdata(fig, 'stop') && getappdata(fig, 'stop')));
end

function throwOnErrors(ff, what)
e = ff.errors();
if ~isempty(e), error('ff_instrument_record:scpi', '%s: the FieldFox reported %s', what, strjoin(e, '; ')); end
end

function tf = sessionOpen(ff)
% STAT:OPER:SAM:COND? = 64 while a Record/Playback session is open (SA mode only).
v = queryNumSafe(ff, "STAT:OPER:SAM:COND?");
ff.errors();
tf = isfinite(v) && bitand(uint32(v), 64) > 0;
end

%% ---- files on the FieldFox -------------------------------------------------

function roots = gpsaRoots(storage)
% Where the FieldFox keeps sessions. Keysight: "USERDATA\RecordPlayback\GPSA";
% "\UserData" is an old alias of "[INTERNAL]:". Both spellings are tried.
switch storage
    case "SD",  d = "[SDCARD]:";
    case "USB", d = "[USBDISK]:";
    otherwise,  d = "[INTERNAL]:";
end
roots = d + ["\RecordPlayback\GPSA", "\USERDATA\RecordPlayback\GPSA"];
end

function S = listSessions(ff, storage)
% All session folders on the drive: struct array root, name.
S = struct('root', {}, 'name', {});
for root = gpsaRoots(storage)
    [E, ok] = listDir(ff, root);
    if ~ok, continue, end
    for e = E(:).'
        if e.isDir, S(end + 1) = struct('root', root, 'name', e.name); end %#ok<AGROW>
    end
end
end

function [root, name] = newSession(after, before)
% The session that RECP:SESS:NEW created = the folder that was not there before.
key = @(S) lower(string({S.root}) + "\" + string({S.name}));
if isempty(after), root = ""; name = ""; return, end
isNew = ~ismember(key(after), key(before));
if ~any(isNew)
    warning('ff_instrument_record:newsession', 'Could not see which session folder is new; using the newest.');
    isNew = true(size(after));
end
c = after(isNew);
[~, i] = max(sessionRank(string({c.name})));
root = c(i).root; name = c(i).name;
end

function [root, name] = findRemoteSession(ff, storage, want)
S = listSessions(ff, storage);
root = ""; name = "";
if isempty(S), return, end
names = string({S.name});
if strcmpi(want, "latest")
    [~, i] = max(sessionRank(names));
else
    i = find(strcmpi(names, want), 1);
    if isempty(i)
        error('ff_instrument_record:session', 'No session "%s" on the FieldFox. Sessions: %s', want, strjoin(names, ', '));
    end
end
root = S(i).root; name = S(i).name;
fprintf('Session %s\\%s\n', root, name);
end

function r = sessionRank(names)
% "AutoSession12" ranks by its number; other names after them, alphabetically.
n = str2double(regexprep(names, '^AutoSession\s*(\d+)$', '$1', 'ignorecase'));
[~, ~, alpha] = unique(lower(names));
r = n(:);
r(isnan(r)) = 1e6 + alpha(isnan(r));
end

function [E, ok] = listDir(ff, path)
% MMEM:CAT? "<folder>" -> struct array name, isDir. ok = false if the folder is not there.
r = ff.query("MMEM:CAT? """ + path + """", 30);
ok = isempty(ff.errors());
E = parseCatalog(r);
end

function E = parseCatalog(r)
% A.12.45 replies <used>,<free>,"name,,size","[Folder],,",... (checked on the
% unit; an empty folder gives 0,0,""). Folders are shown as [Name]; a plain
% SCPI "name,DIR,size" entry is also accepted. Without the <used>,<free>
% header the reply is read as a plain comma-separated list of names.
r = strtrim(string(r));
E = struct('name', {}, 'isDir', {});
if ~isempty(regexp(r, '^[+-]?\d+\s*,\s*[+-]?\d+', 'once'))
    tok = regexp(r, '"([^"]*)"', 'tokens');
    for k = 1:numel(tok)
        p = strtrim(split(string(tok{k}{1}), ","));
        n = p(1);
        if strlength(n) == 0, continue, end
        br = startsWith(n, "[") && endsWith(n, "]");
        if br, n = extractBetween(n, 2, strlength(n) - 1); end
        E(end + 1) = struct('name', n, 'isDir', br || (numel(p) > 1 && upper(p(2)) == "DIR")); %#ok<AGROW>
    end
    return
end
names = strtrim(split(erase(r, '"'), ","));
for n = names(strlength(names) > 0).'
    br = startsWith(n, "[") && endsWith(n, "]");
    if br, n = extractBetween(n, 2, strlength(n) - 1); end %#ok<FXSET>
    E(end + 1) = struct('name', n, 'isDir', br || ~contains(n, ".")); %#ok<AGROW>
end
end

function [sessDir, files] = downloadSession(ff, job, ip, opt, root, name)
% Copy <root>\<name> into OutDir\RecordPlayback\GPSA\<name>, same layout.
sessDir = string(fullfile(opt.OutDir, "RecordPlayback", "GPSA", name));
if opt.Transfer ~= "ftp"
    try
        fprintf('Downloading %s\\%s with MMEM:DATA? ...\n', root, name);
        files = getTree(ff, job, root + "\" + name, sessDir);
        fprintf('Downloaded %d files (%.1f kB) to %s\n', numel(files), totalKB(files), sessDir);
        return
    catch err
        if opt.Transfer == "scpi", rethrow(err); end
        fprintf('SCPI download failed (%s); trying FTP ...\n', err.message);
        flush(ff.Client, "input"); ff.errors();
    end
end
files = ftpDownload(ip, opt.Storage, name, fileparts(sessDir));
fprintf('Downloaded %d files by FTP to %s\n', numel(files), sessDir);
end

function files = getTree(ff, job, remote, localDir)
[E, ok] = listDir(ff, remote);
if ~ok, error('ff_instrument_record:list', 'Cannot list %s on the FieldFox.', remote); end
if ~exist(localDir, 'dir'), mkdir(localDir); end
files = strings(0, 1);
for e = E(:).'
    if e.isDir
        files = [files; getTree(ff, job, remote + "\" + e.name, fullfile(localDir, e.name))]; %#ok<AGROW>
    else
        bytes = getRemoteFile(ff, job, remote, e.name);
        out = string(fullfile(localDir, e.name));
        fid = fopen(out, 'w'); fwrite(fid, bytes, 'uint8'); fclose(fid);
        files(end + 1, 1) = out; %#ok<AGROW>
    end
end
end

function bytes = getRemoteFile(ff, job, folder, name)
% MMEM:DATA? with the full path; if the FieldFox refuses that, change into the
% folder (MMEM:CDIR) and ask for the bare name. The start-up folder is put back
% before release(), which reloads its backup .sta from there.
first = job('get') == "";                               % short timeouts until a transfer has worked
tmo = 10 + 50*~first;
if job('get') ~= "cdir"
    [bytes, ok] = tryGet(ff, folder + "\" + name, tmo);
    if ok, job('get') = "full"; return, end %#ok<NASGU>
    if ~first, error('ff_instrument_record:get', 'MMEM:DATA? failed for %s\\%s.', folder, name); end
    fprintf('  MMEM:DATA? with a full path was refused; trying MMEM:CDIR + file name.\n');
    if strlength(job('cdir')) == 0, job('cdir') = erase(ff.query("MMEM:CDIR?"), '"'); end
    job('get') = "cdir"; %#ok<NASGU>
end
ff.write("MMEM:CDIR """ + folder + """");
[bytes, ok] = tryGet(ff, name, tmo);
if ~ok, error('ff_instrument_record:get', 'MMEM:DATA? failed for %s\\%s.', folder, name); end
end

function [bytes, ok] = tryGet(ff, name, timeout)
% MMEM:DATA? "<name>" -> file bytes. The same transfer as FieldFox.getFile,
% with a timeout of our own: a refused name gets no reply at all.
c = ff.Client; old = c.Timeout; c.Timeout = timeout;
wid = 'transportlib:client:ReadbinblockWarning';        % what a read timeout gives (no error)
ws = warning('off', wid);
restore = onCleanup(@() restoreRead(c, old, ws));
try
    ff.write("MMEM:DATA? """ + name + """");
    lastwarn('');
    bytes = uint8(readbinblock(c, "uint8")); bytes = bytes(:);
    [~, w] = lastwarn;
    ok = ~strcmp(w, wid);
    if ok, read(c, 1, "uint8"); ok = isempty(ff.errors()); end   % the trailing LF
catch
    bytes = uint8([]); ok = false;
end
if ~ok, flush(c, "input"); ff.errors(); end
end

function files = ftpDownload(ip, storage, name, localGPSA)
% Fallback: anonymous FTP (MATLAB ftp + mget, which copies folders recursively).
% The FTP folder names on the FieldFox are not documented - several are tried.
f = ftp(ip);
c = onCleanup(@() close(f));
switch storage
    case "SD",  tops = ["SDCARD", "Storage Card", "SD"];
    case "USB", tops = ["USBDISK", "Hard Disk", "USB"];
    otherwise,  tops = ["", "USERDATA", "UserData", "INTERNAL"];
end
for top = tops
    p = strjoin([top "RecordPlayback" "GPSA"], "/");
    if top == "", p = "RecordPlayback/GPSA"; end
    try
        cd(f, "/"); cd(f, p);
    catch
        continue
    end
    if ~exist(localGPSA, 'dir'), mkdir(localGPSA); end
    files = string(mget(f, name, localGPSA));
    files = files(isfile(files));
    return
end
error('ff_instrument_record:ftp', 'FTP: no RecordPlayback/GPSA folder found on %s.', ip);
end

function restoreRead(c, timeout, ws)
warning(ws);
if isvalid(c), c.Timeout = timeout; end
end

function s = intervalText(x)
if x == 0, s = 'sweep (as fast as it sweeps)'; else, s = sprintf('%g s', x); end
end

function kb = totalKB(files)
kb = 0;
for k = 1:numel(files), d = dir(files(k)); kb = kb + d.bytes/1024; end
end

%% ---- session files on the PC (layout as in analyze_fieldfox.m) ------------

function [sessDir, name] = findLocalSession(folder, want)
folder = string(folder);
if ~isfolder(folder), error('ff_instrument_record:parse', 'Folder not found: %s', folder); end
if ~isempty(dir(fullfile(folder, 'dSTA*')))
    sessDir = folder;
else
    d = dir(fullfile(folder, '**', 'dSTA*'));
    d = d([d.isdir]);
    if isempty(d), error('ff_instrument_record:parse', 'No Record/Playback session (dSTA* folder) under %s.', folder); end
    cand = unique(string({d.folder}));
    [~, names] = fileparts(cand);
    if strlength(want) > 0 && ~strcmpi(want, "latest")
        i = find(strcmpi(names, want), 1);
        if isempty(i), error('ff_instrument_record:parse', 'No session "%s" in %s.', want, folder); end
    else
        [~, i] = max(sessionRank(names));
        if numel(cand) > 1, fprintf('Sessions found: %s (using %s)\n', strjoin(names, ', '), names(i)); end
    end
    sessDir = cand(i);
end
[~, name] = fileparts(sessDir);
end

function [f, P, t, recNo, hdr] = readSession(sessDir)
sta = dir(fullfile(sessDir, 'dSTA*', '*.STA'));
if isempty(sta), error('ff_instrument_record:parse', 'No dSTA*\\*.STA file in %s.', sessDir); end
hdr = readRecordHeader(fullfile(sta(1).folder, sta(1).name));
if numel(sta) > 1
    warning('ff_instrument_record:states', ['The session has %d setting states (settings changed while ' ...
        'recording); only the first one is used for the frequency axis.'], numel(sta));
end
[P, t, recNo] = readRecordTraces(dir(fullfile(sessDir, 'dTRC*', '*.TRC')), hdr.points);
f = linspace(hdr.center_Hz - hdr.span_Hz/2, hdr.center_Hz + hdr.span_Hz/2, hdr.points);
end

function h = readRecordHeader(staFile)
% Byte layout of the recording's .STA file (little-endian):
%   offset  8  uint16  number of points
%   offset 16  uint64  timestamp, .NET ticks (100 ns since 1-Jan-0001)
%   offset 24  double  center frequency (Hz)
%   offset 32  double  span (Hz)
b = readBytes(staFile);
if numel(b) < 40, error('ff_instrument_record:parse', '%s is too short (%d bytes).', staFile, numel(b)); end
h.points    = double(typecast(b(9:10), 'uint16'));
h.time      = datetime(typecast(b(17:24), 'uint64'), 'ConvertFrom', '.net');
h.center_Hz = typecast(b(25:32), 'double');
h.span_Hz   = typecast(b(33:40), 'double');
end

function [P, t, recNo] = readRecordTraces(files, nPts)
% Each record = 44-byte header + nPts float32 values (dBm) + 20 zero bytes.
%   offset 4   uint32  record number
%   offset 8   uint64  timestamp, .NET ticks
%   offset 44  single  trace values (dBm)
% One .TRC file can hold several records back to back; a session can have
% several dTRC folders. Records are sorted by their number.
recLen = 44 + 4*nPts + 20;
[pk, tk, rk] = deal(cell(numel(files), 1));
for k = 1:numel(files)
    b = readBytes(fullfile(files(k).folder, files(k).name));
    n = floor(numel(b)/recLen);
    if mod(numel(b), recLen) ~= 0
        warning('ff_instrument_record:size', '%s: size %d is not a multiple of %d bytes (%d points).', ...
            files(k).name, numel(b), recLen, nPts);
    end
    B = reshape(b(1:n*recLen), recLen, n);
    rk{k} = double(typecast(reshape(B(5:8, :), [], 1), 'uint32'));
    tk{k} = typecast(reshape(B(9:16, :), [], 1), 'uint64');
    pk{k} = reshape(double(typecast(reshape(B(45:44 + 4*nPts, :), [], 1), 'single')), nPts, n).';
end
recNo = vertcat(zeros(0, 1), rk{:});
t = datetime(vertcat(zeros(0, 1, 'uint64'), tk{:}), 'ConvertFrom', '.net');
P = vertcat(zeros(0, nPts), pk{:});
[recNo, order] = sort(recNo);
t = t(order); P = P(order, :);
end

function b = readBytes(file)
fid = fopen(file, 'r');
if fid < 0, error('Cannot open %s', file); end
b = fread(fid, Inf, '*uint8');
fclose(fid);
end

function T = levelVsTime(f, P, t, rbw)
% Power in the whole span per sweep. With the RBW: bandpower of the PSD
% (point power / noise bandwidth, Gaussian RBW filter: 1.065 x RBW).
% Without it, the power average of the points (dBm per point).
sxx = 10.^(P/10);
if isfinite(rbw) && numel(f) > 1
    lvl = 10*log10(bandpower(sxx.'/(1.065*rbw), f(:), 'psd')).';
    name = 'SpanPower_dBm';
else
    lvl = 10*log10(mean(sxx, 2));
    name = 'MeanLevel_dBm';
end
T = table(t, lvl(:), max(P, [], 2), 'VariableNames', {'Time', name, 'Peak_dBm'});
end

%% ---- small helpers ----------------------------------------------------------

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
    args = [{ip} args]; ip = "";                        % called as ff_instrument_record('Session', x)
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

function s = mergeStruct(s, add)
for n = string(fieldnames(add)).', s.(n) = add.(n); end
end

function v = getfieldOr(S, name, default)
if isfield(S, name), v = S.(name); else, v = default; end
end

function v = defaultTo(v, d)
if isempty(v), v = d; end
end

function v = finiteOr(v, d)
if ~isfinite(v), v = d; end
end

function [scale, unit] = freqUnit(fmax)
if fmax >= 1e9,     scale = 1e9; unit = "GHz";
elseif fmax >= 1e6, scale = 1e6; unit = "MHz";
elseif fmax >= 1e3, scale = 1e3; unit = "kHz";
else,               scale = 1;   unit = "Hz";
end
end

function [scale, unit] = timeUnit(smax)
if smax >= 3*3600,  scale = 3600; unit = "h";
elseif smax >= 180, scale = 60;   unit = "min";
else,               scale = 1;    unit = "s";
end
end

function s = fmtHz(x)
if ~isfinite(x), s = '-'; return, end
[scale, unit] = freqUnit(abs(x));
s = sprintf('%.7g %s', x/scale, unit);
end

function s = fmtDur(sec)
if ~isfinite(sec), s = 'until stopped'; return, end
if sec >= 3600, s = sprintf('%.2f h', sec/3600);
elseif sec >= 120, s = sprintf('%.1f min', sec/60);
else, s = sprintf('%.0f s', sec);
end
end
