function h = fieldfox_studio(varargin)
%FIELDFOX_STUDIO  One window for the whole FieldFox suite: connect, quick actions, launch every tool.
%
%   fieldfox_studio                          open the studio (it does NOT connect by itself)
%   fieldfox_studio('IP', "10.90.36.200")
%   h = fieldfox_studio('Visible', 'off')    build it hidden (scripts, tests): h = handles + functions
%
% Window
%   Connection     IP and port, Connect / Disconnect, the instrument's IDN, its options
%                  and the current mode (Refresh; "SA mode" switches to the spectrum
%                  analyzer). Nothing is sent to the FieldFox before you press Connect.
%   Quick actions  with a connection:
%                    Screenshot        the FieldFox screen, shown here and saved as PNG
%                    Identify signals  ff_signal_id on the trace the FieldFox shows now
%                                      (passive, SA mode): table + labelled plot + recording
%                    Save state        all settings to <State file> on the FieldFox, plus a copy on the PC
%                    Restore state     recall <State file> on the FieldFox (every setting changes)
%                  always: FieldFox Data Link (Keysight's PC program), results folder,
%                  README, copy the last command, clear the log, and the MATLAB apps
%                  that go with the suite (Satellite Link Budget Analyzer, Wireless
%                  Waveform Generator, RF Budget Analyzer, Antenna Designer, ...; only
%                  the apps installed here are listed).
%   Tool cards     one card per tool, grouped by purpose (the card title): Live
%                  monitoring, Inverter EMI, Signal ID, VNA / antenna, Meters,
%                  On-instrument recording, GNSS, Satellite / CubeSat, 5G NTN, Offline
%                  analysis. Each card: what it does, the MATLAB toolboxes it combines,
%                  the FieldFox option it needs, a few key settings, the data source
%                  (FieldFox / demo / file) and Launch. "Command" prints the exact MATLAB
%                  line and copies it: paste it in a script to repeat the run without
%                  the studio. Empty number boxes = keep the FieldFox setting (or the
%                  tool's default). "More" (satellite and 5G cards) takes extra
%                  Name-Value pairs in MATLAB syntax, e.g. 'PolLoss', 3 (see help <tool>).
%   Right side     FieldFox screen, the main result table of the last run, and a log.
%
% Satellite / CubeSat and 5G (Satlab SRS-4 S-band downlink, 2200-2290 MHz, 1 MBd):
%   ff_srs4_analyzer  measures the downlink with the FieldFox (listens on port 2 only)
%   ff_cubesat_link   pass plan + link budget + Doppler (TLE from celestrak.org, 2 h cache)
%   ff_ccsds_linksim  Monte-Carlo BER/PER of the CCSDS frame chain; 'Measured' / 'Pass' modes
%   ff_5g_ntn         5G NR NTN from a LEO pass + coexistence with a measured SRS-4 spectrum
%   Demo chains       the tools connected offline (analyzer -> link sim, pass plan -> link
%                     sim, analyzer -> 5G coexistence)
%   SRS-4 SAFETY: it delivers up to +33 dBm, the FieldFox input is damaged above +27 dBm.
%   Over the air keep >= 1 m between the antennas; conducted use a 30-40 dB pad (first pad
%   >= 10 W) and enter it as Ext gain. Never let the SRS-4 transmit into an open connector.
%   No I/Q on this FieldFox: modulation order and EVM cannot be measured; ppm-level
%   frequency claims need an external 10 MHz GPSDO (FieldFox +-0.7 ppm, SRS-4 spec 0.5 ppm).
%
% Launching a tool
%   - The studio DISCONNECTS first: every tool opens its own connection and only one
%     program should drive the FieldFox at a time. It reconnects afterwards if it was
%     connected, so you can see that the tool put the FieldFox back as it was.
%   - The studio waits until the tool has finished; close the tool's window to stop it
%     (computations such as ff_ccsds_linksim or ff_5g_ntn: Ctrl+C in the Command Window;
%     tools that took control still restore the FieldFox).
%     Questions (inverter OFF/ON, confirmations) come as a dialog or in the MATLAB
%     Command Window: keep an eye on it.
%   - Results go to the base workspace (R_spectrum_pro, T_signal_id, ...) and files to
%     <OutDir>\<tool>\.
%   - RF OUT (red) marks the tools that transmit: ff_vna (port 1, and port 2 for a
%     2-port measurement) and ff_vvm_logger (port 1). Disconnect a LISN, limiter or
%     receiver first. Everything else only listens on port 2 (max +27 dBm, 50 V DC)
%     or uses no instrument at all ("no instrument": the studio stays connected).
%
% Name-Value options (defaults):
%   'IP'        "10.90.36.200"   shown in the IP box
%   'Port'      5025             SCPI socket port
%   'OutDir'    <this folder>\results   tool results, screenshots and state copies
%                                (created only when something is saved)
%   'Visible'   'on'             'off' = hidden window; the tools then run without windows
%   'Prompt'    true             false = no dialogs; the tools get 'Prompt', false
%   'DataLink'  "C:\Program Files (x86)\Keysight\FieldFox Data Link 7.x\FFoxDataLink.exe"
%
% Output h (scripts and tests): Figure, Cards (per tool: Fields, Source, File, Browse,
% Launch, Command, Needs, Panel, Spec), Apps (the MATLAB apps listed) and the functions
% connect(), disconnect(), refresh(), toSA(), setSource(key, id[, files]), command(key),
% launch(key), screenshot(), identify(), saveState(), restoreState(), openDataLink(dry),
% openResults(dry), openReadme(dry), openApp(name, dry), log(), instrument(), close().
% Tool keys: fieldfox_live, ff_spectrum_pro, ff_emi_study, ff_signal_id, ff_vna,
% ff_power_meter, ff_vvm_logger, ff_instrument_record, ff_gnss_margin, ff_srs4_analyzer,
% ff_cubesat_link, ff_ccsds_linksim, ff_5g_ntn, sat_demo, analyze_fieldfox, file_check.
% A launch returns the tool's outputs (R, T, info; the demo chains also C, L, N) and puts
% them in the base workspace as R_<tool>, e.g. R_srs4_analyzer, L_sat_demo.
% Example without the instrument:
%   h = fieldfox_studio('Visible', 'off', 'Prompt', false);
%   h.setSource("ff_vna", "demoF");  out = h.launch("ff_vna");  out.R.metrics
%   h.setSource("ff_srs4_analyzer", "demo");  out = h.launch("ff_srs4_analyzer");  out.R.limits
%   h.close();
%
% Good to know: there is no quasi-peak/CISPR detector (option 361), so EMI results
% are an indicative pre-scan; no I/Q (351), so signal identification is a hypothesis
% from the spectrum shape; below ~30 MHz an antenna near the inverter is in its near
% field; the GPS almanac shipped with MATLAB is from 2020. With nothing connected
% to the ports the tools measure the FieldFox's own noise. The satellite and 5G tools are
% predictions or models (TLE age, assumed antennas, 3GPP study values) except the
% spectrum measurements of ff_srs4_analyzer. See README_FieldFox_Suite.md (section 9).
%
% See also FieldFox, fieldfox_live, ff_spectrum_pro, ff_emi_study, ff_signal_id, ff_vna,
% ff_srs4_analyzer, ff_cubesat_link, ff_ccsds_linksim, ff_5g_ntn.

here = string(fileparts(mfilename('fullpath')));
opt = struct('IP', "10.90.36.200", 'Port', 5025, 'OutDir', "", 'Visible', 'on', 'Prompt', true, ...
    'DataLink', "C:\Program Files (x86)\Keysight\FieldFox Data Link 7.x\FFoxDataLink.exe");
if mod(numel(varargin), 2), error('Options must be Name, Value pairs.'); end
names = fieldnames(opt);
for ia = 1:2:numel(varargin)
    jo = find(strcmpi(names, varargin{ia}), 1);
    if isempty(jo), error('Unknown option "%s".', string(varargin{ia})); end
    opt.(names{jo}) = varargin{ia + 1};
end
if strlength(string(opt.OutDir)) == 0, opt.OutDir = fullfile(here, "results"); end
opt.OutDir = string(opt.OutDir);
opt.Prompt = logical(opt.Prompt);
show = any(strcmpi(string(opt.Visible), ["on" "true" "1"]));
vis = char(matlab.lang.OnOffSwitchState(show));

ff = [];            % FieldFox object, only between Connect and Disconnect
fmt0 = "";          % FORM:DATA at connect (put back after Identify signals)
busy = false;
lastCmd = "";
specs = toolCards();
cards = struct();
apps = appList();
apps = apps(arrayfun(@(a) exist(a.fn) > 0, apps));     %#ok<EXIST> only the apps installed here

%% Window
sz = [1580 920];
ss = get(groot, 'ScreenSize');
if ss(3) > 800, sz = min(sz, ss(3:4) - [40 90]); end       % laptop screens: the tool cards scroll
fig = uifigure('Name', 'FieldFox Studio', 'Position', [30 40 sz], 'Visible', 'off', ...
    'HandleVisibility', 'off', 'CloseRequestFcn', @(~, ~) closeStudio(), 'DeleteFcn', @(~, ~) dropConnection());
main = uigridlayout(fig, [3 2], 'ColumnWidth', {'1x', 440}, 'RowHeight', {30, 214, '1x'}, ...
    'Padding', 8, 'RowSpacing', 6, 'ColumnSpacing', 8);

hdr = uigridlayout(main, [1 2], 'ColumnWidth', {'1x', 'fit'}, 'Padding', 0);
hdr.Layout.Row = 1; hdr.Layout.Column = [1 2];
uilabel(hdr, 'FontSize', 16, 'FontWeight', 'bold', 'Text', ...
    "FieldFox Studio  -  Keysight FieldFox N9914A + MATLAB toolboxes  (inverter EMI, CubeSat SRS-4, 5G NTN)");
statusLbl = uilabel(hdr, 'Text', 'Not connected', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');

top = uigridlayout(main, [1 2], 'ColumnWidth', {'1x', '1x'}, 'Padding', 0, 'ColumnSpacing', 8);
top.Layout.Row = 2; top.Layout.Column = 1;

% Connection
gc = uigridlayout(uipanel(top, 'Title', 'Connection', 'FontWeight', 'bold'), [4 6], ...
    'ColumnWidth', {74, '1x', 32, 62, 86, 86}, 'RowHeight', {24, 22, '1x', 24}, 'Padding', 6, 'RowSpacing', 5);
place(uilabel(gc, 'Text', 'IP address'), 1, 1);
ipEdit = place(uieditfield(gc, 'text', 'Value', string(opt.IP)), 1, 2);
place(uilabel(gc, 'Text', 'Port', 'HorizontalAlignment', 'right'), 1, 3);
portEdit = place(uieditfield(gc, 'numeric', 'Value', opt.Port, 'Limits', [1 65535], ...
    'RoundFractionalValues', 'on', 'ValueDisplayFormat', '%d'), 1, 4);
btnConnect = place(uibutton(gc, 'Text', 'Connect', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~, ~) connect()), 1, 5);
btnDisconnect = place(uibutton(gc, 'Text', 'Disconnect', 'Enable', 'off', 'ButtonPushedFcn', @(~, ~) disconnect()), 1, 6);
place(uilabel(gc, 'Text', 'Instrument'), 2, 1);
idnLbl = place(uilabel(gc, 'Text', '-'), 2, [2 6]);
place(uilabel(gc, 'Text', 'Options', 'VerticalAlignment', 'top'), 3, 1);
optLbl = place(uilabel(gc, 'Text', '-', 'WordWrap', 'on', 'VerticalAlignment', 'top', 'FontSize', 11), 3, [2 6]);
place(uilabel(gc, 'Text', 'Mode'), 4, 1);
modeLbl = place(uilabel(gc, 'Text', '-', 'FontWeight', 'bold'), 4, [2 4]);
btnRefresh = place(uibutton(gc, 'Text', 'Refresh', 'Enable', 'off', 'ButtonPushedFcn', @(~, ~) refresh()), 4, 5);
btnSA = place(uibutton(gc, 'Text', 'SA mode', 'Enable', 'off', 'ButtonPushedFcn', @(~, ~) toSA(), ...
    'Tooltip', 'Switch the FieldFox to the spectrum analyzer (fieldfox_live and Identify signals need it)'), 4, 6);

% Quick actions
gq = uigridlayout(uipanel(top, 'Title', 'Quick actions', 'FontWeight', 'bold'), [5 4], ...
    'RowHeight', {26, 24, 26, 24, '1x'}, 'ColumnWidth', {'1x', '1x', '1x', '1x'}, 'Padding', 6, 'RowSpacing', 6);
btnShot = place(uibutton(gq, 'Text', 'Screenshot', 'Enable', 'off', 'ButtonPushedFcn', @(~, ~) screenshot(), ...
    'Tooltip', 'FieldFox screen shown on the right and saved as PNG in <OutDir>\screenshots'), 1, 1);
btnIdent = place(uibutton(gq, 'Text', 'Identify signals', 'Enable', 'off', 'ButtonPushedFcn', @(~, ~) identify(), ...
    'Tooltip', 'ff_signal_id on the trace the FieldFox shows now (SA mode, passive, ~12 s the first time)'), 1, 2);
btnSave = place(uibutton(gq, 'Text', 'Save state', 'Enable', 'off', 'ButtonPushedFcn', @(~, ~) saveState(), ...
    'Tooltip', 'MMEM:STOR:STAT <State file> on the FieldFox + a copy in <OutDir>\states'), 1, 3);
btnLoad = place(uibutton(gq, 'Text', 'Restore state', 'Enable', 'off', 'ButtonPushedFcn', @(~, ~) restoreState(), ...
    'Tooltip', 'MMEM:LOAD:STAT <State file>: all FieldFox settings change'), 1, 4);
place(uilabel(gq, 'Text', 'State file on FieldFox', 'HorizontalAlignment', 'right'), 2, 1);
stateEdit = place(uieditfield(gq, 'text', 'Value', 'studio_state.sta'), 2, [2 3]);
place(uibutton(gq, 'Text', 'FieldFox Data Link', 'ButtonPushedFcn', @(~, ~) openDataLink(false), ...
    'Tooltip', 'Keysight''s PC program (remote display, files). Do not run a suite tool at the same time.'), 2, 4);
place(uibutton(gq, 'Text', 'Results folder', 'ButtonPushedFcn', @(~, ~) openResults(false)), 3, 1);
place(uibutton(gq, 'Text', 'README', 'ButtonPushedFcn', @(~, ~) openReadme(false)), 3, 2);
place(uibutton(gq, 'Text', 'Copy last command', 'ButtonPushedFcn', @(~, ~) copyCommand(lastCmd)), 3, 3);
place(uibutton(gq, 'Text', 'Clear log', 'ButtonPushedFcn', @(~, ~) set(logArea, 'Value', "")), 3, 4);
place(uilabel(gq, 'Text', 'MATLAB app', 'HorizontalAlignment', 'right'), 4, 1);
if isempty(apps), appItems = "(no app installed)"; appData = ""; else, appItems = [apps.label]; appData = [apps.fn]; end
appDrop = place(uidropdown(gq, 'Items', appItems, 'ItemsData', appData, ...
    'ValueChangedFcn', @(src, ~) set(src, 'Tooltip', appTip(src.Value))), 4, [2 3]);
appDrop.Tooltip = appTip(appDrop.Value);
btnApp = place(uibutton(gq, 'Text', 'Open app', 'ButtonPushedFcn', @(~, ~) openApp(appDrop.Value, false)), 4, 4);
if isempty(apps), appDrop.Enable = 'off'; btnApp.Enable = 'off'; end
place(uilabel(gq, 'FontSize', 11, 'WordWrap', 'on', 'FontColor', [0.7 0 0], 'VerticalAlignment', 'top', 'Text', ...
    "RF OUT only from ff_vna and ff_vvm_logger (port 1; port 2 too for 2-port NA). " + ...
    "All other tools only listen on port 2: max +27 dBm, 50 V DC. SRS-4 (up to +33 dBm): " + ...
    ">= 1 m over the air or a 30-40 dB pad, never direct; never let it transmit into an open connector."), 5, [1 4]);

% Right side: screen, result table, log
right = uigridlayout(main, [3 1], 'RowHeight', {300, '1x', '1x'}, 'Padding', 0, 'RowSpacing', 6);
right.Layout.Row = [2 3]; right.Layout.Column = 2;
screenImg = uiimage(uigridlayout(uipanel(right, 'Title', 'FieldFox screen', 'FontWeight', 'bold'), [1 1], ...
    'Padding', 4), 'ScaleMethod', 'fit', 'Visible', 'off', 'Tooltip', 'Quick actions > Screenshot');
resPanel = uipanel(right, 'Title', 'Results', 'FontWeight', 'bold');
resTbl = uitable(uigridlayout(resPanel, [1 1], 'Padding', 4), 'Data', table(), 'RowName', {});
logArea = uitextarea(uigridlayout(uipanel(right, 'Title', 'Log', 'FontWeight', 'bold'), [1 1], 'Padding', 4), ...
    'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 11, 'Value', "");

% Tool cards
pt = uipanel(main, 'FontWeight', 'bold', 'Title', ...
    'Tools  -  Launch runs the tool (the studio waits), Command prints the MATLAB line');
pt.Layout.Row = 3; pt.Layout.Column = 1;
nRow = ceil(numel(specs)/3);
rowH = zeros(1, nRow);                              % each row as tall as its tallest card
for ic = 1:numel(specs)
    rowH(ceil(ic/3)) = max(rowH(ceil(ic/3)), 241 + 25*fieldRows(specs(ic).fields));
end
gt = uigridlayout(pt, [nRow 3], 'RowHeight', num2cell(rowH), 'ColumnWidth', {'1x', '1x', '1x'}, ...
    'Scrollable', 'on', 'Padding', 6, 'RowSpacing', 8, 'ColumnSpacing', 8);
for i = 1:numel(specs)
    cards.(specs(i).key) = makeCard(gt, specs(i));
    sourceChanged(specs(i).key);
end
updateNeeds();

h = struct('Figure', fig, 'Cards', cards, 'Apps', apps, 'connect', @connect, 'disconnect', @disconnect, ...
    'refresh', @refresh, 'toSA', @toSA, 'setSource', @setSource, 'command', @commandFor, 'launch', @launch, ...
    'screenshot', @screenshot, 'identify', @identify, 'saveState', @saveState, ...
    'restoreState', @restoreState, 'openDataLink', @openDataLink, 'openResults', @openResults, ...
    'openReadme', @openReadme, 'openApp', @openApp, 'log', @() string(logArea.Value), ...
    'instrument', @instrumentNow, 'close', @closeStudio);
logMsg("FieldFox Studio ready. Not connected: press Connect to talk to the FieldFox at " + opt.IP + ".");
fig.Visible = vis;
if nargout == 0, clear h, end

%% ===== nested functions (share the variables above) ========================

    function c = makeCard(parent, s)
        p = uipanel(parent, 'Title', s.group, 'FontWeight', 'bold');
        [nfr, fr, fc] = fieldRows(s.fields);
        g = uigridlayout(p, [8 + nfr 4], 'ColumnWidth', {'fit', '1x', 'fit', '1x'}, ...
            'RowHeight', [{18, 44, 44, 16}, repmat({22}, 1, nfr), {'1x', 22, 22, 26}], ...
            'Padding', [6 4 6 6], 'RowSpacing', 3, 'ColumnSpacing', 4);
        place(uilabel(g, 'Text', s.title + "  (" + s.func + ")", 'FontWeight', 'bold'), 1, [1 3]);
        if any([s.sources.rf])
            place(uilabel(g, 'Text', 'RF OUT', 'FontWeight', 'bold', 'FontColor', 'w', ...
                'BackgroundColor', [0.8 0.1 0.1], 'HorizontalAlignment', 'center'), 1, 4);
        elseif any([s.sources.live])
            place(uilabel(g, 'Text', 'listens only', 'FontColor', [0 0.5 0], 'HorizontalAlignment', 'center'), 1, 4);
        else
            place(uilabel(g, 'Text', 'no instrument', 'FontColor', [0.4 0.4 0.4], 'HorizontalAlignment', 'center'), 1, 4);
        end
        place(uilabel(g, 'Text', s.desc, 'WordWrap', 'on', 'VerticalAlignment', 'top', 'FontSize', 11), 2, [1 4]);
        place(uilabel(g, 'Text', "Toolboxes: " + s.tbx, 'WordWrap', 'on', 'VerticalAlignment', 'top', ...
            'FontSize', 10, 'FontColor', [0.25 0.3 0.55]), 3, [1 4]);
        needs = place(uilabel(g, 'Text', '', 'FontSize', 10), 4, [1 4]);
        F = struct();
        for k = 1:numel(s.fields)
            d = s.fields(k);
            place(uilabel(g, 'Text', d.label, 'HorizontalAlignment', 'right', 'Tooltip', d.tip), 4 + fr(k), fc(k));
            span = fc(k) + 1;
            if d.type == "code", span = [2 4]; end      % "More": the whole row
            F.(d.name) = place(makeField(g, d), 4 + fr(k), span);
        end
        r = 6 + nfr;
        place(uilabel(g, 'Text', 'Data source', 'HorizontalAlignment', 'right'), r, 1);
        sd = place(uidropdown(g, 'Items', [s.sources.label], 'ItemsData', [s.sources.id], ...
            'ValueChangedFcn', @(~, ~) sourceChanged(s.key)), r, [2 4]);
        fe = place(uieditfield(g, 'text', 'Placeholder', 'file or folder (Browse...)'), r + 1, [1 3]);
        fb = place(uibutton(g, 'Text', 'Browse...', 'ButtonPushedFcn', @(~, ~) browse(s.key)), r + 1, 4);
        lb = place(uibutton(g, 'Text', 'Launch', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~, ~) launch(s.key)), r + 2, [1 2]);
        cb = place(uibutton(g, 'Text', 'Command', 'ButtonPushedFcn', @(~, ~) showCommand(s.key), ...
            'Tooltip', 'Print the MATLAB line for these settings and copy it'), r + 2, [3 4]);
        c = struct('Panel', p, 'Fields', F, 'Source', sd, 'File', fe, 'Browse', fb, 'Launch', lb, ...
            'Command', cb, 'Needs', needs, 'Spec', s);
    end

    function sourceChanged(key)
        c = cards.(key);
        si = c.Spec.sources([c.Spec.sources.id] == string(c.Source.Value));
        useFile = strlength(si.file) > 0;
        c.File.Enable = useFile; c.Browse.Enable = useFile;
        c.File.Placeholder = filePlaceholder(si.file);
        for d = c.Spec.fields
            c.Fields.(d.name).Enable = isempty(d.srcs) || any(d.srcs == si.id);
        end
    end

    function setSource(key, id, file)
        % Scripts: choose a card's data source (and file), as the dropdown would
        cards.(key).Source.Value = string(id);
        sourceChanged(key);
        if nargin > 2, cards.(key).File.Value = strjoin(string(file), "; "); end
    end

    function browse(key)
        c = cards.(key);
        si = c.Spec.sources([c.Spec.sources.id] == string(c.Source.Value));
        start = char(opt.OutDir);
        if ~isfolder(start), start = char(here); end
        if si.file == "dir"
            p = uigetdir(start, 'Select the folder');
            if isequal(p, 0), return, end
            P = string(p);
        else
            [flt, multi] = fileFilter(si.file);
            [n, d] = uigetfile(flt, char("Select: " + si.label), start, 'MultiSelect', multi);
            if isequal(n, 0), return, end
            P = fullfile(string(d), string(n));
        end
        c.File.Value = strjoin(P, "; ");
        try figure(fig); catch, end
    end

    %% ---- connection ------------------------------------------------------
    function ok = connect()
        ok = false;
        if busy, return, end
        if ~isempty(ff), disconnect(); end
        ip = string(strtrim(ipEdit.Value)); port = portEdit.Value;
        logMsg(sprintf('Connecting to %s:%d ...', ip, port));
        statusLbl.Text = 'Connecting ...'; drawnow;
        try
            ff = FieldFox(ip, port);
            fmt0 = ff.query("FORM:DATA?");
        catch err
            ff = [];
            msg = regexprep(string(err.message), '<[^>]+>', '');    % no HTML links in the log
            logMsg("Could not connect: " + msg);
            statusLbl.Text = 'Not connected';
            alert("Could not connect to " + ip + ": " + msg, 'Connection');
            return
        end
        ok = true;
        idnLbl.Text = ff.IDN;
        optLbl.Text = optionText(ff.Options);
        logMsg("Connected: " + ff.IDN);
        setConnected(true);
        refresh();
    end

    function disconnect()
        if isempty(ff), return, end
        try delete(ff); catch, end      % the studio never takes control: nothing to restore
        ff = [];
        setConnected(false);
        logMsg('Disconnected from the FieldFox.');
    end

    function refresh()
        if ~needConnection(), return, end
        try
            m = ff.mode();
            modeLbl.Text = m + "  (" + modeName(m) + ")";
        catch err
            logMsg("Lost the connection: " + err.message);
            disconnect();
        end
    end

    function toSA()
        if ~needConnection(), return, end
        try
            ff.setMode("SA");
            logMsg('Switched the FieldFox to SA mode.');
        catch err
            logMsg("Could not switch to SA: " + err.message);
        end
        refresh();
    end

    function setConnected(tf)
        set([btnDisconnect btnRefresh btnSA btnShot btnIdent btnSave btnLoad], 'Enable', tf);
        if tf
            statusLbl.Text = "Connected to " + ff.IP;
            statusLbl.FontColor = [0 0.5 0];
        else
            statusLbl.Text = 'Not connected'; statusLbl.FontColor = [0 0 0];
            idnLbl.Text = '-'; optLbl.Text = '-'; modeLbl.Text = '-';
        end
        updateNeeds();
    end

    function updateNeeds()
        for s = specs
            L = cards.(s.key).Needs;
            if isempty(s.needs)
                L.Text = "No FieldFox needed" + s.note; L.FontColor = [0.4 0.4 0.4];
            elseif isempty(ff)
                L.Text = "Needs FieldFox option " + strjoin(s.needs, " + "); L.FontColor = [0.3 0.3 0.3];
            else
                miss = s.needs(~ismember(s.needs, ff.Options));
                if isempty(miss)
                    L.Text = "Option " + strjoin(s.needs, " + ") + ": installed"; L.FontColor = [0 0.5 0];
                else
                    L.Text = "Option " + strjoin(miss, " + ") + " NOT installed on this FieldFox"; L.FontColor = [0.8 0 0];
                end
            end
        end
    end

    function x = instrumentNow()
        x = ff;         % the studio's FieldFox object ([] = not connected)
    end

    function ok = needConnection()
        ok = ~isempty(ff);
        if ~ok, logMsg('Press Connect first.'); end
    end

    %% ---- quick actions ---------------------------------------------------
    function file = screenshot()
        file = "";
        if ~needConnection() || busy, return, end
        f = fullfile(opt.OutDir, "screenshots", "fieldfox_screen_" + stamp() + ".png");
        ensureDir(fileparts(f));
        try
            ff.screenshot(f);
            screenImg.ImageSource = char(f); screenImg.Visible = 'on';
            file = f;
            logMsg("Screenshot saved: " + f);
        catch err
            logMsg("Screenshot failed: " + err.message);
        end
    end

    function T = identify()
        T = table();
        if ~needConnection() || busy, return, end
        d = fullfile(opt.OutDir, "ff_signal_id"); ensureDir(d);
        b = fullfile(d, "ffsid_" + stamp());
        logMsg(">> [T, info] = ff_signal_id(ff, 'Plot', true, 'SaveTo', """ + b + ".mat"")   % ff = the studio's connection");
        setBusy(true, "Identify signals");
        done = onCleanup(@() setBusy(false, ""));
        try
            [T, info] = ff_signal_id(ff, 'Plot', true, 'Visible', vis, 'Prompt', false, 'SaveTo', b + ".mat");
        catch err
            logMsg("Identify signals: " + err.message);
            putFormatBack();
            return
        end
        putFormatBack();
        publish("ff_signal_id", struct('T', T, 'info', info), b + ".png");
        if height(T) == 0, logMsg('No signals above the noise floor.'); end
    end

    function putFormatBack()
        % ff_signal_id leaves FORM:DATA on REAL,32 when it borrows our connection
        try
            if ff.query("FORM:DATA?") ~= fmt0
                ff.write("FORM:DATA " + fmt0);
                ff.syncFormat();
            end
        catch err
            logMsg("Could not put FORM:DATA back to " + fmt0 + ": " + err.message);
        end
    end

    function ok = saveState()
        ok = false;
        if ~needConnection() || busy, return, end
        name = stateName();
        if name == "", return, end
        try
            ff.saveState(name);
            bytes = ff.getFile(name);
            local = fullfile(opt.OutDir, "states", name);
            ensureDir(fileparts(local));
            fid = fopen(local, 'w');
            try fwrite(fid, bytes, 'uint8'); catch e, if fid >= 0, fclose(fid); end, rethrow(e); end
            fclose(fid);
            logMsg(sprintf('State saved on the FieldFox as %s (copy on the PC: %s, %d bytes).', name, local, numel(bytes)));
            ok = true;
        catch err
            logMsg("Save state failed: " + err.message);
        end
        reportErrors();
    end

    function ok = restoreState()
        ok = false;
        if ~needConnection() || busy, return, end
        name = stateName();
        if name == "", return, end
        if show && opt.Prompt
            q = "Recall " + name + " on the FieldFox? Every setting changes. If the state was saved in " + ...
                "NA, CAT or VVM mode the FieldFox transmits from port 1 afterwards.";
            if ~strcmp(uiconfirm(fig, q, 'Restore state', 'Options', {'Recall', 'Cancel'}, ...
                    'DefaultOption', 2, 'CancelOption', 2, 'Icon', 'warning'), 'Recall')
                return
            end
        end
        try
            ff.loadState(name);
            logMsg("State " + name + " recalled on the FieldFox.");
            ok = true;
        catch err
            logMsg("Restore state failed: " + err.message);
        end
        reportErrors();
        refresh();
    end

    function name = stateName()
        name = strtrim(string(stateEdit.Value));
        if ~endsWith(lower(name), ".sta"), name = name + ".sta"; end
        if isempty(regexp(name, '^[\w\-]+\.sta$', 'once', 'ignorecase'))
            logMsg('State file: use letters, digits, _ and - only (e.g. studio_state.sta).');
            name = "";
        end
    end

    function reportErrors()
        try
            e = ff.errors();
            if ~isempty(e), logMsg("FieldFox reported: " + strjoin(e, "; ")); end
        catch
        end
    end

    function [exe, ok] = openDataLink(dry)
        exe = string(opt.DataLink);
        if ~isfile(exe)
            d = dir("C:\Program Files (x86)\Keysight\FieldFox Data Link*\FFoxDataLink.exe");
            if ~isempty(d), exe = string(fullfile(d(1).folder, d(1).name)); end
        end
        ok = isfile(exe);
        if ~ok
            logMsg("FieldFox Data Link not found at " + exe + ". Install it (keysight.com) or pass 'DataLink', <path to FFoxDataLink.exe>.");
        elseif ~dry
            logMsg("Starting FieldFox Data Link: " + exe + " (it talks to the FieldFox too: do not run a suite tool at the same time).");
            winopen(exe);
        end
    end

    function d = openResults(dry)
        d = opt.OutDir;
        if dry, return, end
        ensureDir(d);
        winopen(d);
    end

    function [f, ok] = openReadme(dry)
        f = fullfile(here, "README_FieldFox_Suite.md");
        ok = isfile(f);
        if ~ok, logMsg("README not found: " + f); elseif ~dry, edit(f); end
    end

    function [fn, ok] = openApp(name, dry)
        % One of the listed MATLAB apps (dry = only check that it is there)
        if nargin < 1 || strlength(string(name)) == 0, name = appDrop.Value; end
        if nargin < 2, dry = false; end
        fn = string(name);
        a = apps([apps.fn] == fn);
        ok = isscalar(a);
        if ~ok, logMsg("MATLAB app not available here: " + fn); return, end
        if dry, return, end
        logMsg(">> " + fn + "   % " + a.label + ": " + a.use);
        try
            feval(fn);
        catch err
            ok = false;
            logMsg("Could not open " + fn + ": " + err.message);
        end
    end

    function t = appTip(fn)
        a = apps([apps.fn] == string(fn));
        t = "MATLAB apps that go with the suite (only the installed ones are listed)";
        if isscalar(a), t = a.label + " (" + a.fn + ", " + a.tbx + "): " + a.use; end
    end

    %% ---- tools -------------------------------------------------------------
    function [cmd, F] = commandFor(key, forRun)
        % forRun: the line the studio runs (ff_vna then asks through its dialogs);
        % otherwise the line to paste into a script
        if nargin < 2, forRun = false; end
        c = cards.(key); s = c.Spec;
        si = s.sources([s.sources.id] == string(c.Source.Value));
        F = strip(split(string(c.File.Value), ";")).';
        F = F(strlength(F) > 0);
        if strlength(si.file) > 0
            if isempty(F) && si.file ~= "imat"      % "imat": the file is optional
                error('studio:file', '%s: choose a file or folder first (Browse...).', si.label);
            end
            bad = F(~isfile(F) & ~isfolder(F));
            if ~isempty(bad), error('studio:file', 'Not found: %s', strjoin(bad, ", ")); end
        else
            F = strings(1, 0);                      % ignore a stale path from another source
        end
        x = struct('ip', string(strtrim(ipEdit.Value)), 'port', portEdit.Value, 'vis', vis, 'show', show, ...
            'prompt', opt.Prompt, 'ask', show && opt.Prompt && forRun, 'out', opt.OutDir, 'stamp', stamp(), 'here', here);
        if si.live && strlength(x.ip) == 0, error('studio:ip', 'Enter the FieldFox IP address.'); end
        cmd = s.build(fieldValues(c), si.id, F, x);
    end

    function showCommand(key)
        try
            cmd = commandFor(key);
        catch err
            logMsg(err.message); return
        end
        logMsg(cmd);
        lastCmd = cmd;
        copyCommand(cmd);
    end

    function out = launch(key)
        out = struct();
        if busy, logMsg('A tool is running: wait until it has finished.'); return, end
        c = cards.(key); s = c.Spec;
        si = s.sources([s.sources.id] == string(c.Source.Value));
        try
            [cmd, F] = commandFor(key, true);
            shown = commandFor(key);
        catch err
            logMsg(err.message); alert(err.message, s.title); return
        end
        wasConnected = ~isempty(ff);
        if si.live && wasConnected
            logMsg("Disconnecting the studio: " + s.func + " opens its own connection.");
            disconnect();
        end
        if si.rf
            logMsg("*** RF OUTPUT *** " + s.func + " transmits from the FieldFox test port(s). " + ...
                "Disconnect anything that must not be driven (LISN, limiter, receiver).");
        end
        if si.live && strlength(s.warn) > 0, logMsg(s.warn); end
        dOut = fullfile(opt.OutDir, s.key);
        before = strings(0, 1);
        if any(s.key == ["analyze_fieldfox" "file_check"])
            dOut = "";
        else
            ensureDir(dOut);
            before = string({dir(dOut).name});
        end
        lastCmd = shown;
        logMsg(">> " + shown);
        if cmd ~= shown, logMsg("   (the studio shows ff_vna's instructions as dialogs)"); end
        setBusy(true, s.func);
        idle = onCleanup(@() unlock());            % also after Ctrl+C
        t0 = tic;
        try
            if s.key == "file_check"
                out.T = fileCheck(F(1));
            else
                out = runCommand(cmd, @askDialog);
            end
        catch err
            setBusy(false, "");
            fprintf(2, '%s\n', getReport(err, 'extended', 'hyperlinks', 'off'));    % where it failed
            msg = strtrim(regexprep(string(err.message), '(See )?<a [^>]*>[^<]*</a>( for [^.\n]*\.)?|<[^>]+>', ''));
            logMsg("ERROR in " + s.func + ": " + msg);
            alert(msg, s.func + " stopped");
            out.error = err;
            if si.live && wasConnected, connect(); end
            return
        end
        setBusy(false, "");
        logMsg(sprintf('%s finished (%.0f s).', s.func, toc(t0)));
        png = "";
        if dOut ~= "", png = fullfile(dOut, "ffsid_" + stamp() + ".png"); end
        publish(key, out, png);
        listed = false;
        for n = intersect(["C" "R" "L" "N"], string(fieldnames(out)).', 'stable')
            o = out.(n);
            if ~(isstruct(o) && isscalar(o) && isfield(o, 'files') && ~isempty(o.files)), continue, end
            fl = regexprep(string(o.files(:)).', '^.*[\\/]', '');
            inDir = "";
            if isfield(o, 'outDir'), inDir = " in " + string(o.outDir); end
            etc = ""; if numel(fl) > 8, etc = ", ..."; end
            logMsg(sprintf('Files of %s (%d)%s: %s%s', n, numel(fl), inDir, strjoin(fl(1:min(8, end)), ", "), etc));
            listed = true;
        end
        if ~listed && dOut ~= ""
            new = setdiff(string({dir(dOut).name}), [before "." ".."]);
            if ~isempty(new), logMsg("New in " + dOut + ": " + strjoin(new(1:min(10, end)), ", ")); end
        end
        if s.key == "analyze_fieldfox" && si.id == "run"
            logMsg("analyze_fieldfox writes its figures and CSV files to " + fullfile(here, "results"));
        end
        if si.live && wasConnected, connect(); end
    end

    function publish(key, out, png)
        % Outputs to the base workspace, main table to the Results panel
        suffix = erase(key, "ff_");
        done = strings(0, 1);
        for n = string(fieldnames(out)).'
            if n == "error", continue, end
            assignin('base', n + "_" + suffix, out.(n));
            done(end + 1) = n + "_" + suffix; %#ok<AGROW>
        end
        if ~isempty(done), logMsg("In the base workspace: " + strjoin(done, ", ")); end
        [T, from] = mainTable(out, cards.(key).Spec.table);
        if strlength(from) > 0
            try
                resTbl.Data = T;
            catch
                resTbl.Data = varfun(@string, T);
            end
            name = regexprep(from, '^(\w+)', "$1_" + suffix);      % R.metrics -> R_vna.metrics
            resPanel.Title = "Results: " + cards.(key).Spec.func + "  (" + name + " in the workspace)";
        end
        if isfield(out, 'info') && isstruct(out.info) && isfield(out.info, 'Figure') ...
                && isgraphics(out.info.Figure) && strlength(png) > 0
            exportgraphics(out.info.Figure, png, 'Resolution', 120);
            logMsg("Labelled spectrum saved: " + png);
            if ~show, delete(out.info.Figure); end      % hidden studio: nobody could see or close it
        end
    end

    function r = askDialog(msg)
        % ff_vna's instructions (QuickCal steps, RF warning) as a dialog
        msg = string(msg);
        if contains(msg, "type n")                  % the RF-output confirmation
            msg = "*** RF OUTPUT *** ff_vna is about to transmit from port 1 (a 2-port " + ...
                "measurement: port 2 too). Disconnect LISNs, limiters and receivers that must " + ...
                "not be driven. Port damage level +27 dBm, 50 V DC.";
            o = {'Continue (RF on)', 'Cancel'}; reply = ["" "n"];
        elseif contains(msg, "type s")
            o = {'Done', 'Skip'}; reply = ["" "s"];
        else
            o = {'Done'}; reply = "";
        end
        msg = replace(msg, ["press Enter, or type s + Enter to skip" "press Enter"], ["press Done (or Skip)" "press Done"]);
        sel = uiconfirm(fig, msg, 'ff_vna', 'Options', o, 'DefaultOption', 1, ...
            'CancelOption', numel(o), 'Icon', 'warning');
        r = reply(strcmp(o, sel));
    end

    function setBusy(tf, name)
        busy = tf;
        for s = specs, cards.(s.key).Launch.Enable = ~tf; end
        set([btnConnect btnShot btnIdent btnSave btnLoad], 'Enable', ~tf);
        if ~tf, setConnected(~isempty(ff)); return, end
        set([btnDisconnect btnRefresh btnSA], 'Enable', false);
        statusLbl.Text = "Running " + name + " ...  stop: close its window or Ctrl+C in the Command Window; answer questions there or in the dialog";
        statusLbl.FontColor = [0.75 0.35 0];
        drawnow;
    end

    %% ---- log and closing -----------------------------------------------------
    function logMsg(msg)
        msg = string(msg);
        fprintf('%s\n', msg);
        if ~isvalid(logArea), return, end
        v = string(logArea.Value);
        v = v(strlength(v) > 0 | (1:numel(v)).' > 1);
        lines = splitlines(msg);
        lines(1) = string(datetime('now', 'Format', 'HH:mm:ss')) + "  " + lines(1);
        v = [v(:); lines(:)];
        logArea.Value = v(max(1, end - 499):end);
        try scroll(logArea, 'bottom'); catch, end
    end

    function alert(msg, ttl)
        if show && opt.Prompt && isvalid(fig), uialert(fig, msg, ttl); end
    end

    function closeStudio()
        if busy
            logMsg('A tool is still running: stop it first (close its window), then close the studio.');
            return
        end
        disconnect();
        delete(fig);
    end

    function dropConnection()
        % The window is gone, also without the close box (delete(h.Figure), close all force):
        % never keep the socket open while a script still holds h
        if ~isempty(ff), try delete(ff); catch, end, ff = []; end
    end

    function unlock()
        % Ctrl+C in a tool skips launch's catch: the studio must not stay busy (and unclosable)
        if busy && isvalid(fig), setBusy(false, ""); end
    end
end

%% ===== local functions =======================================================

function out = runCommand(cmd, ask) %#ok<INUSD> 'ask' is the dialog handle used in ff_vna commands
% Runs the generated MATLAB line in its own workspace and returns R, T, info
% (and C, L, N of the demo chains).
eval(cmd);
out = struct();
for n = intersect(["R" "T" "info" "C" "L" "N"], string(who()).', 'stable')
    out.(n) = eval(n);
end
end

function u = place(u, row, col)
u.Layout.Row = row; u.Layout.Column = col;
end

function u = makeField(g, d)
switch d.type
    case "bool"
        u = uicheckbox(g, 'Text', '', 'Value', d.default, 'Tooltip', d.tip);
    case "list"
        if isempty(d.data)
            u = uidropdown(g, 'Items', d.items, 'Value', d.default, 'Tooltip', d.tip);
        else
            u = uidropdown(g, 'Items', d.items, 'ItemsData', d.data, 'Value', d.default, 'Tooltip', d.tip);
        end
    case "code"
        u = uieditfield(g, 'text', 'Value', d.default, 'Tooltip', d.tip, 'FontName', 'Consolas', ...
            'Placeholder', 'e.g. ''PolLoss'', 3   (Name-Value pairs, MATLAB syntax)');
    otherwise                                       % num, vec, text
        u = uieditfield(g, 'text', 'Value', d.default, 'Tooltip', d.tip);
end
end

function [n, rows, cols] = fieldRows(fields)
% Card layout: two fields per row, a "code" field (More) on a row of its own at the end
normal = string({fields.type}) ~= "code";
k = cumsum(normal);
rows = ceil(k/2); cols = 1 + 2*mod(k - 1, 2);
nn = ceil(nnz(normal)/2);
rows(~normal) = nn + (1:nnz(~normal)); cols(~normal) = 1;
n = nn + nnz(~normal);
if isempty(fields), n = 0; end
end

function v = fieldValues(c)
% Current card settings: numbers scaled to SI units, [] = empty box or "keep"
v = struct();
for d = c.Spec.fields
    u = c.Fields.(d.name);
    switch d.type
        case "num"
            t = strtrim(string(u.Value));
            val = [];
            if t ~= ""
                val = str2double(t);
                if isnan(val), error('studio:value', '%s: "%s" is not a number.', d.label, t); end
                val = val * d.scale;
                if d.scale > 1 && isfinite(val), val = round(val, 3); end
            end
        case "vec"
            val = sscanf(char(replace(string(u.Value), ",", " ")), '%f').';
        case "bool"
            val = logical(u.Value);
        case "list"
            val = string(u.Value);
            if val == "keep", val = []; end
        case "code"
            val = codeText(u.Value, d.label);
        otherwise
            val = string(strtrim(u.Value));
    end
    v.(d.name) = val;
end
end

function t = codeText(t, label)
% "More": extra Name-Value pairs typed as MATLAB text, checked for syntax (not run here)
t = strtrim(string(t));
t = regexprep(t, '^,\s*|,\s*$', '');                % tolerate a leading/trailing comma
if t == "", return, end
ok = ~isempty(regexp(t, '^[''"]', 'once'));         % starts with a name in quotes
if ok
    try
        str2func(char("@() f(" + t + ")"));         % parses the text, runs nothing
    catch
        ok = false;
    end
end
if ~ok
    error('studio:value', '%s: "%s" is not a list of Name-Value pairs in MATLAB syntax, e.g. ''PolLoss'', 3', label, t);
end
end

function [T, from] = mainTable(out, paths)
% First non-empty table among e.g. "R.metrics", "T", "L.measured.table"
T = table(); from = "";
for p = paths
    x = out; ok = true;
    for q = split(p, ".").'
        if ~(isstruct(x) && isscalar(x) && isfield(x, q)), ok = false; break, end
        x = x.(q);
    end
    if ok && istable(x) && (height(x) > 0 || p == paths(end))
        T = x; from = p; return
    end
end
end

%% ---- tool cards ----------------------------------------------------------------

function S = toolCards()
S = [
    card("fieldfox_live", "Live monitoring", "Live monitor", ...
        "Watches the SA trace on the FieldFox: live/average/max hold, waterfall, peaks with band names, channel power, inverter OFF baseline. Read-only: changes no setting.", ...
        "base MATLAB (tcpclient, no VISA)", "233", ...
        src("live", "FieldFox (live, read-only; SA mode)", "", true, false), ...
        [fld("Record", "Record", "bool", false, 'Tip', 'Save every sweep (fieldfox_rec_<time>.mat)'), ...
        fld("AutoBaseline", "Auto baseline", "bool", false, 'Tip', 'The first 20 sweeps (inverter OFF) become the baseline'), ...
        fld("Threshold", "Threshold dB", "num", "6", 'Tip', 'dB above the baseline to call an emission new'), ...
        fld("MaxSweeps", "Max sweeps", "num", "Inf", 'Tip', 'Stop by itself (Inf = until you close the window)')], ...
        @bLive, strings(1, 0))
    card("ff_spectrum_pro", "Live monitoring", "Spectrum pro", ...
        "Controlled sweeps in MATLAB's spectrumAnalyzer: inverter-OFF baseline + mask, channel power / OBW / ACPR, cross-checked with the FieldFox's own channel power and frequency counter.", ...
        "DSP System (spectrumAnalyzer), Signal Processing (bandpower, obw) + ff_signal_id labels", "233", ...
        [src("live", "FieldFox (takes control, restores after)", "", true, false), ...
        src("demo", "Demo: synthetic inverter switching on", "", false, false), ...
        src("file", "Replay a recording (.mat)", "mat", false, false)], ...
        [fld("Center", "Center MHz", "num", "", 'Scale', 1e6, 'Sources', "live", 'Tip', 'Empty = keep the FieldFox setting'), ...
        fld("Span", "Span MHz", "num", "", 'Scale', 1e6, 'Sources', "live", 'Tip', 'Empty = keep the FieldFox setting'), ...
        fld("RBW", "RBW kHz", "num", "", 'Scale', 1e3, 'Sources', "live", 'Tip', 'Empty = keep the FieldFox setting'), ...
        fld("Detector", "Detector", "list", "AVER", 'Items', ["keep", "AVER (RMS)", "POS (peak)", "SAMP", "NEG"], ...
        'Data', ["keep" "AVER" "POS" "SAMP" "NEG"], 'Sources', "live", 'Tip', 'AVER (RMS) gives correct channel power'), ...
        fld("AutoBaseline", "Baseline sweeps", "num", "20", 'Sources', ["live" "file"], ...
        'Tip', 'Inverter OFF: the first N sweeps become the baseline (0 = none). With Prompt you are asked to switch it.'), ...
        fld("MaxSweeps", "Max sweeps", "num", "Inf", 'Tip', 'Stop by itself (Inf = until you close a window)')], ...
        @bSpectrum, ["R.Signals" "R.Log"])
    card("ff_emi_study", "Inverter EMI", "Inverter EMI study", ...
        "OFF-ON-OFF scans 150 kHz-3 GHz + GNSS L1, statistics per frequency, switching-harmonic comb, zero-span bursts, indicative CISPR limits, interference range, PDF report. Pre-scan only (no QP).", ...
        "Statistics & ML (ranksum, bootci, fitlm), Bioinformatics (mafdr), Signal Processing, Communications (propagationModel), Report Generator + ff_signal_id", "233", ...
        [src("live", "FieldFox: OFF-ON-OFF study (dialogs guide you)", "", true, false), ...
        src("demo", "Demo study (synthetic, about 2 min)", "", false, false), ...
        src("file", "Saved study / 1-3 recordings / study folder", "mat3", false, false)], ...
        [fld("Bands", "Bands (1-5)", "vec", "1 2 3 4 5", 'Sources', "live", 'Tip', ...
        '1 Conducted 150k-30M, 2 VHF 30-300M, 3 UHF 300M-1G, 4 L/S 1-3G, 5 GNSS L1'), ...
        fld("Sweeps", "Sweeps", "num", "20", 'Sources', "live", 'Tip', 'Per band, detector and inverter state (at least 3; 20 recommended)'), ...
        fld("Preamp", "Preamp", "bool", false, 'Sources', "live", 'Tip', 'Option 235: lower noise floor, overloads sooner'), ...
        fld("Protocol", "Protocol", "list", "ABA", 'Items', ["ABA", "AB"], 'Sources', "live", 'Tip', 'ABA = OFF, ON, OFF (shows ambient drift)'), ...
        fld("Report", "Report", "list", "pdf", 'Items', ["pdf", "docx", "none"]), ...
        fld("Notes", "Notes", "text", "", 'Tip', 'Site, inverter model, antenna, distance ... (goes into the report)')], ...
        @bEmi, ["R.summary" "R.emissions"])
    card("ff_signal_id", "Signal ID", "Signal identification", ...
        "Finds the signals in a spectrum and guesses the technology (LTE, NR, Wi-Fi, FM, GSM, BLE, DVB-T2, CW ...) with band and channel. No I/Q: a hypothesis from the spectrum shape.", ...
        "5G, LTE, WLAN, Bluetooth, Communications, Satellite Comms (templates), DSP System, Signal Processing, Statistics & ML, Deep Learning (1-D CNN)", "233", ...
        [src("live", "FieldFox screen now (passive, SA mode)", "", true, false), ...
        src("demo", "Demo: synthetic signal", "", false, false), ...
        src("file", "Recording (.mat)", "mat", false, false)], ...
        [fld("Method", "Method", "list", "ensemble", 'Items', ["ensemble", "template", "ml", "cnn", "none"], ...
        'Tip', 'ensemble = templates + bagged trees + CNN; none = detection only'), ...
        fld("Threshold", "Threshold dB", "num", "6", 'Tip', 'dB above the noise floor'), ...
        fld("Detector", "Detector", "list", "auto", 'Items', ["auto", "peak", "average", "sample"], ...
        'Sources', ["live" "file"], 'Tip', 'auto = read from the FieldFox or the recording'), ...
        fld("DemoSignal", "Demo signal", "list", "LTE 10 MHz", 'Items', ["LTE 10 MHz", "NR 20 MHz", ...
        "Wi-Fi 20 MHz", "FM broadcast", "GSM", "BLE 1M", "DVB-T2 8 MHz", "CW"], 'Sources', "demo")], ...
        @bSignal, "T")
    card("ff_vna", "VNA / antenna", "S-parameters / antenna", ...
        "Cables, filters, antennas: return loss, VSWR, Smith chart, TDR, distance to fault, matching network, QuickCal, Touchstone files, comparison with an Antenna Toolbox model.", ...
        "RF (sparameters, smithplot, rational, matchingnetwork), Signal Integrity (tdr), Signal Processing, Antenna", "210", ...
        [src("na", "FieldFox NA mode (RF OUT)", "", true, true), ...
        src("track", "FieldFox SA + tracking generator (RF OUT)", "", true, true), ...
        src("demoA", "Demo antenna (synthetic)", "", false, false), ...
        src("demoF", "Demo filter (synthetic)", "", false, false), ...
        src("file", "Touchstone file (.s1p / .s2p)", "snp", false, false)], ...
        [fld("Ports", "Ports", "list", "1", 'Items', ["1 (S11)", "2 (full 2-port)"], 'Data', ["1" "2"], ...
        'Sources', "na", 'Tip', '2 ports: RF out of port 1 AND port 2 (option 211)'), ...
        fld("Cal", "QuickCal", "list", "none", 'Items', ["none", "quick1", "quick2"], 'Sources', "na", ...
        'Tip', 'Option 112, no cal kit needed; you are told what to connect'), ...
        fld("Start", "Start MHz", "num", "1", 'Scale', 1e6, 'Sources', ["na" "track"]), ...
        fld("Stop", "Stop MHz", "num", "1000", 'Scale', 1e6, 'Sources', ["na" "track"]), ...
        fld("Points", "Points", "num", "401", 'Sources', ["na" "track"]), ...
        fld("Antenna", "Model", "list", "none", 'Items', ["none", "dipole", "monopole", "patchMicrostrip"], ...
        'Tip', 'Antenna Toolbox model compared with the measured S11 (dipole ~20 s, patch ~90 s)')], ...
        @bVna, "R.metrics")
    card("ff_power_meter", "Meters", "Channel power logger", ...
        "Logs the FieldFox's built-in channel power meter over time: strip chart, histogram with fitted distributions, statistics and an SA cross-check.", ...
        "Statistics & ML (fitdist, paramci), Signal Processing (bandpower)", "310", ...
        [src("live", "FieldFox CPM mode (restores after)", "", true, false), ...
        src("demo", "Demo (synthetic readings)", "", false, false), ...
        src("file", "Saved log (cpm_*.mat)", "mat", false, false)], ...
        [fld("Center", "Center MHz", "num", "", 'Scale', 1e6, 'Sources', "live", 'Tip', 'Empty = keep the CPM setting'), ...
        fld("ChannelBW", "Channel MHz", "num", "", 'Scale', 1e6, 'Sources', "live", 'Tip', '10 Hz - 100 MHz; empty = keep'), ...
        fld("Atten", "Atten dB", "list", "keep", 'Items', ["keep", "10", "30"], 'Sources', "live", 'Tip', 'CPM accepts only 10 or 30 dB'), ...
        fld("Duration", "Duration s", "num", "60", 'Sources', "live"), ...
        fld("Interval", "Interval s", "num", "1", 'Sources', "live")], ...
        @bPower, "R.stats")
    card("ff_vvm_logger", "Meters", "Magnitude/phase logger", ...
        "Vector voltmeter at one frequency: magnitude and phase vs time for cable or connector stability, drift fit and Allan deviation.", ...
        "Statistics & ML (fitlm), Navigation (allanvar)", "308", ...
        [src("live", "FieldFox VVM mode (RF OUT, port 1)", "", true, true), ...
        src("demo", "Demo (synthetic readings)", "", false, false), ...
        src("file", "Saved log (vvm_*.mat)", "mat", false, false)], ...
        [fld("Parameter", "Parameter", "list", "S11", 'Items', ["S11", "S21"], 'Sources', "live", 'Tip', 'S21 = port 1 -> device -> port 2'), ...
        fld("Frequency", "Freq MHz", "num", "", 'Scale', 1e6, 'Sources', "live", 'Tip', 'Empty = keep the VVM setting'), ...
        fld("Power", "Power dBm", "num", "", 'Sources', "live", 'Tip', '-45 to +3 dBm; empty = keep'), ...
        fld("Zero", "Zero at start", "bool", false, 'Sources', "live", 'Tip', 'Relative readings: the first is 0 dB / 0 deg'), ...
        fld("Duration", "Duration s", "num", "60", 'Sources', "live"), ...
        fld("Interval", "Interval s", "num", "1", 'Sources', "live")], ...
        @bVvm, "R.adev")
    card("ff_instrument_record", "On-instrument recording", "Recorder", ...
        "Long unattended recording saved on the FieldFox itself (it carries on if the PC sleeps), then download, parse and waterfall. Not gap-free (no real-time option).", ...
        "Signal Processing (bandpower) + ff_signal_id; FTP fallback for the download", "236", ...
        [src("rec", "Record on the FieldFox, then download", "", true, false), ...
        src("fetch", "Download the newest session", "", true, false), ...
        src("parse", "Parse a copied RecordPlayback folder", "dir", false, false)], ...
        [fld("Start", "Start MHz", "num", "", 'Scale', 1e6, 'Sources', "rec", 'Tip', 'Empty = keep the FieldFox setting'), ...
        fld("Stop", "Stop MHz", "num", "", 'Scale', 1e6, 'Sources', "rec", 'Tip', 'Empty = keep the FieldFox setting'), ...
        fld("RBW", "RBW kHz", "num", "", 'Scale', 1e3, 'Sources', ["rec" "parse"], ...
        'Tip', 'Parse: give the RBW used, otherwise ff_signal_id cannot label the result'), ...
        fld("Duration", "Duration min", "num", "10", 'Scale', 60, 'Sources', "rec"), ...
        fld("Interval", "Interval s", "num", "0", 'Sources', "rec", 'Tip', '0 = every sweep')], ...
        @bRecord, ["R.signals" "R.level"])
    card("ff_gnss_margin", "GNSS", "GNSS L1 margin", ...
        "How much the inverter's noise at GPS L1 lowers the C/N0 of every visible satellite: link budget + receiver chain + FieldFox L1 noise, OFF vs ON. The shipped almanac is from 2020.", ...
        "Satellite Comms (satelliteScenario, link), Antenna (patch design), Phased Array, RF (rfbudget), Navigation (semread, lookangles), Signal Processing", ["233" "235"], ...
        [src("live", "FieldFox L1 noise + satellites", "", true, false), ...
        src("sat", "Satellites only (no instrument)", "", false, false), ...
        src("file", "L1 recording(s): OFF [; ON]", "mat2", false, false)], ...
        [fld("Lat", "Latitude", "num", "13.7563", 'Tip', 'deg; default central Bangkok: SET YOUR SITE'), ...
        fld("Lon", "Longitude", "num", "100.5018", 'Tip', 'deg'), ...
        fld("Antenna", "Antenna", "list", "patch", 'Items', ["patch", "isotropic"], 'Tip', 'patch = Antenna Toolbox design at L1 (~20 s)'), ...
        fld("Cable", "Cable", "text", "", 'Tip', '.s2p file (e.g. from ff_vna) or a loss in dB; empty = no cable'), ...
        fld("Mode", "Mode", "list", "offon", 'Items', ["offon", "single"], 'Sources', "live", 'Tip', 'offon = inverter OFF, then ON'), ...
        fld("Sweeps", "Sweeps", "num", "50", 'Sources', "live", 'Tip', 'Averaged per state')], ...
        @bGnss, "R.sats")
    card("ff_srs4_analyzer", "Satellite / CubeSat", "SRS-4 analyzer", ...
        "Finds the SRS-4 downlink (2200-2290 MHz): fc and ppm error, Rs, roll-off, OBW, C/N0, ECSS mask, harmonics, envelope, bursts, EIRP; PDF report. No I/Q: no EVM or constellation.", ...
        "Satellite Comms (ccsdsTMWaveformGenerator), Optimization (lsqcurvefit), Statistics & ML (nlparci, bootstrap), Signal Processing (obw), Phased Array (fspl), Antenna, Report Generator", "233", ...
        [src("live", "FieldFox (takes control, port 2, restores after)", "", true, false), ...
        src("demo", "Demo: synthetic SRS-4 + FieldFox model (~1 min)", "", false, false), ...
        src("file", "Re-analyse srs4_run.mat or a recording (.mat)", "mat", false, false)], ...
        [fld("ExpectedFreq", "Freq MHz", "num", "", 'Scale', 1e6, 'Tip', ...
        'SRS-4 tx.freq if you know it (then the ppm error is judged); empty = unknown. Demo: the synthetic carrier (default 2245 MHz)'), ...
        fld("ExpectedRate", "Rate kBd", "num", "1000", 'Scale', 1e3, 'Tip', 'sl-radio --tx-rate (1000 = 1 MBd)'), ...
        fld("Distance", "Distance m", "num", "", 'Tip', 'Over the air: TX antenna to FieldFox antenna (keep >= 1 m); empty = conducted'), ...
        fld("TxGain", "TX ant. dBi", "text", "", 'Tip', ...
        'SRS-4 antenna gain in dBi (e.g. 6.1 for a Satlab patch), or patch / dipole (Antenna Toolbox design); empty = 0 dBi'), ...
        fld("ExtGain", "Ext gain dB", "num", "", 'Sources', "live", 'Tip', ...
        'Pad + cable between the SRS-4 and port 2, negative = loss (e.g. -40.6 for a 40 dB pad); empty = 0 (over the air)'), ...
        fld("DemoMod", "Demo mod.", "list", "QPSK", 'Items', ["QPSK", "BPSK", "8PSK"], 'Sources', "demo", ...
        'Tip', 'Modulation of the synthetic SRS-4 (the envelope statistics hint at BPSK vs QPSK/8PSK)'), ...
        fld("Report", "Report", "list", "pdf", 'Items', ["pdf", "docx", "none"]), ...
        fld("More", "More", "code", "", 'Tip', ...
        'Extra Name-Value pairs, e.g. ''RxAntennaGain'', 2, ''PolLoss'', 3, ''ExpectedRolloff'', 0.2, ''LineScan'', true (help ff_srs4_analyzer)')], ...
        @bSrs4, "R.limits", 'Warn', "SRS-4 SAFETY: it delivers up to +33 dBm, FieldFox port 2 is damaged above +27 dBm. " + ...
        "Over the air keep >= 1 m between the antennas; conducted only through a 30-40 dB pad (first pad >= 10 W) " + ...
        "entered as Ext gain. Never let the SRS-4 transmit into an open connector.")
    card("ff_cubesat_link", "Satellite / CubeSat", "CubeSat pass planner", ...
        "Passes of a CubeSat (KNACKSAT-2 by default) over your site: time-varying SRS-4 link budget and margin per configuration, data volume, ITU-R P.618 losses, Doppler pre-compensation table; PDF.", ...
        "Satellite Comms (satelliteScenario, access, link, dopplershift, p618, satelliteCNR), Aerospace (tleread), Antenna, Phased Array, Communications (bercoding), Navigation, Mapping, Report Gen.", ...
        strings(1, 0), ...
        [src("online", "CelesTrak TLE (internet, cached 2 h)", "", false, false), ...
        src("offline", "Offline: cached TLE or the sample file", "", false, false), ...
        src("tle", "TLE file (.tle / .txt)", "tle", false, false), ...
        src("orbit", "Circular sun-synchronous orbit (Altitude)", "", false, false)], ...
        [fld("NORAD", "NORAD", "num", "67683", 'Sources', ["online" "offline" "tle"], 'Tip', '67683 KNACKSAT-2, 46320 NAPA-1 (Thai CubeSats)'), ...
        fld("Site", "Site", "vec", "13.7563 100.5018 10", 'Tip', 'lat lon alt_m of YOUR ground station (default: central Bangkok)'), ...
        fld("Duration", "Hours", "num", "24", 'Tip', 'Time window searched for passes'), ...
        fld("MinElevation", "Min elev.", "num", "10", 'Tip', 'deg, elevation mask'), ...
        fld("TxPower", "TX dBm", "num", "33", 'Tip', ...
        'SRS-4 output, 20-33 dBm. --tx-gain=50 is not dBm: measure it (ff_srs4_analyzer, conducted through a pad)'), ...
        fld("Altitude", "Altitude km", "num", "500", 'Sources', "orbit", 'Tip', 'Circular sun-synchronous orbit (LTAN 10:30)'), ...
        fld("Track", "FieldFox track", "bool", false, 'Tip', 'Also write fieldfox_track.csv: centre frequency vs time that follows the Doppler'), ...
        fld("Report", "Report", "list", "pdf", 'Items', ["pdf", "docx", "none"]), ...
        fld("More", "More", "code", "", 'Tip', ...
        'Extra Name-Value pairs, e.g. ''Start'', "2026-10-11 06:00", ''SymbolRate'', 2e6, ''Modulation'', "BPSK", ''DishDiameter'', 3.7 (help ff_cubesat_link)')], ...
        @bCubesat, "R.passes", 'Note', "; TLE from celestrak.org only (2 h cache)")
    card("ff_ccsds_linksim", "Satellite / CubeSat", "CCSDS link simulator", ...
        "Monte-Carlo BER/PER of the SRS-4 CCSDS chain (RS, randomizer, CC, RRC-PSK) with an ideal or realistic receiver, checked against datasheet Fig. 13 / Table 7; PER over passes or at a measured C/N0.", ...
        "Satellite Comms (ccsdsTMWaveformGenerator, ccsdsRSDecode), Communications (vitdec, SymbolSynchronizer, CarrierSynchronizer, berconfint), DSP System, Signal Processing, Statistics", ...
        strings(1, 0), ...
        [src("sim", "BER/PER curves only (~1-2 min)", "", false, false), ...
        src("measured", "+ 'Measured': ff_srs4_analyzer srs4_run.mat", "mat", false, false), ...
        src("pass", "+ 'Pass': ff_cubesat_link cubesat_link.mat", "mat", false, false), ...
        src("both", "+ 'Pass' and 'Measured' (both files)", "pm", false, false), ...
        src("demo", "+ 'Pass' and 'Measured' with demo inputs", "", false, false)], ...
        [fld("Modulation", "Modulation", "list", "keep", 'Items', ["auto", "QPSK", "BPSK", "8PSK"], ...
        'Data', ["keep" "QPSK" "BPSK" "8PSK"], 'Tip', 'auto = from the ''Pass'' input, otherwise QPSK'), ...
        fld("Coding", "Coding", "list", "keep", 'Items', ["auto", "CC + RS", "CC", "RS", "none"], ...
        'Data', ["keep" "concatenated" "CC" "RS" "none"], 'Tip', 'auto = from the ''Pass'' input, otherwise CC + RS (SRS-4 default)'), ...
        fld("Receiver", "Receiver", "list", "ideal", 'Items', ["ideal", "realistic"], ...
        'Tip', 'realistic = AFC + Gardner timing + carrier PLL + ASM frame sync (slower)'), ...
        fld("Impairments", "Impairments", "list", "none", 'Items', ["none", "srs4", "leo"], ...
        'Tip', 'srs4 = 0.5 ppm, assumed phase noise, IQ imbalance; leo = srs4 + 45 kHz Doppler and -800 Hz/s'), ...
        fld("Preset", "Preset", "list", "quick", 'Items', ["quick", "full"], 'Tip', 'full = exam-quality curves (10-30 min)'), ...
        fld("Validate", "Validate", "bool", true, 'Tip', 'Datasheet Figure 13 sweep (realistic receiver) + Table 7 comparison'), ...
        fld("Cache", "Cache", "list", "on", 'Items', ["on", "off", "refresh"], ...
        'Tip', 'support\linksim\linksim_cache.mat keeps simulated points: a second run is fast; off = write nothing'), ...
        fld("More", "More", "code", "", 'Tip', ...
        'Extra Name-Value pairs, e.g. ''SymbolRate'', 2e6, ''Rolloff'', 0.35, ''AFCRange'', 100e3, ''Table7Sweep'', true, ''ImplLoss'', 1.5 (help ff_ccsds_linksim)')], ...
        @bLinksim, ["R.measured.table" "R.pass.table" "R.curves"], 'Note', "; reads analyzer / pass results")
    card("ff_5g_ntn", "5G NTN", "5G NR from a LEO satellite", ...
        "5G NR from a LEO pass: Doppler/delay vs numerology, PDSCH BLER in NTN-TDL channels with/without Doppler pre-compensation, throughput, coexistence with an SRS-4.", ...
        "5G (nrDLSCH, nrPDSCH, nrOFDMModulate, nrTDLChannel NTN, nrChannelEstimate, nrEqualizeMMSE), Satellite Comms (satelliteScenario, dopplershift, p618), Aerospace, Mapping, Report Generator", ...
        strings(1, 0), ...
        [src("sim", "Simulate: circular LEO pass over the site (~2-3 min)", "imat", false, false), ...
        src("tle", "Simulate: real orbit from a TLE (CelesTrak, 2 h cache)", "imat", false, false), ...
        src("reuse", "Reuse the BLER curves of ntn5g_run.mat (~40 s what-if)", "rmat", false, false)], ...
        [fld("MaxElevation", "Max elev.", "num", "90", 'Sources', "sim", 'Tip', 'deg, peak elevation of the constructed pass'), ...
        fld("NORAD", "NORAD", "num", "67683", 'Sources', "tle", 'Tip', ...
        '67683 = KNACKSAT-2 (its orbit is used for geometry only); a TLE file instead: More = ''TLE'', "my.tle"'), ...
        fld("Band", "Band", "list", "n256", 'Items', ["n256", "n255"], 'Tip', 'n256 = S-band DL 2170-2200 MHz (next to the SRS-4 band), n255 = L-band'), ...
        fld("ChannelBW", "Channel", "list", "5", 'Items', ["5 MHz", "10 MHz", "15 MHz", "20 MHz"], 'Data', ["5" "10" "15" "20"]), ...
        fld("EIRPDensity", "EIRP dBW/MHz", "num", "", 'Tip', 'Satellite EIRP density; empty = TR 38.821 Set-1 S-band (34 for LEO-600)'), ...
        fld("Preset", "Preset", "list", "quick", 'Items', ["quick", "full"], 'Sources', ["sim" "tle"], ...
        'Tip', 'full = smooth BLER curves (10-30 min)'), ...
        fld("Interferer", "Interferer", "list", "none", 'Items', ["none", "demo SRS-4", "file (box below)"], ...
        'Data', ["none" "demo" "file"], 'Tip', ...
        'Coexistence: how much an SRS-4 spectrum leaks into the 5G channel. file = an ff_srs4_analyzer srs4_run.mat or a recording'), ...
        fld("Report", "Report", "list", "pdf", 'Items', ["pdf", "docx", "none"]), ...
        fld("More", "More", "code", "", 'Tip', ...
        'Extra Name-Value pairs, e.g. ''Site'', [13.73 100.78 5], ''Channel'', "NTN-TDL-A", ''InterfVictimDistance'', 100 (help ff_5g_ntn)')], ...
        @b5g, ["R.interference.summary" "R.throughput"], 'Note', "; can read an ff_srs4_analyzer run")
    card("sat_demo", "Offline analysis", "Satellite demo chains", ...
        "The satellite tools connected, without instrument or internet: analyzer demo -> link sim 'Measured'; offline pass plan -> link sim 'Pass'; analyzer demo at 2201 MHz -> 5G NTN coexistence.", ...
        "What the four satellite / 5G tools combine: Satellite Comms, 5G, Communications, Signal Processing, Optimization, Statistics & ML, Aerospace, Antenna, Phased Array, Report Generator", ...
        strings(1, 0), ...
        [src("meas", "SRS-4 demo -> link sim 'Measured' (~2-3 min)", "", false, false), ...
        src("pass", "Offline pass plan -> link sim 'Pass' (~3 min)", "", false, false), ...
        src("coex", "SRS-4 demo at 2201 MHz -> 5G NTN coexistence (~4 min)", "", false, false)], ...
        [fld("Report", "Reports", "list", "none", 'Items', ["none", "pdf"], 'Tip', 'PDF reports of the analyzer / pass plan / 5G run (slower)'), ...
        fld("Cache", "Link-sim cache", "list", "on", 'Items', ["on", "off", "refresh"], 'Sources', ["meas" "pass"], ...
        'Tip', 'support\linksim\linksim_cache.mat (off = write nothing)')], ...
        @bSatDemo, ["L.measured.table" "L.pass.table" "N.interference.summary" "R.limits" "C.passes"])
    card("analyze_fieldfox", "Offline analysis", "Copied files", ...
        "Script: OFF vs ON sweeps from files copied off the FieldFox (.sta, Record/Playback). Edit the folder at its top first; it writes to results\ next to it and runs clear, clc, close all.", ...
        "base MATLAB", strings(1, 0), ...
        [src("run", "Run the script", "", false, false), ...
        src("edit", "Open it in the Editor (settings at the top)", "", false, false)], ...
        repmat(fld("x", "x", "num", ""), 1, 0), @bAnalyze, strings(1, 0))
    card("file_check", "Offline analysis", "Which tool opens it?", ...
        "Pick a recording, Touchstone or TLE file, a run file (srs4_run, cubesat_link, ntn5g_run) or a folder: shows what it holds and the commands of the tools that can open it.", ...
        "base MATLAB", strings(1, 0), ...
        src("check", "File or folder to check", "any", false, false), ...
        repmat(fld("x", "x", "num", ""), 1, 0), @(~, ~, F, ~) "% studio file check: " + F(1), "T")
    ];
S = S.';                                            % a row, so that "for s = S" visits each card
end

function c = card(key, group, title, desc, tbx, needs, sources, fields, build, tbl, varargin)
% varargin: 'Note' (added to "No FieldFox needed"), 'Warn' (logged before a live launch)
func = string(key);
if key == "file_check", func = "file check"; end
if key == "sat_demo", func = "demo chains"; end
x = struct('Note', "", 'Warn', "");
for k = 1:2:numel(varargin), x.(varargin{k}) = string(varargin{k + 1}); end
c = struct('key', string(key), 'group', string(group), 'title', string(title), 'func', func, ...
    'desc', string(desc), 'tbx', string(tbx), 'needs', string(needs), 'sources', sources, ...
    'fields', fields, 'build', build, 'table', string(tbl), 'note', x.Note, 'warn', x.Warn);
end

function s = src(id, label, file, live, rf)
s = struct('id', string(id), 'label', string(label), 'file', string(file), 'live', live, 'rf', rf);
end

function f = fld(name, label, type, default, varargin)
f = struct('name', string(name), 'label', string(label), 'type', string(type), 'default', default, ...
    'scale', 1, 'items', strings(1, 0), 'data', strings(1, 0), 'srcs', strings(1, 0), 'tip', "");
key = struct('Scale', 'scale', 'Items', 'items', 'Data', 'data', 'Sources', 'srcs', 'Tip', 'tip');
for k = 1:2:numel(varargin), f.(key.(varargin{k})) = varargin{k + 1}; end
f.srcs = string(f.srcs); f.tip = string(f.tip);
end

%% ---- command builders: (values, source id, files, context) -> MATLAB line -----

function c = bLive(v, ~, ~, x)
mx = v.MaxSweeps;
if ~x.show && isinf(mx), mx = 100; end              % a hidden window cannot be closed by hand
c = call("", "fieldfox_live", lit(x.ip), [{'Record', v.Record, 'AutoBaseline', v.AutoBaseline, ...
    'Threshold', v.Threshold, 'MaxSweeps', mx, 'OutDir', toolDir(x, "fieldfox_live"), 'Visible', x.vis}, ...
    portArg(x, 'Port')]);
end

function c = bSpectrum(v, src, F, x)
nv = {'MaxSweeps', finiteOr(v.MaxSweeps)};
switch src
    case "live"
        first = lit(x.ip);
        nv = [{'Center', v.Center, 'Span', v.Span, 'RBW', v.RBW, 'Detector', v.Detector, ...
            'AutoBaseline', positiveOr(v.AutoBaseline)}, nv, portArg(x, 'Port')];
    case "demo"
        first = '""'; nv = [{'Replay', "demo"}, nv];
    otherwise
        first = '""'; nv = [{'Replay', F(1), 'AutoBaseline', positiveOr(v.AutoBaseline)}, nv];
end
c = call("R", "ff_spectrum_pro", first, [nv, {'OutDir', toolDir(x, "ff_spectrum_pro"), ...
    'Visible', x.vis, 'Prompt', x.prompt}]);
end

function c = bEmi(v, src, F, x)
pre = "";
nv = {'Report', v.Report, 'Notes', replace(v.Notes, newline, " ")};
switch src
    case "live"
        first = lit(x.ip);
        b = unique(round(v.Bands));
        if isempty(b) || any(b < 1 | b > 5), error('studio:value', 'Bands: numbers 1 to 5, e.g. "1 2 5".'); end
        if ~isequal(b, 1:5)
            pre = "P = ff_emi_study([], 'DryRun', true); ";
            nv = [{'Bands', code("P.bands(" + mat2str(b) + ", :)")}, nv];
        end
        nv = [nv, {'Sweeps', v.Sweeps, 'Preamp', v.Preamp, 'Protocol', v.Protocol}, portArg(x, 'Port')];
    case "demo"
        first = "[]"; nv = [{'Data', "demo"}, nv];
    otherwise
        first = "[]"; nv = [{'Data', F}, nv];
end
c = pre + call("R", "ff_emi_study", first, [nv, {'OutDir', fullfile(toolDir(x, "ff_emi_study"), ...
    "emi_study_" + x.stamp), 'Visible', x.vis, 'Prompt', x.prompt}]);
end

function c = bSignal(v, src, F, x)
det = v.Detector;
if det == "auto", det = []; end
nv = {'Method', v.Method, 'Threshold', v.Threshold, 'Detector', det, 'Plot', true, ...
    'Visible', x.vis, 'Prompt', x.prompt};
switch src
    case "live"
        rec = fullfile(toolDir(x, "ff_signal_id"), "ffsid_" + x.stamp + ".mat");
        c = call("[T, info]", "ff_signal_id", lit(x.ip), [nv, {'SaveTo', rec}, portArg(x, 'Port')]);
    case "demo"
        fc = dictionary(["LTE 10 MHz" "NR 20 MHz" "Wi-Fi 20 MHz" "FM broadcast" "GSM" "BLE 1M" "DVB-T2 8 MHz" "CW"], ...
            [1842.5e6 3500e6 2437e6 98.5e6 947.4e6 2426e6 538e6 433.92e6]);
        nv{6} = "peak";
        c = "M = ff_signal_classifier('Train', false); [f, P, rbw] = M.synth(" + lit(v.DemoSignal) + ...
            ", 'Fc', " + numLit(fc(v.DemoSignal)) + ", 'SNR', 30, 'Points', 40, 'Detector', ""peak""); " + ...
            call("[T, info]", "ff_signal_id", "f, P, rbw", nv);
    otherwise
        c = call("[T, info]", "ff_signal_id", lit(F(1)), nv);
end
end

function c = bVna(v, src, F, x)
ant = v.Antenna; if ant == "none", ant = []; end
cal = v.Cal; if cal == "none", cal = []; end
nv = {'Antenna', ant};
switch src
    case "na"
        first = lit(x.ip);
        nv = [{'Ports', str2double(v.Ports), 'Cal', cal, 'Start', v.Start, 'Stop', v.Stop, 'Points', v.Points}, nv, portArg(x, 'TcpPort')];
    case "track"
        first = lit(x.ip);
        nv = [{'Method', "tracking", 'Start', v.Start, 'Stop', v.Stop, 'Points', v.Points}, nv, portArg(x, 'TcpPort')];
    case "demoA"
        first = "[]"; nv = [{'Demo', "antenna"}, nv];
    case "demoF"
        first = "[]"; nv = [{'Demo', "filter"}, nv];
    otherwise
        first = "[]"; nv = [{'File', F(1)}, nv];
end
pr = x.prompt;
if x.ask && any(src == ["na" "track"]), pr = code("ask"); end   % the studio's dialog for the steps
c = call("R", "ff_vna", first, [nv, {'SaveFigures', true, 'OutDir', toolDir(x, "ff_vna"), ...
    'Visible', x.vis, 'Prompt', pr}]);
end

function c = bGnss(v, src, F, x)
cable = []; t = strtrim(v.Cable);
if t ~= "", cable = str2double(t); if isnan(cable), cable = t; end, end
nv = {'Lat', v.Lat, 'Lon', v.Lon, 'Antenna', v.Antenna, 'Cable', cable};
switch src
    case "live"
        first = lit(x.ip); nv = [nv, {'Mode', v.Mode, 'Sweeps', v.Sweeps}, portArg(x, 'Port')];
    case "sat"
        first = '""';
    otherwise
        first = '""'; nv = [nv, {'Spectrum', F(1)}];
        if numel(F) > 1, nv = [nv, {'SpectrumOn', F(2)}]; end
end
c = call("R", "ff_gnss_margin", first, [nv, {'OutDir', toolDir(x, "ff_gnss_margin"), ...
    'Visible', x.vis, 'Prompt', x.prompt}]);
end

function c = bPower(v, src, F, x)
switch src
    case "live"
        first = lit(x.ip);
        if ~isempty(v.Atten), v.Atten = str2double(v.Atten); end    % "keep" -> []
        nv = [{'Center', v.Center, 'ChannelBW', v.ChannelBW, 'Atten', v.Atten, ...
            'Duration', v.Duration, 'Interval', v.Interval}, portArg(x, 'Port')];
    case "demo"
        first = ""; nv = {'Demo', true};
    otherwise
        first = ""; nv = {'File', F(1)};
end
c = call("R", "ff_power_meter", first, [nv, {'OutDir', toolDir(x, "ff_power_meter"), ...
    'Visible', x.vis, 'Prompt', x.prompt}]);
end

function c = bVvm(v, src, F, x)
switch src
    case "live"
        first = lit(x.ip);
        nv = [{'Parameter', v.Parameter, 'Frequency', v.Frequency, 'Power', v.Power, 'Zero', v.Zero, ...
            'Duration', v.Duration, 'Interval', v.Interval}, portArg(x, 'Port')];
    case "demo"
        first = ""; nv = {'Demo', true};
    otherwise
        first = ""; nv = {'File', F(1)};
end
c = call("R", "ff_vvm_logger", first, [nv, {'OutDir', toolDir(x, "ff_vvm_logger"), ...
    'Visible', x.vis, 'Prompt', x.prompt}]);
end

function c = bRecord(v, src, F, x)
switch src
    case "rec"
        first = lit(x.ip);
        nv = [{'Start', v.Start, 'Stop', v.Stop, 'RBW', v.RBW, 'Duration', v.Duration, ...
            'Interval', v.Interval}, portArg(x, 'Port')];
    case "fetch"
        first = lit(x.ip); nv = [{'Session', "latest"}, portArg(x, 'Port')];
    otherwise
        first = "'parse', " + lit(F(1)); nv = {'RBW', v.RBW};
end
c = call("R", "ff_instrument_record", first, [nv, {'OutDir', toolDir(x, "ff_instrument_record"), ...
    'Visible', x.vis, 'Prompt', x.prompt}]);
end

function c = bAnalyze(~, src, ~, x)
f = lit(fullfile(x.here, "analyze_fieldfox.m"));
if src == "edit", c = "edit(" + f + ");"; else, c = "run(" + f + ");"; end
end

function c = bSrs4(v, src, F, x)
nv = {'ExpectedFreq', v.ExpectedFreq, 'ExpectedRate', v.ExpectedRate, 'Distance', v.Distance, ...
    'TxAntennaGain', gainArg(v.TxGain, "TX ant. dBi")};
switch src
    case "live"
        if ~isempty(v.Distance) && v.Distance < 1
            error('studio:value', 'Distance %.2f m: keep the SRS-4 antenna at least 1 m from the FieldFox antenna.', v.Distance);
        end
        first = lit(x.ip); nv = [nv, {'ExtGain', v.ExtGain}, portArg(x, 'Port')];
    case "demo"
        first = "[]"; nv = [{'Demo', true, 'DemoModulation', v.DemoMod}, nv];
    otherwise
        first = "[]"; nv = [{'File', F(1)}, nv];
end
c = withMore(call("R", "ff_srs4_analyzer", first, [nv, {'Report', v.Report, ...
    'OutDir', fullfile(toolDir(x, "ff_srs4_analyzer"), "srs4_" + x.stamp), 'Visible', x.vis, 'Prompt', x.prompt}]), v.More);
end

function c = bCubesat(v, src, F, x)
if numel(v.Site) ~= 3 || abs(v.Site(1)) > 90 || abs(v.Site(2)) > 180
    error('studio:value', 'Site: three numbers, latitude (deg) longitude (deg) altitude (m), e.g. "13.7563 100.5018 10".');
end
switch src
    case "online",  nv = {'NORAD', v.NORAD};
    case "offline", nv = {'NORAD', v.NORAD, 'Offline', true};
    case "tle",     nv = {'TLE', F(1), 'NORAD', v.NORAD};
    otherwise,      nv = {'Altitude', v.Altitude};
end
nv = [nv, {'Site', v.Site, 'Duration', v.Duration, 'MinElevation', v.MinElevation, 'TxPower_dBm', v.TxPower}, ...
    notDefault('FieldFoxTrack', v.Track, false), {'Report', v.Report, ...
    'OutDir', fullfile(toolDir(x, "ff_cubesat_link"), "cubesat_link_" + x.stamp), 'Visible', x.vis, 'Prompt', x.prompt}];
c = withMore(call("R", "ff_cubesat_link", "", nv), v.More);
end

function c = bLinksim(v, src, F, x)
switch src
    case "measured", nv = {'Measured', F(1)};
    case "pass",     nv = {'Pass', F(1)};
    case "both"
        [p, m] = passAndMeasured(F);
        nv = {'Pass', p, 'Measured', m};
    case "demo",     nv = {'Pass', "demo", 'Measured', "demo"};
    otherwise,       nv = {};
end
nv = [nv, {'Modulation', v.Modulation, 'Coding', v.Coding, 'Receiver', v.Receiver, 'Impairments', v.Impairments, ...
    'Preset', v.Preset}, notDefault('Validate', v.Validate, true), cacheArg(v.Cache), ...
    {'OutDir', fullfile(toolDir(x, "ff_ccsds_linksim"), "linksim_" + x.stamp), 'Visible', x.vis, 'Prompt', x.prompt}];
c = withMore(call("R", "ff_ccsds_linksim", "", nv), v.More);
end

function c = b5g(v, src, F, x)
switch src
    case "sim",  nv = {'MaxElevation', v.MaxElevation, 'Preset', v.Preset};
    case "tle",  nv = {'TLE', v.NORAD, 'Preset', v.Preset};
    otherwise,   nv = {'File', F(1)};
end
k = 1 + (src == "reuse");                           % reuse: ntn5g_run.mat; interferer
if v.Interferer ~= "file" && numel(F) >= k
    error('studio:file', 'A file for the interferer is in the box: set Interferer to "file (box below)", or remove it.');
end
switch v.Interferer
    case "demo", nv = [nv, {'Interference', "demo"}];
    case "file"
        if numel(F) < k
            error('studio:file', ['Interferer = file: put an ff_srs4_analyzer srs4_run.mat or a recording ' ...
                'in the file box (after the ntn5g_run.mat, separated by ";", when reusing curves).']);
        end
        nv = [nv, {'Interference', F(k)}];
end
nv = [nv, {'Band', v.Band, 'ChannelBW', str2double(v.ChannelBW)*1e6, 'EIRPDensity', v.EIRPDensity, 'Report', v.Report, ...
    'OutDir', fullfile(toolDir(x, "ff_5g_ntn"), "ntn5g_" + x.stamp), 'Visible', x.vis, 'Prompt', x.prompt}];
c = withMore(call("R", "ff_5g_ntn", "", nv), v.More);
end

function c = bSatDemo(v, src, ~, x)
% Two tools connected through their returned results (R / C -> L, R -> N)
o = fullfile(toolDir(x, "sat_demo"), "chain_" + x.stamp);
tail = {'Visible', x.vis, 'Prompt', x.prompt};
cache = cacheArg(v.Cache);
switch src
    case "meas"
        c = call("R", "ff_srs4_analyzer", "[]", [{'Demo', true, 'ExpectedFreq', 2245e6, 'Report', v.Report, ...
            'OutDir', fullfile(o, "srs4")}, tail]) + " " + ...
            call("L", "ff_ccsds_linksim", "", [{'Measured', code("R"), 'Validate', false}, cache, ...
            {'OutDir', fullfile(o, "linksim")}, tail]);
    case "pass"
        c = call("C", "ff_cubesat_link", "", [{'Offline', true, 'Report', v.Report, 'OutDir', fullfile(o, "cubesat")}, tail]) + ...
            " " + call("L", "ff_ccsds_linksim", "", [{'Pass', code("C"), 'Validate', false}, cache, ...
            {'OutDir', fullfile(o, "linksim")}, tail]);
    otherwise
        c = call("R", "ff_srs4_analyzer", "[]", [{'Demo', true, 'ExpectedFreq', 2201e6, 'Distance', 3, ...
            'TxAntennaGain', 6.1, 'Report', "none", 'OutDir', fullfile(o, "srs4_2201")}, tail]) + " " + ...
            call("N", "ff_5g_ntn", "", [{'Interference', code("R"), 'Report', v.Report, 'OutDir', fullfile(o, "ntn5g")}, tail]);
end
end

function g = gainArg(t, label)
% antenna gain box: a number in dBi, or "patch" / "dipole" (Antenna Toolbox design)
t = strtrim(string(t)); g = [];
if t == "", return, end
g = str2double(t);
if isnan(g)
    g = lower(t);
    if ~any(g == ["patch" "dipole"]), error('studio:value', '%s: a number in dBi, or patch / dipole.', label); end
end
end

function a = cacheArg(s)
% link-sim cache list: on = the tool's default (left out), off = false, refresh
switch string(s)
    case "off",     a = {'Cache', false};
    case "refresh", a = {'Cache', "refresh"};
    otherwise,      a = {};
end
end

function a = notDefault(name, val, default)
a = {};
if ~isequal(val, default), a = {name, val}; end
end

function [p, m] = passAndMeasured(F)
% "both": which file is the ff_cubesat_link result and which the ff_srs4_analyzer run
if numel(F) ~= 2, error('studio:file', 'Pass and Measured: give two files, cubesat_link.mat; srs4_run.mat'); end
isRun = false(1, 2);
for k = 1:2
    w = string({whos('-file', F(k)).name});
    isRun(k) = any(ismember(["run_data" "limits"], w));
end
if isequal(isRun, [true false]), F = F([2 1]);
elseif ~isequal(isRun, [false true])
    error('studio:file', 'Pass and Measured: one cubesat_link.mat (ff_cubesat_link) and one srs4_run.mat (ff_srs4_analyzer).');
end
p = F(1); m = F(2);
end

function c = withMore(c, extra)
% append the card's "More" Name-Value text inside the call: f(..., <extra>); (last value wins)
if strlength(extra) == 0, return, end
c = extractBefore(c, strlength(c) - 1) + ", " + extra + ");";
end

%% ---- MATLAB text of a call -------------------------------------------------------

function s = call(out, func, first, nv)
% 'out = func(first, Name, Value, ...);' skipping empty values
parts = strings(1, 0);
if strlength(first) > 0, parts(end + 1) = first; end
for k = 1:2:numel(nv)
    val = nv{k + 1};
    if isempty(val) || (isstring(val) && isscalar(val) && strlength(val) == 0), continue, end
    parts(end + 1) = "'" + nv{k} + "', " + lit(val); %#ok<AGROW>
end
s = func + "(" + strjoin(parts, ", ") + ");";
if strlength(out) > 0, s = out + " = " + s; end
end

function s = lit(v)
% MATLAB literal of a value
if isstruct(v)
    s = v.code;
elseif ischar(v)
    s = "'" + replace(string(v), "'", "''") + "'";
elseif isstring(v)
    q = """" + replace(v, """", """""") + """";
    s = q; if ~isscalar(q), s = "[" + strjoin(q, " ") + "]"; end
elseif islogical(v)
    t = ["false" "true"]; s = t(v + 1);
else
    s = numLit(v);
    if ~isscalar(v), s = "[" + strjoin(arrayfun(@numLit, v), " ") + "]"; end
end
end

function s = numLit(v)
s = string(sprintf('%.10g', v));
if ~isscalar(v) || ~isfinite(v), return, end
for k = [9 6 3]                                     % 1.8425e9, 95e6, 150e3 (but 5261, 1001)
    m = v / 10^k;
    e = string(sprintf('%.10g', m)) + "e" + k;
    if abs(v) >= 10^k && abs(m - round(m, 4)) < 1e-9*abs(m) && strlength(e) < strlength(s)
        s = e; return
    end
end
end

function c = code(text)
c = struct('code', string(text));
end

function a = portArg(x, name)
a = {};
if x.port ~= 5025, a = {name, x.port}; end
end

function v = finiteOr(v)
if ~isempty(v) && isinf(v), v = []; end
end

function v = positiveOr(v)
if ~isempty(v) && v <= 0, v = []; end
end

function d = toolDir(x, func)
d = fullfile(x.out, func);
end

%% ---- file helpers ----------------------------------------------------------------

function T = fileCheck(p)
% What a file holds and which suite tools can open it (README section 6)
p = string(p);
[~, name, ext] = fileparts(p); ext = lower(ext);
rows = strings(0, 2);                               % Tool, Command
note = "";
if isfolder(p)
    if isfile(fullfile(p, "emi_study.mat"))
        rows(end + 1, :) = ["ff_emi_study", "R = ff_emi_study([], 'Data', " + lit(p) + ");"];
        note = "ff_emi_study folder";
    end
    if ~isempty(dir(fullfile(p, "**", "dSTA00000")))
        rows(end + 1, :) = ["ff_instrument_record", "R = ff_instrument_record('parse', " + lit(p) + ", 'RBW', <RBW in Hz>);"];
        note = "Record/Playback session folder (give the RBW you used)";
    end
    if note == "", note = "Folder without an EMI study or a Record/Playback session"; end
elseif any(ext == [".s1p" ".s2p" ".s3p" ".s4p"])
    rows(end + 1, :) = ["ff_vna", "R = ff_vna([], 'File', " + lit(p) + ");"];
    if ext == ".s2p"
        rows(end + 1, :) = ["ff_gnss_margin", "R = ff_gnss_margin(""10.90.36.200"", 'Cable', " + lit(p) + ");  % as the antenna cable"];
        rows(end + 1, :) = ["ff_vna", "R = ff_vna(""10.90.36.200"", 'Fixture', " + lit(p) + ");  % remove it from a measurement"];
    end
    note = "Touchstone file (" + ext + ")";
elseif ext == ".mat"
    w = whos('-file', p); vars = string({w.name});
    if ismember("study", vars)
        rows(end + 1, :) = ["ff_emi_study", "R = ff_emi_study([], 'Data', " + lit(p) + ");"];
        note = "Saved ff_emi_study (re-analyse with other settings)";
    elseif all(ismember(["f" "P"], vars))
        S = load(p, 'f'); nf = numel(S.f);
        nSw = w(vars == "P").size(1);
        if nf == 1 && ismember("phase_deg", vars)
            rows(end + 1, :) = ["ff_vvm_logger", "R = ff_vvm_logger('File', " + lit(p) + ");"];
            note = sprintf("VVM log: %d readings at %s", nSw, fmtHz(S.f));
        elseif nf == 1
            rows(end + 1, :) = ["ff_power_meter", "R = ff_power_meter('File', " + lit(p) + ");"];
            note = sprintf("Power-meter log: %d readings at %s (not a spectrum)", nSw, fmtHz(S.f));
        else
            rbw = NaN;
            if ismember("settings", vars)
                st = load(p, 'settings');
                if isstruct(st.settings) && isfield(st.settings, 'rbw'), rbw = st.settings.rbw; end
            end
            note = sprintf("Spectrum recording: %d sweeps x %d points, %s - %s, RBW %s", nSw, nf, ...
                fmtHz(S.f(1)), fmtHz(S.f(end)), fmtHz(rbw));
            rows(end + 1, :) = ["ff_spectrum_pro", "R = ff_spectrum_pro("""", 'Replay', " + lit(p) + ");"];
            rows(end + 1, :) = ["ff_emi_study", "R = ff_emi_study([], 'Data', " + lit(p) + ");  % split OFF/ON at the largest time gap"];
            if isfinite(rbw)
                rows(end + 1, :) = ["ff_signal_id", "T = ff_signal_id(" + lit(p) + ");"];
            else
                rows(end + 1, :) = ["ff_signal_id", "S = load(" + lit(p) + "); T = ff_signal_id(S.f, S.P, <RBW in Hz>);  % no RBW stored"];
            end
            if S.f(1) <= 1575.42e6 - 1e6 && S.f(end) >= 1575.42e6 + 1e6 && S.f(end) - S.f(1) <= 60e6
                rows(end + 1, :) = ["ff_gnss_margin", "R = ff_gnss_margin("""", 'Spectrum', " + lit(p) + ");"];
            end
            if S.f(1) < 2290e6 && S.f(end) > 2200e6         % overlaps the SRS-4 downlink band
                rows(end + 1, :) = ["ff_srs4_analyzer", "R = ff_srs4_analyzer([], 'File', " + lit(p) + ");  % spectrum analysis of an SRS-4 recording"];
                rows(end + 1, :) = ["ff_5g_ntn", "R = ff_5g_ntn('Interference', " + lit(p) + ");  % coexistence with 5G NTN band n256"];
            end
            rows(end + 1, :) = ["fieldfox_live", "fieldfox_live(""10.90.36.200"", 'BaselineFile', " + lit(p) + ");  % same span on the FieldFox"];
            if ismember("S21_dB", vars), rows(end + 1, :) = ["ff_vna", "R = ff_vna([], 'File', " + lit(p) + ");"]; end
        end
    elseif ismember("sats", vars)
        note = "ff_gnss_margin results (load it; the recordings are the gnss_l1_*.mat files)";
    elseif ismember("MaskFailures", vars)
        note = "ff_spectrum_pro log (load it; the recording is the matching _rec.mat)";
    elseif all(ismember(["run_data" "results"], vars))
        note = "ff_srs4_analyzer run (raw sweeps + results; re-analyse it with new options)";
        rows(end + 1, :) = ["ff_srs4_analyzer", "R = ff_srs4_analyzer([], 'File', " + lit(p) + ", 'ExpectedFreq', <tx.freq in Hz>);"];
        rows(end + 1, :) = ["ff_ccsds_linksim", "L = ff_ccsds_linksim('Measured', " + lit(p) + ");  % BER/PER at the measured C/N0"];
        rows(end + 1, :) = ["ff_5g_ntn", "N = ff_5g_ntn('Interference', " + lit(p) + ");  % 5G NTN coexistence"];
    elseif isscalar(vars) && any(vars == ["R" "results"])
        S = load(p, vars); r = S.(vars);
        fn = strings(1, 0);
        if isstruct(r), fn = string(fieldnames(r)).'; end
        if isstruct(r) && all(ismember(["passes" "orbit"], fn))
            nm = "the satellite";
            if isstruct(r.orbit) && isfield(r.orbit, 'name'), nm = string(r.orbit.name); end
            note = sprintf("ff_cubesat_link result: %d passes of %s", height(r.passes), nm);
            rows(end + 1, :) = ["ff_ccsds_linksim", "L = ff_ccsds_linksim('Pass', " + lit(p) + ");  % PER and throughput over these passes"];
        elseif isstruct(r) && all(ismember(["curves" "summary"], fn))
            note = "ff_ccsds_linksim result (load it: R.curves, R.summary, R.pass, R.measured)";
        elseif isstruct(r) && all(ismember(["linkLevel" "timeline"], fn))
            note = "ff_5g_ntn run (BLER curves + pass); reuse the curves for fast what-if budgets";
            rows(end + 1, :) = ["ff_5g_ntn", "N = ff_5g_ntn('File', " + lit(p) + ", 'EIRPDensity', 28);  % what-if"];
            rows(end + 1, :) = ["ff_5g_ntn", "N = ff_5g_ntn('File', " + lit(p) + ", 'Interference', ""demo"");  % + coexistence"];
        else
            note = "MAT file with the variable " + vars + " (fields: " + strjoin(fn(1:min(6, end)), ", ") + "): load it";
        end
    else
        note = "Not a suite recording (variables: " + strjoin(vars, ", ") + ")";
    end
elseif ext == ".sta"
    note = "FieldFox state file: recall it on the FieldFox (Restore state), or set stateFile at the top of analyze_fieldfox.m to read its CAT trace and SA settings";
elseif any(ext == [".tle" ".txt"])
    L = strtrim(readlines(p));
    L1 = L(startsWith(L, "1 ") & strlength(L) >= 69);
    if isempty(L1)
        note = "Text file without TLE lines";
    else
        ids = unique(strtrim(extractBetween(L1, 3, 7)), 'stable');
        note = sprintf("TLE file: %d satellite(s), NORAD %s", numel(ids), strjoin(ids, ", "));
        rows(end + 1, :) = ["ff_cubesat_link", "R = ff_cubesat_link('TLE', " + lit(p) + ", 'NORAD', " + ids(1) + ");"];
        rows(end + 1, :) = ["ff_5g_ntn", "R = ff_5g_ntn('TLE', " + lit(p) + ");  % geometry of a 5G NTN pass on this orbit"];
    end
else
    note = "Unknown file type " + ext;
end
fprintf('%s: %s\n', name + ext, note);
if isempty(rows), rows = ["-" "-"]; end
T = table(rows(:, 1), rows(:, 2), repmat(note, size(rows, 1), 1), 'VariableNames', ["Tool" "Command" "Note"]);
end

function [flt, multi] = fileFilter(kind)
mat = {'*.mat', 'Recordings (*.mat)'};
switch kind
    case "snp", flt = {'*.s1p;*.s2p;*.s3p;*.s4p', 'Touchstone (*.s1p, *.s2p, ...)'};
    case "tle", flt = {'*.tle;*.txt', 'TLE files (*.tle, *.txt)'; '*.*', 'All files'};
    case "any", flt = {'*.mat;*.s1p;*.s2p;*.sta;*.tle', 'Suite files (*.mat, *.sNp, *.sta, *.tle)'; '*.*', 'All files'};
    otherwise, flt = mat;
end
multi = 'off';
if any(kind == ["mat2" "mat3" "pm" "rmat"]), multi = 'on'; end
end

function t = filePlaceholder(kind)
switch kind
    case "mat",  t = "recording .mat (Browse...)";
    case "mat2", t = "OFF recording [; ON recording]";
    case "mat3", t = "emi_study.mat, a study folder, or 1-3 recordings (OFF; ON; OFF)";
    case "snp",  t = "Touchstone .s1p / .s2p";
    case "dir",  t = "folder copied from the FieldFox (RecordPlayback)";
    case "any",  t = "file or folder to check";
    case "tle",  t = "TLE file, e.g. support\tle\thai_cubesats.tle";
    case "pm",   t = "cubesat_link.mat; srs4_run.mat (pick both)";
    case "imat", t = "Interferer = file: srs4_run.mat or a recording (else leave empty)";
    case "rmat", t = "ntn5g_run.mat [; srs4_run.mat or a recording for Interferer = file]";
    otherwise,   t = "(no file needed)";
end
end

function ensureDir(d)
if ~isfolder(d), mkdir(d); end
end

function copyCommand(cmd)
if strlength(cmd) == 0, return, end
try
    clipboard('copy', char(cmd));
catch
end
end

function s = optionText(opts)
known = ["112" "QuickCal"; "210" "VNA T/R"; "211" "2-port VNA"; "233" "spectrum analyzer"; ...
    "235" "preamp"; "236" "interference analyzer + Record/Playback"; "308" "vector voltmeter"; ...
    "310" "channel power meter"; "350" "real-time SA"; "351" "I/Q"; "361" "EMI/CISPR detectors"];
parts = strings(1, numel(opts));
for k = 1:numel(opts)
    i = find(known(:, 1) == opts(k), 1);
    parts(k) = opts(k);
    if ~isempty(i), parts(k) = opts(k) + " " + known(i, 2); end
end
s = strjoin(parts, ", ");
miss = known(ismember(known(:, 1), ["350" "351" "361"]) & ~ismember(known(:, 1), opts), :);
if ~isempty(miss), s = s + "   |   not installed: " + strjoin(miss(:, 1) + " " + miss(:, 2), ", "); end
end

function A = appList()
% MATLAB apps that go with the suite (launcher names checked in R2026a with which/exist)
a = @(fn, label, tbx, use) struct('fn', string(fn), 'label', string(label), 'tbx', string(tbx), 'use', string(use));
A = [a("satelliteLinkBudgetAnalyzer", "Satellite Link Budget Analyzer", "Satellite Comms", ...
        "uplink/downlink budget tables; compare with ff_cubesat_link R.budget")
    a("wirelessWaveformGenerator", "Wireless Waveform Generator", "Communications / 5G / LTE / WLAN", ...
        "5G NR, LTE, WLAN, DVB-S2, PSK/QAM waveforms with RRC filters and impairments, spectrum + constellation")
    a("rfBudgetAnalyzer", "RF Budget Analyzer", "RF", ...
        "cascade gain / noise figure / IP3 of a receive chain (e.g. antenna, LNA, cable, FieldFox)")
    a("antennaDesigner", "Antenna Designer", "Antenna", ...
        "design an S-band patch or a GNSS antenna; pattern, impedance, gain")
    a("matchingNetworkDesigner", "Matching Network Designer", "RF", ...
        "L/Pi/T matching from a measured S11 (ff_vna Touchstone file)")
    a("signalAnalyzer", "Signal Analyzer", "Signal Processing", ...
        "spectra, spectrograms and persistence of signals in the workspace (e.g. zero-span envelopes)")
    a("bertool", "Bit Error Rate Analysis (bertool)", "Communications", ...
        "theoretical BER curves (PSK, coding) next to ff_ccsds_linksim results")
    a("sensorArrayAnalyzer", "Sensor Array Analyzer", "Phased Array", ...
        "array patterns and steering, e.g. a satellite or ground-station array")];
end

function n = modeName(m)
names = dictionary(["SA" "NA" "CAT" "VVM" "CPM" "PM" "IA"], ["spectrum analyzer" "network analyzer" ...
    "cable and antenna test" "vector voltmeter" "channel power meter" "power meter" "interference analyzer"]);
n = "";
if isKey(names, m), n = names(m); end
end

function s = fmtHz(x)
if ~isfinite(x), s = "-"; return, end
u = ["Hz" "kHz" "MHz" "GHz"]; e = min(3, max(0, floor(log10(max(abs(x), 1))/3)));
s = string(sprintf('%.6g %s', x/1000^e, u(e + 1)));
end

function s = stamp()
s = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
end
