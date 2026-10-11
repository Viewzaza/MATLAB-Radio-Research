function R = ff_spectrum_pro(ip, varargin)
%FF_SPECTRUM_PRO  Controlled spectrum analysis of a FieldFox with MATLAB toolboxes.
%
%   R = ff_spectrum_pro                          FieldFox at 10.90.36.200, its current SA settings
%   R = ff_spectrum_pro(ip, Name, Value, ...)
%   R = ff_spectrum_pro("", 'Replay', "rec.mat")  same analysis on a saved recording (no instrument)
%   R = ff_spectrum_pro("", 'Replay', "demo")     synthetic example: an "inverter" switches on halfway
%
% What it does
%   Takes control of the FieldFox (its state is backed up, single-sweep mode),
%   applies the SA settings you give, then triggers every sweep itself
%   (INIT:IMM;*OPC?) and reads it in binary, so every trace is one complete, fresh
%   sweep. When it stops (window closed, MaxSweeps, Duration, error or Ctrl+C) the
%   instrument is put back exactly as it was. It never switches on the RF source.
%   Stop it by closing a window: then everything is saved. After Ctrl+C only the
%   recording is kept (it is saved every 2 minutes while running).
%
% Windows
%   Spectrum Analyzer (DSP System Toolbox spectrumAnalyzer): spectrum + spectrogram,
%     max hold, peak finder, channel measurement (occupied BW or ACPR) on 'Channel',
%     the baseline as a 2nd trace and the spectral mask with pass/fail.
%   Companion figure: wall-clock waterfall, labelled signals, measurement log.
%   Close either window to stop.
%
% Every CheckEvery sweeps
%   1. Signals: peaks of the recent average, labelled by ff_signal_id when it is
%      installed (band / technology guess), otherwise the scope's peak finder.
%   2. MATLAB check: bandpower / obw (Signal Processing Toolbox) on the same trace
%      the scope measured.
%   3. Instrument cross-check (live only): the FieldFox's own channel power or ACPR
%      (SENS:MEAS:CHAN, CALC:MEAS:DATA?) and its frequency counter on the strongest
%      peak (precise carrier frequency). Costs one extra sweep. A command the
%      FieldFox rejects (SYST:ERR?) is skipped from then on; a timeout switches
%      the cross-check off. The sweep loop pauses for the check (R.Log.Check_s).
%   All three go to R.Log and the log CSV.
%
% EMI baseline (e.g. solar inverter OFF, then ON)
%   'AutoBaseline', 20   the first 20 sweeps (inverter OFF) become the baseline.
%                        With Prompt=true you are asked to switch it OFF, then ON.
%   'BaselineFile', f    use the average of a recording made earlier (same span).
%   The mask is baseline + Threshold dB, or your own 'Limit'. The scope shows
%   pass/fail; every sweep that crosses the mask is logged in R.MaskFailures with
%   its wall-clock time, worst frequency and excess (dB) and the scope's own count
%   (getSpectralMaskStatus). R.MaskEvents holds the times of the scope's
%   MaskTestFailed events: they come in bursts, fewer than the failing sweeps.
%
% Name-Value options (defaults)
%   SA settings ([] = keep what the FieldFox has):
%   'Center' 'Span' 'Start' 'Stop'   []   Hz
%   'RBW'           []      Hz
%   'Points'        []      trace points (e.g. 401, 1001)
%   'Preamp'        []      true / false
%   'Atten'         []      dB
%   'Detector'      []      "POS", "AVER" (RMS, best for channel power), "SAMP", "NEG", "NORM", "AUTO"
%   'RefLevel'      []      dBm (top of the display)
%   'Average'       []      >1 = FieldFox trace average; each reading is a full set (slower)
%   Measurements:
%   'ChannelType'   "obw"   "obw" (occupied bandwidth) or "acpr"
%   'Channel'       []      [center width] in Hz; [] = SA center, span/4
%   'ACPROffset'    []      Hz from channel center to the adjacent channels; [] = width
%   'PercentOBW'    99
%   'NumPeaks'      8
%   'CheckEvery'    20      sweeps between signal labelling / log / cross-check
%   'CrossCheck'    true    run the instrument-side measurements
%   Baseline and mask:
%   'AutoBaseline'  0       sweeps averaged into the baseline (0 = none)
%   'BaselineFile'  ""      recording .mat (f, P) used as the baseline
%   'Threshold'     10      dB above the baseline average for the mask
%   'Limit'         []      own upper limit: dBm (scalar) or [f_Hz level_dBm] rows
%   Run:
%   'MaxSweeps'     Inf
%   'Duration'      Inf     seconds
%   'History'       300     sweeps in the waterfall
%   'Replay'        ""      recording .mat (or "demo") instead of the instrument
%   'ReplaySpeed'   Inf     1 = as recorded, Inf = as fast as possible
%   'Record'        true    save the sweeps (same format as fieldfox_live)
%   'Pictures'      true    save _report.png (spectrum, max hold, baseline, mask, peaks,
%                           waterfall) and, in an interactive session, the scope (_scope.png)
%   'OutDir'        pwd
%   'Visible'       'on'    'off' = no windows (live without MaxSweeps/Duration: 200 sweeps)
%   'Prompt'        true    false = never wait for the user
%   'Port'          5025
%
% Output R: f, P, t, settings, idn (the sweeps, as in a recording), Baseline, Mask,
%   MaskStatus, MaskFailures, MaskEvents, Log, Signals, SignalLog, Waterfall,
%   files, nSweeps, rate, stopReason.
%   R.Log columns: Scope_* = spectrumAnalyzer (ScopeSynced = it had shown this
%   sweep), ML_* = bandpower/obw on the same sweep, Inst_* = FieldFox's own result
%   (on the extra sweep, with ML_InstChanP_dBm computed from that sweep's trace),
%   Counter_Hz = frequency counter, Peak_minus_Counter_Hz = how far the trace point
%   is from the true carrier, Overload = ADC over-range (1: lower the input or add
%   attenuation), Note = instrument messages, SignalID_s / Check_s = time taken.
%
% Good to know
%   - Channel power from a trace is accurate with the AVER (RMS) detector and points
%     spaced no wider than the RBW. MATLAB treats the RBW as a rectangular noise
%     bandwidth, the FieldFox corrects for its Gaussian filter, so a few tenths of
%     a dB difference in the cross-check are normal.
%   - spectrumAnalyzer integrates a channel with the trapezoid rule (edge points
%     half weight), bandpower counts every point inside fully: for a flat channel
%     of n points they differ by 10*log10((n-1)/n) dB, but a carrier sitting on a
%     channel edge counts up to 3 dB less in the scope. Keep carriers inside.
%   - The FieldFox channel measurement is centred on the SA centre frequency, so the
%     cross-check uses SA centre +- width/2 even when 'Channel' is elsewhere.
%   - No quasi-peak / CISPR detectors (option 361 is not installed): treat results
%     as an indicative EMI pre-scan, not a compliance test. Below ~30 MHz an antenna
%     close to the inverter measures its near field.
%   - The spectrogram inside spectrumAnalyzer runs on simulated time; the companion
%     waterfall uses the real clock.
%
% Examples
%   R = ff_spectrum_pro("", 'Replay', "demo");            % try it without the FieldFox
%   R = ff_spectrum_pro("10.90.36.200", 'Center', 95e6, 'Span', 20e6, 'RBW', 100e3, ...
%       'Detector', "AVER", 'AutoBaseline', 20, 'Channel', [91e6 6e6]);   % inverter OFF, then ON
%   R2 = ff_spectrum_pro("", 'Replay', R.files(1));        % same analysis on that recording
%
% See also fieldfox_live, FieldFox, spectrumAnalyzer, bandpower, obw, ff_signal_id.

if nargin < 1, ip = ""; end
opt = struct('Center', [], 'Span', [], 'Start', [], 'Stop', [], 'RBW', [], 'Points', [], ...
    'Preamp', [], 'Atten', [], 'Detector', [], 'RefLevel', [], 'Average', [], ...
    'ChannelType', "obw", 'Channel', [], 'ACPROffset', [], 'PercentOBW', 99, 'NumPeaks', 8, ...
    'CheckEvery', 20, 'CrossCheck', true, 'AutoBaseline', 0, 'BaselineFile', "", 'Threshold', 10, ...
    'Limit', [], 'MaxSweeps', Inf, 'Duration', Inf, 'History', 300, 'Replay', "", 'ReplaySpeed', Inf, ...
    'Record', true, 'Pictures', true, 'OutDir', pwd, 'Visible', 'on', 'Prompt', true, 'Port', 5025);
if mod(numel(varargin), 2), error('Options come in Name, Value pairs.'); end
names = fieldnames(opt);
for k = 1:2:numel(varargin)
    hit = find(strcmpi(names, varargin{k}), 1);
    if isempty(hit), error('Unknown option "%s".', varargin{k}); end
    opt.(names{hit}) = varargin{k + 1};
end
opt.ChannelType = lower(string(opt.ChannelType));
if ~ismember(opt.ChannelType, ["obw" "acpr"]), error('ChannelType must be "obw" or "acpr".'); end
if ~isempty(opt.Channel) && numel(opt.Channel) ~= 2, error('Channel must be [center width] in Hz.'); end
opt.CheckEvery = max(1, round(opt.CheckEvery));
show = any(strcmpi(string(opt.Visible), ["on" "true" "1"]));
opt.Visible = matlab.lang.OnOffSwitchState(show);
offline = strlength(string(opt.Replay)) > 0;
if ~show && ~offline && isinf(opt.MaxSweeps) && isinf(opt.Duration)
    opt.MaxSweeps = 200;
    fprintf('Visible is off and no MaxSweeps/Duration given: stopping after 200 sweeps.\n');
end
if ~exist(opt.OutDir, 'dir'), mkdir(opt.OutDir); end

%% Source of the sweeps: the instrument or a recording jjjjjj
ff = []; ctrl = []; %#ok<NASGU>
if offline
    if strcmpi(string(opt.Replay), "demo")
        S = demoRecording();
        if opt.AutoBaseline == 0 && strlength(string(opt.BaselineFile)) == 0 && isempty(opt.Limit)
            opt.AutoBaseline = 60;
            fprintf('Demo: baseline = first 60 sweeps; the "inverter" switches on at sweep 121.\n');
        end
    else
        S = loadRecording(opt.Replay);
    end
    s = settingsFromRecording(S);
    idn = string(S.idn);
    repP = double(S.P); repT = S.t(:); iRep = 0;
    if numel(repT) ~= size(repP, 1), repT = NaT(size(repP, 1), 1); end
    fprintf('Replaying %d sweeps of %s\n', size(repP, 1), string(opt.Replay));
else
    if strlength(string(ip)) == 0, ip = "10.90.36.200"; end
    fprintf('Connecting to %s:%d ...\n', ip, opt.Port);
    ff = FieldFox(ip, opt.Port);
    fprintf('Connected: %s\n', ff.IDN);
    ff.takeControl();                          % backup .sta + single sweep
    ctrl = onCleanup(@() releaseInstrument(ff));
    cfg = saArgs(opt);
    if ~isempty(cfg), ff.saConfigure(cfg{:}); elseif ff.mode() ~= "SA", ff.setMode("SA"); end
    s = ff.saSettings();
    idn = ff.IDN;
end
s.f = s.f(:).';
s.binStep = (s.stop - s.start) / max(s.points - 1, 1);
s.center = (s.start + s.stop) / 2;
useAvg = ~offline && ~isempty(opt.Average) && opt.Average > 1;
if ~isempty(opt.Average), s.average = opt.Average; end
p = split(idn, ","); model = strtrim(p(min(2, numel(p))));
[fScale, fUnit] = freqUnit(s.stop);
ch = channelSetup();
hasSigId = exist('ff_signal_id', 'file') > 0;
crossOK = ~offline && opt.CrossCheck;
xChan = crossOK; xFco = crossOK; xOvl = crossOK;    % instrument features still in use

%% State
N        = s.points;
H        = nan(opt.History, N);          % waterfall, newest sweep in row 1
Ht       = NaT(opt.History, 1);
nH       = 0;
nSweeps  = 0;
rec      = struct('P', zeros(0, N), 't', NaT(0, 1), 'n', 0);
base     = struct('collecting', false, 'ready', false, 'acc', zeros(0, N), 'avg', [], 'max', [], 'source', "");
mask     = struct('table', [], 'lim', [], 'source', "");
maskFails0 = 0; scopeFails = 0; evPrev = 0;
failRows = struct([]); logRows = struct([]);
Signals = table(); SignalLog = table();
sigWarned = false; spanWarned = false;
rate = NaN; running = true; stopReason = ""; finished = false; seenVisible = false;
files = strings(0, 1); R = struct();
runStamp = stamp();
recFile = fullfile(opt.OutDir, "ff_spectrum_pro_" + runStamp + "_rec.mat");

if hasSigId     % its first call loads models (~10 s): do it now, not in the middle of the run
    fprintf('Loading ff_signal_id ...\n');
    try [~] = ff_signal_id(s.f, -100 + randn(1, N), s.rbw); catch, end
end

%% Spectrum Analyzer scope
% No nested-function callbacks on the windows, the listener or an onCleanup: they
% would keep this workspace alive after an error or Ctrl+C, so the onCleanup that
% restores the FieldFox would not run until MATLAB exits.
sa = makeScope();
evStore = containers.Map('KeyType', 'double', 'ValueType', 'any');   % MaskTestFailed times
lst = addlistener(sa.SpectralMask, 'MaskTestFailed', @(~, ~) logMaskEvent(evStore));

%% Companion figure: waterfall, signals, log (closing it stops the run)
fig = figure('Name', "FieldFox spectrum pro - " + model, 'NumberTitle', 'off', 'Color', 'w', ...
    'Position', [80 60 1250 780], 'Visible', opt.Visible);
axW = axes(fig, 'Position', [0.06 0.50 0.86 0.45]);
hW = imagesc(axW, s.f([1 end])/fScale, [1 opt.History], H, 'AlphaData', 0);
colormap(axW, turbo); cb = colorbar(axW); cb.Label.String = 'dBm';
xlabel(axW, "Frequency (" + fUnit + ")"); ylabel(axW, 'Time (newest at top)');
title(axW, 'Waterfall (wall clock)');
lims = yLimits();
clim(axW, lims);
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.01 0.385 0.48 0.03], ...
    'String', sprintf('Signals (average of the last %d sweeps)', opt.CheckEvery), ...
    'BackgroundColor', 'w', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
tblS = uitable(fig, 'Units', 'normalized', 'Position', [0.01 0.06 0.48 0.32], 'RowName', [], ...
    'ColumnName', {'Freq (MHz)', 'dBm', 'vs base', 'OBW (kHz)', 'Band', 'Technology', 'Score'}, ...
    'ColumnWidth', {75, 50, 55, 65, 130, 120, 45}, 'Data', cell(0, 7));
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.51 0.385 0.48 0.03], ...
    'String', 'Measurement log: scope vs MATLAB vs FieldFox', ...
    'BackgroundColor', 'w', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
tblL = uitable(fig, 'Units', 'normalized', 'Position', [0.51 0.06 0.48 0.32], 'RowName', [], ...
    'ColumnName', {'Time', 'Sweep', 'Scope dBm', 'MATLAB dBm', 'FieldFox dBm', 'FF-ML dB', ...
    'OBW (kHz)', 'Counter (MHz)'}, 'ColumnWidth', {60, 45, 65, 70, 75, 55, 65, 95}, 'Data', cell(0, 8));
status = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.01 0.005 0.98 0.04], ...
    'BackgroundColor', 'w', 'HorizontalAlignment', 'left', 'FontSize', 9, 'String', '');

try
    %% Baseline and mask
    if ~isempty(opt.Limit), setMask(opt.Limit, "user limit"); end
    if strlength(string(opt.BaselineFile)) > 0, loadBaseline(opt.BaselineFile); end
    if opt.AutoBaseline > 0 && ~base.ready
        if ~offline && ~ask(sprintf(['Baseline: set up the reference condition now (e.g. inverter OFF).\n' ...
                'The next %d sweeps are averaged into the baseline.'], opt.AutoBaseline))
            running = false; stopReason = "Cancelled";
        end
        base.collecting = true;
    end

    %% Sweep loop
    tStart = tic; tLast = tic; tSlow = tic; tSave = tic; saveEvery = 120;
    while running
        [y, tNow] = acquire();
        if isempty(y), stopReason = "End of replay"; break, end
        if numel(y) ~= N
            stopReason = "Frequency settings changed on the instrument"; break
        end
        processSweep(y, tNow);
        if mod(nSweeps, opt.CheckEvery) == 0, periodicCheck(y, tNow); end

        dt = toc(tLast); tLast = tic;
        if isnan(rate), rate = 1/dt; else, rate = 0.9*rate + 0.1/dt; end
        if toc(tSlow) > 0.25
            updateSlow();
            tSlow = tic;
        end
        if toc(tSave) > saveEvery                % so that Ctrl+C keeps what was recorded
            tSave = tic; saveRecording(); saveEvery = max(120, 30*toc(tSave)); tSave = tic;
        end
        drawnow limitrate                        % also delivers the MaskTestFailed events
        if show && ~isvalid(fig), stopReason = "Window closed"; running = false; end
        if nSweeps >= opt.MaxSweeps, stopReason = "MaxSweeps"; running = false; end
        if toc(tStart) >= opt.Duration, stopReason = "Duration"; running = false; end
    end
    if nSweeps > 0 && isempty(logRows), periodicCheck(H(1, :), Ht(1)); end
catch runErr
    stopReason = "Error: " + runErr.message;
    finish();                                    % save what we have; onCleanup restores the FieldFox
    rethrow(runErr);
end
finish();

%% ===== nested helpers (share the variables above) ==========================

    function [y, tNow] = acquire()
        if offline
            iRep = iRep + 1;
            if iRep > size(repP, 1), y = []; tNow = NaT; return, end
            if isfinite(opt.ReplaySpeed) && iRep > 1 && ~isnat(repT(iRep))
                pause(max(0, seconds(repT(iRep) - repT(iRep - 1))/opt.ReplaySpeed - toc(tLast)));
            end
            y = double(repP(iRep, :)); tNow = repT(iRep);
            if isnat(tNow), tNow = datetime('now'); end
        else
            if useAvg, ff.restartAverage(); else, ff.sweep(); end
            y = reshape(ff.saTrace(), 1, []);         % row, like a recording
            tNow = datetime('now');
        end
    end

    function processSweep(y, tNow)
        nSweeps = nSweeps + 1;
        H = [y; H(1:end - 1, :)];
        Ht = [tNow; Ht(1:end - 1)];
        nH = min(nH + 1, opt.History);
        if opt.Record && ~offline
            if rec.n == size(rec.P, 1)          % grow geometrically (long runs)
                m = size(rec.P, 1) + 1000; rec.P(2*m, N) = 0; rec.t(2*m, 1) = NaT;
            end
            rec.n = rec.n + 1; rec.P(rec.n, :) = y; rec.t(rec.n) = tNow;
        end
        if nSweeps == 1
            sa.PeakFinder.MinHeight = median(y) + 10;   % peaks 10 dB above the noise floor
        end
        if base.collecting
            base.acc(end + 1, :) = y;
            if size(base.acc, 1) >= opt.AutoBaseline, finishBaseline("first " + opt.AutoBaseline + " sweeps"); end
        end

        if base.ready, sa([y(:) base.avg(:)]); else, sa([y(:) nan(N, 1)]); end %#ok<NOEFF> (steps the scope)
        checkMask(y, tNow);
    end

    function checkMask(y, tNow)
        % MATLAB compares every sweep with the same mask (exact, cheap) and logs
        % each failing sweep. The scope's own count (getSpectralMaskStatus, ~10 ms,
        % can be one sweep behind) is polled on failures and every 10 sweeps.
        if isempty(mask.table), return, end
        over = y - mask.lim;
        mlFail = any(over > 0);
        if mlFail || mod(nSweeps, 10) == 0
            st = getSpectralMaskStatus(sa);
            if isfield(st, 'NumTotalTests')   % empty until the scope has tested the new mask
                scopeFails = st.NumTotalTests - st.NumPassedTests - maskFails0;
            end
        end
        if mlFail
            [ex, i] = max(over);
            row = struct('Time', tNow, 'Sweep', nSweeps, 'WorstFreq_Hz', s.f(i), 'Level_dBm', y(i), ...
                'Limit_dBm', mask.lim(i), 'Excess_dB', ex, 'BinsOver', nnz(over > 0), ...
                'ScopeFailCount', scopeFails, 'EventFired', evStore.Count > evPrev, 'Mask', mask.source);
            if isempty(failRows), failRows = row; else, failRows(end + 1) = row; end
        end
        evPrev = evStore.Count;
    end

    %% Measurements every CheckEvery sweeps
    function periodicCheck(y, tNow)
        tCheck = tic;
        row = struct('Time', tNow, 'Sweep', nSweeps, 'Chan_Hz', ch.center, 'ChanBW_Hz', ch.bw);
        % scope (DSP System Toolbox): it works asynchronously, so first wait until
        % it shows this sweep (typically < 0.1 s); then its results belong to it
        row.ScopeSynced = waitScope(y);
        md = getMeasurementsData(sa);
        c = cellOrValue(md, 'ChannelMeasurements');
        row.Scope_ChanP_dBm = getf(c, 'ChannelPower');
        row.Scope_OBW_Hz = getf(c, 'OccupiedBW');
        row.Scope_ACPRL_dB = getf(c, 'ACPRLower');
        row.Scope_ACPRU_dB = getf(c, 'ACPRUpper');
        % MATLAB (Signal Processing Toolbox) on the same sweep
        [row.ML_ChanP_dBm, row.ML_OBW_Hz, row.ML_ACPRL_dB, row.ML_ACPRU_dB] = mlChannel(y, ch.center, ch.bw);
        [row.Peak_dBm, iPk] = max(y);
        row.Peak_Hz = s.f(iPk);
        pk = cellOrValue(md, 'PeakFinder');
        row.ScopePeak_Hz = NaN;
        if isstruct(pk) && ~isempty(pk.Frequency), row.ScopePeak_Hz = pk.Frequency(1); end
        row = crossCheck(row);
        row.SignalID_s = labelSignals(pk, tNow);
        row.Check_s = toc(tCheck);                 % how long the sweep loop paused
        if isempty(logRows), logRows = row; else, logRows(end + 1) = row; end
        showLog();
    end

    function ok = waitScope(y)
        ok = false; t0 = tic;
        while toc(t0) < 0.5
            sd = getSpectrumData(sa);
            sp = sd.Spectrum{1};
            if size(sp, 1) == N && max(abs(sp(:, 1) - y(:))) < 1e-6, ok = true; return, end
            pause(0.005);
        end
    end

    function row = crossCheck(row)
        % FieldFox's own channel power / ACPR and frequency counter, one extra sweep.
        row.Inst_Center_Hz = NaN; row.Inst_IBW_Hz = NaN; row.Inst_ChanP_dBm = NaN;
        row.Inst_PSD_dBmHz = NaN; row.Inst_ACPRL_dBc = NaN; row.Inst_ACPRU_dBc = NaN;
        row.ML_InstChanP_dBm = NaN; row.Inst_minus_ML_dB = NaN;
        row.Marker_Hz = NaN; row.Marker_dBm = NaN; row.Counter_Hz = NaN; row.Peak_minus_Counter_Hz = NaN;
        row.Overload = NaN; row.Note = "";
        if ~crossOK, return, end
        ibw = min(ch.bw, s.stop - s.start);
        % Checked on this unit (A.12.45): SENS:MEAS:CHAN, SENS:CME:IBW, SENS:ACP:OFFS1:*,
        % CALC:MARK1:STAT NORM (ON is rejected), CALC:MARK1:FUNC:MAX, CALC:MARK:FCO ON,
        % CALC:MARK:FCO:X?. CALC:MEAS:DATA? follows FORM (binary block when REAL,32):
        % CHP -> [power dBm, PSD dBm/Hz], ACPR -> 21 values, OBW -> [OBW Hz, power dBm].
        % Each group is still checked with SYST:ERR? and switched off alone if
        % rejected; a timeout switches the cross-check off.
        step = "setup";
        try
            if xFco
                ff.write("CALC:MARK1:STAT NORM");
                ff.write("CALC:MARK1:FUNC:MAX");       % marker 1 to the peak of the last sweep
                ff.write("CALC:MARK:FCO ON");          % ...and counted during the next sweep
                xFco = accepted("frequency counter", "CALC:MARK:AOFF");
            end
            if xChan
                ff.write("SENS:MEAS:CHAN " + upper(replace(opt.ChannelType, "obw", "chp")));
                ff.write(sprintf("SENS:CME:IBW %.10g", ibw));
                if opt.ChannelType == "acpr"
                    ff.write(sprintf("SENS:ACP:OFFS1:FREQ %.10g", ch.off));
                    ff.write(sprintf("SENS:ACP:OFFS1:BWID %.10g", ibw));
                    ff.write("SENS:ACP:OFFS1:STAT ON");
                end
                xChan = accepted("channel measurement", "SENS:MEAS:CHAN NONE");
                if xChan && ~keepSettings()        % the measurement must not move the span
                    ff.write("SENS:MEAS:CHAN NONE"); keepSettings(); xChan = false;
                    row.Note = "channel measurement changes span/RBW: skipped from now on";
                    warning('ff_spectrum_pro:crossCheck', 'The FieldFox channel measurement changes the span or RBW: skipped from now on.');
                end
            end
            if ~xFco && ~xChan, crossOK = false; return, end

            step = "extra sweep";
            if useAvg, ff.restartAverage(); else, ff.sweep(); end
            yc = reshape(ff.saTrace(), 1, []);
            if numel(yc) ~= N, error('extra sweep has %d points instead of %d', numel(yc), N); end
            if xChan
                step = "CALC:MEAS:DATA?";
                v = ff.queryBlock("CALC:MEAS:DATA?", "single", 5);
                v(v < -9e4) = NaN;                 % undefined offsets: -9.87e4 placeholders
                if numel(v) >= 2
                    row.Inst_ChanP_dBm = v(1); row.Inst_PSD_dBmHz = v(2);
                end
                if numel(v) >= 9                   % ACPR: offset 1 lower/upper rel dBc
                    row.Inst_ACPRL_dBc = v(6); row.Inst_ACPRU_dBc = v(9);
                elseif opt.ChannelType == "acpr"
                    row.Note = "CALC:MEAS:DATA? returned " + numel(v) + " values (21 expected for ACPR)";
                end
                row.Inst_Center_Hz = s.center; row.Inst_IBW_Hz = ibw;
                row.ML_InstChanP_dBm = mlChannel(yc, s.center, ibw);
                row.Inst_minus_ML_dB = row.Inst_ChanP_dBm - row.ML_InstChanP_dBm;
                ff.write("SENS:MEAS:CHAN NONE");
            end
            if xFco
                step = "marker / counter query";
                row.Marker_Hz = scalarNum(ff.queryNum("CALC:MARK1:X?", 5));
                row.Marker_dBm = scalarNum(ff.queryNum("CALC:MARK1:Y?", 5));
                row.Counter_Hz = scalarNum(ff.queryNum("CALC:MARK:FCO:X?", 5));
                [~, iPk] = max(yc);
                row.Peak_minus_Counter_Hz = s.f(iPk) - row.Counter_Hz;
                ff.write("CALC:MARK:FCO OFF");
                ff.write("CALC:MARK:AOFF");
            end
            if xOvl
                try
                    row.Overload = scalarNum(ff.queryNum("CALC:MEAS:WAOR?", 5));   % 1 = ADC over-range: reduce input
                catch
                    row.Overload = NaN; resync(ff);
                end
                if isnan(row.Overload), xOvl = false; end
            end
            step = "*OPC? / SYST:ERR?";
            ff.opc(30);
            e = ff.errors();
            if ~isempty(e)
                row.Note = strjoin([row.Note; e(:)], "; ");
                warning('ff_spectrum_pro:crossCheck', 'Instrument cross-check reported: %s', strjoin(e, "; "));
            end
            if xChan, keepSettings(); end
        catch err
            crossOK = false;
            msg = err.message;
            if contains(msg, "Input should be a string"), msg = "no reply (timeout)"; end   % empty readline
            row.Note = "cross-check failed at " + step + ": " + msg;
            warning('ff_spectrum_pro:crossCheck', 'Instrument cross-check disabled (%s: %s).', step, msg);
            resync(ff);                            % a late reply must not shift later answers
            try ff.write("CALC:MARK:AOFF"); ff.write("SENS:MEAS:CHAN NONE"); ff.errors(); keepSettings(); catch, end
        end

        function ok = accepted(what, undo)
            e = ff.errors();
            ok = isempty(e);
            if ~ok
                ff.write(undo); ff.errors();
                row.Note = strjoin([row.Note; what + ": " + e(:)], "; ");
                row.Note = regexprep(row.Note, "^; ", "");
                warning('ff_spectrum_pro:crossCheck', 'FieldFox rejected the %s commands (%s): skipped from now on.', what, strjoin(e, '; '));
            end
        end
    end

    function ok = keepSettings()
        % If switching the channel measurement moved span/points/RBW, put them back.
        moved = @(n) abs(n.start - s.start) > s.binStep/2 || abs(n.stop - s.stop) > s.binStep/2 || ...
            n.points ~= N || abs(n.rbw - s.rbw) > 1e-3*s.rbw;
        ok = ~moved(ff.saSettings());
        if ok, return, end
        if ~spanWarned
            warning('ff_spectrum_pro:settings', 'Channel measurement changed span/points/RBW; restoring them.');
            spanWarned = true;
        end
        ff.saConfigure(Start=s.start, Stop=s.stop, Points=N, RBW=s.rbw);
        ok = ~moved(ff.saSettings());
    end

    function tId = labelSignals(pk, tNow)
        tId = NaN;
        nA = max(1, min(opt.CheckEvery, nH));
        Pavg = 10*log10(mean(10.^(H(1:nA, :)/10), 1));    % power average, not dB
        T = table();
        if hasSigId
            try
                t0 = tic;
                idDet = "auto";                         % the template views depend on the detector
                if isfield(s, 'detector') && strlength(string(s.detector)) > 0, idDet = string(s.detector); end
                T = ff_signal_id(s.f, Pavg, s.rbw, 'Detector', idDet);
                tId = toc(t0);
            catch err
                if ~sigWarned
                    warning('ff_spectrum_pro:signalId', 'ff_signal_id failed (%s); using the peak finder.', err.message);
                    sigWarned = true;
                end
                T = table();
            end
        end
        if isempty(T) || ~istable(T) || ~ismember('Freq_Hz', T.Properties.VariableNames)
            T = peakTable(pk);
        end
        T.VsBaseline_dB = nan(height(T), 1);
        if base.ready && height(T) > 0
            T.VsBaseline_dB = T.Level_dBm - interp1(s.f, base.avg, T.Freq_Hz, 'nearest', NaN);
        end
        Signals = T;
        keep = intersect({'Freq_Hz', 'Level_dBm', 'OBW_Hz', 'BW26_Hz', 'Band', 'Channel', 'Technology', ...
            'Score', 'Top3', 'Method', 'VsBaseline_dB'}, T.Properties.VariableNames, 'stable');
        L = addvars(T(:, keep), repmat(tNow, height(T), 1), repmat(nSweeps, height(T), 1), ...
            'Before', 1, 'NewVariableNames', {'Time', 'Sweep'});
        try
            SignalLog = [SignalLog; L];
        catch err
            warning('ff_spectrum_pro:signalLog', 'Signal log row skipped: %s', err.message);
        end
        if ~isvalid(fig), return, end
        rows = cell(height(T), 7);
        for i = 1:height(T)
            rows(i, :) = {sprintf('%.4f', T.Freq_Hz(i)/1e6), sprintf('%.1f', T.Level_dBm(i)), ...
                sprintf('%+.1f', T.VsBaseline_dB(i)), sprintf('%.1f', T.OBW_Hz(i)/1e3), ...
                cellText(T.Band(i)), cellText(T.Technology(i)), sprintf('%.2f', T.Score(i))};
        end
        tblS.Data = rows;
    end

    function T = peakTable(pk)
        % Same columns as ff_signal_id, from the scope's peak finder.
        if isstruct(pk) && ~isempty(pk.Frequency), fr = pk.Frequency(:); lv = pk.Value(:);
        else, fr = zeros(0, 1); lv = zeros(0, 1);
        end
        n = numel(fr); e = strings(n, 1);
        T = table(fr, lv, nan(n, 1), nan(n, 1), e, e, repmat("Unknown", n, 1), nan(n, 1), e, ...
            repmat("spectrumAnalyzer PeakFinder", n, 1), 'VariableNames', ...
            {'Freq_Hz', 'Level_dBm', 'OBW_Hz', 'BW26_Hz', 'Band', 'Channel', 'Technology', 'Score', 'Top3', 'Method'});
    end

    function showLog()
        if ~isvalid(fig), return, end
        L = logRows(end:-1:max(1, end - 14));          % newest first
        rows = cell(numel(L), 8);
        for i = 1:numel(L)
            obwHz = L(i).Scope_OBW_Hz; if isnan(obwHz), obwHz = L(i).ML_OBW_Hz; end
            rows(i, :) = {char(string(L(i).Time, 'HH:mm:ss')), L(i).Sweep, sprintf('%.2f', L(i).Scope_ChanP_dBm), ...
                sprintf('%.2f', L(i).ML_ChanP_dBm), sprintf('%.2f', L(i).Inst_ChanP_dBm), ...
                sprintf('%+.2f', L(i).Inst_minus_ML_dB), sprintf('%.1f', obwHz/1e3), sprintf('%.6f', L(i).Counter_Hz/1e6)};
        end
        tblL.Data = rows;
    end

    function [chanP, obwHz, acprL, acprU] = mlChannel(y, fc, bw)
        % Channel power, occupied BW and ACPR with Signal Processing Toolbox functions.
        % sxx = power per trace point (mW); PSD = sxx/RBW (RBW as rectangular noise
        % BW). A channel = the trace points inside it, each counted fully
        % (spectrumAnalyzer gives its two edge points half weight: trapezoid rule).
        f = s.f(:); sxx = 10.^(y(:)/10);
        bp = @(r) 10*log10(bandpower(sxx/s.rbw, f, r, "psd"));
        r = binRange(fc, bw);
        chanP = NaN; obwHz = NaN; acprL = NaN; acprU = NaN;
        if isempty(r), return, end
        chanP = bp(r);
        obwHz = obw(sxx, f, max(s.rbw, s.binStep), r, opt.PercentOBW);   % obw rejects RBW < point spacing (wide spans)
        if fc - ch.off - bw/2 >= f(1), rl = binRange(fc - ch.off, bw); if ~isempty(rl), acprL = bp(rl) - chanP; end, end
        if fc + ch.off + bw/2 <= f(end), ru = binRange(fc + ch.off, bw); if ~isempty(ru), acprU = bp(ru) - chanP; end, end
    end

    function r = binRange(c, w)
        % first and last trace frequency inside c +- w/2
        tol = 1e-6*s.binStep;
        i1 = find(s.f >= c - w/2 - tol, 1); i2 = find(s.f <= c + w/2 + tol, 1, 'last');
        if isempty(i1) || isempty(i2) || i2 <= i1, r = []; else, r = s.f([i1 i2]); end
    end

    %% Baseline and mask
    function finishBaseline(src)
        base.avg = 10*log10(mean(10.^(base.acc/10), 1));
        base.max = max(base.acc, [], 1);
        base.collecting = false; base.ready = true; base.source = src;
        sa.LineVisible = {'on', 'on'};               % 2nd scope channel: baseline
        sa.ShowLegend = true;
        if isempty(opt.Limit), setMask([s.f(:) base.avg(:) + opt.Threshold], "baseline + " + opt.Threshold + " dB"); end
        fprintf('Baseline ready (%s).\n', src);
        if ~offline && startsWith(src, "first") && ~ask('Baseline done. Now switch the inverter (device under test) ON.')
            running = false; stopReason = "Cancelled";
        end
    end

    function loadBaseline(file)
        B = loadRecording(file);
        if numel(B.f) ~= N || abs(B.f(1) - s.start) > s.binStep/2 || abs(B.f(end) - s.stop) > s.binStep/2
            warning('ff_spectrum_pro:baseline', ['Baseline file %s has a different frequency axis ' ...
                '(%s - %s, %d points); not used.'], file, fmtHz(B.f(1)), fmtHz(B.f(end)), numel(B.f));
            return
        end
        base.acc = double(B.P);
        [~, nm, ext] = fileparts(string(file));
        finishBaseline("file " + nm + ext + " (" + size(B.P, 1) + " sweeps)");
    end

    function setMask(M, src)
        if isscalar(M)
            lim = repmat(M, 1, N);
        else
            if size(M, 2) ~= 2 || any(diff(M(:, 1)) < 0), error('Limit must be a scalar or [f_Hz level_dBm] rows with rising f.'); end
            lim = maskLevels(M, s.f);
        end
        sa.SpectralMask.EnabledMasks = "upper";
        sa.SpectralMask.ReferenceLevel = "custom";
        sa.SpectralMask.CustomReferenceLevel = 0;
        sa.SpectralMask.SelectedChannel = 1;       % test the live trace only
        sa.SpectralMask.UpperMask = M;
        mask.table = M; mask.lim = lim; mask.source = src;
        try                                        % scope failures before this mask
            st = getSpectralMaskStatus(sa);
            maskFails0 = st.NumTotalTests - st.NumPassedTests;
        catch
            maskFails0 = 0;
        end
        scopeFails = 0;
    end

    %% Display
    function sa = makeScope()
        df = s.binStep;
        % Exact bin mapping for odd and even N (verified in R2026a):
        % SampleRate = N*df, FrequencyOffset = start + floor((N-1)/2)*df.
        sa = spectrumAnalyzer(InputDomain="frequency", SampleRate=N*df, ...
            FrequencyOffset=s.start + floor((N - 1)/2)*df, ...
            InputUnits="dBm", SpectrumUnits="dBm", ViewType="spectrum-and-spectrogram", ...
            PlotMaxHoldTrace=true, ChannelNames={'Live', 'Baseline'}, LineVisible={'on', 'off'}, YLimits=yLimits(), ...
            Title=sprintf('%s  %s - %s  RBW %s', model, fmtHz(s.start), fmtHz(s.stop), fmtHz(s.rbw)));
        % Always 2 channels [live baseline]; the baseline is NaN (hidden) until it
        % exists, because changing the channel count needs release() (~12 s).
        sa.Visible = show;
        sa.Name = "FieldFox spectrum pro";
        sa.AxesScaling = "manual";                 % keep the FieldFox display range
        sa.ColorLimits = yLimits();                % spectrogram colours (default -80..20 dBm)
        if isfinite(s.rbw) && s.rbw > 0 && s.rbw <= N*df/2
            sa.RBWSource = "property";           % channel power scaled by bin spacing / RBW
            sa.RBW = s.rbw;
        else
            warning('ff_spectrum_pro:rbw', 'RBW %g Hz does not fit the span; scope channel power is a plain bin sum.', s.rbw);
        end
        sa.PeakFinder.Enabled = true;
        sa.PeakFinder.NumPeaks = opt.NumPeaks;
        sa.PeakFinder.LabelPeaks = true;
        sa.ChannelMeasurements.Enabled = true;
        sa.ChannelMeasurements.Type = replace(opt.ChannelType, "obw", "occupied-bandwidth");
        sa.ChannelMeasurements.CenterFrequency = ch.center;
        sa.ChannelMeasurements.Span = ch.bw;
        if opt.ChannelType == "obw"
            sa.ChannelMeasurements.PercentOccupiedBW = opt.PercentOBW;
        else
            sa.ChannelMeasurements.NumOffsets = 1;
            sa.ChannelMeasurements.ACPROffsets = ch.off;
            sa.ChannelMeasurements.AdjacentBW = ch.bw;
        end
    end

    function lim = yLimits()
        if isfield(s, 'ref') && isfield(s, 'pdiv') && all(isfinite([s.ref s.pdiv])) && s.pdiv > 0
            lim = [s.ref - 10*s.pdiv, s.ref];
        elseif offline
            lim = [floor(min(repP(:))) - 5, ceil(max(repP(:))) + 5];
        else
            lim = [-120 0];
        end
    end

    function c = channelSetup()
        c.center = s.center; c.bw = (s.stop - s.start)/4;
        if ~isempty(opt.Channel), c.center = opt.Channel(1); c.bw = opt.Channel(2); end
        c.off = c.bw;
        if ~isempty(opt.ACPROffset), c.off = opt.ACPROffset; end
        if c.center - c.bw/2 < s.start || c.center + c.bw/2 > s.stop
            error('Channel %s +- %s is outside the sweep %s - %s.', fmtHz(c.center), fmtHz(c.bw/2), fmtHz(s.start), fmtHz(s.stop));
        end
        if opt.ChannelType == "acpr" && (c.center - c.off - c.bw/2 < s.start || c.center + c.off + c.bw/2 > s.stop)
            warning('ff_spectrum_pro:acpr', 'Adjacent channels reach outside the sweep: ACPR will be partly undefined.');
        end
    end

    function updateSlow()
        if show && isvalid(sa)
            v = isVisible(sa);
            if v, seenVisible = true; elseif seenVisible, stopReason = "Scope closed"; running = false; end
        end
        if ~isvalid(fig), return, end
        set(hW, 'CData', H, 'AlphaData', ~isnan(H));
        timeTicks(axW);
        txt = sprintf('%.1f sweeps/s   #%d', rate, nSweeps);
        if base.collecting, txt = sprintf('%s   BASELINE %d/%d', txt, size(base.acc, 1), opt.AutoBaseline); end
        if base.ready, txt = sprintf('%s   baseline: %s', txt, base.source); end
        if ~isempty(mask.table), txt = sprintf('%s   mask (%s): %d failing sweeps', txt, mask.source, numel(failRows)); end
        if ~crossOK && ~offline && opt.CrossCheck, txt = [txt '   cross-check OFF']; end
        status.String = txt;
    end

    function timeTicks(ax)
        % wall-clock labels on the waterfall rows
        idx = unique(round(linspace(1, max(nH, 1), 6)));
        idx = idx(~isnat(Ht(idx)));
        if ~isempty(idx), set(ax, 'YTick', idx, 'YTickLabel', cellstr(string(Ht(idx), 'HH:mm:ss'))); end
    end

    function reportFigure(F)
        % Picture for a report from the scope's own data (getSpectrumData also works
        % without a display) plus the wall-clock waterfall, drawn into the hidden F.
        tl = tiledlayout(F, 2, 1, 'TileSpacing', 'compact');
        ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
        sd = getSpectrumData(sa);
        fv = sd.FrequencyVector{1}/fScale;
        if ismember('MaxHoldTrace', sd.Properties.VariableNames)
            plot(ax, fv, sd.MaxHoldTrace{1}(:, 1), 'Color', [0.95 0.6 0.6], 'DisplayName', 'Max hold');
        end
        plot(ax, fv, sd.Spectrum{1}(:, 1), 'Color', [0 0.45 0.74], 'DisplayName', 'Last sweep');
        if base.ready, plot(ax, s.f/fScale, base.avg, 'k', 'DisplayName', "Baseline: " + base.source); end
        if ~isempty(mask.table)
            plot(ax, s.f/fScale, mask.lim, 'r--', 'LineWidth', 1.2, 'DisplayName', "Mask: " + mask.source);
        end
        xline(ax, (ch.center + [-1 1]*ch.bw/2)/fScale, ':', 'Color', [0.2 0.6 0.2], 'HandleVisibility', 'off');
        pk = cellOrValue(getMeasurementsData(sa), 'PeakFinder');
        if isstruct(pk) && ~isempty(pk.Frequency)
            plot(ax, pk.Frequency/fScale, pk.Value, 'kv', 'MarkerFaceColor', 'y', 'DisplayName', 'Peaks');
        end
        xlim(ax, s.f([1 end])/fScale); ylim(ax, yLimits());
        xlabel(ax, "Frequency (" + fUnit + ")"); ylabel(ax, 'dBm'); legend(ax, 'Location', 'northeast');
        title(ax, sprintf('%s   RBW %s   %d sweeps   %d failing the mask   (%s)', model, fmtHz(s.rbw), ...
            nSweeps, numel(failRows), stopReason), 'Interpreter', 'none');
        ax = nexttile(tl);
        imagesc(ax, s.f([1 end])/fScale, [1 opt.History], H, 'AlphaData', ~isnan(H));
        colormap(ax, turbo); clim(ax, yLimits()); c = colorbar(ax); c.Label.String = 'dBm';
        timeTicks(ax);
        xlabel(ax, "Frequency (" + fUnit + ")"); title(ax, 'Waterfall (wall clock, newest at top)');
    end

    function ok = ask(msg)
        ok = true;
        fprintf('%s\n', msg);
        if ~opt.Prompt, return, end
        if show
            ok = strcmp(questdlg(msg, 'ff_spectrum_pro', 'Continue', 'Stop', 'Continue'), 'Continue');
        else
            ok = ~strcmpi(strtrim(input('Enter = continue, s = stop: ', 's')), 's');
        end
    end

    %% Saving and result
    function saveRecording()
        % The sweeps in the fieldfox_live format (f, P, t, settings, idn).
        if offline || ~opt.Record || rec.n == 0, return, end
        settings = rmfield(s, 'f'); settings.tool = "ff_spectrum_pro";
        out = struct('f', s.f, 'P', rec.P(1:rec.n, :), 't', rec.t(1:rec.n), 'settings', settings, 'idn', idn);
        if rec.n*N*8 > 1.5e9, save(recFile, '-struct', 'out', '-v7.3'); else, save(recFile, '-struct', 'out'); end
    end

    function finish()
        if finished, return, end
        finished = true;
        try delete(lst); catch, end
        if stopReason == "", stopReason = "Interrupted"; end
        base.acc = [];
        settings = rmfield(s, 'f'); settings.tool = "ff_spectrum_pro";
        Log = struct2table(logRows, 'AsArray', true);
        MaskFailures = struct2table(failRows, 'AsArray', true);
        if ~isempty(failRows), MaskFailures.Mask = string(MaskFailures.Mask); end
        try maskStatus = getSpectralMaskStatus(sa); catch, maskStatus = struct(); end
        maskEvents = values(evStore); maskEvents = vertcat(NaT(0, 1), maskEvents{:});
        Waterfall = struct('P', H(1:nH, :), 't', Ht(1:nH));
        Baseline = rmfield(base, {'collecting', 'acc'}); Baseline.f = s.f;
        b = fullfile(opt.OutDir, "ff_spectrum_pro_" + runStamp);
        try
            saveRecording();
            if isfile(recFile), files(end + 1, 1) = recFile; end
            Mask = mask.table; Options = opt;
            save(b + "_log.mat", 'Log', 'MaskFailures', 'maskEvents', 'Signals', 'SignalLog', ...
                'Baseline', 'Mask', 'maskStatus', 'settings', 'idn', 'Options', 'stopReason');
            files(end + 1, 1) = b + "_log.mat";
            if ~isempty(logRows), writetable(Log, b + "_log.csv"); files(end + 1, 1) = b + "_log.csv"; end
            if ~isempty(failRows), writetable(MaskFailures, b + "_maskfail.csv"); files(end + 1, 1) = b + "_maskfail.csv"; end
            if height(SignalLog) > 0, writetable(SignalLog, b + "_signals.csv"); files(end + 1, 1) = b + "_signals.csv"; end
        catch err
            warning('ff_spectrum_pro:save', 'Saving failed: %s', err.message);
        end
        F = [];
        try                                        % pictures for a report
            if opt.Pictures && nSweeps > 0
                F = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 780]);
                reportFigure(F);
                exportgraphics(F, b + "_report.png", 'Resolution', 120); close(F);
                files(end + 1, 1) = b + "_report.png";
                % printToFigure copies what the scope has drawn: empty in batch/headless runs
                if show && ~batchStartupOptionUsed && isvalid(sa) && isVisible(sa)
                    F = printToFigure(sa, Visible=false);
                    if ~isempty(F), exportgraphics(F, b + "_scope.png", 'Resolution', 120); close(F); files(end + 1, 1) = b + "_scope.png"; end
                end
            end
        catch err
            warning('ff_spectrum_pro:export', 'Picture export failed: %s', err.message);
            if ~isempty(F) && isvalid(F), delete(F); end    % no hidden figure left behind
        end
        if offline, P = repP(1:min(iRep, end), :); t = repT(1:min(iRep, end));
        else, P = rec.P(1:rec.n, :); t = rec.t(1:rec.n);
        end
        R = struct('f', s.f, 'P', P, 't', t, 'settings', settings, 'idn', idn, ...
            'Baseline', Baseline, 'Mask', mask.table, 'MaskStatus', maskStatus, ...
            'MaskFailures', MaskFailures, 'MaskEvents', maskEvents, 'Log', Log, ...
            'Signals', Signals, 'SignalLog', SignalLog, 'Waterfall', Waterfall, ...
            'files', files, 'nSweeps', nSweeps, 'rate', rate, 'stopReason', stopReason);
        if ~show
            try delete(sa); catch, end
            if isvalid(fig), delete(fig); end
        elseif isvalid(fig)
            updateSlow();
        end
        fprintf('Stopped (%s) after %d sweeps, %.1f sweeps/s.\n', stopReason, nSweeps, rate);
        if ~isempty(mask.table) && isfield(maskStatus, 'NumTotalTests')
            fprintf('Mask "%s": %d failing sweeps (MATLAB), %d (spectrumAnalyzer), listener events %d.\n', ...
                mask.source, numel(failRows), maskStatus.NumTotalTests - maskStatus.NumPassedTests - maskFails0, evStore.Count);
        end
        for i = 1:numel(files), fprintf('  saved %s\n', files(i)); end
    end
end

%% ===== local functions =====================================================

function c = saArgs(opt)
% Name-value list for FieldFox.saConfigure with only the options given.
c = {};
for nm = ["Start" "Stop" "Center" "Span" "Points" "RBW" "RefLevel" "Atten" "Preamp" "Detector" "Average"]
    v = opt.(nm);
    if ~isempty(v), c(end + 1:end + 2) = {char(nm), v}; end
end
if ~isempty(opt.Average)
    if opt.Average > 1, c(end + 1:end + 2) = {'TraceType', "AVG"}; else, c(end + 1:end + 2) = {'TraceType', "CLRW"}; end
end
end

function x = scalarNum(v)
% First number of a query reply, NaN if there is none.
if isempty(v), x = NaN; else, x = v(1); end
end

function resync(ff)
% After a timeout the reply may still arrive later and would be read as the answer
% to the next query. Drop pending input, then read until our own *IDN? comes back.
try
    pause(1);
    flush(ff.Client, "input");
    ff.write("*IDN?");
    t0 = tic;
    while toc(t0) < 20
        r = readline(ff.Client);
        if ~isempty(r) && strtrim(r) == ff.IDN, return, end
    end
    warning('ff_spectrum_pro:resync', 'No reply from the FieldFox after a timeout.');
catch err
    warning('ff_spectrum_pro:resync', 'Could not resynchronise with the FieldFox: %s', err.message);
end
end

function logMaskEvent(store)
% MaskTestFailed listener: wall-clock time of each event (store is a handle).
store(store.Count + 1) = datetime('now'); %#ok<NASGU> (handle object)
end

function releaseInstrument(ff)
try
    ff.release();
    fprintf('FieldFox restored to its previous state.\n');
catch err
    warning('ff_spectrum_pro:release', 'Could not restore the FieldFox: %s', err.message);
end
end

function S = loadRecording(file)
% Recording .mat (f, P, t, settings, idn); only f and P are required.
w = warning('off', 'MATLAB:load:variableNotFound');
S = load(file, 'f', 'P', 't', 'settings', 'idn');
warning(w);
if ~isfield(S, 'f') || ~isfield(S, 'P'), error('%s has no f and P (not a recording).', file); end
if numel(S.f) < 5
    error(['%s is not a swept spectrum (%d frequency point): power meter and VVM logs are ' ...
        're-analysed with their own tool (''File'').'], file, numel(S.f));
end
if size(S.P, 2) ~= numel(S.f), S.P = S.P.'; end
if size(S.P, 2) ~= numel(S.f), error('%s: P does not match f (%d points).', file, numel(S.f)); end
if ~isfield(S, 't') || ~isdatetime(S.t), S.t = NaT(size(S.P, 1), 1); end
if ~isfield(S, 'idn'), S.idn = "unknown instrument"; end
end

function s = settingsFromRecording(S)
f = double(S.f(:).');
s = struct('start', f(1), 'stop', f(end), 'points', numel(f), 'rbw', NaN, 'vbw', NaN, ...
    'ref', NaN, 'pdiv', NaN, 'f', f);
if isfield(S, 'settings') && isstruct(S.settings)
    for nm = ["rbw" "vbw" "ref" "pdiv" "atten" "preamp" "detector"]
        if isfield(S.settings, nm), s.(nm) = S.settings.(nm); end
    end
end
if ~isfinite(s.rbw)
    s.rbw = (f(end) - f(1))/max(numel(f) - 1, 1);
    warning('ff_spectrum_pro:rbw', 'Recording has no RBW; using the point spacing (%g Hz).', s.rbw);
end
end

function lim = maskLevels(M, f)
% Upper-limit level at each frequency of a [f level] table (steps allowed).
lim = nan(size(f));
for k = 1:size(M, 1) - 1
    if M(k + 1, 1) <= M(k, 1), continue, end
    in = f >= M(k, 1) & f <= M(k + 1, 1);
    v = M(k, 2) + (f(in) - M(k, 1)) * (M(k + 1, 2) - M(k, 2)) / (M(k + 1, 1) - M(k, 1));
    lim(in) = max(lim(in), v);
end
end

function v = cellOrValue(md, name)
v = [];
if ismember(name, md.Properties.VariableNames)
    v = md.(name);
    if iscell(v), v = v{1}; end
end
end

function x = getf(c, name)
x = NaN;
if isstruct(c) && isfield(c, name) && ~isempty(c.(name)), x = double(c.(name)(1)); end
end

function S = demoRecording()
% Synthetic FM-band sweeps: carriers + a broadband "inverter" that starts at sweep 121.
f = linspace(88e6, 108e6, 401); n = 240; rbw = 100e3; sg = rbw/2.3548;
rs = RandStream('mt19937ar', 'Seed', 1);
cw = @(fc, p) 10.^(p/10) .* exp(-(f - fc(:)).^2/(2*sg^2));
band = @(fc, p, bw) 10^(p/10)/bw*1.065*rbw*0.5*(erf((fc + bw/2 - f)/(sqrt(2)*sg)) - erf((fc - bw/2 - f)/(sqrt(2)*sg)));
sig = cw(98.5e6, -35) + cw(94e6, -55) + band(103e6, -50, 2e6);
inv = band(91e6, -62, 6e6) + sum(cw(89e6 + (0:7)*0.5e6, -72), 1);
floorMw = 10^((-165 + 10*log10(1.065*rbw))/10);
P = zeros(n, numel(f));
for k = 1:n
    P(k, :) = 10*log10(floorMw*10.^(randn(rs, 1, numel(f))/10) + sig + (k > 120)*inv);
end
S.f = f; S.P = P;
S.t = datetime('now') + seconds((0:n - 1).'*0.1);
S.settings = struct('start', f(1), 'stop', f(end), 'points', numel(f), 'rbw', rbw, 'vbw', rbw, ...
    'ref', -20, 'pdiv', 10, 'binStep', f(2) - f(1));
S.idn = "Synthetic,DEMO,0,0";
end

function c = cellText(x)
% char for a table cell; missing -> ''
x = string(x);
if isempty(x) || ismissing(x(1)), c = ''; else, c = char(x(1)); end
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
