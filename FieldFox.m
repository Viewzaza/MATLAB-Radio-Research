classdef FieldFox < handle
%FIELDFOX  SCPI connection to a Keysight FieldFox over LAN (socket 5025).
%
%   ff = FieldFox("10.90.36.200");        connect (base MATLAB tcpclient, no VISA)
%   ff.takeControl();                      back up the instrument state, single-sweep mode
%   ff.saConfigure(Center=100e6, Span=20e6, RBW=100e3, Preamp=true);
%   ff.sweep();                            one synchronized sweep
%   [y, f] = ff.saTrace();                 dBm trace + frequency axis (binary transfer)
%   ff.release();                          restore the instrument exactly as it was
%   clear ff                               (release also runs automatically)
%
% Mode helpers:  ff.setMode("NA"), ff.naFrequencies(), ff.naSData(trace)
% Files:         ff.screenshot(), ff.getFile(name), ff.saveState(name), ff.loadState(name)
% Low level:     ff.write(cmd), ff.query(cmd), ff.queryNum(cmd), ff.queryBlock(cmd, type)
%
% Verified on N9914A firmware A.12.45 (options 112,210,211,233,235,236,308,310).

    properties (SetAccess = private)
        IP        string
        Port      double
        IDN       string
        Model     string
        Serial    string
        Firmware  string
        Options   string        % e.g. ["112" "210" ...]
        Client                  % tcpclient
        InControl logical = false
    end

    properties
        Verbose logical = false % print every command
        DefaultTimeout double = 10
    end

    properties (Access = private)
        Format   string = ""    % current FORM:DATA as last set ("ASC" / "REAL,32" / "REAL,64")
        Saved    struct = struct()
        BackupFile string = ""
    end

    methods
        %% ---- connection ----------------------------------------------------
        function obj = FieldFox(ip, port)
            arguments
                ip (1,1) string = "10.90.36.200"
                port (1,1) double = 5025
            end
            obj.IP = ip; obj.Port = port;
            obj.Client = tcpclient(ip, port, 'Timeout', obj.DefaultTimeout, 'ConnectTimeout', 5);
            configureTerminator(obj.Client, "LF");
            obj.IDN = obj.query("*IDN?");
            p = split(obj.IDN, ",");
            if numel(p) >= 4
                obj.Model = strtrim(p(2)); obj.Serial = strtrim(p(3)); obj.Firmware = strtrim(p(4));
            end
            opt = erase(obj.query("*OPT?"), '"');
            obj.Options = strtrim(split(opt, ",")).';
            obj.Format = obj.formatNow();
        end

        function delete(obj)
            try
                if obj.InControl, obj.release(); end
            catch err
                warning('FieldFox:release', 'Could not restore the instrument: %s', err.message);
            end
            obj.Client = [];
        end

        %% ---- basic I/O -----------------------------------------------------
        function write(obj, cmd)
            if obj.Verbose, fprintf('>> %s\n', cmd); end
            writeline(obj.Client, cmd);
        end

        function r = query(obj, cmd, timeout)
            if nargin > 2, c = obj.tempTimeout(timeout); end %#ok<NASGU>
            obj.write(cmd);
            w = warning('off', 'transportlib:client:ReadlineWarning');
            r = readline(obj.Client);
            warning(w);
            if isempty(r)
                flush(obj.Client);       % drop a late reply so it cannot answer the next query
                error('FieldFox:timeout', 'No reply to "%s" within %g s.', cmd, obj.Client.Timeout);
            end
            r = strtrim(r);
            if obj.Verbose, fprintf('<< %s\n', r); end
        end

        function v = queryNum(obj, cmd, varargin)
            % Numbers from an ASCII reply (comma-separated lists allowed).
            r = obj.query(cmd, varargin{:});
            v = sscanf(char(erase(r, '"')), '%f,').';
        end

        function v = queryBlock(obj, cmd, precision, timeout)
            % Binary block reply ('#<n><len><bytes>', little endian) as double.
            % precision: "single" (FORM REAL,32) or "double" (FORM REAL,64).
            arguments
                obj; cmd string; precision string = "single"; timeout double = NaN
            end
            if precision == "single", obj.setFormat("REAL,32"); else, obj.setFormat("REAL,64"); end
            if ~isnan(timeout), c = obj.tempTimeout(timeout); end %#ok<NASGU>
            obj.write(cmd);
            w = warning('off', 'transportlib:client:ReadbinblockWarning');
            v = double(readbinblock(obj.Client, precision)).';
            warning(w);
            if isempty(v)
                flush(obj.Client);
                error('FieldFox:timeout', 'No binary reply to "%s" within %g s.', cmd, obj.Client.Timeout);
            end
            read(obj.Client, 1, "uint8");               % the trailing LF
        end

        function opc(obj, timeout)
            % Wait for all pending operations (*OPC?).
            if nargin < 2, timeout = 60; end
            obj.query("*OPC?", timeout);
        end

        function errs = errors(obj)
            % Read and clear the instrument error queue. Empty = no errors.
            errs = strings(0, 1);
            for k = 1:50
                e = obj.query("SYST:ERR?");
                if startsWith(e, "0,") || startsWith(e, "+0,"), break, end
                errs(end + 1, 1) = e; %#ok<AGROW>
            end
        end

        function syncFormat(obj)
            % Re-read FORM:DATA from the instrument (call after MMEM:LOAD:STAT
            % or when another program may have changed it).
            obj.Format = obj.formatNow();
        end

        function ok = hasOption(obj, opt)
            ok = any(obj.Options == string(opt));
        end
%% fix the things that is going onthere 
        %% ---- modes and control ---------------------------------------------
        function m = mode(obj)
            m = erase(obj.query("INST:SEL?"), '"');
        end

        function setMode(obj, m)
            % "SA", "NA", "CAT", "VVM", "CPM", ... (case sensitive on the FieldFox)
            if obj.mode() == string(m), return, end
            obj.query("INST:SEL """ + m + """;*OPC?", 60);
            if obj.InControl, obj.query("INIT:CONT 0;*OPC?", 30); end
        end

        function takeControl(obj, backup)
            % Remember the instrument state and switch to single-sweep mode.
            % backup=true also saves ALL mode settings to a .sta file on the
            % instrument so release() can put everything back exactly. The file
            % name holds the time and this MATLAB's process ID, so a backup left
            % by a MATLAB that was killed is found and restored next time.
            if nargin < 2, backup = true; end
            if obj.InControl, return, end
            obj.Saved.cdir = erase(obj.query("MMEM:CDIR?"), '"');
            obj.write("MMEM:CDIR ""[INTERNAL]:""");          % backups always live here
            obj.recoverCrashedSession();
            obj.Saved.mode = obj.mode();
            obj.Saved.cont = obj.query("INIT:CONT?");
            obj.Saved.format = obj.formatNow();
            if backup
                obj.BackupFile = sprintf("matlab_backup_%s_%d.sta", ...
                    string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')), feature('getpid'));
                obj.query("MMEM:STOR:STAT """ + obj.BackupFile + """;*OPC?", 60);
                e = obj.errors();
                if ~isempty(e)
                    warning('FieldFox:backup', 'State backup failed (%s): release() can only restore mode, format and sweep mode.', strjoin(e, '; '));
                    obj.BackupFile = "";
                end
            end
            obj.write("MMEM:CDIR """ + obj.Saved.cdir + """");
            obj.query("INIT:CONT 0;*OPC?", 30);
            obj.InControl = true;
            obj.errors();     % clear the queue
        end

        function release(obj)
            % Restore what takeControl() saved. Every step runs even if an
            % earlier one fails; the first error is reported as a warning.
            if ~obj.InControl, return, end
            obj.InControl = false;
            msg = strings(0, 1);
            step = @(f) tryStep(f);
            if strlength(obj.BackupFile) > 0
                step(@() obj.write("MMEM:CDIR ""[INTERNAL]:"""));
                step(@() obj.query("MMEM:LOAD:STAT """ + obj.BackupFile + """;*OPC?", 120));
                step(@() obj.write("MMEM:DEL """ + obj.BackupFile + """"));
                step(@() obj.write("MMEM:CDIR """ + obj.Saved.cdir + """"));
                obj.BackupFile = "";
            end
            step(@() obj.syncFormat());
            step(@() restoreMode());
            step(@() obj.setFormat(obj.Saved.format));
            step(@() obj.write("INIT:CONT " + obj.Saved.cont));
            if ~isempty(msg)
                warning('FieldFox:release', 'Restore was incomplete: %s', strjoin(msg, '; '));
            end

            function tryStep(f)
                try
                    f();
                catch err
                    msg(end + 1, 1) = err.message;
                end
            end

            function restoreMode()
                if obj.mode() ~= obj.Saved.mode
                    obj.query("INST:SEL """ + obj.Saved.mode + """;*OPC?", 60);
                end
            end
        end

        function recoverCrashedSession(obj)
            % A MATLAB that was killed while in control leaves its backup .sta
            % behind (and the FieldFox in single sweep). Restore the newest such
            % backup whose MATLAB process is no longer running on this PC.
            r = obj.query("MMEM:CAT? ""[INTERNAL]:""", 30);
            tok = regexp(char(r), '"(matlab_backup[^",]*\.sta)', 'tokens');
            if isempty(tok), return, end
            names = string(cellfun(@(c) c{1}, tok, 'UniformOutput', false));
            pid = str2double(regexprep(names, '^matlab_backup_\d{8}_\d{6}_(\d+)\.sta$', '$1'));
            live = false(size(names));
            for k = find(isfinite(pid) & pid ~= feature('getpid'))   % this MATLAB is not in control yet
                [~, out] = system(sprintf('tasklist /FI "PID eq %d" /NH', pid(k)));
                live(k) = contains(lower(out), "matlab");
            end
            stale = sort(names(~live));       % sorts by time; the old fixed name sorts first
            if isempty(stale), return, end
            newest = stale(end);
            warning('FieldFox:recover', ['A previous MATLAB session did not restore the FieldFox. ' ...
                'Restoring the state it saved (%s) first.'], newest);
            obj.query("MMEM:LOAD:STAT """ + newest + """;*OPC?", 120);
            for n = stale(:).'
                obj.write("MMEM:DEL """ + n + """");
            end
            obj.write("FORM:DATA ASC,0");
            obj.write("INIT:CONT 1");
            obj.syncFormat();
            obj.errors();
        end

        function t = sweep(obj, timeout)
            % One synchronized sweep (INIT:IMM;*OPC?). Returns the sweep time (s).
            % In free-run mode (not in control) it just waits for the current sweep.
            if nargin < 2, timeout = 120; end
            tic;
            if obj.InControl
                obj.query("INIT:IMM;*OPC?", timeout);
            else
                obj.query("*OPC?", timeout);
            end
            t = toc;
        end

        function restartAverage(obj, timeout)
            % Complete set of averages in single mode (INIT:REST;*OPC?).
            if nargin < 2, timeout = 300; end
            obj.query("INIT:REST;*OPC?", timeout);
        end

        %% ---- spectrum analyzer -----------------------------------------------
        function s = saSettings(obj)
            s.start  = obj.queryNum("SENS:FREQ:STAR?");
            s.stop   = obj.queryNum("SENS:FREQ:STOP?");
            s.points = obj.queryNum("SENS:SWE:POIN?");
            s.rbw    = obj.queryNum("SENS:BAND:RES?");
            s.vbw    = obj.queryNum("SENS:BAND:VID?");
            s.ref    = obj.queryNum("DISP:WIND:TRAC:Y:RLEV?");
            s.pdiv   = obj.queryNum("DISP:WIND:TRAC:Y:PDIV?");
            s.atten  = obj.queryNum("SENS:POW:RF:ATT?");
            s.preamp = obj.queryNum("SENS:POW:RF:GAIN?") == 1;
            s.detector = erase(obj.query("SENS:DET:FUNC?"), '"');
            s.f      = linspace(s.start, s.stop, s.points);
        end

        function saConfigure(obj, opts)
            % Set SA parameters. Only the ones given are changed.
            arguments
                obj
                opts.Center double = []
                opts.Span double = []
                opts.Start double = []
                opts.Stop double = []
                opts.RBW double = []          % Hz, or "auto" via RBWAuto=true
                opts.RBWAuto logical = logical([])
                opts.VBW double = []
                opts.Points double = []
                opts.RefLevel double = []
                opts.Atten double = []        % dB
                opts.AttenAuto logical = logical([])
                opts.Preamp logical = logical([])
                opts.Detector string = string([])   % AUTO|NORM|POS|NEG|SAMP|AVER
                opts.Average double = []       % sweeps (1 = off)
                opts.TraceType string = string([]) % CLRW|AVG|MAXH|MINH
            end
            if obj.mode() ~= "SA", obj.setMode("SA"); end
            set = @(cmd, v) obj.write(sprintf(cmd, v));
            if ~isempty(opts.Start),  set("SENS:FREQ:STAR %.10g", opts.Start); end
            if ~isempty(opts.Stop),   set("SENS:FREQ:STOP %.10g", opts.Stop); end
            if ~isempty(opts.Center), set("SENS:FREQ:CENT %.10g", opts.Center); end
            if ~isempty(opts.Span),   set("SENS:FREQ:SPAN %.10g", opts.Span); end
            if ~isempty(opts.Points), set("SENS:SWE:POIN %d", opts.Points); end
            if ~isempty(opts.RBWAuto), set("SENS:BAND:RES:AUTO %d", opts.RBWAuto); end
            if ~isempty(opts.RBW),    set("SENS:BAND:RES %.10g", opts.RBW); end
            if ~isempty(opts.VBW),    set("SENS:BAND:VID %.10g", opts.VBW); end
            if ~isempty(opts.AttenAuto), set("SENS:POW:RF:ATT:AUTO %d", opts.AttenAuto); end
            if ~isempty(opts.Atten),  set("SENS:POW:RF:ATT %g", opts.Atten); end
            if ~isempty(opts.Preamp), set("SENS:POW:RF:GAIN %d", opts.Preamp); end
            if ~isempty(opts.RefLevel), set("DISP:WIND:TRAC:Y:RLEV %g", opts.RefLevel); end
            if ~isempty(opts.Detector), obj.write("SENS:DET:FUNC " + opts.Detector); end
            if ~isempty(opts.TraceType), obj.write("TRAC1:TYPE " + opts.TraceType); end
            if ~isempty(opts.Average), set("SENS:AVER:COUN %d", opts.Average); end
            obj.opc(30);
            e = obj.errors();
            if ~isempty(e), warning('FieldFox:config', 'Instrument reported: %s', strjoin(e, '; ')); end
        end

        function [y, f] = saTrace(obj, trace)
            % SA trace in dBm (binary transfer) and its frequency axis in Hz.
            if nargin < 2, trace = 1; end
            y = obj.queryBlock(sprintf("TRAC%d:DATA?", trace), "single");
            if nargout > 1
                f = linspace(obj.queryNum("SENS:FREQ:STAR?"), obj.queryNum("SENS:FREQ:STOP?"), numel(y));
            end
        end

        %% ---- network analyzer ------------------------------------------------
        function f = naFrequencies(obj)
            f = obj.queryBlock("SENS:FREQ:DATA?", "double");
        end

        function s = naSData(obj, trace)
            % Complex, error-corrected data of NA trace <trace> (CALC:DATA:SDATA?).
            obj.write(sprintf("CALC:PAR%d:SEL", trace));
            v = obj.queryBlock("CALC:DATA:SDATA?", "double");
            s = complex(v(1:2:end), v(2:2:end));
        end

        %% ---- files and state ---------------------------------------------------
        function saveState(obj, name)
            obj.query("MMEM:STOR:STAT """ + name + """;*OPC?", 60);
        end

        function loadState(obj, name)
            obj.query("MMEM:LOAD:STAT """ + name + """;*OPC?", 120);
            obj.syncFormat();
            if obj.InControl, obj.query("INIT:CONT 0;*OPC?", 30); end
        end

        function data = getFile(obj, name)
            % Bytes of a file in the instrument's current directory (MMEM:DATA?).
            c = obj.tempTimeout(60); %#ok<NASGU>
            obj.write("MMEM:DATA? """ + name + """");
            data = readbinblock(obj.Client, "uint8");
            read(obj.Client, 1, "uint8");
            data = uint8(data(:));
        end

        function img = screenshot(obj, file)
            % Screen image of the FieldFox. Saves to <file> (.png) if given.
            tmp = "matlab_screen.png";
            obj.query("MMEM:STOR:IMAG """ + tmp + """;*OPC?", 30);
            bytes = obj.getFile(tmp);
            obj.write("MMEM:DEL """ + tmp + """");
            if nargin < 2 || strlength(string(file)) == 0
                file = [tempname '.png'];
                del = onCleanup(@() delete(file));
            end
            fid = fopen(file, 'w'); fwrite(fid, bytes, 'uint8'); fclose(fid);
            img = imread(file);
        end
    end

    methods (Access = private)
        function f = formatNow(obj)
            r = obj.query("FORM:DATA?");
            if startsWith(r, "ASC"), f = "ASC";
            elseif contains(r, "64"), f = "REAL,64";
            else, f = "REAL,32";
            end
        end

        function setFormat(obj, fmt)
            if obj.Format == fmt, return, end
            if fmt == "ASC", obj.write("FORM:DATA ASC,0"); else, obj.write("FORM:DATA " + fmt); end
            obj.Format = fmt;
        end

        function c = tempTimeout(obj, t)
            old = obj.Client.Timeout;
            obj.Client.Timeout = t;
            c = onCleanup(@() restoreTimeout(obj.Client, old));
        end
    end
end

function restoreTimeout(client, t)
if ~isempty(client) && isvalid(client), client.Timeout = t; end
end
