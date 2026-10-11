%% analyze_fieldfox.m
% Solar inverter EMI check: compares spectrum sweeps recorded with the
% inverter OFF and with it ON, using files saved by a Keysight FieldFox N9914A.
%
%   RecordPlayback\GPSA\AutoSession01
%                          SA "Record/Playback" session (binary files)
%                          - dSTA00000\*.STA  sweep settings
%                          - dTRC00000\*.TRC  recorded spectrum sweeps
%   CALIBRATE50OHM.sta     Saved instrument state. It is really a ZIP file
%                          holding XML/CSV files: a Cable & Antenna Test (CAT)
%                          return-loss trace and the spectrum analyzer settings.
%
% Outputs, written to the "results" folder next to this script:
%   fig1_off_vs_on.png        average spectrum OFF vs ON, and the difference
%   fig2_off_vs_on_zoom.png   same, zoomed to low frequencies
%   fig3_waterfall.png        every recorded sweep as a color map
%   fig4_return_loss.png      return loss from CALIBRATE50OHM.sta
%   band_noise_floor.csv      noise floor per frequency band, OFF vs ON
%   detected_signals.csv      signals above the noise floor, OFF vs ON
%   recorded_traces.csv       all sweeps as numbers (open in Excel)

clear; clc; close all;

%% ---- Settings --------------------------------------------------------------
dataDir = fullfile(getenv('USERPROFILE'), 'Downloads', ...
    'drive-download-20261008T090937Z-1-001');
stateFile  = fullfile(dataDir, 'CALIBRATE50OHM.sta');
sessionDir = fullfile(dataDir, 'RecordPlayback', 'GPSA', 'AutoSession01');

% Record numbers taken with the inverter OFF and ON, e.g. offRecords = 0:2;
% Leave both empty to split automatically at the largest time gap.
offRecords = [];
onRecords  = [];
excludeRecords = [];    % e.g. 4 to leave out an odd sweep

bandEdges_Hz  = [10e6 30e6 300e6 1e9 3e9 5e9];  % bands for the noise-floor table
zoomMax_Hz    = 500e6;  % upper limit of the zoomed plot
aboveFloor_dB = 10;     % a signal must be this far above the noise floor
dcCutoff_Hz   = 10e6;   % ignore the analyzer's own LO feedthrough near 0 Hz

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
outDir = fullfile(here, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---- 1. Load the recorded sweeps -------------------------------------------
st  = readFieldFoxState(stateFile);
hdr = readRecordHeader(fullfile(sessionDir, 'dSTA00000', '00000000.STA'));
[P, tRec, recNo] = readRecordTraces(fullfile(sessionDir, 'dTRC00000'), hdr.points);
nRec = size(P, 1);

f = linspace(hdr.center_Hz - hdr.span_Hz/2, hdr.center_Hz + hdr.span_Hz/2, hdr.points);
binStep = f(2) - f(1);
tSec  = seconds(tRec - tRec(1));
valid = f > dcCutoff_Hz;

fprintf('Instrument : %s, firmware %s\n', st.model, st.firmware);
fprintf('Span       : %.2f MHz to %.3f GHz, %d points (%.2f MHz per point)\n', ...
    f(1)/1e6, f(end)/1e9, hdr.points, binStep/1e6);
fprintf('Settings   : RBW %.0f kHz, VBW %.0f kHz, atten %g dB, preamp %s\n\n', ...
    st.sa.rbw_Hz/1e3, st.sa.vbw_Hz/1e3, st.sa.atten_dB, st.sa.preamp);

%% ---- 2. Split sweeps into inverter OFF / ON --------------------------------
if isempty(offRecords) && isempty(onRecords)
    [gap, iGap] = max(diff(tSec));
    isOff = (1:nRec).' <= iGap;
    isOn  = ~isOff;
    fprintf('OFF/ON split at the largest time gap: %.0f s between record %d and %d.\n', ...
        gap, recNo(iGap), recNo(iGap + 1));
else
    isOff = ismember(recNo, offRecords);
    isOn  = ismember(recNo, onRecords);
end
isOff = isOff & ~ismember(recNo, excludeRecords);
isOn  = isOn  & ~ismember(recNo, excludeRecords);
if ~any(isOff) || ~any(isOn)
    error('Need at least one OFF sweep and one ON sweep.');
end
%% the things must be in the value that can be fix form now accordking to the desk 

group = repmat("OFF", nRec, 1); group(isOn) = "ON";
group(~isOff & ~isOn) = "unused";
floorEach = median(P(:, valid), 2);
sweeps = table(recNo, string(datestr(tRec, 'HH:MM:SS')), tSec, group, floorEach, ...
    'VariableNames', {'Record', 'Time', 't_s', 'Inverter', 'NoiseFloor_dBm'});
disp(sweeps);

% Average in linear power (mW), not in dB
offAvg  = 10*log10(mean(10.^(P(isOff, :)/10), 1));
onAvg   = 10*log10(mean(10.^(P(isOn,  :)/10), 1));
delta   = onAvg - offAvg;

% A change is only convincing if it shows in EVERY sweep:
% "always higher" = the quietest ON sweep is above the loudest OFF sweep.
upMargin   = min(P(isOn, :), [], 1) - max(P(isOff, :), [], 1);
downMargin = min(P(isOff, :), [], 1) - max(P(isOn, :), [], 1);
alwaysUp   = valid & upMargin > 0;
alwaysDown = valid & downMargin > 0;

%% ---- 3. Noise floor per band -----------------------------------------------
% The median of the bins in a band ignores the few bins that hold strong
% carriers, so it tracks the noise floor that inverter noise would raise.
nb = numel(bandEdges_Hz) - 1;
bandName = strings(nb, 1); nBins = zeros(nb, 1);
offFloor = zeros(nb, 1); onFloor = zeros(nb, 1); nUp = zeros(nb, 1); nDown = zeros(nb, 1);
for b = 1:nb
    in = valid & f >= bandEdges_Hz(b) & f < bandEdges_Hz(b + 1);
    bandName(b) = sprintf('%s - %s', fmtHz(bandEdges_Hz(b)), fmtHz(bandEdges_Hz(b + 1)));
    nBins(b)    = nnz(in);
    offFloor(b) = median(offAvg(in));
    onFloor(b)  = median(onAvg(in));
    nUp(b)      = nnz(alwaysUp & in);
    nDown(b)    = nnz(alwaysDown & in);
end
bands = table(bandName, nBins, offFloor, onFloor, onFloor - offFloor, nUp, nDown, ...
    'VariableNames', {'Band', 'Points', 'OFF_dBm', 'ON_dBm', 'Change_dB', ...
                      'PointsAlwaysHigher', 'PointsAlwaysLower'});
fprintf('\n--- Noise floor per band (median of average trace) ---\n');
disp(bands);
fprintf(['PointsAlwaysHigher/Lower: points where every ON sweep is above (below)\n' ...
         'every OFF sweep. If "higher" and "lower" counts are similar, the change\n' ...
         'is normal sweep-to-sweep variation, not the inverter.\n\n']);

%% ---- 4. Signals: OFF vs ON -------------------------------------------------
maxHold  = max(P(isOff | isOn, :), [], 1);
thresh   = median(maxHold(valid)) + aboveFloor_dB;
isPk = islocalmax(maxHold, 'MinProminence', aboveFloor_dB/2) & valid & maxHold > thresh;
pk = find(isPk);

sig = table(f(pk).'/1e6, offAvg(pk).', onAvg(pk).', delta(pk).', ...
    alwaysUp(pk).' | alwaysDown(pk).', strings(numel(pk), 1), ...
    'VariableNames', {'Freq_MHz', 'OFF_dBm', 'ON_dBm', 'Change_dB', 'InEverySweep', 'PossibleSource'});
for k = 1:numel(pk)
    sig.PossibleSource(k) = guessBand(f(pk(k)), binStep/2);
end
fprintf('--- Signals above %.1f dBm (frequency accuracy +/- %.1f MHz) ---\n', thresh, binStep/2e6);
disp(sig);

newUp = find(alwaysUp & delta > 6);
fprintf('Points where ON is always higher AND the average rose more than 6 dB: %d\n', numel(newUp));
for k = newUp
    fprintf('  %8.1f MHz : OFF %6.1f dBm -> ON %6.1f dBm (+%.1f dB)\n', ...
        f(k)/1e6, offAvg(k), onAvg(k), delta(k));
end
fprintf('\n');

%% ---- 5. Figures ------------------------------------------------------------
plotOffOn(f, offAvg, onAvg, delta, alwaysUp, [0 f(end)], ...
    sprintf('Inverter OFF (%d sweeps) vs ON (%d sweeps), RBW %.0f MHz', ...
    nnz(isOff), nnz(isOn), st.sa.rbw_Hz/1e6));
exportgraphics(gcf, fullfile(outDir, 'fig1_off_vs_on.png'), 'Resolution', 150);

plotOffOn(f, offAvg, onAvg, delta, alwaysUp, [0 zoomMax_Hz], ...
    sprintf('Zoom 0 - %s', fmtHz(zoomMax_Hz)));
exportgraphics(gcf, fullfile(outDir, 'fig2_off_vs_on_zoom.png'), 'Resolution', 150);

% Waterfall: one row per sweep, line between the OFF and ON groups
figure('Name', 'Waterfall', 'Position', [100 100 1100 450]);
Pshow = P; Pshow(:, ~valid) = NaN;
imagesc(f/1e9, 1:nRec, Pshow, 'AlphaData', ~isnan(Pshow));
axis xy; colormap(jet); cb = colorbar; cb.Label.String = 'Power (dBm)';
yticks(1:nRec);
yticklabels(compose('#%d %s (%s)', recNo, string(datestr(tRec, 'HH:MM:SS')), group));
if ~isempty(find(isOn, 1))
    yline(find(isOn, 1) - 0.5, 'w-', 'LineWidth', 2);
end
xlabel('Frequency (GHz)'); ylabel('Sweep');
title('Waterfall of recorded sweeps');
exportgraphics(gcf, fullfile(outDir, 'fig3_waterfall.png'), 'Resolution', 150);

%% ---- 6. Return loss from CALIBRATE50OHM.sta --------------------------------
% RL = -20*log10(|Gamma|). Matched 50-ohm load: RL > 20 dB.
% Open or short (or nothing connected): RL ~ 0 dB.
rl = st.cat.rl_dB;
fprintf('--- Return loss in CALIBRATE50OHM.sta (calibrated %s) ---\n', st.cat.calTime);
fprintf('RL min/mean/max : %.2f / %.2f / %.2f dB, %.3f MHz to %.3f GHz\n', ...
    min(rl), mean(rl), max(rl), st.cat.freq_Hz(1)/1e6, st.cat.freq_Hz(end)/1e9);
if median(rl) < 3
    fprintf('Result: RL ~ 0 dB, the port saw an open/short (or nothing), not a 50-ohm load.\n\n');
elseif median(rl) > 20
    fprintf('Result: well matched to 50 ohm.\n\n');
else
    fprintf('Result: partly matched.\n\n');
end
figure('Name', 'Return loss', 'Position', [100 100 800 350]);
plot(st.cat.freq_Hz/1e9, rl, 'LineWidth', 1.5); grid on;
xlabel('Frequency (GHz)'); ylabel('Return loss (dB)');
title('CAT return loss, CALIBRATE50OHM.sta');
exportgraphics(gcf, fullfile(outDir, 'fig4_return_loss.png'), 'Resolution', 150);

%% ---- 7. Save numbers for Excel / your report -------------------------------
T = array2table(P.', 'VariableNames', compose('Rec%d_%s_dBm', recNo, group));
T = [table(f.', 'VariableNames', {'Freq_Hz'}), T, ...
     table(offAvg.', onAvg.', delta.', alwaysUp.', ...
     'VariableNames', {'OFF_avg_dBm', 'ON_avg_dBm', 'Change_dB', 'ON_always_higher'})];
writetable(T, fullfile(outDir, 'recorded_traces.csv'));
writetable(bands, fullfile(outDir, 'band_noise_floor.csv'));
writetable(sig, fullfile(outDir, 'detected_signals.csv'));
fprintf('Saved figures and CSV files to:\n  %s\n', outDir);


%% ===== Local functions ========================================================

function plotOffOn(f, offAvg, onAvg, delta, alwaysUp, xl, ttl)
figure('Position', [100 100 1100 650]);
ax1 = subplot(3, 1, [1 2]); hold on; grid on;
plot(f/1e9, offAvg, 'b', 'LineWidth', 1.2);
plot(f/1e9, onAvg, 'r', 'LineWidth', 1.2);
plot(f(alwaysUp)/1e9, onAvg(alwaysUp), 'ko', 'MarkerSize', 4);
ylabel('Power (dBm)'); title(ttl);
legend('Inverter OFF (average)', 'Inverter ON (average)', 'ON higher in every sweep', ...
    'Location', 'northeast');
ax2 = subplot(3, 1, 3); hold on; grid on;
stem(f/1e9, delta, 'Color', [0.6 0.6 0.6], 'Marker', 'none');
stem(f(alwaysUp)/1e9, delta(alwaysUp), 'r', 'Marker', 'none', 'LineWidth', 1.5);
yline(0, 'k');
xlabel('Frequency (GHz)'); ylabel('ON - OFF (dB)');
linkaxes([ax1 ax2], 'x'); xlim(ax1, xl/1e9);
end

function s = fmtHz(x)
if x >= 1e9, s = sprintf('%g GHz', x/1e9); else, s = sprintf('%g MHz', x/1e6); end
end

function st = readFieldFoxState(staFile)
% A FieldFox .sta file is a ZIP archive of XML and CSV files.
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));
unzip(staFile, tmp);

main = fileread(findFile(tmp, 'MainForm'));
st.model    = xmlValue(main, 'm_ModelNumber');
st.firmware = xmlValue(main, 'm_FirmwareVersion');
st.savedOn  = xmlValue(main, 'm_DateTime');

% Spectrum analyzer settings
sa = fileread(findFile(tmp, 'SaState.xml'));
st.sa.start_Hz = str2double(xmlValue(sa, 'StartFrequency'));
st.sa.stop_Hz  = str2double(xmlValue(sa, 'StopFrequency'));
st.sa.rbw_Hz   = str2double(xmlValue(sa, 'ResolutionBandwidth'));
st.sa.vbw_Hz   = str2double(xmlValue(sa, 'VideoBandwidth'));
st.sa.points   = str2double(xmlValue(sa, 'NumberOfPoints'));
st.sa.atten_dB = str2double(xmlValue(sa, 'RfAttenuation'));
st.sa.preamp   = xmlValue(sa, 'PreampOnState');
tr = fileread(findFile(tmp, 'SaTrace1.xml'));
st.sa.trace_dBm = sscanf(xmlValue(tr, 'Data'), '%f,').';

% Cable & Antenna Test trace, stored as CSV between BEGIN and END
csv  = fileread(findFile(tmp, 'DAT_File.dat.tmp'));
body = regexp(csv, 'BEGIN\s*(.*?)\s*END', 'tokens', 'once');
v    = reshape(sscanf(body{1}, '%f,%f'), 2, []).';
st.cat.freq_Hz = v(:, 1);
st.cat.rl_dB   = v(:, 2);
cal = regexp(csv, '! CORR TIMESTAMP1 ([^\r\n]*)', 'tokens', 'once');
if isempty(cal), st.cat.calTime = 'none'; else, st.cat.calTime = strtrim(cal{1}); end
end

function h = readRecordHeader(staFile)
% Byte layout of the recording's .STA file (little-endian):
%   offset  8  uint16  number of points
%   offset 16  uint64  timestamp, .NET ticks (100 ns since 1-Jan-0001)
%   offset 24  double  center frequency (Hz)
%   offset 32  double  span (Hz)
b = readBytes(staFile);
h.points    = double(typecast(b(9:10), 'uint16'));
h.time      = datetime(typecast(b(17:24), 'uint64'), 'ConvertFrom', '.net');
h.center_Hz = typecast(b(25:32), 'double');
h.span_Hz   = typecast(b(33:40), 'double');
end

function [P, t, recNo] = readRecordTraces(trcDir, nPts)
% Each record = 44-byte header + nPts float32 values (dBm) + 20 zero bytes.
%   offset 4   uint32  record number
%   offset 8   uint64  timestamp, .NET ticks
%   offset 44  single  trace values (dBm)
% One .TRC file can hold several records back to back.
recLen = 44 + 4*nPts + 20;
files  = dir(fullfile(trcDir, '*.TRC'));
P = zeros(0, nPts); t = NaT(0, 1); recNo = zeros(0, 1);
for k = 1:numel(files)
    b = readBytes(fullfile(files(k).folder, files(k).name));
    if mod(numel(b), recLen) ~= 0
        warning('%s: size %d is not a multiple of %d bytes.', files(k).name, numel(b), recLen);
    end
    for r = 0:floor(numel(b)/recLen) - 1
        rec = b(r*recLen + (1:recLen));
        recNo(end+1, 1) = double(typecast(rec(5:8), 'uint32')); %#ok<AGROW>
        t(end+1, 1)     = datetime(typecast(rec(9:16), 'uint64'), 'ConvertFrom', '.net'); %#ok<AGROW>
        P(end+1, :)     = double(typecast(rec(45:44 + 4*nPts), 'single')); %#ok<AGROW>
    end
end
[recNo, order] = sort(recNo);
t = t(order); P = P(order, :);
end

function name = guessBand(fHz, tol)
% Common services in Thailand. A match only says the frequency fits a band.
bands = {
    'FM radio',                 87.5e6,  108e6
    'Digital TV (UHF)',          470e6,  694e6
    'LTE/5G 700 (B28)',          703e6,  803e6
    'GSM/LTE 850 (B5)',          824e6,  894e6
    'GSM/LTE 900 (B8)',          880e6,  960e6
    'LTE 1800 (B3)',            1710e6, 1880e6
    '3G/LTE 2100 (B1)',         1920e6, 2170e6
    'LTE TDD 2300 (B40)',       2300e6, 2400e6
    'Wi-Fi/Bluetooth 2.4 GHz',  2400e6, 2483.5e6
    '5G NR 2600 (n41)',         2496e6, 2690e6
    '5G NR 3500 (n78)',         3300e6, 3800e6
    'Wi-Fi 5 GHz',              5150e6, 5850e6 };
hit = fHz >= [bands{:, 2}] - tol & fHz <= [bands{:, 3}] + tol;
if any(hit)
    name = string(strjoin(bands(hit, 1), ' / '));
else
    name = "unknown";
end
end

function p = findFile(root, name)
d = dir(fullfile(root, '**', name));
if isempty(d), error('"%s" was not found inside the state file.', name); end
p = fullfile(d(1).folder, d(1).name);
end

function v = xmlValue(txt, tag)
tok = regexp(txt, ['<' tag '(?:\s[^>]*)?>(.*?)</' tag '>'], 'tokens', 'once');
if isempty(tok), v = ''; else, v = strtrim(tok{1}); end
end

function b = readBytes(file)
fid = fopen(file, 'r');
if fid < 0, error('Cannot open %s', file); end
b = fread(fid, Inf, '*uint8');
fclose(fid);
end
