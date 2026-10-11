function TPL = ff_signal_templates(varargin)
%FF_SIGNAL_TEMPLATES  Reference power-spectrum templates of common radio signals.
%
%   TPL = ff_signal_templates                 load the cached templates (build if missing)
%   TPL = ff_signal_templates('Rebuild', true) synthesise them again (about 1.5 min)
%   ff_signal_templates('Plot', true)          show all template spectra
%
% A FieldFox gives power spectra only (no I/Q), so a signal can only be
% recognised by the SHAPE and WIDTH of its spectrum. This function makes
% the reference shapes: it generates real standard-compliant baseband
% waveforms with the MATLAB wireless toolboxes and turns each into a
% normalised power spectral density (pwelch):
%
%   5G NR DL      nrDLCarrierConfig + nrWaveformGenerator   5/10/20/40 MHz @ 15 and 30 kHz, 100 MHz @ 30 kHz
%   LTE DL        lteRMCDL + lteRMCDLTool (all RBs + OCNG)  1.4/3/5/10/15/20 MHz
%   Wi-Fi         wlanNonHTConfig/wlanHTConfig/wlanVHTConfig/wlanHESUConfig + wlanWaveformGenerator, 20/40/80 MHz
%   Bluetooth     bleWaveformGenerator (LE1M, LE2M), bluetoothWaveformGenerator (BR, EDR2M)
%   GSM           gsmDownlinkConfig + gsmFrame (all 8 timeslots busy)
%   FM broadcast  comm.FMBroadcastModulator (stereo + RDS, dsp.ColoredNoise "music")
%   GPS L1 C/A    gpsWaveformGenerator (Satellite Communications Toolbox)
%   DVB-T2 8 MHz  comm.OFDMModulator (8K FFT, 6817 carriers; generic OFDM, no pilots)
%   CW            a pure carrier (shows as the RBW filter shape)
%
% Each template also stores "views": the spectrum as a FieldFox would draw
% it through its Gaussian RBW filter, pre-computed for RBW = 10 Hz...6 MHz.
%
% Assumptions: fully loaded carriers (every resource block busy) and ideal
% transmitters. A lightly loaded LTE/NR cell or a bursty Wi-Fi link looks
% weaker and more ragged on the analyzer than its template.
%
% Name-Value options (defaults):
%   'Rebuild'  false   build again even if the cache exists
%   'File'     ff_signal_templates.mat next to this file
%   'Plot'     false   figure of all templates
%   'Visible'  'on'    figure visibility ('off' for batch use)
%   'Prompt'   true    accepted for suite compatibility (no dialogs are used)
%   'Verbose'  true
%
% Output struct TPL:
%   .Templates  struct array: Name, Class, Technology, Generator, NominalBW_Hz,
%               OBW_Hz (99 %), BW26_Hz (-26 dB), Fs_Hz, f (offset Hz), psd
%               (1/Hz, unit total power), psd_dB (0 dB peak), Views
%   .Classes    the class labels (e.g. "LTE 10 MHz") used by ff_signal_classifier
%   .RBWGrid    RBW values (Hz) of the views
%
% See also ff_signal_classifier, ff_signal_id.

opt = struct('Rebuild', false, 'File', "", 'Plot', false, 'Visible', 'on', 'Prompt', true, 'Verbose', true);
for k = 1:2:numel(varargin)
    if ~isfield(opt, varargin{k}), error('Unknown option "%s".', varargin{k}); end
    opt.(varargin{k}) = varargin{k + 1};
end
if strlength(string(opt.File)) == 0
    opt.File = fullfile(fileparts(mfilename('fullpath')), 'ff_signal_templates.mat');
end

persistent cache cacheFile
if ~opt.Rebuild && ~isempty(cache) && cacheFile == string(opt.File)
    TPL = cache;
elseif ~opt.Rebuild && isfile(opt.File)
    S = load(opt.File, 'TPL');
    TPL = S.TPL;
else
    TPL = buildAll(opt.Verbose);
    TPL.File = string(opt.File);
    save(opt.File, 'TPL', '-v7');
    if opt.Verbose, fprintf('Templates saved to %s\n', opt.File); end
end
cache = TPL; cacheFile = string(opt.File);
if opt.Plot, plotTemplates(TPL, opt.Visible); end
end

%% ===== building ============================================================

function TPL = buildAll(verbose)
rng(1);                                    % reproducible templates
tAll = tic;
list = {};
add = @(name, cls, tech, gen, bw, fn) {name, cls, tech, gen, bw, fn};

% 5G NR: max. resource blocks per TS 38.101-1 Table 5.3.2-1
nr = [5 15 25; 5 30 11; 10 15 52; 10 30 24; 20 15 106; 20 30 51; 40 15 216; 40 30 106; 100 30 273];
for r = 1:size(nr, 1)
    b = nr(r, 1); scs = nr(r, 2); nrb = nr(r, 3);
    list(end + 1, :) = add(sprintf('NR %g MHz %d kHz SCS', b, scs), sprintf('NR %g MHz', b), "NR", ...
        "nrWaveformGenerator", b*1e6, @() genNR(b, scs, nrb)); %#ok<AGROW>
end
lte = [1.4 6; 3 15; 5 25; 10 50; 15 75; 20 100];
for r = 1:size(lte, 1)
    list(end + 1, :) = add(sprintf('LTE %g MHz (%d RB)', lte(r, 1), lte(r, 2)), sprintf('LTE %g MHz', lte(r, 1)), "LTE", ...
        "lteRMCDLTool", lte(r, 1)*1e6, @() genLTE(lte(r, 2))); %#ok<AGROW>
end
wl = {'Wi-Fi 802.11a/g 20 MHz', 'nonht', 20; 'Wi-Fi 802.11n 20 MHz', 'ht', 20; 'Wi-Fi 802.11ax 20 MHz', 'he', 20
      'Wi-Fi 802.11n 40 MHz', 'ht', 40;     'Wi-Fi 802.11ax 40 MHz', 'he', 40
      'Wi-Fi 802.11ac 80 MHz', 'vht', 80;   'Wi-Fi 802.11ax 80 MHz', 'he', 80};
for r = 1:size(wl, 1)
    list(end + 1, :) = add(wl{r, 1}, sprintf('Wi-Fi %d MHz', wl{r, 3}), "Wi-Fi", "wlanWaveformGenerator", ...
        wl{r, 3}*1e6, @() genWLAN(wl{r, 2}, wl{r, 3})); %#ok<AGROW>
end
list(end + 1, :) = add('BLE 1M', 'BLE 1M', "BLE", "bleWaveformGenerator", 2e6, @() genBLE('LE1M'));
list(end + 1, :) = add('BLE 2M', 'BLE 2M', "BLE", "bleWaveformGenerator", 2e6, @() genBLE('LE2M'));
list(end + 1, :) = add('Bluetooth BR (GFSK)', 'Bluetooth BR/EDR', "Bluetooth", "bluetoothWaveformGenerator", 1e6, @() genBR('BR'));
list(end + 1, :) = add('Bluetooth EDR 2M', 'Bluetooth BR/EDR', "Bluetooth", "bluetoothWaveformGenerator", 1e6, @() genBR('EDR2M'));
list(end + 1, :) = add('GSM (GMSK)', 'GSM', "GSM", "gsmFrame", 200e3, @genGSM);
list(end + 1, :) = add('FM broadcast stereo+RDS', 'FM broadcast', "FM", "comm.FMBroadcastModulator", 200e3, @genFM);
list(end + 1, :) = add('GPS L1 C/A', 'GPS L1 C/A', "GNSS", "gpsWaveformGenerator", 2.046e6, @genGPS);
list(end + 1, :) = add('DVB-T2 8 MHz (8K)', 'DVB-T2 8 MHz', "DVB-T2", "comm.OFDMModulator", 8e6, @genDVBT2);
list(end + 1, :) = add('CW carrier', 'CW', "CW", "(analytic)", 0, []);

rbwGrid = 10.^(1:0.2:6.8);                % 10 Hz ... 6.3 MHz, 5 per decade
T = struct([]);
for k = 1:size(list, 1)
    t0 = tic;
    [name, cls, tech, gen, bw, fn] = list{k, :};
    if isempty(fn)                         % CW: a single line
        f = 0; psd = 1; fs = NaN; o99 = 0; b26 = 0;
    else
        [x, fs] = fn();
        x = x(movmax(abs(x), 65) > 1e-9 * max(abs(x)));   % drop idle gaps between bursts
        nfft = 2^nextpow2(fs / (max(bw, 200e3) / 400));
        nfft = min(max(nfft, 512), 16384);
        nfft = min(nfft, 2^floor(log2(numel(x) / 4)));
        [p, f] = pwelch(x, hann(nfft), nfft/2, nfft, fs, 'centered');
        f = f(:).'; df = f(2) - f(1);
        psd = p(:).' / (sum(p) * df);     % unit total power
        o99 = obw(psd, f);
        b26 = powerbw(psd, f, [], 26);
    end
    T(k).Name = string(name); T(k).Class = string(cls); T(k).Technology = tech;
    T(k).Generator = gen; T(k).NominalBW_Hz = bw; T(k).OBW_Hz = o99; T(k).BW26_Hz = b26;
    T(k).Fs_Hz = fs; T(k).f = single(f); T(k).psd = single(psd);
    T(k).psd_dB = single(10*log10(psd / max(psd) + 1e-30));
    T(k).Views = renderViews(f, psd, rbwGrid);
    if verbose
        fprintf('  %-28s OBW %8.3f MHz  -26dB BW %8.3f MHz  (%.1f s)\n', name, o99/1e6, b26/1e6, toc(t0));
    end
end
TPL.Templates = T;
TPL.Classes = unique([T.Class], 'stable');
TPL.RBWGrid = rbwGrid;
TPL.Created = datetime('now');
TPL.MATLAB = string(version);
if verbose, fprintf('%d templates, %d classes built in %.0f s\n', numel(T), numel(TPL.Classes), toc(tAll)); end
end

function V = renderViews(f, psd, rbwGrid)
% The trace a swept analyzer would show through a Gaussian RBW filter
% (3 dB width = rbw, peak gain 1), normalised to a 0 dB peak.
V.rbw = rbwGrid;
V.f = cell(size(rbwGrid)); V.s = cell(size(rbwGrid)); V.gain = zeros(size(rbwGrid));
for i = 1:numel(rbwGrid)
    sig = rbwGrid(i) / sqrt(8*log(2));                 % Gaussian std (power response)
    if isscalar(f)                                     % CW: the filter shape itself
        x = linspace(-5*sig, 5*sig, 201);
        y = exp(-x.^2 / (2*sig^2)); g = 1;
    else
        df = f(2) - f(1);
        pw = psd * df;                                 % power per frequency bin
        if sig < df/2                                  % RBW finer than the template grid
            x = f; y = psd * 1.0645 * rbwGrid(i);      % PSD x noise bandwidth
        else
            pad = ceil(5*sig / df);
            x = f(1) - pad*df : df : f(end) + pad*df + df/2;
            x = x(1:numel(f) + 2*pad);
            h = exp(-((-pad:pad)*df).^2 / (2*sig^2));
            y = conv([zeros(1, pad) pw zeros(1, pad)], h, 'same');
        end
        g = max(y);                                    % displayed peak / total power
        keep = find(y > g * 1e-8);                     % down to -80 dB
        x = x(keep(1):keep(end)); y = y(keep(1):keep(end));
        step = max(1, floor(numel(x) / 2048));
        x = x(1:step:end); y = y(1:step:end);
    end
    V.f{i} = single(x); V.s{i} = single(y / max(y)); V.gain(i) = g;
end
end

%% ===== waveform generators ===============================================

function [x, fs] = genNR(bwMHz, scs, nrb)
cfg = nrDLCarrierConfig('FrequencyRange', 'FR1', 'ChannelBandwidth', bwMHz, 'NumSubframes', 2);
cfg.SCSCarriers{1}.SubcarrierSpacing = scs; cfg.SCSCarriers{1}.NSizeGrid = nrb;
cfg.BandwidthParts{1}.SubcarrierSpacing = scs; cfg.BandwidthParts{1}.NSizeBWP = nrb;
cfg.SSBurst.Enable = false;            % keep the full-load PDSCH shape only
cfg.PDCCH{1}.Enable = false; cfg.CSIRS{1}.Enable = false;
cfg.PDSCH{1}.PRBSet = 0:nrb - 1;
cfg.PDSCH{1}.Modulation = '64QAM';
[x, info] = nrWaveformGenerator(cfg);
[x, fs] = txFilter(x, info.ResourceGrids(1).Info.SampleRate, nrb*12*scs*1e3, bwMHz*1e6/2);
end

function [x, fs] = genLTE(nrb)
% No 'RC' field: every subframe scheduled; OCNG fills unused resource elements.
rmc = struct('NDLRB', nrb, 'CellRefP', 1, 'TotSubframes', 4, ...
    'OCNGPDSCHEnable', 'On', 'OCNGPDCCHEnable', 'On');
rmc.PDSCH.PRBSet = (0:nrb - 1).';
rmc.PDSCH.Modulation = '64QAM';
[x, ~, info] = lteRMCDLTool(rmc, randi([0 1], 5000, 1));
bw = [1.4 3 5 10 15 20]; bw = bw([6 15 25 50 75 100] == nrb);
[x, fs] = txFilter(x, info.SamplingRate, nrb*180e3, bw*1e6/2);
end

function [x, fs] = genWLAN(kind, bwMHz)
cbw = sprintf('CBW%d', bwMHz);
switch kind
    case 'nonht', cfg = wlanNonHTConfig('ChannelBandwidth', cbw, 'MCS', 4, 'PSDULength', 1500);
    case 'ht',    cfg = wlanHTConfig('ChannelBandwidth', cbw, 'MCS', 4, 'PSDULength', 1500);
    case 'vht',   cfg = wlanVHTConfig('ChannelBandwidth', cbw, 'MCS', 4, 'APEPLength', 1500);
    case 'he',    cfg = wlanHESUConfig('ChannelBandwidth', cbw, 'MCS', 4, 'APEPLength', 1500);
end
if kind == "he", n = getPSDULength(cfg); else, n = cfg.PSDULength; end
osf = 2;                                   % room for the spectral skirts
x = wlanWaveformGenerator(randi([0 1], 8*n, 1), cfg, 'NumPackets', 4, 'IdleTime', 0, 'OversamplingFactor', osf);
fs = wlanSampleRate(cfg) * osf;
end

function [x, fs] = genBLE(mode)
sps = 8; rsym = 1e6 * (1 + strcmp(mode, 'LE2M'));
x = [];
for k = 1:12
    x = [x; bleWaveformGenerator(randi([0 1], 2000, 1), 'Mode', mode, 'SamplesPerSymbol', sps, ...
        'ChannelIndex', randi([0 36]))]; %#ok<AGROW>
end
fs = sps * rsym;
end

function [x, fs] = genBR(mode)
cfg = bluetoothWaveformConfig('Mode', mode, 'PacketType', 'DH5', 'SamplesPerSymbol', 8);
cfg.PayloadLength = 300;
x = bluetoothWaveformGenerator(randi([0 1], getPayloadLength(cfg)*8*6, 1), cfg);
fs = 8e6;
end

function [x, fs] = genGSM()
cfg = gsmDownlinkConfig;                   % 8 normal bursts per frame, 16 samples/symbol
x = gsmFrame(cfg, 12);
fs = gsmInfo(cfg).SampleRate;
end

function [x, fs] = genFM()
% Stereo programme ("music" = pink noise, 15 kHz audio band) + RDS at 57 kHz.
sps = 10; groups = 8;
fsA = 40 * 1187.5;                         % 47.5 kHz audio rate (RDS-friendly)
fs = 8 * 57e3;                             % 456 kHz output rate
rbds = comm.RBDSWaveformGenerator('GroupsPerFrame', groups, 'SamplesPerSymbol', sps);
nA = 104 * sps * groups * fsA / (1187.5 * sps);
noise = dsp.ColoredNoise('Color', 'pink', 'SamplesPerFrame', nA, 'NumChannels', 2);
mod = comm.FMBroadcastModulator('AudioSampleRate', fsA, 'SampleRate', fs, ...
    'Stereo', true, 'RBDS', true, 'RBDSSamplesPerSymbol', sps, 'FilterTimeConstant', 50e-6);
a = lowpass(noise(), 15e3, fsA);
a = max(min(a / (3 * std(a(:))), 1), -1);  % peaks at the 75 kHz deviation limit
x = mod(a, rbds());
end

function [x, fs] = genGPS()
fs = 10.23e6;
g = gpsWaveformGenerator('SampleRate', fs, 'PRNID', 7);
x = g(randi([0 1], 3, 1));                 % 3 navigation bits = 60 ms
end

function [x, fs] = genDVBT2()
% DVB-T2 8 MHz, 8K mode: T = 7/64 us, 6817 active carriers, guard interval 1/8.
nsym = 24;
m = comm.OFDMModulator('FFTLength', 8192, 'NumGuardBandCarriers', [688; 687], ...
    'InsertDCNull', false, 'CyclicPrefixLength', 1024, 'NumSymbols', nsym);
x = m(qammod(randi([0 63], 6817, nsym), 64, 'UnitAveragePower', true));
[x, fs] = txFilter(x, 64e6/7, 6817*64e6/7/8192, 4.2e6);   % mask: -50 dB at +-4.2 MHz
end

function [x, fs] = txFilter(x, fs, occupied, edge)
% Oversample x2 and apply a transmitter channel filter: flat over the
% occupied band, >= 55 dB down at the channel edge (real base stations
% meet ACLR ~45 dB; a raw OFDM waveform has a -35 dB pedestal instead).
x = resample(x, 2, 1); fs = 2*fs;
steep = 1 - (edge - occupied/2) / (fs/2 - occupied/2);
x = lowpass(x, occupied/2, fs, 'ImpulseResponse', 'fir', 'Steepness', min(steep, 0.995), ...
    'StopbandAttenuation', 60);
end

%% ===== plot ==============================================================

function plotTemplates(TPL, vis)
T = TPL.Templates;
fig = figure('Name', 'FieldFox signal templates', 'Color', 'w', 'Visible', vis, 'Position', [80 80 1300 800]);
tl = tiledlayout(fig, 'flow', 'TileSpacing', 'compact');
for k = 1:numel(T)
    ax = nexttile(tl);
    if isscalar(T(k).f)
        stem(ax, 0, 0, 'BaseValue', -80, 'LineWidth', 1.5);
        xlim(ax, [-1 1]);
    else
        plot(ax, T(k).f/1e6, T(k).psd_dB, 'LineWidth', 1);
        xlim(ax, T(k).f([1 end])/1e6);
    end
    ylim(ax, [-80 3]); grid(ax, 'on');
    title(ax, T(k).Name, 'FontSize', 8, 'FontWeight', 'normal');
    if T(k).OBW_Hz > 0, xregion(ax, [-1 1]*T(k).OBW_Hz/2e6, 'FaceAlpha', 0.1); end
end
xlabel(tl, 'Offset from centre (MHz)'); ylabel(tl, 'PSD (dB rel. peak)');
title(tl, 'Reference spectra (shaded = 99 % occupied bandwidth)');
end
