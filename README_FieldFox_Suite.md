# FieldFox Suite: MATLAB tools for a Keysight FieldFox N9914A

These tools help you study the RF interference from a solar-cell inverter. Each
one connects the FieldFox (serial MY61272405, firmware A.12.45, IP 10.90.36.200)
to MATLAB toolboxes. All of them also run **without the instrument**: on a demo,
on a saved recording, or on a Touchstone file. Try the offline version first.

Four further tools cover satellite communication: an SRS-4 CubeSat S-band downlink
measured with the FieldFox, pass planning, a CCSDS link simulator and 5G NTN. See
**section 9**.

```matlab
cd('C:\Users\WICHARIT\Documents\MATLAB')     % run everything from this folder
help ff_emi_study                            % every tool has full help text
```

No Keysight software or VISA is needed. The tools talk SCPI over the LAN
(socket 5025) with MATLAB's `tcpclient`, through the class `FieldFox.m`.

---

## 1. The tools at a glance

| Tool | What it does | FieldFox options | MATLAB toolboxes it combines | RF out? |
|---|---|---|---|---|
| `fieldfox_live` | Watches the SA trace live and records it. Read-only: it does not change any settings. | 233 | base MATLAB | no |
| `ff_signal_id` | Finds the signals in a spectrum and guesses what they are (LTE, NR, Wi-Fi, FM, GSM, BLE, GPS, CW...), with band and channel. | 233 (live read) | 5G, LTE, WLAN, Bluetooth, Communications, Satellite Comms, DSP System, Signal Processing, Statistics & ML, Deep Learning | no |
| `ff_spectrum_pro` | Controlled sweeps shown in MATLAB's `spectrumAnalyzer`. Inverter-OFF baseline, mask test, channel power / OBW / ACPR, and a cross-check against the FieldFox's own channel power and frequency counter. | 233 (+235) | DSP System, Signal Processing (+ `ff_signal_id`) | no |
| `ff_emi_study` | The main inverter study: OFF-ON-OFF scans from 150 kHz to 3 GHz plus GNSS L1, statistics per frequency, harmonic comb, zero-span bursts, indicative CISPR limits, interference range, PDF/DOCX report. | 233 (+235) | Statistics & ML, Bioinformatics (`mafdr`), Signal Processing, Communications (propagation models), Report Generator (+ `ff_signal_id`) | no |
| `ff_gnss_margin` | Shows how much the inverter's noise at GPS L1 (1575.42 MHz) lowers the C/N0 of every visible satellite, and how much margin is left. | 233 + 235 | Satellite Comms, Antenna, Phased Array, RF, Navigation, Signal Processing (+ `ff_signal_id`) | no (it switches the source off) |
| `ff_vna` | S-parameters of cables, filters and antennas: return loss, Smith chart, TDR, distance to fault, matching network, de-embedding, antenna model comparison. | 210, 211, 112 (QuickCal); 210+233 tracking generator | RF, Signal Integrity, Signal Processing, Antenna | **YES** |
| `ff_power_meter` | Logs the built-in channel power meter over time, with statistics and an SA cross-check. | 310 (+233) | Statistics & ML, Signal Processing | no |
| `ff_vvm_logger` | Logs magnitude and phase over time at one frequency (cable or phase stability): drift and Allan deviation. | 308 | Statistics & ML, Navigation (`allanvar`) | **YES** |
| `ff_instrument_record` | Long unattended recordings saved on the FieldFox itself (it keeps going if the PC sleeps), then downloads and parses them. | 236 (+233) | Signal Processing (+ `ff_signal_id`) | no |
| `analyze_fieldfox` | Script: offline OFF/ON comparison of files copied from the FieldFox (`.sta`, Record/Playback). Edit the settings at the top, then run it. | none (files) | base MATLAB | no |

Helpers used by the tools (you rarely call these yourself):
- `FieldFox.m`: the connection class.
- `ff_signal_templates.m` and `ff_signal_classifier.m`: the reference spectra and the trained models used by `ff_signal_id`.
- `ff_signal_templates.mat` and `ff_signal_models.mat`: cached results. **Do not delete them.** Rebuilding them takes about 1.5 min plus about 4 min.

Not used by the suite: the I/Q option (351, not installed), real-time SA (350, not installed), and the EMI/CISPR detectors (361, not installed). The on-instrument interference-analyzer screens and the field-strength and limit-line menus are not used either, because MATLAB does those jobs.

---

## 2. Safety: read this before you connect anything

- **Two tools transmit RF:**
  - `ff_vna`, default -10 dBm. NA mode transmits from port 1, and also from port 2 with `'Ports', 2`. The tracking generator transmits from port 1.
  - `ff_vvm_logger`, a CW signal from port 1.

  Both print a `*** RF OUTPUT ***` banner first. With `Prompt` on, they also ask you to confirm.
- **During those tests, disconnect anything that must not be driven:** a LISN, a limiter, a GNSS receiver or LNA, or an antenna you may not transmit with. Only the cable or device under test should be connected.
- **All the other tools only listen.** They take signals in on **port 2**, the SA input.
- **Damage levels:** the FieldFox ports are damaged by more than +27 dBm or 50 V DC. With a LISN, always put a transient limiter and a DC block in front of the FieldFox.
- **Switching the inverter:** use its own AC breaker and DC isolator. Never pull DC (PV) connectors under load.
- **Tools that take control** restore the FieldFox when they finish: `ff_spectrum_pro`, `ff_emi_study`, `ff_gnss_margin`, `ff_vna`, `ff_power_meter`, `ff_vvm_logger`, `ff_instrument_record` and `ff_srs4_analyzer`.
  - At the start they save the instrument's full state to `matlab_backup.sta` on the FieldFox.
  - When they finish, they load it back, and also after an error or Ctrl+C.
  - If MATLAB itself crashes, recall that state on the FieldFox: *Save/Recall > Recall State*.
- **Run one tool at a time.** Close `fieldfox_live` before you start a tool that takes control.

---

## 3. Quick start without the instrument

```matlab
R = ff_spectrum_pro("", 'Replay', "demo");             % spectrum + baseline mask on a synthetic "inverter"
R = ff_emi_study([], 'Data', "demo", 'Prompt', false); % full OFF-ON-OFF study + PDF report (~2 min)
T = ff_signal_id("fieldfox_rec_20261009_101500.mat")   % identify signals in any recording
R = ff_gnss_margin("", 'Lat', 13.73, 'Lon', 100.78);   % satellites + C/N0 margins (set YOUR site)
R = ff_vna([], 'Demo', "antenna");                     % S11 analysis on synthetic data
R = ff_power_meter('Demo', true);   R = ff_vvm_logger('Demo', true);
R = ff_instrument_record('parse', "<copied RecordPlayback folder>", 'RBW', 3e6);
```

For scripts and batch runs, add `'Visible', 'off', 'Prompt', false`. Then no window
opens and no dialog waits for you. Figures are still saved as PNG.

---

## 4. One-line live commands (FieldFox at 10.90.36.200)

| Goal | Command |
|---|---|
| Watch and record (read-only) | `fieldfox_live("10.90.36.200", 'Record', true)` |
| What signals are on screen now? (read-only) | `T = ff_signal_id("10.90.36.200", 'Plot', true)` |
| Inverter OFF then ON, live mask and measurements | `R = ff_spectrum_pro("10.90.36.200", 'Center', 95e6, 'Span', 20e6, 'RBW', 100e3, 'Detector', "AVER", 'AutoBaseline', 20)` |
| Full EMI study (dialogs guide OFF -> ON -> OFF) | `R = ff_emi_study("10.90.36.200", 'Sweeps', 20, 'Notes', "roof, 3 m, biconical")` |
| Edit the scan plan first | `P = ff_emi_study([], 'DryRun', true); B = P.bands([1 2 5], :); R = ff_emi_study("10.90.36.200", 'Bands', B)` |
| GPS L1 margin (antenna -> LNA -> port 2) | `R = ff_gnss_margin("10.90.36.200", 'Lat', 13.73, 'Lon', 100.78, 'Cable', "cable.s2p")` |
| Cable or filter S-parameters, 2-port QuickCal (**RF out**) | `R = ff_vna("10.90.36.200", 'Ports', 2, 'Cal', "quick2", 'Start', 1e9, 'Stop', 2e9)` |
| Antenna S11 vs an Antenna Toolbox model (**RF out**) | `R = ff_vna("10.90.36.200", 'Antenna', "dipole", 'DesignFrequency', 433e6)` |
| Channel power over 10 min | `R = ff_power_meter("10.90.36.200", 'Center', 100e6, 'ChannelBW', 1e6, 'Duration', 600)` |
| Cable phase stability for 30 min (**RF out**) | `R = ff_vvm_logger("10.90.36.200", 'Parameter', "S21", 'Frequency', 100e6, 'Duration', 1800)` |
| Overnight recording on the FieldFox | `R = ff_instrument_record("10.90.36.200", 'Start', 150e3, 'Stop', 30e6, 'RBW', 10e3, 'Duration', 8*3600, 'Interval', 10)` |
| Fetch a recording made earlier | `R = ff_instrument_record("10.90.36.200", 'Session', "latest")` |

---

## 5. A suggested workflow for the inverter study

1. **Survey.** Use `fieldfox_live` or `ff_signal_id` to see the ambient signals: FM, LTE/NR, TV, Wi-Fi.
2. **Characterise your set-up** (RF out, inverter off, nothing else connected):
   - measure your antenna cable with `ff_vna(..., 'Ports', 2)`;
   - save the `.s2p`, then use it later as `'Fixture'` in `ff_vna` or `'Cable'` in `ff_gnss_margin`.
3. **EMI study.** Run `ff_emi_study` with the inverter OFF -> ON -> OFF.
   - It flags the frequencies that rise significantly and finds the switching-harmonic comb.
   - Its zero span shows mains-related bursts.
   - It writes a report.
   - Add `'Transducer', "af.csv"` (antenna factor) for dBuV/m, or a LISN table for dBuV.
4. **Watch a band in detail.** Use `ff_spectrum_pro` with `'AutoBaseline'` (OFF) and then the inverter ON.
5. **Day-long behaviour.** Run `ff_instrument_record` overnight, then replay the result in `ff_spectrum_pro` or `ff_emi_study`.
6. **GNSS impact.** Use `ff_gnss_margin` (OFF, then ON) with your antenna, LNA and cable.
7. **Stability checks** if needed: `ff_power_meter` (channel power over time) or `ff_vvm_logger` (cable phase).

---

## 6. Sharing data between tools

Every tool that records sweeps writes a MAT file with the same variables as
`fieldfox_live`:

- `f`: frequencies in Hz (1 x N)
- `P`: levels in dBm (sweeps x N)
- `t`: datetime
- `settings`: at least `rbw`
- `idn`

So recordings move freely between tools. This was tested end to end on 2026-10-09.

| File | Made by | Can be replayed in |
|---|---|---|
| `fieldfox_rec_*.mat` | `fieldfox_live` | `ff_spectrum_pro('Replay' / 'BaselineFile')`, `ff_emi_study('Data', 1-3 files)`, `ff_signal_id(file)`, `ff_gnss_margin('Spectrum', for a ~20 MHz span around L1)` |
| `ff_spectrum_pro_*_rec.mat` | `ff_spectrum_pro` | same as above |
| `emi_study_*\rec\rec_<state>_b<band>_<det>.mat`, `emi_study.mat` | `ff_emi_study` | same; the GNSS-band `AVER` files also work in `ff_gnss_margin`; the whole folder works in `ff_emi_study('Data', folder)` |
| `gnss_l1_*_off.mat` / `_on.mat` | `ff_gnss_margin` | `ff_gnss_margin('Spectrum', off, 'SpectrumOn', on)` and the spectrum tools |
| `ffrec_<session>.mat` | `ff_instrument_record` | the spectrum tools. Give `'RBW'` when you parse, otherwise `ff_signal_id` cannot be used. |
| `ff_vna_*.s1p` / `.s2p` | `ff_vna` | `ff_vna('File' / 'Fixture')`, `ff_gnss_margin('Cable', file or the ff_vna result R)`, any RF Toolbox function |
| `cpm_*.mat`, `vvm_*.mat` | `ff_power_meter`, `ff_vvm_logger` | only their own tool (`'File'`). They hold one frequency, not a spectrum. |
| `srs4_*\srs4_run.mat`, `srs4_*\rec\rec_<step>.mat` | `ff_srs4_analyzer` | `ff_srs4_analyzer('File')`, `ff_ccsds_linksim('Measured')`, `ff_5g_ntn('Interference')`. The `rec_*.mat` files also work in the spectrum tools. `rec_env` / `rec_burst` are zero span, so they are not spectra. |
| `cubesat_link_*\cubesat_link.mat` | `ff_cubesat_link` | `ff_ccsds_linksim('Pass')` |
| `ntn5g_*\ntn5g_run.mat` | `ff_5g_ntn` | `ff_5g_ntn('File')`: reuses the BLER curves |

How recordings are split into inverter states:
- `ff_emi_study` with **one** recording splits it at the largest time gap: before the gap is OFF, after it is ON.
- With **two** recordings: OFF, ON.
- With **three**: OFF, ON, OFF again.

---

## 7. Honest limitations: what the numbers can and cannot tell you

- **No quasi-peak or CISPR-average detector** (option 361 is not installed). Every limit comparison is an **indicative pre-scan**, not a compliance test.
  - The RBW is a 3 dB filter; CISPR uses 6 dB.
  - A roof or garden is not a test site.
  - A "pass" means something. A "fail" needs an accredited EMC lab.
- **No I/Q** (option 351 is not installed), so nothing can be demodulated.
  - `ff_signal_id` guesses a technology from the **width and shape** of the power spectrum only. Treat its answer as a hypothesis.
  - LTE and NR of the same bandwidth look almost identical.
  - Its accuracy (about 87-91 %) was measured on **synthetic** traces, so it is an upper bound.
  - "Unresolved" means too few points across the signal. "Unknown" means it matches none of the 22 known classes, which is typical for inverter noise.
- **Near field below about 30 MHz.** Closer than wavelength/2pi (1.6 m at 30 MHz, 320 m at 150 kHz), an antenna measures the near field. There the noise travels mainly on the cables, so use a LISN. The interference-range estimate is only made in the far field.
- **GPS almanac.** The almanac shipped with MATLAB is from **2020-01-11** (about 6.7 years old), so satellite positions and PRNs are only indicative. For the real sky, download a current SEM almanac (navcen.uscg.gov), a GPS TLE file (celestrak.org) or a RINEX navigation file, and pass it with `'Almanac'`.
- **Satellite C/N0 is a model.** It uses a typical EIRP, an ideal patch antenna and typical thresholds. Only the OFF -> ON change is blamed on the inverter.
- **With nothing connected** you measure only the FieldFox's own noise. The tools print notes that say so (for example "open port", "nothing connected", "less noise than the antenna chain must produce").
- **Recordings are not gap-free.** There is no real-time option (350), so short bursts between sweeps can be missed. Zero span in `ff_emi_study` shows bursts up to about RBW/2 only.
- **Absolute level accuracy** is about +-1.5 dB. Channel power from a trace is accurate only with the AVER (RMS) detector and a point spacing no wider than the RBW.
- **Firmware commands still to be checked on this unit.** Several SCPI commands are documented by Keysight but have not yet been run on firmware A.12.45: zero-span sweep time, VNA averaging restart, QuickCal steps, Record/Playback session commands, the CPM trigger and the VVM data format.
  - The tools check `SYST:ERR?` after each one.
  - If a command is rejected, they fall back to another method or skip that feature, and they tell you.
  - The tools were tested against simulated instruments ("mocks") and synthetic data. The first real-instrument runs follow the live test plan.
- **QuickCal is temporary.** A QuickCal made inside `ff_vna` is used only for that run, because the restore brings back the old state. To keep a cal for several DUTs, hold control yourself (see `help ff_vna`).

---

## 8. Troubleshooting

| Problem | What it means / what to do |
|---|---|
| `The FieldFox is in NA mode. Switch it to SA ...` | The read-only tools (`fieldfox_live`, live `ff_signal_id`) never change the mode. Select SA on the FieldFox. |
| `No reply to "..." within 10 s` (`FieldFox:timeout`) | Check the LAN cable and the IP address, and make sure no other program holds the connection. The instrument is restored anyway. |
| The first `ff_signal_id` call takes about 12 s | It loads the trained models (Deep Learning runtime). Later calls take well under 1 s. |
| `... is not a swept spectrum (1 frequency point)` | You gave a power-meter or VVM log to a spectrum tool. Open it with `ff_power_meter('File', f)` or `ff_vvm_logger('File', f)`. |
| `... has no settings.rbw` | The recording does not store its RBW. Use `ff_signal_id(f, P, rbw)`, or parse with `ff_instrument_record('parse', folder, 'RBW', ...)`. |
| The FieldFox looks different after a crash | Recall `matlab_backup.sta` on the FieldFox (Save/Recall), or run `ff = FieldFox("10.90.36.200"); ff.loadState("matlab_backup.sta")`. |
| Something to check offline before a measurement | Run the tool with `'Data', "demo"`, `'Replay', "demo"`, `'Demo', ...` or `'File', ...`. Every tool has an offline mode. |

---

## 9. Satellite / CubeSat / 5G tools (Satlab SRS-4 S-band downlink)

These four tools are about a CubeSat radio, the **Satlab SRS-4**:
- TX 2200-2290 MHz, BPSK/QPSK/8PSK, RRC roll-off 0.2/0.25/0.35, 20-33 dBm;
- CCSDS framing: ASM, RS(255,223), randomizer, convolutional code r 1/2 K 7.

Your setting `sl-radio --radio=srs4 --tx-rate=1000000 --rx-rate=1000000 --tx-gain=50` means 1 MBd.
One QPSK RS+CC frame is 2072 symbols = 2.072 ms. `--tx-gain=50` is a gain setting, not dBm.
Only `ff_srs4_analyzer` uses the FieldFox, and it only listens. The others are pure computation.

| Tool | What it does | FieldFox | MATLAB toolboxes it combines |
|---|---|---|---|
| `ff_srs4_analyzer` | Finds the SRS-4 downlink in 2200-2290 MHz. It measures:<ul><li>fc and frequency error (ppm), Rs, roll-off: raised-cosine fit with confidence intervals and bootstrap;</li><li>OBW 99 / 99.9 %, channel power, C/N0 and Eb/N0;</li><li>the ECSS OBW limit and the SFCG/ECSS mask;</li><li>2nd harmonic and spurious (SFDR 60 dBc);</li><li>zero-span envelope statistics (a BPSK vs QPSK/8PSK hint);</li><li>burst / frame timing and idle-frame lines;</li><li>EIRP and TX power from Friis.</li></ul>It writes a PDF report. | takes control, SA on Port 2, restores | Satellite Comms (`ccsdsTMWaveformGenerator` reference), Signal Processing, Optimization (`lsqcurvefit`), Statistics & ML (`nlparci`, bootstrap), Phased Array (`fspl`), Antenna (`patchMicrostrip`, `dipole`), Communications (Rapp PA), Report Generator |
| `ff_cubesat_link` | Passes of a CubeSat over your site, from a CelesTrak TLE (KNACKSAT-2 by default) or an orbit. Per pass it gives a time-varying link budget and margin for every SRS-4 configuration (requirements from datasheet Table 7), data volume, ITU-R P.618 losses, and a Doppler pre-compensation table. It can also write a FieldFox Doppler-tracking table. | none | Satellite Comms (`satelliteScenario`, `access`, `link`, `ebno`, `dopplershift`, `p618PropagationLosses`, `satelliteCNR`), Aerospace (`tleread`), Antenna, Phased Array, Communications (`berawgn`, `bercoding`), Statistics, Navigation (`skyplot`), Mapping, Report Generator |
| `ff_ccsds_linksim` | Monte-Carlo BER/PER of the real SRS-4 frame chain. It has an ideal receiver and a "realistic" one (AFC, Gardner timing, PLL, ASM sync), Doppler, phase noise and IQ impairments. It checks the result against datasheet Fig. 13 and Table 7. Two modes connect it to the other tools: `'Pass'` (PER and throughput over `ff_cubesat_link` passes) and `'Measured'` (BER/PER at the C/N0 measured by `ff_srs4_analyzer`). | none | Satellite Comms (`ccsdsTMWaveformGenerator`, `ccsdsRSDecode`), Communications (`vitdec`, `comm.SymbolSynchronizer`, `comm.CarrierSynchronizer`, `comm.CoarseFrequencyCompensator`, `berconfint`), DSP System, Signal Processing, Statistics, Parallel (optional) |
| `ff_5g_ntn` | 5G NR NTN downlink from a LEO pass over Bangkok:<ul><li>Doppler, delay and timing advance vs numerology;</li><li>PDSCH BLER through NTN-TDL/CDL channels, with and without UE Doppler pre-compensation;</li><li>link budget and throughput per pass.</li></ul>With `'Interference'` it also estimates **coexistence**: how much a measured SRS-4 spectrum leaks into the n256 channel next to it, as I/N, separation distance and throughput loss. | none (reads analyzer results / recordings) | 5G (full PDSCH chain, `nrTDLChannel`/`nrCDLChannel` NTN profiles), Satellite Comms, Aerospace (`lla2eci`), Mapping (`reckon`), Communications, Signal Processing, Parallel (optional), Report Generator |

### Offline first (no instrument; about 1-3 min each)

```matlab
q = {'Visible', 'off', 'Prompt', false};
R = ff_srs4_analyzer([], 'Demo', true, q{:});             % synthetic SRS-4 + FieldFox model, PDF report
C = ff_cubesat_link('Offline', true, q{:});               % cached / sample TLE, no internet
L = ff_ccsds_linksim('Pass', C, 'Measured', R, q{:});     % BER/PER curves + pass PER + measured margin
N = ff_5g_ntn('Interference', R, q{:});                   % 5G NTN pass + coexistence with that SRS-4 spectrum
N2 = ff_5g_ntn('File', N, 'EIRPDensity', 28, q{:});       % fast what-if: reuses the BLER curves (~40 s)
```

How the tools connect (all through saved files or the returned `R`):
- `ff_srs4_analyzer`: its `srs4_run.mat` (or its `R`) goes to `ff_ccsds_linksim('Measured')` and to `ff_5g_ntn('Interference')`.
- `ff_cubesat_link`: its `cubesat_link.mat` (or its `R`) goes to `ff_ccsds_linksim('Pass')`. The configuration (modulation, FEC, symbol rate) comes with it.
- The analyzer's `rec\rec_*.mat` files use the shared format. Replay them in `ff_spectrum_pro` or give them to `ff_signal_id`.
  - `ff_signal_id` has no CCSDS class: it calls the SRS-4 "LTE 1.4 MHz" or "Unresolved". Use `ff_srs4_analyzer` for this signal.

### Live (FieldFox at 10.90.36.200, SRS-4 transmitting)

```matlab
R = ff_srs4_analyzer("10.90.36.200", 'Distance', 3, 'TxAntennaGain', 6.1, 'RxAntennaGain', 0, ...
      'Prompt', false, 'Visible', 'off', 'OutDir', fullfile(pwd, 'srs4_ota'));      % over the air
R = ff_srs4_analyzer("10.90.36.200", 'ExtGain', -40.6, 'Prompt', false);           % conducted, 40 dB pad
C = ff_cubesat_link('Site', [13.7290 100.7750 5], 'FieldFoxTrack', true);          % live TLE (2 h cache)
```

Add `'ExpectedFreq', <tx.freq in Hz>` when you know the SRS-4 channel. Then the frequency error is judged against 0.5 / 2.5 ppm.

### SRS-4 measurement setup and SAFETY

- **Over the air:**
  - SRS-4 antenna, then **at least 1 m** of air (3 m is better), then an S-band antenna, then a short cable to FieldFox **Port 2**.
  - At 1 m, 33 dBm + 6 dBi - 39.5 dB free-space loss gives about 0 dBm at Port 2.
  - Enter the distance, the antenna gains and the polarisation loss.
  - Expect at least +-3 dB on EIRP (multipath, antennas).
- **Conducted:**
  - SRS-4 TX connector, then a **30-40 dB attenuator**, then a cable, then Port 2.
  - The first pad must be rated **>= 10 W** (the SRS-4 gives up to 2 W continuous).
  - Measure pad + cable S21 with `ff_vna` first and pass it as `'ExtGain'` (negative = loss). The SRS-4 must be disconnected during that step: `ff_vna` transmits.
  - **Never** connect the SRS-4 directly to the FieldFox. +33 dBm is above the +27 dBm damage level. The analyzer stops above +20 dBm at its input.
- **Never let the SRS-4 transmit into an open connector** (no antenna, no load). The reflected power can damage its power amplifier. Connect the antenna or the pad first, then enable TX.
- Transmitting in 2200-2290 MHz over the air needs a frequency authorisation. Follow your lab's rules (NBTC in Thailand), keep tests short, or use the conducted set-up.
- The FieldFox never transmits in these four tools.
- **Frequency claims:** the FieldFox internal reference is +-0.7 ppm, which is worse than the SRS-4 spec of 0.5 ppm. For ppm-level verdicts, lock the FieldFox to an external 10 MHz GPSDO. The tool reads `ROSC:SOUR?`.
- **No I/Q** (option 351). EVM, the constellation and the modulation order cannot be measured. The envelope statistics are only a hint, and `ff_ccsds_linksim('Measured')` *predicts* BER/PER from C/N0.

### Data packages (installed in `support\`)

| Folder | Used by | Notes |
|---|---|---|
| `support\p618\iturp618v14Maps.mat` (ITU-R P.618 maps, 337 MB) | `ff_cubesat_link`, `ff_5g_ntn` | Put on the path by `startup.m`. If it is missing, a fixed loss is used and a note says so. |
| `support\tle\<NORAD>.tle` (+ `.err`), `thai_cubesats.tle` | `ff_cubesat_link`, `ff_5g_ntn` | Shared CelesTrak cache, refreshed after 2 h. After a failed download there is no retry for 2 h. The sample file holds KNACKSAT-2 (67683) and NAPA-1 (46320) for `'Offline'`. |
| `support\linksim\linksim_cache.mat` | `ff_ccsds_linksim` | BER/PER points, created on the first run. Use `'Cache', "refresh"` to recompute. |
| aeroDataPackage ephemerides | `ff_cubesat_link('Eclipse', true)` | **Not installed.** The eclipse option fails cleanly in its own scenario; the rest of the run is unaffected. |

Limitations to keep in mind:
- `ff_cubesat_link` and `ff_5g_ntn` are **predictions**. They depend on TLE age, assumed antennas and G/T, and 3GPP TR 38.821 study values.
- `ff_ccsds_linksim` RS columns: datasheet Table 7 is 1-2.6 dB worse than the model. That gap is unexplained (open question).
- Coexistence with an over-the-air FieldFox spectrum can rarely prove I/N <= -6 dB. The analyzer floor is far above a handheld UE's noise density, and the tool says so.
