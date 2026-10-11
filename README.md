# MATLAB-Radio-Research

MATLAB tools for RF, wireless and satellite-communication research with a
**Keysight FieldFox N9914A** handheld analyzer (SA, 2-port VNA, QuickCal,
vector voltmeter, channel power meter, Record/Playback), combined with MATLAB's
high-level toolboxes: Satellite Communications, 5G, LTE, WLAN, Bluetooth,
Communications, DSP, Signal Processing, RF, Antenna, Signal Integrity,
Statistics and Machine Learning, Deep Learning and Report Generator.

The full user guide is in [README_FieldFox_Suite.md](README_FieldFox_Suite.md).

## Quick start

```matlab
cd <this folder>
fieldfox_studio            % launcher with every tool (connects only when you press Connect)
```

Every tool also works offline (`'Demo'`, `'File'` or `'Replay'` options), so it can be
used without the instrument. Tools that control the analyzer save its full state
first and restore it afterwards (including after a crash or a killed MATLAB).

## Tools

| Area | Tools |
|---|---|
| Core | `FieldFox.m` (SCPI over LAN, no VISA needed), `fieldfox_studio.m` (launcher) |
| Spectrum monitoring | `fieldfox_live`, `ff_spectrum_pro` |
| Automatic survey and EMI | `ff_auto_survey`, `ff_emi_study` (solar-inverter OFF/ON study with statistics and PDF report) |
| Signal identification | `ff_signal_id` (5G/LTE/WLAN/Bluetooth/GSM/FM templates + ML classifier), `ff_5g_analyzer` (NR/LTE carriers, TDD pattern from zero span, ACLR) |
| Satellite / CubeSat | `ff_satsig` (signal library incl. OrigamiSat-1/2, KNACKSAT-2; overhead prediction; automatic Doppler-based identification), `ff_sat_tracker` (Doppler-tracking capture), `ff_srs4_analyzer` (Satlab SRS-4 S-band downlink vs datasheet), `ff_cubesat_link` (passes, link budget, Doppler table), `ff_ccsds_linksim` (CCSDS link-level Monte Carlo), `ff_gnss_margin` (GPS L1 interference margin) |
| 5G NTN | `ff_5g_ntn` (5G NR from LEO: geometry, NTN channel link level, coexistence) |
| VNA / antennas / meters | `ff_vna` (S-parameters, Smith chart, TDR, matching network, Antenna Toolbox comparison), `ff_power_meter`, `ff_vvm_logger`, `ff_instrument_record` |
| Offline analysis | `analyze_fieldfox` (FieldFox .sta / Record-Playback files) |

## Requirements

- MATLAB R2026a with the toolboxes listed above (each tool checks what it needs).
- Keysight FieldFox on the LAN (default `10.90.36.200`, SCPI socket 5025). No VISA library is required.
- Optional data: the ITU-R P.618 digital maps (download from MathWorks into `support/p618`, added to the path by `startup.m`).

## Safety

The FieldFox input (Port 2) is damaged above +27 dBm or 50 V DC. Use attenuators for
conducted measurements of transmitters (for example the SRS-4, up to +33 dBm) and never
let a transmitter work into an open connector. VNA, vector voltmeter and tracking
generator modes transmit RF from the analyzer ports.

## Limits

The N9914A has no I/Q capture option, so signals are identified and characterised from
power spectra and zero-span envelopes; nothing is demodulated, and EVM cannot be measured.
The EMI study is an indicative pre-scan (no CISPR quasi-peak detector).

## License

MIT (see [LICENSE](LICENSE)).
