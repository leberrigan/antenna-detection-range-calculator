Here's a self-contained reference note — copy this to the other agent as-is.

---

## Receiver Sensitivity & Detection Threshold — Reference Note

**Context:** Motus-style wildlife tracking receiver, hardware = RFM69HCW (Semtech SX1231H core). Tags: CTT 434 MHz, ~-15.3 dBm EIRP, 433.95 MHz, 32-bit FSK @ 25 kbps, 124.2 kHz occupied bandwidth, transmit every 2-20 s.

### 1. Core relationship

```
SNR (dB) = RSSI (dBm) - NoiseFloor (dBm)
NoiseFloor (dBm) = -174 + 10·log10(B_Hz) + NF(dB)
```

- `-174 dBm/Hz` = thermal noise floor at 290 K (a physical constant, not hardware-specific)
- `B_Hz` = receiver's channel/filter bandwidth (NOT the transmitter's occupied bandwidth — these are often different and easy to conflate)
- `NF` = receiver noise figure (hardware-specific, from datasheet)

### 2. RFM69HCW-specific values (from the actual datasheet, not estimated)

| Parameter | Value | Source |
|---|---|---|
| NF (highest-gain LNA stage, G1) | 7 dB | Table 13, Receiver Performance Summary |
| **DemodSnr** (SNR the chip's own demodulator needs) | **8 dB** | §3.4.3.2, AGC Reference formula — chip vendor's own documented figure |
| Published sensitivity, 1.2 kbps / 5 kHz FDA | -118 dBm typ | Table 6 |
| Published sensitivity, 4.8 kbps | -114 dBm typ | Table 6 |
| Published sensitivity, 38.4 kbps / 40 kHz FDA | -105 dBm typ | Table 6 |
| With SensitivityBoost (RegTestLna=0x2D) | -120 dBm | Table 6 |
| RSSI dynamic range (AGC on) | -115 to 0 dBm | Table 6 |
| Bandwidth constraint | BitRate < 2 × RxBw | §3.4.6 — for 25 kbps CTT signal, RxBw must be ≥ 12.5 kHz |

**The headline number:** the chip's own datasheet defines its detection threshold as **8 dB SNR** in whatever channel bandwidth it's configured for. That's a vendor-documented constant, not something to derive — use it directly if you need a "does this look like a detectable signal" heuristic.

### 3. The chip's own sensitivity formula (useful for a bandwidth-dependent threshold)

```
Sensitivity (dBm) = -174 + NF + DemodSnr + 10·log10(2·RxBw) + FadingMargin
                   = -174 + 7 + 8 + 10·log10(2·RxBw) + 5
```

`RxBw` is the single-sideband receiver channel filter width (a config register on the chip, one of a fixed discrete set — not free to pick any value). Since RxBw is the main unknown, tabulate a few options rather than assuming one:

| RxBw (SSB) | Sensitivity, no margin | Sensitivity, +5 dB fading margin |
|---|---|---|
| 125 kHz (≈ matches tag's 124.2 kHz occupied BW) | -105.0 dBm | -100.0 dBm |
| 62.5 kHz | -108.0 dBm | -103.0 dBm |
| 25 kHz | -112.0 dBm | -107.0 dBm |
| 12.5 kHz (tightest legal, per BitRate<2×RxBw) | -115.0 dBm | -110.0 dBm |

### 4. Important caveat for a real Motus/SensorGnome-style station

Many of these receivers don't actually use the RFM69's on-chip FSK packet demodulator — they run the radio in raw/continuous RSSI or I/Q mode and do custom software correlation against the known tag pulse pattern in the host software. If that's the case here:
- the 8 dB DemodSnr figure doesn't directly apply (custom correlation can beat it via processing gain)
- the empirically-observed RSSI floor in real field data is the better ground truth than any datasheet number
- **check field-observed RSSI floor against the tightest-bandwidth row above (12.5 kHz, ~-115 dBm) before trusting the datasheet numbers wholesale** — in this project's data that row matched almost exactly, suggesting either a narrow effective filter or correlation gain roughly compensating for the wider one

### 5. What NOT to compute from raw RSSI logs alone

RSSI is normally only logged for packets that successfully decoded — the distribution is **left-censored**. A KDE/histogram of logged RSSI will show the population of successful detections, not the sensitivity floor directly. The left tail edge is a reasonable empirical proxy for the floor, but a rigorous number needs paired success/miss data (known tag transmit schedule vs. what was actually logged) or a bench sweep.

### 6. Suggested inputs to expose in a UI doing this calculation

- Receiver channel bandwidth (RxBw) — dropdown of the RFM69's discrete legal values, not a free-text field
- Noise figure — default 7 dB (G1), let user override if a different gain stage is relevant
- Demodulator SNR requirement — default 8 dB (chip spec), override for custom software demod
- Fading margin — default 5 dB, or 0 for "raw" chip sensitivity
- Show both the *computed* number (formula above) and the *empirical* number (from field RSSI data) side by side — they answer different questions and shouldn't be collapsed into one value

---

Source: RFM69HCW-V1.1 datasheet (SparkFun/HopeRF), Table 6, Table 13, Table 14, §3.4.3.2, §3.4.6.