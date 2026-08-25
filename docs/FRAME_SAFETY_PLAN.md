# Frame-Safety Architecture: Manifest & Integration Plan

## The Problem

Every bug we've encountered maps to a missing frame distinction:

| Bug | Missing frame distinction |
|---|---|
| GLONASS time offset | Temporal: GLONASST ≠ GPST (−10782 s = −3 h + 18 leap s) |
| BeiDou time offset | Temporal: BDT ≠ GPST (+14 s) |
| λ_IF confusion | Signal-frame: IF wavelength (10.7 cm) ≠ wide-lane (86.2 cm) |
| Galileo band-2 E5a/E5b | Frequency-frame: same band number, different physical signal |
| CAPO vertical bias | Antenna-frame: ARP ≠ phase center (39.7 mm Up) |
| UNR truth offset | Reference-frame: IGS20 ≠ broadcast WGS84 realization |
| OHLN coastal v_p50 | Physical-model: ocean tide loading absent |
| Gradient sign | Azimuthal-frame: az from North CW vs from East CCW |
| Relational coupling (PCV) | Cross-station zenith angles passed as independent floats, allowing physically impossible geometry; caught in minutes by review, not weeks by benchmark anomaly (3 orders of magnitude cheaper detection) |

## Three-Track Solution (parallel implementation)

### Track A: Time-system types (`gnss_time.rs`)
- `TimeSystem` enum: Gps / Glonass / Bdt / Gst
- `GnssTime` with `sys` tag + `to_gpst()` as the ONLY cross-system path
- Replaces ad-hoc `+18.0 - 10800.0` and `+14.0` scattered in parsers
- **Red-team focus**: leap-second table management, sign conventions, week boundaries

### Track B: Reference-frame positions (`frames.rs`)
- `ReferenceFrame` trait + concrete frames (ITRF2014, IGS20, WGS84)
- `EcefPos<F>` newtype with frame-safe operations
- Helmert conversions between frames, epoch-propagated
- **Red-team focus**: Helmert parameter signs, epoch propagation formula, sub-mm precision

### Track C: Signal-frequency registry (`frequencies.rs`)
- `Signal` enum identifying physical signals unambiguously
- `(Constellation, RINEX_type)` → `Signal` mapping (resolves E5a/E5b once)
- GLONASS FDMA k-dependence built into frequency lookup
- **Red-team focus**: ICD frequency values, band-numbering conventions, FDMA edge cases

## Integration Plan (after all three tracks complete)

1. Replace ad-hoc corrections in parsers with Track A types
2. Replace raw Vector3<f64> position passing with Track B types at API boundaries
3. Replace get_frequency() calls with Track C lookups
4. End-to-end test: parse → SPP → DD → AR → output, asserting frame consistency

## What This Does NOT Cover (future phases)

- Attitude/body frames (satellite PCO orientation) — needs orbit-frame computation
- Atmospheric mapping function selection (NMF/GMF/VMF1) — separate abstraction
- Ocean tide loading models — external data dependency
- Tropo gradient convention standardization (currently cot(el), could be CH97)

## Success Criteria

- Zero ad-hoc numeric offsets outside of tested conversion functions
- Compiler prevents mixing incompatible frames
- All frequency/band lookups go through one authoritative table
- Every conversion function is independently verified against published standards
