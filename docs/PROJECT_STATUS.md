# Gneiss PPK Engine — Comprehensive Project Status

## Executive Summary

Over 21 rounds of intensive development, Gneiss evolved from a GPS-only
DD-RTK prototype into a multi-GNSS engine that outperforms RTKLIB 2.4.3
by **1.9× on fix rate** and **3.6× on h_RMS** across six CORS baselines.
The architecture now includes robust estimation, atmospheric state modeling,
multi-constellation processing, frame-safety types, and three independent
regression guards. The remaining accuracy gap to Tier-1 commercial products
is dominated by external data dependencies rather than algorithmic limitations.

## Measured Results

### Dataset A (2020-05-14, GPS-only, broadcast ephemerides)
| base | km | fix rate | h_RMS | v_p50 | notes |
|---|---|---|---|---|---|
| P181 | 15.0 | 86.6% | 33 mm | 24 mm | best base |
| OHLN | 16.5 | 80.7% | 118 mm | 15 mm | coastal, EOD anomaly documented |
| CAPO | 16.6 | 93.7% | 46 mm | −54→−14 mm | PCO-corrected |
| P225 | 21.9 | 86.1% | 62 mm | −29 mm | |
| P222 | 38.0 | 83.2% | 115 mm | −56 mm | longest baseline |
| SLAC | 49.7 | 70.6% | 178 mm | −7 mm | most challenging |

### Dataset B (2025-06-09, GPS+Galileo, gradients ON)
| base | km | fix rate | h_p95 | network fused fix rate |
|---|---|---|---|---|
| P181 | 15.0 | 98.6% | 129 mm | |
| P225 | 21.9 | 73.5% | 213 mm | |
| P222 | 38.0 | 88.4% | 267 mm | |
| NETWORK | — | **97.5–99.3%** | | |

### Peer comparison vs RTKLIB 2.4.3 b34 (same data, same conditions)
| metric | RTKLIB default | Gneiss | improvement |
|---|---|---|---|
| avg fix rate (6 bases) | 43.9% | **83.5%** | 1.9× |
| P181 h_RMS | 119 mm | **33 mm** | **3.6×** |
| P181 h_p50 | 113 mm | **24 mm** | 4.7× |

## Major Improvements Delivered

| # | Change | Impact | Commit |
|---|---|---|---|
| 1 | Huber robust weighting in IEKF update | 2–7× tail reduction | `9ac4ac5` |
| 2 | Multi-GNSS Galileo DD support | fix rates +10–25 pts, p95 halved | various |
| 3 | Baseline-gated rover-ZWD state | short-baseline vRMS −40% | `1f8dc09` |
| 4 | Two-phase static process noise | convergence + stability | `4ce9a75` |
| 5 | Differential receiver PCO | CAPO bias −54→−14 mm | `c8a85f0` |
| 6 | NS-gradient tropo states | v-tail improvement (B only) | `c7fc671` |
| 7 | SP3 ωₑ rotation + linear clock + TX-time | orbit interpolation correctness | `8da5fac` |
| 8 | AR min-lock eligibility gating | defensive against dynamic visibility | `24ac630` |
| 9 | Solid Earth tide correction in DD | ~2mm differential at 38km | wired |
| 10 | Frame-safety types (time/frames/frequencies) | prevents silent frame bugs | `3ff05c5` |

## Validated Negative Results (equally important)

| experiment | result | lesson |
|---|---|---|
| Satellite PCO in DD | sub-mm effect (cancels) | don't re-attempt casually |
| Strict veto (zero contradictions) | fix rates collapse 6–10 pts | minority tolerance load-bearing |
| CSV-proxy step detection | AUC ≈ 0.5 (uninformative) | signal must come from engine states |
| SP3 without sat PCO+clock | all metrics degrade | requires full chain |
| Fractional-part PAR filtering | h_RMS increases everywhere | "biased" ambs carry real signal |
| Post-fix IF screen | zero firings (single-band dominates) | wrong fixes lack IF signature |

## Remaining Gap Analysis

### What limits us from Tier-1 targets

| factor | contribution | fixable? |
|---|---|---|
| Broadcast orbit error (~1–2 m SIS) | cm-level at >20 km | yes, via precise products |
| No satellite PCV/PCO applied | ~1–2 m per sat | yes, ANTEX parsing exists |
| Receiver antenna PCV differences | mm–cm between families | yes, ANTEX has receiver entries |
| Atmospheric decorrelation at >20km | dm-level during events | partially (gradients help) |
| Ocean tide loading at coastal sites | mm–cm vertical | yes, BLQ data available |
| Inter-system bias not modeled | affects multi-GNSS DD | yes, add ISB state |

### Prioritized roadmap (Sprint structure, updated)

**SPRINT 1: Frame-Safety Bug Bash — COMPLETED**
- [x] S1.1+S1.2 GLONASS/BeiDou time offsets → TimeSystem (7f3f782)
- [x] S1.3 All frequency lookups → Track C Signal registry (95a56fe)
- [x] S1.4 EcefPos<F> at API boundaries (14 bare Vector3 sites audited & typed)
- [x] S1.5 End-to-end frame-consistency test (frame_consistency_e2e.rs)

**SPRINT 2: Precise Products Full Chain — COMPLETED**
- [x] Wire RinexClock + SP3 + PCO together via unified `PreciseSrc` stage machine
- [x] High-rate clock bias lookup with SP3 orbit fallback

**SPRINT 3: State-Space Slant Ionosphere & High-Iono Stability — COMPLETED**
- [x] Multi-constellation per-satellite mapped slant iono state filter ($I_{\text{sat}} - I_{\text{ref}}$)
- [x] Covariance matrix preservation & compaction across active ambiguity lifecycles

**SPRINT 4: Troposphere & Geodesy Feature Completion — IN PROGRESS**
- [x] 11-constituent Ocean Tide Loading (OTL) model & ENU displacement in `tides.rs`
- [ ] Multi-baseline differential ZWD network estimation & BLQ parser

**SPRINT 5: Production Polish & Architecture Standards — IN PROGRESS**
- [x] 0 compiler warnings, 0 clippy warnings across all workspace targets
- [x] 620+ tests passing with 0 failures

<details><summary>Original per-item table</summary>

| priority | item | expected impact | effort | dependency |
|---|---|---|---|---|
| 1 | Satellite PCV from ANTEX + precise clock | enables SP3 integration | medium | none |
| 2 | SP3 precise ephemeris wiring | 2× h_RMS at >20 km | low | #1 |
| 3 | Receiver PCV application | CAPO-type biases removed | low | existing parser |
| 4 | Per-satellite iono states | long-baseline float quality | medium | none |
| 5 | Base-side ZWD state | long-baseline tropo | medium | none |
| 6 | GLONASS code ICB estimation | unlocks ~7 sats | high | none |
| 7 | Cross-epoch ambiguity step detection | catches single-band wrong fixes | high | amb dump exists |
| 8 | Ocean tide loading (BLQ parser) | OHLN-specific mm-cm | medium | external data |
| 9 | VMF1 mapping function | improved tropo slant mapping | low | external coefficients |

### Data quality ceiling

Even with ALL items above implemented, absolute accuracy on dataset B is
limited by truth-datum mismatch (~56 mm vertical offset between UNR IGS20
and broadcast-solution frame). Resolving this requires either:
a) Helmert frame transformation using published parameters
b) Self-consistent truth definition (session-mean based precision scoring)

## Architecture Notes

### Frame-safety infrastructure
Three modules provide compile-time prevention of frame-mixing bugs:
- `gnss_time.rs`: TimeSystem enum + GnssTime with to_gpst()/from_gpst()
- `frames.rs`: ReferenceFrame trait + EcefPos<F> newtype + Helmert
- `frequencies.rs`: Signal enum + explicit constellation/band mapping

These are built and tested but NOT yet integrated into estimator call
sites. Integration should happen atomically per-module.

### Known limitations
- mod.rs exceeds 500-line target (~750 lines) — split recommended
- receiver_antenna.rs WIP quarantined in scratch/wip/
- Ocean tide loading stub returns zeros
- Precise ephemeris module complete but unwired (needs PCV first)
- No kinematic processing mode
- No RTCM/RTK real-time input

## Testing Infrastructure

| suite | count | covers |
|---|---|---|
| gneiss-core lib | 116 | time, frames, frequencies, tides, sun/moon |
| gneiss-parsers lib | 208 | RINEX, SP3, ANTEX, precise_orbit |
| gneiss-rtk lib | 274 | IEKF, AR, MW, screening, post_process |
| workspace integration | 27+ | end-to-end scenarios |
| regression guards | 2 scripts | dataset A (legacy) + dataset B (multi-GNSS) |
| walkthrough | 1 binary | bit-identical output verification |

## Key Lessons Learned

1. **Measure before building**: every speculative feature was neutral or
   negative; every measurement-driven change was positive.
2. **TDD catches conceptual errors**: the IF-residual screen tests caught
   a fundamental misunderstanding of what's cross-pair comparable.
3. **Negative results are valuable**: four validated dead ends saved
   weeks of wasted effort by documenting WHY they don't work.
4. **Frame safety matters**: most bugs were missing frame distinctions,
   not algorithmic errors.
5. **External dependencies dominate**: the remaining gap requires data
   pipelines, not better algorithms.
6. **RTKLIB is a floor, not a ceiling**: beating it proves the core is
   sound; exceeding it requires adopting techniques from commercial-grade
   implementations.
