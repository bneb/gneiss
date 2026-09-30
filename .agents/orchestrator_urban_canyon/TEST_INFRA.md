# Test Infrastructure: Urban Canyon Fix Rate Expansion and Multipath Mitigation

**Document Version**: 1.0.0  
**Target Architecture**: Gneiss GNSS/INS Positioning Engine (`crates/gneiss-rtk`)  
**Scope**: 4-Tier Opaque-Box E2E Integration Test Suite for R1–R4  
**Date**: 2026-09-24  

---

## 1. Overview & Objectives

In dense urban canyons (such as Tokyo Shinjuku and Hong Kong Whampoa), GNSS signals suffer from severe pseudorange multipath (5m–20m code steps), Non-Line-of-Sight (NLOS) reflections, low signal-to-noise ratios (C/N0 < 25 dB-Hz), and rapid cycle slips caused by building and overpass occlusions.

The purpose of this E2E test infrastructure is to provide an uncompromising, mathematically grounded test suite that verifies:
1. **R1: Adaptive C/N0 (SNR) and Elevation Observation Covariance Weighting** — Smooth, monotonic covariance scaling preventing Kalman filter gain chatter.
2. **R2: Code-Minus-Carrier (CMC) Multipath Detection & Down-Weighting** — Decoupled pseudorange blunder suppression preserving clean millimeter carrier phase measurements and shielding wide-lane tracking.
3. **R3: Doppler-Assisted Cycle Slip Detection & Phase Continuity Validation** — Instantaneous detection of 1-cycle and half-cycle slips under simulated vehicle dynamics, re-seeding ambiguity states and resetting lock counters.
4. **R4: C/N0- and Elevation-Prioritized Partial Ambiguity Resolution (PAR)** — Physical quality metric ranking of candidate subsets with DOP geometry guards and positive semi-definite covariance conditioning ($Q_{aa} \succeq 0$, $P \succ 0$).

---

## 2. Four-Tier Test Suite Architecture

The test suite is structured into four distinct verification tiers:

```
┌────────────────────────────────────────────────────────────────────────┐
│                        4-TIER TEST SUITE MATRIX                        │
├────────────────────────────────────────────────────────────────────────┤
│ Tier 1: Feature Coverage (>=5 tests per feature, >=20 tests total)     │
│   • R1: Elevation scaling, SNR attenuation, DD covariance definiteness │
│   • R2: CMC arc residual tracking, code inflation, phase preservation  │
│   • R3: Doppler 1-cycle slip, half-cycle slip, time gap detection      │
│   • R4: PAR subset ordering, DOP guard, covariance definiteness        │
├────────────────────────────────────────────────────────────────────────┤
│ Tier 2: Boundary & Corner Cases (>=5 tests per feature, >=20 tests)    │
│   • R1: Horizon limit (el < 5°), extreme low SNR (< 15 dB-Hz), bounds  │
│   • R2: Extreme 20m code steps, zero multipath noise, blunder caps     │
│   • R3: Dynamic acceleration, high Doppler noise, short/long cadence   │
│   • R4: Minimal subset size (k=4), near-collinear satellites, DOP edge │
├────────────────────────────────────────────────────────────────────────┤
│ Tier 3: Cross-Feature Pairwise Interactions (>=6 tests)                │
│   • R1 + R2: SNR attenuation combined with CMC multipath step          │
│   • R2 + R3: Code multipath step vs Doppler cycle slip discrimination  │
│   • R2 + R4: CMC-de-weighted satellites omitted from PAR fixing        │
│   • R3 + R4: Slipped satellites lock-reset and excluded from PAR       │
│   • R1 + R4: Multi-constellation DD formation and subset selection     │
│   • R1 + R2 + R3 + R4: Full-pipeline innovation and conditioning cycle │
├────────────────────────────────────────────────────────────────────────┤
│ Tier 4: Real-World Urban Canyon Mission Scenarios (>=5 scenarios)      │
│   • Scenario 1: Tokyo Shinjuku Skyscraper Canyon (deep multipath)      │
│   • Scenario 2: Hong Kong Whampoa High-Rise Urban Canyon (rapid slips) │
│   • Scenario 3: Highway Overpass Outage & Rapid Reacquisition          │
│   • Scenario 4: Asymmetric CORS Base vs Urban Canyon Rover Noise       │
│   • Scenario 5: Collinear Street Canyon Geometry & DOP Degeneracy      │
└────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Feature Inventory & Mapping

| Feature ID | Requirement | Primary Module | Mathematical Invariant / Contract |
|---|---|---|---|
| **F-R1-1** | Elevation Variance Scaling | `gneiss_core::variance` | $\frac{\partial \sigma^2}{\partial \theta} \le 0$, finite at $\theta \to 0$ |
| **F-R1-2** | C/N0 Exponential Scaling | `gneiss_core::variance` | $\frac{\partial \sigma^2}{\partial S} \le 0$, $f_{\text{SNR}}(S) \le f_{\max} = 1000$ |
| **F-R1-3** | Smooth DD Covariance | `rtk_iekf::formation_cov` | $R_{DD} = D R_{\text{undiff}} D^T \succ 0$, $\lambda_{\min}(R_{DD}) > 0$ |
| **F-R1-4** | Continuous Derivatives | `gneiss_core::variance` | $C^1$-smooth transition across threshold ($40\text{ dB-Hz}$) |
| **F-R1-5** | Float SNR Representation | `gneiss_core::obs` | Fractional C/N0 eliminates $1\text{ dB}$ quantization step |
| **F-R2-1** | CMC Residual Tracking | `rtk_iekf::formation` | $CMC = P - \lambda \Phi - 2I$, baseline frozen on step |
| **F-R2-2** | Decoupled PR Screening | `rtk_iekf::screen` | Prefit $|y_P| > 15\text{m}$ suppresses code row, retains phase row |
| **F-R2-3** | Adaptive Code Inflation | `rtk_iekf::formation_cov` | $\sigma_P^2 \leftarrow \sigma_{P,0}^2 + \sigma_{\text{mp}}^2$, $\sigma_\Phi^2$ unchanged |
| **F-R2-4** | MW Tracker Shielding | `rtk_iekf::mw` | Code multipath step does not trigger false MW arc reset |
| **F-R2-5** | Ambiguity Seed Protection | `rtk_iekf::state` | Slipped/new ambiguity seed variance inflated under code error |
| **F-R3-1** | Doppler 1-Cycle Slip | `post_process::screening` | $|\Delta\Phi + f_D \Delta t| > 0.30\text{ cyc} \implies \text{Slip}$ |
| **F-R3-2** | Doppler Half-Cycle Slip | `post_process::screening` | $0.50\text{ cyc}$ slip detected at $6\sigma$ clearance |
| **F-R3-3** | Multi-Band Slip Screening | `post_process::screening` | Checks bands 1, 2, 5, 6, 7 across GPS, Gal, BDS |
| **F-R3-4** | Unconditional Base Slip | `rtk_iekf::formation` | Base slips evaluated unconditionally and propagated to DD |
| **F-R3-5** | Tracking Lock Reset | `rtk_iekf::state` | Lock epoch counter reset to 0 upon cycle slip |
| **F-R4-1** | Multi-Metric PAR Ranking | `ambiguity::par`, `ar_subsets` | Candidates sorted by Composite Quality Metric (CQM) |
| **F-R4-2** | Geometry & DOP Guard | `ambiguity::par`, `dop` | Subsets with $\text{PDOP} > 6.0$ or $k < 4$ rejected |
| **F-R4-3** | Sub-Covariance Definiteness | `ambiguity::par` | $Q_{\text{sub}} \succ 0$, $\lambda_{\min}(Q_{\text{sub}}) \ge 10^{-9}$ |
| **F-R4-4** | Zero Off-Diagonal Leakage | `rtk_iekf::ar` | Cross-covariances $P_{ia} = 0, P_{ai} = 0$ for fixed $a$ |
| **F-R4-5** | Conditional State Definiteness | `rtk_iekf::ar` | $\lambda_{\min}(P) \ge 10^{-6}$ maintained after integer conditioning |

---

## 4. Test Suite Code Layout

To comply strictly with `AGENTS.md` (all files < 500 LOC, all functions <= 32 LOC, nesting < 3 levels), the test suite is partitioned into focused modular files:

```
crates/gneiss-rtk/tests/
├── test_urban_canyon_e2e.rs         # Test runner entry point (< 100 LOC)
└── urban_canyon/
    ├── mod.rs                        # Module re-exports (< 50 LOC)
    ├── common.rs                     # Synthetic observation fixtures (< 350 LOC)
    ├── tier1_features.rs             # Tier 1 tests: R1, R2, R3, R4 (< 450 LOC)
    ├── tier2_boundaries.rs           # Tier 2 tests: boundary & corner cases (< 450 LOC)
    ├── tier3_interactions.rs         # Tier 3 tests: cross-feature pairwise (< 400 LOC)
    └── tier4_scenarios.rs            # Tier 4 tests: realistic mission scenarios (< 450 LOC)
```

---

## 5. Scenario Definitions (Tier 4)

### Scenario 1: Tokyo Shinjuku Skyscraper Canyon
- **Environment**: Narrow street surrounded by 200m skyscraper glass facades.
- **Signal Conditions**:
  - 2 zenith satellites: direct LOS, $\theta > 65^\circ$, $S = 46\text{ dB-Hz}$, clean code & phase.
  - 2 mid-elevation satellites: reflected from skyscraper, $\theta = 35^\circ$, $S = 24\text{ dB-Hz}$, $+14.0\text{ m}$ code multipath step.
  - 2 low-elevation satellites: heavily obstructed, $\theta = 12^\circ$, $S = 18\text{ dB-Hz}$.
- **Expected Outcome**:
  - CMC detector identifies $+14.0\text{ m}$ code jump and inflates code variance.
  - Prefit screening suppresses contaminated code rows without deleting carrier phases.
  - PAR fixes ambiguities on high-quality zenith and clean satellites first.
  - Zero false fixes, positioning error $< 1.0\text{ m}$.

### Scenario 2: Hong Kong Whampoa High-Rise Urban Canyon
- **Environment**: High-density residential towers with balcony overhangs and street foliage.
- **Signal Conditions**: Frequent signal interruptions, half-cycle ($0.5$ cyc) and 1-cycle slips on alternating satellites every 3–5 epochs.
- **Expected Outcome**:
  - Doppler detector flags half-cycle slips instantaneously ($0.50$ cyc discrepancy).
  - Ambiguity variances on slipped satellites are immediately re-seeded.
  - Lock duration counters reset to 0, preventing premature fixed subset inclusion.
  - Continuous tracking maintained on clean satellites with zero false integer fixes.

### Scenario 3: Highway Overpass Outage & Rapid Reacquisition
- **Environment**: Elevated road / rail underpass resulting in a complete 3-second GNSS blackout.
- **Signal Conditions**: 0 satellites visible for $\Delta t = 3.0\text{ s}$. Signals re-emerge with low initial C/N0 ($28\text{ dB-Hz}$) and high initial code variance.
- **Expected Outcome**:
  - Cycle slip detector triggers time-gap reset ($\Delta t > 2\text{ s}$).
  - Ambiguities are cleanly initialized without stale covariance leakage.
  - Filter safely re-converges to float solution before smoothly locking integers.

### Scenario 4: Asymmetric CORS Base vs Urban Canyon Rover Noise
- **Environment**: Rover in deep canyon, base station on open-sky geodetic monument.
- **Signal Conditions**:
  - Base station: all satellites $S > 45\text{ dB-Hz}$, $\sigma_{\text{code}} = 0.20\text{ m}$.
  - Rover: satellites attenuated down to $20\text{--}25\text{ dB-Hz}$.
- **Expected Outcome**:
  - Error propagation evaluates rover and base single-difference noise independently.
  - Base station reference covariance contribution is not artificially corrupted by rover attenuation.
  - Double-difference covariance matrix $R_{DD}$ is strictly positive definite.

### Scenario 5: Collinear Street Canyon Geometry & DOP Degeneracy Guard
- **Environment**: Narrow urban canyon running strictly North-South; buildings block East-West satellites.
- **Signal Conditions**: All 5 visible satellites lie along azimuths $0^\circ$ and $180^\circ$ (collinear geometry, East DOP $> 12.0$).
- **Expected Outcome**:
  - PAR candidate subset geometry check evaluates PDOP and directional DOP.
  - Geometry guard rejects ill-conditioned subset fixing.
  - Engine remains safely in float mode, preventing geometric magnification of residual errors.

---

## 6. Coverage Thresholds & Quality Invariants

- **Total Test Count**: $\ge 50$ tests across Tiers 1–4.
- **Pass Rate**: 100% (zero failures, zero ignored).
- **Compiler & Clippy Warnings**: Exactly 0 under `cargo clippy --workspace --all-targets -- -D warnings`.
- **Production Code Modularity**: File size < 500 LOC, function size <= 32 LOC, nesting depth < 3 levels.
- **Execution Performance**: Full test suite executes in $< 10$ seconds.

---

## 7. Execution Commands

```bash
# Compile and run the Urban Canyon E2E test suite
cargo test -p gneiss-rtk --test test_urban_canyon_e2e

# Run with full output logging
cargo test -p gneiss-rtk --test test_urban_canyon_e2e -- --nocapture

# Run clippy verification across the test target
cargo clippy -p gneiss-rtk --test test_urban_canyon_e2e --all-targets -- -D warnings
```
