# Project: Urban Canyon Fix Rate Expansion and Multipath Mitigation

## Architecture
This project enhances the Gneiss RTK/PPK double-difference engine and tightly coupled INS estimator to dramatically expand ambiguity fix rates and suppress tail errors in severe urban canyons (Tokyo Shinjuku and Hong Kong Whampoa) while maintaining zero false integer fixes and passing all CI regression smoke guards.

### Core Architecture Components:
1. **Physical Observation Covariance**:
   - `crates/gneiss-core/src/obs.rs`: Fractional C/N0 retrieval (`get_snr_f64`).
   - `crates/gneiss-core/src/variance.rs`: Physically grounded $C^1$-smooth SIGMA-SNR elevation and noise model:
     $$\sigma^2(\theta, S) = \left( a^2 + \frac{b^2}{\sin^2 \theta + \sin^2 \theta_0} \right) \cdot f_{\text{SNR}}(S)$$
     with smooth logistic activation ($\tau = 1.5\text{ dB-Hz}$) and algebraic saturation ceiling ($f_{\max} = 1000.0$).
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`: Smooth double-difference covariance formulation with exact rover and base satellite noise tracking ($R_{DD} \succ 0$).

2. **Code-Minus-Carrier (CMC) Multipath Mitigation**:
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`: Arc-level dual-frequency ($MP$) and single-frequency ($CMC = P - \Phi - 2I$) tracking.
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/screen.rs`: Decoupled code blunder suppression preventing the deletion of valid carrier phase measurements.
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`: Shielding the Melbourne-Wübbena running filter from false cycle slip resets during code multipath steps.

3. **Doppler-Assisted Cycle Slip Validation**:
   - `crates/gneiss-rtk/src/post_process/screening.rs`: Multi-band Doppler phase continuity validation across GPS, Galileo (E1, E5a, E5b), and BeiDou (B1, B2a, B3I) with tuned adaptive thresholds ($0.28\text{--}0.35$ cyc) to detect half-cycle and 1-cycle slips.
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`: Unconditional base-side slip monitoring, propagation of detected slips to wide-lane trackers (`mw.rs`, `pw_tracker`), and tracking arc epoch resets.

4. **C/N0- and Elevation-Prioritized Partial Ambiguity Resolution (PAR)**:
   - `crates/gneiss-rtk/src/ambiguity/ar_subsets.rs` & `crates/gneiss-rtk/src/ambiguity/par.rs`: Composite Quality Metric (CQM) ranking candidate ambiguities by elevation, SNR, continuous lock duration, and CMC multipath variance.
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs`: Geometric DOP guard (minimum 4 satellites) and strict covariance definiteness enforcement ($P \succ 0, Q_{aa} \succeq 0$) eliminating off-diagonal covariance leakage upon integer conditioning.

---

## Feature Inventory
| # | Feature | Description | Milestone | Source |
|---|---------|-------------|-----------|--------|
| 1 | Fractional C/N0 API | Add `get_snr_f64` in `obs.rs` to eliminate 1 dB integer quantization | M1 | Survey R1 |
| 2 | Unified SIGMA-SNR Model | Implement $C^1$-smooth monotonic variance function in `variance.rs` | M1 | Survey R1 |
| 3 | Smooth DD Covariance | Implement continuous DD code/phase covariance scaling in `formation_cov.rs` | M1 | Survey R1 |
| 4 | Decoupled PR Gross Error Screen | Modify `screen.rs` to zero code weight without deleting carrier phase | M2 | Survey R2 |
| 5 | CMC Arc Multipath Tracking | Track $MP$ / $CMC$ residuals per satellite arc in `formation.rs` | M2 | Survey R2 |
| 6 | Adaptive Code De-Weighting | Inflate pseudorange variance on code multipath without affecting phase | M2 | Survey R2 |
| 7 | MW Multipath Shielding | Prevent code multipath jumps from tripping false slip resets in `mw.rs` | M2 | Survey R2 |
| 8 | Multi-band Doppler Slip Detection | Extend Doppler slip detector to bands 1, 2, 5, 6, 7 with 0.30 cyc threshold | M3 | Survey R3 |
| 9 | Unconditional Base Slip & Leak Plug | Check base slips unconditionally and propagate flags to `mw.rs` & `pw_tracker` | M3 | Survey R3 |
| 10 | CQM Candidate Ambiguity Ranking | Order PAR candidates by elevation, C/N0, lock time, and CMC variance | M4 | Survey R3 |
| 11 | PAR DOP Guard & Subset Constraints | Add DOP threshold and enforce minimum subset size $\ge 4$ in `ar_subsets.rs` | M4 | Survey R3 |
| 12 | Positive-Definite Ambiguity Conditioning | Zero cross-covariances and enforce $\lambda_{\min}(P) \ge 10^{-6}$ in `ar.rs` | M4 | Survey R3 |
| 13 | E2E Benchmark & UrbanNav Validation | Pass 100% E2E test suite, expand fix rates in Shinjuku/Whampoa, pass smoke guards | M5 | Prompt |

---

## Milestones
| # | Name | Scope | Dependencies | Status |
|---|------|-------|-------------|--------|
| M1 | Adaptive SNR & Elevation Observation Covariance | F1, F2, F3 (`obs.rs`, `variance.rs`, `formation_cov.rs`) | none | DONE |
| M2 | Code-Minus-Carrier Multipath Screening & Mitigation | F4, F5, F6, F7 (`screen.rs`, `formation.rs`, `mw.rs`) | M1 | DONE |
| M3 | Doppler-Assisted Cycle Slip Detection & Continuity | F8, F9 (`screening.rs`, `formation.rs`, `mw.rs`) | none (independent) | DONE |
| M4 | C/N0- and Elevation-Prioritized PAR Engine | F10, F11, F12 (`par.rs`, `ar_subsets.rs`, `ar.rs`) | M1, M2 | DONE |
| M5 | Final E2E Test Pass & Urban Canyon Benchmark Hardening | F13 (Full E2E test suite, UrbanNav benchmarks, CI smoke guards) | M1, M2, M3, M4, TEST_READY | DONE |

---

## Interface Contracts

### 1. `obs.rs` ↔ `variance.rs` & `formation_cov.rs`
- `SatObs::get_snr_f64(&self, freq_band: u8) -> Option<f64>`: Returns raw floating-point C/N0 in dB-Hz without integer quantization.
- Existing `get_snr(&self, freq_band: u8) -> Option<u8>` is preserved for backward compatibility.

### 2. `variance.rs` ↔ `formation_cov.rs`
- `pub fn sigma_snr_variance(snr_dbhz: f64, el_rad: f64, is_phase: bool) -> f64`
- Inputs: `snr_dbhz: f64`, `el_rad: f64`, `is_phase: bool`.
- Output: strictly positive, finite variance $\sigma^2 > 0$ with continuous first derivatives $\frac{\partial \sigma^2}{\partial \theta} \le 0$ and $\frac{\partial \sigma^2}{\partial S} \le 0$.

### 3. `formation.rs` ↔ `screen.rs` & `formation_cov.rs`
- When code gross error is detected on satellite pair `k`, `screen_gross_pr_errors` marks the code observation for de-weighting / suppression (`pr_suppressed = true` or `pr_var_m2 = 1e6`) rather than dropping the `DoubleDiffMeasurement`.
- The corresponding carrier phase `dd_cp_cycles` remains active with nominal variance.

### 4. `screening.rs` ↔ `formation.rs` & `mw.rs`
- `CycleSlipDetector::check_epoch` checks bands `[1, 2, 5, 6, 7]` with adaptive threshold `(0.30 * dt).clamp(0.28, 1.0)`.
- Returns set of slipped `SatelliteId`s.
- `formation.rs` marks rover and base slips into `DoubleDiffMeasurement::slip`, feeds `slip: true` into `WidelaneTracker::update_tracker_from_obs`, and resets `pair_epochs` to 0.

### 5. `ar_subsets.rs` ↔ `ar.rs`
- `select_par_candidates` computes Composite Quality Metric (CQM):
  $$CQM_i = w_{el} \sin\theta_i + w_{snr} \frac{S_i - 20}{30} + w_{lock} \min\left(1, \frac{t_{lock}}{30}\right) - w_{cmc} \frac{\sigma_{cmc}}{2.0} - w_{var} \sigma_{a_i}$$
  and sorts candidate ambiguities in descending order of quality.
- `condition_state_on_integers`: For all fixed indices $k$, zeros cross-covariances $P_{ki} = 0, P_{ik} = 0$ for $i \ne k$, sets $P_{kk} = 10^{-4}$, and verifies $\lambda_{\min}(P) \ge 10^{-6}$.

---

## Code Layout & Write Boundaries
- `crates/gneiss-core/src/obs.rs`: M1 Worker
- `crates/gneiss-core/src/variance.rs`: M1 Worker
- `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`: M1 Worker
- `crates/gneiss-rtk/src/estimators/rtk_iekf/screen.rs`: M2 Worker
- `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`: M2 Worker / M3 Worker (sequential dependency)
- `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`: M2 Worker / M3 Worker
- `crates/gneiss-rtk/src/post_process/screening.rs`: M3 Worker
- `crates/gneiss-rtk/src/ambiguity/par.rs`: M4 Worker
- `crates/gneiss-rtk/src/ambiguity/ar_subsets.rs`: M4 Worker
- `crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs`: M4 Worker
- `crates/gneiss-rtk/tests/`: E2E Test Writer (exclusive ownership of dedicated E2E test files)

## Verification Invariants
- All code strictly complies with `AGENTS.md`: LOC < 500, func < 32, nest < 3, 0 unwrap in prod.
- 0 compiler warnings, 0 clippy warnings (`cargo clippy --workspace --all-targets -- -D warnings`).
- `cargo test --workspace` passes 100%.
- Both smoke scripts pass: `python3 scripts/check_network_benchmark.py --smoke`, `python3 scripts/check_multignss_benchmark.py --smoke`.
