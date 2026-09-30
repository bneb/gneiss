# Handoff Report — Milestone 2: Code-Minus-Carrier (CMC) Multipath Detection & Down-Weighting

**Agent**: Worker M2 (`worker_m2_urban_canyon`)  
**Parent**: Orchestrator Urban Canyon Gen 2 (`5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0`)  
**Date**: 2026-09-24  
**Type**: Hard Handoff (Task Complete)

---

## 1. Observation

1. **Pre-existing gross pseudorange screen dropped entire DD measurements**:
   In `crates/gneiss-rtk/src/estimators/rtk_iekf/screen.rs:46`, whenever the prefit DD pseudorange residual exceeded `GROSS_PR_ERROR_THRESHOLD_M` (20.0 m), the measurement was unconditionally dropped via `measurements.remove(idx)` even when valid, millimeter-level carrier phase (`dd_cp_cycles.is_some()`) was present.
2. **E2E test suite expectation in `tier1_features.rs:165`**:
   `test_r2_screen_gross_error_preserves_carrier_phase` previously asserted `assert_eq!(meas.len(), 1)` based on blunder deletion. Updating `screen.rs` to retain measurements with active carrier phase caused `meas.len() == 2`. Under orchestrator authorization, line 165 was updated to:
   ```rust
   assert_eq!(meas.len(), 2, "Both measurements retained: carrier phase preserved");
   let screened = meas.iter().find(|m| m.key.sat == 2).expect("corrupted pair retained");
   assert!(screened.dd_cp_cycles.is_some(), "Carrier phase must remain active on screened pair");
   assert!(screened.pr_var_m2 >= 1.0e6, "Code variance must be inflated/suppressed");
   assert_eq!(screened.cp_var_cycles2, 0.0005, "Carrier variance must remain nominal");
   ```
3. **Target file LOC constraints (`AGENTS.md`)**:
   `crates/gneiss-rtk/src/estimators/rtk_iekf/mod.rs` was at 482 LOC (near the 500 LOC ceiling). Adding `cmc_tracker` to `GnssRtkIekf` directly in `mod.rs` would have breached 500 LOC. Instead, `CmcTracker` was integrated inside `WidelaneTracker` in `mw.rs`.
   Post-implementation line counts for all owned/modified files:
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/screen.rs`: 187 LOC (< 500)
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`: 452 LOC (< 500)
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/update/robust.rs`: 389 LOC (< 500)
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`: 485 LOC (< 500)
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`: 235 LOC (< 500)
4. **Code standard verification checks**:
   - Zero `unwrap()` calls in production code across all modified files.
   - All functions in our implementation strictly <= 32 LOC with nesting depth < 3 levels.
   - `cargo clippy --workspace --all-targets -- -D warnings`: passed with 0 warnings.
   - `cargo test -p gneiss-rtk --lib`: passed with 441 passed, 0 failed.
   - `cargo test --test test_urban_canyon_e2e`: passed with 51 passed, 0 failed.
   - `python3 scripts/check_network_benchmark.py --smoke`: ALL CHECKS PASSED (OHLN fix rate 99.1% >= 79.0%, P181 99.7% >= 85.0%, SLAC 71.1% >= 60.0%).
   - `python3 scripts/check_multignss_benchmark.py --smoke`: ALL CHECKS PASSED (P181 fix rate 97.7% >= 97.5%, P225 92.6% >= 71.0%, P222 99.4% >= 86.0%, network fused 99.8% >= 96.5%).

---

## 2. Logic Chain

1. **Decoupling Pseudorange Screening (F4)**:
   - *Premise*: Severe code multipath or diffraction blunders (> 20.0 m) frequently contaminate code pseudorange in urban canyons while the direct Line-of-Sight carrier phase remains intact.
   - *Mechanism*: In `screen.rs`, `screen_gross_pr_errors` checks `measurements[idx].dd_cp_cycles.is_some()`. If true, the pair is retained in `measurements`, its `pr_var_m2` is inflated to `GROSS_PR_DEWEIGHT_VAR_M2` ($1.0 \times 10^8\text{ m}^2$), its carrier phase variance is untouched, and its key is recorded in `rejected` to prevent redundant screening. If `dd_cp_cycles.is_none()`, the pair is dropped with `measurements.remove(idx)`.
   - *Result*: Kalman update effectively zeroes out the gain on the corrupted pseudorange ($K \propto H^T R^{-1} \to 0$) while retaining millimeter-level phase constraints.

2. **Code-Minus-Carrier (CMC) Tracking and Baseline Freezing (F5)**:
   - *Premise*: Geometry-free $CMC = \rho_{DD} - \lambda \cdot \Phi_{DD} = 2 I_{DD} - \lambda N_{DD} + M_{P} + \epsilon$. Over short spans without cycle slips, ionosphere and integer ambiguities vary smoothly, making deviations $|\Delta CMC| > 2.5\text{ m}$ a direct proxy for code multipath $M_P$.
   - *Mechanism*: `CmcTrack` in `update/robust.rs` maintains a baseline estimate across continuous tracking arcs. During the initial 5 epochs (`CMC_WARMUP_EPOCHS`), it accumulates an unweighted running mean. Once warmed up, if $|CMC - \text{baseline}| > 2.5\text{ m}$, the baseline is frozen to avoid contaminating the filter's memory, and the deviation is stored as `multipath_m`. If a cycle slip occurs (`slip == true`), the arc resets immediately.
   - *Result*: Robust, outlier-resistant detection of sudden multipath jumps.

3. **Adaptive Pseudorange Down-Weighting (F6)**:
   - *Premise*: Code multipath should increase code variance $R_{PP} \leftarrow R_{PP} + \sigma_{\text{mp}}^2$ without degrading carrier phase variance $R_{\Phi\Phi}$.
   - *Mechanism*: In `formation.rs`, `apply_cmc_deweight` feeds observations to `CmcTracker`. If multipath is detected, `apply_cmc_downweighting` computes $R_{PP} + \sigma_{\text{mp}}^2$, updating `meas.pr_var_m2` while keeping `meas.cp_var_cycles2` intact.
   - *Result*: The filter automatically de-weights contaminated code measurements proportionally to multipath severity.

4. **Melbourne-Wübbena Multipath Shielding (F7)**:
   - *Premise*: Sudden code multipath jumps can cause wide-lane innovation jumps $> 1.0\text{ cycle}$ in Melbourne-Wübbena combinations. If unshielded, the wide-lane tracker misinterprets the code jump as a cycle slip, resetting the running average and discarding accumulated fix history.
   - *Mechanism*: In `mw.rs`, `WidelaneTracker` tracks CMC state alongside wide-lane arcs. When `tracker.is_cmc_multipath(&key)` is true along a continuous carrier arc, `MwTrack::absorb(&mut self, x: f64, shielded: bool)` shields the track: large innovations ($> 1.0\text{ cycle}$) do not reset the track and do not corrupt the running mean.
   - *Result*: Prevents false cycle-slip resets and protects ambiguity fix continuity.

---

## 3. Caveats

1. **Cycle Slip Detection Dependency**:
   CMC arc tracking relies on cycle slip flags (`lli_slip`). If an unflagged carrier cycle slip coincides with a code multipath jump, CMC tracking will treat the step as code multipath until Doppler cycle-slip screening (Milestone 3 / R3) or LLI detects the slip.
2. **Satellite Geometry Reset**:
   When a satellite sets or is dropped from active tracking, `retain_active(&ddm.active_keys)` purges inactive CMC tracks to prevent memory leakage across long runs.
3. **No Caveats Beyond Above**:
   All requirements from the dispatch, `PROJECT.md`, and survey reports are fully addressed and verified.

---

## 4. Conclusion

Milestone 2 (R2) is complete, robust, and verified:
- **F4**: Prefit PR screening decoupled from carrier phase in `screen.rs`.
- **F5 & F6**: CMC arc tracking and adaptive code down-weighting implemented in `formation.rs` and `robust.rs`.
- **F7**: Melbourne-Wübbena multipath shielding implemented in `mw.rs`.
- Zero compiler/clippy warnings under `-D warnings`.
- Zero `unwrap()` calls in production code.
- All files strictly < 500 LOC, functions <= 32 LOC, nesting < 3 levels.
- 441 lib tests passed, 51 e2e tests passed, both CI smoke benchmark guards passed.

---

## 5. Verification Method

To independently reproduce and verify this work, run:

```bash
# 1. Clippy lint check (must produce 0 warnings)
cargo clippy --workspace --all-targets -- -D warnings

# 2. Gneiss RTK unit tests (441 tests must pass)
cargo test -p gneiss-rtk --lib

# 3. Urban Canyon E2E and tier1 integration tests (51 tests must pass)
cargo test --test test_urban_canyon_e2e

# 4. CI Benchmark smoke guards (both must report ALL CHECKS PASSED)
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke

# 5. Production unwrap and LOC audit
wc -l crates/gneiss-rtk/src/estimators/rtk_iekf/screen.rs \
      crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs \
      crates/gneiss-rtk/src/estimators/rtk_iekf/update/robust.rs \
      crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs \
      crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs
```
