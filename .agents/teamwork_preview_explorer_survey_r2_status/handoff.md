# Handoff Report — Frontier R2: Integer PPP-AR & eval_ppp Survey

**Agent**: Explorer R2 (`teamwork_preview_explorer_survey_r2_status`)  
**Date**: 2026-09-13  
**Handoff Type**: Hard (Investigation & Survey Complete)  
**Target Recipient**: Orchestrator / Sub-Orchestrator M2 (`a6307386-3f81-4920-9a31-a6d124a2f8d6`)  

---

## 1. Observation

1. **Unit Tests**:
   - `cargo test -p gneiss-parsers --lib sinex_bia`:
     ```
     running 2 tests
     test sinex_bia::tests::test_sinex_bias_parsing ... ok
     test sinex_bia::fallback_tests::test_sinex_bias_fallback ... ok
     test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 265 filtered out; finished in 0.00s
     ```
   - `cargo test -p gneiss-rtk --lib ambiguity::ppp_ar`:
     ```
     running 6 tests
     test ambiguity::ppp_ar::tests::test_between_satellite_single_difference_wide_lane_fixing ... ok
     test ambiguity::ppp_ar::tests::test_wide_lane_rounding_with_fractional_biases ... ok
     test ambiguity::ppp_ar::tests::test_sd_form_and_backsubstitute ... ok
     test ambiguity::ppp_ar::tests::test_resolve_sd_with_fixed_wl ... ok
     test ambiguity::ppp_ar::tests::test_fix_single_diff_ambiguities ... ok
     test ambiguity::ppp_ar::tests::test_ppp_mw_tracker_and_fixed_wl ... ok
     test result: ok. 6 passed; 0 failed; 0 ignored; 0 measured; 379 filtered out; finished in 0.00s
     ```
2. **Benchmark Execution (`PPP_ONLY=f9p cargo run --release --bin eval_ppp`)**:
   - Gneiss PPP vs CSRS-PPP ($N=600$):
     * Horizontal: $\text{p50} = 3.2210\text{ m}$, $\text{p68} = 3.2257\text{ m}$, $\text{p95} = 3.2379\text{ m}$, $\text{RMS} = 3.2216\text{ m}$, $\max = 3.2498\text{ m}$
     * 3D Position Error: $\text{p50} = 3.3594\text{ m}$, $\text{p68} = 3.3631\text{ m}$, $\text{p95} = 3.3723\text{ m}$, $\text{RMS} = 3.3582\text{ m}$, $\max = 3.3799\text{ m}$
     * East: $\text{mean} = +2.9951\text{ m}$, $\sigma = 0.0102\text{ m}$, $\text{RMS} = 2.9951\text{ m}$
     * North: $\text{mean} = +1.1866\text{ m}$, $\sigma = 0.0084\text{ m}$, $\text{RMS} = 1.1866\text{ m}$
     * Up: $\text{mean} = -0.9479\text{ m}$, $\sigma = 0.0211\text{ m}$, $\text{RMS} = 0.9481\text{ m}$
   - Gneiss PPP vs RTK Truth ($N=583$ matched):
     * Horizontal: $\text{p50} = 2.9265\text{ m}$, $\text{p68} = 2.9289\text{ m}$, $\text{p95} = 2.9317\text{ m}$, $\text{RMS} = 2.9260\text{ m}$, $\max = 2.9330\text{ m}$
     * 3D Position Error: $\text{p50} = 3.0470\text{ m}$, $\text{p68} = 3.0498\text{ m}$, $\text{p95} = 3.0547\text{ m}$, $\text{RMS} = 3.0478\text{ m}$, $\max = 3.0569\text{ m}$
     * East: $\text{mean} = +2.7168\text{ m}$, $\sigma = 0.0043\text{ m}$, $\text{RMS} = 2.7168\text{ m}$
     * North: $\text{mean} = +1.0866\text{ m}$, $\sigma = 0.0029\text{ m}$, $\text{RMS} = 1.0866\text{ m}$
     * Up: $\text{mean} = -0.8528\text{ m}$, $\sigma = 0.0053\text{ m}$, $\text{RMS} = 0.8529\text{ m}$
   - Commercial Baseline: CSRS-PPP vs RTK Truth ($N=583$ matched):
     * Horizontal: $\text{p50} = 0.2959\text{ m}$, $\text{p68} = 0.3011\text{ m}$, $\text{p95} = 0.3139\text{ m}$, $\text{RMS} = 0.2962\text{ m}$, $\max = 0.3285\text{ m}$
     * 3D Position Error: $\text{p50} = 0.3141\text{ m}$, $\text{p68} = 0.3194\text{ m}$, $\text{p95} = 0.3282\text{ m}$, $\text{RMS} = 0.3118\text{ m}$, $\max = 0.3366\text{ m}$
     * East: $\text{mean} = -0.2782\text{ m}$, $\sigma = 0.0116\text{ m}$, $\text{RMS} = 0.2785\text{ m}$
     * North: $\text{mean} = -0.1007\text{ m}$, $\sigma = 0.0081\text{ m}$, $\text{RMS} = 0.1011\text{ m}$
     * Up: $\text{mean} = +0.0949\text{ m}$, $\sigma = 0.0215\text{ m}$, $\text{RMS} = 0.0973\text{ m}$
   - Ambiguity Fix Rate: **0.0%** ($0$ out of $583$ epochs fixed).
3. **Configuration & Code Inspection**:
   - `crates/gneiss-rtk/src/bin/eval_ppp.rs:436`: `is_kinematic: false`.
   - `crates/gneiss-rtk/src/swfg/engine/mod.rs:187-190`: `if self.is_ppp && !self.is_kinematic { let pose_id = self.solver.ensure_static_pose(); ... }`. Trajectory solution `position_ecef` was identically `[-1276969.350, -4717194.952, 4087249.058]` for all 600 epochs.
   - `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs:20`: `tracker.update((constellation_id, satellite), mw_cycles, 0, slip);` — epoch is hardcoded to `0`.
   - `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs:174` & `ar_integration.rs:194`: `inject_fixed_priors` adds `PriorFactor` ($10^8$ hard constraint) directly to un-differenced ambiguity states, locking the receiver clock state to the float estimate.
   - `crates/gneiss-rtk/src/swfg/engine/mod.rs:382`: `SwfgSolution` always returns `error: None`. In `forward.rs:262`, `let is_fix = sol.error.is_some_and(|e| e < 0.05);` is always false, forcing `quality = 4` (Float).
   - `crates/gneiss-parsers/src/receiver_pcv/mod.rs:195-217`: Standalone `ReceiverPcv` methods `pco_correction_m` and `pcv_correction_m` exist, but are not invoked in `epoch.rs`.
4. **AGENTS.md Code Metrics**:
   - File size: `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs` is **506 LOC** (exceeds 500 limit).
   - Function sizes: `resolve_sd_with_fixed_wl` is 61 LOC; `execute_ar_step` is 74 LOC; `execute_ppp_ar_step` is 74 LOC; `apply_validated_fix` is 57 LOC (exceed 32 LOC limit).
   - Nesting depth: `ar_handler.rs` and `sinex_bia.rs` contain nesting levels of 4 and 5 (exceed 3 levels limit).
   - Production unwrap: Exactly 0 `unwrap()` in production code.
   - Warnings: 0 warnings under `cargo check --workspace`.

---

## 2. Logic Chain

1. From Observation 3 (`eval_ppp.rs:436` having `is_kinematic: false` and `mod.rs:187` creating a single static pose), the current evaluation run did not perform kinematic trajectory estimation; it fitted a single static 3D point to a moving vehicle.
2. From Observation 2, the vehicle trajectory differs from this static point by $2.926\text{ m}$ Horizontal RMS and $3.048\text{ m}$ 3D RMS vs RTK truth ($3.222\text{ m}$ Horizontal and $3.358\text{ m}$ 3D RMS vs CSRS-PPP).
3. From Observation 2, the acceptance criterion specifically demands: "Resolves integer ambiguities on the F9P kinematic drive, closing discrepancy vs CSRS-PPP (0.296 m RMS) to sub-meter kinematic accuracy."
4. Since $3.358\text{ m} \gg 1.0\text{ m}$ and $3.048\text{ m} \gg 1.0\text{ m}$, and ambiguity fix rate is $0.0\%$, the acceptance criterion is mathematically not met.
5. From Observation 3 (`tracker.update` passing epoch 0, `inject_fixed_priors` locking receiver clock, and `sol.error` being `None`), integer ambiguity resolution cannot succeed or report fixed quality in the factor graph.
6. From Observation 4, `ppp_ar.rs` violates the 500 LOC ceiling (506 LOC), and several core functions violate the 32 LOC ceiling and nesting depth limit.
7. Therefore, achieving the acceptance criteria requires setting `is_kinematic: true`, passing the true epoch index to the MW tracker, constraining ambiguities via relative between-satellite single-difference factors (preserving receiver clock freedom), projecting receiver antenna PCO/PCV, propagating fix quality into `FilteredEpoch`, and modularizing `ppp_ar.rs` into `mw_tracker.rs`.

---

## 3. Caveats

1. The current run was strictly read-only per Explorer archetype constraints; no production or test source files were modified during this investigation.
2. The benchmark run evaluated 600 epochs of the u-blox ZED-F9P dataset (`max_epochs: 600`), matching the standard evaluation harness limit.
3. The CODE MGEX products used were `com21374` from `datasets/wtzr_ppp_1224/` which corresponds to the same GPS day (week 2137 day 4) as `f9p_ppp_1224/rover.obs`.

---

## 4. Conclusion

The Frontier R2 PPP-AR implementation has operational unit test suites and foundational components (SINEX BIA indexed parser, ANTEX nadir PCV, SD Wide-Lane/Narrow-Lane math), but **does NOT meet the sub-meter kinematic accuracy target** or the **AGENTS.md file/function size limits**.

Exact numerical status:
- **Ambiguity Fix Rate**: $0.0\%$ ($0 / 583$ epochs)
- **3D Error vs RTK Truth**: $\text{p50} = 3.047\text{ m}$, $\text{RMS} = 3.048\text{ m}$
- **Discrepancy vs CSRS-PPP**: $\text{p50} = 3.359\text{ m}$, $\text{RMS} = 3.358\text{ m}$
- **Target**: Sub-meter kinematic accuracy ($< 1.0\text{ m}$ RMS) and integer ambiguity resolution.

The exact blueprint of required code modifications is documented in `/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r2_status/report.md`.

---

## 5. Verification Method

To independently verify these findings:
1. Run parser and AR unit tests:
   ```bash
   cargo test -p gneiss-parsers --lib sinex_bia
   cargo test -p gneiss-rtk --lib ambiguity::ppp_ar
   ```
2. Run the evaluation benchmark:
   ```bash
   PPP_ONLY=f9p cargo run --release --bin eval_ppp
   ```
3. Run the ENU analysis script:
   ```bash
   python3 .agents/teamwork_preview_explorer_survey_r2_status/compute_enu.py
   ```
4. Verify AGENTS.md code metrics:
   ```bash
   wc -l crates/gneiss-rtk/src/ambiguity/ppp_ar.rs
   ```
   (Outputs `506`, confirming violation of $< 500$ LOC limit).
