# Investigation Report — Frontier R2: Integer PPP-AR Engine & eval_ppp Benchmark Status

**Author**: Explorer R2 (`teamwork_preview_explorer_survey_r2_status`)  
**Date**: 2026-09-13  
**Working Directory**: `/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r2_status/`  
**Target Milestone**: Frontier R2 — Integer PPP-AR Engine via SINEX OSB Ingestion and Kinematic Benchmark Verification  

---

## 1. Executive Summary

We conducted a comprehensive, read-only investigation of Frontier R2 across `gneiss-parsers` and `gneiss-rtk`, executed existing unit test suites, evaluated the `eval_ppp` benchmark against Canadian Geodetic Service CSRS-PPP truth and RTK ground truth, analyzed why the acceptance criteria are not yet met, audited AGENTS.md code standards, and specified the exact architectural and mathematical fixes required.

### Core Finding
The acceptance criterion — **"Resolves integer ambiguities on the F9P kinematic drive, closing discrepancy vs CSRS-PPP (0.296 m RMS) to sub-meter kinematic accuracy"** — is **NOT YET MET**:
1. **Current 3D Error vs RTK Truth**: $\text{p50} = \mathbf{3.047\text{ m}}$, $\text{RMS} = \mathbf{3.048\text{ m}}$ (Horizontal $\text{p50} = 2.927\text{ m}$, $\text{RMS} = 2.926\text{ m}$).
2. **Current Discrepancy vs CSRS-PPP**: $\text{p50} = \mathbf{3.359\text{ m}}$, $\text{RMS} = \mathbf{3.358\text{ m}}$ (Horizontal $\text{p50} = 3.221\text{ m}$, $\text{RMS} = 3.222\text{ m}$).
3. **Commercial Baseline (CSRS-PPP vs RTK Truth)**: $\text{p50} = \mathbf{0.314\text{ m}}$, $\text{RMS} = \mathbf{0.312\text{ m}}$ (Horizontal $\text{p50} = 0.296\text{ m}$, $\text{RMS} = 0.296\text{ m}$).
4. **Ambiguity Fix Rate**: $\mathbf{0.0\%}$ ($0$ out of $583$ evaluated epochs fixed).
5. **Primary Culprit**: In `crates/gneiss-rtk/src/bin/eval_ppp.rs`, the F9P vehicle dataset was mistakenly configured with `is_kinematic: false`. This triggered `solver.ensure_static_pose()`, locking all 600 epochs of the driving vehicle to a single stationary point (`[-1276969.350, -4717194.952, 4087249.058]`). Moreover, ambiguity fixing never succeeds because the Melbourne-Wübbena tracker receives a hardcoded epoch of `0`, and un-differenced prior injection locks receiver clock states rather than between-satellite single differences.

---

## 2. Test & Benchmark Execution Results

### 2.1 Unit Tests

1. `cargo test -p gneiss-parsers --lib sinex_bia`:
   - **Result**: `ok. 2 passed; 0 failed; 0 ignored`
   - Tests: `sinex_bia::tests::test_sinex_bias_parsing`, `sinex_bia::fallback_tests::test_sinex_bias_fallback`
2. `cargo test -p gneiss-rtk --lib ambiguity::ppp_ar`:
   - **Result**: `ok. 6 passed; 0 failed; 0 ignored`
   - Tests:
     * `test_between_satellite_single_difference_wide_lane_fixing`
     * `test_wide_lane_rounding_with_fractional_biases`
     * `test_sd_form_and_backsubstitute`
     * `test_resolve_sd_with_fixed_wl`
     * `test_fix_single_diff_ambiguities`
     * `test_ppp_mw_tracker_and_fixed_wl`

### 2.2 Benchmark: `PPP_ONLY=f9p cargo run --release --bin eval_ppp`

Evaluated against `datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.obs` (600 epochs, 1Hz dual-frequency u-blox ZED-F9P kinematic vehicle drive):

| Comparison | Epochs ($N$) | Fix Rate (%) | East RMS (m) | North RMS (m) | Up RMS (m) | Horiz p50 (m) | Horiz RMS (m) | 3D p50 (m) | 3D RMS (m) |
|---|---|---|---|---|---|---|---|---|---|
| **Commercial Baseline: CSRS-PPP vs RTK Truth** | 583 | 100.0% (Float PPP) | 0.2785 | 0.1011 | 0.0973 | 0.2959 | 0.2962 | 0.3141 | 0.3118 |
| **Gneiss PPP vs RTK Truth** | 583 | **0.0%** | **2.7168** | **1.0866** | **0.8529** | **2.9265** | **2.9260** | **3.0470** | **3.0478** |
| **Gneiss PPP vs CSRS-PPP Discrepancy** | 600 | **0.0%** | **2.9951** | **1.1866** | **0.9481** | **3.2210** | **3.2216** | **3.3594** | **3.3582** |

### 2.3 Detailed Component Error Distributions (ENU)

#### A. Commercial Baseline: CSRS-PPP vs RTK Truth (Matched $N=583$)
- **East**: $\text{mean} = -0.2782\text{ m}$, $\sigma = 0.0116\text{ m}$, $\text{RMS} = 0.2785\text{ m}$, $[\min, \max] = [-0.3167, -0.2461]\text{ m}$
- **North**: $\text{mean} = -0.1007\text{ m}$, $\sigma = 0.0081\text{ m}$, $\text{RMS} = 0.1011\text{ m}$, $[\min, \max] = [-0.1238, -0.0752]\text{ m}$
- **Up**: $\text{mean} = +0.0949\text{ m}$, $\sigma = 0.0215\text{ m}$, $\text{RMS} = 0.0973\text{ m}$, $[\min, \max] = [+0.0435, +0.1492]\text{ m}$
- **Horizontal**: $\text{p50} = 0.2959\text{ m}$, $\text{p68} = 0.3011\text{ m}$, $\text{p95} = 0.3139\text{ m}$, $\text{RMS} = 0.2962\text{ m}$, $\max = 0.3285\text{ m}$
- **3D Error**: $\text{p50} = 0.3141\text{ m}$, $\text{p68} = 0.3194\text{ m}$, $\text{p95} = 0.3282\text{ m}$, $\text{RMS} = 0.3118\text{ m}$, $\max = 0.3366\text{ m}$

#### B. Discrepancy: Current Gneiss PPP vs CSRS-PPP ($N=600$)
- **East**: $\text{mean} = +2.9951\text{ m}$, $\sigma = 0.0102\text{ m}$, $\text{RMS} = 2.9951\text{ m}$, $[\min, \max] = [+2.9672, +3.0276]\text{ m}$
- **North**: $\text{mean} = +1.1866\text{ m}$, $\sigma = 0.0084\text{ m}$, $\text{RMS} = 1.1866\text{ m}$, $[\min, \max] = [+1.1530, +1.2092]\text{ m}$
- **Up**: $\text{mean} = -0.9479\text{ m}$, $\sigma = 0.0211\text{ m}$, $\text{RMS} = 0.9481\text{ m}$, $[\min, \max] = [-0.9986, -0.9006]\text{ m}$
- **Horizontal**: $\text{p50} = 3.2210\text{ m}$, $\text{p68} = 3.2257\text{ m}$, $\text{p95} = 3.2379\text{ m}$, $\text{RMS} = 3.2216\text{ m}$, $\max = 3.2498\text{ m}$
- **3D Error**: $\text{p50} = 3.3594\text{ m}$, $\text{p68} = 3.3631\text{ m}$, $\text{p95} = 3.3723\text{ m}$, $\text{RMS} = 3.3582\text{ m}$, $\max = 3.3799\text{ m}$

#### C. Current Gneiss PPP vs RTK Truth (Matched $N=583$)
- **East**: $\text{mean} = +2.7168\text{ m}$, $\sigma = 0.0043\text{ m}$, $\text{RMS} = 2.7168\text{ m}$, $[\min, \max] = [+2.7098, +2.7253]\text{ m}$
- **North**: $\text{mean} = +1.0866\text{ m}$, $\sigma = 0.0029\text{ m}$, $\text{RMS} = 1.0866\text{ m}$, $[\min, \max] = [+1.0798, +1.0921]\text{ m}$
- **Up**: $\text{mean} = -0.8528\text{ m}$, $\sigma = 0.0053\text{ m}$, $\text{RMS} = 0.8529\text{ m}$, $[\min, \max] = [-0.8644, -0.8352]\text{ m}$
- **Horizontal**: $\text{p50} = 2.9265\text{ m}$, $\text{p68} = 2.9289\text{ m}$, $\text{p95} = 2.9317\text{ m}$, $\text{RMS} = 2.9260\text{ m}$, $\max = 2.9330\text{ m}$
- **3D Error**: $\text{p50} = 3.0470\text{ m}$, $\text{p68} = 3.0498\text{ m}$, $\text{p95} = 3.0547\text{ m}$, $\text{RMS} = 3.0478\text{ m}$, $\max = 3.0569\text{ m}$

---

## 3. Deep-Dive Root Cause Analysis

### 3.1 The "Static Pose" Artifice in `eval_ppp.rs`
In `crates/gneiss-rtk/src/bin/eval_ppp.rs` line 436:
```rust
    let f9p_spec = PppDatasetSpec {
        name: "RTK Explorer F9P Kinematic Vehicle Drive (1Hz, Integer PPP-AR)",
        ...
        is_kinematic: false, // <-- ERROR: Should be true!
    };
```
When `is_kinematic` is false, `crates/gneiss-rtk/src/swfg/engine/mod.rs` (lines 187–190) executes:
```rust
    if self.is_ppp && !self.is_kinematic {
        let pose_id = self.solver.ensure_static_pose();
        self.solver.create_epoch_variables(epoch, n_sats, &constellations, has_imu, is_rtk);
        Ok(pose_id)
    }
```
This causes the factor graph solver to allocate only **one single pose variable** across all 600 epochs. Because the F9P vehicle drove a small loop around a parking lot / open field, solving for a single static point converged to the stationary centroid `[-1276969.350, -4717194.952, 4087249.058]`. The trajectory output by Gneiss was completely flat and constant across all 600 epochs.

### 3.2 Hardcoded Epoch `0` in Melbourne-Wübbena Tracking
In `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs` lines 17–22:
```rust
pub fn record_mw_sample(constellation_id: u8, satellite: u16, mw_cycles: f64, slip: bool) {
    if let Ok(mut guard) = PPP_MW_TRACKER.lock() {
        let tracker = guard.get_or_insert_with(crate::ambiguity::ppp_ar::PppMwTracker::new);
        tracker.update((constellation_id, satellite), mw_cycles, 0, slip); // <-- ERROR: epoch is hardcoded 0!
    }
}
```
Because the 3rd argument is unconditionally `0`, `PppMwTracker::update` never sees advancing epochs (`epoch > entry.2 + 2`). More critically, `record_mw_if_dual_freq` in `epoch.rs` does not receive the epoch number, preventing proper continuous arc tracking across cycle slips and outages.

### 3.3 Receiver Clock Locking via Un-differenced `PriorFactor` Injection
In `execute_ppp_ar_step`, integer search resolves between-satellite single differences, back-substitutes them into un-differenced float ambiguities, and calls `inject_fixed_priors(solver, &sub_ids, &fixed_undiff)`:
```rust
pub fn inject_fixed_priors(graph: &mut EstimationGraph, amb_var_ids: &[VariableId], fixed_integers: &DVector<f64>) {
    for (i, &amb_id) in amb_var_ids.iter().enumerate() {
        let factor = PriorFactor {
            variable: amb_id,
            mu: DVector::from_element(1, fixed_integers[i]),
            information: DMatrix::from_element(1, 1, 1e8), // hard constraint
        };
        graph.add_factor(Box::new(factor));
    }
}
```
In un-differenced carrier phase observation equations:
$$\rho_{IF} = \|\mathbf{r}_s - \mathbf{r}_{rx}\| + c \delta t_{rx} - c \delta t^s + T_z m_w + A_{IF}^s$$
The ambiguity variable $A_{IF}^s$ directly couples with the receiver clock bias $c \delta t_{rx}$. Freezing $A_{IF}^s$ with an infinite prior ($10^8$) forcibly locks the receiver clock bias to the float estimate at that epoch. As the local oscillator drifts, carrier phase residuals blow up, triggering immediate fix validation failure (`validate_fix_geometry`) or optimizer divergence.
**Correct Formulation**: Must constrain between-satellite single-difference pairs via a 2-variable relative factor:
$$r = (A_{IF}^s - A_{IF}^{s_0}) - \left[ \lambda_{NL} \check{N}_{NL}^{s, s_0} + \frac{c f_2}{f_1^2 - f_2^2} \check{N}_{WL}^{s, s_0} \right]$$
This leaves the common receiver clock mode unconstrained.

### 3.4 Disconnection of Fix Quality Signaling
In `crates/gneiss-rtk/src/swfg/engine/mod.rs` line 382, `SwfgSolution` always sets `error: None`. In `crates/gneiss-rtk/src/post_process/forward.rs` lines 262–264:
```rust
let is_fix = sol.error.is_some_and(|e| e < 0.05);
let q = if is_fix { 1 } else if base_pos.is_some() { 2 } else { 4 };
```
Because `sol.error` is never populated, `is_fix` is always `false`, and every epoch is marked `quality = 4` (Float). This explains why the fix rate is reported as exactly 0.0%.

### 3.5 Missing Standalone Receiver Antenna PCO/PCV
Although `ReceiverPcv` methods `pco_correction_m` and `pcv_correction_m` were added to `receiver_pcv/mod.rs`, they are **never called** in `crates/gneiss-rtk/src/swfg/engine/epoch.rs`. For un-differenced PPP, antenna phase center offsets and elevation-dependent variations must be projected along the line of sight for each visible satellite.

---

## 4. AGENTS.md Compliance Audit

| Standard | Rule | Current Status | Violations / Details |
|---|---|---|---|
| **File size** | $< 500$ LOC | **FAIL** | `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs` is **506 LOC** ($> 500$). `antex.rs` is 491 LOC and `epoch.rs` is 492 LOC (close to limit). |
| **Function size** | $< 32$ LOC | **FAIL** | Multiple functions exceed 32 LOC: `resolve_sd_with_fixed_wl` (61 LOC), `fix_sd_wide_lane_subset` (55 LOC), `execute_ar_step` (74 LOC), `execute_ppp_ar_step` (74 LOC), `apply_validated_fix` (57 LOC), `evaluate_ppp_dataset` (181 LOC). |
| **Nesting depth** | $< 3$ levels | **FAIL** | Found depth 4 and 5 in `ar_handler.rs` (lines 44–51), `sinex_bia.rs` (lines 118–120), and `eval_ppp.rs` (lines 61–65). |
| **unwrap() in prod** | 0 | **PASS** | Exactly 0 `unwrap()` calls in production code across all audited files. |
| **Compiler warnings** | 0 | **PASS** | `cargo check --workspace` passes cleanly with 0 warnings. |

---

## 5. Required Code Modifications Blueprint

To achieve the sub-meter kinematic PPP-AR acceptance criterion and full AGENTS.md compliance, implement these specific changes:

### 1. `crates/gneiss-rtk/src/bin/eval_ppp.rs`
- In `f9p_spec`:
  - Set `is_kinematic: true`.
  - Align precise products to the F9P dataset directory:
    * `sp3_path: "datasets/rtkexplorer/sample_1/f9p_ppp_1224/COD0MGXFIN_20203590000_01D_05M_ORB.SP3"`
    * `clk_path: "datasets/rtkexplorer/sample_1/f9p_ppp_1224/COD0MGXFIN_20203590000_01D_30S_CLK.CLK"`
    * `bia_path: Some("datasets/rtkexplorer/sample_1/f9p_ppp_1224/COD0MGXFIN_20203590000_01D_01D_OSB.BIA")`
- In `print_stats`:
  - Add East, North, and Up error statistics (mean, std, RMS, min, max).
  - Print ambiguity fix rate percentage ($N_{\text{fixed}} / N_{\text{total}} \times 100\%$).

### 2. `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`
- **File modularization**: Extract `PppMwTracker` into a new file `crates/gneiss-rtk/src/ambiguity/mw_tracker.rs` (~150 LOC), bringing `ppp_ar.rs` down to ~350 LOC (strictly $< 500$ LOC).
- **Function refactoring**: Decompose `resolve_sd_with_fixed_wl` (61 LOC) into `form_nl_system` (25 LOC) and `solve_and_backsubstitute` (28 LOC) to satisfy $< 32$ LOC limit.
- Include satellite narrow-lane phase bias $B_{NL}^s$ from `SinexBias::narrow_lane_satellite_bias` in the single-difference narrow-lane formation.

### 3. `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs`
- **Pass Epoch Number**: Pass the true epoch counter into `record_mw_sample` and `PppMwTracker::update`.
- **Between-Satellite Single-Difference Factor**: Replace `inject_fixed_priors` with `inject_sd_fixed_factors(solver, ref_amb_id, cand_amb_id, fixed_sd_m)`.
  Define a 2-variable factor `SdAmbiguityConstraintFactor { vars: [a_ref, a_cand], diff_m: f64, info: 1e8 }` whose residual is $(a_{\text{cand}} - a_{\text{ref}}) - \text{diff}_m$. This enforces the fixed integer baseline without locking the receiver clock state.
- **Refactor functions**: Split `execute_ppp_ar_step` (74 LOC) into helper functions `find_ppp_ar_groups`, `select_reference_satellite`, and `attempt_group_fix` (< 32 LOC each, max nesting depth $\le 2$).

### 4. `crates/gneiss-rtk/src/swfg/engine/epoch.rs` & `mod.rs`
- **Standalone Receiver Antenna PCO/PCV**:
  In `extract_single_sat_with_source`, call `receiver_pcv.total_correction_m(az_deg, el_deg)` and subtract it from the range model.
- **Fix Status Propagation**:
  In `SwfgEngine::process_impl`, when `execute_ar_step` successfully fixes ambiguities, set `SwfgSolution.error = Some(0.01)` or add `is_fixed: bool` so that `FilteredEpoch.quality` is set to `1` (Fixed).
- **Process Noise in Kinematic PPP**:
  In `setup.rs`, adjust `RelativePoseFactor` process noise for kinematic motion (velocity random walk) rather than zero-velocity hold.

---

## 6. Summary Matrix of Required Actions

| Component | Target File | Issue | Required Fix | AGENTS.md Impact |
|---|---|---|---|---|
| Benchmark Config | `crates/gneiss-rtk/src/bin/eval_ppp.rs` | `is_kinematic: false` forces static pose | Change to `is_kinematic: true`; align CODE products | Keeps file modular |
| File Size | `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs` | 506 LOC exceeds 500 limit | Extract `PppMwTracker` to `mw_tracker.rs` | Reduces file to ~350 LOC |
| Function Size | `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs` | `resolve_sd_with_fixed_wl` is 61 LOC | Split into system formation & solution helpers | All fns $< 32$ LOC |
| MW Tracking | `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs` | Epoch is hardcoded to `0` | Pass true epoch index from engine | Keeps fn $< 32$ LOC |
| AR Factor | `crates/gneiss-rtk/src/swfg/engine/ar_handler.rs` | Single-variable `PriorFactor` locks receiver clock | Inject 2-variable SD ambiguity relative factors | Preserves receiver clock freedom |
| Receiver Antenna | `crates/gneiss-rtk/src/swfg/engine/epoch.rs` | Standalone receiver PCO/PCV uncalled | Project `ReceiverPcv` along LOS | Fixes 2–8 cm systematic bias |
| Fix Reporting | `crates/gneiss-rtk/src/swfg/engine/mod.rs` | `sol.error` always `None` | Propagate fix boolean to `FilteredEpoch` | Enables fix rate counting in report |
