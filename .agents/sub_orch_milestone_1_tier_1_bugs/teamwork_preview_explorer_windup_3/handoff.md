# Handoff Report — Phase Wind-Up Correction Sign Bug

## 1. Observation
We observed the application of the phase wind-up correction in the following locations:

- **RTK Double-Differenced correction (`crates/gneiss-rtk/src/engine/measurement.rs`):**
  Lines 42–49:
  ```rust
  fn apply_windup_to_obs(obs: &mut DdObservation, windup: f64) {
      if let Some(cp) = &mut obs.cp_l1 {
          *cp += windup;
      }
      if let Some(cp2) = &mut obs.cp_l2 {
          *cp2 += windup;
      }
  }
  ```

- **PPP Pre-alignment (`crates/gneiss-rtk/src/engine/ppp.rs`):**
  Lines 495–504:
  ```rust
  let l_meas = if sat.is_iono_free && sat.cp2.is_some() {
      crate::engine::ppp_math::compute_iono_free(
          (cp1 + wup) * sat.lam1,
          (sat.cp2.unwrap() + wup) * sat.lam2,
          sat.f1,
          sat.f2,
      )
  } else {
      (cp1 + wup) * sat.lam1
  };
  ```
  Lines 545–546:
  ```rust
  let l1_m = (cp1 + wup) * sat.lam1;
  let l2_m = (sat.cp2.unwrap() + wup) * sat.lam2;
  ```
  Lines 588–589:
  ```rust
  let l1_meas = (cp1 + wup) * sat.lam1;
  let l2_meas = (sat.cp2.unwrap() + wup) * sat.lam2;
  ```

- **PPP IEKF solver (`crates/gneiss-rtk/src/engine/ppp_iekf.rs`):**
  Lines 897–898:
  ```rust
  let windup = *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0);
  let l_meas = (cp1 + windup) * sat.lam1;
  ```
  Lines 1026:
  ```rust
  let res_l1 = (sat.cp1.unwrap() + windup) * sat.lam1 - (expected_base - i1 + n1);
  ```
  Line 1043:
  ```rust
  let res_l2 = (sat.cp2.unwrap() + windup) * sat.lam2 - (expected_base - gamma * i1 + n2);
  ```

In all cases, the computed wind-up value `wup`/`windup` is added directly to the raw carrier phase `cp1`/`cp2` values.

---

## 2. Logic Chain
1. According to the electromagnetic and geometric equations of GNSS antenna rotation (Wu et al. 1992), phase wind-up represents a phase shift added to the received carrier phase signal during the transmission/reception process:
   $$\phi_{\text{meas}} = \phi_{\text{geom}} + N + \phi_{\text{windup}}$$
2. Therefore, to compare the observation against the geometric model in a filter (such as an EKF or factor graph), the wind-up correction must be subtracted from the observed carrier phase to obtain the corrected phase matching the geometric propagation model:
   $$\phi_{\text{corrected}} = \phi_{\text{meas}} - \phi_{\text{windup}}$$
3. The observations show that `gneiss-rtk` consistently adds the correction: `cp + windup`.
4. As a result, the sign of the phase wind-up correction applied to the observations is opposite to the correct physical model.

---

## 3. Caveats
- We assume that the computed wind-up value from `gneiss_core::windup::phase_windup` has a positive physical sign matching the rotation angle (direction of the right-handed circularly polarized wave). A review of `gneiss_core::windup::phase_windup` shows that it follows the standard Wu et al. (1992) convention where positive wind-up corresponds to positive rotation. Therefore, the sign discrepancy resides in the correction application step in `gneiss-rtk` rather than the core calculation in `gneiss-core`.
- The factor graph solver in `crates/gneiss-rtk/src/estimators/factor_graph/gnss_factors.rs` currently implements a placeholder/testing model for carrier phase that does not incorporate wind-up. Once wind-up is added to the factor graph solver, it should similarly be subtracted.

---

## 4. Conclusion
The phase wind-up correction in `gneiss-rtk` has an opposite sign because it is added to the raw carrier phase observations (`cp + windup` / `cp += windup`), whereas the correct physical model requires subtracting it (`cp - windup` / `cp -= windup`).

Actionable fix:
- In `crates/gneiss-rtk/src/engine/measurement.rs`: Replace `*cp += windup` with `*cp -= windup`.
- In `crates/gneiss-rtk/src/engine/ppp.rs`: Replace `+ wup` with `- wup` when applying it to carrier phase.
- In `crates/gneiss-rtk/src/engine/ppp_iekf.rs`: Replace `+ windup` with `- windup` in the carrier phase residual and measurement prediction logic.

---

## 5. Verification Method
- **Command**: Run `cargo test` to verify that the existing tests pass.
- **Inspect**: Inspect the files changed (`measurement.rs`, `ppp.rs`, and `ppp_iekf.rs`) to confirm that all instances of adding phase wind-up (`+ wup` or `+= windup`) have been changed to subtraction (`- wup` or `-= windup`).
- **Regression Test**: Implement the unit test `test_phase_windup_correction_sign_rtk` (detailed in `analysis.md`) in `crates/gneiss-rtk/src/engine/measurement.rs` and verify it passes. The test is validated when a positive wind-up value decreases the corrected phase.
