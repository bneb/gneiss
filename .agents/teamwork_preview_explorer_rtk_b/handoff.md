# Test Suite Static Analysis Audit Report — crates/gneiss-rtk/src/engine

This report presents findings from a read-only static analysis audit of the test suite under the `crates/gneiss-rtk/src/engine` directory.

---

## 1. Observations

### Observation A: Missing `#[test]` attributes in `updater_math.rs`
In `crates/gneiss-rtk/src/engine/updater_math.rs`, multiple functions that start with `test_` and contain critical assertions are defined inside `#[cfg(test)]` modules (or in the parent module) but do not have the `#[test]` attribute. Consequently, they are completely omitted by the test runner.

* **File Path**: `crates/gneiss-rtk/src/engine/updater_math.rs`
* **Unrun Functions**:
  1. **Line 180**:
     ```rust
     fn test_evaluate_post_fit_outliers() {
     ```
  2. **Line 222** (inside `mod threshold_tests`):
     ```rust
     fn test_get_pre_fit_threshold() {
     ```
  3. **Line 237** (inside `mod threshold_tests`):
     ```rust
     fn test_compute_scalar_thresholds() {
     ```
  4. **Line 277** (inside `mod loose_coupling_tests`):
     ```rust
     fn test_compute_loose_coupling_innovations() {
     ```
  5. **Line 321** (inside `mod filter_tests`):
     ```rust
     fn test_filter_pre_fit_residuals() {
     ```
  6. **Line 384** (inside `mod check_pre_fit_tests`):
     ```rust
     fn test_check_pre_fit_residual() {
     ```
  7. **Line 422** (inside `mod missed_mutant_tests`):
     ```rust
     fn test_evaluate_post_fit_outliers_exact_abs_thresh() {
     ```
  8. **Line 440** (inside `mod missed_mutant_tests`):
     ```rust
     fn test_evaluate_post_fit_outliers_meas_type_3_abs_outlier() {
     ```
  9. **Line 458** (inside `mod missed_mutant_tests`):
     ```rust
     fn test_evaluate_post_fit_outliers_equal_ratio() {
     ```
  10. **Line 477** (inside `mod missed_mutant_tests`):
     ```rust
     fn test_evaluate_post_fit_outliers_exact_ratio_1() {
     ```
  11. **Line 494** (inside `mod missed_mutant_tests`):
     ```rust
     fn test_evaluate_post_fit_outliers_exact_valid_count_4() {
     ```
  12. **Line 511** (inside `mod missed_mutant_tests`):
     ```rust
     fn test_populate_loosely_coupled_jacobian() {
     ```

* **Tool Outputs & Verification**:
  * Running `cargo test --package gneiss-rtk` showed 111 passed tests, none of which corresponded to `engine::updater_math::*`.
  * Running `cargo build` generated a compiler warning confirming that the outer test function is dead code:
    ```
    warning: function `test_evaluate_post_fit_outliers` is never used
       --> crates/gneiss-rtk/src/engine/updater_math.rs:180:8
        |
    180 |     fn test_evaluate_post_fit_outliers() {
        |        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ```

### Observation B: Empty Test with no assertions in `measurement.rs`
An empty test function decorated with `#[test]` and `#[ignore]` exists in `measurement.rs`, verifying nothing.

* **File Path**: `crates/gneiss-rtk/src/engine/measurement.rs`
* **Line 825-827**:
  ```rust
  #[test]
  #[ignore]
  fn test_compute_dd_carrier_phase() {}
  ```

### Observation C: Assertions checking `is_ok()` / `is_err()`
Some tests verify if a function returns an `Ok` or `Err` result without asserting the exact contents of the result. However, for functions returning `Result<(), EkfError>`, there is no contained value to assert when it is `Ok`.

1. **File Path**: `crates/gneiss-rtk/src/engine/spp_tight.rs`
   * **Line 297**:
     ```rust
     assert!(process_spp_tightly_coupled(&mut engine, &obs).is_err());
     ```
     *Checks that SPP tightly coupled returns an error on startup with no observations, but does not verify the specific error variant.*
2. **File Path**: `crates/gneiss-rtk/src/engine/tests_updater.rs`
   * **Line 191** (`test_fix_and_hold_updates_imu_states`):
     ```rust
     assert!(res.is_ok());
     ```
   * **Line 252** (`test_loosely_coupled_jacobian`):
     ```rust
     assert!(res.is_ok());
     ```
   * **Line 276** (`test_update_loosely_coupled_huber`):
     ```rust
     assert!(res.is_ok(), "Huber scaling should prevent rejection of the huge error");
     ```
     *These verify that EKF updates succeed. Since the return type of these operations is `Result<(), EkfError>`, there is no contained success value to inspect.*

### Observation D: Tolerances in Approximations
Several tests use looser tolerances (e.g. > 0.01) for float comparisons:

1. **File Path**: `crates/gneiss-rtk/src/engine/tests_predictor.rs`
   * **Line 24**:
     ```rust
     assert!((yaw.abs() - core::f64::consts::FRAC_PI_2).abs() < 0.1);
     ```
     *Used to verify Euler yaw rotation over 100 IMU prediction epochs. 0.1 radians is reasonable since IMU integration propagates over time.*
2. **File Path**: `crates/gneiss-rtk/src/engine/jacobian_verify.rs`
   * **Line 267**:
     ```rust
     assert!(trace.abs() / diag_scale < 0.1, ...);
     ```
     *Used to verify that the trace of the numerical gravity Jacobian is close to 0 (since gravitational Laplacian in free space is 0). Earth's J2 oblateness terms prevent it from being exactly 0; a 10% relative tolerance is physically justified.*
   * **Line 317**:
     ```rust
     assert!(err < 0.01, ...);
     ```
     *Used to compare numerical vs analytical state transition matrix Φ vel->pos block. The tolerance of 1% is due to numerical differentiation approximations.*
   * **Line 368**:
     ```rust
     assert!(diff < 1e-3, ...);
     ```
     *Checks the GNSS lever arm Jacobian against a non-linear attitude perturbation. 1e-3 is reasonable because finite differences introduce O(eps) truncation and rounding errors.*

### Observation E: Redundant / Empty Attribute in `ppp_fg.rs`
There is a redundant `#[test]` attribute before `test_snr_scale` in `ppp_fg.rs`.

* **File Path**: `crates/gneiss-rtk/src/engine/ppp_fg.rs`
* **Line 444-446**:
  ```rust
  #[test]

  #[test]
  fn test_snr_scale() {
  ```

---

## 2. Logic Chain

1. In Rust, a function inside a test module or test file must have the `#[test]` (or similar test framework) attribute to be picked up by `cargo test`.
2. Static inspection of `crates/gneiss-rtk/src/engine/updater_math.rs` reveals 12 functions starting with `test_` that do not have `#[test]` attributes (Observation A).
3. The lack of test execution for these 12 functions was confirmed by running `cargo test --package gneiss-rtk` (where no `engine::updater_math::*` tests are listed) and `cargo build` (which flagged `test_evaluate_post_fit_outliers` as dead/unused code).
4. Similarly, `test_compute_dd_carrier_phase` in `crates/gneiss-rtk/src/engine/measurement.rs:827` has an empty body `{}` and contains no assertions (Observation B).
5. These unrun/empty tests result in a significant gap in test coverage, meaning critical math code (pre-fit residual filtering, post-fit outlier evaluation, lever-arm Jacobian computation, carrier phase double-difference calculation) is not being verified.
6. Other assertions checking `.is_ok()`/`.is_err()` (Observation C) and approximation tolerances (Observation D) are physically/mathematically justified and do not hide bugs.

---

## 3. Caveats

* We assumed that there are no external integrations or custom test runners (e.g. integration scripts) that manually call the missing `test_*` functions in `updater_math.rs`. However, since they are private to the module (declared as `fn test_...` without `pub`), they cannot be accessed outside the module, confirming they are unused.
* We assumed that the empty test `test_compute_dd_carrier_phase` was intended to be implemented but left as a stub (marked `#[ignore]`). It is possible it was intentionally left as a placeholder.

---

## 4. Conclusion

1. **Critical Defect**: The test suite has a major coverage gap due to 12 missing `#[test]` attributes in `crates/gneiss-rtk/src/engine/updater_math.rs`. These functions test crucial GNSS/INS EKF updater math (ratio thresholds, pre-fit residual filtering, outlier detection, and loosely coupled jacobians) but are completely bypassed by `cargo test`.
2. **Defect**: `test_compute_dd_carrier_phase` in `crates/gneiss-rtk/src/engine/measurement.rs:827` is a silent test with an empty body, resulting in zero verification of double-differenced carrier phase measurements.
3. **Style Issue**: There is a redundant duplicate `#[test]` attribute at lines 444-446 in `crates/gneiss-rtk/src/engine/ppp_fg.rs`.

---

## 5. Verification Method

To verify these findings:
1. Run `cargo test --package gneiss-rtk`. Inspect the list of executed tests; verify that no tests from `engine::updater_math` are executed.
2. Run `cargo build --package gneiss-rtk`. Inspect compiler warnings and verify that `test_evaluate_post_fit_outliers` is flagged as an unused function.
3. Open `crates/gneiss-rtk/src/engine/updater_math.rs` and inspect lines 180, 222, 237, 277, 321, 384, 422, 440, 458, 477, 494, and 511 to confirm the absence of `#[test]` attributes.
4. Open `crates/gneiss-rtk/src/engine/measurement.rs` at line 827 to confirm `test_compute_dd_carrier_phase` has an empty body `{}`.
