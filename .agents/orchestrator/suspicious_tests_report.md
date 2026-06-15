# Suspicious Tests Audit Report — Gneiss Workspace

This report details the findings of a comprehensive static analysis audit of the test suite in the `gneiss` workspace. The goal of this audit was to identify mistakes in test assertions, unrun tests, silent validations, and overly loose tolerances that could potentially hide actual bugs in the GNSS/INS processing codebase.

---

## Executive Summary
The audit across the `gneiss` workspace member crates (`gneiss-core`, `gneiss-rtk`, `gneiss-parsers`, `gneiss-geodesy`, `gneiss-ntrip`, `gneiss-fetch`) and root integration tests revealed:
1. **16 Unrun / Dead Test Functions**: 12 critical math unit tests in `gneiss-rtk` are never run because of missing `#[test]` annotations, and 4 test files are unregistered or uncompiled.
2. **5 Silent Tests (No Verification)**: Tests that execute complex parser or EKF logic but discard the output or perform zero assertions.
3. **3 Trivial / Weak Assertions**: Assertions that are tautologies (always true) or check trivial properties (like matrix dimensions) rather than actual values.
4. **3 Overly Loose Tolerances**: High-precision calculations checked with tolerances large enough to mask incorrect physical constants or coordinate errors.
5. **Logical & Attribute Discrepancies**: Inconsistent comparison operators and redundant test decorators.

---

## 1. Unrun / Dead Test Files and Modules

### A. Dead Integration Test File
- **File**: `/Users/kevin/projects/gneiss/tests/src/ppp_integration.rs` (specifically `test_ppp_skeleton` starting at line 7)
- **Status**: Completely uncompiled/unrun.
- **Cause**: The module is never declared in the library entry point (`/Users/kevin/projects/gneiss/tests/src/lib.rs`). 
- **Impact**: Any bugs or regressions in the PPP integration framework will go completely unnoticed.

### B. Unregistered Scratch/Test Files
- **Files**:
  - `/Users/kevin/projects/gneiss/tests/test_rotation.rs` (lines 1-7)
  - `/Users/kevin/projects/gneiss/tests/test_size.rs` (lines 1-7)
  - `/Users/kevin/projects/gneiss/crates/gneiss-core/tests_rotation.rs` (specifically `test_rot` starting at line 3)
- **Status**: Ignored by Cargo.
- **Cause**: These files define a `main()` entry point instead of `#[test]` attributes (or in the case of `tests_rotation.rs`, define a `#[test]` function but reside outside `src/` or `tests/` target folders) and are not registered in their respective `Cargo.toml` targets.
- **Impact**: Rotation math and layout size checks remain unverified.

---

## 2. Missing `#[test]` Attributes (Critical Unrun Unit Tests)

In `/Users/kevin/projects/gneiss/crates/gneiss-rtk/src/engine/updater_math.rs`, **12 unit tests** containing critical assertions are defined inside `#[cfg(test)]` blocks but are missing the `#[test]` decorator. They are completely ignored by `cargo test`.

| Line | Function Name | Missing Verification |
|---|---|---|
| 180 | `test_evaluate_post_fit_outliers` | Evaluation of EKF state post-fit residuals |
| 222 | `threshold_tests::test_get_pre_fit_threshold` | Pre-fit threshold retrieval logic |
| 237 | `threshold_tests::test_compute_scalar_thresholds` | Scalar threshold computation |
| 277 | `loose_coupling_tests::test_compute_loose_coupling_innovations` | Loose coupling innovation logic |
| 321 | `filter_tests::test_filter_pre_fit_residuals` | Residual pre-filtering |
| 384 | `check_pre_fit_tests::test_check_pre_fit_residual` | Pre-fit residual validator |
| 422 | `missed_mutant_tests::test_evaluate_post_fit_outliers_exact_abs_thresh` | Boundary absolute threshold check |
| 440 | `missed_mutant_tests::test_evaluate_post_fit_outliers_meas_type_3_abs_outlier` | Measurement type outlier check |
| 458 | `missed_mutant_tests::test_evaluate_post_fit_outliers_equal_ratio` | Ratio outlier check |
| 477 | `missed_mutant_tests::test_evaluate_post_fit_outliers_exact_ratio_1` | Boundary ratio check |
| 494 | `missed_mutant_tests::test_evaluate_post_fit_outliers_exact_valid_count_4` | Minimum valid measurement check |
| 511 | `missed_mutant_tests::test_populate_loosely_coupled_jacobian` | Loosely coupled Jacobian structure |

*Note: Running `cargo build` confirms the compiler warning: `warning: function ... is never used`.*

---

## 3. Silent Tests (Zero Assertions or Discarded Outputs)

### A. Ignored Double-Difference Clock Elimination
- **File**: `/Users/kevin/projects/gneiss/crates/gneiss-rtk/src/estimators/ekf/filter.rs` (lines 491-507, `test_double_difference_eliminates_clocks`)
- **Issue**: Calls `compute_double_difference` at line 506 with complex simulated observations to test receiver/satellite clock elimination, but assigns the result to `_dd` and does **not** assert anything about it.
- **Impact**: Swapped variables, calculation bugs, or incorrect frequency constants in the clock-elimination math will pass silently.

### B. Empty EKF Stability Test
- **File**: `/Users/kevin/projects/gneiss/crates/gneiss-rtk/src/tests_ekf.rs` (lines 8-27, `test_ekf_update_stability`)
- **Issue**: Declares state variables and matrices but contains **no assertions** and never calls any EKF updater or filter routines.
- **Impact**: It acts as a dead stub and provides zero validation.

### C. Ignored GPS Ephemeris Verification
- **File**: `/Users/kevin/projects/gneiss/crates/gneiss-parsers/tests/integration_test.rs` (lines 37-39, 66-72)
- **Issue**: In the RTCM3 parser integration test, the result of decoding message 1019 (`parse_1019(frame.payload)`) is assigned to `_eph` at line 37 and discarded. There is no assertion verifying that `eph_frames > 0` or checking the fields inside `_eph`.
- **Impact**: If the GPS ephemeris parser fails to decode correctly, the test will still pass because it only asserts MSM frame count.

### D. No-assertion Skeleton Replay Test
- **File**: `/Users/kevin/projects/gneiss/tests/src/urbannav_integration.rs` (lines 8-39, `test_urbannav_tst_replay_skeleton`)
- **Issue**: Replays a UBX dataset file in a `while let Ok((rem, _frame)) = parse_ubx_frame(remaining)` loop. It increments a parsed counter but has no assertion checking that any frames were successfully decoded or that the dataset was actually parsed.
- **Impact**: If the dataset path doesn't exist, it returns early. If the dataset does exist but the parser immediately fails on the first byte, the loop terminates, `parsed_count` is 0, and the test passes silently.

### E. Ignored Stub Test
- **File**: `/Users/kevin/projects/gneiss/crates/gneiss-rtk/src/engine/measurement.rs` (line 827, `test_compute_dd_carrier_phase`)
- **Issue**: The test body is completely empty (`{}`) and marked `#[ignore]`.
- **Impact**: Zero verification of carrier phase double-differencing.

---

## 4. Trivial / Weak Assertions

### A. Tautological Assertion
- **File**: `/Users/kevin/projects/gneiss/tests/src/ppp_integration.rs` (line 64)
- **Issue**: `assert!(result.is_err() || result.is_ok());`
- **Impact**: A `Result` is always either `Ok` or `Err`. This statement is a tautology (evaluates to `assert!(true)`) and verifies absolutely nothing about the output of `process_epoch`.

### B. Dimension-only Robust Inversion Assertion
- **File**: `/Users/kevin/projects/gneiss/crates/gneiss-rtk/src/math/inversion.rs` (line 37, `test_invert_matrix_robust`)
- **Issue**: Checks robust matrix inversion on a singular matrix by checking the dimensions only: `assert!(inv_singular.nrows() == 3);`.
- **Impact**: It fails to check if the robust fallback logic returned correct mathematical values (like a regularized identity or pseudo-inverse elements). If the code returns numerical garbage or zeros, the test still passes.

### C. Overly Broad Configuration Assertion
- **File**: `/Users/kevin/projects/gneiss/tests/src/lib.rs` (line 52)
- **Issue**: `assert!(matches!(engine.config.mode, EngineMode::RtkIns | EngineMode::SppIns | EngineMode::PppIns));`
- **Impact**: The test explicitly sets `config.mode = EngineMode::RtkIns`. Checking that the mode is one of three valid modes is too weak and would fail to detect a bug that incorrectly sets it to `SppIns` or `PppIns`.

---

## 5. Overly Loose Numerical Tolerances

### A. Huge Tropospheric Delay Range
- **File**: `/Users/kevin/projects/gneiss/crates/gneiss-core/src/atmosphere.rs` (line 205, `test_tropo_delay`)
- **Issue**: `assert!(delay > 2.0 && delay < 10.0);`
- **Impact**: The exact Saastamoinen tropospheric delay under the test input conditions is `~4.94025` meters. A check of `(2.0, 10.0)` represents an error range of -59% to +102%. Major bugs in mathematical scale factors or model coefficients could occur, and the test would still pass.

### B. Wide Wavelength Tolerance Masking Constellation Mismatch
- **File**: `/Users/kevin/projects/gneiss/crates/gneiss-core/src/signal.rs` (line 75, `test_glonass_fdma_wavelengths`)
- **Issue**: `assert!((w1 - 0.1873).abs() < 0.01);`
- **Impact**: The exact GLONASS L1 channel -4 wavelength is `~0.187399` m. GPS L1 wavelength is `~0.19029` m. The difference is `0.00299` m, which is well within the `0.01` tolerance. If the code mistakenly calculates the GPS wavelength instead of the GLONASS channel wavelength, the assertion will still pass.

### C. Translation-only Helmert Tolerance
- **File**: `/Users/kevin/projects/gneiss/crates/gneiss-geodesy/src/helmert.rs` (lines 136-138)
- **Issue**: Uses a tolerance of `1e-4` (0.1 mm) when checking translation results under zero rotation and scale parameters.
- **Impact**: Translation-only transformations represent simple floating-point addition and are mathematically exact. Tight tolerances of `1e-9` (as used in other Helmert tests) should be enforced. A tolerance of `1e-4` could mask sub-millimeter offsets or rounding errors.

---

## 6. Assertion Logic and Other Discrepancies

### A. Strict Inequality Check in GDOP vs PDOP
- **File**: `/Users/kevin/projects/gneiss/crates/gneiss-core/src/dop.rs` (line 113)
- **Issue**: `assert!(dop.gdop > dop.pdop, "GDOP must be >= PDOP");`
- **Impact**: The text assertion correctly states `GDOP >= PDOP`, but the code checks for strict inequality `>`. If TDOP is exactly 0, GDOP equals PDOP, and the test will crash on valid physical outputs.

### B. Missing Epoch TOW Verification in RINEX NAV Parser
- **File**: `/Users/kevin/projects/gneiss/crates/gneiss-parsers/src/rinex.rs` (line 634)
- **Issue**: `test_parse_rinex_3_nav_date` verifies that `eph.toe().week == 2137` but omits checking the Time of Week (TOW) field.
- **Impact**: Slicing logic for RINEX 3 seconds field parsing could truncate decimal values or compute incorrect TOW without failing the test.

### C. Redundant `#[test]` Decorator
- **File**: `/Users/kevin/projects/gneiss/crates/gneiss-rtk/src/engine/ppp_fg.rs` (lines 444-446)
- **Issue**: Contains duplicate `#[test]` attributes:
  ```rust
  #[test]

  #[test]
  fn test_snr_scale() {
  ```

### D. Verbatim Test Duplication
- **File**: `/Users/kevin/projects/gneiss/crates/gneiss-rtk/src/measurements/doppler.rs`
- **Issue**: The tests `test_doppler_exact_range_mutant` (lines 366 and 573), `test_doppler_short_range_continue_mutant` (lines 404 and 611), and `test_missing_ephemeris_first_mutant` (lines 445 and 652) are duplicated verbatim in both `mod tests` and `mod missing_eph_tests`.

---

## Conclusion & Action Plan
This audit highlights several gaps where tests are completely unrun or assertion logic is weak/trivial. 

To improve test coverage and accuracy:
1. **Register the missing files** in `/Users/kevin/projects/gneiss/tests/src/lib.rs` and crate `Cargo.toml` configurations.
2. **Add missing `#[test]` decorators** to the 12 updater math functions.
3. **Tighten tolerances** in geodesy transformations, tropospheric delay bounds, and GLONASS wavelength checks.
4. **Implement robust assertions** for silent tests to verify values and parameters (e.g., verifying RTCM3 GPS ephemeris decoded values, EKF update calculations, and double-difference clock cancellation results).
