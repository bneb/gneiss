# Audit Report: GNSS-RTK Crate Test Suite Static Analysis

This report documents the findings of a read-only static analysis audit of the test suite in `crates/gneiss-rtk` for non-engine directories.

---

## 1. Observation

During the static analysis audit, the following suspicious test files and code blocks were directly observed:

### Observation A: Silent Test in `crates/gneiss-rtk/src/estimators/ekf/filter.rs`
- **File Path**: `crates/gneiss-rtk/src/estimators/ekf/filter.rs`
- **Line Numbers**: 491–507
- **Test Function**: `test_double_difference_eliminates_clocks`
- **Verbatim Code**:
```rust
    #[test]
    fn test_double_difference_eliminates_clocks() {
        let sat_ref = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat_a = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let true_r_rover_ref = 20_000_000.0;
        let true_r_rover_a   = 21_000_000.0;
        let true_r_base_ref  = 20_005_000.0;
        let true_r_base_a    = 21_004_000.0;
        let rover_clk = 300.0; let base_clk = -150.0;
        let sat_ref_clk = 1000.0; let sat_a_clk = -500.0;
        let rover_ref_obs = DdObservation { sat: sat_ref, pr_l1: true_r_rover_ref + rover_clk - sat_ref_clk, pr_l2: None, cp_l1: Some(0.0), cp_l2: None, doppler: 0.0, snr: 45.0, locktime: Some(1000) };
        let rover_a_obs = DdObservation { sat: sat_a, pr_l1: true_r_rover_a + rover_clk - sat_a_clk, pr_l2: None, cp_l1: Some(0.0), cp_l2: None, doppler: 0.0, snr: 45.0, locktime: Some(1000) };
        let base_ref_obs = DdObservation { sat: sat_ref, pr_l1: true_r_base_ref + base_clk - sat_ref_clk, pr_l2: None, cp_l1: Some(0.0), cp_l2: None, doppler: 0.0, snr: 45.0, locktime: Some(1000) };
        let base_a_obs = DdObservation { sat: sat_a, pr_l1: true_r_base_a + base_clk - sat_a_clk, pr_l2: None, cp_l1: Some(0.0), cp_l2: None, doppler: 0.0, snr: 45.0, locktime: Some(1000) };
        let f1 = 1575.42e6;
        let f2 = 1227.60e6;
        let _dd = compute_double_difference(&rover_ref_obs, &rover_a_obs, &base_ref_obs, &base_a_obs, f1, f2, f1, f2);
    }
```

### Observation B: Empty Silent Test in `crates/gneiss-rtk/src/tests_ekf.rs`
- **File Path**: `crates/gneiss-rtk/src/tests_ekf.rs`
- **Line Numbers**: 8–27
- **Test Function**: `test_ekf_update_stability`
- **Verbatim Code**:
```rust
    #[test]
    fn test_ekf_update_stability() {
        let initial_pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0));
        let mut _state = RtkState::new(GpsTime::new(0, 0.0), initial_pos, 10.0);
        
        // Z = [1.0, 2.0], H = [I | 0] (first 2 states are X and Y)
        let _z = DVector::from_vec(vec![1.0, 2.0]);
        let mut h = DMatrix::zeros(2, 18); 
        h[(0, 0)] = 1.0; 
        h[(1, 1)] = 1.0;
        let _r = DMatrix::from_diagonal(&DVector::from_vec(vec![0.1, 0.1]));
        
        // Small subset of states to test update logic
        // P initial = diag(10, 10)
        // K = P * H^T * (H P H^T + R)^-1
        // K = 10 * 1 * (10 + 0.1)^-1 = 10 / 10.1 = 0.990099
        // dx = K * Z = 0.99 * 1.0 = 0.990099
        
    }
```

### Observation C: Weak Assertion in `crates/gneiss-rtk/src/math/inversion.rs`
- **File Path**: `crates/gneiss-rtk/src/math/inversion.rs`
- **Line Numbers**: 29–38
- **Test Function**: `test_invert_matrix_robust`
- **Verbatim Code (relevant snippet)**:
```rust
        let mut m_singular = CovMatrix::zeros(3, 3);
        m_singular[(0, 0)] = 1.0;
        let inv_singular = invert_matrix_robust(&m_singular);
        assert!(inv_singular.nrows() == 3);
```

### Observation D: Verbatim Test Duplication in `crates/gneiss-rtk/src/measurements/doppler.rs`
- **File Path**: `crates/gneiss-rtk/src/measurements/doppler.rs`
- **Test Functions**: `test_doppler_exact_range_mutant`, `test_doppler_short_range_continue_mutant`, and `test_missing_ephemeris_first_mutant`
- **Verbatim Code (first occurrence, in `mod tests` lines 365–478)**:
```rust
    #[test]
    fn test_doppler_exact_range_mutant() { ... }

    #[test]
    fn test_doppler_short_range_continue_mutant() { ... }

    #[test]
    fn test_missing_ephemeris_first_mutant() { ... }
```
- **Verbatim Code (second occurrence, in `mod missing_eph_tests` lines 572–685)**:
```rust
    #[test]
    fn test_doppler_exact_range_mutant() { ... }

    #[test]
    fn test_doppler_short_range_continue_mutant() { ... }

    #[test]
    fn test_missing_ephemeris_first_mutant() { ... }
```

---

## 2. Logic Chain

### For Observation A (`test_double_difference_eliminates_clocks`):
1. **Fact**: The test setups reference observations (`rover_ref_obs`, `rover_a_obs`, `base_ref_obs`, `base_a_obs`) with intentionally added receiver clock biases (`rover_clk = 300.0`, `base_clk = -150.0`) and satellite clock biases (`sat_ref_clk = 1000.0`, `sat_a_clk = -500.0`).
2. **Fact**: The test invokes `compute_double_difference` to verify that double differencing eliminates these clock biases.
3. **Fact**: The return value `_dd` of the function is completely ignored and no assertions check `_dd` (e.g., verifying that the pseudorange DD equals the geometric range DD, which is exactly `1000.0`).
4. **Conclusion**: This is a silent test with zero verification. It will pass even if the double-difference math incorrectly leaves clock biases or introduces regressions.

### For Observation B (`test_ekf_update_stability`):
1. **Fact**: The test initiates variables (`_state`, `_z`, `h`, `_r`) but never calls the EKF state correction/update function (`updater::update` or similar).
2. **Fact**: There are no assertions in the test body at all.
3. **Conclusion**: This is a completely empty silent test that verifies absolutely nothing about EKF update stability.

### For Observation C (`test_invert_matrix_robust`):
1. **Fact**: The robust inversion function `invert_matrix_robust` handles singular matrices by falling back to SVD pseudo-inverse or regularized identity multiplication.
2. **Fact**: The test checks singular matrix handling by asserting `assert!(inv_singular.nrows() == 3)`.
3. **Fact**: Checking only the row dimension of a matrix does not verify any of the matrix element values. If the function returned an all-zeros matrix, an identity matrix, or numerical garbage, the assertion would still pass.
4. **Conclusion**: This is a trivial/weak assertion that fails to verify the correctness of the robust matrix inversion fallback values.

### For Observation D (`doppler.rs` duplicated tests):
1. **Fact**: The three mutant tests are duplicated word-for-word in two separate modules in the same file.
2. **Conclusion**: This is redundant, bloats compile times, and increases maintenance effort since any test updates must be synchronised across both modules.

---

## 3. Caveats

- The engine tests under `crates/gneiss-rtk/src/engine` were explicitly excluded from this audit as per the instructions, so any potential issues there have not been inspected.
- No source code files were modified; this report serves as a read-only assessment of the current state of the test suite.

---

## 4. Conclusion

The audit identified two critical silent tests with no assertions (`test_double_difference_eliminates_clocks` and `test_ekf_update_stability`), one weak/trivial dimension assertion (`test_invert_matrix_robust`), and three verbatim duplicated tests. 

- **Actionable Steps**:
  1. Add assertions to `test_double_difference_eliminates_clocks` to verify that `_dd` correctly eliminates receiver/satellite clocks (e.g. asserting `(dd.pr_l1 - 1000.0).abs() < 1e-6`).
  2. Implement EKF update logic invocation and subsequent state correction verification in `test_ekf_update_stability`.
  3. Replace the dimension assertion in `test_invert_matrix_robust` with checks on actual elements of `inv_singular` to verify the singular matrix fallback values are mathematically correct.
  4. Remove the duplicate tests from `doppler.rs` (either from `mod tests` or `mod missing_eph_tests`).

---

## 5. Verification Method

To independently verify these findings:

1. **Commands to Run**:
   Run the project test suite using cargo:
   ```bash
   cargo test --package gneiss-rtk --lib
   ```
   All tests will pass, confirming that the silent and weak assertions do not fail.

2. **Files to Inspect**:
   - `crates/gneiss-rtk/src/estimators/ekf/filter.rs` (lines 491–507)
   - `crates/gneiss-rtk/src/tests_ekf.rs` (lines 8–27)
   - `crates/gneiss-rtk/src/math/inversion.rs` (lines 29–38)
   - `crates/gneiss-rtk/src/measurements/doppler.rs` (lines 365–478 and lines 572–685)
