# Handoff Report — Bug 2 Velocity-Attitude Transition Sign Mismatch Audit

## 1. Observation
- **File Checked**: `crates/gneiss-rtk/src/engine/predictor.rs` at line 86:
  ```rust
  let vel_att = f_e_skew * dt;
  ```
- **File Checked**: `crates/gneiss-rtk/src/engine/tests_predictor.rs` at lines 159–205:
  ```rust
  #[test]
  fn test_transition_matrix_velocity_attitude_coupling() {
      let time = GpsTime::new(2000, 0.0);
      let pos = Coordinate::new(
          Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
          Datum::WGS84,
          Frame::ECEF,
          time,
      );
      let mut state = RtkState::new(time, pos, 1.0);
      state.attitude = UnitQuaternion::identity();
      state.accel_bias = Vector3::zeros();

      let dt = 0.5;
      let accel = Vector3::new(1.0, 2.0, 3.0);
      let imu_meas = ImuMeasurement::new(0, accel, Vector3::zeros());

      let phi = predictor::compute_transition_matrix(&state, dt, &[imu_meas]);

      // Compute the expected coupling block: skew(f_e) * dt
      // Since state.attitude is identity and accel_bias is zero, f_e = accel = [1.0, 2.0, 3.0]
      // skew(f_e) is:
      // [ 0.0, -3.0,  2.0]
      // [ 3.0,  0.0, -1.0]
      // [-2.0,  1.0,  0.0]
      //
      // Multiplying by dt (0.5):
      // [ 0.0, -1.5,  1.0]
      // [ 1.5,  0.0, -0.5]
      // [-1.0,  0.5,  0.0]
      let expected_vel_att = nalgebra::Matrix3::new(
          0.0, -1.5, 1.0,
          1.5, 0.0, -0.5,
          -1.0, 0.5, 0.0,
      );

      for r in 0..3 {
          for c in 0..3 {
              let actual_val = phi[(3 + r, 6 + c)];
              let expected_val = expected_vel_att[(r, c)];
              assert!(
                  (actual_val - expected_val).abs() < 1e-10,
                  "Mismatch at ({}, {}): expected {}, got {}",
                  r, c, expected_val, actual_val
              );
          }
      }
  }
  ```
- **Tool Commands Run**:
  - `cargo test --package gneiss-rtk --lib -- engine::tests_predictor`
    - Result: `5 passed; 0 failed; 0 ignored;`
  - `cargo test --workspace`
    - Result: `258 passed; 0 failed; 1 ignored;`

## 2. Logic Chain
1. The EKF perturbation equation for velocity in ECEF is given by $\delta \dot{v} = + [f_e \times] \psi$, where $f_e$ is specific force in ECEF and $\psi$ is the attitude error vector. This yields a positive transition matrix coupling term $+ f_e\_skew * dt$.
2. In `predictor.rs`, line 86 defines `let vel_att = f_e_skew * dt;` with a positive sign (no leading minus sign). This matches the correct physics model.
3. In `tests_predictor.rs`, `test_transition_matrix_velocity_attitude_coupling` sets up `accel = [1, 2, 3]` and asserts that the transition matrix's velocity-attitude block matches the manually derived matrix `expected_vel_att` with a positive coupling sign. If the negative sign was present, the actual values would have inverted signs and would fail the assertion.
4. Running `cargo test --workspace` compiles and passes all tests successfully, demonstrating no regressions.
5. Analysis of the source code shows no hardcoded test shortcuts, bypasses, or facade implementations.

## 3. Caveats
- Assumes that the EKF operates at a high enough rate (e.g. >= 100Hz) such that the first-order transition matrix approximation $\Phi \approx I + F dt$ remains valid. This is standard in all practical GNSS/INS integrations.

## 4. Conclusion

## Forensic Audit Report

**Work Product**: crates/gneiss-rtk/src/engine/predictor.rs and crates/gneiss-rtk/src/engine/tests_predictor.rs
**Profile**: General Project
**Verdict**: CLEAN

### Phase Results
- **Hardcoded output detection**: PASS — No hardcoded test results or bypass strings.
- **Facade detection**: PASS — Full dynamic implementation of transition matrix computing.
- **Pre-populated artifact detection**: PASS — No fabricated test logs or artifacts exist.
- **Build and run**: PASS — Cargo workspace builds and test suite passes.
- **Output verification**: PASS — Correct positive coupling sign verified mathematically and behavioral regression test asserted.
- **Dependency audit**: PASS — Uses standard library/core crates, no external solver or blackbox delegation for EKF prediction logic.

## 5. Verification Method
To independently verify the audit:
1. View the source file `crates/gneiss-rtk/src/engine/predictor.rs` at line 86 to verify that `let vel_att = f_e_skew * dt;` is defined with no leading negative sign.
2. View the test file `crates/gneiss-rtk/src/engine/tests_predictor.rs` starting at line 159 to check the implementation of `test_transition_matrix_velocity_attitude_coupling`.
3. Run the specific test using:
   ```bash
   cargo test --package gneiss-rtk --lib -- engine::tests_predictor::tests::test_transition_matrix_velocity_attitude_coupling
   ```
