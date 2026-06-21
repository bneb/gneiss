# Changes - Bug 2: Velocity-Attitude Transition Sign Mismatch

## Modified Files

### `crates/gneiss-rtk/src/engine/predictor.rs`
- Negated the attitude-to-velocity coupling block computation inside the `compute_transition_matrix` function:
  ```rust
  let vel_att = -f_e_skew * dt;
  ```

### `crates/gneiss-rtk/src/engine/tests_predictor.rs`
- Added the unit test `test_transition_matrix_velocity_attitude_coupling` to verify the sign and correctness of the coupling block in the state transition matrix:
  ```rust
  #[test]
  fn test_transition_matrix_velocity_attitude_coupling() {
      // ... test setup ...
      let phi = predictor::compute_transition_matrix(&state, dt, &[imu_meas]);
      let expected_vel_att = nalgebra::Matrix3::new(
          0.0, 1.5, -1.0,
          -1.5, 0.0, 0.5,
          1.0, -0.5, 0.0,
      );
      // assert equality ...
  }
  ```

## Verification Details
1. **Red Phase (Buggy Code)**:
   - Command: `cargo test -p gneiss-rtk test_transition_matrix_velocity_attitude_coupling`
   - Output: Mismatch failure at `(0, 1)`, where expected was `1.5` but got `-1.5` due to the incorrect positive sign of the coupling matrix block.
2. **Green Phase (Fixed Code)**:
   - Command: `cargo test -p gneiss-rtk test_transition_matrix_velocity_attitude_coupling`
   - Output: Success (`test result: ok. 1 passed; 0 failed`).
3. **Workspace Verification**:
   - Command: `cargo test --workspace`
   - Output: All tests passed across all crates (including `gneiss-rtk` and integration tests).
