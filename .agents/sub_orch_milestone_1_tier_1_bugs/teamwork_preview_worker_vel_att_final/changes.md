# Changes made

## 1. Velocity-Attitude Transition Sign Mismatch Fix
- **File**: `crates/gneiss-rtk/src/engine/predictor.rs`
- **Details**: Changed `let vel_att = f_e_skew * dt;` to `let vel_att = -f_e_skew * dt;` at line 86. This aligns the transition matrix block with the negative skew symmetric coupling of the specific force.

## 2. Unit Test Update
- **File**: `crates/gneiss-rtk/src/engine/tests_predictor.rs`
- **Details**: Updated `test_transition_matrix_velocity_attitude_coupling` to assert the correct negative coupling:
  ```rust
  let expected_vel_att = nalgebra::Matrix3::new(
      0.0, 1.5, -1.0,
      -1.5, 0.0, 0.5,
      1.0, -0.5, 0.0,
  );
  ```

## 3. Formatting
- **File**: `crates/gneiss-rtk/src/engine/tests_predictor.rs` and `crates/gneiss-rtk/src/engine/predictor.rs`
- **Command**: `cargo fmt` and `cargo fmt --check`
- **Details**: Formatted modified files and verified that formatting is compliant across the workspace with no violations.

## Verification Results
- **Command**: `cargo test -p gneiss-rtk`
  - **Result**: Passed (255 tests passed, 0 failed, 1 ignored)
- **Command**: `cargo test --workspace`
  - **Result**: Passed (all tests passed)
