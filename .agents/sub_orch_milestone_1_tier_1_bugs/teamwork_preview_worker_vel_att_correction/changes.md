# Changes — Velocity-Attitude Transition Sign Mismatch Correction

## File: crates/gneiss-rtk/src/engine/tests_predictor.rs

- Modified `test_transition_matrix_velocity_attitude_coupling` to assert positive coupling sign (`skew(f_e) * dt`) instead of the incorrect negative sign.
- Match updated matrix values:
  - Row 0: `0.0, -1.5, 1.0`
  - Row 1: `1.5, 0.0, -0.5`
  - Row 2: `-1.0, 0.5, 0.0`

## File: crates/gneiss-rtk/src/engine/predictor.rs

- Verified that `crates/gneiss-rtk/src/engine/predictor.rs` at line 86 correctly uses the positive sign:
  ```rust
  let vel_att = f_e_skew * dt;
  ```
