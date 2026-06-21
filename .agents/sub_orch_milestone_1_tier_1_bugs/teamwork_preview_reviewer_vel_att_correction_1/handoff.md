# Handoff Report — Bug 2 Velocity-Attitude Transition Sign Mismatch Verification

## 1. Observation

We directly observed the following files and tool execution results in the project directory `/Users/kevin/projects/gneiss`:

1. **Source Code Implementation** (`crates/gneiss-rtk/src/engine/predictor.rs` lines 86-91):
   ```rust
   let vel_att = f_e_skew * dt;
   for r in 0..3 {
       for c in 0..3 {
           phi[(3 + r, 6 + c)] = vel_att[(r, c)];
       }
   }
   ```
   Positive sign coupling `let vel_att = f_e_skew * dt;` is implemented.

2. **Unit Test Implementation** (`crates/gneiss-rtk/src/engine/tests_predictor.rs` lines 159-205):
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
   This unit test explicitly asserts that the transition matrix's velocity-attitude coupling matches the manual derivation of `f_e_skew * dt`.

3. **Test Executions**:
   - `cargo test --package gneiss-rtk --lib -- engine::tests_predictor`
     - Command completed successfully with output:
       ```
       running 5 tests
       test engine::tests_predictor::tests::test_coupling_lever_arm_to_attitude ... ok
       test engine::tests_predictor::tests::test_imu_prediction_rotation ... ok
       test engine::tests_predictor::tests::test_transition_matrix_velocity_attitude_coupling ... ok
       test engine::tests_predictor::tests::test_physics_centrifugal_cancellation ... ok
       test engine::tests_predictor::tests::test_physics_stationary_gravity ... ok
       test result: ok. 5 passed; 0 failed; 0 ignored; 0 measured; 251 filtered out
       ```
   - `cargo test --workspace`
     - Command completed successfully: `255 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.26s` for `gneiss-rtk` tests.

4. **Code Formatting Check**:
   - `cargo fmt --check`
     - Command failed with exit code 1 due to formatting layout differences in `crates/gneiss-rtk/src/engine/tests_predictor.rs`:
       ```
       Diff in /Users/kevin/projects/gneiss/crates/gneiss-rtk/src/engine/tests_predictor.rs:185:
       -        let expected_vel_att = nalgebra::Matrix3::new(
       -            0.0, -1.5, 1.0,
       -            1.5, 0.0, -0.5,
       -            -1.0, 0.5, 0.0,
       -        );
       +        let expected_vel_att =
       +            nalgebra::Matrix3::new(0.0, -1.5, 1.0, 1.5, 0.0, -0.5, -1.0, 0.5, 0.0);
       ```

5. **Clippy Static Analysis**:
   - `cargo clippy --all-targets`
     - Command completed successfully with zero clippy errors (only code warnings unrelated to this issue).

---

## 2. Logic Chain

1. **Bug Correction Identification**: In `crates/gneiss-rtk/src/engine/predictor.rs`, the velocity-attitude transition block `vel_att` is calculated as `f_e_skew * dt`. Comparing this to the previous commit `da013e27a4be9319e98e5389afa792140f4b49f4` history reveals that the line was indeed updated from `let vel_att = -f_e_skew * dt;` to `let vel_att = f_e_skew * dt;`, implementing the positive sign coupling as requested.
2. **Behavioral Test Verification**: The unit test `test_transition_matrix_velocity_attitude_coupling` in `crates/gneiss-rtk/src/engine/tests_predictor.rs` manually constructs the expected matrix block using the positive coupling formulation (`f_e_skew * dt`) and compares it against the calculated transition matrix `phi`. Since `cargo test` succeeds for this test, the implementation correctly computes a positive coupling sign, and any future sign change would fail this regression test.
3. **Workspace Integrity**: Running `cargo test --workspace` passes all tests without regressions, meaning the sign correction is compatible with the entire navigation filter suite.
4. **Code Layout and Formatting**: Running `cargo fmt --check` fails because `expected_vel_att` matrix initialization in the unit test `test_transition_matrix_velocity_attitude_coupling` is formatted over multiple lines instead of a single line, violating the project formatting guidelines. Because the formatting check fails, the verdict must be `REQUEST_CHANGES`.

---

## 3. Caveats

* **Local vs Global Attitude Definition Controversy**: The EKF attitude error is left-multiplied in the global ECEF frame (`R_true = (I + [δψ×]) R_est`). In this framework, the actual perturbation of the specific force ECEF vector is $R_{true} f^b \approx f^e - [f^e\times] \delta\psi$. This would theoretically result in a transition matrix coupling term of $-f^e_{skew} dt$. However, the user request specifically mandates ensuring positive sign coupling `vel_att = f_e_skew * dt;` is implemented, which we verified is the case. We accept this positive sign constraint as it is explicitly requested by the mission.
* No other caveats.

---

## 4. Conclusion

### Review Summary

**Verdict**: **REQUEST_CHANGES**

### Findings

#### [Minor] Finding 1: Code Formatting Non-Conformance

- **What**: `cargo fmt --check` fails due to formatting issues in `tests_predictor.rs`.
- **Where**: `crates/gneiss-rtk/src/engine/tests_predictor.rs` at line 185-189.
- **Why**: Matrix3 definition is not formatted using standard rustfmt layout, breaking CI/CD verification.
- **Suggestion**: Format the file using `cargo fmt` to align with the workspace-wide formatting rules.

### Verified Claims

- Positive sign coupling `vel_att = f_e_skew * dt;` is implemented → verified via `view_file` on `predictor.rs` → PASS
- Unit test `test_transition_matrix_velocity_attitude_coupling` asserts positive sign coupling → verified via `view_file` on `tests_predictor.rs` → PASS
- Build and tests pass cleanly → verified via `cargo test` command execution → PASS

### Coverage Gaps

- None. The unit test covers the exact state transition matrix components changed.

---

### Challenge Summary

**Overall risk assessment**: **LOW**

### Challenges

#### [Low] Challenge 1: Code Formatting Breaking CI

- **Assumption challenged**: That the code is ready for merge.
- **Attack scenario**: Pull request merged with formatting violations, triggering CI build failures.
- **Blast radius**: Prevents developers from merging other changes or block pipeline.
- **Mitigation**: Requiring formatting fix before approval (`REQUEST_CHANGES`).

---

## 5. Verification Method

To independently verify the review findings:

1. **Check Transition Coupling Sign**:
   Open `crates/gneiss-rtk/src/engine/predictor.rs` and verify line 86 reads:
   ```rust
   let vel_att = f_e_skew * dt;
   ```

2. **Verify Formatting Violations**:
   Run the following command in the project root:
   ```bash
   cargo fmt --check
   ```
   Verify that it reports a formatting diff in `crates/gneiss-rtk/src/engine/tests_predictor.rs`.

3. **Verify Tests**:
   Run the test suite:
   ```bash
   cargo test --package gneiss-rtk --lib -- engine::tests_predictor::tests::test_transition_matrix_velocity_attitude_coupling
   ```
