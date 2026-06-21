# Handoff Report - Bug 2: Velocity-Attitude Transition Sign Mismatch

## 1. Observation
- In `crates/gneiss-rtk/src/engine/predictor.rs`, the attitude-to-velocity coupling block of the state transition matrix was computed at line 86:
  ```rust
  let vel_att = f_e_skew * dt;
  ```
- Added a regression test `test_transition_matrix_velocity_attitude_coupling` in `crates/gneiss-rtk/src/engine/tests_predictor.rs` at line 157.
- Running the new test on the unmodified codebase using `cargo test -p gneiss-rtk test_transition_matrix_velocity_attitude_coupling` failed with:
  ```
  thread 'engine::tests_predictor::tests::test_transition_matrix_velocity_attitude_coupling' panicked at crates/gneiss-rtk/src/engine/tests_predictor.rs:198:17:
  Mismatch at (0, 1): expected 1.5, got -1.5
  ```
- After negating the coupling block in `crates/gneiss-rtk/src/engine/predictor.rs`:
  ```rust
  let vel_att = -f_e_skew * dt;
  ```
- Running the same test passed successfully:
  ```
  running 1 test
  test engine::tests_predictor::tests::test_transition_matrix_velocity_attitude_coupling ... ok
  ```
- Running the workspace tests with `cargo test --workspace` succeeded with:
  ```
  test result: ok. 255 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.20s
  ...
  test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
  ```

## 2. Logic Chain
- **Observation 1**: The transition matrix equation uses positive `f_e_skew * dt` for the attitude-to-velocity coupling.
- **Observation 2**: Setting up a standard EKF state transition matrix with positive acceleration inputs produces an attitude-to-velocity coupling that has a positive sign mismatch relative to the theoretical analytical transition matrix $-[f_e \times] dt$.
- **Observation 3**: Creating a regression test with a known non-zero acceleration vector, an identity attitude, and zero accelerometer bias allows us to isolate this coupling block.
- **Observation 4**: The regression test fails on the original code, proving it acts as a robust test for this specific sign bug.
- **Observation 5**: Negating `f_e_skew` to get `let vel_att = -f_e_skew * dt;` correctly aligns the code with the theoretical EKF transition matrix.
- **Observation 6**: Running the regression test after applying the change results in a pass, and running workspace-wide tests verifies no regressions.

## 3. Caveats
- No caveats.

## 4. Conclusion
- The velocity-attitude transition sign mismatch has been resolved by negating the `vel_att` computation in `crates/gneiss-rtk/src/engine/predictor.rs`.
- The regression test `test_transition_matrix_velocity_attitude_coupling` in `crates/gneiss-rtk/src/engine/tests_predictor.rs` correctly verifies the fix and prevents future sign regressions.

## 5. Verification Method
1. Check the modified files:
   - `crates/gneiss-rtk/src/engine/predictor.rs`
   - `crates/gneiss-rtk/src/engine/tests_predictor.rs`
2. Run the specific regression test:
   - `cargo test -p gneiss-rtk test_transition_matrix_velocity_attitude_coupling`
3. Run all tests in the workspace to verify stability:
   - `cargo test --workspace`
