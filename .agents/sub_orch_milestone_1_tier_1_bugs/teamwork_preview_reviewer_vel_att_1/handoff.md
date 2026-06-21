# Handoff Report: Bug 2 Review (Velocity-Attitude Transition Sign Mismatch)

## 1. Observation
- File changes verified in `crates/gneiss-rtk/src/engine/predictor.rs` at line 86 (originally line 86-88):
  ```rust
-        let vel_att = f_e_skew * dt;
+        let vel_att = -f_e_skew * dt;
  ```
- File changes verified in `crates/gneiss-rtk/src/engine/tests_predictor.rs` where the regression test `test_transition_matrix_velocity_attitude_coupling` was added from line 159 to 206.
- The command `cargo test -p gneiss-rtk test_transition_matrix_velocity_attitude_coupling` was executed, returning:
  ```
running 1 test
test engine::tests_predictor::tests::test_transition_matrix_velocity_attitude_coupling ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 255 filtered out; finished in 0.01s
  ```
- The command `cargo test --workspace` was executed, returning:
  ```
test result: ok. 255 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.24s
...
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
  ```
- EKF attitude correction is implemented in `crates/gneiss-rtk/src/engine/updater.rs` line 32-41:
  ```rust
  fn apply_attitude_correction(state: &mut RtkState, dx: &DVector<f64>) {
      let d_theta = Vector3::new(dx[6], dx[7], dx[8]);
      if d_theta.norm() <= 1e-10 {
          return;
      }
      let dq =
          UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_normalize(d_theta), d_theta.norm());
      state.attitude = dq * state.attitude;
      state.attitude.renormalize();
  }
  ```

## 2. Logic Chain
- In `crates/gneiss-rtk/src/engine/updater.rs`, the EKF attitude correction applies an additive update in the ECEF frame as $\hat{C}_b^e = \text{Rot}(d\theta) C_b^e \approx (I + [d\theta \times]) C_b^e$.
- In the error-state formulation, the estimated ECEF acceleration is $\hat{f}^e = \hat{C}_b^e \hat{f}^b \approx (I + [\psi \times]) C_b^e (f^b + \delta f^b) \approx f^e + C_b^e \delta f^b + [\psi \times] f^e = f^e + C_b^e \delta f^b - [f^e \times] \psi$, where $\psi$ is the attitude error.
- The derivative of ECEF velocity error $\delta v$ with respect to attitude error $\psi$ is therefore $\frac{\partial \dot{v}}{\partial \psi} = - [f^e \times]$.
- In the state transition matrix $\Phi$, this coupling term over a step $\Delta t$ is approximated as $\Phi_{v \psi} = - [f^e \times] \Delta t$.
- In `crates/gneiss-rtk/src/engine/predictor.rs`, the code prior to the fix computed `vel_att = f_e_skew * dt`, which is $+ [f^e \times] \Delta t$.
- The fix correctly flipped the sign to `-f_e_skew * dt`, which corresponds to $- [f^e \times] \Delta t$, aligning the transition matrix with the EKF attitude error correction definition.
- A new test in `crates/gneiss-rtk/src/engine/tests_predictor.rs` verifies that the computed transition matrix coupling block matches the expected $- [f^e \times] \Delta t$ matrix for a known acceleration.
- Running `cargo test --workspace` validates that all tests in the workspace pass, with no regressions.

## 3. Caveats
- The state transition matrix is a first-order Taylor approximation ($\Phi \approx I + F \Delta t$). Under extreme acceleration changes or very large time steps, this approximation may degrade, but this is a standard limitation of linearized EKF propagation rather than an implementation bug.

## 4. Conclusion
- Final verdict: **APPROVE**.
- The fix for Bug 2 correctly resolves the sign mismatch in the velocity-attitude coupling term of the transition matrix. The code layout complies with the project conventions (tests are co-located in the `src/engine/tests_predictor.rs` module).

## 5. Verification Method
- Execute the specific unit test:
  `cargo test -p gneiss-rtk test_transition_matrix_velocity_attitude_coupling`
- Run the entire test suite to ensure no regressions:
  `cargo test --workspace`
