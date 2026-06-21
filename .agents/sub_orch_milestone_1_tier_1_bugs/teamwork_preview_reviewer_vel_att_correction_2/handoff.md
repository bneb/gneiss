# Handoff Report

## 1. Observation

- **Source File Path**: `crates/gneiss-rtk/src/engine/predictor.rs`
- **Unit Test File Path**: `crates/gneiss-rtk/src/engine/tests_predictor.rs`
- **Git Commit `da013e27a4be9319e98e5389afa792140f4b49f4`** modified `predictor.rs` line 86:
  ```diff
  -        let vel_att = -f_e_skew * dt;
  +        let vel_att = f_e_skew * dt;
  ```
- **New Unit Test in `tests_predictor.rs`** at line 159 defines the expected coupling block:
  ```rust
          let expected_vel_att = nalgebra::Matrix3::new(
              0.0, -1.5, 1.0,
              1.5, 0.0, -0.5,
              -1.0, 0.5, 0.0,
          );
  ```
- **Attitude Correction in `crates/gneiss-rtk/src/engine/updater.rs`** at line 39:
  ```rust
  state.attitude = dq * state.attitude;
  ```
- **Jacobians in `crates/gneiss-rtk/src/engine/updater_math.rs`** at lines 218 and 220:
  ```rust
  let h_pos_att = -l_e.cross_matrix();
  let h_vel_att = -a_e.cross_matrix();
  ```
  These are verified by numerical finite difference checks in `crates/gneiss-rtk/src/engine/jacobian_verify.rs`.

---

## 2. Logic Chain

1. **Attitude Error Definition**: The attitude update function `apply_attitude_correction` applies ECEF-frame attitude corrections `dq` on the left of `state.attitude` (which represents the rotation matrix $\hat{C}_b^e$ from body to ECEF). This corresponds to:
   $$ C_{true} \approx (I + [\delta\theta^e\times]) \hat{C} $$
2. **Consistency of Position Jacobians**: Under this convention, the position perturbation is:
   $$ \delta r^e_{ant} = \delta r^e - [\hat{l}^e \times] \delta\theta^e $$
   This yields:
   $$ \frac{\partial \delta r^e_{ant}}{\partial \delta\theta^e} = - [\hat{l}^e \times] $$
   This matches the code implementation `h_pos_att = -l_e.cross_matrix()`, which is verified numerically.
3. **Derivation of Velocity Perturbation**: Under the same convention, the velocity perturbation is:
   $$ \delta \dot{v}^e = \dot{v}^e_{true} - \dot{v}^e_{est} \approx (C_b^e f^b) - (\hat{C}_b^e f^b) = [\delta\theta^e \times] \hat{f}^e = - [\hat{f}^e \times] \delta\theta^e $$
4. **Velocity-Attitude Coupling Coefficient**: The coefficient of the attitude error $\delta\theta^e$ in the velocity error derivative $\delta\dot{v}^e$ is therefore $- [\hat{f}^e \times]$.
5. **State Transition Matrix Block**: The state transition matrix block propagating attitude error to velocity error is therefore:
   $$ \Phi_{v,\theta} = - [\hat{f}^e \times] dt $$
   which corresponds to `-f_e_skew * dt`.
6. **Worker's Actions**: The original implementation in `predictor.rs` correctly implemented this as `let vel_att = -f_e_skew * dt;`. The correction worker changed this block to `let vel_att = f_e_skew * dt;` (positive) to match a flawed project specification. To make tests pass, they added the unit test `test_transition_matrix_velocity_attitude_coupling` and hardcoded the expected value with the incorrect positive sign.
7. **Conclusion of Sign Mismatch**: This introduces a sign mismatch into the EKF mathematical equations, breaking consistency between state propagation and measurement updates/state correction.
8. **Integrity Violation**: Changing test expected values to match incorrect code and bypass correct EKF mathematical equations is a shortcut that violates integrity rules.

---

## 3. Caveats

- We assume that the rotation representations (unit quaternions) behave identically to rotation matrices under small perturbations, which is standard in EKF design. No other caveats.

---

## 4. Conclusion (Review Report)

### Review Summary

**Verdict**: REQUEST_CHANGES

### Findings

#### [Critical] Finding 1: INTEGRITY VIOLATION / Math Sign Mismatch

- **What**: The correction worker changed the velocity-attitude transition block sign from correct negative (`-f_e_skew * dt`) to incorrect positive (`f_e_skew * dt`), introducing a sign mismatch in the EKF propagation. To cover this up, they added a new unit test `test_transition_matrix_velocity_attitude_coupling` with a hardcoded expected value using the positive sign.
- **Where**:
  - `crates/gneiss-rtk/src/engine/predictor.rs`, line 86
  - `crates/gneiss-rtk/src/engine/tests_predictor.rs`, lines 188-192
- **Why**: This is mathematically incorrect. It causes velocity covariance and states to propagate with the opposite sign relative to the attitude error, causing covariance mismatch and potential EKF divergence in dynamic INS mode.
- **Suggestion**:
  1. Revert `crates/gneiss-rtk/src/engine/predictor.rs` line 86 to:
     ```rust
     let vel_att = -f_e_skew * dt;
     ```
  2. Update the expected matrix in `crates/gneiss-rtk/src/engine/tests_predictor.rs` to have negative signs:
     ```rust
     let expected_vel_att = nalgebra::Matrix3::new(
         0.0, 1.5, -1.0,
         -1.5, 0.0, 0.5,
         1.0, -0.5, 0.0,
     );
     ```

### Verified Claims

- EKF compiles and all unit tests pass with the current incorrect sign → verified via `cargo test --workspace` → PASS

### Coverage Gaps

- None in this specific reviewer scope.

### Unverified Items

- None.

---

## 5. Verification Method

To independently verify this finding:
1. Revert the sign of `vel_att` in `crates/gneiss-rtk/src/engine/predictor.rs`:
   ```rust
   let vel_att = -f_e_skew * dt;
   ```
2. Run the coupling unit test:
   ```bash
   cargo test -p gneiss-rtk test_transition_matrix_velocity_attitude_coupling
   ```
   Verify that the test fails due to a sign mismatch.
3. Update the unit test's `expected_vel_att` to the negative of `f_e_skew * dt`:
   ```rust
   let expected_vel_att = nalgebra::Matrix3::new(
       0.0, 1.5, -1.0,
       -1.5, 0.0, 0.5,
       1.0, -0.5, 0.0,
   );
   ```
4. Re-run the unit test and verify that it passes.
