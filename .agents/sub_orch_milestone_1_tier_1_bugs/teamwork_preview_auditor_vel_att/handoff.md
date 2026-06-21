# Handoff Report — Bug 2 Fix Forensic Audit

This handoff report documents the forensic integrity audit of the fix for **Bug 2: Velocity-Attitude Transition Sign Mismatch** in the `gneiss` GNSS/INS EKF system.

---

## 1. Observation

We directly observed the following modifications and test executions in the codebase:

### Code Modifications
1. **In `crates/gneiss-rtk/src/engine/predictor.rs` (lines 83-91):**
   ```rust
   for i in 0..3 {
       phi[(i, 3 + i)] = dt;
   }

   let vel_att = -f_e_skew * dt;
   for r in 0..3 {
       for c in 0..3 {
           phi[(3 + r, 6 + c)] = vel_att[(r, c)];
       }
   }
   ```
   *Source Reference*: The working tree diff reveals `let vel_att = f_e_skew * dt;` was replaced by `let vel_att = -f_e_skew * dt;`.

2. **In `crates/gneiss-rtk/src/engine/tests_predictor.rs` (lines 159-205):**
   A new unit test `test_transition_matrix_velocity_attitude_coupling` was added:
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

       // Compute the expected coupling block: -skew(f_e) * dt
       // Since state.attitude is identity and accel_bias is zero, f_e = accel = [1.0, 2.0, 3.0]
       // skew(f_e) is:
       // [ 0.0, -3.0,  2.0]
       // [ 3.0,  0.0, -1.0]
       // [-2.0,  1.0,  0.0]
       //
       // Negating it and multiplying by dt (0.5):
       // [ 0.0,  1.5, -1.0]
       // [-1.5,  0.0,  0.5]
       // [ 1.0, -0.5,  0.0]
       let expected_vel_att = nalgebra::Matrix3::new(
           0.0, 1.5, -1.0,
           -1.5, 0.0, 0.5,
           1.0, -0.5, 0.0,
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

### Test Execution Results
Running `cargo test --package gneiss-rtk` outputs:
```
test engine::tests_predictor::tests::test_transition_matrix_velocity_attitude_coupling ... ok
test engine::tests_predictor::tests::test_coupling_lever_arm_to_attitude ... ok
test engine::tests_predictor::tests::test_physics_centrifugal_cancellation ... ok
test engine::tests_predictor::tests::test_physics_stationary_gravity ... ok
test engine::tests_predictor::tests::test_imu_prediction_rotation ... ok
```
We also observed one unrelated failure in `engine::ppp_iekf::mutant_killer_tests::test_resolve_widelane_ar_insufficient` due to a different bug under active development.

---

## 2. Logic Chain

1. **Perturbation Framework**: In `crates/gneiss-rtk/src/engine/updater.rs`, `apply_attitude_correction` applies corrections via a left-multiplied ECEF perturbation:
   $$R_{true} = R_{pert} R_{est} = \exp([\delta \psi \times]) R_{est} \approx (I + [\delta \psi \times]) R_{est}$$
   This choice uses a **positive** sign ($+\delta \psi$) for the attitude error rotation vector.

2. **Specific Force Dynamics**:
   The ECEF specific force is $f_e = R f_b$. Perturbing it gives:
   $$f_{e, true} = R_{true} f_b \approx (I + [\delta \psi \times]) f_{e, est} = f_{e, est} + \delta \psi \times f_{e, est}$$
   Using the cross product's anti-symmetric property:
   $$\delta \psi \times f_{e, est} = - (f_{e, est} \times \delta \psi) = - [f_{e, est} \times] \delta \psi$$
   Therefore:
   $$\frac{\partial f_e}{\partial \delta \psi} = - [f_e \times] = - \text{f\_e\_skew}$$

3. **Transition Matrix Coupling**:
   The velocity dynamics propagate attitude corrections over time interval $\Delta t$ into velocity corrections:
   $$\delta v_{k+1} \approx \delta v_k + \frac{\partial f_e}{\partial \delta \psi} \Delta t \delta \psi_k = \delta v_k - [f_e \times] \Delta t \delta \psi_k$$
   Thus, the velocity-attitude transition coupling block $\Phi_{v,\psi}$ must equal $- f_e\_skew \times \Delta t$.

4. **Code Correctness**:
   * The codebase originally had $+ f_e\_skew \times \Delta t$, which was mathematically incorrect.
   * The fix correctly introduces the negative sign (`-f_e_skew * dt`) and aligns it with the rest of the EKF's left-multiplied ECEF attitude error definition.
   * The newly added test `test_transition_matrix_velocity_attitude_coupling` enforces the correct sign, and would fail if the sign were flipped.

---

## 3. Caveats

* **Unrelated Test Failures**: There is a failing test `engine::ppp_iekf::mutant_killer_tests::test_resolve_widelane_ar_insufficient` in the workspace, which is related to a different PPP ambiguity resolution bug and does not impact or invalidate this predictor sign correction.
* No other caveats.

---

## 4. Conclusion

### Forensic Audit Report

**Work Product**: `crates/gneiss-rtk/src/engine/predictor.rs` and `crates/gneiss-rtk/src/engine/tests_predictor.rs`  
**Profile**: General Project  
**Verdict**: **CLEAN**

#### Phase Results
- **Hardcoded output detection**: PASS — No hardcoded mock bypasses or hardcoded test returns were found in the implementation.
- **Facade detection**: PASS — The implementation of `compute_transition_matrix` dynamically calculates and constructs the transition matrix using actual IMU data and attitude states.
- **Pre-populated artifact detection**: PASS — No pre-populated test result logs or fake attestation files exist in the workspace for this bug.
- **Build and run**: PASS — The `gneiss-rtk` package successfully compiles, and the target predictor tests pass cleanly.
- **Dependency audit**: PASS — Core mechanization and transition matrix logic are built directly within the crate and do not delegate work to unauthorized third-party libraries.

---

### Adversarial Review / Challenge Report

#### Challenge Summary
**Overall risk assessment**: **LOW**

#### Challenges

##### [Low] Challenge 1: Local vs. Global Perturbation Definition
* **Assumption challenged**: Whether the attitude error state in this EKF is defined globally or locally.
* **Attack scenario**: If the error state was defined in the local/body frame, the transition Jacobian would be $-[f_e \times] R_b^e$.
* **Blast radius**: If this assumption were wrong, the transition matrix would lack the $R_b^e$ rotation, leading to EKF divergence.
* **Mitigation**: Confirmed via `apply_attitude_correction` in `updater.rs` that the correction is left-multiplied (`state.attitude = dq * state.attitude`), proving the error is globally (ECEF) defined.

---

## 5. Verification Method

To verify the audit findings and fix correctness independently:

1. **Verify the transition sign in source code**:
   Inspect line 86 of `crates/gneiss-rtk/src/engine/predictor.rs`:
   ```rust
   let vel_att = -f_e_skew * dt;
   ```
2. **Execute the predictor tests**:
   Run the following terminal command in the project root directory:
   ```bash
   cargo test --package gneiss-rtk --lib -- engine::tests_predictor::tests
   ```
   All five tests under `engine::tests_predictor::tests` must compile and pass successfully.
