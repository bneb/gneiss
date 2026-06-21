# Handoff Report: Bug 2 — Velocity-Attitude Transition Sign Mismatch

## 1. Observation
In `crates/gneiss-rtk/src/engine/predictor.rs`, lines 86–91:
```rust
        let vel_att = f_e_skew * dt;
        for r in 0..3 {
            for c in 0..3 {
                phi[(3 + r, 6 + c)] = vel_att[(r, c)];
            }
        }
```
where the variables `f_e` and `f_e_skew` are defined in the same file at lines 78–79:
```rust
        let f_e = state.attitude * (imu_buffer.last().unwrap().accel - state.accel_bias);
        let f_e_skew = skew_symmetric(&f_e);
```

In `crates/gneiss-rtk/src/engine/updater.rs`, lines 32–41:
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
1. In `updater.rs`, attitude corrections are applied via left-multiplying the attitude quaternion:
   $$R_b^e \leftarrow \delta R \cdot R_b^e$$
   where $\delta R$ represents a small rotation vector $\psi^e$ in the ECEF (global) frame.
2. For small error angles, this representation yields:
   $$R_b^e \approx (I + [\psi^e \times]) \hat{R}_b^e$$
3. In `predictor.rs`, the specific force projection into the ECEF frame is:
   $$f^e = R_b^e f^b \approx (I + [\psi^e \times]) \hat{R}_b^e (\hat{f}^b - \delta b_a)$$
   Ignoring second-order terms:
   $$f^e \approx \hat{f}^e - \hat{R}_b^e \delta b_a + [\psi^e \times] \hat{f}^e$$
   where $\hat{f}^e = \hat{R}_b^e \hat{f}^b$.
4. Under the anti-commutative property of the cross product, we have:
   $$[\psi^e \times] \hat{f}^e = - [\hat{f}^e \times] \psi^e$$
5. Consequently, the velocity error dynamics is:
   $$\dot{\delta v}^e \approx - [\hat{f}^e \times] \psi^e - \hat{R}_b^e \delta b_a - 2 [\omega_{ie}^e \times] \delta v^e$$
6. Therefore, the transition matrix block mapping attitude error $\psi^e$ to velocity error $\delta v^e$ (at row index 3..6 and column index 6..9) must be:
   $$\Phi_{v, \psi} = - [\hat{f}^e \times] dt$$
7. In `predictor.rs`, `vel_att` is computed as:
   `let vel_att = f_e_skew * dt;`
   Since `f_e_skew` is defined as the skew-symmetric matrix of the estimated ECEF specific force $[\hat{f}^e \times]$, this evaluates to $+ [\hat{f}^e \times] dt$. This is an exact sign mismatch (positive instead of negative).

## 3. Caveats
No caveats. The coordinate systems and error definitions are consistent across the filter, and the sign inversion is mathematically verified.

## 4. Conclusion
There is a sign mismatch in `crates/gneiss-rtk/src/engine/predictor.rs` at line 86.
The transition matrix term for the attitude-to-velocity coupling block `vel_att` must be negated.

### Proposed Code Change:
In `crates/gneiss-rtk/src/engine/predictor.rs` at line 86:
```rust
<<<<
        let vel_att = f_e_skew * dt;
====
        let vel_att = -f_e_skew * dt;
>>>>
```

## 5. Verification Method

### 1. Test Execution Command
Run the existing tests to ensure compilation and basic prediction stability:
```bash
cargo test -p gneiss-rtk
```

### 2. New Verification Test Case
To protect against regressions and verify this fix, add the following test case to `crates/gneiss-rtk/src/engine/tests_predictor.rs`:

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

        // Specific force pointing down in ECEF/body frame
        let accel = Vector3::new(0.0, 0.0, -9.81);
        let imu_meas = ImuMeasurement::new(0, accel, Vector3::zeros());

        let dt = 0.1;
        let phi = predictor::compute_transition_matrix(&state, dt, &[imu_meas]);

        // Specific force f_e = [0.0, 0.0, -9.81].
        // For a small rotation delta_theta around Y-axis (psi_y):
        // R_b^e = [  1   0   delta_theta ]
        //         [  0   1        0      ]
        //         [-delta_theta   0   1  ]
        //
        // f_e_true = R_b^e * f_b = [-9.81 * delta_theta, 0.0, -9.81]^T
        // f_e_est = [0.0, 0.0, -9.81]^T
        // v_dot_error = f_e_true - f_e_est = [-9.81 * delta_theta, 0.0, 0.0]^T
        //
        // Thus, the velocity error in X changes by -9.81 * delta_theta * dt.
        // The transition matrix element phi[(3, 7)] (mapping psi_y to delta_v_x)
        // must be -9.81 * dt = -0.981.
        let expected_coupling = -9.81 * dt;
        let computed_coupling = phi[(3, 7)];

        assert!(
            (computed_coupling - expected_coupling).abs() < 1e-10,
            "Sign mismatch in velocity-attitude coupling: expected {}, got {}",
            expected_coupling,
            computed_coupling
        );
    }
```
*Note: Under the current incorrect implementation, this test fails because `computed_coupling` equals `+0.981` instead of `-0.981`.*
