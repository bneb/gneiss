# Handoff Report — Review of Bug 2: Velocity-Attitude Transition Sign Mismatch

## 1. Observation
The worker modified the following files in the project workspace `/Users/kevin/projects/gneiss`:
* **`crates/gneiss-rtk/src/engine/predictor.rs`**:
  * Line 86: Replaced `let vel_att = f_e_skew * dt;` with `let vel_att = -f_e_skew * dt;`.
* **`crates/gneiss-rtk/src/engine/tests_predictor.rs`**:
  * Added the test `test_transition_matrix_velocity_attitude_coupling` (lines 159-205) to verify the transition matrix's velocity-attitude coupling block against analytical calculations.

Test execution in `/Users/kevin/projects/gneiss`:
* `cargo test -p gneiss-rtk test_transition_matrix_velocity_attitude_coupling` returns:
  `test engine::tests_predictor::tests::test_transition_matrix_velocity_attitude_coupling ... ok`
* `cargo test --workspace` returns:
  `test result: ok. 255 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.21s`

---

## 2. Logic Chain
1. Under a standard left-multiplicative ECEF-frame perturbation model, the estimated attitude rotation matrix is represented as:
   $$\hat{R}_b^e = (I + [\delta \theta \times]) R_b^e$$
   where $\delta \theta$ is the attitude error vector in the ECEF frame.
2. The true and estimated specific force projections in ECEF are $f^e = R_b^e f^b$ and $\hat{f}^e = \hat{R}_b^e \hat{f}^b \approx (I + [\delta \theta \times]) R_b^e (f^b - \delta f_b^b)$ respectively (neglecting second-order error terms).
3. The derivative of the velocity error vector $\delta v^e = \hat{v}^e - v^e$ includes the specific force coupling term:
   $$\dot{\delta v}^e = \hat{f}^e - f^e \approx [\delta \theta \times] f^e - R_b^e \delta f_b^b = \delta \theta \times f^e - R_b^e \delta f_b^b$$
4. Applying the cross-product skew-symmetric property ($\delta \theta \times f^e = - f^e \times \delta \theta = - [f^e\times] \delta \theta$) yields:
   $$\dot{\delta v}^e \supset - [f^e\times] \delta \theta$$
5. Over a time step $dt$, the first-order Euler approximation of the state transition matrix ($\Phi$) block mapping the attitude error state $\delta \theta$ to the velocity error state $\delta v^e$ is given by:
   $$\Phi_{v, \theta} = - [f^e\times] dt$$
6. The original implementation used $\Phi_{v, \theta} = + [f^e\times] dt$. The worker corrected this to $- [f^e\times] dt$ via `let vel_att = -f_e_skew * dt;`.
7. The new unit test validates the corrected sign by asserting that $\Phi_{v, \theta}$ matches the matrix element-by-element for a non-zero specific force vector.

---

## 3. Caveats
* **First-Order Propagation**: The transition matrix uses a first-order Euler step. In high-attitude-rate environments or for large $dt$, higher-order integration terms might be necessary.
* **Single specific force sample**: The transition matrix evaluates the specific force $f^e$ using only the last IMU measurement in the prediction buffer. For large time steps with varying acceleration, this can introduce integration errors.

---

## 4. Conclusion
The worker's fix correctly resolves the sign mismatch in the velocity-attitude transition block. The co-located unit test successfully reproduces and verifies the correct mathematical coupling. Verdict is **APPROVE**.

---

## 5. Verification Method
Run the following commands in the workspace `/Users/kevin/projects/gneiss`:
```bash
cargo test -p gneiss-rtk test_transition_matrix_velocity_attitude_coupling
cargo test --workspace
```
Check that the output compiles cleanly and all tests pass.

---

## Quality Review Report

### Review Summary
* **Verdict**: APPROVE

### Findings
* None. The implementation is correct, conforms to Rust standards, and is properly integrated.

### Verified Claims
* `let vel_att = -f_e_skew * dt;` correctly models the negative skew coupling $\rightarrow$ Verified via analytical derivation and the unit test `test_transition_matrix_velocity_attitude_coupling` $\rightarrow$ PASS.
* Code complies with layout rules $\rightarrow$ Verified that tests are co-located in the library crate under `crates/gneiss-rtk/src/engine/tests_predictor.rs` and no source/test code is written to `.agents/` $\rightarrow$ PASS.

### Coverage Gaps
* None.

### Unverified Items
* None.

---

## Challenge (Adversarial Review) Report

### Challenge Summary
* **Overall risk assessment**: LOW

### Challenges

#### [Low] Challenge 1: Discretization Error under Changing Dynamics
* **Assumption challenged**: The specific force $f^e$ is constant over the prediction interval $dt$ and can be approximated by the last IMU measurement in the buffer.
* **Attack scenario**: High-frequency vibration or rapid acceleration changes (e.g. vehicle collision or rough terrain traversal) where the specific force changes dramatically within the $dt$ window.
* **Blast radius**: Suboptimal covariance propagation in the EKF, leading to minor transient estimation errors.
* **Mitigation**: For typical GNSS/INS EKF systems running at IMU rate, $dt$ is very small (10ms), rendering this discretization error negligible.

#### [Medium] Challenge 2: Small-Angle Assumption Failure
* **Assumption challenged**: The attitude error $\delta \theta$ is small, ensuring that the linear approximation of the transition matrix holds.
* **Attack scenario**: Filter startup or loss of GNSS signal under high dynamics, where the attitude uncertainty is large (e.g. $> 10$ degrees).
* **Blast radius**: The state transition matrix linear prediction may diverge or cause filter instability due to unmodeled non-linear coupling.
* **Mitigation**: Large attitude uncertainties are represented by large values in the covariance matrix $P$. EKF measurement updates will dominate and correct the state, but users should be aware that the filter relies on quick convergence to a small-angle regime.

### Stress Test Results
* **Scenario**: Stationary state with gravity cancelation $\rightarrow$ Stable velocity and position $\rightarrow$ PASS.
* **Scenario**: IMU rotation prediction $\rightarrow$ Correct yaw tracking $\rightarrow$ PASS.
* **Scenario**: Lever-arm coupling test $\rightarrow$ Corrections propagate properly $\rightarrow$ PASS.

### Unchallenged Areas
* **Area**: WGS84 Gravity formula implementation detail $\rightarrow$ Reason: Out of scope for this specific velocity-attitude transition bug review.
