# RTK Verification and Testing Skill

This skill governs the development life cycle of GNSS/INS RTK systems. It ensures that critical EKF math (such as analytical Jacobians) and file parsers (such as RINEX/UBX/RTCM3) are thoroughly tested, preventing silent regressions and EKF divergence.

## 1. Numerical Jacobian Verification

Analytical Jacobians in EKFs (e.g., measurement models, motion constraints, and process transition matrices) are highly prone to:
* Sign flips (especially w.r.t. attitude error or coordinate transformations).
* Incorrect scaling (e.g., using degrees instead of radians, or neglecting frequency factors).
* Mismatched state indices.

### The Rule of Numerical Verification
Whenever you modify or add a Jacobian in the filter or constraints, **you MUST write a unit test that verifies the analytical Jacobian against a numerical Jacobian computed using central finite differences.**

#### Mathematical Definition
For a state vector $x$ and a measurement function $h(x)$, the analytical Jacobian is $H = \frac{\partial h}{\partial x}$.
The numerical Jacobian $H_{num}$ column $j$ is computed as:
$$
H_{num}[:, j] = \frac{h(x + \epsilon e_j) - h(x - \epsilon e_j)}{2 \epsilon}
$$
where $e_j$ is a unit vector along the $j$-th state dimension, and $\epsilon$ is a small perturbation (typically $10^{-6}$ to $10^{-8}$).

#### Rust Implementation Example
For attitude error states, the perturbation is applied as a small rotation:
```rust
fn compute_numerical_jacobian(
    state: &RtkState,
    h_func: impl Fn(&RtkState) -> DVector<f64>,
    epsilon: f64
) -> DMatrix<f64> {
    let base_meas = h_func(state);
    let m = base_meas.len();
    let n = state.covariance.nrows(); // total state size
    let mut h_num = DMatrix::zeros(m, n);

    for j in 0..n {
        let mut state_pos = state.clone();
        let mut state_neg = state.clone();
        
        // Perturb state element j
        if j >= 6 && j <= 8 {
            // Attitude error perturbation (small-angle rotation in ECEF)
            let mut dpsi_pos = Vector3::zeros();
            dpsi_pos[j - 6] = epsilon;
            let dq_pos = UnitQuaternion::from_scaled_axis(dpsi_pos);
            state_pos.attitude = dq_pos * state_pos.attitude;

            let mut dpsi_neg = Vector3::zeros();
            dpsi_neg[j - 6] = -epsilon;
            let dq_neg = UnitQuaternion::from_scaled_axis(dpsi_neg);
            state_neg.attitude = dq_neg * state_neg.attitude;
        } else if j >= 0 && j < 3 {
            // Position perturbation
            state_pos.position.vector[j] += epsilon;
            state_neg.position.vector[j] -= epsilon;
        } else if j >= 3 && j < 6 {
            // Velocity perturbation
            state_pos.velocity[j - 3] += epsilon;
            state_neg.velocity[j - 3] -= epsilon;
        } else {
            // Other scalar states (biases, ambiguities)
            // ...
        }

        let meas_pos = h_func(&state_pos);
        let meas_neg = h_func(&state_neg);
        let col = (meas_pos - meas_neg) / (2.0 * epsilon);
        h_num.set_column(j, &col);
    }
    h_num
}
```

Compare the analytical $H_{ana}$ with $H_{num}$ using an element-wise tolerance check:
```rust
assert!((H_ana - H_num).abs().max() < 1e-5, "Jacobian verification failed!");
```

---

## 2. Parser Edge-Case Testing

GNSS parsing bugs (e.g., ignoring LLI flags, parsing exponents incorrectly, or failing to handle missing frequency channels) directly lead to corrupted measurement models and EKF failures.

### The Rule of Parser Validation
Every GNSS parser (RINEX, UBX, RTCM) must have dedicated test coverage demonstrating:
1. **Flag Parsing**: Verification of cycle-slip indicators (LLI bits), lock-time counters, and quality flags.
2. **Malformed Input Safety**: Graceful error handling or exclusion of corrupt data lines instead of panic/crash.
3. **Empty Field Handling**: Correct identification of empty or missing observation values (e.g., whitespace spacing in RINEX, or missing fields in UBX).

---

## 3. Red-to-Green Bug Isolation

Before applying a fix for EKF divergence or parser failures:
1. **Reproduce the bug in a unit test first** (e.g., a test that asserts the incorrect Jacobian sign or the ignored field).
2. Verify the test fails (Red).
3. Apply the fix.
4. Verify the test passes (Green).
