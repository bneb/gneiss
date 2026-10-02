#[cfg(test)]
use super::*;
use std::collections::BTreeMap;
use crate::swfg::variables::{VariableId, VariableKind, VariableNode};

fn build_test_values(
pose_val: &[f64; 6],
vel_val: &[f64; 3],
) -> (VariableId, VariableId, VariableValues) {
let mut vars = BTreeMap::new();
let id_pose = VariableId::new(1);
let id_vel = VariableId::new(2);

let mut node_pose = VariableNode::new(id_pose, VariableKind::Pose { epoch: 1 });
node_pose.set_value(pose_val);
vars.insert(id_pose, node_pose);

let mut node_vel = VariableNode::new(id_vel, VariableKind::Velocity { epoch: 1 });
node_vel.set_value(vel_val);
vars.insert(id_vel, node_vel);

let values = VariableValues::build(&vars);
(id_pose, id_vel, values)
}

/// Reproduce `SlidingWindowSolver::apply_delta` (swfg/solver/mod.rs:336) for a
/// pose: the attitude increment is applied on the RIGHT, `q_new = q_old *
/// exp(delta)`, and the rotation vector is re-extracted with `scaled_axis()`.
fn pose_after_solver_delta(pose: &[f64; 6], delta: Vector3<f64>) -> [f64; 6] {
let q_old = UnitQuaternion::from_scaled_axis(Vector3::new(pose[3], pose[4], pose[5]));
let rot = (q_old * UnitQuaternion::from_scaled_axis(delta)).scaled_axis();
[pose[0], pose[1], pose[2], rot.x, rot.y, rot.z]
}

// --------------------------------------------------------------- Odometer

#[test]
fn test_odometer_velocity_factor_residual_zero_at_truth() {
let pose = [100.0, 200.0, 300.0, 0.0, 0.0, 0.0];
let vel = [10.0, 0.0, 0.0];
let (id_pose, id_vel, values) = build_test_values(&pose, &vel);

let factor = OdometerVelocityFactor::new(
id_pose,
id_vel,
Vector3::new(10.0, 0.0, 0.0),
Vector3::new(0.01, 0.01, 0.01),
);

let r = factor.residual(&values);
assert!(r.norm() < 1e-12);
assert_eq!(factor.variables(), &[id_pose, id_vel]);
assert_eq!(factor.robust_threshold(), Some(3.0));
}

#[test]
fn test_odometer_velocity_factor_numerical_jacobian() {
let pose = [100.0, 200.0, 300.0, 0.1, -0.2, 0.3];
let vel = [12.0, -3.0, 5.0];
let (id_pose, id_vel, values) = build_test_values(&pose, &vel);

let factor = OdometerVelocityFactor::new(
id_pose,
id_vel,
Vector3::new(10.0, 0.0, 0.0),
Vector3::new(0.01, 0.01, 0.01),
);

let j_analytic = factor.jacobian(&values);
let eps = 1e-7;

for comp in 0..3 {
let mut vel_p = vel;
vel_p[comp] += eps;
let (_, _, v_p) = build_test_values(&pose, &vel_p);
let r_p = factor.residual(&v_p);

let mut vel_m = vel;
vel_m[comp] -= eps;
let (_, _, v_m) = build_test_values(&pose, &vel_m);
let r_m = factor.residual(&v_m);

let j_num = (r_p - r_m) / (2.0 * eps);
let (s_vel, _) = values.index_of(id_vel).expect("index of vel");
for row in 0..3 {
let diff = (j_analytic[(row, s_vel + comp)] - j_num[row]).abs();
assert!(diff < 1e-5, "vel jacobian mismatch comp={comp}, row={row}, diff={diff}");
}
}
}

/// The attitude block was previously untested, and it is exactly where the
/// engine's two rotation parameterisations diverge: `-R^T skew(v)` for a
/// left/global retraction versus `+skew(R^T v)` for the right retraction the
/// solver actually applies. Finite-differencing through `apply_delta` pins it.
#[test]
fn odometer_attitude_jacobian_matches_the_solver_retraction() {
let pose = [100.0, 200.0, 300.0, 0.25, -0.4, 0.55];
let vel = [12.0, -3.0, 5.0];
let (id_pose, id_vel, values) = build_test_values(&pose, &vel);
let factor = OdometerVelocityFactor::new(
id_pose,
id_vel,
Vector3::new(10.0, 0.0, 0.0),
Vector3::new(0.01, 0.01, 0.01),
);
let (s_pose, _) = values.index_of(id_pose).expect("index of pose");
let j_analytic = factor.jacobian(&values);
let eps = 1e-7;
for row in 0..3 {
for comp in 0..3 {
let mut d = Vector3::zeros();
d[comp] = eps;
let up = pose_after_solver_delta(&pose, d);
let dn = pose_after_solver_delta(&pose, -d);
let r_up = factor.residual(&build_test_values(&up, &vel).2)[row];
let r_dn = factor.residual(&build_test_values(&dn, &vel).2)[row];
let j_num = (r_up - r_dn) / (2.0 * eps);
let got = j_analytic[(row, s_pose + 3 + comp)];
assert!(
(got - j_num).abs() < 1e-5,
"attitude J[{row},{comp}]: analytic={got}, solver-path finite difference={j_num}"
);
}
}
}

#[test]
fn odometer_jacobian_blocks_match_the_closed_form() {
// dv_body/dv_ecef = R_ecef->body = q^-1, and dv_body/d(dtheta) = skew(v_body)
// under the solver's right retraction. Translation columns stay zero.
let pose = [0.0, 0.0, 0.0, 0.0, 0.0, std::f64::consts::FRAC_PI_2];
let vel = [10.0, 0.0, 0.0];
let (id_pose, id_vel, values) = build_test_values(&pose, &vel);
let factor = OdometerVelocityFactor::new(
id_pose,
id_vel,
Vector3::new(0.0, -10.0, 0.0),
Vector3::new(0.01, 0.01, 0.01),
);
let (s_pose, _) = values.index_of(id_pose).expect("pose offset");
let (s_vel, _) = values.index_of(id_vel).expect("vel offset");
let j = factor.jacobian(&values);
let near = |got: f64, want: f64| assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
// Yaw 90 deg: q^-1 = [[0,1,0],[-1,0,0],[0,0,1]].
near(j[(0, s_vel)], 0.0);
near(j[(0, s_vel + 1)], 1.0);
near(j[(0, s_vel + 2)], 0.0);
near(j[(1, s_vel)], -1.0);
near(j[(1, s_vel + 1)], 0.0);
near(j[(1, s_vel + 2)], 0.0);
near(j[(2, s_vel)], 0.0);
near(j[(2, s_vel + 1)], 0.0);
near(j[(2, s_vel + 2)], 1.0);
// v_body = q^-1 [10,0,0] = [0,-10,0]; skew(v_body) = [[0,0,-10],[0,0,0],[10,0,0]].
near(j[(0, s_pose + 3)], 0.0);
near(j[(0, s_pose + 4)], 0.0);
near(j[(0, s_pose + 5)], -10.0);
near(j[(1, s_pose + 3)], 0.0);
near(j[(1, s_pose + 4)], 0.0);
near(j[(1, s_pose + 5)], 0.0);
near(j[(2, s_pose + 3)], 10.0);
near(j[(2, s_pose + 4)], 0.0);
near(j[(2, s_pose + 5)], 0.0);
// Translation never enters a body-frame velocity constraint.
for c in 0..3 {
near(j[(0, c)], 0.0);
near(j[(1, c)], 0.0);
near(j[(2, c)], 0.0);
}
// The residual really is zero at this configuration (measurement = [0,-10,0]).
assert!(factor.residual(&values).norm() < 1e-12);
}

#[test]
fn odometer_information_is_the_inverse_of_each_axis_variance() {
// 1/0.04 = 25, 1/0.09 = 1/0.09, 1/0.16 = 6.25 — diagonals only.
let f = OdometerVelocityFactor::new(
VariableId::new(1),
VariableId::new(2),
Vector3::zeros(),
Vector3::new(0.04, 0.09, 0.16),
);
let info = f.information();
assert_eq!(info[(0, 0)], 25.0);
assert_eq!(info[(1, 1)], 1.0 / 0.09);
assert_eq!(info[(2, 2)], 6.25);
for r in 0..3 {
for c in 0..3 {
if r != c {
assert_eq!(info[(r, c)], 0.0);
}
}
}
}

#[test]
fn sanitize_variance_branches() {
// Finite and positive: floored at 1e-6, passed through otherwise.
assert_eq!(sanitize_variance(4.0), 4.0);
assert_eq!(sanitize_variance(1e-9), 1e-6);
// Not finite or not positive: replaced by the 1e-4 default.
assert_eq!(sanitize_variance(0.0), 1e-4);
assert_eq!(sanitize_variance(-1.0), 1e-4);
assert_eq!(sanitize_variance(f64::NAN), 1e-4);
assert_eq!(sanitize_variance(f64::INFINITY), 1e-4);
}

#[test]
fn odometer_information_stays_finite_for_degenerate_variances() {
let f = OdometerVelocityFactor::new(
VariableId::new(1),
VariableId::new(2),
Vector3::zeros(),
Vector3::new(f64::NAN, -1.0, f64::INFINITY),
);
let info = f.information();
assert_eq!(info, DMatrix::from_diagonal(&DVector::from_vec(vec![1e4, 1e4, 1e4])));
}

#[test]
fn test_aux_factors_nan_and_degenerate_handling() {
let pose_nan = [f64::NAN, 200.0, 300.0, 0.0, 0.0, 0.0];
let vel = [10.0, 0.0, 0.0];
let (id_pose, id_vel, values_nan) = build_test_values(&pose_nan, &vel);

let odo = OdometerVelocityFactor::new(
id_pose,
id_vel,
Vector3::new(10.0, 0.0, 0.0),
Vector3::new(f64::NAN, -1.0, 0.0),
);

assert_eq!(odo.residual(&values_nan), DVector::zeros(3));
assert_eq!(odo.jacobian(&values_nan), DMatrix::zeros(3, values_nan.total_dim()));

let info = odo.information();
assert!(info[(0, 0)].is_finite());
assert!(info[(1, 1)].is_finite());
assert!(info[(2, 2)].is_finite());

// A NaN in the *velocity* must also disarm the factor, not just a NaN pose.
let (_, _, values_nan_vel) = build_test_values(&[0.; 6], &[f64::NAN, 0.0, 0.0]);
assert_eq!(odo.residual(&values_nan_vel), DVector::zeros(3));
assert_eq!(odo.jacobian(&values_nan_vel), DMatrix::zeros(3, 9));
}

#[test]
fn odometer_degrades_to_zero_when_a_variable_is_absent() {
// Only the velocity is in the graph: neither the `get` nor the `index_of`
// lookup succeeds, and the documented fallback is a zero residual / Jacobian.
let id_vel = VariableId::new(2);
let mut vars = BTreeMap::new();
let mut n = VariableNode::new(id_vel, VariableKind::Velocity { epoch: 1 });
n.set_value(&[10.0, 0.0, 0.0]);
vars.insert(id_vel, n);
let values = VariableValues::build(&vars);

let f = OdometerVelocityFactor::new(
VariableId::new(1),
id_vel,
Vector3::new(10.0, 0.0, 0.0),
Vector3::new(0.01, 0.01, 0.01),
);
assert_eq!(f.residual(&values), DVector::zeros(3));
assert_eq!(f.jacobian(&values), DMatrix::zeros(3, 3));
}

// ---------------------------------------------------------- Dual antenna

#[test]
fn test_dual_antenna_heading_factor_residual_and_jacobian() {
let pose = [100.0, 200.0, 300.0, 0.0, 0.0, 0.0];
let vel = [0.0, 0.0, 0.0];
let (id_pose, _, values) = build_test_values(&pose, &vel);

let q = UnitQuaternion::from_scaled_axis(Vector3::new(pose[3], pose[4], pose[5]));
let baseline_body = Vector3::new(1.0, 0.0, 0.0);
let measured_ecef = q * baseline_body;

let factor = DualAntennaHeadingFactor::new(
id_pose,
baseline_body,
measured_ecef,
Vector3::new(0.005, 0.005, 0.005),
);

let r = factor.residual(&values);
assert!(r.norm() < 1e-12);
assert_eq!(factor.variables(), &[id_pose]);
assert_eq!(factor.robust_threshold(), Some(3.0));

let j_analytic = factor.jacobian(&values);
let eps = 1e-7;

for comp in 0..3 {
let mut pose_p = pose;
pose_p[3 + comp] += eps;
let (_, _, v_p) = build_test_values(&pose_p, &vel);
let r_p = factor.residual(&v_p);

let mut pose_m = pose;
pose_m[3 + comp] -= eps;
let (_, _, v_m) = build_test_values(&pose_m, &vel);
let r_m = factor.residual(&v_m);

let j_num = (r_p - r_m) / (2.0 * eps);
let (s_pose, _) = values.index_of(id_pose).expect("index of pose");
for row in 0..3 {
let diff = (j_analytic[(row, s_pose + 3 + comp)] - j_num[row]).abs();
assert!(diff < 1e-5, "attitude jacobian mismatch comp={comp}, row={row}, diff={diff}");
}
}
}

/// The test above perturbs the *rotation vector*; at zero attitude that is
/// indistinguishable from the solver's right retraction. Re-run the comparison
/// at a non-zero attitude, where the two parameterisations genuinely differ.
#[test]
fn dual_antenna_attitude_jacobian_holds_at_nonzero_attitude() {
let pose = [0.0, 0.0, 0.0, 0.35, -0.5, 0.2];
let vel = [0.0; 3];
let (id_pose, _, values) = build_test_values(&pose, &vel);
let baseline_body = Vector3::new(0.30, -0.10, 0.05);
let q = UnitQuaternion::from_scaled_axis(Vector3::new(pose[3], pose[4], pose[5]));
let factor = DualAntennaHeadingFactor::new(
id_pose,
baseline_body,
q * baseline_body,
Vector3::new(0.005, 0.005, 0.005),
);
assert!(factor.residual(&values).norm() < 1e-12);
let (s_pose, _) = values.index_of(id_pose).expect("index of pose");
let j_analytic = factor.jacobian(&values);
let eps = 1e-7;
for row in 0..3 {
for comp in 0..3 {
let mut d = Vector3::zeros();
d[comp] = eps;
let up = pose_after_solver_delta(&pose, d);
let dn = pose_after_solver_delta(&pose, -d);
let j_num = (factor.residual(&build_test_values(&up, &vel).2)[row]
- factor.residual(&build_test_values(&dn, &vel).2)[row])
/ (2.0 * eps);
let got = j_analytic[(row, s_pose + 3 + comp)];
assert!(
(got - j_num).abs() < 1e-5,
"dual-antenna J[{row},{comp}]: analytic={got}, solver-path finite difference={j_num}"
);
}
}
}

#[test]
fn dual_antenna_jacobian_is_minus_r_times_skew_of_the_body_baseline() {
// Hand table at a yaw of exactly 90 degrees. `R = Rz(pi/2) =
// [[0,-1,0],[1,0,0],[0,0,1]]` and `b = [L, 0, 0]`, so
//   skew(b) = [[0,0,0],[0,0,-L],[0,L,0]]
//   R skew(b) = [[0, 0, L], [0, 0, 0], [0, L, 0]]
//   J = -R skew(b) = [[0, 0, -L], [0, 0, 0], [0, -L, 0]]
// Compare with the (incorrect) `-skew(R b)`, which at this attitude would be
// [[0, 0, 0], [0, 0, 0], [0, -L, 0]] — column 0 would vanish.
let pose = [0.0, 0.0, 0.0, 0.0, 0.0, std::f64::consts::FRAC_PI_2];
let vel = [0.0; 3];
let (id_pose, _, values) = build_test_values(&pose, &vel);
let l = 1.5_f64;
let baseline_body = Vector3::new(l, 0.0, 0.0);
let q = UnitQuaternion::from_scaled_axis(Vector3::new(pose[3], pose[4], pose[5]));
let factor = DualAntennaHeadingFactor::new(
id_pose,
baseline_body,
q * baseline_body,
Vector3::new(0.005, 0.005, 0.005),
);
let (s_pose, _) = values.index_of(id_pose).expect("index of pose");
let j = factor.jacobian(&values);
let near = |got: f64, want: f64| assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
let want = [
[0.0, 0.0, -l],
[0.0, 0.0, 0.0],
[0.0, -l, 0.0],
];
for row in 0..3 {
for comp in 0..3 {
near(j[(row, s_pose + 3 + comp)], want[row][comp]);
}
}
// Translation columns are identically zero: the baseline lives in the body
// frame, so only attitude is observable from a single pose.
for c in 0..3 {
for row in 0..3 {
near(j[(row, c)], 0.0);
}
}
}

#[test]
fn dual_antenna_information_is_diagonal_inverse_variance() {
let f = DualAntennaHeadingFactor::new(
VariableId::new(1),
Vector3::new(1.0, 0.0, 0.0),
Vector3::new(0.0, 1.0, 0.0),
Vector3::new(0.01, 0.04, 0.09),
);
let info = f.information();
assert_eq!(info, DMatrix::from_diagonal(&DVector::from_vec(vec![100.0, 25.0, 1.0 / 0.09])));
}

#[test]
fn dual_antenna_degrades_to_zero_for_absent_or_non_finite_pose() {
let factor = DualAntennaHeadingFactor::new(
VariableId::new(1),
Vector3::new(1.0, 0.0, 0.0),
Vector3::new(0.0, 1.0, 0.0),
Vector3::new(0.01, 0.01, 0.01),
);
// Pose entirely absent.
let mut vars = BTreeMap::new();
let id_vel = VariableId::new(2);
let mut n = VariableNode::new(id_vel, VariableKind::Velocity { epoch: 1 });
n.set_value(&[0.0; 3]);
vars.insert(id_vel, n);
let values = VariableValues::build(&vars);
assert_eq!(factor.residual(&values), DVector::zeros(3));
assert_eq!(factor.jacobian(&values), DMatrix::zeros(3, 3));
// Pose present but non-finite.
let (_, _, bad) = build_test_values(&[f64::INFINITY, 0., 0., 0., 0., 0.], &[0.0; 3]);
assert_eq!(factor.residual(&bad), DVector::zeros(3));
assert_eq!(factor.jacobian(&bad), DMatrix::zeros(3, 9));
}
