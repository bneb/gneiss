//! Vehicle motion and kinematic constraints (NHC, Doppler dynamics).

use nalgebra::{DMatrix, DVector, Matrix3, UnitQuaternion, Vector3};

use crate::swfg::factor::Factor;
use crate::swfg::variables::{VariableId, VariableValues};

/// Non-Holonomic Constraint (NHC) Factor for land vehicles.
/// Constrains the lateral (y) and vertical (z) velocity in the body frame to be near zero.
#[derive(Clone, Debug)]
pub struct NhcFactor {
    pub var_pose: VariableId,
    pub var_vel: VariableId,
    pub variance_y: f64,
    pub variance_z: f64,
    pub variables: Vec<VariableId>,
}

impl Factor for NhcFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let pose = match values.get(self.var_pose) {
            Some(p) => p,
            None => return DVector::zeros(2),
        };
        let vel = match values.get(self.var_vel) {
            Some(v) => v,
            None => return DVector::zeros(2),
        };

        let q = UnitQuaternion::from_scaled_axis(Vector3::new(pose[3], pose[4], pose[5]));
        let v_ecef = Vector3::new(vel[0], vel[1], vel[2]);
        let v_body = q.inverse_transform_vector(&v_ecef);

        DVector::from_vec(vec![v_body.y, v_body.z])
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(2, total_dim);

        let start_pose = match values.index_of(self.var_pose) {
            Some((s, _)) => s,
            None => return j,
        };
        let start_vel = match values.index_of(self.var_vel) {
            Some((s, _)) => s,
            None => return j,
        };

        let pose = match values.get(self.var_pose) {
            Some(p) => p,
            None => return j,
        };
        let vel = match values.get(self.var_vel) {
            Some(v) => v,
            None => return j,
        };

        let q = UnitQuaternion::from_scaled_axis(Vector3::new(pose[3], pose[4], pose[5]));
        let r_b2e = q.to_rotation_matrix().into_inner();
        let r_e2b = r_b2e.transpose();
        let v_ecef = Vector3::new(vel[0], vel[1], vel[2]);

        let dv_dvel = r_e2b;
        j[(0, start_vel)] = dv_dvel[(1, 0)];
        j[(0, start_vel + 1)] = dv_dvel[(1, 1)];
        j[(0, start_vel + 2)] = dv_dvel[(1, 2)];
        j[(1, start_vel)] = dv_dvel[(2, 0)];
        j[(1, start_vel + 1)] = dv_dvel[(2, 1)];
        j[(1, start_vel + 2)] = dv_dvel[(2, 2)];

        // `SlidingWindowSolver::apply_delta` (swfg/solver/mod.rs:336) retracts
        // the attitude on the RIGHT: `q_new = q_old * exp(delta)`, i.e.
        // `R_new = R * Exp(dtheta)`.  Then `R_new^T v = Exp(-dtheta) R^T v`,
        // whose first-order variation is `+skew(R^T v) dtheta` — the skew of
        // the *body-frame* velocity, not `-R^T skew(v)` (that form belongs to a
        // left/global retraction, `R_new = Exp(dtheta) * R`).  This also matches
        // `OdometerVelocityFactor` in `pipeline/aux_factors.rs`, which residuals
        // the same `q^-1 * v_ecef` and therefore must share this Jacobian.
        let v_body = r_e2b * v_ecef;
        let dv_drot = Matrix3::new(
            0.0, -v_body.z, v_body.y,
            v_body.z, 0.0, -v_body.x,
            -v_body.y, v_body.x, 0.0,
        );

        j[(0, start_pose + 3)] = dv_drot[(1, 0)];
        j[(0, start_pose + 4)] = dv_drot[(1, 1)];
        j[(0, start_pose + 5)] = dv_drot[(1, 2)];
        j[(1, start_pose + 3)] = dv_drot[(2, 0)];
        j[(1, start_pose + 4)] = dv_drot[(2, 1)];
        j[(1, start_pose + 5)] = dv_drot[(2, 2)];

        j
    }

    fn information(&self) -> DMatrix<f64> {
        let mut info = DMatrix::zeros(2, 2);
        info[(0, 0)] = 1.0 / self.variance_y.max(1e-4);
        info[(1, 1)] = 1.0 / self.variance_z.max(1e-4);
        info
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(3.0)
    }
}

/// Between-epoch dynamics constraint using Doppler-derived velocity.
#[derive(Clone)]
pub struct DopplerVelocityFactor {
    pub var_pose_prev: VariableId,
    pub var_pose_curr: VariableId,
    pub velocity_ecef: Vector3<f64>,
    pub dt: f64,
    pub variance_m2: f64,
    pub variables: [VariableId; 2],
}

impl std::fmt::Debug for DopplerVelocityFactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DopplerVelocityFactor")
            .field("dt", &self.dt)
            .finish()
    }
}

impl Factor for DopplerVelocityFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let prev = match values.get(self.var_pose_prev) {
            Some(p) => p,
            None => return DVector::zeros(3),
        };
        let curr = match values.get(self.var_pose_curr) {
            Some(c) => c,
            None => return DVector::zeros(3),
        };
        let predicted_delta = self.velocity_ecef * self.dt;
        let actual_delta = Vector3::new(
            curr[0] - prev[0],
            curr[1] - prev[1],
            curr[2] - prev[2],
        );
        let r = actual_delta - predicted_delta;
        DVector::from_column_slice(r.as_slice())
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(3, total_dim);
        if let Some((s_prev, _)) = values.index_of(self.var_pose_prev) {
            for k in 0..3 {
                j[(k, s_prev + k)] = -1.0;
            }
        }
        if let Some((s_curr, _)) = values.index_of(self.var_pose_curr) {
            for k in 0..3 {
                j[(k, s_curr + k)] = 1.0;
            }
        }
        j
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::identity(3, 3) / self.variance_m2.max(1e-4)
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(10.0)
    }
}

#[cfg(test)]
mod tests {
use super::*;
use crate::swfg::variables::{VariableKind, VariableNode};
use std::collections::BTreeMap;

/// Pack `(id, kind, value)` triples in `VariableId` order.
fn values_of(nodes: &[(VariableId, VariableKind, Vec<f64>)]) -> VariableValues {
let mut vars = BTreeMap::new();
for (id, kind, val) in nodes {
vars.insert(*id, VariableNode { id: *id, kind: *kind, value: DVector::from_vec(val.clone()) });
}
VariableValues::build(&vars)
}

/// Reproduce `SlidingWindowSolver::apply_delta` (swfg/solver/mod.rs:336) for a
/// pose variable: the attitude increment is applied on the RIGHT of the nominal
/// attitude, `q_new = q_old * exp(delta)`, and the stored rotation vector is
/// re-extracted with `scaled_axis()`. Every attitude Jacobian here must be the
/// derivative with respect to *that* parameterisation.
fn pose_after_solver_delta(pose: &[f64; 6], delta: Vector3<f64>) -> [f64; 6] {
let q_old = UnitQuaternion::from_scaled_axis(Vector3::new(pose[3], pose[4], pose[5]));
let dq = UnitQuaternion::from_scaled_axis(delta);
let rot = (q_old * dq).scaled_axis();
[pose[0], pose[1], pose[2], rot.x, rot.y, rot.z]
}

fn nhc_at(pose: VariableId, vel: VariableId, vy: f64, vz: f64) -> NhcFactor {
NhcFactor { var_pose: pose, var_vel: vel, variance_y: vy, variance_z: vz, variables: vec![] }
}

fn doppler(prev: VariableId, curr: VariableId, v: [f64; 3], dt: f64, var: f64) -> DopplerVelocityFactor {
DopplerVelocityFactor {
var_pose_prev: prev,
var_pose_curr: curr,
velocity_ecef: Vector3::new(v[0], v[1], v[2]),
dt,
variance_m2: var,
variables: [prev, curr],
}
}

// ------------------------------------------------------------------- NHC

#[test]
fn nhc_variables_and_robust_threshold_are_advertised() {
let p = VariableId::new(1);
let v = VariableId::new(2);
let nhc = NhcFactor { var_pose: p, var_vel: v, variance_y: 0.01, variance_z: 0.01, variables: vec![p, v] };
assert_eq!(nhc.variables(), &[p, v]);
assert_eq!(nhc.robust_threshold(), Some(3.0));
}

#[test]
fn nhc_information_is_the_inverse_of_the_lateral_and_vertical_variances() {
// 1/0.04 = 25 and 1/0.16 = 6.25, both exact in binary.
let info = nhc_at(VariableId::new(1), VariableId::new(2), 0.04, 0.16).information();
assert_eq!(info, DMatrix::from_diagonal(&DVector::from_vec(vec![25.0, 6.25])));
// Lateral and vertical slip are uncoupled: no off-diagonal.
assert_eq!(info[(0, 1)], 0.0);
assert_eq!(info[(1, 0)], 0.0);
}

#[test]
fn nhc_information_floors_a_nonpositive_variance() {
// `variance.max(1e-4)` caps the weight at 1e4 instead of dividing by zero.
let info = nhc_at(VariableId::new(1), VariableId::new(2), 0.0, -1.0).information();
assert_eq!(info[(0, 0)], 1e4);
assert_eq!(info[(1, 1)], 1e4);
assert!(info.iter().all(|v| v.is_finite()));
}

/// Hand computation at a yaw of exactly 90 degrees.
///
/// `q = from_scaled_axis([0, 0, pi/2])` is `Rz(pi/2)`, so
///   R_body->ecef = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
///   R_ecef->body = R^T    = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]
/// With `v_ecef = [10, 0, 0]`:
///   v_body = R^T v = [0, -10, 0]  and  r = [v_body.y, v_body.z] = [-10, 0]`
///   dr/dv_ecef  = rows 1 and 2 of R^T      = [-1, 0, 0] and [0, 0, 1]
/// `dr/d(dtheta)` for `r = (R^T v)[1..3]` under the solver's retraction
/// `R_new = R * Exp(dtheta)` (solver/mod.rs:336) is rows 1 and 2 of
/// `skew(R^T v)`: `R_new^T v = Exp(-dtheta) R^T v`, whose first-order variation
/// is `+skew(R^T v) dtheta`. With `v_body = [0, -10, 0]`,
///   skew(v_body) = [[0, -v.z, v.y], [v.z, 0, -v.x], [-v.y, v.x, 0]]
///                = [[0, 0, -10], [0, 0, 0], [10, 0, 0]]
/// so the attitude block is row 1 = `[0, 0, 0]` for `r[0] = v_body.y` and
/// row 2 = `[10, 0, 0]` for `r[1] = v_body.z`.
///
/// Physical reading of the `10`: rolling by `dtheta_x` maps the body-frame
/// velocity by `v_body - dtheta_x e_x x v_body = v_body + (0, 0, 10 dtheta_x)`,
/// so a vertical NHC residual picks up `+10` per radian of roll.
fn nhc_yaw_90_fixture() -> (VariableId, VariableId, VariableValues) {
let pose = VariableId::new(1);
let vel = VariableId::new(2);
let half_pi = std::f64::consts::FRAC_PI_2;
let values = values_of(&[
(pose, VariableKind::Pose { epoch: 0 }, vec![1., 2., 3., 0., 0., half_pi]),
(vel, VariableKind::Velocity { epoch: 0 }, vec![10., 0., 0.]),
]);
(pose, vel, values)
}

#[test]
fn nhc_residual_is_the_body_frame_lateral_and_vertical_velocity() {
let (pose, vel, values) = nhc_yaw_90_fixture();
assert_eq!(nhc_at(pose, vel, 0.01, 0.01).residual(&values), DVector::from_vec(vec![-10.0, 0.0]));
}

#[test]
fn nhc_jacobian_blocks_match_hand_derivation_at_yaw_90() {
let (pose, vel, values) = nhc_yaw_90_fixture();
let j = nhc_at(pose, vel, 0.01, 0.01).jacobian(&values);
let (s_pose, _) = values.index_of(pose).unwrap();
let (s_vel, _) = values.index_of(vel).unwrap();
assert_eq!((s_pose, s_vel), (0, 6));
assert_eq!((j.nrows(), j.ncols()), (2, 9));
// cos(pi/2) is 6.1e-17 in binary, not 0, so compare to the exact closed form
// with an absolute tolerance of 1e-12.
let near = |got: f64, want: f64| {
assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
};
// Velocity block: rows 1 and 2 of R_ecef->body.
near(j[(0, s_vel)], -1.0);
near(j[(0, s_vel + 1)], 0.0);
near(j[(0, s_vel + 2)], 0.0);
near(j[(1, s_vel)], 0.0);
near(j[(1, s_vel + 1)], 0.0);
near(j[(1, s_vel + 2)], 1.0);
// Attitude block: row 1 of skew(v_body) for the lateral residual, row 2 for
// the vertical one. v_body = [0, -10, 0].
near(j[(0, s_pose + 3)], 0.0);
near(j[(0, s_pose + 4)], 0.0);
near(j[(0, s_pose + 5)], 0.0);
near(j[(1, s_pose + 3)], 10.0);
near(j[(1, s_pose + 4)], 0.0);
near(j[(1, s_pose + 5)], 0.0);
// Translation columns are untouched: the NHC constrains velocity expressed in
// a body frame that depends only on attitude.
for c in 0..3 {
near(j[(0, c)], 0.0);
near(j[(1, c)], 0.0);
}
}

#[test]
fn nhc_attitude_jacobian_matches_the_solver_retraction() {
let pose = VariableId::new(1);
let vel = VariableId::new(2);
let pose_val = [1.0, 2.0, 3.0, 0.3, -0.2, 0.4];
let vel_val = [12.0, -3.0, 5.0];
let nhc = nhc_at(pose, vel, 0.01, 0.01);
let values = values_of(&[
(pose, VariableKind::Pose { epoch: 0 }, pose_val.to_vec()),
(vel, VariableKind::Velocity { epoch: 0 }, vel_val.to_vec()),
]);
let (s_pose, _) = values.index_of(pose).unwrap();
let analytic = nhc.jacobian(&values);
let eps = 1e-7;
let residual_row = |row: usize, p: [f64; 6]| {
let v = values_of(&[
(pose, VariableKind::Pose { epoch: 0 }, p.to_vec()),
(vel, VariableKind::Velocity { epoch: 0 }, vel_val.to_vec()),
]);
nhc.residual(&v)[row]
};
for row in 0..2 {
for comp in 0..3 {
let mut d = Vector3::zeros();
d[comp] = eps;
let numeric = (residual_row(row, pose_after_solver_delta(&pose_val, d))
- residual_row(row, pose_after_solver_delta(&pose_val, -d)))
/ (2.0 * eps);
let got = analytic[(row, s_pose + 3 + comp)];
assert!(
(got - numeric).abs() < 1e-5,
"NHC attitude J[{row},{comp}]: analytic={got}, solver-path finite difference={numeric}"
);
}
}
}

#[test]
fn nhc_zero_residual_for_body_forward_motion_at_any_attitude() {
// If v_ecef is exactly the body x-axis rotated into ECEF, then v_body = [v,0,0]
// and both constrained components vanish whatever the attitude.
for (rx, ry, rz) in [(0.0, 0.0, 0.0), (0.7, -0.4, 1.1), (2.9, 1.3, -2.2)] {
let pose = VariableId::new(1);
let vel = VariableId::new(2);
let q = UnitQuaternion::from_scaled_axis(Vector3::new(rx, ry, rz));
let v_ecef = q * Vector3::new(10.0, 0.0, 0.0);
let values = values_of(&[
(pose, VariableKind::Pose { epoch: 0 }, vec![0., 0., 0., rx, ry, rz]),
(vel, VariableKind::Velocity { epoch: 0 }, vec![v_ecef.x, v_ecef.y, v_ecef.z]),
]);
assert!(
nhc_at(pose, vel, 0.01, 0.01).residual(&values).norm() < 1e-12,
"attitude ({rx},{ry},{rz}) left a lateral/vertical component"
);
}
}

#[test]
fn nhc_degrades_to_zero_when_a_variable_is_absent() {
// Window management can drop a variable out from under a live factor; the
// documented behaviour is a zero residual, never a panic or a stale read.
let pose = VariableId::new(1);
let other = VariableId::new(2);
let values = values_of(&[(other, VariableKind::Velocity { epoch: 0 }, vec![1., 2., 3.])]);
assert_eq!(nhc_at(pose, other, 0.01, 0.01).residual(&values), DVector::zeros(2));
assert_eq!(nhc_at(pose, other, 0.01, 0.01).jacobian(&values), DMatrix::zeros(2, 3));
assert_eq!(nhc_at(pose, VariableId::new(9), 0.01, 0.01).residual(&values), DVector::zeros(2));
}

// ------------------------------------------------- DopplerVelocityFactor

#[test]
fn doppler_residual_is_zero_when_the_pose_follows_the_predicted_velocity() {
// v*dt = [10,0,0]*0.5 = [5,0,0]; the displacement matches exactly.
let prev = VariableId::new(1);
let curr = VariableId::new(2);
let values = values_of(&[
(prev, VariableKind::Pose { epoch: 0 }, vec![0.; 6]),
(curr, VariableKind::Pose { epoch: 1 }, vec![5., 0., 0., 0., 0., 0.]),
]);
let f = doppler(prev, curr, [10.0, 0.0, 0.0], 0.5, 1.0);
assert_eq!(f.residual(&values), DVector::zeros(3));
assert_eq!(f.variables(), &[prev, curr]);
assert_eq!(f.robust_threshold(), Some(10.0));
assert!(format!("{f:?}").contains("0.5"), "Debug must surface dt");
}

#[test]
fn doppler_residual_is_displacement_minus_predicted() {
// Displacement [4,0,0] minus predicted [5,0,0] = [-1,0,0]. Swapping the sign
// convention would give +1, so this pins the ordering.
let prev = VariableId::new(1);
let curr = VariableId::new(2);
let values = values_of(&[
(prev, VariableKind::Pose { epoch: 0 }, vec![0.; 6]),
(curr, VariableKind::Pose { epoch: 1 }, vec![4., 0., 0., 0., 0., 0.]),
]);
assert_eq!(
doppler(prev, curr, [10.0, 0.0, 0.0], 0.5, 1.0).residual(&values),
DVector::from_vec(vec![-1.0, 0.0, 0.0])
);
// Reversing the predicted velocity to [-5, 0, 0] makes the residual
// [4 - (-5), 0, 0] = [9, 0, 0] — pinning both the sign and the `v * dt` scaling.
assert_eq!(
doppler(prev, curr, [-10.0, 0.0, 0.0], 0.5, 1.0).residual(&values),
DVector::from_vec(vec![9.0, 0.0, 0.0])
);
// Halving dt predicts [2.5, 0, 0] and leaves [1.5, 0, 0].
assert_eq!(
doppler(prev, curr, [10.0, 0.0, 0.0], 0.25, 1.0).residual(&values),
DVector::from_vec(vec![1.5, 0.0, 0.0])
);
}

#[test]
fn doppler_jacobian_is_minus_identity_then_identity() {
let prev = VariableId::new(1);
let curr = VariableId::new(2);
let values = values_of(&[
(prev, VariableKind::Pose { epoch: 0 }, vec![0.; 6]),
(curr, VariableKind::Pose { epoch: 1 }, vec![0.; 6]),
]);
let j = doppler(prev, curr, [1.0, 2.0, 3.0], 0.5, 1.0).jacobian(&values);
assert_eq!((j.nrows(), j.ncols()), (3, 12));
for r in 0..3 {
for c in 0..12 {
let expect = match c {
x if x == r => -1.0,
x if x == 6 + r => 1.0,
_ => 0.0,
};
assert_eq!(j[(r, c)], expect, "J[{r},{c}]");
}
}
}

#[test]
fn doppler_information_is_inverse_variance_with_a_floor() {
let prev = VariableId::new(1);
let curr = VariableId::new(2);
assert_eq!(doppler(prev, curr, [0.0; 3], 1.0, 0.25).information(), DMatrix::<f64>::identity(3, 3) * 4.0);
assert_eq!(doppler(prev, curr, [0.0; 3], 1.0, 0.0).information(), DMatrix::<f64>::identity(3, 3) * 1e4);
}

#[test]
fn doppler_degrades_to_zero_when_a_pose_is_absent() {
let prev = VariableId::new(1);
let curr = VariableId::new(2);
let values = values_of(&[(curr, VariableKind::Pose { epoch: 1 }, vec![0.; 6])]);
let f = doppler(prev, curr, [1.0, 0.0, 0.0], 1.0, 1.0);
assert_eq!(f.residual(&values), DVector::zeros(3));
// The missing variable contributes no columns at all: no mass may leak into
// the block of a variable that is not in the graph.
let j = f.jacobian(&values);
assert_eq!((j.nrows(), j.ncols()), (3, 6));
for c in 0..6 {
let expected = if c < 3 { 1.0 } else { 0.0 };
for r in 0..3 {
assert_eq!(j[(r, c)], if r == c { expected } else { 0.0 }, "J[{r},{c}]");
}
}
}
}
