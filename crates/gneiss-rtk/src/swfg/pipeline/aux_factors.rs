//! Auxiliary sensor factors for the sliding-window factor graph.
//!
//! Includes body-frame velocity constraints (odometer / DVL) and dual-antenna
//! GNSS heading baseline constraints.

use nalgebra::{DMatrix, DVector, Matrix3, UnitQuaternion, Vector3};

use crate::swfg::factor::Factor;
use crate::swfg::variables::{VariableId, VariableValues};

/// Wheel Odometry / DVL 3D Velocity Factor in the body frame.
///
/// Constrains the body-frame forward velocity $v_x^b$ measured by an odometer or DVL,
/// and enforces Non-Holonomic Constraints ($v_y^b \approx 0, v_z^b \approx 0$).
///
/// Residual:
///   r(pose, vel) = R_b^e(q)^T * v^e - v_measured^b
#[derive(Clone, Debug)]
pub struct OdometerVelocityFactor {
    pub var_pose: VariableId,
    pub var_vel: VariableId,
    pub measured_v_body: Vector3<f64>,
    pub variances: Vector3<f64>,
    pub variables: Vec<VariableId>,
}

impl OdometerVelocityFactor {
    /// Create a new odometer/DVL velocity factor.
    pub fn new(
        var_pose: VariableId,
        var_vel: VariableId,
        measured_v_body: Vector3<f64>,
        variances: Vector3<f64>,
    ) -> Self {
        Self {
            var_pose,
            var_vel,
            measured_v_body,
            variances,
            variables: vec![var_pose, var_vel],
        }
    }
}

impl Factor for OdometerVelocityFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let (Some(p), Some(v)) = (values.get(self.var_pose), values.get(self.var_vel)) else {
            return DVector::zeros(3);
        };
        if !p.iter().all(|x| x.is_finite()) || !v.iter().all(|x| x.is_finite()) {
            return DVector::zeros(3);
        }

        let q = UnitQuaternion::from_scaled_axis(Vector3::new(p[3], p[4], p[5]));
        let v_ecef = Vector3::new(v[0], v[1], v[2]);
        let v_body_pred = q.inverse() * v_ecef;
        let diff = v_body_pred - self.measured_v_body;
        DVector::from_vec(vec![diff.x, diff.y, diff.z])
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(3, total_dim);
        let (Some((s_pose, _)), Some((s_vel, _))) =
            (values.index_of(self.var_pose), values.index_of(self.var_vel))
        else {
            return j;
        };
        let (Some(p), Some(v)) = (values.get(self.var_pose), values.get(self.var_vel)) else {
            return j;
        };
        if !p.iter().all(|x| x.is_finite()) || !v.iter().all(|x| x.is_finite()) {
            return j;
        }

        let q = UnitQuaternion::from_scaled_axis(Vector3::new(p[3], p[4], p[5]));
        let v_ecef = Vector3::new(v[0], v[1], v[2]);
        let v_body = q.inverse() * v_ecef;
        let r_t = q.inverse().to_rotation_matrix().into_inner();

        let skew_v = Matrix3::new(
            0.0, -v_body.z, v_body.y,
            v_body.z, 0.0, -v_body.x,
            -v_body.y, v_body.x, 0.0,
        );

        for k in 0..3 {
            for row in 0..3 {
                j[(row, s_vel + k)] = r_t[(row, k)];
                j[(row, s_pose + 3 + k)] = skew_v[(row, k)];
            }
        }
        j
    }

    fn information(&self) -> DMatrix<f64> {
        let mut info = DMatrix::zeros(3, 3);
        let vx = sanitize_variance(self.variances.x);
        let vy = sanitize_variance(self.variances.y);
        let vz = sanitize_variance(self.variances.z);
        info[(0, 0)] = 1.0 / vx;
        info[(1, 1)] = 1.0 / vy;
        info[(2, 2)] = 1.0 / vz;
        info
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(3.0)
    }
}

/// Dual-Antenna GNSS Baseline / Heading Factor.
///
/// Constrains the baseline vector b^b between primary and secondary antennas
/// in the body frame relative to the estimated ECEF baseline dp^e = p_sec^e - p_prim^e.
///
/// Residual:
///   r(pose) = R_b^e(q) * b^b - dp_measured^e
#[derive(Clone, Debug)]
pub struct DualAntennaHeadingFactor {
    pub var_pose: VariableId,
    pub baseline_body: Vector3<f64>,
    pub measured_baseline_ecef: Vector3<f64>,
    pub variances: Vector3<f64>,
    pub variables: Vec<VariableId>,
}

impl DualAntennaHeadingFactor {
    /// Create a new dual-antenna heading factor.
    pub fn new(
        var_pose: VariableId,
        baseline_body: Vector3<f64>,
        measured_baseline_ecef: Vector3<f64>,
        variances: Vector3<f64>,
    ) -> Self {
        Self {
            var_pose,
            baseline_body,
            measured_baseline_ecef,
            variances,
            variables: vec![var_pose],
        }
    }
}

impl Factor for DualAntennaHeadingFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let Some(p) = values.get(self.var_pose) else {
            return DVector::zeros(3);
        };
        if !p.iter().all(|x| x.is_finite()) {
            return DVector::zeros(3);
        }

        let q = UnitQuaternion::from_scaled_axis(Vector3::new(p[3], p[4], p[5]));
        let b_ecef_pred = q * self.baseline_body;
        let diff = b_ecef_pred - self.measured_baseline_ecef;
        DVector::from_vec(vec![diff.x, diff.y, diff.z])
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(3, total_dim);
        let Some((s_pose, _)) = values.index_of(self.var_pose) else {
            return j;
        };
        let Some(p) = values.get(self.var_pose) else {
            return j;
        };
        if !p.iter().all(|x| x.is_finite()) {
            return j;
        }

        let q = UnitQuaternion::from_scaled_axis(Vector3::new(p[3], p[4], p[5]));
        let b_ecef = q * self.baseline_body;

        let skew_b_ecef = Matrix3::new(
            0.0, -b_ecef.z, b_ecef.y,
            b_ecef.z, 0.0, -b_ecef.x,
            -b_ecef.y, b_ecef.x, 0.0,
        );
        let j_att = -skew_b_ecef;

        for k in 0..3 {
            for row in 0..3 {
                j[(row, s_pose + 3 + k)] = j_att[(row, k)];
            }
        }
        j
    }

    fn information(&self) -> DMatrix<f64> {
        let mut info = DMatrix::zeros(3, 3);
        let vx = sanitize_variance(self.variances.x);
        let vy = sanitize_variance(self.variances.y);
        let vz = sanitize_variance(self.variances.z);
        info[(0, 0)] = 1.0 / vx;
        info[(1, 1)] = 1.0 / vy;
        info[(2, 2)] = 1.0 / vz;
        info
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(3.0)
    }
}

#[inline]
fn sanitize_variance(v: f64) -> f64 {
    if v.is_finite() && v > 0.0 {
        v.max(1e-6)
    } else {
        1e-4
    }
}

#[cfg(test)]
mod tests {
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
    }
}
