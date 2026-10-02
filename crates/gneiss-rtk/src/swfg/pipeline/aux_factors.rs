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
        let r_b2e = q.to_rotation_matrix().into_inner();

        // `d(R b)/d(dtheta)` under the solver's retraction `R_new = R *
        // Exp(dtheta)` (swfg/solver/mod.rs:336) is `R (dtheta x b) =
        // -R skew(b)`.  This is *not* `-skew(R b)`: rotation invariance of the
        // cross product gives `skew(R b) = R skew(b) R^T`, so the two differ by
        // a factor of `R^T` and coincide only at the identity attitude.
        let b = &self.baseline_body;
        let skew_b_body = Matrix3::new(
            0.0, -b.z, b.y,
            b.z, 0.0, -b.x,
            -b.y, b.x, 0.0,
        );
        let j_att = -r_b2e * skew_b_body;

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
#[path = "aux_factors_tests.rs"]
mod tests;
