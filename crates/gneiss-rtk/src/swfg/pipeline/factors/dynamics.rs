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

        let v_skew = Matrix3::new(
            0.0, -v_ecef.z, v_ecef.y,
            v_ecef.z, 0.0, -v_ecef.x,
            -v_ecef.y, v_ecef.x, 0.0,
        );
        let dv_drot = -r_e2b * v_skew;

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
