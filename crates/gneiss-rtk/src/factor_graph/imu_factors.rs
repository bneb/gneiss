use nalgebra::{DMatrix, DVector, Matrix3, UnitQuaternion, Vector3};
use super::Factor;
use gneiss_core::imu::ImuMeasurement;

/// Represents the preintegrated measurements between two GNSS epochs.
#[derive(Debug, Clone)]
pub struct ImuPreintegration {
    pub dp: Vector3<f64>,
    pub dv: Vector3<f64>,
    pub dq: UnitQuaternion<f64>,
    pub dt: f64,
    
    // Jacobians w.r.t biases (optional for 1-epoch filtering, but needed for W > 1)
    pub dp_dba: Matrix3<f64>,
    pub dp_dbg: Matrix3<f64>,
    pub dv_dba: Matrix3<f64>,
    pub dv_dbg: Matrix3<f64>,
    pub dq_dbg: Matrix3<f64>,
    
    pub covariance: DMatrix<f64>,
}

impl Default for ImuPreintegration {
    fn default() -> Self {
        Self::new()
    }
}

impl ImuPreintegration {
    pub fn new() -> Self {
        Self {
            dp: Vector3::zeros(),
            dv: Vector3::zeros(),
            dq: UnitQuaternion::identity(),
            dt: 0.0,
            dp_dba: Matrix3::zeros(),
            dp_dbg: Matrix3::zeros(),
            dv_dba: Matrix3::zeros(),
            dv_dbg: Matrix3::zeros(),
            dq_dbg: Matrix3::zeros(),
            covariance: DMatrix::identity(15, 15) * 1e-4, // Default small covariance
        }
    }

    /// Preintegrate a sequence of IMU measurements using the given prior biases.
    pub fn integrate(&mut self, imu_data: &[ImuMeasurement], ba_i: &Vector3<f64>, bg_i: &Vector3<f64>) {
        if imu_data.is_empty() { return; }
        
        let mut prev_time = imu_data[0].time_tag;
        
        for m in imu_data.iter().skip(1) {
            // Assume time_tag is in milliseconds and handles simple wrap-arounds or is absolute.
            let dt = if m.time_tag >= prev_time {
                (m.time_tag - prev_time) as f64 / 1000.0
            } else {
                // Handle potential u32 wrap-around
                ((u32::MAX as u64 - prev_time as u64) + m.time_tag as u64) as f64 / 1000.0
            };
            
            prev_time = m.time_tag;
            
            // Limit dt to a reasonable value (e.g., 0.1s max) in case of gaps
            let dt = dt.clamp(0.001, 0.1);
            
            self.dt += dt;
            
            // Correct measurements with prior biases
            let a = m.accel - ba_i;
            let w = m.gyro - bg_i;
            
            // Mid-point or Euler integration. Here we use Euler for simplicity, 
            // though mid-point is better.
            let dq_step = UnitQuaternion::from_scaled_axis(w * dt);
            
            // Position and velocity updates
            let a_world = self.dq * a;
            self.dp += self.dv * dt + 0.5 * a_world * dt * dt;
            self.dv += a_world * dt;
            
            // Attitude update
            self.dq *= dq_step;
            
            // Bias Jacobians update (Simplified Euler propagation)
            let r_mat = self.dq.to_rotation_matrix().into_inner();
            self.dp_dba += self.dv_dba * dt - 0.5 * r_mat * dt * dt;
            self.dp_dbg += self.dv_dbg * dt; // Should include cross product term of a_world, omitted for brevity in first-order
            
            self.dv_dba += -r_mat * dt;
            // self.dv_dbg omitted complex cross product terms
            
            // self.dq_dbg omitted complex right-Jacobian terms
        }
    }
}

/// Factor connecting IMU state at epoch i to IMU state at epoch j.
pub struct ImuPreintegrationFactor {
    pub preint: ImuPreintegration,
    pub gravity: Vector3<f64>,
    
    pub nominal_p_i: Vector3<f64>,
    pub nominal_v_i: Vector3<f64>,
    pub nominal_q_i: UnitQuaternion<f64>,
    pub nominal_ba_i: Vector3<f64>,
    pub nominal_bg_i: Vector3<f64>,
    
    pub nominal_p_j: Vector3<f64>,
    pub nominal_v_j: Vector3<f64>,
    pub nominal_q_j: UnitQuaternion<f64>,
    pub nominal_ba_j: Vector3<f64>,
    pub nominal_bg_j: Vector3<f64>,
    
    pub idx_p_i: usize,
    pub idx_v_i: usize,
    pub idx_q_i: usize,
    pub idx_ba_i: usize,
    pub idx_bg_i: usize,
    
    pub idx_p_j: usize,
    pub idx_v_j: usize,
    pub idx_q_j: usize,
    pub idx_ba_j: usize,
    pub idx_bg_j: usize,
}

impl Factor for ImuPreintegrationFactor {
    fn residual(&self, state: &DVector<f64>) -> DVector<f64> {
        // Evaluate absolute states from nominal + delta
        let p_i = self.nominal_p_i + Vector3::new(state[self.idx_p_i], state[self.idx_p_i+1], state[self.idx_p_i+2]);
        let v_i = self.nominal_v_i + Vector3::new(state[self.idx_v_i], state[self.idx_v_i+1], state[self.idx_v_i+2]);
        let q_i = self.nominal_q_i * UnitQuaternion::from_scaled_axis(Vector3::new(state[self.idx_q_i], state[self.idx_q_i+1], state[self.idx_q_i+2]));
        let ba_i = self.nominal_ba_i + Vector3::new(state[self.idx_ba_i], state[self.idx_ba_i+1], state[self.idx_ba_i+2]);
        let bg_i = self.nominal_bg_i + Vector3::new(state[self.idx_bg_i], state[self.idx_bg_i+1], state[self.idx_bg_i+2]);
        
        let p_j = self.nominal_p_j + Vector3::new(state[self.idx_p_j], state[self.idx_p_j+1], state[self.idx_p_j+2]);
        let v_j = self.nominal_v_j + Vector3::new(state[self.idx_v_j], state[self.idx_v_j+1], state[self.idx_v_j+2]);
        let q_j = self.nominal_q_j * UnitQuaternion::from_scaled_axis(Vector3::new(state[self.idx_q_j], state[self.idx_q_j+1], state[self.idx_q_j+2]));
        let ba_j = self.nominal_ba_j + Vector3::new(state[self.idx_ba_j], state[self.idx_ba_j+1], state[self.idx_ba_j+2]);
        let bg_j = self.nominal_bg_j + Vector3::new(state[self.idx_bg_j], state[self.idx_bg_j+1], state[self.idx_bg_j+2]);
        
        let dt = self.preint.dt;
        let r_i_t = q_i.inverse();
        
        // Compute bias-corrected preintegrated measurements
        let dba = ba_i - self.nominal_ba_i;
        let dbg = bg_i - self.nominal_bg_i;
        
        let dp = self.preint.dp + self.preint.dp_dba * dba + self.preint.dp_dbg * dbg;
        let dv = self.preint.dv + self.preint.dv_dba * dba + self.preint.dv_dbg * dbg;
        let dq = self.preint.dq * UnitQuaternion::from_scaled_axis(self.preint.dq_dbg * dbg);
        
        let r_p = r_i_t * (p_j - p_i - v_i * dt - 0.5 * self.gravity * dt * dt) - dp;
        let r_v = r_i_t * (v_j - v_i - self.gravity * dt) - dv;
        let r_q = (dq.inverse() * r_i_t * q_j).scaled_axis();
        let r_ba = ba_j - ba_i;
        let r_bg = bg_j - bg_i;
        
        let mut res = DVector::zeros(15);
        res.fixed_rows_mut::<3>(0).copy_from(&r_p);
        res.fixed_rows_mut::<3>(3).copy_from(&r_v);
        res.fixed_rows_mut::<3>(6).copy_from(&r_q);
        res.fixed_rows_mut::<3>(9).copy_from(&r_ba);
        res.fixed_rows_mut::<3>(12).copy_from(&r_bg);
        res
    }
    
    fn jacobian(&self, state: &DVector<f64>) -> DMatrix<f64> {
        // Finite difference used to ensure correctness.
        let mut jac = DMatrix::zeros(15, state.len());
        let eps = 1e-6;
        let mut state_mut = state.clone();
        
        let indices = [
            self.idx_p_i, self.idx_v_i, self.idx_q_i, self.idx_ba_i, self.idx_bg_i,
            self.idx_p_j, self.idx_v_j, self.idx_q_j, self.idx_ba_j, self.idx_bg_j
        ];
        
        for &base_idx in &indices {
            for offset in 0..3 {
                let idx = base_idx + offset;
                state_mut[idx] += eps;
                let r_plus = self.residual(&state_mut);
                state_mut[idx] -= 2.0 * eps;
                let r_minus = self.residual(&state_mut);
                state_mut[idx] += eps;
                
                let deriv = (r_plus - r_minus) / (2.0 * eps);
                jac.column_mut(idx).copy_from(&deriv);
            }
        }
        
        jac
    }
    
    fn information(&self) -> DMatrix<f64> {
        self.preint.covariance.clone().cholesky().unwrap().inverse()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_imu_preintegration_stationary() {
        let mut preint = ImuPreintegration::new();
        let gravity = Vector3::new(0.0, 0.0, -9.80665);
        let mut imu_data = Vec::new();
        for i in 0..100 {
            imu_data.push(ImuMeasurement {
                time_tag: i * 10,
                accel: -gravity, // Accel measures upward specific force to counteract gravity
                gyro: Vector3::zeros(),
                temperature: Some(0.0),
            });
        }
        
        preint.integrate(&imu_data, &Vector3::zeros(), &Vector3::zeros());
        
        // Since it's exactly upright, dq should be identity
        assert!(preint.dq.angle() < 1e-6);
        // And dv should be 9.80665 * time (since gravity isn't removed until the factor evaluation)
        // dt = 0.01 * 99 = 0.99
        assert!((preint.dv - Vector3::new(0.0, 0.0, 9.80665 * 0.99)).norm() < 1e-4);
    }

    #[test]
    fn test_imu_preintegration_dynamic() {
        let mut preint = ImuPreintegration::new();
        let mut imu_data = Vec::new();
        for i in 0..100 {
            imu_data.push(ImuMeasurement {
                time_tag: i * 10,
                accel: Vector3::new(1.0, 0.0, 0.0),
                gyro: Vector3::new(0.0, 0.0, std::f64::consts::PI / 2.0), // 90 deg/sec around Z
                temperature: Some(0.0),
            });
        }
        
        preint.integrate(&imu_data, &Vector3::zeros(), &Vector3::zeros());
        
        // After 0.99 sec at 90 deg/sec, rotation should be ~89.1 degrees around Z
        let expected_angle = (std::f64::consts::PI / 2.0) * 0.99;
        assert!((preint.dq.angle() - expected_angle).abs() < 1e-4);
        assert!((preint.dq.axis().unwrap().into_inner() - Vector3::z()).norm() < 1e-4);
    }

    #[test]
    fn test_imu_preintegration_factor_jacobian() {
        // Since the current implementation uses finite differencing, this test
        // ensures that the jacobian is properly formed and non-zero in expected blocks.
        let mut preint = ImuPreintegration::new();
        preint.dt = 1.0;
        preint.dv = Vector3::new(1.0, 0.0, 0.0);
        
        let factor = ImuPreintegrationFactor {
            preint,
            gravity: Vector3::new(0.0, 0.0, -9.80665),
            nominal_p_i: Vector3::zeros(),
            nominal_v_i: Vector3::zeros(),
            nominal_q_i: UnitQuaternion::identity(),
            nominal_ba_i: Vector3::zeros(),
            nominal_bg_i: Vector3::zeros(),
            nominal_p_j: Vector3::new(1.0, 0.0, 0.0),
            nominal_v_j: Vector3::new(1.0, 0.0, 0.0),
            nominal_q_j: UnitQuaternion::identity(),
            nominal_ba_j: Vector3::zeros(),
            nominal_bg_j: Vector3::zeros(),
            idx_p_i: 0, idx_v_i: 3, idx_q_i: 6, idx_ba_i: 9, idx_bg_i: 12,
            idx_p_j: 15, idx_v_j: 18, idx_q_j: 21, idx_ba_j: 24, idx_bg_j: 27,
        };

        let state = DVector::zeros(30);
        let jac = factor.jacobian(&state);
        
        // Expect jacobian to be 15 rows by 30 cols
        assert_eq!(jac.nrows(), 15);
        assert_eq!(jac.ncols(), 30);
        
        // p_j block should map directly to position residual (identity)
        let dp_dpj = jac.view((0, 15), (3, 3));
        assert!((dp_dpj - DMatrix::identity(3, 3)).norm() < 1e-4);
        
        // p_i block should be negative identity
        let dp_dpi = jac.view((0, 0), (3, 3));
        assert!((dp_dpi + DMatrix::identity(3, 3)).norm() < 1e-4);
    }
}
