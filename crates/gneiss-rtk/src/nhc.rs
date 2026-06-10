use crate::filter::RtkState;
use crate::engine::updater;
use nalgebra::{DMatrix, DVector, Matrix3};

/// Applies Non-Holonomic Constraints (NHC) to the EKF state.
/// This assumes the vehicle's lateral and vertical velocity in the body frame is zero.
pub fn apply_nhc(
    state: &mut RtkState, 
    sigma_lateral: f64, 
    sigma_vertical: f64,
    imu_mounting_angles: &Option<[f64; 3]>,
    imu_to_nhc_lever_arm: &[f64; 3],
    omega_b: &nalgebra::Vector3<f64>
) -> Result<(), &'static str> {
    let r_e_b = state.attitude.inverse().to_rotation_matrix();
    let v_b_imu = r_e_b * state.velocity;
    
    let l_b = nalgebra::Vector3::new(imu_to_nhc_lever_arm[0], imu_to_nhc_lever_arm[1], imu_to_nhc_lever_arm[2]);
    let v_b_axle = v_b_imu + omega_b.cross(&l_b);

    let r_b_v = if let Some(angles) = imu_mounting_angles {
        // Angles: Roll, Pitch, Yaw
        let roll = angles[0];
        let pitch = angles[1];
        let yaw = angles[2];
        nalgebra::Rotation3::from_euler_angles(roll, pitch, yaw)
    } else {
        nalgebra::Rotation3::identity()
    };
    
    let v_v = r_b_v * v_b_axle;
    
    // Measurement z: we want lateral (y) and vertical (z) velocity in vehicle frame to be 0
    let z = DVector::from_column_slice(&[-v_v.y, -v_v.z]);
    
    let n = state.covariance.nrows();
    let mut h = DMatrix::<f64>::zeros(2, n);
    
    // Jacobian w.r.t velocity in ECEF
    // v_v = R_b_v * v_b_axle = R_b_v * (R_e_b * v_e + omega x l)
    // d(v_v)/dv_e = R_b_v * R_e_b
    let rot = r_b_v * r_e_b;
    let dr_dv = rot.matrix();
    h[(0, 3)] = dr_dv[(1, 0)]; h[(0, 4)] = dr_dv[(1, 1)]; h[(0, 5)] = dr_dv[(1, 2)];
    h[(1, 3)] = dr_dv[(2, 0)]; h[(1, 4)] = dr_dv[(2, 1)]; h[(1, 5)] = dr_dv[(2, 2)];
    
    // Jacobian w.r.t attitude
    // d(v_v)/d(psi) = R_b_v * [v_b_imu x] * R_e^b
    let v_b_skew = Matrix3::new(
        0.0, -v_b_imu.z,  v_b_imu.y,
        v_b_imu.z,  0.0, -v_b_imu.x,
       -v_b_imu.y,  v_b_imu.x,  0.0
    );
    let dr_dpsi = r_b_v.matrix() * v_b_skew * r_e_b.matrix();
    h[(0, 6)] = dr_dpsi[(1, 0)]; h[(0, 7)] = dr_dpsi[(1, 1)]; h[(0, 8)] = dr_dpsi[(1, 2)];
    h[(1, 6)] = dr_dpsi[(2, 0)]; h[(1, 7)] = dr_dpsi[(2, 1)]; h[(1, 8)] = dr_dpsi[(2, 2)];
    
    // Jacobian w.r.t gyro bias
    // v_b_axle = ... + (gyro - bg) x l_b = ... + l_b x bg
    // d(v_v)/dbg = R_b_v * [l_b x]
    let l_b_skew = Matrix3::new(
        0.0, -l_b.z,  l_b.y,
        l_b.z,  0.0, -l_b.x,
       -l_b.y,  l_b.x,  0.0
    );
    let dr_dbg = r_b_v.matrix() * l_b_skew;
    h[(0, 12)] = dr_dbg[(1, 0)]; h[(0, 13)] = dr_dbg[(1, 1)]; h[(0, 14)] = dr_dbg[(1, 2)];
    h[(1, 12)] = dr_dbg[(2, 0)]; h[(1, 13)] = dr_dbg[(2, 1)]; h[(1, 14)] = dr_dbg[(2, 2)];

    let r = DMatrix::from_diagonal(&DVector::from_column_slice(&[
        sigma_lateral * sigma_lateral,
        sigma_vertical * sigma_vertical
    ]));

    updater::update(state, &z, &h, &r, 1e9, None, true).map_err(|_| "NHC update failed")?;
    Ok(())
}

pub fn apply_zupt(state: &mut RtkState, sigma: f64) -> Result<(), &'static str> {
    let z = DVector::from_column_slice(&[-state.velocity.x, -state.velocity.y, -state.velocity.z]);
    let n = state.covariance.nrows();
    let mut h = DMatrix::<f64>::zeros(3, n);
    for i in 0..3 { h[(i, 3 + i)] = 1.0; }
    
    let r = DMatrix::from_diagonal(&DVector::from_element(3, sigma * sigma));
    updater::update(state, &z, &h, &r, 1e9, None, true).map_err(|_| "ZUPT update failed")?;
    Ok(())}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::RtkState;
    use gneiss_core::time::GpsTime;
    use nalgebra::{Vector3, UnitQuaternion};

    fn compute_numerical_jacobian_nhc(
        state: &RtkState,
        imu_mounting_angles: &Option<[f64; 3]>,
        imu_to_nhc_lever_arm: &[f64; 3],
        omega_b: &Vector3<f64>,
        epsilon: f64
    ) -> DMatrix<f64> {
        let n = state.covariance.nrows();
        let mut h_num = DMatrix::zeros(2, n);

        for j in 0..n {
            let mut state_pos = state.clone();
            let mut state_neg = state.clone();

            let mut omega_b_pos = *omega_b;
            let mut omega_b_neg = *omega_b;

            if j >= 6 && j <= 8 {
                let mut dpsi_pos = Vector3::zeros();
                dpsi_pos[j - 6] = epsilon;
                let dq_pos = UnitQuaternion::from_scaled_axis(dpsi_pos);
                state_pos.attitude = dq_pos * state_pos.attitude;

                let mut dpsi_neg = Vector3::zeros();
                dpsi_neg[j - 6] = -epsilon;
                let dq_neg = UnitQuaternion::from_scaled_axis(dpsi_neg);
                state_neg.attitude = dq_neg * state_neg.attitude;
            } else if j >= 3 && j < 6 {
                state_pos.velocity[j - 3] += epsilon;
                state_neg.velocity[j - 3] -= epsilon;
            } else if j >= 12 && j < 15 {
                state_pos.gyro_bias[j - 12] += epsilon;
                state_neg.gyro_bias[j - 12] -= epsilon;
                // omega_b = gyro - bg
                omega_b_pos[j - 12] -= epsilon;
                omega_b_neg[j - 12] += epsilon;
            }

            let get_z = |s: &RtkState, ob: &Vector3<f64>| -> DVector<f64> {
                let r_e_b = s.attitude.inverse().to_rotation_matrix();
                let v_b_imu = r_e_b * s.velocity;
                let l_b = Vector3::new(imu_to_nhc_lever_arm[0], imu_to_nhc_lever_arm[1], imu_to_nhc_lever_arm[2]);
                let v_b_axle = v_b_imu + ob.cross(&l_b);
                let r_b_v = if let Some(angles) = imu_mounting_angles {
                    nalgebra::Rotation3::from_euler_angles(angles[0], angles[1], angles[2])
                } else {
                    nalgebra::Rotation3::identity()
                };
                let v_v = r_b_v * v_b_axle;
                DVector::from_column_slice(&[v_v.y, v_v.z])
            };

            let meas_pos = get_z(&state_pos, &omega_b_pos);
            let meas_neg = get_z(&state_neg, &omega_b_neg);
            let col = (meas_pos - meas_neg) / (2.0 * epsilon);
            h_num.set_column(j, &col);
        }
        h_num
    }

    #[test]
    fn test_nhc_jacobian() {
        let mut state = RtkState::new(GpsTime::new(0, 0.0), gneiss_core::coords::Coordinate::new(Vector3::new(1.0, 2.0, 3.0), gneiss_core::coords::Datum::WGS84, gneiss_core::coords::Frame::ECEF, GpsTime::new(0, 0.0)), 1.0);
        state.velocity = Vector3::new(10.0, 20.0, 30.0);
        state.attitude = UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3);
        
        let imu_mounting_angles = Some([0.05, -0.05, 0.1]);
        let imu_to_nhc_lever_arm = [1.5, 0.2, -0.5];
        let omega_b = Vector3::new(0.01, -0.02, 0.05);

        // Run apply_nhc to get analytical H
        // We will intercept the H matrix by making a dummy updater update.
        // Actually apply_nhc calls updater::update. It's easier to just copy the analytical Jacobian code here.
        let r_e_b = state.attitude.inverse().to_rotation_matrix();
        let v_b_imu = r_e_b * state.velocity;
        let r_b_v = nalgebra::Rotation3::from_euler_angles(imu_mounting_angles.unwrap()[0], imu_mounting_angles.unwrap()[1], imu_mounting_angles.unwrap()[2]);
        
        let n = state.covariance.nrows();
        let mut h_ana = DMatrix::<f64>::zeros(2, n);
        let rot = r_b_v * r_e_b;
        let dr_dv = rot.matrix();
        h_ana[(0, 3)] = dr_dv[(1, 0)]; h_ana[(0, 4)] = dr_dv[(1, 1)]; h_ana[(0, 5)] = dr_dv[(1, 2)];
        h_ana[(1, 3)] = dr_dv[(2, 0)]; h_ana[(1, 4)] = dr_dv[(2, 1)]; h_ana[(1, 5)] = dr_dv[(2, 2)];
        
        let v_b_skew = nalgebra::Matrix3::new(
            0.0, -v_b_imu.z,  v_b_imu.y,
            v_b_imu.z,  0.0, -v_b_imu.x,
           -v_b_imu.y,  v_b_imu.x,  0.0
        );
        let dr_dpsi = r_b_v.matrix() * v_b_skew * r_e_b.matrix();
        h_ana[(0, 6)] = dr_dpsi[(1, 0)]; h_ana[(0, 7)] = dr_dpsi[(1, 1)]; h_ana[(0, 8)] = dr_dpsi[(1, 2)];
        h_ana[(1, 6)] = dr_dpsi[(2, 0)]; h_ana[(1, 7)] = dr_dpsi[(2, 1)]; h_ana[(1, 8)] = dr_dpsi[(2, 2)];
        
        let l_b = Vector3::new(imu_to_nhc_lever_arm[0], imu_to_nhc_lever_arm[1], imu_to_nhc_lever_arm[2]);
        let l_b_skew = nalgebra::Matrix3::new(
            0.0, -l_b.z,  l_b.y,
            l_b.z,  0.0, -l_b.x,
           -l_b.y,  l_b.x,  0.0
        );
        let dr_dbg = r_b_v.matrix() * l_b_skew;
        h_ana[(0, 12)] = dr_dbg[(1, 0)]; h_ana[(0, 13)] = dr_dbg[(1, 1)]; h_ana[(0, 14)] = dr_dbg[(1, 2)];
        h_ana[(1, 12)] = dr_dbg[(2, 0)]; h_ana[(1, 13)] = dr_dbg[(2, 1)]; h_ana[(1, 14)] = dr_dbg[(2, 2)];

        let h_num = compute_numerical_jacobian_nhc(&state, &imu_mounting_angles, &imu_to_nhc_lever_arm, &omega_b, 1e-6);

        let diff = (h_ana.clone() - h_num.clone()).abs().max();
        println!("Max diff: {}", diff);
        assert!(diff < 1e-5, "NHC Jacobian verification failed!");
    }
}
