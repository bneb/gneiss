use crate::engine::fgo::factors::skew_symmetric;
use crate::engine::updater;
use crate::filter::RtkState;
use nalgebra::{DMatrix, DVector, Vector3};

/// Assigns a 3x3 matrix block to a 2xN measurement Jacobian matrix
fn assign_jacobian_block(h: &mut DMatrix<f64>, col: usize, m: nalgebra::Matrix3<f64>) {
    for i in 0..3 {
        h[(0, col + i)] = m[(1, i)];
        h[(1, col + i)] = m[(2, i)];
    }
}

/// Applies Non-Holonomic Constraints (NHC) to the EKF state.
/// This assumes the vehicle's lateral and vertical velocity in the body frame is zero.
pub fn apply_nhc(
    state: &mut RtkState,
    sigma_lateral: f64,
    sigma_vertical: f64,
    imu_to_nhc_lever_arm: &[f64; 3],
    omega_b: &Vector3<f64>,
    tuning: &crate::engine::config::EkfTuningConfig,
) -> Result<(), &'static str> {
    let r_e_b = state.attitude.to_rotation_matrix().transpose();

    // Convert IMU velocity and NHC lever arm to vehicle frame
    let v_b_imu = r_e_b * state.velocity;
    let l_b = Vector3::from_column_slice(imu_to_nhc_lever_arm);

    // Velocity at NHC point in vehicle frame (v_v = v_b_imu + omega x l_b)
    let v_v = v_b_imu + omega_b.cross(&l_b);
    let z = DVector::from_column_slice(&[-v_v.y, -v_v.z]);

    let mut h = DMatrix::<f64>::zeros(2, state.covariance.nrows());

    // Jacobian w.r.t velocity in ECEF (dr_dv = R_e_b)
    assign_jacobian_block(&mut h, 3, *r_e_b.matrix());

    // Jacobian w.r.t attitude (dr_dpsi = [v_b_imu x] * R_e_b)
    let v_b_skew = skew_symmetric(&v_b_imu);
    assign_jacobian_block(&mut h, 6, v_b_skew * r_e_b.matrix());

    // Jacobian w.r.t gyro bias (dr_dbg = [l_b x])
    let l_b_skew = skew_symmetric(&l_b);
    assign_jacobian_block(&mut h, 12, l_b_skew);

    let r = DMatrix::from_diagonal(&DVector::from_column_slice(&[
        sigma_lateral * sigma_lateral,
        sigma_vertical * sigma_vertical,
    ]));

    updater::update::<crate::engine::updater_math::TightCoupling>(
        state, &z, &h, &r, 1e9, None, tuning,
    )
    .map_err(|_| "NHC update failed")?;
    Ok(())
}

pub fn apply_zupt(
    state: &mut RtkState,
    sigma: f64,
    tuning: &crate::engine::config::EkfTuningConfig,
) -> Result<(), &'static str> {
    let z = DVector::from_column_slice(&[-state.velocity.x, -state.velocity.y, -state.velocity.z]);
    let n = state.covariance.nrows();
    let mut h = DMatrix::<f64>::zeros(3, n);
    for i in 0..3 {
        h[(i, 3 + i)] = 1.0;
    }

    let r = DMatrix::from_diagonal(&DVector::from_element(3, sigma * sigma));
    updater::update::<crate::engine::updater_math::TightCoupling>(
        state, &z, &h, &r, 1e9, None, tuning,
    )
    .map_err(|_| "ZUPT update failed")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::RtkState;
    use gneiss_core::time::GpsTime;
    use nalgebra::{UnitQuaternion, Vector3};

    fn compute_numerical_jacobian_nhc(
        state: &RtkState,
        _imu_mounting_angles: &Option<[f64; 3]>,
        imu_to_nhc_lever_arm: &[f64; 3],
        omega_b: &Vector3<f64>,
        epsilon: f64,
    ) -> DMatrix<f64> {
        let n = state.covariance.nrows();
        let mut h_num = DMatrix::zeros(2, n);

        for j in 0..n {
            let mut state_pos = state.clone();
            let mut state_neg = state.clone();

            let mut omega_b_pos = *omega_b;
            let mut omega_b_neg = *omega_b;

            if (6..=8).contains(&j) {
                let mut dpsi_pos = Vector3::zeros();
                dpsi_pos[j - 6] = epsilon;
                let dq_pos = UnitQuaternion::from_scaled_axis(dpsi_pos);
                state_pos.attitude = dq_pos * state_pos.attitude;

                let mut dpsi_neg = Vector3::zeros();
                dpsi_neg[j - 6] = -epsilon;
                let dq_neg = UnitQuaternion::from_scaled_axis(dpsi_neg);
                state_neg.attitude = dq_neg * state_neg.attitude;
            } else if (3..6).contains(&j) {
                state_pos.velocity[j - 3] += epsilon;
                state_neg.velocity[j - 3] -= epsilon;
            } else if (12..15).contains(&j) {
                state_pos.gyro_bias[j - 12] += epsilon;
                state_neg.gyro_bias[j - 12] -= epsilon;
                // omega_b = gyro - bg
                omega_b_pos[j - 12] -= epsilon;
                omega_b_neg[j - 12] += epsilon;
            }

            let get_z = |s: &RtkState, ob: &Vector3<f64>| -> DVector<f64> {
                let r_e_b = s.attitude.inverse().to_rotation_matrix();
                let v_b_imu = r_e_b * s.velocity;
                let l_b = Vector3::new(
                    imu_to_nhc_lever_arm[0],
                    imu_to_nhc_lever_arm[1],
                    imu_to_nhc_lever_arm[2],
                );
                let v_v = v_b_imu + ob.cross(&l_b);
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
        let mut state = RtkState::new(
            GpsTime::new(0, 0.0),
            gneiss_core::coords::Coordinate::new(
                Vector3::new(1.0, 2.0, 3.0),
                gneiss_core::coords::Datum::WGS84,
                gneiss_core::coords::Frame::ECEF,
                GpsTime::new(0, 0.0),
            ),
            1.0,
        );
        state.velocity = Vector3::new(10.0, 20.0, 30.0);
        state.attitude = UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3);

        let imu_to_nhc_lever_arm = [1.5, 0.2, -0.5];
        let omega_b = Vector3::new(0.01, -0.02, 0.05);

        let r_e_b = state.attitude.to_rotation_matrix().transpose();
        let v_b_imu = r_e_b * state.velocity;

        let n = state.covariance.nrows();
        let mut h_ana = DMatrix::<f64>::zeros(2, n);
        let dr_dv = r_e_b.matrix();
        h_ana[(0, 3)] = dr_dv[(1, 0)];
        h_ana[(0, 4)] = dr_dv[(1, 1)];
        h_ana[(0, 5)] = dr_dv[(1, 2)];
        h_ana[(1, 3)] = dr_dv[(2, 0)];
        h_ana[(1, 4)] = dr_dv[(2, 1)];
        h_ana[(1, 5)] = dr_dv[(2, 2)];

        let v_b_skew = nalgebra::Matrix3::new(
            0.0, -v_b_imu.z, v_b_imu.y, v_b_imu.z, 0.0, -v_b_imu.x, -v_b_imu.y, v_b_imu.x, 0.0,
        );
        let dr_dpsi = v_b_skew * r_e_b.matrix();
        h_ana[(0, 6)] = dr_dpsi[(1, 0)];
        h_ana[(0, 7)] = dr_dpsi[(1, 1)];
        h_ana[(0, 8)] = dr_dpsi[(1, 2)];
        h_ana[(1, 6)] = dr_dpsi[(2, 0)];
        h_ana[(1, 7)] = dr_dpsi[(2, 1)];
        h_ana[(1, 8)] = dr_dpsi[(2, 2)];

        let l_b = Vector3::new(
            imu_to_nhc_lever_arm[0],
            imu_to_nhc_lever_arm[1],
            imu_to_nhc_lever_arm[2],
        );
        let l_b_skew =
            nalgebra::Matrix3::new(0.0, -l_b.z, l_b.y, l_b.z, 0.0, -l_b.x, -l_b.y, l_b.x, 0.0);
        let dr_dbg = l_b_skew;
        h_ana[(0, 12)] = dr_dbg[(1, 0)];
        h_ana[(0, 13)] = dr_dbg[(1, 1)];
        h_ana[(0, 14)] = dr_dbg[(1, 2)];
        h_ana[(1, 12)] = dr_dbg[(2, 0)];
        h_ana[(1, 13)] = dr_dbg[(2, 1)];
        h_ana[(1, 14)] = dr_dbg[(2, 2)];

        let h_num =
            compute_numerical_jacobian_nhc(&state, &None, &imu_to_nhc_lever_arm, &omega_b, 1e-6);

        let diff = (h_ana.clone() - h_num.clone()).abs().max();
        println!("Max diff: {}", diff);
        assert!(diff < 1e-5, "NHC Jacobian verification failed!");
    }

    #[test]
    fn test_apply_nhc_updates_state() {
        let mut state = RtkState::new(
            GpsTime::new(0, 0.0),
            gneiss_core::coords::Coordinate::new(
                Vector3::new(1.0, 2.0, 3.0),
                gneiss_core::coords::Datum::WGS84,
                gneiss_core::coords::Frame::ECEF,
                GpsTime::new(0, 0.0),
            ),
            1.0,
        );
        state.velocity = Vector3::new(10.0, 5.0, -2.0);
        state.attitude = UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3);
        let tuning = crate::engine::config::EkfTuningConfig::default();

        // Use distinct initial covariance values to catch indexing mutants
        for i in 0..state.covariance.nrows() {
            state.covariance[(i, i)] = (i as f64 + 1.0) * 0.1;
        }

        let initial_cov = state.covariance.clone();
        let imu_to_nhc = [1.0, 2.0, 3.0];
        let omega = Vector3::new(0.01, -0.02, 0.05);

        apply_nhc(&mut state, 0.1, 0.1, &imu_to_nhc, &omega, &tuning).unwrap();

        assert!(state.covariance[(3, 3)] < initial_cov[(3, 3)]);
        assert!(state.covariance[(4, 4)] < initial_cov[(4, 4)]);
        assert!(state.covariance[(5, 5)] < initial_cov[(5, 5)]);
        assert!(state.covariance[(6, 6)] < initial_cov[(6, 6)]);
        assert!(state.covariance[(12, 12)] < initial_cov[(12, 12)]); // gyro bias
    }

    #[test]
    fn test_apply_zupt() {
        let mut state = RtkState::new(
            GpsTime::new(0, 0.0),
            gneiss_core::coords::Coordinate::new(
                Vector3::new(1.0, 2.0, 3.0),
                gneiss_core::coords::Datum::WGS84,
                gneiss_core::coords::Frame::ECEF,
                GpsTime::new(0, 0.0),
            ),
            1.0,
        );
        state.velocity = Vector3::new(10.0, 5.0, -2.0);
        let tuning = crate::engine::config::EkfTuningConfig::default();
        let initial_cov = state.covariance.clone();

        apply_zupt(&mut state, 0.1, &tuning).unwrap();

        assert!(state.covariance[(3, 3)] < initial_cov[(3, 3)]);
        assert!(state.covariance[(4, 4)] < initial_cov[(4, 4)]);
        assert!(state.covariance[(5, 5)] < initial_cov[(5, 5)]);
    }
}
