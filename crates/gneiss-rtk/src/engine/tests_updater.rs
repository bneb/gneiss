#[cfg(test)]
mod tests {
    use crate::filter::RtkState;
    use crate::engine::updater::UpdateError;
    use gneiss_core::time::GpsTime;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use nalgebra::{Vector3, DVector, DMatrix};


    #[test]
    fn test_apply_state_correction() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.5);
        
        // Construct dx with a unique value at every core state index.
        // State layout: [0-2]=pos, [3-5]=vel, [6-8]=att, [9-11]=accel_bias,
        // [12-14]=gyro_bias, [15]=clk_bias, [16]=isb_glo, [17]=isb_gal,
        // [18]=isb_bds, [19]=clk_drift, [20]=zwd
        let mut dx = DVector::zeros(state.covariance.nrows());
        dx[0] = 1.0; dx[1] = -2.0; dx[2] = 3.0;     // position
        dx[3] = 0.1; dx[4] = 0.2; dx[5] = -0.3;      // velocity
        // dx[6..8] = 0 (skip attitude for this test)
        dx[9] = 0.01; dx[10] = 0.02; dx[11] = 0.03;  // accel bias
        dx[12] = 0.04; dx[13] = 0.05; dx[14] = 0.06; // gyro bias
        dx[15] = 100.0;                                // clock bias
        dx[16] = 10.0;                                 // ISB GLONASS
        dx[17] = 20.0;                                 // ISB Galileo
        dx[18] = 30.0;                                 // ISB BeiDou
        dx[19] = 1.5;                                  // clock drift
        dx[20] = 0.05;                                 // ZWD
        
        crate::engine::updater::apply_state_correction(&mut state, &dx);
        
        // Position
        assert_eq!(state.position.vector.x, 1.0);
        assert_eq!(state.position.vector.y, -2.0);
        assert_eq!(state.position.vector.z, 3.0);
        // Velocity
        assert_eq!(state.velocity.x, 0.1);
        assert_eq!(state.velocity.y, 0.2);
        assert_eq!(state.velocity.z, -0.3);
        // IMU biases
        assert_eq!(state.accel_bias.x, 0.01);
        assert_eq!(state.accel_bias.y, 0.02);
        assert_eq!(state.accel_bias.z, 0.03);
        assert_eq!(state.gyro_bias.x, 0.04);
        assert_eq!(state.gyro_bias.y, 0.05);
        assert_eq!(state.gyro_bias.z, 0.06);
        // Clock and ISBs
        assert_eq!(state.rcv_clk_bias, 100.0);
        assert_eq!(state.isb_glo, 10.0);
        assert_eq!(state.isb_gal, 20.0);
        assert_eq!(state.isb_bds, 30.0);
        assert_eq!(state.rcv_clk_drift, 1.5);
        assert!((state.zwd - 0.15).abs() < 1e-14); // 0.1 initial + 0.05 correction
    }

    #[test]
    fn test_apply_state_correction_attitude_global_frame() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.5);
        
        // Initial attitude: 90 degrees around Z
        let q_init = nalgebra::UnitQuaternion::from_axis_angle(&nalgebra::Vector3::z_axis(), core::f64::consts::FRAC_PI_2);
        state.attitude = q_init;
        
        let mut dx = DVector::zeros(state.covariance.nrows());
        // Global frame error state: rotation around X axis in ECEF frame
        dx[6] = 0.1; 
        
        crate::engine::updater::apply_state_correction(&mut state, &dx);
        
        // The rotation should be applied in the global (ECEF) frame.
        // q_new = q_roll_global * q_init
        
        let v_b = Vector3::new(0.0, 1.0, 0.0);
        // q_init rotates (0,1,0) to (-1,0,0)
        // Then roll around X by 0.1 leaves (-1,0,0) unchanged!
        let v_e = state.attitude * v_b;
        
        assert!((v_e.x - -1.0).abs() < 1e-6);
        assert!((v_e.y - 0.0).abs() < 1e-6);
        assert!((v_e.z - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_apply_joseph_covariance_update() {
        let state_cov = DMatrix::from_diagonal(&DVector::from_element(3, 10.0));
        let k = DMatrix::from_element(3, 1, 0.5);
        let h = DMatrix::from_element(1, 3, 1.0);
        let r = DMatrix::from_element(1, 1, 1.0);
        
        let p_new = crate::engine::updater::apply_joseph_covariance_update(&state_cov, &k, &h, &r);
        
        // Ensure symmetric
        assert_eq!(p_new[(0, 1)], p_new[(1, 0)]);
        assert_eq!(p_new[(1, 2)], p_new[(2, 1)]);
        assert_eq!(p_new[(0, 2)], p_new[(2, 0)]);
    }

    #[test]
    fn test_filter_pre_fit_residuals() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 2.5);
        
        let mut z = DVector::zeros(3);
        let h = DMatrix::from_fn(3, state.covariance.ncols(), |r, c| if r == c { 1.0 } else { 0.0 });
        let r = DMatrix::from_diagonal(&DVector::from_element(3, 1.0));
        
        // s_ii = P(0,0) + R(0,0) = 2.5 + 1.0 = 3.5 for all measurements
        
        // PR (type 0): z=10.0, threshold = max_inn^2 = 225, nu^2/s_ii = 28.5 → Valid
        z[0] = 10.0;
        // Phase (type 1): z=3.0, threshold = CP_CHI2 = 100, nu^2/s_ii = 2.57 → Valid
        z[1] = 3.0;
        // Doppler (type 3): z=3.0, threshold = DOP_CHI2 = 50, nu^2/s_ii = 2.57 → Valid
        z[2] = 3.0;
        
        use gneiss_core::sat::{SatelliteId, Constellation};
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let meas_types = [(sat1, 0), (sat1, 1), (sat1, 3)];
        
        let _valid_idx = crate::engine::updater::filter_pre_fit_residuals(&z, &h, &r, &state.covariance, 15.0, Some(&meas_types));
        
        // Now make phase invalid: z=20.0, nu^2/s_ii = 400/3.5 = 114 > 100 → Invalid
        z[1] = 20.0;
        let valid_idx = crate::engine::updater::filter_pre_fit_residuals(&z, &h, &r, &state.covariance, 15.0, Some(&meas_types));
        
        assert!(valid_idx.contains(&0), "PR should pass");
        assert!(!valid_idx.contains(&1), "Phase should be rejected");
        assert!(valid_idx.contains(&2), "Doppler should pass");
    }

    #[test]
    fn test_update_loosely_coupled_huber() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos.clone(), 2.5);
        let mut gnss_state = RtkState::new(time, pos.clone(), 2.5);
        
        // Give gnss_state a position error
        gnss_state.position.vector.x += 15.0; // Huge error
        
        // High Mahalanobis sq will be generated.
        let tuning = crate::engine::config::EkfTuningConfig {
            loosely_coupled_mahalanobis_sq: 10.0,
            huber_threshold_loosely: 3.0,
            ..Default::default()
        };
        
        let lever_arm = Vector3::zeros();
        let omega_b = Vector3::zeros();
        
        let res = crate::engine::updater::update_loosely_coupled(&mut state, &gnss_state, lever_arm, omega_b, &tuning);
        assert!(res.is_ok(), "Huber scaling should prevent rejection of the huge error");
        
        // state position should be updated but not fully 15m due to variance inflation
        assert!(state.position.vector.x > 0.0 && state.position.vector.x < 15.0);
    }

}
