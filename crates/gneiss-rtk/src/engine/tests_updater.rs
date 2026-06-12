#[cfg(test)]
mod tests {
    use crate::filter::RtkState;
    use crate::engine::updater::UpdateError;
    use gneiss_core::time::GpsTime;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use nalgebra::{Vector3, DVector, DMatrix};

    #[test]
    fn test_robust_outlier_rejection_multipath() {
        // Setup a state with a relatively constrained position and unconstrained clock
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.5); // P_pos = 2.5 (SPP variance)
        state.covariance[(15, 15)] = 1e6; // Clock bias highly uncertain

        // Create 6 satellites: 5 good, 1 bad (500m multipath)
        let mut z = DVector::zeros(6);
        let mut h = DMatrix::zeros(6, state.covariance.ncols());
        let mut r = DMatrix::zeros(6, 6);

        // Good satellites (residuals near 0)
        for i in 0..5 {
            z[i] = 1.0; // Small noise
            h[(i, 0)] = 0.5; // Arbitrary geometry
            h[(i, 1)] = 0.5;
            h[(i, 2)] = 0.5;
            h[(i, 15)] = 1.0; // Clock bias
            r[(i, i)] = 25.0; // PR variance
        }

        // Bad satellite (500m multipath)
        z[5] = 500.0;
        h[(5, 0)] = -0.5;
        h[(5, 1)] = -0.5;
        h[(5, 2)] = -0.5;
        h[(5, 15)] = 1.0;
        r[(5, 5)] = 25.0;

        let result = crate::engine::updater::update(
            &mut state, 
            &z, 
            &h, 
            &r, 
            15.0, // max_innovation
            None,
            false
        );

        assert!(result.is_ok());
        let final_valid = result.unwrap();

        // The worst outlier (index 5) should be rejected.
        assert!(!final_valid.contains(&5usize), "Outlier was not rejected! final_valid: {:?}", final_valid);
        assert_eq!(final_valid.len(), 5usize, "Good measurements were rejected!");
    }

    #[test]
    fn test_robust_outlier_rejection_clock_jump() {
        // Setup a state with a relatively constrained position and unconstrained clock
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.5);
        state.covariance[(15, 15)] = 1e6; // Clock bias highly uncertain

        // Create 6 satellites: ALL jump by 10,000m (true clock jump)
        let mut z = DVector::zeros(6);
        let mut h = DMatrix::zeros(6, state.covariance.ncols());
        let mut r = DMatrix::zeros(6, 6);

        for i in 0..6 {
            z[i] = 10000.0; // Large clock jump
            h[(i, 0)] = 0.5;
            h[(i, 1)] = 0.5;
            h[(i, 2)] = 0.5;
            h[(i, 15)] = 1.0; // Clock bias
            r[(i, i)] = 25.0; // PR variance
        }

        let result = crate::engine::updater::update(
            &mut state, 
            &z, 
            &h, 
            &r, 
            15.0, 
            None,
            false
        );

        assert!(result.is_ok());
        let final_valid = result.unwrap();

        // No measurements should be rejected because they are consistent with a clock jump.
        assert_eq!(final_valid.len(), 6usize, "True clock jump was rejected!");
    }

    #[test]
    fn test_apply_state_correction() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.5);
        
        let mut dx = DVector::zeros(state.covariance.nrows());
        dx[0] = 1.0; dx[1] = -2.0; dx[2] = 3.0; // pos
        dx[3] = 0.1; dx[4] = 0.2; dx[5] = -0.3; // vel
        dx[15] = 100.0; dx[16] = 1.5; // clock
        
        crate::engine::updater::apply_state_correction(&mut state, &dx);
        
        assert_eq!(state.position.vector.x, 1.0);
        assert_eq!(state.position.vector.y, -2.0);
        assert_eq!(state.position.vector.z, 3.0);
        assert_eq!(state.velocity.x, 0.1);
        assert_eq!(state.velocity.y, 0.2);
        assert_eq!(state.velocity.z, -0.3);
        assert_eq!(state.rcv_clk_bias, 100.0);
        assert_eq!(state.rcv_clk_drift, 1.5);
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
        let mut h = DMatrix::zeros(3, state.covariance.ncols());
        let mut r = DMatrix::from_diagonal(&DVector::from_element(3, 1.0));
        
        // 0: PR - valid (nu=10.0, var=1.0)
        z[0] = 10.0;
        h[(0, 0)] = 1.0;
        
        // 1: Phase - invalid (nu=1000.0) -> threshold is max_inn * 10000, so it might pass if max_inn=15
        z[1] = 500.0; 
        h[(1, 0)] = 1.0;
        
        // 2: Doppler - valid (nu=50.0) -> threshold is max_inn * 1000
        z[2] = 50.0;
        h[(2, 0)] = 1.0;
        
        use gneiss_core::sat::{SatelliteId, Constellation};
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let meas_types = [(sat1, 0), (sat1, 1), (sat1, 3)];
        
        let valid_idx = crate::engine::updater::filter_pre_fit_residuals(&z, &h, &r, &state.covariance, 15.0, Some(&meas_types));
        
        // PR threshold: 15*15 = 225. nu*nu / s_ii = 100 / (2.5 + 1) = 28.5. Valid.
        // Phase threshold: 15*10000 = 150000. nu*nu = 250000. 250000 / 3.5 = 71428. Valid.
        // Let's make phase explicitly invalid by increasing z[1]
        
        // Adjust z[1] to 1000.0 so nu^2 / s_ii = 1e6 / 3.5 = 285k > 150k.
        z[1] = 1000.0;
        let valid_idx = crate::engine::updater::filter_pre_fit_residuals(&z, &h, &r, &state.covariance, 15.0, Some(&meas_types));
        
        assert!(valid_idx.contains(&0));
        assert!(!valid_idx.contains(&1));
        assert!(valid_idx.contains(&2));
    }

    #[test]
    fn test_evaluate_post_fit_outliers() {
        let v = DVector::from_vec(vec![10.0, 50.0, 10.0]); // Res: PR=10m, Phase=50m, Dopp=10m/s
        let s = DMatrix::from_diagonal(&DVector::from_vec(vec![1.0, 1.0, 1.0])); // Var=1
        let current_z = DVector::from_vec(vec![10.0, 50.0, 10.0]);
        let current_valid = vec![0, 1, 2];
        
        use gneiss_core::sat::{SatelliteId, Constellation};
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let meas_types = [(sat1, 0), (sat1, 1), (sat1, 3)]; // PR, Phase, Doppler
        
        // threshold PR=15, Phase=5.0, Doppler=30
        let (worst_idx, max_ratio) = crate::engine::updater::evaluate_post_fit_outliers(
            &v, &s, &current_z, &current_valid, Some(&meas_types), 15.0, false
        );
        
        // v[0] = 10, thresh=15 -> OK
        // v[1] = 50, thresh=5 -> OUTLIER. Ratio = 50/1 = 50
        // v[2] = 10, thresh=30 -> OK
        
        assert_eq!(worst_idx, Some(1));
        assert_eq!(max_ratio, 50.0);
    }
}
