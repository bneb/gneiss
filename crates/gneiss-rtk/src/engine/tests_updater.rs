#[cfg(test)]
mod tests {
    use crate::filter::RtkState;

    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, DVector, Vector3};

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
        dx[0] = 1.0;
        dx[1] = -2.0;
        dx[2] = 3.0; // position
        dx[3] = 0.1;
        dx[4] = 0.2;
        dx[5] = -0.3; // velocity
                      // dx[6..8] = 0 (skip attitude for this test)
        dx[9] = 0.01;
        dx[10] = 0.02;
        dx[11] = 0.03; // accel bias
        dx[12] = 0.04;
        dx[13] = 0.05;
        dx[14] = 0.06; // gyro bias
        dx[15] = 100.0; // clock bias
        dx[16] = 10.0; // ISB GLONASS
        dx[17] = 20.0; // ISB Galileo
        dx[18] = 30.0; // ISB BeiDou
        dx[19] = 1.5; // clock drift
        dx[20] = 0.05; // ZWD

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
        let q_init = nalgebra::UnitQuaternion::from_axis_angle(
            &nalgebra::Vector3::z_axis(),
            core::f64::consts::FRAC_PI_2,
        );
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
    fn test_apply_state_correction_mutants_bounds() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.5);

        // Add an ambiguity so state.ambiguities.len() > 0
        state.add_ambiguity(
            gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
            1,
            5.3,
            1.0,
        );

        // Pass a dx of EXACTLY CORE_STATE_SIZE.
        // If the code had `dx.len() >= CORE_STATE_SIZE` it would enter the loop and try to access dx[21] which panics.
        let dx = DVector::zeros(crate::filter::CORE_STATE_SIZE);

        // This should run without panicking
        crate::engine::updater::apply_state_correction(&mut state, &dx);

        // Also test the apply_imu_and_clock_correction bound which requires dx.len() >= 15.
        // If it was changed to > 15, then length 15 dx wouldn't trigger it. Oh wait, apply_imu_and_clock_correction has `if dx.len() > 15`.
        // To kill the `if dx.len() > 15` mutated to `>=`, we need dx.len() == 15 exactly.
        let dx15 = DVector::zeros(15);
        crate::engine::updater::apply_state_correction(&mut state, &dx15);

        let mut dx_full = DVector::zeros(crate::filter::CORE_STATE_SIZE + 1);
        dx_full[crate::filter::CORE_STATE_SIZE] = 1.5;
        crate::engine::updater::apply_state_correction(&mut state, &dx_full);
        assert!((state.ambiguities[0] - 6.8).abs() < 1e-10); // 5.3 + 1.5 = 6.8

        // Test d_theta.norm() > 1e-10 mutants
        let mut dx_rot = DVector::zeros(crate::filter::CORE_STATE_SIZE);
        dx_rot[6] = 2e-10; // norm > 1e-10
        let orig_att = state.attitude.clone();
        crate::engine::updater::apply_state_correction(&mut state, &dx_rot);
        assert!((state.attitude.i - orig_att.i).abs() > 0.0); // should change

        let mut dx_rot_small = DVector::zeros(crate::filter::CORE_STATE_SIZE);
        dx_rot_small[6] = 0.5e-10; // norm < 1e-10
        let orig_att = state.attitude.clone();
        crate::engine::updater::apply_state_correction(&mut state, &dx_rot_small);
        assert_eq!(state.attitude, orig_att); // should not change
    }

    #[test]
    fn test_apply_joseph_covariance_update() {
        let state_cov = DMatrix::from_diagonal(&DVector::from_element(3, 10.0));
        let k = DMatrix::from_element(3, 1, 0.5);
        let h = DMatrix::from_element(1, 3, 1.0);
        let r = DMatrix::from_element(1, 1, 1.0);

        let p_new = crate::math::covariance::apply_joseph_covariance_update(&state_cov, &k, &h, &r);

        // Ensure symmetric
        assert_eq!(p_new[(0, 1)], p_new[(1, 0)]);
        assert_eq!(p_new[(1, 2)], p_new[(2, 1)]);
        assert_eq!(p_new[(0, 2)], p_new[(2, 0)]);

        // Ensure correct values to kill mutants.
        // P = 10 * I. KH is all 0.5. I - KH has 0.5 on diag, -0.5 elsewhere.
        // (I - KH) * P * (I - KH)^T:
        // row 0: [0.5, -0.5, -0.5] * 10 = [5, -5, -5]
        // dot product with itself: 25 + 25 + 25 = 75
        // off-diagonals: [5, -5, -5] dot [-5, 0.5*10, -5] = -25 - 25 + 25 = -25.
        // K R K^T = 0.25 everywhere.
        // So diag = 7.75, off-diag = -2.25.
        assert!((p_new[(0, 0)] - 7.75).abs() < 1e-6);
        assert!((p_new[(0, 1)] - (-2.25)).abs() < 1e-6);
    }

    #[test]
    fn test_filter_pre_fit_residuals() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 2.5);

        let mut z = DVector::zeros(3);
        let h = DMatrix::from_fn(
            3,
            state.covariance.ncols(),
            |r, c| if r == c { 1.0 } else { 0.0 },
        );
        let r = DMatrix::from_diagonal(&DVector::from_element(3, 1.0));

        // s_ii = P(0,0) + R(0,0) = 2.5 + 1.0 = 3.5 for all measurements

        // PR (type 0): z=10.0, threshold = max_inn^2 = 225, nu^2/s_ii = 28.5 → Valid
        z[0] = 10.0;
        // Phase (type 1): z=3.0, threshold = CP_CHI2 = 100, nu^2/s_ii = 2.57 → Valid
        z[1] = 3.0;
        // Doppler (type 3): z=3.0, threshold = DOP_CHI2 = 50, nu^2/s_ii = 2.57 → Valid
        z[2] = 3.0;

        use gneiss_core::sat::{Constellation, SatelliteId};
        let sat1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let meas_types = [(sat1, 0), (sat1, 1), (sat1, 3)];

        let _valid_idx = crate::engine::updater_math::filter_pre_fit_residuals::<
            crate::engine::updater_math::LooseCoupling,
        >(&z, &h, &r, &state.covariance, 15.0, Some(&meas_types));

        // Now make phase invalid: z=20.0, nu^2/s_ii = 400/3.5 = 114 > 100 → Invalid
        let mut z = z.clone();
        z[1] = 20.0;
        let valid_idx = crate::engine::updater_math::filter_pre_fit_residuals::<
            crate::engine::updater_math::LooseCoupling,
        >(&mut z, &h, &r, &state.covariance, 15.0, Some(&meas_types));

        assert!(valid_idx.contains(&0), "PR should pass");
        assert!(!valid_idx.contains(&1), "Phase should be rejected");
        assert!(valid_idx.contains(&2), "Doppler should pass");
    }

    #[test]
    fn test_fix_and_hold_updates_imu_states() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);

        // Add two ambiguities with known float values
        state.add_ambiguity(
            gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
            1,
            5.3,
            1.0,
        );
        state.add_ambiguity(
            gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 2,
            },
            1,
            3.8,
            1.0,
        );

        // Set non-trivial cross-covariance between position and ambiguity states
        let amb_start = crate::filter::CORE_STATE_SIZE;
        state.covariance[(0, amb_start)] = 0.5;
        state.covariance[(amb_start, 0)] = 0.5;
        // Cross-covariance between attitude and ambiguity
        state.covariance[(6, amb_start)] = 0.2;
        state.covariance[(amb_start, 6)] = 0.2;
        // Cross-covariance between gyro bias and ambiguity
        state.covariance[(12, amb_start)] = 0.1;
        state.covariance[(amb_start, 12)] = 0.1;

        let initial_att = state.attitude;
        let initial_gyro_bias = state.gyro_bias;

        // Integer ambiguities: 5.0 and 4.0
        let z_dd = DVector::from_vec(vec![5.0 - 4.0]); // DD = amb[0] - amb[1] = 1.0
                                                       // D matrix maps DD to SD: row selects amb[0] - amb[1]
        let n_cols = state.covariance.ncols();
        let mut d_full = DMatrix::zeros(1, n_cols);
        d_full[(0, amb_start)] = 1.0;
        d_full[(0, amb_start + 1)] = -1.0;

        let var = 0.001; // Tight fix variance

        let res = crate::engine::updater::apply_fix_and_hold(&mut state, &z_dd, &d_full, var);
        assert!(res.is_ok());

        // The key assertion: attitude and gyro bias SHOULD be updated
        // because the covariance has cross-correlation between these states
        // and the ambiguity states. Before the fix, these were zeroed.
        let att_changed = (state.attitude.quaternion() - initial_att.quaternion()).norm() > 1e-15
            || state.gyro_bias != initial_gyro_bias;

        // With cross-covariance, the fix should propagate to IMU states
        assert!(
            att_changed,
            "Fix-and-hold must update IMU states via cross-covariance"
        );

        // D_full is [1, -1] for ambiguities.
        // H_P = D * P. P has 0.5 for pos-amb cross, and 0.2 for att-amb cross.
        // D * P will yield some specific H_P vector.
        // By checking the exact output, we kill math operator mutants.
        // E.g. k[(i, j)] *= 0.1 damping.
        assert!((state.gyro_bias.x - initial_gyro_bias.x).abs() > 1e-10);

        // dx = k * v = (0.1 / 2.001) * 0.5? Wait, v = 1.0.
        // No, v = z_dd - d_full * a_sd = 1.0 - (5.3 - 3.8) = 1.0 - 1.5 = -0.5
        // dx = (0.1 / 2.001) * (-0.5) = -0.024987506246876564
        let expected_gyro_x = -0.024987506246876564;
        assert!((state.gyro_bias.x - expected_gyro_x).abs() < 1e-5);
    }

    #[test]
    fn test_fix_and_hold_covariance_reduces() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);

        state.add_ambiguity(
            gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
            1,
            5.3,
            2.0,
        );
        state.add_ambiguity(
            gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 2,
            },
            1,
            3.8,
            2.0,
        );

        let amb_start = crate::filter::CORE_STATE_SIZE;
        let pre_var_0 = state.covariance[(amb_start, amb_start)];
        let pre_var_1 = state.covariance[(amb_start + 1, amb_start + 1)];

        let z_dd = DVector::from_vec(vec![5.0 - 4.0]);
        let n_cols = state.covariance.ncols();
        let mut d_full = DMatrix::zeros(1, n_cols);
        d_full[(0, amb_start)] = 1.0;
        d_full[(0, amb_start + 1)] = -1.0;

        crate::engine::updater::apply_fix_and_hold(&mut state, &z_dd, &d_full, 0.001).unwrap();

        // Ambiguity variance should decrease after fix
        let post_var_0 = state.covariance[(amb_start, amb_start)];
        let post_var_1 = state.covariance[(amb_start + 1, amb_start + 1)];

        assert!(
            post_var_0 < pre_var_0,
            "Ambiguity variance should decrease after fix: {} >= {}",
            post_var_0,
            pre_var_0
        );
        assert!(
            post_var_1 < pre_var_1,
            "Ambiguity variance should decrease after fix: {} >= {}",
            post_var_1,
            pre_var_1
        );
    }

    #[test]
    fn test_loosely_coupled_jacobian() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos.clone(), 1.0);
        let gnss_state = RtkState::new(time, pos, 1.0);

        let lever_arm = Vector3::new(0.5, 0.3, 0.1);
        let omega_b = Vector3::new(0.01, 0.02, 0.03);
        let tuning = crate::engine::config::EkfTuningConfig::default();

        let initial_pos = state.position.vector;
        let initial_vel = state.velocity;

        let res = crate::engine::updater::update_loosely_coupled(
            &mut state,
            &gnss_state,
            lever_arm,
            omega_b,
            &tuning,
        );
        assert!(res.is_ok());

        // Loosely coupled innovation has 6 elements (3 position + 3 velocity),
        // and the update should have modified the state from its initial zeros.
        assert_ne!(
            state.position.vector, initial_pos,
            "position should be updated by loosely coupled correction"
        );
        assert_ne!(
            state.velocity, initial_vel,
            "velocity should be updated by loosely coupled correction"
        );
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

        let res = crate::engine::updater::update_loosely_coupled(
            &mut state,
            &gnss_state,
            lever_arm,
            omega_b,
            &tuning,
        );
        assert!(
            res.is_ok(),
            "Huber scaling should prevent rejection of the huge error"
        );

        assert!(state.position.vector.x > 0.1 && state.position.vector.x < 5.0);
    }

    #[test]
    fn test_update_loosely_coupled_exact_math() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos.clone(), 2.0);
        let mut gnss_state = RtkState::new(time, pos.clone(), 2.0);

        // P = 4.0, R = 4.0, S = 8.0. K = P/S = 0.5.
        gnss_state.position.vector.x += 10.0;
        gnss_state.position.vector.y += 10.0;

        let tuning = crate::engine::config::EkfTuningConfig {
            loosely_coupled_mahalanobis_sq: 1000.0, // no rejection
            huber_threshold_loosely: 1000.0,        // no huber scaling
            ..Default::default()
        };

        let lever_arm = Vector3::new(1.0, 0.0, 0.0);
        let omega_b = Vector3::new(0.0, 1.0, 0.0);

        let res = crate::engine::updater::update_loosely_coupled(
            &mut state,
            &gnss_state,
            lever_arm,
            omega_b,
            &tuning,
        );
        assert!(res.is_ok());

        // s = 8.0. If mutant makes it p*r = 16.0, K=4/16=0.25, dx=2.5.
        // True dx = 0.5 * 9.0 = 4.5
        assert!((state.position.vector.x - pos.vector.x).abs() > 1.0);

        // Check jacobian was populated (if state.covariance.nrows() >= 21)
        // With lever arm, the jacobian H will have non-zero attitude mapping, so attitude will be updated.
        assert!(
            state.attitude.vector().norm() > 0.0,
            "Attitude should be updated via jacobian"
        );
    }

    #[test]
    fn test_enforce_symmetry_operator() {
        let mut p = DMatrix::from_diagonal(&DVector::from_element(2, 10.0));
        p[(0, 1)] = 2.0;
        p[(1, 0)] = 4.0;
        // apply_joseph_covariance_update creates p_new = i_kh * p * i_kh.T + k * r * k.T
        // With K=0, H=0, p_new = p.
        let k = DMatrix::zeros(2, 1);
        let h = DMatrix::zeros(1, 2);
        let r = DMatrix::zeros(1, 1);
        let p_new = crate::math::covariance::apply_joseph_covariance_update(&p, &k, &h, &r);

        // (2.0 + 4.0) * 0.5 = 3.0
        assert_eq!(p_new[(0, 1)], 3.0, "Symmetry operator + or *0.5 mutated");
        assert_eq!(p_new[(1, 0)], 3.0, "Symmetry operator + or *0.5 mutated");
    }

    #[test]
    fn test_compute_update_iteration_math() {
        let state_cov = DMatrix::from_diagonal(&DVector::from_element(2, 2.0));
        let h = DMatrix::from_element(1, 2, 1.0);
        let r = DMatrix::from_element(1, 1, 1.0);
        let z = DVector::from_element(1, 5.0);

        let tuning = crate::engine::config::EkfTuningConfig::default();
        let valid = vec![0];

        // s = H P H^T + R = [1 1] [2 0; 0 2] [1; 1] + 1 = 2 + 2 + 1 = 5.
        // S_inv = 1/5 = 0.2.
        // K = P H^T S_inv = [2 0; 0 2] [1; 1] * 0.2 = [0.4; 0.4].
        // dx = K z = [0.4; 0.4] * 5.0 = [2.0; 2.0].
        // v = z - H dx = 5.0 - [1 1] [2.0; 2.0] = 5.0 - 4.0 = 1.0.

        let res = crate::engine::updater::compute_update_iteration::<
            crate::engine::updater::TightCoupling,
        >(&state_cov, &z, &h, &r, &valid, None, 10.0, &tuning)
        .unwrap();

        assert!((res.dx[0] - 2.0).abs() < 1e-10);
        assert!((res.dx[1] - 2.0).abs() < 1e-10);

        // evaluate_post_fit_outliers takes `v`, if v was z + H dx = 5.0 + 4.0 = 9.0
        // then the outlier ratio would be 9.0 / sqrt(5.0) vs 1.0 / sqrt(5.0)
        assert!(
            res.max_outlier_ratio < 1.0,
            "v should be 1.0, ratio 1.0/sqrt(5.0) ~ 0.44. If mutant +, v=9.0, ratio ~4.0"
        );
    }

    #[test]
    fn test_compute_update_iteration_singular_s() {
        let state_cov = DMatrix::from_diagonal(&DVector::from_element(1, f64::INFINITY));
        let h = DMatrix::from_element(1, 1, 1.0);
        let r = DMatrix::from_element(1, 1, 1.0);
        let z = DVector::from_element(1, 5.0);
        let valid = vec![0];
        let tuning = crate::engine::config::EkfTuningConfig::default();

        let res = crate::engine::updater::compute_update_iteration::<
            crate::engine::updater::TightCoupling,
        >(&state_cov, &z, &h, &r, &valid, None, 10.0, &tuning);
        assert!(res.is_err(), "Should return Err for Inf in S");

        let state_cov_nan = DMatrix::from_diagonal(&DVector::from_element(1, f64::NAN));
        let res2 = crate::engine::updater::compute_update_iteration::<
            crate::engine::updater::TightCoupling,
        >(&state_cov_nan, &z, &h, &r, &valid, None, 10.0, &tuning);
        assert!(res2.is_err(), "Should return Err for NaN in S");

        let state_cov_large = DMatrix::from_diagonal(&DVector::from_element(1, 2e15));
        let res3 = crate::engine::updater::compute_update_iteration::<
            crate::engine::updater::TightCoupling,
        >(&state_cov_large, &z, &h, &r, &valid, None, 10.0, &tuning);
        assert!(res3.is_err(), "Should return Err for > 1e15 in S");
    }

    #[test]
    fn test_evaluate_post_fit_outliers_math_operators() {
        let v = DVector::from_element(1, 10.0);
        let s = DMatrix::from_element(1, 1, 4.0); // sqrt is 2.0
        let current_z = DVector::from_element(1, 20.0);
        let current_valid = vec![0];
        use gneiss_core::sat::{Constellation, SatelliteId};
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let meas_types = vec![(sat, 3)]; // Doppler (type 3)
        let mut tuning = crate::engine::config::EkfTuningConfig::default();
        tuning.doppler_outlier_ratio_mult = 2.0;
        tuning.dop_abs_thresh = 50.0;
        let max_innovation = 6.0;

        // thresh = max_innovation * mult = 6.0 * 2.0 = 12.0
        // ratio = v.abs() / s.sqrt() = 10.0 / 2.0 = 5.0
        // Since 5.0 <= 12.0, it should NOT be an outlier.
        let (worst_idx, _, _) = crate::engine::updater_math::evaluate_post_fit_outliers::<
            crate::engine::updater_math::LooseCoupling,
        >(
            &v,
            &s,
            &current_z,
            &current_valid,
            Some(&meas_types),
            max_innovation,
            &tuning,
        );
        assert_eq!(
            worst_idx, None,
            "Doppler mult operator * mutated to / (would yield thresh 3.0, causing outlier)"
        );

        // If we increase ratio to 15.0 by setting v=30.0, it should be an outlier.
        let v2 = DVector::from_element(1, 30.0);
        let (worst_idx2, _, _) = crate::engine::updater_math::evaluate_post_fit_outliers::<
            crate::engine::updater_math::LooseCoupling,
        >(
            &v2,
            &s,
            &current_z,
            &current_valid,
            Some(&meas_types),
            max_innovation,
            &tuning,
        );
        assert_eq!(worst_idx2, Some(0));
    }

    #[test]
    fn test_apply_state_correction_ambiguity_bound() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.5);
        state.add_ambiguity(
            gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
            1,
            5.3,
            1.0,
        );
        let dx = DVector::zeros(crate::filter::CORE_STATE_SIZE);
        crate::engine::updater::apply_state_correction(&mut state, &dx);
        assert!((state.ambiguities[0] - 5.3).abs() < 1e-14);
    }

    #[test]
    fn test_apply_imu_and_clock_correction_bounds() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.5);
        let mut dx = DVector::zeros(15);
        dx[6] = 1e-10;
        let original_att = state.attitude;
        crate::engine::updater::apply_state_correction(&mut state, &dx);
        assert_eq!(state.attitude, original_att);
    }

    #[test]
    fn test_compute_update_iteration_s_threshold() {
        let n_state = 3;
        let n_obs = 1;
        let state_cov = DMatrix::identity(n_state, n_state);
        let current_z = DVector::from_element(n_obs, 1.0);
        let current_h = DMatrix::zeros(n_obs, n_state);
        let mut current_r = DMatrix::identity(n_obs, n_obs);
        current_r[(0, 0)] = 1e15;

        let tuning = crate::engine::config::EkfTuningConfig::default();
        let res = crate::engine::updater::compute_update_iteration::<
            crate::engine::updater_math::LooseCoupling,
        >(
            &state_cov,
            &current_z,
            &current_h,
            &current_r,
            &[0],
            None,
            3.0,
            &tuning,
        );
        assert!(res.is_ok());
        let result = res.unwrap();
        assert_eq!(
            result.dx.len(),
            3,
            "dx should have 3 elements matching n_state"
        );
    }

    #[test]
    fn test_apply_fix_and_hold_mutants() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 2.5);

        let sat1 = gneiss_core::sat::SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: 1,
        };
        let sat2 = gneiss_core::sat::SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: 2,
        };
        state.add_ambiguity(sat1, 1, 10.0, 1.0);
        state.add_ambiguity(sat2, 1, 20.0, 1.0);

        let z_dd = DVector::from_element(1, 1.5);
        let mut d_full = DMatrix::zeros(1, crate::filter::CORE_STATE_SIZE + 2);
        d_full[(0, crate::filter::CORE_STATE_SIZE)] = 1.0;
        d_full[(0, crate::filter::CORE_STATE_SIZE + 1)] = -1.0;

        let var = 0.5;
        let original_cov = state.covariance.clone();
        crate::engine::updater::apply_fix_and_hold(&mut state, &z_dd, &d_full, var).unwrap();

        assert!((state.ambiguities[0] - 10.0).abs() > 1e-5);
        assert!((state.ambiguities[1] - 20.0).abs() > 1e-5);
        assert_ne!(state.covariance, original_cov);
    }
}
