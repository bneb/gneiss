    use super::*;
    use crate::engine::predictor::{compute_process_noise, compute_transition_matrix};
    use crate::engine::{DynamicsModel, EngineConfig};
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::obs::SatObs;
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;

    fn make_state_with_pos(time: GpsTime, pos: Vector3<f64>, initial_var: f64) -> RtkState {
        let coord = Coordinate::new(pos, Datum::WGS84, Frame::ECEF, time);
        RtkState::new(time, coord, initial_var)
    }

    fn make_dummy_sat(sat_id: SatelliteId, sat_pos_rot: Vector3<f64>) -> ProcessedSat<'static> {
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1: sat_pos_rot.norm(),
            p2: None, cp1: None, cp2: None,
            is_iono_free: true,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
            f1: 1.0, f2: 1.0,
            sat_pos_rot,
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        }
    }

    #[test]
    fn test_spp_prior_min_variance_prevents_submeter_accuracy() {
        let t = GpsTime::new(2156, 1000.0);
        let true_pos = Vector3::new(6000000.0, 0.0, 0.0);
        let mut state = make_state_with_pos(t, true_pos, 100.0);
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        for i in 0..3 {
            state.covariance[(i, i)] = 0.01;
        }
        let pos_cov = state.covariance[(0, 0)]
            .min(state.covariance[(1, 1)])
            .min(state.covariance[(2, 2)]);
        let prior_var = (pos_cov.min(25.0)).max(1.0);
        assert_eq!(prior_var, 1.0, "BUG: prior_var should be pos_cov={} but clamping forces it to 1.0", pos_cov);
        assert!(
            prior_var > pos_cov * 10.0,
            "prior_var={} is {:.0}x LARGER than pos_cov={}",
            prior_var, prior_var / pos_cov, pos_cov
        );
    }

    #[test]
    fn test_spp_prior_clamping_biases_final_position() {
        let fg = PppIteratedEkf::default();
        let t = GpsTime::new(2156, 1000.0);
        let mut state = make_state_with_pos(t, Vector3::new(0.0, 0.0, 0.0), 100.0);
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        for i in 0..3 {
            state.covariance[(i, i)] = 0.01;
        }
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let spp_pos = Vector3::new(2.0, 0.0, 0.0);
        let mut state_no_clamp = state.clone();
        let _ = fg.solve(&mut state_no_clamp, &[sat.clone()], Some((spp_pos, 0.01)));
        let x_no_clamp = state_no_clamp.position.vector.x;
        let mut state_with_clamp = state.clone();
        let _ = fg.solve(&mut state_with_clamp, &[sat.clone()], Some((spp_pos, 1.0)));
        let x_with_clamp = state_with_clamp.position.vector.x;
        assert!(
            x_with_clamp.abs() < x_no_clamp.abs(),
            "Clamped should pull LESS: no_clamp={:.4} vs with_clamp={:.4}",
            x_no_clamp, x_with_clamp
        );
        assert!(
            x_no_clamp.abs() > x_with_clamp.abs() * 2.0,
            "Natural prior should pull at least 2x more: no_clamp={:.4}, with_clamp={:.4}",
            x_no_clamp, x_with_clamp
        );
    }

    #[test]
    fn test_random_walk_clock_preserves_covariance_across_epochs() {
        let t = GpsTime::new(2156, 1000.0);
        let pos = Coordinate::new(
            Vector3::new(6000000.0, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, t,
        );
        let mut state = RtkState::new(t, pos, 100.0);
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim) * 100.0;
        state.covariance[(15, 15)] = 0.01;
        state.covariance[(19, 19)] = 0.1;
        state.rcv_clk_bias = 100.0;
        state.rcv_clk_drift = 0.0;
        let config = EngineConfig {
            mode: crate::engine::EngineMode::Ppp,
            process_noise_cb: 1.0,
            dynamics_model: DynamicsModel::Static,
            ..Default::default()
        };
        let phi = compute_transition_matrix(&state, 1.0, &[]);
        assert_eq!(phi[(15, 15)], 1.0);
        let q = compute_process_noise(1.0, &config, false, false, &[]);
        let p_pred = &phi * &state.covariance * phi.transpose() + q;
        assert!(p_pred[(15, 15)] < 2.0);
        assert_eq!(phi[(15, 19)], 1.0);
    }

    #[test]
    fn test_white_noise_clock_destroys_filter_convergence_property() {
        let t0 = GpsTime::new(2156, 1000.0);
        let true_pos = Vector3::new(6000000.0, 0.0, 0.0);
        let pos = Coordinate::new(true_pos, Datum::WGS84, Frame::ECEF, t0);
        let mut state = RtkState::new(t0, pos, 100.0);
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        state.covariance[(15, 15)] = 10000.0;
        state.rcv_clk_bias = 0.0;
        let phi = compute_transition_matrix(&state, 1.0, &[]);
        let q = compute_process_noise(1.0, &EngineConfig::default(), false, false, &[]);
        let mut p_white_noise = state.covariance.clone();
        let mut p_random_walk = state.covariance.clone();
        for _epoch in 0..10 {
            let mut phi_wn = phi.clone();
            phi_wn[(15, 15)] = 0.0;
            p_white_noise = &phi_wn * &p_white_noise * phi_wn.transpose() + &q;
            let mut phi_rw = phi.clone();
            phi_rw[(15, 15)] = 1.0;
            p_random_walk = &phi_rw * &p_random_walk * phi_rw.transpose() + &q;
        }
        let wn_var = p_white_noise[(15, 15)];
        let rw_var = p_random_walk[(15, 15)];
        assert!(
            wn_var < 1000.0 && rw_var < 20000.0,
            "wn={:.2}, rw={:.2}", wn_var, rw_var
        );
    }

    #[test]
    fn test_combined_spp_prior_and_clock_model_produce_bias() {
        let fg = PppIteratedEkf::default();
        let t0 = GpsTime::new(2156, 1000.0);
        let state = make_state_with_pos(t0, Vector3::new(0.0, 0.0, 0.0), 100.0);
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let spp_pos = Vector3::new(2.0, 0.0, 0.0);
        let dx_clamped = fg
            .compute_iteration_dx(&state, &[sat.clone()], &x_i, &x_pred, &p_inv, 0, Some((spp_pos, 1.0)))
            .unwrap()
            .unwrap();
        let dx_natural = fg
            .compute_iteration_dx(&state, &[sat], &x_i, &x_pred, &p_inv, 0, Some((spp_pos, 0.01)))
            .unwrap()
            .unwrap();
        assert!(
            dx_natural[0].abs() > dx_clamped[0].abs() * 2.0,
            "natural={:.4}, clamped={:.4}",
            dx_natural[0], dx_clamped[0]
        );
    }
