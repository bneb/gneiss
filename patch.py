import re

with open('crates/gneiss-rtk/src/engine/ppp_fg.rs', 'r') as f:
    content = f.read()

test_code = """
    #[test]
    fn test_solve_outlier_rejection_limit() {
        let fg = PppFactorGraph::default();
        let mut state = dummy_rtk_state();
        let mut sats = Vec::new();

        let obs1 = SatObs { sat: SatelliteId { constellation: Constellation::Gps, prn: 1 }, observations: vec![] };
        let obs2 = SatObs { sat: SatelliteId { constellation: Constellation::Gps, prn: 2 }, observations: vec![] };
        let obs3 = SatObs { sat: SatelliteId { constellation: Constellation::Gps, prn: 3 }, observations: vec![] };
        
        let obs1_ref = Box::leak(Box::new(obs1));
        let obs2_ref = Box::leak(Box::new(obs2));
        let obs3_ref = Box::leak(Box::new(obs3));

        let mut add_sat = |obs: &'static SatObs, bad: f64| {
            state.ambiguity_keys.push((obs.sat, 0));
            state.ambiguities.push(0.0);
            sats.push(ProcessedSat {
                sat_obs: obs,
                dt_sat_m: 0.0, p1: 0.0, p2: None, cp1: Some(bad), cp2: None,
                is_iono_free: false, osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
                los: Vector3::zeros(), dist: 0.0, el: 15.01_f64.to_radians(), snr: 45.0, doppler: 0.0,
                lam1: 0.19, lam2: 0.24, tropo_dry: 0.0, map_wet: 0.0, iono_delay: 5.0,
                f1: 1.0, f2: 1.0, sat_pos_rot: Vector3::zeros(), sat_vel: Vector3::zeros(),
                sat_clock_drift: 0.0, rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0
            });
        };

        add_sat(obs1_ref, 10000.0);
        add_sat(obs2_ref, 10000.0);
        add_sat(obs3_ref, 10000.0);

        state.covariance = DMatrix::identity(CORE_STATE_SIZE + 3, CORE_STATE_SIZE + 3);

        assert_eq!(fg.solve(&mut state, &sats), Err(EngineError::MaxIterationsReached));
    }
"""

with open('crates/gneiss-rtk/src/engine/ppp_fg.rs', 'w') as f:
    f.write(content.replace('fn test_solve_empty_sats() {\n        let fg = PppFactorGraph::default();', test_code + '\n    #[test]\n    fn test_solve_empty_sats() {\n        let fg = PppFactorGraph::default();'))

