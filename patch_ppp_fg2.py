import re

with open('crates/gneiss-rtk/src/engine/ppp_fg.rs', 'r') as f:
    content = f.read()

tests = """
    #[test]
    fn test_find_ar_candidates() {
        let fg = PppFactorGraph::default();
        let mut state = dummy_rtk_state();
        let sat_id1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.ambiguity_keys.push((sat_id1, 1));
        state.ambiguity_keys.push((sat_id1, 2));
        state.ambiguities.push(0.0);
        state.ambiguities.push(0.0);

        let obs = SatObs { sat: sat_id1, observations: vec![] };
        let mut sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0, p1: 0.0, p2: None, cp1: Some(0.0), cp2: Some(0.0),
            is_iono_free: false, osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0, el: 15.01_f64.to_radians(), snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24, tropo_dry: 0.0, map_wet: 0.0, iono_delay: 5.0,
            f1: 1.0, f2: 1.0, sat_pos_rot: Vector3::zeros(), sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0, rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0
        };
        let cands = fg.find_ar_candidates(&state, &[sat.clone()]);
        assert_eq!(cands.len(), 1);

        // Test `!s.is_iono_free` mutant
        sat.is_iono_free = true;
        assert_eq!(fg.find_ar_candidates(&state, &[sat.clone()]).len(), 0);
        sat.is_iono_free = false;

        // Test `s.cp2.is_some()` mutant
        sat.cp2 = None;
        assert_eq!(fg.find_ar_candidates(&state, &[sat.clone()]).len(), 0);
        sat.cp2 = Some(0.0);

        // Test constellation
        let sat_id_glo = SatelliteId { constellation: Constellation::Glonass, prn: 1 };
        let obs_glo = SatObs { sat: sat_id_glo, observations: vec![] };
        sat.sat_obs = &obs_glo;
        assert_eq!(fg.find_ar_candidates(&state, &[sat.clone()]).len(), 0);
    }

    #[test]
    fn test_resolve_cascade_ar_bounds() {
        let fg = PppFactorGraph::default();
        let mut state = dummy_rtk_state();
        let mut sats = Vec::new();

        for i in 1..=4 {
            let sat_id = SatelliteId { constellation: Constellation::Gps, prn: i };
            state.ambiguity_keys.push((sat_id, 1));
            state.ambiguity_keys.push((sat_id, 2));
            state.ambiguities.push(0.0);
            state.ambiguities.push(0.0);
            let obs = SatObs { sat: sat_id, observations: vec![] };
            // Need to box or somehow make `sat_obs` live long enough. 
            // In Rust tests we can't easily reference a local loop variable inside a Vec.
            // We can just create them and leak or collect.
        }
    }
"""

content = content.replace("    fn test_find_worst_outlier() {", tests + "\n    fn test_find_worst_outlier() {")

with open('crates/gneiss-rtk/src/engine/ppp_fg.rs', 'w') as f:
    f.write(content)
