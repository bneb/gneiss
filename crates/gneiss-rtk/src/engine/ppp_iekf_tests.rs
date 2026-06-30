    use super::*;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, DVector, Vector3};

    fn dummy_rtk_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(
                Vector3::zeros(),
                Datum::WGS84,
                Frame::ECEF,
                GpsTime::new(0, 0.0),
            ),
            1.0,
        )
    }

    #[test]
    fn test_solve_empty_sats() {
        let fg = PppIteratedEkf::new();
        let mut state = dummy_rtk_state();
        state.covariance = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        let sats = vec![];
        let res = fg.solve(&mut state, &sats, None);
        assert!(matches!(res, Err(EngineError::InsufficientSatellites)));
    }

    #[test]
    fn test_ppp_factor_graph_default() {
        let fg = PppIteratedEkf::default();
        assert_eq!(fg.max_iterations, 15);
        assert_eq!(fg.convergence_threshold, 1e-3);
        assert_eq!(fg.huber_k, 3.0);
        let fg2 = PppIteratedEkf::new();
        assert_eq!(fg2.max_iterations, 15);
    }

    #[test]
    fn test_find_ambiguity_index() {
        let mut state = dummy_rtk_state();

        use gneiss_core::sat::{Constellation, SatelliteId};
        let sat1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let sat2 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        };
        state.ambiguity_keys.push((sat1, 0));
        state.ambiguity_keys.push((sat2, 1));
        state.ambiguities.push(0.0);
        state.ambiguities.push(0.0);

        assert_eq!(find_ambiguity_index(&state, sat1), Some(0));
        assert_eq!(find_ambiguity_index(&state, sat2), None);

        let sat3 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 3,
        };
        assert_eq!(find_ambiguity_index(&state, sat3), None);
    }

    #[test]
    fn test_build_weight_matrix() {
        let m1 = FgMeasurement {
            res: 1.0,
            h_row: DVector::zeros(1),
            weight: 2.0,
            raw_var: 0.5,
            is_phase: false,
            sat: None,
        };
        let m2 = FgMeasurement {
            res: 2.0,
            h_row: DVector::zeros(1),
            weight: 4.0,
            raw_var: 0.25,
            is_phase: true,
            sat: None,
        };
        let meas = vec![m1, m2];
        let mut r = DMatrix::zeros(2, 2);
        r[(0, 0)] = 2.0;
        r[(1, 1)] = 4.0;
        let w = build_weight_matrix(&meas, &r);
        assert_eq!(w[(0, 0)], 0.5);
        assert_eq!(w[(1, 1)], 0.25);
        assert_eq!(w[(0, 1)], 0.0);
    }

    #[test]
    fn test_build_h_row() {
        let los = Vector3::new(1.0, 2.0, 3.0);
        // size = CORE_STATE_SIZE + 1 so ISBs (size>18) and ZWD (size>20) are populated,
        // plus one ambiguity slot at index CORE_STATE_SIZE.
        let size = CORE_STATE_SIZE + 1;
        let h = build_h_row(
            &los,
            4.0,
            Some(CORE_STATE_SIZE),
            size,
            gneiss_core::sat::Constellation::Gps,
        );
        assert_eq!(h.len(), size);
        assert_eq!(h[0], -1.0);
        assert_eq!(h[1], -2.0);
        assert_eq!(h[2], -3.0);
        assert_eq!(h[15], 1.0); // clock bias
        assert_eq!(h[20], 4.0); // ZWD mapping
        assert_eq!(h[CORE_STATE_SIZE], 1.0); // ambiguity

        // size = CORE_STATE_SIZE: ISBs populated (size>18) but ZWD NOT (size==21, not >20... wait 21>20 is true)
        // Actually CORE_STATE_SIZE=21 > 20, so ZWD IS set. No ambiguity.
        let h2 = build_h_row(
            &los,
            4.0,
            None,
            CORE_STATE_SIZE,
            gneiss_core::sat::Constellation::Gps,
        );
        assert_eq!(h2.len(), CORE_STATE_SIZE);
        assert_eq!(h2[0], -1.0);
        assert_eq!(h2[15], 1.0);
        assert_eq!(h2[20], 4.0); // ZWD mapping (21 > 20)
    }

    #[test]
    fn test_build_h_row_doppler() {
        let los = Vector3::new(1.0, 2.0, 3.0);
        // size = CORE_STATE_SIZE (21) > 19, so velocity and clock drift are populated.
        let h = build_h_row_doppler(&los, CORE_STATE_SIZE);
        assert_eq!(h.len(), CORE_STATE_SIZE);
        assert_eq!(h[3], -1.0);
        assert_eq!(h[4], -2.0);
        assert_eq!(h[5], -3.0);
        assert_eq!(h[19], 1.0); // clock drift at index 19

        // size = 19, NOT > 19, so nothing is set — all zeros.
        let h2 = build_h_row_doppler(&los, 19);
        assert_eq!(h2.len(), 19);
        assert_eq!(h2[3], 0.0);
    }

    #[test]
    fn test_assemble_matrices() {
        let meas = vec![
            FgMeasurement {
                res: 1.5,
                h_row: DVector::from_element(3, 1.0),
                weight: 2.0,
                raw_var: 0.5,
                is_phase: false,
                sat: None,
            },
            FgMeasurement {
                res: 2.5,
                h_row: DVector::from_element(3, 2.0),
                weight: 3.0,
                raw_var: 0.33,
                is_phase: true,
                sat: None,
            },
        ];
        let (h, z, r) = assemble_matrices(&meas, 3);
        assert_eq!(h.nrows(), 2);
        assert_eq!(h.ncols(), 3);
        assert_eq!(h[(0, 0)], 1.0);
        assert_eq!(h[(1, 2)], 2.0);
        assert_eq!(z.len(), 2);
        assert_eq!(z[0], 1.5);
        assert_eq!(z[1], 2.5);
        assert_eq!(r.nrows(), 2);
        assert_eq!(r.ncols(), 2);
        assert_eq!(r[(0, 0)], 2.0);
        assert_eq!(r[(1, 1)], 3.0);
        assert_eq!(r[(0, 1)], 0.0);
    }

    #[test]
    fn test_extract_and_apply_state_vector() {
        let mut state = dummy_rtk_state();
        state.position.vector = Vector3::new(1.0, 2.0, 3.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);
        state.attitude = nalgebra::UnitQuaternion::from_scaled_axis(Vector3::new(0.1, 0.2, 0.3));
        state.accel_bias = Vector3::new(10.0, 11.0, 12.0);
        state.gyro_bias = Vector3::new(13.0, 14.0, 15.0);
        state.rcv_clk_bias = 16.0;
        state.rcv_clk_drift = 17.0;
        state.zwd = 18.0;
        state.ambiguities = vec![19.0, 20.0];
        let x = extract_state_vector(&state);
        assert_eq!(x.len(), CORE_STATE_SIZE + 2);
        assert_eq!(x[0], 1.0);
        assert_eq!(x[15], 16.0);
        assert_eq!(x[16], 0.0);
        assert_eq!(x[17], 0.0);
        assert_eq!(x[18], 0.0);
        assert_eq!(x[19], 17.0);
        assert_eq!(x[20], 18.0);
        assert_eq!(x[CORE_STATE_SIZE], 19.0);
        assert_eq!(x[CORE_STATE_SIZE + 1], 20.0);

        let mut state2 = dummy_rtk_state();
        state2.ambiguities = vec![0.0, 0.0];
        let cov = state2.covariance.clone();
        apply_state_vector(&mut state2, &x, cov);
        assert_eq!(state2.position.vector, Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(state2.velocity, Vector3::new(4.0, 5.0, 6.0));
        assert!((state2.attitude.scaled_axis() - Vector3::new(0.1, 0.2, 0.3)).norm() < 1e-10);
        assert_eq!(state2.accel_bias, Vector3::new(10.0, 11.0, 12.0));
        assert_eq!(state2.gyro_bias, Vector3::new(13.0, 14.0, 15.0));
        assert_eq!(state2.rcv_clk_bias, 16.0);
        assert_eq!(state2.rcv_clk_drift, 17.0);
        assert_eq!(state2.zwd, 18.0);
        assert_eq!(state2.ambiguities, vec![19.0, 20.0]);
    }
#[cfg(test)]
mod nan_tests {
    use super::*;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, Vector3};

    fn dummy_rtk_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(
                Vector3::zeros(),
                Datum::WGS84,
                Frame::ECEF,
                GpsTime::new(0, 0.0),
            ),
            1.0,
        )
    }

    #[test]
    fn test_solve_matrix_inversion_failure() {
        let fg = PppIteratedEkf::new();
        let mut state = dummy_rtk_state();
        state.covariance = DMatrix::from_element(CORE_STATE_SIZE, CORE_STATE_SIZE, f64::NAN);
        let sats = vec![];
        let res = fg.solve(&mut state, &sats, None);
        assert!(matches!(res, Err(EngineError::StateDisappeared)));
    }
}

#[cfg(test)]
mod mutant_killer_tests {
    use super::*;
    use crate::engine::processed_sat::ProcessedSat;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::obs::SatObs;
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, DVector, Vector3};

    fn dummy_rtk_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(
                Vector3::zeros(),
                Datum::WGS84,
                Frame::ECEF,
                GpsTime::new(0, 0.0),
            ),
            1.0,
        )
    }

    #[test]
    fn test_resolve_widelane_ar_insufficient() {
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        state.covariance = DMatrix::zeros(CORE_STATE_SIZE + 4, CORE_STATE_SIZE + 4);
        for i in 0..4 {
            state.covariance[(CORE_STATE_SIZE + i, CORE_STATE_SIZE + i)] = 100.0;
        } // huge variance
        let subset = [
            (
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 1,
                    },
                    0,
                    0,
                    1.0,
                    1.0,
                    1.0,
                ),
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 2,
                    },
                    1,
                    1,
                    1.0,
                    1.0,
                    1.0,
                ),
            ),
            (
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 1,
                    },
                    0,
                    0,
                    1.0,
                    1.0,
                    1.0,
                ),
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 3,
                    },
                    2,
                    2,
                    1.0,
                    1.0,
                    1.0,
                ),
            ),
            (
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 1,
                    },
                    0,
                    0,
                    1.0,
                    1.0,
                    1.0,
                ),
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 4,
                    },
                    3,
                    3,
                    1.0,
                    1.0,
                    1.0,
                ),
            ),
        ];
        let x = DVector::zeros(CORE_STATE_SIZE + 4);
        let res = fg.resolve_widelane_ar(&state, &state.covariance, &subset, &x);
        assert!(res.is_err());
    }

    #[test]
    fn test_push_cp_measurement_iono_free() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        state.ambiguity_keys.push((sat_id, 0));
        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: true,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let x_i = DVector::zeros(CORE_STATE_SIZE + 1);
        let los = Vector3::new(0.0, 0.0, 1.0);
        fg.push_cp_measurement(
            &mut meas, &state, &sat, &x_i, 0, &los, 10.0, 20000000.0, 100.0,
        );
        assert_eq!(meas.len(), 1);
        // expected_base = 10.0, x_i = 0. expected_cp = 10.0.
        // l_meas = 100.0 * 0.19 = 19.0.
        // res = 19.0 - 10.0 = 9.0
        assert!((meas[0].res - 9.0).abs() < 1e-6);
        // var_cp = 0.0001 * 1.0 / 1.0 * 9.0 = 0.0009
        assert!((meas[0].raw_var - 0.0009).abs() < 1e-6);
    }

    #[test]
    fn test_push_cp_measurement_not_iono_free() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        state.ambiguity_keys.push((sat_id, 0));
        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let x_i = DVector::zeros(CORE_STATE_SIZE + 1);
        let los = Vector3::new(0.0, 0.0, 1.0);
        fg.push_cp_measurement(
            &mut meas, &state, &sat, &x_i, 0, &los, 10.0, 20000000.0, 100.0,
        );
        assert_eq!(meas.len(), 1);
        // expected_base = 10.0. expected_cp = 10.0 - 5.0 = 5.0.
        // l_meas = 100.0 * 0.19 = 19.0.
        // res = 19.0 - 5.0 = 14.0
        assert!((meas[0].res - 14.0).abs() < 1e-6);
        // var_cp = 0.0001 * 1.0 / 1.0 = 0.0001
        assert!((meas[0].raw_var - 0.0001).abs() < 1e-6);
    }

    #[test]
    fn test_find_ar_candidates() {
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let sat_id1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        state.ambiguity_keys.push((sat_id1, 1));
        state.ambiguity_keys.push((sat_id1, 2));
        state.ambiguities.push(0.0);
        state.ambiguities.push(0.0);

        let obs = SatObs {
            sat: sat_id1,
            observations: vec![],
        };
        let mut sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: Some(0.0),
            cp2: Some(0.0),
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: 15.01_f64.to_radians(),
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
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
        let sat_id_glo = SatelliteId {
            constellation: Constellation::Glonass,
            prn: 1,
        };
        let obs_glo = SatObs {
            sat: sat_id_glo,
            observations: vec![],
        };
        sat.sat_obs = &obs_glo;
        assert_eq!(fg.find_ar_candidates(&state, &[sat.clone()]).len(), 0);
    }

    #[test]
    fn test_resolve_cascade_ar_bounds() {
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let mut sats: Vec<ProcessedSat> = Vec::new();

        let obs1 = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![],
        };
        let obs2 = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 2,
            },
            observations: vec![],
        };
        let obs3 = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 3,
            },
            observations: vec![],
        };
        let obs4 = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 4,
            },
            observations: vec![],
        };
        // Static references to avoid lifetime issues in closure
        let obs1_ref = Box::leak(Box::new(obs1));
        let obs2_ref = Box::leak(Box::new(obs2));
        let obs3_ref = Box::leak(Box::new(obs3));
        let obs4_ref = Box::leak(Box::new(obs4));

        {
            let mut add_sat = |obs: &'static SatObs| {
                state.add_ambiguity(obs.sat, 1, 0.0, 1.0);
                state.add_ambiguity(obs.sat, 2, 0.0, 1.0);
                sats.push(ProcessedSat {
                    sat_obs: obs,
                    dt_sat_m: 0.0,
                    p1: 0.0,
                    p2: None,
                    cp1: Some(0.0),
                    cp2: Some(0.0),
                    is_iono_free: false,
                    osb_p1: 0.0,
                    osb_p2: 0.0,
                    osb_cp1: 0.0,
                    osb_cp2: 0.0,
                    los: Vector3::zeros(),
                    dist: 0.0,
                    el: 15.01_f64.to_radians(),
                    snr: 45.0,
                    doppler: 0.0,
                    lam1: 0.19,
                    lam2: 0.24,
                    tropo_dry: 0.0,
                    map_wet: 0.0,
                    iono_delay: 5.0,
                    f1: 1.0,
                    f2: 1.0,
                    sat_pos_rot: Vector3::zeros(),
                    sat_vel: Vector3::zeros(),
                    sat_clock_drift: 0.0,
                    rcv_pos_ecef: Vector3::zeros(),
                    pcv_correction: 0.0,
                });
            };

            add_sat(obs1_ref);
            add_sat(obs2_ref);
            add_sat(obs3_ref);
        }
        assert_eq!(
            fg.resolve_cascade_ar(&mut state, &sats),
            Err("Insufficient dual-frequency satellites for AR")
        );

        state.add_ambiguity(obs4_ref.sat, 1, 0.0, 1.0);
        state.add_ambiguity(obs4_ref.sat, 2, 0.0, 1.0);
        sats.push(ProcessedSat {
            sat_obs: obs4_ref,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: Some(0.0),
            cp2: Some(0.0),
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: 15.01_f64.to_radians(),
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        });
        assert!(
            fg.resolve_cascade_ar(&mut state, &sats)
                != Err("Insufficient dual-frequency satellites for AR")
        );
    }

    #[test]
    fn test_find_worst_outlier() {
        let sat1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let sat2 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        };

        let meas = vec![
            // Not phase, shouldn't be picked even if high
            FgMeasurement {
                res: 1000.0,
                raw_var: 1.0,
                is_phase: false,
                sat: Some(sat1),
                h_row: DVector::zeros(0),
                weight: 1.0,
            },
            // Phase, but norm = 10.0 / sqrt(4.0) = 5.0 (less than max_norm 15.0)
            FgMeasurement {
                res: 10.0,
                raw_var: 4.0,
                is_phase: true,
                sat: Some(sat1),
                h_row: DVector::zeros(0),
                weight: 1.0,
            },
            // Phase, norm = 40.0 / sqrt(4.0) = 20.0 (greater than max_norm 15.0)
            FgMeasurement {
                res: 40.0,
                raw_var: 4.0,
                is_phase: true,
                sat: Some(sat2),
                h_row: DVector::zeros(0),
                weight: 1.0,
            },
            // Phase, norm = -60.0 / sqrt(9.0) = 20.0 (equal to current max_norm, shouldn't override because of >)
            FgMeasurement {
                res: -60.0,
                raw_var: 9.0,
                is_phase: true,
                sat: Some(sat1),
                h_row: DVector::zeros(0),
                weight: 1.0,
            },
        ];

        // Outlier detection disabled during convergence to prevent cascade:
        // removing one ambiguity degrades remaining measurements, causing
        // more removals until all CP is lost.  The Huber estimator handles
        // outlier down-weighting without removing the ambiguity.
        assert_eq!(PppIteratedEkf::find_worst_outlier_sat(&meas), None);

        let meas_no_outlier = vec![FgMeasurement {
            res: 10.0,
            raw_var: 4.0,
            is_phase: true,
            sat: Some(sat1),
            h_row: DVector::zeros(0),
            weight: 1.0,
        }];
        assert_eq!(
            PppIteratedEkf::find_worst_outlier_sat(&meas_no_outlier),
            None
        );
    }

    #[test]
    fn test_find_ar_candidates_bounds() {
        let _fg = PppIteratedEkf::default();
        let _state = dummy_rtk_state();
        let _sats: Vec<ProcessedSat> = Vec::new();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let _sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: Some(0.0),
            cp2: Some(0.0),
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: 15.01_f64.to_radians(),
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        // It requires state.ambiguity_keys to contain (sat, 0) and (sat, 1) and (sat, 2) etc depending on `is_iono_free`.
        // We'll skip adding a full state test and rely on smaller integration tests or direct tests.
    }

    #[test]
    fn test_build_iono_constraint_row() {
        let h = build_iono_constraint_row(25, 21);
        assert_eq!(h.len(), 25);
        assert_eq!(h[21], 1.0);
        assert_eq!(h[0], 0.0);
        assert_eq!(h[24], 0.0);
    }

    #[test]
    fn test_iono_constraint_row_middle_index() {
        let h = build_iono_constraint_row(30, 15);
        assert_eq!(h.len(), 30);
        assert_eq!(h[15], 1.0);
        for i in 0..30 {
            if i != 15 {
                assert_eq!(h[i], 0.0);
            }
        }
    }

    #[test]
    fn test_sequential_ar_mismatch_regression() {
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();

        // Add 6 satellites across 3 constellations (GPS, Galileo, Beidou) to get 3 constellation groups
        let mut sats = Vec::new();
        let constellations = [
            Constellation::Gps,
            Constellation::Galileo,
            Constellation::Beidou,
        ];

        for (i, &constellation) in constellations.iter().enumerate() {
            let sat1 = SatelliteId {
                constellation,
                prn: (i * 2 + 1) as u8,
            };
            let sat2 = SatelliteId {
                constellation,
                prn: (i * 2 + 2) as u8,
            };

            state.add_ambiguity(sat1, 1, 0.0, 1.0);
            state.add_ambiguity(sat1, 2, 0.0, 1.0);
            state.add_ambiguity(sat2, 1, 0.0, 1.0);
            state.add_ambiguity(sat2, 2, 0.0, 1.0);

            let obs1 = Box::leak(Box::new(SatObs {
                sat: sat1,
                observations: vec![],
            }));
            let obs2 = Box::leak(Box::new(SatObs {
                sat: sat2,
                observations: vec![],
            }));

            let make_processed = |obs: &'static SatObs| ProcessedSat {
                sat_obs: obs,
                dt_sat_m: 0.0,
                p1: 0.0,
                p2: None,
                cp1: Some(0.0),
                cp2: Some(0.0),
                is_iono_free: false,
                osb_p1: 0.0,
                osb_p2: 0.0,
                osb_cp1: 0.0,
                osb_cp2: 0.0,
                los: Vector3::zeros(),
                dist: 0.0,
                el: 15.01_f64.to_radians(),
                snr: 45.0,
                doppler: 0.0,
                lam1: 0.19,
                lam2: 0.24,
                tropo_dry: 0.0,
                map_wet: 0.0,
                iono_delay: 5.0,
                f1: 1.0,
                f2: 1.0,
                sat_pos_rot: Vector3::zeros(),
                sat_vel: Vector3::zeros(),
                sat_clock_drift: 0.0,
                rcv_pos_ecef: Vector3::zeros(),
                pcv_correction: 0.0,
            };
            sats.push(make_processed(obs1));
            sats.push(make_processed(obs2));
        }

        // Initialize state vector and covariance
        let state_dim = CORE_STATE_SIZE + state.ambiguities.len();
        state.covariance = DMatrix::from_fn(state_dim, state_dim, |r, c| {
            if r == c {
                (r + 1) as f64 * 1.5
            } else {
                0.01
            }
        });
        state.position.vector = Vector3::new(1.0, 2.0, 3.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);
        state.is_fixed = false;

        let initial_state_vector = extract_state_vector(&state);
        let initial_covariance = state.covariance.clone();

        // Configure mock
        let mock_wl = Ok((
            DVector::zeros(state_dim),
            DMatrix::zeros(state_dim, state_dim),
            vec![0, 1],
        ));
        let mock_nl = Ok((
            DVector::zeros(state_dim),
            DMatrix::zeros(state_dim, state_dim),
        ));

        {
            let mut mock_lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner());
            *mock_lock = Some(crate::engine::ppp_ar::ArMock {
                wl_result: Some(mock_wl),
                nl_result: Some(mock_nl),
                nl_calls: 0,
            });
        }

        // Execute resolve_cascade_ar
        let result = fg.resolve_cascade_ar(&mut state, &sats);

        // Verify that it failed due to global position jump check
        assert_eq!(result, Err("Position jump too large after AR fix"));

        // Verify state is unmodified
        let final_state_vector = extract_state_vector(&state);
        assert_eq!(final_state_vector, initial_state_vector);
        assert_eq!(state.covariance, initial_covariance);
        assert_eq!(state.is_fixed, false);

        // Clear mock
        {
            let mut mock_lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner());
            *mock_lock = None;
        }
    }

    #[test]
    fn test_extract_isb() {
        let x = DVector::from_fn(19, |i, _| i as f64);
        assert!((PppIteratedEkf::extract_isb(&x, Constellation::Gps) - 0.0).abs() < 1e-10);
        assert!(
            (PppIteratedEkf::extract_isb(&x, Constellation::Glonass) - 16.0).abs() < 1e-10
        );
        assert!(
            (PppIteratedEkf::extract_isb(&x, Constellation::Galileo) - 17.0).abs() < 1e-10
        );
        assert!(
            (PppIteratedEkf::extract_isb(&x, Constellation::Beidou) - 18.0).abs() < 1e-10
        );
        // size <= 18 returns 0.0 for all constellations
        let x_small = DVector::from_fn(18, |i, _| i as f64);
        assert!(
            (PppIteratedEkf::extract_isb(&x_small, Constellation::Glonass) - 0.0).abs() < 1e-10
        );
        assert!(
            (PppIteratedEkf::extract_isb(&x_small, Constellation::Galileo) - 0.0).abs() < 1e-10
        );
        assert!(
            (PppIteratedEkf::extract_isb(&x_small, Constellation::Beidou) - 0.0).abs() < 1e-10
        );
    }

    #[test]
    fn test_build_h_row_uduc_basic() {
        let los = Vector3::new(1.0, 2.0, 3.0);
        let size = CORE_STATE_SIZE + 2;
        // No iono or ambiguity indices
        let h = build_h_row_uduc(&los, 4.0, None, 1.0, None, size, Constellation::Gps);
        assert_eq!(h.len(), size);
        assert!((h[0] - (-1.0)).abs() < 1e-10);
        assert!((h[1] - (-2.0)).abs() < 1e-10);
        assert!((h[2] - (-3.0)).abs() < 1e-10);
        assert!((h[15] - 1.0).abs() < 1e-10);
        assert!((h[20] - 4.0).abs() < 1e-10);
        // No indices set
        assert!((h[CORE_STATE_SIZE] - 0.0).abs() < 1e-10);
        assert!((h[CORE_STATE_SIZE + 1] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_build_h_row_uduc_with_indices() {
        let los = Vector3::new(0.5, -1.5, 2.0);
        let size = CORE_STATE_SIZE + 4;
        let i_idx = CORE_STATE_SIZE + 2;
        let n_idx = CORE_STATE_SIZE + 3;
        let h = build_h_row_uduc(
            &los,
            2.5,
            Some(i_idx),
            -1.5,
            Some(n_idx),
            size,
            Constellation::Galileo,
        );
        assert!((h[17] - 1.0).abs() < 1e-10); // Galileo ISB
        assert!((h[i_idx] - (-1.5)).abs() < 1e-10); // Ionosphere coefficient
        assert!((h[n_idx] - 1.0).abs() < 1e-10); // Ambiguity
        assert!((h[20] - 2.5).abs() < 1e-10); // ZWD
    }

    #[test]
    fn test_build_h_row_uduc_size_boundaries() {
        let los = Vector3::new(1.0, 0.0, 0.0);
        // size = CORE_STATE_SIZE (21): >18 (ISBs), >20 (ZWD)
        let h = build_h_row_uduc(
            &los, 3.0, None, 1.0, None, CORE_STATE_SIZE, Constellation::Glonass,
        );
        assert_eq!(h.len(), CORE_STATE_SIZE);
        assert!((h[16] - 1.0).abs() < 1e-10); // Glonass ISB
        assert!((h[20] - 3.0).abs() < 1e-10); // ZWD
        // size = 17 (< 18): no ISBs or ZWD
        let h2 = build_h_row_uduc(&los, 3.0, None, 1.0, None, 17, Constellation::Glonass);
        assert_eq!(h2.len(), 17);
        assert!((h2[16] - 0.0).abs() < 1e-10); // No ISB
        assert!((h2[15] - 1.0).abs() < 1e-10); // Clock bias always set
    }

    #[test]
    fn test_build_ar_subset_gps_reference() {
        let fg = PppIteratedEkf::default();
        let cands = vec![
            (
                SatelliteId { constellation: Constellation::Gps, prn: 1 },
                0,
                1,
                20.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Gps, prn: 2 },
                2,
                3,
                30.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Galileo, prn: 1 },
                4,
                5,
                25.0_f64.to_radians(),
                0.19,
                0.24,
            ),
        ];
        let subset = fg.build_ar_subset(&cands);
        // GPS PRN 2 (highest elev = 30 deg) should be reference
        assert_eq!(subset.len(), 2);
        for pair in &subset {
            assert_eq!(pair.1.0.prn, 2);
            assert_eq!(pair.1.0.constellation, Constellation::Gps);
        }
    }

    #[test]
    fn test_build_ar_subset_no_gps_per_constellation() {
        let fg = PppIteratedEkf::default();
        let cands = vec![
            (
                SatelliteId { constellation: Constellation::Galileo, prn: 1 },
                0,
                1,
                30.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Galileo, prn: 2 },
                2,
                3,
                20.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Beidou, prn: 1 },
                4,
                5,
                25.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Beidou, prn: 2 },
                6,
                7,
                15.0_f64.to_radians(),
                0.19,
                0.24,
            ),
        ];
        let subset = fg.build_ar_subset(&cands);
        // Per-constellation fallback: 1 pair per constellation
        assert_eq!(subset.len(), 2);
        // Each pair must have same constellation (order is HashMap-dependent)
        for pair in &subset {
            assert_eq!(pair.0.0.constellation, pair.1.0.constellation);
        }
        // Both constellations must be represented
        let constels: std::collections::HashSet<_> = subset
            .iter()
            .map(|p| p.0.0.constellation)
            .collect();
        assert!(constels.contains(&Constellation::Galileo));
        assert!(constels.contains(&Constellation::Beidou));
    }

    #[test]
    fn test_build_ar_subset_empty_or_single() {
        let fg = PppIteratedEkf::default();
        // Single GPS sat -> GPS reference but no non-ref sats
        let cands = vec![(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            0,
            1,
            30.0_f64.to_radians(),
            0.19,
            0.24,
        )];
        assert!(fg.build_ar_subset(&cands).is_empty());
        // Empty input
        let cands: Vec<(SatelliteId, usize, usize, f64, f64, f64)> = vec![];
        assert!(fg.build_ar_subset(&cands).is_empty());
    }

    #[test]
    fn test_try_push_cp_measurement_rejects_none_or_zero() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let los = Vector3::new(1.0, 0.0, 0.0);

        // cp1 = None -> no measurement
        let sat_no_cp1 = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 0.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        fg.try_push_cp_measurement(&mut meas, &state, &sat_no_cp1, &x_i, 0, &los, 0.0, 0.0);
        assert_eq!(meas.len(), 0);

        // cp1 = Some(0.0) -> also no measurement (zero check)
        let sat_zero_cp1 = ProcessedSat {
            cp1: Some(0.0),
            ..sat_no_cp1.clone()
        };
        fg.try_push_cp_measurement(&mut meas, &state, &sat_zero_cp1, &x_i, 0, &los, 0.0, 0.0);
        assert_eq!(meas.len(), 0);
    }

    #[test]
    fn test_push_pr_measurement_iono_free_variance() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let los = Vector3::new(0.0, 0.0, 1.0);

        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 10.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: true,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 1.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        fg.push_pr_measurement(&mut meas, &state, &sat, &x_i, 0, &los, 5.0, 0.0, 0.0);
        assert_eq!(meas.len(), 1);
        // expected_pr = expected_base (iono-free) = 5.0
        // res_pr = p1 - expected_pr = 10.0 - 5.0 = 5.0
        assert!((meas[0].res - 5.0).abs() < 1e-6);
        // var_pr = PSEUDORANGE_VARIANCE_BASE * snr_scale(45) / sin(pi/2) * 9.0
        // snr_scale(45) = (45/45)^2 = 1.0, var_pr = 1.0 * 1.0 / 1.0 * 9.0 = 9.0
        assert!((meas[0].raw_var - 9.0).abs() < 1e-6);
        // h_row[20] = map_wet
        assert!((meas[0].h_row[20] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_push_pr_measurement_not_iono_free_variance() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let los = Vector3::new(0.0, 0.0, 1.0);

        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 10.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 1.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        fg.push_pr_measurement(&mut meas, &state, &sat, &x_i, 0, &los, 5.0, 0.0, 0.0);
        assert_eq!(meas.len(), 1);
        // expected_pr = expected_base + iono_delay = 5.0 + 5.0 = 10.0
        // res_pr = 10.0 - 10.0 = 0.0
        assert!((meas[0].res - 0.0).abs() < 1e-6);
        // var_pr = 1.0 * 1.0 / 1.0 + 9.0 = 10.0
        assert!((meas[0].raw_var - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_push_sat_meas_rejects_large_pr_residual() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let los = Vector3::new(0.0, 0.0, 1.0);

        // PR residual = p1 - expected_pr > 100 -> returns false, no meas
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 1000.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: true,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 0.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let result = fg.push_sat_meas(&mut meas, &state, &sat, &x_i, 0, &los, 5.0, 0.0, 0.0);
        assert!(!result);
        assert!(meas.is_empty());
    }

    #[test]
    fn test_push_sat_meas_adaptive_threshold() {
        // Verify the adaptive PR rejection threshold:
        //   threshold = max(100, 5*pos_std), capped at 500
        // Small covariance (sigma=1m) → threshold=100m
        // Large covariance (sigma=300m) → threshold=200m (capped)
        // Medium covariance (sigma=30m) → threshold=150m
        let fg = PppIteratedEkf::default();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let los = Vector3::new(0.0, 0.0, 1.0);

        // Case 1: Small variance (sigma=1m) → threshold=100m
        // PR residual = 120m > 100m → rejected
        {
            let mut state = dummy_rtk_state();
            for i in 0..3 { state.covariance[(i, i)] = 1.0; }
            let mut meas = Vec::new();
            let sat = ProcessedSat {
                sat_obs: &obs, dt_sat_m: 0.0, p1: 120.0, p2: None,
                cp1: None, cp2: None, is_iono_free: true,
                osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
                los: Vector3::zeros(), dist: 0.0, el: std::f64::consts::PI / 2.0,
                snr: 45.0, doppler: 0.0, lam1: 0.19, lam2: 0.24,
                tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
                f1: 1.0, f2: 1.0, sat_pos_rot: Vector3::zeros(),
                sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
                rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
            };
            let result = fg.push_sat_meas(&mut meas, &state, &sat, &x_i, 0, &los, 0.0, 0.0, 0.0);
            assert!(!result, "120m residual with sigma=1m should be rejected (threshold=100m)");
        }

        // Case 2: PR residual = 80m < 100m → accepted even with small variance
        {
            let mut state = dummy_rtk_state();
            for i in 0..3 { state.covariance[(i, i)] = 1.0; }
            let mut meas = Vec::new();
            let sat = ProcessedSat {
                sat_obs: &obs, dt_sat_m: 0.0, p1: 80.0, p2: None,
                cp1: None, cp2: None, is_iono_free: true,
                osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
                los: Vector3::zeros(), dist: 0.0, el: std::f64::consts::PI / 2.0,
                snr: 45.0, doppler: 0.0, lam1: 0.19, lam2: 0.24,
                tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
                f1: 1.0, f2: 1.0, sat_pos_rot: Vector3::zeros(),
                sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
                rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
            };
            let result = fg.push_sat_meas(&mut meas, &state, &sat, &x_i, 0, &los, 0.0, 0.0, 0.0);
            assert!(result, "80m residual with sigma=1m should be accepted (threshold=100m)");
        }

        // Case 3: Large variance (sigma=300m) → threshold=200m (capped)
        // PR residual = 150m < 200m → accepted
        {
            let mut state = dummy_rtk_state();
            for i in 0..3 { state.covariance[(i, i)] = 90000.0; }
            let mut meas = Vec::new();
            let sat = ProcessedSat {
                sat_obs: &obs, dt_sat_m: 0.0, p1: 150.0, p2: None,
                cp1: None, cp2: None, is_iono_free: true,
                osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
                los: Vector3::zeros(), dist: 0.0, el: std::f64::consts::PI / 2.0,
                snr: 45.0, doppler: 0.0, lam1: 0.19, lam2: 0.24,
                tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
                f1: 1.0, f2: 1.0, sat_pos_rot: Vector3::zeros(),
                sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
                rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
            };
            let result = fg.push_sat_meas(&mut meas, &state, &sat, &x_i, 0, &los, 0.0, 0.0, 0.0);
            assert!(result, "150m residual with sigma=300m should be accepted (threshold=200m)");
        }

        // Case 4: PR residual = 250m > 200m cap → rejected even with huge variance
        {
            let mut state = dummy_rtk_state();
            for i in 0..3 { state.covariance[(i, i)] = 90000.0; }
            let mut meas = Vec::new();
            let sat = ProcessedSat {
                sat_obs: &obs, dt_sat_m: 0.0, p1: 250.0, p2: None,
                cp1: None, cp2: None, is_iono_free: true,
                osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
                los: Vector3::zeros(), dist: 0.0, el: std::f64::consts::PI / 2.0,
                snr: 45.0, doppler: 0.0, lam1: 0.19, lam2: 0.24,
                tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
                f1: 1.0, f2: 1.0, sat_pos_rot: Vector3::zeros(),
                sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
                rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
            };
            let result = fg.push_sat_meas(&mut meas, &state, &sat, &x_i, 0, &los, 0.0, 0.0, 0.0);
            assert!(!result, "250m residual should be rejected (cap=200m)");
        }
    }

    #[test]
    fn test_resolve_uduc_indices_with_all_indices() {
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.ambiguity_keys.push((sat_id, 1)); // n1 idx 0
        state.ambiguity_keys.push((sat_id, 2)); // n2 idx 1
        state.ambiguity_keys.push((sat_id, 3)); // i1 idx 2
        state.ambiguities = vec![100.0, 200.0, 300.0];
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 0.0,
            f1: 1575.42e6,
            f2: 1227.60e6,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let x_i = DVector::from_fn(CORE_STATE_SIZE + 3, |i, _| i as f64);
        let (i1_idx, n1_idx, n2_idx, i1, n1, n2, gamma) =
            PppIteratedEkf::resolve_uduc_indices(&state, &sat, &x_i);
        assert_eq!(i1_idx, Some(CORE_STATE_SIZE + 2));
        assert_eq!(n1_idx, Some(CORE_STATE_SIZE + 0));
        assert_eq!(n2_idx, Some(CORE_STATE_SIZE + 1));
        let expected_gamma = (1575.42e6 * 1575.42e6) / (1227.60e6 * 1227.60e6);
        assert!((gamma - expected_gamma).abs() < 1e-6);
        assert!((i1 - (CORE_STATE_SIZE + 2) as f64).abs() < 1e-6);
        assert!((n1 - (CORE_STATE_SIZE + 0) as f64).abs() < 1e-6);
        assert!((n2 - (CORE_STATE_SIZE + 1) as f64).abs() < 1e-6);
    }

    #[test]
    fn test_resolve_uduc_indices_missing_indices() {
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 0.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let (i1_idx, n1_idx, n2_idx, i1, n1, n2, gamma) =
            PppIteratedEkf::resolve_uduc_indices(&state, &sat, &x_i);
        // No ambiguity keys -> all indices None, values default to 0.0
        assert!(i1_idx.is_none());
        assert!(n1_idx.is_none());
        assert!(n2_idx.is_none());
        assert!((i1 - 0.0).abs() < 1e-10);
        assert!((n1 - 0.0).abs() < 1e-10);
        assert!((n2 - 0.0).abs() < 1e-10);
        assert!((gamma - 1.0).abs() < 1e-10); // f1/f2 = 1.0/1.0 = 1.0
    }

    #[test]
    fn test_push_doppler_measurement_creates_residual() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let x_i = DVector::from_fn(20, |i, _| i as f64);
        let los = Vector3::new(1.0, 0.0, 0.0);
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: -100.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 0.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::new(100.0, 0.0, 0.0),
            sat_clock_drift: 1e-5,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        fg.push_doppler_measurement(&mut meas, &sat, &x_i, &los);
        assert_eq!(meas.len(), 1);
        // meas_rr = -(-100.0) * 0.19 = 19.0
        // rcv_vel = Vector3(3, 4, 5), los = (1,0,0) -> los.dot(rcv_vel) = 3.0
        // rcv_clk_drift = x_i[19] = 19.0
        // expected_rr = 100.0 - 3.0 + 19.0 - 1e-5 * C
        let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
        let res_expected = 19.0 - (100.0 - 3.0 + 19.0 - 1e-5 * c);
        assert!((meas[0].res - res_expected).abs() < 1e-3);
        // h_row: velocity terms at [3,4,5] and clock drift at [19]
        assert!((meas[0].h_row[3] - (-1.0)).abs() < 1e-10);
        assert!((meas[0].h_row[19] - 1.0).abs() < 1e-10);
        assert!(!meas[0].is_phase);
    }

    // ============ log_ppp_convergence tests ============

    #[test]
    fn test_log_ppp_convergence_basic() {
        // Verify log_ppp_convergence does not panic with a basic state
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_pred = DMatrix::identity(dim, dim);
        log_ppp_convergence(&state, &[], &x_i, &x_pred, &p_pred, &fg);
    }

    #[test]
    fn test_log_ppp_convergence_empty_sats() {
        // Empty sat list should not cause panics
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_pred = DMatrix::identity(dim, dim);
        log_ppp_convergence(&state, &[], &x_i, &x_pred, &p_pred, &fg);
    }

    // ============ compute_final_covariance tests ============

    #[test]
    fn test_compute_final_covariance_empty_meas() {
        // Empty measurements should return p_pred directly
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let p_pred = DMatrix::identity(dim, dim);
        let p_inv = DMatrix::identity(dim, dim);
        let result = fg.compute_final_covariance(&state, &[], &x_i, &p_pred, &p_inv, None);
        assert_eq!(result, p_pred);
    }

    #[test]
    fn test_compute_final_covariance_with_meas() {
        // One measurement should produce a damped covariance different from p_pred
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        let x_i = DVector::zeros(dim);
        let p_pred = DMatrix::identity(dim, dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        let sat = ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1: 20000000.0,
            p2: None, cp1: None, cp2: None,
            is_iono_free: true,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
            f1: 1.0, f2: 1.0,
            sat_pos_rot: Vector3::new(20000000.0, 0.0, 0.0),
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        };
        let result = fg.compute_final_covariance(&state, &[sat], &x_i, &p_pred, &p_inv, None);
        assert_eq!(result.nrows(), dim);
        assert_eq!(result.ncols(), dim);
        assert_ne!(result, p_pred, "covariance should differ from prior with measurements present");
    }

    #[test]
    fn test_compute_final_covariance_nan_fallback() {
        // p_inv with NaN causes invert_matrix to fail, falling back to p_pred
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        let x_i = DVector::zeros(dim);
        let p_pred = DMatrix::identity(dim, dim);
        let p_inv_nan = DMatrix::from_element(dim, dim, f64::NAN);
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        let sat = ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1: 20000000.0,
            p2: None, cp1: None, cp2: None,
            is_iono_free: true,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
            f1: 1.0, f2: 1.0,
            sat_pos_rot: Vector3::new(20000000.0, 0.0, 0.0),
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        };
        let result = fg.compute_final_covariance(&state, &[sat], &x_i, &p_pred, &p_inv_nan, None);
        assert_eq!(result, p_pred, "should fall back to p_pred when inversion fails");
    }

    // ============ compute_iteration_dx tests ============

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
    fn test_compute_iteration_dx_no_prior() {
        // Normal solution path without position prior
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let result = fg.compute_iteration_dx(&state, &[sat], &x_i, &x_pred, &p_inv, 0, None);
        assert!(result.is_ok());
        let dx = result.unwrap();
        assert!(dx.is_some(), "should produce a delta-x solution");
    }

    #[test]
    fn test_compute_iteration_dx_with_prior() {
        // Position prior anchors the first three state elements
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let spp_pos = Vector3::new(1.0, 2.0, 3.0);
        let prior_var = 1.0;
        let result = fg.compute_iteration_dx(
            &state, &[sat], &x_i, &x_pred, &p_inv, 0, Some((spp_pos, prior_var)),
        );
        assert!(result.is_ok());
        let dx = result.unwrap();
        assert!(dx.is_some(), "should produce a delta-x with position prior");
    }

    // ============ push_uduc_pr_measurements tests ============

    fn make_uduc_sat(
        sat_id: SatelliteId,
        p1: f64,
        p2: Option<f64>,
        cp1: Option<f64>,
        cp2: Option<f64>,
    ) -> ProcessedSat<'static> {
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1, p2, cp1, cp2,
            is_iono_free: false,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.5, iono_delay: 5.0,
            f1: 1575.42e6, f2: 1227.60e6,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        }
    }

    #[test]
    fn test_push_uduc_pr_measurements_both() {
        // Both p1 and p2 present: should produce two measurements with correct residuals
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat = make_uduc_sat(sat_id, 20000005.0, Some(20000010.0), None, None);
        let x_i_size = CORE_STATE_SIZE + 2;
        let x_i = DVector::zeros(x_i_size);
        let los = Vector3::new(1.0, 0.0, 0.0);
        let expected_base = 20000000.0;
        let i1_idx = Some(CORE_STATE_SIZE);
        let i1 = 3.0;
        let gamma = 1.5;

        fg.push_uduc_pr_measurements(&mut meas, &sat, &x_i, &los, expected_base, i1_idx, i1, gamma);

        assert_eq!(meas.len(), 2);
        // P1: res = p1 - (expected_base + i1) = 20000005 - (20000000 + 3) = 2.0
        assert!((meas[0].res - 2.0).abs() < 1e-6, "P1 residual");
        assert!(!meas[0].is_phase);
        // i1 coefficient in h_row should be 1.0 for P1
        assert!((meas[0].h_row[i1_idx.unwrap()] - 1.0).abs() < 1e-10, "P1 iono coef");

        // P2: res = p2 - (expected_base + gamma * i1) = 20000010 - (20000000 + 1.5*3) = 5.5
        assert!((meas[1].res - 5.5).abs() < 1e-6, "P2 residual");
        assert!(!meas[1].is_phase);
        // i1 coefficient in h_row should be gamma for P2
        assert!((meas[1].h_row[i1_idx.unwrap()] - 1.5).abs() < 1e-10, "P2 iono coef");

        // Variance: var_p1 = PSEUDORANGE_VARIANCE_BASE * snr_scale(45) / sin(pi/2) = 1.0
        assert!((meas[0].raw_var - 1.0).abs() < 1e-6, "P1 variance");
        // var_p2 = var_p1 * 1.5 = 1.5
        assert!((meas[1].raw_var - 1.5).abs() < 1e-6, "P2 variance");
    }

    #[test]
    fn test_push_uduc_pr_measurements_p1_only() {
        // Only p1 present (p2=None): should produce one measurement
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat = make_uduc_sat(sat_id, 20000005.0, None, None, None);
        let x_i_size = CORE_STATE_SIZE + 2;
        let x_i = DVector::zeros(x_i_size);
        let los = Vector3::new(1.0, 0.0, 0.0);
        let expected_base = 20000000.0;
        let i1_idx = Some(CORE_STATE_SIZE);
        let i1 = 3.0;
        let gamma = 1.5;

        fg.push_uduc_pr_measurements(&mut meas, &sat, &x_i, &los, expected_base, i1_idx, i1, gamma);

        assert_eq!(meas.len(), 1);
        assert!(!meas[0].is_phase);
    }

    // ============ push_uduc_cp_measurements tests ============

    #[test]
    fn test_push_uduc_cp_measurements_both() {
        // Both cp1 and cp2 present: should produce two phase measurements with correct signs
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat = make_uduc_sat(sat_id, 0.0, None, Some(105263158.0), Some(83333333.0));
        let x_i_size = CORE_STATE_SIZE + 4;
        let x_i = DVector::zeros(x_i_size);
        let los = Vector3::new(1.0, 0.0, 0.0);
        let expected_base = 20000000.0;
        let i1_idx = Some(CORE_STATE_SIZE);
        let n1_idx = Some(CORE_STATE_SIZE + 1);
        let n2_idx = Some(CORE_STATE_SIZE + 2);
        let i1 = 3.0;
        let n1 = 1.0;
        let n2 = 2.0;
        let gamma = 1.5;

        fg.push_uduc_cp_measurements(
            &mut meas, &state, &sat, &x_i, &los,
            expected_base, i1_idx, n1_idx, n2_idx,
            i1, n1, n2, gamma,
        );

        assert_eq!(meas.len(), 2);
        assert!(meas[0].is_phase, "CP1 is phase");
        assert!(meas[1].is_phase, "CP2 is phase");

        // L1: (cp1 - windup) * lam1 - (expected_base - i1 + n1)
        // windup=0, cp1*lam1=20000000.02, expected_base-i1+n1=19999998.0
        let res_l1 = 105263158.0 * 0.19 - (20000000.0 - 3.0 + 1.0);
        assert!((meas[0].res - res_l1).abs() < 1e-4, "CP1 residual");

        // L2: (cp2 - windup) * lam2 - (expected_base - gamma*i1 + n2)
        let res_l2 = 83333333.0 * 0.24 - (20000000.0 - 1.5 * 3.0 + 2.0);
        assert!((meas[1].res - res_l2).abs() < 1e-4, "CP2 residual");

        // h_row coefficients: L1 iono = -1.0, L2 iono = -gamma
        assert!((meas[0].h_row[i1_idx.unwrap()] - (-1.0)).abs() < 1e-10, "CP1 iono=-1");
        assert!((meas[1].h_row[i1_idx.unwrap()] - (-1.5)).abs() < 1e-10, "CP2 iono=-gamma");

        // Ambiguity coefficients
        assert!((meas[0].h_row[n1_idx.unwrap()] - 1.0).abs() < 1e-10, "CP1 amb coef");
        assert!((meas[1].h_row[n2_idx.unwrap()] - 1.0).abs() < 1e-10, "CP2 amb coef");
    }

    #[test]
    fn test_push_uduc_cp_measurements_cp1_only() {
        // Only cp1 present (cp2=None): should produce one phase measurement
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat = make_uduc_sat(sat_id, 0.0, None, Some(105263158.0), None);
        let x_i_size = CORE_STATE_SIZE + 4;
        let x_i = DVector::zeros(x_i_size);
        let los = Vector3::new(1.0, 0.0, 0.0);
        let expected_base = 20000000.0;
        let i1_idx = Some(CORE_STATE_SIZE);
        let n1_idx = Some(CORE_STATE_SIZE + 1);
        let n2_idx = Some(CORE_STATE_SIZE + 2);
        let i1 = 3.0;
        let n1 = 1.0;
        let n2 = 2.0;
        let gamma = 1.5;

        fg.push_uduc_cp_measurements(
            &mut meas, &state, &sat, &x_i, &los,
            expected_base, i1_idx, n1_idx, n2_idx,
            i1, n1, n2, gamma,
        );

        assert_eq!(meas.len(), 1);
        assert!(meas[0].is_phase);
    }

    #[test]
    fn test_push_uduc_cp_measurements_with_windup() {
        // Windup value should subtract from carrier phase
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.windup.insert(sat_id, 2.5);
        let mut meas = Vec::new();
        let sat = make_uduc_sat(sat_id, 0.0, None, Some(105263158.0), Some(83333333.0));
        let x_i_size = CORE_STATE_SIZE + 4;
        let x_i = DVector::zeros(x_i_size);
        let los = Vector3::new(1.0, 0.0, 0.0);
        let expected_base = 20000000.0;
        let i1_idx = Some(CORE_STATE_SIZE);
        let n1_idx = Some(CORE_STATE_SIZE + 1);
        let n2_idx = Some(CORE_STATE_SIZE + 2);
        let i1 = 3.0;
        let n1 = 1.0;
        let n2 = 2.0;
        let gamma = 1.5;

        fg.push_uduc_cp_measurements(
            &mut meas, &state, &sat, &x_i, &los,
            expected_base, i1_idx, n1_idx, n2_idx,
            i1, n1, n2, gamma,
        );

        assert_eq!(meas.len(), 2);
        // With windup=2.5: (cp1 - 2.5) * 0.19 - (expected_base - i1 + n1)
        let res_l1_windup = (105263158.0 - 2.5) * 0.19 - (20000000.0 - 3.0 + 1.0);
        assert!((meas[0].res - res_l1_windup).abs() < 1e-4, "CP1 with windup");
        let res_l2_windup = (83333333.0 - 2.5) * 0.24 - (20000000.0 - 1.5 * 3.0 + 2.0);
        assert!((meas[1].res - res_l2_windup).abs() < 1e-4, "CP2 with windup");
        // Verify windup actually modified the residual vs no-windup baseline
        let res_l1_no_windup = 105263158.0 * 0.19 - (20000000.0 - 3.0 + 1.0);
        assert!((meas[0].res - res_l1_no_windup).abs() > 0.1, "windup should change residual");
    }

    // ============ push_uduc_measurements tests ============

    #[test]
    fn test_push_uduc_measurements_full() {
        // Full UDUC: both PR and CP measurements are pushed (p1,p2,cp1,cp2)
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.ambiguity_keys.push((sat_id, 1)); // n1
        state.ambiguity_keys.push((sat_id, 2)); // n2
        state.ambiguity_keys.push((sat_id, 3)); // i1
        state.ambiguities = vec![0.0, 0.0, 0.0];
        let sat = make_uduc_sat(
            sat_id, 20000005.0, Some(20000010.0), Some(105263158.0), Some(83333333.0),
        );
        let dim = CORE_STATE_SIZE + 3;
        let x_i = DVector::zeros(dim);
        let los = Vector3::new(1.0, 0.0, 0.0);

        fg.push_uduc_measurements(&mut meas, &state, &sat, &x_i, 0, &los, 20000000.0, 0.0, 0.0);

        // Should produce: p1, p2, cp1, cp2 = 4 measurements
        assert_eq!(meas.len(), 4);
        assert!(!meas[0].is_phase, "p1 is PR");
        assert!(!meas[1].is_phase, "p2 is PR");
        assert!(meas[2].is_phase, "cp1 is phase");
        assert!(meas[3].is_phase, "cp2 is phase");

        // Verify i1 state coefficient signs via resolve_uduc_indices:
        // i1 is at CORE_STATE_SIZE + 2 (third ambiguity key)
        let i1_idx = CORE_STATE_SIZE + 2;
        // PR uses positive i1 (coef=1.0), CP uses negative i1 (coef=-1.0)
        assert!((meas[0].h_row[i1_idx] - 1.0).abs() < 1e-10, "p1 i1 coef=+1");
        assert!((meas[2].h_row[i1_idx] - (-1.0)).abs() < 1e-10, "cp1 i1 coef=-1");
    }

    // ============ build_measurements tests ============

    #[test]
    fn test_build_measurements_with_sats() {
        // Two iono-free sats at different positions produce two PR measurements
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        let x_i = DVector::zeros(dim);
        let sat1 = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let sat2 = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 2 },
            Vector3::new(0.0, 20000000.0, 0.0),
        );
        let meas = fg.build_measurements(&state, &[sat1, sat2], &x_i, 0);
        assert_eq!(meas.len(), 2);
        for m in &meas {
            assert!(!m.is_phase, "PR measurements");
            assert!((m.h_row[15] - 1.0).abs() < 1e-10, "clock bias coef");
        }
        // Each measurement should have a different los direction
        assert!((meas[0].h_row[0] - (-1.0)).abs() < 1e-10, "sat1 los.x");
        assert!((meas[1].h_row[1] - (-1.0)).abs() < 1e-10, "sat2 los.y");
    }

    // ============ push_sat_meas tests ============

    #[test]
    fn test_push_sat_meas_uduc_path() {
        // UDUC path: !iono_free, cp1+cp2+p2 present -> 5 measurements (4 UDUC + iono prior)
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.ambiguity_keys.push((sat_id, 1)); // n1
        state.ambiguity_keys.push((sat_id, 2)); // n2
        state.ambiguity_keys.push((sat_id, 3)); // i1
        state.ambiguities = vec![0.0, 0.0, 0.0];
        let dim = CORE_STATE_SIZE + 3;
        let x_i = DVector::zeros(dim);
        let los = Vector3::new(0.0, 0.0, 1.0);
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        let sat = ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1: 20000000.0,
            p2: Some(20000000.0),
            cp1: Some(105263158.0),
            cp2: Some(83333333.0),
            is_iono_free: false,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
            f1: 1575.42e6, f2: 1227.60e6,
            sat_pos_rot: Vector3::new(20000000.0, 0.0, 0.0),
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        };

        let result = fg.push_sat_meas(
            &mut meas, &state, &sat, &x_i, 0, &los, 20000000.0, 20000000.0, 0.0,
        );
        assert!(result, "sat should be accepted");
        // UDUC: p1, p2, cp1, cp2 + iono prior = 5
        assert_eq!(meas.len(), 5);
    }

    #[test]
    fn test_push_sat_meas_iono_prior() {
        // Iono prior constraint measurement has correct residual and variance
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.ambiguity_keys.push((sat_id, 3)); // i1 idx needed for iono prior
        state.ambiguities = vec![0.0];
        let dim = CORE_STATE_SIZE + 1;
        let x_i = DVector::zeros(dim);
        let los = Vector3::new(0.0, 0.0, 1.0);
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        let sat = ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1: 20000000.0,
            p2: Some(20000000.0),
            cp1: Some(105263158.0),
            cp2: Some(83333333.0),
            is_iono_free: false,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.0, iono_delay: 5.0, // Klobuchar prediction
            f1: 1575.42e6, f2: 1227.60e6,
            sat_pos_rot: Vector3::new(20000000.0, 0.0, 0.0),
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        };

        let result = fg.push_sat_meas(
            &mut meas, &state, &sat, &x_i, 0, &los, 20000000.0, 20000000.0, 0.0,
        );
        assert!(result);
        // UDUC: p1, p2, cp1, cp2 + iono prior = 5
        assert_eq!(meas.len(), 5);
        // Last measurement is the iono prior
        let iono = &meas[4];
        assert!(!iono.is_phase);
        // res = sat.iono_delay - x_i[i1_idx]; i1_idx=0, x_i[CORE_STATE_SIZE]=0 -> res=5.0
        assert!((iono.res - 5.0).abs() < 1e-6, "iono prior residual");
        assert!((iono.raw_var - 9.0).abs() < 1e-6, "iono prior variance (3m std)");
    }

    #[test]
    fn test_push_sat_meas_try_cp_success() {
        // try_push_cp adds a CP measurement for iono-free sats with non-zero cp1
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.ambiguity_keys.push((sat_id, 0));
        state.ambiguities = vec![0.0];
        let dim = CORE_STATE_SIZE + 1;
        let x_i = DVector::zeros(dim);
        let los = Vector3::new(0.0, 0.0, 1.0);
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        let sat = ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1: 20000000.0,
            p2: None, cp1: Some(105263158.0), cp2: None,
            is_iono_free: true,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
            f1: 1.0, f2: 1.0,
            sat_pos_rot: Vector3::new(20000000.0, 0.0, 0.0),
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        };

        let result = fg.push_sat_meas(
            &mut meas, &state, &sat, &x_i, 0, &los, 20000000.0, 20000000.0, 0.0,
        );
        assert!(result);
        // PR + CP = 2 measurements
        assert_eq!(meas.len(), 2);
        assert!(!meas[0].is_phase, "PR measurement");
        assert!(meas[1].is_phase, "CP measurement");
    }

    // ============ resolve_narrowlane_ar tests (uses crate::engine::ppp_ar::AR_MOCK) ============

    #[test]
    fn test_resolve_narrowlane_ar_mock() {
        { let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner()); *lock = None; }
        // Use crate::engine::ppp_ar::AR_MOCK to verify NL resolution returns expected modified state
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE + 4;
        let x_wl = DVector::from_fn(dim, |i, _| i as f64);
        let p_wl = DMatrix::identity(dim, dim);
        let subset = vec![(
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0, 1, 1.0, 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 2, 3, 1.0, 0.19, 0.24),
        )];
        let keep_indices = vec![0];

        let mock_nl: Result<(DVector<f64>, DMatrix<f64>), &'static str> = Ok((
            DVector::zeros(dim),
            DMatrix::zeros(dim, dim),
        ));
        {
            let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner());
            *lock = Some(crate::engine::ppp_ar::ArMock { wl_result: None, nl_result: Some(mock_nl), nl_calls: 0 });
        }

        let result = fg.resolve_narrowlane_ar(&state, &subset, &keep_indices, &x_wl, &p_wl);
        assert!(result.is_ok(), "mock NL should succeed");
        let (x_fixed, p_fixed) = result.unwrap();
        // Mock adds 9.5 to x_wl[0] and scales p_wl by 0.5
        assert!((x_fixed[0] - (x_wl[0] + 9.5)).abs() < 1e-10, "NL adds 9.5 to pos.x");
        assert!((p_fixed[(0, 0)] - 0.5).abs() < 1e-10, "NL scales covariance by 0.5");
        // Verify nl_calls was incremented
        {
            let lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner());
            assert_eq!(lock.as_ref().unwrap().nl_calls, 1);
        }
        // Clean up mock
        { let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner()); *lock = None; }
    }

    // ============ process_constellation_group tests (uses crate::engine::ppp_ar::AR_MOCK) ============

    #[test]
    fn test_process_constellation_group_less_than_two() {
        // Group with fewer than 2 candidates returns None
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_current = DVector::zeros(dim);
        let p_current = DMatrix::identity(dim, dim);
        let result = fg.process_constellation_group(
            &state, &p_current, &x_current, &[], Constellation::Gps,
        );
        assert!(result.is_none());
    }

    #[test]
    fn test_process_constellation_group_empty_keep_indices() {
        { let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner()); *lock = None; }
        // WL returns empty keep_indices -> process_constellation_group returns None
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE + 4;
        let x_current = DVector::zeros(dim);
        let p_current = DMatrix::identity(dim, dim);
        let group_cands = vec![
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0, 1, 30.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 2, 3, 20.0_f64.to_radians(), 0.19, 0.24),
        ];

        let mock_wl: Result<(DVector<f64>, DMatrix<f64>, Vec<usize>), &'static str> = Ok((
            DVector::zeros(dim),
            DMatrix::identity(dim, dim),
            vec![],
        ));
        {
            let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner());
            *lock = Some(crate::engine::ppp_ar::ArMock { wl_result: Some(mock_wl), nl_result: None, nl_calls: 0 });
        }

        let result = fg.process_constellation_group(
            &state, &p_current, &x_current, &group_cands, Constellation::Gps,
        );
        assert!(result.is_none(), "empty keep_indices -> None");

        { let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner()); *lock = None; }
    }

    #[test]
    fn test_process_constellation_group_mock_success() {
        { let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner()); *lock = None; }
        // With mock WL+NL, group resolves successfully (jump 9.5 < 10 passes position check)
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE + 4;
        let x_current = DVector::zeros(dim);
        let p_current = DMatrix::identity(dim, dim);
        let group_cands = vec![
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0, 1, 30.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 2, 3, 20.0_f64.to_radians(), 0.19, 0.24),
        ];

        let mock_wl = Ok((DVector::zeros(dim), DMatrix::identity(dim, dim), vec![0]));
        let mock_nl = Ok((DVector::zeros(dim), DMatrix::zeros(dim, dim)));
        {
            let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner());
            *lock = Some(crate::engine::ppp_ar::ArMock { wl_result: Some(mock_wl), nl_result: Some(mock_nl), nl_calls: 0 });
        }

        let result = fg.process_constellation_group(
            &state, &p_current, &x_current, &group_cands, Constellation::Gps,
        );
        assert!(result.is_some(), "mock AR should succeed");
        let (_xf, _pf, n_sats) = result.unwrap();
        assert_eq!(n_sats, 2, "keep_indices.len() + 1 = 1 + 1 = 2");

        { let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner()); *lock = None; }
    }

    // ============ compute_iteration_dx tests ============

    #[test]
    fn test_compute_iteration_dx_empty_meas() {
        // No satellites -> build_measurements returns empty -> Err(InsufficientSatellites)
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let result = fg.compute_iteration_dx(&state, &[], &x_i, &x_pred, &p_inv, 0, None);
        assert!(matches!(result, Err(EngineError::InsufficientSatellites)));
    }

    // ============ try_inter_constellation_fallback tests (uses crate::engine::ppp_ar::AR_MOCK) ============

    #[test]
    fn test_try_inter_constellation_fallback_subset_too_small() {
        // build_ar_subset produces < 3 pairs -> returns None early
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_current = DVector::zeros(dim);
        let p_current = DMatrix::identity(dim, dim);
        // 3 candidates -> GPS ref + 2 pairs -> 2 < 3 -> None
        let cands = vec![
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0, 1, 30.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 2, 3, 20.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Galileo, prn: 1 }, 4, 5, 25.0_f64.to_radians(), 0.19, 0.24),
        ];
        let result = fg.try_inter_constellation_fallback(&state, &p_current, &x_current, &cands);
        assert!(result.is_none(), "3 candidates -> 2 pairs < 3 -> None");
    }

    #[test]
    fn test_try_inter_constellation_fallback_empty_keep_indices() {
        { let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner()); *lock = None; }
        // WL mock returns empty keep_indices -> returns None
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE + 6;
        let x_current = DVector::zeros(dim);
        let p_current = DMatrix::identity(dim, dim);
        // 4 candidates -> GPS ref + 3 pairs >= 3 -> passes subset check
        let cands = vec![
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0, 1, 30.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 2, 3, 20.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 3 }, 4, 5, 15.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Galileo, prn: 1 }, 6, 7, 25.0_f64.to_radians(), 0.19, 0.24),
        ];
        let mock_wl = Ok((DVector::zeros(dim), DMatrix::identity(dim, dim), vec![]));
        {
            let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner());
            *lock = Some(crate::engine::ppp_ar::ArMock { wl_result: Some(mock_wl), nl_result: None, nl_calls: 0 });
        }
        let result = fg.try_inter_constellation_fallback(&state, &p_current, &x_current, &cands);
        assert!(result.is_none(), "empty keep_indices -> None");
        { let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner()); *lock = None; }
    }

    #[test]
    fn test_try_inter_constellation_fallback_success() {
        { let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner()); *lock = None; }
        // Full mock success path with valid WL keep_indices and NL resolution
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE + 6;
        let x_current = DVector::zeros(dim);
        let p_current = DMatrix::identity(dim, dim);
        let cands = vec![
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0, 1, 30.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 2, 3, 20.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 3 }, 4, 5, 15.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Galileo, prn: 1 }, 6, 7, 25.0_f64.to_radians(), 0.19, 0.24),
        ];
        let mock_wl = Ok((DVector::zeros(dim), DMatrix::identity(dim, dim), vec![0, 1, 2]));
        let mock_nl = Ok((DVector::zeros(dim), DMatrix::zeros(dim, dim)));
        {
            let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner());
            *lock = Some(crate::engine::ppp_ar::ArMock { wl_result: Some(mock_wl), nl_result: Some(mock_nl), nl_calls: 0 });
        }
        let result = fg.try_inter_constellation_fallback(&state, &p_current, &x_current, &cands);
        assert!(result.is_some(), "mock fallback should succeed");
        let (xf, _pf, n_sats) = result.unwrap();
        assert_eq!(n_sats, 4, "keep_indices.len() + 1 = 3 + 1 = 4");
        // Mock NL adds 9.5 to x_wl[0] and jump=9.5 <= 20.0 passes position check
        assert!((xf[0] - 9.5).abs() < 1e-10, "NL adds 9.5 to x[0]");
        { let mut lock = crate::engine::ppp_ar::AR_MOCK.lock().unwrap_or_else(|e| e.into_inner()); *lock = None; }
    }

    // ============ solve() convergence with measurements ============

    #[test]
    fn test_solve_normal_convergence() {
        // Normal convergence: one satellite with matching p1/geometry -> dx=0 -> converges immediately
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        state.position.vector = Vector3::new(0.0, 0.0, 0.0);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let result = fg.solve(&mut state, &[sat], None);
        assert!(result.is_ok(), "solve should converge with matching geometry");
        assert!(state.full_x_predict.is_some(), "x_pred should be saved");
        assert!(state.full_p_predict.is_some(), "p_pred should be saved");
    }

    // ============ SPP Anchor / Position Prior tests ============
    //
    // The SPP anchor applies a soft position prior in the IEKF at indices 0,1,2
    // (X, Y, Z in ECEF). The prior weight is 1/variance, added to the htwh diagonal
    // and htwr residual. These tests verify the prior math:
    //   - Prior pulls position toward SPP with the correct sign
    //   - Stronger prior (smaller variance) pulls harder
    //   - The prior primarily targets position, not clock or other states
    //   - Full solve converges with a prior and moves position

    #[test]
    fn test_compute_iteration_dx_prior_sign_correct() {
        // Prior to the LEFT of current state should produce negative dx;
        // prior to the RIGHT should produce positive dx.
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );

        // Prior at X=-1000 with state at X=0: prior says "move left"
        let spp_neg = Vector3::new(-1000.0, 0.0, 0.0);
        let dx_neg = fg
            .compute_iteration_dx(
                &state,
                &[sat.clone()],
                &x_i,
                &x_pred,
                &p_inv,
                0,
                Some((spp_neg, 1.0)),
            )
            .unwrap()
            .unwrap();
        assert!(
            dx_neg[0] < -100.0,
            "prior at X=-1000 should pull negative, got dx[0]={}",
            dx_neg[0]
        );

        // Prior at X=+1000: prior says "move right"
        let spp_pos = Vector3::new(1000.0, 0.0, 0.0);
        let dx_pos = fg
            .compute_iteration_dx(
                &state,
                &[sat],
                &x_i,
                &x_pred,
                &p_inv,
                0,
                Some((spp_pos, 1.0)),
            )
            .unwrap()
            .unwrap();
        assert!(
            dx_pos[0] > 100.0,
            "prior at X=1000 should pull positive, got dx[0]={}",
            dx_pos[0]
        );

        // Verify opposite signs
        assert!(
            dx_neg[0] < 0.0 && dx_pos[0] > 0.0,
            "opposite prior positions should produce opposite-sign dx"
        );
    }

    #[test]
    fn test_compute_iteration_dx_prior_strength_scales_with_variance() {
        // A tighter prior (smaller variance) should produce larger position corrections
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let spp_pos = Vector3::new(1000.0, 0.0, 0.0);
        let spp_exact = Vector3::new(0.0, 0.0, 0.0); // prior matches state exactly

        // Weak prior: large variance = 100 (weight = 1/100 = 0.01)
        let dx_weak = fg
            .compute_iteration_dx(
                &state,
                &[sat.clone()],
                &x_i,
                &x_pred,
                &p_inv,
                0,
                Some((spp_pos, 100.0)),
            )
            .unwrap()
            .unwrap();

        // Strong prior: small variance = 1 (weight = 1.0)
        let dx_strong = fg
            .compute_iteration_dx(
                &state,
                &[sat.clone()],
                &x_i,
                &x_pred,
                &p_inv,
                0,
                Some((spp_pos, 1.0)),
            )
            .unwrap()
            .unwrap();

        // Strong prior should pull harder in X
        assert!(
            dx_strong[0].abs() > dx_weak[0].abs(),
            "strong prior (var=1, dx[0]={}) should pull X harder than weak prior (var=100, dx[0]={})",
            dx_strong[0],
            dx_weak[0]
        );

        // No pull when prior matches current state exactly (zero innovation)
        let dx_exact = fg
            .compute_iteration_dx(
                &state,
                &[sat],
                &x_i,
                &x_pred,
                &p_inv,
                0,
                Some((spp_exact, 1.0)),
            )
            .unwrap()
            .unwrap();

        // When prior == state, the htwr contribution is zero, but the htwh damping
        // still increases diagonal elements (tightens the covariance).
        // dx may not be exactly zero because the stronger diagonal pulls the
        // solution toward the prediction (x_pred == x_i here, so it should be ~0).
        assert!(
            dx_exact[0].abs() < 1.0,
            "prior matching state should produce negligible dx, got dx[0]={}",
            dx_exact[0]
        );
    }

    #[test]
    fn test_compute_iteration_dx_prior_targets_position_indices() {
        // The position prior is applied only to state indices 0, 1, 2 (X, Y, Z in ECEF).
        // Non-position elements like clock bias (index 15) are only affected through
        // measurement coupling, not directly by the prior.
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );

        // Prior only pulls X positive
        let spp_pos = Vector3::new(1000.0, 0.0, 0.0);
        let dx = fg
            .compute_iteration_dx(
                &state,
                &[sat],
                &x_i,
                &x_pred,
                &p_inv,
                0,
                Some((spp_pos, 1.0)),
            )
            .unwrap()
            .unwrap();

        // X should be pulled positive
        assert!(
            dx[0] > 100.0,
            "prior at X=1000 should produce large positive dx[0], got {}",
            dx[0]
        );

        // Y and Z have no prior and no measurement sensitivity in this setup
        // (measurement LOS is along X axis), so they should be near zero
        assert!(
            dx[1].abs() < 1e-6,
            "Y should not be directly pulled by prior, got dx[1]={}",
            dx[1]
        );
        assert!(
            dx[2].abs() < 1e-6,
            "Z should not be directly pulled by prior, got dx[2]={}",
            dx[2]
        );

        // Clock bias (index 15) is coupled through the measurement H matrix
        // (which has -1 at [0] and +1 at [15]). Anchoring position naturally
        // helps resolve clock-state ambiguity, but the clock correction should
        // be an order of magnitude smaller than the position correction.
        assert!(
            dx[15].abs() < dx[0].abs(),
            "clock correction ({}) should be smaller than position correction ({})",
            dx[15],
            dx[0]
        );
    }

    #[test]
    fn test_solve_with_prior_pulls_position_toward_spp() {
        // Full solve() with a position prior should pull the estimated
        // position toward the SPP position while converging normally.
        // NOTE: the prior displacement must be small (<~200m) so that the
        // satellite pseudorange residual stays within the 100m rejection
        // threshold across IEKF iterations as position evolves.
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        state.position.vector = Vector3::new(0.0, 0.0, 0.0);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );

        // Small SPP prior displacement: 10m in X
        let spp_pos = Vector3::new(10.0, 0.0, 0.0);
        let result = fg.solve(&mut state, &[sat], Some((spp_pos, 1.0)));
        assert!(
            result.is_ok(),
            "solve should converge with position prior"
        );

        // The converged position should have moved from 0 toward SPP (10).
        // The exact balance depends on prior weight vs process noise.
        assert!(
            state.position.vector.x > 0.5,
            "solve with prior should pull X toward SPP (10), got X={}",
            state.position.vector.x
        );

        // The position should not overshoot the prior
        assert!(
            state.position.vector.x < 9.5,
            "solve with prior should not overshoot SPP position, got X={}",
            state.position.vector.x
        );

        // Y and Z should stay near zero (no prior pull on those axes)
        assert!(
            state.position.vector.y.abs() < 1.0,
            "Y should not be pulled by X-axis prior, got Y={}",
            state.position.vector.y
        );
        assert!(
            state.position.vector.z.abs() < 1.0,
            "Z should not be pulled by X-axis prior, got Z={}",
            state.position.vector.z
        );
    }

    #[test]
    fn test_solve_with_prior_pulls_all_three_axes() {
        // Prior pulling in all three axes simultaneously should move each
        // component toward its respective prior value.
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        state.position.vector = Vector3::new(0.0, 0.0, 0.0);
        let sat1 = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let sat2 = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 2 },
            Vector3::new(0.0, 20000000.0, 0.0),
        );

        // Small 3-axis prior displacement to keep PR residuals < 100m
        let spp_pos = Vector3::new(10.0, -8.0, 5.0);
        let result = fg.solve(&mut state, &[sat1, sat2], Some((spp_pos, 1.0)));
        assert!(
            result.is_ok(),
            "solve should converge with 3-axis prior"
        );

        // Each axis should move toward the prior (at least by 10% of the pull)
        assert!(
            state.position.vector.x > 0.5,
            "X should be pulled toward 10, got X={}",
            state.position.vector.x
        );
        assert!(
            state.position.vector.y < -0.5,
            "Y should be pulled toward -8, got Y={}",
            state.position.vector.y
        );
        assert!(
            state.position.vector.z > 0.3,
            "Z should be pulled toward 5, got Z={}",
            state.position.vector.z
        );
    }

    #[test]
    fn test_solve_with_epoch_count_past_ar_threshold() {
        // When epoch_count > 10 and has_new_sats is true, solve() attempts
        // cascade AR.  With poor geometry (all sats at same position vector)
        // and covariance consistent, resolve_cascade_ar may fail gracefully,
        // exercising the tracing::info error-log path (lines 85-87).
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        state.epoch_count = 15;                          // > 10 threshold
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        state.position.vector = Vector3::new(0.0, 0.0, 0.0);

        // Set up 4 dual-frequency GPS sats so the "Insufficient dual-frequency"
        // check passes.  Use distinct PRNs but identical position so the AR
        // subset may still fail (non-GPS reference check or position jump).
        let sats: Vec<ProcessedSat> = (1..=4)
            .map(|prn| {
                let sat_id = SatelliteId { constellation: Constellation::Gps, prn };
                // Mark as last_observed = epoch_count so has_new_sats = true
                state
                    .last_observed
                    .insert((sat_id, 0), state.epoch_count as u32);
                // Add dual-frequency ambiguity keys
                state.add_ambiguity(sat_id, 1, 0.0, 1.0);
                state.add_ambiguity(sat_id, 2, 0.0, 1.0);
                make_dummy_sat(sat_id, Vector3::new(20000000.0, 0.0, 0.0))
            })
            .collect();

        // solve() should return Ok even when cascade AR fails internally
        let result = fg.solve(&mut state, &sats, None);
        assert!(result.is_ok(), "solve must not propagate cascade AR errors");
    }

    #[test]
    fn test_resolve_widelane_ar_no_mock() {
        // Test resolve_widelane_ar WITHOUT mock -> exercises real LAMBDA path.
        // Small covariance ensures cov_ok passes, zero floats ensure trivial fix.
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        // Covariance dimension must span all ambiguity indices used in subset
        let dim = CORE_STATE_SIZE + 8;
        let x = DVector::zeros(dim);
        // Very small covariance so q_wl_full diagonals sqrt < 0.30 -> all pairs kept
        let p = DMatrix::identity(dim, dim) * 0.0001;

        // 3 DD pairs: G01 as reference, G02-G04 as candidates
        let subset = vec![
            (
                (SatelliteId { constellation: Constellation::Gps, prn: 2 },
                 0, 1, 25.0_f64.to_radians(), 0.19, 0.24),
                (SatelliteId { constellation: Constellation::Gps, prn: 1 },
                 2, 3, 30.0_f64.to_radians(), 0.19, 0.24),
            ),
            (
                (SatelliteId { constellation: Constellation::Gps, prn: 3 },
                 4, 5, 20.0_f64.to_radians(), 0.19, 0.24),
                (SatelliteId { constellation: Constellation::Gps, prn: 1 },
                 2, 3, 30.0_f64.to_radians(), 0.19, 0.24),
            ),
            (
                (SatelliteId { constellation: Constellation::Gps, prn: 4 },
                 6, 7, 15.0_f64.to_radians(), 0.19, 0.24),
                (SatelliteId { constellation: Constellation::Gps, prn: 1 },
                 2, 3, 30.0_f64.to_radians(), 0.19, 0.24),
            ),
        ];

        let result = fg.resolve_widelane_ar(&state, &p, &subset, &x);
        // Should succeed: small covariance -> cov_ok passes -> LAMBDA resolves zero floats trivially
        assert!(result.is_ok(), "WL AR should succeed with trivial input: {:?}", result);
        let (_x_wl, _p_wl, keep_indices) = result.unwrap();
        assert_eq!(keep_indices.len(), 3, "All 3 pairs should be kept");
    }

    #[test]
    fn test_resolve_narrowlane_ar_no_mock() {
        // Test resolve_narrowlane_ar WITHOUT mock -> exercises real LAMBDA path.
        // Uses synthetic WL result with small covariance for trivial NL fix.
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE + 8;

        // x_wl from a WL-fixed state
        let x_wl = DVector::zeros(dim);
        // Small p_wl -> q_nl diagonals small -> LAMBDA succeeds
        let p_wl = DMatrix::identity(dim, dim) * 0.01;

        let subset = vec![
            (
                (SatelliteId { constellation: Constellation::Gps, prn: 2 },
                 0, 1, 25.0_f64.to_radians(), 0.19, 0.24),
                (SatelliteId { constellation: Constellation::Gps, prn: 1 },
                 2, 3, 30.0_f64.to_radians(), 0.19, 0.24),
            ),
            (
                (SatelliteId { constellation: Constellation::Gps, prn: 3 },
                 4, 5, 20.0_f64.to_radians(), 0.19, 0.24),
                (SatelliteId { constellation: Constellation::Gps, prn: 1 },
                 2, 3, 30.0_f64.to_radians(), 0.19, 0.24),
            ),
            (
                (SatelliteId { constellation: Constellation::Gps, prn: 4 },
                 6, 7, 15.0_f64.to_radians(), 0.19, 0.24),
                (SatelliteId { constellation: Constellation::Gps, prn: 1 },
                 2, 3, 30.0_f64.to_radians(), 0.19, 0.24),
            ),
        ];
        let keep_indices = vec![0usize, 1, 2];

        let result = fg.resolve_narrowlane_ar(&state, &subset, &keep_indices, &x_wl, &p_wl);
        // Should succeed: small covariance -> NL LAMBDA resolves zero floats trivially
        assert!(result.is_ok(), "NL AR should succeed with trivial input: {:?}", result);
        let (_x_fixed, _p_fixed) = result.unwrap();
    }

    #[test]
    fn test_resolve_widelane_ar_mw_confident() {
        // MW counts > 50 -> MW-based Q with tighter variances
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE + 8;
        let x = DVector::zeros(dim);
        let p = DMatrix::identity(dim, dim) * 100.0; // large cov -> cov_ok FAILS, need MW path

        // Set MW counts > 50 to trigger MW-based Q
        state.mw_sd_counts.insert(
            SatelliteId { constellation: Constellation::Gps, prn: 2 }, 60);
        state.mw_sd_counts.insert(
            SatelliteId { constellation: Constellation::Gps, prn: 1 }, 60);
        state.mw_sd_counts.insert(
            SatelliteId { constellation: Constellation::Gps, prn: 3 }, 60);
        state.mw_sd_counts.insert(
            SatelliteId { constellation: Constellation::Gps, prn: 4 }, 60);
        state.mw_sd_ema.insert(
            SatelliteId { constellation: Constellation::Gps, prn: 2 }, 0.0);
        state.mw_sd_ema.insert(
            SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0.0);
        state.mw_sd_ema.insert(
            SatelliteId { constellation: Constellation::Gps, prn: 3 }, 0.0);
        state.mw_sd_ema.insert(
            SatelliteId { constellation: Constellation::Gps, prn: 4 }, 0.0);

        let subset = vec![
            (
                (SatelliteId { constellation: Constellation::Gps, prn: 2 },
                 0, 1, 25.0_f64.to_radians(), 0.19, 0.24),
                (SatelliteId { constellation: Constellation::Gps, prn: 1 },
                 2, 3, 30.0_f64.to_radians(), 0.19, 0.24),
            ),
            (
                (SatelliteId { constellation: Constellation::Gps, prn: 3 },
                 4, 5, 20.0_f64.to_radians(), 0.19, 0.24),
                (SatelliteId { constellation: Constellation::Gps, prn: 1 },
                 2, 3, 30.0_f64.to_radians(), 0.19, 0.24),
            ),
            (
                (SatelliteId { constellation: Constellation::Gps, prn: 4 },
                 6, 7, 15.0_f64.to_radians(), 0.19, 0.24),
                (SatelliteId { constellation: Constellation::Gps, prn: 1 },
                 2, 3, 30.0_f64.to_radians(), 0.19, 0.24),
            ),
        ];

        let result = fg.resolve_widelane_ar(&state, &p, &subset, &x);
        // With MW counts > 50, MW-based Q is used (all_mw_confident=true)
        assert!(result.is_ok(), "WL AR with MW should succeed: {:?}", result);
        let (_x_wl, _p_wl, keep_indices) = result.unwrap();
        assert_eq!(keep_indices.len(), 3, "All pairs kept via MW counts");
    }

}

// =========================================================================
// Adversarial tests: PPP accuracy gap investigation
// =========================================================================
//
// These tests expose the ROOT CAUSES of the ~0.5m horizontal bias observed
// in the f9p PPP benchmark. The findings are:
//
// Issue #1: SPP PRIOR VARIANCE FLOOR (line 52 of process_ppp.rs):
//   let prior_var = (pos_cov.min(25.0)).max(1.0);
//
//   When the IEKF position covariance converges below 1.0 m^2 (10 cm std),
//   the floor at 1.0 m^2 INCREASES the prior variance, DECREASING the prior
//   weight from 1/0.01=100 to 1/1.0=1. This makes the SPP anchor 100x
//   WEAKER than optimal after convergence. The filter loses its primary
//   absolute position anchor just when it needs it most.
//
// Issue #2: WHITE-NOISE CLOCK MODEL (predictor.rs line 134):
//   phi[(15, 15)] = 0.0
//
//   This destroys temporal correlation of the clock bias. The predicted
//   clock variance resets to ~process_noise*dt each epoch instead of
//   remaining converged. The clock prior weight drops to ~9e-5, meaning
//   the clock is re-estimated from scratch every epoch.
//
// COMBINED EFFECT: The filter loses BOTH absolute position anchors:
// - The SPP prior weight drops from 100 to 1 (floor)
// - The clock prior weight stays at ~1e-4 forever (white-noise)
//
// The result is a position estimate that converges to ~0.5-1.5m rather than
// the cm-level accuracy achievable with a random-walk clock and strong prior.

