#![allow(clippy::unwrap_used)]

    use super::*;
    use crate::swfg::config::SppConfig;
    use gneiss_core::ephemeris::Ephemeris;
    use gneiss_core::sat::{Constellation, SatelliteId};

    fn make_gps_eph(sat: SatelliteId, toc: GpsTime, m0: f64) -> Ephemeris {
        use gneiss_core::ephemeris::GpsEphemeris;
        Ephemeris::Gps(GpsEphemeris {
            sat, toe: toc, toc, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0, e: 0.001, sqrt_a: 26_560_000.0_f64.sqrt(),
            delta_n: 0.0, omega0: 0.0, omega_dot: 0.0, i0: 0.96, idot: 0.0,
            omega: 0.0, tgd: 0.0, iode: 0, iodc: 0,
        })
    }

    fn make_epoch(time: GpsTime, n_sats: usize) -> EpochObs {
        use gneiss_core::obs::{ObsCode, ObsType, Observation, SignalCode};
        let mut obs = EpochObs { time, satellites: Vec::new() };
        for prn in 1..=(n_sats as u8) {
            let sat = SatelliteId { constellation: Constellation::Gps, prn };
            let sat_obs = gneiss_core::obs::SatObs {
                sat,
                observations: vec![Observation {
                    code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 1, attribute: 'C' } },
                    value: 20_000_000.0 + (prn as f64) * 100.0,
                    lock_time: None, lli: None,
                }],
            };
            obs.satellites.push(sat_obs);
        }
        obs
    }

    #[test]
    fn epoch_to_raw_obs_returns_observations() {
        let time = GpsTime::new(2200, 100.0);
        let ephs: Vec<Ephemeris> = (1..=20).map(|prn| {
            let m0 = (prn as f64 - 1.0) * std::f64::consts::TAU / 20.0;
            make_gps_eph(SatelliteId { constellation: Constellation::Gps, prn }, time, m0)
        }).collect();

        let config = EngineConfig::Spp(SppConfig::default());
        let engine = SwfgEngine::new(&config, ephs);
        let rover = make_epoch(time, 16);
        let raw = engine.epoch_to_raw_obs(&rover).unwrap();
        assert!(raw.len() >= 4);
        assert!(raw[0].pr_l1 > 0.0);
    }

    #[test]
    fn engine_processes_single_epoch() {
        let time = GpsTime::new(2200, 100.0);
        let ephs: Vec<Ephemeris> = (1..=20).map(|prn| {
            let m0 = (prn as f64 - 1.0) * std::f64::consts::TAU / 20.0;
            make_gps_eph(SatelliteId { constellation: Constellation::Gps, prn }, time, m0)
        }).collect();

        let r = gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M;
        let config = EngineConfig::Spp(SppConfig { initial_position: Some([r, 0.0, 0.0]), ..SppConfig::default() });
        let mut engine = SwfgEngine::new(&config, ephs);
        let rover = make_epoch(time, 16);
        let sol = engine.process_epoch(&rover).unwrap();
        assert!(sol.n_satellites >= 4);
        assert!(sol.position_ecef.norm() > 1e6);
    }

    /// Regression test for the IfbGlonass orphan-variable bug
    /// (docs/SOLVER_MODE_MATRIX.md): a GLONASS satellite tracked in the
    /// observation file with NO matching ephemeris (missing/partial nav
    /// data, e.g. a GPS-only broadcast file against multi-GNSS
    /// observations) must not crash rover-only/PPP processing. Before the
    /// fix, the IfbGlonass variable was created from raw satellite
    /// tracking regardless of whether the satellite ever reached a
    /// processed observation, so it got zero factors and every solve in
    /// the session failed with OrphanVariable, forever.
    #[test]
    fn engine_processes_epoch_with_untracked_glonass_satellite() {
        let time = GpsTime::new(2200, 100.0);
        // Only GPS ephemerides -- deliberately no GLONASS entry, so the
        // GLONASS satellite below cannot survive ephemeris matching.
        let ephs: Vec<Ephemeris> = (1..=20).map(|prn| {
            let m0 = (prn as f64 - 1.0) * std::f64::consts::TAU / 20.0;
            make_gps_eph(SatelliteId { constellation: Constellation::Gps, prn }, time, m0)
        }).collect();

        let r = gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M;
        let config = EngineConfig::Spp(SppConfig { initial_position: Some([r, 0.0, 0.0]), ..SppConfig::default() });
        let mut engine = SwfgEngine::new(&config, ephs);

        let mut rover = make_epoch(time, 16);
        use gneiss_core::obs::{ObsCode, ObsType, Observation, SignalCode};
        rover.satellites.push(gneiss_core::obs::SatObs {
            sat: SatelliteId { constellation: Constellation::Glonass, prn: 1 },
            observations: vec![Observation {
                code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 1, attribute: 'C' } },
                value: 20_000_100.0,
                lock_time: None, lli: None,
            }],
        });

        let sol = engine.process_epoch(&rover);
        assert!(sol.is_ok(), "untracked-ephemeris GLONASS satellite must not break the solve: {sol:?}");
    }

    #[test]
    fn engine_processes_rtk_epoch() {
        use crate::swfg::config::RtkConfig;
        let time = GpsTime::new(2200, 100.0);
        let ephs: Vec<Ephemeris> = (1..=20).map(|prn| {
            let m0 = (prn as f64 - 1.0) * std::f64::consts::TAU / 20.0;
            make_gps_eph(SatelliteId { constellation: Constellation::Gps, prn }, time, m0)
        }).collect();

        let r = gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M;
        let base_pos = Vector3::new(r, 0.0, 0.0);
        let config = EngineConfig::Rtk(RtkConfig { initial_position: Some([r + 10.0, 10.0, 0.0]), ..RtkConfig::default() });
        let mut engine = SwfgEngine::new(&config, ephs);
        let rover = make_epoch(time, 16);
        let base = make_epoch(time, 16);
        let sol = engine.process_rtk_epoch(&rover, &base, base_pos).unwrap();
        assert!(sol.n_satellites >= 4);
    }

    #[test]
    fn engine_processes_ppp_epoch() {
        use crate::swfg::config::PppConfig;
        let time = GpsTime::new(2200, 100.0);
        let ephs: Vec<Ephemeris> = (1..=20).map(|prn| {
            let m0 = (prn as f64 - 1.0) * std::f64::consts::TAU / 20.0;
            make_gps_eph(SatelliteId { constellation: Constellation::Gps, prn }, time, m0)
        }).collect();

        let r = gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M;
        let config = EngineConfig::Ppp(PppConfig { initial_position: Some([r, 0.0, 0.0]), is_kinematic: false, ..PppConfig::default() });
        let mut engine = SwfgEngine::new(&config, ephs);
        let rover = make_epoch(time, 16);
        let sol = engine.process_epoch(&rover);
        println!("PPP test sol: {:?}", sol);
        assert!(sol.is_ok(), "PPP solve failed: {:?}", sol.err());
    }

    /// 20 GPS ephemerides spread around the orbit at `time`.
    fn gps_scene(time: GpsTime) -> Vec<Ephemeris> {
        (1..=20)
            .map(|prn| {
                let m0 = (prn as f64 - 1.0) * std::f64::consts::TAU / 20.0;
                make_gps_eph(SatelliteId { constellation: Constellation::Gps, prn }, time, m0)
            })
            .collect()
    }

    fn spp_engine(time: GpsTime) -> SwfgEngine {
        let r = gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M;
        let config = EngineConfig::Spp(SppConfig { initial_position: Some([r, 0.0, 0.0]), ..SppConfig::default() });
        SwfgEngine::new(&config, gps_scene(time))
    }

    #[test]
    fn current_position_defaults_to_the_earth_equator_before_any_epoch() {
        let engine = spp_engine(GpsTime::new(2200, 100.0));
        assert_eq!(
            engine.get_current_position(),
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0)
        );
    }

    #[test]
    fn is_ppp_follows_the_selected_mode() {
        let time = GpsTime::new(2200, 100.0);
        let spp = SwfgEngine::new(&EngineConfig::Spp(SppConfig::default()), gps_scene(time));
        assert!(!spp.is_ppp());
        use crate::swfg::config::PppConfig;
        let ppp = SwfgEngine::new(&EngineConfig::Ppp(PppConfig::default()), gps_scene(time));
        assert!(ppp.is_ppp());
    }

    #[test]
    fn engine_solves_a_sequence_of_epochs_and_keeps_the_solution_finite() {
        // SPP runs a one-epoch window, so every epoch after the first also
        // exercises marginalization of the epoch that just left the window.
        let t0 = GpsTime::new(2200, 100.0);
        let mut engine = spp_engine(t0);
        for k in 0..12u32 {
            let t = GpsTime::new(2200, 100.0 + f64::from(k));
            let sol = engine.process_epoch(&make_epoch(t, 16)).unwrap_or_else(|e| panic!("epoch {k} failed: {e}"));
            assert_eq!(sol.time, t, "solution must be stamped with the epoch time");
            assert!(sol.n_satellites >= 4, "epoch {k} kept only {} satellites", sol.n_satellites);
            assert!(sol.position_ecef.iter().all(|v| v.is_finite()), "epoch {k} produced {}", sol.position_ecef);
            assert!(
                (5.0e6..7.0e6).contains(&sol.position_ecef.norm()),
                "epoch {k} position {} is not on Earth", sol.position_ecef.norm()
            );
        }
        let final_pos = engine.get_current_position();
        assert!((5.0e6..7.0e6).contains(&final_pos.norm()));
    }

    #[test]
    fn rtk_engine_solves_a_sequence_of_epochs() {
        use crate::swfg::config::RtkConfig;
        let t0 = GpsTime::new(2200, 100.0);
        let r = gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M;
        let base_pos = Vector3::new(r, 0.0, 0.0);
        let config = EngineConfig::Rtk(RtkConfig { initial_position: Some([r + 10.0, 10.0, 0.0]), ..RtkConfig::default() });
        let mut engine = SwfgEngine::new(&config, gps_scene(t0));
        for k in 0..8u32 {
            let t = GpsTime::new(2200, 100.0 + f64::from(k));
            let sol = engine
                .process_rtk_epoch(&make_epoch(t, 16), &make_epoch(t, 16), base_pos)
                .unwrap_or_else(|e| panic!("RTK epoch {k} failed: {e}"));
            assert!(sol.position_ecef.iter().all(|v| v.is_finite()));
            assert!((5.0e6..7.0e6).contains(&sol.position_ecef.norm()), "epoch {k}: {}", sol.position_ecef.norm());
        }
    }

    #[test]
    fn engine_couples_preintegrated_imu_into_an_rtk_epoch() {
        use crate::swfg::config::RtkConfig;
        use crate::swfg::imu_preintegration::ImuPreintegration;
        let t0 = GpsTime::new(2200, 100.0);
        let r = gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M;
        let base_pos = Vector3::new(r, 0.0, 0.0);
        let config = EngineConfig::RtkIns(crate::swfg::config::RtkInsConfig {
            rtk: RtkConfig { initial_position: Some([r + 10.0, 10.0, 0.0]), ..RtkConfig::default() },
            imu: crate::swfg::config::ImuConfig::default(),
        });
        let mut engine = SwfgEngine::new(&config, gps_scene(t0));

        // One second of rest: the accelerometer senses the specific force that
        // cancels gravity, so a stationary platform integrates to zero.
        let mut preint = ImuPreintegration::new();
        preint.dt = 1.0;
        preint.is_stationary = true;
        for k in 0..4u32 {
            let t = GpsTime::new(2200, 100.0 + f64::from(k));
            let sol = engine
                .process_rtk_epoch_with_imu(&make_epoch(t, 16), &make_epoch(t, 16), base_pos, Some(preint.clone()))
                .unwrap_or_else(|e| panic!("INS epoch {k} failed: {e}"));
            assert!(sol.position_ecef.iter().all(|v| v.is_finite()), "epoch {k}: {}", sol.position_ecef);
            assert!((5.0e6..7.0e6).contains(&sol.position_ecef.norm()), "epoch {k}: {}", sol.position_ecef.norm());
        }
    }

    #[test]
    fn engine_rejects_an_epoch_with_too_few_usable_observations() {
        // Fewer than four satellites cannot determine a position, and without
        // IMU there is nothing to fall back on.
        let t = GpsTime::new(2200, 100.0);
        let mut engine = spp_engine(t);
        let err = engine.process_epoch(&make_epoch(t, 2)).expect_err("two satellites must be rejected");
        assert!(err.contains("too few observations"), "unexpected error: {err}");
    }
