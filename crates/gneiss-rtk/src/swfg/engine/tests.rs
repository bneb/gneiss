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
