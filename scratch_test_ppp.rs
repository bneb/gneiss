#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{EngineConfig, EngineMode, ProcessingEngine, EngineError};
    use crate::filter::RtkState;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use gneiss_core::sat::{SatelliteId, Constellation};
    use gneiss_core::ephemeris::{Ephemeris, GpsEphemeris};
    use gneiss_core::obs::{EpochObs, SatObs, Observation, ObsCode, SignalCode, ObsType};
    use nalgebra::{Vector3, DMatrix, DVector};

    #[test]
    fn test_snr_scale_mutants() {
        assert_eq!(snr_scale(45), 1.0);
        assert_eq!(snr_scale(35), 10.0);
        assert_eq!(snr_scale(25), 100.0);
    }

    fn create_mock_obs() -> EpochObs {
        EpochObs {
            time: GpsTime::new(2000, 0.0),
            satellites: vec![],
        }
    }

    #[test]
    fn test_process_ppp_valid_pos_fallback() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::Ppp;
        let mut engine = ProcessingEngine::new(config);
        let obs = create_mock_obs();

        let res = process_ppp(&mut engine, &obs);
        assert!(res.is_err()); 

        let time = GpsTime::new(2000, 0.0);
        engine.current_state = Some(RtkState::new(time, Coordinate::new(Vector3::new(10.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time), 1.0));
        let res2 = process_ppp(&mut engine, &obs);
        assert!(res2.is_ok());

        engine.current_state = Some(RtkState::new(time, Coordinate::new(Vector3::new(f64::NAN, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time), 1.0));
        let res3 = process_ppp(&mut engine, &obs);
        assert!(res3.is_ok());

        engine.current_state = Some(RtkState::new(time, Coordinate::new(Vector3::new(2000.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time), 1.0));
        let res4 = process_ppp(&mut engine, &obs);
        assert_eq!(res4.unwrap_err(), EngineError::InsufficientSatellites);
    }

    fn mock_eph(sat: SatelliteId) -> Ephemeris {
        Ephemeris::Gps(GpsEphemeris {
            sat,
            toe: GpsTime::new(2000, 0.0),
            toc: GpsTime::new(2000, 0.0),
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.0, sqrt_a: 5153.0, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.0, idot: 0.0, omega: 0.0, tgd: 0.0,
            iode: 0, iodc: 0,
        })
    }

    fn pseudo_obs(freq_band: u8, val: f64) -> Observation {
        Observation {
            code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band, attribute: 'C' } },
            value: val,
            lock_time: None,
        }
    }
    
    fn phase_obs(freq_band: u8, val: f64) -> Observation {
        Observation {
            code: ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band, attribute: 'C' } },
            value: val,
            lock_time: Some(10),
        }
    }

    #[test]
    fn test_process_ppp_sat_skipping() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::Ppp;
        let mut engine = ProcessingEngine::new(config);
        
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(6378000.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        engine.current_state = Some(RtkState::new(time, pos, 1.0));
        
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        engine.ephemerides.push(mock_eph(sat1));

        let mut obs = create_mock_obs();
        obs.satellites.push(SatObs {
            sat: sat1,
            observations: vec![
                pseudo_obs(1, 20000000.0),
            ]
        });
        assert_eq!(process_ppp(&mut engine, &obs).unwrap_err(), EngineError::InsufficientSatellites);

        let sat2 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        obs.satellites[0].sat = sat2; // no eph
        obs.satellites[0].observations.push(pseudo_obs(2, 20000000.0));
        assert_eq!(process_ppp(&mut engine, &obs).unwrap_err(), EngineError::InsufficientSatellites);

        obs.satellites[0].sat = sat1; // valid eph
        engine.current_state.as_mut().unwrap().position.vector = Vector3::new(-6378000.0, 0.0, 0.0); // el < 15
        assert_eq!(process_ppp(&mut engine, &obs).unwrap_err(), EngineError::InsufficientSatellites);
    }

    #[test]
    fn test_process_ppp_full_observation() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::Ppp;
        let mut engine = ProcessingEngine::new(config);
        
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(6378000.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        engine.current_state = Some(RtkState::new(time, pos, 1.0));
        
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        engine.ephemerides.push(mock_eph(sat1));

        let mut obs = create_mock_obs();
        obs.satellites.push(SatObs {
            sat: sat1,
            observations: vec![
                pseudo_obs(1, 40000000.0),
                pseudo_obs(2, 40000000.0),
                phase_obs(1, 100000000.0),
                phase_obs(2, 70000000.0),
            ]
        });

        engine.current_state.as_mut().unwrap().gf_values.insert(sat1, 10.0);
        
        let res = process_ppp(&mut engine, &obs);
        assert!(res.is_ok());
        
        let state = res.unwrap();
        assert!(state.ambiguity_keys.contains(&(sat1, 0)));
        assert_eq!(*state.locktimes.get(&(sat1, 1)).unwrap(), 10);
    }

    #[test]
    fn test_process_ppp_clock_jump() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::Ppp;
        let mut engine = ProcessingEngine::new(config);
        
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(6378000.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        engine.current_state = Some(RtkState::new(time, pos, 1.0));
        
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        engine.ephemerides.push(mock_eph(sat1));

        let mut obs = create_mock_obs();
        obs.satellites.push(SatObs {
            sat: sat1,
            observations: vec![
                pseudo_obs(1, 40000000.0),
                pseudo_obs(2, 40000000.0),
                phase_obs(1, 100000000.0),
                phase_obs(2, 70000000.0),
            ]
        });

        let res = process_ppp(&mut engine, &obs);
        assert!(res.is_ok());
        let state = res.unwrap();
        assert!(state.rcv_clk_bias > 10_000_000.0);
    }

    #[test]
    fn test_process_ppp_slip_boundary() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::Ppp;
        let mut engine = ProcessingEngine::new(config);
        
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(Vector3::new(6378000.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time);
        engine.current_state = Some(RtkState::new(time, pos, 1.0));
        
        let sat1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        engine.ephemerides.push(mock_eph(sat1));

        let mut obs = create_mock_obs();
        obs.satellites.push(SatObs {
            sat: sat1,
            observations: vec![
                pseudo_obs(1, 40000000.0),
                pseudo_obs(2, 40000000.0),
                phase_obs(1, 100.0),
                phase_obs(2, 100.0),
            ]
        });

        let res1 = process_ppp(&mut engine, &obs);
        assert!(res1.is_ok());
        
        obs.satellites[0].observations[2].lock_time = Some(11);
        obs.satellites[0].observations[3].lock_time = Some(11); 

        obs.satellites[0].observations[2].value = 100.0 + 0.04 / 0.19029367279836488;
        let res2 = process_ppp(&mut engine, &obs);
        assert!(res2.is_ok());
        assert!(engine.current_state.as_ref().unwrap().ambiguity_keys.contains(&(sat1, 0)), "Should not have slipped");

        obs.satellites[0].observations[2].value = 100.0 + 0.06 / 0.19029367279836488;
        obs.satellites[0].observations[2].lock_time = Some(12);
        obs.satellites[0].observations[3].lock_time = Some(12);
        
        let res3 = process_ppp(&mut engine, &obs);
        assert!(res3.is_ok());
        assert!(!engine.current_state.as_ref().unwrap().ambiguity_keys.contains(&(sat1, 0)), "Should have slipped");
    }
}
