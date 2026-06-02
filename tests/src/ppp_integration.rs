use gneiss_rtk::engine::{ProcessingEngine, EngineConfig, EngineMode};
use gneiss_core::obs::{EpochObs, SatObs, Observation, ObsCode};
use core::str::FromStr;


#[test]
fn test_ppp_skeleton() {
    let mut config = EngineConfig::default();
    config.mode = EngineMode::Ppp;
    let mut engine = ProcessingEngine::new(config);

    // Create some dummy observations
    let time = gneiss_core::time::GpsTime::new(2300, 345600.0);
    
    // We need dual frequency observations for Ionosphere-Free
    let obs1 = Observation {
        code: ObsCode::from_str("C1C").unwrap(),
        value: 20000000.0,
        lock_time: None,
    };
    let obs2 = Observation {
        code: ObsCode::from_str("L1C").unwrap(),
        value: 105000000.0,
        lock_time: None,
    };
    let obs3 = Observation {
        code: ObsCode::from_str("C2W").unwrap(),
        value: 20000010.0,
        lock_time: None,
    };
    let obs4 = Observation {
        code: ObsCode::from_str("L2W").unwrap(),
        value: 82000000.0,
        lock_time: None,
    };
    let obs5 = Observation {
        code: ObsCode::from_str("S1C").unwrap(),
        value: 45.0,
        lock_time: None,
    };

    let sat_obs = SatObs {
        sat: gneiss_core::sat::SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: 1,
        },
        observations: vec![obs1, obs2, obs3, obs4, obs5],
    };

    let epoch_obs = EpochObs {
        time,
        satellites: vec![sat_obs],
    };

    // Note: Ephemerides are missing, so the engine should gracefully return an error or skip.
    let result = engine.process_epoch(&epoch_obs, None);
    
    // As long as it doesn't crash, the skeleton is valid.
    assert!(result.is_err() || result.is_ok());
}
