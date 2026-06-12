use gneiss_core::obs::EpochObs;
use crate::filter::RtkState;
use crate::engine::{EngineError, ProcessingEngine};

pub fn process_tight_fg<'a>(engine: &'a mut ProcessingEngine, rover_obs: &EpochObs) -> Result<&'a RtkState, EngineError> {
    engine.process_spp(rover_obs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::obs::EpochObs;
    use crate::engine::ProcessingEngine;
    use crate::engine::EngineConfig;
    use gneiss_core::time::GpsTime;

    #[test]
    fn test_process_tight_fg() {
        let config = EngineConfig::default();
        let mut engine = ProcessingEngine::new(config);
        
        // This is a minimal valid test for the dummy method to catch the `Err` mutation
        let obs = EpochObs {
            time: GpsTime::new(0, 0.0),
            satellites: Vec::new(),
        };
        let res = process_tight_fg(&mut engine, &obs);
        assert!(res.is_err());
    }
}
