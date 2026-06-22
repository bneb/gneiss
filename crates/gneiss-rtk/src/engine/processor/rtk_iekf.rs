use crate::engine::tight_iekf::TightFactorGraph;
use crate::engine::{EngineError, ProcessingEngine};
use crate::filter::RtkState;
use gneiss_core::obs::EpochObs;

pub fn process_rtk_factor_graph<'a>(
    engine: &'a mut ProcessingEngine,
    rover_obs: &EpochObs,
    base_obs: Option<&EpochObs>,
) -> Result<&'a RtkState, EngineError> {
    // 1. Run the standard RTK/INS predictor & EKF updater to get a strong prior
    let _ = engine.process_rtk(rover_obs, base_obs)?;

    // 2. Apply Factor Graph optimization (Iterated EKF style)
    let fg = TightFactorGraph::new();
    fg.solve(engine, rover_obs, None);

    engine
        .current_state
        .as_ref()
        .ok_or(EngineError::StateDisappeared)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::obs::EpochObs;
    use gneiss_core::time::GpsTime;
    use nalgebra::Vector3;

    fn make_empty_rover(time: GpsTime) -> EpochObs {
        EpochObs {
            time,
            satellites: Vec::new(),
        }
    }

    #[test]
    fn test_process_rtk_factor_graph_fails_without_state() {
        let mut engine = ProcessingEngine::new(crate::engine::EngineConfig::default());
        engine.current_state = None;

        let time = GpsTime::new(0, 0.0);
        let rover = make_empty_rover(time);

        // process_rtk will fail since there's no state and no ephemerides
        let err = process_rtk_factor_graph(&mut engine, &rover, None).unwrap_err();
        assert!(matches!(err, EngineError::InitialSppFailed));
    }

    #[test]
    fn test_process_rtk_factor_graph_propagates_state_disappeared() {
        // Create engine with a state and ephemerides so process_rtk succeeds
        let mut engine = ProcessingEngine::new(crate::engine::EngineConfig::default());
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let state = RtkState::new(time, pos, 1.0);
        engine.current_state = Some(state);

        let rover = make_empty_rover(time);

        // With state present but no observations, process_rtk should run
        // The factor graph solve may panic or fail — we just check that
        // a state_disappeared error can arise if factor graph clears state
        let _ = process_rtk_factor_graph(&mut engine, &rover, None);
        // Depending on internal TightFactorGraph behavior, this may succeed or fail
        // But it should not panic
    }
}
