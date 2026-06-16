use crate::engine::{EngineError, ProcessingEngine};
use gneiss_core::obs::EpochObs;
use crate::engine::tight_fg::TightFactorGraph;
use crate::filter::RtkState;

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
    
    engine.current_state.as_ref().ok_or(EngineError::StateDisappeared)
}
