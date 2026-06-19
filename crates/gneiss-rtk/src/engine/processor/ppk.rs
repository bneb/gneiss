use super::ProcessingEngine;
use crate::engine::EngineError;
use crate::filter::RtkState;

impl ProcessingEngine {
    pub fn run_combined_ppk(&mut self) -> Result<Vec<RtkState>, EngineError> {
        crate::engine::smoother::run_combined_ppk(self)
    }
}
