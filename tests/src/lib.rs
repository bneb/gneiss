#![cfg_attr(test, allow(clippy::unwrap_used))]

#[cfg(test)]
pub mod benchmark_matrix;
pub mod post_process_simulation;

#[cfg(test)]
mod integration {
    use gneiss_rtk::swfg::config::EngineConfig;
    use gneiss_rtk::swfg::engine::SwfgEngine;

    #[test]
    fn test_swfg_engine_creation() {
        let config = EngineConfig::Spp(Default::default());
        let engine = SwfgEngine::new(&config, vec![]);
        let pos = engine.get_current_position();
        // Default position is on the equator at prime meridian
        assert!(pos.norm() > 1e6);
    }
}
