//! Qinertia-Grade Offline Post-Processing Engine for RTK and PPK.
//!
//! Orchestrates the complete 4-pass offline pipeline:
//! 1. Screening & Quality Control (GF/MW/Doppler cycle slip detection, ZUPT, base refinement)
//! 2. Forward Trajectory Filtering + AR
//! 3. Backward Trajectory Filtering + AR
//! 4. Optimal Bidirectional Fusion + Covariance Intersection Smoothing

pub mod backward;
pub mod combiner;
pub mod forward;
pub mod quality;
pub mod screening;

use nalgebra::Vector3;

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;

use crate::swfg::config::EngineConfig;
use crate::swfg::imu_preintegration::ImuSample;

pub use combiner::SmoothedEpoch;
pub use quality::QualityReport;
pub use screening::ScreeningReport;

/// Complete result of a post-processing run.
#[derive(Debug, Clone)]
pub struct PostProcessResult {
    /// Screen report from Pass 1.
    pub screening: ScreeningReport,
    /// Final smoothed trajectory from Pass 4.
    pub trajectory: Vec<SmoothedEpoch>,
    /// Statistical quality assessment.
    pub quality: QualityReport,
}

/// Post-processing configuration options.
#[derive(Debug, Clone, Default)]
pub struct PostProcessOptions {
    pub enable_bidirectional: bool,
    pub base_position: Option<Vector3<f64>>,
    pub initial_rover_position: Option<Vector3<f64>>,
    pub klobuchar_alpha: Option<[f64; 4]>,
    pub klobuchar_beta: Option<[f64; 4]>,
}

/// Execute the complete offline post-processing pipeline.
pub fn execute_post_process(
    config: &EngineConfig,
    ephemerides: &[Ephemeris],
    rover_epochs: &[EpochObs],
    base_epochs: Option<&[EpochObs]>,
    imu_samples: Option<&[ImuSample]>,
    options: &PostProcessOptions,
) -> Result<PostProcessResult, String> {
    if rover_epochs.is_empty() {
        return Err("No rover epochs provided".to_string());
    }

    // Pass 1: Screening & Quality Control
    let screening = screening::screen_dataset(rover_epochs, base_epochs, options.base_position);
    let base_pos = options.base_position.or(screening.refined_base_pos);

    let klob = match (options.klobuchar_alpha, options.klobuchar_beta) {
        (Some(a), Some(b)) => Some((a, b)),
        _ => None,
    };

    // Pass 2: Forward Pass
    let forward_traj = forward::run_forward_pass(
        config, ephemerides, klob, rover_epochs, base_epochs, base_pos, imu_samples, options.initial_rover_position,
    );

    // Pass 3: Backward Pass (if enabled)
    let backward_map = if options.enable_bidirectional {
        let initial_rover_pos = forward_traj.last().map(|e| e.position_ecef);
        backward::run_backward_pass(
            config, ephemerides, klob, rover_epochs, base_epochs, base_pos, imu_samples, initial_rover_pos,
        )
    } else {
        std::collections::BTreeMap::new()
    };

    // Pass 4: Optimal Bidirectional Fusion
    let smoothed_traj = combiner::combine_trajectories(&forward_traj, &backward_map);
    let quality_rep = quality::generate_quality_report(&smoothed_traj);

    Ok(PostProcessResult {
        screening,
        trajectory: smoothed_traj,
        quality: quality_rep,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::time::GpsTime;

    #[test]
    fn test_post_process_empty_error() {
        let config = EngineConfig::Spp(Default::default());
        let res = execute_post_process(&config, &[], &[], None, None, &PostProcessOptions::default());
        assert!(res.is_err());
    }

    #[test]
    fn test_post_process_minimal() {
        let config = EngineConfig::Spp(Default::default());
        let rover = vec![EpochObs { time: GpsTime::new(2000, 100.0), satellites: Vec::new() }];
        let options = PostProcessOptions { enable_bidirectional: false, ..Default::default() };
        let res = execute_post_process(&config, &[], &rover, None, None, &options);
        assert!(res.is_ok());
    }
}
