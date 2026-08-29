//! Qinertia-Grade Offline Post-Processing Engine for RTK and PPK.
//!
//! Orchestrates the complete 4-pass offline pipeline:
//! 1. Screening & Quality Control (GF/MW/Doppler cycle slip detection, ZUPT, base refinement)
//! 2. Forward Trajectory Filtering + AR
//! 3. Backward Trajectory Filtering + AR
//! 4. Optimal Bidirectional Fusion + Covariance Intersection Smoothing

pub mod antenna;
pub mod backward;
pub mod combiner;
pub mod dynamics;
pub mod forward;
pub mod iekf_pass;
pub mod network;
pub mod quality;
pub mod screening;
pub mod sidereal;

use nalgebra::Vector3;

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;

use crate::swfg::config::EngineConfig;
use crate::swfg::imu_preintegration::ImuSample;

pub use combiner::SmoothedEpoch;
pub use dynamics::ProcessingDynamics;
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
    /// Acceleration random-walk process noise (m/s^2). Defaults to 1.0
    /// (kinematic). Static datasets must pass a small value (e.g. 1e-6):
    /// at 30 s epochs the position Q scales as dt^3/3 * q, so 1.0
    /// re-randomizes position ~55 m per epoch and prevents the float
    /// solution from converging.
    pub q_accel: Option<f64>,
    /// Opt-in Melbourne–Wübbena wide-lane cascade AR for long baselines.
    /// Default off (bit-identical legacy path); when on, it only affects
    /// epochs the joint FAR/PAR left float, and requires six confidently
    /// fixed wide-lane pairs before claiming a fix.
    pub widelane_ar: bool,
    /// Opt-in tropospheric N/E gradient states. Default off; the
    /// multi-GNSS profile enables them (measured: network fused
    /// 96.5 -> 97.5% with v-tail improvements, no regressions).
    pub tropo_gradients: bool,
    /// Network-solved satellite wide-lane UPDs (cycles): consumed by the
    /// MW tracker when `widelane_ar` is on. Produced by a Phase-A
    /// pre-pass (`mw::solve_network_upd`) over all bases.
    pub network_sat_upd: Option<std::collections::HashMap<u16, f64>>,
    /// Opt-in receiver antenna PCV correction `(rover, base)`. Default off
    /// (bit-identical legacy path); when set, the elevation-dependent
    /// differential receiver PCV is removed from DD carrier phase before
    /// ambiguity estimation. Arc-wrapped so cloning options stays cheap.
    pub receiver_pcv: Option<std::sync::Arc<ReceiverPcvPair>>,
    /// Rover motion model. Default [`ProcessingDynamics::Static`] keeps
    /// every legacy behaviour byte-identical; [`ProcessingDynamics::
    /// Kinematic`] selects mobile Q, no monument lock, widened
    /// innovation gates, and sigma-scaled combiner thresholds.
    pub dynamics: ProcessingDynamics,
    /// Opt-in GLONASS participation in DD formation and AR. Default off:
    /// code inter-channel biases don't cancel between receivers, so phase
    /// ambiguities absorb per-satellite constants but code ICBs still
    /// leak via float seeding (docs/NETWORK_RTK_NEXT_STEPS.md). Applied
    /// unconditionally regardless of baseline length or `widelane_ar` --
    /// previously this was set independently inside each of the forward
    /// and backward passes, nested inside unrelated baseline-length
    /// gating that differed between the two, so a >25km baseline with
    /// widelane_ar on could silently run forward and backward with
    /// DIFFERENT constellation sets. This field is now the single source
    /// of truth for both passes.
    pub enable_glonass: bool,
    /// Opt-in temporal honesty gate on each base's fused trajectory
    /// (`network::apply_continuity_gate_dynamics`): downgrades a fixed
    /// claim to float when it jumps farther than the motion model allows
    /// between adjacent epochs, the only independent signal against a
    /// forward/backward pair that confidently agrees on the same wrong
    /// integer (docs/NETWORK_RTK_NEXT_STEPS.md, "Combiner trust model").
    /// Default off and must stay opt-in: flipping it for a consumer with
    /// an existing bit-identical walkthrough reference is a deliberate,
    /// signed-off change (see "Walkthrough reference refresh"), not a
    /// side effect of wiring it up.
    pub continuity_gate: bool,
}

/// Paired receiver antenna PCV models consumed by the DD engine.
#[derive(Debug, Clone)]
pub struct ReceiverPcvPair {
    pub rover: std::sync::Arc<gneiss_parsers::receiver_antenna::ReceiverAntenna>,
    pub base: std::sync::Arc<gneiss_parsers::receiver_antenna::ReceiverAntenna>,
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
        config, ephemerides, klob, rover_epochs, base_epochs, base_pos, imu_samples, options.initial_rover_position, options.q_accel, options.dynamics, options.widelane_ar, options.tropo_gradients,
        options.network_sat_upd.clone(),
        options.receiver_pcv.clone(),
        options.enable_glonass,
    );

    // Pass 3: Backward Pass (if enabled)
    let backward_map = if options.enable_bidirectional {
        let initial_rover_pos = forward_traj.last().map(|e| e.position_ecef);
        backward::run_backward_pass(
            config, ephemerides, klob, rover_epochs, base_epochs, base_pos, imu_samples, initial_rover_pos, options.q_accel, options.dynamics, options.widelane_ar, options.tropo_gradients,
            options.network_sat_upd.clone(),
            options.receiver_pcv.clone(),
            options.enable_glonass,
        )
    } else {
        std::collections::BTreeMap::new()
    };

    // Pass 4: Optimal Bidirectional Fusion
    let smoothed_traj = combiner::combine_trajectories(
        &forward_traj,
        &backward_map,
        options.widelane_ar,
        options.dynamics,
    );
    let smoothed_traj = if options.continuity_gate {
        network::apply_continuity_gate_dynamics(
            smoothed_traj,
            network::CONTINUITY_MAX_DT_S,
            options.dynamics,
        )
    } else {
        smoothed_traj
    };
    let quality_rep = quality::generate_quality_report(&smoothed_traj);

    Ok(PostProcessResult {
        screening,
        trajectory: smoothed_traj,
        quality: quality_rep,
    })
}

/// Nearest-rank percentile of a pre-sorted slice (`q` in `[0.0, 1.0]`).
/// Shared by `quality.rs` and `sidereal/mod.rs`, which both need exactly
/// this convention -- not `gneiss_core::metrics::compute_statistics`'s
/// `ceil`-based p95/p99, a different (also valid) convention that
/// callers here were never using in the first place.
pub(crate) fn percentile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((sorted.len() as f64 * q).floor() as usize).min(sorted.len() - 1);
    sorted[idx]
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
