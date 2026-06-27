//! IF→UDUC mode switching for the RTKLIB PPP solver.
//!
//! The hybrid approach: IF mode for coarse convergence + AR (~50 epochs),
//! then warm-handoff to UDUC for sub-meter accuracy. UDUC separates
//! ionosphere from geometry, eliminating the 7m vertical IF code bias.

use std::collections::HashMap;

use gneiss_core::sat::SatelliteId;
use nalgebra::Vector3;

/// PPP solver mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PppSolverMode {
    /// IF mode — default
    If,
    /// UDUC mode — active after successful handoff
    Uduc,
}

/// Saved IF state for UDUC handoff or fallback recovery.
#[derive(Clone)]
pub struct IfSnapshot {
    pub position: Vector3<f64>,
    pub clock_bias: f64,
    pub zwd: f64,
    /// IF bias values (meters) for AR-fixed satellites
    pub bias_values: Vec<f64>,
    /// Diagonal of IF bias covariance
    pub bias_vars: Vec<f64>,
    /// Satellite IDs for bias slots
    pub sat_ids: Vec<SatelliteId>,
    /// MW widelane EMA
    pub mw_ema: HashMap<SatelliteId, (u32, f64)>,
    pub epoch: u32,
}

/// Tracks IF→UDUC mode state and saved snapshots.
#[derive(Clone)]
pub struct PppModeManager {
    pub mode: PppSolverMode,
    pub if_snapshot: Option<IfSnapshot>,
    pub uduc_epoch_count: u32,
    pub fallback_count: u32,
    pub if_epoch_count: u32,
}

impl Default for PppModeManager {
    fn default() -> Self {
        Self {
            mode: PppSolverMode::If,
            if_snapshot: None,
            uduc_epoch_count: 0,
            fallback_count: 0,
            if_epoch_count: 0,
        }
    }
}
