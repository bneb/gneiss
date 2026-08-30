//! Single Point Positioning (SPP) Weighted Non-Linear Least Squares Estimator.

pub mod measurements;
pub mod raim;
pub mod solver;

#[cfg(test)]
mod tests;

use gneiss_core::atmosphere::KlobucharParams;
use gneiss_core::coords::Coordinate;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_core::time::GpsTime;
use serde::{Deserialize, Serialize};

pub use measurements::build_measurements;
pub use solver::spp_wnlls_step;
use self::raim::{apply_raim, seed_initial_state};
use self::solver::solve_spp_iteratively;

pub(crate) const LIGHT_SPEED: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
pub(crate) const OMEGA_E: f64 = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S;
pub(crate) const MIN_WEIGHT: f64 = 1e-10;
pub(crate) const HEIGHT_CONSTRAINT_VAR: f64 = 100.0;
pub(crate) const MIN_EARTH_RADIUS_M: f64 = 6_000_000.0;
pub(crate) const MIN_ATMOSPHERE_ELEVATION_RAD: f64 = 5.0 * core::f64::consts::PI / 180.0;

/// Represents the current estimated state of the receiver.
#[derive(Debug, Clone, PartialEq)]
pub struct SppState {
    /// Receiver position in a specific Datum and Frame.
    pub position: Coordinate,
    /// Receiver clock bias in meters (c * dt) for GPS.
    pub cdt: f64,
    /// Receiver clock bias in meters for Galileo.
    pub cdt_gal: f64,
    /// Receiver clock bias in meters for BeiDou.
    pub cdt_bds: f64,
    /// Receiver clock bias in meters for GLONASS.
    pub cdt_glo: f64,
}

impl SppState {
    pub fn new(position: Coordinate, cdt: f64, cdt_gal: f64, cdt_bds: f64, cdt_glo: f64) -> Self {
        Self {
            position,
            cdt,
            cdt_gal,
            cdt_bds,
            cdt_glo,
        }
    }

    pub fn get_cdt(&self, constellation: gneiss_core::sat::Constellation) -> f64 {
        match constellation {
            gneiss_core::sat::Constellation::Galileo => self.cdt_gal,
            gneiss_core::sat::Constellation::Beidou => self.cdt_bds,
            gneiss_core::sat::Constellation::Glonass => self.cdt_glo,
            _ => self.cdt,
        }
    }
}

/// A single satellite measurement for use in the SPP estimator.
#[derive(Debug, Clone)]
pub struct SppMeasurement {
    pub constellation: gneiss_core::sat::Constellation,
    pub raw_pr: f64,
    pub snr: f64,
    pub doppler: f64,
    pub time: GpsTime,
    pub eph: Ephemeris,
    pub is_iono_free: bool,
    pub freq_band: u8,
}

/// Errors that can occur during an SPP WNLLS step.
#[derive(Debug, Clone, PartialEq)]
pub enum SppError {
    /// Not enough measurements to solve the state (need at least 4).
    NotEnoughMeasurements,
    /// The matrix inversion failed (singular matrix, poor geometry).
    MatrixInversionFailed,
    /// Maximum iterations reached without convergence.
    ConvergenceFailed,
    /// Geometry exploded or was too poor to use.
    PoorGeometry,
}

/// Configuration for the Single Point Positioning solver.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SppConfig {
    pub max_iterations: usize,
    pub convergence_threshold: f64,
    pub geometry_variance_threshold: f64,
    pub enable_sagnac: bool,
    pub enable_tropo: bool,
    pub enable_iono: bool,
    pub raim_outlier_m: f64,
    pub snr_a: f64,
    pub snr_b: f64,
    pub min_measurements_init: usize,
    pub raim_mad_multiplier: f64,
    pub elevation_mask_rad: f64,
}

impl Default for SppConfig {
    fn default() -> Self {
        Self {
            max_iterations: 15,
            convergence_threshold: 1e-4,
            geometry_variance_threshold: 10000.0,
            enable_sagnac: false,
            enable_tropo: false,
            enable_iono: false,
            raim_outlier_m: 50.0,
            snr_a: 1.0,
            snr_b: 150.0,
            min_measurements_init: 3,
            raim_mad_multiplier: 7.413,
            elevation_mask_rad: 0.261799,
        }
    }
}

pub fn compute_spp(
    epoch: &EpochObs,
    ephemerides: &[Ephemeris],
    iono_params: Option<&KlobucharParams>,
    config: &SppConfig,
    prev_state: Option<&SppState>,
) -> Result<SppState, SppError> {
    let measurements = build_measurements(epoch, ephemerides, config);
    if measurements.len() < config.min_measurements_init {
        return Err(SppError::NotEnoughMeasurements);
    }
    let seed_state = seed_initial_state(&measurements, prev_state);
    let state = solve_spp_iteratively(seed_state.clone(), &measurements, iono_params, config)?;
    apply_raim(state, seed_state, &measurements, iono_params, config)
}
