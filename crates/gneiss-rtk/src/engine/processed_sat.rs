use gneiss_core::obs::SatObs;
use nalgebra::Vector3;

/// Encapsulates all intermediate calculated geometry and observables for a single satellite.
/// This prevents monolithic factor graph functions from juggling 18+ local variables in loops.
#[derive(Clone)]
pub struct ProcessedSat<'a> {
    pub sat_obs: &'a SatObs,
    pub dt_sat_m: f64,
    pub p_meas: f64,
    pub is_iono_free: bool,
    pub cp1: Option<f64>,
    pub cp2: Option<f64>,
    pub los: Vector3<f64>,
    pub dist: f64,
    pub el: f64,
    pub snr: f64,
    pub doppler: f64,
    pub lam1: f64,
    pub lam2: f64,
    pub tropo_dry: f64,
    pub map_wet: f64,
    pub iono_delay: f64,
    pub f1: f64,
    pub f2: f64,
    pub sat_pos_rot: Vector3<f64>,
    pub sat_vel: Vector3<f64>,
    pub sat_clock_drift: f64,
    pub rcv_pos_ecef: Vector3<f64>,
    pub pcv_correction: f64,
}
