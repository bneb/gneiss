//! Dataset Generation and Loss Functions for ML
//!
//! This module provides functions to export GNSS measurements and FGO post-fit residuals
//! to a CSV file for training the RAIM GNN. It also implements numerically stable
//! Heteroscedastic Negative Log-Likelihood (NLL) loss functions using `candle_core`.

use candle_core::{Result, Tensor};
use std::collections::HashMap;
use std::fs::OpenOptions;
use std::io::Write;

/// Maximum absolute bounds for predicted log-variance to prevent exponential detonation
pub const MAX_LOG_VAR: f32 = 10.0;
/// Epsilon to prevent divide-by-zero during Softplus loss calculation
pub const EPSILON: f32 = 1e-6;

/// Exports an epoch of GNSS observations and their post-fit squared residuals for training the GNN RAIM model.
#[allow(clippy::too_many_arguments)]
pub fn export_epoch_to_csv(
    path: &str,
    epoch_num: u32,
    matched_obs: &[(crate::filter::DdObservation, crate::filter::DdObservation)],
    rov_llh: nalgebra::Vector3<f64>,
    pos_apc: nalgebra::Vector3<f64>,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    state_time: gneiss_core::time::GpsTime,
    innovations: &nalgebra::DVector<f64>,
    h_matrix: &nalgebra::DMatrix<f64>,
    valid_indices: &[usize], // Indices into innovations/H corresponding to the accepted measurements
    measurement_types: &[(gneiss_core::sat::SatelliteId, u8, f64)],
    dx: &nalgebra::DVector<f64>,
) {
    let sat_to_residual =
        extract_post_fit_residuals(innovations, h_matrix, dx, valid_indices, measurement_types);
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .unwrap();

    for (rov, _) in matched_obs {
        if let Some(&r_sq) = sat_to_residual.get(&rov.sat) {
            let (az, el) = calculate_az_el(rov, rov_llh, pos_apc, ephemerides, state_time);
            writeln!(
                file,
                "{},{},{:.2},{:.2},{:.2},{:.4},{:.0},{:.6}",
                epoch_num,
                rov.sat,
                rov.snr,
                el,
                az,
                rov.doppler,
                rov.locktime.unwrap_or(0),
                r_sq
            )
            .unwrap();
        }
    }
}

/// Computes the post-fit residuals and maps them by SatelliteId.
fn extract_post_fit_residuals(
    innovations: &nalgebra::DVector<f64>,
    h_matrix: &nalgebra::DMatrix<f64>,
    dx: &nalgebra::DVector<f64>,
    valid_indices: &[usize],
    measurement_types: &[(gneiss_core::sat::SatelliteId, u8, f64)],
) -> HashMap<gneiss_core::sat::SatelliteId, f64> {
    let post_fit = innovations - h_matrix * dx;
    let mut sat_to_residual = HashMap::new();
    for (i, &idx) in valid_indices.iter().enumerate() {
        let (sat_id, type_code, _) = measurement_types[idx];
        if type_code == 0 {
            // Pseudorange
            sat_to_residual.insert(sat_id, post_fit[i].powi(2));
        }
    }
    sat_to_residual
}

/// Calculates the azimuth and elevation for a given observation.
fn calculate_az_el(
    rov: &crate::filter::DdObservation,
    rov_llh: nalgebra::Vector3<f64>,
    pos_apc: nalgebra::Vector3<f64>,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    state_time: gneiss_core::time::GpsTime,
) -> (f64, f64) {
    if let Some(eph) = ephemerides.iter().find(|e| e.sat() == rov.sat) {
        let (sat_pos, _) =
            crate::engine::measurement_math::get_sat_state(eph, rov.pr_l1, state_time, pos_apc);
        let (a, e) = gneiss_core::coords::az_el(rov_llh, pos_apc, sat_pos);
        (a.to_degrees(), e.to_degrees())
    } else {
        (0.0, 0.0)
    }
}

/// NLL Loss (Heteroscedastic) where GNN predicts `log_sigma_sq` directly
pub fn nll_loss_logvar(log_sigma_sq: &Tensor, r_sq: &Tensor, mask: &Tensor) -> Result<Tensor> {
    let log_sigma_sq = log_sigma_sq.clamp(-MAX_LOG_VAR, MAX_LOG_VAR)?;
    // exp(-log_var) is equivalent to 1 / sigma^2 but strictly stable
    let inv_sigma_sq = log_sigma_sq.neg()?.exp()?;
    let term2 = r_sq.broadcast_mul(&inv_sigma_sq)?;

    let unmasked_loss = log_sigma_sq.broadcast_add(&term2)?;
    let masked_loss = unmasked_loss.broadcast_mul(mask)?;

    let valid_count = mask.sum_all()?;
    masked_loss.sum_all()?.broadcast_div(&valid_count)
}

/// NLL Loss for Softplus output (current GNN architecture)
pub fn nll_loss_softplus(sigma_sq: &Tensor, r_sq: &Tensor, mask: &Tensor) -> Result<Tensor> {
    let eps = Tensor::new(EPSILON, sigma_sq.device())?;
    let sigma_sq_safe = sigma_sq.broadcast_add(&eps)?;

    let log_sigma_sq = sigma_sq_safe.log()?;
    let term2 = r_sq.broadcast_div(&sigma_sq_safe)?;

    let unmasked_loss = log_sigma_sq.broadcast_add(&term2)?;

    let masked_loss = unmasked_loss.broadcast_mul(mask)?;
    let valid_count = mask.sum_all()?;
    let total_loss = masked_loss.sum_all()?;

    total_loss.broadcast_div(&valid_count)
}

/// Normalizes a single raw observation into neural network feature space [0, 1] or [-1, 1].
pub fn normalize_features(
    snr: f32,
    elevation_deg: f32,
    azimuth_deg: f32,
    doppler: f32,
    locktime_ms: f32,
) -> [f32; 5] {
    let norm_snr = snr / 50.0;
    let norm_el = elevation_deg / 90.0;
    let norm_az = azimuth_deg / 360.0;
    let norm_dop = (doppler / 5000.0).clamp(-1.0, 1.0);
    let norm_lock = (locktime_ms / 1000.0).clamp(0.0, 1.0);

    [norm_snr, norm_el, norm_az, norm_dop, norm_lock]
}

#[cfg(test)]
mod tests {
    use super::*;
    use candle_core::Device;

    #[test]
    fn test_normalize_features() {
        // Normal bounds
        let f1 = normalize_features(45.0, 45.0, 180.0, 1000.0, 500.0);
        assert!((f1[0] - 0.9).abs() < 1e-5);
        assert!((f1[1] - 0.5).abs() < 1e-5);
        assert!((f1[2] - 0.5).abs() < 1e-5);
        assert!((f1[3] - 0.2).abs() < 1e-5);
        assert!((f1[4] - 0.5).abs() < 1e-5);

        // Clamping bounds
        let f2 = normalize_features(0.0, 0.0, 0.0, 6000.0, 2000.0);
        assert_eq!(f2[3], 1.0); // Doppler clamped to 1.0
        assert_eq!(f2[4], 1.0); // Locktime clamped to 1.0

        let f3 = normalize_features(0.0, 0.0, 0.0, -10000.0, -500.0);
        assert_eq!(f3[3], -1.0); // Doppler clamped to -1.0
        assert_eq!(f3[4], 0.0); // Locktime clamped to 0.0
    }

    #[test]
    fn test_nll_loss_logvar_stability() -> Result<()> {
        let device = Device::Cpu;
        // Test with extreme values
        // log_sigma_sq contains very small, very large, and normal values
        let log_sigma_sq = Tensor::new(&[[-100.0f32, 100.0, 0.0, -10.0, 10.0]], &device)?;
        // r_sq contains very large residuals and 0
        let r_sq = Tensor::new(&[[1e6f32, 1e6, 1.0, 0.0, 1e6]], &device)?;
        // All valid
        let mask = Tensor::new(&[[1.0f32, 1.0, 1.0, 1.0, 1.0]], &device)?;

        let loss = nll_loss_logvar(&log_sigma_sq, &r_sq, &mask)?;
        let loss_val = loss.to_vec0::<f32>()?;

        // Loss should not be NaN or Infinity
        assert!(!loss_val.is_nan(), "Loss is NaN");
        assert!(!loss_val.is_infinite(), "Loss is infinite");

        Ok(())
    }
}
