//! Robust weighting, outlier screening, and scalar residual filters.

use std::collections::HashMap;
use nalgebra::Vector3;
use crate::estimators::rtk_iekf::state::{DoubleDiffKey, RtkState};
use super::system::compute_tropo_dd;
use super::{pcv_corrected_cp, DoubleDiffMeasurement};

/// Robust estimation: inflate measurement variance when its normalized
/// innovation squared exceeds this threshold (~3-sigma point).
pub const ROBUST_INNOVATION_THRESHOLD: f64 = 9.0;

/// Gradient state seed variance (m^2): ~2 mm of horizon slant delay each.
pub const GRAD_INIT_VAR_M2: f64 = 4.0e-6;
/// Gradient random-walk rate (m^2/s).
pub const GRAD_RW_M2_PER_S: f64 = 5.0e-11;
/// Elevation floor (rad) inside the cot(el) gradient mapping.
pub const GRAD_MIN_SIN_EL: f64 = 0.17;

/// Seed variance of the rover ZWD residual state (m^2): ~15 cm zenith.
pub const ZWD_INIT_VAR_M2: f64 = 0.0225;
/// Random-walk variance rate of the rover ZWD residual (m^2/s).
pub const ZWD_RW_M2_PER_S: f64 = 3e-7;
/// Per-epoch saturation bound on the ZWD correction (m).
pub const MAX_ZWD_STEP_M: f64 = 0.05;

/// Phase-innovation cycle-slip gate (cycles).
pub const PHASE_INNOVATION_GATE_CYCLES: f64 = 500.0;

/// Huber-style variance inflation: keep nominal R below the normalized
/// innovation threshold, decay weight smoothly above it.
pub fn robust_inflate(innovation: f64, variance: f64, gate_scale: f64) -> f64 {
    let threshold = ROBUST_INNOVATION_THRESHOLD * gate_scale;
    let nis = innovation * innovation / variance;
    if nis > threshold {
        variance * (nis / threshold)
    } else {
        variance
    }
}

/// Scalar random-walk Kalman update for the rover ZWD residual with step saturation.
pub fn update_zwd_scalar(
    zwd: f64,
    var: f64,
    rw_m2_per_s: f64,
    dt_s: f64,
    pairs: &[(f64, f64, f64)],
) -> (f64, f64) {
    let prior_var = var + rw_m2_per_s * dt_s.abs().max(1e-3);
    let mut denom = 1.0 / prior_var;
    let mut num = 0.0;
    for &(h, y, r) in pairs {
        if r <= 0.0 {
            continue;
        }
        denom += h * h / r;
        num += h * y / r;
    }
    let post_var = 1.0 / denom;
    let delta = (post_var * num).clamp(-MAX_ZWD_STEP_M, MAX_ZWD_STEP_M);
    (zwd + delta, post_var)
}

/// Keys whose DD phase innovation exceeds the slip gate under the current state.
pub fn phase_innovation_outliers(
    state: &RtkState,
    measurements: &[DoubleDiffMeasurement],
    max_cycles: f64,
) -> Vec<DoubleDiffKey> {
    let mut out = Vec::new();
    for m in measurements {
        let Some(cp_obs) = pcv_corrected_cp(m) else { continue };
        let Some(amb_idx) = state.get_amb_idx(&m.key) else { continue };
        let cur_pos = state.pos_ecef;
        let r_sat = (m.sat_pos - cur_pos).norm();
        let r_ref = (m.ref_pos - cur_pos).norm();
        let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
        let geom_dd = (r_sat - r_ref) - base_dd
            + compute_tropo_dd(m.sat_pos, m.ref_pos, m.base_pos, cur_pos);
        let pred = geom_dd / m.lambda + state.to_dvector()[amb_idx];
        if (cp_obs - pred).abs() > max_cycles {
            out.push(m.key);
        }
    }
    out
}

/// Per-pair iono-free residual outlier screen over a fixed ambiguity set.
pub fn if_residual_outliers(
    pos: Vector3<f64>,
    measurements: &[DoubleDiffMeasurement],
    fixed_n1: &HashMap<DoubleDiffKey, f64>,
    fixed_n2: &HashMap<DoubleDiffKey, f64>,
) -> Vec<DoubleDiffKey> {
    const MAX_DEV_M: f64 = 0.05;

    let mut b2: HashMap<DoubleDiffKey, &DoubleDiffMeasurement> = HashMap::new();
    for m in measurements {
        if m.key.freq_band == 2 && m.dd_cp_cycles.is_some() {
            b2.insert(DoubleDiffKey { freq_band: 1, ..m.key }, m);
        }
    }

    let mut items: Vec<(DoubleDiffKey, f64)> = Vec::new();
    for m in measurements {
        if m.key.freq_band != 1 {
            continue;
        }
        let Some(cp1) = pcv_corrected_cp(m) else { continue };
        let Some(m2) = b2.get(&m.key) else { continue };
        let Some(cp2) = pcv_corrected_cp(m2) else { continue };
        let k2 = DoubleDiffKey { freq_band: 2, ..m.key };
        let (Some(n1), Some(n2)) = (fixed_n1.get(&m.key), fixed_n2.get(&k2)) else {
            continue;
        };

        let r_sat = (m.sat_pos - pos).norm();
        let r_ref = (m.ref_pos - pos).norm();
        let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
        let geom_m = (r_sat - r_ref) - base_dd;

        let d1 = cp1 * m.lambda - geom_m - n1 * m.lambda;
        let d2 = cp2 * m2.lambda - geom_m - n2 * m2.lambda;

        let res = (d1 / (m.lambda * m.lambda) - d2 / (m2.lambda * m2.lambda))
            / (1.0 / (m.lambda * m.lambda) - 1.0 / (m2.lambda * m2.lambda));
        items.push((m.key, res));
    }
    if items.len() < 3 {
        return Vec::new();
    }
    let mut rs: Vec<f64> = items.iter().map(|(_, v)| *v).collect();
    rs.sort_by(|a, b| a.total_cmp(b));
    let med = rs[rs.len() / 2];
    items.into_iter()
        .filter(|(_, v)| (v - med).abs() > MAX_DEV_M)
        .map(|(k, _)| k)
        .collect()
}
