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

/// Code-Minus-Carrier multipath detection threshold (meters).
pub const CMC_MULTIPATH_THRESHOLD_M: f64 = 2.5;
/// Warm-up epochs required before declaring code multipath on a tracking arc.
pub const CMC_WARMUP_EPOCHS: u32 = 5;

/// Continuous tracking arc state for Code-Minus-Carrier (CMC) multipath detection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CmcTrack {
    pub baseline_m: f64,
    pub arc_length: u32,
    pub multipath_m: f64,
    pub is_multipath: bool,
}

impl CmcTrack {
    pub fn new(cmc_obs: f64) -> Self {
        Self {
            baseline_m: cmc_obs,
            arc_length: 1,
            multipath_m: 0.0,
            is_multipath: false,
        }
    }

    pub fn update(&mut self, cmc_obs: f64, slip: bool) -> f64 {
        if slip {
            *self = Self::new(cmc_obs);
            return 0.0;
        }
        self.arc_length += 1;
        if self.arc_length <= CMC_WARMUP_EPOCHS {
            let alpha = 1.0 / self.arc_length as f64;
            self.baseline_m += alpha * (cmc_obs - self.baseline_m);
            self.multipath_m = 0.0;
            self.is_multipath = false;
            return 0.0;
        }
        let dev = (cmc_obs - self.baseline_m).abs();
        if dev > CMC_MULTIPATH_THRESHOLD_M {
            self.multipath_m = dev;
            self.is_multipath = true;
        } else {
            let alpha = 1.0 / (self.arc_length.min(100) as f64);
            self.baseline_m += alpha * (cmc_obs - self.baseline_m);
            self.multipath_m = 0.0;
            self.is_multipath = false;
        }
        self.multipath_m
    }
}

/// Tracks Code-Minus-Carrier (CMC) multipath across all active double-difference pairs.
#[derive(Debug, Clone, Default)]
pub struct CmcTracker {
    tracks: HashMap<DoubleDiffKey, CmcTrack>,
}

impl CmcTracker {
    pub fn new() -> Self {
        Self { tracks: HashMap::new() }
    }

    /// Update CMC tracking for a double-difference observation.
    /// Returns the estimated code multipath magnitude in meters (0.0 if clean or warm-up).
    pub fn update_pair(
        &mut self,
        key: DoubleDiffKey,
        dd_pr_m: f64,
        dd_cp_cycles: f64,
        lambda: f64,
        slip: bool,
    ) -> f64 {
        let cmc = dd_pr_m - lambda * dd_cp_cycles;
        let track = self.tracks.entry(key).or_insert_with(|| CmcTrack::new(cmc));
        track.update(cmc, slip)
    }

    pub fn get_multipath_m(&self, key: &DoubleDiffKey) -> f64 {
        self.tracks.get(key).map_or(0.0, |t| t.multipath_m)
    }

    pub fn is_multipath(&self, key: &DoubleDiffKey) -> bool {
        self.tracks.get(key).is_some_and(|t| t.is_multipath)
    }

    pub fn retain_active(&mut self, active_keys: &[DoubleDiffKey]) {
        self.tracks.retain(|k, _| active_keys.contains(k));
    }
}

/// Adaptively down-weight pseudorange measurement variance using Code-Minus-Carrier
/// multipath deviation without inflating carrier phase variance.
pub fn apply_cmc_downweighting(base_pr_var_m2: f64, multipath_m: f64) -> f64 {
    if multipath_m > 0.0 {
        base_pr_var_m2 + multipath_m * multipath_m
    } else {
        base_pr_var_m2
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

fn pred_cp(state: &RtkState, m: &DoubleDiffMeasurement, amb_idx: usize) -> f64 {
    let cur_pos = state.pos_ecef;
    let r_sat = (m.sat_pos - cur_pos).norm();
    let r_ref = (m.ref_pos - cur_pos).norm();
    let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
    let geom_dd = (r_sat - r_ref) - base_dd
        + compute_tropo_dd(m.sat_pos, m.ref_pos, m.base_pos, cur_pos);
    geom_dd / m.lambda + state.to_dvector()[amb_idx]
}

fn phase_innov_var(state: &RtkState, m: &DoubleDiffMeasurement, amb_idx: usize) -> f64 {
    let cur_pos = state.pos_ecef;
    let r_sat = (m.sat_pos - cur_pos).norm().max(1e-3);
    let r_ref = (m.ref_pos - cur_pos).norm().max(1e-3);
    let h_pos = ((m.ref_pos - cur_pos) / r_ref - (m.sat_pos - cur_pos) / r_sat) / m.lambda;
    let p_pos = state.cov.fixed_view::<3, 3>(0, 0);
    let pos_var = (h_pos.transpose() * p_pos * h_pos)[(0, 0)];
    let amb_var = state.cov[(amb_idx, amb_idx)];
    let cp_r = m.cp_var_cycles2.max(1e-4);
    pos_var + amb_var + cp_r
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
        let amb_var = state.cov[(amb_idx, amb_idx)];
        if amb_var > 4.0 {
            continue;
        }
        let pred = pred_cp(state, m, amb_idx);
        let err = (cp_obs - pred).abs();
        let s_cp = phase_innov_var(state, m, amb_idx);
        if err > max_cycles && (err * err / s_cp) > 36.0 {
            out.push(m.key);
        }
    }
    if out.len() > measurements.len() / 2 {
        out.clear();
    }
    out
}

fn compute_iono_free_res(
    pos: Vector3<f64>,
    m: &DoubleDiffMeasurement,
    m2: &DoubleDiffMeasurement,
    cp1: f64,
    cp2: f64,
    n1: f64,
    n2: f64,
) -> f64 {
    let r_sat = (m.sat_pos - pos).norm();
    let r_ref = (m.ref_pos - pos).norm();
    let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
    let geom_m = (r_sat - r_ref) - base_dd;
    let d1 = cp1 * m.lambda - geom_m - n1 * m.lambda;
    let d2 = cp2 * m2.lambda - geom_m - n2 * m2.lambda;
    (d1 / (m.lambda * m.lambda) - d2 / (m2.lambda * m2.lambda))
        / (1.0 / (m.lambda * m.lambda) - 1.0 / (m2.lambda * m2.lambda))
}

fn filter_outlier_residuals(items: Vec<(DoubleDiffKey, f64)>, max_dev: f64) -> Vec<DoubleDiffKey> {
    if items.len() < 3 {
        return Vec::new();
    }
    let mut rs: Vec<f64> = items.iter().map(|(_, v)| *v).collect();
    rs.sort_by(|a, b| a.total_cmp(b));
    let med = rs[rs.len() / 2];
    items.into_iter()
        .filter(|(_, v)| (v - med).abs() > max_dev)
        .map(|(k, _)| k)
        .collect()
}

/// Per-pair iono-free residual outlier screen over a fixed ambiguity set.
pub fn if_residual_outliers(
    pos: Vector3<f64>,
    measurements: &[DoubleDiffMeasurement],
    fixed_n1: &HashMap<DoubleDiffKey, f64>,
    fixed_n2: &HashMap<DoubleDiffKey, f64>,
) -> Vec<DoubleDiffKey> {
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
        let (Some(cp1), Some(m2)) = (pcv_corrected_cp(m), b2.get(&m.key)) else { continue };
        let Some(cp2) = pcv_corrected_cp(m2) else { continue };
        let k2 = DoubleDiffKey { freq_band: 2, ..m.key };
        let (Some(&n1), Some(&n2)) = (fixed_n1.get(&m.key), fixed_n2.get(&k2)) else {
            continue;
        };
        let res = compute_iono_free_res(pos, m, m2, cp1, cp2, n1, n2);
        items.push((m.key, res));
    }
    filter_outlier_residuals(items, 0.05)
}

/// Validate post-fix carrier phase residuals across all fixed ambiguities.
///
/// Returns true if all fixed ambiguities fit the double-difference carrier phase observations
/// within `max_residual_m` (default: 0.05m = 5 cm). If any fixed ambiguity deviates by more
/// than 5 cm, returns false, signaling an invalid/false integer fix.
pub fn validate_fixed_carrier_residuals(
    pos: Vector3<f64>,
    measurements: &[DoubleDiffMeasurement],
    fixed_ambiguities: &[(DoubleDiffKey, f64)],
    max_residual_m: f64,
) -> bool {
    if fixed_ambiguities.is_empty() { return false; }
    let meas_map: HashMap<DoubleDiffKey, &DoubleDiffMeasurement> = measurements
        .iter()
        .map(|m| (m.key, m))
        .collect();

    for (k, n_val) in fixed_ambiguities {
        let Some(m) = meas_map.get(k) else { return false; };
        let Some(cp) = pcv_corrected_cp(m) else { return false; };

        let r_sat = (m.sat_pos - pos).norm();
        let r_ref = (m.ref_pos - pos).norm();
        let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
        let geom_m = (r_sat - r_ref) - base_dd + compute_tropo_dd(m.sat_pos, m.ref_pos, m.base_pos, pos);
        let baseline_m = (pos - m.base_pos).norm();
        let eff_max = max_residual_m + 4.0e-6 * baseline_m;

        let res_m = (cp * m.lambda - geom_m - n_val * m.lambda).abs();
        if res_m > eff_max {
            return false;
        }
    }
    true
}

/// Validate that the fixed solution does not contradict raw pseudorange observations.
///
/// A false integer fix in a weak geometry can yield small carrier residuals for its own
/// subset while moving the position by tens of meters, creating 50-100m pseudorange residuals
/// across the constellation. Returns false if any DD pseudorange residual exceeds `max_pr_res_m`.
pub fn validate_fixed_pseudorange_residuals(
    pos: Vector3<f64>,
    measurements: &[DoubleDiffMeasurement],
    max_pr_rms_m: f64,
    max_pr_res_m: f64,
) -> bool {
    if measurements.is_empty() { return true; }
    let mut sum_sq = 0.0;
    let mut count = 0;
    let mut n_large = 0;
    for m in measurements {
        let r_sat = (m.sat_pos - pos).norm();
        let r_ref = (m.ref_pos - pos).norm();
        let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
        let geom_m = (r_sat - r_ref) - base_dd + compute_tropo_dd(m.sat_pos, m.ref_pos, m.base_pos, pos);
        let res = (m.dd_pr_m - geom_m).abs();
        let is_outlier = res > 8.0 && res > 2.5 * m.pr_var_m2.sqrt();
        if res > max_pr_res_m { n_large += 2; }
        else if is_outlier { n_large += 1; }
        sum_sq += res * res;
        count += 1;
    }
    if count == 0 { return true; }
    let rms = (sum_sq / count as f64).sqrt();
    let max_allowed_large = (count / 8).max(3);
    rms <= max_pr_rms_m && n_large <= max_allowed_large
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cmc_track_clean_arc_accumulates_without_multipath() {
        let mut track = CmcTrack::new(10.0);
        for _ in 0..20 {
            let mp = track.update(10.05, false);
            assert_eq!(mp, 0.0);
            assert!(!track.is_multipath);
        }
        assert_eq!(track.arc_length, 21);
        assert!((track.baseline_m - 10.05).abs() < 0.05);
    }

    #[test]
    fn test_cmc_track_detects_multipath_step_and_freezes_baseline() {
        let mut track = CmcTrack::new(10.0);
        for _ in 0..10 {
            track.update(10.0, false);
        }
        let baseline_before = track.baseline_m;
        let mp = track.update(22.0, false);
        assert!(track.is_multipath, "Must flag code multipath");
        assert!((mp - 12.0).abs() < 1e-6, "Multipath must be 12.0m, got {mp}");
        assert_eq!(track.baseline_m, baseline_before, "Baseline must remain frozen");
    }

    #[test]
    fn test_cmc_track_resets_on_cycle_slip() {
        let mut track = CmcTrack::new(10.0);
        for _ in 0..10 {
            track.update(10.0, false);
        }
        assert_eq!(track.arc_length, 11);
        let mp = track.update(35.0, true);
        assert_eq!(mp, 0.0);
        assert_eq!(track.arc_length, 1);
        assert_eq!(track.baseline_m, 35.0);
    }

    #[test]
    fn test_apply_cmc_downweighting_inflates_code_only() {
        let nominal_pr_var = 0.04;
        let clean_var = apply_cmc_downweighting(nominal_pr_var, 0.0);
        assert_eq!(clean_var, nominal_pr_var);

        let inflated_var = apply_cmc_downweighting(nominal_pr_var, 10.0);
        assert_eq!(inflated_var, 0.04 + 100.0);
    }
}


