//! Integer Ambiguity Resolution (LAMBDA + FFRT + PAR) for RTK.

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use crate::ambiguity::{ffrt, lambda};
use super::ar_subsets::{
    extract_subset, generate_omission_subsets, generate_two_omission_subsets,
    partition_constellation_subsets, select_par_candidates,
};
use super::state::{DoubleDiffKey, RtkState};
use super::update;
use super::DoubleDiffMeasurement;

/// Result of ambiguity resolution attempt for an epoch.
#[derive(Debug, Clone)]
pub struct ArResult {
    pub position_ecef: Vector3<f64>,
    pub cov_position: Matrix3<f64>,
    pub ratio: f64,
    pub is_fixed: bool,
    pub num_ambiguities: usize,
    /// Integer values of the fixed ambiguity subset, keyed by DD pair.
    pub fixed_ambiguities: Vec<(DoubleDiffKey, f64)>,
}

/// Attempt integer ambiguity resolution (FAR + PAR) on the current float state.
///
/// Resolution runs jointly over all bands: per-band LAMBDA with the looser
/// FFRT thresholds accepts integer vectors whose elements sit 0.4-2+ cycles
/// off their float values (the float DD ambiguities carry unmodeled bias),
/// which biases the conditional position by decimetres. The joint threshold
/// is conservative but only ever fixes near-integer subsets.
pub fn resolve_ambiguities(
    state: &RtkState,
    min_ambiguities: usize,
    target_pf: f64,
    is_kinematic: bool,
) -> ArResult {
    resolve_ambiguities_screened(state, min_ambiguities, target_pf, is_kinematic, None)
}

/// Attempt integer ambiguity resolution with pre-acceptance carrier residual screening.
pub fn resolve_ambiguities_screened(
    state: &RtkState,
    min_ambiguities: usize,
    target_pf: f64,
    is_kinematic: bool,
    dd_meas: Option<&[DoubleDiffMeasurement]>,
) -> ArResult {
    let (a_float, q_amb) = state.extract_amb_block();
    let n_amb = a_float.len();
    let (pos, cov) = (state.pos_ecef, state.extract_pos_cov());
    let float_trace = cov[(0, 0)] + cov[(1, 1)] + cov[(2, 2)];
    let max_float_trace = if is_kinematic { 3.50 } else { 25.0 };
    if n_amb < min_ambiguities || n_amb < 3 || float_trace > max_float_trace {
        return build_float_result(pos, cov, 0.0, n_amb);
    }
    if let Some(far_res) = try_full_ar(state, &a_float, &q_amb, target_pf, is_kinematic, dd_meas) {
        return far_res;
    }
    if let Some(par_res) = try_partial_ar(state, &a_float, &q_amb, min_ambiguities, target_pf, is_kinematic, dd_meas) {
        return par_res;
    }
    build_float_result(pos, cov, 0.0, n_amb)
}

fn try_full_ar(
    state: &RtkState,
    a_float: &DVector<f64>,
    q_amb: &DMatrix<f64>,
    target_pf: f64,
    is_kinematic: bool,
    dd_meas: Option<&[DoubleDiffMeasurement]>,
) -> Option<ArResult> {
    let n_amb = a_float.len();
    if is_kinematic && (0..n_amb).any(|i| q_amb[(i, i)] > 1.0) {
        return None;
    }
    let l_res = lambda::resolve_lambda(a_float, q_amb).ok()?;
    let thresh = ffrt::calculate_threshold(n_amb, target_pf);
    let min_ratio = if is_kinematic {
        thresh.max(2.0)
    } else {
        thresh
    };
    if l_res.ratio < min_ratio {
        return None;
    }
    let full_idx: Vec<usize> = (0..n_amb).collect();
    let (pos, cov) = project_subset_fixed(state, a_float, &l_res.best_integers, q_amb, &full_idx)?;
    let fixed = fixed_subset_ambiguities(state, &full_idx, &l_res.best_integers);
    if let Some(dd) = dd_meas {
        if !update::validate_fixed_carrier_residuals(pos, dd, &fixed, 0.05) {
            return None;
        }
    }
    Some(ArResult {
        position_ecef: pos,
        cov_position: cov,
        ratio: l_res.ratio,
        is_fixed: true,
        num_ambiguities: n_amb,
        fixed_ambiguities: fixed,
    })
}

pub(crate) fn float_result(state: &RtkState) -> ArResult {
    build_float_result(state.pos_ecef, state.extract_pos_cov(), 0.0, state.ambiguities.len())
}

fn build_float_result(pos: Vector3<f64>, cov: Matrix3<f64>, ratio: f64, n_amb: usize) -> ArResult {
    ArResult {
        position_ecef: pos,
        cov_position: cov,
        ratio,
        is_fixed: false,
        num_ambiguities: n_amb,
        fixed_ambiguities: Vec::new(),
    }
}

fn fixed_subset_ambiguities(state: &RtkState, indices: &[usize], integers: &DVector<f64>) -> Vec<(DoubleDiffKey, f64)> {
    indices.iter().enumerate()
        .map(|(i, &idx)| (state.ambiguities[idx].0, integers[i]))
        .collect()
}

struct ParContext<'a> {
    state: &'a RtkState,
    a_float: &'a DVector<f64>,
    q_amb: &'a DMatrix<f64>,
    target_pf: f64,
    is_kinematic: bool,
    dd_meas: Option<&'a [DoubleDiffMeasurement]>,
}

fn try_partial_ar(
    state: &RtkState,
    a_float: &DVector<f64>,
    q_amb: &DMatrix<f64>,
    min_ambs: usize,
    target_pf: f64,
    is_kinematic: bool,
    dd_meas: Option<&[DoubleDiffMeasurement]>,
) -> Option<ArResult> {
    let n_amb = a_float.len();
    let float_trace = state.cov[(0, 0)] + state.cov[(1, 1)] + state.cov[(2, 2)];
    let max_float_trace = if is_kinematic { 3.50 } else { 25.0 };
    if n_amb <= 4 || float_trace > max_float_trace {
        return None;
    }
    let (sorted_indices, max_k) = select_par_candidates(state, a_float, q_amb, min_ambs, is_kinematic);
    let min_k = if is_kinematic { min_ambs.max(6) } else { min_ambs.max(4) };
    if max_k < min_k {
        return None;
    }
    let ctx = ParContext { state, a_float, q_amb, target_pf, is_kinematic, dd_meas };
    if let Some(res) = eval_par_prefix_subsets(&ctx, &sorted_indices, min_k, max_k) {
        return Some(res);
    }
    if is_kinematic {
        if let Some(res) = eval_par_partition_subsets(&ctx, &sorted_indices, min_k) {
            return Some(res);
        }
    }
    let pool_len = max_k.min(min_k + 4);
    if pool_len > min_k {
        eval_par_omission_subsets(&ctx, &sorted_indices, pool_len, min_k)
    } else {
        None
    }
}

fn eval_par_prefix_subsets(
    ctx: &ParContext<'_>,
    sorted_indices: &[usize],
    min_k: usize,
    max_k: usize,
) -> Option<ArResult> {
    for k in (min_k..=max_k).rev() {
        if let Some(res) = eval_par_subset(ctx, &sorted_indices[..k]) {
            return Some(res);
        }
    }
    None
}

fn eval_par_partition_subsets(
    ctx: &ParContext<'_>,
    sorted_indices: &[usize],
    min_k: usize,
) -> Option<ArResult> {
    let clusters = partition_constellation_subsets(ctx.state, sorted_indices, min_k);
    for sub in clusters {
        if let Some(res) = eval_par_subset(ctx, sub.as_slice()) {
            return Some(res);
        }
    }
    None
}

fn eval_par_omission_subsets(
    ctx: &ParContext<'_>,
    sorted_indices: &[usize],
    pool_len: usize,
    min_k: usize,
) -> Option<ArResult> {
    let subsets1 = generate_omission_subsets(sorted_indices, pool_len);
    for sub in subsets1 {
        if sub.len >= min_k {
            if let Some(res) = eval_par_subset(ctx, sub.as_slice()) {
                return Some(res);
            }
        }
    }
    if ctx.is_kinematic && pool_len >= min_k + 2 {
        let subsets2 = generate_two_omission_subsets(sorted_indices, pool_len, min_k);
        for sub in subsets2 {
            if let Some(res) = eval_par_subset(ctx, sub.as_slice()) {
                return Some(res);
            }
        }
    }
    None
}

fn eval_par_subset(
    ctx: &ParContext<'_>,
    subset_idx: &[usize],
) -> Option<ArResult> {
    let k = subset_idx.len();
    let (sub_a, sub_q) = extract_subset(ctx.a_float, ctx.q_amb, subset_idx);
    let l_res = lambda::resolve_lambda(&sub_a, &sub_q).ok()?;
    let thresh = ffrt::calculate_threshold(k, ctx.target_pf);
    let min_ratio = if ctx.is_kinematic {
        thresh.max(2.0)
    } else {
        thresh
    };
    if l_res.ratio < min_ratio {
        return None;
    }
    let (pos, cov) = project_subset_fixed(ctx.state, &sub_a, &l_res.best_integers, &sub_q, subset_idx)?;
    let fixed = fixed_subset_ambiguities(ctx.state, subset_idx, &l_res.best_integers);
    if let Some(dd) = ctx.dd_meas {
        if !update::validate_fixed_carrier_residuals(pos, dd, &fixed, 0.05) {
            return None;
        }
    }
    Some(ArResult {
        position_ecef: pos,
        cov_position: cov,
        ratio: l_res.ratio,
        is_fixed: true,
        num_ambiguities: k,
        fixed_ambiguities: fixed,
    })
}

pub(crate) fn project_subset_fixed(
    state: &RtkState,
    sub_a_float: &DVector<f64>,
    sub_a_fixed: &DVector<f64>,
    sub_q_amb: &DMatrix<f64>,
    indices: &[usize],
) -> Option<(Vector3<f64>, Matrix3<f64>)> {
    let q_inv = sub_q_amb.clone().try_inverse()?;
    let k = indices.len();

    let off = state.amb_offset();
    let mut p_xa = DMatrix::zeros(3, k);
    for r in 0..3 {
        for (c, &idx) in indices.iter().enumerate() {
            p_xa[(r, c)] = state.cov[(r, off + idx)];
        }
    }

    let da = sub_a_float - sub_a_fixed;
    let dx = &p_xa * &q_inv * &da;
    let dx_norm = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
    let float_3d_std = (state.cov[(0, 0)] + state.cov[(1, 1)] + state.cov[(2, 2)]).sqrt();

    // Reject fix if jump exceeds 3-sigma or maximum physical bound
    if dx_norm > (3.0 * float_3d_std).max(0.50) || dx_norm > 2.0 {
        return None;
    }

    let fix_pos = state.pos_ecef - Vector3::new(dx[0], dx[1], dx[2]);

    let cov_reduction = &p_xa * &q_inv * &p_xa.transpose();
    let mut p_xx = state.extract_pos_cov() - cov_reduction;
    p_xx = 0.5 * (p_xx + p_xx.transpose());
    for i in 0..3 {
        p_xx[(i, i)] = p_xx[(i, i)].max(1e-6);
    }

    Some((fix_pos, p_xx))
}

fn compute_integer_conditioning(
    state: &RtkState,
    indices: &[usize],
    a_fix: &[f64],
) -> Option<(DVector<f64>, DMatrix<f64>)> {
    let k = indices.len();
    let mut q_aa = DMatrix::zeros(k, k);
    for (r, &i) in indices.iter().enumerate() {
        for (c, &j) in indices.iter().enumerate() {
            q_aa[(r, c)] = state.cov[(i, j)];
        }
    }
    let q_inv = q_aa.try_inverse()?;
    let dim = state.dim();
    let mut p_all_a = DMatrix::zeros(dim, k);
    for r in 0..dim {
        for (c, &idx) in indices.iter().enumerate() {
            p_all_a[(r, c)] = state.cov[(r, idx)];
        }
    }
    let x_vec = state.to_dvector();
    let da = DVector::from_iterator(k, indices.iter().enumerate().map(|(i, &idx)| x_vec[idx] - a_fix[i]));
    let dx = &p_all_a * &q_inv * &da;
    let cov_red = &p_all_a * &q_inv * &p_all_a.transpose();
    Some((dx, cov_red))
}

/// Condition state and covariance on fixed integer ambiguities (fix-and-hold).
pub fn condition_state_on_integers(
    state: &mut RtkState,
    fixed_ambiguities: &[(DoubleDiffKey, f64)],
) -> bool {
    let mut indices = Vec::with_capacity(fixed_ambiguities.len());
    let mut a_fix = Vec::with_capacity(fixed_ambiguities.len());
    for (key, val) in fixed_ambiguities {
        if let Some(idx) = state.get_amb_idx(key) {
            indices.push(idx);
            a_fix.push(*val);
        }
    }
    if indices.is_empty() { return false; }
    let Some((dx, cov_red)) = compute_integer_conditioning(state, &indices, &a_fix) else { return false; };
    let dx_pos = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
    let float_pos_std = (state.cov[(0, 0)] + state.cov[(1, 1)] + state.cov[(2, 2)]).sqrt();
    if dx_pos > (3.0 * float_pos_std).max(0.50) || dx_pos > 2.0 { return false; }

    let mut x_vec = state.to_dvector() - dx;
    for (i, &idx) in indices.iter().enumerate() { x_vec[idx] = a_fix[i]; }
    state.update_from_dvector(&x_vec);

    let mut new_cov = &state.cov - cov_red;
    new_cov = 0.5 * (&new_cov + &new_cov.transpose());
    for i in 0..state.dim() { new_cov[(i, i)] = new_cov[(i, i)].max(1e-6); }
    for &idx in &indices { new_cov[(idx, idx)] = 1e-4; }
    state.cov = new_cov;
    true
}

/// Check candidate integers consistency across consecutive epochs.
/// Returns true if candidate integers match previous epoch, incrementing `consecutive`.
/// Resets to 1 if integers changed or first epoch.
pub fn update_fix_hysteresis(
    consecutive: &mut u32,
    last_ambs: &mut Vec<(DoubleDiffKey, f64)>,
    current_ambs: &[(DoubleDiffKey, f64)],
) {
    let matched = !last_ambs.is_empty()
        && current_ambs.iter().all(|(k, v)| {
            last_ambs.iter().find(|(pk, _)| pk == k).is_none_or(|(_, pv)| (v - pv).abs() < 1e-3)
        });
    *consecutive = if matched { *consecutive + 1 } else { 1 };
    *last_ambs = current_ambs.to_vec();
}

/// Reset the fix hysteresis counter and previous ambiguities.
pub fn reset_fix_hysteresis(
    consecutive: &mut u32,
    last_ambs: &mut Vec<(DoubleDiffKey, f64)>,
) {
    *consecutive = 0;
    last_ambs.clear();
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::time::GpsTime;
    use crate::estimators::rtk_iekf::state::DoubleDiffKey;

    #[test]
    fn test_ar_insufficient_ambiguities_returns_float() {
        let mut state = RtkState::new(Vector3::new(1.0, 2.0, 3.0), GpsTime::new(2000, 100.0));
        let k1 = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
        state.ensure_ambiguity(k1, 10.2, 0.05);

        let res = resolve_ambiguities(&state, 4, 0.001, false);
        assert!(!res.is_fixed);
        assert_eq!(res.num_ambiguities, 1);
    }

    #[test]
    fn test_partial_ambiguity_resolution_fixes_clean_subset() {
        let mut state = RtkState::new(Vector3::new(100.0, 200.0, 300.0), GpsTime::new(2000, 100.0));
        state.cov[(0, 0)] = 0.04;
        state.cov[(1, 1)] = 0.04;
        state.cov[(2, 2)] = 0.04;
        // Add 4 clean ambiguities with tiny variance and 1 dirty ambiguity with huge variance
        for i in 2..=5 {
            let key = DoubleDiffKey { constellation_id: 0, sat: i, ref_sat: 1, freq_band: 1 };
            state.ensure_ambiguity(key, (i * 10) as f64 + 0.001, 0.0005);
        }
        // Dirty ambiguity on sat 6
        let k_dirty = DoubleDiffKey { constellation_id: 0, sat: 6, ref_sat: 1, freq_band: 1 };
        state.ensure_ambiguity(k_dirty, 45.5, 100.0);

        let res = resolve_ambiguities(&state, 4, 0.001, false);
        assert!(res.is_fixed, "PAR should fix the 4 clean ambiguities");
        assert_eq!(res.num_ambiguities, 4);
    }

    #[test]
    fn test_condition_state_on_integers_updates_variance_and_position() {
        let mut state = RtkState::new(Vector3::new(100.0, 200.0, 300.0), GpsTime::new(2000, 100.0));
        let k1 = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
        state.ensure_ambiguity(k1, 10.02, 0.04);
        let idx = state.get_amb_idx(&k1).unwrap();
        state.cov[(0, idx)] = 0.01;
        state.cov[(idx, 0)] = 0.01;

        let ok = condition_state_on_integers(&mut state, &[(k1, 10.0)]);
        assert!(ok);
        assert!((state.ambiguities[0].1 - 10.0).abs() < 1e-9);
        assert_eq!(state.cov[(idx, idx)], 1e-4);
    }

    #[test]
    fn test_fix_hysteresis_tracking() {
        let mut consecutive = 0;
        let mut last_ambs = Vec::new();
        let k1 = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
        let k2 = DoubleDiffKey { constellation_id: 0, sat: 3, ref_sat: 1, freq_band: 1 };

        // Epoch 1: first fix
        update_fix_hysteresis(&mut consecutive, &mut last_ambs, &[(k1, 10.0), (k2, -5.0)]);
        assert_eq!(consecutive, 1);

        // Epoch 2: matching fix
        update_fix_hysteresis(&mut consecutive, &mut last_ambs, &[(k1, 10.0), (k2, -5.0)]);
        assert_eq!(consecutive, 2);

        // Epoch 3: matching fix
        update_fix_hysteresis(&mut consecutive, &mut last_ambs, &[(k1, 10.0), (k2, -5.0)]);
        assert_eq!(consecutive, 3);

        // Epoch 4: integer jump on k2 -> resets to 1
        update_fix_hysteresis(&mut consecutive, &mut last_ambs, &[(k1, 10.0), (k2, -4.0)]);
        assert_eq!(consecutive, 1);

        // Reset on float
        reset_fix_hysteresis(&mut consecutive, &mut last_ambs);
        assert_eq!(consecutive, 0);
        assert!(last_ambs.is_empty());
    }
}
