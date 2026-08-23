//! Integer Ambiguity Resolution (LAMBDA + FFRT + PAR) for RTK.

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use crate::ambiguity::{ffrt, lambda};
use super::state::{DoubleDiffKey, RtkState};

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
pub fn resolve_ambiguities(state: &RtkState, min_ambiguities: usize, target_pf: f64) -> ArResult {
    let (a_float, q_amb) = state.extract_amb_block();
    let n_amb = a_float.len();
    let float_pos = state.pos_ecef;
    let float_cov = state.extract_pos_cov();

    if n_amb < min_ambiguities || n_amb < 3 {
        return build_float_result(float_pos, float_cov, 0.0, n_amb);
    }

    // 1. Full Ambiguity Resolution (FAR)
    if let Ok(l_res) = lambda::resolve_lambda(&a_float, &q_amb) {
        let thresh = ffrt::calculate_threshold(n_amb, target_pf);
        if l_res.ratio >= thresh {
            let full_idx: Vec<usize> = (0..n_amb).collect();
            let fixed = fixed_subset_ambiguities(state, &full_idx, &l_res.best_integers);
            if let Some((pos, cov)) = project_subset_fixed(state, &a_float, &l_res.best_integers, &q_amb, &full_idx) {
                return ArResult { position_ecef: pos, cov_position: cov, ratio: l_res.ratio, is_fixed: true, num_ambiguities: n_amb, fixed_ambiguities: fixed };
            }
        }
    }

    // 2. Partial Ambiguity Resolution (PAR)
    if let Some(par_res) = try_partial_ar(state, &a_float, &q_amb, min_ambiguities, target_pf) {
        return par_res;
    }

    build_float_result(float_pos, float_cov, 0.0, n_amb)
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

fn try_partial_ar(
    state: &RtkState,
    a_float: &DVector<f64>,
    q_amb: &DMatrix<f64>,
    min_ambs: usize,
    target_pf: f64,
) -> Option<ArResult> {
    let n_amb = a_float.len();
    let float_trace = state.cov[(0, 0)] + state.cov[(1, 1)] + state.cov[(2, 2)];
    if n_amb <= 4 || float_trace > 2.0 {
        return None;
    }

    // Sort ambiguity indices by diagonal variance ascending
    let mut sorted_indices: Vec<usize> = (0..n_amb).collect();
    sorted_indices.sort_by(|&i, &j| q_amb[(i, i)].total_cmp(&q_amb[(j, j)]));

    for k in (min_ambs.max(4)..n_amb).rev() {
        let subset_idx = &sorted_indices[0..k];
        let (sub_a, sub_q) = extract_subset(a_float, q_amb, subset_idx);
        if let Ok(l_res) = lambda::resolve_lambda(&sub_a, &sub_q) {
            let thresh = ffrt::calculate_threshold(k, target_pf);
            if l_res.ratio >= thresh {
                if let Some((pos, cov)) = project_subset_fixed(state, &sub_a, &l_res.best_integers, &sub_q, subset_idx) {
                    let fixed = fixed_subset_ambiguities(state, subset_idx, &l_res.best_integers);
                    return Some(ArResult {
                        position_ecef: pos,
                        cov_position: cov,
                        ratio: l_res.ratio,
                        is_fixed: true,
                        num_ambiguities: k,
                        fixed_ambiguities: fixed,
                    });
                }
            }
        }
    }
    None
}

fn extract_subset(
    a_float: &DVector<f64>,
    q_amb: &DMatrix<f64>,
    indices: &[usize],
) -> (DVector<f64>, DMatrix<f64>) {
    let k = indices.len();
    let mut sub_a = DVector::zeros(k);
    let mut sub_q = DMatrix::zeros(k, k);

    for (r, &i) in indices.iter().enumerate() {
        sub_a[r] = a_float[i];
        for (c, &j) in indices.iter().enumerate() {
            sub_q[(r, c)] = q_amb[(i, j)];
        }
    }
    (sub_a, sub_q)
}

fn project_subset_fixed(
    state: &RtkState,
    sub_a_float: &DVector<f64>,
    sub_a_fixed: &DVector<f64>,
    sub_q_amb: &DMatrix<f64>,
    indices: &[usize],
) -> Option<(Vector3<f64>, Matrix3<f64>)> {
    let q_inv = sub_q_amb.clone().try_inverse()?;
    let k = indices.len();

    let mut p_xa = DMatrix::zeros(3, k);
    for r in 0..3 {
        for (c, &idx) in indices.iter().enumerate() {
            p_xa[(r, c)] = state.cov[(r, 6 + idx)];
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

        let res = resolve_ambiguities(&state, 4, 0.001);
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

        let res = resolve_ambiguities(&state, 4, 0.001);
        assert!(res.is_fixed, "PAR should fix the 4 clean ambiguities");
        assert_eq!(res.num_ambiguities, 4);
    }
}
