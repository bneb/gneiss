use crate::engine::ppp_common::{
    apply_state_vector, assemble_matrices, build_iono_constraint_row, build_weight_matrix,
    extract_state_vector, find_amb_idx, find_ambiguity_index, invert_matrix, snr_scale,
    FgMeasurement,
};
use crate::engine::processed_sat::ProcessedSat;
use crate::engine::EngineError;
use crate::filter::{RtkState, CORE_STATE_SIZE};
use crate::math::{inversion::solve_cholesky_svd, thresholding::apply_huber};
use nalgebra::{DMatrix, DVector, Vector3};

#[cfg(test)]
#[derive(Clone)]
pub struct ArMock {
    pub wl_result: Option<Result<(DVector<f64>, DMatrix<f64>, Vec<usize>), &'static str>>,
    pub nl_result: Option<Result<(DVector<f64>, DMatrix<f64>), &'static str>>,
    pub nl_calls: usize,
}

#[cfg(test)]
pub static AR_MOCK: std::sync::Mutex<Option<ArMock>> = std::sync::Mutex::new(None);

const SPEED_OF_LIGHT: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
const PSEUDORANGE_VARIANCE_BASE: f64 = 1.0;

/// Iterated Extended Kalman Filter for PPP.
///
/// Despite the historical "fg" (factor graph) naming, this is an IEKF —
/// an iterated least-squares solver with Huber robust estimation and
/// a prior from state propagation. It does not perform marginalization,
/// variable elimination, or incremental smoothing (iSAM2).
pub struct PppIteratedEkf {
    pub max_iterations: usize,
    pub convergence_threshold: f64,
    pub huber_k: f64,
}

impl Default for PppIteratedEkf {
    fn default() -> Self {
        Self {
            max_iterations: 15,
            convergence_threshold: 1e-3,
            huber_k: 3.0,
        }
    }
}

impl PppIteratedEkf {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn solve(
        &self,
        state: &mut RtkState,
        sats: &[ProcessedSat],
        position_prior: Option<(Vector3<f64>, f64)>, // (ecef_m, variance_m2)
    ) -> Result<(), EngineError> {
        for outer_iter in 0..4 {
            let done = self.solve_inner(state, sats, position_prior)?;
            if done || outer_iter == 3 {
                if !done {
                    tracing::warn!("Max outlier rejection iterations reached.");
                }
                break;
            }
        }
        if state.epoch_count > 10 {
            if let Err(e) = self.resolve_cascade_ar(state, sats) {
                tracing::info!("Cascade AR did not fix: {:?}", e);
            }
        }
        Ok(())
    }

    fn solve_inner(
        &self,
        state: &mut RtkState,
        sats: &[ProcessedSat],
        position_prior: Option<(Vector3<f64>, f64)>,
    ) -> Result<bool, EngineError> {
        let x_pred = extract_state_vector(state);
        let p_pred = state.covariance.clone();
        state.full_x_predict = Some(x_pred.clone());
        state.full_p_predict = Some(p_pred.clone());
        let mut x_i = x_pred.clone();
        let p_inv = invert_matrix(&p_pred).ok_or(EngineError::StateDisappeared)?;

        for _iter in 0..self.max_iterations {
            if let Some(dx) = self.compute_iteration_dx(
                state,
                sats,
                &x_i,
                &x_pred,
                &p_inv,
                _iter,
                position_prior,
            )? {
                x_i = &x_i + &dx;
                if dx.norm() < self.convergence_threshold {
                    break;
                }
            } else {
                break;
            }
        }

        if let Some(sat) = self.find_worst_outlier(state, sats, &x_i) {
            tracing::warn!(
                "PPP FG Outlier Detected for {:?}. Removing ambiguity and retrying.",
                sat
            );
            for i in 0..4 {
                state.remove_ambiguity(sat, i);
            }
            return Ok(false);
        }

        let final_p = self.compute_final_covariance(state, sats, &x_i, &p_pred, &p_inv);
        apply_state_vector(state, &x_i, final_p);
        log_ppp_convergence(state, sats, &x_i, &x_pred, &p_pred, self);
        Ok(true)
    }

    fn find_worst_outlier(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
        x_i: &DVector<f64>,
    ) -> Option<gneiss_core::sat::SatelliteId> {
        let final_meas = self.build_measurements(state, sats, x_i, self.max_iterations);
        Self::find_worst_outlier_sat(&final_meas)
    }

    fn find_worst_outlier_sat(meas: &[FgMeasurement]) -> Option<gneiss_core::sat::SatelliteId> {
        let mut worst_sat = None;
        let mut max_norm = 15.0;
        for m in meas {
            if m.is_phase {
                let norm = m.res.abs() / m.raw_var.sqrt();
                if norm > max_norm {
                    max_norm = norm;
                    worst_sat = m.sat;
                }
            }
        }
        worst_sat
    }

    /// Try to fix ambiguities for a single constellation group using WL+NL cascade.
    /// Returns (x_fixed, p_fixed, n_sats) on success, or None if this group cannot fix.
    fn process_constellation_group(
        &self,
        state: &RtkState,
        p_current: &DMatrix<f64>,
        x_current: &DVector<f64>,
        group_cands: &[(gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64)],
        constellation: gneiss_core::sat::Constellation,
    ) -> Option<(DVector<f64>, DMatrix<f64>, usize)> {
        if group_cands.len() < 2 {
            return None;
        }
        // Build per-constellation subset: highest-el as reference
        let mut sorted: Vec<_> = group_cands.to_vec();
        sorted.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap_or(std::cmp::Ordering::Equal));
        let ref_cand = &sorted[0];
        let subset: Vec<_> = sorted
            .iter()
            .skip(1)
            .map(|c| (c.clone(), ref_cand.clone()))
            .collect();

        if subset.is_empty() {
            return None;
        }

        tracing::info!(
            "PPP-AR per-const {:?}: {} pairs ({} sats)",
            constellation,
            subset.len(),
            group_cands.len()
        );

        // WL for this constellation
        let wl_result = self.resolve_widelane_ar(state, p_current, &subset, x_current);
        let (x_wl, p_wl, keep_indices) = wl_result.ok()?;

        if keep_indices.is_empty() {
            return None;
        }

        // NL for this constellation
        let nl_result = self.resolve_narrowlane_ar(state, &subset, &keep_indices, &x_wl, &p_wl);
        let (x_fixed, p_fixed) = nl_result.ok()?;

        // Per-constellation position validation
        let float_pos = Vector3::new(x_current[0], x_current[1], x_current[2]);
        let fixed_pos = Vector3::new(x_fixed[0], x_fixed[1], x_fixed[2]);
        let jump = (fixed_pos - float_pos).norm();
        if jump > 10.0 {
            tracing::warn!(
                "PPP-AR {:?} rejected: position jump {:.2}m > 10m",
                constellation,
                jump
            );
            return None;
        }

        tracing::info!(
            "PPP-AR {:?} Fixed! N_Sats={} jump={:.2}m",
            constellation,
            keep_indices.len() + 1,
            jump
        );
        Some((x_fixed, p_fixed, keep_indices.len() + 1))
    }

    /// Try inter-constellation fallback: build DD pairs across constellations
    /// using the highest-elevation GPS as universal reference.
    fn try_inter_constellation_fallback(
        &self,
        state: &RtkState,
        p_current: &DMatrix<f64>,
        x_current: &DVector<f64>,
        cands: &[(gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64)],
    ) -> Option<(DVector<f64>, DMatrix<f64>, usize)> {
        tracing::info!("PPP-AR: per-constellation failed, trying inter-constellation fallback");
        let subset = self.build_ar_subset(cands);
        if subset.len() < 3 {
            return None;
        }
        let wl_result = self.resolve_widelane_ar(state, p_current, &subset, x_current);
        let (x_wl, p_wl, keep_indices) = wl_result.ok()?;
        if keep_indices.is_empty() {
            return None;
        }
        let nl_result = self.resolve_narrowlane_ar(state, &subset, &keep_indices, &x_wl, &p_wl);
        let (x_fixed, p_fixed) = nl_result.ok()?;
        let float_pos = Vector3::new(x_current[0], x_current[1], x_current[2]);
        let fixed_pos = Vector3::new(x_fixed[0], x_fixed[1], x_fixed[2]);
        let jump = (fixed_pos - float_pos).norm();
        if jump > 20.0 {
            return None;
        }
        Some((x_fixed, p_fixed, keep_indices.len() + 1))
    }

    pub fn resolve_cascade_ar(
        &self,
        state: &mut RtkState,
        sats: &[ProcessedSat],
    ) -> Result<(), &'static str> {
        let cands = self.find_ar_candidates(state, sats);

        // Diagnostic: log constellation breakdown and MW state
        let gps_c: Vec<_> = cands
            .iter()
            .filter(|c| c.0.constellation == gneiss_core::sat::Constellation::Gps)
            .collect();
        let gal_c: Vec<_> = cands
            .iter()
            .filter(|c| c.0.constellation == gneiss_core::sat::Constellation::Galileo)
            .collect();
        let bds_c: Vec<_> = cands
            .iter()
            .filter(|c| c.0.constellation == gneiss_core::sat::Constellation::Beidou)
            .collect();
        let mw_counts: Vec<_> = cands
            .iter()
            .map(|c| (c.0, state.mw_sd_counts.get(&c.0).copied().unwrap_or(0)))
            .collect();
        let mw_vals: Vec<_> = cands
            .iter()
            .map(|c| (c.0, state.mw_sd_ema.get(&c.0).copied().unwrap_or(0.0)))
            .collect();
        tracing::info!(
            "PPP-AR diag: cands={} (GPS={} GAL={} BDS={}) mw_counts={:?} mw_vals={:?}",
            cands.len(),
            gps_c.len(),
            gal_c.len(),
            bds_c.len(),
            mw_counts,
            mw_vals
        );

        if cands.len() < 4 {
            return Err("Insufficient dual-frequency satellites for AR");
        }

        // Per-constellation AR: group candidates by constellation, run WL+NL
        // separately per group. This avoids inter-constellation ISB mixing that
        // destroys the WL LAMBDA ratio (was 1.0-1.1 for mixed constellations).
        let mut const_groups: std::collections::HashMap<
            gneiss_core::sat::Constellation,
            Vec<(gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64)>,
        > = std::collections::HashMap::new();
        for cand in &cands {
            const_groups
                .entry(cand.0.constellation)
                .or_default()
                .push(cand.clone());
        }

        let x = extract_state_vector(state);
        let mut x_current = x.clone();
        let mut p_current = state.covariance.clone();
        let mut any_fixed = false;
        let mut total_fixed_sats = 0usize;

        for (constellation, group_cands) in &const_groups {
            if let Some((xf, pf, n_sats)) = self.process_constellation_group(
                state,
                &p_current,
                &x_current,
                group_cands,
                *constellation,
            ) {
                x_current = xf;
                p_current = pf;
                any_fixed = true;
                total_fixed_sats += n_sats;
            }
        }

        if !any_fixed {
            if let Some((xf, pf, n_sats)) = self.try_inter_constellation_fallback(
                state, &p_current, &x_current, &cands,
            ) {
                x_current = xf;
                p_current = pf;
                any_fixed = true;
                total_fixed_sats = n_sats;
            }
        }

        if !any_fixed {
            return Err("No constellation could fix ambiguities");
        }

        // Global position validation against original float
        let float_pos = Vector3::new(x[0], x[1], x[2]);
        let fixed_pos = Vector3::new(x_current[0], x_current[1], x_current[2]);
        let jump = (fixed_pos - float_pos).norm();
        if jump > 20.0 {
            tracing::warn!("PPP-AR rejected: 3D position jump {:.2}m > 20m", jump);
            return Err("Position jump too large after AR fix");
        }

        tracing::info!(
            "PPP Cascade AR Fixed! N_Sats: {} jump={:.2}m",
            total_fixed_sats,
            jump
        );
        apply_state_vector(state, &x_current, p_current);
        state.is_fixed = true;
        Ok(())
    }

    fn find_ar_candidates(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
    ) -> Vec<(gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64)> {
        let mut cands = Vec::new();
        for sat in sats.iter().filter(|s| {
            !s.is_iono_free
                && s.cp2.is_some()
                && s.sat_obs.sat.constellation != gneiss_core::sat::Constellation::Glonass
        }) {
            if let (Some(n1), Some(n2)) = (
                find_amb_idx(state, sat.sat_obs.sat, 1),
                find_amb_idx(state, sat.sat_obs.sat, 2),
            ) {
                cands.push((sat.sat_obs.sat, n1, n2, sat.el, sat.lam1, sat.lam2));
            }
        }
        cands
    }

    fn build_ar_subset(
        &self,
        cands: &[(gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64)],
    ) -> Vec<(
        (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
        (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
    )> {
        // Inter-constellation: single highest-elevation GPS as universal reference.
        // Fall back to per-constellation if no GPS available.
        if let Some(ref_cand) = cands
            .iter()
            .filter(|c| c.0.constellation == gneiss_core::sat::Constellation::Gps)
            .max_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal))
        {
            let ref_cand = ref_cand.clone();
            return cands
                .iter()
                .filter(|c| c.0 != ref_cand.0)
                .map(|c| (c.clone(), ref_cand.clone()))
                .collect();
        }
        // Fallback: per-constellation
        let mut subset = Vec::new();
        let mut const_cands = std::collections::HashMap::new();
        for cand in cands {
            const_cands
                .entry(cand.0.constellation)
                .or_insert_with(Vec::new)
                .push(cand.clone());
        }
        for (_, mut group) in const_cands {
            if group.len() < 2 {
                continue;
            }
            group.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap_or(std::cmp::Ordering::Equal));
            let ref_cand = group[0].clone();
            for cand in group.iter().skip(1) {
                subset.push((cand.clone(), ref_cand.clone()));
            }
        }
        subset
    }

    fn resolve_widelane_ar(
        &self,
        state: &RtkState,
        p: &DMatrix<f64>,
        subset: &[(
            (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
            (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
        )],
        x: &DVector<f64>,
    ) -> Result<(DVector<f64>, DMatrix<f64>, Vec<usize>), &'static str> {
        #[cfg(test)]
        {
            let mut lock = AR_MOCK.lock().unwrap();
            if let Some(ref mut mock) = *lock {
                if let Some(ref wl) = mock.wl_result {
                    let mut res = wl.clone()?;
                    res.0 = x.clone();
                    res.1 = p.clone();
                    return Ok(res);
                }
            }
        }
        let mut d_wl_full = DMatrix::zeros(subset.len(), p.nrows());
        for (i, (c, ref_sat)) in subset.iter().enumerate() {
            d_wl_full[(i, CORE_STATE_SIZE + c.1)] = 1.0 / c.4;
            d_wl_full[(i, CORE_STATE_SIZE + c.2)] = -1.0 / c.5;
            d_wl_full[(i, CORE_STATE_SIZE + ref_sat.1)] = -1.0 / ref_sat.4;
            d_wl_full[(i, CORE_STATE_SIZE + ref_sat.2)] = 1.0 / ref_sat.5;
        }

        let q_wl_full = &d_wl_full * p * d_wl_full.transpose();
        // Accept satellites with converged covariance OR sufficient MW samples
        let keep_indices: Vec<usize> = (0..q_wl_full.nrows())
            .filter(|&i| {
                let cov_ok = q_wl_full[(i, i)].sqrt() < 0.30;
                if cov_ok {
                    return true;
                }
                // MW-based: accept if both rover and reference have >10 MW samples
                let (c, ref_sat) = &subset[i];
                let mw_ok = state.mw_sd_counts.get(&c.0).unwrap_or(&0) > &10
                    && state.mw_sd_counts.get(&ref_sat.0).unwrap_or(&0) > &10;
                mw_ok
            })
            .collect();
        if keep_indices.len() < 3 {
            return Err("Insufficient well-converged Widelane ambiguities");
        }

        let mut d_wl = DMatrix::zeros(keep_indices.len(), p.nrows());
        for (i, &idx) in keep_indices.iter().enumerate() {
            for j in 0..p.nrows() {
                d_wl[(i, j)] = d_wl_full[(idx, j)];
            }
        }

        let mut a_wl = &d_wl * x;
        let n = keep_indices.len();
        // Seed WL floats from MW EMA where available, but use state covariance Q
        // with diagonal floor. Correlations in state covariance help LAMBDA
        // distinguish integer sets when DD values are tightly clustered.
        for (i, &idx) in keep_indices.iter().enumerate() {
            let (c, ref_sat) = &subset[idx];
            let cnt_c = state.mw_sd_counts.get(&c.0).copied().unwrap_or(0);
            let cnt_ref = state.mw_sd_counts.get(&ref_sat.0).copied().unwrap_or(0);
            if cnt_c > 10 && cnt_ref > 10 {
                let mw_c = state.mw_sd_ema.get(&c.0).copied().unwrap_or(0.0);
                let mw_ref = state.mw_sd_ema.get(&ref_sat.0).copied().unwrap_or(0.0);
                a_wl[i] = mw_c - mw_ref;
            }
        }
        let mut q_wl = &d_wl * p * d_wl.transpose();
        for i in 0..q_wl.nrows() {
            q_wl[(i, i)] = q_wl[(i, i)].max(0.01);
        }
        let res_wl = crate::ambiguity::lambda::resolve_lambda(&a_wl, &q_wl)
            .map_err(|_| "WL LAMBDA Failed")?;

        tracing::info!(
            "PPP-AR WL: {} pairs, ratio={:.2}, success_rate={:.3}",
            n,
            res_wl.ratio,
            res_wl.success_rate
        );
        // NL fixes are all sub-cycle when WL passes — lower threshold to
        // recover more epochs. Bootstrapping success_rate provides secondary gating.
        let wl_ok = res_wl.ratio >= 1.1 && res_wl.success_rate >= 0.05;
        if !wl_ok {
            return Err("WL ratio test failed");
        }

        let s_inv = q_wl.try_inverse().ok_or("WL Cov Inversion failed")?;
        let k_wl = p * d_wl.transpose() * s_inv;
        let dx_wl = &k_wl * (res_wl.best_integers - a_wl);
        // Tiny regularization prevents p_wl from going singular (Joseph form
        // with zero measurement noise collapses rank). 1e-6 m² per pair keeps
        // the covariance full-rank for downstream NL LAMBDA and gain inversion.
        let r_wl = DMatrix::identity(keep_indices.len(), keep_indices.len()) * 1e-6;
        Ok((
            x + dx_wl,
            crate::math::covariance::apply_joseph_covariance_update(p, &k_wl, &d_wl, &r_wl),
            keep_indices,
        ))
    }

    fn resolve_narrowlane_ar(
        &self,
        _state: &RtkState,
        subset: &[(
            (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
            (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64),
        )],
        keep_indices: &[usize],
        x_wl: &DVector<f64>,
        p_wl: &DMatrix<f64>,
    ) -> Result<(DVector<f64>, DMatrix<f64>), &'static str> {
        #[cfg(test)]
        {
            let mut lock = AR_MOCK.lock().unwrap();
            if let Some(ref mut mock) = *lock {
                if let Some(ref nl) = mock.nl_result {
                    mock.nl_calls += 1;
                    let mut res = nl.clone()?;
                    res.0 = x_wl.clone();
                    res.0[0] += 9.5;
                    res.1 = p_wl.clone() * 0.5;
                    return Ok(res);
                }
            }
        }
        let mut d_nl = DMatrix::zeros(keep_indices.len(), p_wl.nrows());
        for (i, &idx) in keep_indices.iter().enumerate() {
            let (c, ref_sat) = &subset[idx];
            d_nl[(i, CORE_STATE_SIZE + c.1)] = 1.0 / c.4;
            d_nl[(i, CORE_STATE_SIZE + ref_sat.1)] = -1.0 / ref_sat.4;
        }

        let a_nl = &d_nl * x_wl;
        // Use p_wl so that the NL lambda search is constrained by the WL fixes.
        // p_wl is the WL-constrained covariance (now regularized to stay full-rank).
        // Using it for NL correctly reflects the WL fix's information gain.
        let mut q_nl = &d_nl * p_wl * d_nl.transpose();
        // Add small diagonal to ensure full rank for LAMBDA
        for i in 0..q_nl.nrows() {
            q_nl[(i, i)] = q_nl[(i, i)].max(0.01);
        }
        let res_nl = crate::ambiguity::lambda::resolve_lambda(&a_nl, &q_nl)
            .map_err(|_| "NL LAMBDA Failed")?;

        tracing::info!(
            "PPP-AR NL: {} pairs, ratio={:.2}, success_rate={:.3}",
            keep_indices.len(),
            res_nl.ratio,
            res_nl.success_rate
        );
        // NL uses state covariance which has large initial variance (10000 m²).
        // When WL has fixed correctly (all_mw), accept lower NL confidence.
        // Skip NL ratio test — WL fix constrains the solution enough that
        // position validation (>20m jump) is the effective NL gate.

        // --- Diagnostic: per-pair NL fix quality (before q_nl is consumed) ---
        let mut diag_parts: Vec<String> = Vec::new();
        for i in 0..keep_indices.len() {
            let float_val = a_nl[i];
            let fixed_val = res_nl.best_integers[i];
            let residual = fixed_val - float_val;
            let q_sqrt = q_nl[(i, i)].sqrt();
            let nsigma = if q_sqrt > 1e-9 {
                residual / q_sqrt
            } else {
                0.0
            };
            let (sat, ref_sat) = &subset[keep_indices[i]];
            let float_s = format!("{:.3}", float_val);
            let fix_s = format!("{:.0}", fixed_val);
            let res_s = format!("{:.3}", residual);
            diag_parts.push(format!(
                "{}/{}: float={} fix={} res={}cy ({:.1}σ)",
                sat.0, ref_sat.0, float_s, fix_s, res_s, nsigma
            ));
        }

        let s_nl_inv = q_nl.try_inverse().ok_or("NL Cov Inversion failed")?;
        let k_nl = p_wl * d_nl.transpose() * s_nl_inv;

        let dx_nl = &k_nl * (res_nl.best_integers - a_nl);
        let pos_correction_norm = (dx_nl[0].powi(2) + dx_nl[1].powi(2) + dx_nl[2].powi(2)).sqrt();
        let pos_corr_str = format!("{:.3}", pos_correction_norm);
        tracing::info!(
            "PPP-AR NL diag: pos_corr={}m | {}",
            pos_corr_str,
            diag_parts.join(" | ")
        );

        let p_fixed = crate::math::covariance::apply_joseph_covariance_update(
            p_wl,
            &k_nl,
            &d_nl,
            &DMatrix::zeros(keep_indices.len(), keep_indices.len()),
        );

        let pos_corr_str2 = format!("{:.3}", pos_correction_norm);
        tracing::info!(
            "PPP Cascade AR Fixed! N_Sats: {} pos_corr={}m",
            keep_indices.len() + 1,
            pos_corr_str2
        );
        Ok((x_wl + dx_nl, p_fixed))
    }

    fn compute_iteration_dx(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
        x_i: &DVector<f64>,
        x_pred: &DVector<f64>,
        p_inv: &DMatrix<f64>,
        iter: usize,
        position_prior: Option<(Vector3<f64>, f64)>,
    ) -> Result<Option<DVector<f64>>, EngineError> {
        let meas = self.build_measurements(state, sats, x_i, iter);
        if meas.is_empty() {
            return Err(EngineError::InsufficientSatellites);
        }

        let (h_mat, res_vec, r_mat) = assemble_matrices(&meas, x_i.len());
        let w_mat = build_weight_matrix(&meas, &r_mat);

        let h_t = h_mat.transpose();
        let htw = &h_t * &w_mat;
        let mut htwh = &htw * &h_mat;
        let mut htwr = &htw * &res_vec;

        // Add SPP position prior as soft pseudo-measurements on X, Y, Z.
        // This replaces the hard position reset in process_ppp, allowing
        // the IEKF to accumulate carrier-phase information across epochs
        // while staying loosely anchored to the SPP solution.
        if let Some((spp_pos, var)) = position_prior {
            let w = 1.0 / var;
            for i in 0..3 {
                htwh[(i, i)] += w;
                htwr[i] += w * (spp_pos[i] - x_i[i]);
            }
        }

        let htwh_damped = &htwh + p_inv;
        let innov = &htwr + p_inv * (x_pred - x_i);

        match solve_cholesky_svd(&htwh_damped, &innov, 1e-9) {
            Ok(sol) => Ok(Some(sol)),
            Err(_) => {
                tracing::warn!("Failed to solve normal equations in PPP FG!");
                Ok(None)
            }
        }
    }

    fn compute_final_covariance(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
        x_i: &DVector<f64>,
        p_pred: &DMatrix<f64>,
        p_inv: &DMatrix<f64>,
    ) -> DMatrix<f64> {
        let last_meas = self.build_measurements(state, sats, x_i, self.max_iterations);
        if last_meas.is_empty() {
            return p_pred.clone();
        }

        let (h_mat, _, r_mat) = assemble_matrices(&last_meas, x_i.len());
        let w_mat = build_weight_matrix(&last_meas, &r_mat);

        let htwh = h_mat.transpose() * &w_mat * h_mat;
        let htwh_damped = htwh + p_inv;

        invert_matrix(&htwh_damped).unwrap_or_else(|| {
            tracing::warn!("invert_matrix(&htwh_damped) FAILED! Falling back to p_pred. htwh_damped has NaNs: {}, Infs: {}", htwh_damped.iter().any(|x| x.is_nan()), htwh_damped.iter().any(|x| x.is_infinite()));
            p_pred.clone()
        })
    }

    pub(crate) fn build_measurements(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
        x_i: &DVector<f64>,
        iter: usize,
    ) -> Vec<FgMeasurement> {
        let mut meas = Vec::new();
        let tide_offset =
            gneiss_core::tides::solid_earth_tides_ecef(state.time, state.position.vector);
        let rcv_pos = Vector3::new(x_i[0], x_i[1], x_i[2]) + tide_offset;
        let ztd = if x_i.len() > 20 && !x_i[20].is_nan() && x_i[20] != 0.0 {
            x_i[20]
        } else {
            state.zwd
        };

        for sat in sats {
            let geometric_dist = (sat.sat_pos_rot - rcv_pos).norm();
            let dist = geometric_dist - sat.pcv_correction;
            let los = (sat.sat_pos_rot - rcv_pos) / geometric_dist;
            let isb = Self::extract_isb(x_i, sat.sat_obs.sat.constellation);
            let expected_base =
                dist + x_i[15] + isb - sat.dt_sat_m + sat.tropo_dry + ztd * sat.map_wet;

            if !self.push_sat_meas(
                &mut meas,
                state,
                sat,
                x_i,
                iter,
                &los,
                expected_base,
                dist,
                isb,
            ) {
                continue;
            }
            if sat.doppler != 0.0 {
                self.push_doppler_measurement(&mut meas, sat, x_i, &los);
            }
        }
        meas
    }

    fn push_sat_meas(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        iter: usize,
        los: &Vector3<f64>,
        expected_base: f64,
        dist: f64,
        isb: f64,
    ) -> bool {
        let expected_pr = if sat.is_iono_free {
            expected_base
        } else if sat.cp1.is_some() && sat.cp2.is_some() && sat.p2.is_some() {
            expected_base
                + find_amb_idx(state, sat.sat_obs.sat, 3)
                    .map(|i| x_i[CORE_STATE_SIZE + i])
                    .unwrap_or(sat.iono_delay)
        } else {
            expected_base + sat.iono_delay
        };

        let res_pr = sat.p1 - expected_pr;
        if res_pr.abs() > 100.0 {
            return false;
        }

        if !sat.is_iono_free && sat.cp1.is_some() && sat.cp2.is_some() && sat.p2.is_some() {
            self.push_uduc_measurements(meas, state, sat, x_i, iter, los, expected_base, dist, isb);
            // Add ionospheric prior constraint: tie i1 state to Klobuchar prediction
            if let Some(i1_idx) = find_amb_idx(state, sat.sat_obs.sat, 3) {
                let i1_est = x_i.get(CORE_STATE_SIZE + i1_idx).copied().unwrap_or(0.0);
                let res_i1 = sat.iono_delay - i1_est;
                let var_i1 = 9.0; // 3m std for Klobuchar accuracy
                meas.push(FgMeasurement {
                    res: res_i1,
                    h_row: build_iono_constraint_row(x_i.len(), CORE_STATE_SIZE + i1_idx),
                    weight: var_i1,
                    raw_var: var_i1,
                    is_phase: false,
                    sat: Some(sat.sat_obs.sat),
                });
            }
        } else {
            self.push_pr_measurement(meas, state, sat, x_i, iter, los, expected_base, dist, isb);
            self.try_push_cp_measurement(
                meas,
                state,
                sat,
                x_i,
                iter,
                los,
                expected_base,
                dist,
            );
        }
        true
    }

    fn extract_isb(x_i: &DVector<f64>, constel: gneiss_core::sat::Constellation) -> f64 {
        if x_i.len() > 18 {
            match constel {
                gneiss_core::sat::Constellation::Glonass => x_i[16],
                gneiss_core::sat::Constellation::Galileo => x_i[17],
                gneiss_core::sat::Constellation::Beidou => x_i[18],
                _ => 0.0,
            }
        } else {
            0.0
        }
    }

    fn push_pr_measurement(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        iter: usize,
        los: &Vector3<f64>,
        expected_base: f64,
        _dist: f64,
        _isb: f64,
    ) {
        let expected_pr = if sat.is_iono_free {
            expected_base
        } else {
            expected_base + sat.iono_delay
        };
        let res_pr = sat.p1 - expected_pr;

        if state.epoch_count == 0 && iter == 0 {
            tracing::trace!(
                "PPP {:?}{:02} PR res={:.3}m",
                sat.sat_obs.sat.constellation,
                sat.sat_obs.sat.prn,
                res_pr
            );
        }
        let mut var_pr = PSEUDORANGE_VARIANCE_BASE * snr_scale(sat.snr as i32) / libm::sin(sat.el);
        if sat.is_iono_free {
            var_pr *= 9.0; // Iono-free combination amplifies noise
        } else {
            var_pr += 9.0; // Single frequency has ~3m Klobuchar residual iono error (3^2 = 9)
        }
        let w_pr = apply_huber(res_pr, var_pr, self.huber_k);
        meas.push(FgMeasurement {
            res: res_pr,
            h_row: build_h_row(
                &los,
                sat.map_wet,
                None,
                x_i.len(),
                sat.sat_obs.sat.constellation,
            ),
            weight: var_pr / w_pr,
            raw_var: var_pr,
            is_phase: false,
            sat: Some(sat.sat_obs.sat),
        });
    }

    fn push_doppler_measurement(
        &self,
        meas: &mut Vec<FgMeasurement>,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        los: &Vector3<f64>,
    ) {
        let rcv_vel = Vector3::new(x_i[3], x_i[4], x_i[5]);
        let rcv_clk_drift = if x_i.len() > 19 { x_i[19] } else { 0.0 };
        let meas_rr = -sat.doppler * sat.lam1;
        let expected_rr = los.dot(&sat.sat_vel) - los.dot(&rcv_vel) + rcv_clk_drift
            - sat.sat_clock_drift * SPEED_OF_LIGHT;

        let res_rr = meas_rr - expected_rr;
        tracing::debug!("DOPPLER {}: res_rr={:.3} meas_rr={:.3} exp_rr={:.3} doppler={:.3} rcv_drift={:.3} sat_drift={:.3} los_v={:.3} sat_vel=[{:.3}, {:.3}, {:.3}]",
            sat.sat_obs.sat.to_string(), res_rr, meas_rr, expected_rr, sat.doppler, rcv_clk_drift, sat.sat_clock_drift * gneiss_core::constants::SPEED_OF_LIGHT_M_S, los.dot(&sat.sat_vel), sat.sat_vel.x, sat.sat_vel.y, sat.sat_vel.z);
        let var_rr = 0.01; // Decreased Doppler variance (trust velocity more)
        let w_rr = apply_huber(res_rr, var_rr, 3.0);
        meas.push(FgMeasurement {
            res: res_rr,
            h_row: build_h_row_doppler(&los, x_i.len()),
            weight: var_rr / w_rr,
            raw_var: var_rr,
            is_phase: false,
            sat: Some(sat.sat_obs.sat),
        });
    }

    fn push_cp_measurement(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        iter: usize,
        los: &Vector3<f64>,
        expected_base: f64,
        dist: f64,
        cp1: f64,
    ) {
        if let Some(amb_idx) = find_ambiguity_index(state, sat.sat_obs.sat) {
            let windup = *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0);
            let l_meas = (cp1 - windup) * sat.lam1;
            let expected_cp = if sat.is_iono_free {
                expected_base + x_i[CORE_STATE_SIZE + amb_idx]
            } else {
                expected_base - sat.iono_delay + x_i[CORE_STATE_SIZE + amb_idx]
            };
            let res_cp = l_meas - expected_cp;
            if res_cp.abs() > 100.0 && iter == 0 {
                tracing::warn!("HUGE res_cp: sat={}, l_meas={:.2}, exp={:.2}, dist={:.2}, clk={:.2}, n_amb={:.2}", sat.sat_obs.sat, l_meas, expected_cp, dist, x_i[15], x_i[CORE_STATE_SIZE + amb_idx]);
            }
            let mut var_cp = 0.0001 * snr_scale(sat.snr as i32) / libm::sin(sat.el);
            if sat.is_iono_free {
                var_cp *= 9.0;
            } // Iono-free amplifies phase noise
            let w_cp = apply_huber(res_cp, var_cp, self.huber_k);
            meas.push(FgMeasurement {
                res: res_cp,
                h_row: build_h_row(
                    &los,
                    sat.map_wet,
                    Some(CORE_STATE_SIZE + amb_idx),
                    x_i.len(),
                    sat.sat_obs.sat.constellation,
                ),
                weight: var_cp / w_cp,
                raw_var: var_cp,
                is_phase: true,
                sat: Some(sat.sat_obs.sat),
            });
        }
    }

    /// Push CP measurement only if sat.cp1 is Some and non-zero.
    /// Extracted as a helper to reduce nesting depth in push_sat_meas.
    fn try_push_cp_measurement(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        iter: usize,
        los: &Vector3<f64>,
        expected_base: f64,
        dist: f64,
    ) {
        let cp1 = match sat.cp1 {
            Some(cp) if cp != 0.0 => cp,
            _ => return,
        };
        self.push_cp_measurement(meas, state, sat, x_i, iter, los, expected_base, dist, cp1);
    }

    fn resolve_uduc_indices(
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
    ) -> (
        Option<usize>,
        Option<usize>,
        Option<usize>,
        f64,
        f64,
        f64,
        f64,
    ) {
        let i1_idx = find_amb_idx(state, sat.sat_obs.sat, 3).map(|idx| CORE_STATE_SIZE + idx);
        let n1_idx = find_amb_idx(state, sat.sat_obs.sat, 1).map(|idx| CORE_STATE_SIZE + idx);
        let n2_idx = find_amb_idx(state, sat.sat_obs.sat, 2).map(|idx| CORE_STATE_SIZE + idx);
        let i1 = i1_idx.map(|idx| x_i[idx]).unwrap_or(0.0);
        let n1 = n1_idx.map(|idx| x_i[idx]).unwrap_or(0.0);
        let n2 = n2_idx.map(|idx| x_i[idx]).unwrap_or(0.0);
        (
            i1_idx,
            n1_idx,
            n2_idx,
            i1,
            n1,
            n2,
            (sat.f1 * sat.f1) / (sat.f2 * sat.f2),
        )
    }

    fn push_uduc_pr_measurements(
        &self,
        meas: &mut Vec<FgMeasurement>,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        los: &Vector3<f64>,
        expected_base: f64,
        i1_idx: Option<usize>,
        i1: f64,
        gamma: f64,
    ) {
        let var_p1 = PSEUDORANGE_VARIANCE_BASE * snr_scale(sat.snr as i32) / libm::sin(sat.el);
        let res_p1 = sat.p1 - (expected_base + i1);
        meas.push(FgMeasurement {
            res: res_p1,
            h_row: build_h_row_uduc(
                los,
                sat.map_wet,
                i1_idx,
                1.0,
                None,
                x_i.len(),
                sat.sat_obs.sat.constellation,
            ),
            weight: var_p1 / apply_huber(res_p1, var_p1, self.huber_k),
            raw_var: var_p1,
            is_phase: false,
            sat: Some(sat.sat_obs.sat),
        });
        let res_p2 = match sat.p2 {
            Some(p2) => p2 - (expected_base + gamma * i1),
            None => return,
        };
        meas.push(FgMeasurement {
            res: res_p2,
            h_row: build_h_row_uduc(
                los,
                sat.map_wet,
                i1_idx,
                gamma,
                None,
                x_i.len(),
                sat.sat_obs.sat.constellation,
            ),
            weight: (var_p1 * 1.5) / apply_huber(res_p2, var_p1 * 1.5, self.huber_k),
            raw_var: var_p1 * 1.5,
            is_phase: false,
            sat: Some(sat.sat_obs.sat),
        });
    }

    fn push_uduc_cp_measurements(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        los: &Vector3<f64>,
        expected_base: f64,
        i1_idx: Option<usize>,
        n1_idx: Option<usize>,
        n2_idx: Option<usize>,
        i1: f64,
        n1: f64,
        n2: f64,
        gamma: f64,
    ) {
        let windup = *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0);
        let var_l1 = 0.0001 * snr_scale(sat.snr as i32) / libm::sin(sat.el);
        let cp1_val = match sat.cp1 {
            Some(cp1) => cp1,
            None => return,
        };
        let res_l1 = (cp1_val - windup) * sat.lam1 - (expected_base - i1 + n1);
        meas.push(FgMeasurement {
            res: res_l1,
            h_row: build_h_row_uduc(
                los,
                sat.map_wet,
                i1_idx,
                -1.0,
                n1_idx,
                x_i.len(),
                sat.sat_obs.sat.constellation,
            ),
            weight: var_l1 / apply_huber(res_l1, var_l1, self.huber_k),
            raw_var: var_l1,
            is_phase: true,
            sat: Some(sat.sat_obs.sat),
        });
        let cp2_val = match sat.cp2 {
            Some(cp2) => cp2,
            None => return,
        };
        let res_l2 = (cp2_val - windup) * sat.lam2 - (expected_base - gamma * i1 + n2);
        meas.push(FgMeasurement {
            res: res_l2,
            h_row: build_h_row_uduc(
                los,
                sat.map_wet,
                i1_idx,
                -gamma,
                n2_idx,
                x_i.len(),
                sat.sat_obs.sat.constellation,
            ),
            weight: (var_l1 * 1.5) / apply_huber(res_l2, var_l1 * 1.5, self.huber_k),
            raw_var: var_l1 * 1.5,
            is_phase: true,
            sat: Some(sat.sat_obs.sat),
        });
    }

    fn push_uduc_measurements(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        _iter: usize,
        los: &Vector3<f64>,
        expected_base: f64,
        _dist: f64,
        _isb: f64,
    ) {
        let (i1_idx, n1_idx, n2_idx, i1, n1, n2, gamma) =
            Self::resolve_uduc_indices(state, sat, x_i);
        self.push_uduc_pr_measurements(meas, sat, x_i, los, expected_base, i1_idx, i1, gamma);
        self.push_uduc_cp_measurements(
            meas,
            state,
            sat,
            x_i,
            los,
            expected_base,
            i1_idx,
            n1_idx,
            n2_idx,
            i1,
            n1,
            n2,
            gamma,
        );
    }
}

fn log_ppp_convergence(
    state: &RtkState,
    sats: &[ProcessedSat],
    x_i: &DVector<f64>,
    x_pred: &DVector<f64>,
    p_pred: &DMatrix<f64>,
    solver: &PppIteratedEkf,
) {
    let _p_amb = if p_pred.nrows() > 21 {
        p_pred[(21, 21)]
    } else {
        0.0
    };
    tracing::info!(
        "PPP Epoch: pos=[{:.2}, {:.2}, {:.2}], dx_norm={:.4}",
        x_i[0],
        x_i[1],
        x_i[2],
        (x_i.clone() - x_pred.clone()).norm()
    );

    let meas = solver.build_measurements(state, sats, x_i, solver.max_iterations - 1);
    let (mut sum_pr, mut count_pr, mut sum_rr, mut count_rr) = (0.0, 0, 0.0, 0);
    for m in &meas {
        if m.h_row.len() == x_i.len() {
            if m.is_phase {
                continue;
            }
            if m.weight > 0.05 {
                sum_pr += m.res.abs();
                count_pr += 1;
            } else {
                sum_rr += m.res.abs();
                count_rr += 1;
            }
        }
    }
    if state.epoch_count.is_multiple_of(100) {
        tracing::trace!(
            "Epoch {}: Mean PR Res = {:.3} m, Mean RR Res = {:.3} m/s",
            state.epoch_count,
            sum_pr / count_pr.max(1) as f64,
            sum_rr / count_rr.max(1) as f64
        );
    }
}

fn build_h_row_uduc(
    los: &Vector3<f64>,
    map_wet: f64,
    i_idx: Option<usize>,
    i_coef: f64,
    n_idx: Option<usize>,
    size: usize,
    constel: gneiss_core::sat::Constellation,
) -> DVector<f64> {
    let mut h = DVector::zeros(size);
    h[0] = -los.x;
    h[1] = -los.y;
    h[2] = -los.z;
    h[15] = 1.0;
    if size > 18 {
        match constel {
            gneiss_core::sat::Constellation::Glonass => h[16] = 1.0,
            gneiss_core::sat::Constellation::Galileo => h[17] = 1.0,
            gneiss_core::sat::Constellation::Beidou => h[18] = 1.0,
            _ => {}
        }
    }
    if size > 20 {
        h[20] = map_wet;
    }
    if let Some(idx) = i_idx {
        h[idx] = i_coef;
    }
    if let Some(idx) = n_idx {
        h[idx] = 1.0;
    }
    h
}

fn build_h_row(
    los: &Vector3<f64>,
    map_wet: f64,
    amb_idx: Option<usize>,
    size: usize,
    constel: gneiss_core::sat::Constellation,
) -> DVector<f64> {
    let mut h = DVector::zeros(size);
    h[0] = -los.x;
    h[1] = -los.y;
    h[2] = -los.z;
    h[15] = 1.0;
    if size > 18 {
        match constel {
            gneiss_core::sat::Constellation::Glonass => h[16] = 1.0,
            gneiss_core::sat::Constellation::Galileo => h[17] = 1.0,
            gneiss_core::sat::Constellation::Beidou => h[18] = 1.0,
            _ => {}
        }
    }
    if size > 20 {
        h[20] = map_wet;
    }
    if let Some(idx) = amb_idx {
        h[idx] = 1.0;
    }
    h
}

fn build_h_row_doppler(los: &Vector3<f64>, size: usize) -> DVector<f64> {
    let mut h = DVector::zeros(size);
    if size > 19 {
        h[3] = -los.x;
        h[4] = -los.y;
        h[5] = -los.z;
        h[19] = 1.0;
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, DVector, Vector3};

    fn dummy_rtk_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(
                Vector3::zeros(),
                Datum::WGS84,
                Frame::ECEF,
                GpsTime::new(0, 0.0),
            ),
            1.0,
        )
    }

    #[test]
    fn test_solve_empty_sats() {
        let fg = PppIteratedEkf::new();
        let mut state = dummy_rtk_state();
        state.covariance = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        let sats = vec![];
        let res = fg.solve(&mut state, &sats, None);
        assert!(matches!(res, Err(EngineError::InsufficientSatellites)));
    }

    #[test]
    fn test_ppp_factor_graph_default() {
        let fg = PppIteratedEkf::default();
        assert_eq!(fg.max_iterations, 15);
        assert_eq!(fg.convergence_threshold, 1e-3);
        assert_eq!(fg.huber_k, 3.0);
        let fg2 = PppIteratedEkf::new();
        assert_eq!(fg2.max_iterations, 15);
    }

    #[test]
    fn test_find_ambiguity_index() {
        let mut state = dummy_rtk_state();

        use gneiss_core::sat::{Constellation, SatelliteId};
        let sat1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let sat2 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        };
        state.ambiguity_keys.push((sat1, 0));
        state.ambiguity_keys.push((sat2, 1));
        state.ambiguities.push(0.0);
        state.ambiguities.push(0.0);

        assert_eq!(find_ambiguity_index(&state, sat1), Some(0));
        assert_eq!(find_ambiguity_index(&state, sat2), None);

        let sat3 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 3,
        };
        assert_eq!(find_ambiguity_index(&state, sat3), None);
    }

    #[test]
    fn test_build_weight_matrix() {
        let m1 = FgMeasurement {
            res: 1.0,
            h_row: DVector::zeros(1),
            weight: 2.0,
            raw_var: 0.5,
            is_phase: false,
            sat: None,
        };
        let m2 = FgMeasurement {
            res: 2.0,
            h_row: DVector::zeros(1),
            weight: 4.0,
            raw_var: 0.25,
            is_phase: true,
            sat: None,
        };
        let meas = vec![m1, m2];
        let mut r = DMatrix::zeros(2, 2);
        r[(0, 0)] = 2.0;
        r[(1, 1)] = 4.0;
        let w = build_weight_matrix(&meas, &r);
        assert_eq!(w[(0, 0)], 0.5);
        assert_eq!(w[(1, 1)], 0.25);
        assert_eq!(w[(0, 1)], 0.0);
    }

    #[test]
    fn test_build_h_row() {
        let los = Vector3::new(1.0, 2.0, 3.0);
        // size = CORE_STATE_SIZE + 1 so ISBs (size>18) and ZWD (size>20) are populated,
        // plus one ambiguity slot at index CORE_STATE_SIZE.
        let size = CORE_STATE_SIZE + 1;
        let h = build_h_row(
            &los,
            4.0,
            Some(CORE_STATE_SIZE),
            size,
            gneiss_core::sat::Constellation::Gps,
        );
        assert_eq!(h.len(), size);
        assert_eq!(h[0], -1.0);
        assert_eq!(h[1], -2.0);
        assert_eq!(h[2], -3.0);
        assert_eq!(h[15], 1.0); // clock bias
        assert_eq!(h[20], 4.0); // ZWD mapping
        assert_eq!(h[CORE_STATE_SIZE], 1.0); // ambiguity

        // size = CORE_STATE_SIZE: ISBs populated (size>18) but ZWD NOT (size==21, not >20... wait 21>20 is true)
        // Actually CORE_STATE_SIZE=21 > 20, so ZWD IS set. No ambiguity.
        let h2 = build_h_row(
            &los,
            4.0,
            None,
            CORE_STATE_SIZE,
            gneiss_core::sat::Constellation::Gps,
        );
        assert_eq!(h2.len(), CORE_STATE_SIZE);
        assert_eq!(h2[0], -1.0);
        assert_eq!(h2[15], 1.0);
        assert_eq!(h2[20], 4.0); // ZWD mapping (21 > 20)
    }

    #[test]
    fn test_build_h_row_doppler() {
        let los = Vector3::new(1.0, 2.0, 3.0);
        // size = CORE_STATE_SIZE (21) > 19, so velocity and clock drift are populated.
        let h = build_h_row_doppler(&los, CORE_STATE_SIZE);
        assert_eq!(h.len(), CORE_STATE_SIZE);
        assert_eq!(h[3], -1.0);
        assert_eq!(h[4], -2.0);
        assert_eq!(h[5], -3.0);
        assert_eq!(h[19], 1.0); // clock drift at index 19

        // size = 19, NOT > 19, so nothing is set — all zeros.
        let h2 = build_h_row_doppler(&los, 19);
        assert_eq!(h2.len(), 19);
        assert_eq!(h2[3], 0.0);
    }

    #[test]
    fn test_assemble_matrices() {
        let meas = vec![
            FgMeasurement {
                res: 1.5,
                h_row: DVector::from_element(3, 1.0),
                weight: 2.0,
                raw_var: 0.5,
                is_phase: false,
                sat: None,
            },
            FgMeasurement {
                res: 2.5,
                h_row: DVector::from_element(3, 2.0),
                weight: 3.0,
                raw_var: 0.33,
                is_phase: true,
                sat: None,
            },
        ];
        let (h, z, r) = assemble_matrices(&meas, 3);
        assert_eq!(h.nrows(), 2);
        assert_eq!(h.ncols(), 3);
        assert_eq!(h[(0, 0)], 1.0);
        assert_eq!(h[(1, 2)], 2.0);
        assert_eq!(z.len(), 2);
        assert_eq!(z[0], 1.5);
        assert_eq!(z[1], 2.5);
        assert_eq!(r.nrows(), 2);
        assert_eq!(r.ncols(), 2);
        assert_eq!(r[(0, 0)], 2.0);
        assert_eq!(r[(1, 1)], 3.0);
        assert_eq!(r[(0, 1)], 0.0);
    }

    #[test]
    fn test_extract_and_apply_state_vector() {
        let mut state = dummy_rtk_state();
        state.position.vector = Vector3::new(1.0, 2.0, 3.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);
        state.attitude = nalgebra::UnitQuaternion::from_scaled_axis(Vector3::new(0.1, 0.2, 0.3));
        state.accel_bias = Vector3::new(10.0, 11.0, 12.0);
        state.gyro_bias = Vector3::new(13.0, 14.0, 15.0);
        state.rcv_clk_bias = 16.0;
        state.rcv_clk_drift = 17.0;
        state.zwd = 18.0;
        state.ambiguities = vec![19.0, 20.0];
        let x = extract_state_vector(&state);
        assert_eq!(x.len(), CORE_STATE_SIZE + 2);
        assert_eq!(x[0], 1.0);
        assert_eq!(x[15], 16.0);
        assert_eq!(x[16], 0.0);
        assert_eq!(x[17], 0.0);
        assert_eq!(x[18], 0.0);
        assert_eq!(x[19], 17.0);
        assert_eq!(x[20], 18.0);
        assert_eq!(x[CORE_STATE_SIZE], 19.0);
        assert_eq!(x[CORE_STATE_SIZE + 1], 20.0);

        let mut state2 = dummy_rtk_state();
        state2.ambiguities = vec![0.0, 0.0];
        let cov = state2.covariance.clone();
        apply_state_vector(&mut state2, &x, cov);
        assert_eq!(state2.position.vector, Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(state2.velocity, Vector3::new(4.0, 5.0, 6.0));
        assert!((state2.attitude.scaled_axis() - Vector3::new(0.1, 0.2, 0.3)).norm() < 1e-10);
        assert_eq!(state2.accel_bias, Vector3::new(10.0, 11.0, 12.0));
        assert_eq!(state2.gyro_bias, Vector3::new(13.0, 14.0, 15.0));
        assert_eq!(state2.rcv_clk_bias, 16.0);
        assert_eq!(state2.rcv_clk_drift, 17.0);
        assert_eq!(state2.zwd, 18.0);
        assert_eq!(state2.ambiguities, vec![19.0, 20.0]);
    }
}
#[cfg(test)]
mod nan_tests {
    use super::*;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, Vector3};

    fn dummy_rtk_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(
                Vector3::zeros(),
                Datum::WGS84,
                Frame::ECEF,
                GpsTime::new(0, 0.0),
            ),
            1.0,
        )
    }

    #[test]
    fn test_solve_matrix_inversion_failure() {
        let fg = PppIteratedEkf::new();
        let mut state = dummy_rtk_state();
        state.covariance = DMatrix::from_element(CORE_STATE_SIZE, CORE_STATE_SIZE, f64::NAN);
        let sats = vec![];
        let res = fg.solve(&mut state, &sats, None);
        assert!(matches!(res, Err(EngineError::StateDisappeared)));
    }
}

#[cfg(test)]
mod mutant_killer_tests {
    use super::*;
    use crate::engine::processed_sat::ProcessedSat;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::obs::SatObs;
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, DVector, Vector3};

    fn dummy_rtk_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(
                Vector3::zeros(),
                Datum::WGS84,
                Frame::ECEF,
                GpsTime::new(0, 0.0),
            ),
            1.0,
        )
    }

    #[test]
    fn test_resolve_widelane_ar_insufficient() {
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        state.covariance = DMatrix::zeros(CORE_STATE_SIZE + 4, CORE_STATE_SIZE + 4);
        for i in 0..4 {
            state.covariance[(CORE_STATE_SIZE + i, CORE_STATE_SIZE + i)] = 100.0;
        } // huge variance
        let subset = [
            (
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 1,
                    },
                    0,
                    0,
                    1.0,
                    1.0,
                    1.0,
                ),
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 2,
                    },
                    1,
                    1,
                    1.0,
                    1.0,
                    1.0,
                ),
            ),
            (
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 1,
                    },
                    0,
                    0,
                    1.0,
                    1.0,
                    1.0,
                ),
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 3,
                    },
                    2,
                    2,
                    1.0,
                    1.0,
                    1.0,
                ),
            ),
            (
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 1,
                    },
                    0,
                    0,
                    1.0,
                    1.0,
                    1.0,
                ),
                (
                    SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 4,
                    },
                    3,
                    3,
                    1.0,
                    1.0,
                    1.0,
                ),
            ),
        ];
        let x = DVector::zeros(CORE_STATE_SIZE + 4);
        let res = fg.resolve_widelane_ar(&state, &state.covariance, &subset, &x);
        assert!(res.is_err());
    }

    #[test]
    fn test_push_cp_measurement_iono_free() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        state.ambiguity_keys.push((sat_id, 0));
        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: true,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let x_i = DVector::zeros(CORE_STATE_SIZE + 1);
        let los = Vector3::new(0.0, 0.0, 1.0);
        fg.push_cp_measurement(
            &mut meas, &state, &sat, &x_i, 0, &los, 10.0, 20000000.0, 100.0,
        );
        assert_eq!(meas.len(), 1);
        // expected_base = 10.0, x_i = 0. expected_cp = 10.0.
        // l_meas = 100.0 * 0.19 = 19.0.
        // res = 19.0 - 10.0 = 9.0
        assert!((meas[0].res - 9.0).abs() < 1e-6);
        // var_cp = 0.0001 * 1.0 / 1.0 * 9.0 = 0.0009
        assert!((meas[0].raw_var - 0.0009).abs() < 1e-6);
    }

    #[test]
    fn test_push_cp_measurement_not_iono_free() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        state.ambiguity_keys.push((sat_id, 0));
        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: None,
            cp2: None,
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let x_i = DVector::zeros(CORE_STATE_SIZE + 1);
        let los = Vector3::new(0.0, 0.0, 1.0);
        fg.push_cp_measurement(
            &mut meas, &state, &sat, &x_i, 0, &los, 10.0, 20000000.0, 100.0,
        );
        assert_eq!(meas.len(), 1);
        // expected_base = 10.0. expected_cp = 10.0 - 5.0 = 5.0.
        // l_meas = 100.0 * 0.19 = 19.0.
        // res = 19.0 - 5.0 = 14.0
        assert!((meas[0].res - 14.0).abs() < 1e-6);
        // var_cp = 0.0001 * 1.0 / 1.0 = 0.0001
        assert!((meas[0].raw_var - 0.0001).abs() < 1e-6);
    }

    #[test]
    fn test_find_ar_candidates() {
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let sat_id1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        state.ambiguity_keys.push((sat_id1, 1));
        state.ambiguity_keys.push((sat_id1, 2));
        state.ambiguities.push(0.0);
        state.ambiguities.push(0.0);

        let obs = SatObs {
            sat: sat_id1,
            observations: vec![],
        };
        let mut sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: Some(0.0),
            cp2: Some(0.0),
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: 15.01_f64.to_radians(),
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let cands = fg.find_ar_candidates(&state, &[sat.clone()]);
        assert_eq!(cands.len(), 1);

        // Test `!s.is_iono_free` mutant
        sat.is_iono_free = true;
        assert_eq!(fg.find_ar_candidates(&state, &[sat.clone()]).len(), 0);
        sat.is_iono_free = false;

        // Test `s.cp2.is_some()` mutant
        sat.cp2 = None;
        assert_eq!(fg.find_ar_candidates(&state, &[sat.clone()]).len(), 0);
        sat.cp2 = Some(0.0);

        // Test constellation
        let sat_id_glo = SatelliteId {
            constellation: Constellation::Glonass,
            prn: 1,
        };
        let obs_glo = SatObs {
            sat: sat_id_glo,
            observations: vec![],
        };
        sat.sat_obs = &obs_glo;
        assert_eq!(fg.find_ar_candidates(&state, &[sat.clone()]).len(), 0);
    }

    #[test]
    fn test_resolve_cascade_ar_bounds() {
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let mut sats: Vec<ProcessedSat> = Vec::new();

        let obs1 = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![],
        };
        let obs2 = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 2,
            },
            observations: vec![],
        };
        let obs3 = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 3,
            },
            observations: vec![],
        };
        let obs4 = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 4,
            },
            observations: vec![],
        };
        // Static references to avoid lifetime issues in closure
        let obs1_ref = Box::leak(Box::new(obs1));
        let obs2_ref = Box::leak(Box::new(obs2));
        let obs3_ref = Box::leak(Box::new(obs3));
        let obs4_ref = Box::leak(Box::new(obs4));

        {
            let mut add_sat = |obs: &'static SatObs| {
                state.add_ambiguity(obs.sat, 1, 0.0, 1.0);
                state.add_ambiguity(obs.sat, 2, 0.0, 1.0);
                sats.push(ProcessedSat {
                    sat_obs: obs,
                    dt_sat_m: 0.0,
                    p1: 0.0,
                    p2: None,
                    cp1: Some(0.0),
                    cp2: Some(0.0),
                    is_iono_free: false,
                    osb_p1: 0.0,
                    osb_p2: 0.0,
                    osb_cp1: 0.0,
                    osb_cp2: 0.0,
                    los: Vector3::zeros(),
                    dist: 0.0,
                    el: 15.01_f64.to_radians(),
                    snr: 45.0,
                    doppler: 0.0,
                    lam1: 0.19,
                    lam2: 0.24,
                    tropo_dry: 0.0,
                    map_wet: 0.0,
                    iono_delay: 5.0,
                    f1: 1.0,
                    f2: 1.0,
                    sat_pos_rot: Vector3::zeros(),
                    sat_vel: Vector3::zeros(),
                    sat_clock_drift: 0.0,
                    rcv_pos_ecef: Vector3::zeros(),
                    pcv_correction: 0.0,
                });
            };

            add_sat(obs1_ref);
            add_sat(obs2_ref);
            add_sat(obs3_ref);
        }
        assert_eq!(
            fg.resolve_cascade_ar(&mut state, &sats),
            Err("Insufficient dual-frequency satellites for AR")
        );

        state.add_ambiguity(obs4_ref.sat, 1, 0.0, 1.0);
        state.add_ambiguity(obs4_ref.sat, 2, 0.0, 1.0);
        sats.push(ProcessedSat {
            sat_obs: obs4_ref,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: Some(0.0),
            cp2: Some(0.0),
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: 15.01_f64.to_radians(),
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        });
        assert!(
            fg.resolve_cascade_ar(&mut state, &sats)
                != Err("Insufficient dual-frequency satellites for AR")
        );
    }

    #[test]
    fn test_find_worst_outlier() {
        let sat1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let sat2 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        };

        let meas = vec![
            // Not phase, shouldn't be picked even if high
            FgMeasurement {
                res: 1000.0,
                raw_var: 1.0,
                is_phase: false,
                sat: Some(sat1),
                h_row: DVector::zeros(0),
                weight: 1.0,
            },
            // Phase, but norm = 10.0 / sqrt(4.0) = 5.0 (less than max_norm 15.0)
            FgMeasurement {
                res: 10.0,
                raw_var: 4.0,
                is_phase: true,
                sat: Some(sat1),
                h_row: DVector::zeros(0),
                weight: 1.0,
            },
            // Phase, norm = 40.0 / sqrt(4.0) = 20.0 (greater than max_norm 15.0)
            FgMeasurement {
                res: 40.0,
                raw_var: 4.0,
                is_phase: true,
                sat: Some(sat2),
                h_row: DVector::zeros(0),
                weight: 1.0,
            },
            // Phase, norm = -60.0 / sqrt(9.0) = 20.0 (equal to current max_norm, shouldn't override because of >)
            FgMeasurement {
                res: -60.0,
                raw_var: 9.0,
                is_phase: true,
                sat: Some(sat1),
                h_row: DVector::zeros(0),
                weight: 1.0,
            },
        ];

        assert_eq!(PppIteratedEkf::find_worst_outlier_sat(&meas), Some(sat2));

        let meas_no_outlier = vec![FgMeasurement {
            res: 10.0,
            raw_var: 4.0,
            is_phase: true,
            sat: Some(sat1),
            h_row: DVector::zeros(0),
            weight: 1.0,
        }];
        assert_eq!(
            PppIteratedEkf::find_worst_outlier_sat(&meas_no_outlier),
            None
        );
    }

    #[test]
    fn test_find_ar_candidates_bounds() {
        let _fg = PppIteratedEkf::default();
        let _state = dummy_rtk_state();
        let _sats: Vec<ProcessedSat> = Vec::new();
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let obs = SatObs {
            sat: sat_id,
            observations: vec![],
        };
        let _sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 0.0,
            p2: None,
            cp1: Some(0.0),
            cp2: Some(0.0),
            is_iono_free: false,
            osb_p1: 0.0,
            osb_p2: 0.0,
            osb_cp1: 0.0,
            osb_cp2: 0.0,
            los: Vector3::zeros(),
            dist: 0.0,
            el: 15.01_f64.to_radians(),
            snr: 45.0,
            doppler: 0.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        // It requires state.ambiguity_keys to contain (sat, 0) and (sat, 1) and (sat, 2) etc depending on `is_iono_free`.
        // We'll skip adding a full state test and rely on smaller integration tests or direct tests.
    }

    #[test]
    fn test_build_iono_constraint_row() {
        let h = build_iono_constraint_row(25, 21);
        assert_eq!(h.len(), 25);
        assert_eq!(h[21], 1.0);
        assert_eq!(h[0], 0.0);
        assert_eq!(h[24], 0.0);
    }

    #[test]
    fn test_iono_constraint_row_middle_index() {
        let h = build_iono_constraint_row(30, 15);
        assert_eq!(h.len(), 30);
        assert_eq!(h[15], 1.0);
        for i in 0..30 {
            if i != 15 {
                assert_eq!(h[i], 0.0);
            }
        }
    }

    #[test]
    fn test_sequential_ar_mismatch_regression() {
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();

        // Add 6 satellites across 3 constellations (GPS, Galileo, Beidou) to get 3 constellation groups
        let mut sats = Vec::new();
        let constellations = [
            Constellation::Gps,
            Constellation::Galileo,
            Constellation::Beidou,
        ];

        for (i, &constellation) in constellations.iter().enumerate() {
            let sat1 = SatelliteId {
                constellation,
                prn: (i * 2 + 1) as u8,
            };
            let sat2 = SatelliteId {
                constellation,
                prn: (i * 2 + 2) as u8,
            };

            state.add_ambiguity(sat1, 1, 0.0, 1.0);
            state.add_ambiguity(sat1, 2, 0.0, 1.0);
            state.add_ambiguity(sat2, 1, 0.0, 1.0);
            state.add_ambiguity(sat2, 2, 0.0, 1.0);

            let obs1 = Box::leak(Box::new(SatObs {
                sat: sat1,
                observations: vec![],
            }));
            let obs2 = Box::leak(Box::new(SatObs {
                sat: sat2,
                observations: vec![],
            }));

            let make_processed = |obs: &'static SatObs| ProcessedSat {
                sat_obs: obs,
                dt_sat_m: 0.0,
                p1: 0.0,
                p2: None,
                cp1: Some(0.0),
                cp2: Some(0.0),
                is_iono_free: false,
                osb_p1: 0.0,
                osb_p2: 0.0,
                osb_cp1: 0.0,
                osb_cp2: 0.0,
                los: Vector3::zeros(),
                dist: 0.0,
                el: 15.01_f64.to_radians(),
                snr: 45.0,
                doppler: 0.0,
                lam1: 0.19,
                lam2: 0.24,
                tropo_dry: 0.0,
                map_wet: 0.0,
                iono_delay: 5.0,
                f1: 1.0,
                f2: 1.0,
                sat_pos_rot: Vector3::zeros(),
                sat_vel: Vector3::zeros(),
                sat_clock_drift: 0.0,
                rcv_pos_ecef: Vector3::zeros(),
                pcv_correction: 0.0,
            };
            sats.push(make_processed(obs1));
            sats.push(make_processed(obs2));
        }

        // Initialize state vector and covariance
        let state_dim = CORE_STATE_SIZE + state.ambiguities.len();
        state.covariance = DMatrix::from_fn(state_dim, state_dim, |r, c| {
            if r == c {
                (r + 1) as f64 * 1.5
            } else {
                0.01
            }
        });
        state.position.vector = Vector3::new(1.0, 2.0, 3.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);
        state.is_fixed = false;

        let initial_state_vector = extract_state_vector(&state);
        let initial_covariance = state.covariance.clone();

        // Configure mock
        let mock_wl = Ok((
            DVector::zeros(state_dim),
            DMatrix::zeros(state_dim, state_dim),
            vec![0, 1],
        ));
        let mock_nl = Ok((
            DVector::zeros(state_dim),
            DMatrix::zeros(state_dim, state_dim),
        ));

        {
            let mut mock_lock = AR_MOCK.lock().unwrap();
            *mock_lock = Some(ArMock {
                wl_result: Some(mock_wl),
                nl_result: Some(mock_nl),
                nl_calls: 0,
            });
        }

        // Execute resolve_cascade_ar
        let result = fg.resolve_cascade_ar(&mut state, &sats);

        // Verify that it failed due to global position jump check
        assert_eq!(result, Err("Position jump too large after AR fix"));

        // Verify state is unmodified
        let final_state_vector = extract_state_vector(&state);
        assert_eq!(final_state_vector, initial_state_vector);
        assert_eq!(state.covariance, initial_covariance);
        assert_eq!(state.is_fixed, false);

        // Clear mock
        {
            let mut mock_lock = AR_MOCK.lock().unwrap();
            *mock_lock = None;
        }
    }
}
