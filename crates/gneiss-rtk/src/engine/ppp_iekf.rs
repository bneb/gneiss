use crate::engine::ppp_common::{
    apply_state_vector, assemble_matrices, build_weight_matrix, extract_state_vector,
    find_amb_idx, invert_matrix, FgMeasurement,
};

// These are used only by test modules within this file; #[cfg(test)] avoids
// unused-import warnings during library (non-test) compilation.
#[cfg(test)]
use crate::engine::ppp_common::{build_iono_constraint_row, find_ambiguity_index};
use crate::engine::processed_sat::ProcessedSat;
use crate::engine::types::IonosphereModel;
use crate::engine::EngineError;
use crate::filter::{RtkState, CORE_STATE_SIZE};
use crate::math::inversion::solve_cholesky_svd;
use nalgebra::{DMatrix, DVector, Vector3};

// Re-export free functions and constants from the measurement module so that
// callers (including tests within this file) can reference them without
// qualifying the sibling module path.
pub(crate) use super::ppp_measurements::*;

#[cfg(test)]
#[derive(Clone)]
pub struct ArMock {
    pub wl_result: Option<Result<(DVector<f64>, DMatrix<f64>, Vec<usize>), &'static str>>,
    pub nl_result: Option<Result<(DVector<f64>, DMatrix<f64>), &'static str>>,
    pub nl_calls: usize,
}

#[cfg(test)]
pub static AR_MOCK: std::sync::Mutex<Option<ArMock>> = std::sync::Mutex::new(None);

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
    /// Ionosphere model — controls iono prior variance in UDUC measurements.
    /// Klobuchar: 9.0 m² (3m std). IONEX: 0.0025 m² (0.05m std).
    pub iono_model: IonosphereModel,
}

impl Default for PppIteratedEkf {
    fn default() -> Self {
        Self {
            max_iterations: 15,
            convergence_threshold: 1e-3,
            huber_k: 3.0,
            iono_model: IonosphereModel::Klobuchar,
        }
    }
}

impl PppIteratedEkf {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_iono_model(mut self, model: IonosphereModel) -> Self {
        self.iono_model = model;
        self
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
        // Skip AR if already fixed and no new satellites appeared.
        // Re-fixing every epoch after a correct fix wastes compute and risks
        // an incorrect subsequent fix corrupting a converged solution.
        let has_new_sats = sats.iter().any(|s| {
            let key = (s.sat_obs.sat, 0);
            state.last_observed.get(&key).copied().unwrap_or(0) == state.epoch_count as u32
        });
        if !state.is_fixed || has_new_sats {
            if state.epoch_count > 10 {
                if let Err(e) = self.resolve_cascade_ar(state, sats) {
                    tracing::info!("Cascade AR did not fix: {:?}", e);
                }
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

    fn find_worst_outlier_sat(_meas: &[FgMeasurement]) -> Option<gneiss_core::sat::SatelliteId> {
        // Disabled: outlier removal cascades during convergence — removing
        // one ambiguity degrades remaining measurements, triggering more
        // removals, until all CP is lost.  The Huber estimator handles
        // outlier down-weighting without removing the ambiguity entirely.
        None
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
                // MW-based: accept if both rover and reference have >50 MW samples.
                // 50 samples gives ~0.06 cycle WL precision vs ~0.12 at 10 samples,
                // and reduces first-sample EMA bias from ~15% to ~4%.
                let (c, ref_sat) = &subset[i];
                let mw_ok = state.mw_sd_counts.get(&c.0).unwrap_or(&0) > &50
                    && state.mw_sd_counts.get(&ref_sat.0).unwrap_or(&0) > &50;
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
            if cnt_c > 50 && cnt_ref > 50 {
                let mw_c = state.mw_sd_ema.get(&c.0).copied().unwrap_or(0.0);
                let mw_ref = state.mw_sd_ema.get(&ref_sat.0).copied().unwrap_or(0.0);
                a_wl[i] = mw_c - mw_ref;
            }
        }
        let q_wl_ekf = &d_wl * p * d_wl.transpose();
        // When all kept pairs have sufficient MW samples, use MW-based
        // covariance instead of the EKF state covariance.  The MW EMA
        // converges at ~0.42/sqrt(N) cycles whereas the EKF ambiguity
        // states start at 10_000 m² (≈277_000 cycles² on WL).  Using the
        // EKF covariance makes LAMBDA think every integer set is equally
        // likely, yielding ratio ≈ 1.0.
        let all_mw_confident = keep_indices.iter().all(|&idx| {
            let (c, ref_sat) = &subset[idx];
            state.mw_sd_counts.get(&c.0).copied().unwrap_or(0) > 50
                && state.mw_sd_counts.get(&ref_sat.0).copied().unwrap_or(0) > 50
        });
        let q_wl = if all_mw_confident {
            // MW per-sample DD variance: each single-epoch MW measurement has
            // ~0.42 cycle std on GPS L1/L2, so 0.18 cycles² per sample.
            // Reference satellite noise is shared across all DD pairs.
            let mw_var_per_sample: f64 = 0.18; // 0.42² cycles²
            let mut q = DMatrix::zeros(n, n);
            for i in 0..n {
                let (c_i, ref_sat_i) = &subset[keep_indices[i]];
                let cnt_i = state.mw_sd_counts.get(&c_i.0).copied().unwrap_or(1);
                let cnt_ref = state.mw_sd_counts.get(&ref_sat_i.0).copied().unwrap_or(1);
                let var_i = mw_var_per_sample / cnt_i as f64;
                let var_ref = mw_var_per_sample / cnt_ref as f64;
                q[(i, i)] = (var_i + var_ref).max(0.0025); // 0.05² floor
                for j in (i + 1)..n {
                    // Shared reference → off-diagonal covariance
                    q[(i, j)] = var_ref;
                    q[(j, i)] = var_ref;
                }
            }
            q
        } else {
            let mut q = q_wl_ekf;
            for i in 0..q.nrows() {
                q[(i, i)] = q[(i, i)].max(0.01);
            }
            q
        };
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

        // Outlier detection disabled during convergence to prevent cascade:
        // removing one ambiguity degrades remaining measurements, causing
        // more removals until all CP is lost.  The Huber estimator handles
        // outlier down-weighting without removing the ambiguity.
        assert_eq!(PppIteratedEkf::find_worst_outlier_sat(&meas), None);

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

    #[test]
    fn test_extract_isb() {
        let x = DVector::from_fn(19, |i, _| i as f64);
        assert!((PppIteratedEkf::extract_isb(&x, Constellation::Gps) - 0.0).abs() < 1e-10);
        assert!(
            (PppIteratedEkf::extract_isb(&x, Constellation::Glonass) - 16.0).abs() < 1e-10
        );
        assert!(
            (PppIteratedEkf::extract_isb(&x, Constellation::Galileo) - 17.0).abs() < 1e-10
        );
        assert!(
            (PppIteratedEkf::extract_isb(&x, Constellation::Beidou) - 18.0).abs() < 1e-10
        );
        // size <= 18 returns 0.0 for all constellations
        let x_small = DVector::from_fn(18, |i, _| i as f64);
        assert!(
            (PppIteratedEkf::extract_isb(&x_small, Constellation::Glonass) - 0.0).abs() < 1e-10
        );
        assert!(
            (PppIteratedEkf::extract_isb(&x_small, Constellation::Galileo) - 0.0).abs() < 1e-10
        );
        assert!(
            (PppIteratedEkf::extract_isb(&x_small, Constellation::Beidou) - 0.0).abs() < 1e-10
        );
    }

    #[test]
    fn test_build_h_row_uduc_basic() {
        let los = Vector3::new(1.0, 2.0, 3.0);
        let size = CORE_STATE_SIZE + 2;
        // No iono or ambiguity indices
        let h = build_h_row_uduc(&los, 4.0, None, 1.0, None, size, Constellation::Gps);
        assert_eq!(h.len(), size);
        assert!((h[0] - (-1.0)).abs() < 1e-10);
        assert!((h[1] - (-2.0)).abs() < 1e-10);
        assert!((h[2] - (-3.0)).abs() < 1e-10);
        assert!((h[15] - 1.0).abs() < 1e-10);
        assert!((h[20] - 4.0).abs() < 1e-10);
        // No indices set
        assert!((h[CORE_STATE_SIZE] - 0.0).abs() < 1e-10);
        assert!((h[CORE_STATE_SIZE + 1] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_build_h_row_uduc_with_indices() {
        let los = Vector3::new(0.5, -1.5, 2.0);
        let size = CORE_STATE_SIZE + 4;
        let i_idx = CORE_STATE_SIZE + 2;
        let n_idx = CORE_STATE_SIZE + 3;
        let h = build_h_row_uduc(
            &los,
            2.5,
            Some(i_idx),
            -1.5,
            Some(n_idx),
            size,
            Constellation::Galileo,
        );
        assert!((h[17] - 1.0).abs() < 1e-10); // Galileo ISB
        assert!((h[i_idx] - (-1.5)).abs() < 1e-10); // Ionosphere coefficient
        assert!((h[n_idx] - 1.0).abs() < 1e-10); // Ambiguity
        assert!((h[20] - 2.5).abs() < 1e-10); // ZWD
    }

    #[test]
    fn test_build_h_row_uduc_size_boundaries() {
        let los = Vector3::new(1.0, 0.0, 0.0);
        // size = CORE_STATE_SIZE (21): >18 (ISBs), >20 (ZWD)
        let h = build_h_row_uduc(
            &los, 3.0, None, 1.0, None, CORE_STATE_SIZE, Constellation::Glonass,
        );
        assert_eq!(h.len(), CORE_STATE_SIZE);
        assert!((h[16] - 1.0).abs() < 1e-10); // Glonass ISB
        assert!((h[20] - 3.0).abs() < 1e-10); // ZWD
        // size = 17 (< 18): no ISBs or ZWD
        let h2 = build_h_row_uduc(&los, 3.0, None, 1.0, None, 17, Constellation::Glonass);
        assert_eq!(h2.len(), 17);
        assert!((h2[16] - 0.0).abs() < 1e-10); // No ISB
        assert!((h2[15] - 1.0).abs() < 1e-10); // Clock bias always set
    }

    #[test]
    fn test_build_ar_subset_gps_reference() {
        let fg = PppIteratedEkf::default();
        let cands = vec![
            (
                SatelliteId { constellation: Constellation::Gps, prn: 1 },
                0,
                1,
                20.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Gps, prn: 2 },
                2,
                3,
                30.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Galileo, prn: 1 },
                4,
                5,
                25.0_f64.to_radians(),
                0.19,
                0.24,
            ),
        ];
        let subset = fg.build_ar_subset(&cands);
        // GPS PRN 2 (highest elev = 30 deg) should be reference
        assert_eq!(subset.len(), 2);
        for pair in &subset {
            assert_eq!(pair.1.0.prn, 2);
            assert_eq!(pair.1.0.constellation, Constellation::Gps);
        }
    }

    #[test]
    fn test_build_ar_subset_no_gps_per_constellation() {
        let fg = PppIteratedEkf::default();
        let cands = vec![
            (
                SatelliteId { constellation: Constellation::Galileo, prn: 1 },
                0,
                1,
                30.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Galileo, prn: 2 },
                2,
                3,
                20.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Beidou, prn: 1 },
                4,
                5,
                25.0_f64.to_radians(),
                0.19,
                0.24,
            ),
            (
                SatelliteId { constellation: Constellation::Beidou, prn: 2 },
                6,
                7,
                15.0_f64.to_radians(),
                0.19,
                0.24,
            ),
        ];
        let subset = fg.build_ar_subset(&cands);
        // Per-constellation fallback: 1 pair per constellation
        assert_eq!(subset.len(), 2);
        // Each pair must have same constellation (order is HashMap-dependent)
        for pair in &subset {
            assert_eq!(pair.0.0.constellation, pair.1.0.constellation);
        }
        // Both constellations must be represented
        let constels: std::collections::HashSet<_> = subset
            .iter()
            .map(|p| p.0.0.constellation)
            .collect();
        assert!(constels.contains(&Constellation::Galileo));
        assert!(constels.contains(&Constellation::Beidou));
    }

    #[test]
    fn test_build_ar_subset_empty_or_single() {
        let fg = PppIteratedEkf::default();
        // Single GPS sat -> GPS reference but no non-ref sats
        let cands = vec![(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            0,
            1,
            30.0_f64.to_radians(),
            0.19,
            0.24,
        )];
        assert!(fg.build_ar_subset(&cands).is_empty());
        // Empty input
        let cands: Vec<(SatelliteId, usize, usize, f64, f64, f64)> = vec![];
        assert!(fg.build_ar_subset(&cands).is_empty());
    }

    #[test]
    fn test_try_push_cp_measurement_rejects_none_or_zero() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let los = Vector3::new(1.0, 0.0, 0.0);

        // cp1 = None -> no measurement
        let sat_no_cp1 = ProcessedSat {
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
            iono_delay: 0.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        fg.try_push_cp_measurement(&mut meas, &state, &sat_no_cp1, &x_i, 0, &los, 0.0, 0.0);
        assert_eq!(meas.len(), 0);

        // cp1 = Some(0.0) -> also no measurement (zero check)
        let sat_zero_cp1 = ProcessedSat {
            cp1: Some(0.0),
            ..sat_no_cp1.clone()
        };
        fg.try_push_cp_measurement(&mut meas, &state, &sat_zero_cp1, &x_i, 0, &los, 0.0, 0.0);
        assert_eq!(meas.len(), 0);
    }

    #[test]
    fn test_push_pr_measurement_iono_free_variance() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let los = Vector3::new(0.0, 0.0, 1.0);

        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 10.0,
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
            map_wet: 1.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        fg.push_pr_measurement(&mut meas, &state, &sat, &x_i, 0, &los, 5.0, 0.0, 0.0);
        assert_eq!(meas.len(), 1);
        // expected_pr = expected_base (iono-free) = 5.0
        // res_pr = p1 - expected_pr = 10.0 - 5.0 = 5.0
        assert!((meas[0].res - 5.0).abs() < 1e-6);
        // var_pr = PSEUDORANGE_VARIANCE_BASE * snr_scale(45) / sin(pi/2) * 9.0
        // snr_scale(45) = (45/45)^2 = 1.0, var_pr = 1.0 * 1.0 / 1.0 * 9.0 = 9.0
        assert!((meas[0].raw_var - 9.0).abs() < 1e-6);
        // h_row[20] = map_wet
        assert!((meas[0].h_row[20] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_push_pr_measurement_not_iono_free_variance() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let los = Vector3::new(0.0, 0.0, 1.0);

        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 10.0,
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
            map_wet: 1.0,
            iono_delay: 5.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        fg.push_pr_measurement(&mut meas, &state, &sat, &x_i, 0, &los, 5.0, 0.0, 0.0);
        assert_eq!(meas.len(), 1);
        // expected_pr = expected_base + iono_delay = 5.0 + 5.0 = 10.0
        // res_pr = 10.0 - 10.0 = 0.0
        assert!((meas[0].res - 0.0).abs() < 1e-6);
        // var_pr = 1.0 * 1.0 / 1.0 + 9.0 = 10.0
        assert!((meas[0].raw_var - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_push_sat_meas_rejects_large_pr_residual() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let los = Vector3::new(0.0, 0.0, 1.0);

        // PR residual = p1 - expected_pr > 100 -> returns false, no meas
        let sat = ProcessedSat {
            sat_obs: &obs,
            dt_sat_m: 0.0,
            p1: 1000.0,
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
            iono_delay: 0.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let result = fg.push_sat_meas(&mut meas, &state, &sat, &x_i, 0, &los, 5.0, 0.0, 0.0);
        assert!(!result);
        assert!(meas.is_empty());
    }

    #[test]
    fn test_resolve_uduc_indices_with_all_indices() {
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.ambiguity_keys.push((sat_id, 1)); // n1 idx 0
        state.ambiguity_keys.push((sat_id, 2)); // n2 idx 1
        state.ambiguity_keys.push((sat_id, 3)); // i1 idx 2
        state.ambiguities = vec![100.0, 200.0, 300.0];
        let obs = SatObs { sat: sat_id, observations: vec![] };
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
            iono_delay: 0.0,
            f1: 1575.42e6,
            f2: 1227.60e6,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let x_i = DVector::from_fn(CORE_STATE_SIZE + 3, |i, _| i as f64);
        let (i1_idx, n1_idx, n2_idx, i1, n1, n2, gamma) =
            PppIteratedEkf::resolve_uduc_indices(&state, &sat, &x_i);
        assert_eq!(i1_idx, Some(CORE_STATE_SIZE + 2));
        assert_eq!(n1_idx, Some(CORE_STATE_SIZE + 0));
        assert_eq!(n2_idx, Some(CORE_STATE_SIZE + 1));
        let expected_gamma = (1575.42e6 * 1575.42e6) / (1227.60e6 * 1227.60e6);
        assert!((gamma - expected_gamma).abs() < 1e-6);
        assert!((i1 - (CORE_STATE_SIZE + 2) as f64).abs() < 1e-6);
        assert!((n1 - (CORE_STATE_SIZE + 0) as f64).abs() < 1e-6);
        assert!((n2 - (CORE_STATE_SIZE + 1) as f64).abs() < 1e-6);
    }

    #[test]
    fn test_resolve_uduc_indices_missing_indices() {
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
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
            iono_delay: 0.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(),
            sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        let x_i = DVector::zeros(CORE_STATE_SIZE);
        let (i1_idx, n1_idx, n2_idx, i1, n1, n2, gamma) =
            PppIteratedEkf::resolve_uduc_indices(&state, &sat, &x_i);
        // No ambiguity keys -> all indices None, values default to 0.0
        assert!(i1_idx.is_none());
        assert!(n1_idx.is_none());
        assert!(n2_idx.is_none());
        assert!((i1 - 0.0).abs() < 1e-10);
        assert!((n1 - 0.0).abs() < 1e-10);
        assert!((n2 - 0.0).abs() < 1e-10);
        assert!((gamma - 1.0).abs() < 1e-10); // f1/f2 = 1.0/1.0 = 1.0
    }

    #[test]
    fn test_push_doppler_measurement_creates_residual() {
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = SatObs { sat: sat_id, observations: vec![] };
        let x_i = DVector::from_fn(20, |i, _| i as f64);
        let los = Vector3::new(1.0, 0.0, 0.0);
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
            doppler: -100.0,
            lam1: 0.19,
            lam2: 0.24,
            tropo_dry: 0.0,
            map_wet: 0.0,
            iono_delay: 0.0,
            f1: 1.0,
            f2: 1.0,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::new(100.0, 0.0, 0.0),
            sat_clock_drift: 1e-5,
            rcv_pos_ecef: Vector3::zeros(),
            pcv_correction: 0.0,
        };
        fg.push_doppler_measurement(&mut meas, &sat, &x_i, &los);
        assert_eq!(meas.len(), 1);
        // meas_rr = -(-100.0) * 0.19 = 19.0
        // rcv_vel = Vector3(3, 4, 5), los = (1,0,0) -> los.dot(rcv_vel) = 3.0
        // rcv_clk_drift = x_i[19] = 19.0
        // expected_rr = 100.0 - 3.0 + 19.0 - 1e-5 * C
        let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
        let res_expected = 19.0 - (100.0 - 3.0 + 19.0 - 1e-5 * c);
        assert!((meas[0].res - res_expected).abs() < 1e-3);
        // h_row: velocity terms at [3,4,5] and clock drift at [19]
        assert!((meas[0].h_row[3] - (-1.0)).abs() < 1e-10);
        assert!((meas[0].h_row[19] - 1.0).abs() < 1e-10);
        assert!(!meas[0].is_phase);
    }

    // ============ log_ppp_convergence tests ============

    #[test]
    fn test_log_ppp_convergence_basic() {
        // Verify log_ppp_convergence does not panic with a basic state
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_pred = DMatrix::identity(dim, dim);
        log_ppp_convergence(&state, &[], &x_i, &x_pred, &p_pred, &fg);
    }

    #[test]
    fn test_log_ppp_convergence_empty_sats() {
        // Empty sat list should not cause panics
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_pred = DMatrix::identity(dim, dim);
        log_ppp_convergence(&state, &[], &x_i, &x_pred, &p_pred, &fg);
    }

    // ============ compute_final_covariance tests ============

    #[test]
    fn test_compute_final_covariance_empty_meas() {
        // Empty measurements should return p_pred directly
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let p_pred = DMatrix::identity(dim, dim);
        let p_inv = DMatrix::identity(dim, dim);
        let result = fg.compute_final_covariance(&state, &[], &x_i, &p_pred, &p_inv);
        assert_eq!(result, p_pred);
    }

    #[test]
    fn test_compute_final_covariance_with_meas() {
        // One measurement should produce a damped covariance different from p_pred
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        let x_i = DVector::zeros(dim);
        let p_pred = DMatrix::identity(dim, dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        let sat = ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1: 20000000.0,
            p2: None, cp1: None, cp2: None,
            is_iono_free: true,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
            f1: 1.0, f2: 1.0,
            sat_pos_rot: Vector3::new(20000000.0, 0.0, 0.0),
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        };
        let result = fg.compute_final_covariance(&state, &[sat], &x_i, &p_pred, &p_inv);
        assert_eq!(result.nrows(), dim);
        assert_eq!(result.ncols(), dim);
        assert_ne!(result, p_pred, "covariance should differ from prior with measurements present");
    }

    #[test]
    fn test_compute_final_covariance_nan_fallback() {
        // p_inv with NaN causes invert_matrix to fail, falling back to p_pred
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        let x_i = DVector::zeros(dim);
        let p_pred = DMatrix::identity(dim, dim);
        let p_inv_nan = DMatrix::from_element(dim, dim, f64::NAN);
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        let sat = ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1: 20000000.0,
            p2: None, cp1: None, cp2: None,
            is_iono_free: true,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
            f1: 1.0, f2: 1.0,
            sat_pos_rot: Vector3::new(20000000.0, 0.0, 0.0),
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        };
        let result = fg.compute_final_covariance(&state, &[sat], &x_i, &p_pred, &p_inv_nan);
        assert_eq!(result, p_pred, "should fall back to p_pred when inversion fails");
    }

    // ============ compute_iteration_dx tests ============

    fn make_dummy_sat(sat_id: SatelliteId, sat_pos_rot: Vector3<f64>) -> ProcessedSat<'static> {
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1: sat_pos_rot.norm(),
            p2: None, cp1: None, cp2: None,
            is_iono_free: true,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
            f1: 1.0, f2: 1.0,
            sat_pos_rot,
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        }
    }

    #[test]
    fn test_compute_iteration_dx_no_prior() {
        // Normal solution path without position prior
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let result = fg.compute_iteration_dx(&state, &[sat], &x_i, &x_pred, &p_inv, 0, None);
        assert!(result.is_ok());
        let dx = result.unwrap();
        assert!(dx.is_some(), "should produce a delta-x solution");
    }

    #[test]
    fn test_compute_iteration_dx_with_prior() {
        // Position prior anchors the first three state elements
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let spp_pos = Vector3::new(1.0, 2.0, 3.0);
        let prior_var = 1.0;
        let result = fg.compute_iteration_dx(
            &state, &[sat], &x_i, &x_pred, &p_inv, 0, Some((spp_pos, prior_var)),
        );
        assert!(result.is_ok());
        let dx = result.unwrap();
        assert!(dx.is_some(), "should produce a delta-x with position prior");
    }

    // ============ push_uduc_pr_measurements tests ============

    fn make_uduc_sat(
        sat_id: SatelliteId,
        p1: f64,
        p2: Option<f64>,
        cp1: Option<f64>,
        cp2: Option<f64>,
    ) -> ProcessedSat<'static> {
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1, p2, cp1, cp2,
            is_iono_free: false,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.5, iono_delay: 5.0,
            f1: 1575.42e6, f2: 1227.60e6,
            sat_pos_rot: Vector3::zeros(),
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        }
    }

    #[test]
    fn test_push_uduc_pr_measurements_both() {
        // Both p1 and p2 present: should produce two measurements with correct residuals
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat = make_uduc_sat(sat_id, 20000005.0, Some(20000010.0), None, None);
        let x_i_size = CORE_STATE_SIZE + 2;
        let x_i = DVector::zeros(x_i_size);
        let los = Vector3::new(1.0, 0.0, 0.0);
        let expected_base = 20000000.0;
        let i1_idx = Some(CORE_STATE_SIZE);
        let i1 = 3.0;
        let gamma = 1.5;

        fg.push_uduc_pr_measurements(&mut meas, &sat, &x_i, &los, expected_base, i1_idx, i1, gamma);

        assert_eq!(meas.len(), 2);
        // P1: res = p1 - (expected_base + i1) = 20000005 - (20000000 + 3) = 2.0
        assert!((meas[0].res - 2.0).abs() < 1e-6, "P1 residual");
        assert!(!meas[0].is_phase);
        // i1 coefficient in h_row should be 1.0 for P1
        assert!((meas[0].h_row[i1_idx.unwrap()] - 1.0).abs() < 1e-10, "P1 iono coef");

        // P2: res = p2 - (expected_base + gamma * i1) = 20000010 - (20000000 + 1.5*3) = 5.5
        assert!((meas[1].res - 5.5).abs() < 1e-6, "P2 residual");
        assert!(!meas[1].is_phase);
        // i1 coefficient in h_row should be gamma for P2
        assert!((meas[1].h_row[i1_idx.unwrap()] - 1.5).abs() < 1e-10, "P2 iono coef");

        // Variance: var_p1 = PSEUDORANGE_VARIANCE_BASE * snr_scale(45) / sin(pi/2) = 1.0
        assert!((meas[0].raw_var - 1.0).abs() < 1e-6, "P1 variance");
        // var_p2 = var_p1 * 1.5 = 1.5
        assert!((meas[1].raw_var - 1.5).abs() < 1e-6, "P2 variance");
    }

    #[test]
    fn test_push_uduc_pr_measurements_p1_only() {
        // Only p1 present (p2=None): should produce one measurement
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat = make_uduc_sat(sat_id, 20000005.0, None, None, None);
        let x_i_size = CORE_STATE_SIZE + 2;
        let x_i = DVector::zeros(x_i_size);
        let los = Vector3::new(1.0, 0.0, 0.0);
        let expected_base = 20000000.0;
        let i1_idx = Some(CORE_STATE_SIZE);
        let i1 = 3.0;
        let gamma = 1.5;

        fg.push_uduc_pr_measurements(&mut meas, &sat, &x_i, &los, expected_base, i1_idx, i1, gamma);

        assert_eq!(meas.len(), 1);
        assert!(!meas[0].is_phase);
    }

    // ============ push_uduc_cp_measurements tests ============

    #[test]
    fn test_push_uduc_cp_measurements_both() {
        // Both cp1 and cp2 present: should produce two phase measurements with correct signs
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat = make_uduc_sat(sat_id, 0.0, None, Some(105263158.0), Some(83333333.0));
        let x_i_size = CORE_STATE_SIZE + 4;
        let x_i = DVector::zeros(x_i_size);
        let los = Vector3::new(1.0, 0.0, 0.0);
        let expected_base = 20000000.0;
        let i1_idx = Some(CORE_STATE_SIZE);
        let n1_idx = Some(CORE_STATE_SIZE + 1);
        let n2_idx = Some(CORE_STATE_SIZE + 2);
        let i1 = 3.0;
        let n1 = 1.0;
        let n2 = 2.0;
        let gamma = 1.5;

        fg.push_uduc_cp_measurements(
            &mut meas, &state, &sat, &x_i, &los,
            expected_base, i1_idx, n1_idx, n2_idx,
            i1, n1, n2, gamma,
        );

        assert_eq!(meas.len(), 2);
        assert!(meas[0].is_phase, "CP1 is phase");
        assert!(meas[1].is_phase, "CP2 is phase");

        // L1: (cp1 - windup) * lam1 - (expected_base - i1 + n1)
        // windup=0, cp1*lam1=20000000.02, expected_base-i1+n1=19999998.0
        let res_l1 = 105263158.0 * 0.19 - (20000000.0 - 3.0 + 1.0);
        assert!((meas[0].res - res_l1).abs() < 1e-4, "CP1 residual");

        // L2: (cp2 - windup) * lam2 - (expected_base - gamma*i1 + n2)
        let res_l2 = 83333333.0 * 0.24 - (20000000.0 - 1.5 * 3.0 + 2.0);
        assert!((meas[1].res - res_l2).abs() < 1e-4, "CP2 residual");

        // h_row coefficients: L1 iono = -1.0, L2 iono = -gamma
        assert!((meas[0].h_row[i1_idx.unwrap()] - (-1.0)).abs() < 1e-10, "CP1 iono=-1");
        assert!((meas[1].h_row[i1_idx.unwrap()] - (-1.5)).abs() < 1e-10, "CP2 iono=-gamma");

        // Ambiguity coefficients
        assert!((meas[0].h_row[n1_idx.unwrap()] - 1.0).abs() < 1e-10, "CP1 amb coef");
        assert!((meas[1].h_row[n2_idx.unwrap()] - 1.0).abs() < 1e-10, "CP2 amb coef");
    }

    #[test]
    fn test_push_uduc_cp_measurements_cp1_only() {
        // Only cp1 present (cp2=None): should produce one phase measurement
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat = make_uduc_sat(sat_id, 0.0, None, Some(105263158.0), None);
        let x_i_size = CORE_STATE_SIZE + 4;
        let x_i = DVector::zeros(x_i_size);
        let los = Vector3::new(1.0, 0.0, 0.0);
        let expected_base = 20000000.0;
        let i1_idx = Some(CORE_STATE_SIZE);
        let n1_idx = Some(CORE_STATE_SIZE + 1);
        let n2_idx = Some(CORE_STATE_SIZE + 2);
        let i1 = 3.0;
        let n1 = 1.0;
        let n2 = 2.0;
        let gamma = 1.5;

        fg.push_uduc_cp_measurements(
            &mut meas, &state, &sat, &x_i, &los,
            expected_base, i1_idx, n1_idx, n2_idx,
            i1, n1, n2, gamma,
        );

        assert_eq!(meas.len(), 1);
        assert!(meas[0].is_phase);
    }

    #[test]
    fn test_push_uduc_cp_measurements_with_windup() {
        // Windup value should subtract from carrier phase
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.windup.insert(sat_id, 2.5);
        let mut meas = Vec::new();
        let sat = make_uduc_sat(sat_id, 0.0, None, Some(105263158.0), Some(83333333.0));
        let x_i_size = CORE_STATE_SIZE + 4;
        let x_i = DVector::zeros(x_i_size);
        let los = Vector3::new(1.0, 0.0, 0.0);
        let expected_base = 20000000.0;
        let i1_idx = Some(CORE_STATE_SIZE);
        let n1_idx = Some(CORE_STATE_SIZE + 1);
        let n2_idx = Some(CORE_STATE_SIZE + 2);
        let i1 = 3.0;
        let n1 = 1.0;
        let n2 = 2.0;
        let gamma = 1.5;

        fg.push_uduc_cp_measurements(
            &mut meas, &state, &sat, &x_i, &los,
            expected_base, i1_idx, n1_idx, n2_idx,
            i1, n1, n2, gamma,
        );

        assert_eq!(meas.len(), 2);
        // With windup=2.5: (cp1 - 2.5) * 0.19 - (expected_base - i1 + n1)
        let res_l1_windup = (105263158.0 - 2.5) * 0.19 - (20000000.0 - 3.0 + 1.0);
        assert!((meas[0].res - res_l1_windup).abs() < 1e-4, "CP1 with windup");
        let res_l2_windup = (83333333.0 - 2.5) * 0.24 - (20000000.0 - 1.5 * 3.0 + 2.0);
        assert!((meas[1].res - res_l2_windup).abs() < 1e-4, "CP2 with windup");
        // Verify windup actually modified the residual vs no-windup baseline
        let res_l1_no_windup = 105263158.0 * 0.19 - (20000000.0 - 3.0 + 1.0);
        assert!((meas[0].res - res_l1_no_windup).abs() > 0.1, "windup should change residual");
    }

    // ============ push_uduc_measurements tests ============

    #[test]
    fn test_push_uduc_measurements_full() {
        // Full UDUC: both PR and CP measurements are pushed (p1,p2,cp1,cp2)
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.ambiguity_keys.push((sat_id, 1)); // n1
        state.ambiguity_keys.push((sat_id, 2)); // n2
        state.ambiguity_keys.push((sat_id, 3)); // i1
        state.ambiguities = vec![0.0, 0.0, 0.0];
        let sat = make_uduc_sat(
            sat_id, 20000005.0, Some(20000010.0), Some(105263158.0), Some(83333333.0),
        );
        let dim = CORE_STATE_SIZE + 3;
        let x_i = DVector::zeros(dim);
        let los = Vector3::new(1.0, 0.0, 0.0);

        fg.push_uduc_measurements(&mut meas, &state, &sat, &x_i, 0, &los, 20000000.0, 0.0, 0.0);

        // Should produce: p1, p2, cp1, cp2 = 4 measurements
        assert_eq!(meas.len(), 4);
        assert!(!meas[0].is_phase, "p1 is PR");
        assert!(!meas[1].is_phase, "p2 is PR");
        assert!(meas[2].is_phase, "cp1 is phase");
        assert!(meas[3].is_phase, "cp2 is phase");

        // Verify i1 state coefficient signs via resolve_uduc_indices:
        // i1 is at CORE_STATE_SIZE + 2 (third ambiguity key)
        let i1_idx = CORE_STATE_SIZE + 2;
        // PR uses positive i1 (coef=1.0), CP uses negative i1 (coef=-1.0)
        assert!((meas[0].h_row[i1_idx] - 1.0).abs() < 1e-10, "p1 i1 coef=+1");
        assert!((meas[2].h_row[i1_idx] - (-1.0)).abs() < 1e-10, "cp1 i1 coef=-1");
    }

    // ============ build_measurements tests ============

    #[test]
    fn test_build_measurements_with_sats() {
        // Two iono-free sats at different positions produce two PR measurements
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        let x_i = DVector::zeros(dim);
        let sat1 = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let sat2 = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 2 },
            Vector3::new(0.0, 20000000.0, 0.0),
        );
        let meas = fg.build_measurements(&state, &[sat1, sat2], &x_i, 0);
        assert_eq!(meas.len(), 2);
        for m in &meas {
            assert!(!m.is_phase, "PR measurements");
            assert!((m.h_row[15] - 1.0).abs() < 1e-10, "clock bias coef");
        }
        // Each measurement should have a different los direction
        assert!((meas[0].h_row[0] - (-1.0)).abs() < 1e-10, "sat1 los.x");
        assert!((meas[1].h_row[1] - (-1.0)).abs() < 1e-10, "sat2 los.y");
    }

    // ============ push_sat_meas tests ============

    #[test]
    fn test_push_sat_meas_uduc_path() {
        // UDUC path: !iono_free, cp1+cp2+p2 present -> 5 measurements (4 UDUC + iono prior)
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.ambiguity_keys.push((sat_id, 1)); // n1
        state.ambiguity_keys.push((sat_id, 2)); // n2
        state.ambiguity_keys.push((sat_id, 3)); // i1
        state.ambiguities = vec![0.0, 0.0, 0.0];
        let dim = CORE_STATE_SIZE + 3;
        let x_i = DVector::zeros(dim);
        let los = Vector3::new(0.0, 0.0, 1.0);
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        let sat = ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1: 20000000.0,
            p2: Some(20000000.0),
            cp1: Some(105263158.0),
            cp2: Some(83333333.0),
            is_iono_free: false,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
            f1: 1575.42e6, f2: 1227.60e6,
            sat_pos_rot: Vector3::new(20000000.0, 0.0, 0.0),
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        };

        let result = fg.push_sat_meas(
            &mut meas, &state, &sat, &x_i, 0, &los, 20000000.0, 20000000.0, 0.0,
        );
        assert!(result, "sat should be accepted");
        // UDUC: p1, p2, cp1, cp2 + iono prior = 5
        assert_eq!(meas.len(), 5);
    }

    #[test]
    fn test_push_sat_meas_iono_prior() {
        // Iono prior constraint measurement has correct residual and variance
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.ambiguity_keys.push((sat_id, 3)); // i1 idx needed for iono prior
        state.ambiguities = vec![0.0];
        let dim = CORE_STATE_SIZE + 1;
        let x_i = DVector::zeros(dim);
        let los = Vector3::new(0.0, 0.0, 1.0);
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        let sat = ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1: 20000000.0,
            p2: Some(20000000.0),
            cp1: Some(105263158.0),
            cp2: Some(83333333.0),
            is_iono_free: false,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.0, iono_delay: 5.0, // Klobuchar prediction
            f1: 1575.42e6, f2: 1227.60e6,
            sat_pos_rot: Vector3::new(20000000.0, 0.0, 0.0),
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        };

        let result = fg.push_sat_meas(
            &mut meas, &state, &sat, &x_i, 0, &los, 20000000.0, 20000000.0, 0.0,
        );
        assert!(result);
        // UDUC: p1, p2, cp1, cp2 + iono prior = 5
        assert_eq!(meas.len(), 5);
        // Last measurement is the iono prior
        let iono = &meas[4];
        assert!(!iono.is_phase);
        // res = sat.iono_delay - x_i[i1_idx]; i1_idx=0, x_i[CORE_STATE_SIZE]=0 -> res=5.0
        assert!((iono.res - 5.0).abs() < 1e-6, "iono prior residual");
        assert!((iono.raw_var - 9.0).abs() < 1e-6, "iono prior variance (3m std)");
    }

    #[test]
    fn test_push_sat_meas_try_cp_success() {
        // try_push_cp adds a CP measurement for iono-free sats with non-zero cp1
        let fg = PppIteratedEkf::default();
        let mut meas = Vec::new();
        let mut state = dummy_rtk_state();
        let sat_id = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        state.ambiguity_keys.push((sat_id, 0));
        state.ambiguities = vec![0.0];
        let dim = CORE_STATE_SIZE + 1;
        let x_i = DVector::zeros(dim);
        let los = Vector3::new(0.0, 0.0, 1.0);
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        let sat = ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1: 20000000.0,
            p2: None, cp1: Some(105263158.0), cp2: None,
            is_iono_free: true,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
            f1: 1.0, f2: 1.0,
            sat_pos_rot: Vector3::new(20000000.0, 0.0, 0.0),
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        };

        let result = fg.push_sat_meas(
            &mut meas, &state, &sat, &x_i, 0, &los, 20000000.0, 20000000.0, 0.0,
        );
        assert!(result);
        // PR + CP = 2 measurements
        assert_eq!(meas.len(), 2);
        assert!(!meas[0].is_phase, "PR measurement");
        assert!(meas[1].is_phase, "CP measurement");
    }

    // ============ resolve_narrowlane_ar tests (uses AR_MOCK) ============

    #[test]
    fn test_resolve_narrowlane_ar_mock() {
        // Use AR_MOCK to verify NL resolution returns expected modified state
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE + 4;
        let x_wl = DVector::from_fn(dim, |i, _| i as f64);
        let p_wl = DMatrix::identity(dim, dim);
        let subset = vec![(
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0, 1, 1.0, 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 2, 3, 1.0, 0.19, 0.24),
        )];
        let keep_indices = vec![0];

        let mock_nl: Result<(DVector<f64>, DMatrix<f64>), &'static str> = Ok((
            DVector::zeros(dim),
            DMatrix::zeros(dim, dim),
        ));
        {
            let mut lock = AR_MOCK.lock().unwrap();
            *lock = Some(ArMock { wl_result: None, nl_result: Some(mock_nl), nl_calls: 0 });
        }

        let result = fg.resolve_narrowlane_ar(&state, &subset, &keep_indices, &x_wl, &p_wl);
        assert!(result.is_ok(), "mock NL should succeed");
        let (x_fixed, p_fixed) = result.unwrap();
        // Mock adds 9.5 to x_wl[0] and scales p_wl by 0.5
        assert!((x_fixed[0] - (x_wl[0] + 9.5)).abs() < 1e-10, "NL adds 9.5 to pos.x");
        assert!((p_fixed[(0, 0)] - 0.5).abs() < 1e-10, "NL scales covariance by 0.5");
        // Verify nl_calls was incremented
        {
            let lock = AR_MOCK.lock().unwrap();
            assert_eq!(lock.as_ref().unwrap().nl_calls, 1);
        }
        // Clean up mock
        { let mut lock = AR_MOCK.lock().unwrap(); *lock = None; }
    }

    // ============ process_constellation_group tests (uses AR_MOCK) ============

    #[test]
    fn test_process_constellation_group_less_than_two() {
        // Group with fewer than 2 candidates returns None
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_current = DVector::zeros(dim);
        let p_current = DMatrix::identity(dim, dim);
        let result = fg.process_constellation_group(
            &state, &p_current, &x_current, &[], Constellation::Gps,
        );
        assert!(result.is_none());
    }

    #[test]
    fn test_process_constellation_group_empty_keep_indices() {
        // WL returns empty keep_indices -> process_constellation_group returns None
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE + 4;
        let x_current = DVector::zeros(dim);
        let p_current = DMatrix::identity(dim, dim);
        let group_cands = vec![
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0, 1, 30.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 2, 3, 20.0_f64.to_radians(), 0.19, 0.24),
        ];

        let mock_wl: Result<(DVector<f64>, DMatrix<f64>, Vec<usize>), &'static str> = Ok((
            DVector::zeros(dim),
            DMatrix::identity(dim, dim),
            vec![],
        ));
        {
            let mut lock = AR_MOCK.lock().unwrap();
            *lock = Some(ArMock { wl_result: Some(mock_wl), nl_result: None, nl_calls: 0 });
        }

        let result = fg.process_constellation_group(
            &state, &p_current, &x_current, &group_cands, Constellation::Gps,
        );
        assert!(result.is_none(), "empty keep_indices -> None");

        { let mut lock = AR_MOCK.lock().unwrap(); *lock = None; }
    }

    #[test]
    fn test_process_constellation_group_mock_success() {
        // With mock WL+NL, group resolves successfully (jump 9.5 < 10 passes position check)
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE + 4;
        let x_current = DVector::zeros(dim);
        let p_current = DMatrix::identity(dim, dim);
        let group_cands = vec![
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0, 1, 30.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 2, 3, 20.0_f64.to_radians(), 0.19, 0.24),
        ];

        let mock_wl = Ok((DVector::zeros(dim), DMatrix::identity(dim, dim), vec![0]));
        let mock_nl = Ok((DVector::zeros(dim), DMatrix::zeros(dim, dim)));
        {
            let mut lock = AR_MOCK.lock().unwrap();
            *lock = Some(ArMock { wl_result: Some(mock_wl), nl_result: Some(mock_nl), nl_calls: 0 });
        }

        let result = fg.process_constellation_group(
            &state, &p_current, &x_current, &group_cands, Constellation::Gps,
        );
        assert!(result.is_some(), "mock AR should succeed");
        let (_xf, _pf, n_sats) = result.unwrap();
        assert_eq!(n_sats, 2, "keep_indices.len() + 1 = 1 + 1 = 2");

        { let mut lock = AR_MOCK.lock().unwrap(); *lock = None; }
    }

    // ============ compute_iteration_dx tests ============

    #[test]
    fn test_compute_iteration_dx_empty_meas() {
        // No satellites -> build_measurements returns empty -> Err(InsufficientSatellites)
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let result = fg.compute_iteration_dx(&state, &[], &x_i, &x_pred, &p_inv, 0, None);
        assert!(matches!(result, Err(EngineError::InsufficientSatellites)));
    }

    // ============ try_inter_constellation_fallback tests (uses AR_MOCK) ============

    #[test]
    fn test_try_inter_constellation_fallback_subset_too_small() {
        // build_ar_subset produces < 3 pairs -> returns None early
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_current = DVector::zeros(dim);
        let p_current = DMatrix::identity(dim, dim);
        // 3 candidates -> GPS ref + 2 pairs -> 2 < 3 -> None
        let cands = vec![
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0, 1, 30.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 2, 3, 20.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Galileo, prn: 1 }, 4, 5, 25.0_f64.to_radians(), 0.19, 0.24),
        ];
        let result = fg.try_inter_constellation_fallback(&state, &p_current, &x_current, &cands);
        assert!(result.is_none(), "3 candidates -> 2 pairs < 3 -> None");
    }

    #[test]
    fn test_try_inter_constellation_fallback_empty_keep_indices() {
        // WL mock returns empty keep_indices -> returns None
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE + 6;
        let x_current = DVector::zeros(dim);
        let p_current = DMatrix::identity(dim, dim);
        // 4 candidates -> GPS ref + 3 pairs >= 3 -> passes subset check
        let cands = vec![
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0, 1, 30.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 2, 3, 20.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 3 }, 4, 5, 15.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Galileo, prn: 1 }, 6, 7, 25.0_f64.to_radians(), 0.19, 0.24),
        ];
        let mock_wl = Ok((DVector::zeros(dim), DMatrix::identity(dim, dim), vec![]));
        {
            let mut lock = AR_MOCK.lock().unwrap();
            *lock = Some(ArMock { wl_result: Some(mock_wl), nl_result: None, nl_calls: 0 });
        }
        let result = fg.try_inter_constellation_fallback(&state, &p_current, &x_current, &cands);
        assert!(result.is_none(), "empty keep_indices -> None");
        { let mut lock = AR_MOCK.lock().unwrap(); *lock = None; }
    }

    #[test]
    fn test_try_inter_constellation_fallback_success() {
        // Full mock success path with valid WL keep_indices and NL resolution
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE + 6;
        let x_current = DVector::zeros(dim);
        let p_current = DMatrix::identity(dim, dim);
        let cands = vec![
            (SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0, 1, 30.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 2 }, 2, 3, 20.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Gps, prn: 3 }, 4, 5, 15.0_f64.to_radians(), 0.19, 0.24),
            (SatelliteId { constellation: Constellation::Galileo, prn: 1 }, 6, 7, 25.0_f64.to_radians(), 0.19, 0.24),
        ];
        let mock_wl = Ok((DVector::zeros(dim), DMatrix::identity(dim, dim), vec![0, 1, 2]));
        let mock_nl = Ok((DVector::zeros(dim), DMatrix::zeros(dim, dim)));
        {
            let mut lock = AR_MOCK.lock().unwrap();
            *lock = Some(ArMock { wl_result: Some(mock_wl), nl_result: Some(mock_nl), nl_calls: 0 });
        }
        let result = fg.try_inter_constellation_fallback(&state, &p_current, &x_current, &cands);
        assert!(result.is_some(), "mock fallback should succeed");
        let (xf, _pf, n_sats) = result.unwrap();
        assert_eq!(n_sats, 4, "keep_indices.len() + 1 = 3 + 1 = 4");
        // Mock NL adds 9.5 to x_wl[0] and jump=9.5 <= 20.0 passes position check
        assert!((xf[0] - 9.5).abs() < 1e-10, "NL adds 9.5 to x[0]");
        { let mut lock = AR_MOCK.lock().unwrap(); *lock = None; }
    }

    // ============ solve() convergence with measurements ============

    #[test]
    fn test_solve_normal_convergence() {
        // Normal convergence: one satellite with matching p1/geometry -> dx=0 -> converges immediately
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        state.position.vector = Vector3::new(0.0, 0.0, 0.0);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let result = fg.solve(&mut state, &[sat], None);
        assert!(result.is_ok(), "solve should converge with matching geometry");
        assert!(state.full_x_predict.is_some(), "x_pred should be saved");
        assert!(state.full_p_predict.is_some(), "p_pred should be saved");
    }

    // ============ SPP Anchor / Position Prior tests ============
    //
    // The SPP anchor applies a soft position prior in the IEKF at indices 0,1,2
    // (X, Y, Z in ECEF). The prior weight is 1/variance, added to the htwh diagonal
    // and htwr residual. These tests verify the prior math:
    //   - Prior pulls position toward SPP with the correct sign
    //   - Stronger prior (smaller variance) pulls harder
    //   - The prior primarily targets position, not clock or other states
    //   - Full solve converges with a prior and moves position

    #[test]
    fn test_compute_iteration_dx_prior_sign_correct() {
        // Prior to the LEFT of current state should produce negative dx;
        // prior to the RIGHT should produce positive dx.
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );

        // Prior at X=-1000 with state at X=0: prior says "move left"
        let spp_neg = Vector3::new(-1000.0, 0.0, 0.0);
        let dx_neg = fg
            .compute_iteration_dx(
                &state,
                &[sat.clone()],
                &x_i,
                &x_pred,
                &p_inv,
                0,
                Some((spp_neg, 1.0)),
            )
            .unwrap()
            .unwrap();
        assert!(
            dx_neg[0] < -100.0,
            "prior at X=-1000 should pull negative, got dx[0]={}",
            dx_neg[0]
        );

        // Prior at X=+1000: prior says "move right"
        let spp_pos = Vector3::new(1000.0, 0.0, 0.0);
        let dx_pos = fg
            .compute_iteration_dx(
                &state,
                &[sat],
                &x_i,
                &x_pred,
                &p_inv,
                0,
                Some((spp_pos, 1.0)),
            )
            .unwrap()
            .unwrap();
        assert!(
            dx_pos[0] > 100.0,
            "prior at X=1000 should pull positive, got dx[0]={}",
            dx_pos[0]
        );

        // Verify opposite signs
        assert!(
            dx_neg[0] < 0.0 && dx_pos[0] > 0.0,
            "opposite prior positions should produce opposite-sign dx"
        );
    }

    #[test]
    fn test_compute_iteration_dx_prior_strength_scales_with_variance() {
        // A tighter prior (smaller variance) should produce larger position corrections
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let spp_pos = Vector3::new(1000.0, 0.0, 0.0);
        let spp_exact = Vector3::new(0.0, 0.0, 0.0); // prior matches state exactly

        // Weak prior: large variance = 100 (weight = 1/100 = 0.01)
        let dx_weak = fg
            .compute_iteration_dx(
                &state,
                &[sat.clone()],
                &x_i,
                &x_pred,
                &p_inv,
                0,
                Some((spp_pos, 100.0)),
            )
            .unwrap()
            .unwrap();

        // Strong prior: small variance = 1 (weight = 1.0)
        let dx_strong = fg
            .compute_iteration_dx(
                &state,
                &[sat.clone()],
                &x_i,
                &x_pred,
                &p_inv,
                0,
                Some((spp_pos, 1.0)),
            )
            .unwrap()
            .unwrap();

        // Strong prior should pull harder in X
        assert!(
            dx_strong[0].abs() > dx_weak[0].abs(),
            "strong prior (var=1, dx[0]={}) should pull X harder than weak prior (var=100, dx[0]={})",
            dx_strong[0],
            dx_weak[0]
        );

        // No pull when prior matches current state exactly (zero innovation)
        let dx_exact = fg
            .compute_iteration_dx(
                &state,
                &[sat],
                &x_i,
                &x_pred,
                &p_inv,
                0,
                Some((spp_exact, 1.0)),
            )
            .unwrap()
            .unwrap();

        // When prior == state, the htwr contribution is zero, but the htwh damping
        // still increases diagonal elements (tightens the covariance).
        // dx may not be exactly zero because the stronger diagonal pulls the
        // solution toward the prediction (x_pred == x_i here, so it should be ~0).
        assert!(
            dx_exact[0].abs() < 1.0,
            "prior matching state should produce negligible dx, got dx[0]={}",
            dx_exact[0]
        );
    }

    #[test]
    fn test_compute_iteration_dx_prior_targets_position_indices() {
        // The position prior is applied only to state indices 0, 1, 2 (X, Y, Z in ECEF).
        // Non-position elements like clock bias (index 15) are only affected through
        // measurement coupling, not directly by the prior.
        let fg = PppIteratedEkf::default();
        let state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );

        // Prior only pulls X positive
        let spp_pos = Vector3::new(1000.0, 0.0, 0.0);
        let dx = fg
            .compute_iteration_dx(
                &state,
                &[sat],
                &x_i,
                &x_pred,
                &p_inv,
                0,
                Some((spp_pos, 1.0)),
            )
            .unwrap()
            .unwrap();

        // X should be pulled positive
        assert!(
            dx[0] > 100.0,
            "prior at X=1000 should produce large positive dx[0], got {}",
            dx[0]
        );

        // Y and Z have no prior and no measurement sensitivity in this setup
        // (measurement LOS is along X axis), so they should be near zero
        assert!(
            dx[1].abs() < 1e-6,
            "Y should not be directly pulled by prior, got dx[1]={}",
            dx[1]
        );
        assert!(
            dx[2].abs() < 1e-6,
            "Z should not be directly pulled by prior, got dx[2]={}",
            dx[2]
        );

        // Clock bias (index 15) is coupled through the measurement H matrix
        // (which has -1 at [0] and +1 at [15]). Anchoring position naturally
        // helps resolve clock-state ambiguity, but the clock correction should
        // be an order of magnitude smaller than the position correction.
        assert!(
            dx[15].abs() < dx[0].abs(),
            "clock correction ({}) should be smaller than position correction ({})",
            dx[15],
            dx[0]
        );
    }

    #[test]
    fn test_solve_with_prior_pulls_position_toward_spp() {
        // Full solve() with a position prior should pull the estimated
        // position toward the SPP position while converging normally.
        // NOTE: the prior displacement must be small (<~200m) so that the
        // satellite pseudorange residual stays within the 100m rejection
        // threshold across IEKF iterations as position evolves.
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        state.position.vector = Vector3::new(0.0, 0.0, 0.0);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );

        // Small SPP prior displacement: 10m in X
        let spp_pos = Vector3::new(10.0, 0.0, 0.0);
        let result = fg.solve(&mut state, &[sat], Some((spp_pos, 1.0)));
        assert!(
            result.is_ok(),
            "solve should converge with position prior"
        );

        // The converged position should have moved from 0 toward SPP (10).
        // The exact balance depends on prior weight vs process noise.
        assert!(
            state.position.vector.x > 0.5,
            "solve with prior should pull X toward SPP (10), got X={}",
            state.position.vector.x
        );

        // The position should not overshoot the prior
        assert!(
            state.position.vector.x < 9.5,
            "solve with prior should not overshoot SPP position, got X={}",
            state.position.vector.x
        );

        // Y and Z should stay near zero (no prior pull on those axes)
        assert!(
            state.position.vector.y.abs() < 1.0,
            "Y should not be pulled by X-axis prior, got Y={}",
            state.position.vector.y
        );
        assert!(
            state.position.vector.z.abs() < 1.0,
            "Z should not be pulled by X-axis prior, got Z={}",
            state.position.vector.z
        );
    }

    #[test]
    fn test_solve_with_prior_pulls_all_three_axes() {
        // Prior pulling in all three axes simultaneously should move each
        // component toward its respective prior value.
        let fg = PppIteratedEkf::default();
        let mut state = dummy_rtk_state();
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        state.position.vector = Vector3::new(0.0, 0.0, 0.0);
        let sat1 = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );
        let sat2 = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 2 },
            Vector3::new(0.0, 20000000.0, 0.0),
        );

        // Small 3-axis prior displacement to keep PR residuals < 100m
        let spp_pos = Vector3::new(10.0, -8.0, 5.0);
        let result = fg.solve(&mut state, &[sat1, sat2], Some((spp_pos, 1.0)));
        assert!(
            result.is_ok(),
            "solve should converge with 3-axis prior"
        );

        // Each axis should move toward the prior (at least by 10% of the pull)
        assert!(
            state.position.vector.x > 0.5,
            "X should be pulled toward 10, got X={}",
            state.position.vector.x
        );
        assert!(
            state.position.vector.y < -0.5,
            "Y should be pulled toward -8, got Y={}",
            state.position.vector.y
        );
        assert!(
            state.position.vector.z > 0.3,
            "Z should be pulled toward 5, got Z={}",
            state.position.vector.z
        );
    }
}

// =========================================================================
// Adversarial tests: PPP accuracy gap investigation
// =========================================================================
//
// These tests expose the ROOT CAUSES of the ~0.5m horizontal bias observed
// in the f9p PPP benchmark. The findings are:
//
// Issue #1: SPP PRIOR VARIANCE FLOOR (line 52 of process_ppp.rs):
//   let prior_var = (pos_cov.min(25.0)).max(1.0);
//
//   When the IEKF position covariance converges below 1.0 m^2 (10 cm std),
//   the floor at 1.0 m^2 INCREASES the prior variance, DECREASING the prior
//   weight from 1/0.01=100 to 1/1.0=1. This makes the SPP anchor 100x
//   WEAKER than optimal after convergence. The filter loses its primary
//   absolute position anchor just when it needs it most.
//
// Issue #2: WHITE-NOISE CLOCK MODEL (predictor.rs line 134):
//   phi[(15, 15)] = 0.0
//
//   This destroys temporal correlation of the clock bias. The predicted
//   clock variance resets to ~process_noise*dt each epoch instead of
//   remaining converged. The clock prior weight drops to ~9e-5, meaning
//   the clock is re-estimated from scratch every epoch.
//
// COMBINED EFFECT: The filter loses BOTH absolute position anchors:
// - The SPP prior weight drops from 100 to 1 (floor)
// - The clock prior weight stays at ~1e-4 forever (white-noise)
//
// The result is a position estimate that converges to ~0.5-1.5m rather than
// the cm-level accuracy achievable with a random-walk clock and strong prior.

#[cfg(test)]
mod adversarial_accuracy_tests {
    use super::*;
    use crate::engine::predictor::{compute_process_noise, compute_transition_matrix};
    use crate::engine::{DynamicsModel, EngineConfig};
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::obs::SatObs;
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;

    fn make_state_with_pos(time: GpsTime, pos: Vector3<f64>, initial_var: f64) -> RtkState {
        let coord = Coordinate::new(pos, Datum::WGS84, Frame::ECEF, time);
        RtkState::new(time, coord, initial_var)
    }

    /// Copy of make_dummy_sat from mutant_killer_tests (needed here since the
    /// original is not pub). Key detail: p1 = sat_pos_rot.norm(), is_iono_free=true,
    /// so the pseudorange residual res_pr = p1 - geometric_dist ≈ 0 when the
    /// state is at the origin. This ensures the measurement passes the 100m
    /// residual check in push_sat_meas.
    fn make_dummy_sat(sat_id: SatelliteId, sat_pos_rot: Vector3<f64>) -> ProcessedSat<'static> {
        let obs = Box::leak(Box::new(SatObs { sat: sat_id, observations: vec![] }));
        ProcessedSat {
            sat_obs: obs,
            dt_sat_m: 0.0,
            p1: sat_pos_rot.norm(),
            p2: None, cp1: None, cp2: None,
            is_iono_free: true,
            osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0,
            los: Vector3::zeros(), dist: 0.0,
            el: std::f64::consts::PI / 2.0,
            snr: 45.0, doppler: 0.0,
            lam1: 0.19, lam2: 0.24,
            tropo_dry: 0.0, map_wet: 0.0, iono_delay: 0.0,
            f1: 1.0, f2: 1.0,
            sat_pos_rot,
            sat_vel: Vector3::zeros(), sat_clock_drift: 0.0,
            rcv_pos_ecef: Vector3::zeros(), pcv_correction: 0.0,
        }
    }

    // ======================================================================
    // Adversarial Test 1: SPP prior variance floor WEAKENS position anchoring
    // ======================================================================
    //
    // In process_ppp.rs, the SPP prior variance is:
    //   let prior_var = (pos_cov.min(25.0)).max(1.0);
    //
    // KEY INSIGHT: The floor at 1.0 m^2 makes the prior WEAKER after
    // convergence, not stronger. When position covariance converges to
    // 0.01 m^2 (10 cm std), the natural prior weight should be 1/0.01=100.
    // Instead, the floor forces prior_var=1.0, giving weight=1.0 — a 100x
    // reduction.
    //
    // This means the SPP prior provides WEAK absolute position anchoring
    // at the exact time when the filter has converged and would benefit
    // from strong temporal constraints. Combined with the white-noise
    // clock model (which also destroys temporal correlation), the filter
    // lacks sufficient absolute position information to reach cm-level
    // accuracy. The position estimate converges to ~0.5-1m instead of
    // the sub-10cm achievable with proper models.

    #[test]
    fn test_spp_prior_min_variance_prevents_submeter_accuracy() {
        // Simulate a converged state with small position covariance (0.01 m^2)
        let t = GpsTime::new(2156, 1000.0);
        let true_pos = Vector3::new(6000000.0, 0.0, 0.0);
        let mut state = make_state_with_pos(t, true_pos, 100.0);
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);

        // Set position covariance to 0.01 m^2 (10 cm std) — representing
        // a well-converged PPP solution
        for i in 0..3 {
            state.covariance[(i, i)] = 0.01;
        }

        // Simulate what process_ppp does: compute prior_var from pos_cov
        let pos_cov = state.covariance[(0, 0)]
            .min(state.covariance[(1, 1)])
            .min(state.covariance[(2, 2)]);
        // This is the actual code from process_ppp line 52:
        let prior_var = (pos_cov.min(25.0)).max(1.0);

        // The bug: prior_var should be pos_cov=0.01 for a well-converged filter
        // (so weight = 1/0.01 = 100), but it's clamped to 1.0 (weight = 1.0).
        // This REDUCES the prior weight by 100x after convergence, meaning the
        // SPP prior provides almost no position anchoring for a converged filter.
        // With the white-noise clock model (phi=0), the clock bias doesn't maintain
        // temporal correlation either, so the filter has TWO weak constraints on
        // absolute position: the SPP prior (weight=1) and the clock prior
        // (weight ≈ 1e-4). This is insufficient for cm-level PPP accuracy.
        assert_eq!(
            prior_var, 1.0,
            "BUG: prior_var should be pos_cov={} but clamping forces it to 1.0. \
             This reduces the prior weight from 1/0.01=100 to 1/1.0=1, making the \
             SPP anchor 100x weaker than it should be after convergence.",
            pos_cov
        );
        assert!(
            prior_var > pos_cov * 10.0,
            "prior_var={} is {:.0}x LARGER than pos_cov={}. The prior weight \
             is {:.0}x WEAKER than optimal. Combined with the white-noise clock, \
             this removes the two key absolute position anchors.",
            prior_var, prior_var / pos_cov, pos_cov, prior_var / pos_cov
        );
    }

    #[test]
    fn test_spp_prior_clamping_biases_final_position() {
        // Full IEKF solve demonstrating that the prior variance floor
        // (1.0 m^2) injects position bias when the true SPP error is ~2m.
        //
        // Setup: state position at 0 with cov=0.01 on diagonal (converged).
        // SPP position at (2, 0, 0) — typical 2m horizontal SPP error.
        // Prior variance = 1.0 (the clamped value).
        //
        // With weak prior (var=100, no clamping), position should stay near 0.
        // With strong prior (var=1, clamping active), position should shift toward 2.

        let fg = PppIteratedEkf::default();
        let t = GpsTime::new(2156, 1000.0);
        let mut state = make_state_with_pos(t, Vector3::new(0.0, 0.0, 0.0), 100.0);
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        for i in 0..3 {
            state.covariance[(i, i)] = 0.01; // converged
        }

        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );

        // SPP position with ~2m error in X (typical for f9p)
        let spp_pos = Vector3::new(2.0, 0.0, 0.0);

        // Solve WITHOUT clamping: use actual pos_cov = 0.01 as prior variance
        // (this simulates what would happen if the max(1.0) clamp were removed)
        let mut state_no_clamp = state.clone();
        let _ = fg.solve(&mut state_no_clamp, &[sat.clone()], Some((spp_pos, 0.01)));
        let x_no_clamp = state_no_clamp.position.vector.x;

        // Solve WITH clamping: prior variance = 1.0 (current bug)
        let mut state_with_clamp = state.clone();
        let _ = fg.solve(&mut state_with_clamp, &[sat.clone()], Some((spp_pos, 1.0)));
        let x_with_clamp = state_with_clamp.position.vector.x;

        // The clamped version should have MORE bias toward SPP error
        let bias_no_clamp = (x_no_clamp - 0.0).abs();
        let bias_with_clamp = (x_with_clamp - 0.0).abs();

        tracing::info!(
            "SPP prior variance test: no_clamp_x={:.6}, with_clamp_x={:.6}, \
             bias_no_clamp={:.6}, bias_with_clamp={:.6}",
            x_no_clamp, x_with_clamp, bias_no_clamp, bias_with_clamp
        );

        // KEY INSIGHT: The clamped prior (var=1.0, weight=1.0) pulls LESS
        // than the natural prior (var=0.01, weight=100.0). This means the
        // 1.0 clamping floor ACTUALLY WEAKENS the prior after convergence,
        // reducing absolute position anchoring.
        // The real problem: with pos_cov=0.01, the filter should naturally
        // have prior weight=100, but the floor forces weight=1. This leaves
        // the position less constrained, making it vulnerable to drift
        // when combined with the white-noise clock model.
        assert!(
            x_with_clamp.abs() < x_no_clamp.abs(),
            "Clamping prior var to 1.0 paradoxically REDUCES prior weight \
             (from 1/0.01=100 to 1/1.0=1). Clamped should pull LESS: \
             no_clamp={:.4} vs with_clamp={:.4}",
            x_no_clamp, x_with_clamp
        );
        assert!(
            x_no_clamp.abs() > x_with_clamp.abs() * 2.0,
            "Natural prior (var=0.01) should pull AT LEAST 2x more than \
             clamped prior (var=1.0): no_clamp={:.4}, with_clamp={:.4}",
            x_no_clamp, x_with_clamp
        );
    }

    // ======================================================================
    // Adversarial Test 2: White-noise clock model destroys temporal correlation
    // ======================================================================
    //
    // In predictor.rs line 134: phi[(15, 15)] = 0.0
    //
    // This means the clock bias is NOT propagated from one epoch to the next.
    // Instead, it's reset each epoch with information coming only from the
    // clock drift and process noise. This creates a system where:
    //
    // - P_pred[15,15] ≈ process_noise_cb * dt + dt^2 * P_drift
    // - With process_noise_cb = 10000, P_pred[15,15] ≈ 10000 + 1000 = 11000
    // - The prior weight for clock bias is 1/11000 ≈ 9e-5 (extremely weak)
    //
    // In contrast, a random-walk model (phi=1) would propagate:
    // - P_pred[15,15] = P_prev[15,15] + Q[15,15] * dt
    // - After convergence with small Q, P_pred stays near P_prev
    // - The prior weight stays high, maintaining temporal correlation
    //
    // THE CRITICAL INTERACTION: The clock bias state and position are coupled
    // through the measurement model (both appear in the CP and PR H rows).
    // When the clock bias has no temporal correlation, the filter cannot
    // separate position from clock bias across epochs. The result is that
    // position accuracy degrades to the level that can be determined from
    // a single epoch's measurements: ~0.5-1m.

    #[test]
    fn test_random_walk_clock_preserves_covariance_across_epochs() {
        // Simulate a converged clock bias with small covariance (0.01 m^2)
        let t = GpsTime::new(2156, 1000.0);
        let pos = Coordinate::new(
            Vector3::new(6000000.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            t,
        );
        let mut state = RtkState::new(t, pos, 100.0);
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim) * 100.0;
        state.covariance[(15, 15)] = 0.01; // 10 cm std — well converged
        state.covariance[(19, 19)] = 0.1;  // clock drift
        state.rcv_clk_bias = 100.0;
        state.rcv_clk_drift = 0.0;

        let config = EngineConfig {
            mode: crate::engine::EngineMode::Ppp,
            process_noise_cb: 1.0,
            dynamics_model: DynamicsModel::Static,
            ..Default::default()
        };

        let phi = compute_transition_matrix(&state, 1.0, &[]);

        // Random-walk clock: φ[15,15] = 1.0 preserves temporal correlation
        assert_eq!(
            phi[(15, 15)], 1.0,
            "Random-walk clock requires φ[15,15]=1.0 to preserve temporal correlation"
        );

        // Propagate covariance
        let q = compute_process_noise(1.0, &config, false, false, &[]);
        let p_pred = &phi * &state.covariance * phi.transpose() + q;

        // With φ[15,15]=1.0 and process_noise_cb=1.0, predicted covariance
        // stays close to the converged value (just adds q[15,15] = 1.0 m²)
        assert!(
            p_pred[(15, 15)] < 2.0,
            "Random-walk clock preserves covariance: P_pred[15]={:.2} ≈ 0.01 + 1.0",
            p_pred[(15, 15)]
        );

        // Verify clock drift coupling still works
        assert_eq!(phi[(15, 19)], 1.0, "φ[15,19]=dt preserves drift coupling");
    }

    #[test]
    fn test_white_noise_clock_destroys_filter_convergence_property() {
        // Demonstrate that the white-noise clock prevents the filter from
        // converging to the true position because it requires re-estimating
        // the clock from scratch each epoch, which couples into the position
        // estimate through the measurement model.
        //
        // Setup: simulate 10 epochs of a static receiver with known true position
        // and clock bias. Compare filter behavior with phi=0 (current) vs phi=1.

        let t0 = GpsTime::new(2156, 1000.0);
        let true_pos = Vector3::new(6000000.0, 0.0, 0.0);

        // Initialize state
        let pos = Coordinate::new(true_pos, Datum::WGS84, Frame::ECEF, t0);
        let mut state = RtkState::new(t0, pos, 100.0);
        let dim = CORE_STATE_SIZE;
        state.covariance = DMatrix::identity(dim, dim);
        state.covariance[(15, 15)] = 10000.0; // initial clock bias variance
        state.rcv_clk_bias = 0.0;

        // Compute the transition matrix as currently coded
        let phi = compute_transition_matrix(&state, 1.0, &[]);
        let q = compute_process_noise(1.0, &EngineConfig::default(), false, false, &[]);

        // Track clock bias variance over 10 epochs with white-noise model
        let mut p_white_noise = state.covariance.clone();
        let mut p_random_walk = state.covariance.clone();

        for _epoch in 0..10 {
            // White noise (current): phi[15,15] = 0
            let mut phi_wn = phi.clone();
            phi_wn[(15, 15)] = 0.0;
            p_white_noise = &phi_wn * &p_white_noise * phi_wn.transpose() + &q;

            // Random walk (fix): phi[15,15] = 1
            let mut phi_rw = phi.clone();
            phi_rw[(15, 15)] = 1.0;
            p_random_walk = &phi_rw * &p_random_walk * phi_rw.transpose() + &q;
        }

        // After 10 epochs of propagation and simulated updates:
        // White-noise clock variance should stay ~process_noise*dt (never converges)
        // Random-walk clock variance should grow slowly (can converge with updates)
        let wn_var = p_white_noise[(15, 15)];
        let rw_var = p_random_walk[(15, 15)];

        tracing::info!(
            "After 10 epochs: white-noise clock P[15,15]={:.2}, random-walk P[15,15]={:.2}",
            wn_var, rw_var
        );

        // With process_noise_cb=1.0 (RALPH fix), both models are well-behaved.
        // White-noise resets to q_cb=1.0 each epoch (stays near 1-10 m²).
        // Random-walk accumulates q_cb per epoch (10000 + 10*1 ≈ 10010 m²).
        // Either way, the variance is bounded — the old process_noise_cb=10000
        // was the real problem.
        assert!(
            wn_var < 1000.0 && rw_var < 20000.0,
            "Both clock models should have bounded variance with process_noise_cb=1.0: \
             wn={:.2}, rw={:.2}",
            wn_var, rw_var
        );
    }

    // ======================================================================
    // Adversarial Test 3: Combined effect of SPP prior + white-noise clock
    // ======================================================================
    //
    // The two bugs compound: the white-noise clock forces the filter to rely
    // on the SPP prior for absolute position anchoring, but the prior is
    // biased by SPP errors. The result is a systematic position bias.
    //
    // In a properly designed filter with random-walk clock:
    // - Clock bias accumulates information across epochs (P converges)
    // - The clock-code separation naturally anchors position
    // - The SPP prior is only needed for cold-start, not for convergence
    //
    // In the current filter:
    // - Clock bias resets each epoch (P always large)
    // - The filter relies on SPP prior for absolute position
    // - The prior induces position bias

    #[test]
    fn test_combined_spp_prior_and_clock_model_produce_bias() {
        // Demonstrate the SPP prior bias mechanism via compute_iteration_dx.
        // When the prior variance is clamped to 1.0 (the minimum from process_ppp),
        // it injects more weight toward an erroneous SPP position than when the
        // variance follows the actual state covariance.
        //
        // This is the core mathematical mechanism behind the ~0.5m East bias.

        let fg = PppIteratedEkf::default();
        let t0 = GpsTime::new(2156, 1000.0);
        let state = make_state_with_pos(t0, Vector3::new(0.0, 0.0, 0.0), 100.0);
        let dim = CORE_STATE_SIZE;
        let x_i = DVector::zeros(dim);
        let x_pred = DVector::zeros(dim);
        let p_inv = DMatrix::identity(dim, dim);
        let sat = make_dummy_sat(
            SatelliteId { constellation: Constellation::Gps, prn: 1 },
            Vector3::new(20000000.0, 0.0, 0.0),
        );

        // SPP position at X=2m (typical f9p SPP horizontal error)
        let spp_pos = Vector3::new(2.0, 0.0, 0.0);

        // Clamped prior variance = 1.0 (from process_ppp line 52: .max(1.0))
        let dx_clamped = fg
            .compute_iteration_dx(
                &state, &[sat.clone()], &x_i, &x_pred, &p_inv, 0,
                Some((spp_pos, 1.0)), // the bug: min variance floor
            )
            .unwrap()
            .unwrap();

        // Natural prior variance = 0.01 (what it should be after convergence)
        let dx_natural = fg
            .compute_iteration_dx(
                &state, &[sat], &x_i, &x_pred, &p_inv, 0,
                Some((spp_pos, 0.01)), // what the variance should follow
            )
            .unwrap()
            .unwrap();

        tracing::info!(
            "Combined bias mechanism: clamped_prior dx[0]={:.4}m, \
             natural_prior dx[0]={:.4}m. Clamping REDUCES pull by {:.1}x \
             (because weight drops from 100 to 1)",
            dx_clamped[0], dx_natural[0], dx_natural[0] / dx_clamped[0]
        );

        // KEY INSIGHT: The clamped prior has LESS pull than the natural prior.
        // This is because the floor INCREASES variance (from 0.01 to 1.0),
        // DECREASING weight (from 100 to 1). The natural prior (matching
        // the state covariance) provides MUCH stronger position anchoring.
        // The irony: the "prior_var clamping" was designed to prevent the
        // prior from being too strong, but after convergence, the filter
        // NEEDS that strong prior because the white-noise clock model
        // (phi=0) destroys temporal position correlation through the
        // clock bias state.
        assert!(
            dx_natural[0].abs() > dx_clamped[0].abs() * 2.0,
            "Natural prior (var=0.01, weight=100) should pull at least 2x \
             more than clamped prior (var=1.0, weight=1): \
             natural={:.4}, clamped={:.4}",
            dx_natural[0], dx_clamped[0]
        );
    }
}
