use crate::engine::ppp_common::{apply_state_vector, extract_state_vector, find_amb_idx};
use crate::engine::processed_sat::ProcessedSat;
use crate::filter::{RtkState, CORE_STATE_SIZE};
use nalgebra::{DMatrix, DVector, Vector3};

// Type aliases for AR candidate tuples used throughout this module.
// (sat_id, n1_idx, n2_idx, el_rad, lam1, lam2)
type ArCandidate = (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64);
// WL result: (state_vector, covariance, kept_indices)
type WlResult = (DVector<f64>, DMatrix<f64>, Vec<usize>);
// NL result: (state_vector, covariance)
type NlResult = (DVector<f64>, DMatrix<f64>);

#[cfg(test)]
#[derive(Clone)]
pub struct ArMock {
    pub wl_result: Option<Result<WlResult, &'static str>>,
    pub nl_result: Option<Result<NlResult, &'static str>>,
    pub nl_calls: usize,
}

#[cfg(test)]
pub static AR_MOCK: std::sync::Mutex<Option<ArMock>> = std::sync::Mutex::new(None);

use super::ppp_iekf::PppIteratedEkf;

impl PppIteratedEkf {
    /// Try to fix ambiguities for a single constellation group using WL+NL cascade.
    /// Returns (x_fixed, p_fixed, n_sats) on success, or None if this group cannot fix.
    pub(super) fn process_constellation_group(
        &self,
        state: &RtkState,
        p_current: &DMatrix<f64>,
        x_current: &DVector<f64>,
        group_cands: &[ArCandidate],
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
            .map(|c| (*c, *ref_cand))
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
    pub(super) fn try_inter_constellation_fallback(
        &self,
        state: &RtkState,
        p_current: &DMatrix<f64>,
        x_current: &DVector<f64>,
        cands: &[ArCandidate],
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
            Vec<ArCandidate>,
        > = std::collections::HashMap::new();
        for cand in &cands {
            const_groups
                .entry(cand.0.constellation)
                .or_default()
                .push(*cand);
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

    pub(super) fn find_ar_candidates(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
    ) -> Vec<ArCandidate> {
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

    pub(super) fn build_ar_subset(
        &self,
        cands: &[ArCandidate],
    ) -> Vec<(
        ArCandidate,
        ArCandidate,
    )> {
        // Inter-constellation: single highest-elevation GPS as universal reference.
        // Fall back to per-constellation if no GPS available.
        if let Some(ref_cand) = cands
            .iter()
            .filter(|c| c.0.constellation == gneiss_core::sat::Constellation::Gps)
            .max_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal))
        {
            let ref_cand = *ref_cand;
            return cands
                .iter()
                .filter(|c| c.0 != ref_cand.0)
                .map(|c| (*c, ref_cand))
                .collect();
        }
        // Fallback: per-constellation
        let mut subset = Vec::new();
        let mut const_cands = std::collections::HashMap::new();
        for cand in cands {
            const_cands
                .entry(cand.0.constellation)
                .or_insert_with(Vec::new)
                .push(*cand);
        }
        for (_, mut group) in const_cands {
            if group.len() < 2 {
                continue;
            }
            group.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap_or(std::cmp::Ordering::Equal));
            let ref_cand = group[0];
            for cand in group.iter().skip(1) {
                subset.push((*cand, ref_cand));
            }
        }
        subset
    }

    pub(super) fn resolve_widelane_ar(
        &self,
        state: &RtkState,
        p: &DMatrix<f64>,
        subset: &[(
            ArCandidate,
            ArCandidate,
        )],
        x: &DVector<f64>,
    ) -> Result<WlResult, &'static str> {
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
                let cov_ok = q_wl_full[(i, i)].sqrt() < 0.30; // cycles — safe threshold
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
        let all_mw_confident = keep_indices.iter().all(|&idx| {
            let (c, ref_sat) = &subset[idx];
            state.mw_sd_counts.get(&c.0).copied().unwrap_or(0) > 50
                && state.mw_sd_counts.get(&ref_sat.0).copied().unwrap_or(0) > 50
        });
        let q_wl = if all_mw_confident {
            // MW-based Q: tight (~0.18/N cycles²)
            // MW per-sample DD variance: each single-epoch MW measurement has
            // ~0.42 cycle std on GPS L1/L2, so 0.18 cycles² per sample.
            // Reference satellite noise is shared across all DD pairs.
            let mw_var_per_sample: f64 = 0.18;
            let mut q = DMatrix::zeros(n, n);
            for i in 0..n {
                let (c_i, ref_sat_i) = &subset[keep_indices[i]];
                let cnt_i = state.mw_sd_counts.get(&c_i.0).copied().unwrap_or(1);
                let cnt_ref = state.mw_sd_counts.get(&ref_sat_i.0).copied().unwrap_or(1);
                let var_i = mw_var_per_sample / cnt_i as f64;
                let var_ref = mw_var_per_sample / cnt_ref as f64;
                q[(i, i)] = (var_i + var_ref).max(0.0025);
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
        // Use the configurable AR ratio threshold (default 1.5, recommend 3.0
        // for stations with weak geometry).  Bootstrapping success_rate
        // provides secondary gating.  Hardcoded 1.1 was too lenient and
        // allowed wrong integer fixes on 3 of 7 IGS stations.
        let wl_ok = res_wl.ratio >= self.lambda_min_ratio && res_wl.success_rate >= 0.05;
        if !wl_ok {
            return Err("WL ratio test failed");
        }

        let s_inv = q_wl.try_inverse().ok_or("WL Cov Inversion failed")?;
        let k_wl = p * d_wl.transpose() * s_inv;
        let dx_wl = &k_wl * (res_wl.best_integers - a_wl);
        // Tiny regularization prevents p_wl from going singular (Joseph form
        // with zero measurement noise collapses rank). 1e-6 m² per pair keeps
        // the covariance full-rank for downstream NL LAMBDA and gain inversion.
        let r_wl = DMatrix::identity(keep_indices.len(), keep_indices.len()) * 0.01; // σ=10cm soft lock
        Ok((
            x + dx_wl,
            crate::math::covariance::apply_joseph_covariance_update(p, &k_wl, &d_wl, &r_wl),
            keep_indices,
        ))
    }

    pub(super) fn resolve_narrowlane_ar(
        &self,
        _state: &RtkState,
        subset: &[(
            ArCandidate,
            ArCandidate,
        )],
        keep_indices: &[usize],
        x_wl: &DVector<f64>,
        p_wl: &DMatrix<f64>,
    ) -> Result<NlResult, &'static str> {
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

        // --- Diagnostic: per-pair NL fix quality ---
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
            &(DMatrix::identity(keep_indices.len(), keep_indices.len()) * 0.01), // σ=10cm soft lock
        );

        let pos_corr_str2 = format!("{:.3}", pos_correction_norm);
        tracing::info!(
            "PPP Cascade AR Fixed! N_Sats: {} pos_corr={}m",
            keep_indices.len() + 1,
            pos_corr_str2
        );
        Ok((x_wl + dx_nl, p_fixed))
    }
}
