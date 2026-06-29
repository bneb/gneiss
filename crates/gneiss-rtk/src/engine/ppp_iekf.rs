use crate::engine::ppp_common::{
    apply_state_vector, assemble_matrices, build_weight_matrix, extract_state_vector,
    invert_matrix, FgMeasurement,
};

// These are used only by test modules within this file; #[cfg(test)] avoids
// unused-import warnings during library (non-test) compilation.
#[cfg(test)]
use crate::engine::ppp_common::{build_iono_constraint_row, find_ambiguity_index};
use crate::engine::processed_sat::ProcessedSat;
use crate::engine::types::IonosphereModel;
use crate::engine::EngineError;
use crate::filter::RtkState;

#[cfg(test)]
use crate::filter::CORE_STATE_SIZE;
use crate::math::inversion::solve_cholesky_svd;
use nalgebra::{DMatrix, DVector, Vector3};

// Re-export free functions and constants from the measurement module so that
// callers (including tests within this file) can reference them without
// qualifying the sibling module path.
pub(crate) use super::ppp_measurements::*;

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
    /// LAMBDA AR minimum ratio threshold (default 3.0 for safety).
    /// Hardcoded 1.1 allowed wrong fixes on 3 of 7 IGS stations.
    pub lambda_min_ratio: f64,
}

impl Default for PppIteratedEkf {
    fn default() -> Self {
        Self {
            max_iterations: 15,
            convergence_threshold: 1e-3,
            huber_k: 3.0,
            iono_model: IonosphereModel::Klobuchar,
            lambda_min_ratio: 3.0, // safe threshold — avoids wrong WL fixes
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

    pub fn with_lambda_min_ratio(mut self, ratio: f64) -> Self {
        self.lambda_min_ratio = ratio;
        self
    }

    pub fn solve(
        &self,
        state: &mut RtkState,
        sats: &[ProcessedSat],
        position_prior: Option<(Vector3<f64>, f64)>, // (ecef_m, variance_m2)
    ) -> Result<(), EngineError> {
        self.solve_with_fixed_amb(state, sats, position_prior, &[])
    }

    /// Solve with externally-fixed ambiguity constraints.
    /// `fixed_amb` is a list of `(ambiguity_index, value_meters)` pairs.
    /// Each is added as a pseudo-measurement with σ=1mm.
    pub fn solve_with_fixed_amb(
        &self,
        state: &mut RtkState,
        sats: &[ProcessedSat],
        position_prior: Option<(Vector3<f64>, f64)>,
        fixed_amb: &[(usize, f64)],
    ) -> Result<(), EngineError> {
        // Store target values for the pseudo-measurements.
        // These are (ambiguity_index, target_meters).
        let targets: Vec<(usize, f64)> = fixed_amb.to_vec();
        let has_constraints = !targets.is_empty();
        let saved_prior = position_prior; // Clone before inner solve consumes it
        for outer_iter in 0..4 {
            let done = if has_constraints {
                self.solve_inner_fixed(state, sats, position_prior, &targets)?
            } else {
                self.solve_inner(state, sats, position_prior)?
            };
            if done || outer_iter == 3 { break; }
        }
        let has_new_sats = sats.iter().any(|s| {
            let key = (s.sat_obs.sat, 0);
            state.last_observed.get(&key).copied().unwrap_or(0) == state.epoch_count as u32
        });
        if (!state.is_fixed || has_new_sats) && state.epoch_count > 10 {
            if let Err(e) = self.resolve_cascade_ar(state, sats) {
                tracing::info!("Cascade AR did not fix: {:?}", e);
            }
        }
        // Re-apply tight position prior after AR to prevent wrong fixes
        // from cascading into large position errors. With σ=1cm prior,
        // the IEKF measurement update will reconcile the fixed ambiguities
        // with the known position.
        if let Some((prior_pos, prior_var)) = saved_prior {
            if prior_var < 0.01 {
                // Tight prior: gently pull position toward known coordinates
                for i in 0..3 {
                    let innovation = prior_pos[i] - state.position.vector[i];
                    let weight = 0.1; // Soft correction: 10% toward prior per epoch
                    state.position.vector[i] += innovation * weight;
                }
            }
        }
        Ok(())
    }

    fn solve_inner_fixed(
        &self,
        state: &mut RtkState,
        sats: &[ProcessedSat],
        position_prior: Option<(Vector3<f64>, f64)>,
        targets: &[(usize, f64)],
    ) -> Result<bool, EngineError> {
        let x_pred = extract_state_vector(state);
        let p_pred = state.covariance.clone();
        state.full_x_predict = Some(x_pred.clone());
        state.full_p_predict = Some(p_pred.clone());
        let mut x_i = x_pred.clone();
        let p_inv = invert_matrix(&p_pred).ok_or(EngineError::StateDisappeared)?;
        for _iter in 0..self.max_iterations {
            let mut meas = self.build_measurements(state, sats, &x_i, _iter);
            if meas.is_empty() { return Err(EngineError::InsufficientSatellites); }
            // Add fixed ambiguity constraints: target = n_if, residual = n_if - x_i[idx]
            for &(amb_idx, target) in targets {
                let i = crate::filter::CORE_STATE_SIZE + amb_idx;
                if i < x_i.len() {
                    let mut h = DVector::zeros(x_i.len());
                    h[i] = 1.0;
                    meas.push(FgMeasurement { res: target - x_i[i], h_row: h, weight: 1e6, raw_var: 1e-6, is_phase: false, sat: None });
                }
            }
            let (h_mat, res_vec, r_mat) = assemble_matrices(&meas, x_i.len());
            let w_mat = build_weight_matrix(&meas, &r_mat);
            let mut htwh = h_mat.transpose() * &w_mat * &h_mat;
            let mut htwr = h_mat.transpose() * &w_mat * &res_vec;
            if let Some((spp_pos, var)) = position_prior {
                let w = 1.0 / var;
                for j in 0..3 { htwh[(j, j)] += w; htwr[j] += w * (spp_pos[j] - x_i[j]); }
            }
            let htwh_damped = &htwh + &p_inv;
            let innov = &htwr + &p_inv * (&x_pred - &x_i);
            match solve_cholesky_svd(&htwh_damped, &innov, 1e-6) {
                Ok(sol) => { if sol.norm() < self.convergence_threshold { break; } x_i += sol; }
                Err(_) => { tracing::warn!("Failed to solve normal equations in PPP FG!"); break; }
            }
        }
        let final_p = self.compute_final_covariance(state, sats, &x_i, &p_pred, &p_inv, position_prior);
        apply_state_vector(state, &x_i, final_p);
        Ok(true)
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

        let final_p = self.compute_final_covariance(state, sats, &x_i, &p_pred, &p_inv, position_prior);
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

    #[allow(clippy::too_many_arguments)]
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

        match solve_cholesky_svd(&htwh_damped, &innov, 1e-6) {
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
        position_prior: Option<(Vector3<f64>, f64)>,
    ) -> DMatrix<f64> {
        let last_meas = self.build_measurements(state, sats, x_i, self.max_iterations);
        if last_meas.is_empty() {
            return p_pred.clone();
        }

        let (h_mat, _, r_mat) = assemble_matrices(&last_meas, x_i.len());
        let w_mat = build_weight_matrix(&last_meas, &r_mat);

        let mut htwh = h_mat.transpose() * &w_mat * h_mat;

        // Include the position prior in the posterior covariance.
        // Without this, the prior decays via process noise each epoch:
        // after 100 epochs at 3e-5/epoch, σ_pos grows from 1cm to 5.6cm.
        // The prior is re-applied here so the posterior reflects ongoing
        // constraint from the known position.
        if let Some((_spp_pos, var)) = position_prior {
            let w = 1.0 / var;
            for i in 0..3 {
                htwh[(i, i)] += w;
            }
        }

        let htwh_damped = htwh + p_inv;

        let mut final_p = invert_matrix(&htwh_damped).unwrap_or_else(|| {
            tracing::warn!("invert_matrix(&htwh_damped) FAILED! Falling back to p_pred. htwh_damped has NaNs: {}, Infs: {}", htwh_damped.iter().any(|x| x.is_nan()), htwh_damped.iter().any(|x| x.is_infinite()));
            p_pred.clone()
        });

        // Clamp diagonal elements to safe range. Tight position priors
        // (σ=1cm) can make the normal equations nearly singular, producing
        // negative or exploding variances after inversion.
        for i in 0..final_p.nrows() {
            final_p[(i, i)] = final_p[(i, i)].clamp(1e-12, 1e8);
        }

        final_p
    }
}

#[cfg(test)]
#[path = "ppp_iekf_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "ppp_ar_tests.rs"]
mod adversarial_tests;
