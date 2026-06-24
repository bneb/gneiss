//! Multi-epoch sliding-window factor graph PPP optimizer.
//!
//! Replaces the single-epoch IEKF with a joint optimization across N epochs
//! to break the ~5m architectural accuracy floor. Between-epoch dynamics
//! constraints allow carrier phase to contribute relative position information
//! at mm precision without integer AR.
//!
//! Architecture (from SPRINT_PLAN.md):
//!   Phase A1: 2-epoch joint optimization — validate concept
//!   Phase A2: N-epoch sliding window (5-10 epochs)
//!   Phase A3: Shared ambiguities across window

use nalgebra::{DMatrix, DVector, Vector3};

use crate::engine::ppp_common::{apply_state_vector, extract_state_vector};
use crate::engine::processed_sat::ProcessedSat;
use crate::engine::EngineError;
use crate::filter::{RtkState, CORE_STATE_SIZE};

/// Multi-epoch sliding-window factor graph optimizer for PPP.
///
/// Jointly optimizes position, clock, tropo, ISB, and ambiguities across a
/// window of N epochs. Between-epoch dynamics constraints (position, clock,
/// tropo) allow carrier phase to contribute relative information.
pub struct MultiEpochOptimizer {
    /// Window size in epochs (2 = Phase A1, 5-10 = Phase A2)
    pub window_size: usize,
    /// Maximum LM iterations per window
    pub max_iterations: usize,
    /// Convergence threshold (norm of dx)
    pub convergence_threshold: f64,
    /// Huber loss parameter for robust estimation
    pub huber_k: f64,
}

impl Default for MultiEpochOptimizer {
    fn default() -> Self {
        Self {
            window_size: 2,
            max_iterations: 15,
            convergence_threshold: 1e-3,
            huber_k: 3.0,
        }
    }
}

impl MultiEpochOptimizer {
    pub fn new(window_size: usize) -> Self {
        Self {
            window_size,
            ..Default::default()
        }
    }

    /// Solve the multi-epoch PPP optimization.
    ///
    /// Phase A1 (current): 2-epoch joint optimization.
    /// State vector: [x_{k-1} (CORE_STATE+M), x_k (CORE_STATE+M)]
    ///
    /// Factors:
    ///   1. PR/CP/Doppler factors for epoch k (from `sats`)
    ///   2. Dynamics factor: position constrained between epochs
    ///   3. Clock constraint: random-walk between epochs
    ///   4. SPP position prior on current epoch
    ///   5. Weak prior on x_{k-1} from predicted covariance
    pub fn solve(
        &self,
        state: &mut RtkState,
        sats: &[ProcessedSat],
        position_prior: Option<(Vector3<f64>, f64)>,
    ) -> Result<(), EngineError> {
        if self.window_size < 2 {
            return Err(EngineError::StateDisappeared);
        }

        // ── Phase A1: 2-epoch joint optimization ──────────────────────
        let x_original = extract_state_vector(state);
        let n_state = x_original.len();

        // If no history available, fall back to IEKF
        let _x_prev = match &state.full_x_predict {
            Some(xp) => xp.clone(),
            None => return Err(EngineError::StateDisappeared),
        };

        let p_prev = match &state.full_p_predict {
            Some(pp) => pp.clone(),
            None => return Err(EngineError::StateDisappeared),
        };

        let dim = 2 * n_state;

        // --- Pre-compute linear factors (these don't change with iteration) ---
        let mut linear_rows: Vec<(DVector<f64>, f64, f64)> = Vec::new();

        // Factor 2: Dynamics prior — position stays constant (static)
        for i in 0..3 {
            let var = 3.0; // 0.1 m²/s * 30s
            let mut h_row = DVector::zeros(dim);
            h_row[i] = -1.0;
            h_row[n_state + i] = 1.0;
            linear_rows.push((h_row, 0.0, 1.0 / var));
        }

        // Factor 3: Clock constraint (random walk)
        if n_state > 15 {
            let var_clk = 30.0; // 1.0 m²/s * 30s
            let mut h_row = DVector::zeros(dim);
            h_row[15] = -1.0;
            h_row[n_state + 15] = 1.0;
            linear_rows.push((h_row, 0.0, 1.0 / var_clk));
        }

        // Factor 5: Prior on x_prev from predicted covariance P^{-1}.
        // Every element of x_prev must be constrained to avoid rank deficiency.
        // p_prev covers CORE_STATE_SIZE elements; ambiguities use state.covariance.
        let core_dim = CORE_STATE_SIZE.min(n_state);
        for i in 0..core_dim {
            let raw_var = p_prev[(i, i)];
            if raw_var <= 1e-12 || !raw_var.is_finite() {
                continue;
            }
            let prior_var = raw_var.clamp(0.01, 1e4);
            let mut h_row = DVector::zeros(dim);
            h_row[i] = 1.0;
            linear_rows.push((h_row, 0.0, 1.0 / prior_var));
        }
        // Ambiguity states: use current posterior covariance from state
        for i in core_dim..n_state {
            let raw_var = state.covariance[(i, i)];
            if raw_var <= 1e-12 || !raw_var.is_finite() {
                continue;
            }
            let prior_var = raw_var.clamp(0.01, 1e4);
            let mut h_row = DVector::zeros(dim);
            h_row[i] = 1.0;
            linear_rows.push((h_row, 0.0, 1.0 / prior_var));
        }

        // --- LM iteration ---
        let iektf = crate::engine::ppp_iekf::PppIteratedEkf::new();
        let max_iter = 3;
        let mut x_curr = x_original.clone();
        let mut lambda = 1.0; // LM damping

        for _iter in 0..max_iter {
            // Factor 1: Measurements at epoch k (re-linearized each iteration)
            let meas = iektf.build_measurements(state, sats, &x_curr, 0);
            if meas.is_empty() {
                return Err(EngineError::InsufficientSatellites);
            }

            // Factor 4: SPP position prior on current epoch
            let mut spp_rows: Vec<(DVector<f64>, f64, f64)> = Vec::new();
            if let Some((prior_pos, prior_var)) = position_prior {
                for i in 0..3 {
                    let mut h_row = DVector::zeros(dim);
                    h_row[n_state + i] = 1.0;
                    let res = x_curr[i] - prior_pos[i];
                    spp_rows.push((h_row, res, 1.0 / prior_var.max(1.0)));
                }
            }

            // Build normal equations
            let mut lhs = DMatrix::zeros(dim, dim);
            let mut rhs = DVector::zeros(dim);

            // Measurement rows
            for m in &meas {
                let mut h_row = DVector::zeros(dim);
                h_row
                    .rows_mut(n_state, m.h_row.len())
                    .copy_from(&m.h_row);
                let weight = 1.0 / m.raw_var.max(1e-12);
                for i in 0..dim {
                    let h_i_w = h_row[i] * weight;
                    for j in 0..dim {
                        lhs[(i, j)] += h_i_w * h_row[j];
                    }
                    rhs[i] += h_i_w * m.res;
                }
            }

            // SPP prior rows (recomputed each iter since residual depends on x_curr)
            for (h, z, w) in &spp_rows {
                for i in 0..dim {
                    let h_i_w = h[i] * w;
                    for j in 0..dim {
                        lhs[(i, j)] += h_i_w * h[j];
                    }
                    rhs[i] += h_i_w * z;
                }
            }

            // Linear rows (constant across iterations)
            for (h, z, w) in &linear_rows {
                for i in 0..dim {
                    let h_i_w = h[i] * w;
                    for j in 0..dim {
                        lhs[(i, j)] += h_i_w * h[j];
                    }
                    rhs[i] += h_i_w * z;
                }
            }

            // LM damping
            for i in 0..dim {
                lhs[(i, i)] += lambda;
            }

            // Solve
            let dx = match crate::math::inversion::solve_cholesky_svd(&lhs, &rhs, 1e-8) {
                Ok(dx) => dx,
                Err(e) => {
                    tracing::warn!("MultiEpoch solve failed: {:?}", e);
                    return Err(EngineError::StateDisappeared);
                }
            };

            // Update x_curr for current epoch block only
            let dx_curr = dx.rows(n_state, n_state).clone_owned();
            x_curr += dx_curr;

            // Reduce LM damping on convergence
            let pos_corr = dx.rows(n_state, 3).norm();
            if pos_corr < self.convergence_threshold {
                lambda *= 0.1;
                if pos_corr < 1e-4 {
                    break;
                }
            } else if pos_corr > 100.0 {
                lambda *= 10.0; // Diverge guard: increase damping
                // If diverging badly, fall back to IEKF and bail
                if pos_corr > 1000.0 {
                    tracing::warn!(
                        "MultiEpoch diverging (pos_corr={:.1}m), falling back to IEKF",
                        pos_corr
                    );
                    return Err(EngineError::StateDisappeared);
                }
            } else {
                lambda *= 0.5;
            }
        }

        // Extract posterior covariance for current block.
        // Use the (2,2) block of lhs^{-1} from the last iteration (without LM damping).
        // Rebuild lhs_final without regularization to get the information matrix,
        // then invert the (2,2) block.
        let meas_final = iektf.build_measurements(state, sats, &x_curr, 0);
        let mut lhs_final = DMatrix::zeros(dim, dim);
        for m in &meas_final {
            let mut h_row = DVector::zeros(dim);
            h_row
                .rows_mut(n_state, m.h_row.len())
                .copy_from(&m.h_row);
            let weight = 1.0 / m.raw_var.max(1e-12);
            for i in 0..dim {
                let h_i_w = h_row[i] * weight;
                for j in 0..dim {
                    lhs_final[(i, j)] += h_i_w * h_row[j];
                }
            }
        }
        // Add SPP prior to covariance info matrix
        if position_prior.is_some() {
            for i in 0..3 {
                let mut h_row = DVector::zeros(dim);
                h_row[n_state + i] = 1.0;
                let w = 1.0; // SPP prior variance clamped to ~1.0
                for k in 0..dim {
                    let h_k_w = h_row[k] * w;
                    for j in 0..dim {
                        lhs_final[(k, j)] += h_k_w * h_row[j];
                    }
                }
            }
        }
        for (h, _z, w) in &linear_rows {
            for i in 0..dim {
                let h_i_w = h[i] * w;
                for j in 0..dim {
                    lhs_final[(i, j)] += h_i_w * h[j];
                }
            }
        }
        for i in 0..dim {
            lhs_final[(i, i)] += 1e-8;
        }

        // Compute marginal covariance of x_k via Schur complement:
        //   P_k = (I_22 - I_21 * I_11^{-1} * I_12)^{-1}
        // Using I_22^{-1} alone is anti-conservative (underestimates variance),
        // which causes filter overconfidence and eventual divergence.
        let i11 = lhs_final.view((0, 0), (n_state, n_state)).into_owned();
        let i12 = lhs_final.view((0, n_state), (n_state, n_state)).into_owned();
        let i21 = lhs_final.view((n_state, 0), (n_state, n_state)).into_owned();
        let i22 = lhs_final.view((n_state, n_state), (n_state, n_state)).into_owned();

        let p_opt = if self.window_size == 1 {
            // Single epoch: no Schur complement needed
            crate::math::inversion::invert_matrix_robust(&i22)
        } else {
            // Schur complement with regularization
            let inv_i11 = crate::math::inversion::invert_matrix_robust(&i11);
            let s = &i21 * inv_i11 * &i12;
            let eff_info = i22 - s;
            // Add regularization
            let mut eff_reg = eff_info;
            for k in 0..eff_reg.nrows() {
                eff_reg[(k, k)] += 1e-8;
            }
            crate::math::inversion::invert_matrix_robust(&eff_reg)
        };

        apply_state_vector(state, &x_curr, p_opt);

        let pos_corr = ((x_curr[0] - x_original[0]).powi(2)
            + (x_curr[1] - x_original[1]).powi(2)
            + (x_curr[2] - x_original[2]).powi(2))
        .sqrt();
        tracing::info!(
            "MultiEpoch PPP: win={} pos_corr={:.3}m",
            self.window_size,
            pos_corr
        );
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multi_epoch_constructor() {
        let opt = MultiEpochOptimizer::new(2);
        assert_eq!(opt.window_size, 2);
        assert_eq!(opt.max_iterations, 15);
    }

    #[test]
    fn test_multi_epoch_rejects_single_epoch_window() {
        let opt = MultiEpochOptimizer::new(1);
        let time = gneiss_core::time::GpsTime::new(2082, 0.0);
        let coord = gneiss_core::coords::Coordinate::new(
            nalgebra::Vector3::zeros(),
            gneiss_core::coords::Datum::WGS84,
            gneiss_core::coords::Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, coord, 10.0);
        let sats: Vec<ProcessedSat> = vec![];
        let result = opt.solve(&mut state, &sats, None);
        assert!(result.is_err());
    }
}
