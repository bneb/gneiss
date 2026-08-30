pub mod normal_eq;
#[cfg(test)]
mod tests;

use normal_eq::{accumulate_factor_normal_equations, solve_linear_system};
use nalgebra::{DMatrix, DVector};

use crate::swfg::config::EngineConfig;
use crate::swfg::graph::EstimationGraph;
use crate::swfg::variables::{VariableId, VariableKind, VariableValues};

/// Error types for the sliding-window solver.
#[derive(Debug)]
pub enum SolveError {
    /// The normal equations are singular or near-singular.
    SingularSystem,
    /// Maximum iterations exceeded without convergence.
    DidNotConverge { final_error: f64, iterations: usize },
    /// A factor referenced a variable not in the graph.
    MissingVariable(VariableId),
    /// Schur complement marginalization failed (singular block).
    MarginalizationFailed(String),
    /// Unconnected variable detected in graph structure.
    OrphanVariable(String),
}

/// The sliding-window factor graph solver.
///
/// Maintains a window of the most recent N epochs.  Each call to
/// `process_epoch` adds GNSS + IMU factors for one new epoch,
/// runs Levenberg-Marquardt optimization, and marginalizes the
/// oldest epoch if the window is full.
pub struct SlidingWindowSolver {
    pub graph: EstimationGraph,
    /// Maximum number of epochs in the sliding window.
    window_size: usize,
    /// Current epoch counter (monotonically increasing).
    current_epoch: u32,
    /// Maximum LM iterations per solve.
    max_lm_iterations: usize,
    /// Convergence tolerance (norm of delta).
    convergence_tol: f64,
}

impl SlidingWindowSolver {
    /// Create a new solver from the engine configuration.
    pub fn new(config: &EngineConfig) -> Self {
        let window_size = match config {
            EngineConfig::Spp(_) => 1,
            EngineConfig::Ppp(ppp) => ppp.window_size.max(1),
            EngineConfig::Rtk(rtk) => rtk.window_size.max(1),
            EngineConfig::RtkIns(rtk_ins) => rtk_ins.rtk.window_size.max(1),
            EngineConfig::PppIns(ppp_ins) => ppp_ins.ppp.window_size.max(1),
        };
        Self {
            graph: EstimationGraph::new(),
            window_size,
            current_epoch: 0,
            max_lm_iterations: 10,
            convergence_tol: 1e-4,
        }
    }

    /// Create per-epoch variables for a new epoch.  Returns the IDs of
    /// all variables created, for later marginalization.
    ///
    /// When `is_dd_mode` is true (DD-RTK), clock biases, ZWD, and IFB are
    /// skipped because double-differencing eliminates receiver clocks and
    /// common-mode atmospheric errors.
    pub fn create_epoch_variables(
        &mut self,
        epoch: u32,
        _n_satellites: usize,
        _constellations: &[u8],
        has_imu: bool,
        is_dd_mode: bool,
    ) -> Vec<VariableId> {
        let mut vars = Vec::new();

        // Every epoch gets Pose
        let pose_id = self.graph.add_variable(VariableKind::Pose { epoch });
        vars.push(pose_id);

        // Velocity only if IMU is present
        if has_imu {
            let vel_id = self.graph.add_variable(VariableKind::Velocity { epoch });
            vars.push(vel_id);
        }

        if !is_dd_mode {
            // Clock biases per constellation (eliminated by DD)
            for &constellation_id in _constellations {
                let clk_id = self.graph.add_variable(VariableKind::ClockBias {
                    epoch,
                    constellation_id,
                });
                vars.push(clk_id);
            }

            // Troposphere ZWD (common-mode canceled by DD for short baselines)
            let zwd_id = self.graph.add_variable(VariableKind::TropoZwd { epoch });
            vars.push(zwd_id);
        }

        // IMU bias (one per session, created as soon as IMU is present)
        if has_imu && !self.graph.variables.values().any(|n| matches!(n.kind, VariableKind::ImuBias)) {
            let imu_id = self.graph.add_variable(VariableKind::ImuBias);
            vars.push(imu_id);
            let bias_prior = crate::swfg::pipeline::bias_factor::ImuBiasPriorFactor::new(
                imu_id,
                nalgebra::Vector3::zeros(),
                nalgebra::Vector3::zeros(),
                0.04,
                0.0001,
            );
            self.graph.add_factor(Box::new(bias_prior));
        }

        if !has_imu {
            // GNSS-only mode: the orientation part of Pose (indices 3, 4, 5) is completely unobserved.
            // Add a prior to constrain it, making the Hessian positive definite so that `try_inverse` succeeds.
            let mut info = nalgebra::DMatrix::zeros(6, 6);
            info[(3, 3)] = 1e8;
            info[(4, 4)] = 1e8;
            info[(5, 5)] = 1e8;
            let orient_prior = crate::swfg::factor::PriorFactor::new(
                pose_id,
                nalgebra::DVector::zeros(6),
                1.0, // Scale is overridden by the custom info matrix below
            );
            // Overwrite with our custom info matrix
            let mut orient_prior = orient_prior;
            orient_prior.information = info;
            self.graph.add_factor(Box::new(orient_prior));
        }

        vars
    }

    /// Lazily find-or-create the session's single GLONASS inter-frequency-
    /// bias variable. Deliberately NOT created ahead of time from raw
    /// satellite tracking (that was the bug: a GLONASS satellite can be
    /// tracked in the observation file with no matching ephemeris —
    /// missing/partial nav data — never survive to a processed
    /// observation, and never get a factor, leaving a pre-created
    /// variable permanently orphaned and every subsequent solve failing).
    /// Call this only where a factor referencing the result is about to
    /// be built, mirroring [`Self::ensure_ambiguity`]'s pattern: existence
    /// and factor coverage can then never disagree.
    pub fn ensure_ifb_glonass(&mut self) -> VariableId {
        for node in self.graph.variables.values() {
            if matches!(node.kind, VariableKind::IfbGlonass) {
                return node.id;
            }
        }
        self.graph.add_variable(VariableKind::IfbGlonass)
    }

    /// Create an ambiguity variable for a satellite-frequency pair.
    /// These persist across ALL epochs — they are NOT per-epoch.
    pub fn ensure_ambiguity(&mut self, satellite: u16, frequency: u8) -> VariableId {
        // Check if already exists (by scanning — could be optimized with a
        // separate lookup map).
        for node in self.graph.variables.values() {
            if let VariableKind::Ambiguity {
                satellite: s,
                frequency: f,
            } = node.kind
            {
                if s == satellite && f == frequency {
                    return node.id;
                }
            }
        }
        self.graph
            .add_variable(VariableKind::Ambiguity { satellite, frequency })
    }

    /// Create a double-differenced ambiguity variable for a (sat, ref_sat) pair.
    /// These persist across ALL epochs in RTK mode.
    pub fn ensure_dd_ambiguity(
        &mut self,
        constellation_id: u8,
        satellite: u16,
        ref_satellite: u16,
        frequency: u8,
        arc: u32,
        initial_value: f64,
    ) -> (VariableId, bool) {
        for node in self.graph.variables.values() {
            if let VariableKind::DdAmbiguity {
                constellation_id: c,
                satellite: s,
                ref_satellite: r,
                frequency: f,
                arc: a,
            } = node.kind
            {
                if c == constellation_id && s == satellite && r == ref_satellite && f == frequency && a == arc {
                    return (node.id, false);
                }
            }
        }
        let id = self.graph.add_variable(VariableKind::DdAmbiguity {
            constellation_id,
            satellite,
            ref_satellite,
            frequency,
            arc,
        });
        if let Some(node) = self.graph.variables.get_mut(&id) {
            node.set_value(&[initial_value]);
        }
        (id, true)
    }

    /// Build the normal equations (J^T W J, J^T W r) from all factors
    /// and the marginal prior.
    pub fn build_normal_equations(
        &self,
        values: &VariableValues,
    ) -> (DMatrix<f64>, DVector<f64>, f64) {
        let total_dim = values.total_dim();
        let mut jtj = DMatrix::zeros(total_dim, total_dim);
        let mut jtr = DVector::zeros(total_dim);
        let mut total_error = 0.0_f64;

        // Accumulate from factors
        for factor in &self.graph.factors {
            total_error += accumulate_factor_normal_equations(factor.as_ref(), values, &mut jtj, &mut jtr);
        }

        // Add marginal prior if present.  The prior's Hessian and gradient
        // are in the subspace of kept variables — we expand them to the
        // full state vector by mapping each variable to its current index.
        if let Some(ref prior) = self.graph.marginal_prior {
            let mut kept_indices: Vec<(usize, usize)> = Vec::new(); // (prior_idx, sys_idx)
            let mut current_x = prior.x0.clone();
            
            let mut prior_idx = 0;
            for &(var_id, dim) in &prior.variables {
                if let (Some((start, _)), Some(val)) = (values.index_of(var_id), values.get(var_id)) {
                    for k in 0..dim {
                        current_x[prior_idx + k] = val[k];
                        kept_indices.push((prior_idx + k, start + k));
                    }
                }
                prior_idx += dim;
            }
            
            let delta_x = &current_x - &prior.x0;
            let h_delta = &prior.hessian * &delta_x;
            
            // Prior cost = 1/2 * delta_x^T * H * delta_x + g^T * delta_x
            let prior_cost = 0.5 * (&delta_x.transpose() * &h_delta)[(0, 0)] + (&prior.gradient.transpose() * &delta_x)[(0, 0)];
            total_error += prior_cost;
            
            for &(p_i, sys_i) in &kept_indices {
                for &(p_j, sys_j) in &kept_indices {
                    jtj[(sys_i, sys_j)] += prior.hessian[(p_i, p_j)];
                }
                jtr[sys_i] += prior.gradient[p_i] + h_delta[p_i];
            }
        }

        (jtj, jtr, total_error)
    }

    /// Build normal equations ONLY for factors that are connected to the marginalized variables,
    /// PLUS the existing marginal prior (to re-center and merge it).
    /// This prevents double-counting active factors that remain in the graph!
    fn build_marginalization_equations(
        &self,
        values: &VariableValues,
        marg_ids: &[VariableId],
    ) -> (DMatrix<f64>, DVector<f64>) {
        let total_dim = values.total_dim();
        let mut jtj = DMatrix::zeros(total_dim, total_dim);
        let mut jtr = DVector::zeros(total_dim);

        // 1. Accumulate from factors connected to at least one marginalized variable
        for factor in &self.graph.factors {
            let connected = factor.variables().iter().any(|id| marg_ids.contains(id));
            if connected {
                accumulate_factor_normal_equations(factor.as_ref(), values, &mut jtj, &mut jtr);
            }
        }

        // 2. Add existing marginal prior (to re-center and merge it)
        if let Some(ref prior) = self.graph.marginal_prior {
            let mut kept_indices: Vec<(usize, usize)> = Vec::new();
            let mut current_x = prior.x0.clone();
            
            let mut prior_idx = 0;
            for &(var_id, dim) in &prior.variables {
                if let (Some((start, _)), Some(val)) = (values.index_of(var_id), values.get(var_id)) {
                    for k in 0..dim {
                        current_x[prior_idx + k] = val[k];
                        kept_indices.push((prior_idx + k, start + k));
                    }
                }
                prior_idx += dim;
            }
            
            let delta_x = &current_x - &prior.x0;
            let h_delta = &prior.hessian * &delta_x;
            
            for &(p_i, sys_i) in &kept_indices {
                for &(p_j, sys_j) in &kept_indices {
                    jtj[(sys_i, sys_j)] += prior.hessian[(p_i, p_j)];
                }
                jtr[sys_i] += prior.gradient[p_i] + h_delta[p_i];
            }
        }

        (jtj, jtr)
    }

    pub fn apply_delta(&mut self, delta: &DVector<f64>) {
        let mut offset = 0;
        for node in self.graph.variables.values_mut() {
            let dim = node.value.len();
            if let VariableKind::Pose { .. } = node.kind {
                // Position update is linear
                node.value[0] += delta[offset];
                node.value[1] += delta[offset + 1];
                node.value[2] += delta[offset + 2];
                // Attitude update on SO(3) manifold: q_new = q_old * exp(delta)
                let q_old = nalgebra::UnitQuaternion::from_scaled_axis(
                    nalgebra::Vector3::new(node.value[3], node.value[4], node.value[5])
                );
                let dq = nalgebra::UnitQuaternion::from_scaled_axis(
                    nalgebra::Vector3::new(delta[offset + 3], delta[offset + 4], delta[offset + 5])
                );
                let q_new = q_old * dq;
                let new_rot_vec = q_new.scaled_axis();
                node.value[3] = new_rot_vec.x;
                node.value[4] = new_rot_vec.y;
                node.value[5] = new_rot_vec.z;
            } else {
                node.value += delta.rows(offset, dim);
            }
            offset += dim;
        }
    }

    /// Extract the current state vector (packed variable values).
    pub fn extract_state_vector(&self) -> DVector<f64> {
        let total_dim = self.graph.total_dim();
        let mut state = DVector::zeros(total_dim);
        let mut offset = 0;
        for node in self.graph.variables.values() {
            let dim = node.value.len();
            state.rows_mut(offset, dim).copy_from(&node.value);
            offset += dim;
        }
        state
    }

        /// Run the full Levenberg-Marquardt solver on the current graph.
    ///
    /// Returns the optimized state vector, or an error if the solve fails.
    pub fn solve(&mut self) -> Result<DVector<f64>, SolveError> {
        self.graph.validate_graph_structure().map_err(SolveError::OrphanVariable)?;
        let mut lambda = 1e-6;
        let mut prev_error = f64::MAX;
        let mut consecutive_rejections = 0;

        let values = VariableValues::build(&self.graph.variables);
        let (mut jtj, mut jtr, mut current_error) = self.build_normal_equations(&values);

        for _iteration in 0..self.max_lm_iterations {
            let total_dim = jtj.nrows();
            
            // Damped normal equations: (J^T W J + λ diag(J^T W J) + 1e-4 I) Δx = -J^T W r
            // Marquardt scale-invariant damping with Tikhonov regularization guarantees positive definiteness.
            let mut damped = jtj.clone();
            for i in 0..total_dim {
                damped[(i, i)] += lambda * damped[(i, i)] + 1e-4;
            }

            let delta = match solve_linear_system(&damped, &(-&jtr)) {
                Some(d) => d,
                None => return Err(SolveError::SingularSystem),
            };

            if delta.norm() < self.convergence_tol {
                break; // converged
            }
            
            self.apply_delta(&delta);
            let new_values = VariableValues::build(&self.graph.variables);
            let new_error = self.evaluate_total_error(&new_values);

            if new_error < current_error {
                // Accept step: compute new Jacobians and normal equations
                consecutive_rejections = 0;
                lambda = (lambda * 0.5).max(1e-7);
                let (new_jtj, new_jtr, _) = self.build_normal_equations(&new_values);
                jtj = new_jtj;
                jtr = new_jtr;
                // Stop if relative error reduction is tiny
                if (prev_error - new_error).abs() / current_error.max(1.0) < 1e-4 {
                    break; // converged
                }
                prev_error = current_error;
                current_error = new_error;
            } else {
                // Reject step: revert state vector and increase damping without computing Jacobians
                self.apply_delta(&(-&delta)); // Revert
                lambda *= 10.0;
                consecutive_rejections += 1;
                if consecutive_rejections >= 3 {
                    break;
                }
            }
        }

        Ok(self.extract_state_vector())
    }

    /// Compute total robust error across all factors without calculating Jacobians.
    pub fn evaluate_total_error(&self, values: &VariableValues) -> f64 {
        let mut total = 0.0;
        for factor in &self.graph.factors {
            let r = factor.residual(values);
            let w = factor.information();
            let r_norm = r.norm();
            let raw_s = (&r.transpose() * &w * &r)[(0, 0)];
            let s = factor.robust_threshold().map_or(raw_s, |k| {
                if r_norm < 1e-12 { raw_s } else { raw_s * crate::swfg::factor::compute_robust_error(r_norm, k, factor.use_cauchy()) / (r_norm * r_norm) }
            });
            total += s;
        }
        total
    }

    /// Advance to the next epoch.
    pub fn advance_epoch(&mut self) {
        self.current_epoch += 1;
    }

    pub fn current_epoch(&self) -> u32 {
        self.current_epoch
    }

    pub fn window_size(&self) -> usize {
        self.window_size
    }

    /// Marginalize the oldest epoch when the window exceeds `window_size`.
    pub fn marginalize_oldest_epoch(&mut self, oldest_epoch: u32) -> Result<(), SolveError> {
        let (marg_ids, kept_ids) = crate::swfg::marginalization::partition_variables(&self.graph, oldest_epoch);
        if marg_ids.is_empty() {
            return Ok(());
        }
        let values = VariableValues::build(&self.graph.variables);
        let (jtj, jtr) = self.build_marginalization_equations(&values, &marg_ids);
        
        if let Some(prior) = crate::swfg::marginalization::marginalize(
            &self.graph, &jtj, &jtr, &marg_ids, &kept_ids
        ) {
            crate::swfg::marginalization::apply_marginalization(&mut self.graph, prior, &marg_ids);
        } else {
            // If marginalization fails (singular block), discard oldest epoch variables so they do not linger as orphans
            for &id in &marg_ids {
                self.graph.remove_variable(id);
            }
        }
        Ok(())
    }
}
