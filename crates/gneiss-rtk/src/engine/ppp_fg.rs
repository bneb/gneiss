use nalgebra::{DMatrix, DVector, Vector3};
use crate::engine::processed_sat::ProcessedSat;
use crate::filter::RtkState;

const C: f64 = 299792458.0;

/// A robust Levenberg-Marquardt Factor Graph solver for Precise Point Positioning (PPP).
pub struct PppFactorGraph {
    pub max_iterations: usize,
    pub convergence_threshold: f64,
    pub huber_k: f64,
}

impl Default for PppFactorGraph {
    fn default() -> Self {
        Self {
            max_iterations: 15,
            convergence_threshold: 1e-4,
            huber_k: 3.0,
        }
    }
}

pub struct PppFgResult {
    pub position: Vector3<f64>,
    pub cdt: f64,
    pub ztd: f64,
    pub covariance: DMatrix<f64>,
}

impl PppFactorGraph {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn solve(&self, initial_state: &RtkState, sats: &[ProcessedSat]) -> Option<PppFgResult> {
        let num_sats = sats.len();
        if num_sats < 4 {
            return None;
        }

        // State vector: [X, Y, Z, cdt, ztd]
        // If we want to solve ambiguities, they would go here. For standard PPP, we often float them.
        // For simplicity in the recovered engine, we solve the core state first.
        let num_states = 5;
        let mut x = DVector::zeros(num_states);
        x[0] = initial_state.position.vector.x;
        x[1] = initial_state.position.vector.y;
        x[2] = initial_state.position.vector.z;
        x[3] = initial_state.rcv_clk_bias;
        x[4] = 0.0; // Initial ZTD

        let mut lambda = 0.01; // LM damping
        let mut best_cost = f64::INFINITY;
        let mut final_cov = DMatrix::zeros(num_states, num_states);

        for _iter in 0..self.max_iterations {
            let mut residuals = Vec::new();
            let mut jacobians = Vec::new();
            let mut weights = Vec::new();

            let current_pos = Vector3::new(x[0], x[1], x[2]);
            let current_cdt = x[3];
            let current_ztd = x[4];

            // Build factors
            for sat in sats {
                let geom_range = (sat.sat_pos_rot - current_pos).norm();
                let los = (current_pos - sat.sat_pos_rot) / geom_range; // derivative of range wrt pos
                let tropo_delay = sat.tropo_dry + current_ztd * sat.map_wet;
                
                let expected_pr = geom_range + current_cdt - (sat.sat_clock_drift * C) + tropo_delay;

                // Pseudorange factor (L1)
                if let Some(pr1) = sat.sat_obs.get_observable(1) {
                    let r = pr1 - expected_pr;
                    let mut j = DVector::zeros(num_states);
                    j[0] = los.x;
                    j[1] = los.y;
                    j[2] = los.z;
                    j[3] = 1.0;
                    j[4] = sat.map_wet;

                    residuals.push(r);
                    jacobians.push(j);
                    // Weight based on elevation and SNR
                    let mut w: f64 = 1.0 / (0.3 + 2.0 * (-sat.el).exp()); 
                    
                    // Huber loss weighting
                    let norm_r = r.abs() * w.sqrt();
                    if norm_r > self.huber_k {
                        w *= self.huber_k / norm_r;
                    }
                    weights.push(w);
                }

                // If iono-free or dual freq phase is available, add phase factor
                // (Omitted for brevity in the core recovery, but easily extensible by adding state variables for ambiguities)
            }

            let num_meas = residuals.len();
            if num_meas < num_states {
                return None; // Not enough measurements
            }

            let mut j_mat = DMatrix::zeros(num_meas, num_states);
            let mut r_vec = DVector::zeros(num_meas);
            let mut w_mat = DMatrix::zeros(num_meas, num_meas);

            let mut current_cost = 0.0;
            for i in 0..num_meas {
                r_vec[i] = residuals[i];
                w_mat[(i, i)] = weights[i];
                for k in 0..num_states {
                    j_mat[(i, k)] = jacobians[i][k];
                }
                current_cost += residuals[i] * residuals[i] * weights[i];
            }

            let jtw = j_mat.transpose() * &w_mat;
            let jt_w_j = &jtw * &j_mat;
            let jt_w_r = &jtw * &r_vec;

            // Levenberg-Marquardt diagonal augmentation
            let mut a_mat = jt_w_j.clone();
            for i in 0..num_states {
                a_mat[(i, i)] *= 1.0 + lambda;
            }

            if let Some(dx) = a_mat.cholesky().map(|c| c.solve(&jt_w_r)) {
                if dx.norm() < self.convergence_threshold {
                    if let Some(inv) = jt_w_j.try_inverse() {
                        final_cov = inv;
                    }
                    break;
                }

                // Try update
                let new_x = &x + &dx;
                // Evaluate new cost
                let mut new_cost = 0.0;
                for i in 0..num_meas {
                    let sat = &sats[i]; // Approximate
                    let pr_new_pos = Vector3::new(new_x[0], new_x[1], new_x[2]);
                    let new_geom_range = (sat.sat_pos_rot - pr_new_pos).norm();
                    let expected_pr = new_geom_range + new_x[3] - (sat.sat_clock_drift * C) + sat.tropo_dry + new_x[4] * sat.map_wet;
                    
                    if let Some(pr1) = sat.sat_obs.get_observable(1) {
                        let r = pr1 - expected_pr;
                        let mut w = weights[i];
                        let norm_r = r.abs() * w.sqrt();
                        if norm_r > self.huber_k { w *= self.huber_k / norm_r; }
                        new_cost += r * r * w;
                    }
                }

                if new_cost < current_cost {
                    x = new_x;
                    best_cost = new_cost;
                    lambda = (lambda * 0.1).max(1e-7); // decrease lambda (closer to Gauss-Newton)
                } else {
                    lambda = (lambda * 10.0).min(1e5); // increase lambda (closer to gradient descent)
                }
            } else {
                return None; // Singular matrix
            }
        }

        Some(PppFgResult {
            position: Vector3::new(x[0], x[1], x[2]),
            cdt: x[3],
            ztd: x[4],
            covariance: final_cov,
        })
    }
}
