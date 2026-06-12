pub mod gnss_factors;
pub mod imu_factors;

use nalgebra::{DMatrix, DVector};

/// A generic interface for a Factor in the graph.
pub trait Factor {
    /// Compute the residual error vector.
    fn residual(&self, state: &DVector<f64>) -> DVector<f64>;
    
    /// Compute the Jacobian of the residual with respect to the state.
    fn jacobian(&self, state: &DVector<f64>) -> DMatrix<f64>;
    
    /// Information matrix (inverse covariance) of the measurement.
    fn information(&self) -> DMatrix<f64>;
    
    /// Optional robust Huber loss threshold (k). If None, uses pure L2 loss.
    fn robust_threshold(&self) -> Option<f64> {
        None
    }
    
    fn is_cauchy_rejectable(&self) -> bool {
        false
    }
}

/// A generic Prior Factor on the entire state vector.
pub struct PriorFactor {
    pub information: DMatrix<f64>,
}

impl Factor for PriorFactor {
    fn residual(&self, state: &DVector<f64>) -> DVector<f64> {
        // Since state is the error state delta_x, and the prior mean is the nominal state (delta_x = 0),
        // the prior error is observed (0) - predicted (state) = -state.
        -state.clone()
    }
    
    fn jacobian(&self, state: &DVector<f64>) -> DMatrix<f64> {
        -DMatrix::identity(state.len(), state.len())
    }
    
    fn information(&self) -> DMatrix<f64> {
        self.information.clone()
    }
}

/// A simple Levenberg-Marquardt Optimizer for Factor Graphs.
pub struct FactorGraphOptimizer {
    pub factors: Vec<Box<dyn Factor>>,
}

impl Default for FactorGraphOptimizer {
    fn default() -> Self {
        Self::new()
    }
}

impl FactorGraphOptimizer {
    pub fn new() -> Self {
        Self { factors: Vec::new() }
    }
    
    pub fn add_factor(&mut self, factor: Box<dyn Factor>) {
        self.factors.push(factor);
    }
    
    /// Optimize the state vector using Levenberg-Marquardt.
    pub fn optimize(&self, initial_state: &DVector<f64>, max_iters: usize, tol: f64) -> (DVector<f64>, DMatrix<f64>) {
        let mut state = initial_state.clone();
        let mut lambda = 1e-3;
        
        for factor in &self.factors {
            let info = factor.information();
            let res = factor.residual(&state);
            let jac = factor.jacobian(&state);
            if info.nrows() != res.nrows() || jac.nrows() != res.nrows() || jac.ncols() != state.nrows() {
                println!("DIMENSION MISMATCH! info {}x{}, res {}x1, jac {}x{}", info.nrows(), info.ncols(), res.nrows(), jac.nrows(), jac.ncols());
            }
        }
        
        for _iter in 0..max_iters {
            let mut h = DMatrix::zeros(state.len(), state.len());
            let mut b = DVector::zeros(state.len());
            let mut current_error = 0.0;
            
            // Build the normal equations
            for factor in &self.factors {
                let res = factor.residual(&state);
                let jac = factor.jacobian(&state);
                let mut info = factor.information();
                
                let maha_sq = (res.transpose() * &info * &res)[0];
                let mut cost = 0.5 * maha_sq;
                
                if let Some(k) = factor.robust_threshold() {
                    let e = maha_sq.sqrt();
                    if _iter == 0 {
                        tracing::debug!("Iter 0 factor: res={:.2}, e={:.2}, k={:.2}, info={:.4}", res[0], e, k, info[(0,0)]);
                    }
                    if e > k * 3.0 && factor.is_cauchy_rejectable() {
                        // Extreme outlier -> Cauchy-like aggressive downweighting or hard reject
                        // We use a smooth aggressive downweighting: weight = (k/e)^3
                        cost = k * (e - 0.5 * k); // keep cost monotonic for line search if we had one
                        let weight = (k / e) * (k / e) * (k / e);
                        info *= weight;
                    } else if e > k {
                        cost = k * (e - 0.5 * k);
                        let weight = k / e;
                        info *= weight;
                    }
                }
                
                let j_t_info = jac.transpose() * &info;
                h += &j_t_info * &jac;
                b += &j_t_info * &res;
                
                current_error += cost;
            }
            
            // Add Levenberg-Marquardt damping
            for i in 0..state.len() {
                h[(i, i)] += lambda * h[(i, i)].max(1e-9) + 1e-6; // Add small constant to guarantee positive definiteness
            }
            
            if _iter == 0 {
                let mut max_h_diag = 0.0;
                let mut max_h_idx = 0;
                for i in 0..state.len() {
                    if h[(i, i)] > max_h_diag {
                        max_h_diag = h[(i, i)];
                        max_h_idx = i;
                    }
                }
                tracing::debug!("Iter 0: b_norm={}, H_trace={}", b.norm(), h.trace());
                tracing::debug!("Iter 0: state norm = {}", state.norm());
                let (max_i, _) = b.argmax();
                tracing::debug!("Iter 0: max b is at index {}, value {}", max_i, b[max_i]);
                tracing::debug!("Iter 0: max H diag is at index {}, value {}", max_h_idx, max_h_diag);
            }
            
            let mut delta = match h.clone().cholesky() {
                Some(d) => d.solve(&b),
                None => {
                    // Fallback to SVD if Cholesky fails
                    match h.clone().svd(true, true).solve(&b, 1e-14) {
                        Ok(sol) => sol,
                        Err(_) => {
                            tracing::warn!("H matrix SVD failed, breaking optimization early.");
                            break;
                        }
                    }
                }
            };
            
            // Trust region removed
            
            if delta.norm() < tol {
                break;
            }
            
            let mut new_state = state.clone();
            new_state -= &delta;
            
            // Compute new error to accept/reject step
            let mut new_error = 0.0;
            for factor in &self.factors {
                let res = factor.residual(&new_state);
                let info = factor.information();
                let maha_sq = (res.transpose() * &info * &res)[0];
                let mut cost = 0.5 * maha_sq;
                
                if let Some(k) = factor.robust_threshold() {
                    let e = maha_sq.sqrt();
                    if e > k {
                        cost = k * (e - 0.5 * k);
                    }
                }
                new_error += cost;
            }
            
            // Allow accepting slightly worse steps due to non-linearity, but not catastrophic
            if new_error < current_error * 1.5 {
                tracing::debug!("LM Iter {}: ACCEPTED step norm {:.2}, error {:.2} -> {:.2}", _iter, delta.norm(), current_error, new_error);
                state = new_state;
                if new_error < current_error {
                    lambda /= 10.0;
                }
            } else {
                tracing::debug!("LM Iter {}: REJECTED step norm {:.2}, error {:.2} -> {:.2}", _iter, delta.norm(), current_error, new_error);
                lambda *= 10.0;
            }
        }
        
        // Compute final covariance (inverse Hessian)
        let mut h = DMatrix::zeros(state.len(), state.len());
        for factor in &self.factors {
            let jac = factor.jacobian(&state);
            let info = factor.information();
            let res = factor.residual(&state);
            
            let maha_sq = (res.transpose() * &info * &res)[0];
            let mut weight = 1.0;
            
            // Apply robust loss weighting (IRLS)
            if let Some(k) = factor.robust_threshold() {
                let e_over_k_sq = maha_sq / (k * k);
                weight = 1.0 / (1.0 + e_over_k_sq);
            }
            
            let w_info = info * weight;
            h += jac.transpose() * &w_info * jac;
        }
        let cov = match h.clone().cholesky() {
            Some(chol) => {
                let inv = chol.inverse();
                tracing::debug!("Final H diag[12]: {}, cov diag[12]: {}", h[(12, 12)], inv[(12, 12)]);
                inv
            },
            None => {
                tracing::debug!("Cholesky decomposition failed for final covariance.");
                let pinv = h.clone().pseudo_inverse(1e-9).unwrap_or_else(|_| DMatrix::identity(state.len(), state.len()) * 1e-6);
                tracing::debug!("Final H diag[12]: {}, pinv cov diag[12]: {}", h[(12, 12)], pinv[(12, 12)]);
                pinv
            }
        };
        
        (state, cov)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{DMatrix, DVector};

    struct MockFactor {
        target: DVector<f64>,
    }

    impl Factor for MockFactor {
        fn residual(&self, state: &DVector<f64>) -> DVector<f64> {
            state - &self.target
        }
        
        fn jacobian(&self, _state: &DVector<f64>) -> DMatrix<f64> {
            DMatrix::identity(self.target.len(), self.target.len())
        }
        
        fn information(&self) -> DMatrix<f64> {
            DMatrix::identity(self.target.len(), self.target.len())
        }
    }

    #[test]
    fn test_factor_graph_optimizer_convergence() {
        let mut optimizer = FactorGraphOptimizer::new();
        let target = DVector::from_vec(vec![5.0, -3.0, 42.0]);
        optimizer.add_factor(Box::new(MockFactor { target: target.clone() }));
        
        let initial_state = DVector::from_vec(vec![0.0, 0.0, 0.0]);
        let (optimized, cov) = optimizer.optimize(&initial_state, 10, 1e-4);
        
        assert!((optimized - target).norm() < 1e-3);
        assert!((cov - DMatrix::identity(3, 3)).norm() < 1e-6);
    }

    #[test]
    fn test_prior_factor() {
        let info = DMatrix::from_diagonal(&DVector::from_vec(vec![2.0, 3.0]));
        let factor = PriorFactor { information: info.clone() };
        
        let state = DVector::from_vec(vec![1.5, -2.5]);
        let res = factor.residual(&state);
        assert_eq!(res, -state.clone());
        
        let jac = factor.jacobian(&state);
        assert_eq!(jac, -DMatrix::identity(2, 2));
        
        let info_out = factor.information();
        assert_eq!(info_out, info);
    }
}
pub fn foo(a: i32, b: i32) -> i32 { a + b }
