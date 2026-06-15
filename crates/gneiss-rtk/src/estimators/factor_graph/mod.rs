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
    
    pub fn optimize(&self, initial_state: &DVector<f64>, max_iters: usize, tol: f64) -> (DVector<f64>, DMatrix<f64>) {
        let mut state = initial_state.clone();
        let mut lambda = 1e-3;
        
        for iter in 0..max_iters {
            let (current_error, h, b) = self.build_normal_equations(&state, iter);
            if iter == 0 { tracing::debug!("Iter 0: b_norm={}, H_trace={}", b.norm(), h.trace()); }
            
            let delta = match Self::solve_normal_equations(h, &b, lambda) {
                Some(d) => d,
                None => { tracing::warn!("H matrix solve failed, breaking."); break; }
            };
            
            if delta.norm() < tol { break; }
            
            let new_state = &state - &delta;
            let new_error = self.compute_error(&new_state);
            
            if new_error < current_error * 1.5 {
                state = new_state;
                if new_error < current_error { lambda /= 10.0; }
            } else {
                lambda *= 10.0;
            }
        }
        (state.clone(), self.compute_final_covariance(&state))
    }

    fn build_normal_equations(&self, state: &DVector<f64>, iter: usize) -> (f64, DMatrix<f64>, DVector<f64>) {
        let mut h = DMatrix::zeros(state.len(), state.len());
        let mut b = DVector::zeros(state.len());
        let mut current_error = 0.0;
        
        for factor in &self.factors {
            let (mut info, res, jac) = (factor.information(), factor.residual(state), factor.jacobian(state));
            let maha_sq = (res.transpose() * &info * &res)[0];
            let mut cost = 0.5 * maha_sq;
            
            if let Some(k) = factor.robust_threshold() {
                let e = maha_sq.sqrt();
                if iter == 0 { tracing::debug!("Iter 0 factor: res={:.2}, e={:.2}, k={:.2}", res[0], e, k); }
                if e > k * 3.0 && factor.is_cauchy_rejectable() {
                    cost = k * (e - 0.5 * k);
                    info *= (k / e).powi(3);
                } else if e > k {
                    cost = k * (e - 0.5 * k);
                    info *= k / e;
                }
            }
            
            let j_t_info = jac.transpose() * &info;
            h += &j_t_info * &jac;
            b += &j_t_info * &res;
            current_error += cost;
        }
        (current_error, h, b)
    }

    fn solve_normal_equations(mut h: DMatrix<f64>, b: &DVector<f64>, lambda: f64) -> Option<DVector<f64>> {
        for i in 0..h.nrows() { h[(i, i)] += lambda * h[(i, i)].max(1e-9) + 1e-6; }
        h.clone().cholesky().map(|d| d.solve(b)).or_else(|| h.svd(true, true).solve(b, 1e-14).ok())
    }

    fn compute_error(&self, state: &DVector<f64>) -> f64 {
        self.factors.iter().map(|f| {
            let res = f.residual(state);
            let maha_sq = (res.transpose() * f.information() * &res)[0];
            match f.robust_threshold() {
                Some(k) if maha_sq.sqrt() > k => k * (maha_sq.sqrt() - 0.5 * k),
                _ => 0.5 * maha_sq,
            }
        }).sum()
    }

    fn compute_final_covariance(&self, state: &DVector<f64>) -> DMatrix<f64> {
        let mut h = DMatrix::zeros(state.len(), state.len());
        for factor in &self.factors {
            let (info, res, jac) = (factor.information(), factor.residual(state), factor.jacobian(state));
            let maha_sq = (res.transpose() * &info * &res)[0];
            let weight = factor.robust_threshold().map_or(1.0, |k| 1.0 / (1.0 + maha_sq / (k * k)));
            h += jac.transpose() * &(info * weight) * jac;
        }
        h.clone().cholesky().map(|c| c.inverse()).unwrap_or_else(|| h.pseudo_inverse(1e-9).unwrap_or_else(|_| DMatrix::identity(state.len(), state.len()) * 1e-6))
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
