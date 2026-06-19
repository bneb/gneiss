use crate::engine::fgo::graph::FactorGraph;
use nalgebra::DMatrix;

/// Configuration for the Levenberg-Marquardt solver.
#[derive(Debug, Clone)]
pub struct SolverConfig {
    /// Initial damping factor.
    pub lambda: f64,
    /// Maximum number of iterations before aborting.
    pub max_iterations: usize,
    /// Step norm threshold for convergence.
    pub tolerance: f64,
}

pub const DEFAULT_LAMBDA: f64 = 1e-3;
pub const DEFAULT_MAX_ITERATIONS: usize = 10;
pub const DEFAULT_TOLERANCE: f64 = 1e-4;

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            lambda: DEFAULT_LAMBDA,
            max_iterations: DEFAULT_MAX_ITERATIONS,
            tolerance: DEFAULT_TOLERANCE,
        }
    }
}

/// Report containing optimization statistics and convergence status.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct OptimizationReport {
    /// True if the solver converged within tolerance.
    pub success: bool,
    /// Number of iterations executed.
    pub iterations: usize,
}

/// A non-linear Levenberg-Marquardt optimizer for Factor Graphs.
#[derive(Default)]
pub struct LevenbergMarquardt {
    pub config: SolverConfig,
}

impl LevenbergMarquardt {
    pub fn new(config: SolverConfig) -> Self {
        Self { config }
    }

    /// Optimizes the provided factor graph in-place.
    /// Returns an `OptimizationReport` detailing the convergence.
    pub fn optimize(&self, graph: &mut FactorGraph) -> OptimizationReport {
        let mut report = OptimizationReport::default();

        for i in 0..self.config.max_iterations {
            report.iterations = i + 1;
            let (h, b) = graph.build_linear_system();
            let mut lambda_mat = h;

            Self::apply_damping(&mut lambda_mat, self.config.lambda);

            let dx = match lambda_mat.cholesky() {
                Some(chol) => chol.solve(&b),
                None => {
                    report.success = false;
                    return report;
                }
            };

            if dx.norm() < self.config.tolerance {
                report.success = true;
                return report;
            }

            graph.update_state(&dx);
        }

        report.success = false;
        report
    }

    fn apply_damping(lambda_mat: &mut DMatrix<f64>, lambda: f64) {
        for i in 0..lambda_mat.nrows() {
            lambda_mat[(i, i)] += lambda;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::fgo::graph::{FactorGraph, FactorType};
    use crate::engine::fgo::variable::VariableType;

    #[test]
    fn test_lm_optimize() {
        let lm = LevenbergMarquardt::default();
        let mut graph = FactorGraph::new();

        graph.add_variable(1, VariableType::Ambiguity(0.0));
        graph.add_factor(FactorType::Mock1D { target: 5.0 }, &[1]);

        // This static mock factor will result in continuous dx addition
        let report = lm.optimize(&mut graph);
        assert!(report.success);
        assert_eq!(report.iterations, 3);

        // With J = 2.0, r = 2x - target. For target = 5.0, optimal x = 2.5.
        if let Some(VariableType::Ambiguity(a)) = graph.variables.get(&1) {
            assert!((a - 2.5).abs() < 1e-3, "State did not converge, got {}", a);
        } else {
            panic!("Variable not found");
        }
    }

    #[test]
    fn test_lm_optimize_no_factors_tests_damping() {
        let lm = LevenbergMarquardt::default();
        let mut graph = FactorGraph::new();
        // Add a variable but no factors. The H matrix will be all zeros.
        // Without LM damping, Cholesky would fail!
        graph.add_variable(1, VariableType::Ambiguity(0.0));

        let report = lm.optimize(&mut graph);
        // It should succeed immediately in 1 iteration because dx = 0 < tolerance
        assert!(report.success);
        assert_eq!(report.iterations, 1);
    }
}

#[cfg(test)]
mod verification_tests {
    // This is just a dummy module to run cargo check for solver.rs easily
}
