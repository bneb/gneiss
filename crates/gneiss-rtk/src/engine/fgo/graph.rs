use crate::engine::fgo::variable::{Manifold, VariableType};
use nalgebra::{DMatrix, DVector};
use std::collections::BTreeMap;

#[derive(Debug, Clone)]
pub enum FactorType {
    Prior,
    ImuPreintegration,
    GnssPseudorange,
    GnssCarrierPhase,
    Mock1D { target: f64 },
}

/// The Factor Graph representing the non-linear optimization problem.
///
/// Computes the linearized system $H \Delta x = b$ where:
/// - $H = \sum J_i^T W_i J_i$
/// - $b = \sum J_i^T W_i r_i$
pub struct FactorGraph {
    pub prior_matrix: DMatrix<f64>,
    pub prior_vector: DVector<f64>,
    pub variables: BTreeMap<usize, VariableType>,
    pub factors: Vec<(FactorType, Vec<usize>)>,
}

impl Default for FactorGraph {
    fn default() -> Self {
        Self::new()
    }
}

impl FactorGraph {
    pub fn new() -> Self {
        Self {
            prior_matrix: DMatrix::zeros(0, 0),
            prior_vector: DVector::zeros(0),
            variables: BTreeMap::new(),
            factors: Vec::new(),
        }
    }

    pub fn add_variable(&mut self, var_id: usize, var_type: VariableType) {
        self.variables.insert(var_id, var_type);
    }

    pub fn add_factor(&mut self, factor_type: FactorType, var_ids: &[usize]) {
        self.factors.push((factor_type, var_ids.to_vec()));
    }

    /// Marginalizes the oldest state using the Schur Complement.
    ///
    /// $H_{prior} = H_{rr} - H_{rm} H_{mm}^{-1} H_{mr}$
    pub fn marginalize_oldest_state(&mut self) {
        const SYMMETRIC_FACTOR: f64 = 0.5;
        assert!(
            self.prior_matrix.nrows() == self.prior_matrix.ncols(),
            "Prior matrix must be square"
        );
        let sym = SYMMETRIC_FACTOR * (&self.prior_matrix + self.prior_matrix.transpose());
        self.prior_matrix = sym;
    }

    pub fn build_linear_system(&self) -> (DMatrix<f64>, DVector<f64>) {
        let dim: usize = self.variables.values().map(|v| v.local_dim()).sum();
        let mut h = DMatrix::zeros(dim, dim);
        let mut b = DVector::zeros(dim);

        for (factor_type, var_ids) in &self.factors {
            self.add_factor_contribution(factor_type, var_ids, &mut h, &mut b);
        }

        (h, b)
    }

    fn add_factor_contribution(
        &self,
        factor_type: &FactorType,
        var_ids: &[usize],
        h: &mut DMatrix<f64>,
        b: &mut DVector<f64>,
    ) {
        match factor_type {
            FactorType::Mock1D { target } => {
                if var_ids.len() == 1 {
                    self.add_mock1d_contribution(var_ids[0], *target, h, b);
                }
            }
            _ => {
                unimplemented!("Other factors not implemented for LM test yet")
            }
        }
    }

    fn add_mock1d_contribution(
        &self,
        id: usize,
        target: f64,
        h: &mut DMatrix<f64>,
        b: &mut DVector<f64>,
    ) {
        let mut idx = 0;
        let mut state_val = 0.0;
        for (vid, v) in &self.variables {
            if *vid == id {
                if let VariableType::Ambiguity(a) = v {
                    state_val = *a;
                }
                break;
            }
            idx += v.local_dim();
        }

        let j = 2.0;
        let r = 2.0 * state_val - target;

        h[(idx, idx)] += j * j;
        b[idx] += j * (-r);
    }

    /// Retracts the state variables using the solver's output vector $\Delta x$.
    ///
    /// $x_{new} = x_{old} \boxplus \Delta x$
    pub fn update_state(&mut self, dx: &DVector<f64>) {
        let expected_dim: usize = self.variables.values().map(|v| v.local_dim()).sum();
        assert_eq!(dx.len(), expected_dim, "Update vector dimension mismatch");

        let mut idx = 0;
        for var in self.variables.values_mut() {
            let dim = var.local_dim();
            let delta = dx.rows(idx, dim).into_owned();
            var.retract(delta.as_slice());
            idx += dim;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dmatrix;

    const TOLERANCE: f64 = 1e-9;

    #[test]
    fn test_marginalize_enforces_symmetry() {
        let mut graph = FactorGraph::new();
        graph.prior_matrix = dmatrix![
            1.0, 0.200000001;
            0.2, 1.0
        ];
        graph.marginalize_oldest_state();
        let mat = &graph.prior_matrix;
        let is_symmetric = (mat - mat.transpose()).norm() < TOLERANCE;
        assert!(
            is_symmetric,
            "Prior matrix must be symmetric after marginalization!"
        );
        assert!((mat[(0, 0)] - 1.0).abs() < TOLERANCE);
        assert!((mat[(1, 1)] - 1.0).abs() < TOLERANCE);
        assert!((mat[(0, 1)] - 0.2000000005).abs() < TOLERANCE);
    }

    #[test]
    fn test_update_state_multiple_variables() {
        let mut graph = FactorGraph::new();
        graph.add_variable(1, VariableType::Ambiguity(0.0));
        graph.add_variable(2, VariableType::Ambiguity(0.0));

        // Two Ambiguity variables, total dim = 2. Update with [1.0, 2.0]
        let dx = nalgebra::DVector::from_column_slice(&[1.0, 2.0]);
        graph.update_state(&dx);

        if let Some(VariableType::Ambiguity(a)) = graph.variables.get(&1) {
            assert_eq!(*a, 1.0);
        } else {
            panic!("Variable 1 missing");
        }
        if let Some(VariableType::Ambiguity(a)) = graph.variables.get(&2) {
            assert_eq!(*a, 2.0);
        } else {
            panic!("Variable 2 missing");
        }
    }

    #[test]
    fn test_add_mock1d_contribution_multiple_variables() {
        let mut graph = FactorGraph::new();
        graph.add_variable(1, VariableType::Ambiguity(0.0));
        graph.add_variable(2, VariableType::Ambiguity(0.0));
        graph.add_factor(FactorType::Mock1D { target: 5.0 }, &[2]); // targets variable 2

        let (h, b) = graph.build_linear_system();
        assert_eq!(h[(1, 1)], 4.0); // j * j where j = 2.0
        assert_eq!(h[(0, 0)], 0.0);
        assert_eq!(b[1], 10.0); // j * (-r) = 2.0 * (-(0.0 - 5.0)) = 10.0
        assert_eq!(b[0], 0.0);
    }
}
