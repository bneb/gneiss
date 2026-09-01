use nalgebra::{DMatrix, DVector};
use std::collections::{BTreeMap, VecDeque};

use crate::swfg::factor::Factor;
use crate::swfg::variables::{VariableId, VariableKind, VariableNode};

/// Metadata for one epoch in the sliding window.
#[derive(Debug, Clone)]
pub struct EpochMetadata {
    /// GPS time of this epoch.
    pub time: f64,
    /// Variables created at this epoch.  Used during marginalization to
    /// identify which variables to remove when this epoch slides out.
    pub variables: Vec<VariableId>,
    /// Number of satellite observations in this epoch (for diagnostics).
    pub n_satellites: usize,
}

/// A dense prior factor produced by Schur-complement marginalization
/// of old epochs.  Encodes all information from marginalized variables
/// as a quadratic cost on the remaining (active) variables.
///
/// Cost: 1/2 * x^T H x - g^T x + const
/// where H = hessian, g = gradient.
#[derive(Debug, Clone)]
pub struct MarginalPriorFactor {
    /// Variables this prior connects to (subset of active variables).
    /// Stores both the VariableId and its dimension.
    pub variables: Vec<(VariableId, usize)>,
    /// Hessian (information) matrix.  Dimension: total_dim × total_dim
    /// where total_dim = sum of dimensions of `variables`.
    pub hessian: DMatrix<f64>,
    /// Information-weighted residual vector (gradient term).
    /// Dimension: total_dim × 1.
    pub gradient: DVector<f64>,
    /// The state of the variables at the time of marginalization.
    /// Used to compute the gradient penalty H * (x - x0)
    pub x0: DVector<f64>,
}

/// The core factor graph data structure.
///
/// Variables are stored in a `BTreeMap` ordered by `VariableId`, which
/// provides deterministic iteration order — important for reproducible
/// solves and tests.
///
/// The graph supports adding variables and factors, running LM
/// optimization, and Schur-complement marginalization of old epochs.
#[derive(Debug)]
pub struct EstimationGraph {
    /// All active variables, ordered by VariableId.
    pub variables: BTreeMap<VariableId, VariableNode>,
    /// All active factors.
    pub factors: Vec<Box<dyn Factor>>,
    /// Dense prior from marginalized epochs (None if no marginalization
    /// has occurred yet).
    pub marginal_prior: Option<MarginalPriorFactor>,
    /// Epoch metadata for the current window (oldest first).
    pub window: VecDeque<EpochMetadata>,
    /// Next VariableId to assign (monotonically increasing).
    next_id: u64,
}

impl EstimationGraph {
    pub fn new() -> Self {
        Self {
            variables: BTreeMap::new(),
            factors: Vec::new(),
            marginal_prior: None,
            window: VecDeque::new(),
            next_id: 0,
        }
    }

    /// Allocate a new variable with the given kind.  Returns its
    /// unique `VariableId`.  The variable starts with a zero value;
    /// the caller should set the initial estimate if known (e.g.,
    /// from SPP for position, from previous epoch for ambiguities).
    pub fn add_variable(&mut self, kind: VariableKind) -> VariableId {
        let id = VariableId::new(self.next_id);
        self.next_id += 1;
        self.variables.insert(id, VariableNode::new(id, kind));
        id
    }

    /// Add a factor to the graph.
    pub fn add_factor(&mut self, factor: Box<dyn Factor>) {
        self.factors.push(factor);
    }

    /// Remove a variable and all factors that reference it.  Used when
    /// an epoch slides out of the window and its per-epoch variables
    /// are either marginalized or discarded.
    pub fn remove_variable(&mut self, id: VariableId) {
        self.variables.remove(&id);
        self.factors.retain(|f| !f.variables().contains(&id));
    }

    /// Set a variable's value.  Used to initialize variables from external
    /// sources (e.g., SPP position, known base station coordinates).
    pub fn set_value(&mut self, id: VariableId, value: &[f64]) {
        if let Some(node) = self.variables.get_mut(&id) {
            node.set_value(value);
        }
    }

    /// Clear all factors (e.g., before rebuilding for a new epoch).
    pub fn clear_factors(&mut self) {
        self.factors.clear();
    }

    /// Number of active variables.
    pub fn n_variables(&self) -> usize {
        self.variables.len()
    }

    /// Number of active factors.
    pub fn n_factors(&self) -> usize {
        self.factors.len()
    }

    /// Total dimension of the state vector.
    pub fn total_dim(&self) -> usize {
        self.variables.values().map(|v| v.value.len()).sum()
    }

    /// Validate graph integrity: returns an error if any active variable is unconnected (orphan).
    pub fn validate_graph_structure(&self) -> Result<(), String> {
        let mut connected = std::collections::HashSet::new();
        for factor in &self.factors {
            for &v in factor.variables() {
                connected.insert(v);
            }
        }
        if let Some(ref prior) = self.marginal_prior {
            for &(v, _) in &prior.variables {
                connected.insert(v);
            }
        }
        for (&id, node) in &self.variables {
            if !connected.contains(&id) {
                return Err(format!("Orphan variable in graph: ID {} ({:?})", id.as_u64(), node.kind));
            }
        }
        Ok(())
    }
}

impl Default for EstimationGraph {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::swfg::factor::PriorFactor;

    #[test]
    fn graph_add_variable_assigns_sequential_ids() {
        let mut g = EstimationGraph::new();
        let id0 = g.add_variable(VariableKind::Pose { epoch: 0 });
        let id1 = g.add_variable(VariableKind::Velocity { epoch: 0 });
        assert_eq!(id0.as_u64(), 0);
        assert_eq!(id1.as_u64(), 1);
        assert_eq!(g.n_variables(), 2);
    }

    #[test]
    fn graph_add_factor_increments_count() {
        let mut g = EstimationGraph::new();
        let id = g.add_variable(VariableKind::TropoZwd { epoch: 0 });
        let factor = PriorFactor {
            variable: id,
            mu: DVector::from_element(1, 0.1),
            information: DMatrix::from_element(1, 1, 100.0),
        };
        g.add_factor(Box::new(factor));
        assert_eq!(g.n_factors(), 1);
    }

    #[test]
    fn graph_remove_variable_removes_factors() {
        let mut g = EstimationGraph::new();
        let id = g.add_variable(VariableKind::Ambiguity { satellite: 1, frequency: 1, arc: 0 });
        let factor = PriorFactor {
            variable: id,
            mu: DVector::from_element(1, 0.0),
            information: DMatrix::from_element(1, 1, 1.0),
        };
        g.add_factor(Box::new(factor));
        g.remove_variable(id);
        assert_eq!(g.n_variables(), 0);
        assert_eq!(g.n_factors(), 0);
    }

    #[test]
    fn graph_total_dim_sums_variable_dims() {
        let mut g = EstimationGraph::new();
        g.add_variable(VariableKind::Pose { epoch: 0 });        // 6
        g.add_variable(VariableKind::Velocity { epoch: 0 });    // 3
        g.add_variable(VariableKind::Ambiguity { satellite: 1, frequency: 1, arc: 0 }); // 1
        assert_eq!(g.total_dim(), 10);
    }

    #[test]
    fn graph_validate_structure_detects_orphan_variable() {
        let mut g = EstimationGraph::new();
        let id0 = g.add_variable(VariableKind::Pose { epoch: 0 });
        let _id1 = g.add_variable(VariableKind::Velocity { epoch: 0 }); // Orphan!
        g.add_factor(Box::new(PriorFactor {
            variable: id0,
            mu: DVector::from_element(6, 0.0),
            information: DMatrix::identity(6, 6),
        }));
        assert!(g.validate_graph_structure().is_err());
    }
}
