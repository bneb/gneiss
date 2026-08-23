use nalgebra::{DMatrix, DVector};

use crate::swfg::variables::{VariableId, VariableValues};

/// A factor (cost term) in the nonlinear least-squares problem.
///
/// The total cost is:  E(x) = Σ ||r_i(x)||²_{W_i}
/// where r_i = residual and W_i = information matrix (inverse covariance).
///
/// The LM solver linearizes around the current estimate:
///   J_i = ∂r_i/∂x  (Jacobian)
///   Builds normal equations: J^T W J Δx = -J^T W r
///
/// Each implementation is responsible for computing its own residual,
/// Jacobian, and information matrix from the variable values.
pub trait Factor: std::fmt::Debug {
    /// The variables this factor connects to, in order.
    /// The order must match the column order in the Jacobian:
    ///   J = [∂r/∂v_0, ∂r/∂v_1, ..., ∂r/∂v_k]
    fn variables(&self) -> &[VariableId];

    /// Residual r(x) = h(x) - z.
    /// Dimension: measurement_dim × 1.
    fn residual(&self, values: &VariableValues) -> DVector<f64>;

    /// Jacobian J = ∂r/∂x evaluated at the current estimate.
    /// Dimension: measurement_dim × total_dim, where total_dim is the
    /// sum of dimensions of all active variables.
    ///
    /// The Jacobian is sparse in the variable dimension: columns for
    /// variables NOT in `self.variables()` are all zeros.
    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64>;

    /// Information matrix (inverse measurement covariance).
    /// Dimension: measurement_dim × measurement_dim.
    fn information(&self) -> DMatrix<f64>;

    /// Optional robust loss threshold.
    /// - `None`: quadratic (L2) loss
    /// - `Some(k)`: Huber loss with threshold k, or Cauchy with scale k
    fn robust_threshold(&self) -> Option<f64> {
        None
    }

    /// Whether the Cauchy loss function should be used instead of Huber.
    /// Cauchy is more aggressive at rejecting large outliers.
    fn use_cauchy(&self) -> bool {
        false
    }
}

/// A factor that penalizes deviation from a fixed prior value.
/// Cost: 1/2 ||x - mu||²_W  where W = information matrix.
#[derive(Debug, Clone)]
pub struct PriorFactor {
    pub variable: VariableId,
    pub mu: DVector<f64>,
    pub information: DMatrix<f64>,
}

impl PriorFactor {
    pub fn new(variable: VariableId, mu: DVector<f64>, variance: f64) -> Self {
        let dim = mu.len();
        let information = DMatrix::identity(dim, dim) / variance.max(1e-12);
        Self {
            variable,
            mu,
            information,
        }
    }
}

impl Factor for PriorFactor {
    fn variables(&self) -> &[VariableId] {
        std::slice::from_ref(&self.variable)
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let x = values.get(self.variable).expect("variable in graph");
        x.into_owned() - &self.mu
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        if let Some((start, dim)) = values.index_of(self.variable) {
            let mut j = DMatrix::zeros(dim, total_dim);
            for i in 0..dim {
                j[(i, start + i)] = 1.0;
            }
            j
        } else {
            DMatrix::zeros(0, total_dim)
        }
    }

    fn information(&self) -> DMatrix<f64> {
        self.information.clone()
    }
}

/// A factor connecting two pose variables p1 and p2 with a relative motion prior.
/// Residual: r = x_{p2} - x_{p1}.
/// Jacobian: ∂r/∂x_{p1} = -I, ∂r/∂x_{p2} = +I.
#[derive(Debug, Clone)]
pub struct RelativePoseFactor {
    pub vars: [VariableId; 2],
    pub information: DMatrix<f64>,
}

impl RelativePoseFactor {
    pub fn new(p1: VariableId, p2: VariableId, variance: f64) -> Self {
        let information = DMatrix::identity(6, 6) / variance.max(1e-12);
        Self {
            vars: [p1, p2],
            information,
        }
    }
}

impl Factor for RelativePoseFactor {
    fn variables(&self) -> &[VariableId] {
        &self.vars
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let v1 = values.get(self.vars[0]).expect("p1 in graph");
        let v2 = values.get(self.vars[1]).expect("p2 in graph");
        v2.into_owned() - v1.into_owned()
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(6, total_dim);
        if let Some((start1, dim1)) = values.index_of(self.vars[0]) {
            for i in 0..dim1.min(6) {
                j[(i, start1 + i)] = -1.0;
            }
        }
        if let Some((start2, dim2)) = values.index_of(self.vars[1]) {
            for i in 0..dim2.min(6) {
                j[(i, start2 + i)] = 1.0;
            }
        }
        j
    }

    fn information(&self) -> DMatrix<f64> {
        self.information.clone()
    }
}

/// A factor that constrains the attitude components of a Pose to zero.
/// This prevents solver explosions when IMU is disabled but the Pose variable is still 6-DOF.
#[derive(Debug, Clone)]
pub struct AttitudePriorFactor {
    pub var_pose: VariableId,
    pub information: DMatrix<f64>,
}

impl AttitudePriorFactor {
    pub fn new(var_pose: VariableId, variance: f64) -> Self {
        Self {
            var_pose,
            information: DMatrix::identity(3, 3) / variance.max(1e-12),
        }
    }
}

impl Factor for AttitudePriorFactor {
    fn variables(&self) -> &[VariableId] {
        std::slice::from_ref(&self.var_pose)
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let pose = values.get(self.var_pose).expect("pose in graph");
        DVector::from_row_slice(&[pose[3], pose[4], pose[5]])
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(3, total_dim);
        let (start, _) = values.index_of(self.var_pose).expect("pose in graph");
        
        j[(0, start + 3)] = 1.0;
        j[(1, start + 4)] = 1.0;
        j[(2, start + 5)] = 1.0;
        
        j
    }

    fn information(&self) -> DMatrix<f64> {
        self.information.clone()
    }
}

/// Helper to compute M-estimator robust weight for a residual norm.
pub fn compute_robust_weight(r_norm: f64, k: f64, use_cauchy: bool) -> f64 {
    if use_cauchy {
        1.0 / (1.0 + (r_norm / k).powi(2))
    } else if r_norm <= k {
        1.0
    } else {
        k / r_norm
    }
}

/// Helper to compute M-estimator robust loss/error for a residual norm.
pub fn compute_robust_error(r_norm: f64, k: f64, use_cauchy: bool) -> f64 {
    let s = r_norm * r_norm;
    if use_cauchy {
        let k2 = k * k;
        k2 * (1.0 + s / k2).ln()
    } else if r_norm <= k {
        s
    } else {
        2.0 * k * r_norm - k * k
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::swfg::variables::{VariableKind, VariableNode};
    use std::collections::BTreeMap;

    #[test]
    fn robust_weight_huber_within_threshold_is_one() {
        assert_eq!(compute_robust_weight(2.0, 3.0, false), 1.0);
    }

    #[test]
    fn robust_weight_huber_above_threshold_downweights() {
        let w = compute_robust_weight(6.0, 3.0, false);
        assert!((w - 0.5).abs() < 1e-12);
    }

    #[test]
    fn robust_weight_cauchy_downweights_squared() {
        // Cauchy with k=3 at r=3 => 1 / (1 + 1) = 0.5
        let w = compute_robust_weight(3.0, 3.0, true);
        assert!((w - 0.5).abs() < 1e-12);
    }

    #[test]
    fn prior_factor_residual_is_zero_at_mu() {
        let kind = VariableKind::Pose { epoch: 0 };
        let id = VariableId::new(0);
        let mut vars = BTreeMap::new();
        let mut node = VariableNode::new(id, kind);
        node.value.copy_from(&DVector::from_vec(vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0]));
        vars.insert(id, node);

        let values = VariableValues::build(&vars);
        let factor = PriorFactor {
            variable: id,
            mu: DVector::from_vec(vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0]),
            information: DMatrix::identity(6, 6),
        };

        let r = factor.residual(&values);
        assert!(r.norm() < 1e-12, "residual at prior mean should be zero");
    }

    #[test]
    fn prior_factor_residual_is_nonzero_away_from_mu() {
        let kind = VariableKind::Velocity { epoch: 0 };
        let id = VariableId::new(0);
        let mut vars = BTreeMap::new();
        let mut node = VariableNode::new(id, kind);
        node.value.copy_from(&DVector::from_vec(vec![0.0, 0.0, 0.0]));
        vars.insert(id, node);

        let values = VariableValues::build(&vars);
        let factor = PriorFactor {
            variable: id,
            mu: DVector::from_vec(vec![1.0, 0.0, 0.0]),
            information: DMatrix::identity(3, 3),
        };

        let r = factor.residual(&values);
        assert!((r.norm() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn relative_pose_factor_residual_and_jacobian() {
        let p1 = VariableId::new(0);
        let p2 = VariableId::new(1);
        let mut vars = BTreeMap::new();

        let mut node1 = VariableNode::new(p1, VariableKind::Pose { epoch: 0 });
        node1.value.copy_from(&DVector::from_vec(vec![10.0, 20.0, 30.0, 0.0, 0.0, 0.0]));
        vars.insert(p1, node1);

        let mut node2 = VariableNode::new(p2, VariableKind::Pose { epoch: 1 });
        node2.value.copy_from(&DVector::from_vec(vec![12.0, 20.0, 30.0, 0.0, 0.0, 0.0]));
        vars.insert(p2, node2);

        let values = VariableValues::build(&vars);
        let rel_factor = RelativePoseFactor::new(p1, p2, 1.0);

        let res = rel_factor.residual(&values);
        assert_eq!(res.len(), 6);
        assert!((res[0] - 2.0).abs() < 1e-12);

        let j = rel_factor.jacobian(&values);
        assert_eq!(j.nrows(), 6);
        assert_eq!(j.ncols(), 12);
        assert_eq!(j[(0, 0)], -1.0);
        assert_eq!(j[(0, 6)], 1.0);
    }
}
