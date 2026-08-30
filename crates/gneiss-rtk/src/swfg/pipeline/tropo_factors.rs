//! Tropospheric Zenith Wet Delay (ZWD) and Horizontal Gradient Factor Graph Constraints.

use nalgebra::{DMatrix, DVector, Vector2};
use crate::swfg::factor::Factor;
use crate::swfg::variables::{VariableId, VariableValues};

/// Random-walk temporal factor between consecutive ZWD states:
///   r = \text{ZWD}_{k+1} - \text{ZWD}_k \sim \mathcal{N}(0, \sigma^2)
#[derive(Debug, Clone)]
pub struct ZwdRandomWalkFactor {
    pub var_zwd_prev: VariableId,
    pub var_zwd_curr: VariableId,
    pub variance: f64, // sigma^2
    pub variables: Vec<VariableId>,
}

impl ZwdRandomWalkFactor {
    pub fn new(var_zwd_prev: VariableId, var_zwd_curr: VariableId, dt_s: f64, q_zwd_m_per_sqrt_s: f64) -> Self {
        let sigma = (q_zwd_m_per_sqrt_s * libm::sqrt(dt_s.max(0.1))).max(1e-5);
        let variance = sigma * sigma;
        Self {
            var_zwd_prev,
            var_zwd_curr,
            variance,
            variables: vec![var_zwd_prev, var_zwd_curr],
        }
    }
}

impl Factor for ZwdRandomWalkFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let zwd_prev = values.get(self.var_zwd_prev).map(|v| v[0]).unwrap_or(0.0);
        let zwd_curr = values.get(self.var_zwd_curr).map(|v| v[0]).unwrap_or(0.0);
        DVector::from_element(1, zwd_curr - zwd_prev)
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(1, total_dim);

        if let Some((s_prev, _)) = values.index_of(self.var_zwd_prev) {
            j[(0, s_prev)] = -1.0;
        }
        if let Some((s_curr, _)) = values.index_of(self.var_zwd_curr) {
            j[(0, s_curr)] = 1.0;
        }

        j
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::from_element(1, 1, 1.0 / self.variance)
    }
}

/// Computes the Chen & Herring (1992) gradient mapping function:
///   m_{\text{grad}}(e) = \frac{1}{\sin(e) \tan(e) + 0.0032}
#[inline]
pub fn gradient_mapping_function(elev_rad: f64) -> f64 {
    let sin_e = libm::sin(elev_rad.max(0.05));
    let tan_e = libm::tan(elev_rad.max(0.05));
    1.0 / (sin_e * tan_e + 0.0032)
}

/// Computes the slant tropospheric gradient delay:
///   \Delta \tau = m_{\text{grad}}(e) \cdot (G_N \cos A + G_E \sin A)
#[inline]
pub fn slant_gradient_delay(elev_rad: f64, azim_rad: f64, gradients_ne: Vector2<f64>) -> f64 {
    let m_g = gradient_mapping_function(elev_rad);
    let dir = gradients_ne.x * libm::cos(azim_rad) + gradients_ne.y * libm::sin(azim_rad);
    m_g * dir
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;
    use crate::swfg::variables::{VariableKind, VariableNode};

    #[test]
    fn test_zwd_random_walk_factor_residual_and_jacobian() {
        let var1 = VariableId::new(0);
        let var2 = VariableId::new(1);
        let factor = ZwdRandomWalkFactor::new(var1, var2, 30.0, 1e-4);

        let mut map = BTreeMap::new();
        let mut node1 = VariableNode::new(var1, VariableKind::TropoZwd { epoch: 0 });
        node1.value.copy_from(&DVector::from_vec(vec![0.150]));
        map.insert(var1, node1);

        let mut node2 = VariableNode::new(var2, VariableKind::TropoZwd { epoch: 1 });
        node2.value.copy_from(&DVector::from_vec(vec![0.152]));
        map.insert(var2, node2);

        let values = VariableValues::build(&map);

        let res = factor.residual(&values);
        assert!((res[0] - 0.002).abs() < 1e-9);

        let j = factor.jacobian(&values);
        assert_eq!(j.nrows(), 1);
        assert_eq!(j[(0, 0)], -1.0);
        assert_eq!(j[(0, 1)], 1.0);
    }
}
