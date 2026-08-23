//! IMU bias factors for sliding-window factor graph optimization.

use nalgebra::{DMatrix, DVector, Vector3};

use crate::swfg::factor::Factor;
use crate::swfg::variables::{VariableId, VariableValues};

/// Prior and random walk factor for 6-DOF IMU bias [ba_x, ba_y, ba_z, bg_x, bg_y, bg_z].
///
/// Residual:
///   r(bias) = bias - nominal_bias
#[derive(Clone, Debug)]
pub struct ImuBiasPriorFactor {
    pub var_bias: VariableId,
    pub nominal_bias: Vector3<f64>,
    pub nominal_gyro_bias: Vector3<f64>,
    pub accel_variance: f64,
    pub gyro_variance: f64,
    pub variables: Vec<VariableId>,
}

impl ImuBiasPriorFactor {
    /// Create a new IMU bias prior factor.
    pub fn new(
        var_bias: VariableId,
        nominal_bias: Vector3<f64>,
        nominal_gyro_bias: Vector3<f64>,
        accel_variance: f64,
        gyro_variance: f64,
    ) -> Self {
        Self {
            var_bias,
            nominal_bias,
            nominal_gyro_bias,
            accel_variance,
            gyro_variance,
            variables: vec![var_bias],
        }
    }
}

impl Factor for ImuBiasPriorFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let Some(bias) = values.get(self.var_bias) else {
            return DVector::zeros(6);
        };
        if bias.len() < 6 || !bias.iter().all(|x| x.is_finite()) {
            return DVector::zeros(6);
        }

        let ba = Vector3::new(bias[0], bias[1], bias[2]);
        let bg = Vector3::new(bias[3], bias[4], bias[5]);
        let d_ba = ba - self.nominal_bias;
        let d_bg = bg - self.nominal_gyro_bias;

        DVector::from_vec(vec![d_ba.x, d_ba.y, d_ba.z, d_bg.x, d_bg.y, d_bg.z])
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(6, total_dim);
        let Some((s_bias, _)) = values.index_of(self.var_bias) else {
            return j;
        };

        for i in 0..6 {
            j[(i, s_bias + i)] = 1.0;
        }
        j
    }

    fn information(&self) -> DMatrix<f64> {
        let mut info = DMatrix::zeros(6, 6);
        let var_a = if self.accel_variance > 1e-12 { self.accel_variance } else { 1e-4 };
        let var_g = if self.gyro_variance > 1e-12 { self.gyro_variance } else { 1e-6 };

        for i in 0..3 {
            info[(i, i)] = 1.0 / var_a;
        }
        for i in 3..6 {
            info[(i, i)] = 1.0 / var_g;
        }
        info
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::swfg::graph::EstimationGraph;
    use crate::swfg::variables::VariableKind;

    #[test]
    fn test_imu_bias_prior_factor_residual_and_jacobian() {
        let mut graph = EstimationGraph::new();
        let var_bias = graph.add_variable(VariableKind::ImuBias);
        graph.set_value(var_bias, &[0.05, -0.02, 0.01, 0.001, -0.002, 0.003]);

        let factor = ImuBiasPriorFactor::new(
            var_bias,
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0),
            0.04,
            0.0001,
        );

        let values = VariableValues::build(&graph.variables);
        let res = factor.residual(&values);
        assert_eq!(res.len(), 6);
        assert!((res[0] - 0.05).abs() < 1e-6);
        assert!((res[1] - (-0.02)).abs() < 1e-6);
        assert!((res[3] - 0.001).abs() < 1e-6);

        let jac = factor.jacobian(&values);
        assert_eq!(jac.nrows(), 6);
        assert_eq!(jac.ncols(), 6);
        for i in 0..6 {
            assert!((jac[(i, i)] - 1.0).abs() < 1e-6);
        }

        let info = factor.information();
        assert!((info[(0, 0)] - 1.0 / 0.04).abs() < 1e-6);
        assert!((info[(3, 3)] - 1.0 / 0.0001).abs() < 1e-6);
    }
}
