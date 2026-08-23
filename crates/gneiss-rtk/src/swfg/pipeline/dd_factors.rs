//! Double-differenced RTK factor graph factors.

use nalgebra::{DMatrix, DVector, Vector3};

use crate::swfg::factor::Factor;
use crate::swfg::variables::{VariableId, VariableValues};

/// Compute elevation-dependent pseudorange variance.
/// Model: σ = a + b / sin(el), matching RTKLIB's approach.
pub fn elevation_pr_variance(el_rad: f64) -> f64 {
    let sin_el = el_rad.sin().max(0.087); // clamp at ~5°
    let sigma = 0.3 + 0.3 / sin_el;
    sigma * sigma
}

/// Compute elevation-dependent carrier phase variance.
pub fn elevation_cp_variance(el_rad: f64) -> f64 {
    let sin_el = el_rad.sin().max(0.087);
    let sigma = 0.003 + 0.003 / sin_el;
    sigma * sigma
}

/// Double-differenced pseudorange factor for RTK.
///
/// Residual: DD_PR_obs - (DD_geometric_range)
/// where DD = (rover_sat - rover_ref) - (base_sat - base_ref)
#[derive(Clone)]
pub struct DdPseudorangeFactor {
    pub var_pose: VariableId,
    /// DD pseudorange observation (meters).
    pub dd_pr_obs: f64,
    /// Satellite position ECEF.
    pub sat_pos: Vector3<f64>,
    /// Reference satellite position ECEF.
    pub ref_pos: Vector3<f64>,
    /// Base station position ECEF.
    pub base_pos: Vector3<f64>,
    /// Base-to-sat range minus base-to-ref range (precomputed).
    pub base_dd_range: f64,
    /// Measurement variance (m²).
    pub variance_m2: f64,
    /// Satellite elevation (radians).
    pub elevation_rad: f64,
    /// Reference satellite elevation (radians).
    pub ref_elevation_rad: f64,
    pub variables: Vec<VariableId>,
}

impl std::fmt::Debug for DdPseudorangeFactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DdPseudorangeFactor").finish()
    }
}

impl Factor for DdPseudorangeFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let pose = match values.get(self.var_pose) {
            Some(p) => p,
            None => return DVector::zeros(1),
        };
        let rx_pos = Vector3::new(pose[0], pose[1], pose[2]);
        let rover_sat_range = (self.sat_pos - rx_pos).norm();
        let rover_ref_range = (self.ref_pos - rx_pos).norm();
        let rover_dd_range = rover_sat_range - rover_ref_range;
        let predicted = rover_dd_range - self.base_dd_range;
        DVector::from_element(1, self.dd_pr_obs - predicted)
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(1, total_dim);
        let start = match values.index_of(self.var_pose) {
            Some((s, _)) => s,
            None => return j,
        };
        let pose = match values.get(self.var_pose) {
            Some(p) => p,
            None => return j,
        };
        let rx_pos = Vector3::new(pose[0], pose[1], pose[2]);
        let los_sat = (self.sat_pos - rx_pos).normalize();
        let los_ref = (self.ref_pos - rx_pos).normalize();
        for k in 0..3 {
            j[(0, start + k)] = los_sat[k] - los_ref[k];
        }
        j
    }

    fn information(&self) -> DMatrix<f64> {
        let var_dd = elevation_pr_variance(self.elevation_rad)
            + elevation_pr_variance(self.ref_elevation_rad);
        DMatrix::from_element(1, 1, 1.0 / var_dd.max(1e-4))
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(10.0)
    }

    fn use_cauchy(&self) -> bool {
        true
    }
}

/// Double-differenced carrier phase factor for RTK.
///
/// Residual: DD_CP_obs - (DD_geometric_range + lambda * DD_N)
#[derive(Clone)]
pub struct DdCarrierPhaseFactor {
    pub var_pose: VariableId,
    pub var_amb: VariableId,
    /// DD carrier phase observation (meters).
    pub dd_cp_obs_m: f64,
    /// Satellite position ECEF.
    pub sat_pos: Vector3<f64>,
    /// Reference satellite position ECEF.
    pub ref_pos: Vector3<f64>,
    /// Base station position ECEF.
    pub base_pos: Vector3<f64>,
    /// Base DD range (precomputed).
    pub base_dd_range: f64,
    /// Carrier wavelength (meters).
    pub lambda: f64,
    /// Measurement variance (m²).
    pub variance_m2: f64,
    /// Satellite elevation (radians).
    pub elevation_rad: f64,
    /// Reference satellite elevation (radians).
    pub ref_elevation_rad: f64,
    pub is_new_amb: bool,
    pub variables: Vec<VariableId>,
}

impl std::fmt::Debug for DdCarrierPhaseFactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DdCarrierPhaseFactor").finish()
    }
}

impl Factor for DdCarrierPhaseFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let pose = match values.get(self.var_pose) {
            Some(p) => p,
            None => return DVector::zeros(1),
        };
        let amb = match values.get(self.var_amb) {
            Some(a) => a[0],
            None => return DVector::zeros(1),
        };
        let rx_pos = Vector3::new(pose[0], pose[1], pose[2]);
        let rover_dd = (self.sat_pos - rx_pos).norm() - (self.ref_pos - rx_pos).norm();
        let predicted = rover_dd - self.base_dd_range + self.lambda * amb;
        DVector::from_element(1, self.dd_cp_obs_m - predicted)
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(1, total_dim);
        let s_pose = match values.index_of(self.var_pose) {
            Some((s, _)) => s,
            None => return j,
        };
        let s_amb = match values.index_of(self.var_amb) {
            Some((s, _)) => s,
            None => return j,
        };
        let pose = match values.get(self.var_pose) {
            Some(p) => p,
            None => return j,
        };
        let rx_pos = Vector3::new(pose[0], pose[1], pose[2]);
        let los_sat = (self.sat_pos - rx_pos).normalize();
        let los_ref = (self.ref_pos - rx_pos).normalize();
        for k in 0..3 {
            j[(0, s_pose + k)] = los_sat[k] - los_ref[k];
        }
        j[(0, s_amb)] = -self.lambda;
        j
    }

    fn information(&self) -> DMatrix<f64> {
        let var_dd = elevation_cp_variance(self.elevation_rad)
            + elevation_cp_variance(self.ref_elevation_rad);
        DMatrix::from_element(1, 1, 1.0 / var_dd.max(1e-4))
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(0.50)
    }

    fn use_cauchy(&self) -> bool {
        true
    }
}

/// Double-differenced Doppler velocity factor between consecutive poses.
///
/// Residual: (P_curr - P_prev) . (e_sat - e_ref) - DD_Doppler_m_s * dt
#[derive(Clone)]
pub struct DdDopplerFactor {
    pub var_pose_prev: VariableId,
    pub var_pose_curr: VariableId,
    /// DD Doppler velocity (m/s).
    pub dd_doppler_m_s: f64,
    /// Time delta between epochs (seconds).
    pub dt: f64,
    /// Satellite position ECEF.
    pub sat_pos: Vector3<f64>,
    /// Reference satellite position ECEF.
    pub ref_pos: Vector3<f64>,
    /// Measurement variance (m²).
    pub variance_m2: f64,
    pub variables: Vec<VariableId>,
}

impl std::fmt::Debug for DdDopplerFactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DdDopplerFactor").finish()
    }
}

impl Factor for DdDopplerFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let pose_prev = match values.get(self.var_pose_prev) {
            Some(p) => p,
            None => return DVector::zeros(1),
        };
        let pose_curr = match values.get(self.var_pose_curr) {
            Some(p) => p,
            None => return DVector::zeros(1),
        };
        let pos_prev = Vector3::new(pose_prev[0], pose_prev[1], pose_prev[2]);
        let pos_curr = Vector3::new(pose_curr[0], pose_curr[1], pose_curr[2]);

        let delta_pos = pos_curr - pos_prev;
        let los_sat = (self.sat_pos - pos_curr).normalize();
        let los_ref = (self.ref_pos - pos_curr).normalize();
        let dd_los = los_sat - los_ref;

        let proj_disp = delta_pos.dot(&dd_los);
        let expected_disp = self.dd_doppler_m_s * self.dt;

        DVector::from_element(1, proj_disp - expected_disp)
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(1, total_dim);
        let s_prev = match values.index_of(self.var_pose_prev) {
            Some((s, _)) => s,
            None => return j,
        };
        let s_curr = match values.index_of(self.var_pose_curr) {
            Some((s, _)) => s,
            None => return j,
        };

        let pose_curr = match values.get(self.var_pose_curr) {
            Some(p) => p,
            None => return j,
        };
        let pos_curr = Vector3::new(pose_curr[0], pose_curr[1], pose_curr[2]);
        let los_sat = (self.sat_pos - pos_curr).normalize();
        let los_ref = (self.ref_pos - pos_curr).normalize();
        let dd_los = los_sat - los_ref;

        for k in 0..3 {
            j[(0, s_prev + k)] = -dd_los[k];
            j[(0, s_curr + k)] = dd_los[k];
        }
        j
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::from_element(1, 1, 1.0 / self.variance_m2.max(1e-4))
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(0.20)
    }

    fn use_cauchy(&self) -> bool {
        true
    }
}

/// A constraint factor to enforce N1 - N2 = N_WL during widelane bootstrapping.
pub struct WidelaneConstraintFactor {
    pub var_amb1: VariableId,
    pub var_amb2: VariableId,
    pub fixed_n_wl: f64,
    pub variance: f64,
    pub variables: Vec<VariableId>,
}

impl std::fmt::Debug for WidelaneConstraintFactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("WidelaneConstraintFactor")
            .field("fixed_n_wl", &self.fixed_n_wl)
            .finish()
    }
}

impl Factor for WidelaneConstraintFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let n1 = match values.get(self.var_amb1) {
            Some(v) => v[0],
            None => return DVector::zeros(1),
        };
        let n2 = match values.get(self.var_amb2) {
            Some(v) => v[0],
            None => return DVector::zeros(1),
        };
        DVector::from_element(1, n1 - n2 - self.fixed_n_wl)
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(1, total_dim);
        
        let start_n1 = match values.index_of(self.var_amb1) {
            Some((s, _)) => s,
            None => return j,
        };
        let start_n2 = match values.index_of(self.var_amb2) {
            Some((s, _)) => s,
            None => return j,
        };
        
        j[(0, start_n1)] = 1.0;
        j[(0, start_n2)] = -1.0;
        j
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::from_element(1, 1, 1.0 / self.variance)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_elevation_variances() {
        let el = std::f64::consts::FRAC_PI_2;
        let v_pr = elevation_pr_variance(el);
        let v_cp = elevation_cp_variance(el);
        assert!(v_pr > 0.0);
        assert!(v_cp > 0.0);
        assert!(v_cp < v_pr);
    }
}
