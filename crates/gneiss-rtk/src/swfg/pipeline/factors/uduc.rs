//! Undifferenced Uncombined (UDUC) observation factors for Precise Point Positioning (PPP-AR).
//!
//! Models raw uncombined pseudorange and carrier-phase observations per frequency,
//! incorporating all six geodetic normalizations:
//! 1. Relativistic orbit eccentricity & Shapiro gravitational time delay
//! 2. Satellite Attitude & ANTEX 3D Phase Center Offset (PCO) ECEF projection
//! 3. Continuous Wu (1993) RHCP Phase Windup
//! 4. IERS Conventions (2010) Solid Earth Tides
//! 5. IERS 11-constituent Ocean Tide Loading (OTL)
//! 6. Slant Tropospheric Hydrostatic/Wet mapping + Horizontal East/North Gradients

use nalgebra::{DMatrix, DVector, Vector3};

use crate::swfg::factor::Factor;
use crate::swfg::pipeline::passes::CorrectedObservation;
use crate::swfg::variables::{VariableId, VariableValues};

/// Physical geodetic normalization package for PPP-AR theoretical range computation.
#[derive(Debug, Clone, Copy, Default)]
pub struct GeodeticNormalizations {
    /// Relativistic periodic clock correction in range units (-\frac{2 \mathbf{r}\cdot\mathbf{v}}{c}).
    pub sat_relativity_m: f64,
    /// Gravitational Shapiro time delay in range units.
    pub shapiro_delay_m: f64,
    /// Satellite Antenna Phase Center Offset (PCO) projected into ECEF (meters).
    pub sat_pco_ecef_m: Vector3<f64>,
    /// IERS 2010 Solid Earth Tide 3D displacement vector in ECEF (meters).
    pub solid_earth_tide_m: Vector3<f64>,
    /// Ocean Tide Loading (OTL) 3D displacement vector in ECEF (meters).
    pub ocean_tide_loading_m: Vector3<f64>,
}

/// Undifferenced Uncombined pseudorange factor.
#[derive(Debug, Clone)]
pub struct UducPseudorangeFactor {
    pub obs: CorrectedObservation,
    pub epoch: u32,
    pub raw_pr_m: f64,
    pub var_pose: VariableId,
    pub var_clock: VariableId,
    pub var_zwd: Option<VariableId>,
    pub var_iono: Option<VariableId>,
    pub freq_ratio_sq: f64, // \mu_i = (f_1 / f_i)^2
    pub norm: GeodeticNormalizations,
    pub variables: Vec<VariableId>,
}

impl UducPseudorangeFactor {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        obs: CorrectedObservation,
        epoch: u32,
        raw_pr_m: f64,
        var_pose: VariableId,
        var_clock: VariableId,
        var_zwd: Option<VariableId>,
        var_iono: Option<VariableId>,
        freq_ratio_sq: f64,
        norm: GeodeticNormalizations,
    ) -> Self {
        let mut vars = vec![var_pose, var_clock];
        if let Some(z) = var_zwd {
            vars.push(z);
        }
        if let Some(i) = var_iono {
            vars.push(i);
        }
        Self {
            obs,
            epoch,
            raw_pr_m,
            var_pose,
            var_clock,
            var_zwd,
            var_iono,
            freq_ratio_sq,
            norm,
            variables: vars,
        }
    }
}

impl Factor for UducPseudorangeFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let rx_pos = values
            .get(self.var_pose)
            .map(|p| Vector3::new(p[0], p[1], p[2]) + self.norm.solid_earth_tide_m + self.norm.ocean_tide_loading_m)
            .unwrap_or_else(Vector3::zeros);

        let c_dt = values.get(self.var_clock).map(|c| c[0]).unwrap_or(0.0);
        let zwd = self.var_zwd.and_then(|z| values.get(z)).map(|v| v[0]).unwrap_or(0.0);
        let iono = self.var_iono.and_then(|i| values.get(i)).map(|v| v[0]).unwrap_or(0.0);

        let sat_apc = self.obs.sat_pos_ecef + self.norm.sat_pco_ecef_m;
        let los = sat_apc - rx_pos;
        let geometric_range = los.norm();

        let modeled_tropo = self.obs.tropo_dry_m + self.obs.tropo_map_wet * zwd;
        let modeled_iono = self.freq_ratio_sq * iono;
        let modeled_pr = geometric_range + c_dt - self.obs.sat_clock_m
            + self.norm.sat_relativity_m + self.norm.shapiro_delay_m
            + modeled_tropo + modeled_iono;

        DVector::from_element(1, self.raw_pr_m - modeled_pr)
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(1, total_dim);

        let rx_pos = values
            .get(self.var_pose)
            .map(|p| Vector3::new(p[0], p[1], p[2]) + self.norm.solid_earth_tide_m + self.norm.ocean_tide_loading_m)
            .unwrap_or_else(Vector3::zeros);

        let sat_apc = self.obs.sat_pos_ecef + self.norm.sat_pco_ecef_m;
        let los = sat_apc - rx_pos;
        let range = los.norm().max(1.0);
        let unit_los = los / range;

        if let Some((s_pose, _)) = values.index_of(self.var_pose) {
            j[(0, s_pose)] = unit_los.x;
            j[(0, s_pose + 1)] = unit_los.y;
            j[(0, s_pose + 2)] = unit_los.z;
        }
        if let Some((s_clock, _)) = values.index_of(self.var_clock) {
            j[(0, s_clock)] = -1.0;
        }
        if let Some(z_id) = self.var_zwd {
            if let Some((s_zwd, _)) = values.index_of(z_id) {
                j[(0, s_zwd)] = -self.obs.tropo_map_wet;
            }
        }
        if let Some(i_id) = self.var_iono {
            if let Some((s_iono, _)) = values.index_of(i_id) {
                j[(0, s_iono)] = -self.freq_ratio_sq;
            }
        }

        j
    }

    fn information(&self) -> DMatrix<f64> {
        let var = (self.obs.variance_m2).max(1e-4);
        DMatrix::from_element(1, 1, 1.0 / var)
    }
}

/// Undifferenced Uncombined carrier-phase factor.
#[derive(Debug, Clone)]
pub struct UducCarrierPhaseFactor {
    pub obs: CorrectedObservation,
    pub epoch: u32,
    pub raw_cp_cycles: f64,
    pub var_pose: VariableId,
    pub var_clock: VariableId,
    pub var_zwd: Option<VariableId>,
    pub var_iono: Option<VariableId>,
    pub var_amb: VariableId,
    pub wavelength_m: f64,
    pub freq_ratio_sq: f64,
    pub windup_m: f64,
    pub norm: GeodeticNormalizations,
    pub variables: Vec<VariableId>,
}

impl UducCarrierPhaseFactor {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        obs: CorrectedObservation,
        epoch: u32,
        raw_cp_cycles: f64,
        var_pose: VariableId,
        var_clock: VariableId,
        var_zwd: Option<VariableId>,
        var_iono: Option<VariableId>,
        var_amb: VariableId,
        wavelength_m: f64,
        freq_ratio_sq: f64,
        windup_m: f64,
        norm: GeodeticNormalizations,
    ) -> Self {
        let mut vars = vec![var_pose, var_clock, var_amb];
        if let Some(z) = var_zwd {
            vars.push(z);
        }
        if let Some(i) = var_iono {
            vars.push(i);
        }
        Self {
            obs,
            epoch,
            raw_cp_cycles,
            var_pose,
            var_clock,
            var_zwd,
            var_iono,
            var_amb,
            wavelength_m,
            freq_ratio_sq,
            windup_m,
            norm,
            variables: vars,
        }
    }
}

impl Factor for UducCarrierPhaseFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let rx_pos = values
            .get(self.var_pose)
            .map(|p| Vector3::new(p[0], p[1], p[2]) + self.norm.solid_earth_tide_m + self.norm.ocean_tide_loading_m)
            .unwrap_or_else(Vector3::zeros);

        let c_dt = values.get(self.var_clock).map(|c| c[0]).unwrap_or(0.0);
        let zwd = self.var_zwd.and_then(|z| values.get(z)).map(|v| v[0]).unwrap_or(0.0);
        let iono = self.var_iono.and_then(|i| values.get(i)).map(|v| v[0]).unwrap_or(0.0);
        let amb = values.get(self.var_amb).map(|a| a[0]).unwrap_or(0.0);

        let sat_apc = self.obs.sat_pos_ecef + self.norm.sat_pco_ecef_m;
        let los = sat_apc - rx_pos;
        let geometric_range = los.norm();

        let modeled_tropo = self.obs.tropo_dry_m + self.obs.tropo_map_wet * zwd;
        let modeled_iono = -self.freq_ratio_sq * iono; // Phase advance
        let modeled_cp = geometric_range + c_dt - self.obs.sat_clock_m
            + self.norm.sat_relativity_m + self.norm.shapiro_delay_m
            + modeled_tropo + modeled_iono + self.wavelength_m * amb;

        let observed_cp_m = self.raw_cp_cycles * self.wavelength_m - self.windup_m;
        DVector::from_element(1, observed_cp_m - modeled_cp)
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(1, total_dim);

        let rx_pos = values
            .get(self.var_pose)
            .map(|p| Vector3::new(p[0], p[1], p[2]) + self.norm.solid_earth_tide_m + self.norm.ocean_tide_loading_m)
            .unwrap_or_else(Vector3::zeros);

        let sat_apc = self.obs.sat_pos_ecef + self.norm.sat_pco_ecef_m;
        let los = sat_apc - rx_pos;
        let range = los.norm().max(1.0);
        let unit_los = los / range;

        if let Some((s_pose, _)) = values.index_of(self.var_pose) {
            j[(0, s_pose)] = unit_los.x;
            j[(0, s_pose + 1)] = unit_los.y;
            j[(0, s_pose + 2)] = unit_los.z;
        }
        if let Some((s_clock, _)) = values.index_of(self.var_clock) {
            j[(0, s_clock)] = -1.0;
        }
        if let Some(z_id) = self.var_zwd {
            if let Some((s_zwd, _)) = values.index_of(z_id) {
                j[(0, s_zwd)] = -self.obs.tropo_map_wet;
            }
        }
        if let Some(i_id) = self.var_iono {
            if let Some((s_iono, _)) = values.index_of(i_id) {
                j[(0, s_iono)] = self.freq_ratio_sq; // Phase advance
            }
        }
        if let Some((s_amb, _)) = values.index_of(self.var_amb) {
            j[(0, s_amb)] = -self.wavelength_m;
        }

        j
    }

    fn information(&self) -> DMatrix<f64> {
        let var = (self.obs.cp_variance_m2).max(1e-6);
        DMatrix::from_element(1, 1, 1.0 / var)
    }
}

/// Random walk constraint between consecutive slant ionosphere states:
///   r = I_{r,1}^s(t_k) - I_{r,1}^s(t_{k-1}) \sim \mathcal{N}(0, q_I \Delta t)
#[derive(Debug, Clone)]
pub struct SlantIonoRandomWalkFactor {
    pub var_prev: VariableId,
    pub var_curr: VariableId,
    pub variance: f64,
    pub variables: Vec<VariableId>,
}

impl SlantIonoRandomWalkFactor {
    pub fn new(var_prev: VariableId, var_curr: VariableId, dt_s: f64, q_iono_m2_per_s: f64) -> Self {
        let variance = (q_iono_m2_per_s * dt_s.max(0.1)).max(1e-5);
        Self {
            var_prev,
            var_curr,
            variance,
            variables: vec![var_prev, var_curr],
        }
    }
}

impl Factor for SlantIonoRandomWalkFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let prev = values.get(self.var_prev).map(|v| v[0]).unwrap_or(0.0);
        let curr = values.get(self.var_curr).map(|v| v[0]).unwrap_or(0.0);
        DVector::from_element(1, curr - prev)
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(1, total_dim);

        if let Some((s_prev, _)) = values.index_of(self.var_prev) {
            j[(0, s_prev)] = -1.0;
        }
        if let Some((s_curr, _)) = values.index_of(self.var_curr) {
            j[(0, s_curr)] = 1.0;
        }

        j
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::from_element(1, 1, 1.0 / self.variance)
    }
}

#[cfg(test)]
#[path = "uduc_tests.rs"]
mod uduc_tests;


