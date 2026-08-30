//! Factor graph observation factors.
//!
//! Undifferenced: PseudorangeFactor, CarrierPhaseFactor (SPP/PPP)
//! Dynamics: DopplerVelocityFactor (between-epoch constraint)

pub mod dynamics;
#[cfg(test)]
mod tests;

pub use dynamics::{DopplerVelocityFactor, NhcFactor};

use nalgebra::{DMatrix, DVector, Vector3};

use crate::swfg::factor::Factor;
use crate::swfg::pipeline::passes::CorrectedObservation;
use crate::swfg::variables::{VariableId, VariableValues};

/// Build a pseudorange factor from a corrected observation.
pub fn build_pseudorange_factor(
    obs: &CorrectedObservation,
    epoch: u32,
    var_pose: VariableId,
    var_clock: Option<VariableId>,
    var_zwd: Option<VariableId>,
    var_ifb: Option<VariableId>,
) -> Box<dyn Factor> {
    Box::new(PseudorangeFactor {
        obs: obs.clone(),
        epoch,
        var_pose,
        var_clock,
        var_zwd,
        var_ifb,
        variables: build_var_list(var_pose, var_clock, var_zwd, var_ifb),
    })
}

/// Build a carrier-phase factor from a corrected observation.
#[allow(clippy::too_many_arguments)]
pub fn build_carrier_phase_factor(
    obs: &CorrectedObservation,
    epoch: u32,
    var_pose: VariableId,
    var_clock: Option<VariableId>,
    var_zwd: Option<VariableId>,
    var_amb: VariableId,
    var_ifb: Option<VariableId>,
    windup_m: f64,
) -> Box<dyn Factor> {
    Box::new(CarrierPhaseFactor {
        obs: obs.clone(),
        epoch,
        var_pose,
        var_clock,
        var_zwd,
        var_amb,
        var_ifb,
        windup_m,
        variables: build_var_list_cp(var_pose, var_clock, var_zwd, var_amb, var_ifb),
    })
}

fn build_var_list(
    pose: VariableId,
    clock: Option<VariableId>,
    zwd: Option<VariableId>,
    ifb: Option<VariableId>,
) -> Vec<VariableId> {
    let mut v = vec![pose];
    if let Some(z) = zwd {
        v.push(z);
    }
    if let Some(c) = clock {
        v.push(c);
    }
    if let Some(i) = ifb {
        v.push(i);
    }
    v
}

fn build_var_list_cp(
    pose: VariableId,
    clock: Option<VariableId>,
    zwd: Option<VariableId>,
    amb: VariableId,
    ifb: Option<VariableId>,
) -> Vec<VariableId> {
    let mut v = vec![pose, amb];
    if let Some(z) = zwd {
        v.push(z);
    }
    if let Some(c) = clock {
        v.push(c);
    }
    if let Some(i) = ifb {
        v.push(i);
    }
    v
}

/// Carrier-phase factor: r = (cp_obs - windup) - (range + clock + tropo - iono + lambda*N + ifb).
pub struct CarrierPhaseFactor {
    pub obs: CorrectedObservation,
    pub epoch: u32,
    pub var_pose: VariableId,
    pub var_clock: Option<VariableId>,
    pub var_zwd: Option<VariableId>,
    pub var_amb: VariableId,
    pub var_ifb: Option<VariableId>,
    pub windup_m: f64,
    pub variables: Vec<VariableId>,
}

impl std::fmt::Debug for CarrierPhaseFactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CarrierPhaseFactor")
            .field("sat", &self.obs.satellite)
            .field("epoch", &self.epoch)
            .finish()
    }
}

impl Factor for CarrierPhaseFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let pose = match values.get(self.var_pose) {
            Some(p) => p,
            None => return DVector::zeros(1),
        };
        let rx_pos = Vector3::new(pose[0], pose[1], pose[2]);
        let geometric_range = (self.obs.sat_pos_ecef - rx_pos).norm();
        let zwd = self
            .var_zwd
            .and_then(|z| values.get(z))
            .map_or(0.0, |v| v[0]);
        let amb_cycles = match values.get(self.var_amb) {
            Some(v) => v[0],
            None => 0.0,
        };
        let lambda = gneiss_core::constants::SPEED_OF_LIGHT_M_S / self.obs.f1.max(1.0);
        let rx_clk = self
            .var_clock
            .and_then(|c| values.get(c))
            .map_or(0.0, |v| v[0]);
        let ifb_term = self
            .var_ifb
            .and_then(|i| values.get(i))
            .map_or(0.0, |v| v[0] * self.obs.freq_num as f64);

        let predicted = geometric_range
            - self.obs.sat_clock_m
            + rx_clk
            + self.obs.tropo_dry_m
            + zwd * self.obs.tropo_map_wet
            - self.obs.iono_l1_m
            + lambda * amb_cycles
            + ifb_term;

        let cp_raw_m = self.obs.cp_l1.unwrap_or(0.0) * lambda;
        // Phase windup is ALWAYS subtracted from carrier phase per AGENTS.md
        let cp_corr_m = cp_raw_m - self.windup_m;

        DVector::from_element(1, cp_corr_m - predicted)
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(1, total_dim);
        let start_pose = match values.index_of(self.var_pose) {
            Some((s, _)) => s,
            None => return j,
        };
        let start_amb = match values.index_of(self.var_amb) {
            Some((s, _)) => s,
            None => return j,
        };

        let pose = match values.get(self.var_pose) {
            Some(p) => p,
            None => return j,
        };
        let rx_pos = Vector3::new(pose[0], pose[1], pose[2]);
        let los = (self.obs.sat_pos_ecef - rx_pos).normalize();

        for k in 0..3 {
            j[(0, start_pose + k)] = los[k];
        }
        if let Some(zwd_id) = self.var_zwd {
            if let Some((start_zwd, _)) = values.index_of(zwd_id) {
                j[(0, start_zwd)] = -self.obs.tropo_map_wet;
            }
        }
        let lambda = gneiss_core::constants::SPEED_OF_LIGHT_M_S / self.obs.f1.max(1.0);
        j[(0, start_amb)] = -lambda;

        if let Some(clk_id) = self.var_clock {
            if let Some((start_clk, _)) = values.index_of(clk_id) {
                j[(0, start_clk)] = -1.0;
            }
        }
        if let Some(ifb_id) = self.var_ifb {
            if let Some((start_ifb, _)) = values.index_of(ifb_id) {
                j[(0, start_ifb)] = -self.obs.freq_num as f64;
            }
        }

        j
    }

    fn information(&self) -> DMatrix<f64> {
        DMatrix::from_element(1, 1, 1.0 / self.obs.cp_variance_m2.max(1e-9))
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(0.05)
    }

    fn use_cauchy(&self) -> bool {
        false
    }
}

/// Pseudorange factor: r = pr_obs - (range + clock + tropo + iono + ifb).
pub struct PseudorangeFactor {
    pub obs: CorrectedObservation,
    pub epoch: u32,
    pub var_pose: VariableId,
    pub var_clock: Option<VariableId>,
    pub var_zwd: Option<VariableId>,
    pub var_ifb: Option<VariableId>,
    pub variables: Vec<VariableId>,
}

impl std::fmt::Debug for PseudorangeFactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PseudorangeFactor")
            .field("sat", &self.obs.satellite)
            .field("epoch", &self.epoch)
            .finish()
    }
}

impl Factor for PseudorangeFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        let pose = match values.get(self.var_pose) {
            Some(p) => p,
            None => return DVector::zeros(1),
        };
        let rx_pos = Vector3::new(pose[0], pose[1], pose[2]);
        let geometric_range = (self.obs.sat_pos_ecef - rx_pos).norm();
        let zwd = self
            .var_zwd
            .and_then(|z| values.get(z))
            .map_or(0.0, |v| v[0]);
        let rx_clk = self
            .var_clock
            .and_then(|c| values.get(c))
            .map_or(0.0, |v| v[0]);
        let ifb_term = self
            .var_ifb
            .and_then(|i| values.get(i))
            .map_or(0.0, |v| v[0] * self.obs.freq_num as f64);

        if let Some(pr2) = self.obs.pr_l2 {
            let f1_sq = self.obs.f1 * self.obs.f1;
            let f2_sq = self.obs.f2.max(1.0) * self.obs.f2.max(1.0);
            let gamma = f1_sq / f2_sq;
            if (gamma - 1.0).abs() > 0.1 {
                let pr_if = (gamma * self.obs.pr_l1 - pr2) / (gamma - 1.0);
                let predicted = geometric_range
                    - self.obs.sat_clock_m
                    + rx_clk
                    + self.obs.tropo_dry_m
                    + zwd * self.obs.tropo_map_wet
                    + ifb_term;
                return DVector::from_element(1, pr_if - predicted);
            }
        }

        let predicted = geometric_range
            - self.obs.sat_clock_m
            + rx_clk
            + self.obs.tropo_dry_m
            + zwd * self.obs.tropo_map_wet
            + self.obs.iono_l1_m
            + ifb_term;

        DVector::from_element(1, self.obs.pr_l1 - predicted)
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(1, total_dim);
        let start_pose = match values.index_of(self.var_pose) {
            Some((s, _)) => s,
            None => return j,
        };

        let pose = match values.get(self.var_pose) {
            Some(p) => p,
            None => return j,
        };
        let rx_pos = Vector3::new(pose[0], pose[1], pose[2]);
        let los = (self.obs.sat_pos_ecef - rx_pos).normalize();

        for k in 0..3 {
            j[(0, start_pose + k)] = los[k];
        }
        if let Some(zwd_id) = self.var_zwd {
            if let Some((start_zwd, _)) = values.index_of(zwd_id) {
                j[(0, start_zwd)] = -self.obs.tropo_map_wet;
            }
        }

        if let Some(clk_id) = self.var_clock {
            if let Some((start_clk, _)) = values.index_of(clk_id) {
                j[(0, start_clk)] = -1.0;
            }
        }
        if let Some(ifb_id) = self.var_ifb {
            if let Some((start_ifb, _)) = values.index_of(ifb_id) {
                j[(0, start_ifb)] = -self.obs.freq_num as f64;
            }
        }

        j
    }

    fn information(&self) -> DMatrix<f64> {
        let var = if self.obs.pr_l2.is_some() {
            let f1_sq = self.obs.f1 * self.obs.f1;
            let f2_sq = self.obs.f2.max(1.0) * self.obs.f2.max(1.0);
            let gamma = f1_sq / f2_sq;
            if (gamma - 1.0).abs() > 0.1 {
                let factor = (gamma * gamma + 1.0) / ((gamma - 1.0) * (gamma - 1.0));
                self.obs.variance_m2 * factor
            } else {
                self.obs.variance_m2
            }
        } else {
            self.obs.variance_m2
        };
        DMatrix::from_element(1, 1, 1.0 / var.max(1e-4))
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(6.0)
    }

    fn use_cauchy(&self) -> bool {
        false
    }
}
