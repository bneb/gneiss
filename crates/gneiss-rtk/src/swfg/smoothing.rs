//! Bi-directional batch factor graph smoothing for GNSS/INS trajectory post-processing.
//!
//! Spans all $T$ epochs of a trajectory in a unified Levenberg-Marquardt factor graph
//! to achieve zero initial PPP convergence delay and instant ambiguity re-lock.

use nalgebra::{DVector, Vector3};

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_core::time::GpsTime;

use crate::swfg::ar_integration::{attempt_ar_fix, collect_ambiguity_variables, extract_ambiguity_state, inject_fixed_priors, ArResult};
use crate::swfg::factor::PriorFactor;
use crate::swfg::graph::EstimationGraph;
use crate::swfg::pipeline::{build_carrier_phase_factor, build_pseudorange_factor, MeasurementPipeline};
use crate::swfg::solver::SlidingWindowSolver;
use crate::swfg::variables::{VariableId, VariableKind, VariableValues};


/// Configuration for batch factor graph trajectory smoothing.
#[derive(Debug, Clone)]
pub struct BatchSmoothingConfig {
    pub max_iterations: usize,
    pub convergence_tol: f64,
    pub enable_ar: bool,
    pub ar_ratio_threshold: f64,
}

impl Default for BatchSmoothingConfig {
    fn default() -> Self {
        Self {
            max_iterations: 20,
            convergence_tol: 1e-4,
            enable_ar: true,
            ar_ratio_threshold: 2.0,
        }
    }
}

/// A batch factor graph spanning an entire trajectory.
pub struct BatchFactorGraph {
    pub config: BatchSmoothingConfig,
    pub graph: EstimationGraph,
    pub pipeline: MeasurementPipeline,
    pub ephemerides: Vec<Ephemeris>,
    pub epochs: Vec<GpsTime>,
    pub pose_ids: Vec<VariableId>,
}

impl BatchFactorGraph {
    pub fn new(config: BatchSmoothingConfig, ephemerides: Vec<Ephemeris>) -> Self {
        Self {
            config,
            graph: EstimationGraph::new(),
            pipeline: MeasurementPipeline::new(),
            ephemerides,
            epochs: Vec::new(),
            pose_ids: Vec::new(),
        }
    }

    /// Add an epoch of observations to the batch factor graph.
    pub fn add_epoch(&mut self, epoch_idx: u32, obs: &EpochObs, base: Option<(&EpochObs, Vector3<f64>)>, initial_pos: Option<Vector3<f64>>) -> VariableId {
        let init = initial_pos.unwrap_or_else(|| Vector3::new(-3963427.0, 3350882.0, 3694866.0));
        let pose_id = self.graph.add_variable(VariableKind::Pose { epoch: epoch_idx });
        self.graph.set_value(pose_id, &[init.x, init.y, init.z, 0.0, 0.0, 0.0]);

        if epoch_idx == 0 {
            // Anchor first epoch with a weak prior if no base position is supplied
            let prior = PriorFactor::new(
                pose_id,
                DVector::from_row_slice(&[init.x, init.y, init.z, 0.0, 0.0, 0.0]),
                100.0,
            );
            self.graph.add_factor(Box::new(prior));
        }

        let clock_id = self.graph.add_variable(VariableKind::ClockBias { epoch: epoch_idx, constellation_id: 0 });

        let mut zwd_id = None;
        if base.is_none() {
            let id = self.graph.add_variable(VariableKind::TropoZwd { epoch: epoch_idx });
            self.graph.set_value(id, &[0.1]);
            zwd_id = Some(id);
        }

        let rx_state = crate::swfg::pipeline::ReceiverState {
            position_ecef: init,
            llh_rad: gneiss_core::coords::ecef_to_llh(init),
            clock_bias_m: vec![0.0],
            zwd_m: 0.1,
            ifb_glo: 0.0,
            time: obs.time,
        };

        let mut raw_obs = crate::swfg::engine::epoch::extract_raw_observations(obs, &self.ephemerides, Some(init))
            .unwrap_or_default();

        if let Some((base_obs, base_pos)) = base {
            let base_raw = crate::swfg::engine::epoch::extract_raw_observations(base_obs, &self.ephemerides, Some(base_pos))
                .unwrap_or_default();
            raw_obs = raw_obs.into_iter().filter_map(|rov| {
                base_raw.iter().find(|b| b.satellite == rov.satellite).map(|bas| {
                    let range_b = (bas.sat_pos_ecef - base_pos).norm();
                    let l1_wavelength = gneiss_core::constants::SPEED_OF_LIGHT_M_S / bas.f1;
                    let l2_wavelength = gneiss_core::constants::SPEED_OF_LIGHT_M_S / bas.f2;
                    crate::swfg::pipeline::RawObservation {
                        pr_l1: (rov.pr_l1 - bas.pr_l1) + range_b,
                        pr_l2: rov.pr_l2.zip(bas.pr_l2).map(|(r, b)| (r - b) + range_b),
                        cp_l1: rov.cp_l1.zip(bas.cp_l1).map(|(r, b)| (r - b) + range_b / l1_wavelength),
                        cp_l2: rov.cp_l2.zip(bas.cp_l2).map(|(r, b)| (r - b) + range_b / l2_wavelength),
                        sat_clock_m: 0.0,
                        tropo_dry_m: 0.0,
                        iono_l1_m: 0.0,
                        ..rov
                    }
                })
            }).collect();
        }

        let corrected = self.pipeline.process(&raw_obs, &rx_state);

        let mut clock_offsets = Vec::new();
        for c_obs in &corrected {
            let geometric_range = (c_obs.sat_pos_ecef - init).norm();
            let predicted_no_clk = geometric_range 
                - c_obs.sat_clock_m 
                + c_obs.tropo_dry_m 
                + c_obs.iono_l1_m;
            clock_offsets.push(c_obs.pr_l1 - predicted_no_clk);
        }
        if !clock_offsets.is_empty() {
            clock_offsets.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let median_clk = clock_offsets[clock_offsets.len() / 2];
            
            // Seed the clock variable to prevent robust estimators from completely downweighting everything
            if let Some(node) = self.graph.variables.get_mut(&clock_id) {
                node.set_value(&[median_clk, 0.0, 0.0]);
            }
        }

        for c_obs in &corrected {
            let pr_factor = build_pseudorange_factor(
                c_obs, epoch_idx, pose_id, Some(clock_id), zwd_id, None,
            );
            self.graph.add_factor(pr_factor);

            if c_obs.cp_l1.is_some() {
                let amb_id = self.ensure_ambiguity(c_obs.satellite, 1);
                let cp_factor = build_carrier_phase_factor(
                    c_obs, epoch_idx, pose_id, Some(clock_id), zwd_id, amb_id, None, 0.0,
                );
                self.graph.add_factor(cp_factor);
            }
        }

        self.epochs.push(obs.time);
        self.pose_ids.push(pose_id);
        pose_id
    }

    fn ensure_ambiguity(&mut self, satellite: u16, frequency: u8) -> VariableId {
        let existing = self.graph.variables.iter().find_map(|(id, n)| match n.kind {
            VariableKind::Ambiguity { satellite: s, frequency: f } if s == satellite && f == frequency => Some(*id),
            _ => None,
        });
        if let Some(id) = existing {
            id
        } else {
            self.graph.add_variable(VariableKind::Ambiguity { satellite, frequency })
        }
    }

    /// Add a smooth kinematic motion prior between consecutive pose variables.
    pub fn add_smoothness_factors(&mut self, max_velocity_m_s: f64, dt_s: f64) {
        let var_m = (max_velocity_m_s * dt_s).powi(2);
        for i in 0..(self.pose_ids.len() - 1) {
            let p1 = self.pose_ids[i];
            let p2 = self.pose_ids[i + 1];
            let factor = crate::swfg::factor::RelativePoseFactor::new(p1, p2, var_m);
            self.graph.add_factor(Box::new(factor));
        }
    }

    /// Run full-batch bi-directional optimization over all trajectory epochs.
    pub fn solve(&mut self) -> Result<Vec<(GpsTime, Vector3<f64>)>, String> {
        let engine_cfg = crate::swfg::config::EngineConfig::Spp(crate::swfg::config::SppConfig::default());
        let mut solver = SlidingWindowSolver::new(&engine_cfg);
        solver.graph = std::mem::take(&mut self.graph);

        solver.solve().map_err(|e| format!("Batch solve failed: {:?}", e))?;

        if self.config.enable_ar {
            let amb_ids = collect_ambiguity_variables(&solver.graph);
            if amb_ids.len() >= 4 {
                let values = VariableValues::build(&solver.graph.variables);
                let (jtj, _, _) = solver.build_normal_equations(&values);
                let (float_amb, amb_cov) = extract_ambiguity_state(&solver.graph, &jtj, &amb_ids);
                if let ArResult::Fixed { integers, .. } = attempt_ar_fix(&float_amb, &amb_cov, self.config.ar_ratio_threshold) {
                    for (i, &amb_id) in amb_ids.iter().enumerate() {
                        solver.graph.set_value(amb_id, &[integers[i]]);
                    }
                    inject_fixed_priors(&mut solver.graph, &amb_ids, &integers);
                    let _ = solver.solve();
                }
            }
        }

        let values = VariableValues::build(&solver.graph.variables);
        let mut trajectory = Vec::new();

        for (i, &pose_id) in self.pose_ids.iter().enumerate() {
            if let Some(pos) = values.get(pose_id) {
                trajectory.push((self.epochs[i], Vector3::new(pos[0], pos[1], pos[2])));
            }
        }

        self.graph = std::mem::take(&mut solver.graph);
        Ok(trajectory)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::sat::{Constellation, SatelliteId};

    fn make_test_obs(time: GpsTime) -> EpochObs {
        use gneiss_core::obs::{ObsCode, ObsType, Observation, SatObs, SignalCode};
        let mut obs = EpochObs {
            time,
            satellites: Vec::new(),
        };
        for prn in 1..=5 {
            let sat = SatelliteId { constellation: Constellation::Gps, prn };
            obs.satellites.push(SatObs {
                sat,
                observations: vec![Observation {
                    code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 1, attribute: 'C' } },
                    value: 20_000_000.0 + (prn as f64) * 100.0,
                    lock_time: None,
                    lli: None,
                }],
            });
        }
        obs
    }

    #[test]
    fn test_batch_smoothing_scaffolding() {
        let config = BatchSmoothingConfig::default();
        let mut batch = BatchFactorGraph::new(config, Vec::new());
        let t0 = GpsTime::new(2000, 100.0);
        let obs0 = make_test_obs(t0);
        let p0 = batch.add_epoch(0, &obs0, None, None);
        assert_eq!(batch.pose_ids.len(), 1);
        assert_eq!(batch.pose_ids[0], p0);
    }
}
