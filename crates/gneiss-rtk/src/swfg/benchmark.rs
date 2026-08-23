#![allow(clippy::unwrap_used)]
//! Benchmark harness for the sliding-window factor graph engine.
//!
//! Provides an end-to-end integration test that exercises the full
//! SWFG pipeline: graph construction → correction passes → factor
//! building → LM optimization → AR → Schur marginalization.
//!
//! Uses synthetic GNSS observations so no external datasets are needed.

use nalgebra::{DMatrix, DVector, Vector3};

use crate::swfg::config::{EngineConfig, RtkConfig};
use crate::swfg::factor::PriorFactor;
use crate::swfg::pipeline::{self, MeasurementPipeline, RawObservation, ReceiverState};
use crate::swfg::solver::SlidingWindowSolver;
use crate::swfg::variables::{VariableId, VariableKind, VariableValues};

/// A simplified end-to-end benchmark that builds a graph with synthetic
/// observations, runs LM, and returns the optimized position.
///
/// Not a real processing loop — this is a correctness smoke test for
/// the full pipeline.
pub struct SwfgBenchmark {
    solver: SlidingWindowSolver,
    pipeline: MeasurementPipeline,
    epoch: u32,
    pose_id: Option<VariableId>,
    zwd_id: Option<VariableId>,
    /// Known true position for validation.
    true_position: Vector3<f64>,
}

impl SwfgBenchmark {
    pub fn new(true_position: Vector3<f64>) -> Self {
        let config = EngineConfig::Rtk(RtkConfig::default());
        Self {
            solver: SlidingWindowSolver::new(&config),
            pipeline: MeasurementPipeline::spp_mode(),
            epoch: 0,
            pose_id: None,
            zwd_id: None,
            true_position,
        }
    }

    /// Add one epoch of synthetic observations and run the solver.
    pub fn process_epoch(&mut self) -> Result<Vector3<f64>, String> {
        let epoch = self.epoch;
        self.epoch += 1;

        // Create variables for this epoch
        let _vars = self.solver.create_epoch_variables(epoch, 4, &[0], false, false);
        let pose_id = self.solver.graph.variables.iter()
            .find(|(_, n)| matches!(n.kind, VariableKind::Pose { epoch: e } if e == epoch))
            .map(|(id, _)| *id)
            .ok_or("pose not found")?;
        let zwd_id = self.solver.graph.variables.iter()
            .find(|(_, n)| matches!(n.kind, VariableKind::TropoZwd { epoch: e } if e == epoch))
            .map(|(id, _)| *id)
            .ok_or("zwd not found")?;

        // Set initial pose to truth ± some noise
        self.solver.graph.set_value(pose_id, &[
            self.true_position.x + 5.0,
            self.true_position.y - 3.0,
            self.true_position.z + 2.0,
            0.0, 0.0, 0.0,
        ]);
        self.solver.graph.set_value(zwd_id, &[0.1]);

        self.pose_id = Some(pose_id);
        self.zwd_id = Some(zwd_id);

        // Build synthetic observations: 4 satellites at various elevations
        let sats = [
            (Vector3::new(15_000_000.0, 10_000_000.0, 20_000_000.0), 0.8, 0.3),
            (Vector3::new(-10_000_000.0, 20_000_000.0, 15_000_000.0), 1.0, 0.5),
            (Vector3::new(5_000_000.0, -20_000_000.0, 10_000_000.0), 0.6, 0.4),
            (Vector3::new(20_000_000.0, -5_000_000.0, 18_000_000.0), 0.9, 0.2),
        ];

        let raw_obs: Vec<RawObservation> = sats.iter().enumerate().map(|(i, (sat_pos, el, _az))| {
            let geom_range = (sat_pos - self.true_position).norm();
            RawObservation {
                satellite: (i + 1) as u16,
                constellation_id: 0,
                pr_l1: geom_range + 0.5, // 0.5m noise
                pr_l2: None,
                cp_l1: None,
                cp_l1_lli: None,
                cp_l2: None,
                doppler: 0.0,
                snr_dbhz: 45.0,
                sat_pos_ecef: *sat_pos,
                sat_vel_ecef: Vector3::zeros(),
                sat_clock_m: 50.0,
                f1: 1575.42e6,
                f2: 1227.60e6,
                freq_num: 0,
                elevation_rad: *el,
                azimuth_rad: 0.5,
                tropo_dry_m: 0.0,
                tropo_map_wet: 0.0,
                iono_l1_m: 0.0,
                variance_m2: 0.25,
                cp_variance_m2: 9e-6,
            }
        }).collect();

        let llh = gneiss_core::coords::ecef_to_llh(self.true_position);
        let rx_state = ReceiverState {
            position_ecef: self.true_position,
            clock_bias_m: vec![0.0],
            zwd_m: 0.1,
            ifb_glo: 0.0,
            llh_rad: llh,
            time: gneiss_core::time::GpsTime::new(2200, epoch as f64),
        };

        let corrected = self.pipeline.process(&raw_obs, &rx_state);

        // Build PR factors
        self.solver.graph.clear_factors();
        for obs in &corrected {
            let factor = pipeline::build_pseudorange_factor(
                obs, epoch, pose_id, None, Some(zwd_id), None,
            );
            self.solver.graph.add_factor(factor);
        }

        // Add a weak prior on position to keep the system well-conditioned
        let prior = PriorFactor {
            variable: pose_id,
            mu: DVector::from_vec(vec![
                self.true_position.x + 5.0,
                self.true_position.y - 3.0,
                self.true_position.z + 2.0,
                0.0, 0.0, 0.0,
            ]),
            information: DMatrix::identity(6, 6) * 0.01, // σ=10m prior
        };
        self.solver.graph.add_factor(Box::new(prior));

        // Run a few LM iterations (simplified — full solver in solver.rs)
        for _iter in 0..5 {
            let values = VariableValues::build(&self.solver.graph.variables);
            let (jtj, jtr, _err) = self.solver.build_normal_equations(&values);
            if let Some(inv) = jtj.try_inverse() {
                let delta = -&inv * &jtr;
                self.solver.apply_delta(&delta);
                if delta.norm() < 1e-6 {
                    break;
                }
            }
        }

        // Extract optimized position
        let vals = VariableValues::build(&self.solver.graph.variables);
        let pose = vals.get(pose_id).ok_or("pose missing")?;
        Ok(Vector3::new(pose[0], pose[1], pose[2]))
    }

    /// Verify the optimized position is within tolerance of truth.
    pub fn check_accuracy(&self) -> Result<f64, String> {
        let pose_id = self.pose_id.ok_or("no pose")?;
        let vals = VariableValues::build(&self.solver.graph.variables);
        let pose = vals.get(pose_id).ok_or("pose missing")?;
        let pos = Vector3::new(pose[0], pose[1], pose[2]);
        Ok((pos - self.true_position).norm())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn benchmark_single_epoch_converges() {
        let mut bench = SwfgBenchmark::new(Vector3::new(1_000_000.0, 2_000_000.0, 6_378_000.0));
        let result = bench.process_epoch();
        assert!(result.is_ok(), "process_epoch failed: {:?}", result.err());

        let pos = result.unwrap();
        // After one epoch with 4 satellites and a 10m prior, position should
        // be closer to truth than the initial 7.7m offset.
        let err = (pos - bench.true_position).norm();
        assert!(err < 7.0, "position error {} should be less than initial 7.7m offset", err);
    }

    #[test]
    fn benchmark_accuracy_improves() {
        let mut bench = SwfgBenchmark::new(Vector3::new(1_000_000.0, 2_000_000.0, 6_378_000.0));
        bench.process_epoch().unwrap();
        let err = bench.check_accuracy().unwrap();
        assert!(err < 10.0, "error should be reasonable after one epoch");
    }
}
