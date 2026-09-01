//! SWFG processing engine — bridges the factor graph to real GNSS data.

pub mod accumulator;
pub mod ar_handler;
pub mod builder;
pub mod epoch;
pub mod uduc_builder;

use std::collections::{BTreeSet, HashMap, HashSet};
use nalgebra::Vector3;

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_core::time::GpsTime;

use crate::swfg::config::EngineConfig;
use crate::swfg::engine::builder::RtkFactorContext;
use crate::swfg::imu_preintegration::ImuPreintegration;
use crate::swfg::pipeline::{MeasurementPipeline, RawObservation, ReceiverState};
use crate::swfg::solver::SlidingWindowSolver;
use crate::swfg::variables::{VariableId, VariableKind, VariableValues};

/// A solution produced by the SWFG engine for one epoch.
#[derive(Debug, Clone)]
pub struct SwfgSolution {
    pub time: GpsTime,
    pub position_ecef: Vector3<f64>,
    pub clock_bias_m: f64,
    pub n_satellites: usize,
    pub solver_iterations: usize,
    pub error: Option<f64>,
}

/// The SWFG processing engine.
pub struct SwfgEngine {
    solver: SlidingWindowSolver,
    pipeline: MeasurementPipeline,
    ephemerides: Vec<Ephemeris>,
    epoch: u32,
    pub(crate) current_pose: Option<VariableId>,
    prev_position: Option<Vector3<f64>>,
    prev_prev_position: Option<Vector3<f64>>,
    prev_time: Option<GpsTime>,
    current_attitude: Option<nalgebra::UnitQuaternion<f64>>,
    initial_position: Option<Vector3<f64>>,
    slip_counts: HashMap<u16, u32>,
    ref_sat_per_constellation: HashMap<u8, u16>,
    elevation_mask_rad: f64,
    mw_accumulator: accumulator::DdPseudorangeAccumulator,
    is_ppp: bool,
    is_kinematic: bool,
    prev_zwd: Option<VariableId>,
    windup_trackers: HashMap<u16, gneiss_geodesy::windup::PhaseWindupTracker>,
    pub precise_orbits: Option<std::sync::Arc<gneiss_parsers::precise_orbit::PreciseOrbit>>,
    pub precise_clocks: Option<std::sync::Arc<gneiss_parsers::rinex_clk::RinexClock>>,
    pub sinex_bias: Option<std::sync::Arc<gneiss_parsers::sinex_bia::SinexBias>>,
}

impl SwfgEngine {
    pub fn new(config: &EngineConfig, ephemerides: Vec<Ephemeris>) -> Self {
        let (pipeline, initial_position, elevation_mask_deg, is_ppp, is_kinematic) = match config {
            EngineConfig::Ppp(c) => (MeasurementPipeline::ppp_mode(), c.initial_position.map(|p| Vector3::new(p[0], p[1], p[2])), c.elevation_mask_deg, true, c.is_kinematic),
            EngineConfig::Rtk(c) => (MeasurementPipeline::rtk_mode(), c.initial_position.map(|p| Vector3::new(p[0], p[1], p[2])), c.elevation_mask_deg, false, false),
            EngineConfig::Spp(c) => (MeasurementPipeline::spp_mode(), c.initial_position.map(|p| Vector3::new(p[0], p[1], p[2])), c.elevation_mask_deg, false, false),
            EngineConfig::RtkIns(c) => (MeasurementPipeline::rtk_mode(), c.rtk.initial_position.map(|p| Vector3::new(p[0], p[1], p[2])), c.rtk.elevation_mask_deg, false, true),
            EngineConfig::PppIns(c) => (MeasurementPipeline::ppp_mode(), None, c.ppp.elevation_mask_deg, true, c.ppp.is_kinematic),
        };
        Self {
            solver: SlidingWindowSolver::new(config),
            pipeline, ephemerides, epoch: 0, current_pose: None,
            prev_position: None, prev_prev_position: None, prev_time: None,
            current_attitude: None, initial_position, slip_counts: HashMap::new(),
            ref_sat_per_constellation: HashMap::new(),
            elevation_mask_rad: elevation_mask_deg.to_radians(),
            mw_accumulator: accumulator::DdPseudorangeAccumulator::new(10),
            is_ppp, is_kinematic, prev_zwd: None, windup_trackers: HashMap::new(),
            precise_orbits: None, precise_clocks: None, sinex_bias: None,
        }
    }

    pub fn set_precise_products(&mut self, orbits: std::sync::Arc<gneiss_parsers::precise_orbit::PreciseOrbit>, clocks: Option<std::sync::Arc<gneiss_parsers::rinex_clk::RinexClock>>) {
        self.precise_orbits = Some(orbits);
        self.precise_clocks = clocks;
    }

    pub fn set_sinex_bias(&mut self, bias: std::sync::Arc<gneiss_parsers::sinex_bia::SinexBias>) {
        self.sinex_bias = Some(bias);
    }

    fn get_initial_position(&self, rover: &EpochObs) -> Vector3<f64> {
        let default_eq = Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        if let Some(prev) = self.prev_position { return prev; }
        if let Some(init) = self.initial_position { return init; }
        if let Some(spp) = epoch::compute_spp_seeding(rover, &self.ephemerides) {
            if spp.norm() > 1e6 && (spp - default_eq).norm() > 1000.0 { return spp; }
        }
        default_eq
    }

    pub fn process_epoch(&mut self, rover: &EpochObs) -> Result<SwfgSolution, String> {
        self.process_impl(rover, None, None)
    }

    pub fn process_rtk_epoch(
        &mut self,
        rover: &EpochObs,
        base: &EpochObs,
        base_position: Vector3<f64>,
    ) -> Result<SwfgSolution, String> {
        self.process_impl(rover, Some((base, base_position)), None)
    }

    pub fn process_rtk_epoch_with_imu(
        &mut self,
        rover: &EpochObs,
        base: &EpochObs,
        base_position: Vector3<f64>,
        imu_preint: Option<ImuPreintegration>,
    ) -> Result<SwfgSolution, String> {
        self.process_impl(rover, Some((base, base_position)), imu_preint)
    }

    fn process_impl(
        &mut self,
        rover: &EpochObs,
        base: Option<(&EpochObs, Vector3<f64>)>,
        imu_preint: Option<ImuPreintegration>,
    ) -> Result<SwfgSolution, String> {
        let is_rtk = base.is_some();
        let has_imu = imu_preint.is_some();
        let init_pos = self.get_initial_position(rover);
        let epoch = self.epoch;
        self.epoch += 1;

        let pose_id = self.setup_epoch_variables(rover, epoch, has_imu, is_rtk)?;
        let prev_pose_id = self.current_pose;
        self.current_pose = Some(pose_id);

        self.setup_imu_and_rel_factors(epoch, pose_id, prev_pose_id, init_pos, &imu_preint);
        self.setup_priors(epoch, pose_id, prev_pose_id, init_pos, has_imu);
        let dt_sec = self.prev_time.map_or(1.0, |t| (rover.time.tow - t.tow).abs());
        let zwd_id = self.setup_zwd(epoch, is_rtk, dt_sec);

        let (raw_obs, base_raw) = self.filter_observations(rover, base, init_pos)?;
        if raw_obs.len() < 4 && !has_imu {
            return Err(format!("too few observations: {}", raw_obs.len()));
        }

        let corrected = self.apply_corrections(rover, &raw_obs, init_pos);
        self.seed_clocks(epoch, &corrected, init_pos, is_rtk);

        let cp_records = self.build_observation_factors(
            epoch, pose_id, prev_pose_id, init_pos, rover.time, base,
            &corrected, &base_raw, zwd_id,
        )?;

        self.prune_unobserved_variables();
        self.solver.solve().map_err(|e| format!("solve failed: {:?}", e))?;

        ar_handler::execute_ar_step(&mut self.solver, pose_id, init_pos, is_rtk, epoch);
        let sol = self.extract_and_finalize_solution(
            pose_id, init_pos, rover.time, corrected.len(), has_imu, is_rtk, &cp_records, epoch,
        )?;
        Ok(sol)
    }

    fn setup_epoch_variables(
        &mut self, rover: &EpochObs, epoch: u32, has_imu: bool, is_rtk: bool,
    ) -> Result<VariableId, String> {
        let n_sats = rover.satellites.len();
        let constellations: Vec<u8> = rover.satellites.iter()
            .map(|s| s.sat.constellation as u8)
            .collect::<BTreeSet<_>>().into_iter().collect();
        if self.is_ppp && !self.is_kinematic {
            let pose_id = self.solver.ensure_static_pose();
            self.solver.create_epoch_variables(epoch, n_sats, &constellations, has_imu, is_rtk);
            Ok(pose_id)
        } else {
            self.solver.create_epoch_variables(epoch, n_sats, &constellations, has_imu, is_rtk);
            self.solver.graph.variables.iter()
                .find(|(_, n)| matches!(n.kind, VariableKind::Pose { epoch: e } if e == epoch))
                .map(|(id, _)| *id)
                .ok_or_else(|| "pose not found after create_epoch_variables".to_string())
        }
    }

    fn setup_imu_and_rel_factors(
        &mut self,
        epoch: u32,
        pose_id: VariableId,
        prev_pose_id: Option<VariableId>,
        init_pos: Vector3<f64>,
        imu_preint: &Option<ImuPreintegration>,
    ) {
        if imu_preint.is_some() && self.current_attitude.is_none() {
            let llh = gneiss_core::coords::ecef_to_llh(init_pos);
            let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(llh).transpose();
            let rot = nalgebra::Rotation3::from_matrix_unchecked(ned_to_ecef);
            self.current_attitude = Some(nalgebra::UnitQuaternion::from_rotation_matrix(&rot));
        }
        if let (Some(preint), Some(prev_p)) = (imu_preint, prev_pose_id) {
            self.add_imu_preintegration_factors(epoch, pose_id, prev_p, init_pos, preint);
            self.current_attitude = self.current_attitude.map(|q| q * preint.dq);
        } else if let Some(prev_p) = prev_pose_id {
            if prev_p != pose_id && self.solver.graph.variables.contains_key(&prev_p) {
                let mut rel_info = nalgebra::DMatrix::zeros(6, 6);
                let q_pos = if self.is_ppp && !self.is_kinematic { 1.0 / 1e-4 } else { 1.0 / 25.0 };
                for i in 0..3 { rel_info[(i, i)] = q_pos; }
                for i in 3..6 { rel_info[(i, i)] = 1.0; }
                let rel = crate::swfg::factor::RelativePoseFactor { vars: [prev_p, pose_id], information: rel_info };
                self.solver.graph.add_factor(Box::new(rel));
            }
        }
    }

    fn add_imu_preintegration_factors(
        &mut self, epoch: u32, pose_id: VariableId, prev_p: VariableId,
        init_pos: Vector3<f64>, preint: &ImuPreintegration,
    ) {
        let vel_i = self.solver.graph.variables.iter()
            .find(|(_, n)| matches!(n.kind, VariableKind::Velocity { epoch: e } if e == epoch - 1))
            .map(|(id, _)| *id);
        let vel_j = self.solver.graph.variables.iter()
            .find(|(_, n)| matches!(n.kind, VariableKind::Velocity { epoch: e } if e == epoch))
            .map(|(id, _)| *id);
        let bias = self.solver.graph.variables.iter()
            .find(|(_, n)| matches!(n.kind, VariableKind::ImuBias))
            .map(|(id, _)| *id);

        if let (Some(vi), Some(vj), Some(b)) = (vel_i, vel_j, bias) {
            let att_i = self.current_attitude.unwrap_or_else(nalgebra::UnitQuaternion::identity);
            let att_j = att_i * preint.dq;
            let grav = -9.80665 * init_pos.normalize();
            let fac = crate::swfg::imu_preintegration::ImuPreintegrationFactor::new(
                preint.clone(), grav, Vector3::zeros(), Vector3::zeros(),
                att_i, Vector3::zeros(), Vector3::zeros(),
                att_j, Vector3::zeros(), Vector3::zeros(),
                prev_p, vi, pose_id, vj, b,
            );
            self.solver.graph.add_factor(Box::new(fac));
            let nhc = crate::swfg::pipeline::OdometerVelocityFactor::new(
                pose_id, vj, Vector3::zeros(), Vector3::new(25.0, 0.0025, 0.0025),
            );
            self.solver.graph.add_factor(Box::new(nhc));
            let speed = preint.dp.norm() / preint.dt.max(1e-3);
            if preint.dt > 0.05 && speed < 0.15 {
                let zupt = crate::swfg::pipeline::OdometerVelocityFactor::new(
                    pose_id, vj, Vector3::zeros(), Vector3::new(0.0001, 0.0001, 0.0001),
                );
                self.solver.graph.add_factor(Box::new(zupt));
            }
        }
    }

    fn setup_priors(
        &mut self, epoch: u32, pose_id: VariableId, prev_pose_id: Option<VariableId>,
        init_pos: Vector3<f64>, has_imu: bool,
    ) {
        let rot_axis = if has_imu {
            self.current_attitude.map(|q| q.scaled_axis()).unwrap_or_else(Vector3::zeros)
        } else { Vector3::zeros() };
        if epoch == 0 || prev_pose_id.is_none() || prev_pose_id != Some(pose_id) {
            self.solver.graph.set_value(pose_id, &[init_pos.x, init_pos.y, init_pos.z, rot_axis.x, rot_axis.y, rot_axis.z]);
        }

        if epoch == 0 || prev_pose_id.is_none() {
            let mut prior_info = nalgebra::DMatrix::zeros(6, 6);
            for i in 0..3 { prior_info[(i, i)] = 1.0 / 100_000.0; }
            if !has_imu { for i in 3..6 { prior_info[(i, i)] = 1.0; } }
            let pos_prior = crate::swfg::factor::PriorFactor {
                variable: pose_id,
                mu: nalgebra::DVector::from_row_slice(&[init_pos.x, init_pos.y, init_pos.z, 0.0, 0.0, 0.0]),
                information: prior_info,
            };
            self.solver.graph.add_factor(Box::new(pos_prior));
            if has_imu {
                if let Some(vel_id) = self.solver.graph.variables.iter()
                    .find(|(_, n)| matches!(n.kind, VariableKind::Velocity { epoch: e } if e == epoch))
                    .map(|(id, _)| *id)
                {
                    let vprior = crate::swfg::factor::PriorFactor::new(vel_id, nalgebra::DVector::zeros(3), 25.0);
                    self.solver.graph.add_factor(Box::new(vprior));
                }
            }
        } else if !has_imu && prev_pose_id != Some(pose_id) {
            let mut att_info = nalgebra::DMatrix::zeros(6, 6);
            for i in 3..6 { att_info[(i, i)] = 1.0; }
            let att_prior = crate::swfg::factor::PriorFactor {
                variable: pose_id,
                mu: nalgebra::DVector::from_row_slice(&[init_pos.x, init_pos.y, init_pos.z, 0.0, 0.0, 0.0]),
                information: att_info,
            };
            self.solver.graph.add_factor(Box::new(att_prior));
        }
    }

    fn setup_zwd(&mut self, epoch: u32, is_rtk: bool, dt_sec: f64) -> Option<VariableId> {
        if is_rtk { return None; }
        let id = self.solver.graph.variables.iter()
            .find(|(_, n)| matches!(n.kind, VariableKind::TropoZwd { epoch: e } if e == epoch))
            .map(|(id, _)| *id)?;
        self.solver.graph.set_value(id, &[0.1]);
        if let Some(prev_id) = self.prev_zwd {
            let rw_factor = crate::swfg::pipeline::tropo_factors::ZwdRandomWalkFactor::new(
                prev_id, id, dt_sec.max(1.0), 1e-8,
            );
            self.solver.graph.add_factor(Box::new(rw_factor));
        } else {
            let zwd_prior = crate::swfg::factor::PriorFactor::new(id, nalgebra::DVector::from_element(1, 0.1), 0.04);
            self.solver.graph.add_factor(Box::new(zwd_prior));
        }
        self.prev_zwd = Some(id);
        Some(id)
    }

    fn extract_obs(&self, obs: &EpochObs, pos: Vector3<f64>) -> Result<Vec<RawObservation>, String> {
        if let Some(ref orbits) = self.precise_orbits {
            let src = crate::estimators::rtk_iekf::satpos::PreciseSrc {
                orbits, clocks: self.precise_clocks.as_deref(),
            };
            epoch::extract_raw_observations_with_source(obs, &src, Some(pos), self.sinex_bias.as_deref())
        } else {
            epoch::extract_raw_observations(obs, &self.ephemerides, Some(pos))
        }
    }

    fn filter_observations(
        &self, rover: &EpochObs, base: Option<(&EpochObs, Vector3<f64>)>, init_pos: Vector3<f64>,
    ) -> Result<(Vec<RawObservation>, Option<Vec<RawObservation>>), String> {
        let mut raw_obs = self.extract_obs(rover, init_pos)?;
        let el_mask = self.elevation_mask_rad;
        raw_obs.retain(|r| r.elevation_rad >= el_mask && r.snr_dbhz >= 25.0);
        for r in &mut raw_obs { if r.snr_dbhz < 30.0 { r.cp_l1 = None; r.cp_l2 = None; } }

        let base_raw = if let Some((base_obs, base_pos)) = base {
            let mut b_raw = self.extract_obs(base_obs, base_pos)?;
            b_raw.retain(|r| r.elevation_rad >= el_mask && r.snr_dbhz >= 25.0);
            for r in &mut b_raw { if r.snr_dbhz < 30.0 { r.cp_l1 = None; r.cp_l2 = None; } }
            let base_sats: HashSet<(u8, u16)> = b_raw.iter().map(|b| (b.constellation_id, b.satellite)).collect();
            raw_obs.retain(|r| base_sats.contains(&(r.constellation_id, r.satellite)));
            Some(b_raw)
        } else { None };

        Ok((raw_obs, base_raw))
    }

    fn apply_corrections(&self, rover: &EpochObs, raw_obs: &[RawObservation], pos: Vector3<f64>) -> Vec<crate::swfg::pipeline::passes::CorrectedObservation> {
        let rx_state = ReceiverState {
            position_ecef: pos, clock_bias_m: vec![0.0], zwd_m: 0.1, ifb_glo: 0.0,
            llh_rad: gneiss_core::coords::ecef_to_llh(pos), time: rover.time,
        };
        let mut corrected = self.pipeline.process(raw_obs, &rx_state);
        corrected.retain(|o| o.pr_l1 > 0.0);
        corrected
    }

    fn seed_clocks(&mut self, epoch: u32, corrected: &[crate::swfg::pipeline::passes::CorrectedObservation], pos: Vector3<f64>, is_rtk: bool) {
        if is_rtk { return; }
        for c in corrected.iter().map(|o| o.constellation_id).collect::<BTreeSet<u8>>() {
            if let Some((c_id, _)) = self.solver.graph.variables.iter().find(|(_, n)| matches!(n.kind, VariableKind::ClockBias { epoch: e, constellation_id } if e == epoch && constellation_id == c)) {
                let c_id = *c_id;
                let mut offsets: Vec<f64> = corrected.iter().filter(|o| o.constellation_id == c).map(|obs| {
                    let pr = if let Some(pr2) = obs.pr_l2 {
                        let gamma = (obs.f1 / obs.f2.max(1.0)).powi(2);
                        (gamma * obs.pr_l1 - pr2) / (gamma - 1.0)
                    } else {
                        obs.pr_l1 - obs.iono_l1_m
                    };
                    pr - ((obs.sat_pos_ecef - pos).norm() - obs.sat_clock_m + obs.tropo_dry_m)
                }).collect();
                if !offsets.is_empty() {
                    offsets.sort_by(|a, b| a.total_cmp(b));
                    let med = offsets[offsets.len() / 2];
                    self.solver.graph.set_value(c_id, &[med]);
                    let prior = crate::swfg::factor::PriorFactor::new(c_id, nalgebra::DVector::from_element(1, med), 10_000.0);
                    self.solver.graph.add_factor(Box::new(prior));
                }
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn build_observation_factors(
        &mut self,
        epoch: u32,
        pose_id: VariableId,
        prev_pose_id: Option<VariableId>,
        init_pos: Vector3<f64>,
        rover_time: GpsTime,
        base: Option<(&EpochObs, Vector3<f64>)>,
        corrected: &[crate::swfg::pipeline::passes::CorrectedObservation],
        base_raw: &Option<Vec<RawObservation>>,
        zwd_id: Option<VariableId>,
    ) -> Result<Vec<builder::CpMeasurementRecord>, String> {
        let dt_sec = self.prev_time.map_or(1.0, |t| (rover_time.tow - t.tow).abs());
        if let Some((_base_obs, base_pos)) = base {
            let base_r = base_raw.as_ref().ok_or("base_raw missing")?;
            let ctx = RtkFactorContext {
                epoch, pose_id, prev_pose_id, init_pos,
                prev_position: self.prev_position, dt_sec, base_pos,
                corrected, base_raw: base_r,
            };
            let records = builder::build_rtk_dd_factors(
                &mut self.solver, &ctx, &mut self.ref_sat_per_constellation,
                &mut self.slip_counts, &mut self.mw_accumulator,
            );
            Ok(records)
        } else {
            builder::build_undifferenced_factors(
                &mut self.solver, corrected, epoch, pose_id, zwd_id, &mut self.slip_counts,
                &mut self.windup_trackers, rover_time, init_pos,
            );
            Ok(Vec::new())
        }
    }

    fn prune_unobserved_variables(&mut self) {
        let active_ids: HashSet<VariableId> = self.solver.graph.factors.iter()
            .flat_map(|f| f.variables()).copied().collect();
        let unref_ids: Vec<VariableId> = self.solver.graph.variables.iter()
            .filter(|(id, node)| {
                matches!(node.kind, VariableKind::Ambiguity { .. } | VariableKind::DdAmbiguity { .. } | VariableKind::ClockBias { .. })
                    && !active_ids.contains(id)
            })
            .map(|(id, _)| *id).collect();
        for id in unref_ids {
            self.solver.graph.remove_variable(id);
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn extract_and_finalize_solution(
        &mut self,
        pose_id: VariableId,
        init_pos: Vector3<f64>,
        time: GpsTime,
        n_sats: usize,
        has_imu: bool,
        is_rtk: bool,
        cp_records: &[builder::CpMeasurementRecord],
        epoch: u32,
    ) -> Result<SwfgSolution, String> {
        let vals = VariableValues::build(&self.solver.graph.variables);
        let pose = vals.get(pose_id).ok_or("pose missing after solve")?;
        let mut pos_ecef = Vector3::new(pose[0], pose[1], pose[2]);

        if is_rtk {
            ar_handler::check_postfit_cycle_slips(
                &mut self.solver, pos_ecef, cp_records, &mut self.slip_counts,
            );
        }
        if pos_ecef.x.is_nan() || pos_ecef.y.is_nan() || pos_ecef.z.is_nan() || pos_ecef.norm() < 1e6 {
            pos_ecef = self.prev_position.unwrap_or(init_pos);
        }

        if let Some(pose_val) = self.solver.graph.variables.get(&pose_id).map(|n| &n.value) {
            if pose_val.len() >= 6 && has_imu {
                let rot = Vector3::new(pose_val[3], pose_val[4], pose_val[5]);
                if rot.norm() > 1e-4 {
                    self.current_attitude = Some(nalgebra::UnitQuaternion::from_scaled_axis(rot));
                }
            }
        }

        self.prev_prev_position = self.prev_position;
        self.prev_position = Some(pos_ecef);
        self.prev_time = Some(time);

        let window_size = self.solver.window_size() as u32;
        if epoch >= window_size {
            let oldest = epoch - window_size;
            let _ = self.solver.marginalize_oldest_epoch(oldest);
        }

        Ok(SwfgSolution {
            time, position_ecef: pos_ecef, clock_bias_m: 0.0,
            n_satellites: n_sats, solver_iterations: 1, error: None,
        })
    }

    pub fn get_current_position(&self) -> Vector3<f64> {
        self.prev_position.unwrap_or_else(|| {
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0)
        })
    }

    pub fn set_klobuchar(&mut self, alpha: [f64; 4], beta: [f64; 4]) {
        self.pipeline.set_klobuchar(alpha, beta);
    }

    pub fn epoch_to_raw_obs(&self, rover: &EpochObs) -> Result<Vec<RawObservation>, String> {
        epoch::extract_raw_observations(rover, &self.ephemerides, self.initial_position)
    }
}


#[cfg(test)]
mod tests;
