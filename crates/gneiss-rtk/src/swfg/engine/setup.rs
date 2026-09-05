//! Epoch variable, prior, and IMU factor setup routines for SWFG engine.

use nalgebra::Vector3;
use crate::swfg::imu_preintegration::ImuPreintegration;
use crate::swfg::solver::SlidingWindowSolver;
use crate::swfg::variables::{VariableId, VariableKind};

fn init_attitude_from_pos(current: &mut Option<nalgebra::UnitQuaternion<f64>>, pos: Vector3<f64>) {
    let llh = gneiss_core::coords::ecef_to_llh(pos);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(llh).transpose();
    let rot = nalgebra::Rotation3::from_matrix_unchecked(ned_to_ecef);
    *current = Some(nalgebra::UnitQuaternion::from_rotation_matrix(&rot));
}

#[allow(clippy::too_many_arguments)]
pub fn setup_imu_and_rel_factors(
    solver: &mut SlidingWindowSolver,
    current_attitude: &mut Option<nalgebra::UnitQuaternion<f64>>,
    epoch: u32,
    pose_id: VariableId,
    prev_pose_id: Option<VariableId>,
    init_pos: Vector3<f64>,
    imu_preint: &Option<ImuPreintegration>,
    is_ppp: bool,
    is_kinematic: bool,
) {
    if imu_preint.is_some() && current_attitude.is_none() {
        init_attitude_from_pos(current_attitude, init_pos);
    }
    if let (Some(preint), Some(prev_p)) = (imu_preint, prev_pose_id) {
        add_imu_preintegration_factors(solver, current_attitude, epoch, pose_id, prev_p, init_pos, preint);
        *current_attitude = current_attitude.map(|q| q * preint.dq);
    } else if let Some(prev_p) = prev_pose_id {
        add_relative_pose_factor(solver, prev_p, pose_id, is_ppp, is_kinematic);
    }
}

fn add_relative_pose_factor(
    solver: &mut SlidingWindowSolver,
    prev_p: VariableId,
    pose_id: VariableId,
    is_ppp: bool,
    is_kinematic: bool,
) {
    if prev_p == pose_id || !solver.graph.variables.contains_key(&prev_p) {
        return;
    }
    let mut rel_info = nalgebra::DMatrix::zeros(6, 6);
    let q_pos = if is_ppp && !is_kinematic { 1.0 / 1e-4 } else { 1.0 / 25.0 };
    for i in 0..3 { rel_info[(i, i)] = q_pos; }
    for i in 3..6 { rel_info[(i, i)] = 1.0; }
    let rel = crate::swfg::factor::RelativePoseFactor { vars: [prev_p, pose_id], information: rel_info };
    solver.graph.add_factor(Box::new(rel));
}

fn add_imu_preintegration_factors(
    solver: &mut SlidingWindowSolver,
    current_attitude: &Option<nalgebra::UnitQuaternion<f64>>,
    epoch: u32,
    pose_id: VariableId,
    prev_p: VariableId,
    init_pos: Vector3<f64>,
    preint: &ImuPreintegration,
) {
    let vel_i = solver.graph.variables.iter()
        .find(|(_, n)| matches!(n.kind, VariableKind::Velocity { epoch: e } if e == epoch - 1))
        .map(|(id, _)| *id);
    let vel_j = solver.graph.variables.iter()
        .find(|(_, n)| matches!(n.kind, VariableKind::Velocity { epoch: e } if e == epoch))
        .map(|(id, _)| *id);
    let bias = solver.graph.variables.iter()
        .find(|(_, n)| matches!(n.kind, VariableKind::ImuBias))
        .map(|(id, _)| *id);

    if let (Some(vi), Some(vj), Some(b)) = (vel_i, vel_j, bias) {
        let att_i = current_attitude.unwrap_or_else(nalgebra::UnitQuaternion::identity);
        let att_j = att_i * preint.dq;
        let grav = -9.80665 * init_pos.normalize();
        let fac = crate::swfg::imu_preintegration::ImuPreintegrationFactor::new(
            preint.clone(), grav, Vector3::zeros(), Vector3::zeros(),
            att_i, Vector3::zeros(), Vector3::zeros(),
            att_j, Vector3::zeros(), Vector3::zeros(),
            prev_p, vi, pose_id, vj, b,
        );
        solver.graph.add_factor(Box::new(fac));
        add_motion_constraints(solver, pose_id, vj, preint);
    }
}

fn add_motion_constraints(
    solver: &mut SlidingWindowSolver,
    pose_id: VariableId,
    vj: VariableId,
    preint: &ImuPreintegration,
) {
    let nhc = crate::swfg::pipeline::OdometerVelocityFactor::new(
        pose_id, vj, Vector3::zeros(), Vector3::new(25.0, 0.0025, 0.0025),
    );
    solver.graph.add_factor(Box::new(nhc));
    let speed = preint.dp.norm() / preint.dt.max(1e-3);
    if preint.dt > 0.05 && speed < 0.15 {
        let zupt = crate::swfg::pipeline::OdometerVelocityFactor::new(
            pose_id, vj, Vector3::zeros(), Vector3::new(0.0001, 0.0001, 0.0001),
        );
        solver.graph.add_factor(Box::new(zupt));
    }
}

pub fn setup_priors(
    solver: &mut SlidingWindowSolver,
    current_attitude: &Option<nalgebra::UnitQuaternion<f64>>,
    epoch: u32,
    pose_id: VariableId,
    prev_pose_id: Option<VariableId>,
    init_pos: Vector3<f64>,
    has_imu: bool,
) {
    let rot_axis = if has_imu {
        current_attitude.map(|q| q.scaled_axis()).unwrap_or_else(Vector3::zeros)
    } else { Vector3::zeros() };
    if epoch == 0 || prev_pose_id.is_none() || prev_pose_id != Some(pose_id) {
        solver.graph.set_value(pose_id, &[init_pos.x, init_pos.y, init_pos.z, rot_axis.x, rot_axis.y, rot_axis.z]);
    }

    if epoch == 0 || prev_pose_id.is_none() {
        setup_initial_priors(solver, epoch, pose_id, init_pos, has_imu);
    } else if !has_imu && prev_pose_id != Some(pose_id) {
        add_attitude_prior(solver, pose_id, init_pos);
    }
}

fn setup_initial_priors(
    solver: &mut SlidingWindowSolver,
    epoch: u32,
    pose_id: VariableId,
    init_pos: Vector3<f64>,
    has_imu: bool,
) {
    let mut prior_info = nalgebra::DMatrix::zeros(6, 6);
    for i in 0..3 { prior_info[(i, i)] = 1.0 / 100_000.0; }
    if !has_imu { for i in 3..6 { prior_info[(i, i)] = 1.0; } }
    let pos_prior = crate::swfg::factor::PriorFactor {
        variable: pose_id,
        mu: nalgebra::DVector::from_row_slice(&[init_pos.x, init_pos.y, init_pos.z, 0.0, 0.0, 0.0]),
        information: prior_info,
    };
    solver.graph.add_factor(Box::new(pos_prior));
    if has_imu {
        add_initial_velocity_prior(solver, epoch);
    }
}

fn add_initial_velocity_prior(solver: &mut SlidingWindowSolver, epoch: u32) {
    if let Some(vel_id) = solver.graph.variables.iter()
        .find(|(_, n)| matches!(n.kind, VariableKind::Velocity { epoch: e } if e == epoch))
        .map(|(id, _)| *id)
    {
        let vprior = crate::swfg::factor::PriorFactor::new(vel_id, nalgebra::DVector::zeros(3), 25.0);
        solver.graph.add_factor(Box::new(vprior));
    }
}

fn add_attitude_prior(solver: &mut SlidingWindowSolver, pose_id: VariableId, init_pos: Vector3<f64>) {
    let mut att_info = nalgebra::DMatrix::zeros(6, 6);
    for i in 3..6 { att_info[(i, i)] = 1.0; }
    let att_prior = crate::swfg::factor::PriorFactor {
        variable: pose_id,
        mu: nalgebra::DVector::from_row_slice(&[init_pos.x, init_pos.y, init_pos.z, 0.0, 0.0, 0.0]),
        information: att_info,
    };
    solver.graph.add_factor(Box::new(att_prior));
}
