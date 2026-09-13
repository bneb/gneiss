//! N-Pass Iterative Calibration Architecture for Extrinsics, Intrinsics, Biases & Offsets.
//!
//! Iteratively executes the post-processing engine over passes 1..N:
//! 1. Runs forward-backward filtering and smoothing.
//! 2. Extracts dynamic maneuvering intervals (for antenna-to-IMU lever arm).
//! 3. Extracts stationary intervals (for IMU accelerometer and gyro biases).
//! 4. Extracts fixed-epoch residual offsets (for antenna phase center body offsets).
//! 5. Evaluates convergence criteria against tolerances.
//! 6. Updates calibration parameters and re-processes until convergence.

use std::collections::BTreeMap;
use nalgebra::{Matrix3, Vector3};

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use crate::post_process::lever_arm::{LeverArmEstimator, LeverArmObservation};
use crate::post_process::{
    execute_post_process, PostProcessOptions, PostProcessResult, SmoothedEpoch,
};
use crate::swfg::config::EngineConfig;
use crate::swfg::imu_preintegration::stationary::compute_stationary_metrics;
use crate::swfg::imu_preintegration::ImuSample;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct CalibrationParameters {
    pub lever_arm_body: Option<Vector3<f64>>,
    pub lever_arm_std: Option<Vector3<f64>>,
    pub imu_accel_bias: Option<Vector3<f64>>,
    pub imu_gyro_bias: Option<Vector3<f64>>,
    pub antenna_body_offset: Option<Vector3<f64>>,
    pub boresight_rpy_rad: Option<Vector3<f64>>,
}

/// Convergence criteria for the multi-pass calibration loop.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CalibrationConvergenceCriteria {
    pub max_iterations: usize,
    pub lever_arm_tol_m: f64,
    pub accel_bias_tol_mps2: f64,
    pub gyro_bias_tol_radps: f64,
}

impl Default for CalibrationConvergenceCriteria {
    fn default() -> Self {
        Self {
            max_iterations: 3,
            lever_arm_tol_m: 0.005,
            accel_bias_tol_mps2: 0.005,
            gyro_bias_tol_radps: 0.0005,
        }
    }
}

/// Map from timestamp in milliseconds to (ECEF ground position, [roll, pitch, yaw] in radians).
pub type ReferencePointMap = BTreeMap<i64, (Vector3<f64>, [f64; 3])>;

/// Options controlling the multi-pass calibration engine.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiPassCalibrationOptions {
    pub criteria: CalibrationConvergenceCriteria,
    pub estimate_lever_arm: bool,
    pub estimate_imu_biases: bool,
    pub estimate_antenna_offset: bool,
    pub reference_points: Option<ReferencePointMap>,
}

impl Default for MultiPassCalibrationOptions {
    fn default() -> Self {
        Self {
            criteria: CalibrationConvergenceCriteria::default(),
            estimate_lever_arm: true,
            estimate_imu_biases: true,
            estimate_antenna_offset: true,
            reference_points: None,
        }
    }
}

/// Record of one iteration pass during calibration.
#[derive(Debug, Clone, PartialEq)]
pub struct CalibrationIterationRecord {
    pub pass_number: usize,
    pub parameters: CalibrationParameters,
    pub lever_arm_delta_m: f64,
    pub accel_bias_delta_mps2: f64,
    pub gyro_bias_delta_radps: f64,
    pub fixed_epochs: usize,
    pub total_epochs: usize,
}

/// Comprehensive report of the multi-pass calibration procedure.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiPassCalibrationReport {
    pub passes_executed: usize,
    pub converged: bool,
    pub initial_parameters: CalibrationParameters,
    pub final_parameters: CalibrationParameters,
    pub history: Vec<CalibrationIterationRecord>,
}

/// Pre-estimate static IMU biases from stationary segments in the raw IMU data.
pub fn estimate_static_imu_biases(
    samples: &[ImuSample],
) -> (Option<Vector3<f64>>, Option<Vector3<f64>>) {
    if samples.len() < 20 {
        return (None, None);
    }
    let mut gyro_sum = Vector3::zeros();
    let mut accel_norm_sum = 0.0;
    let mut count = 0usize;
    for chunk in samples.chunks(20) {
        if let Some(m) = compute_stationary_metrics(chunk) {
            if m.is_stationary() {
                for s in chunk {
                    gyro_sum += s.gyro;
                    accel_norm_sum += s.accel.norm();
                    count += 1;
                }
            }
        }
    }
    if count < 20 {
        return (None, None);
    }
    let c_f = count as f64;
    let gyro_bias = gyro_sum / c_f;
    let accel_bias_z = accel_norm_sum / c_f - 9.80665;
    (Some(Vector3::new(0.0, 0.0, accel_bias_z)), Some(gyro_bias))
}

/// Build DCM from vehicle velocity (course-over-ground) assuming flat/level driving.
pub fn heading_pitch_to_dcm(vel_ecef: Vector3<f64>, pos_ecef: Vector3<f64>) -> Option<Matrix3<f64>> {
    let speed = vel_ecef.norm();
    if speed < 1.0 {
        return None;
    }
    let llh = gneiss_core::coords::ecef_to_llh(pos_ecef);
    let e_to_n = gneiss_core::coords::ecef_to_ned_matrix(llh);
    let v_ned = e_to_n * vel_ecef;
    let heading = libm::atan2(v_ned.y, v_ned.x);
    let pitch = libm::atan2(-v_ned.z, libm::sqrt(v_ned.x * v_ned.x + v_ned.y * v_ned.y));
    let (sh, ch) = (libm::sin(heading), libm::cos(heading));
    let (sp, cp) = (libm::sin(pitch), libm::cos(pitch));
    let r_b_n = Matrix3::new(
        cp * ch, -sh, sp * ch,
        cp * sh,  ch, sp * sh,
        -sp,     0.0, cp,
    );
    Some(e_to_n.transpose() * r_b_n)
}

/// Compute single lever arm observation from trajectory and IMU window.
fn build_lever_arm_obs(
    ep_prev: &SmoothedEpoch,
    ep_curr: &SmoothedEpoch,
    ep_next: &SmoothedEpoch,
    imu_samples: &[ImuSample],
) -> Option<LeverArmObservation> {
    let dt = ep_next.time.tow - ep_prev.time.tow;
    if dt <= 0.01 || dt > 2.0 {
        return None;
    }
    let v_prev = ep_prev.velocity_ecef?;
    let v_next = ep_next.velocity_ecef?;
    let a_gnss_ecef = (v_next - v_prev) / dt;
    let r_b2e = ep_curr.attitude.map(|q| q.to_rotation_matrix().into_inner())
        .or_else(|| heading_pitch_to_dcm(v_next, ep_curr.position_ecef))?;
    let g_ecef = -9.80665 * ep_curr.position_ecef.normalize();
    let a_gnss_body = r_b2e.transpose() * (a_gnss_ecef - g_ecef);

    let (omega, alpha, accel_imu) = summarize_imu_window(ep_curr.time.tow, imu_samples)?;
    Some(LeverArmObservation {
        omega_body: omega,
        alpha_body: alpha,
        accel_imu_body: accel_imu,
        accel_gnss_body: a_gnss_body,
    })
}

/// Extract IMU angular rate, angular acceleration, and specific force near tow.
fn summarize_imu_window(
    tow: f64,
    samples: &[ImuSample],
) -> Option<(Vector3<f64>, Vector3<f64>, Vector3<f64>)> {
    let tow_us = (tow * 1_000_000.0).round() as u64;
    let idx = samples.binary_search_by_key(&tow_us, |s| s.time_us)
        .unwrap_or_else(|i| i);
    if idx < 2 || idx + 2 >= samples.len() {
        return None;
    }
    let s_prev = &samples[idx - 1];
    let s_curr = &samples[idx];
    let s_next = &samples[idx + 1];
    let dt_us = (s_next.time_us as f64 - s_prev.time_us as f64) * 1e-6;
    if dt_us <= 1e-4 {
        return None;
    }
    let alpha = (s_next.gyro - s_prev.gyro) / dt_us;
    Some((s_curr.gyro, alpha, s_curr.accel))
}

/// Extract dynamic lever arm observations across turning maneuvers in the trajectory.
pub fn extract_lever_arm_observations(
    trajectory: &[SmoothedEpoch],
    imu_samples: &[ImuSample],
) -> Vec<LeverArmObservation> {
    if trajectory.len() < 3 || imu_samples.is_empty() {
        return Vec::new();
    }
    let mut obs = Vec::new();
    for i in 1..(trajectory.len() - 1) {
        if let Some(o) = build_lever_arm_obs(&trajectory[i - 1], &trajectory[i], &trajectory[i + 1], imu_samples) {
            if o.omega_body.norm() > 0.04 || o.alpha_body.norm() > 0.04 {
                obs.push(o);
            }
        }
    }
    obs
}

/// Estimate residual body-frame antenna offset from fixed RTK epochs against reference.
pub fn estimate_body_antenna_offset(
    trajectory: &[SmoothedEpoch],
    reference_points: &ReferencePointMap,
) -> Option<Vector3<f64>> {
    let mut sum_body = Vector3::zeros();
    let mut count = 0usize;
    for ep in trajectory {
        if ep.quality != 1 {
            continue;
        }
        let tow_ms = (ep.time.tow * 1000.0).round() as i64;
        if let Some((_, &(ref_pos, rpy))) = reference_points.range((tow_ms - 150)..=(tow_ms + 150)).next() {
            let r_b_e = rpy_to_rbe(ref_pos, rpy[0], rpy[1], rpy[2]);
            let diff_e = ep.position_ecef - ref_pos;
            sum_body += r_b_e.transpose() * diff_e;
            count += 1;
        }
    }
    if count >= 10 {
        Some(sum_body / (count as f64))
    } else {
        None
    }
}

fn rpy_to_rbe(pos: Vector3<f64>, roll: f64, pitch: f64, yaw: f64) -> Matrix3<f64> {
    let (sr, cr) = (libm::sin(roll), libm::cos(roll));
    let (sp, cp) = (libm::sin(pitch), libm::cos(pitch));
    let (sy, cy) = (libm::sin(yaw), libm::cos(yaw));
    let r_b_n = Matrix3::new(
        cp * cy, sr * sp * cy - cr * sy, cr * sp * cy + sr * sy,
        cp * sy, sr * sp * sy + cr * cy, cr * sp * sy - sr * cy,
        -sp,     sr * cp,                cr * cp,
    );
    let llh = gneiss_core::coords::ecef_to_llh(pos);
    let e_to_n = gneiss_core::coords::ecef_to_ned_matrix(llh);
    e_to_n.transpose() * r_b_n
}

/// Apply calibrated lever arm offset to output trajectory.
pub fn apply_calibration_to_trajectory(
    trajectory: &mut [SmoothedEpoch],
    offset_body: Vector3<f64>,
) {
    if offset_body.norm() < 1e-5 {
        return;
    }
    for ep in trajectory.iter_mut() {
        let r_b2e = ep.attitude.map(|q| q.to_rotation_matrix().into_inner())
            .or_else(|| {
                ep.velocity_ecef.and_then(|v| heading_pitch_to_dcm(v, ep.position_ecef))
            });
        if let Some(r) = r_b2e {
            ep.position_ecef -= r * offset_body;
        }
    }
}

/// Check parameter convergence across successive passes.
fn check_convergence(
    prev: &CalibrationParameters,
    curr: &CalibrationParameters,
    crit: &CalibrationConvergenceCriteria,
) -> (bool, f64, f64, f64) {
    let d_arm = match (prev.lever_arm_body, curr.lever_arm_body) {
        (Some(p), Some(c)) => (c - p).norm(),
        _ => 0.0,
    };
    let d_accel = match (prev.imu_accel_bias, curr.imu_accel_bias) {
        (Some(p), Some(c)) => (c - p).norm(),
        _ => 0.0,
    };
    let d_gyro = match (prev.imu_gyro_bias, curr.imu_gyro_bias) {
        (Some(p), Some(c)) => (c - p).norm(),
        _ => 0.0,
    };
    let conv = d_arm <= crit.lever_arm_tol_m
        && d_accel <= crit.accel_bias_tol_mps2
        && d_gyro <= crit.gyro_bias_tol_radps;
    (conv, d_arm, d_accel, d_gyro)
}

/// Perform one calibration analysis step on a completed post-process result.
fn analyze_iteration_result(
    result: &PostProcessResult,
    imu_samples: Option<&[ImuSample]>,
    opts: &MultiPassCalibrationOptions,
    current: &mut CalibrationParameters,
) {
    if opts.estimate_lever_arm {
        if let Some(samples) = imu_samples {
            let obs = extract_lever_arm_observations(&result.trajectory, samples);
            if let Some(est) = LeverArmEstimator::estimate(&obs) {
                current.lever_arm_body = Some(est.lever_arm_body);
                current.lever_arm_std = Some(est.std_body);
            }
        }
    }
    if opts.estimate_antenna_offset {
        if let Some(ref refs) = opts.reference_points {
            if let Some(off) = estimate_body_antenna_offset(&result.trajectory, refs) {
                current.antenna_body_offset = Some(off);
                if current.lever_arm_body.is_none() {
                    current.lever_arm_body = Some(off);
                }
            }
        }
    }
}

fn init_calibration_params(
    calib_opts: &MultiPassCalibrationOptions,
    imu_samples: Option<&[ImuSample]>,
) -> CalibrationParameters {
    let mut params = CalibrationParameters::default();
    if calib_opts.estimate_imu_biases {
        if let Some(samples) = imu_samples {
            let (ba, bg) = estimate_static_imu_biases(samples);
            params.imu_accel_bias = ba;
            params.imu_gyro_bias = bg;
        }
    }
    params
}

fn record_iteration(
    pass: usize,
    curr: &CalibrationParameters,
    prev: &CalibrationParameters,
    crit: &CalibrationConvergenceCriteria,
    res: &PostProcessResult,
) -> (CalibrationIterationRecord, bool) {
    let fixed_count = res.trajectory.iter().filter(|e| e.quality == 1).count();
    let (conv, d_arm, d_acc, d_gyr) = check_convergence(prev, curr, crit);
    let rec = CalibrationIterationRecord {
        pass_number: pass,
        parameters: curr.clone(),
        lever_arm_delta_m: d_arm,
        accel_bias_delta_mps2: d_acc,
        gyro_bias_delta_radps: d_gyr,
        fixed_epochs: fixed_count,
        total_epochs: res.trajectory.len(),
    };
    (rec, conv)
}

#[allow(clippy::too_many_arguments)]
fn run_calibration_loop(
    config: &EngineConfig,
    ephems: &[Ephemeris],
    rov: &[EpochObs],
    base: Option<&[EpochObs]>,
    imu: Option<&[ImuSample]>,
    opts: &PostProcessOptions,
    calib: &MultiPassCalibrationOptions,
    curr: &mut CalibrationParameters,
) -> Result<(PostProcessResult, Vec<CalibrationIterationRecord>, bool), String> {
    let mut history = Vec::new();
    let mut prev = curr.clone();
    let mut last_res = None;
    let mut converged = false;
    for pass in 1..=calib.criteria.max_iterations.max(1) {
        let mut pass_opts = opts.clone();
        pass_opts.calibration = Some(curr.clone());
        let res = execute_post_process(config, ephems, rov, base, imu, &pass_opts)?;
        analyze_iteration_result(&res, imu, calib, curr);
        let (rec, conv) = record_iteration(pass, curr, &prev, &calib.criteria, &res);
        history.push(rec);
        prev = curr.clone();
        last_res = Some(res);
        if pass > 1 && conv {
            converged = true;
            break;
        }
    }
    let res = last_res.ok_or_else(|| "No passes executed".to_string())?;
    Ok((res, history, converged))
}

/// Execute the multi-pass calibrated post-processing pipeline.
pub fn execute_calibrated_post_process(
    config: &EngineConfig,
    ephemerides: &[Ephemeris],
    rover_epochs: &[EpochObs],
    base_epochs: Option<&[EpochObs]>,
    imu_samples: Option<&[ImuSample]>,
    base_options: &PostProcessOptions,
    calib_opts: &MultiPassCalibrationOptions,
) -> Result<(PostProcessResult, MultiPassCalibrationReport), String> {
    let mut curr_params = init_calibration_params(calib_opts, imu_samples);
    let initial_params = curr_params.clone();
    let (mut res, history, conv) = run_calibration_loop(
        config, ephemerides, rover_epochs, base_epochs, imu_samples,
        base_options, calib_opts, &mut curr_params,
    )?;
    if let Some(offset) = curr_params.antenna_body_offset.or(curr_params.lever_arm_body) {
        apply_calibration_to_trajectory(&mut res.trajectory, offset);
    }
    let report = MultiPassCalibrationReport {
        passes_executed: history.len(),
        converged: conv,
        initial_parameters: initial_params,
        final_parameters: curr_params,
        history,
    };
    Ok((res, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::UnitQuaternion;

    #[test]
    fn test_static_imu_bias_estimation() {
        let mut samples = Vec::new();
        for i in 0..60 {
            samples.push(ImuSample {
                accel: Vector3::new(0.01, -0.02, 9.85),
                gyro: Vector3::new(0.005, -0.002, 0.001),
                time_us: i * 50_000,
            });
        }
        let (ba, bg) = estimate_static_imu_biases(&samples);
        let gyro = bg.expect("gyro bias estimated");
        assert!((gyro.x - 0.005).abs() < 1e-4);
        assert!((gyro.y + 0.002).abs() < 1e-4);
        let accel = ba.expect("accel bias estimated");
        assert!((accel.z - (9.85 - 9.80665)).abs() < 1e-3);
    }

    #[test]
    fn test_convergence_criteria() {
        let crit = CalibrationConvergenceCriteria::default();
        let p1 = CalibrationParameters {
            lever_arm_body: Some(Vector3::new(0.50, 0.20, -0.10)),
            ..Default::default()
        };
        let p2 = CalibrationParameters {
            lever_arm_body: Some(Vector3::new(0.502, 0.201, -0.101)),
            ..Default::default()
        };
        let (conv, d_arm, _, _) = check_convergence(&p1, &p2, &crit);
        assert!(conv);
        assert!(d_arm < 0.005);
    }

    #[test]
    fn test_apply_calibration_to_trajectory() {
        let mut traj = vec![SmoothedEpoch {
            time: gneiss_core::time::GpsTime::new(2000, 100.0),
            position_ecef: Vector3::new(100.0, 200.0, 300.0),
            velocity_ecef: Some(Vector3::new(10.0, 0.0, 0.0)),
            attitude: Some(UnitQuaternion::identity()),
            cov_position: Matrix3::identity(),
            std_east: 0.01,
            std_north: 0.01,
            std_up: 0.02,
            separation_3d: 0.01,
            quality: 1,
            n_satellites: 8,
        }];
        apply_calibration_to_trajectory(&mut traj, Vector3::new(1.0, 0.5, -0.2));
        assert!((traj[0].position_ecef.x - 99.0).abs() < 1e-4);
        assert!((traj[0].position_ecef.y - 199.5).abs() < 1e-4);
        assert!((traj[0].position_ecef.z - 300.2).abs() < 1e-4);
    }
}
