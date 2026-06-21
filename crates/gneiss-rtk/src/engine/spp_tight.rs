use gneiss_core::obs::EpochObs;
use nalgebra::{DMatrix, DVector, Vector3};

use crate::engine::{EngineError, ProcessingEngine};
use crate::filter::RtkState;
use crate::spp::{build_measurements, SppConfig, SppMeasurement};
use gneiss_core::atmosphere::{AtmosphereModel, TropoParams};
use gneiss_core::constants::{EARTH_ROTATION_RATE_RAD_S, SPEED_OF_LIGHT_M_S};
use gneiss_core::coords::{az_el, ecef_to_llh};
use gneiss_core::sat::Constellation;

const MIN_GEOM_RANGE: f64 = 1e-6;
const MIN_ELEVATION_RAD: f64 = 5.0 * core::f64::consts::PI / 180.0;
const MAX_CONSECUTIVE_REJECTIONS: usize = 10;
const CLK_BIAS_VAR_RESET: f64 = 1e6;
const COV_POS_RESET: f64 = 100.0;
const COV_VEL_RESET: f64 = 10.0;
const COV_CLK_RESET: f64 = 100000.0;
const COV_DRIFT_RESET: f64 = 1000.0;
const MIN_VALID_MEASUREMENTS: usize = 3;

fn compute_sagnac_correction(sat_ecef: Vector3<f64>, geometric_pr: f64) -> Vector3<f64> {
    let tof = geometric_pr / SPEED_OF_LIGHT_M_S;
    let theta = EARTH_ROTATION_RATE_RAD_S * tof;
    let cos_t = f64::cos(theta);
    let sin_t = f64::sin(theta);
    Vector3::new(
        sat_ecef.x * cos_t + sat_ecef.y * sin_t,
        -sat_ecef.x * sin_t + sat_ecef.y * cos_t,
        sat_ecef.z,
    )
}

fn compute_sagnac_velocity_correction(sat_vel: Vector3<f64>, geometric_pr: f64) -> Vector3<f64> {
    let tof = geometric_pr / SPEED_OF_LIGHT_M_S;
    let theta = EARTH_ROTATION_RATE_RAD_S * tof;
    let cos_t = f64::cos(theta);
    let sin_t = f64::sin(theta);
    Vector3::new(
        sat_vel.x * cos_t + sat_vel.y * sin_t,
        -sat_vel.x * sin_t + sat_vel.y * cos_t,
        sat_vel.z,
    )
}

pub fn process_spp_tightly_coupled<'a>(
    engine: &'a mut ProcessingEngine,
    rover_obs: &'a EpochObs,
) -> Result<&'a RtkState, EngineError> {
    if engine.current_state.is_none() {
        return engine.process_spp(rover_obs);
    }
    predict_and_align_state(engine, rover_obs);
    let measurements: Vec<_> =
        build_measurements(rover_obs, &engine.ephemerides, &SppConfig::default())
            .into_iter()
            .filter(|m| m.snr >= 25.0)
            .collect();
    if measurements.is_empty() {
        return Err(EngineError::NoObservations);
    }
    bootstrap_clock_bias(engine, rover_obs);
    let (z_vec, h_mat, r_mat, meas_types) = build_ekf_matrices(engine, &measurements);
    let rejected = update_ekf(engine, &z_vec, &h_mat, &r_mat, &meas_types);
    handle_rejection(engine, rover_obs, rejected)?;
    finalize_epoch(engine, rover_obs);
    Ok(engine.current_state.as_ref().unwrap())
}

fn predict_and_align_state(engine: &mut ProcessingEngine, rover_obs: &EpochObs) {
    let dt = rover_obs.time.tow - engine.current_state.as_ref().unwrap().time.tow;
    engine.predict_state(dt);
    let state = engine.current_state.as_mut().unwrap();
    state.time = rover_obs.time;
    state.position.epoch = rover_obs.time;
}

fn bootstrap_clock_bias(engine: &mut ProcessingEngine, rover_obs: &EpochObs) {
    if let Ok(spp_res) = crate::spp::compute_spp(
        rover_obs,
        &engine.ephemerides,
        engine.klobuchar_params.as_ref(),
        &SppConfig::default(),
        None,
    ) {
        let state = engine.current_state.as_mut().unwrap();
        state.rcv_clk_bias = spp_res.cdt;
        state.covariance[(15, 15)] = CLK_BIAS_VAR_RESET;
    }
}

fn get_apc_kinematics(
    engine: &ProcessingEngine,
) -> (
    Vector3<f64>,
    Vector3<f64>,
    nalgebra::Rotation3<f64>,
    Vector3<f64>,
    Vector3<f64>,
) {
    let state = engine.current_state.as_ref().unwrap();
    let r_b_e = state.attitude.to_rotation_matrix();
    let lever_arm = if state.ins_aligned {
        Vector3::from_column_slice(&engine.config.imu_to_antenna_lever_arm)
    } else {
        Vector3::zeros()
    };
    let l_e = r_b_e * lever_arm;
    let pos_apc = state.position.vector + l_e;

    let omega_b = engine
        .imu_history
        .last()
        .and_then(|buf: &Vec<gneiss_core::imu::ImuMeasurement>| buf.last())
        .map(|last_imu| last_imu.gyro - state.gyro_bias)
        .unwrap_or_else(Vector3::zeros);

    let omega_ie_e = Vector3::new(0.0, 0.0, EARTH_ROTATION_RATE_RAD_S);
    let omega_eb_b = omega_b - r_b_e.transpose() * omega_ie_e;
    let v_apc = state.velocity + r_b_e * omega_eb_b.cross(&lever_arm);
    (pos_apc, v_apc, r_b_e, lever_arm, omega_eb_b)
}

struct EkfContext<'a> {
    pos_apc: Vector3<f64>,
    v_apc: Vector3<f64>,
    r_b_e: nalgebra::Rotation3<f64>,
    lever_arm: Vector3<f64>,
    omega_eb_b: Vector3<f64>,
    rec_llh: Vector3<f64>,
    state: &'a RtkState,
    engine: &'a ProcessingEngine,
    n_cols: usize,
}

struct MatrixTarget<'a> {
    z: &'a mut DVector<f64>,
    h: &'a mut DMatrix<f64>,
    r: &'a mut DMatrix<f64>,
    types: &'a mut Vec<(gneiss_core::sat::SatelliteId, u8)>,
    row: &'a mut usize,
}

struct SatGeometry {
    geom_r: f64,
    los: Vector3<f64>,
    el: f64,
    az: f64,
    cdt_rx: f64,
    sat_clk: f64,
    sat_vel: Vector3<f64>,
    sat_drift: f64,
}

fn build_ekf_matrices(
    engine: &ProcessingEngine,
    measurements: &[SppMeasurement],
) -> (
    DVector<f64>,
    DMatrix<f64>,
    DMatrix<f64>,
    Vec<(gneiss_core::sat::SatelliteId, u8)>,
) {
    let (pos_apc, v_apc, r_b_e, lever_arm, omega_eb_b) = get_apc_kinematics(engine);
    let ctx = EkfContext {
        pos_apc,
        v_apc,
        r_b_e,
        lever_arm,
        omega_eb_b,
        rec_llh: ecef_to_llh(pos_apc),
        state: engine.current_state.as_ref().unwrap(),
        engine,
        n_cols: engine.current_state.as_ref().unwrap().covariance.ncols(),
    };

    let num_dop = measurements.iter().filter(|m| m.doppler != 0.0).count();
    let total_rows = measurements.len() + num_dop;

    let mut z_vec = DVector::zeros(total_rows);
    let mut h_mat = DMatrix::zeros(total_rows, ctx.n_cols);
    let mut r_mat = DMatrix::zeros(total_rows, total_rows);
    let mut m_types = Vec::with_capacity(total_rows);

    let mut row_idx = 0;
    for m in measurements {
        let mut target = MatrixTarget {
            z: &mut z_vec,
            h: &mut h_mat,
            r: &mut r_mat,
            types: &mut m_types,
            row: &mut row_idx,
        };
        process_measurement(&ctx, m, &mut target);
    }
    (z_vec, h_mat, r_mat, m_types)
}

fn get_cdt_rx(state: &RtkState, constellation: Constellation) -> f64 {
    match constellation {
        Constellation::Gps => state.rcv_clk_bias,
        Constellation::Galileo => state.rcv_clk_bias + state.isb_gal,
        Constellation::Beidou => state.rcv_clk_bias + state.isb_bds,
        Constellation::Glonass => state.rcv_clk_bias + state.isb_glo,
        _ => state.rcv_clk_bias,
    }
}

fn process_measurement(ctx: &EkfContext, m: &SppMeasurement, target: &mut MatrixTarget) {
    let cdt_rx = get_cdt_rx(ctx.state, m.constellation);
    let t_rcv_true = m.time.tow - cdt_rx / SPEED_OF_LIGHT_M_S;
    let t_tx_nom =
        gneiss_core::time::GpsTime::new(m.time.week, t_rcv_true - m.raw_pr / SPEED_OF_LIGHT_M_S);
    let (_, _, dt_s, _) = m.eph.position(t_tx_nom);
    let t_tx_true = gneiss_core::time::GpsTime::new(
        m.time.week,
        t_rcv_true - m.raw_pr / SPEED_OF_LIGHT_M_S - dt_s,
    );
    let (sat_pos_raw, sat_vel_raw, sat_clk, sat_drift) = m.eph.position(t_tx_true);

    let sat_ecef = compute_sagnac_correction(sat_pos_raw, m.raw_pr - cdt_rx);
    let sat_vel = compute_sagnac_velocity_correction(sat_vel_raw, m.raw_pr - cdt_rx);

    let dx = ctx.pos_apc.x - sat_ecef.x;
    let dy = ctx.pos_apc.y - sat_ecef.y;
    let dz = ctx.pos_apc.z - sat_ecef.z;
    let geom_r = f64::sqrt(dx * dx + dy * dy + dz * dz).max(MIN_GEOM_RANGE);
    let los = Vector3::new(dx / geom_r, dy / geom_r, dz / geom_r);
    let (az, el) = az_el(ctx.rec_llh, ctx.pos_apc, sat_ecef);

    let geom = SatGeometry {
        geom_r,
        los,
        el,
        az,
        cdt_rx,
        sat_clk,
        sat_vel,
        sat_drift,
    };
    process_pseudorange(ctx, m, target, &geom);
    if m.doppler != 0.0 {
        process_doppler(ctx, m, target, &geom);
    }
}

fn populate_pr_attitude_jacobian(
    ctx: &EkfContext,
    los: &Vector3<f64>,
    h: &mut DMatrix<f64>,
    r_idx: usize,
) {
    let h_pos_att = -(ctx.r_b_e * ctx.lever_arm).cross_matrix();
    let pr_h_att = los.transpose() * h_pos_att;
    h[(r_idx, 6)] = pr_h_att[0];
    h[(r_idx, 7)] = pr_h_att[1];
    h[(r_idx, 8)] = pr_h_att[2];
}

fn populate_pr_clock_jacobian(h: &mut DMatrix<f64>, r_idx: usize, constel: Constellation) {
    h[(r_idx, 15)] = 1.0;
    match constel {
        Constellation::Glonass => h[(r_idx, 16)] = 1.0,
        Constellation::Galileo => h[(r_idx, 17)] = 1.0,
        Constellation::Beidou => h[(r_idx, 18)] = 1.0,
        _ => {}
    }
}

fn process_pseudorange(
    ctx: &EkfContext,
    m: &SppMeasurement,
    target: &mut MatrixTarget,
    geom: &SatGeometry,
) {
    let safe_el = geom.el.max(MIN_ELEVATION_RAD);
    let tropo = AtmosphereModel::tropo_nmf(&TropoParams::default(), ctx.rec_llh, safe_el, m.time);
    let iono = ctx.engine.klobuchar_params.as_ref().map_or(0.0, |p| {
        AtmosphereModel::iono_klobuchar(p, ctx.rec_llh, geom.az, safe_el, m.time)
    });
    let expected_pr = geom.geom_r + geom.cdt_rx - geom.sat_clk * SPEED_OF_LIGHT_M_S + tropo + iono;
    let r_idx = *target.row;
    target.z[r_idx] = m.raw_pr - expected_pr;
    for i in 0..3 {
        target.h[(r_idx, i)] = geom.los[i];
    }
    if ctx.n_cols > 15 {
        populate_pr_attitude_jacobian(ctx, &geom.los, target.h, r_idx);
        populate_pr_clock_jacobian(target.h, r_idx, m.constellation);
    }
    let v_scale = gneiss_core::variance::observation_variance(
        m.snr,
        geom.el,
        ctx.engine.config.tuning.snr_a,
        ctx.engine.config.tuning.snr_b,
    );
    target.r[(r_idx, r_idx)] = ctx.engine.config.tuning.pr_base_var * v_scale;
    target.types.push((m.eph.sat(), 0));
    *target.row += 1;
}

fn populate_dop_jacobian(ctx: &EkfContext, los: &Vector3<f64>, h: &mut DMatrix<f64>, r_idx: usize) {
    let a_0 = ctx.r_b_e * ctx.omega_eb_b.cross(&ctx.lever_arm);
    let h_vel_att = -a_0.cross_matrix();
    let h_vel_bg = ctx.r_b_e.matrix() * ctx.lever_arm.cross_matrix();
    let dop_h_att = los.transpose() * h_vel_att;
    let dop_h_bg = los.transpose() * h_vel_bg;
    h[(r_idx, 6)] = dop_h_att[0];
    h[(r_idx, 7)] = dop_h_att[1];
    h[(r_idx, 8)] = dop_h_att[2];
    h[(r_idx, 12)] = dop_h_bg[0];
    h[(r_idx, 13)] = dop_h_bg[1];
    h[(r_idx, 14)] = dop_h_bg[2];
    h[(r_idx, 19)] = 1.0;
}

fn process_doppler(
    ctx: &EkfContext,
    m: &SppMeasurement,
    target: &mut MatrixTarget,
    geom: &SatGeometry,
) {
    let rel_vel = ctx.v_apc - geom.sat_vel;
    let expected_dop =
        geom.los.dot(&rel_vel) + ctx.state.rcv_clk_drift - geom.sat_drift * SPEED_OF_LIGHT_M_S;
    let f1 = gneiss_core::signal::satellite_frequencies(m.eph.sat(), m.eph.freq_num()).0;
    let measured_dop_ms = -m.doppler * (SPEED_OF_LIGHT_M_S / f1);
    let r_idx = *target.row;
    target.z[r_idx] = measured_dop_ms - expected_dop;
    target.h[(r_idx, 3)] = geom.los.x;
    target.h[(r_idx, 4)] = geom.los.y;
    target.h[(r_idx, 5)] = geom.los.z;
    if ctx.n_cols > 19 {
        populate_dop_jacobian(ctx, &geom.los, target.h, r_idx);
    }
    let v_scale = gneiss_core::variance::observation_variance(
        m.snr,
        geom.el,
        ctx.engine.config.tuning.snr_a,
        ctx.engine.config.tuning.snr_b,
    );
    target.r[(r_idx, r_idx)] = ctx.engine.config.tuning.dop_base_var * v_scale;
    target.types.push((m.eph.sat(), 3));
    *target.row += 1;
}

fn update_ekf(
    engine: &mut ProcessingEngine,
    z: &DVector<f64>,
    h: &DMatrix<f64>,
    r: &DMatrix<f64>,
    types: &[(gneiss_core::sat::SatelliteId, u8)],
) -> bool {
    let state = engine.current_state.as_mut().unwrap();
    let res = crate::engine::updater::update::<crate::engine::updater_math::TightCoupling>(
        state,
        z,
        h,
        r,
        engine.config.spp_consistency_threshold_m,
        Some(types),
        &engine.config.tuning,
    );
    if let Ok((valid_indices, _)) = res {
        valid_indices.len() < MIN_VALID_MEASUREMENTS
    } else {
        true
    }
}

fn handle_rejection(
    engine: &mut ProcessingEngine,
    rover_obs: &EpochObs,
    rejected: bool,
) -> Result<(), EngineError> {
    let state = engine.current_state.as_mut().unwrap();
    if rejected {
        state.consecutive_rejections += 1;
        if state.consecutive_rejections > MAX_CONSECUTIVE_REJECTIONS {
            reset_state_from_spp(engine, rover_obs)?;
        }
    } else {
        state.consecutive_rejections = 0;
    }
    Ok(())
}

fn reset_state_from_spp(
    engine: &mut ProcessingEngine,
    rover_obs: &EpochObs,
) -> Result<(), EngineError> {
    let spp_res = crate::spp::compute_spp(
        rover_obs,
        &engine.ephemerides,
        engine.klobuchar_params.as_ref(),
        &SppConfig::default(),
        None,
    )
    .map_err(|_| EngineError::NoObservations)?;
    let state = engine.current_state.as_mut().unwrap();
    state.position = spp_res.position;
    state.velocity = Vector3::zeros();
    state.rcv_clk_bias = spp_res.cdt;
    state.rcv_clk_drift = 0.0;
    state.decouple_position();
    for i in 0..3 {
        state.covariance[(i, i)] = COV_POS_RESET;
    }
    for i in 3..6 {
        state.covariance[(i, i)] = COV_VEL_RESET;
    }
    if state.covariance.nrows() > 15 {
        state.covariance[(15, 15)] = COV_CLK_RESET;
        state.covariance[(19, 19)] = COV_DRIFT_RESET;
    }
    state.is_reset = true;
    state.consecutive_rejections = 0;
    Ok(())
}

fn finalize_epoch(engine: &mut ProcessingEngine, rover_obs: &EpochObs) {
    engine.attempt_kinematic_alignment();
    let final_state = engine.current_state.as_ref().unwrap().clone();
    engine.state_history.push(final_state);
    engine.obs_history.push((rover_obs.clone(), None));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{EngineConfig, EngineMode};
    use gneiss_core::time::GpsTime;

    #[test]
    fn test_process_spp_tightly_coupled_startup() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::SppIns;
        let obs = EpochObs {
            time: GpsTime::new(2000, 1000.0),
            satellites: vec![],
        };
        assert!(process_spp_tightly_coupled(&mut engine, &obs).is_err());
    }
}
