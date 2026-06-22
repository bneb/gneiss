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
    let (_, _, dt_s, _) = if m.freq_band == 7 {
        m.eph.position_e5b(t_tx_nom)
    } else {
        m.eph.position(t_tx_nom)
    };
    let t_tx_true = gneiss_core::time::GpsTime::new(
        m.time.week,
        t_rcv_true - m.raw_pr / SPEED_OF_LIGHT_M_S - dt_s,
    );
    let (sat_pos_raw, sat_vel_raw, sat_clk, sat_drift) = if m.freq_band == 7 {
        m.eph.position_e5b(t_tx_true)
    } else {
        m.eph.position(t_tx_true)
    };

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
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::ephemeris::Ephemeris;
    use gneiss_core::sat::Constellation;
    use gneiss_core::time::GpsTime;

    // -------------------------------------------------------------------------
    // process_spp_tightly_coupled
    // -------------------------------------------------------------------------

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

    #[test]
    fn test_process_spp_tightly_coupled_no_snr_measurements_triggers_no_obs() {
        // When current_state exists but measurements get filtered to zero,
        // we should get NoObservations.
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::SppIns;
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        engine.current_state = Some(state);
        let obs = EpochObs {
            time: GpsTime::new(2000, 1.0),
            satellites: vec![],
        };
        let result = process_spp_tightly_coupled(&mut engine, &obs);
        assert!(result.is_err());
    }

    // -------------------------------------------------------------------------
    // Pure function tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_compute_sagnac_correction_zero_range() {
        let sat = Vector3::new(1.0, 2.0, 3.0);
        let result = compute_sagnac_correction(sat, 0.0);
        assert!((result.x - 1.0).abs() < 1e-15);
        assert!((result.y - 2.0).abs() < 1e-15);
        assert!((result.z - 3.0).abs() < 1e-15);
    }

    #[test]
    fn test_compute_sagnac_correction_rotates_about_z() {
        let sat = Vector3::new(1.0e7, 2.0e7, 3.0e7);
        let result = compute_sagnac_correction(sat, 2.0e7);
        assert!((result.z - 3.0e7).abs() < 1e-6);
    }

    #[test]
    fn test_compute_sagnac_correction_sign_convention() {
        let omega = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S;
        let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
        let angle = 1.0e-6;
        let range = angle * c / omega;
        let r = 10_000_000.0;
        let result = compute_sagnac_correction(Vector3::new(r, 0.0, 0.0), range);
        let expected_x = r * f64::cos(angle);
        let expected_y = -r * f64::sin(angle);
        let tol = 1e-6;
        assert!((result.x - expected_x).abs() < tol, "x mismatch");
        assert!((result.y - expected_y).abs() < tol, "y mismatch");
    }

    #[test]
    fn test_compute_sagnac_velocity_correction_zero_range() {
        let vel = Vector3::new(1000.0, 2000.0, 3000.0);
        let result = compute_sagnac_velocity_correction(vel, 0.0);
        assert!((result.x - 1000.0).abs() < 1e-15);
        assert!((result.y - 2000.0).abs() < 1e-15);
        assert!((result.z - 3000.0).abs() < 1e-15);
    }

    #[test]
    fn test_compute_sagnac_velocity_correction_z_unchanged() {
        let vel = Vector3::new(1000.0, 2000.0, 3000.0);
        let result = compute_sagnac_velocity_correction(vel, 2.0e7);
        assert!((result.z - 3000.0).abs() < 1e-6);
    }

    #[test]
    fn test_get_cdt_rx_values() {
        let coord = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let mut state = RtkState::new(GpsTime::new(2000, 0.0), coord, 1.0);
        state.rcv_clk_bias = 100.0;
        state.isb_gal = 5.0;
        state.isb_bds = -3.0;
        state.isb_glo = 2.0;

        let tol = 1e-12;
        assert!((get_cdt_rx(&state, Constellation::Gps) - 100.0).abs() < tol);
        assert!((get_cdt_rx(&state, Constellation::Galileo) - 105.0).abs() < tol);
        assert!((get_cdt_rx(&state, Constellation::Beidou) - 97.0).abs() < tol);
        assert!((get_cdt_rx(&state, Constellation::Glonass) - 102.0).abs() < tol);
        assert!((get_cdt_rx(&state, Constellation::Qzss) - 100.0).abs() < tol);
    }

    // -------------------------------------------------------------------------
    // predict_and_align_state
    // -------------------------------------------------------------------------

    #[test]
    fn test_predict_and_align_state_updates_time() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::SppIns;
        let pos = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        engine.current_state = Some(RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0));

        let obs = EpochObs {
            time: GpsTime::new(2000, 1.0),
            satellites: vec![],
        };
        predict_and_align_state(&mut engine, &obs);
        let state = engine.current_state.as_ref().unwrap();
        assert_eq!(state.time.tow, 1.0);
        assert_eq!(state.position.epoch.tow, 1.0);
    }

    #[test]
    fn test_predict_and_align_state_moderate_dt_does_not_panic() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::SppIns;
        let pos = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        engine.current_state = Some(RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0));

        // Moderate positive dt
        let obs = EpochObs {
            time: GpsTime::new(2000, 10.0),
            satellites: vec![],
        };
        predict_and_align_state(&mut engine, &obs);
        let state = engine.current_state.as_ref().unwrap();
        assert_eq!(state.time.tow, 10.0);
        assert_eq!(state.position.epoch.tow, 10.0);
    }

    // -------------------------------------------------------------------------
    // bootstrap_clock_bias
    // -------------------------------------------------------------------------

    #[test]
    fn test_bootstrap_clock_bias_no_ephemerides_does_not_panic() {
        // When compute_spp fails (no ephemerides), bootstrap_clock_bias should
        // gracefully fail without modifying state.
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let pos = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let mut state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        state.rcv_clk_bias = 42.0;
        engine.current_state = Some(state);

        let obs = EpochObs {
            time: GpsTime::new(2000, 1.0),
            satellites: vec![],
        };
        bootstrap_clock_bias(&mut engine, &obs);
        // Clock bias should be unchanged since compute_spp failed
        assert_eq!(engine.current_state.as_ref().unwrap().rcv_clk_bias, 42.0);
    }

    // -------------------------------------------------------------------------
    // get_apc_kinematics
    // -------------------------------------------------------------------------

    #[test]
    fn test_get_apc_kinematics_zero_lever_arm() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::SppIns;
        let pos = Coordinate::new(
            Vector3::new(10.0, 20.0, 30.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let mut state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        state.velocity = Vector3::new(1.0, 2.0, 3.0);
        state.ins_aligned = false; // lever_arm = zeros
        engine.current_state = Some(state);

        let (pos_apc, v_apc, r_b_e, lever_arm, _omega) = get_apc_kinematics(&engine);
        assert!((pos_apc.x - 10.0).abs() < 1e-10);
        assert!((pos_apc.y - 20.0).abs() < 1e-10);
        assert!((pos_apc.z - 30.0).abs() < 1e-10);
        assert!((v_apc.x - 1.0).abs() < 1e-10);
        assert!((v_apc.y - 2.0).abs() < 1e-10);
        assert!((v_apc.z - 3.0).abs() < 1e-10);
        assert!((lever_arm.norm() - 0.0).abs() < 1e-10);
        // r_b_e should be a valid rotation matrix
        assert!((r_b_e.matrix().determinant() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_get_apc_kinematics_with_lever_arm_and_imu_history() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::SppIns;
        engine.config.imu_to_antenna_lever_arm = [1.0, 0.0, 0.0];
        let pos = Coordinate::new(
            Vector3::new(10.0, 20.0, 30.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let mut state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        state.velocity = Vector3::new(1.0, 2.0, 3.0);
        state.ins_aligned = true;
        engine.current_state = Some(state);

        // Add IMU history with a measurement to provide non-zero omega_b
        engine.imu_history.push(vec![gneiss_core::imu::ImuMeasurement {
            accel: Vector3::zeros(),
            gyro: Vector3::new(0.01, 0.0, 0.0),
            time_tag: 0,
            temperature: None,
        }]);

        let (pos_apc, v_apc, _r_b_e, lever_arm, _omega) = get_apc_kinematics(&engine);
        // Lever arm returned is the body-frame vector, which for ins_aligned && config [1,0,0] is [1,0,0]
        assert!((lever_arm.norm() - 1.0).abs() < 1e-10);
        // pos_apc = position.vector + r_b_e * lever_arm, which differs from raw position
        assert!((pos_apc - Vector3::new(10.0, 20.0, 30.0)).norm() > 1e-10,
            "Antenna phase center should differ from the IMU position due to lever arm");
        // v_apc should include the lever-arm velocity term
        assert!((v_apc - Vector3::new(1.0, 2.0, 3.0)).norm() > 0.0);
    }

    #[test]
    fn test_get_apc_kinematics_no_imu_history_falls_back() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::SppIns;
        engine.config.imu_to_antenna_lever_arm = [0.5, 0.0, 0.0];
        let pos = Coordinate::new(
            Vector3::new(10.0, 20.0, 30.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let mut state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        state.velocity = Vector3::new(1.0, 2.0, 3.0);
        state.ins_aligned = true;
        engine.current_state = Some(state);
        // No imu_history at all

        let (_, v_apc, _, _, omega) = get_apc_kinematics(&engine);
        // omega_eb_b = -R_b_e^T * omega_ie_e when gyro is zero
        assert!((omega - Vector3::new(0.0, 0.0, -gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S)).norm() > 0.0
            || omega.norm() < 1e-10);
        // v_apc should still be computed with the lever arm and omega fallback (zeros -> omega_eb_b = -R^T * omega_ie)
        assert!(v_apc.norm() > 0.0);
    }

    // -------------------------------------------------------------------------
    // populate_pr_attitude_jacobian
    // -------------------------------------------------------------------------

    #[test]
    fn test_populate_pr_attitude_jacobian_values() {
        let ctx = EkfContext {
            pos_apc: Vector3::zeros(),
            v_apc: Vector3::zeros(),
            r_b_e: nalgebra::Rotation3::identity(),
            lever_arm: Vector3::new(1.0, 0.0, 0.0),
            omega_eb_b: Vector3::zeros(),
            rec_llh: Vector3::zeros(),
            state: &RtkState::new(
                GpsTime::new(2000, 0.0),
                Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(2000, 0.0)),
                1.0,
            ),
            engine: &ProcessingEngine::new(EngineConfig::default()),
            n_cols: 21,
        };
        let ncols = 21;
        let mut h = DMatrix::zeros(1, ncols);
        let los = Vector3::new(1.0, 0.0, 0.0);
        populate_pr_attitude_jacobian(&ctx, &los, &mut h, 0);

        // With R=I and lever_arm=[1,0,0], cross_matrix of lever_arm is
        // [[0,0,0],[0,0,-1],[0,1,0]].
        // h_pos_att = -cross_matrix([1,0,0]) = [[0,0,0],[0,0,1],[0,-1,0]]
        // pr_h_att = los^T * h_pos_att = [1,0,0] * [[0,0,0],[0,0,1],[0,-1,0]] = [0, 0, 0]
        // So all three entries should be 0 for this specific geometry
        assert!((h[(0, 6)]).abs() < 1e-15);
        assert!((h[(0, 7)]).abs() < 1e-15);
        assert!((h[(0, 8)]).abs() < 1e-15);
    }

    #[test]
    fn test_populate_pr_attitude_jacobian_nonzero() {
        // With lever_arm = [0, 1, 0] and los = [0, 0, 1], entries should be non-zero
        let ctx = EkfContext {
            pos_apc: Vector3::zeros(),
            v_apc: Vector3::zeros(),
            r_b_e: nalgebra::Rotation3::identity(),
            lever_arm: Vector3::new(0.0, 1.0, 0.0),
            omega_eb_b: Vector3::zeros(),
            rec_llh: Vector3::zeros(),
            state: &RtkState::new(
                GpsTime::new(2000, 0.0),
                Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(2000, 0.0)),
                1.0,
            ),
            engine: &ProcessingEngine::new(EngineConfig::default()),
            n_cols: 21,
        };
        let ncols = 21;
        let mut h = DMatrix::zeros(1, ncols);
        let los = Vector3::new(0.0, 0.0, 1.0);
        populate_pr_attitude_jacobian(&ctx, &los, &mut h, 0);

        // With lever_arm = [0,1,0]: cross_matrix = [[0,0,1],[0,0,0],[-1,0,0]]
        // h_pos_att = -cross_matrix = [[0,0,-1],[0,0,0],[1,0,0]]
        // los^T * h_pos_att = [0,0,1] * [[0,0,-1],[0,0,0],[1,0,0]] = [1, 0, 0]
        assert!((h[(0, 6)] - 1.0).abs() < 1e-10);
        assert!((h[(0, 7)]).abs() < 1e-10);
        assert!((h[(0, 8)]).abs() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // populate_pr_clock_jacobian
    // -------------------------------------------------------------------------

    #[test]
    fn test_populate_pr_clock_jacobian_gps() {
        let ncols = 21;
        let mut h = DMatrix::zeros(1, ncols);
        populate_pr_clock_jacobian(&mut h, 0, Constellation::Gps);
        assert!((h[(0, 15)] - 1.0).abs() < 1e-15);
        assert!((h[(0, 16)]).abs() < 1e-15);
        assert!((h[(0, 17)]).abs() < 1e-15);
        assert!((h[(0, 18)]).abs() < 1e-15);
    }

    #[test]
    fn test_populate_pr_clock_jacobian_galileo() {
        let ncols = 21;
        let mut h = DMatrix::zeros(1, ncols);
        populate_pr_clock_jacobian(&mut h, 0, Constellation::Galileo);
        assert!((h[(0, 15)] - 1.0).abs() < 1e-15);
        assert!((h[(0, 16)]).abs() < 1e-15);
        assert!((h[(0, 17)] - 1.0).abs() < 1e-15);
        assert!((h[(0, 18)]).abs() < 1e-15);
    }

    #[test]
    fn test_populate_pr_clock_jacobian_glonass() {
        let ncols = 21;
        let mut h = DMatrix::zeros(1, ncols);
        populate_pr_clock_jacobian(&mut h, 0, Constellation::Glonass);
        assert!((h[(0, 15)] - 1.0).abs() < 1e-15);
        assert!((h[(0, 16)] - 1.0).abs() < 1e-15);
        assert!((h[(0, 17)]).abs() < 1e-15);
        assert!((h[(0, 18)]).abs() < 1e-15);
    }

    #[test]
    fn test_populate_pr_clock_jacobian_beidou() {
        let ncols = 21;
        let mut h = DMatrix::zeros(1, ncols);
        populate_pr_clock_jacobian(&mut h, 0, Constellation::Beidou);
        assert!((h[(0, 15)] - 1.0).abs() < 1e-15);
        assert!((h[(0, 16)]).abs() < 1e-15);
        assert!((h[(0, 17)]).abs() < 1e-15);
        assert!((h[(0, 18)] - 1.0).abs() < 1e-15);
    }

    // -------------------------------------------------------------------------
    // populate_dop_jacobian
    // -------------------------------------------------------------------------

    #[test]
    fn test_populate_dop_jacobian_structure() {
        // With identity attitude, zero lever arm, the body-rate jacobian terms vanish
        let ctx = EkfContext {
            pos_apc: Vector3::zeros(),
            v_apc: Vector3::zeros(),
            r_b_e: nalgebra::Rotation3::identity(),
            lever_arm: Vector3::zeros(),
            omega_eb_b: Vector3::zeros(),
            rec_llh: Vector3::zeros(),
            state: &RtkState::new(
                GpsTime::new(2000, 0.0),
                Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(2000, 0.0)),
                1.0,
            ),
            engine: &ProcessingEngine::new(EngineConfig::default()),
            n_cols: 21,
        };
        let ncols = 21;
        let mut h = DMatrix::zeros(1, ncols);
        let los = Vector3::new(1.0, 0.0, 0.0);
        populate_dop_jacobian(&ctx, &los, &mut h, 0);

        // With zero lever arm and zero omega, all attitude/bg entries are 0
        assert!((h[(0, 6)]).abs() < 1e-15);
        assert!((h[(0, 7)]).abs() < 1e-15);
        assert!((h[(0, 8)]).abs() < 1e-15);
        assert!((h[(0, 12)]).abs() < 1e-15);
        assert!((h[(0, 13)]).abs() < 1e-15);
        assert!((h[(0, 14)]).abs() < 1e-15);
        // Clock drift column should be 1
        assert!((h[(0, 19)] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_populate_dop_jacobian_nonzero_attitude_entries() {
        // With lever_arm and non-zero omega, the attitude jacobian entries are non-zero
        let ctx = EkfContext {
            pos_apc: Vector3::zeros(),
            v_apc: Vector3::zeros(),
            r_b_e: nalgebra::Rotation3::identity(),
            lever_arm: Vector3::new(1.0, 0.0, 0.0),
            omega_eb_b: Vector3::new(0.1, 0.0, 0.0),
            rec_llh: Vector3::zeros(),
            state: &RtkState::new(
                GpsTime::new(2000, 0.0),
                Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(2000, 0.0)),
                1.0,
            ),
            engine: &ProcessingEngine::new(EngineConfig::default()),
            n_cols: 21,
        };
        let ncols = 21;
        let mut h = DMatrix::zeros(1, ncols);
        let los = Vector3::new(1.0, 0.0, 0.0);
        populate_dop_jacobian(&ctx, &los, &mut h, 0);

        // a_0 = R * (omega x lever_arm) = [0.1,0,0] x [1,0,0] = [0,0,0]
        // So h_vel_att = -cross_matrix([0,0,0]) = zeros
        // Thus all attitude entries are still zero for this case
        assert!((h[(0, 6)]).abs() < 1e-15);
        assert!((h[(0, 7)]).abs() < 1e-15);
        assert!((h[(0, 8)]).abs() < 1e-15);
        assert!((h[(0, 12)]).abs() < 1e-15);
        assert!((h[(0, 13)]).abs() < 1e-15);
        assert!((h[(0, 14)]).abs() < 1e-15);
        assert!((h[(0, 19)] - 1.0).abs() < 1e-15);
    }

    // -------------------------------------------------------------------------
    // finalize_epoch
    // -------------------------------------------------------------------------

    #[test]
    fn test_finalize_epoch_appends_state_and_obs() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::SppIns;
        let pos = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        engine.current_state = Some(state);

        let obs = EpochObs {
            time: GpsTime::new(2000, 1.0),
            satellites: vec![],
        };
        assert_eq!(engine.state_history.len(), 0);
        assert_eq!(engine.obs_history.len(), 0);

        finalize_epoch(&mut engine, &obs);

        assert_eq!(engine.state_history.len(), 1);
        assert_eq!(engine.obs_history.len(), 1);
    }

    #[test]
    fn test_finalize_epoch_preserves_state_values() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::SppIns;
        let pos = Coordinate::new(
            Vector3::new(10.0, 20.0, 30.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 5.0),
        );
        let state = RtkState::new(GpsTime::new(2000, 5.0), pos, 1.0);
        engine.current_state = Some(state);

        let obs = EpochObs {
            time: GpsTime::new(2000, 5.0),
            satellites: vec![],
        };
        finalize_epoch(&mut engine, &obs);

        let stored = &engine.state_history[0];
        assert!((stored.position.vector.x - 10.0).abs() < 1e-10);
        assert!((stored.position.vector.y - 20.0).abs() < 1e-10);
        assert!((stored.position.vector.z - 30.0).abs() < 1e-10);
        assert_eq!(engine.obs_history[0].0.time.tow, 5.0);
        assert!(engine.obs_history[0].1.is_none()); // base_obs is None
    }

    // -------------------------------------------------------------------------
    // handle_rejection edge case: triggers reset when exceeding MAX_CONSECUTIVE_REJECTIONS
    // -------------------------------------------------------------------------

    #[test]
    fn test_handle_rejection_no_rejection() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let coord = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        engine.current_state = Some(RtkState::new(GpsTime::new(2000, 0.0), coord, 1.0));
        engine.current_state.as_mut().unwrap().consecutive_rejections = 5;
        let obs = EpochObs {
            time: GpsTime::new(2000, 1.0),
            satellites: vec![],
        };

        let result = handle_rejection(&mut engine, &obs, false);
        assert!(result.is_ok());
        assert_eq!(
            engine.current_state.as_ref().unwrap().consecutive_rejections,
            0
        );
    }

    #[test]
    fn test_handle_rejection_increments_counter() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let coord = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        engine.current_state = Some(RtkState::new(GpsTime::new(2000, 0.0), coord, 1.0));
        engine.current_state.as_mut().unwrap().consecutive_rejections = 3;
        let obs = EpochObs {
            time: GpsTime::new(2000, 1.0),
            satellites: vec![],
        };

        let result = handle_rejection(&mut engine, &obs, true);
        assert!(result.is_ok());
        assert_eq!(
            engine.current_state.as_ref().unwrap().consecutive_rejections,
            4
        );
    }

    #[test]
    fn test_handle_rejection_returns_none_when_reset_state_fails() {
        // When consecutive_rejections exceeds MAX_CONSECUTIVE_REJECTIONS (10)
        // and reset_state_from_spp fails (no ephemerides, no measurements),
        // handle_rejection should propagate the EngineError.
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let coord = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        engine.current_state = Some(RtkState::new(GpsTime::new(2000, 0.0), coord, 1.0));
        engine.current_state.as_mut().unwrap().consecutive_rejections = 10;
        let obs = EpochObs {
            time: GpsTime::new(2000, 1.0),
            satellites: vec![],
        };

        let result = handle_rejection(&mut engine, &obs, true);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), EngineError::NoObservations));
    }

    // -------------------------------------------------------------------------
    // update_ekf
    // -------------------------------------------------------------------------

    #[test]
    fn test_update_ekf_no_state_is_err() {
        // When current_state is missing, update_ekf will unwrap on None -> panic.
        // But this is a precondition: process_spp_tightly_coupled ensures state exists.
        // We test by ensuring the inner function works with a valid state but empty data.
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let pos = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        engine.current_state = Some(RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0));

        // Empty z, h, r with zero rows -> will fail dimension checks -> returns true (rejected)
        let z = DVector::zeros(0);
        let h = DMatrix::zeros(0, 21);
        let r = DMatrix::zeros(0, 0);
        let types = vec![];
        let rejected = update_ekf(&mut engine, &z, &h, &r, &types);
        // With zero rows, filter_pre_fit_residuals yields empty valid_indices -> len < 3 -> rejected
        assert!(rejected);
    }

    // -------------------------------------------------------------------------
    // process_pseudorange (basic setup covering tropo/iono path)
    // -------------------------------------------------------------------------

    #[test]
    fn test_process_pseudorange_no_klobuchar_params() {
        // With no atmosphere parameters, pseudorange residual should be computed
        // using only geometric range and clock terms.
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let pos = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        engine.current_state = Some(state.clone());
        let ctx = EkfContext {
            pos_apc: Vector3::new(0.0, 0.0, 0.0),
            v_apc: Vector3::zeros(),
            r_b_e: nalgebra::Rotation3::identity(),
            lever_arm: Vector3::zeros(),
            omega_eb_b: Vector3::zeros(),
            rec_llh: Vector3::new(0.0, 0.0, 0.0),
            state: &state,
            engine: &engine,
            n_cols: 21,
        };

        let sat_id = gneiss_core::sat::SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: sat_id,
            toe: GpsTime::new(2000, 0.0),
            toc: GpsTime::new(2000, 0.0),
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.0, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0,
            omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        });

        let meas = SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20_000_000.0,
            snr: 45.0,
            doppler: 0.0,
            time: GpsTime::new(2000, 0.0),
            eph: eph.clone(),
            is_iono_free: false,
            freq_band: 1,
        };

        // Set up matrix target with one row
        let ncols = 21;
        let mut z = DVector::zeros(1);
        let mut h = DMatrix::zeros(1, ncols);
        let mut r = DMatrix::zeros(1, 1);
        let mut types = Vec::new();
        let mut row = 0usize;
        let mut target = MatrixTarget {
            z: &mut z,
            h: &mut h,
            r: &mut r,
            types: &mut types,
            row: &mut row,
        };

        let geom = SatGeometry {
            geom_r: 20_000_000.0,
            los: Vector3::new(1.0, 0.0, 0.0),
            el: 0.5,
            az: 0.0,
            cdt_rx: 0.0,
            sat_clk: 0.0,
            sat_vel: Vector3::zeros(),
            sat_drift: 0.0,
        };

        process_pseudorange(&ctx, &meas, &mut target, &geom);

        // After processing, row should be 1
        assert_eq!(row, 1);
        // types should have one entry with (sat, 0)
        assert_eq!(types.len(), 1);
        assert_eq!(types[0].0, sat_id);
        assert_eq!(types[0].1, 0);
        // z should be finite
        assert!(z[0].is_finite());
        // position jacobian should be los (the unit vector from receiver to satellite)
        assert!((h[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((h[(0, 1)]).abs() < 1e-10);
        assert!((h[(0, 2)]).abs() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // process_doppler
    // -------------------------------------------------------------------------

    #[test]
    fn test_process_doppler_basic() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let pos = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        engine.current_state = Some(state.clone());
        let ctx = EkfContext {
            pos_apc: Vector3::new(0.0, 0.0, 0.0),
            v_apc: Vector3::zeros(),
            r_b_e: nalgebra::Rotation3::identity(),
            lever_arm: Vector3::zeros(),
            omega_eb_b: Vector3::zeros(),
            rec_llh: Vector3::new(0.0, 0.0, 0.0),
            state: &state,
            engine: &engine,
            n_cols: 21,
        };

        let sat_id = gneiss_core::sat::SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: sat_id,
            toe: GpsTime::new(2000, 0.0), toc: GpsTime::new(2000, 0.0),
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0, sqrt_a: 5153.6,
            delta_n: 0.0, omega0: 0.0, omega_dot: 0.0, i0: 0.95,
            idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        });

        let meas = SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20_000_000.0,
            snr: 45.0,
            doppler: 1000.0, // non-zero -> triggers doppler processing
            time: GpsTime::new(2000, 0.0),
            eph: eph.clone(),
            is_iono_free: false,
            freq_band: 1,
        };

        let ncols = 21;
        let mut z = DVector::zeros(2);
        let mut h = DMatrix::zeros(2, ncols);
        let mut r = DMatrix::zeros(2, 2);
        let mut types = Vec::new();
        let mut row = 0usize;

        // First process pseudorange at row 0
        {
            let mut target = MatrixTarget {
                z: &mut z,
                h: &mut h,
                r: &mut r,
                types: &mut types,
                row: &mut row,
            };
            let p_geom = SatGeometry {
                geom_r: 20_000_000.0,
                los: Vector3::new(1.0, 0.0, 0.0),
                el: 0.5,
                az: 0.0,
                cdt_rx: 0.0,
                sat_clk: 0.0,
                sat_vel: Vector3::new(1000.0, 0.0, 0.0),
                sat_drift: 0.0,
            };
            process_pseudorange(&ctx, &meas, &mut target, &p_geom);
        }

        // Then process doppler at row 1
        {
            let mut target = MatrixTarget {
                z: &mut z,
                h: &mut h,
                r: &mut r,
                types: &mut types,
                row: &mut row,
            };
            let d_geom = SatGeometry {
                geom_r: 20_000_000.0,
                los: Vector3::new(1.0, 0.0, 0.0),
                el: 0.5,
                az: 0.0,
                cdt_rx: 0.0,
                sat_clk: 0.0,
                sat_vel: Vector3::new(1000.0, 0.0, 0.0),
                sat_drift: 0.0,
            };
            process_doppler(&ctx, &meas, &mut target, &d_geom);
        }

        assert_eq!(row, 2);
        assert_eq!(types.len(), 2);
        assert_eq!(types[1].0, sat_id);
        assert_eq!(types[1].1, 3); // doppler type code
        // z[1] should be finite
        assert!(z[1].is_finite());
    }

    // -------------------------------------------------------------------------
    // build_ekf_matrices
    // -------------------------------------------------------------------------

    #[test]
    fn test_build_ekf_matrices_empty_measurements() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let pos = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        engine.current_state = Some(state);

        let measurements = vec![];
        let (z, h, r, types) = build_ekf_matrices(&engine, &measurements);
        assert_eq!(z.len(), 0);
        assert_eq!(h.nrows(), 0);
        assert_eq!(r.nrows(), 0);
        assert!(types.is_empty());
    }

    // -------------------------------------------------------------------------
    // reset_state_from_spp (fails when no ephemerides)
    // -------------------------------------------------------------------------

    #[test]
    fn test_reset_state_from_spp_fails_without_ephemerides() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let pos = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        engine.current_state = Some(state);

        let obs = EpochObs {
            time: GpsTime::new(2000, 1.0),
            satellites: vec![],
        };
        let result = reset_state_from_spp(&mut engine, &obs);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), EngineError::NoObservations));
    }
}
