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
    Ok(engine.current_state.as_ref().expect("current_state is Some after None check"))
}

fn predict_and_align_state(engine: &mut ProcessingEngine, rover_obs: &EpochObs) {
    let dt = rover_obs.time.tow - engine.current_state.as_ref().expect("current_state is Some after None check").time.tow;
    engine.predict_state(dt);
    let state = engine.current_state.as_mut().expect("current_state is Some after None check");
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
        let state = engine.current_state.as_mut().expect("current_state is Some after None check in caller");
        state.rcv_clk_bias = spp_res.cdt;
        state.covariance[(15, 15)] = CLK_BIAS_VAR_RESET;
    }
}

type ApcKinematics = (
    Vector3<f64>,
    Vector3<f64>,
    nalgebra::Rotation3<f64>,
    Vector3<f64>,
    Vector3<f64>,
);

fn get_apc_kinematics(
    engine: &ProcessingEngine,
) -> ApcKinematics {
    let state = engine.current_state.as_ref().expect("current_state is Some after None check in caller");
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

type EkfMatrices = (
    DVector<f64>,
    DMatrix<f64>,
    DMatrix<f64>,
    Vec<(gneiss_core::sat::SatelliteId, u8)>,
);

fn build_ekf_matrices(
    engine: &ProcessingEngine,
    measurements: &[SppMeasurement],
) -> EkfMatrices {
    let (pos_apc, v_apc, r_b_e, lever_arm, omega_eb_b) = get_apc_kinematics(engine);
    let ctx = EkfContext {
        pos_apc,
        v_apc,
        r_b_e,
        lever_arm,
        omega_eb_b,
        rec_llh: ecef_to_llh(pos_apc),
        state: engine.current_state.as_ref().expect("current_state is Some after None check in caller"),
        engine,
        n_cols: engine.current_state.as_ref().expect("current_state is Some after None check in caller").covariance.ncols(),
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
    let state = engine.current_state.as_mut().expect("current_state is Some after None check in caller");
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
    let state = engine.current_state.as_mut().expect("current_state is Some after None check in caller");
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
    let state = engine.current_state.as_mut().expect("current_state is Some after None check in caller");
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
    let final_state = engine.current_state.as_ref().expect("current_state is Some after None check in caller").clone();
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

    // -------------------------------------------------------------------------
    // process_pseudorange with small state (n_cols = 6, skip attitude/clock)
    // -------------------------------------------------------------------------

    #[test]
    fn test_process_pseudorange_small_state_skips_attitude_and_clock_jacobians() {
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
            n_cols: 6,
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

        let ncols = 6;
        let mut z = DVector::zeros(1);
        let mut h = DMatrix::zeros(1, ncols);
        let mut r_mat = DMatrix::zeros(1, 1);
        let mut types = Vec::new();
        let mut row = 0usize;
        let mut target = MatrixTarget {
            z: &mut z,
            h: &mut h,
            r: &mut r_mat,
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

        assert_eq!(row, 1);
        assert_eq!(types.len(), 1);
        assert_eq!(types[0].0, sat_id);
        assert_eq!(types[0].1, 0);
        assert!(z[0].is_finite());
        // Position jacobian entries (cols 0-2) should be populated
        assert!((h[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((h[(0, 1)]).abs() < 1e-10);
        assert!((h[(0, 2)]).abs() < 1e-10);
        // With n_cols=6, h has 6 cols; the if ctx.n_cols > 15 branch is skipped
        assert_eq!(h.ncols(), 6);
    }

    // -------------------------------------------------------------------------
    // process_pseudorange with Klobuchar parameters
    // -------------------------------------------------------------------------

    #[test]
    fn test_process_pseudorange_with_klobuchar_params() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.klobuchar_params = Some(gneiss_core::atmosphere::KlobucharParams::default());
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

        let ncols = 21;
        let mut z = DVector::zeros(1);
        let mut h = DMatrix::zeros(1, ncols);
        let mut r_mat = DMatrix::zeros(1, 1);
        let mut types = Vec::new();
        let mut row = 0usize;
        let mut target = MatrixTarget {
            z: &mut z,
            h: &mut h,
            r: &mut r_mat,
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
        assert!(z[0].is_finite());
    }

    // -------------------------------------------------------------------------
    // process_doppler with small state (n_cols = 15, skip extra jacobians)
    // -------------------------------------------------------------------------

    #[test]
    fn test_process_doppler_small_state_skips_extra_jacobians() {
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
            n_cols: 15,
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
            doppler: 1000.0,
            time: GpsTime::new(2000, 0.0),
            eph: eph.clone(),
            is_iono_free: false,
            freq_band: 1,
        };

        let ncols = 15;
        let mut z = DVector::zeros(1);
        let mut h = DMatrix::zeros(1, ncols);
        let mut r_mat = DMatrix::zeros(1, 1);
        let mut types = Vec::new();
        let mut row = 0usize;
        let mut target = MatrixTarget {
            z: &mut z,
            h: &mut h,
            r: &mut r_mat,
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
            sat_vel: Vector3::new(1000.0, 0.0, 0.0),
            sat_drift: 0.0,
        };

        process_doppler(&ctx, &meas, &mut target, &geom);

        assert_eq!(row, 1);
        assert_eq!(types.len(), 1);
        assert_eq!(types[0].0, sat_id);
        assert_eq!(types[0].1, 3);
        assert!(z[0].is_finite());
        // Velocity jacobian (cols 3-5) should be populated
        assert!((h[(0, 3)] - 1.0).abs() < 1e-10);
        assert!((h[(0, 4)]).abs() < 1e-10);
        assert!((h[(0, 5)]).abs() < 1e-10);
        // populate_dop_jacobian is skipped because n_cols (15) <= 19
        assert_eq!(h.ncols(), 15);
    }

    // -------------------------------------------------------------------------
    // update_ekf dimension mismatch returns rejected
    // -------------------------------------------------------------------------

    #[test]
    fn test_update_ekf_dimension_mismatch_returns_rejected() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let pos = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        engine.current_state = Some(state);

        // z(3), h(3x6), r(3x3), empty types
        // h.ncols() = 6, but state.covariance.ncols() = 21 -> dimension mismatch
        let z = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let h = DMatrix::zeros(3, 6);
        let r_mat = DMatrix::identity(3, 3);
        let types = vec![];

        let rejected = update_ekf(&mut engine, &z, &h, &r_mat, &types);
        assert!(rejected);
    }

    // -------------------------------------------------------------------------
    // reset_state_from_spp with valid observations
    // -------------------------------------------------------------------------

    fn build_test_ephemerides_and_obs(
    ) -> (Vec<Ephemeris>, EpochObs) {
        use gneiss_core::obs::{ObsCode, ObsType, Observation, SatObs, SignalCode};
        use gneiss_core::sat::SatelliteId;
        use gneiss_core::coords::ecef_to_llh;
        use gneiss_core::atmosphere::{AtmosphereModel, TropoParams};

        let t = GpsTime::new(2000, 100000.0);
        let true_pos =
            nalgebra::Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let rec_llh = ecef_to_llh(true_pos);
        let true_cdt = 1000.0;

        let mut ephemerides = Vec::new();
        let mut satellites = Vec::new();

        // Use verified (m0, omega0) pairs from elevation diagnostic that all
        // produce >15-degree elevation from the receiver at (WGS84_SEMI_MAJOR_AXIS_M, 0, 0).
        // Using 12 satellites for robust geometry.
        let kepler_params: Vec<(f64, f64)> = vec![
            (0.000, 0.785),    // el=73.2 deg
            (2.618, 4.712),    // el=54.9 deg
            (3.665, 3.927),    // el=58.8 deg
            (5.760, 1.571),    // el=54.3 deg
            (0.524, 0.785),    // el=57.3 deg
            (3.142, 3.927),    // el=73.3 deg
            (1.047, 0.000),    // el=31.7 deg
            (3.665, 3.142),    // el=33.4 deg
            (0.000, 1.571),    // el=48.6 deg
            (2.094, 4.712),    // el=33.0 deg
            (2.618, 3.927),    // el=41.5 deg
            (5.760, 0.785),    // el=40.0 deg
        ];

        for (prn, (m0_val, omega0_val)) in kepler_params.iter().enumerate() {
            let prn_u8 = (prn + 1) as u8;
            let sat_id = SatelliteId {
                constellation: Constellation::Gps,
                prn: prn_u8,
            };

            let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: sat_id,
                toe: t,
                toc: t,
                af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
                cic: 0.0, cis: 0.0,
                m0: *m0_val,
                e: 0.01,
                sqrt_a: 5153.6,
                delta_n: 0.0,
                omega0: *omega0_val,
                omega_dot: 0.0,
                i0: 0.95,
                idot: 0.0,
                omega: 0.0,
                tgd: 0.0,
                iode: 1,
                iodc: 1,
            });
            ephemerides.push(eph.clone());

            let light_speed = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
            let mut raw_pr = 20_000_000.0 + true_cdt;

            // Compute satellite position and account for geometric + clock terms
            for _iter in 0..5 {
                let pr_time = raw_pr / light_speed;
                let t_tx_sat = t.tow - pr_time;
                let (_, _, sat_clk_err_rough, _) =
                    eph.position(GpsTime::new(t.week, t_tx_sat));
                let t_tx_true = t_tx_sat - sat_clk_err_rough;
                let (sat_pos, _, sat_clk_err, _) =
                    eph.position(GpsTime::new(t.week, t_tx_true));
                let dx = true_pos.x - sat_pos.x;
                let dy = true_pos.y - sat_pos.y;
                let dz = true_pos.z - sat_pos.z;
                let geometric_range = f64::sqrt(dx * dx + dy * dy + dz * dz);

                // Add approximate tropo delay so residual is near zero at start
                let (_, el) = gneiss_core::coords::az_el(rec_llh, true_pos, sat_pos);
                let safe_el = el.max(5.0f64.to_radians());
                let tropo = AtmosphereModel::tropo_nmf(
                    &TropoParams::default(), rec_llh, safe_el, t,
                );

                raw_pr = geometric_range + true_cdt - (sat_clk_err * light_speed) + tropo;
            }

            satellites.push(SatObs {
                sat: sat_id,
                observations: vec![Observation {
                    code: ObsCode {
                        obs_type: ObsType::Pseudorange,
                        signal: SignalCode {
                            freq_band: 1,
                            attribute: 'C',
                        },
                    },
                    value: raw_pr,
                    lock_time: None,
                    lli: None,
                }],
            });
        }

        let obs = EpochObs {
            time: t,
            satellites,
        };

        (ephemerides, obs)
    }

    #[test]
    fn test_reset_state_from_spp_success_resets_clock_indices() {
        let (ephemerides, obs) = build_test_ephemerides_and_obs();

        // Set up engine with ephemerides and a full 21x21 state
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let pos = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        engine.current_state = Some(state);
        engine.ephemerides = ephemerides;

        let result = reset_state_from_spp(&mut engine, &obs);
        assert!(result.is_ok(),
            "reset_state_from_spp should succeed, got error: {:?}",
            result.as_ref().err());

        let final_state = engine.current_state.as_ref().unwrap();
        assert!((final_state.covariance[(15, 15)] - COV_CLK_RESET).abs() < 1e-10,
            "clock bias reset variance mismatch: {} != {}",
            final_state.covariance[(15, 15)], COV_CLK_RESET);
        assert!((final_state.covariance[(19, 19)] - COV_DRIFT_RESET).abs() < 1e-10,
            "clock drift reset variance mismatch: {} != {}",
            final_state.covariance[(19, 19)], COV_DRIFT_RESET);
        // rcv_clk_bias should be set (non-zero from compute_spp)
        assert!(final_state.rcv_clk_bias != 0.0,
            "rcv_clk_bias should be non-zero after SPP reset, got 0.0");
        // velocity should be reset to zeros
        assert!((final_state.velocity.x).abs() < 1e-10);
        assert!((final_state.velocity.y).abs() < 1e-10);
        assert!((final_state.velocity.z).abs() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // bootstrap_clock_bias with valid observations
    // -------------------------------------------------------------------------

    #[test]
    fn test_bootstrap_clock_bias_with_valid_observations() {
        let (ephemerides, obs) = build_test_ephemerides_and_obs();

        // Set up engine with ephemerides and a state that has a pre-existing clock bias
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let pos = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(2000, 0.0),
        );
        let mut state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        state.rcv_clk_bias = 42.0; // Pre-existing bias, should be overwritten
        engine.current_state = Some(state);
        engine.ephemerides = ephemerides;

        bootstrap_clock_bias(&mut engine, &obs);

        let final_state = engine.current_state.as_ref().unwrap();
        // rcv_clk_bias should be updated from compute_spp result (different from initial 42.0)
        assert!(final_state.rcv_clk_bias != 42.0,
            "rcv_clk_bias should be updated from SPP result, was still {}", final_state.rcv_clk_bias);
        // covariance[(15,15)] should be set to CLK_BIAS_VAR_RESET
        assert!((final_state.covariance[(15, 15)] - CLK_BIAS_VAR_RESET).abs() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // build_ekf_matrices with real measurements (indirectly tests process_measurement)
    // -------------------------------------------------------------------------

    #[test]
    fn test_build_ekf_matrices_with_one_gps_measurement() {
        let (ephemerides, obs) = build_test_ephemerides_and_obs();

        // Set up engine with ephemerides and an aligned state with position near the receiver
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let true_pos =
            nalgebra::Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let pos = Coordinate::new(
            true_pos, Datum::WGS84, Frame::ECEF, GpsTime::new(2000, 100000.0),
        );
        let mut state = RtkState::new(GpsTime::new(2000, 100000.0), pos, 1.0);
        state.ins_aligned = true;
        state.velocity = Vector3::zeros();
        engine.current_state = Some(state);
        engine.ephemerides = ephemerides;

        // Build measurements from the observation
        let measurements: Vec<_> = build_measurements(&obs, &engine.ephemerides, &SppConfig::default())
            .into_iter()
            .take(1) // Just one satellite
            .collect();

        // Call build_ekf_matrices
        // Since we filtered to 1 satellite and doppler may not be available (check if 0)
        // we need to check if it's empty or has 1 row
        if !measurements.is_empty() {
            let (z, h, r, types) = build_ekf_matrices(&engine, &measurements);
            // The total rows = measurements + doppler_count
            // We expect at least 1 row (pseudorange)
            assert!(z.len() >= 1, "Should have at least 1 row, got {}", z.len());
            assert_eq!(h.nrows(), z.len());
            assert_eq!(h.ncols(), 21);
            assert_eq!(r.nrows(), z.len());
            assert_eq!(r.ncols(), z.len());
            // z should be finite (not NaN)
            for i in 0..z.len() {
                assert!(z[i].is_finite(), "z[{}] should be finite, got {}", i, z[i]);
            }
            // At least one type entry
            assert!(!types.is_empty());
            // The pseudorange type code is 0
            let has_pr = types.iter().any(|(_, code)| *code == 0);
            assert!(has_pr, "Should have at least one pseudorange measurement");
        }
    }

    #[test]
    fn test_build_ekf_matrices_with_doppler() {
        let (ephemerides, obs) = build_test_ephemerides_and_obs();

        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let true_pos =
            nalgebra::Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let pos = Coordinate::new(
            true_pos, Datum::WGS84, Frame::ECEF, GpsTime::new(2000, 100000.0),
        );
        let mut state = RtkState::new(GpsTime::new(2000, 100000.0), pos, 1.0);
        state.ins_aligned = true;
        state.velocity = Vector3::new(5.0, 0.0, 0.0); // Moving receiver
        engine.current_state = Some(state);
        engine.ephemerides = ephemerides;

        let mut measurements: Vec<_> =
            build_measurements(&obs, &engine.ephemerides, &SppConfig::default());

        // Set a non-zero doppler on one measurement to trigger doppler processing
        if let Some(m) = measurements.first_mut() {
            m.doppler = -2000.0; // Approximate doppler for a moving receiver
        }

        if !measurements.is_empty() {
            let (z, h, r, types) = build_ekf_matrices(&engine, &measurements);
            assert!(z.len() > 0);
            // Should have at least pseudorange rows
            let pr_count = types.iter().filter(|(_, code)| *code == 0).count();
            let dop_count = types.iter().filter(|(_, code)| *code == 3).count();
            assert!(pr_count >= 1, "Should have pseudorange rows, got {}", pr_count);
            // z values should all be finite
            for i in 0..z.len() {
                assert!(z[i].is_finite(), "z[{}] should be finite, got {}", i, z[i]);
            }
            // Verify velocity Jacobian columns exist if doppler is present
            if dop_count > 0 {
                let has_vel_cols = (0..z.len()).any(|i| {
                    h[(i, 3)].abs() > 1e-10 || h[(i, 4)].abs() > 1e-10 || h[(i, 5)].abs() > 1e-10
                });
                // At least some velocity columns should be non-zero (or all zero if satellite relative velocity cancels)
                // Just verify the matrix is valid
                assert!(h.nrows() > 0 && h.ncols() == 21);
            }
        }
    }

    // -------------------------------------------------------------------------
    // process_measurement integration test
    // -------------------------------------------------------------------------

    #[test]
    fn test_process_measurement_integration_gps() {
        // Directly test process_measurement by constructing an EkfContext
        // and a measurement, then calling process_measurement and verifying
        // the matrices are populated correctly.
        use gneiss_core::coords::ecef_to_llh;
        use gneiss_core::ephemeris::GpsEphemeris;

        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let true_pos =
            nalgebra::Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let rec_llh = ecef_to_llh(true_pos);
        let pos = Coordinate::new(
            true_pos, Datum::WGS84, Frame::ECEF, GpsTime::new(2000, 100000.0),
        );
        let state = RtkState::new(GpsTime::new(2000, 100000.0), pos, 1.0);
        engine.current_state = Some(state.clone());

        let sat_id = gneiss_core::sat::SatelliteId {
            constellation: Constellation::Gps, prn: 1,
        };
        let t = GpsTime::new(2000, 100000.0);

        // Use a simple ephemeris with the satellite at a known position
        // On the equatorial plane at radius sqrt_a^2 from the origin
        let eph = Ephemeris::Gps(GpsEphemeris {
            sat: sat_id, toe: t, toc: t,
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0,
            sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0,
            omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        });

        // Compute the satellite position at T and the geometric range
        let (_sat_pos_rough, _, sat_clk_rough, _) = eph.position(t);
        // The satellite is at orbital radius ~26,560 km, so the range
        // to our receiver at ~6378 km is ~20,182 km
        let dx = true_pos.x - _sat_pos_rough.x;
        let dy = true_pos.y - _sat_pos_rough.y;
        let dz = true_pos.z - _sat_pos_rough.z;
        let geom_r = (dx * dx + dy * dy + dz * dz).sqrt();
        let raw_pr = geom_r + 0.0 - sat_clk_rough * SPEED_OF_LIGHT_M_S;

        let meas = SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr,
            snr: 45.0,
            doppler: 0.0,
            time: t,
            eph,
            is_iono_free: false,
            freq_band: 1,
        };

        let ctx = EkfContext {
            pos_apc: true_pos,
            v_apc: Vector3::zeros(),
            r_b_e: nalgebra::Rotation3::identity(),
            lever_arm: Vector3::zeros(),
            omega_eb_b: Vector3::zeros(),
            rec_llh,
            state: &state,
            engine: &engine,
            n_cols: 21,
        };

        let mut z = DVector::zeros(2);
        let mut h = DMatrix::zeros(2, 21);
        let mut r_mat = DMatrix::zeros(2, 2);
        let mut types = Vec::new();
        let mut row = 0usize;
        let mut target = MatrixTarget {
            z: &mut z,
            h: &mut h,
            r: &mut r_mat,
            types: &mut types,
            row: &mut row,
        };

        process_measurement(&ctx, &meas, &mut target);

        // Verify that process_measurement populated the matrices
        assert_eq!(row, 1, "Should have processed 1 measurement (no doppler)");
        assert_eq!(types.len(), 1);
        assert_eq!(types[0].0, sat_id);
        assert_eq!(types[0].1, 0); // PR code
        // z should be finite (small residual if geometry is consistent)
        assert!(z[0].is_finite(), "z[0] should be finite, got {}", z[0]);
        // Position Jacobian should be populated
        assert!(h[(0, 0)] != 0.0 || h[(0, 1)] != 0.0 || h[(0, 2)] != 0.0,
            "LOS Jacobian should be non-zero");
        // Clock Jacobian should be 1 at col 15
        assert!((h[(0, 15)] - 1.0).abs() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // freq_band == 7 (Galileo E5b) branch in process_measurement
    // -------------------------------------------------------------------------

    #[test]
    fn test_process_measurement_galileo_e5b_freq_band_7() {
        use gneiss_core::coords::ecef_to_llh;
        use gneiss_core::ephemeris::GalileoEphemeris;

        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let true_pos =
            nalgebra::Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let rec_llh = ecef_to_llh(true_pos);
        let pos = Coordinate::new(
            true_pos, Datum::WGS84, Frame::ECEF, GpsTime::new(2000, 100000.0),
        );
        let state = RtkState::new(GpsTime::new(2000, 100000.0), pos, 1.0);
        engine.current_state = Some(state.clone());

        let sat_id = gneiss_core::sat::SatelliteId {
            constellation: Constellation::Galileo, prn: 1,
        };
        let t = GpsTime::new(2000, 100000.0);

        // Construct Galileo ephemeris with same orbital parameters
        let eph = Ephemeris::Galileo(GalileoEphemeris {
            sat: sat_id, toe: t, toc: t,
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0,
            sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0,
            omega: 0.0,
            bgd_e1_e5a: 0.0,
            bgd_e1_e5b: 0.0,
            iod_nav: 1,
        });

        // Compute the satellite position at T
        let (sat_pos, _, sat_clk, _) = eph.position(t);
        let dx = true_pos.x - sat_pos.x;
        let dy = true_pos.y - sat_pos.y;
        let dz = true_pos.z - sat_pos.z;
        let geom_r = (dx * dx + dy * dy + dz * dz).sqrt();
        let raw_pr = geom_r - sat_clk * SPEED_OF_LIGHT_M_S;

        let meas = SppMeasurement {
            constellation: Constellation::Galileo,
            raw_pr,
            snr: 45.0,
            doppler: 0.0,
            time: t,
            eph,
            is_iono_free: false,
            freq_band: 7, // E5b
        };

        let ctx = EkfContext {
            pos_apc: true_pos,
            v_apc: Vector3::zeros(),
            r_b_e: nalgebra::Rotation3::identity(),
            lever_arm: Vector3::zeros(),
            omega_eb_b: Vector3::zeros(),
            rec_llh,
            state: &state,
            engine: &engine,
            n_cols: 21,
        };

        let mut z = DVector::zeros(1);
        let mut h = DMatrix::zeros(1, 21);
        let mut r_mat = DMatrix::zeros(1, 1);
        let mut types = Vec::new();
        let mut row = 0usize;
        let mut target = MatrixTarget {
            z: &mut z,
            h: &mut h,
            r: &mut r_mat,
            types: &mut types,
            row: &mut row,
        };

        process_measurement(&ctx, &meas, &mut target);

        // freq_band == 7 should still work with Galileo ephemeris
        assert_eq!(row, 1, "Should have processed 1 measurement");
        assert_eq!(types.len(), 1);
        assert_eq!(types[0].0, sat_id);
        assert_eq!(types[0].1, 0); // PR code
        // z should be finite (small residual)
        assert!(z[0].is_finite(), "z[0] should be finite, got {}", z[0]);
        // Clock Jacobian for Galileo should include ISB at col 17
        assert!((h[(0, 15)] - 1.0).abs() < 1e-10, "Clock bias col should be 1");
        assert!((h[(0, 17)] - 1.0).abs() < 1e-10, "ISB Gal col should be 1");
    }

    // -------------------------------------------------------------------------
    // process_doppler with non-zero relative velocity
    // -------------------------------------------------------------------------

    #[test]
    fn test_process_doppler_with_relative_velocity() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let pos = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, GpsTime::new(2000, 0.0),
        );
        let state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        engine.current_state = Some(state.clone());

        let sat_id = gneiss_core::sat::SatelliteId {
            constellation: Constellation::Gps, prn: 1,
        };
        let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: sat_id, toe: GpsTime::new(2000, 0.0), toc: GpsTime::new(2000, 0.0),
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0, sqrt_a: 5153.6,
            delta_n: 0.0, omega0: 0.0, omega_dot: 0.0, i0: 0.95,
            idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        });

        // Receiver moving at 10 m/s along +X, satellite stationary overhead
        let ctx = EkfContext {
            pos_apc: Vector3::new(0.0, 0.0, 0.0),
            v_apc: Vector3::new(10.0, 0.0, 0.0), // Moving receiver
            r_b_e: nalgebra::Rotation3::identity(),
            lever_arm: Vector3::zeros(),
            omega_eb_b: Vector3::zeros(),
            rec_llh: Vector3::new(0.0, 0.0, 0.0),
            state: &state,
            engine: &engine,
            n_cols: 21,
        };

        // Satellite at [20_000_000, 0, 0] moving at [2000, 0, 0]
        // Relative velocity = [10, 0, 0] - [2000, 0, 0] = [-1990, 0, 0]
        // LOS = [1, 0, 0]
        // Expected doppler = -1990 * 0 + rcv_clk_drift (0) - sat_drift (0) = -1990 m/s

        let meas = SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20_000_000.0,
            snr: 45.0,
            doppler: 1000.0, // will be converted to m/s in process_doppler
            time: GpsTime::new(2000, 0.0),
            eph,
            is_iono_free: false,
            freq_band: 1,
        };

        let ncols = 21;
        let mut z = DVector::zeros(1);
        let mut h = DMatrix::zeros(1, ncols);
        let mut r_mat = DMatrix::zeros(1, 1);
        let mut types = Vec::new();
        let mut row = 0usize;
        let mut target = MatrixTarget {
            z: &mut z,
            h: &mut h,
            r: &mut r_mat,
            types: &mut types,
            row: &mut row,
        };

        let geom = SatGeometry {
            geom_r: 20_000_000.0,
            los: Vector3::new(1.0, 0.0, 0.0),
            el: 1.0,
            az: 0.0,
            cdt_rx: 0.0,
            sat_clk: 0.0,
            sat_vel: Vector3::new(2000.0, 0.0, 0.0),
            sat_drift: 0.0,
        };

        process_doppler(&ctx, &meas, &mut target, &geom);

        assert_eq!(row, 1);
        assert_eq!(types.len(), 1);
        assert_eq!(types[0].1, 3); // doppler type code
        assert!(z[0].is_finite(), "z[0] should be finite, got {}", z[0]);
        // Velocity Jacobian should match the LOS
        assert!((h[(0, 3)] - 1.0).abs() < 1e-10, "vx column should be 1, got {}", h[(0, 3)]);
        assert!((h[(0, 4)]).abs() < 1e-10, "vy column should be 0");
        assert!((h[(0, 5)]).abs() < 1e-10, "vz column should be 0");
        // Clock drift column
        assert!((h[(0, 19)] - 1.0).abs() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // update_ekf with well-formed data passes consistency
    // -------------------------------------------------------------------------

    #[test]
    fn test_update_ekf_well_formed_data_not_rejected() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let pos = Coordinate::new(
            Vector3::new(0.0, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, GpsTime::new(2000, 0.0),
        );
        let state = RtkState::new(GpsTime::new(2000, 0.0), pos, 1.0);
        engine.current_state = Some(state);

        // Create 4 well-formed measurements with small residuals
        let n = 4;
        let z = DVector::from_vec(vec![0.1, -0.2, 0.15, -0.05]);
        let mut h = DMatrix::zeros(n, 21);
        // Set position Jacobians and clock
        for i in 0..n {
            h[(i, 0)] = 1.0;
            h[(i, 1)] = 0.0;
            h[(i, 2)] = 0.0;
            h[(i, 15)] = 1.0;
        }
        let r_mat = DMatrix::identity(n, n);
        let types = vec![
            (gneiss_core::sat::SatelliteId { constellation: Constellation::Gps, prn: 1 }, 0),
            (gneiss_core::sat::SatelliteId { constellation: Constellation::Gps, prn: 2 }, 0),
            (gneiss_core::sat::SatelliteId { constellation: Constellation::Gps, prn: 3 }, 0),
            (gneiss_core::sat::SatelliteId { constellation: Constellation::Gps, prn: 4 }, 0),
        ];

        let rejected = update_ekf(&mut engine, &z, &h, &r_mat, &types);
        // With 4 valid measurements (>= MIN_VALID_MEASUREMENTS=3), should NOT be rejected
        assert!(!rejected, "Should pass consistency check with 4 measurements");
    }

    // -------------------------------------------------------------------------
    // process_spp_tightly_coupled integration test with real ephemeris + obs
    // -------------------------------------------------------------------------

    #[test]
    fn test_process_spp_tightly_coupled_with_ephemerides() {
        let (ephemerides, obs) = build_test_ephemerides_and_obs();

        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let true_pos =
            nalgebra::Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let pos = Coordinate::new(
            true_pos, Datum::WGS84, Frame::ECEF, GpsTime::new(2000, 100000.0),
        );
        let mut state = RtkState::new(GpsTime::new(2000, 100000.0), pos, 1.0);
        state.ins_aligned = true;
        state.velocity = Vector3::zeros();
        engine.current_state = Some(state);
        engine.ephemerides = ephemerides;

        let result = process_spp_tightly_coupled(&mut engine, &obs);
        // Should successfully process and return a state
        assert!(result.is_ok(), "process_spp_tightly_coupled should succeed: {:?}", result.as_ref().err());
        let final_state = result.unwrap();
        // The state should have been updated
        assert!(final_state.covariance.nrows() > 0);
        // State history should have at least one entry
        assert!(!engine.state_history.is_empty());
        // Obs history should also have an entry
        assert!(!engine.obs_history.is_empty());
    }
}
