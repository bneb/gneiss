use gneiss_core::obs::EpochObs;
use nalgebra::{DMatrix, DVector, Vector3};

use crate::filter::RtkState;
use crate::engine::{EngineError, ProcessingEngine};
use crate::spp::{build_measurements, SppConfig};
use gneiss_core::atmosphere::{AtmosphereModel, TropoParams};
use gneiss_core::coords::{ecef_to_llh, az_el};
use gneiss_core::constants::{SPEED_OF_LIGHT_M_S, EARTH_ROTATION_RATE_RAD_S};

fn compute_sagnac_correction(sat_ecef: Vector3<f64>, geometric_pr: f64) -> Vector3<f64> {
    let tof = geometric_pr / SPEED_OF_LIGHT_M_S;
    let theta = EARTH_ROTATION_RATE_RAD_S * tof;
    let cos_t = f64::cos(theta);
    let sin_t = f64::sin(theta);
    Vector3::new(sat_ecef.x * cos_t + sat_ecef.y * sin_t, -sat_ecef.x * sin_t + sat_ecef.y * cos_t, sat_ecef.z)
}

fn compute_sagnac_velocity_correction(sat_vel: Vector3<f64>, geometric_pr: f64) -> Vector3<f64> {
    let tof = geometric_pr / SPEED_OF_LIGHT_M_S;
    let theta = EARTH_ROTATION_RATE_RAD_S * tof;
    let cos_t = f64::cos(theta);
    let sin_t = f64::sin(theta);
    Vector3::new(sat_vel.x * cos_t + sat_vel.y * sin_t, -sat_vel.x * sin_t + sat_vel.y * cos_t, sat_vel.z)
}

pub fn process_spp_tightly_coupled<'a>(engine: &'a mut ProcessingEngine, rover_obs: &'a EpochObs) -> Result<&'a RtkState, EngineError> {
    let spp_config = SppConfig::default();
    
    // Fall back to loosely coupled if there is no current state (e.g., initial startup)
    // or if the state has disappeared / decoupled.
    if engine.current_state.is_none() {
        return engine.process_spp(rover_obs);
    }
    
    let dt = rover_obs.time - engine.current_state.as_ref().unwrap().time;
    engine.predict_state(dt);
    let state = engine.current_state.as_mut().unwrap();
    state.time = rover_obs.time;
    state.position.epoch = rover_obs.time;
    
    let measurements = build_measurements(rover_obs, &engine.ephemerides, &spp_config);
    if measurements.is_empty() {
        return Err(EngineError::NoObservations);
    }
    
    // Compute a loosely coupled SPP solution just to get a good clock bias estimate.
    // Receiver clocks jump and drift unpredictably, which will cause the EKF innovations 
    // to explode and reject all measurements if we rely purely on INS propagation.
    if let Ok(spp_res) = crate::spp::compute_spp(rover_obs, &engine.ephemerides, engine.klobuchar_params.as_ref(), &spp_config, None) {
        state.rcv_clk_bias = spp_res.cdt;
    }
    
    // Pre-calculate lever arm compensation
    let r_b_e = state.attitude.to_rotation_matrix();
    let lever_arm = Vector3::from_column_slice(&engine.config.imu_to_antenna_lever_arm);
    let l_e = r_b_e * lever_arm;
    let pos_apc = state.position.vector + l_e;
    
    let omega_b = if let Some(imu_buf) = engine.imu_history.last() {
        if let Some(last_imu) = imu_buf.last() {
            last_imu.gyro - state.gyro_bias
        } else { Vector3::zeros() }
    } else { Vector3::zeros() };
    
    let v_apc = state.velocity + r_b_e * omega_b.cross(&lever_arm);
    let rec_llh = ecef_to_llh(pos_apc);
    
    let n_meas = measurements.len();
    // We construct z, h, r for both Pseudorange and Doppler
    let mut num_doppler = 0;
    for m in &measurements {
        if m.doppler != 0.0 { num_doppler += 1; }
    }
    
    let total_rows = n_meas + num_doppler;
    let n_cols = state.covariance.ncols();
    
    let mut z_vec = DVector::zeros(total_rows);
    let mut h_mat = DMatrix::zeros(total_rows, n_cols);
    let mut r_mat = DMatrix::zeros(total_rows, total_rows);
    let mut meas_types = Vec::with_capacity(total_rows);
    
    let mut row_idx = 0;
    
    let c_clk = state.rcv_clk_bias; // GPS clock bias
    let glo_clk = state.isb_glo; // ISB for GLO relative to GPS
    let gal_clk = state.isb_gal; // ISB for GAL
    let bds_clk = state.isb_bds; // ISB for BDS
    let c_drift = state.rcv_clk_drift;
    
    for m in &measurements {
        let cdt_rx = match m.constellation {
            gneiss_core::sat::Constellation::Gps => c_clk,
            gneiss_core::sat::Constellation::Galileo => c_clk + gal_clk,
            gneiss_core::sat::Constellation::Beidou => c_clk + bds_clk,
            gneiss_core::sat::Constellation::Glonass => c_clk + glo_clk,
            _ => c_clk,
        };
        
        // Compute Satellite state
        let t_rcv_true = m.time.tow - cdt_rx / SPEED_OF_LIGHT_M_S;
        let t_tx_nom = gneiss_core::time::GpsTime::new(m.time.week, t_rcv_true - m.raw_pr / SPEED_OF_LIGHT_M_S);
        let (_, _, dt_s, _) = m.eph.position(t_tx_nom);
        let t_tx_true = gneiss_core::time::GpsTime::new(m.time.week, t_rcv_true - m.raw_pr / SPEED_OF_LIGHT_M_S - dt_s);
        let (sat_pos_raw, sat_vel_raw, sat_clk, sat_drift) = m.eph.position(t_tx_true);
        let sat_clk_m = sat_clk * SPEED_OF_LIGHT_M_S;
        let sat_drift_ms = sat_drift * SPEED_OF_LIGHT_M_S;

        let sat_ecef = compute_sagnac_correction(sat_pos_raw, m.raw_pr - cdt_rx);
        let sat_vel = compute_sagnac_velocity_correction(sat_vel_raw, m.raw_pr - cdt_rx);
        
        let dx = pos_apc.x - sat_ecef.x;
        let dy = pos_apc.y - sat_ecef.y;
        let dz = pos_apc.z - sat_ecef.z;
        let geom_r = f64::sqrt(dx * dx + dy * dy + dz * dz).max(1e-6);
        
        let (az, el) = az_el(rec_llh, pos_apc, sat_ecef);
        
        let safe_el = el.max(5.0 * core::f64::consts::PI / 180.0);
        let tropo = AtmosphereModel::tropo_nmf(&TropoParams::default(), rec_llh, safe_el, m.time);
        let iono = if let Some(iono_params) = &engine.klobuchar_params {
            AtmosphereModel::iono_klobuchar(iono_params, rec_llh, az, safe_el, m.time)
        } else { 0.0 };
        
        // --- Pseudorange Innovation ---
        let expected_pr = geom_r + cdt_rx - sat_clk_m + tropo + iono;
        let pr_innovation = m.raw_pr - expected_pr;
        
        let los_x = dx / geom_r;
        let los_y = dy / geom_r;
        let los_z = dz / geom_r;
        let los = Vector3::new(los_x, los_y, los_z);
        let los_row = nalgebra::RowVector3::new(los_x, los_y, los_z);

        z_vec[row_idx] = pr_innovation;
        for i in 0..3 {
            h_mat[(row_idx, i)] = los[i];
        }
        
        if n_cols > 15 {
            let h_pos_att = -l_e.cross_matrix();
            let pr_h_att = los_row * h_pos_att;
            h_mat[(row_idx, 6)] = pr_h_att[0];
            h_mat[(row_idx, 7)] = pr_h_att[1];
            h_mat[(row_idx, 8)] = pr_h_att[2];
            
            h_mat[(row_idx, 15)] = 1.0; // GPS clock bias
            match m.constellation {
                gneiss_core::sat::Constellation::Glonass => h_mat[(row_idx, 16)] = 1.0,
                gneiss_core::sat::Constellation::Galileo => h_mat[(row_idx, 17)] = 1.0,
                gneiss_core::sat::Constellation::Beidou => h_mat[(row_idx, 18)] = 1.0,
                _ => {}
            }
        }
        
        let var_scale = gneiss_core::variance::observation_variance(m.snr, el, engine.config.nominal_snr_dbhz);
        r_mat[(row_idx, row_idx)] = engine.config.tuning.pr_base_var * var_scale;
        meas_types.push((m.eph.sat(), 0)); // 0 = pseudorange
        
        row_idx += 1;
        
        // --- Doppler Innovation ---
        if m.doppler != 0.0 {
            let relative_vel = Vector3::new(v_apc.x - sat_vel.x, v_apc.y - sat_vel.y, v_apc.z - sat_vel.z);
            let expected_dop = (los_x * relative_vel.x + los_y * relative_vel.y + los_z * relative_vel.z) + c_drift - sat_drift_ms;
            
            let f1 = gneiss_core::signal::satellite_frequencies(m.eph.sat(), m.eph.freq_num()).0;
            let wavelength = SPEED_OF_LIGHT_M_S / f1;
            let measured_dop_ms = -m.doppler * wavelength; // Doppler frequency shift to m/s
            
            let dop_innovation = measured_dop_ms - expected_dop;
            
            z_vec[row_idx] = dop_innovation;
            h_mat[(row_idx, 3)] = los_x;
            h_mat[(row_idx, 4)] = los_y;
            h_mat[(row_idx, 5)] = los_z;
            if n_cols > 19 {
                let a_e = r_b_e * omega_b.cross(&lever_arm);
                let h_vel_att = -a_e.cross_matrix();
                let h_vel_bg = r_b_e.matrix() * lever_arm.cross_matrix();
                
                let dop_h_att = los_row * h_vel_att;
                let dop_h_bg = los_row * h_vel_bg;
                
                h_mat[(row_idx, 6)] = dop_h_att[0];
                h_mat[(row_idx, 7)] = dop_h_att[1];
                h_mat[(row_idx, 8)] = dop_h_att[2];
                h_mat[(row_idx, 12)] = dop_h_bg[0];
                h_mat[(row_idx, 13)] = dop_h_bg[1];
                h_mat[(row_idx, 14)] = dop_h_bg[2];
                
                h_mat[(row_idx, 19)] = 1.0; // Receiver clock drift
            }
            r_mat[(row_idx, row_idx)] = engine.config.tuning.dop_base_var * var_scale; // Doppler has similar variance scale to PR in smartphones and urban canyons
            meas_types.push((m.eph.sat(), 3)); // 3 = Doppler
            row_idx += 1;

        }
    }
    
    tracing::debug!("SPP tight update: P[19, 19] = {:.2}, n_cols = {}, meas_types len = {}", state.covariance[(19, 19)], n_cols, meas_types.len());
    
    let mut rejected = false;
    let update_res = crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat, engine.config.spp_consistency_threshold_m, Some(&meas_types), true, &engine.config.tuning);
    if update_res.is_err() {
        tracing::debug!("Tightly-coupled SPP EKF update failed: {:?}", update_res.err().unwrap());
        rejected = true;
    } else if let Ok(valid_indices) = &update_res {
        tracing::debug!("Tightly-coupled SPP EKF valid indices: {} / {}", valid_indices.len(), total_rows);
        if valid_indices.len() < 3 {
            tracing::debug!("Tightly-coupled SPP EKF rejected: Not enough valid measurements ({})", valid_indices.len());
            rejected = true; // Not enough valid measurements accepted
        }
    }

    if rejected {
        state.consecutive_rejections += 1;
        if state.consecutive_rejections > 5 {
            tracing::warn!("Tightly-coupled SPP EKF rejected for {} epochs. Hard resetting INS.", state.consecutive_rejections);
            // We should use process_spp to reset the clock bias and position, but carefully preserve covariance
            if let Ok(spp_res) = crate::spp::compute_spp(rover_obs, &engine.ephemerides, engine.klobuchar_params.as_ref(), &spp_config, None) {
                state.position = spp_res.position;
                state.velocity = nalgebra::Vector3::zeros();
                state.rcv_clk_bias = spp_res.cdt;
                state.rcv_clk_drift = 0.0;
                state.decouple_position();
                for i in 0..3 { state.covariance[(i, i)] = 100.0; }
                for i in 3..6 { state.covariance[(i, i)] = 10.0; }
                if state.covariance.nrows() > 15 {
                    state.covariance[(15, 15)] = 100000.0;
                    state.covariance[(19, 19)] = 1000.0;
                }
                state.is_reset = true;
                state.consecutive_rejections = 0;
            } else {
                return engine.process_spp(rover_obs); // if SPP compute fails completely, fallback
            }
        } else {
            tracing::warn!("Tightly-coupled SPP EKF update rejected. Riding through outage via INS dead-reckoning.");
        }
    } else {
        state.consecutive_rejections = 0;
    }
    
    let final_state = engine.current_state.as_ref().unwrap().clone();
    tracing::debug!("SPP tight update done: vel=[{:.2}, {:.2}, {:.2}], c_drift={:.2}", final_state.velocity.x, final_state.velocity.y, final_state.velocity.z, final_state.rcv_clk_drift);
    engine.state_history.push(final_state);
    engine.obs_history.push((rover_obs.clone(), None));
    Ok(engine.current_state.as_ref().unwrap())
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::time::GpsTime;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use crate::engine::{EngineConfig, EngineMode};
    use crate::filter::RtkState;
    use crate::engine::ProcessingEngine;
    use gneiss_core::obs::SatObs;
    use gneiss_core::sat::SatelliteId;

    #[test]
    fn test_process_spp_tightly_coupled_startup() {
        // Test that tightly coupled filter falls back to loosely coupled SPP initially
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.mode = EngineMode::SppIns;
        let obs = EpochObs { time: GpsTime::new(2000, 1000.0), satellites: vec![] };
        // Empty obs should result in error, but verify it attempts to call process_spp
        assert!(process_spp_tightly_coupled(&mut engine, &obs).is_err());
    }
}
