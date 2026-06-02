use gneiss_core::obs::EpochObs;
use crate::filter::RtkState;
use crate::engine::{EngineError, ProcessingEngine};
use nalgebra::{Vector3, DMatrix, DVector};

const LIGHT_SPEED: f64 = 299792458.0;

fn snr_scale(snr: i32) -> f64 {
    (10.0f64).powf((45.0 - snr as f64) / 10.0)
}

pub fn process_ppp<'a>(engine: &'a mut ProcessingEngine, rover_obs: &EpochObs) -> Result<&'a RtkState, EngineError> {
    let valid_pos = {
        if engine.current_state.is_none() {
            return engine.process_spp(rover_obs);
        }
        let state = engine.current_state.as_ref().unwrap();
        state.position.vector.norm().is_normal() && state.position.vector.norm() >= 1000.0
    };
    if !valid_pos {
        return engine.process_spp(rover_obs);
    }
    let state = engine.current_state.as_mut().unwrap();
    
    let mut h_rows = Vec::new();
    let mut z_vec = Vec::new();
    let mut r_vec = Vec::new();
    let mut sat_vec = Vec::new();

    let mut sats_to_process = Vec::new();

    for (_idx, sat_obs) in rover_obs.satellites.iter().enumerate() {
        let pr1 = sat_obs.get_observable(1);
        let pr2 = sat_obs.get_observable(2);
        let cp1 = sat_obs.get_observable_phase(1);
        let cp2 = sat_obs.get_observable_phase(2);
        
        if pr1.is_none() || pr2.is_none() { continue; }
        let eph = match engine.ephemerides.iter().find(|e| e.sat() == sat_obs.sat) {
            Some(e) => e,
            None => continue,
        };
        
        let pr1 = sat_obs.get_observable(1);
        let pr2 = sat_obs.get_observable(2);
        let cp1 = sat_obs.get_observable_phase(1);
        let cp2 = sat_obs.get_observable_phase(2);

        if pr1.is_none() || pr2.is_none() || cp1.is_none() || cp2.is_none() {
            continue;
        }

        let f1 = gneiss_core::signal::satellite_frequencies(sat_obs.sat, eph.freq_num()).0;
        let f2 = gneiss_core::signal::satellite_frequencies(sat_obs.sat, eph.freq_num()).1;
        let p_if = crate::combinations::iono_free(pr1.unwrap(), pr2.unwrap(), f1, f2);

        let mut rcv_pos_ecef = Vector3::new(state.position.vector.x, state.position.vector.y, state.position.vector.z);
        let set_disp = gneiss_core::tides::solid_earth_tides_ecef(rover_obs.time, rcv_pos_ecef);
        rcv_pos_ecef += set_disp;
        let _rcv_pos_llh = gneiss_core::coords::ecef_to_llh(rcv_pos_ecef);
        
        let (sat_pos, _sat_vel, sat_clk, _) = eph.position(rover_obs.time);
        let sat_pos_rot = sat_pos;

        let dist = (sat_pos_rot - rcv_pos_ecef).norm();
        let rcv_pos_llh = gneiss_core::coords::ecef_to_llh(rcv_pos_ecef);
        let (_az, el) = gneiss_core::coords::az_el(rcv_pos_llh, rcv_pos_ecef, sat_pos_rot);

        if el < libm::asin(0.261799) { continue; }

        let dt_sat_m = sat_clk * LIGHT_SPEED;
        let snr = sat_obs.get_snr(1).unwrap_or(45);

        let los = (sat_pos_rot - rcv_pos_ecef) / dist;
        let tropo_params = gneiss_core::atmosphere::TropoParams::default();
        let z_dry = 0.0022768 * tropo_params.press_hpa / (1.0 - 0.00266 * libm::cos(2.0 * rcv_pos_llh.x) - 0.00028 * rcv_pos_llh.z / 1000.0);
        let map_wet = 1.0 / libm::sin(el);
        let tropo_dry = z_dry / libm::sin(el);

        let lam1 = LIGHT_SPEED / gneiss_core::signal::satellite_frequencies(sat_obs.sat, eph.freq_num()).0;
        let lam2 = LIGHT_SPEED / gneiss_core::signal::satellite_frequencies(sat_obs.sat, eph.freq_num()).1;

        sats_to_process.push((sat_obs, dt_sat_m, p_if, cp1, cp2, los, dist, el, snr, lam1, lam2, tropo_dry, map_wet, f1, f2, sat_pos_rot, rcv_pos_ecef));
    }

    // PASS 1: State modifications
    for (sat_obs, dt_sat_m, p_if, cp1, cp2, _los, dist, _el, _snr, lam1, lam2, tropo_dry, map_wet, f1, f2, sat_pos_rot, rcv_pos_ecef) in &sats_to_process {
        let expected_p = dist + state.rcv_clk_bias - dt_sat_m + tropo_dry + state.zwd * map_wet;
        let res_p = p_if - expected_p;

        if let (Some(cp1_cyc), Some(cp2_cyc)) = (cp1, cp2) {
            if *cp1_cyc != 0.0 && *cp2_cyc != 0.0 {
                let windup = gneiss_core::windup::phase_windup(*sat_pos_rot, gneiss_core::sun::sun_position_ecef(rover_obs.time), *rcv_pos_ecef, *state.windup.get(&sat_obs.sat).unwrap_or(&0.0));
                state.windup.insert(sat_obs.sat, windup);

                let l1_m = (cp1_cyc + windup) * lam1;
                let l2_m = (cp2_cyc + windup) * lam2;
                let l_if = crate::combinations::iono_free(l1_m, l2_m, *f1, *f2);
                let l_gf = l1_m - l2_m;
                
                let mut slip = false;
                if let Some(&prev_gf) = state.gf_values.get(&sat_obs.sat) {
                    if (l_gf - prev_gf).abs() > 0.05 { slip = true; }
                }
                state.gf_values.insert(sat_obs.sat, l_gf);
                
                let l1_lock = sat_obs.get_locktime(1);
                let prev_l1 = *state.locktimes.get(&(sat_obs.sat, 1)).unwrap_or(&0);
                let mut new_l1 = prev_l1.saturating_add(1);
                if let Some(lk) = l1_lock {
                    if lk == 0 || (lk as u16) < prev_l1 { slip = true; new_l1 = lk as u16; } else { new_l1 = lk as u16; }
                } else if new_l1 == 0 && state.locktimes.contains_key(&(sat_obs.sat, 1)) { slip = true; }
                state.locktimes.insert((sat_obs.sat, 1), new_l1);
                
                let l2_lock = sat_obs.get_locktime(2);
                let prev_l2 = *state.locktimes.get(&(sat_obs.sat, 2)).unwrap_or(&0);
                let mut new_l2 = prev_l2.saturating_add(1);
                if let Some(lk) = l2_lock {
                    if lk == 0 || (lk as u16) < prev_l2 { slip = true; new_l2 = lk as u16; } else { new_l2 = lk as u16; }
                } else if new_l2 == 0 && state.locktimes.contains_key(&(sat_obs.sat, 2)) { slip = true; }
                state.locktimes.insert((sat_obs.sat, 2), new_l2);

                if slip { state.remove_ambiguity(sat_obs.sat, 0); }

                let res_l = l_if - expected_p;
                if !state.ambiguity_keys.contains(&(sat_obs.sat, 0)) {
                    let amb_est = res_l - res_p;
                    state.add_ambiguity(sat_obs.sat, 0, amb_est, 100.0);
                }
                state.last_observed.insert((sat_obs.sat, 0), state.epoch_count as u32);
            }
        }
    }

    // PASS 2: Measurement generation
    for (sat_obs, dt_sat_m, p_if, cp1, cp2, los, dist, el, snr, lam1, lam2, tropo_dry, map_wet, f1, f2, _sat_pos_rot, _rcv_pos_ecef) in &sats_to_process {
        let expected_p = dist + state.rcv_clk_bias - dt_sat_m + tropo_dry + state.zwd * map_wet;
        let res_p = p_if - expected_p;

        let mut h_row_p = vec![0.0; state.covariance.ncols()];
        h_row_p[0] = -los.x;
        h_row_p[1] = -los.y;
        h_row_p[2] = -los.z;
        h_row_p[17] = *map_wet;
        h_row_p[15] = 1.0;

        h_rows.push(h_row_p);
        z_vec.push(res_p);
        sat_vec.push(sat_obs.sat);
        
        let var_p = 9.0 * snr_scale(*snr as i32) / libm::sin(*el);
        r_vec.push(var_p);

        if let (Some(cp1_cyc), Some(cp2_cyc)) = (cp1, cp2) {
            if *cp1_cyc != 0.0 && *cp2_cyc != 0.0 {
                let windup = *state.windup.get(&sat_obs.sat).unwrap_or(&0.0);
                let l1_m = (cp1_cyc + windup) * lam1;
                let l2_m = (cp2_cyc + windup) * lam2;
                let l_if = crate::combinations::iono_free(l1_m, l2_m, *f1, *f2);

                let amb_idx = state.ambiguity_keys.iter().position(|&(s, f)| s == sat_obs.sat && f == 0).unwrap();
                let n_if_est = state.ambiguities[amb_idx];

                let mut h_row_l = vec![0.0; state.covariance.ncols()];
                h_row_l[0] = -los.x;
                h_row_l[1] = -los.y;
                h_row_l[2] = -los.z;
                h_row_l[17] = *map_wet;
                h_row_l[15] = 1.0;
                h_row_l[crate::filter::CORE_STATE_SIZE + amb_idx] = 1.0;

                let expected_l = expected_p + n_if_est;
                let res_l = l_if - expected_l;

                h_rows.push(h_row_l);
                z_vec.push(res_l);
                sat_vec.push(sat_obs.sat);
                
                let f1_4 = f1 * f1 * f1 * f1;
                let f2_4 = f2 * f2 * f2 * f2;
                let if_factor = (f1_4 + f2_4) / (f1 * f1 - f2 * f2).powi(2);
                let var_l = if_factor * 0.0001 * snr_scale(*snr as i32) / libm::sin(*el);
                r_vec.push(var_l);
            }
        }
    }

    if z_vec.is_empty() {
        return Err(EngineError::InsufficientSatellites);
    }


    // Clock jump detection
    let mut pr_residuals = Vec::new();
    for i in 0..z_vec.len() {
        // Pseudorange rows have 1.0 in h_rows[i][15] and 0.0 in ambiguity states
        let is_pr = h_rows[i][15] == 1.0 && !h_rows[i][crate::filter::CORE_STATE_SIZE..].iter().any(|&x| x != 0.0);
        if is_pr {
            pr_residuals.push(z_vec[i]);
        }
    }
    
    if !pr_residuals.is_empty() {
        pr_residuals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median_res = pr_residuals[pr_residuals.len() / 2];
        
        if median_res.abs() > 100_000.0 {
            tracing::warn!("Clock jump detected! Median residual = {:.2}m", median_res);
            
            // Shift clock bias
            state.rcv_clk_bias += median_res;
            
            // Shift all ambiguities to counteract clock shift in phase
            for i in 0..state.ambiguities.len() {
                state.ambiguities[i] -= median_res;
            }
            
            // Expand clock covariance to allow re-convergence
            state.covariance[(15, 15)] += 1e6;
            
            // Apply shift to pseudorange residuals
            for i in 0..z_vec.len() {
                let is_pr = h_rows[i][15] == 1.0 && !h_rows[i][crate::filter::CORE_STATE_SIZE..].iter().any(|&x| x != 0.0);
                if is_pr {
                    z_vec[i] -= median_res;
                }
            }
        }
    }

    let h = DMatrix::from_fn(z_vec.len(), state.covariance.ncols(), |r, c| h_rows[r][c]);

    let z = DVector::from_vec(z_vec);
    let mut r_mat = DMatrix::zeros(r_vec.len(), r_vec.len());
    for i in 0..r_vec.len() {
        r_mat[(i, i)] = r_vec[i];
    }

    crate::engine::updater::update(state, &z, &h, &r_mat).map_err(|e| {
        tracing::error!("PPP Update Error: {:?}", e);
        EngineError::StateDisappeared
    })?;

    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
    Ok(state)
}
