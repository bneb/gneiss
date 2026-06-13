use gneiss_core::obs::EpochObs;
use crate::filter::RtkState;
use crate::engine::{EngineError, ProcessingEngine};
use crate::engine::processed_sat::ProcessedSat;
use crate::engine::ppp_fg::PppFactorGraph;
use nalgebra::Vector3;

const LIGHT_SPEED: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

pub fn process_ppp<'a>(engine: &'a mut ProcessingEngine, rover_obs: &'a EpochObs) -> Result<&'a RtkState, EngineError> {
    if !valid_pos(engine) {
        return engine.process_spp(rover_obs);
    }
    
    let dt = rover_obs.time - engine.current_state.as_ref().unwrap().time;
    engine.predict_state(dt);
    let state = engine.current_state.as_mut().unwrap();
    state.time = rover_obs.time;
    state.position.epoch = rover_obs.time;

    let sats = build_sats(engine, rover_obs);
    if sats.is_empty() {
        return Err(EngineError::InsufficientSatellites);
    }

    let state = engine.current_state.as_mut().unwrap();
    update_phase_ambiguities(state, &sats, rover_obs.time);
    
    let fg = PppFactorGraph::new();
    fg.solve(state, &sats)?;
    
    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
    
    // We need to clone state BEFORE taking a reference again
    let final_state = engine.current_state.as_ref().unwrap().clone();
    engine.state_history.push(final_state);
    engine.obs_history.push((rover_obs.clone(), None));
    
    Ok(engine.current_state.as_ref().unwrap())
}

fn valid_pos(engine: &ProcessingEngine) -> bool {
    if let Some(state) = &engine.current_state {
        state.position.vector.norm().is_normal() && state.position.vector.norm() >= 1000.0
    } else {
        false
    }
}

fn build_sats<'a>(engine: &ProcessingEngine, rover_obs: &'a EpochObs) -> Vec<ProcessedSat<'a>> {
    let mut sats = Vec::new();
    let state = engine.current_state.as_ref().unwrap();
    let mut rcv_pos_ecef = Vector3::new(state.position.vector.x, state.position.vector.y, state.position.vector.z);
    rcv_pos_ecef += gneiss_core::tides::solid_earth_tides_ecef(rover_obs.time, rcv_pos_ecef);
    let rcv_pos_llh = gneiss_core::coords::ecef_to_llh(rcv_pos_ecef);

    for sat_obs in &rover_obs.satellites {
        let pr1 = match sat_obs.get_observable(1) {
            Some(p1) => p1,
            None => continue,
        };
        let eph = match engine.ephemerides.iter().find(|e| e.sat() == sat_obs.sat) {
            Some(e) => e, None => continue,
        };

        let f1 = gneiss_core::signal::satellite_frequencies(sat_obs.sat, eph.freq_num()).0;
        
        let tau_pr = pr1 / LIGHT_SPEED;
        let t_tx_nom = gneiss_core::time::GpsTime::new(rover_obs.time.week, rover_obs.time.tow - tau_pr);
        let (_, _, dt_s, _) = eph.position(t_tx_nom);
        let t_tx_true = gneiss_core::time::GpsTime::new(rover_obs.time.week, rover_obs.time.tow - tau_pr - dt_s);
        let (raw_vec, raw_vel, sat_clk, sat_drift) = eph.position(t_tx_true);
        
        let mut sat_pos = raw_vec;
        let mut sat_vel = raw_vel;
        for _ in 0..2 {
            let geometric_range = (sat_pos - rcv_pos_ecef).norm();
            let true_tau = geometric_range / LIGHT_SPEED;
            let theta = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * true_tau;
            let cos_t = libm::cos(theta);
            let sin_t = libm::sin(theta);
            sat_pos = nalgebra::Vector3::new(
                raw_vec.x * cos_t + raw_vec.y * sin_t,
                -raw_vec.x * sin_t + raw_vec.y * cos_t,
                raw_vec.z
            );
            sat_vel = nalgebra::Vector3::new(
                raw_vel.x * cos_t + raw_vel.y * sin_t,
                -raw_vel.x * sin_t + raw_vel.y * cos_t,
                raw_vel.z
            );
        }
        
        let dist = (sat_pos - rcv_pos_ecef).norm();
        let (az, el) = gneiss_core::coords::az_el(rcv_pos_llh, rcv_pos_ecef, sat_pos);
        if el < libm::asin(0.261799) { continue; }

        let tropo_params = gneiss_core::atmosphere::TropoParams::default();
        let z_dry = 0.0022768 * tropo_params.press_hpa / (1.0 - 0.00266 * libm::cos(2.0 * rcv_pos_llh.x) - 0.00028 * rcv_pos_llh.z / 1000.0);
        let tropo_dry = z_dry / libm::sin(el);

        let klobuchar = engine.klobuchar_params.unwrap_or_default();
        let iono_delay = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(&klobuchar, rcv_pos_llh, az, el, rover_obs.time);

        sats.push(ProcessedSat {
            sat_obs, dt_sat_m: sat_clk * LIGHT_SPEED, p_meas: pr1, is_iono_free: false,
            cp1: sat_obs.get_observable_phase(1), cp2: None,
            los: (sat_pos - rcv_pos_ecef) / dist, dist, el, snr: sat_obs.get_snr(1).unwrap_or(45) as f64,
            doppler: sat_obs.get_doppler(1).unwrap_or(0.0), lam1: LIGHT_SPEED / f1, lam2: LIGHT_SPEED / f1,
            tropo_dry, map_wet: 1.0 / libm::sin(el), iono_delay,
            f1, f2: f1, sat_pos_rot: sat_pos, sat_vel, sat_clock_drift: sat_drift,
            rcv_pos_ecef, pcv_correction: 0.0,
        });
    }
    sats
}

fn update_phase_ambiguities(state: &mut RtkState, sats: &[ProcessedSat], time: gneiss_core::time::GpsTime) {
    for sat in sats {
        if let Some(cp1_cyc) = sat.cp1 {
            if cp1_cyc == 0.0 { continue; }
            let windup = gneiss_core::windup::phase_windup(sat.sat_pos_rot, gneiss_core::sun::sun_position_ecef(time), sat.rcv_pos_ecef, *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0));
            state.windup.insert(sat.sat_obs.sat, windup);

            let l1_m = (cp1_cyc + windup) * sat.lam1;
            
            let mut slip = false;
            
            let lk_obs = sat.sat_obs.get_locktime(1);
            let prev = *state.locktimes.get(&(sat.sat_obs.sat, 1)).unwrap_or(&0);
            let mut new_lk = prev.saturating_add(1);
            if let Some(lk) = lk_obs {
                if lk == 0 || lk < prev { slip = true; new_lk = lk; } else { new_lk = lk; }
            } else if new_lk == 0 && state.locktimes.contains_key(&(sat.sat_obs.sat, 1)) {
                slip = true;
            }
            state.locktimes.insert((sat.sat_obs.sat, 1), new_lk);

            if slip { state.remove_ambiguity(sat.sat_obs.sat, 0); }

            let isb = match sat.sat_obs.sat.constellation {
                gneiss_core::sat::Constellation::Glonass => state.isb_glo,
                gneiss_core::sat::Constellation::Galileo => state.isb_gal,
                gneiss_core::sat::Constellation::Beidou => state.isb_bds,
                _ => 0.0,
            };

            let expected_p = sat.dist + state.rcv_clk_bias + isb - sat.dt_sat_m + sat.tropo_dry + state.zwd * sat.map_wet;
            let expected_with_iono = expected_p - sat.iono_delay;
            if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 0)) {
                state.add_ambiguity(sat.sat_obs.sat, 0, l1_m - expected_with_iono, 10000.0);
            }
            state.last_observed.insert((sat.sat_obs.sat, 0), state.epoch_count as u32);
        }
    }
}
