use gneiss_core::obs::EpochObs;
use crate::filter::RtkState;
use crate::engine::{EngineError, ProcessingEngine};
use crate::engine::processed_sat::ProcessedSat;
use crate::engine::ppp_fg::PppFactorGraph;
use nalgebra::Vector3;

const LIGHT_SPEED: f64 = 299792458.0;

pub fn process_ppp<'a>(engine: &'a mut ProcessingEngine, rover_obs: &'a EpochObs) -> Result<&'a RtkState, EngineError> {
    if !valid_pos(engine) {
        return engine.process_spp(rover_obs);
    }
    
    let sats = build_sats(engine, rover_obs);
    if sats.is_empty() {
        return Err(EngineError::InsufficientSatellites);
    }

    let state = engine.current_state.as_mut().unwrap();
    update_phase_ambiguities(state, &sats, rover_obs.time);
    
    let fg = PppFactorGraph::new();
    fg.solve(state, &sats)?;
    
    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
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
        let (pr1, pr2) = match (sat_obs.get_observable(1), sat_obs.get_observable(2)) {
            (Some(p1), Some(p2)) => (p1, p2),
            _ => continue,
        };
        let eph = match engine.ephemerides.iter().find(|e| e.sat() == sat_obs.sat) {
            Some(e) => e, None => continue,
        };

        let f1 = gneiss_core::signal::satellite_frequencies(sat_obs.sat, eph.freq_num()).0;
        let f2 = gneiss_core::signal::satellite_frequencies(sat_obs.sat, eph.freq_num()).1;
        let p_if = crate::combinations::iono_free(pr1, pr2, f1, f2);
        
        let (sat_pos, sat_vel, sat_clk, sat_drift) = eph.position(rover_obs.time);
        let dist = (sat_pos - rcv_pos_ecef).norm();
        let (_az, el) = gneiss_core::coords::az_el(rcv_pos_llh, rcv_pos_ecef, sat_pos);
        if el < libm::asin(0.261799) { continue; }

        let tropo_params = gneiss_core::atmosphere::TropoParams::default();
        let z_dry = 0.0022768 * tropo_params.press_hpa / (1.0 - 0.00266 * libm::cos(2.0 * rcv_pos_llh.x) - 0.00028 * rcv_pos_llh.z / 1000.0);
        let tropo_dry = z_dry / libm::sin(el);

        sats.push(ProcessedSat {
            sat_obs, dt_sat_m: sat_clk * LIGHT_SPEED, p_meas: p_if, is_iono_free: true,
            cp1: sat_obs.get_observable_phase(1), cp2: sat_obs.get_observable_phase(2),
            los: (sat_pos - rcv_pos_ecef) / dist, dist, el, snr: sat_obs.get_snr(1).unwrap_or(45) as f64,
            doppler: 0.0, lam1: LIGHT_SPEED / f1, lam2: LIGHT_SPEED / f2,
            tropo_dry, map_wet: 1.0 / libm::sin(el), iono_delay: 0.0,
            f1, f2, sat_pos_rot: sat_pos, sat_vel, sat_clock_drift: sat_drift,
            rcv_pos_ecef, pcv_correction: 0.0,
        });
    }
    sats
}

fn update_phase_ambiguities(state: &mut RtkState, sats: &[ProcessedSat], time: gneiss_core::time::GpsTime) {
    for sat in sats {
        if let (Some(cp1_cyc), Some(cp2_cyc)) = (sat.cp1, sat.cp2) {
            if cp1_cyc == 0.0 || cp2_cyc == 0.0 { continue; }
            let windup = gneiss_core::windup::phase_windup(sat.sat_pos_rot, gneiss_core::sun::sun_position_ecef(time), sat.rcv_pos_ecef, *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0));
            state.windup.insert(sat.sat_obs.sat, windup);

            let l1_m = (cp1_cyc + windup) * sat.lam1;
            let l2_m = (cp2_cyc + windup) * sat.lam2;
            let l_if = crate::combinations::iono_free(l1_m, l2_m, sat.f1, sat.f2);
            let l_gf = l1_m - l2_m;
            
            let mut slip = false;
            if let Some(&prev_gf) = state.gf_values.get(&sat.sat_obs.sat) {
                if (l_gf - prev_gf).abs() > 0.05 { slip = true; }
            }
            state.gf_values.insert(sat.sat_obs.sat, l_gf);
            
            if check_locktime(state, sat.sat_obs, 1) || check_locktime(state, sat.sat_obs, 2) {
                slip = true;
            }

            if slip { state.remove_ambiguity(sat.sat_obs.sat, 0); }

            let expected_p = sat.dist + state.rcv_clk_bias - sat.dt_sat_m + sat.tropo_dry + state.zwd * sat.map_wet;
            if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 0)) {
                state.add_ambiguity(sat.sat_obs.sat, 0, (l_if - expected_p) - (sat.p_meas - expected_p), 100.0);
            }
            state.last_observed.insert((sat.sat_obs.sat, 0), state.epoch_count as u32);
        }
    }
}

fn check_locktime(state: &mut RtkState, obs: &gneiss_core::obs::SatObs, band: u8) -> bool {
    let lk_obs = obs.get_locktime(band);
    let prev = *state.locktimes.get(&(obs.sat, band)).unwrap_or(&0);
    let mut new_lk = prev.saturating_add(1);
    let mut slip = false;
    if let Some(lk) = lk_obs {
        if lk == 0 || lk < prev { slip = true; new_lk = lk; } else { new_lk = lk; }
    } else if new_lk == 0 && state.locktimes.contains_key(&(obs.sat, band)) {
        slip = true;
    }
    state.locktimes.insert((obs.sat, band), new_lk);
    slip
}
