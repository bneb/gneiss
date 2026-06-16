use gneiss_core::obs::EpochObs;
use crate::filter::RtkState;
use crate::engine::{EngineError, ProcessingEngine};
use crate::engine::processed_sat::ProcessedSat;
use crate::engine::ppp_fg::PppFactorGraph;
use nalgebra::Vector3;
use chrono::TimeZone;

const LIGHT_SPEED: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

pub fn process_ppp<'a>(engine: &'a mut ProcessingEngine, rover_obs: &'a EpochObs) -> Result<&'a RtkState, EngineError> {
    if !valid_pos(engine) {
        return engine.process_spp(rover_obs);
    }
    
    let dt = rover_obs.time.tow - engine.current_state.as_ref().unwrap().time.tow;
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
    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
    
    let fg = PppFactorGraph::new();
    fg.solve(state, &sats)?;
    
    state.epoch_count += 1;
    
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
        let eph = match engine.ephemerides.iter().find(|e| e.sat() == sat_obs.sat) {
            Some(e) => e, None => continue,
        };
        let freqs = gneiss_core::signal::satellite_frequencies(sat_obs.sat, eph.freq_num());
        let f1 = freqs.0;
        let mut f2 = freqs.1;
        if f2 == 0.0 { f2 = f1; }

        let mut p1_opt = match sat_obs.sat.constellation {
            gneiss_core::sat::Constellation::Beidou => sat_obs.get_observable(2),
            _ => sat_obs.get_observable(1),
        };
        let mut p2_opt = match sat_obs.sat.constellation {
            gneiss_core::sat::Constellation::Galileo => sat_obs.get_observable(7).or(sat_obs.get_observable(5)),
            gneiss_core::sat::Constellation::Beidou => sat_obs.get_observable(7).or(sat_obs.get_observable(6)),
            _ => sat_obs.get_observable(2),
        };
        
        if let Some(p1) = p1_opt.as_mut() {
            if let Some(&dcb) = engine.dcbs.get(&(sat_obs.sat, "P1C1".to_string())) {
                *p1 -= dcb * 1e-9 * LIGHT_SPEED;
            }
        }
        if let Some(p2) = p2_opt.as_mut() {
            if let Some(&dcb) = engine.dcbs.get(&(sat_obs.sat, "P2C2".to_string())) {
                *p2 -= dcb * 1e-9 * LIGHT_SPEED;
            }
        }
        
        let mut pr1 = match p1_opt {
            Some(p) => p,
            None => continue,
        };
        let cp1 = match sat_obs.sat.constellation {
            gneiss_core::sat::Constellation::Beidou => sat_obs.get_observable_phase(2),
            _ => sat_obs.get_observable_phase(1),
        };
        let cp2 = match sat_obs.sat.constellation {
            gneiss_core::sat::Constellation::Galileo => sat_obs.get_observable_phase(7).or(sat_obs.get_observable_phase(5)),
            gneiss_core::sat::Constellation::Beidou => sat_obs.get_observable_phase(7).or(sat_obs.get_observable_phase(6)),
            _ => sat_obs.get_observable_phase(2),
        };
        let mut is_iono_free = false;

        if !engine.sp3_epochs.is_empty() || engine.clk_data.is_some() {
            if let (Some(p1), Some(p2)) = (p1_opt, p2_opt) {
                let gamma = (f1 * f1) / (f2 * f2);
                pr1 = (gamma * p1 - p2) / (gamma - 1.0);
                is_iono_free = true;
            }
        }
        
        // Allow single-frequency satellites in PPP by using Klobuchar ionosphere correction.
        // This is critical for consumer devices (like Pixel 4) where most satellites are single-frequency.
        
        let tau_pr = pr1 / LIGHT_SPEED;
        let t_tx_nom = gneiss_core::time::GpsTime::new(rover_obs.time.week, rover_obs.time.tow - tau_pr);
        
        let mut dt_s = 0.0;
        if let Some(clk_data) = &engine.clk_data {
            if let Some(bias) = clk_data.get_clock_bias(sat_obs.sat, t_tx_nom) {
                dt_s = bias;
            }
        }
        
        let brdc_clk = eph.position(t_tx_nom).2;
        let diff_clk = dt_s - brdc_clk;
        if diff_clk.abs() > 1e-8 { // > 3 meters
            tracing::debug!("DEBUG CLK DIFF: sat={}, sp3_clk={:.9}, brdc_clk={:.9}, diff={:.9}s ({:.3}m)", sat_obs.sat, dt_s, brdc_clk, diff_clk, diff_clk * LIGHT_SPEED);
        }
        
        let mut clk_found = false;
        if dt_s != 0.0 { clk_found = true; }
        
        let precise_mode = !engine.sp3_epochs.is_empty() || engine.clk_data.is_some();
        if !precise_mode {
            dt_s = brdc_clk;
            clk_found = true;
        }
        if precise_mode && !clk_found {
            if let Some((_, sp3_clk)) = crate::engine::ssr::get_precise_orbit(&engine.sp3_epochs, sat_obs.sat, t_tx_nom, 10) {
                if !sp3_clk.is_nan() && sp3_clk != 0.0 {
                    dt_s = sp3_clk;
                    clk_found = true;
                }
            }
        }

        if precise_mode && !clk_found {
            continue; // Skip satellite if precise clock is missing in PPP
        }

        let t_tx_true = gneiss_core::time::GpsTime::new(rover_obs.time.week, rover_obs.time.tow - tau_pr - dt_s);
        let (mut raw_vec, raw_vel, _, sat_drift) = eph.position(t_tx_true);
        let mut sat_clk = if precise_mode { dt_s } else { brdc_clk };
        
        if precise_mode {
            // Relativistic clock correction for precise orbits (broadcast already includes this via af0, af1, af2)
            let dt_rel = -2.0 * raw_vec.dot(&raw_vel) / (LIGHT_SPEED * LIGHT_SPEED);
            sat_clk += dt_rel;
            let mut sp3_found = false;
            let mut sp3_pos_val = nalgebra::Vector3::zeros();
            if precise_mode {
                if let Some((sp3_pos, _)) = crate::engine::ssr::get_precise_orbit(&engine.sp3_epochs, sat_obs.sat, t_tx_true, 10) {
                    raw_vec = sp3_pos;
                    sp3_pos_val = sp3_pos;
                    sp3_found = true;
                }
            }
            let brdc_pos = eph.position(t_tx_true).0;
            if !sp3_found {
                raw_vec = brdc_pos;
            } else {
                let diff = sp3_pos_val - brdc_pos;
                tracing::debug!("DEBUG POS DIFF: sat={}, sp3=[{:.3}, {:.3}, {:.3}], brdc=[{:.3}, {:.3}, {:.3}], diff=[{:.3}, {:.3}, {:.3}], norm={:.3}", 
                    sat_obs.sat, 
                    sp3_pos_val.x, sp3_pos_val.y, sp3_pos_val.z, 
                    brdc_pos.x, brdc_pos.y, brdc_pos.z, 
                    diff.x, diff.y, diff.z, 
                    diff.norm());
            }
            
            if !sp3_found {
                continue; // Skip satellite if precise orbit is missing in PPP
            }
        }
        
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
        if el < 0.261799 { continue; } // 15° elevation mask in radians

        let tropo_params = gneiss_core::atmosphere::TropoParams::default();
        let p_z = tropo_params.press_hpa * libm::pow(1.0 - 0.0000226 * rcv_pos_llh.z, 5.225);
        let z_dry = 0.0022768 * p_z / (1.0 - 0.00266 * libm::cos(2.0 * rcv_pos_llh.x) - 0.00028 * rcv_pos_llh.z / 1000.0);
        let (m_h, m_w) = gneiss_core::atmosphere::AtmosphereModel::nmf_mapping_functions(rcv_pos_llh, el, rover_obs.time);
        let tropo_dry = z_dry * m_h;

        let klobuchar = engine.klobuchar_params.unwrap_or_default();
        let mut iono_delay = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(&klobuchar, rcv_pos_llh, az, el, rover_obs.time);
        if is_iono_free {
            iono_delay = 0.0;
        }

        let mut pcv_correction = 0.0;
        if precise_mode {
            if let Some(db) = &engine.antex {
                let gps_epoch = chrono::Utc.with_ymd_and_hms(1980, 1, 6, 0, 0, 0).unwrap();
                let dt_seconds = rover_obs.time.week as i64 * 604800 + rover_obs.time.tow as i64;
                let utc_time = gps_epoch + chrono::Duration::seconds(dt_seconds);
                
                if let Some(ant) = db.find_satellite(&sat_obs.sat.to_string(), utc_time) {
                    let freq_code1 = match sat_obs.sat.constellation {
                        gneiss_core::sat::Constellation::Gps => "G01",
                        gneiss_core::sat::Constellation::Glonass => "R01",
                        gneiss_core::sat::Constellation::Galileo => "E01",
                        gneiss_core::sat::Constellation::Beidou => "C02",
                        _ => "G01",
                    };
                    let freq_code2 = match sat_obs.sat.constellation {
                        gneiss_core::sat::Constellation::Gps => "G02",
                        gneiss_core::sat::Constellation::Glonass => "R02",
                        gneiss_core::sat::Constellation::Galileo => "E05",
                        gneiss_core::sat::Constellation::Beidou => "C07",
                        _ => "G02",
                    };
                    
                    let get_pco = |code: &str| -> nalgebra::Vector3<f64> {
                        if let Some(pcv) = ant.frequencies.get(code) {
                            pcv.pco / 1000.0 // millimeters to meters
                        } else if let Some(pcv) = ant.frequencies.values().next() {
                            pcv.pco / 1000.0
                        } else {
                            nalgebra::Vector3::zeros()
                        }
                    };
                    
                    let pco1 = get_pco(freq_code1);
                    
                    let pco = if is_iono_free && sat_obs.get_observable_phase(2).is_some() {
                        let pco2 = get_pco(freq_code2);
                        let gamma = (f1 * f1) / (f2 * f2);
                        (pco1 * gamma - pco2) / (gamma - 1.0)
                    } else {
                        pco1
                    };
                    
                    let k = (rcv_pos_ecef - sat_pos).normalize();
                    let sat_z = -sat_pos.normalize();
                    let e_sun = (gneiss_core::sun::sun_position_ecef(rover_obs.time) - sat_pos).normalize();
                    let sat_y = sat_z.cross(&e_sun).normalize();
                    let sat_x = sat_y.cross(&sat_z).normalize();
                    
                    // In ANTEX, North=X, East=Y, Up=Z in satellite body frame
                    let pco_ecef = sat_x * pco.x + sat_y * pco.y + sat_z * pco.z;
                    
                    // The distance to the antenna phase center is |dist_vec - pco_ecef|
                    // which is roughly |dist_vec| - k.dot(&pco_ecef)
                    pcv_correction = -k.dot(&pco_ecef);
                }
            }
        }
        
        tracing::trace!("PPP sat={}, iono_free={}, dist={:.1}, clk_m={:.3}, pcv={:.4}", sat_obs.sat, is_iono_free, dist, sat_clk * LIGHT_SPEED, pcv_correction);

        sats.push(ProcessedSat {
            sat_obs, dt_sat_m: sat_clk * LIGHT_SPEED, p_meas: pr1, is_iono_free,
            cp1, cp2,
            los: (sat_pos - rcv_pos_ecef) / dist, dist: dist - pcv_correction, el, snr: sat_obs.get_snr(1).unwrap_or(45) as f64,
            doppler: sat_obs.get_doppler(1).unwrap_or(0.0), lam1: LIGHT_SPEED / f1, lam2: LIGHT_SPEED / f2,
            tropo_dry, map_wet: m_w, iono_delay,
            f1, f2, sat_pos_rot: sat_pos, sat_vel, sat_clock_drift: sat_drift,
            rcv_pos_ecef, pcv_correction,
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

            let l_meas = if sat.is_iono_free && sat.cp2.is_some() {
                crate::combinations::iono_free((cp1_cyc + windup) * sat.lam1, (sat.cp2.unwrap() + windup) * sat.lam2, sat.f1, sat.f2)
            } else {
                (cp1_cyc + windup) * sat.lam1
            };
            
            let mut slip = false;
            
            let prev = *state.locktimes.get(&(sat.sat_obs.sat, 1)).unwrap_or(&0);
            let mut new_lk = prev.saturating_add(1);

            for obs in &sat.sat_obs.observations {
                if obs.code.obs_type == gneiss_core::obs::ObsType::CarrierPhase {
                    if let Some(lli) = obs.lli {
                        if (lli & 1) != 0 || (lli & 2) != 0 { slip = true; new_lk = 0; }
                    } else if let Some(lk) = obs.lock_time {
                        if lk == 0 || lk < prev { slip = true; new_lk = lk; } else { new_lk = new_lk.min(lk); }
                    }
                }
            }
            if new_lk == 0 && state.locktimes.contains_key(&(sat.sat_obs.sat, 1)) {
                slip = true;
            }
            state.locktimes.insert((sat.sat_obs.sat, 1), new_lk);

            if slip { 
                tracing::error!("Cycle slip detected for {} (gap > 1 or lli)", sat.sat_obs.sat);
                state.remove_ambiguity(sat.sat_obs.sat, 0); 
            }

            let isb = match sat.sat_obs.sat.constellation {
                gneiss_core::sat::Constellation::Glonass => state.isb_glo,
                gneiss_core::sat::Constellation::Galileo => state.isb_gal,
                gneiss_core::sat::Constellation::Beidou => state.isb_bds,
                _ => 0.0,
            };

            let expected_p = sat.dist + state.rcv_clk_bias + isb - sat.dt_sat_m + sat.tropo_dry + state.zwd * sat.map_wet;
            let expected_with_iono = if sat.is_iono_free && sat.cp2.is_some() {
                expected_p
            } else {
                expected_p - sat.iono_delay
            };
            if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 0)) {
                state.add_ambiguity(sat.sat_obs.sat, 0, l_meas - expected_with_iono, 10000.0);
            }
            state.last_observed.insert((sat.sat_obs.sat, 0), state.epoch_count as u32);
        }
    }
}
