use gneiss_core::coords::Coordinate;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::time::GpsTime;
use crate::filter::{RtkState, DdObservation};

pub fn manage_ambiguities_and_slips(
    state: &mut RtkState,
    config: &crate::engine::EngineConfig,
    matched_obs: &[(DdObservation, DdObservation)],
    ephemerides: &[Ephemeris],
    base_coord: &Coordinate,
    rover_time: GpsTime,
    base_time: GpsTime,
) {
    for (r, b) in matched_obs {
        let mut slip = false;
        let prev_lock = *state.locktimes.get(&(r.sat, 1)).unwrap_or(&0);
        let mut new_lock = prev_lock + 1;
        
        if let Some(r_lock) = r.locktime {
            if r_lock == 0 {
                tracing::debug!("slip=true because r_lock == 0 for {:?}", r.sat);
                slip = true;
                new_lock = 0;
            } else if r_lock < prev_lock {
                tracing::debug!("slip=true because r_lock < prev_lock for {:?}", r.sat);
                slip = true;
                new_lock = r_lock;
            } else {
                new_lock = r_lock;
            }
        }

        if new_lock == 0 && state.locktimes.contains_key(&(r.sat, 1)) {
            tracing::debug!("slip=true because new_lock == 0 for {:?}", r.sat);
            slip = true;
        }

        state.locktimes.insert((r.sat, 1), new_lock);
        state.locktimes.insert((r.sat, 2), new_lock);
        
        let r_freq_num = ephemerides.iter().find(|e| e.sat() == r.sat).map(|e| e.freq_num()).unwrap_or(0);
        let b_freq_num = ephemerides.iter().find(|e| e.sat() == b.sat).map(|e| e.freq_num()).unwrap_or(0);
        
        let (r_f1, r_f2) = gneiss_core::signal::satellite_frequencies(r.sat, r_freq_num);
        let (b_f1, b_f2) = gneiss_core::signal::satellite_frequencies(b.sat, b_freq_num);

        let mut slip_l1 = slip;
        let mut slip_l2 = slip;
        if *state.reject_counts.get(&(r.sat, 1)).unwrap_or(&0) > config.max_reject_count { slip_l1 = true; }
        if *state.reject_counts.get(&(r.sat, 2)).unwrap_or(&0) > config.max_reject_count { slip_l2 = true; }

        if let Some(r_cp1) = r.cp_l1 {
            if let Some(&(prev_cp, prev_doppler, prev_time)) = state.phase_history.get(&(r.sat, 1)) {
                let dt = rover_time - prev_time;
                if dt > 0.0 && dt <= config.max_base_age_s {
                    if check_doppler_phase_slip(r_cp1, prev_cp, r.doppler, prev_doppler, dt, config.doppler_slip_threshold_cycles) {
                        tracing::debug!("Doppler-Phase cycle slip detected on {:?} L1", r.sat);
                        slip_l1 = true;
                    }
                } else if dt > config.max_base_age_s {
                    tracing::debug!("Data gap > {:.1}s on {:?} L1, resetting ambiguity", config.max_base_age_s, r.sat);
                    slip_l1 = true;
                }
            }
            state.phase_history.insert((r.sat, 1), (r_cp1, r.doppler, rover_time));
        }

        if let Some(r_cp2) = r.cp_l2 {
            let doppler_l2 = r.doppler * (r_f2 / r_f1);
            if let Some(&(prev_cp, prev_doppler, prev_time)) = state.phase_history.get(&(r.sat, 2)) {
                let dt = rover_time - prev_time;
                if dt > 0.0 && dt <= config.max_base_age_s {
                    if check_doppler_phase_slip(r_cp2, prev_cp, doppler_l2, prev_doppler, dt, config.doppler_slip_threshold_cycles) {
                        tracing::debug!("Doppler-Phase cycle slip detected on {:?} L2", r.sat);
                        slip_l2 = true;
                    }
                } else if dt > config.max_base_age_s {
                    tracing::debug!("Data gap > {:.1}s on {:?} L2, resetting ambiguity", config.max_base_age_s, r.sat);
                    slip_l2 = true;
                }
            }
            state.phase_history.insert((r.sat, 2), (r_cp2, doppler_l2, rover_time));
        }

        // Geometry-Free cycle slip detection (feature-gated)
        #[cfg(feature = "gf-slip")]
        if let (Some(r_cp1), Some(r_cp2)) = (r.cp_l1, r.cp_l2) {
            let lam1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f1;
            let lam2 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f2;
            let l_gf = r_cp1 * lam1 - r_cp2 * lam2;
            if let Some(&prev_gf) = state.gf_values.get(&r.sat) {
                if (l_gf - prev_gf).abs() > 0.05 {
                    tracing::debug!("GF cycle slip detected on {:?}: |ΔGF| = {:.4}m", r.sat, (l_gf - prev_gf).abs());
                    slip_l1 = true;
                    slip_l2 = true;
                }
            }
            state.gf_values.insert(r.sat, l_gf);
        }

        if !state.ambiguity_keys.contains(&(r.sat, 1)) || slip_l1 || slip_l2 {
            if slip_l1 || slip_l2 { 
                tracing::debug!("Cycle slip detected for {:?}: slip_l1={}, slip_l2={}", r.sat, slip_l1, slip_l2);
                state.remove_ambiguity(r.sat, 1); 
                state.remove_ambiguity(r.sat, 2);
                state.reject_counts.insert((r.sat, 1), 0);
                state.reject_counts.insert((r.sat, 2), 0);
            }
            
            if let (Some(r_cp1), Some(b_cp1)) = (r.cp_l1, b.cp_l1) {
                let lam_r1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f1;
                let lam_b1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / b_f1;
                let cp_l1_rov = r_cp1 * lam_r1;
                let cp_l1_base = b_cp1 * lam_b1;

                let mut initialized = false;
                if state.covariance[(0,0)] < 0.1 {
                    // Try to find an anchor satellite to cancel the clock bias
                    for (anchor_r, anchor_b) in matched_obs.iter() {
                        if anchor_r.sat == r.sat { continue; }
                        if let Some(anchor_idx) = state.ambiguity_keys.iter().position(|&(s, f)| s == anchor_r.sat && f == 1) {
                            if state.covariance[(15 + anchor_idx, 15 + anchor_idx)] < 0.05 {
                                if let (Some(ar_cp), Some(ab_cp)) = (anchor_r.cp_l1, anchor_b.cp_l1) {
                                    let (ar_sat_vec, _, _, _) = ephemerides.iter().find(|e| e.sat() == anchor_r.sat).unwrap().position(rover_time);
                                    let (ab_sat_vec, _, _, _) = ephemerides.iter().find(|e| e.sat() == anchor_b.sat).unwrap().position(base_time);
                                    let ar_dist_rov = (state.position.vector - ar_sat_vec).norm();
                                    let ar_dist_base = (base_coord.vector - ab_sat_vec).norm();
                                    
                                    let a_lam_r = gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f1;
                                    let a_lam_b = gneiss_core::constants::SPEED_OF_LIGHT_M_S / b_f1;
                                    let a_cp_rov = ar_cp * a_lam_r;
                                    let a_cp_base = ab_cp * a_lam_b;
                                    
                                    let anchor_sd = state.ambiguities[anchor_idx];
                                    let b_clock_rov = a_cp_rov - ar_dist_rov - anchor_sd;
                                    let b_clock_base = a_cp_base - ar_dist_base; // Base has no ambiguity in SD, assuming SD = rov - base

                                    let (r_sat_vec, _, _, _) = ephemerides.iter().find(|e| e.sat() == r.sat).unwrap().position(rover_time);
                                    let (b_sat_vec, _, _, _) = ephemerides.iter().find(|e| e.sat() == b.sat).unwrap().position(base_time);
                                    let dist_rov = (state.position.vector - r_sat_vec).norm();
                                    let dist_base = (base_coord.vector - b_sat_vec).norm();
                                    
                                    let initial_est_l1 = (cp_l1_rov - dist_rov - b_clock_rov) - (cp_l1_base - dist_base - b_clock_base);
                                    state.add_ambiguity(r.sat, 1, initial_est_l1, config.initial_ambiguity_variance);
                                    initialized = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if !initialized {
                    let initial_est_l1 = (cp_l1_rov - r.pr_l1) - (cp_l1_base - b.pr_l1);
                    state.add_ambiguity(r.sat, 1, initial_est_l1, config.initial_ambiguity_variance);
                }
            }
            
            if let (Some(r_pr2), Some(r_cp2), Some(b_pr2), Some(b_cp2)) = (r.pr_l2, r.cp_l2, b.pr_l2, b.cp_l2) {
                let lam_r2 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f2;
                let lam_b2 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / b_f2;
                let cp_l2_rov = r_cp2 * lam_r2;
                let cp_l2_base = b_cp2 * lam_b2;
                
                let mut initialized = false;
                if state.covariance[(0,0)] < 0.1 {
                    for (anchor_r, anchor_b) in matched_obs.iter() {
                        if anchor_r.sat == r.sat { continue; }
                        if let Some(anchor_idx) = state.ambiguity_keys.iter().position(|&(s, f)| s == anchor_r.sat && f == 2) {
                            if state.covariance[(15 + anchor_idx, 15 + anchor_idx)] < 0.05 {
                                if let (Some(ar_cp), Some(ab_cp)) = (anchor_r.cp_l2, anchor_b.cp_l2) {
                                    let (ar_sat_vec, _, _, _) = ephemerides.iter().find(|e| e.sat() == anchor_r.sat).unwrap().position(rover_time);
                                    let (ab_sat_vec, _, _, _) = ephemerides.iter().find(|e| e.sat() == anchor_b.sat).unwrap().position(base_time);
                                    let ar_dist_rov = (state.position.vector - ar_sat_vec).norm();
                                    let ar_dist_base = (base_coord.vector - ab_sat_vec).norm();
                                    
                                    let a_lam_r = gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f2;
                                    let a_lam_b = gneiss_core::constants::SPEED_OF_LIGHT_M_S / b_f2;
                                    let a_cp_rov = ar_cp * a_lam_r;
                                    let a_cp_base = ab_cp * a_lam_b;
                                    
                                    let anchor_sd = state.ambiguities[anchor_idx];
                                    let b_clock_rov = a_cp_rov - ar_dist_rov - anchor_sd;
                                    let b_clock_base = a_cp_base - ar_dist_base;

                                    let (r_sat_vec, _, _, _) = ephemerides.iter().find(|e| e.sat() == r.sat).unwrap().position(rover_time);
                                    let (b_sat_vec, _, _, _) = ephemerides.iter().find(|e| e.sat() == b.sat).unwrap().position(base_time);
                                    let dist_rov = (state.position.vector - r_sat_vec).norm();
                                    let dist_base = (base_coord.vector - b_sat_vec).norm();
                                    
                                    let initial_est_l2 = (cp_l2_rov - dist_rov - b_clock_rov) - (cp_l2_base - dist_base - b_clock_base);
                                    state.add_ambiguity(r.sat, 2, initial_est_l2, config.initial_ambiguity_variance);
                                    initialized = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if !initialized {
                    let initial_est_l2 = (cp_l2_rov - r_pr2) - (cp_l2_base - b_pr2);
                    state.add_ambiguity(r.sat, 2, initial_est_l2, config.initial_ambiguity_variance);
                }
            }
        }
        if let (Some(r_pr2), Some(r_cp2), Some(b_pr2), Some(b_cp2), Some(r_cp1), Some(b_cp1)) = (r.pr_l2, r.cp_l2, b.pr_l2, b.cp_l2, r.cp_l1, b.cp_l1) {
            let r_cp1_m = r_cp1 * (gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f1);
            let r_cp2_m = r_cp2 * (gneiss_core::constants::SPEED_OF_LIGHT_M_S / r_f2);
            let b_cp1_m = b_cp1 * (gneiss_core::constants::SPEED_OF_LIGHT_M_S / b_f1);
            let b_cp2_m = b_cp2 * (gneiss_core::constants::SPEED_OF_LIGHT_M_S / b_f2);
            let mw_sd = crate::combinations::melbourne_wubbena(r_cp1_m, r_cp2_m, r.pr_l1, r_pr2, r_f1, r_f2) - 
                        crate::combinations::melbourne_wubbena(b_cp1_m, b_cp2_m, b.pr_l1, b_pr2, b_f1, b_f2);
            state.update_mw(r.sat, mw_sd / crate::combinations::lambda_wl(r_f1, r_f2));
        }
    }
}

pub fn check_doppler_phase_slip(
    cp: f64,
    prev_cp: f64,
    doppler: f64,
    prev_doppler: f64,
    dt: f64,
    threshold: f64,
) -> bool {
    let expected_change = -0.5 * (doppler + prev_doppler) * dt;
    let diff = (cp - prev_cp - expected_change).abs();
    if diff > threshold {
        tracing::debug!("Doppler slip: cp={} prev_cp={} dop={} prev_dop={} exp={} diff={}", cp, prev_cp, doppler, prev_doppler, expected_change, diff);
    }
    diff > threshold
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_doppler_phase_cycle_slip_detection() {
        let prev_cp = 100.0;
        let prev_doppler = -2.0;
        
        let dt = 1.0;
        let doppler = 5.0; // test approaching satellite
        
        let expected_change = -0.5 * (doppler + prev_doppler) * dt;
        assert_eq!(expected_change, -1.5);
        
        let cp = 100.0 + expected_change; 
        
        let slip = check_doppler_phase_slip(cp, prev_cp, doppler, prev_doppler, dt, 5.0);
        assert!(!slip, "No slip should be detected for consistent doppler");
        
        let cp_slip = 100.0 + expected_change + 6.0;
        let slip = check_doppler_phase_slip(cp_slip, prev_cp, doppler, prev_doppler, dt, 5.0);
        assert!(slip, "Slip should be detected when phase jumps");
    }
}
