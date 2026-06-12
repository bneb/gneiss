use crate::filter::DdObservation;
use gneiss_core::obs::{EpochObs, ObsType};

pub fn match_observations(rover_obs: &EpochObs, base_obs: &EpochObs, ephemerides: &[gneiss_core::ephemeris::Ephemeris]) -> Vec<(DdObservation, DdObservation)> {
    let mut matched_obs = Vec::new();
    for r_sat in &rover_obs.satellites {
        if !ephemerides.iter().any(|e| e.sat() == r_sat.sat) {
            continue;
        }
        if let Some(b_sat) = base_obs.satellites.iter().find(|s| s.sat == r_sat.sat) {
            let r_pr_l1 = r_sat.observations.iter().find(|o| o.code.obs_type == ObsType::Pseudorange && o.code.signal.freq_band == 1);
            let r_pr_l2 = r_sat.observations.iter().find(|o| o.code.obs_type == ObsType::Pseudorange && o.code.signal.freq_band == 2);
            let r_cp_l1 = r_sat.observations.iter().find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 1);
            let r_cp_l2 = r_sat.observations.iter().find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 2);
            let r_dop = r_sat.observations.iter().find(|o| o.code.obs_type == ObsType::Doppler && o.code.signal.freq_band == 1).map(|o| o.value).unwrap_or(0.0);
            let b_pr_l1 = b_sat.observations.iter().find(|o| o.code.obs_type == ObsType::Pseudorange && o.code.signal.freq_band == 1);
            let b_pr_l2 = b_sat.observations.iter().find(|o| o.code.obs_type == ObsType::Pseudorange && o.code.signal.freq_band == 2);
            let b_cp_l1 = b_sat.observations.iter().find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 1);
            let b_cp_l2 = b_sat.observations.iter().find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 2);
            let b_dop = b_sat.observations.iter().find(|o| o.code.obs_type == ObsType::Doppler && o.code.signal.freq_band == 1).map(|o| o.value).unwrap_or(0.0);
            let r_snr = r_sat.observations.iter().find(|o| o.code.obs_type == ObsType::Snr && o.code.signal.freq_band == 1).map(|o| o.value).unwrap_or(25.0);
            let r_lock = r_sat.observations.iter().find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 1).and_then(|o| o.lock_time);

            if r_pr_l1.is_none() { tracing::debug!("Sat {:?} missing rover PR1. Rover Obs: {:?}", r_sat.sat, r_sat.observations.iter().map(|o| o.code).collect::<Vec<_>>()); }
            if b_pr_l1.is_none() { tracing::debug!("Sat {:?} missing base PR1. Base Obs: {:?}", b_sat.sat, b_sat.observations.iter().map(|o| o.code).collect::<Vec<_>>()); }
            if r_cp_l1.is_none() { tracing::debug!("Sat {:?} missing rover CP1", r_sat.sat); }
            if b_cp_l1.is_none() { tracing::debug!("Sat {:?} missing base CP1. Base Obs: {:?}", b_sat.sat, b_sat.observations.iter().map(|o| o.code).collect::<Vec<_>>()); }

            if let (Some(r_pr1), Some(b_pr1)) = (r_pr_l1, b_pr_l1) {
                tracing::debug!("Sat {:?} L1 PR: Rover {} (val: {}), Base {} (val: {})", 
                    r_sat.sat, r_pr1.code, r_pr1.value, b_pr1.code, b_pr1.value);
                if let (Some(r_pr2), Some(b_pr2)) = (r_pr_l2, b_pr_l2) {
                    tracing::debug!("Sat {:?} L2 PR: Rover {} (val: {}), Base {} (val: {})", 
                        r_sat.sat, r_pr2.code, r_pr2.value, b_pr2.code, b_pr2.value);
                }
                matched_obs.push((
                    DdObservation { 
                        sat: r_sat.sat, 
                        pr_l1: r_pr1.value, 
                        pr_l2: r_pr_l2.map(|o| o.value), 
                        cp_l1: r_cp_l1.map(|o| o.value), 
                        cp_l2: r_cp_l2.map(|o| o.value), 
                        doppler: r_dop, 
                        snr: r_snr, 
                        locktime: r_lock 
                    },
                    DdObservation { 
                        sat: b_sat.sat, 
                        pr_l1: b_pr1.value, 
                        pr_l2: b_pr_l2.map(|o| o.value), 
                        cp_l1: b_cp_l1.map(|o| o.value), 
                        cp_l2: b_cp_l2.map(|o| o.value), 
                        doppler: b_dop, 
                        snr: 25.0, 
                        locktime: Some(1000) 
                    }
                ));
            }
        }
    }
    matched_obs
}
