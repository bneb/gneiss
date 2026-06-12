use crate::filter::{RtkState, DdObservation};
use crate::combinations::{melbourne_wubbena, lambda_ewl};
use crate::engine::updater::update;
use nalgebra::{DMatrix, DVector};

pub fn apply_tcar(
    state: &mut RtkState,
    rover_obs: &DdObservation,
    base_obs: &DdObservation,
    r_freqs: [f64; 3],
    b_freqs: [f64; 3]
) {
    let [r_f1, r_f2, r_f5] = r_freqs;
    let [_b_f1, b_f2, b_f5] = b_freqs;
    let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

    // 1. Extra-Wide-Lane (EWL) Update: L2 and L5
    if let (Some(idx_2), Some(idx_5)) = (
        state.ambiguity_keys.iter().position(|&(s, f)| s == rover_obs.sat && f == 2),
        state.ambiguity_keys.iter().position(|&(s, f)| s == rover_obs.sat && f == 5)
    ) {
        if let (Some(r_pr2), Some(r_cp2), Some(b_pr2), Some(b_cp2), Some(r_pr5), Some(r_cp5), Some(b_pr5), Some(b_cp5)) = 
            (rover_obs.pr_l2, rover_obs.cp_l2, base_obs.pr_l2, base_obs.cp_l2, rover_obs.pr_l5, rover_obs.cp_l5, base_obs.pr_l5, base_obs.cp_l5) 
        {
            let r_cp2_m = r_cp2 * (c / r_f2);
            let r_cp5_m = r_cp5 * (c / r_f5);
            let b_cp2_m = b_cp2 * (c / b_f2);
            let b_cp5_m = b_cp5 * (c / b_f5);

            let mw_ewl_r = melbourne_wubbena(r_cp2_m, r_cp5_m, r_pr2, r_pr5, r_f2, r_f5);
            let mw_ewl_b = melbourne_wubbena(b_cp2_m, b_cp5_m, b_pr2, b_pr5, b_f2, b_f5);
            let mw_ewl_sd = mw_ewl_r - mw_ewl_b;

            let lam_ewl = lambda_ewl(r_f2, r_f5);
            let n_ewl = (mw_ewl_sd / lam_ewl).round();

            if (mw_ewl_sd / lam_ewl - n_ewl).abs() < 0.25 {
                let mut h = DMatrix::zeros(1, state.covariance.nrows());
                h[(0, crate::filter::CORE_STATE_SIZE + idx_2)] = r_f2 / c;
                h[(0, crate::filter::CORE_STATE_SIZE + idx_5)] = -r_f5 / c;

                let p_ewl = (&h * &state.covariance * h.transpose())[(0,0)];
                if p_ewl > 0.05 {
                    tracing::debug!("TCAR EWL Update: Sat {:?} N_ewl = {}", rover_obs.sat, n_ewl);
                    let z = DVector::from_element(1, n_ewl);
                    let r = DMatrix::from_element(1, 1, 0.05); 
                    let _ = update(state, &z, &h, &r, 10.0, None, false);
                }
            }
        }
    }

    // 2. Wide-Lane (WL) Update: L1 and L2
    if let (Some(idx_1), Some(idx_2)) = (
        state.ambiguity_keys.iter().position(|&(s, f)| s == rover_obs.sat && f == 1),
        state.ambiguity_keys.iter().position(|&(s, f)| s == rover_obs.sat && f == 2)
    ) {
        if let (Some(&ema_wl), Some(&count)) = (state.mw_sd_ema.get(&rover_obs.sat), state.mw_sd_counts.get(&rover_obs.sat)) {
            if count >= 10 {
                let n_wl = ema_wl.round();
                if (ema_wl - n_wl).abs() < 0.15 {
                    let mut h = DMatrix::zeros(1, state.covariance.nrows());
                    h[(0, crate::filter::CORE_STATE_SIZE + idx_1)] = r_f1 / c;
                    h[(0, crate::filter::CORE_STATE_SIZE + idx_2)] = -r_f2 / c;

                    let p_wl = (&h * &state.covariance * h.transpose())[(0,0)];
                    if p_wl > 0.05 {
                        tracing::debug!("TCAR WL Update: Sat {:?} N_wl = {}", rover_obs.sat, n_wl);
                        let z = DVector::from_element(1, n_wl);
                        let r = DMatrix::from_element(1, 1, 0.01); 
                        let _ = update(state, &z, &h, &r, 10.0, None, false);
                    }
                }
            }
        }
    }
}
