import re

with open("crates/gneiss-rtk/src/filter.rs", "r") as f:
    content = f.read()

replacement = """pub fn select_ar_candidates(
    state: &RtkState,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    ar_min_lock: u32,
) -> Vec<(usize, usize, u16, f64)> {
    let constellations = [gneiss_core::sat::Constellation::Gps, gneiss_core::sat::Constellation::Galileo];
    let candidates = filter_by_locktime(state, &constellations, ar_min_lock);
    let mut candidate_vars = compute_candidate_variance(state, ephemerides, &candidates);
    candidate_vars.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));
    candidate_vars
}

fn filter_by_locktime(
    state: &RtkState,
    constellations: &[gneiss_core::sat::Constellation],
    ar_min_lock: u32,
) -> Vec<(usize, usize, u16)> {
    let num_amb = state.ambiguities.len();
    let mut candidates = Vec::new();
    for &constell in constellations {
        let mut best_ref_idx = None;
        let mut max_lock = 0;
        for i in 0..num_amb {
            let (sat, freq) = state.ambiguity_keys[i];
            if sat.constellation != constell || freq != 1 { continue; }
            let lock = *state.locktimes.get(&(sat, freq)).unwrap_or(&0);
            if lock >= ar_min_lock as u16 && lock > max_lock {
                max_lock = lock; best_ref_idx = Some(i);
            }
        }
        if let Some(ref_idx) = best_ref_idx {
            let ref_sat_id = state.ambiguity_keys[ref_idx].0;
            let l2_ref_idx = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_sat_id && f == 2);
            for i in 0..num_amb {
                if i == ref_idx || Some(i) == l2_ref_idx { continue; }
                let (rov_sat, freq) = state.ambiguity_keys[i];
                if rov_sat.constellation != constell { continue; }
                let lock = *state.locktimes.get(&(rov_sat, freq)).unwrap_or(&0);
                if lock >= ar_min_lock as u16 {
                    if freq == 1 { candidates.push((i, ref_idx, lock)); }
                    else if let Some(r2_idx) = l2_ref_idx { candidates.push((i, r2_idx, lock)); }
                }
            }
        }
    }
    candidates
}

fn compute_candidate_variance(
    state: &RtkState,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    candidates: &[(usize, usize, u16)],
) -> Vec<(usize, usize, u16, f64)> {
    let num_amb = state.ambiguities.len();
    let q_sd = state.covariance.view((CORE_STATE_SIZE, CORE_STATE_SIZE), (num_amb, num_amb));
    let mut candidate_vars = Vec::new();
    for &(rov, r_idx, lock) in candidates {
        let (rov_sat_id, freq_band) = state.ambiguity_keys[rov];
        if rov_sat_id.constellation == gneiss_core::sat::Constellation::Glonass { continue; }
        let freq_num = ephemerides.iter().find(|e| e.sat() == rov_sat_id).map(|e| e.freq_num()).unwrap_or(0);
        let (f1, f2) = gneiss_core::signal::satellite_frequencies(rov_sat_id, freq_num);
        let lam = gneiss_core::constants::SPEED_OF_LIGHT_M_S / if freq_band == 1 { f1 } else { f2 };
        let q_dd = q_sd[(rov, rov)] + q_sd[(r_idx, r_idx)] - 2.0 * q_sd[(rov, r_idx)];
        let var_cycles = q_dd / (lam * lam);
        if var_cycles < 10000.0 { candidate_vars.push((rov, r_idx, lock, var_cycles)); }
    }
    candidate_vars
}

pub fn build_lambda_matrices(
    state: &RtkState,
    candidates: &[(usize, usize, u16, f64)],
    subset_size: usize,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
) -> (DMatrix<f64>, DVector<f64>, DMatrix<f64>) {
    let (d_mat, a_cycles) = build_lambda_design_matrix(state, candidates, subset_size, ephemerides);
    let q_cycles = build_lambda_variance_matrix(state, candidates, subset_size, ephemerides);
    (d_mat, a_cycles, q_cycles)
}

fn build_lambda_design_matrix(
    state: &RtkState,
    candidates: &[(usize, usize, u16, f64)],
    subset_size: usize,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
) -> (DMatrix<f64>, DVector<f64>) {
    let num_amb = state.ambiguities.len();
    let a_sd = nalgebra::DVector::from_vec(state.ambiguities.clone());
    let mut d_mat = DMatrix::zeros(subset_size, num_amb);
    let mut a_cycles = DVector::zeros(subset_size);
    for row in 0..subset_size {
        let (rov, r_idx, _, _) = candidates[row];
        let (rov_sat_id, freq_band) = state.ambiguity_keys[rov];
        let freq_num = ephemerides.iter().find(|e| e.sat() == rov_sat_id).map(|e| e.freq_num()).unwrap_or(0);
        let (f1, f2) = gneiss_core::signal::satellite_frequencies(rov_sat_id, freq_num);
        let lam = gneiss_core::constants::SPEED_OF_LIGHT_M_S / if freq_band == 1 { f1 } else { f2 };
        d_mat[(row, rov)] = 1.0;
        d_mat[(row, r_idx)] = -1.0;
        a_cycles[row] = (a_sd[rov] - a_sd[r_idx]) / lam;
    }
    (d_mat, a_cycles)
}

fn build_lambda_variance_matrix(
    state: &RtkState,
    candidates: &[(usize, usize, u16, f64)],
    subset_size: usize,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
) -> DMatrix<f64> {
    let num_amb = state.ambiguities.len();
    let q_sd = state.covariance.view((CORE_STATE_SIZE, CORE_STATE_SIZE), (num_amb, num_amb));
    let mut q_cycles = DMatrix::zeros(subset_size, subset_size);
    for r in 0..subset_size {
        for c in 0..subset_size {
            let (rov_r, ref_r, _, _) = candidates[r];
            let (rov_c, ref_c, _, _) = candidates[c];
            let freq_r = state.ambiguity_keys[rov_r].1;
            let freq_c = state.ambiguity_keys[rov_c].1;
            let (sat_r, _) = state.ambiguity_keys[rov_r];
            let freq_num_r = ephemerides.iter().find(|e| e.sat() == sat_r).map(|e| e.freq_num()).unwrap_or(0);
            let (f1_r, f2_r) = gneiss_core::signal::satellite_frequencies(sat_r, freq_num_r);
            let lam_r = gneiss_core::constants::SPEED_OF_LIGHT_M_S / if freq_r == 1 { f1_r } else { f2_r };
            let (sat_c, _) = state.ambiguity_keys[rov_c];
            let freq_num_c = ephemerides.iter().find(|e| e.sat() == sat_c).map(|e| e.freq_num()).unwrap_or(0);
            let (f1_c, f2_c) = gneiss_core::signal::satellite_frequencies(sat_c, freq_num_c);
            let lam_c = gneiss_core::constants::SPEED_OF_LIGHT_M_S / if freq_c == 1 { f1_c } else { f2_c };
            let q_dd = q_sd[(rov_r, rov_c)] - q_sd[(rov_r, ref_c)] - q_sd[(ref_r, rov_c)] + q_sd[(ref_r, ref_c)];
            q_cycles[(r, c)] = q_dd / (lam_r * lam_c);
        }
    }
    q_cycles
}"""

# Find select_ar_candidates
start_idx = content.find("pub fn select_ar_candidates")
# Find end of file since this is at the end of the file
new_content = content[:start_idx] + replacement + "\n"

with open("crates/gneiss-rtk/src/filter.rs", "w") as f:
    f.write(new_content)

print("Refactored AR candidates and lambda matrices")
