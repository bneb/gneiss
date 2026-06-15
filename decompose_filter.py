import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

old_filter = """#[allow(clippy::too_many_arguments)]
fn filter_innovations_chi_squared(
    state: &mut RtkState,
    state_size: usize,
    chi_square_pr_threshold: f64,
    chi_square_cp_threshold: f64,
    z_all: &[f64],
    h_all: &[Vec<f64>],
    r_all: &[f64],
    type_all: &[(gneiss_core::sat::SatelliteId, u8, f64)],
) -> Vec<usize> {
    let mut safe_indices = Vec::new();
    
    for i in 0..z_all.len() {
        let mut h_row = DMatrix::zeros(1, state_size);
        for c in 0..state_size { h_row[(0, c)] = h_all[i][c]; }
        let s_ii = (&h_row * &state.covariance * h_row.transpose())[(0, 0)] + r_all[i];
        let chi2 = z_all[i] * z_all[i] / s_ii;
        
        let threshold = match type_all[i].1 { 
            0 => chi_square_pr_threshold * chi_square_pr_threshold,  
            1 | 2 => chi_square_cp_threshold * chi_square_cp_threshold,  
            3 => chi_square_pr_threshold * 1000.0, // Doppler relax
            _ => chi_square_pr_threshold * chi_square_pr_threshold   
        };
        
        if chi2 <= threshold { 
            safe_indices.push(i); 
            if z_all[i].abs() > 100000.0 {
                tracing::error!("MASSIVE Z PASSED PRE-FILTER! type: {}, z: {:.1}, chi2: {:.1}, thresh: {:.1}, s_ii: {:.1}, P_pos: {:.1}", type_all[i].1, z_all[i], chi2, threshold, s_ii, state.covariance[(0,0)]);
            }
            if type_all[i].1 == 1 || type_all[i].1 == 2 {
                for c in crate::filter::CORE_STATE_SIZE..state_size {
                    if h_row[(0, c)] > 0.5 {
                        let key = state.ambiguity_keys[c - crate::filter::CORE_STATE_SIZE];
                        state.reject_counts.insert(key, 0);
                    }
                }
            }
        } else {
            tracing::debug!("Rejected meas type {} with inn: {:.3}, chi2: {:.1}, threshold: {:.1}, s_ii: {:.1}", type_all[i].1, z_all[i], chi2, threshold, s_ii);
            if type_all[i].1 == 1 || type_all[i].1 == 2 {
                for c in crate::filter::CORE_STATE_SIZE..state_size {
                    if h_row[(0, c)] > 0.5 {
                        let key = state.ambiguity_keys[c - crate::filter::CORE_STATE_SIZE];
                        let count = *state.reject_counts.get(&key).unwrap_or(&0) + 1;
                        state.reject_counts.insert(key, count);
                    }
                }
            }
        }
    }
    safe_indices
}"""

new_filter = """fn update_reject_counts(state: &mut RtkState, h_row: &DMatrix<f64>, state_size: usize, passed: bool) {
    for c in crate::filter::CORE_STATE_SIZE..state_size {
        if h_row[(0, c)] > 0.5 {
            let key = state.ambiguity_keys[c - crate::filter::CORE_STATE_SIZE];
            if passed {
                state.reject_counts.insert(key, 0);
            } else {
                let count = *state.reject_counts.get(&key).unwrap_or(&0) + 1;
                state.reject_counts.insert(key, count);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn filter_innovations_chi_squared(
    state: &mut RtkState, state_size: usize, chi_pr: f64, chi_cp: f64,
    z_all: &[f64], h_all: &[Vec<f64>], r_all: &[f64], type_all: &[(gneiss_core::sat::SatelliteId, u8, f64)],
) -> Vec<usize> {
    let mut safe_indices = Vec::new();
    for i in 0..z_all.len() {
        let mut h_row = DMatrix::zeros(1, state_size);
        for c in 0..state_size { h_row[(0, c)] = h_all[i][c]; }
        let s_ii = (&h_row * &state.covariance * h_row.transpose())[(0, 0)] + r_all[i];
        let chi2 = z_all[i] * z_all[i] / s_ii;
        
        let threshold = match type_all[i].1 { 
            0 => chi_pr * chi_pr, 1 | 2 => chi_cp * chi_cp, 3 => chi_pr * 1000.0, _ => chi_pr * chi_pr   
        };
        
        let passed = chi2 <= threshold;
        if passed { safe_indices.push(i); }
        else { tracing::debug!("Rejected meas type {} with inn: {:.3}, chi2: {:.1}", type_all[i].1, z_all[i], chi2); }
        
        if type_all[i].1 == 1 || type_all[i].1 == 2 {
            update_reject_counts(state, &h_row, state_size, passed);
        }
    }
    safe_indices
}"""

content = content.replace(old_filter, new_filter)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
