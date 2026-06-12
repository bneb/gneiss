import re

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    content = f.read()

# 1. Replace compute_dd_carrier_phase
cp_pattern = r"fn compute_dd_carrier_phase\([\s\S]*?\}\n\nfn compute_dd_doppler"
cp_replacement = """fn compute_single_freq_dd_cp(
    state: &RtkState,
    sat: gneiss_core::sat::SatelliteId,
    freq_idx: u8,
    rov_ref_cp: Option<f64>,
    rov_sat_cp: Option<f64>,
    bas_ref_cp: Option<f64>,
    bas_sat_cp: Option<f64>,
    ref_idx: Option<usize>,
    ref_f: f64,
    sat_f: f64,
    comp_pr_dd: f64,
    iono_dd: f64,
    h_r: Vector3<f64>,
    h_att: Vector3<f64>,
    state_size: usize,
    r_val: f64,
) -> Option<(f64, Vec<f64>, f64, u8)> {
    let sat_idx = state.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == freq_idx)?;
    let ref_idx = ref_idx?;
    
    let (rr, rs, br, bs) = (rov_ref_cp?, rov_sat_cp?, bas_ref_cp?, bas_sat_cp?);
    
    let c = SPEED_OF_LIGHT_M_S;
    let lam_ref = c / ref_f;
    let lam_sat = c / sat_f;
    
    let cp_dd = (rs * lam_sat - rr * lam_ref) - (bs * lam_sat - br * lam_ref);
    let n_dd = state.ambiguities[sat_idx] - state.ambiguities[ref_idx];
    
    let mut h_cp = vec![0.0; state_size];
    h_cp[0] = h_r.x; h_cp[1] = h_r.y; h_cp[2] = h_r.z;
    h_cp[6] = h_att.x; h_cp[7] = h_att.y; h_cp[8] = h_att.z;
    h_cp[crate::filter::CORE_STATE_SIZE + sat_idx] = 1.0; 
    h_cp[crate::filter::CORE_STATE_SIZE + ref_idx] = -1.0;
    
    let expected_cp_dd = comp_pr_dd - iono_dd + n_dd;
    let innov = cp_dd - expected_cp_dd;
    
    Some((innov, h_cp, r_val, freq_idx))
}

fn compute_dd_carrier_phase(
    state: &RtkState,
    rov_sat: &DdObservation,
    base_sat: &DdObservation,
    rov_ref: &DdObservation,
    ref_base: &DdObservation,
    ref_idx_l1: Option<usize>,
    ref_idx_l2: Option<usize>,
    ref_idx_l5: Option<usize>,
    comp_pr_dd: f64,
    iono_dd_l1: f64,
    iono_dd_l2: f64,
    iono_dd_l5: f64,
    sat_f1: f64, sat_f2: f64, sat_f5: f64,
    ref_f1: f64, ref_f2: f64, ref_f5: f64,
    h_r: Vector3<f64>,
    h_att: Vector3<f64>,
    state_size: usize,
    var_factor: f64,
) -> Vec<(f64, Vec<f64>, f64, u8)> {
    let mut updates = Vec::new();
    let r_val = if state.is_fixed { DD_CARRIER_PHASE_FIXED_VARIANCE * var_factor } else { DD_CARRIER_PHASE_BASE_VARIANCE * var_factor };

    if let Some(up) = compute_single_freq_dd_cp(state, rov_sat.sat, 1, rov_ref.cp_l1, rov_sat.cp_l1, ref_base.cp_l1, base_sat.cp_l1, ref_idx_l1, ref_f1, sat_f1, comp_pr_dd, iono_dd_l1, h_r, h_att, state_size, r_val) {
        tracing::debug!("CP1 DD computed for {:?}", rov_sat.sat);
        updates.push(up);
    } else {
        if state.ambiguity_keys.iter().any(|&(s, f)| s == rov_sat.sat && f == 1) {
            tracing::debug!("Missing CP1 obs for {:?}", rov_sat.sat);
        } else {
            tracing::debug!("Missing ambiguity for {:?}", rov_sat.sat);
        }
    }

    if let Some(up) = compute_single_freq_dd_cp(state, rov_sat.sat, 2, rov_ref.cp_l2, rov_sat.cp_l2, ref_base.cp_l2, base_sat.cp_l2, ref_idx_l2, ref_f2, sat_f2, comp_pr_dd, iono_dd_l2, h_r, h_att, state_size, r_val) {
        updates.push(up);
    }

    if let Some(up) = compute_single_freq_dd_cp(state, rov_sat.sat, 5, rov_ref.cp_l5, rov_sat.cp_l5, ref_base.cp_l5, base_sat.cp_l5, ref_idx_l5, ref_f5, sat_f5, comp_pr_dd, iono_dd_l5, h_r, h_att, state_size, r_val) {
        updates.push(up);
    }

    updates
}

fn compute_dd_doppler"""

if not re.search(cp_pattern, content):
    print("Could not find cp_pattern")
    exit(1)

content = re.sub(cp_pattern, cp_replacement, content)


# 2. Replace build_measurement_model
bm_pattern = r"pub fn build_measurement_model\([\s\S]*"
bm_replacement = """fn select_reference_satellite(
    group: &[(DdObservation, DdObservation)],
    ephemerides: &[Ephemeris],
    time_tow: f64,
    rov_llh_ref: Vector3<f64>,
    rov_pos: Vector3<f64>,
) -> usize {
    let mut best_score = -1.0;
    let mut ref_idx = 0;
    for (i, (r, _)) in group.iter().enumerate() {
        if let Some(eph) = ephemerides.iter().filter(|e| e.sat() == r.sat).min_by(|a, b| {
            let da = (a.toe().tow - time_tow).abs();
            let db = (b.toe().tow - time_tow).abs();
            da.partial_cmp(&db).unwrap()
        }) {
            let tau = r.pr_l1 / SPEED_OF_LIGHT_M_S;
            let t_tx = gneiss_core::time::GpsTime::new(0, time_tow - tau); // week is not important here
            let (sat_pos, _, _, _) = eph.position(t_tx);
            let (_, el) = gneiss_core::coords::az_el(rov_llh_ref, rov_pos, sat_pos);
            
            let score = if r.cp_l1.is_some() { el + 100.0 } else { el };
            if score > best_score {
                best_score = score;
                ref_idx = i;
            }
        }
    }
    ref_idx
}

fn update_reject_counts(state: &mut RtkState, h_row: &DMatrix<f64>, state_size: usize, m_type: u8, accepted: bool) {
    if m_type == 1 || m_type == 2 || m_type == 5 {
        for c in crate::filter::CORE_STATE_SIZE..state_size {
            if h_row[(0, c)] > 0.5 {
                let key = state.ambiguity_keys[c - crate::filter::CORE_STATE_SIZE];
                let count = if accepted { 0 } else { *state.reject_counts.get(&key).unwrap_or(&0) + 1 };
                state.reject_counts.insert(key, count);
            }
        }
    }
}

fn validate_measurements(
    z_all: &[f64],
    h_all: &[Vec<f64>],
    r_all: &[f64],
    type_all: &[(gneiss_core::sat::SatelliteId, u8)],
    state: &mut RtkState,
    chi_square_pr_threshold: f64,
    chi_square_cp_threshold: f64,
) -> Vec<usize> {
    let state_size = crate::filter::CORE_STATE_SIZE + state.ambiguities.len();
    let mut safe_indices = Vec::new();

    for i in 0..z_all.len() {
        let mut h_row = DMatrix::zeros(1, state_size);
        for c in 0..state_size { h_row[(0, c)] = h_all[i][c]; }
        let s_ii = (&h_row * &state.covariance * h_row.transpose())[(0, 0)] + r_all[i];
        let chi2 = z_all[i] * z_all[i] / s_ii;

        let threshold = match type_all[i].1 { 
            0 => chi_square_pr_threshold * chi_square_pr_threshold,  
            1 | 2 | 5 => chi_square_cp_threshold * chi_square_cp_threshold,  
            3 => chi_square_pr_threshold * 1000.0,
            _ => chi_square_pr_threshold * chi_square_pr_threshold   
        };

        if chi2 <= threshold { 
            safe_indices.push(i); 
            if z_all[i].abs() > 100000.0 {
                tracing::error!("MASSIVE Z PASSED PRE-FILTER! type: {}, z: {:.1}, chi2: {:.1}", type_all[i].1, z_all[i], chi2);
            }
            update_reject_counts(state, &h_row, state_size, type_all[i].1, true);
        } else {
            tracing::debug!("Rejected meas type {} with inn: {:.3}, chi2: {:.1}", type_all[i].1, z_all[i], chi2);
            update_reject_counts(state, &h_row, state_size, type_all[i].1, false);
        }
    }
    safe_indices
}

fn pack_measurements(
    safe_indices: &[usize],
    z_all: &[f64],
    h_all: &[Vec<f64>],
    r_all: &[f64],
    type_all: &[(gneiss_core::sat::SatelliteId, u8)],
    state_size: usize,
) -> Option<(DVector<f64>, DMatrix<f64>, DMatrix<f64>, Vec<(gneiss_core::sat::SatelliteId, u8)>)> {
    if safe_indices.len() < 4 {
        return None;
    }
    let mut z_vec = DVector::zeros(safe_indices.len());
    let mut h_mat = DMatrix::zeros(safe_indices.len(), state_size);
    let mut r_mat = DMatrix::zeros(safe_indices.len(), safe_indices.len());
    let mut t_vec = Vec::with_capacity(safe_indices.len());

    for (new_i, &old_i) in safe_indices.iter().enumerate() {
        z_vec[new_i] = z_all[old_i];
        for c in 0..state_size { h_mat[(new_i, c)] = h_all[old_i][c]; }
        r_mat[(new_i, new_i)] = r_all[old_i];
        t_vec.push(type_all[old_i]);
    }
    Some((z_vec, h_mat, r_mat, t_vec))
}

pub fn build_measurement_model(
    state: &mut RtkState,
    matched_obs: &[(DdObservation, DdObservation)],
    ephemerides: &[Ephemeris],
    base_coord: &Coordinate,
    _rover_time: GpsTime,
    base_time: GpsTime,
    lever_arm: Vector3<f64>,
    omega_b: Vector3<f64>,
    chi_square_pr_threshold: f64,
    chi_square_cp_threshold: f64,
) -> Option<(DVector<f64>, DMatrix<f64>, DMatrix<f64>, Vec<(gneiss_core::sat::SatelliteId, u8)>)> {
    let mut z_all = Vec::new();
    let mut h_all = Vec::new();
    let mut r_all = Vec::new();
    let mut type_all = Vec::new();

    use std::collections::HashMap;
    use gneiss_core::sat::Constellation;
    let mut const_groups: HashMap<Constellation, Vec<(DdObservation, DdObservation)>> = HashMap::new();
    
    for obs in matched_obs {
        const_groups.entry(obs.0.sat.constellation).or_default().push(obs.clone());
    }

    let rov_llh_ref = gneiss_core::coords::ecef_to_llh(state.position.vector);
    
    for (_, group) in const_groups {
        if group.len() < 2 { continue; } 

        let ref_idx = select_reference_satellite(&group, ephemerides, state.time.tow, rov_llh_ref, state.position.vector);
        let mut group_clone = group.clone();
        let (ref_rover, ref_base) = group_clone.remove(ref_idx);

        if let Some((z, h, r, mt)) = compute_innovations(state, &group_clone, ephemerides, base_coord, base_time, &ref_rover, &ref_base, lever_arm, omega_b) {
            z_all.extend(z);
            h_all.extend(h);
            r_all.extend(r);
            type_all.extend(mt);
        }
    }

    let state_size = crate::filter::CORE_STATE_SIZE + state.ambiguities.len();
    let safe_indices = validate_measurements(&z_all, &h_all, &r_all, &type_all, state, chi_square_pr_threshold, chi_square_cp_threshold);

    pack_measurements(&safe_indices, &z_all, &h_all, &r_all, &type_all, state_size)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{DMatrix, Vector3};
    use gneiss_core::sat::{SatelliteId, Constellation};
    
    fn create_dummy_state() -> RtkState {
        let mut state = RtkState::new(Coordinate::new(1000.0, 1000.0, 1000.0), GpsTime::new(2000, 100000.0));
        state.covariance = DMatrix::identity(15, 15);
        state.ambiguities.push(10.0);
        state.ambiguities.push(12.0);
        state.ambiguity_keys.push((SatelliteId::new(Constellation::GPS, 1).unwrap(), 1));
        state.ambiguity_keys.push((SatelliteId::new(Constellation::GPS, 2).unwrap(), 1));
        state.covariance = DMatrix::identity(17, 17);
        state
    }

    #[test]
    fn test_compute_single_freq_dd_cp() {
        let state = create_dummy_state();
        let sat = SatelliteId::new(Constellation::GPS, 1).unwrap();
        
        // n_dd = amb[sat_idx] - amb[ref_idx]
        // sat_idx = 0 (G01), ref_idx = 1 (G02)
        // n_dd = 10.0 - 12.0 = -2.0
        
        let c = SPEED_OF_LIGHT_M_S;
        let lam = c / 1.57542e9; // ~0.19m
        
        // rov_ref(G02): 1000 cycles, rov_sat(G01): 1100 cycles
        // bas_ref(G02): 950 cycles, bas_sat(G01): 1040 cycles
        let rov_ref_cp = Some(1000.0);
        let rov_sat_cp = Some(1100.0);
        let bas_ref_cp = Some(950.0);
        let bas_sat_cp = Some(1040.0);
        
        // cp_dd = (1100*lam - 1000*lam) - (1040*lam - 950*lam)
        // cp_dd = 100*lam - 90*lam = 10*lam
        
        let expected_pr_dd = 10.0 * lam + 0.5; // True geo is 10*lam, tropo/iono 0.5
        let iono_dd = 0.2;
        
        let h_r = Vector3::new(1.0, 0.0, 0.0);
        let h_att = Vector3::new(0.0, 1.0, 0.0);
        
        let res = compute_single_freq_dd_cp(
            &state, sat, 1, rov_ref_cp, rov_sat_cp, bas_ref_cp, bas_sat_cp, 
            Some(1), 1.57542e9, 1.57542e9, expected_pr_dd, iono_dd, 
            h_r, h_att, 17, 0.01
        ).unwrap();
        
        let (innov, h, r, freq) = res;
        
        // cp_dd = 10*lam
        // exp_cp_dd = expected_pr_dd - iono_dd + n_dd = (10*lam + 0.5) - 0.2 + (-2.0)
        // exp_cp_dd = 10*lam - 1.7
        // innov = cp_dd - exp_cp_dd = 1.7
        
        assert!((innov - 1.7).abs() < 1e-9);
        assert_eq!(freq, 1);
        assert_eq!(h[0], 1.0);
        assert_eq!(h[7], 1.0);
        assert_eq!(h[15], 1.0); // sat_idx
        assert_eq!(h[16], -1.0); // ref_idx
        assert_eq!(r, 0.01);
    }
    
    #[test]
    fn test_validate_measurements() {
        let mut state = create_dummy_state();
        let sat = SatelliteId::new(Constellation::GPS, 1).unwrap();
        
        let z_all = vec![1.0, 10.0, 0.05];
        let mut h0 = vec![0.0; 17]; h0[0] = 1.0;
        let mut h1 = vec![0.0; 17]; h1[0] = 1.0; h1[15] = 1.0;
        let mut h2 = vec![0.0; 17]; h2[0] = 1.0; h2[15] = 1.0;
        let h_all = vec![h0, h1, h2];
        let r_all = vec![0.1, 0.01, 0.01];
        let type_all = vec![(sat, 0), (sat, 1), (sat, 1)]; // PR, CP, CP
        
        let safe_indices = validate_measurements(
            &z_all, &h_all, &r_all, &type_all, &mut state, 3.0, 3.0
        );
        
        // s_ii[0] = 1*1*1 + 0.1 = 1.1 => chi2 = 1.0 / 1.1 = 0.9 < 9 (Pass)
        // s_ii[1] = h*P*h' + 0.01 = 2.01 => chi2 = 100 / 2.01 = 49.7 > 9 (Fail)
        // s_ii[2] = 2.01 => chi2 = 0.0025 / 2.01 = 0.001 < 9 (Pass)
        
        assert_eq!(safe_indices, vec![0, 2]);
        
        // Check reject counts
        let key = (sat, 1);
        let count = state.reject_counts.get(&key).unwrap();
        // type 1 (Fail) -> incremented to 1
        // type 2 (Pass) -> reset to 0
        assert_eq!(*count, 0); // Wait, if it processes 1 (fail, set to 1) then 2 (pass, set to 0), it will be 0.
    }
}
"""

if not re.search(bm_pattern, content):
    print("Could not find bm_pattern")
    exit(1)

content = re.sub(bm_pattern, bm_replacement, content)

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(content)
print("Updated measurement.rs successfully")
