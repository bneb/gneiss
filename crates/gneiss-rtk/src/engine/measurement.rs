use nalgebra::{DMatrix, DVector, Vector3};
use gneiss_core::coords::Coordinate;
use gneiss_core::time::GpsTime;
use gneiss_core::ephemeris::Ephemeris;
use crate::filter::{RtkState, DdObservation};

pub fn snr_scale(snr: f64) -> f64 {
    let snr_clamped = snr.max(25.0).min(50.0);
    let scale = libm::pow(10.0, (45.0 - snr_clamped) / 10.0);
    if scale > 100.0 { 100.0 } else { scale }
}

pub fn compute_innovations(
    state: &mut RtkState,
    group: &[(DdObservation, DdObservation)],
    ephemerides: &[Ephemeris],
    base_coord: &Coordinate,
    base_time: GpsTime,
    ref_rover: &DdObservation,
    ref_base: &DdObservation,
    lever_arm: Vector3<f64>,
) -> Option<(Vec<f64>, Vec<Vec<f64>>, Vec<f64>, Vec<u8>)> {
    let mut z_vals = Vec::new();
    let mut h_rows = Vec::new();
    let mut r_vals = Vec::new();
    let mut meas_type = Vec::new();

    let state_size = 15 + state.ambiguities.len();
    
    let r_b_e = state.attitude.to_rotation_matrix();
    let pos_apc = state.position.vector + r_b_e * lever_arm;

    let set_rov = gneiss_core::tides::solid_earth_tides_ecef(state.time, pos_apc);
    let set_bas = gneiss_core::tides::solid_earth_tides_ecef(state.time, base_coord.vector);
    
    let pos_apc = pos_apc + set_rov;
    let base_coord_vec = base_coord.vector + set_bas;

    let get_sat_pos = |eph: &Ephemeris, pr: f64, t_rx: GpsTime, rx_pos: Vector3<f64>| -> Vector3<f64> {
        let tau_pr = pr / 299792458.0;
        let t_tx_nom = GpsTime::new(t_rx.week, t_rx.tow - tau_pr);
        let (_, _, dt_s, _) = eph.position(t_tx_nom);
        let t_tx_true = GpsTime::new(t_rx.week, t_rx.tow - tau_pr - dt_s);
        let (raw_vec, _, _, _) = eph.position(t_tx_true);
        
        let mut sat_pos = raw_vec;
        for _ in 0..2 {
            let geometric_range = (sat_pos - rx_pos).norm();
            let true_tau = geometric_range / 299792458.0;
            let theta = 7.2921151467e-5 * true_tau;
            let cos_t = f64::cos(theta);
            let sin_t = f64::sin(theta);
            sat_pos = nalgebra::Vector3::new(
                raw_vec.x * cos_t + raw_vec.y * sin_t,
                -raw_vec.x * sin_t + raw_vec.y * cos_t,
                raw_vec.z
            );
        }
        sat_pos
    };

    let find_eph = |sat| {
        ephemerides.iter()
            .filter(|e| e.sat() == sat)
            .min_by(|a, b| {
                let da = (a.toe().tow - state.time.tow).abs();
                let db = (b.toe().tow - state.time.tow).abs();
                da.partial_cmp(&db).unwrap()
            })
    };

    let ref_eph = find_eph(ref_rover.sat)?;
    
    let ref_sat_vec_rov = get_sat_pos(ref_eph, ref_rover.pr_l1, state.time, pos_apc);
    let ref_sat_vec_bas = get_sat_pos(ref_eph, ref_base.pr_l1, base_time, base_coord_vec);

    let e_ref_rov = (ref_sat_vec_rov - pos_apc).normalize();
    
    let (ref_idx_l1, ref_idx_l2) = {
        let idx1 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover.sat && f == 1)?;
        let idx2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover.sat && f == 2);
        (idx1, idx2)
    };

    let tropo_params = gneiss_core::atmosphere::TropoParams::default();
    let iono_params = gneiss_core::atmosphere::KlobucharParams::default();
    let sun_pos = gneiss_core::sun::sun_position_ecef(state.time);

    for (rover_sat, base_sat) in group {
        if let Some(sat_eph) = find_eph(rover_sat.sat) {
            
            let sat_vec_rov = get_sat_pos(sat_eph, rover_sat.pr_l1, state.time, pos_apc);
            let sat_vec_bas = get_sat_pos(sat_eph, base_sat.pr_l1, base_time, base_coord_vec);

            let e_sat_rov = (sat_vec_rov - pos_apc).normalize();
            let h_r = e_ref_rov - e_sat_rov; 
            
            let lever_ecef = r_b_e * lever_arm;
            let h_att = lever_ecef.cross(&h_r);

            let base_llh = gneiss_core::coords::ecef_to_llh(base_coord_vec);
            let rov_llh = gneiss_core::coords::ecef_to_llh(pos_apc);
            
            let (az_rov_sat, el_rov_sat) = gneiss_core::coords::az_el(rov_llh, pos_apc, sat_vec_rov);
            let (_, el_rov_ref) = gneiss_core::coords::az_el(rov_llh, pos_apc, ref_sat_vec_rov);
            let (az_bas_sat, el_bas_sat) = gneiss_core::coords::az_el(base_llh, base_coord_vec, sat_vec_bas);
            let (_, el_bas_ref) = gneiss_core::coords::az_el(base_llh, base_coord_vec, ref_sat_vec_bas);
            
            let tropo_rov_sat = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&tropo_params, rov_llh, el_rov_sat);
            let tropo_rov_ref = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&tropo_params, rov_llh, el_rov_ref);
            let tropo_bas_sat = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&tropo_params, base_llh, el_bas_sat);
            let tropo_bas_ref = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&tropo_params, base_llh, el_bas_ref);
            
            let tropo_dd = (tropo_rov_sat - tropo_rov_ref) - (tropo_bas_sat - tropo_bas_ref);

            let az_rov_ref = gneiss_core::coords::az_el(rov_llh, pos_apc, ref_sat_vec_rov).0;
            let az_bas_ref = gneiss_core::coords::az_el(base_llh, base_coord_vec, ref_sat_vec_bas).0;
            
            let iono_rov_sat = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(&iono_params, rov_llh, az_rov_sat, el_rov_sat, state.time);
            let iono_rov_ref = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(&iono_params, rov_llh, az_rov_ref, el_rov_ref, state.time);
            let iono_bas_sat = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(&iono_params, base_llh, az_bas_sat, el_bas_sat, state.time);
            let iono_bas_ref = gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(&iono_params, base_llh, az_bas_ref, el_bas_ref, state.time);
            
            let iono_dd_l1 = (iono_rov_sat - iono_rov_ref) - (iono_bas_sat - iono_bas_ref);

            let (ref_f1, ref_f2) = gneiss_core::signal::satellite_frequencies(ref_rover.sat, ref_eph.freq_num());
            let (sat_f1, sat_f2) = gneiss_core::signal::satellite_frequencies(rover_sat.sat, sat_eph.freq_num());
            
            let f_ratio_sat_l2 = (sat_f1 / sat_f2).powi(2);
            let f_ratio_ref_l2 = (ref_f1 / ref_f2).powi(2);
            let iono_dd_l2 = (iono_rov_sat * f_ratio_sat_l2 - iono_rov_ref * f_ratio_ref_l2) - (iono_bas_sat * f_ratio_sat_l2 - iono_bas_ref * f_ratio_ref_l2);

            let prev_w_sat = *state.windup.get(&rover_sat.sat).unwrap_or(&0.0);
            let prev_w_ref = *state.windup.get(&ref_rover.sat).unwrap_or(&0.0);
            
            let w_sat = gneiss_core::windup::phase_windup(sat_vec_rov, sun_pos, pos_apc, prev_w_sat);
            let w_ref = gneiss_core::windup::phase_windup(ref_sat_vec_rov, sun_pos, pos_apc, prev_w_ref);
            
            state.windup.insert(rover_sat.sat, w_sat);
            state.windup.insert(ref_rover.sat, w_ref);

            let mut rov_sat = rover_sat.clone();
            if let Some(cp) = &mut rov_sat.cp_l1 { *cp += w_sat; }
            if let Some(cp2) = &mut rov_sat.cp_l2 { *cp2 += w_sat; }
            
            let mut rov_ref = ref_rover.clone();
            if let Some(cp) = &mut rov_ref.cp_l1 { *cp += w_ref; }
            if let Some(cp2) = &mut rov_ref.cp_l2 { *cp2 += w_ref; }

            let prev_w_bas_sat = *state.windup.get(&base_sat.sat).unwrap_or(&0.0);
            let prev_w_bas_ref = *state.windup.get(&ref_base.sat).unwrap_or(&0.0);
            let w_bas_sat = gneiss_core::windup::phase_windup(sat_vec_bas, sun_pos, base_coord_vec, prev_w_bas_sat);
            let w_bas_ref = gneiss_core::windup::phase_windup(ref_sat_vec_bas, sun_pos, base_coord_vec, prev_w_bas_ref);
            state.windup.insert(base_sat.sat, w_bas_sat);
            state.windup.insert(ref_base.sat, w_bas_ref);

            let mut bas_sat = base_sat.clone();
            if let Some(cp) = &mut bas_sat.cp_l1 { *cp += w_bas_sat; }
            if let Some(cp2) = &mut bas_sat.cp_l2 { *cp2 += w_bas_sat; }

            let mut bas_ref = ref_base.clone();
            if let Some(cp) = &mut bas_ref.cp_l1 { *cp += w_bas_ref; }
            if let Some(cp2) = &mut bas_ref.cp_l2 { *cp2 += w_bas_ref; }

            let scale_sat = snr_scale(rover_sat.snr);
            let scale_ref = snr_scale(ref_rover.snr);

            let sin_el_rov_sat = f64::sin(el_rov_sat).max(0.1);
            let sin_el_rov_ref = f64::sin(el_rov_ref).max(0.1);
            let sin_el_bas_sat = f64::sin(el_bas_sat).max(0.1);
            let sin_el_bas_ref = f64::sin(el_bas_ref).max(0.1);

            let var_factor = (scale_sat / (sin_el_rov_sat * sin_el_rov_sat)) 
                           + (scale_ref / (sin_el_rov_ref * sin_el_rov_ref))
                           + (1.0 / (sin_el_bas_sat * sin_el_bas_sat)) 
                           + (1.0 / (sin_el_bas_ref * sin_el_bas_ref));

            let comp_pr_dd = ( (pos_apc - sat_vec_rov).norm() - (pos_apc - ref_sat_vec_rov).norm() ) - ( (base_coord_vec - sat_vec_bas).norm() - (base_coord_vec - ref_sat_vec_bas).norm() ) + tropo_dd;
            
            let pr_dd_l1 = (rov_sat.pr_l1 - rov_ref.pr_l1) - (base_sat.pr_l1 - ref_base.pr_l1);
            z_vals.push(pr_dd_l1 - (comp_pr_dd + iono_dd_l1));
            let mut h_pr1 = vec![0.0; state_size]; 
            h_pr1[0] = h_r.x; h_pr1[1] = h_r.y; h_pr1[2] = h_r.z;
            h_pr1[6] = h_att.x; h_pr1[7] = h_att.y; h_pr1[8] = h_att.z;
            h_rows.push(h_pr1); r_vals.push(4.0 * var_factor); meas_type.push(0); 

            if let (Some(rr2), Some(rs2), Some(br2), Some(bs2)) = (rov_ref.pr_l2, rov_sat.pr_l2, bas_ref.pr_l2, bas_sat.pr_l2) {
                let pr_dd_l2 = (rs2 - rr2) - (bs2 - br2);
                z_vals.push(pr_dd_l2 - (comp_pr_dd + iono_dd_l2));
                let mut h_pr2 = vec![0.0; state_size]; 
                h_pr2[0] = h_r.x; h_pr2[1] = h_r.y; h_pr2[2] = h_r.z;
                h_pr2[6] = h_att.x; h_pr2[7] = h_att.y; h_pr2[8] = h_att.z;
                h_rows.push(h_pr2); r_vals.push(4.0 * var_factor); meas_type.push(0);
            }

            if let Some(sat_idx) = state.ambiguity_keys.iter().position(|&(s, f)| s == rover_sat.sat && f == 1) {
                if let (Some(rr1), Some(rs1), Some(br1), Some(bs1)) = (rov_ref.cp_l1, rov_sat.cp_l1, bas_ref.cp_l1, bas_sat.cp_l1) {
                    let c = 299792458.0;
                    let lam_ref_1 = c / ref_f1;
                    let lam_sat_1 = c / sat_f1;
                    
                    let cp_l1_rov_ref = rr1 * lam_ref_1;
                    let cp_l1_rov_sat = rs1 * lam_sat_1;
                    let cp_l1_bas_ref = br1 * lam_ref_1;
                    let cp_l1_bas_sat = bs1 * lam_sat_1;
                    
                    let cp_dd_l1 = (cp_l1_rov_sat - cp_l1_rov_ref) - (cp_l1_bas_sat - cp_l1_bas_ref);
                    let n_dd_l1 = state.ambiguities[sat_idx] - state.ambiguities[ref_idx_l1];
                    
                    z_vals.push(cp_dd_l1 - (comp_pr_dd - iono_dd_l1 + n_dd_l1));
                    let mut h_cp1 = vec![0.0; state_size]; 
                    h_cp1[0] = h_r.x; h_cp1[1] = h_r.y; h_cp1[2] = h_r.z;
                    h_cp1[6] = h_att.x; h_cp1[7] = h_att.y; h_cp1[8] = h_att.z;
                    h_cp1[15 + sat_idx] = 1.0; h_cp1[15 + ref_idx_l1] = -1.0;
                    h_rows.push(h_cp1); 
                    r_vals.push(if state.is_fixed { 1e-6 * var_factor } else { 0.0001 * var_factor }); 
                    meas_type.push(1);
                }
            }

            if let (Some(sat_idx), Some(ref_idx)) = (state.ambiguity_keys.iter().position(|&(s, f)| s == rover_sat.sat && f == 2), ref_idx_l2) {
                if let (Some(rr2), Some(rs2), Some(br2), Some(bs2)) = (rov_ref.cp_l2, rov_sat.cp_l2, bas_ref.cp_l2, bas_sat.cp_l2) {
                    let c = 299792458.0;
                    let lam_ref_2 = c / ref_f2;
                    let lam_sat_2 = c / sat_f2;
                    
                    let cp_l2_rov_ref = rr2 * lam_ref_2;
                    let cp_l2_rov_sat = rs2 * lam_sat_2;
                    let cp_l2_bas_ref = br2 * lam_ref_2;
                    let cp_l2_bas_sat = bs2 * lam_sat_2;
                    
                    let cp_dd_l2 = (cp_l2_rov_sat - cp_l2_rov_ref) - (cp_l2_bas_sat - cp_l2_bas_ref);
                    let n_dd_l2 = state.ambiguities[sat_idx] - state.ambiguities[ref_idx];
                    
                    z_vals.push(cp_dd_l2 - (comp_pr_dd - iono_dd_l2 + n_dd_l2));
                    let mut h_cp2 = vec![0.0; state_size]; 
                    h_cp2[0] = h_r.x; h_cp2[1] = h_r.y; h_cp2[2] = h_r.z;
                    h_cp2[6] = h_att.x; h_cp2[7] = h_att.y; h_cp2[8] = h_att.z;
                    h_cp2[15 + sat_idx] = 1.0; h_cp2[15 + ref_idx] = -1.0;
                    h_rows.push(h_cp2); 
                    r_vals.push(if state.is_fixed { 1e-6 * var_factor } else { 0.0001 * var_factor }); 
                    meas_type.push(2);
                }
            }
        }
    }
    
    Some((z_vals, h_rows, r_vals, meas_type))
}

pub fn build_measurement_model(
    state: &mut RtkState,
    matched_obs: &[(DdObservation, DdObservation)],
    ephemerides: &[Ephemeris],
    base_coord: &Coordinate,
    _rover_time: GpsTime,
    base_time: GpsTime,
    lever_arm: Vector3<f64>,
) -> Option<(DVector<f64>, DMatrix<f64>, DMatrix<f64>)> {
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

    for (_, group) in const_groups {
        if group.len() < 2 { continue; } 

        let mut max_snr = -1.0;
        let mut ref_idx = 0;
        for (i, (r, _)) in group.iter().enumerate() {
            if r.snr > max_snr && ephemerides.iter().any(|e| e.sat() == r.sat) {
                max_snr = r.snr;
                ref_idx = i;
            }
        }
        
        let mut group_clone = group.clone();
        let (ref_rover, ref_base) = group_clone.remove(ref_idx);

        if let Some((z, h, r, mt)) = compute_innovations(state, &group_clone, ephemerides, base_coord, base_time, &ref_rover, &ref_base, lever_arm) {
            z_all.extend(z);
            h_all.extend(h);
            r_all.extend(r);
            type_all.extend(mt);
        }
    }

    let state_size = 15 + state.ambiguities.len();
    let mut safe_indices = Vec::new();
    
    for i in 0..z_all.len() {
        let mut h_row = DMatrix::zeros(1, state_size);
        for c in 0..state_size { h_row[(0, c)] = h_all[i][c]; }
        let s_ii = (&h_row * &state.covariance * h_row.transpose())[(0, 0)] + r_all[i];
        let chi2 = z_all[i] * z_all[i] / s_ii;
        
        let threshold = match type_all[i] { 
            0 => 15.0,  
            1 | 2 => 1e7,  
            _ => 15.0   
        };
        
        if chi2 <= threshold { 
            safe_indices.push(i); 
            if type_all[i] == 1 || type_all[i] == 2 {
                for c in 15..state_size {
                    if h_row[(0, c)] > 0.5 {
                        let key = state.ambiguity_keys[c - 15];
                        state.reject_counts.insert(key, 0);
                    }
                }
            }
        } else {
            tracing::debug!("Rejected meas type {} with inn: {:.3}, chi2: {:.1}", type_all[i], z_all[i], chi2);
            if type_all[i] == 1 || type_all[i] == 2 {
                for c in 15..state_size {
                    if h_row[(0, c)] > 0.5 {
                        let key = state.ambiguity_keys[c - 15];
                        let count = *state.reject_counts.get(&key).unwrap_or(&0) + 1;
                        state.reject_counts.insert(key, count);
                    }
                }
            }
        }
    }

    if safe_indices.len() >= 4 {
        let mut z_safe = DVector::zeros(safe_indices.len());
        let mut h_safe = DMatrix::zeros(safe_indices.len(), state_size);
        let mut r_safe = DMatrix::zeros(safe_indices.len(), safe_indices.len());
        for (new_i, &old_i) in safe_indices.iter().enumerate() {
            z_safe[new_i] = z_all[old_i];
            for c in 0..state_size { h_safe[(new_i, c)] = h_all[old_i][c]; }
            r_safe[(new_i, new_i)] = r_all[old_i];
        }
        Some((z_safe, h_safe, r_safe))
    } else {
        None
    }
}
