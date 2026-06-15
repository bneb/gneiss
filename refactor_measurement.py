import re

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    content = f.read()

# Pattern for compute_innovations
compute_inn_pattern = re.compile(r'#\[allow\(clippy::too_many_arguments\)\]\npub fn compute_innovations\(.*?\n    Some\(\(z_vals, h_rows, r_vals, meas_type\)\)\n\}', re.DOTALL)

compute_inn_repl = '''#[allow(clippy::too_many_arguments)]
pub fn compute_single_innovation(
    state: &mut RtkState,
    rover_sat_orig: &DdObservation,
    base_sat_orig: &DdObservation,
    ref_rover_orig: &DdObservation,
    ref_base_orig: &DdObservation,
    env: &MeasurementEnvironment,
    geom: &EkfGeometryContext,
    ref_state: &SatState,
    sat_eph: &Ephemeris,
    e_ref_rov: Vector3<f64>,
    ref_idx_l1: Option<usize>,
    ref_idx_l2: Option<usize>,
    z_vals: &mut Vec<f64>,
    h_rows: &mut Vec<Vec<f64>>,
    r_vals: &mut Vec<f64>,
    meas_type: &mut Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
) {
    let (sat_vec_rov, sat_vel_rov) =
        get_sat_state(sat_eph, rover_sat_orig.pr_l1, state.time, geom.pos_apc);
    let (sat_vec_bas, sat_vel_bas) =
        get_sat_state(sat_eph, base_sat_orig.pr_l1, env.base_time, geom.base_coord_vec);
    let (sat_f1, sat_f2) =
        gneiss_core::signal::satellite_frequencies(rover_sat_orig.sat, sat_eph.freq_num());

    let sat_state = SatState {
        rov_pos: sat_vec_rov,
        rov_vel: sat_vel_rov,
        bas_pos: sat_vec_bas,
        bas_vel: sat_vel_bas,
        f1: sat_f1,
        f2: sat_f2,
    };

    let e_sat_rov = (sat_vec_rov - geom.pos_apc).normalize();
    let h_r = e_ref_rov - e_sat_rov;
    let h_att = geom.compute_attitude_jacobian(&h_r);

    let (tropo_dd, iono_dd_l1, iono_dd_l2) = compute_atmospheric_delays(
        state.time,
        geom.pos_apc,
        geom.base_coord_vec,
        sat_vec_rov,
        ref_state.rov_pos,
        sat_vec_bas,
        ref_state.bas_pos,
        sat_f1,
        sat_f2,
        ref_state.f1,
        ref_state.f2,
    );

    let mut rov_sat = rover_sat_orig.clone();
    let mut bas_sat = base_sat_orig.clone();
    let mut rov_ref = ref_rover_orig.clone();
    let mut bas_ref = ref_base_orig.clone();

    let mut ctx = DdContext {
        rov_sat: &mut rov_sat,
        base_sat: &mut bas_sat,
        rov_ref: &mut rov_ref,
        ref_base: &mut bas_ref,
        sat_state: &sat_state,
        ref_state,
    };

    apply_phase_windup(state, geom.pos_apc, geom.base_coord_vec, &mut ctx);

    let base_llh = gneiss_core::coords::ecef_to_llh(geom.base_coord_vec);
    let rov_llh = gneiss_core::coords::ecef_to_llh(geom.pos_apc);
    let (_, el_rov_sat) = gneiss_core::coords::az_el(rov_llh, geom.pos_apc, sat_vec_rov);
    let (_, el_rov_ref) = gneiss_core::coords::az_el(rov_llh, geom.pos_apc, ref_state.rov_pos);
    let (_, el_bas_sat) = gneiss_core::coords::az_el(base_llh, geom.base_coord_vec, sat_vec_bas);
    let (_, el_bas_ref) = gneiss_core::coords::az_el(base_llh, geom.base_coord_vec, ref_state.bas_pos);

    let (var_factor, ref_var_factor) = compute_variance_factors(&ctx, el_rov_sat, el_rov_ref, el_bas_sat, el_bas_ref, env.tuning.snr_a, env.tuning.snr_b);
    let (h_zwd, zwd_dd) = compute_zwd_mapping(el_rov_sat, el_rov_ref, state.zwd);

    let comp_pr_dd = compute_geometric_dd(
        geom.pos_apc, geom.base_coord_vec, sat_vec_rov, ref_state.rov_pos, sat_vec_bas, ref_state.bas_pos,
    ) + tropo_dd + zwd_dd;

    for update in compute_dd_pseudorange(
        &ctx, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom.state_size,
        var_factor, env, h_zwd, ref_var_factor
    ) {
        z_vals.push(update.0);
        h_rows.push(update.1);
        r_vals.push(update.2);
        meas_type.push((ctx.rov_sat.sat, update.3, update.4));
    }

    for update in compute_dd_carrier_phase(
        state,
        &ctx,
        ref_idx_l1,
        ref_idx_l2,
        comp_pr_dd,
        iono_dd_l1,
        iono_dd_l2,
        h_r,
        h_att,
        geom.state_size,
        var_factor,
        env,
        h_zwd,
        ref_var_factor
    ) {
        z_vals.push(update.0);
        h_rows.push(update.1);
        r_vals.push(update.2);
        meas_type.push((ctx.rov_sat.sat, update.3, update.4));
    }

    if let Some(update) = compute_dd_doppler(
        state,
        &ctx,
        geom.pos_apc,
        geom.base_coord_vec,
        h_r,
        geom.state_size,
        var_factor,
        env,
        ref_var_factor
    ) {
        z_vals.push(update.0);
        h_rows.push(update.1);
        r_vals.push(update.2);
        meas_type.push((ctx.rov_sat.sat, update.3, update.4));
    }
}

pub fn find_ephemeris<'a>(ephemerides: &'a [Ephemeris], sat: gneiss_core::sat::SatelliteId, time_tow: f64) -> Option<&'a Ephemeris> {
    ephemerides
        .iter()
        .filter(|e| e.sat() == sat)
        .min_by(|a, b| {
            let da = (a.toe().tow - time_tow).abs();
            let db = (b.toe().tow - time_tow).abs();
            da.partial_cmp(&db).unwrap()
        })
}

pub fn get_ref_ambiguity_indices(state: &RtkState, sat: gneiss_core::sat::SatelliteId) -> (Option<usize>, Option<usize>) {
    let idx1 = state.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == 1);
    let idx2 = state.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == 2);
    (idx1, idx2)
}

#[allow(clippy::too_many_arguments)]
pub fn compute_innovations(
    state: &mut RtkState,
    group: &[(DdObservation, DdObservation)],
    ref_rover_orig: &DdObservation,
    ref_base_orig: &DdObservation,
    env: &MeasurementEnvironment,
) -> Option<(
    Vec<f64>,
    Vec<Vec<f64>>,
    Vec<f64>,
    Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
)> {
    let mut z_vals = Vec::new();
    let mut h_rows = Vec::new();
    let mut r_vals = Vec::new();
    let mut meas_type = Vec::new();

    let geom = EkfGeometryContext::new(state, env);

    let ref_eph = find_ephemeris(env.ephemerides, ref_rover_orig.sat, state.time.tow)?;
    let (ref_sat_vec_rov, ref_sat_vel_rov) =
        get_sat_state(ref_eph, ref_rover_orig.pr_l1, state.time, geom.pos_apc);
    let (ref_sat_vec_bas, ref_sat_vel_bas) =
        get_sat_state(ref_eph, ref_base_orig.pr_l1, env.base_time, geom.base_coord_vec);
    let (ref_f1, ref_f2) =
        gneiss_core::signal::satellite_frequencies(ref_rover_orig.sat, ref_eph.freq_num());

    let ref_state = SatState {
        rov_pos: ref_sat_vec_rov,
        rov_vel: ref_sat_vel_rov,
        bas_pos: ref_sat_vec_bas,
        bas_vel: ref_sat_vel_bas,
        f1: ref_f1,
        f2: ref_f2,
    };

    let e_ref_rov = (ref_sat_vec_rov - geom.pos_apc).normalize();
    let (ref_idx_l1, ref_idx_l2) = get_ref_ambiguity_indices(state, ref_rover_orig.sat);

    for (rover_sat_orig, base_sat_orig) in group {
        if let Some(sat_eph) = find_ephemeris(env.ephemerides, rover_sat_orig.sat, state.time.tow) {
            compute_single_innovation(
                state, rover_sat_orig, base_sat_orig, ref_rover_orig, ref_base_orig,
                env, &geom, &ref_state, sat_eph, e_ref_rov, ref_idx_l1, ref_idx_l2,
                &mut z_vals, &mut h_rows, &mut r_vals, &mut meas_type
            );
        }
    }

    Some((z_vals, h_rows, r_vals, meas_type))
}'''

content = compute_inn_pattern.sub(compute_inn_repl, content)

build_meas_pattern = re.compile(r'#\[allow\(clippy::type_complexity\)\]\npub fn build_measurement_model\(.*?\n    \} else \{\n        None\n    \}\n\}', re.DOTALL)

build_meas_repl = '''pub fn find_best_reference_satellite(
    state: &RtkState,
    env: &MeasurementEnvironment,
    group: &[(DdObservation, DdObservation)],
) -> Option<usize> {
    let mut best_score = -1.0;
    let mut ref_idx = None;
    let rov_llh_ref = gneiss_core::coords::ecef_to_llh(state.position.vector);

    for (i, (rover_obs, base_obs)) in group.iter().enumerate() {
        let sat_coord = env.ephemerides.iter().find(|e| e.sat() == rover_obs.sat)
            .map(|e| e.position(env.base_time).0);

        let ele = if let Some(sat_pos) = sat_coord {
            gneiss_core::coords::az_el(rov_llh_ref, state.position.vector, sat_pos).1
        } else {
            continue;
        };

        let score = ele * 100.0 +
                    rover_obs.cp_l1.map_or(0.0, |_| 1000.0) +
                    rover_obs.cp_l2.map_or(0.0, |_| 500.0) +
                    base_obs.cp_l1.map_or(0.0, |_| 1000.0) +
                    base_obs.cp_l2.map_or(0.0, |_| 500.0);

        if score > best_score {
            best_score = score;
            ref_idx = Some(i);
        }
    }
    ref_idx
}

#[allow(clippy::too_many_arguments)]
pub fn evaluate_innovation_outliers(
    state: &mut RtkState,
    z_all: &[f64],
    h_all: &[Vec<f64>],
    r_all: &[f64],
    type_all: &[(gneiss_core::sat::SatelliteId, u8, f64)],
    chi_square_pr_threshold: f64,
    chi_square_cp_threshold: f64,
    is_tightly_coupled: bool,
    state_size: usize,
) -> Vec<usize> {
    let mut safe_indices = Vec::new();
    
    for i in 0..z_all.len() {
        let mut h_row = DMatrix::zeros(1, state_size);
        for c in 0..state_size { h_row[(0, c)] = h_all[i][c]; }
        let s_ii = (&h_row * &state.covariance * h_row.transpose())[(0, 0)] + r_all[i];
        let chi2 = z_all[i] * z_all[i] / s_ii;
        
        let mut cp_thresh_sq = chi_square_cp_threshold * chi_square_cp_threshold;
        if is_tightly_coupled {
            cp_thresh_sq = cp_thresh_sq.max(1000000.0);
        }
        
        let threshold = match type_all[i].1 { 
            0 => chi_square_pr_threshold * chi_square_pr_threshold,  
            1 | 2 => cp_thresh_sq,  
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
}

#[allow(clippy::type_complexity)]
pub fn build_measurement_model(
    state: &mut RtkState,
    matched_obs: &[(DdObservation, DdObservation)],
    env: &MeasurementEnvironment,
    chi_square_pr_threshold: f64,
    chi_square_cp_threshold: f64,
    is_tightly_coupled: bool,
) -> Option<(DVector<f64>, DMatrix<f64>, DMatrix<f64>, Vec<(gneiss_core::sat::SatelliteId, u8, f64)>)> {
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

    for (_, mut group) in const_groups {
        if group.len() < 2 { continue; } 

        if let Some(ref_idx) = find_best_reference_satellite(state, env, &group) {
            let (ref_rover, ref_base) = group.remove(ref_idx);

            if let Some((z, h, r, mt)) = compute_innovations(state, &group, &ref_rover, &ref_base, env) {
                z_all.extend(z);
                h_all.extend(h);
                r_all.extend(r);
                type_all.extend(mt);
            }
        }
    }

    let state_size = crate::filter::CORE_STATE_SIZE + state.ambiguities.len();
    
    tracing::trace!("Pre-filter z_all len: {}", z_all.len());
    
    let safe_indices = evaluate_innovation_outliers(
        state, &z_all, &h_all, &r_all, &type_all, 
        chi_square_pr_threshold, chi_square_cp_threshold, 
        is_tightly_coupled, state_size
    );

    if safe_indices.len() >= 4 {
        let mut z_vec = DVector::zeros(safe_indices.len());
        let mut h_mat = DMatrix::zeros(safe_indices.len(), state_size);
        let mut r_diagonals = Vec::new();
        let mut t_vec = Vec::new();

        for (new_i, &old_i) in safe_indices.iter().enumerate() {
            z_vec[new_i] = z_all[old_i];
            for c in 0..state_size { h_mat[(new_i, c)] = h_all[old_i][c]; }
            r_diagonals.push(r_all[old_i]);
            t_vec.push(type_all[old_i]);
        }
        let r_mat = build_dense_covariance_matrix(&r_diagonals, &t_vec);
        Some((z_vec, h_mat, r_mat, t_vec))
    } else {
        None
    }
}'''

content = build_meas_pattern.sub(build_meas_repl, content)

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(content)

