import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

updates_struct = """pub struct EkfUpdates {
    pub z: Vec<f64>,
    pub h: Vec<Vec<f64>>,
    pub r: Vec<f64>,
    pub mt: Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
}

impl EkfUpdates {
    pub fn new() -> Self {
        Self { z: Vec::new(), h: Vec::new(), r: Vec::new(), mt: Vec::new() }
    }
    pub fn push(&mut self, u: (f64, Vec<f64>, f64, u8, f64), sat: gneiss_core::sat::SatelliteId) {
        self.z.push(u.0); self.h.push(u.1); self.r.push(u.2); self.mt.push((sat, u.3, u.4));
    }
    pub fn extend(&mut self, other: Self) {
        self.z.extend(other.z); self.h.extend(other.h); self.r.extend(other.r); self.mt.extend(other.mt);
    }
}
"""

content = content.replace("pub struct DdContext<'a> {", updates_struct + "\npub struct DdContext<'a> {")

old_push = """#[allow(clippy::too_many_arguments)]
fn push_update(
    update: (f64, Vec<f64>, f64, u8, f64), sat: gneiss_core::sat::SatelliteId,
    z: &mut Vec<f64>, h: &mut Vec<Vec<f64>>, r: &mut Vec<f64>,
    mt: &mut Vec<(gneiss_core::sat::SatelliteId, u8, f64)>
) {
    z.push(update.0);
    h.push(update.1);
    r.push(update.2);
    mt.push((sat, update.3, update.4));
}

#[allow(clippy::too_many_arguments)]
fn generate_measurement_updates(
    state: &mut RtkState, ctx: &DdContext, comp_pr_dd: f64, iono_dd_l1: f64, iono_dd_l2: f64,
    h_r: Vector3<f64>, h_att: Vector3<f64>, geom: &EkfGeometryContext, var_factor: f64,
    env: &MeasurementEnvironment, h_zwd: f64, ref_var_factor: f64,
    ref_idx_l1: Option<usize>, ref_idx_l2: Option<usize>,
    z_vals: &mut Vec<f64>, h_rows: &mut Vec<Vec<f64>>, r_vals: &mut Vec<f64>,
    meas_type: &mut Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
) {
    let sat = ctx.rov_sat.sat;
    for u in compute_dd_pseudorange(ctx, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom.state_size, var_factor, env, h_zwd, ref_var_factor) {
        push_update(u, sat, z_vals, h_rows, r_vals, meas_type);
    }

    let sat_idx_l1 = state.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == 1);
    let sat_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == 2);
    for u in compute_dd_carrier_phase(ctx, state.is_fixed, &state.ambiguities, sat_idx_l1, ref_idx_l1, sat_idx_l2, ref_idx_l2, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom.state_size, var_factor, ref_var_factor, env.tuning.cp_base_var, h_zwd) {
        push_update(u, sat, z_vals, h_rows, r_vals, meas_type);
    }

    if let Some(u) = compute_dd_doppler(ctx, geom.pos_apc, geom.base_coord_vec, h_r, &geom.r_b_e, &state.velocity.vector, &state.omega_b, &geom.lever_arm, env.tuning.dop_base_var, geom.state_size, var_factor, ref_var_factor) {
        push_update(u, sat, z_vals, h_rows, r_vals, meas_type);
    }
}"""

new_push = """#[allow(clippy::too_many_arguments)]
fn generate_measurement_updates(
    state: &mut RtkState, ctx: &DdContext, comp_pr_dd: f64, iono_dd_l1: f64, iono_dd_l2: f64,
    h_r: Vector3<f64>, h_att: Vector3<f64>, geom: &EkfGeometryContext, var_factor: f64,
    env: &MeasurementEnvironment, h_zwd: f64, ref_var_factor: f64,
    ref_idx_l1: Option<usize>, ref_idx_l2: Option<usize>, updates: &mut EkfUpdates
) {
    let sat = ctx.rov_sat.sat;
    for u in compute_dd_pseudorange(ctx, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom.state_size, var_factor, env, h_zwd, ref_var_factor) {
        updates.push(u, sat);
    }

    let sat_idx_l1 = state.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == 1);
    let sat_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == sat && f == 2);
    for u in compute_dd_carrier_phase(ctx, state.is_fixed, &state.ambiguities, sat_idx_l1, ref_idx_l1, sat_idx_l2, ref_idx_l2, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom.state_size, var_factor, ref_var_factor, env.tuning.cp_base_var, h_zwd) {
        updates.push(u, sat);
    }

    if let Some(u) = compute_dd_doppler(ctx, geom.pos_apc, geom.base_coord_vec, h_r, &geom.r_b_e, &state.velocity.vector, &state.omega_b, &geom.lever_arm, env.tuning.dop_base_var, geom.state_size, var_factor, ref_var_factor) {
        updates.push(u, sat);
    }
}"""

content = content.replace(old_push, new_push)

old_innov_sig = """pub fn compute_innovations(
    state: &mut RtkState, group: &[(DdObservation, DdObservation)],
    ref_rover_orig: &DdObservation, ref_base_orig: &DdObservation,
    env: &MeasurementEnvironment,
) -> Option<(Vec<f64>, Vec<Vec<f64>>, Vec<f64>, Vec<(gneiss_core::sat::SatelliteId, u8, f64)>)> {"""

new_innov_sig = """pub fn compute_innovations(
    state: &mut RtkState, group: &[(DdObservation, DdObservation)],
    ref_rover_orig: &DdObservation, ref_base_orig: &DdObservation,
    env: &MeasurementEnvironment,
) -> Option<EkfUpdates> {"""

content = content.replace(old_innov_sig, new_innov_sig)

old_innov_body = """    let mut z_vals = Vec::new(); let mut h_rows = Vec::new();
    let mut r_vals = Vec::new(); let mut meas_type = Vec::new();

    let geom = EkfGeometryContext::new(state, env);
    let ref_eph = find_ephemeris(env.ephemerides, ref_rover_orig.sat, state.time.tow)?;
    let ref_state = compute_sat_state(ref_eph, ref_rover_orig, ref_base_orig, state.time, &geom, env.base_time);

    let ref_idx_l1 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 1);
    let ref_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 2);

    for (rover_sat_orig, base_sat_orig) in group {
        if let Some(sat_eph) = find_ephemeris(env.ephemerides, rover_sat_orig.sat, state.time.tow) {
            process_single_satellite_pair(
                state, rover_sat_orig, base_sat_orig, ref_rover_orig, ref_base_orig,
                sat_eph, &ref_state, &geom, env, ref_idx_l1, ref_idx_l2,
                &mut z_vals, &mut h_rows, &mut r_vals, &mut meas_type
            );
        }
    }
    Some((z_vals, h_rows, r_vals, meas_type))"""

new_innov_body = """    let mut updates = EkfUpdates::new();
    let geom = EkfGeometryContext::new(state, env);
    let ref_eph = find_ephemeris(env.ephemerides, ref_rover_orig.sat, state.time.tow)?;
    let ref_state = compute_sat_state(ref_eph, ref_rover_orig, ref_base_orig, state.time, &geom, env.base_time);

    let ref_idx_l1 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 1);
    let ref_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 2);

    for (rover_sat_orig, base_sat_orig) in group {
        if let Some(sat_eph) = find_ephemeris(env.ephemerides, rover_sat_orig.sat, state.time.tow) {
            process_single_satellite_pair(
                state, rover_sat_orig, base_sat_orig, ref_rover_orig, ref_base_orig,
                sat_eph, &ref_state, &geom, env, ref_idx_l1, ref_idx_l2, &mut updates
            );
        }
    }
    Some(updates)"""

content = content.replace(old_innov_body, new_innov_body)

old_pair_sig = """#[allow(clippy::too_many_arguments)]
fn process_single_satellite_pair(
    state: &mut RtkState,
    rover_sat_orig: &DdObservation,
    base_sat_orig: &DdObservation,
    ref_rover_orig: &DdObservation,
    ref_base_orig: &DdObservation,
    sat_eph: &gneiss_core::ephemeris::Ephemeris,
    ref_state: &SatState,
    geom: &EkfGeometryContext,
    env: &MeasurementEnvironment,
    ref_idx_l1: Option<usize>,
    ref_idx_l2: Option<usize>,
    z_vals: &mut Vec<f64>,
    h_rows: &mut Vec<Vec<f64>>,
    r_vals: &mut Vec<f64>,
    meas_type: &mut Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
) {"""

new_pair_sig = """#[allow(clippy::too_many_arguments)]
fn process_single_satellite_pair(
    state: &mut RtkState, rover_sat_orig: &DdObservation, base_sat_orig: &DdObservation,
    ref_rover_orig: &DdObservation, ref_base_orig: &DdObservation,
    sat_eph: &gneiss_core::ephemeris::Ephemeris, ref_state: &SatState,
    geom: &EkfGeometryContext, env: &MeasurementEnvironment,
    ref_idx_l1: Option<usize>, ref_idx_l2: Option<usize>, updates: &mut EkfUpdates,
) {"""

content = content.replace(old_pair_sig, new_pair_sig)

old_pair_call = """generate_measurement_updates(
        state, &ctx, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom, var_factor, env, h_zwd, ref_var_factor, ref_idx_l1, ref_idx_l2,
        z_vals, h_rows, r_vals, meas_type
    );"""

new_pair_call = """generate_measurement_updates(
        state, &ctx, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom, var_factor, env, h_zwd, ref_var_factor, ref_idx_l1, ref_idx_l2, updates
    );"""

content = content.replace(old_pair_call, new_pair_call)

old_build_body = """    let mut z_all = Vec::new();
    let mut h_all = Vec::new();
    let mut r_all = Vec::new();
    let mut type_all = Vec::new();

    let matched_obs = filter_and_match_observations(state, env);
    if matched_obs.len() < 4 { return None; }

    let const_groups = group_measurements_by_constellation(&matched_obs);

    for group in const_groups.values() {
        if group.is_empty() { continue; }
        let ref_idx = select_reference_satellite(group, state, env);
        let mut group_clone = group.clone();
        let (ref_rover, ref_base) = group_clone.remove(ref_idx);

        if let Some((z, h, r, mt)) = compute_innovations(state, &group_clone, &ref_rover, &ref_base, env) {
            z_all.extend(z);
            h_all.extend(h);
            r_all.extend(r);
            type_all.extend(mt);
        }
    }

    let state_size = crate::filter::CORE_STATE_SIZE + state.ambiguities.len();
    tracing::trace!("Pre-filter z_all len: {}", z_all.len());
    
    let safe_indices = filter_innovations_chi_squared(
        state, state_size, chi_square_pr_threshold, chi_square_cp_threshold,
        &z_all, &h_all, &r_all, &type_all
    );

    build_final_measurement_matrices(state_size, safe_indices, &z_all, &h_all, &r_all, &type_all)"""

new_build_body = """    let mut all = EkfUpdates::new();

    let matched_obs = filter_and_match_observations(state, env);
    if matched_obs.len() < 4 { return None; }

    let const_groups = group_measurements_by_constellation(&matched_obs);

    for group in const_groups.values() {
        if group.is_empty() { continue; }
        let ref_idx = select_reference_satellite(group, state, env);
        let mut group_clone = group.clone();
        let (ref_rover, ref_base) = group_clone.remove(ref_idx);

        if let Some(updates) = compute_innovations(state, &group_clone, &ref_rover, &ref_base, env) {
            all.extend(updates);
        }
    }

    let state_size = crate::filter::CORE_STATE_SIZE + state.ambiguities.len();
    tracing::trace!("Pre-filter z_all len: {}", all.z.len());
    
    let safe_indices = filter_innovations_chi_squared(
        state, state_size, chi_square_pr_threshold, chi_square_cp_threshold,
        &all.z, &all.h, &all.r, &all.mt
    );

    build_final_measurement_matrices(state_size, safe_indices, &all.z, &all.h, &all.r, &all.mt)"""

content = content.replace(old_build_body, new_build_body)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
