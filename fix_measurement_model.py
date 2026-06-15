import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

old_gen = """#[allow(clippy::too_many_arguments)]
fn generate_measurement_updates(
    state: &mut RtkState, ctx: &DdContext, comp_pr_dd: f64, iono_dd_l1: f64, iono_dd_l2: f64,
    h_r: Vector3<f64>, h_att: Vector3<f64>, geom: &EkfGeometryContext, var_factor: f64,
    env: &MeasurementEnvironment, h_zwd: f64, ref_var_factor: f64,
    ref_idx_l1: Option<usize>, ref_idx_l2: Option<usize>, updates: &mut EkfUpdates
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

    let r_b_e_rot = state.attitude.to_rotation_matrix();
    if let Some(u) = compute_dd_doppler(ctx, geom.pos_apc, geom.base_coord_vec, h_r, &r_b_e_rot, &state.velocity, &env.omega_b, &env.lever_arm, env.tuning.dop_base_var, geom.state_size, var_factor, ref_var_factor) {
        push_update(u, sat, z_vals, h_rows, r_vals, meas_type);
    }
}"""

new_gen = """#[allow(clippy::too_many_arguments)]
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

    let r_b_e_rot = state.attitude.to_rotation_matrix();
    if let Some(u) = compute_dd_doppler(ctx, geom.pos_apc, geom.base_coord_vec, h_r, &r_b_e_rot, &state.velocity, &env.omega_b, &env.lever_arm, env.tuning.dop_base_var, geom.state_size, var_factor, ref_var_factor) {
        updates.push(u, sat);
    }
}"""

if old_gen in content:
    content = content.replace(old_gen, new_gen)

old_build = """pub fn build_measurement_model(
    state: &mut RtkState,
    matched_obs: &[(DdObservation, DdObservation)],
    env: &MeasurementEnvironment,
    chi_square_pr_threshold: f64,
    chi_square_cp_threshold: f64,
) -> Option<(DVector<f64>, DMatrix<f64>, DMatrix<f64>, Vec<(gneiss_core::sat::SatelliteId, u8, f64)>)> {
    let mut z_all = Vec::new();
    let mut h_all = Vec::new();
    let mut r_all = Vec::new();
    let mut type_all = Vec::new();

    let const_groups = group_measurements_by_constellation(matched_obs);

    for (_, group) in const_groups {
        if group.len() < 2 { continue; } 

        let ref_idx = select_reference_satellite(&group, state, env);
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

    build_final_measurement_matrices(state_size, safe_indices, &z_all, &h_all, &r_all, &type_all)
}"""

new_build = """pub fn build_measurement_model(
    state: &mut RtkState, matched_obs: &[(DdObservation, DdObservation)],
    env: &MeasurementEnvironment, chi_square_pr_threshold: f64, chi_square_cp_threshold: f64,
) -> Option<(DVector<f64>, DMatrix<f64>, DMatrix<f64>, Vec<(gneiss_core::sat::SatelliteId, u8, f64)>)> {
    let mut all = EkfUpdates::new();
    let const_groups = group_measurements_by_constellation(matched_obs);

    for (_, group) in const_groups {
        if group.len() < 2 { continue; } 
        let ref_idx = select_reference_satellite(&group, state, env);
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

    build_final_measurement_matrices(state_size, safe_indices, &all.z, &all.h, &all.r, &all.mt)
}"""

if old_build in content:
    content = content.replace(old_build, new_build)

# Remove push_update fn because it is not used anymore.
old_push_update = """#[allow(clippy::too_many_arguments)]
fn push_update(
    update: (f64, Vec<f64>, f64, u8, f64), sat: gneiss_core::sat::SatelliteId,
    z: &mut Vec<f64>, h: &mut Vec<Vec<f64>>, r: &mut Vec<f64>,
    mt: &mut Vec<(gneiss_core::sat::SatelliteId, u8, f64)>
) {
    z.push(update.0);
    h.push(update.1);
    r.push(update.2);
    mt.push((sat, update.3, update.4));
}"""
content = content.replace(old_push_update, "")

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)

