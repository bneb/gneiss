import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

old_gen = """fn generate_measurement_updates(
    state: &mut RtkState,
    ctx: &DdContext,
    comp_pr_dd: f64,
    iono_dd_l1: f64,
    iono_dd_l2: f64,
    h_r: Vector3<f64>,
    h_att: Vector3<f64>,
    geom: &EkfGeometryContext,
    var_factor: f64,
    env: &MeasurementEnvironment,
    h_zwd: f64,
    ref_var_factor: f64,
    ref_idx_l1: Option<usize>,
    ref_idx_l2: Option<usize>,
    z_vals: &mut Vec<f64>,
    h_rows: &mut Vec<Vec<f64>>,
    r_vals: &mut Vec<f64>,
    meas_type: &mut Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
) {
    for update in compute_dd_pseudorange(
        ctx, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom.state_size,
        var_factor, env, h_zwd, ref_var_factor
    ) {
        z_vals.push(update.0);
        h_rows.push(update.1);
        r_vals.push(update.2);
        meas_type.push((ctx.rov_sat.sat, update.3, update.4));
    }

    let sat_idx_l1 = state.ambiguity_keys.iter().position(|&(s, f)| s == ctx.rov_sat.sat && f == 1);
    let sat_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ctx.rov_sat.sat && f == 2);

    for update in compute_dd_carrier_phase(
        ctx, state.is_fixed, &state.ambiguities, sat_idx_l1, ref_idx_l1, sat_idx_l2, ref_idx_l2,
        comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom.state_size,
        var_factor, ref_var_factor, env.tuning.cp_base_var, h_zwd
    ) {
        z_vals.push(update.0);
        h_rows.push(update.1);
        r_vals.push(update.2);
        meas_type.push((ctx.rov_sat.sat, update.3, update.4));
    }

    let r_b_e_rot = state.attitude.to_rotation_matrix();
    if let Some(update) = compute_dd_doppler(
        ctx, geom.pos_apc, geom.base_coord_vec, h_r, &r_b_e_rot,
        &state.velocity, &env.omega_b, &env.lever_arm, env.tuning.dop_base_var,
        geom.state_size, var_factor, ref_var_factor
    ) {
        z_vals.push(update.0);
        h_rows.push(update.1);
        r_vals.push(update.2);
        meas_type.push((ctx.rov_sat.sat, update.3, update.4));
    }
}"""

new_gen = """fn push_update(
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

    let r_b_e_rot = state.attitude.to_rotation_matrix();
    if let Some(u) = compute_dd_doppler(ctx, geom.pos_apc, geom.base_coord_vec, h_r, &r_b_e_rot, &state.velocity, &env.omega_b, &env.lever_arm, env.tuning.dop_base_var, geom.state_size, var_factor, ref_var_factor) {
        push_update(u, sat, z_vals, h_rows, r_vals, meas_type);
    }
}"""

content = content.replace(old_gen, new_gen)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
