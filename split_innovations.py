import re

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    text = f.read()

start_idx = text.find("pub fn compute_innovations(")
end_idx = text.find("pub fn build_measurement_model(", start_idx)

original = text[start_idx:end_idx]

new_funcs = """#[allow(clippy::too_many_arguments)]
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
    let time_tow = state.time.tow;

    let find_eph = |sat| {
        env.ephemerides
            .iter()
            .filter(|e| e.sat() == sat)
            .min_by(|a, b| {
                let da = (a.toe().tow - time_tow).abs();
                let db = (b.toe().tow - time_tow).abs();
                da.partial_cmp(&db).unwrap()
            })
    };

    let ref_eph = find_eph(ref_rover_orig.sat)?;
    let ref_state = compute_sat_state(ref_eph, ref_rover_orig, ref_base_orig, state.time, &geom, env.base_time);

    let ref_idx_l1 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 1);
    let ref_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 2);

    for (rover_sat_orig, base_sat_orig) in group {
        if let Some(sat_eph) = find_eph(rover_sat_orig.sat) {
            process_single_satellite_pair(
                state, rover_sat_orig, base_sat_orig, ref_rover_orig, ref_base_orig,
                sat_eph, &ref_state, &geom, env, ref_idx_l1, ref_idx_l2,
                &mut z_vals, &mut h_rows, &mut r_vals, &mut meas_type
            );
        }
    }

    Some((z_vals, h_rows, r_vals, meas_type))
}

fn compute_sat_state(
    eph: &gneiss_core::ephemeris::Ephemeris,
    rov_obs: &DdObservation,
    bas_obs: &DdObservation,
    time: gneiss_core::time::GpsTime,
    geom: &EkfGeometryContext,
    base_time: gneiss_core::time::GpsTime,
) -> SatState {
    let (rov_pos, rov_vel) = get_sat_state(eph, rov_obs.pr_l1, time, geom.pos_apc);
    let (bas_pos, bas_vel) = get_sat_state(eph, bas_obs.pr_l1, base_time, geom.base_coord_vec);
    let (f1, f2) = gneiss_core::signal::satellite_frequencies(rov_obs.sat, eph.freq_num());

    SatState { rov_pos, rov_vel, bas_pos, bas_vel, f1, f2 }
}

#[allow(clippy::too_many_arguments)]
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
) {
    let sat_state = compute_sat_state(sat_eph, rover_sat_orig, base_sat_orig, state.time, geom, env.base_time);
    
    let e_ref_rov = (ref_state.rov_pos - geom.pos_apc).normalize();
    let e_sat_rov = (sat_state.rov_pos - geom.pos_apc).normalize();
    let h_r = e_ref_rov - e_sat_rov;
    let h_att = geom.compute_attitude_jacobian(&h_r);

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

    let (comp_pr_dd, iono_dd_l1, iono_dd_l2, var_factor, ref_var_factor, h_zwd) = 
        compute_dd_components(state, geom, env, &ctx);

    generate_measurement_updates(
        state, &ctx, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom, var_factor, env, h_zwd, ref_var_factor, ref_idx_l1, ref_idx_l2,
        z_vals, h_rows, r_vals, meas_type
    );
}

fn compute_dd_components(
    state: &RtkState,
    geom: &EkfGeometryContext,
    env: &MeasurementEnvironment,
    ctx: &DdContext,
) -> (f64, f64, f64, f64, f64, f64) {
    let (tropo_dd, iono_dd_l1, iono_dd_l2) = compute_atmospheric_delays(
        state.time, geom.pos_apc, geom.base_coord_vec,
        ctx.sat_state.rov_pos, ctx.ref_state.rov_pos,
        ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
        ctx.sat_state.f1, ctx.sat_state.f2, ctx.ref_state.f1, ctx.ref_state.f2,
    );

    let base_llh = gneiss_core::coords::ecef_to_llh(geom.base_coord_vec);
    let rov_llh = gneiss_core::coords::ecef_to_llh(geom.pos_apc);
    let (_, el_rov_sat) = gneiss_core::coords::az_el(rov_llh, geom.pos_apc, ctx.sat_state.rov_pos);
    let (_, el_rov_ref) = gneiss_core::coords::az_el(rov_llh, geom.pos_apc, ctx.ref_state.rov_pos);
    let (_, el_bas_sat) = gneiss_core::coords::az_el(base_llh, geom.base_coord_vec, ctx.sat_state.bas_pos);
    let (_, el_bas_ref) = gneiss_core::coords::az_el(base_llh, geom.base_coord_vec, ctx.ref_state.bas_pos);

    let (var_factor, ref_var_factor) = compute_variance_factors(ctx.rov_sat.snr, ctx.rov_ref.snr, el_rov_sat, el_rov_ref, el_bas_sat, el_bas_ref, env.tuning.snr_a, env.tuning.snr_b);
    let (h_zwd, zwd_dd) = compute_zwd_mapping(el_rov_sat, el_rov_ref, state.zwd);

    let comp_pr_dd = compute_geometric_dd(
        geom.pos_apc, geom.base_coord_vec, ctx.sat_state.rov_pos, ctx.ref_state.rov_pos, ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
    ) + tropo_dd + zwd_dd;

    (comp_pr_dd, iono_dd_l1, iono_dd_l2, var_factor, ref_var_factor, h_zwd)
}

#[allow(clippy::too_many_arguments)]
fn generate_measurement_updates(
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

    for update in compute_dd_carrier_phase(
        state, ctx, ref_idx_l1, ref_idx_l2, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom.state_size,
        var_factor, env, h_zwd, ref_var_factor
    ) {
        z_vals.push(update.0);
        h_rows.push(update.1);
        r_vals.push(update.2);
        meas_type.push((ctx.rov_sat.sat, update.3, update.4));
    }

    if let Some(update) = compute_dd_doppler(
        state, ctx, geom.pos_apc, geom.base_coord_vec, h_r, geom.state_size,
        var_factor, env, ref_var_factor
    ) {
        z_vals.push(update.0);
        h_rows.push(update.1);
        r_vals.push(update.2);
        meas_type.push((ctx.rov_sat.sat, update.3, update.4));
    }
}

"""

new_text = text[:start_idx] + new_funcs + text[end_idx:]

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(new_text)

