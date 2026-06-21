use crate::engine::measurement::types::{
    DdContext, EkfUpdates, MeasurementEnvironment, SatState,
};
use crate::engine::measurement::{
    compute_dd_components, find_ephemeris, generate_measurement_updates,
    update_windup_state_and_obs, EkfGeometryContext,
};
use crate::engine::measurement_math;
use crate::filter::{DdObservation, RtkState};

fn compute_sat_state(
    eph: &gneiss_core::ephemeris::Ephemeris,
    rov_obs: &DdObservation,
    bas_obs: &DdObservation,
    time: gneiss_core::time::GpsTime,
    geom: &EkfGeometryContext,
    base_time: gneiss_core::time::GpsTime,
    rcv_clk_bias_m: f64,
) -> SatState {
    let (rov_pos, rov_vel) =
        measurement_math::get_sat_state(eph, rov_obs.pr_l1, rcv_clk_bias_m, time, geom.pos_apc);
    let (bas_pos, bas_vel) = measurement_math::get_sat_state(
        eph,
        bas_obs.pr_l1,
        0.0,
        base_time,
        geom.base_coord_vec,
    );
    let (f1, f2) = gneiss_core::signal::satellite_frequencies(rov_obs.sat, eph.freq_num());

    SatState {
        rov_pos,
        rov_vel,
        bas_pos,
        bas_vel,
        f1,
        f2,
    }
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
    updates: &mut EkfUpdates,
) {
    let sat_state = compute_sat_state(
        sat_eph,
        rover_sat_orig,
        base_sat_orig,
        state.time,
        geom,
        env.base_time,
        state.rcv_clk_bias,
    );

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

    update_windup_state_and_obs(state, &mut ctx, geom);

    let comps = compute_dd_components(state, geom, env, &ctx);
    let ugeom = crate::engine::measurement::types::UpdateGeometry {
        comp_dd: comps.comp_pr_dd,
        h_r,
        h_att,
        h_zwd: comps.h_zwd,
        state_size: geom.state_size,
    };
    let mctx = crate::engine::measurement::types::DdMeasurementContext {
        ctx: &ctx,
        geom: &ugeom,
        comps: &comps,
        env,
    };

    generate_measurement_updates(state, &mctx, geom, ref_idx_l1, ref_idx_l2, updates);
}

pub fn compute_innovations(
    state: &mut RtkState,
    group: &[(DdObservation, DdObservation)],
    ref_rover_orig: &DdObservation,
    ref_base_orig: &DdObservation,
    env: &MeasurementEnvironment,
) -> Option<EkfUpdates> {
    let mut updates = EkfUpdates::new();
    let geom = EkfGeometryContext::new(state, env);
    let ref_eph = find_ephemeris(env.ephemerides, ref_rover_orig.sat, state.time.tow)?;
    let ref_state = compute_sat_state(
        ref_eph,
        ref_rover_orig,
        ref_base_orig,
        state.time,
        &geom,
        env.base_time,
        state.rcv_clk_bias,
    );

    let ref_idx_l1 = state
        .ambiguity_keys
        .iter()
        .position(|&(s, f)| s == ref_rover_orig.sat && f == 1);
    let ref_idx_l2 = state
        .ambiguity_keys
        .iter()
        .position(|&(s, f)| s == ref_rover_orig.sat && f == 2);

    for (rover_sat_orig, base_sat_orig) in group {
        if let Some(sat_eph) = find_ephemeris(env.ephemerides, rover_sat_orig.sat, state.time.tow)
        {
            process_single_satellite_pair(
                state,
                rover_sat_orig,
                base_sat_orig,
                ref_rover_orig,
                ref_base_orig,
                sat_eph,
                &ref_state,
                &geom,
                env,
                ref_idx_l1,
                ref_idx_l2,
                &mut updates,
            );
        }
    }
    Some(updates)
}
