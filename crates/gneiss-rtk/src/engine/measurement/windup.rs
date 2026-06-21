use crate::engine::measurement::EkfGeometryContext;
use crate::engine::measurement::types::DdContext;
use crate::filter::{DdObservation, RtkState};

/// Computes the phase windup correction, updates the state's windup map,
/// and returns the windup values for each satellite-receiver pair.
fn update_windup_state(
    state: &mut RtkState,
    ctx: &mut DdContext,
    geom: &EkfGeometryContext,
) -> (f64, f64, f64, f64) {
    let prev_w_sat = *state.windup.get(&ctx.rov_sat.sat).unwrap_or(&0.0);
    let prev_w_ref = *state.windup.get(&ctx.rov_ref.sat).unwrap_or(&0.0);
    let prev_w_bas_sat = *state.windup.get(&ctx.base_sat.sat).unwrap_or(&0.0);
    let prev_w_bas_ref = *state.windup.get(&ctx.ref_base.sat).unwrap_or(&0.0);

    let sun_pos = gneiss_core::sun::sun_position_ecef(state.time);
    let crate::engine::measurement_math::WindupUpdates {
        w_sat,
        w_ref,
        w_bas_sat,
        w_bas_ref,
    } = crate::engine::measurement_math::compute_phase_windup(
        geom.pos_apc,
        geom.base_coord_vec,
        sun_pos,
        ctx.sat_state.rov_pos,
        ctx.ref_state.rov_pos,
        ctx.sat_state.bas_pos,
        ctx.ref_state.bas_pos,
        prev_w_sat,
        prev_w_ref,
        prev_w_bas_sat,
        prev_w_bas_ref,
    );

    state.windup.insert(ctx.rov_sat.sat, w_sat);
    state.windup.insert(ctx.rov_ref.sat, w_ref);
    state.windup.insert(ctx.base_sat.sat, w_bas_sat);
    state.windup.insert(ctx.ref_base.sat, w_bas_ref);

    (w_sat, w_ref, w_bas_sat, w_bas_ref)
}

pub(crate) fn update_windup_state_and_obs(
    state: &mut RtkState,
    ctx: &mut DdContext,
    geom: &EkfGeometryContext,
) {
    let (w_sat, w_ref, w_bas_sat, w_bas_ref) = update_windup_state(state, ctx, geom);
    apply_windup_to_obs(ctx.rov_sat, w_sat);
    apply_windup_to_obs(ctx.rov_ref, w_ref);
    apply_windup_to_obs(ctx.base_sat, w_bas_sat);
    apply_windup_to_obs(ctx.ref_base, w_bas_ref);
}

pub(crate) fn apply_windup_to_obs(obs: &mut DdObservation, windup: f64) {
    if let Some(cp) = &mut obs.cp_l1 {
        *cp -= windup;
    }
    if let Some(cp2) = &mut obs.cp_l2 {
        *cp2 -= windup;
    }
}
