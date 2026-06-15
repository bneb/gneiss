import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

# 1. apply_phase_windup
# We will make it return the 4 new windup values, and remove `state: &mut RtkState` and `ctx: &mut DdContext`.
old_apply = """pub fn apply_phase_windup(
    state: &mut RtkState,
    pos_apc: Vector3<f64>,
    base_coord_vec: Vector3<f64>,
    ctx: &mut DdContext,
) {
    let sun_pos = gneiss_core::sun::sun_position_ecef(state.time);

    let prev_w_sat = *state.windup.get(&ctx.rov_sat.sat).unwrap_or(&0.0);
    let prev_w_ref = *state.windup.get(&ctx.rov_ref.sat).unwrap_or(&0.0);
    let w_sat = gneiss_core::windup::phase_windup(ctx.sat_state.rov_pos, sun_pos, pos_apc, prev_w_sat);
    let w_ref = gneiss_core::windup::phase_windup(ctx.ref_state.rov_pos, sun_pos, pos_apc, prev_w_ref);
    state.windup.insert(ctx.rov_sat.sat, w_sat);
    state.windup.insert(ctx.rov_ref.sat, w_ref);

    if let Some(cp) = &mut ctx.rov_sat.cp_l1 { *cp += w_sat; }
    if let Some(cp2) = &mut ctx.rov_sat.cp_l2 { *cp2 += w_sat; }
    if let Some(cp) = &mut ctx.rov_ref.cp_l1 { *cp += w_ref; }
    if let Some(cp2) = &mut ctx.rov_ref.cp_l2 { *cp2 += w_ref; }

    let prev_w_bas_sat = *state.windup.get(&ctx.base_sat.sat).unwrap_or(&0.0);
    let prev_w_bas_ref = *state.windup.get(&ctx.ref_base.sat).unwrap_or(&0.0);
    let w_bas_sat = gneiss_core::windup::phase_windup(ctx.sat_state.bas_pos, sun_pos, base_coord_vec, prev_w_bas_sat);
    let w_bas_ref = gneiss_core::windup::phase_windup(ctx.ref_state.bas_pos, sun_pos, base_coord_vec, prev_w_bas_ref);
    state.windup.insert(ctx.base_sat.sat, w_bas_sat);
    state.windup.insert(ctx.ref_base.sat, w_bas_ref);

    if let Some(cp) = &mut ctx.base_sat.cp_l1 { *cp += w_bas_sat; }
    if let Some(cp2) = &mut ctx.base_sat.cp_l2 { *cp2 += w_bas_sat; }
    if let Some(cp) = &mut ctx.ref_base.cp_l1 { *cp += w_bas_ref; }
    if let Some(cp2) = &mut ctx.ref_base.cp_l2 { *cp2 += w_bas_ref; }
}"""

new_apply = """pub fn compute_phase_windup(
    time: GpsTime,
    pos_apc: Vector3<f64>,
    base_coord_vec: Vector3<f64>,
    rov_sat_pos: Vector3<f64>,
    rov_ref_pos: Vector3<f64>,
    bas_sat_pos: Vector3<f64>,
    bas_ref_pos: Vector3<f64>,
    prev_w_sat: f64,
    prev_w_ref: f64,
    prev_w_bas_sat: f64,
    prev_w_bas_ref: f64,
) -> (f64, f64, f64, f64) {
    let sun_pos = gneiss_core::sun::sun_position_ecef(time);
    let w_sat = gneiss_core::windup::phase_windup(rov_sat_pos, sun_pos, pos_apc, prev_w_sat);
    let w_ref = gneiss_core::windup::phase_windup(rov_ref_pos, sun_pos, pos_apc, prev_w_ref);
    let w_bas_sat = gneiss_core::windup::phase_windup(bas_sat_pos, sun_pos, base_coord_vec, prev_w_bas_sat);
    let w_bas_ref = gneiss_core::windup::phase_windup(bas_ref_pos, sun_pos, base_coord_vec, prev_w_bas_ref);
    (w_sat, w_ref, w_bas_sat, w_bas_ref)
}"""

content = content.replace(old_apply, new_apply)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
