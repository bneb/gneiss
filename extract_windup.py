import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

old_proc = """    let mut ctx = DdContext {
        rov_sat: &mut rov_sat,
        base_sat: &mut bas_sat,
        rov_ref: &mut rov_ref,
        ref_base: &mut bas_ref,
        sat_state: &sat_state,
        ref_state,
    };

    let prev_w_sat = *state.windup.get(&ctx.rov_sat.sat).unwrap_or(&0.0);
    let prev_w_ref = *state.windup.get(&ctx.rov_ref.sat).unwrap_or(&0.0);
    let prev_w_bas_sat = *state.windup.get(&ctx.base_sat.sat).unwrap_or(&0.0);
    let prev_w_bas_ref = *state.windup.get(&ctx.ref_base.sat).unwrap_or(&0.0);

    let (w_sat, w_ref, w_bas_sat, w_bas_ref) = compute_phase_windup(
        state.time, geom.pos_apc, geom.base_coord_vec,
        ctx.sat_state.rov_pos, ctx.ref_state.rov_pos,
        ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
        prev_w_sat, prev_w_ref, prev_w_bas_sat, prev_w_bas_ref,
    );

    state.windup.insert(ctx.rov_sat.sat, w_sat);
    state.windup.insert(ctx.rov_ref.sat, w_ref);
    state.windup.insert(ctx.base_sat.sat, w_bas_sat);
    state.windup.insert(ctx.ref_base.sat, w_bas_ref);

    apply_windup_to_obs(ctx.rov_sat, w_sat);
    apply_windup_to_obs(ctx.rov_ref, w_ref);
    apply_windup_to_obs(ctx.base_sat, w_bas_sat);
    apply_windup_to_obs(ctx.ref_base, w_bas_ref);"""

new_proc = """    let mut ctx = DdContext {
        rov_sat: &mut rov_sat,
        base_sat: &mut bas_sat,
        rov_ref: &mut rov_ref,
        ref_base: &mut bas_ref,
        sat_state: &sat_state,
        ref_state,
    };

    update_windup_state_and_obs(state, &mut ctx, geom);"""

helper = """fn update_windup_state_and_obs(state: &mut RtkState, ctx: &mut DdContext, geom: &EkfGeometryContext) {
    let prev_w_sat = *state.windup.get(&ctx.rov_sat.sat).unwrap_or(&0.0);
    let prev_w_ref = *state.windup.get(&ctx.rov_ref.sat).unwrap_or(&0.0);
    let prev_w_bas_sat = *state.windup.get(&ctx.base_sat.sat).unwrap_or(&0.0);
    let prev_w_bas_ref = *state.windup.get(&ctx.ref_base.sat).unwrap_or(&0.0);

    let (w_sat, w_ref, w_bas_sat, w_bas_ref) = compute_phase_windup(
        state.time, geom.pos_apc, geom.base_coord_vec,
        ctx.sat_state.rov_pos, ctx.ref_state.rov_pos,
        ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
        prev_w_sat, prev_w_ref, prev_w_bas_sat, prev_w_bas_ref,
    );

    state.windup.insert(ctx.rov_sat.sat, w_sat);
    state.windup.insert(ctx.rov_ref.sat, w_ref);
    state.windup.insert(ctx.base_sat.sat, w_bas_sat);
    state.windup.insert(ctx.ref_base.sat, w_bas_ref);

    apply_windup_to_obs(ctx.rov_sat, w_sat);
    apply_windup_to_obs(ctx.rov_ref, w_ref);
    apply_windup_to_obs(ctx.base_sat, w_bas_sat);
    apply_windup_to_obs(ctx.ref_base, w_bas_ref);
}

"""

content = content.replace(old_proc, new_proc)
content = content.replace("fn apply_windup_to_obs", helper + "fn apply_windup_to_obs")

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
