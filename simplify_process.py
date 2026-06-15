import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

old_proc = """    let (w_sat, w_ref, w_bas_sat, w_bas_ref) = compute_phase_windup(
        state.time, geom.pos_apc, geom.base_coord_vec,
        ctx.sat_state.rov_pos, ctx.ref_state.rov_pos,
        ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
        prev_w_sat, prev_w_ref, prev_w_bas_sat, prev_w_bas_ref,
    );

    state.windup.insert(ctx.rov_sat.sat, w_sat);
    state.windup.insert(ctx.rov_ref.sat, w_ref);
    state.windup.insert(ctx.base_sat.sat, w_bas_sat);
    state.windup.insert(ctx.ref_base.sat, w_bas_ref);

    if let Some(cp) = &mut ctx.rov_sat.cp_l1 { *cp += w_sat; }
    if let Some(cp2) = &mut ctx.rov_sat.cp_l2 { *cp2 += w_sat; }
    if let Some(cp) = &mut ctx.rov_ref.cp_l1 { *cp += w_ref; }
    if let Some(cp2) = &mut ctx.rov_ref.cp_l2 { *cp2 += w_ref; }
    if let Some(cp) = &mut ctx.base_sat.cp_l1 { *cp += w_bas_sat; }
    if let Some(cp2) = &mut ctx.base_sat.cp_l2 { *cp2 += w_bas_sat; }
    if let Some(cp) = &mut ctx.ref_base.cp_l1 { *cp += w_bas_ref; }
    if let Some(cp2) = &mut ctx.ref_base.cp_l2 { *cp2 += w_bas_ref; }"""

new_proc = """    let (w_sat, w_ref, w_bas_sat, w_bas_ref) = compute_phase_windup(
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

content = content.replace(old_proc, new_proc)

helper = """fn apply_windup_to_obs(obs: &mut DdObservation, windup: f64) {
    if let Some(cp) = &mut obs.cp_l1 { *cp += windup; }
    if let Some(cp2) = &mut obs.cp_l2 { *cp2 += windup; }
}
"""

content = helper + content

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
