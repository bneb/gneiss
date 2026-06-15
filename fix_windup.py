import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

# Replace compute_phase_windup usage
old_call = """    let (w_sat, w_ref, w_bas_sat, w_bas_ref) = crate::engine::measurement_math::compute_phase_windup(
        state.time, geom.pos_apc, geom.base_coord_vec,
        ctx.sat_state.rov_pos, ctx.ref_state.rov_pos,
        ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
        prev_w_sat, prev_w_ref, prev_w_bas_sat, prev_w_bas_ref,
    );"""

new_call = """    let sun_pos = gneiss_core::sun::sun_position_ecef(state.time);
    let crate::engine::measurement_math::WindupUpdates { w_sat, w_ref, w_bas_sat, w_bas_ref } = crate::engine::measurement_math::compute_phase_windup(
        geom.pos_apc, geom.base_coord_vec, sun_pos,
        ctx.sat_state.rov_pos, ctx.ref_state.rov_pos,
        ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
        prev_w_sat, prev_w_ref, prev_w_bas_sat, prev_w_bas_ref,
    );"""

content = content.replace(old_call, new_call)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)

