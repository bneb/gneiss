import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

# Replace compute_phase_windup usage in test
old_call = """        let (w_sat, w_ref, w_bas_sat, w_bas_ref) = crate::engine::measurement_math::compute_phase_windup(
            state.time, Vector3::zeros(), Vector3::zeros(),
            ctx.sat_state.rov_pos, ctx.ref_state.rov_pos,
            ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
            0.0, 0.0, 0.0, 0.0,
        );"""

new_call = """        let sun_pos = gneiss_core::sun::sun_position_ecef(state.time);
        let crate::engine::measurement_math::WindupUpdates { w_sat, w_ref, w_bas_sat, w_bas_ref } = crate::engine::measurement_math::compute_phase_windup(
            Vector3::zeros(), Vector3::zeros(), sun_pos,
            ctx.sat_state.rov_pos, ctx.ref_state.rov_pos,
            ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
            0.0, 0.0, 0.0, 0.0,
        );"""

content = content.replace(old_call, new_call)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)

