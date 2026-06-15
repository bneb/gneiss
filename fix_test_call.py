import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

old_call = """        apply_phase_windup(&mut state, Vector3::zeros(), Vector3::zeros(), &mut ctx);"""
new_call = """        let (w_sat, w_ref, w_bas_sat, w_bas_ref) = compute_phase_windup(
            state.time, Vector3::zeros(), Vector3::zeros(),
            ctx.sat_state.rov_pos, ctx.ref_state.rov_pos,
            ctx.sat_state.bas_pos, ctx.ref_state.bas_pos,
            0.0, 0.0, 0.0, 0.0,
        );
        if let Some(cp) = &mut ctx.rov_sat.cp_l1 { *cp += w_sat; }
        if let Some(cp2) = &mut ctx.rov_sat.cp_l2 { *cp2 += w_sat; }
        if let Some(cp) = &mut ctx.rov_ref.cp_l1 { *cp += w_ref; }
        if let Some(cp2) = &mut ctx.rov_ref.cp_l2 { *cp2 += w_ref; }
        if let Some(cp) = &mut ctx.base_sat.cp_l1 { *cp += w_bas_sat; }
        if let Some(cp2) = &mut ctx.base_sat.cp_l2 { *cp2 += w_bas_sat; }
        if let Some(cp) = &mut ctx.ref_base.cp_l1 { *cp += w_bas_ref; }
        if let Some(cp2) = &mut ctx.ref_base.cp_l2 { *cp2 += w_bas_ref; }"""

content = content.replace(old_call, new_call)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
