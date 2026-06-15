import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

old_cp = """        let updates = compute_dd_carrier_phase(
            &ctx, state.is_fixed, &state.ambiguities, Some(1), Some(3), None, None,
            0.0, 0.0, 0.0, Vector3::new(1.0, 0.0, 0.0), Vector3::zeros(), state_size, 1.0, 1.0, env.tuning.cp_base_var, 0.0
        );"""

new_cp = """        let updates = compute_dd_carrier_phase(
            &ctx, state.is_fixed, &state.ambiguities, Some(0), Some(1), Some(2), Some(3),
            0.0, 0.0, 0.0, Vector3::new(1.0, 0.0, 0.0), Vector3::zeros(), state_size, 1.0, 1.0, env.tuning.cp_base_var, 0.0
        );"""

content = content.replace(old_cp, new_cp)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
