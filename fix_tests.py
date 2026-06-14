import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

# Fix compute_innovations in tests
content = content.replace(
    "let (z, _h, r, _) = super::super::measurement::compute_innovations(&mut state, &matched_obs, &ref_rover, &ref_base, &env).unwrap();",
    "let updates = super::super::measurement::compute_innovations(&mut state, &matched_obs, &ref_rover, &ref_base, &env).unwrap();\n        let z = updates.z;\n        let r = updates.r;"
)

content = content.replace(
    "let (z_vals, h_rows, r_vals, meas_type) = res.unwrap();",
    "let updates = res.unwrap();\n        let (z_vals, h_rows, r_vals, meas_type) = (updates.z, updates.h, updates.r, updates.mt);"
)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
