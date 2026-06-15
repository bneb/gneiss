import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

# Fix compute_dd_doppler return
content = content.replace(
    "return Some((innov, h_dop, dop_base_var * var_factor, 3, dop_base_var * ref_var_factor));",
    "return Some(SingleUpdate { z: innov, h: h_dop, r: dop_base_var * var_factor, type_code: 3, r_ref: dop_base_var * ref_var_factor });"
)

# Fix test assert
content = content.replace("assert!(!u.0.is_nan());", "assert!(!u.z.is_nan());")

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)

