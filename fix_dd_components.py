import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

struct_def = """pub struct DdComponents {
    pub comp_pr_dd: f64,
    pub iono_dd_l1: f64,
    pub iono_dd_l2: f64,
    pub var_factor: f64,
    pub ref_var_factor: f64,
    pub h_zwd: f64,
}

"""

# Insert struct before compute_dd_components
content = content.replace("fn compute_dd_components(", struct_def + "fn compute_dd_components(")

# Fix return type
content = content.replace(") -> (f64, f64, f64, f64, f64, f64) {", ") -> DdComponents {")

# Fix return value
content = content.replace(
    "(comp_pr_dd, iono_dd_l1, iono_dd_l2, var_factor, ref_var_factor, h_zwd)",
    "DdComponents { comp_pr_dd, iono_dd_l1, iono_dd_l2, var_factor, ref_var_factor, h_zwd }"
)

# Fix caller (process_single_satellite_pair)
old_caller = "let (comp_pr_dd, iono_dd_l1, iono_dd_l2, var_factor, ref_var_factor, h_zwd) = compute_dd_components(state, geom, env, ctx);"
new_caller = "let comps = compute_dd_components(state, geom, env, ctx);"
content = content.replace(old_caller, new_caller)

# Fix variables used in generate_measurement_updates
old_call = "generate_measurement_updates(state, ctx, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, geom, var_factor, env, h_zwd, ref_var_factor, ref_idx_l1, ref_idx_l2, updates);"
new_call = "generate_measurement_updates(state, ctx, comps.comp_pr_dd, comps.iono_dd_l1, comps.iono_dd_l2, h_r, h_att, geom, comps.var_factor, env, comps.h_zwd, comps.ref_var_factor, ref_idx_l1, ref_idx_l2, updates);"
content = content.replace(old_call, new_call)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)

