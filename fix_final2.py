with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    text = f.read()

# Update signature
text = text.replace(
    "pub fn compute_variance_factors(\n    ctx: &DdContext,\n    el_rov_sat: f64,\n    el_rov_ref: f64,\n    el_bas_sat: f64,\n    el_bas_ref: f64,\n) -> (f64, f64) {",
    "pub fn compute_variance_factors(\n    ctx: &DdContext,\n    el_rov_sat: f64,\n    el_rov_ref: f64,\n    el_bas_sat: f64,\n    el_bas_ref: f64,\n    snr_a: f64,\n    snr_b: f64,\n) -> (f64, f64) {"
)

# Update caller
text = text.replace(
    "let (var_factor, ref_var_factor) = compute_variance_factors(ctx, el_rov_sat, el_rov_ref, el_bas_sat, el_bas_ref);",
    "let (var_factor, ref_var_factor) = compute_variance_factors(ctx, el_rov_sat, el_rov_ref, el_bas_sat, el_bas_ref, env.tuning.snr_a, env.tuning.snr_b);"
)

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(text)

with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "r") as f:
    text = f.read()

# Revert to old import
text = text.replace("crate::engine::measurement_math::build_dense_covariance_matrix", "crate::engine::measurement::build_dense_covariance_matrix")

with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "w") as f:
    f.write(text)

