with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    text = f.read()

text = text.replace(
    "let (var_factor, ref_var_factor) = compute_variance_factors(&ctx, el_rov_sat, el_rov_ref, el_bas_sat, el_bas_ref);",
    "let (var_factor, ref_var_factor) = compute_variance_factors(ctx.rov_sat.snr, ctx.rov_ref.snr, el_rov_sat, el_rov_ref, el_bas_sat, el_bas_ref, env.tuning.snr_a, env.tuning.snr_b);"
)

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(text)
