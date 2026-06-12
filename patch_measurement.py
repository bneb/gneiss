import sys

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    text = f.read()

# Update compute_atmospheric_delays signature
text = text.replace(
    "sat_f1: f64, sat_f2: f64,\n    ref_f1: f64, ref_f2: f64,\n) -> (f64, f64, f64)",
    "sat_f1: f64, sat_f2: f64, sat_f5: f64,\n    ref_f1: f64, ref_f2: f64, ref_f5: f64,\n) -> (f64, f64, f64, f64)"
)

# Update return of compute_atmospheric_delays
idx = text.find("    (tropo_dd, iono_dd_l1, iono_dd_l2)\n}")
if idx != -1:
    l5_iono = """
    let f_ratio_sat_l5 = (sat_f1 / sat_f5).powi(2);
    let f_ratio_ref_l5 = (ref_f1 / ref_f5).powi(2);
    let iono_dd_l5 = (iono_rov_sat * f_ratio_sat_l5 - iono_rov_ref * f_ratio_ref_l5) - (iono_bas_sat * f_ratio_sat_l5 - iono_bas_ref * f_ratio_ref_l5);
"""
    text = text[:idx] + l5_iono + "    (tropo_dd, iono_dd_l1, iono_dd_l2, iono_dd_l5)\n}" + text[idx + len("    (tropo_dd, iono_dd_l1, iono_dd_l2)\n}"):]

# Duplicate L2 PR block
l2_pr_start = "if let (Some(rr2), Some(rs2), Some(br2), Some(bs2)) = (rov_ref.pr_l2, rov_sat.pr_l2, ref_base.pr_l2, base_sat.pr_l2) {"
l2_pr_end = "updates.push((pr_dd_l2 - (comp_pr_dd + iono_dd_l2), h_pr2, pr_r_val, 2));\n    }"
idx_start = text.find(l2_pr_start)
idx_end = text.find(l2_pr_end)
if idx_start != -1 and idx_end != -1:
    l2_pr_block = text[idx_start:idx_end + len(l2_pr_end)]
    l5_pr_block = l2_pr_block.replace("rr2", "rr5").replace("rs2", "rs5").replace("br2", "br5").replace("bs2", "bs5")
    l5_pr_block = l5_pr_block.replace("pr_l2", "pr_l5").replace("pr_dd_l2", "pr_dd_l5").replace("iono_dd_l2", "iono_dd_l5").replace("h_pr2", "h_pr5")
    l5_pr_block = l5_pr_block.replace("pr_r_val, 2", "pr_r_val, 5")
    
    text = text[:idx_end + len(l2_pr_end)] + "\n\n    " + l5_pr_block + text[idx_end + len(l2_pr_end):]

# Update compute_dd_pseudorange signature
text = text.replace(
    "iono_dd_l2: f64,\n    h_r: Vector3<f64>,",
    "iono_dd_l2: f64,\n    iono_dd_l5: f64,\n    h_r: Vector3<f64>,"
)

# Duplicate L2 CP block
l2_cp_start = "if let (Some(sat_idx), Some(ref_idx)) = (state.ambiguity_keys.iter().position(|&(s, f)| s == rov_sat.sat && f == 2), ref_idx_l2) {"
l2_cp_end = "updates.push((cp_dd_l2 - (comp_pr_dd - iono_dd_l2 + n_dd_l2), h_cp2, r_val, 2));\n        }\n    }"
idx_start = text.find(l2_cp_start)
idx_end = text.find(l2_cp_end)
if idx_start != -1 and idx_end != -1:
    l2_cp_block = text[idx_start:idx_end + len(l2_cp_end)]
    l5_cp_block = l2_cp_block.replace("f == 2", "f == 5").replace("ref_idx_l2", "ref_idx_l5")
    l5_cp_block = l5_cp_block.replace("rr2", "rr5").replace("rs2", "rs5").replace("br2", "br5").replace("bs2", "bs5")
    l5_cp_block = l5_cp_block.replace("cp_l2", "cp_l5").replace("cp_dd_l2", "cp_dd_l5").replace("n_dd_l2", "n_dd_l5")
    l5_cp_block = l5_cp_block.replace("lam_ref_2", "lam_ref_5").replace("lam_sat_2", "lam_sat_5")
    l5_cp_block = l5_cp_block.replace("ref_f2", "ref_f5").replace("sat_f2", "sat_f5")
    l5_cp_block = l5_cp_block.replace("iono_dd_l2", "iono_dd_l5").replace("h_cp2", "h_cp5")
    l5_cp_block = l5_cp_block.replace("r_val, 2", "r_val, 5")
    
    text = text[:idx_end + len(l2_cp_end)] + "\n\n    " + l5_cp_block + text[idx_end + len(l2_cp_end):]

# Update compute_dd_carrier_phase signature
text = text.replace(
    "ref_idx_l2: Option<usize>,\n    comp_pr_dd: f64,\n    iono_dd_l1: f64,\n    iono_dd_l2: f64,\n    sat_f1: f64, sat_f2: f64,\n    ref_f1: f64, ref_f2: f64,",
    "ref_idx_l2: Option<usize>,\n    ref_idx_l5: Option<usize>,\n    comp_pr_dd: f64,\n    iono_dd_l1: f64,\n    iono_dd_l2: f64,\n    iono_dd_l5: f64,\n    sat_f1: f64, sat_f2: f64, sat_f5: f64,\n    ref_f1: f64, ref_f2: f64, ref_f5: f64,"
)

# Update build_measurement_model calls
text = text.replace(
    "let [sat_f1, sat_f2, _] = gneiss_core::signal::satellite_frequencies(rov_sat.sat, sat_freq_num);\n            let [ref_f1, ref_f2, _] = gneiss_core::signal::satellite_frequencies(ref_rover.sat, ref_freq_num);",
    "let [sat_f1, sat_f2, sat_f5] = gneiss_core::signal::satellite_frequencies(rov_sat.sat, sat_freq_num);\n            let [ref_f1, ref_f2, ref_f5] = gneiss_core::signal::satellite_frequencies(ref_rover.sat, ref_freq_num);"
)

text = text.replace(
    "let (tropo_dd, iono_dd_l1, iono_dd_l2) = compute_atmospheric_delays(\n                state.time, pos_apc, base_coord_vec, sat_vec_rov, ref_sat_vec_rov, sat_vec_bas, ref_sat_vec_bas, \n                sat_f1, sat_f2, ref_f1, ref_f2\n            );",
    "let (tropo_dd, iono_dd_l1, iono_dd_l2, iono_dd_l5) = compute_atmospheric_delays(\n                state.time, pos_apc, base_coord_vec, sat_vec_rov, ref_sat_vec_rov, sat_vec_bas, ref_sat_vec_bas, \n                sat_f1, sat_f2, sat_f5, ref_f1, ref_f2, ref_f5\n            );"
)

text = text.replace(
    "let ref_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover.sat && f == 2);",
    "let ref_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover.sat && f == 2);\n            let ref_idx_l5 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover.sat && f == 5);"
)

text = text.replace(
    "compute_dd_pseudorange(&rov_sat, &bas_sat, &rov_ref, &bas_ref, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, state_size, var_factor)",
    "compute_dd_pseudorange(&rov_sat, &bas_sat, &rov_ref, &bas_ref, comp_pr_dd, iono_dd_l1, iono_dd_l2, iono_dd_l5, h_r, h_att, state_size, var_factor)"
)

text = text.replace(
    "compute_dd_carrier_phase(state, &rov_sat, &bas_sat, &rov_ref, &bas_ref, Some(ref_idx_l1), ref_idx_l2, comp_pr_dd, iono_dd_l1, iono_dd_l2, sat_f1, sat_f2, ref_f1, ref_f2, h_r, h_att, state_size, var_factor)",
    "compute_dd_carrier_phase(state, &rov_sat, &bas_sat, &rov_ref, &bas_ref, Some(ref_idx_l1), ref_idx_l2, ref_idx_l5, comp_pr_dd, iono_dd_l1, iono_dd_l2, iono_dd_l5, sat_f1, sat_f2, sat_f5, ref_f1, ref_f2, ref_f5, h_r, h_att, state_size, var_factor)"
)

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(text)
print("Patched measurement.rs")
