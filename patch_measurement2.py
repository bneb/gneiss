import sys

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

# Update compute_atmospheric_delays
sig_start = content.find("fn compute_atmospheric_delays(")
sig_end = content.find("-> (f64, f64, f64) {")
if sig_start != -1 and sig_end != -1:
    sig_block = content[sig_start:sig_end]
    new_sig_block = sig_block.replace("sat_f1: f64, sat_f2: f64,", "sat_f1: f64, sat_f2: f64, sat_f5: f64,")
    new_sig_block = new_sig_block.replace("ref_f1: f64, ref_f2: f64,", "ref_f1: f64, ref_f2: f64, ref_f5: f64,")
    content = content[:sig_start] + new_sig_block + "-> (f64, f64, f64, f64) {\n" + content[sig_end + len("-> (f64, f64, f64) {"):]

ret_start = content.find("(tropo_dd, iono_dd_l1, iono_dd_l2)", sig_start)
if ret_start != -1:
    block_before_ret = content[sig_start:ret_start]
    l5_iono = "\n    let f_ratio_sat_l5 = (sat_f1 / sat_f5).powi(2);\n    let f_ratio_ref_l5 = (ref_f1 / ref_f5).powi(2);\n    let iono_dd_l5 = (iono_rov_sat * f_ratio_sat_l5 - iono_rov_ref * f_ratio_ref_l5) - (iono_bas_sat * f_ratio_sat_l5 - iono_bas_ref * f_ratio_ref_l5);\n\n    "
    content = content[:ret_start] + l5_iono + "(tropo_dd, iono_dd_l1, iono_dd_l2, iono_dd_l5)" + content[ret_start + len("(tropo_dd, iono_dd_l1, iono_dd_l2)"):]

# Duplicate compute_dd_pseudorange block for L2 to create L5
l2_pr_start = content.find("if let (Some(rr2), Some(rs2), Some(br2), Some(bs2)) = (rov_ref.pr_l2, rov_sat.pr_l2, ref_base.pr_l2, base_sat.pr_l2) {")
l2_pr_end = content.find("    updates\n}", l2_pr_start)
if l2_pr_start != -1 and l2_pr_end != -1:
    l2_pr_block = content[l2_pr_start:l2_pr_end]
    l5_pr_block = l2_pr_block.replace("rr2", "rr5").replace("rs2", "rs5").replace("br2", "br5").replace("bs2", "bs5")
    l5_pr_block = l5_pr_block.replace("pr_l2", "pr_l5")
    l5_pr_block = l5_pr_block.replace("pr_dd_l2", "pr_dd_l5")
    l5_pr_block = l5_pr_block.replace("iono_dd_l2", "iono_dd_l5").replace("2))", "5))")
    
    content = content[:l2_pr_end] + l5_pr_block + content[l2_pr_end:]
    
    # Update compute_dd_pseudorange signature
    sig_start = content.find("fn compute_dd_pseudorange(")
    sig_end = content.find("-> Vec<(f64, Vec<f64>, f64, u8)> {", sig_start)
    if sig_start != -1 and sig_end != -1:
        sig_block = content[sig_start:sig_end]
        new_sig_block = sig_block.replace("iono_dd_l2: f64,", "iono_dd_l2: f64,\n    iono_dd_l5: f64,")
        content = content[:sig_start] + new_sig_block + content[sig_end:]

# Update the call sites in build_measurement_model
call_atmo = content.find("let (tropo_dd, iono_dd_l1, iono_dd_l2) = compute_atmospheric_delays(")
if call_atmo != -1:
    content = content[:call_atmo] + "let (tropo_dd, iono_dd_l1, iono_dd_l2, iono_dd_l5) = compute_atmospheric_delays(" + content[call_atmo + len("let (tropo_dd, iono_dd_l1, iono_dd_l2) = compute_atmospheric_delays("):]
    
    # Need to add sat_f5 and ref_f5
    call_atmo_end = content.find(");", call_atmo)
    call_atmo_block = content[call_atmo:call_atmo_end]
    call_atmo_block = call_atmo_block.replace("sat_f1, sat_f2,", "sat_f1, sat_f2, sat_f5,")
    call_atmo_block = call_atmo_block.replace("ref_f1, ref_f2", "ref_f1, ref_f2, ref_f5")
    content = content[:call_atmo] + call_atmo_block + content[call_atmo_end:]

call_pr = content.find("for update in compute_dd_pseudorange(&rov_sat, &bas_sat, &rov_ref, &bas_ref, comp_pr_dd, iono_dd_l1, iono_dd_l2,")
if call_pr != -1:
    content = content.replace("for update in compute_dd_pseudorange(&rov_sat, &bas_sat, &rov_ref, &bas_ref, comp_pr_dd, iono_dd_l1, iono_dd_l2, h_r, h_att, state_size, var_factor) {",
                              "for update in compute_dd_pseudorange(&rov_sat, &bas_sat, &rov_ref, &bas_ref, comp_pr_dd, iono_dd_l1, iono_dd_l2, iono_dd_l5, h_r, h_att, state_size, var_factor) {")

call_cp = content.find("for update in compute_dd_carrier_phase(state, &rov_sat, &bas_sat, &rov_ref, &bas_ref, Some(ref_idx_l1), ref_idx_l2,")
if call_cp != -1:
    content = content.replace("for update in compute_dd_carrier_phase(state, &rov_sat, &bas_sat, &rov_ref, &bas_ref, Some(ref_idx_l1), ref_idx_l2, comp_pr_dd, iono_dd_l1, iono_dd_l2, sat_f1, sat_f2, ref_f1, ref_f2, h_r, h_att, state_size, var_factor) {",
                              "for update in compute_dd_carrier_phase(state, &rov_sat, &bas_sat, &rov_ref, &bas_ref, Some(ref_idx_l1), ref_idx_l2, ref_idx_l5, comp_pr_dd, iono_dd_l1, iono_dd_l2, iono_dd_l5, sat_f1, sat_f2, sat_f5, ref_f1, ref_f2, ref_f5, h_r, h_att, state_size, var_factor) {")

# Find ref_idx_l5 assignment
ref_idx_l2_assign = content.find("let ref_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover.sat && f == 2);")
if ref_idx_l2_assign != -1:
    content = content[:ref_idx_l2_assign] + "let ref_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover.sat && f == 2);\n            let ref_idx_l5 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover.sat && f == 5);" + content[ref_idx_l2_assign + len("let ref_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover.sat && f == 2);"):]

# Find frequencies assignments
freq_assign = content.find("let [sat_f1, sat_f2, _] = gneiss_core::signal::satellite_frequencies(rov_sat.sat, sat_freq_num);")
if freq_assign != -1:
    content = content.replace("let [sat_f1, sat_f2, _] = gneiss_core::signal::satellite_frequencies(rov_sat.sat, sat_freq_num);",
                              "let [sat_f1, sat_f2, sat_f5] = gneiss_core::signal::satellite_frequencies(rov_sat.sat, sat_freq_num);")
    content = content.replace("let [ref_f1, ref_f2, _] = gneiss_core::signal::satellite_frequencies(ref_rover.sat, ref_freq_num);",
                              "let [ref_f1, ref_f2, ref_f5] = gneiss_core::signal::satellite_frequencies(ref_rover.sat, ref_freq_num);")


with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
print("Patched remaining measurement.rs")
