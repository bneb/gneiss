import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

struct_def = """pub struct SingleUpdate {
    pub z: f64,
    pub h: Vec<f64>,
    pub r: f64,
    pub type_code: u8,
    pub r_ref: f64,
}

"""

# We'll put SingleUpdate right before compute_pseudorange_update
content = content.replace("fn compute_pseudorange_update(", struct_def + "fn compute_pseudorange_update(")

content = content.replace("-> (f64, Vec<f64>, f64, u8, f64) {", "-> SingleUpdate {")
content = content.replace("-> Vec<(f64, Vec<f64>, f64, u8, f64)> {", "-> Vec<SingleUpdate> {")

# compute_pseudorange_update returns:
content = content.replace("(pr_dd - (comp_pr_dd + iono_dd), h_pr, r_val, 0, r_ref_val)", "SingleUpdate { z: pr_dd - (comp_pr_dd + iono_dd), h: h_pr, r: r_val, type_code: 0, r_ref: r_ref_val }")

# compute_carrier_phase_update returns:
content = content.replace("(cp_dd - (comp_pr_dd - iono_dd + n_dd * lam_sat), h_cp, r_val, freq_idx, r_ref_val)", "SingleUpdate { z: cp_dd - (comp_pr_dd - iono_dd + n_dd * lam_sat), h: h_cp, r: r_val, type_code: freq_idx, r_ref: r_ref_val }")

# compute_doppler_innovation returns:
content = content.replace("(dd_doppler - comp_doppler, h_dop, r_val, freq_idx, r_ref_val)", "SingleUpdate { z: dd_doppler - comp_doppler, h: h_dop, r: r_val, type_code: freq_idx, r_ref: r_ref_val }")

# In compute_dd_components, where they are consumed:
content = content.replace("for (z, h, r, t, r_ref) in updates {", "for up in updates {")
content = content.replace("out.z_all.push(z);", "out.z_all.push(up.z);")
content = content.replace("out.h_all.push(h);", "out.h_all.push(up.h);")
content = content.replace("out.r_all.push(r);", "out.r_all.push(up.r);")
content = content.replace("out.type_all.push((ctx.rov_sat.sat, t, r_ref));", "out.type_all.push((ctx.rov_sat.sat, up.type_code, up.r_ref));")

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)

