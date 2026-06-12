import re

content = open('crates/gneiss-rtk/src/engine/measurement.rs').read()

# Replace pr_l1 checks
content = content.replace("if rov_sat.pr_l1 > 0.0 && base_sat.pr_l1 > 0.0 && rov_ref.pr_l1 > 0.0 && ref_base.pr_l1 > 0.0 {", "let pr1_valid = [rov_sat.pr_l1, base_sat.pr_l1, rov_ref.pr_l1, ref_base.pr_l1].iter().all(|&x| x > 0.0);\n    if pr1_valid {")

# Replace pr_l2 checks
content = content.replace("if let (Some(rr2), Some(rs2), Some(br2), Some(bs2)) = (rov_ref.pr_l2, rov_sat.pr_l2, ref_base.pr_l2, base_sat.pr_l2) {", "let pr2_vals = [rov_ref.pr_l2, rov_sat.pr_l2, ref_base.pr_l2, base_sat.pr_l2];\n    if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = pr2_vals {")

# Replace cp_l1 checks
content = content.replace("if let (Some(rr1), Some(rs1), Some(br1), Some(bs1)) = (rov_ref.cp_l1, rov_sat.cp_l1, ref_base.cp_l1, base_sat.cp_l1) {", "let cp1_vals = [rov_ref.cp_l1, rov_sat.cp_l1, ref_base.cp_l1, base_sat.cp_l1];\n        if let [Some(rr1), Some(rs1), Some(br1), Some(bs1)] = cp1_vals {")

# Replace cp_l2 checks
content = content.replace("if let (Some(rr2), Some(rs2), Some(br2), Some(bs2)) = (rov_ref.cp_l2, rov_sat.cp_l2, ref_base.cp_l2, base_sat.cp_l2) {", "let cp2_vals = [rov_ref.cp_l2, rov_sat.cp_l2, ref_base.cp_l2, base_sat.cp_l2];\n        if let [Some(rr2), Some(rs2), Some(br2), Some(bs2)] = cp2_vals {")

# Replace doppler checks
content = content.replace("if rov_sat.doppler != 0.0 && rov_ref.doppler != 0.0 && base_sat.doppler != 0.0 && ref_base.doppler != 0.0 {", "let dop_valid = [rov_sat.doppler, rov_ref.doppler, base_sat.doppler, ref_base.doppler].iter().all(|&x| x != 0.0);\n    if dop_valid {")

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
