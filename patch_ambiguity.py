import sys

with open('crates/gneiss-rtk/src/engine/ambiguity.rs', 'r') as f:
    content = f.read()

# Duplicate the block for L2 to create L5
l2_start = content.find("if let (Some(r_pr2), Some(r_cp2), Some(b_pr2), Some(b_cp2)) = (r.pr_l2, r.cp_l2, b.pr_l2, b.cp_l2) {")
l2_end = content.find("if let (Some(r_pr2), Some(r_cp2), Some(b_pr2), Some(b_cp2), Some(r_cp1), Some(b_cp1))")

if l2_start != -1 and l2_end != -1:
    l2_block = content[l2_start:l2_end]
    l5_block = l2_block.replace("r_pr2", "r_pr5").replace("r_cp2", "r_cp5").replace("b_pr2", "b_pr5").replace("b_cp2", "b_cp5")
    l5_block = l5_block.replace("pr_l2", "pr_l5").replace("cp_l2", "cp_l5")
    l5_block = l5_block.replace("lam_r2", "lam_r5").replace("lam_b2", "lam_b5")
    l5_block = l5_block.replace("cp_l2_rov", "cp_l5_rov").replace("cp_l2_base", "cp_l5_base")
    l5_block = l5_block.replace("initial_est_l2", "initial_est_l5")
    l5_block = l5_block.replace("r_f2", "r_f5").replace("b_f2", "b_f5")
    l5_block = l5_block.replace("f == 2", "f == 5")
    l5_block = l5_block.replace("r.sat, 2", "r.sat, 5")
    l5_block = l5_block.replace("slip_l2", "slip_l5")
    
    new_content = content[:l2_end] + l5_block + content[l2_end:]
    
    # Need to add mut slip_l5
    new_content = new_content.replace("let mut slip_l2 = slip;", "let mut slip_l2 = slip;\n        let mut slip_l5 = slip;")
    
    with open('crates/gneiss-rtk/src/engine/ambiguity.rs', 'w') as f:
        f.write(new_content)
    print("Patched ambiguity.rs")
else:
    print("Could not find L2 block")
