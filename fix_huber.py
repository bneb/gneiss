import sys

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'r') as f:
    content = f.read()

# Change r to diag(3.0) so p=2, r=3. p+r = 5. p*r = 6. 
content = content.replace("let r = DMatrix::from_diagonal(&DVector::from_vec(vec![2.0, 2.0]));", "let r = DMatrix::from_diagonal(&DVector::from_vec(vec![3.0, 3.0]));")
content = content.replace("// s_raw = p + r = diag(4.0, 4.0)", "// s_raw = p + r = diag(5.0, 5.0)")
content = content.replace("// s_raw_inv = diag(0.25, 0.25)", "// s_raw_inv = diag(0.2, 0.2)")
# z = [2.0, 2.0] -> z^T * s_inv * z = 4*0.2 + 4*0.2 = 1.6 <= 4.0
content = content.replace("// z = [2.0, 2.0] -> z^T * s_inv * z = 4*0.25 + 4*0.25 = 2.0 <= 4.0", "// z = [2.0, 2.0] -> z^T * s_inv * z = 4*0.2 + 4*0.2 = 1.6 <= 4.0")
content = content.replace("assert!((scaled_r[(0, 0)] - 2.0).abs() < 1e-9); // R remains unscaled", "assert!((scaled_r[(0, 0)] - 3.0).abs() < 1e-9); // R remains unscaled")
# z = [4.0, 4.0] -> z^T * s_inv * z = 16*0.2 + 16*0.2 = 6.4 > 4.0
# scale = 6.4 / 4.0 = 1.6
# R_new = R * 1.6 = diag(4.8, 4.8)
content = content.replace("// z = [4.0, 4.0] -> z^T * s_inv * z = 16*0.25 + 16*0.25 = 8.0 > 4.0", "// z = [4.0, 4.0] -> z^T * s_inv * z = 16*0.2 + 16*0.2 = 6.4 > 4.0")
content = content.replace("// scale = 8.0 / 4.0 = 2.0", "// scale = 6.4 / 4.0 = 1.6")
content = content.replace("// R_new = R * 2.0 = diag(4.0, 4.0)", "// R_new = R * 1.6 = diag(4.8, 4.8)")
content = content.replace("assert!((scaled_r_t2[(0, 0)] - 4.0).abs() < 1e-9);", "assert!((scaled_r_t2[(0, 0)] - 4.8).abs() < 1e-9);")
content = content.replace("assert!((scaled_r_t2[(1, 1)] - 4.0).abs() < 1e-9);", "assert!((scaled_r_t2[(1, 1)] - 4.8).abs() < 1e-9);")

with open('crates/gneiss-rtk/src/engine/updater_math.rs', 'w') as f:
    f.write(content)
