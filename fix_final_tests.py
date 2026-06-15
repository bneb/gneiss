import re

# Fix ppp_fg.rs test
with open('crates/gneiss-rtk/src/engine/ppp_fg.rs', 'r') as f:
    text = f.read()

text = text.replace("assert!((inv2[(0,0)] - 1e6).abs() < 1e-1);", "assert!((inv2[(0,0)] - 1e-6).abs() < 1e-10);")

with open('crates/gneiss-rtk/src/engine/ppp_fg.rs', 'w') as f:
    f.write(text)

# Fix tests_updater.rs (by fixing the math code to actually do symmetry)
with open('crates/gneiss-rtk/src/math/covariance.rs', 'r') as f:
    text = f.read()

text = text.replace("&i_kh * p * i_kh.transpose() + k * r * k.transpose()",
                    "let mut p_new = &i_kh * p * i_kh.transpose() + k * r * k.transpose();\n    p_new.symmetric_part()")
text = text.replace("*p = &i_kh * &*p * i_kh.transpose() + &k_i * r_i * k_i.transpose();",
                    "*p = &i_kh * &*p * i_kh.transpose() + &k_i * r_i * k_i.transpose();\n    *p = p.symmetric_part();")

with open('crates/gneiss-rtk/src/math/covariance.rs', 'w') as f:
    f.write(text)

