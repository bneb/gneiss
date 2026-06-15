import sys

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

bad_str = """            let lhs = j_linearized[(r, c)];
            let rhs = j_numeric[(r, c)];
            assert!(
                (lhs - rhs).abs() < 1e-14,"""

good_str = """            let lhs: f64 = j_linearized[(r, c)];
            let rhs: f64 = j_numeric[(r, c)];
            assert!(
                (lhs - rhs).abs() < 1e-14,"""

content = content.replace(bad_str, good_str)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
