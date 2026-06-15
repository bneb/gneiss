import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

content = content.replace(") -> Option<(f64, Vec<f64>, f64, u8, f64)> {", ") -> Option<SingleUpdate> {")

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)

