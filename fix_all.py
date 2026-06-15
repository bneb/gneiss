import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

# Fix push
content = content.replace(
    "self.z.push(u.0); self.h.push(u.1); self.r.push(u.2); self.mt.push((sat, u.3, u.4));",
    "self.z.push(u.z); self.h.push(u.h); self.r.push(u.r); self.mt.push((sat, u.type_code, u.r_ref));"
)

# Fix test
content = content.replace("u.0.is_nan()", "u.z.is_nan()")

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)

