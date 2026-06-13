import re

with open("crates/gneiss-core/src/sun.rs", "r") as f:
    content = f.read()

# Replace GMST in sun_position_ecef and moon_position_ecef
old_gmst = "let gmst = (4.8949612128230587e-5 * d * 86400.0 + 1.7533685592333) % (2.0 * core::f64::consts::PI);"
new_gmst = "let gmst = ((18.697374558 + 24.06570982441908 * d) * (core::f64::consts::PI / 12.0)) % (2.0 * core::f64::consts::PI);"

if old_gmst in content:
    content = content.replace(old_gmst, new_gmst)
else:
    print("Could not find GMST in sun.rs")

with open("crates/gneiss-core/src/sun.rs", "w") as f:
    f.write(content)
