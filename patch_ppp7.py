import os

with open("crates/gneiss-rtk/src/engine/ppp_math.rs", "r") as f:
    content = f.read()

target = """    let mut p1 = p1_obs.map(|o| o.0);
    let mut p2 = p2_obs.map(|o| o.0);
    let mut cp1 = cp1_obs.map(|o| o.0);
    let mut cp2 = cp2_obs.map(|o| o.0);"""

replacement = """    let mut p1 = p1_obs.map(|o| o.0); let mut p2 = p2_obs.map(|o| o.0); let mut cp1 = cp1_obs.map(|o| o.0); let mut cp2 = cp2_obs.map(|o| o.0);"""

content = content.replace(target, replacement)
with open("crates/gneiss-rtk/src/engine/ppp_math.rs", "w") as f:
    f.write(content)
