import sys

with open("crates/gneiss-rtk/src/engine/tight_fg.rs", "r") as f:
    text = f.read()

# I will write a better replacer
def replace_func(match):
    return "fn compute_pcv(sat_obs: &SatObs, engine: &ProcessingEngine, rover_obs: &EpochObs, sat_pos: Vector3<f64>, b1: u8, rcv_pos: Vector3<f64>, sat_pos_rot: &mut Vector3<f64>) -> f64 { 0.0 }"

import re
text = re.sub(r"fn compute_pcv.*?-> f64 \{.*?^}", replace_func, text, flags=re.DOTALL | re.MULTILINE)

with open("crates/gneiss-rtk/src/engine/tight_fg.rs", "w") as f:
    f.write(text)
