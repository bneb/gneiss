import re

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    text = f.read()

# Since `check_loc` used a naive regex, we can just print lines 212-250
# Actually let's just grep or print it with python
lines = text.split('\n')
for i, line in enumerate(lines):
    if "pub fn compute_dd_doppler" in line:
        for j in range(i, i+50):
            print(f"{j+1}: {lines[j]}")
        break

