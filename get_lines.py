import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    lines = f.readlines()

def count_lines(func_name):
    start = -1
    for i, line in enumerate(lines):
        if line.startswith(f"pub fn {func_name}") or line.startswith(f"fn {func_name}") or line.startswith(f"#[allow(clippy::too_many_arguments)]\nfn {func_name}") or line.startswith(f"#[allow(clippy::too_many_arguments)]\npub fn {func_name}"):
            start = i
            break
        if func_name in line and "fn " in line:
            start = i
            break
    
    if start == -1:
        return -1
    
    braces = 0
    in_func = False
    for i in range(start, len(lines)):
        line = lines[i]
        braces += line.count('{')
        braces -= line.count('}')
        if braces > 0:
            in_func = True
        if in_func and braces == 0:
            return i - start + 1
    return -1

funcs = [
    "compute_innovations",
    "process_single_satellite_pair",
    "build_measurement_model",
    "select_reference_satellite",
    "compute_dd_carrier_phase"
]

for f in funcs:
    print(f"{f}: {count_lines(f)}")

