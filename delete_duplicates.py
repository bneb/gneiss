import re

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    text = f.read()

functions_to_delete = [
    "get_sat_state",
    "compute_atmospheric_delays",
    "range_attitude_jacobian",
    "doppler_attitude_jacobian",
    "compute_variance_factors",
    "compute_zwd_mapping",
    "compute_geometric_dd"
]

lines = text.split('\n')
new_lines = []
skip = False
brace_count = 0
found_first_brace = False

for i, line in enumerate(lines):
    if not skip:
        # Check if line starts one of the functions
        match = False
        for fn in functions_to_delete:
            if line.startswith(f"pub fn {fn}(") or line.startswith(f"pub fn {fn} (") or line.startswith(f"fn {fn}("):
                match = True
                break
            # Also catch the #[allow(clippy...)] right before it
            if line.startswith("#[allow(clippy::too_many_arguments)]"):
                next_line = lines[i+1]
                if next_line.startswith(f"pub fn {fn}") or next_line.startswith(f"fn {fn}"):
                    match = True
                    break
        
        if match:
            skip = True
            brace_count = 0
            found_first_brace = False
            
    if skip:
        if '{' in line:
            found_first_brace = True
        brace_count += line.count('{') - line.count('}')
        if found_first_brace and brace_count <= 0:
            skip = False
    else:
        new_lines.append(line)

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write('\n'.join(new_lines))

