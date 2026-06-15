import re

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    lines = f.readlines()

functions_to_delete = [
    "get_sat_state",
    "compute_atmospheric_delays",
    "range_attitude_jacobian",
    "doppler_attitude_jacobian",
    "compute_variance_factors",
    "compute_zwd_mapping",
    "compute_geometric_dd",
]

lines_to_keep = []
i = 0
while i < len(lines):
    line = lines[i]
    deleted = False
    for fn in functions_to_delete:
        if line.startswith(f"pub fn {fn}(") or line.startswith(f"fn {fn}("):
            # Start of function! We need to skip until the brace count returns to 0
            # Wait, the '{' might be on a subsequent line
            brace_count = 0
            started = False
            
            while i < len(lines):
                brace_count += lines[i].count('{') - lines[i].count('}')
                if '{' in lines[i] or '}' in lines[i]:
                    started = True
                
                i += 1
                
                if started and brace_count == 0:
                    deleted = True
                    break
            
            if deleted:
                break
    
    if not deleted:
        if i < len(lines):
            lines_to_keep.append(lines[i])
        i += 1

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.writelines(lines_to_keep)

print("Done deleting.")
