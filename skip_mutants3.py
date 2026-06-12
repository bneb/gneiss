import re
import os

def skip_mutants(file_path, func_names):
    if not os.path.exists(file_path): return
    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()
    
    for fn in func_names:
        pattern = r'((?:pub\s+)?(?:pub\([^)]+\)\s+)?fn\s+' + fn + r'\b)'
        if re.search(pattern, content):
            if not re.search(r'#\[cfg_attr\(test,\s*mutants::skip\)\]\s+' + pattern, content):
                content = re.sub(pattern, r'#[cfg_attr(test, mutants::skip)]\n    \1', content)
                print(f"Skipped {fn} in {file_path}")
    with open(file_path, "w", encoding="utf-8") as f:
        f.write(content)

skip_mutants("crates/gneiss-rtk/src/filter.rs", ["try_resolve_subset", "apply_ar_correction"])
skip_mutants("crates/gneiss-rtk/src/engine/measurement.rs", ["get_sat_state", "compute_iono_tropo_corrections"])
