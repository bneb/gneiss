import os
import glob

def resolve_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    # For mod.rs files, engine_arch often added structs, we want both or engine_arch's additions.
    # Actually, if we just split by conflict markers:
    parts = []
    lines = content.splitlines(True)
    in_conflict = False
    head_lines = []
    their_lines = []
    state = 0 # 0: normal, 1: HEAD, 2: theirs

    res = []
    for line in lines:
        if line.startswith("<<<<<<< HEAD"):
            state = 1
            head_lines = []
            their_lines = []
        elif line.startswith("======="):
            state = 2
        elif line.startswith(">>>>>>> engine_arch"):
            state = 0
            # Custom logic based on filename
            if "mod.rs" in filepath:
                # Keep both! Just ensure no duplicates
                merged = head_lines + [l for l in their_lines if l not in head_lines]
                res.extend(merged)
            else:
                # For updater.rs, keep HEAD (math_specialist's clean math)
                res.extend(head_lines)
        else:
            if state == 0:
                res.append(line)
            elif state == 1:
                head_lines.append(line)
            elif state == 2:
                their_lines.append(line)

    with open(filepath, 'w') as f:
        f.writelines(res)

for f in ["crates/gneiss-rtk/src/ambiguity/mod.rs",
          "crates/gneiss-rtk/src/estimators/ekf/mod.rs",
          "crates/gneiss-rtk/src/estimators/mod.rs",
          "crates/gneiss-rtk/src/measurements/mod.rs",
          "crates/gneiss-rtk/src/engine/updater.rs"]:
    resolve_file(f)
    os.system(f"git add {f}")

