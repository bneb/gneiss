with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    lines = f.readlines()

for i in range(796, 0, -1):
    if lines[i].startswith("impl ") or lines[i].startswith("pub fn ") or lines[i].startswith("fn "):
        print(f"Parent might be: {lines[i].strip()}")
        break
