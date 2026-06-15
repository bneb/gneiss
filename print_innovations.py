with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    lines = f.readlines()
for i, l in enumerate(lines):
    if "pub fn compute_innovations" in l:
        for j in range(i, i+60):
            print(f"{j+1}: {lines[j]}", end="")
        break
