with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    lines = f.readlines()
for i, l in enumerate(lines):
    if "fn process_single_satellite_pair" in l:
        for j in range(i, i+80):
            print(f"{j+1}: {lines[j]}", end="")
        break
