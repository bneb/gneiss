with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    lines = f.readlines()

for i in range(710, 715):
    print(f"{i+1}: {lines[i].rstrip()}")
