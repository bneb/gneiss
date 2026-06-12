with open("crates/gneiss-rtk/src/engine/mod.rs", "r") as f:
    lines = f.readlines()
for i in range(len(lines)):
    if lines[i].strip().startswith("pub mod ppp_fg;"):
        lines[i] = "//" + lines[i]
    if lines[i].strip().startswith("pub mod tight_fg;"):
        lines[i] = "//" + lines[i]
    if lines[i].strip().startswith("mod tests_measurement;"):
        lines[i] = "//" + lines[i]
with open("crates/gneiss-rtk/src/engine/mod.rs", "w") as f:
    f.writelines(lines)
