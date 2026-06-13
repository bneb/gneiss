with open("crates/gneiss-rtk/src/engine/updater.rs", "r") as f:
    content = f.read()

# Replace body frame rotation with global frame rotation
content = content.replace("state.attitude = state.attitude * dq;", "state.attitude = dq * state.attitude;")

with open("crates/gneiss-rtk/src/engine/updater.rs", "w") as f:
    f.write(content)

print("Fixed apply_state_correction rotation order in updater.rs")
