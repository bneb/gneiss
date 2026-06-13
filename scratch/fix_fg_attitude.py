with open("crates/gneiss-rtk/src/engine/tight_fg.rs", "r") as f:
    content = f.read()

# Replace body frame rotation with global frame rotation
content = content.replace("state.attitude *= nalgebra::UnitQuaternion::from_scaled_axis(rot_vec);", "state.attitude = nalgebra::UnitQuaternion::from_scaled_axis(rot_vec) * state.attitude;")

with open("crates/gneiss-rtk/src/engine/tight_fg.rs", "w") as f:
    f.write(content)

print("Fixed tight_fg.rs rotation order")
