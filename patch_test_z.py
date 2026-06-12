with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "r") as f:
    content = f.read()

content = content.replace("println!(\"Z: {:?}\", z);", "println!(\"Z len: {}, Z: {:?}\", z.len(), z);")

with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "w") as f:
    f.write(content)
