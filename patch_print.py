with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "r") as f:
    content = f.read()

content = content.replace("assert!((tropo - -0.0039775f64).abs() < 1e-4);", "println!(\"tropo: {}, iono1: {}, iono2: {}\", tropo, iono1, iono2); assert!((tropo - -0.0039775f64).abs() < 1e-4);")

with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "w") as f:
    f.write(content)
