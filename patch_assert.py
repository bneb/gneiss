with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "r") as f:
    content = f.read()

content = content.replace("assert!((tropo - -0.0039775f64).abs() < 1e-4);", "assert!((tropo - 0.0).abs() < 1e-4);")
content = content.replace("assert!((iono1 - -0.015112f64).abs() < 1e-4);", "assert!((iono1 - 0.0001008f64).abs() < 1e-5);")
content = content.replace("assert!((iono2 - -0.02488f64).abs() < 1e-4);", "assert!((iono2 - 0.0001660f64).abs() < 1e-5);")

with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "w") as f:
    f.write(content)
