import re

with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "r") as f:
    content = f.read()

content = content.replace("Vector3::new(20", "nalgebra::Vector3::new(20")
content = content.replace("Vector3::new(10", "nalgebra::Vector3::new(10")
content = content.replace("Vector3::new(0", "nalgebra::Vector3::new(0")
content = content.replace("tropo - -0.0039775", "tropo - -0.0039775f64")
content = content.replace("iono1 - -0.015112", "iono1 - -0.015112f64")
content = content.replace("iono2 - -0.02488", "iono2 - -0.02488f64")

with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "w") as f:
    f.write(content)
