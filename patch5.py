with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "r") as f:
    content = f.read()

content = content.replace("    }\n}\n    #[test]", "    }\n    #[test]")

with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "w") as f:
    f.write(content)
