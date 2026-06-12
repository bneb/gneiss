with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "r") as f:
    content = f.read()

# remove trailing }
content = content.replace("    }\n}\n\n    #[test]", "    }\n\n    #[test]")
content += "\n}\n"

with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "w") as f:
    f.write(content)
