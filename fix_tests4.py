import re

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    content = f.read()

# Fix dummy state size 15->21, 17->23
content = content.replace("DMatrix::identity(15, 15)", "DMatrix::identity(21, 21)")
content = content.replace("DMatrix::identity(17, 17)", "DMatrix::identity(23, 23)")
content = content.replace("17, 0.01", "23, 0.01")
content = content.replace("vec![0.0; 17]", "vec![0.0; 23]")
content = content.replace("h[15]", "h[21]")
content = content.replace("h[16]", "h[22]")

# In the test test_validate_measurements
content = content.replace("h1[15] = 1.0;", "h1[21] = 1.0;")
content = content.replace("h2[15] = 1.0;", "h2[21] = 1.0;")

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(content)
print("Fixed test dimensions")
