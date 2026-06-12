with open("crates/gneiss-rtk/src/spp.rs", "r") as f:
    lines = f.readlines()

for i in range(len(lines)-1, -1, -1):
    if lines[i].strip() == "}":
        with open("test_append.rs", "r") as tf:
            test_content = tf.read()
        lines.insert(i, test_content + "\n")
        break

with open("crates/gneiss-rtk/src/spp.rs", "w") as f:
    f.writelines(lines)
