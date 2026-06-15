with open("crates/gneiss-rtk/src/engine/ambiguity.rs", "r") as f:
    text = f.read()
# Let's inspect the file first
with open("ambiguity_debug.txt", "w") as f:
    f.write(text)
