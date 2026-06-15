import ast

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    source = f.read()

# Rust is not Python, so ast won't work.
# Let's write a simple brace matcher that works over characters.

in_string = False
in_comment = False
func_starts = []

lines = source.split("\n")
for i, line in enumerate(lines):
    if line.strip().startswith("pub fn ") or line.strip().startswith("fn "):
        if "test" not in line and "cfg(test)" not in line:
            func_starts.append((i, line.strip().split("(")[0]))

for start_idx, name in func_starts:
    brace_count = 0
    started = False
    for i in range(start_idx, len(lines)):
        line = lines[i]
        brace_count += line.count("{") - line.count("}")
        if "{" in line:
            started = True
        if started and brace_count <= 0:
            length = i - start_idx + 1
            if length > 30:
                print(f"{name}: {length} lines")
            break
