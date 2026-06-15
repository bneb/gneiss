import re

path = "crates/gneiss-rtk/src/engine/updater_math.rs"
with open(path, "r") as f:
    text = f.read()

# Add #[test] before any fn test_ if it is not immediately preceded by #[test]
lines = text.split("\n")
out = []
for i, line in enumerate(lines):
    if "fn test_" in line:
        # Check previous non-empty line
        prev_idx = i - 1
        has_test = False
        while prev_idx >= 0 and lines[prev_idx].strip() == "":
            prev_idx -= 1
        if prev_idx >= 0 and "#[test]" in lines[prev_idx]:
            has_test = True
        if not has_test:
            # maintain indentation
            indent = line[:len(line) - len(line.lstrip())]
            out.append(f"{indent}#[test]")
    out.append(line)

with open(path, "w") as f:
    f.write("\n".join(out))
