import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    lines = f.readlines()

in_func = False
func_name = ""
func_len = 0
bracket_count = 0

for line in lines:
    if line.strip().startswith("pub fn ") or line.strip().startswith("fn "):
        in_func = True
        func_name = line.strip().split("(")[0]
        func_len = 0
        bracket_count = 0
    if in_func:
        func_len += 1
        bracket_count += line.count("{") - line.count("}")
        if bracket_count == 0 and "{" in line:
            # wait, if it's on the same line, e.g. fn foo() {}
            # it will end
            pass
        if bracket_count == 0 and func_len > 1:
            if func_len > 30:
                print(f"{func_name}: {func_len} lines")
            in_func = False

