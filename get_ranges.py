with open("crates/gneiss-rtk/src/engine/updater.rs", "r") as f:
    lines = f.readlines()

for name in ["extract_valid_measurements", "remove_measurement"]:
    for i in range(len(lines)):
        if lines[i].startswith(f"fn {name}(") or lines[i].startswith(f"pub fn {name}("):
            brace_count = 0
            started = False
            for j in range(i, len(lines)):
                brace_count += lines[j].count('{') - lines[j].count('}')
                if '{' in lines[j] or '}' in lines[j]:
                    started = True
                if started and brace_count == 0:
                    print(f"DEL {name}: {i+1} {j+1}")
                    break
            break
