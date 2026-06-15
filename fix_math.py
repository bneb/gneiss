with open("crates/gneiss-rtk/src/engine/updater_math.rs", "r") as f:
    text = f.read()

# Delete the second populate_loosely_coupled_jacobian
lines = text.split('\n')
for i, line in enumerate(lines):
    if line.startswith("pub fn populate_loosely_coupled_jacobian(") and i > 200:
        brace_count = 0
        started = False
        end_idx = -1
        for j in range(i, len(lines)):
            brace_count += lines[j].count('{') - lines[j].count('}')
            if '{' in lines[j] or '}' in lines[j]:
                started = True
            if started and brace_count == 0:
                end_idx = j
                break
        if end_idx != -1:
            del lines[i:end_idx+1]
            break

with open("crates/gneiss-rtk/src/engine/updater_math.rs", "w") as f:
    f.write('\n'.join(lines))
