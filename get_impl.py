with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    text = f.read()

lines = text.split('\n')
for i, line in enumerate(lines):
    if "pub fn compute_innovations" in line:
        brace_count = 0
        for j in range(i, len(lines)):
            brace_count += lines[j].count('{') - lines[j].count('}')
            if brace_count == 0:
                print(f"End of compute_innovations at line {j+1}")
                break
        break
