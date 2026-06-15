with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    lines = f.readlines()

brace_count = 0
for i in range(len(lines)):
    brace_count += lines[i].count('{') - lines[i].count('}')

print(f"Total brace count: {brace_count}")
