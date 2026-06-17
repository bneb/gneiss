import re

def count_loc(file_path):
    with open(file_path, "r") as f:
        content = f.read()

    # Find all function blocks
    blocks = re.split(r'\n(?:pub\s+)?fn\s+', content)
    for i, block in enumerate(blocks):
        if i == 0:
            continue
        
        # Name
        name_match = re.match(r'([a-zA-Z0-9_]+)', block)
        if not name_match:
            continue
        name = name_match.group(1)

        # Count lines until the matching brace closes
        lines = block.split('\n')
        count = 0
        brace_level = 0
        started = False
        
        for line in lines:
            if not started and '{' in line:
                started = True
                brace_level += line.count('{') - line.count('}')
                count += 1
            elif started:
                brace_level += line.count('{') - line.count('}')
                if line.strip() and not line.strip().startswith('//'):
                    count += 1
                if brace_level <= 0:
                    break
        
        print(f"{file_path}: {name} -> {count} LOC")

count_loc("crates/gneiss-rtk/src/engine/ppp.rs")
count_loc("crates/gneiss-rtk/src/engine/ppp_math.rs")
count_loc("crates/gneiss-rtk/src/engine/ppp_fg.rs")
