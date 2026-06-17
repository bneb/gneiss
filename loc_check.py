import sys
import re

def check_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    in_fn = False
    fn_name = ""
    fn_start = 0
    brace_level = 0
    max_brace = 0
    
    for i, line in enumerate(lines):
        line_clean = line.split('//')[0].strip()
        
        if line_clean.startswith('fn ') or line_clean.startswith('pub fn ') or line_clean.startswith('pub(crate) fn '):
            in_fn = True
            fn_name = line_clean.split('fn ')[1].split('(')[0].strip()
            fn_start = i
            brace_level = 0
            max_brace = 0
            
        if in_fn:
            brace_level += line_clean.count('{') - line_clean.count('}')
            if brace_level > max_brace:
                max_brace = brace_level
            
            if brace_level == 0 and '{' in line_clean and '}' in line_clean:
                # One-liner function
                in_fn = False
                continue
                
            if brace_level == 0 and i > fn_start:
                in_fn = False
                loc = i - fn_start + 1
                if loc > 30 or max_brace >= 4:
                    print(f"{filename}:{fn_start+1} {fn_name} - LOC: {loc}, Max Nesting: {max_brace}")

for f in sys.argv[1:]:
    check_file(f)
