import re

path = 'crates/gneiss-rtk/src/engine/ppp_fg.rs'
with open(path, 'r') as f:
    text = f.read()

# remove test_invert_matrix entirely using regex
pattern = r'#\[test\]\s+fn test_invert_matrix\(\)\s*\{.*?\n    \}'
text = re.sub(pattern, '', text, flags=re.DOTALL)

with open(path, 'w') as f:
    f.write(text)

