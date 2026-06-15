import re

path = "crates/gneiss-rtk/src/estimators/ekf/filter_tests.rs"
with open(path, "r") as f:
    text = f.read()

text = re.sub(r'cp_l1:\s*0\.0', 'cp_l1: Some(0.0)', text)
text = re.sub(r'locktime:\s*1000', 'locktime: Some(1000)', text)

with open(path, "w") as f:
    f.write(text)
