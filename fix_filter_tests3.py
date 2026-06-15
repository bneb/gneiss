import re

path = "crates/gneiss-rtk/src/estimators/ekf/filter_tests.rs"
with open(path, "r") as f:
    text = f.read()

text = text.replace("state.ambiguities()", "state.ambiguities")
text = text.replace("state.ambiguity_keys()", "state.ambiguity_keys")

with open(path, "w") as f:
    f.write(text)
