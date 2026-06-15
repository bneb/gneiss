import re

path = "crates/gneiss-rtk/src/estimators/ekf/filter_tests.rs"
with open(path, "r") as f:
    text = f.read()

text = text.replace("state.epoch_count() = 100;", "state.epoch_count = 100;")

with open(path, "w") as f:
    f.write(text)
