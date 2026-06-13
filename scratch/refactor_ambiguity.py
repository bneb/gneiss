import re

with open("crates/gneiss-rtk/src/filter.rs", "r") as f:
    content = f.read()

# We want to move resolve_ambiguities from filter.rs to engine/ambiguity.rs
# But I can do this directly using python script.
