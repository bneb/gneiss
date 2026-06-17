import re

with open("crates/gneiss-rtk/src/engine/ppp.rs", "r") as f:
    orig = f.read()

# We will just write out ppp_math.rs and completely replace ppp.rs.
# To keep this clean, let's use a Python script to define the precise strings.

# Wait, `apply_osb_corrections` and `OsbCorrections` need to be moved to ppp_math.rs.
# The cycle slip detection needs to be moved to ppp_math.rs.
