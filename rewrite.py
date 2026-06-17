import sys

with open("crates/gneiss-rtk/src/engine/ppp.rs", "r") as f:
    ppp_content = f.read()

# Instead of trying to parse it with python, I will just copy the entire ppp.rs content, but chunked into ppp.rs and ppp_math.rs.
# It is actually easier to use standard text editor replacements if I know the line bounds.
# build_sats is from lines 52 to 352
