import re

with open("crates/gneiss-rtk/src/engine/ppp_ins_fg.rs", "r") as f:
    content = f.read()

# E0282: type annotation needed for buf
content = content.replace("|buf| buf.last()", "|buf: &Vec<gneiss_core::imu::ImuMeasurement>| buf.last()")

# E0425: missing lever_arm in solve inside process_ppp_ins_fg
# Ah! The call was `engine.imu_history.as_slice()` etc. It couldn't find engine.config because process_ppp_ins_fg had engine. Let's look at it.
# Actually I need to restore the original backup, do it cleanly in one python script.
