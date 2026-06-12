import re

with open("crates/gneiss-rtk/src/engine/mod.rs", "r") as f:
    text = f.read()

# Add antex back to ProcessingEngine
if "pub antex" not in text:
    text = text.replace("pub precise_clocks: Option<crate::engine::ssr::PreciseClocks>,", "pub precise_clocks: Option<crate::engine::ssr::PreciseClocks>,\n    pub antex: Option<gneiss_parsers::antex::AntexDatabase>,")
    text = text.replace("precise_clocks: None,", "precise_clocks: None,\n            antex: None,")

# Now strip out the methods: process_spp, process_rtk, process_spp_loosely_coupled, process_rtk_loosely_coupled,
# update_ekf, predict_ekf, find_common_satellites.
# Actually, the easiest way to delete them is to find their boundaries.
# Or since they are from lines ~260 to the end of `impl ProcessingEngine`,
# let's just use Python's ast? No, we can just find 'pub fn process_spp_loosely_coupled' and delete everything until 'pub fn get_state'.
# Wait! Does `mod.rs` still have `get_state`? Let's check where the impl ends.

# Let's find exactly what's duplicated.
# In `updater.rs`, we have: update_ekf, update_loosely_coupled
# In `predictor.rs`, we have: predict_ekf, find_common_satellites, predict_imu
# In `spp.rs`, we don't have `spp.rs`, we left `process_spp` in `mod.rs`!
# Wait, did I leave `process_spp` in `mod.rs` or move it to `spp.rs`?
# Let's look at the directory! I didn't see `spp.rs` in `ls -l crates/gneiss-rtk/src/engine`.
# Oh! I didn't create `spp.rs` or `rtk.rs`!
# The duplicates are ONLY `update_ekf` and `predict_ekf`?!
# Let me look at `cargo check` output for E0428.
