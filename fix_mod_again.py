import re

with open("crates/gneiss-rtk/src/engine/mod.rs", "r") as f:
    text = f.read()

# Fix the broken replacement
broken_str = "klobuchar: None,\n            precise_orbits: None,\n            precise_clocks: None,\n            antex: None, "
text = text.replace(broken_str, "")

# It was "self.config," that got replaced!
broken_full = """crate::engine::ambiguity::manage_ambiguities_and_slips(state, &self.config,
            klobuchar: None,
            precise_orbits: None,
            precise_clocks: None,
            antex: None, &matched_obs, &self.ephemerides, &base_coord, rover_obs.time, base.time);"""
fixed_full = """crate::engine::ambiguity::manage_ambiguities_and_slips(state, &self.config, &matched_obs, &self.ephemerides, &base_coord, rover_obs.time, base.time);"""
text = text.replace(broken_full, fixed_full)

# Find where ProcessingEngine is actually instantiated
# ProcessingEngine::new() -> Self
init_pattern = """        Self {
            ephemerides: Vec::new(),
            state_history: Vec::new(),
            obs_history: Vec::new(),
            smooth_history: Vec::new(),
            spp_state_history: Vec::new(),
            current_state: None,
            config: config,"""
init_fixed = """        Self {
            ephemerides: Vec::new(),
            state_history: Vec::new(),
            obs_history: Vec::new(),
            smooth_history: Vec::new(),
            spp_state_history: Vec::new(),
            current_state: None,
            config: config,
            klobuchar: None,
            precise_orbits: None,
            precise_clocks: None,
            antex: None,"""
if "precise_clocks: None" not in text:
    text = text.replace(init_pattern, init_fixed)

# Also fix the `antex` type error
# We don't have Antex parsing yet maybe? No, we had it in gneiss_parsers but it might not be `pub`.
# Let's just use `Option<()>` for antex for now?
# Or just remove it! The `cargo check` said `not found in gneiss_parsers::antex`.
# Wait, let's look at `cargo check` again. "cannot find type Antex in module gneiss_parsers::antex".
# Wait, I didn't add `pub mod antex;` to `crates/gneiss-parsers/src/lib.rs`!
# Let's fix that too.

# Write back
with open("crates/gneiss-rtk/src/engine/mod.rs", "w") as f:
    f.write(text)

with open("crates/gneiss-parsers/src/lib.rs", "r") as f:
    lib_text = f.read()

if "pub mod antex;" not in lib_text:
    lib_text += "\npub mod antex;\n"
    with open("crates/gneiss-parsers/src/lib.rs", "w") as f:
        f.write(lib_text)

print("Fixed mod.rs and lib.rs")
