with open("crates/gneiss-rtk/src/engine/mod.rs", "r") as f:
    text = f.read()

old_block = """            if crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat, self.config.spp_consistency_threshold_m, None, self.config.mode.is_tightly_coupled(), &self.config.tuning).is_err() {
                state.consecutive_rejections += 1;"""

new_block = """            match crate::engine::updater::update(state, &z_vec, &h_mat, &r_mat, self.config.spp_consistency_threshold_m, None, self.config.mode.is_tightly_coupled(), &self.config.tuning) {
                Err(_) | Ok(valid) if valid.len() < z_vec.len() => {
                    state.consecutive_rejections += 1;"""

old_block_end = """                }
            } else {
                state.consecutive_rejections = 0;
            }"""

new_block_end = """                }
                _ => {
                    state.consecutive_rejections = 0;
                }
            }"""

if old_block in text and old_block_end in text:
    print("Found old blocks")
    text = text.replace(old_block, new_block)
    text = text.replace(old_block_end, new_block_end)
    with open("crates/gneiss-rtk/src/engine/mod.rs", "w") as f:
        f.write(text)
else:
    print("Old blocks not found!")
