import sys

filename = "crates/gneiss-rtk/src/engine/mod.rs"
with open(filename, "r") as f:
    content = f.read()

target = """        } else {
            // SPP failed for this epoch.
            if matches!(self.config.mode, EngineMode::Spp) {
                // If pure SPP, do not coast on dead-reckoning. Drop the state.
                self.current_state = None;
                return Err(EngineError::InitialSppFailed);
            }
        }"""

replacement = """        } else {
            // SPP failed for this epoch.
            if matches!(self.config.mode, EngineMode::Spp) {
                // We preserve the state so the next epoch has a good seed, 
                // but we return an error so no output is produced for this epoch.
                if let Some(state) = self.current_state.as_mut() {
                    state.time = rover_obs.time;
                    state.position.epoch = rover_obs.time;
                }
                return Err(EngineError::InitialSppFailed);
            }
        }"""

if target in content:
    content = content.replace(target, replacement, 1)

with open(filename, "w") as f:
    f.write(content)
