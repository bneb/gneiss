import sys

filename = "crates/gneiss-rtk/src/engine/mod.rs"
with open(filename, "r") as f:
    content = f.read()

target = """                // If the jump implies a velocity > 150 m/s (540 km/h), reject it
                if dt > 0.0 && (dist / dt) > 150.0 {
                    tracing::warn!("SPP rejected due to impossible velocity: {:.1} m/s", dist / dt);
                    self.current_state = None;
                    return Err(EngineError::InitialSppFailed);
                }"""

replacement = """                // If the jump implies a velocity > 150 m/s (540 km/h), reject it
                if dt > 0.0 && (dist / dt) > 150.0 {
                    tracing::warn!("SPP rejected due to impossible velocity: {:.1} m/s", dist / dt);
                    // DO NOT wipe the state! Just advance the time so we can coast and seed the next epoch properly.
                    state.time = rover_obs.time;
                    state.position.epoch = rover_obs.time;
                    return Err(EngineError::InitialSppFailed);
                }"""

if target in content:
    content = content.replace(target, replacement, 1)

with open(filename, "w") as f:
    f.write(content)
