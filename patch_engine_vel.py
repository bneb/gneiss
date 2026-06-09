import sys

filename = "crates/gneiss-rtk/src/engine/mod.rs"
with open(filename, "r") as f:
    content = f.read()

target = """        if let Some(pos) = spp_pos {
            if matches!(self.config.mode, EngineMode::Spp) {
                // Pure SPP is an epoch-by-epoch solution. Do not filter.
                let state = self.current_state.as_mut().unwrap();
                state.time = rover_obs.time;
                state.position = pos;
                state.position.epoch = rover_obs.time;"""

replacement = """        if let Some(pos) = spp_pos {
            if matches!(self.config.mode, EngineMode::Spp) {
                // Velocity sanity check for pure SPP
                let state = self.current_state.as_mut().unwrap();
                let dx = pos.vector.x - state.position.vector.x;
                let dy = pos.vector.y - state.position.vector.y;
                let dz = pos.vector.z - state.position.vector.z;
                let dist = f64::sqrt(dx*dx + dy*dy + dz*dz);
                
                // If the jump implies a velocity > 150 m/s (540 km/h), reject it
                if dt > 0.0 && (dist / dt) > 150.0 {
                    tracing::warn!("SPP rejected due to impossible velocity: {:.1} m/s", dist / dt);
                    self.current_state = None;
                    return Err(EngineError::SppFailed);
                }
                
                // Pure SPP is an epoch-by-epoch solution. Do not filter.
                state.time = rover_obs.time;
                state.position = pos;
                state.position.epoch = rover_obs.time;"""

if target in content:
    content = content.replace(target, replacement, 1)

with open(filename, "w") as f:
    f.write(content)
