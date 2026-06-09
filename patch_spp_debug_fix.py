import sys

filename = "crates/gneiss-rtk/src/spp.rs"
with open(filename, "r") as f:
    content = f.read()

target = """            // Check delta from seed (which is prev_state if available)
            if let Some(prev) = prev_state {
                let jx = state.position.vector.x - prev.position.vector.x;
                let jy = state.position.vector.y - prev.position.vector.y;
                let jz = state.position.vector.z - prev.position.vector.z;
                let jump = f64::sqrt(jx*jx + jy*jy + jz*jz);
                if jump > 10000.0 {
                    tracing::warn!("SPP HUGE JUMP: {:.1}m. Iterations: 0", jump);
                }
            }"""

replacement = """            // Check delta from prev_state
            {
                let jx = state.position.vector.x - prev_state.position.vector.x;
                let jy = state.position.vector.y - prev_state.position.vector.y;
                let jz = state.position.vector.z - prev_state.position.vector.z;
                let jump = f64::sqrt(jx*jx + jy*jy + jz*jz);
                if jump > 10000.0 {
                    tracing::warn!("SPP HUGE JUMP: {:.1}m.", jump);
                }
            }"""

if target in content:
    content = content.replace(target, replacement, 1)

with open(filename, "w") as f:
    f.write(content)
