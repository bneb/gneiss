import sys

filename = "crates/gneiss-rtk/src/spp.rs"
with open(filename, "r") as f:
    content = f.read()

target = """        if delta_norm < config.convergence_threshold {
            // Adaptive RAIM: Median Absolute Deviation (MAD)"""

replacement = """        if delta_norm < config.convergence_threshold {
            let llh = gneiss_core::coords::ecef_to_llh(state.position.vector);
            if llh.z < -20000.0 || llh.z > 100000.0 {
                tracing::warn!("SPP converged to invalid altitude: {:.1}m", llh.z);
                return Err(SppError::ConvergenceFailed);
            }

            // Adaptive RAIM: Median Absolute Deviation (MAD)"""

if target in content:
    content = content.replace(target, replacement, 1)

with open(filename, "w") as f:
    f.write(content)
