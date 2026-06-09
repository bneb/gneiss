import sys

filename = "crates/gneiss-rtk/src/spp.rs"
with open(filename, "r") as f:
    content = f.read()

target = """                        if f64::sqrt(c_dx * c_dx + c_dy * c_dy + c_dz * c_dz + c_dcdt * c_dcdt) < config.convergence_threshold {
                            tracing::debug!("compute_spp clean converged: {:?}", clean_state);
                            return Ok(clean_state);
                        }"""

replacement = """                        if f64::sqrt(c_dx * c_dx + c_dy * c_dy + c_dz * c_dz + c_dcdt * c_dcdt) < config.convergence_threshold {
                            let clean_llh = gneiss_core::coords::ecef_to_llh(clean_state.position.vector);
                            if clean_llh.z < -20000.0 || clean_llh.z > 100000.0 {
                                tracing::warn!("SPP clean converged to invalid altitude: {:.1}m", clean_llh.z);
                                return Err(SppError::ConvergenceFailed);
                            }
                            tracing::debug!("compute_spp clean converged: {:?}", clean_state);
                            return Ok(clean_state);
                        }"""

if target in content:
    content = content.replace(target, replacement, 1)

with open(filename, "w") as f:
    f.write(content)
