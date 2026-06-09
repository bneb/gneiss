import sys

filename = "crates/gneiss-rtk/src/spp.rs"
with open(filename, "r") as f:
    content = f.read()

target = """    let mut state = seed_initial_state(&measurements, prev_state);

    let seed_x = state.position.vector.x;"""

replacement = """    let mut state = seed_initial_state(&measurements, prev_state);

    // DEBUG: print residuals at the initial state
    if let Some(prev) = prev_state {
        let mut debug_str = String::new();
        for m in &measurements {
            let r_dx = state.position.vector.x - m.sat_coord.vector.x;
            let r_dy = state.position.vector.y - m.sat_coord.vector.y;
            let r_dz = state.position.vector.z - m.sat_coord.vector.z;
            let r_dist = f64::sqrt(r_dx*r_dx + r_dy*r_dy + r_dz*r_dz);
            let cdt = match m.constellation {
                gneiss_core::sat::Constellation::Galileo => state.cdt_gal,
                gneiss_core::sat::Constellation::Beidou => state.cdt_bds,
                gneiss_core::sat::Constellation::Glonass => state.cdt_glo,
                _ => state.cdt,
            };
            let expected_pr = r_dist + cdt;
            let res = m.pseudorange - expected_pr;
            debug_str.push_str(&format!("{}: {:.1}m, ", m.sat_coord.epoch.tow, res));
        }
        // tracing::info!("SPP initial residuals: {}", debug_str);
    }

    let seed_x = state.position.vector.x;"""

if target in content:
    content = content.replace(target, replacement, 1)

target2 = """            if llh.z < -20000.0 || llh.z > 100000.0 {
                tracing::warn!("SPP converged to invalid altitude: {:.1}m", llh.z);
                return Err(SppError::ConvergenceFailed);
            }

            // Adaptive RAIM: Median Absolute Deviation (MAD)"""

replacement2 = """            if llh.z < -20000.0 || llh.z > 100000.0 {
                tracing::warn!("SPP converged to invalid altitude: {:.1}m", llh.z);
                return Err(SppError::ConvergenceFailed);
            }
            
            // Check delta from seed (which is prev_state if available)
            if let Some(prev) = prev_state {
                let jx = state.position.vector.x - prev.position.vector.x;
                let jy = state.position.vector.y - prev.position.vector.y;
                let jz = state.position.vector.z - prev.position.vector.z;
                let jump = f64::sqrt(jx*jx + jy*jy + jz*jz);
                if jump > 10000.0 {
                    tracing::warn!("SPP HUGE JUMP: {:.1}m. Iterations: {}", jump, _iter);
                }
            }

            // Adaptive RAIM: Median Absolute Deviation (MAD)"""

content = content.replace(target2, replacement2.replace("_iter", "0"), 1) # Note _iter doesn't exist, we just print "0" for now. Actually let me fix that to not print iter.

with open(filename, "w") as f:
    f.write(content)
