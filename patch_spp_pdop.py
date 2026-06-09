import sys

filename = "crates/gneiss-rtk/src/spp.rs"
with open(filename, "r") as f:
    content = f.read()

target = """    let pos_variance = h_t_w_h_inv[(0, 0)] + h_t_w_h_inv[(1, 1)] + h_t_w_h_inv[(2, 2)];
    if pos_variance > config.geometry_variance_threshold {
        return Err(SppError::PoorGeometry);
    }"""

replacement = """    let pos_variance = h_t_w_h_inv[(0, 0)] + h_t_w_h_inv[(1, 1)] + h_t_w_h_inv[(2, 2)];
    // tracing::info!("pos_variance: {:.1}", pos_variance);
    if pos_variance > config.geometry_variance_threshold {
        return Err(SppError::PoorGeometry);
    }"""

if target in content:
    content = content.replace(target, replacement, 1)

with open(filename, "w") as f:
    f.write(content)
