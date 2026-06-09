import sys

filename = "crates/gneiss-rtk/src/spp.rs"
with open(filename, "r") as f:
    content = f.read()

target = """    let dx_vec = h_t_w_h_inv * h_t_w * dz_vector;

    let next_cdt = current_state.cdt + gps_col.map(|c| dx_vec[c]).unwrap_or(0.0);"""

replacement = """    let dx_vec = h_t_w_h_inv * h_t_w * &dz_vector;
    tracing::warn!("dx_vec: {:?}", dx_vec.transpose());
    tracing::warn!("dz_vector (first 5): {:?}", &dz_vector.as_slice()[0..std::cmp::min(5, dz_vector.len())]);
    let next_cdt = current_state.cdt + gps_col.map(|c| dx_vec[c]).unwrap_or(0.0);"""

content = content.replace(target, replacement)

with open(filename, "w") as f:
    f.write(content)
