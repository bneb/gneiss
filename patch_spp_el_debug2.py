import sys

filename = "crates/gneiss-rtk/src/spp.rs"
with open(filename, "r") as f:
    content = f.read()

target = """    for (i, m) in measurements.iter().enumerate() {
        let sat_ecef = Vector3::new(m.sat_coord.vector.x, m.sat_coord.vector.y, m.sat_coord.vector.z);
        let (_, el) = az_el(rec_llh_seed, rec_ecef_seed, sat_ecef);
        tracing::debug!("Sat PRN {} el: {:.2} deg (mask {:.2} deg)", m.sat_coord.vector.x, el * 180.0 / core::f64::consts::PI, el_mask * 180.0 / core::f64::consts::PI);
        if el >= el_mask || rec_ecef_seed.x == 0.0 {
            valid_indices.push(i);
        }
    }
    tracing::warn!("valid_indices.len() = {} / {}", valid_indices.len(), measurements.len());"""

replacement = """    tracing::warn!("rec_ecef_seed: {:?}", rec_ecef_seed);
    tracing::warn!("rec_llh_seed: {:?}", rec_llh_seed);
    for (i, m) in measurements.iter().enumerate() {
        let sat_ecef = Vector3::new(m.sat_coord.vector.x, m.sat_coord.vector.y, m.sat_coord.vector.z);
        let (_, el) = az_el(rec_llh_seed, rec_ecef_seed, sat_ecef);
        if el >= el_mask || rec_ecef_seed.x == 0.0 {
            valid_indices.push(i);
        }
    }"""

content = content.replace(target, replacement)

with open(filename, "w") as f:
    f.write(content)
