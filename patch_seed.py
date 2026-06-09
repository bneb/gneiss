import sys

filename = "crates/gneiss-rtk/src/spp.rs"
with open(filename, "r") as f:
    content = f.read()

target = """    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_z = 0.0;
    for m in measurements {
        sum_x += m.sat_coord.vector.x;
        sum_y += m.sat_coord.vector.y;
        sum_z += m.sat_coord.vector.z;
    }
    let n = measurements.len() as f64;
    let mut seed_x = sum_x / n;
    let mut seed_y = sum_y / n;
    let mut seed_z = sum_z / n;

    // Project onto Earth's surface
    let seed_pos = gneiss_core::coords::llh_to_ecef(Vector3::new(
        gneiss_core::coords::ecef_to_llh(Vector3::new(seed_x, seed_y, seed_z)).x,
        gneiss_core::coords::ecef_to_llh(Vector3::new(seed_x, seed_y, seed_z)).y,
        0.0,
    ));
    seed_x = seed_pos.x;
    seed_y = seed_pos.y;
    seed_z = seed_pos.z;"""

replacement = """    // Just use [0,0,0] as seed
    let seed_x = 0.0;
    let seed_y = 0.0;
    let seed_z = 0.0;"""

if target in content:
    content = content.replace(target, replacement, 1)

with open(filename, "w") as f:
    f.write(content)
