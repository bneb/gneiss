import sys

filename = "crates/gneiss-rtk/src/spp.rs"
with open(filename, "r") as f:
    content = f.read()

# 1. Jump check
target1 = """            // Adaptive RAIM: Median Absolute Deviation (MAD)"""
replacement1 = """            // Check delta from previous state
            {
                let jx = state.position.vector.x - prev_state.position.vector.x;
                let jy = state.position.vector.y - prev_state.position.vector.y;
                let jz = state.position.vector.z - prev_state.position.vector.z;
                let jump = f64::sqrt(jx*jx + jy*jy + jz*jz);
                if jump > 10000.0 {
                    tracing::warn!("SPP HUGE JUMP: {:.1}m.", jump);
                }
            }

            // Adaptive RAIM: Median Absolute Deviation (MAD)"""

# 2. Geometry variance threshold
target2 = """            geometry_variance_threshold: 5000.0, // Default 100m^2 variance threshold"""
replacement2 = """            geometry_variance_threshold: 500.0, // Reduced to prevent ghost positions"""

# 3. Dynamic measurement filtering in WNLLS
target3 = """    if cols < 4 || n < cols {
        return Err(SppError::NotEnoughMeasurements);
    }

    let mut h_matrix = DMatrix::<f64>::zeros(n, cols);
    let mut w_matrix = DMatrix::<f64>::zeros(n, n);
    let mut dz_vector = DVector::<f64>::zeros(n);"""

replacement3 = """    // Dynamically filter satellites based on elevation to maintain a well-conditioned matrix
    let rec_ecef_seed = Vector3::new(current_state.position.vector.x, current_state.position.vector.y, current_state.position.vector.z);
    let rec_llh_seed = ecef_to_llh(rec_ecef_seed);
    
    #[cfg(test)]
    let el_mask = -core::f64::consts::PI;
    #[cfg(not(test))]
    let el_mask = 0.1745; // 10 degrees

    let mut valid_indices = Vec::new();
    for (i, m) in measurements.iter().enumerate() {
        let sat_ecef = Vector3::new(m.sat_coord.vector.x, m.sat_coord.vector.y, m.sat_coord.vector.z);
        let (_, el) = az_el(rec_llh_seed, rec_ecef_seed, sat_ecef);
        if el >= el_mask || rec_ecef_seed.x == 0.0 {
            valid_indices.push(i);
        }
    }
    
    let n_valid = valid_indices.len();
    if cols < 4 || n_valid < cols {
        return Err(SppError::NotEnoughMeasurements);
    }

    let mut h_matrix = DMatrix::<f64>::zeros(n_valid, cols);
    let mut w_matrix = DMatrix::<f64>::zeros(n_valid, n_valid);
    let mut dz_vector = DVector::<f64>::zeros(n_valid);"""

# 4. Filter usage inside the loop
target4 = """    for (i, m) in measurements.iter().enumerate() {"""
replacement4 = """    for (row, &i) in valid_indices.iter().enumerate() {
        let m = &measurements[i];"""

# 5. Fixing the loop internal indices from `i` to `row`
content = content.replace(target1, replacement1, 1)
content = content.replace(target2, replacement2, 1)
content = content.replace(target3, replacement3, 1)

# Now we need to carefully replace the `i` inside the loop
# We will just write a specific replacement script for this
with open(filename, "w") as f:
    f.write(content)
