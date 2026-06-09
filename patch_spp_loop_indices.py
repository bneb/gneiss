import sys

filename = "crates/gneiss-rtk/src/spp.rs"
with open(filename, "r") as f:
    content = f.read()

target1 = """    for (i, m) in measurements.iter().enumerate() {"""
replacement1 = """    for (row, &i) in valid_indices.iter().enumerate() {
        let m = &measurements[i];"""

target2 = """        #[cfg(test)]
        let el_mask = -core::f64::consts::PI; // Disable in unit tests
        #[cfg(not(test))]
        let el_mask = 0.1745; // 10 degrees

        // Exclude satellites below elevation mask
        if el < el_mask && current_state.position.vector.x != 0.0 {
            // Give it an incredibly low weight so it doesn't affect the solution
            w_matrix[(i, i)] = 1e-10;
            continue;
        }"""
replacement2 = """        // Elevation mask is now handled before the loop!"""

content = content.replace(target1, replacement1, 1)
content = content.replace(target2, replacement2, 1)

# Replace all matrix indexing
content = content.replace("h_matrix[(i, 0)] = dx / r_safe;", "h_matrix[(row, 0)] = dx / r_safe;")
content = content.replace("h_matrix[(i, 1)] = dy / r_safe;", "h_matrix[(row, 1)] = dy / r_safe;")
content = content.replace("h_matrix[(i, 2)] = dz / r_safe;", "h_matrix[(row, 2)] = dz / r_safe;")
content = content.replace("h_matrix[(i, c)] = 1.0;", "h_matrix[(row, c)] = 1.0;")
content = content.replace("dz_vector[i] = residual;", "dz_vector[row] = residual;")
content = content.replace("w_matrix[(i, i)] = 1.0 / total_var;", "w_matrix[(row, row)] = 1.0 / total_var;")

with open(filename, "w") as f:
    f.write(content)
