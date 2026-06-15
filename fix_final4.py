with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    text = f.read()

# The block is from line 710 to 729
# Let's extract it
block = """
/// DD measurements sharing a reference satellite have correlated noise equal to
/// the reference satellite's measurement variance.
const DD_CROSS_CORRELATION_SCALE: f64 = 1.0;

pub fn build_dense_covariance_matrix(
    r_diagonals: &[f64],
    meas_types: &[(gneiss_core::sat::SatelliteId, u8, f64)],
) -> DMatrix<f64> {
    let mut r_mat = DMatrix::from_diagonal(&DVector::from_row_slice(r_diagonals));
    for i in 0..meas_types.len() {
        for j in (i + 1)..meas_types.len() {
            if meas_types[i].1 == meas_types[j].1 && meas_types[i].0.constellation == meas_types[j].0.constellation {
                let cov = meas_types[i].2.min(meas_types[j].2) * DD_CROSS_CORRELATION_SCALE;
                r_mat[(i, j)] = cov;
                r_mat[(j, i)] = cov;
            }
        }
    }
    r_mat
}
"""

# Remove it from inside the function
text = text.replace(block, "")

# Append it to the end of the file
text += "\n" + block + "\n"

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(text)

