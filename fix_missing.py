import re

# Add missing functions to updater_math.rs
with open("crates/gneiss-rtk/src/engine/updater_math.rs", "r") as f:
    text = f.read()

missing_funcs = """
pub fn filter_pre_fit_residuals(
    z: &DVector<f64>,
    h: &DMatrix<f64>,
    r: &DMatrix<f64>,
    state_cov: &DMatrix<f64>,
    max_innovation: f64,
    meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>,
) -> Vec<usize> {
    let mut valid_indices = Vec::with_capacity(z.len());
    let hp = h * state_cov;
    
    for i in 0..z.len() {
        let s_ii = hp.row(i).dot(&h.row(i).transpose()) + r[(i, i)];
        let meas_type = meas_types.map_or(0, |m| m[i].1);
        let threshold = get_pre_fit_threshold(meas_type, max_innovation);
        
        if check_pre_fit_residual(z[i], s_ii, r[(i, i)], meas_type, threshold) {
            valid_indices.push(i);
        }
    }
    valid_indices
}

pub fn compute_loose_coupling_innovations(
    r_b_e: &nalgebra::Matrix3<f64>,
    state_pos: &Vector3<f64>,
    state_vel: &Vector3<f64>,
    gnss_pos: &Vector3<f64>,
    gnss_vel: &Vector3<f64>,
    lever_arm: &Vector3<f64>,
    omega_b: &Vector3<f64>,
) -> DVector<f64> {
    let l_e = r_b_e * lever_arm;
    let pos_apc = state_pos + l_e;
    let v_apc = state_vel + r_b_e * omega_b.cross(lever_arm);

    let mut z = DVector::zeros(6);
    z.rows_mut(0, 3).copy_from(&(gnss_pos - pos_apc));
    z.rows_mut(3, 3).copy_from(&(gnss_vel - v_apc));
    z
}
"""

if "pub fn filter_pre_fit_residuals" not in text:
    # insert before the first test
    if "#[test]" in text:
        text = text.replace("#[test]", missing_funcs + "\n#[test]", 1)
    else:
        text += missing_funcs

with open("crates/gneiss-rtk/src/engine/updater_math.rs", "w") as f:
    f.write(text)

# Fix tests_updater.rs
with open("crates/gneiss-rtk/src/engine/tests_updater.rs", "r") as f:
    tests_text = f.read()

tests_text = tests_text.replace("crate::engine::updater::evaluate_post_fit_outliers", "crate::engine::updater_math::evaluate_post_fit_outliers")
tests_text = tests_text.replace("crate::engine::updater::filter_pre_fit_residuals", "crate::engine::updater_math::filter_pre_fit_residuals")
tests_text = tests_text.replace("crate::engine::updater::apply_joseph_covariance_update", "crate::engine::updater_math::apply_joseph_covariance_update")

with open("crates/gneiss-rtk/src/engine/tests_updater.rs", "w") as f:
    f.write(tests_text)

