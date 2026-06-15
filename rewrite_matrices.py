import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

# Replace build_measurement_model return type
old_sig_build = """pub fn build_measurement_model(
    state: &mut RtkState, matched_obs: &[(DdObservation, DdObservation)],
    env: &MeasurementEnvironment, chi_square_pr_threshold: f64, chi_square_cp_threshold: f64,
) -> Option<(DVector<f64>, DMatrix<f64>, DMatrix<f64>, Vec<(gneiss_core::sat::SatelliteId, u8, f64)>)> {"""

new_sig_build = """pub struct EkfMeasurementMatrices {
    pub z: DVector<f64>,
    pub h: DMatrix<f64>,
    pub r: DMatrix<f64>,
    pub mt: Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
}

pub fn build_measurement_model(
    state: &mut RtkState, matched_obs: &[(DdObservation, DdObservation)],
    env: &MeasurementEnvironment, chi_square_pr_threshold: f64, chi_square_cp_threshold: f64,
) -> Option<EkfMeasurementMatrices> {"""

content = content.replace(old_sig_build, new_sig_build)

# Replace build_final_measurement_matrices return type
old_sig_final = """fn build_final_measurement_matrices(
    state_size: usize,
    safe_indices: Vec<usize>,
    z_all: &[f64],
    h_all: &[Vec<f64>],
    r_all: &[f64],
    type_all: &[(gneiss_core::sat::SatelliteId, u8, f64)],
) -> Option<(DVector<f64>, DMatrix<f64>, DMatrix<f64>, Vec<(gneiss_core::sat::SatelliteId, u8, f64)>)> {"""

new_sig_final = """fn build_final_measurement_matrices(
    state_size: usize, safe_indices: Vec<usize>, z_all: &[f64], h_all: &[Vec<f64>], r_all: &[f64], type_all: &[(gneiss_core::sat::SatelliteId, u8, f64)],
) -> Option<EkfMeasurementMatrices> {"""

content = content.replace(old_sig_final, new_sig_final)

# Replace Some((z_vec, h_mat, r_mat, t_vec)) with Some(EkfMeasurementMatrices { ... })
content = content.replace("Some((z_vec, h_mat, r_mat, t_vec))", "Some(EkfMeasurementMatrices { z: z_vec, h: h_mat, r: r_mat, mt: t_vec })")

# Replace compute_phase_windup
old_windup = """#[allow(clippy::too_many_arguments)]
pub fn compute_phase_windup(
    time: GpsTime,
    pos_apc: Vector3<f64>,
    base_coord_vec: Vector3<f64>,
    rov_sat_pos: Vector3<f64>,
    rov_ref_pos: Vector3<f64>,
    bas_sat_pos: Vector3<f64>,
    bas_ref_pos: Vector3<f64>,
    prev_w_sat: f64,
    prev_w_ref: f64,
    prev_w_bas_sat: f64,
    prev_w_bas_ref: f64,
) -> (f64, f64, f64, f64) {
    let sun_pos = gneiss_core::sun::sun_position_ecef(time);
    let w_sat = gneiss_core::windup::phase_windup(rov_sat_pos, sun_pos, pos_apc, prev_w_sat);
    let w_ref = gneiss_core::windup::phase_windup(rov_ref_pos, sun_pos, pos_apc, prev_w_ref);
    let w_bas_sat = gneiss_core::windup::phase_windup(bas_sat_pos, sun_pos, base_coord_vec, prev_w_bas_sat);
    let w_bas_ref = gneiss_core::windup::phase_windup(bas_ref_pos, sun_pos, base_coord_vec, prev_w_bas_ref);
    (w_sat, w_ref, w_bas_sat, w_bas_ref)
}"""
content = content.replace(old_windup, "")

# Replace usage of compute_phase_windup with crate::engine::measurement_math::compute_phase_windup
content = content.replace("compute_phase_windup(", "crate::engine::measurement_math::compute_phase_windup(")
content = content.replace("crate::engine::measurement_math::crate::engine::measurement_math::compute_phase_windup(", "crate::engine::measurement_math::compute_phase_windup(")

# Delete the duplicate test for compute_phase_windup in measurement.rs
# It starts around line 892 and ends around 927. We can just use regex.
import re
content = re.sub(r"#[test]\s*fn test_compute_phase_windup\(\).*?(?=\n    #\[test\])", "", content, flags=re.DOTALL)
content = content.replace("use crate::engine::measurement::{compute_phase_windup, DdContext, SatState};", "use crate::engine::measurement::{DdContext, SatState};")


with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)

