import re

with open("crates/gneiss-rtk/src/engine/processor.rs", "r") as f:
    proc_text = f.read()

proc_text = proc_text.replace("self.config.mode.is_tightly_coupled(),", "")

with open("crates/gneiss-rtk/src/engine/processor.rs", "w") as f:
    f.write(proc_text)

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    meas_text = f.read()

# Replace observation_variance in measurement.rs
meas_text = meas_text.replace(
    "gneiss_core::variance::observation_variance(ctx.rov_ref.snr, el_rov_ref, BASE_SNR_ELEVATION_THRESH_DEG)",
    "gneiss_core::variance::observation_variance(ctx.rov_ref.snr, el_rov_ref, env.tuning.snr_a, env.tuning.snr_b)"
)
meas_text = meas_text.replace(
    "gneiss_core::variance::observation_variance(ctx.rov_sat.snr, el_rov_sat, BASE_SNR_ELEVATION_THRESH_DEG)",
    "gneiss_core::variance::observation_variance(ctx.rov_sat.snr, el_rov_sat, env.tuning.snr_a, env.tuning.snr_b)"
)

# Remove the test boilerplate
test_boilerplate = """
#[derive(Debug, Clone)]
pub struct SingleMeasurementUpdate {
    pub innovation: f64,
    pub h_row: Vec<f64>,
    pub variance: f64,
    pub meas_type: u8,
    pub ref_variance: f64,
}

#[derive(Debug, Clone)]
pub struct MeasurementModel {
    pub z_vec: DVector<f64>,
    pub h_mat: DMatrix<f64>,
    pub r_mat: DMatrix<f64>,
    pub meas_types: Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
}

use crate::engine::measurement_math::*;
"""

# The script that did the replacement might have inserted slightly different indentation or newlines, so let's use regex
pattern = re.compile(r"#\[derive\(Debug, Clone\)\]\s*pub struct SingleMeasurementUpdate \{[^\}]*\}\s*#\[derive\(Debug, Clone\)\]\s*pub struct MeasurementModel \{[^\}]*\}\s*use crate::engine::measurement_math::\*;")
meas_text = pattern.sub("", meas_text)

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(meas_text)
