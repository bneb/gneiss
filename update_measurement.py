import re

with open("crates/gneiss-rtk/src/engine/measurement.rs", "r") as f:
    content = f.read()

# Replace tuples with structs
content = content.replace("Option<(DVector<f64>, DMatrix<f64>, DMatrix<f64>, Vec<(gneiss_core::sat::SatelliteId, u8, f64)>)>", "Option<MeasurementModel>")
content = content.replace("Some((z_vec, h_mat, r_mat, t_vec))", "Some(MeasurementModel { z_vec, h_mat, r_mat, meas_types: t_vec })")

# Replace (f64, Vec<f64>, f64, u8, f64) with SingleMeasurementUpdate
content = content.replace("Vec<(f64, Vec<f64>, f64, u8, f64)>", "Vec<SingleMeasurementUpdate>")
content = content.replace("Option<(f64, Vec<f64>, f64, u8, f64)>", "Option<SingleMeasurementUpdate>")

# Insert structs
structs = """
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

content = content.replace("use crate::filter::{RtkState, DdObservation};", "use crate::filter::{RtkState, DdObservation};\n" + structs)

with open("crates/gneiss-rtk/src/engine/measurement.rs", "w") as f:
    f.write(content)

