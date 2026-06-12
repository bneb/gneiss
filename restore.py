import sys

# I will write it line by line to avoid issues
lines = [
"use gneiss_core::obs::{EpochObs, SatObs};",
"use gneiss_core::coords::Coordinate;",
"use crate::filter::RtkState;",
"use crate::engine::{EngineError, ProcessingEngine};",
"use nalgebra::{Vector3, DMatrix, DVector};",
"use chrono::TimeZone;",
"use crate::engine::processed_sat::ProcessedSat;",
"use crate::factor_graph::{FactorGraphOptimizer, gnss_factors::{ErrorStatePseudorangeFactor, ErrorStateCarrierPhaseFactor}};",
"",
"use gneiss_core::constants::{SPEED_OF_LIGHT_M_S, EARTH_ROTATION_RATE_RAD_S};",
"",
"fn snr_scale(snr: i32) -> f64 {",
"    (10.0f64).powf((45.0 - snr as f64) / 10.0)",
"}",
]
