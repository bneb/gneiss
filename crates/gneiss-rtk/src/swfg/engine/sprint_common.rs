//! Shared constants for the SWFG engine sprint-hunt tests.
#![allow(dead_code)]

use nalgebra::Vector3;

/// CORS base station used by the RTK factor fixtures.
pub(super) const BASE: Vector3<f64> = Vector3::new(-3_961_904.434_1, 3_348_994.266, 3_698_211.706_7);
/// Rover position the synthetic double differences are made consistent with.
pub(super) const TRUE_POS: Vector3<f64> = Vector3::new(-3_960_000.0, 3_349_000.0, 3_698_000.0);
/// GPS L1 carrier wavelength: c / 1575.42 MHz = 0.19029367 m.
pub(super) const LAMBDA_L1: f64 = 299_792_458.0 / 1_575.42e6;