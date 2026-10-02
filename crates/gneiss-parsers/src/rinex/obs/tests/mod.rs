#![allow(clippy::unwrap_used)]

use super::*;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use gneiss_core::obs::ObsType;
use std::collections::HashMap;
use std::io::BufReader;

mod part1;
mod part2;
