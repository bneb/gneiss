use crate::filter::DdObservation;
use crate::engine::measurement::components::DdComponents;
use gneiss_core::coords::Coordinate;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::time::GpsTime;
use nalgebra::Vector3;

pub struct MeasurementEnvironment<'a> {
    pub ephemerides: &'a [Ephemeris],
    pub sp3_epochs: &'a [gneiss_parsers::sp3::Sp3Epoch],
    pub clk_data: Option<&'a gneiss_parsers::rinex_clk::RinexClock>,
    pub base_coord: &'a Coordinate,
    pub base_time: GpsTime,
    pub lever_arm: Vector3<f64>,
    pub omega_b: Vector3<f64>,
    pub tuning: &'a crate::engine::config::EkfTuningConfig,
    pub gnn_variances: std::collections::HashMap<gneiss_core::sat::SatelliteId, f64>,
    pub klobuchar_params: Option<gneiss_core::atmosphere::KlobucharParams>,
}

pub struct SatState {
    pub rov_pos: Vector3<f64>,
    pub rov_vel: Vector3<f64>,
    pub bas_pos: Vector3<f64>,
    pub bas_vel: Vector3<f64>,
    pub f1: f64,
    pub f2: f64,
}

pub struct EkfUpdates {
    pub z: Vec<f64>,
    pub h: Vec<Vec<f64>>,
    pub r: Vec<f64>,
    pub mt: Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
}

impl Default for EkfUpdates {
    fn default() -> Self {
        Self::new()
    }
}

impl EkfUpdates {
    pub fn new() -> Self {
        Self {
            z: Vec::new(),
            h: Vec::new(),
            r: Vec::new(),
            mt: Vec::new(),
        }
    }
    pub fn push(&mut self, u: SingleUpdate, sat: gneiss_core::sat::SatelliteId) {
        self.z.push(u.z);
        self.h.push(u.h);
        self.r.push(u.r);
        self.mt.push((sat, u.type_code, u.r_ref));
    }
    pub fn extend(&mut self, other: Self) {
        self.z.extend(other.z);
        self.h.extend(other.h);
        self.r.extend(other.r);
        self.mt.extend(other.mt);
    }
}

pub struct DdContext<'a> {
    pub rov_sat: &'a mut DdObservation,
    pub base_sat: &'a mut DdObservation,
    pub rov_ref: &'a mut DdObservation,
    pub ref_base: &'a mut DdObservation,
    pub sat_state: &'a SatState,
    pub ref_state: &'a SatState,
}

pub struct SingleUpdate {
    pub z: f64,
    pub h: Vec<f64>,
    pub r: f64,
    pub type_code: u8,
    pub r_ref: f64,
}

pub struct UpdateGeometry {
    pub comp_dd: f64,
    pub h_r: Vector3<f64>,
    pub h_att: Vector3<f64>,
    pub h_zwd: f64,
    pub state_size: usize,
}

pub struct VarianceWeights {
    pub val: f64,
    pub ref_val: f64,
}

pub struct DdMeasurementContext<'a> {
    pub ctx: &'a DdContext<'a>,
    pub geom: &'a UpdateGeometry,
    pub comps: &'a DdComponents,
    pub env: &'a MeasurementEnvironment<'a>,
}

pub struct DdCarrierPhaseParams<'a> {
    pub is_fixed: bool,
    pub ambiguities: &'a [f64],
    pub sat_idx_l1: Option<usize>,
    pub ref_idx_l1: Option<usize>,
    pub sat_idx_l2: Option<usize>,
    pub ref_idx_l2: Option<usize>,
    /// Ionosphere state index and value for the rover satellite (freq band 3)
    pub iono_idx_sat: Option<usize>,
    /// Ionosphere state index and value for the reference satellite (freq band 3)
    pub iono_idx_ref: Option<usize>,
    pub iono_state_vals: &'a [f64], // all ambiguity values (includes iono states)
    pub cp_base_var: f64,
}
