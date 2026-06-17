import os

with open("crates/gneiss-rtk/src/engine/ppp_math.rs", "w") as f:
    f.write("""use gneiss_core::sat::{SatelliteId, Constellation};
use gneiss_core::time::GpsTime;
use gneiss_core::obs::{ObsCode, SatObs};
use nalgebra::Vector3;
use crate::engine::ProcessingEngine;
use chrono::TimeZone;

pub const LIGHT_SPEED: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

pub struct OsbCorrections {
    pub p1: Option<f64>, pub p2: Option<f64>, pub cp1: Option<f64>, pub cp2: Option<f64>,
    pub osb_p1: f64, pub osb_p2: f64, pub osb_cp1: f64, pub osb_cp2: f64,
}

pub fn apply_osb_corrections(
    sinex: &gneiss_parsers::sinex_bia::SinexBias,
    sat: SatelliteId,
    time: GpsTime,
    f1: f64, f2: f64,
    p1_obs: Option<(f64, ObsCode)>,
    p2_obs: Option<(f64, ObsCode)>,
    cp1_obs: Option<(f64, ObsCode)>,
    cp2_obs: Option<(f64, ObsCode)>
) -> OsbCorrections {
    let mut p1 = p1_obs.map(|o| o.0);
    let mut p2 = p2_obs.map(|o| o.0);
    let mut cp1 = cp1_obs.map(|o| o.0);
    let mut cp2 = cp2_obs.map(|o| o.0);
    let mut osb_p1 = 0.0; let mut osb_p2 = 0.0; let mut osb_cp1 = 0.0; let mut osb_cp2 = 0.0;
    let mut is_cp2_fb = false;

    let convert = |b: f64| b * 1e-9 * LIGHT_SPEED;
    if let Some((_, c)) = p1_obs { if let Some(b) = sinex.get_bias(sat, c, time) { osb_p1 = convert(b); } }
    if let Some((_, c)) = p2_obs { if let Some(b) = sinex.get_bias(sat, c, time) { osb_p2 = convert(b); } }
    if let Some((_, c)) = cp1_obs { if let Some(b) = sinex.get_bias(sat, c, time) { osb_cp1 = convert(b); } }
    if let Some((_, c)) = cp2_obs { 
        if let Some(b) = sinex.get_bias(sat, c, time) { 
            osb_cp2 = convert(b); 
            is_cp2_fb = sinex.get_exact_bias(sat, c, time).is_none();
        } 
    }

    if let Some(v) = p1.as_mut() { *v -= osb_p1; }
    if let Some(v) = p2.as_mut() { *v -= osb_p2; }
    if let Some(v) = cp1.as_mut() { *v -= osb_cp1 / (LIGHT_SPEED / f1); }
    if let Some(v) = cp2.as_mut() { 
        *v -= osb_cp2 / (LIGHT_SPEED / f2); 
        if is_cp2_fb && sat.constellation == Constellation::Gps {
            if let Some((_, c)) = cp2_obs {
                let s = c.to_string();
                if s == "L2L" || s == "L2S" || s == "L2X" { *v -= 0.25; }
            }
        }
    }
    OsbCorrections { p1, p2, cp1, cp2, osb_p1, osb_p2, osb_cp1, osb_cp2 }
}

pub fn apply_earth_rotation(
    raw_pos: Vector3<f64>,
    raw_vel: Vector3<f64>,
    rcv_pos: Vector3<f64>,
) -> (Vector3<f64>, Vector3<f64>) {
    let mut sat_pos = raw_pos;
    let mut sat_vel = raw_vel;
    for _ in 0..2 {
        let rng = (sat_pos - rcv_pos).norm();
        let tau = rng / LIGHT_SPEED;
        let th = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * tau;
        let c = libm::cos(th);
        let s = libm::sin(th);
        sat_pos = Vector3::new(raw_pos.x * c + raw_pos.y * s, -raw_pos.x * s + raw_pos.y * c, raw_pos.z);
        sat_vel = Vector3::new(raw_vel.x * c + raw_vel.y * s, -raw_vel.x * s + raw_vel.y * c, raw_vel.z);
    }
    (sat_pos, sat_vel)
}

pub fn detect_cycle_slip(sat_obs: &SatObs, prev: u32) -> (bool, u32) {
    let mut slip = false;
    let mut lk = prev.saturating_add(1);
    for obs in &sat_obs.observations {
        if obs.code.obs_type == gneiss_core::obs::ObsType::CarrierPhase {
            if let Some(lli) = obs.lli {
                if (lli & 1) != 0 || (lli & 2) != 0 { slip = true; lk = 0; }
            } else if let Some(l) = obs.lock_time {
                if l == 0 || l < prev { slip = true; lk = l; } else { lk = lk.min(l); }
            }
        }
    }
    if lk == 0 && prev > 0 { slip = true; }
    (slip, lk)
}

pub fn compute_tropo_dry(rcv_pos_llh: Vector3<f64>, el: f64, t: GpsTime) -> (f64, f64) {
    let tp = gneiss_core::atmosphere::TropoParams::default();
    let p_z = tp.press_hpa * libm::pow(1.0 - 0.0000226 * rcv_pos_llh.z, 5.225);
    let z_dry = 0.0022768 * p_z / (1.0 - 0.00266 * libm::cos(2.0 * rcv_pos_llh.x) - 0.00028 * rcv_pos_llh.z / 1000.0);
    let (mh, mw) = gneiss_core::atmosphere::AtmosphereModel::nmf_mapping_functions(rcv_pos_llh, el, t);
    (z_dry * mh, mw)
}

pub fn compute_iono_free(f1: f64, f2: f64, v1: f64, v2: f64) -> f64 {
    let g = (f1 * f1) / (f2 * f2);
    (g * v1 - v2) / (g - 1.0)
}
""")

mod_rs = "crates/gneiss-rtk/src/engine/mod.rs"
with open(mod_rs, "r") as f:
    mod_content = f.read()
if "pub mod ppp_math;" not in mod_content:
    mod_content = mod_content.replace("pub mod ppp;", "pub mod ppp;\npub mod ppp_math;")
    with open(mod_rs, "w") as f:
        f.write(mod_content)
