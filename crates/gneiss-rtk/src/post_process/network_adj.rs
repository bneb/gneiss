//! Multi-Baseline CORS Network Adjustment Engine.
//!
//! Formulates inter-station baseline graphs and double-difference integer
//! ambiguity network adjustments across regional CORS networks to isolate
//! localized atmospheric delay states (tropospheric ZWD and slant ionosphere).

use std::collections::HashMap;
use nalgebra::Vector3;

use gneiss_core::coords::{az_el, ecef_to_llh};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::sat::SatelliteId;
use gneiss_core::signal::satellite_frequencies;

use crate::spatial::delaunay::EngineError;

/// Regional CORS station metadata and epoch observations.
#[derive(Debug, Clone)]
pub struct CorsStation {
    pub id: String,
    pub pos_ecef: Vector3<f64>,
    pub epochs: Vec<EpochObs>,
}

/// Baseline connecting two CORS stations in the network graph.
#[derive(Debug, Clone)]
pub struct NetworkBaseline {
    pub base_a: String,
    pub base_b: String,
    pub length_m: f64,
}

/// Atmospheric delay state extracted at a CORS station.
#[derive(Debug, Clone, Default)]
pub struct StationAtmosphere {
    pub station_id: String,
    pub zwd_m: f64,
    pub iono_slant_m: HashMap<SatelliteId, f64>,
}

/// Solved network adjustment epoch containing per-station atmospheric states.
#[derive(Debug, Clone)]
pub struct NetworkAdjustmentResult {
    pub master_id: String,
    pub master_pos: Vector3<f64>,
    pub station_atmospheres: HashMap<String, StationAtmosphere>,
    pub fixed_ambiguities: usize,
}

/// Network adjustment solver over a set of regional CORS stations.
pub struct NetworkAdjuster {
    pub master_id: String,
    pub master_pos: Vector3<f64>,
    pub baselines: Vec<NetworkBaseline>,
}

impl NetworkAdjuster {
    /// Construct network adjuster given CORS stations and a nominated master ID.
    pub fn new(stations: &[CorsStation], master_id: &str) -> Result<Self, EngineError> {
        let master = stations.iter().find(|s| s.id == master_id)
            .ok_or_else(|| EngineError::DegenerateMesh(format!("Master station {master_id} not found")))?;
        let mut baselines = Vec::new();
        for s in stations {
            if s.id != master_id {
                let len = (s.pos_ecef - master.pos_ecef).norm();
                baselines.push(NetworkBaseline {
                    base_a: master_id.to_string(),
                    base_b: s.id.clone(),
                    length_m: len,
                });
            }
        }
        Ok(Self {
            master_id: master_id.to_string(),
            master_pos: master.pos_ecef,
            baselines,
        })
    }

    /// Run double-difference network adjustment for a given epoch index.
    pub fn adjust_epoch(
        &self,
        stations: &[CorsStation],
        epoch_idx: usize,
        ephemerides: &[Ephemeris],
    ) -> Result<NetworkAdjustmentResult, EngineError> {
        let master_st = stations.iter().find(|s| s.id == self.master_id)
            .ok_or_else(|| EngineError::DegenerateMesh("Master station not found".into()))?;
        let master_epoch = master_st.epochs.get(epoch_idx)
            .ok_or_else(|| EngineError::InterpolationFailed("Epoch index out of bounds".into()))?;

        let mut atmospheres = HashMap::new();
        let mut total_fixed = 0;
        let mut master_atmo = StationAtmosphere {
            station_id: self.master_id.clone(),
            zwd_m: 0.15, // Nominal a priori ZWD
            iono_slant_m: HashMap::new(),
        };

        for s in stations {
            if s.id == self.master_id { continue; }
            if let Some(sec_epoch) = s.epochs.get(epoch_idx) {
                let (st_atmo, fixed) = adjust_baseline(
                    master_st, master_epoch, s, sec_epoch, ephemerides,
                );
                total_fixed += fixed;
                for (&sat, &val) in &st_atmo.iono_slant_m {
                    master_atmo.iono_slant_m.entry(sat).or_insert(0.0);
                    let _ = val;
                }
                atmospheres.insert(s.id.clone(), st_atmo);
            }
        }
        atmospheres.insert(self.master_id.clone(), master_atmo);

        Ok(NetworkAdjustmentResult {
            master_id: self.master_id.clone(),
            master_pos: self.master_pos,
            station_atmospheres: atmospheres,
            fixed_ambiguities: total_fixed,
        })
    }
}

/// Double-difference carrier phase observations for a satellite pair.
struct DdPairObs {
    sat: SatelliteId,
    dd_rho: f64,
    dd_phi1: f64,
    dd_phi2: f64,
    dd_p1: f64,
    dd_p2: f64,
    el_sec: f64,
}

/// Adjust a single inter-CORS baseline to extract relative atmosphere.
fn adjust_baseline(
    master: &CorsStation,
    m_epoch: &EpochObs,
    sec: &CorsStation,
    s_epoch: &EpochObs,
    ephemerides: &[Ephemeris],
) -> (StationAtmosphere, usize) {
    let pairs = form_dd_pairs(master, m_epoch, sec, s_epoch, ephemerides);
    let mut iono_map = HashMap::new();
    let mut fixed_count = 0;
    let mut zwd_numer = 0.0;
    let mut zwd_denom = 1e-9;

    for p in &pairs {
        if let Some((n1, n2)) = resolve_dd_ambiguities(p) {
            fixed_count += 1;
            let (iono, tropo) = extract_dd_atmosphere(p, n1, n2);
            iono_map.insert(p.sat, iono);
            let mw = 1.0 / p.el_sec.sin().max(0.1);
            zwd_numer += mw * tropo;
            zwd_denom += mw * mw;
        }
    }
    let delta_zwd = if zwd_denom > 1e-8 { zwd_numer / zwd_denom } else { 0.0 };
    let atmo = StationAtmosphere {
        station_id: sec.id.clone(),
        zwd_m: 0.15 + delta_zwd,
        iono_slant_m: iono_map,
    };
    (atmo, fixed_count)
}

/// Context for double-difference baseline adjustment.
struct BaselineContext<'a> {
    master: &'a CorsStation,
    sec: &'a CorsStation,
    m_epoch: &'a EpochObs,
    s_epoch: &'a EpochObs,
    ephemerides: &'a [Ephemeris],
}

/// Reference satellite geometry and observations.
struct RefSatData<'a> {
    m_ref: &'a SatObs,
    s_ref: &'a SatObs,
    rho_ref_m: f64,
    rho_ref_s: f64,
}

/// Form double differences against highest elevation reference satellite.
fn form_dd_pairs(
    master: &CorsStation,
    m_epoch: &EpochObs,
    sec: &CorsStation,
    s_epoch: &EpochObs,
    ephemerides: &[Ephemeris],
) -> Vec<DdPairObs> {
    let ctx = BaselineContext { master, sec, m_epoch, s_epoch, ephemerides };
    let s_llh = ecef_to_llh(sec.pos_ecef);

    let common = find_common_sats(&ctx, s_llh);
    let ref_sat = match common.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)) {
        Some(r) => r.0,
        None => return Vec::new(),
    };

    build_dd_pair_list(&ctx, &common, ref_sat)
}

/// Find common satellites with valid dual-frequency tracking.
fn find_common_sats(
    ctx: &BaselineContext<'_>,
    s_llh: Vector3<f64>,
) -> Vec<(SatelliteId, f64)> {
    let mut common = Vec::new();
    let time = ctx.m_epoch.time;
    for m_sat in &ctx.m_epoch.satellites {
        if ctx.s_epoch.satellites.iter().any(|s| s.sat == m_sat.sat) {
            if let Some(eph) = ctx.ephemerides.iter().find(|e| e.sat() == m_sat.sat) {
                let (sat_p, _, _, _) = eph.position(time);
                let (_, el) = az_el(s_llh, ctx.sec.pos_ecef, sat_p);
                if el > 0.15 {
                    common.push((m_sat.sat, el));
                }
            }
        }
    }
    common
}

/// Build list of double-difference observations across satellites.
fn build_dd_pair_list(
    ctx: &BaselineContext<'_>,
    common: &[(SatelliteId, f64)],
    ref_sat: SatelliteId,
) -> Vec<DdPairObs> {
    let time = ctx.m_epoch.time;
    let (m_ref, s_ref) = match (get_sat(ctx.m_epoch, ref_sat), get_sat(ctx.s_epoch, ref_sat)) {
        (Some(m), Some(s)) => (m, s),
        _ => return Vec::new(),
    };
    let ref_eph = match ctx.ephemerides.iter().find(|e| e.sat() == ref_sat) {
        Some(e) => e,
        None => return Vec::new(),
    };
    let (p_ref, _, _, _) = ref_eph.position(time);
    let ref_data = RefSatData {
        m_ref,
        s_ref,
        rho_ref_m: (p_ref - ctx.master.pos_ecef).norm(),
        rho_ref_s: (p_ref - ctx.sec.pos_ecef).norm(),
    };

    let mut pairs = Vec::new();
    for &(sat, el_sec) in common {
        if sat == ref_sat { continue; }
        if let Some(pair) = make_single_dd_pair(ctx, &ref_data, sat, el_sec) {
            pairs.push(pair);
        }
    }
    pairs
}

/// Form single DD pair observation between satellite and reference.
fn make_single_dd_pair(
    ctx: &BaselineContext<'_>,
    ref_data: &RefSatData<'_>,
    sat: SatelliteId,
    el_sec: f64,
) -> Option<DdPairObs> {
    let (m_s, s_s) = (get_sat(ctx.m_epoch, sat)?, get_sat(ctx.s_epoch, sat)?);
    let sat_eph = ctx.ephemerides.iter().find(|e| e.sat() == sat)?;
    let (p_sat, _, _, _) = sat_eph.position(ctx.m_epoch.time);
    let rho_m = (p_sat - ctx.master.pos_ecef).norm();
    let rho_s = (p_sat - ctx.sec.pos_ecef).norm();
    let dd_rho = (rho_s - rho_m) - (ref_data.rho_ref_s - ref_data.rho_ref_m);

    let (mp1, mp2, mc1, mc2) = extract_l1_l2(m_s)?;
    let (sp1, sp2, sc1, sc2) = extract_l1_l2(s_s)?;
    let (mrp1, mrp2, mrc1, mrc2) = extract_l1_l2(ref_data.m_ref)?;
    let (srp1, srp2, src1, src2) = extract_l1_l2(ref_data.s_ref)?;

    Some(DdPairObs {
        sat,
        dd_rho,
        dd_phi1: (sc1 - mc1) - (src1 - mrc1),
        dd_phi2: (sc2 - mc2) - (src2 - mrc2),
        dd_p1: (sp1 - mp1) - (srp1 - mrp1),
        dd_p2: (sp2 - mp2) - (srp2 - mrp2),
        el_sec,
    })
}

/// Extract dual-frequency carrier phases and pseudoranges.
fn extract_l1_l2(sat_obs: &SatObs) -> Option<(f64, f64, f64, f64)> {
    let mut p1 = None;
    let mut p2 = None;
    let mut c1 = None;
    let mut c2 = None;
    for o in &sat_obs.observations {
        let b = o.code.signal.freq_band;
        match o.code.obs_type {
            gneiss_core::obs::ObsType::Pseudorange => {
                if b == 1 { p1 = Some(o.value); }
                else if b == 2 { p2 = Some(o.value); }
            }
            gneiss_core::obs::ObsType::CarrierPhase => {
                if b == 1 { c1 = Some(o.value); }
                else if b == 2 { c2 = Some(o.value); }
            }
            _ => {}
        }
    }
    Some((p1?, p2?, c1?, c2?))
}

fn get_sat(epoch: &EpochObs, sat: SatelliteId) -> Option<&SatObs> {
    epoch.satellites.iter().find(|s| s.sat == sat)
}

/// Resolve wide-lane and narrow-lane integer ambiguities for a DD pair.
fn resolve_dd_ambiguities(p: &DdPairObs) -> Option<(i32, i32)> {
    let (f1, f2) = satellite_frequencies(p.sat, 0);
    let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    let lam1 = c / f1;
    let lam2 = c / f2;
    let lam_w = c / (f1 - f2);

    let phi_w = p.dd_phi1 - p.dd_phi2;
    let p_nl = (f1 * p.dd_p1 + f2 * p.dd_p2) / (f1 + f2);
    let n_w = (phi_w - p_nl / lam_w).round() as i32;

    let f1_sq = f1 * f1;
    let f2_sq = f2 * f2;
    let phi_if = (f1_sq * p.dd_phi1 * lam1 - f2_sq * p.dd_phi2 * lam2) / (f1_sq - f2_sq);
    let lam_if = c / (f1 + f2);
    let n1_float = (phi_if - p.dd_rho) / lam_if - (f2 / (f1 - f2)) * n_w as f64;
    let n1 = n1_float.round() as i32;
    let n2 = n1 - n_w;
    Some((n1, n2))
}

/// Extract slant ionosphere and tropospheric residual from fixed pair.
fn extract_dd_atmosphere(p: &DdPairObs, n1: i32, n2: i32) -> (f64, f64) {
    let (f1, f2) = satellite_frequencies(p.sat, 0);
    let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    let lam1 = c / f1;
    let lam2 = c / f2;

    let phi_gf = lam1 * (p.dd_phi1 - n1 as f64) - lam2 * (p.dd_phi2 - n2 as f64);
    let gamma = (f1 / f2).powi(2);
    let iono = phi_gf / (gamma - 1.0);

    let f1_sq = f1 * f1;
    let f2_sq = f2 * f2;
    let phi_if_m = (f1_sq * (p.dd_phi1 - n1 as f64) * lam1 - f2_sq * (p.dd_phi2 - n2 as f64) * lam2)
        / (f1_sq - f2_sq);
    let tropo = phi_if_m - p.dd_rho;
    (iono, tropo)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::sat::Constellation;

    #[test]
    fn test_network_adjuster_initialization() {
        let master_id = "P181";
        let stations = vec![
            CorsStation {
                id: "P181".into(),
                pos_ecef: Vector3::new(-2697941.0, -4255089.0, 3898009.0),
                epochs: Vec::new(),
            },
            CorsStation {
                id: "OHLN".into(),
                pos_ecef: Vector3::new(-2686856.0, -4254625.0, 3905990.0),
                epochs: Vec::new(),
            },
        ];
        let adj = NetworkAdjuster::new(&stations, master_id).expect("init failed");
        assert_eq!(adj.baselines.len(), 1);
        assert!(adj.baselines[0].length_m > 10_000.0);
    }

    #[test]
    fn test_resolve_dd_ambiguities_synthetic() {
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let (f1, f2) = satellite_frequencies(sat, 0);
        let c = gneiss_core::constants::SPEED_OF_LIGHT_M_S;
        let lam1 = c / f1;
        let lam2 = c / f2;
        let true_n1 = 5;
        let true_n2 = 4;
        let true_rho = 100.0;

        let pair = DdPairObs {
            sat,
            dd_rho: true_rho,
            dd_phi1: true_rho / lam1 + true_n1 as f64,
            dd_phi2: true_rho / lam2 + true_n2 as f64,
            dd_p1: true_rho,
            dd_p2: true_rho,
            el_sec: 0.8,
        };
        let (n1, n2) = resolve_dd_ambiguities(&pair).expect("ambiguities");
        assert_eq!(n1, true_n1);
        assert_eq!(n2, true_n2);
    }
}
