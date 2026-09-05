//! Satellite position calculation and broadcast orbit extraction.

use nalgebra::Vector3;
use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::time::GpsTime;
use super::sat_pco;

/// GLONASS FDMA frequency channel number for a satellite, from its
/// broadcast ephemeris (0 when unknown -> nominal frequency).
pub(crate) fn glo_freq_num(ephems: &[Ephemeris], sat_id: gneiss_core::sat::SatelliteId) -> i8 {
    for e in ephems {
        if let Ephemeris::Glonass(g) = e {
            if g.sat == sat_id {
                return g.freq_num;
            }
        }
    }
    0
}

pub(crate) fn extract_sat_positions(
    rover: &EpochObs,
    ephems: &[Ephemeris],
    rx_pos: Vector3<f64>,
    min_el: f64,
    precise_orbits: Option<&std::sync::Arc<gneiss_parsers::precise_orbit::PreciseOrbit>>,
) -> Vec<(gneiss_core::sat::SatelliteId, Vector3<f64>)> {
    let mut out = Vec::new();
    let rx_llh = gneiss_core::coords::ecef_to_llh(rx_pos);

    for s in &rover.satellites {
        if let Some(precise) = precise_orbits {
            let sys_char = match s.sat.constellation {
                gneiss_core::sat::Constellation::Gps => 'G',
                gneiss_core::sat::Constellation::Glonass => 'R',
                gneiss_core::sat::Constellation::Galileo => 'E',
                gneiss_core::sat::Constellation::Beidou => 'C',
                _ => 'G',
            };
            let sv_name = format!("{}{:02}", sys_char, s.sat.prn);
            if std::env::var("GNEISS_SP3_PROBE").is_ok() {
                static FIRST: std::sync::atomic::AtomicUsize =
                    std::sync::atomic::AtomicUsize::new(0);
                let n = FIRST.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if n < 4 {
                    let hit = precise.position_at(&sv_name, rover.time).is_some();
                    eprintln!("SP3-PROBE[{}] {} -> {}", n, sv_name, hit);
                }
            }
            if let Some((pos_com, _clk)) =
                precise.position_at_with_hint(&sv_name, rover.time, Some(rx_pos))
            {
                let tau = (rx_pos - pos_com).norm()
                    / gneiss_core::constants::SPEED_OF_LIGHT_M_S;
                let wt = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * tau;
                let (sw, cw) = libm::sincos(wt);
                let rotated = Vector3::new(
                    pos_com.x * cw + pos_com.y * sw,
                    -pos_com.x * sw + pos_com.y * cw,
                    pos_com.z,
                );
                let pos = sat_pco::apply_sat_pco_z(rotated, s.sat.prn as u16);
                let (_az, el) = gneiss_core::coords::az_el(rx_llh, rx_pos, pos);
                if el >= min_el {
                    out.push((s.sat, pos));
                }
                continue;
            }
        }

        let eph_opt = ephems
            .iter()
            .filter(|e| e.sat() == s.sat)
            .min_by(|a, b| {
                (a.toe().tow - rover.time.tow)
                    .abs()
                    .partial_cmp(&(b.toe().tow - rover.time.tow).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        if let Some(eph) = eph_opt {
            let sat_p = compute_signal_sat_pos(s, eph, rover.time);
            let (_az, el) = gneiss_core::coords::az_el(rx_llh, rx_pos, sat_p);
            if el >= min_el {
                out.push((s.sat, sat_p));
            }
        }
    }
    out
}

/// Nearest-TOE broadcast position for a satellite (pipeline delegate).
pub fn broadcast_position_for(
    ephems: &[Ephemeris],
    sv: &gneiss_core::sat::SatelliteId,
    t: GpsTime,
) -> Option<(Vector3<f64>, f64)> {
    let mut best: Option<(&Ephemeris, f64)> = None;
    for cand in ephems {
        if cand.sat().constellation != sv.constellation || cand.sat().prn != sv.prn {
            continue;
        }
        let toe = match cand {
            Ephemeris::Gps(g) => g.toe,
            Ephemeris::Galileo(g) => g.toe,
            Ephemeris::Glonass(_) => return None,
            _ => continue,
        };
        let dt = if toe.week == t.week { (toe.tow - t.tow).abs() } else { f64::INFINITY };
        if best.is_none_or(|(_, d)| dt < d) {
            best = Some((cand, dt));
        }
    }
    let (eph, _) = best?;
    let (p, _, clk, _) = match eph {
        Ephemeris::Gps(g) => g.position(t),
        Ephemeris::Galileo(g) => g.position(t),
        _ => return None,
    };
    Some((p, clk))
}

fn compute_signal_sat_pos(s: &SatObs, eph: &Ephemeris, time: gneiss_core::time::GpsTime) -> Vector3<f64> {
    let pr_m = s.get_observable(1).or_else(|| s.get_observable(2)).unwrap_or(20_000_000.0);
    let tau = pr_m / SPEED_OF_LIGHT_M_S;
    let t_tx = time - tau;
    let (_, _, sat_clk_err_rough, _) = eph.position(t_tx);
    let t_tx_true = t_tx - sat_clk_err_rough;
    let (sat_p, _, _, _) = eph.position(t_tx_true);

    let omega_tau = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * tau;
    let cos_wt = omega_tau.cos();
    let sin_wt = omega_tau.sin();
    Vector3::new(
        sat_p.x * cos_wt + sat_p.y * sin_wt,
        -sat_p.x * sin_wt + sat_p.y * cos_wt,
        sat_p.z,
    )
}
