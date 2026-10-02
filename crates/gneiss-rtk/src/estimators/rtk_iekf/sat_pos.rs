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
                gneiss_core::sat::Constellation::Qzss => 'J',
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
            // GLONASS broadcast ephemerides are NOT Kepler elements: they carry
            // position/velocity/acceleration in PZ-90 and are referenced to
            // GLONASS time (UTC(SU) + 3 h with leap seconds), and require
            // numerical integration rather than the Keplerian propagation used
            // for GPS/Galileo. Evaluating one at GPST misplaces the satellite by
            // tens of kilometres -- measured at ~20 km on RTK Explorer F9P, which
            // poisoned every double-difference fix in the set.
            //
            // GLONASS remains usable from precise (SP3) orbits, handled above.
            // Here we have no correct position, so contribute nothing rather than
            // a wrong one: a bogus position is worse than a missing satellite,
            // because it survives as a large innovation instead of a gap.
            if matches!(eph, Ephemeris::Glonass(_)) {
                continue;
            }
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

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::ephemeris::{BeidouEphemeris, GalileoEphemeris, GpsEphemeris, GlonassEphemeris};
    use gneiss_core::obs::{ObsCode, ObsType, Observation, SignalCode};
    use gneiss_core::sat::{Constellation, SatelliteId};

    /// Deterministic Keplerian record: only `m0` varies, so two records for
    /// the same satellite at different TOEs give visibly different positions.
    fn gps_eph(prn: u8, cons: Constellation, toe_week: u32, toe_tow: f64, m0: f64) -> GpsEphemeris {
        GpsEphemeris {
            sat: SatelliteId { constellation: cons, prn },
            toe: GpsTime::new(toe_week, toe_tow),
            toc: GpsTime::new(toe_week, toe_tow),
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0,
            e: 0.0,
            // sqrt(26 561 000 m): the GPS semi-major axis, in metres.
            sqrt_a: 5153.738,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 0.96,
            idot: 0.0,
            omega: 0.0,
            tgd: 0.0,
            iode: 1,
            iodc: 1,
        }
    }

    /// Galileo twin of [`gps_eph`] with the Galileo-specific tail fields.
    fn galileo_eph(prn: u8, toe_week: u32, toe_tow: f64, m0: f64) -> GalileoEphemeris {
        GalileoEphemeris {
            sat: SatelliteId { constellation: Constellation::Galileo, prn },
            toe: GpsTime::new(toe_week, toe_tow),
            toc: GpsTime::new(toe_week, toe_tow),
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0,
            e: 0.0,
            sqrt_a: 5440.0,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 0.96,
            idot: 0.0,
            omega: 0.0,
            bgd_e1_e5a: 0.0,
            bgd_e1_e5b: 0.0,
            iod_nav: 1,
        }
    }

    fn beidou_eph(prn: u8, toe_week: u32, toe_tow: f64) -> BeidouEphemeris {
        BeidouEphemeris {
            sat: SatelliteId { constellation: Constellation::Beidou, prn },
            toe: GpsTime::new(toe_week, toe_tow),
            toc: GpsTime::new(toe_week, toe_tow),
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 0.0,
            e: 0.0,
            sqrt_a: 5282.0,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 0.96,
            idot: 0.0,
            omega: 0.0,
            tgd1: 0.0,
            tgd2: 0.0,
            aode: 1,
            aodc: 1,
        }
    }

    fn glonass_eph(prn: u8, freq_num: i8) -> GlonassEphemeris {
        GlonassEphemeris {
            sat: SatelliteId { constellation: Constellation::Glonass, prn },
            toe: GpsTime::new(2200, 300_000.0),
            freq_num,
            tau_n: 0.0,
            gamma_n: 0.0,
            delta_tau_n: 0.0,
            x: 15_000.000,
            y: 10_000.000,
            z: 20_000.000,
            vx: 1.0,
            vy: 2.0,
            vz: 3.0,
            ax: 0.0,
            ay: 0.0,
            az: 0.0,
        }
    }

    /// One observation record on `band` of the given type.
    fn obs(band: u8, kind: ObsType, value: f64) -> Observation {
        Observation {
            code: ObsCode { obs_type: kind, signal: SignalCode { freq_band: band, attribute: 'C' } },
            value,
            lock_time: None,
            lli: Some(0),
        }
    }

    fn want9() -> SatelliteId {
        SatelliteId { constellation: Constellation::Gps, prn: 9 }
    }

// ---------------------------------------------------------------------------
// glo_freq_num
// ---------------------------------------------------------------------------

#[test]
fn glo_freq_num_reads_the_channel_from_the_matching_glonass_record() {
    let ephems = vec![
        Ephemeris::Glonass(glonass_eph(8, -3)),
        Ephemeris::Glonass(glonass_eph(9, 1)),
    ];
    let sat = SatelliteId { constellation: Constellation::Glonass, prn: 9 };
    assert_eq!(glo_freq_num(&ephems, sat), 1);
}

#[test]
fn glo_freq_num_is_zero_for_an_absent_or_non_glonass_satellite() {
    let ephems = vec![Ephemeris::Glonass(glonass_eph(8, -3))];
    // Present record is GLONASS R08; asking for R09 finds nothing.
    let absent = SatelliteId { constellation: Constellation::Glonass, prn: 9 };
    assert_eq!(glo_freq_num(&ephems, absent), 0);
    // A GPS request must not be answered by a GLONASS record.
    let gps = SatelliteId { constellation: Constellation::Gps, prn: 8 };
    assert_eq!(glo_freq_num(&ephems, gps), 0);
    assert_eq!(glo_freq_num(&[], absent), 0);
}

// ---------------------------------------------------------------------------
// broadcast_position_for
// ---------------------------------------------------------------------------

#[test]
fn broadcast_position_for_matches_on_constellation_as_well_as_prn() {
    // GPS PRN 5 and Galileo PRN 5 share a PRN. Asking for GPS must return
    // the GPS orbit (m0 = 0 -> satellite near perigee) and not the Galileo
    // one (m0 = 3.0 rad -> a completely different point).
    let ephems = vec![
        Ephemeris::Gps(gps_eph(5, Constellation::Gps, 2200, 300_000.0, 0.0)),
        Ephemeris::Galileo(galileo_eph(5, 2200, 300_000.0, 3.0)),
    ];
    let want = SatelliteId { constellation: Constellation::Gps, prn: 5 };
    let (p_gps, _) = broadcast_position_for(&ephems, &want, GpsTime::new(2200, 300_000.0))
        .expect("GPS orbit must resolve");
    let gal = Ephemeris::Galileo(galileo_eph(5, 2200, 300_000.0, 3.0));
    let (p_gal, _, _, _) = gal.position(GpsTime::new(2200, 300_000.0));
    assert!(
        (p_gps - p_gal).norm() > 1.0e6,
        "returned the Galileo orbit for a GPS query (delta {:.3} m)",
        (p_gps - p_gal).norm()
    );
}

#[test]
fn broadcast_position_for_returns_none_for_a_missing_satellite() {
    let ephems = vec![Ephemeris::Gps(gps_eph(5, Constellation::Gps, 2200, 300_000.0, 0.0))];
    let missing = SatelliteId { constellation: Constellation::Gps, prn: 6 };
    assert!(broadcast_position_for(&ephems, &missing, GpsTime::new(2200, 300_000.0)).is_none());
}

#[test]
fn broadcast_position_for_refuses_glonass_even_when_the_record_matches() {
    // GLONASS broadcast states are PZ-90 integrated positions referenced to
    // GLONASS time; evaluating one through the Keplerian path at GPST would
    // misplace it by tens of km, so the selector must bail out instead.
    let ephems = vec![Ephemeris::Glonass(glonass_eph(8, -3))];
    let want = SatelliteId { constellation: Constellation::Glonass, prn: 8 };
    assert!(broadcast_position_for(&ephems, &want, GpsTime::new(2200, 300_000.0)).is_none());
}

#[test]
fn broadcast_position_for_ignores_keplerless_constellations() {
    // BeiDou/QZSS variants carry no Keplerian toe in this selector's match
    // arm, so a BeiDou-only list resolves to nothing for a BeiDou query.
    let ephems = vec![Ephemeris::Beidou(beidou_eph(6, 2200, 300_000.0))];
    let want = SatelliteId { constellation: Constellation::Beidou, prn: 6 };
    assert!(broadcast_position_for(&ephems, &want, GpsTime::new(2200, 300_000.0)).is_none());
}

#[test]
fn broadcast_position_for_prefers_the_nearest_toe_in_the_same_week() {
    // Two GPS records for PRN 9: TOE 300 000 s (dt = 60 s from the query)
    // and TOE 320 000 s (dt = 20 000 s). The near record must win.
    let near = Ephemeris::Gps(gps_eph(9, Constellation::Gps, 2200, 300_000.0, 0.0));
    let far = Ephemeris::Gps(gps_eph(9, Constellation::Gps, 2200, 320_000.0, 2.5));
    let ephems = vec![far.clone(), near.clone()];
    let t = GpsTime::new(2200, 300_060.0);
    let (got, _) = broadcast_position_for(&ephems, &want9(), t).expect("must resolve");
    let (expect, _, _, _) = near.position(t);
    let (other, _, _, _) = far.position(t);
    assert!((got - expect).norm() < 1e-6, "nearest-TOE record must be selected");
    // The two records differ by m0 = 2.5 rad, i.e. tens of megametres: the
    // returned point must not merely be *some* valid orbit for this satellite.
    assert!(
        (got - other).norm() > 1.0e6,
        "returned the far-TOE record instead of the near one"
    );
}

// ---------------------------------------------------------------------------
// extract_sat_positions
// ---------------------------------------------------------------------------

fn epoch_with(prns: &[(Constellation, u8)]) -> EpochObs {
    EpochObs {
        time: GpsTime::new(2200, 300_000.0),
        satellites: prns.iter().map(|(c, p)| SatObs {
            sat: SatelliteId { constellation: *c, prn: *p },
            observations: vec![obs(1, ObsType::Pseudorange, 22_000_000.0)],
        }).collect(),
    }
}

#[test]
fn extract_sat_positions_drops_glonass_without_a_correct_time_system() {
    // A GLONASS record exists and matches the observation, but the broadcast
    // path has no correct way to evaluate it: contributing nothing is the
    // documented choice ("a bogus position is worse than a missing satellite").
    let rover = epoch_with(&[(Constellation::Glonass, 8), (Constellation::Gps, 5)]);
    let ephems = vec![
        Ephemeris::Glonass(glonass_eph(8, -3)),
        Ephemeris::Gps(gps_eph(5, Constellation::Gps, 2200, 300_000.0, 0.0)),
    ];
    let rx = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);
    let out = extract_sat_positions(&rover, &ephems, rx, 0.0, None);
    assert_eq!(out.len(), 1, "only the GPS satellite may contribute: {out:?}");
    assert_eq!(out[0].0, SatelliteId { constellation: Constellation::Gps, prn: 5 });
}

#[test]
fn extract_sat_positions_honours_the_elevation_cut_off() {
    let rover = epoch_with(&[(Constellation::Gps, 5)]);
    let ephems = vec![Ephemeris::Gps(gps_eph(5, Constellation::Gps, 2200, 300_000.0, 0.0))];
    let rx = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);
    // min_el = 0 (everything kept) vs min_el just under pi/2 (nothing can
    // clear a satellite that is not overhead, so the set must be empty only
    // if it genuinely is low). Pin it down by first reading the elevation.
    let with_all = extract_sat_positions(&rover, &ephems, rx, -1.0, None);
    assert_eq!(with_all.len(), 1);
    let pos = with_all[0].1;
    let llh = gneiss_core::coords::ecef_to_llh(rx);
    let (_, el) = gneiss_core::coords::az_el(llh, rx, pos);
    // An impossible cut-off above this satellite's elevation must reject it.
    let filtered = extract_sat_positions(&rover, &ephems, rx, el + 0.05, None);
    assert!(filtered.is_empty(), "satellite at el={el} must fail a {}-rad mask", el + 0.05);
    // A cut-off below it must keep it.
    assert_eq!(extract_sat_positions(&rover, &ephems, rx, el - 0.05, None).len(), 1);
}

#[test]
fn extract_sat_positions_returns_nothing_for_a_satellite_with_no_ephemeris() {
    let rover = epoch_with(&[(Constellation::Gps, 31)]);
    let ephems = vec![Ephemeris::Gps(gps_eph(5, Constellation::Gps, 2200, 300_000.0, 0.0))];
    let rx = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);
    assert!(extract_sat_positions(&rover, &ephems, rx, 0.0, None).is_empty());
}

// ---------------------------------------------------------------------------
// compute_signal_sat_pos
// ---------------------------------------------------------------------------

#[test]
fn signal_travel_time_uses_the_measured_pseudorange() {
    // tau = P / c. With P = 20 000 km = 2.0e7 m and c = 299 792 458 m/s,
    // tau = 0.0667... s. A 1000 km longer code range must lengthen tau by
    // exactly 1000e3 / c = 3.3356... ms.
    let sv = SatelliteId { constellation: Constellation::Gps, prn: 5 };
    let eph = Ephemeris::Gps(gps_eph(5, Constellation::Gps, 2200, 300_000.0, 0.0));
    let t = GpsTime::new(2200, 300_000.0);
    let mk = |p1: f64| SatObs {
        sat: sv,
        observations: vec![obs(1, ObsType::Pseudorange, p1)],
    };
    let a = compute_signal_sat_pos(&mk(2.0e7), &eph, t);
    let b = compute_signal_sat_pos(&mk(2.1e7), &eph, t);
    // Over 30 minutes of Kepler propagation an extra 3.336 ms of transmit
    // time moves the satellite by roughly v * dt; GPS orbital speed is
    // ~3.87 km/s, so > 10 m of separation is the expected signature.
    assert!(
        (a - b).norm() > 10.0,
        "pseudorange must drive transmit time; separation was {:.3} m",
        (a - b).norm()
    );
    assert!((2.0e7 / SPEED_OF_LIGHT_M_S - 0.06671).abs() < 1e-5);
}

#[test]
fn signal_position_falls_back_to_20_000_km_when_no_code_is_present() {
    // No observable on L1 or L2: the default range is 20 000 000 m, so the
    // computed position must equal the L1-with-20 000 km case exactly.
    let sv = SatelliteId { constellation: Constellation::Gps, prn: 5 };
    let eph = Ephemeris::Gps(gps_eph(5, Constellation::Gps, 2200, 300_000.0, 0.0));
    let t = GpsTime::new(2200, 300_000.0);
    let empty = SatObs { sat: sv, observations: Vec::new() };
    let l1 = SatObs {
        sat: sv,
        observations: vec![obs(1, ObsType::Pseudorange, 20_000_000.0)],
    };
    let a = compute_signal_sat_pos(&empty, &eph, t);
    let b = compute_signal_sat_pos(&l1, &eph, t);
    assert!((a - b).norm() < 1e-6, "default range must match 20 000 km");
}}
