//! Reference satellite selection, elevation hysteresis, and constellation pairing.

use std::collections::HashMap;
use nalgebra::Vector3;
use gneiss_core::coords::{az_el, ecef_to_llh};
use gneiss_core::sat::{Constellation, SatelliteId};

/// Map a [`SatelliteId`] to its unique u16 tracking/DD identifier.
///
/// GPS satellites (PRN 1..32) map to 1..32.
/// QZSS satellites (PRN 1..10) map to 193..202 (standard NMEA/RTCM offset 192 + PRN).
/// Other constellations map directly to their PRN.
#[inline]
pub fn sat_to_prn_u16(sat: SatelliteId) -> u16 {
    match sat.constellation {
        Constellation::Qzss => 192 + sat.prn as u16,
        _ => sat.prn as u16,
    }
}

/// Convert a constellation ID and u16 identifier back into a canonical [`SatelliteId`].
pub fn prn_u16_to_sat(const_id: u8, prn_u16: u16) -> SatelliteId {
    if const_id == Constellation::Gps as u8 && prn_u16 >= 193 {
        SatelliteId {
            constellation: Constellation::Qzss,
            prn: (prn_u16 - 192) as u8,
        }
    } else {
        let constellation = match const_id {
            0 => Constellation::Gps,
            1 => Constellation::Glonass,
            2 => Constellation::Galileo,
            3 => Constellation::Beidou,
            5 => Constellation::Qzss,
            _ => Constellation::Gps,
        };
        SatelliteId {
            constellation,
            prn: prn_u16 as u8,
        }
    }
}

/// Check if an observation [`SatelliteId`] matches a given constellation and u16 identifier.
pub fn sat_matches_id(sat: SatelliteId, const_id: u8, prn_u16: u16) -> bool {
    if const_id == Constellation::Gps as u8 {
        if prn_u16 >= 193 {
            sat.constellation == Constellation::Qzss && sat.prn as u16 == prn_u16 - 192
        } else {
            sat.constellation == Constellation::Gps && sat.prn as u16 == prn_u16
        }
    } else {
        sat.constellation as u8 == const_id && sat.prn as u16 == prn_u16
    }
}

/// Constellations that participate in DD formation, sorted by id.
///
/// QZSS is merged with GPS (constellation 0) because both share identical CDMA frequencies
/// (L1, L2, L5) and receiver tracking hardware. GLONASS requires the opt-in FDMA policy.
pub fn select_constellations(sat_info: &[(SatelliteId, Vector3<f64>)], glo: bool) -> Vec<u8> {
    let mut v = Vec::new();
    let has_gps_or_qzss = sat_info.iter().any(|(s, _)| {
        s.constellation == Constellation::Gps || s.constellation == Constellation::Qzss
    });
    if has_gps_or_qzss {
        v.push(Constellation::Gps as u8);
    }
    for (s, _) in sat_info {
        let c = s.constellation as u8;
        if c == Constellation::Gps as u8 || c == Constellation::Qzss as u8 {
            continue;
        }
        if c == Constellation::Glonass as u8 && !glo {
            continue;
        }
        if !v.contains(&c) {
            v.push(c);
        }
    }
    v.sort_unstable();
    v
}

/// BeiDou Geostationary Orbit (GEO) satellites have PRNs 1-5 and 59-63.
/// Due to near-stationary geometry, they provide minimal line-of-sight angular
/// velocity and should be demoted as reference pivots when IGSO/MEOs are available.
pub fn is_beidou_geo(sat: SatelliteId) -> bool {
    sat.constellation == Constellation::Beidou
        && ((1..=5).contains(&sat.prn) || (59..=63).contains(&sat.prn))
}

fn sat_tier1_eligible(r: &gneiss_core::obs::SatObs, b: &gneiss_core::obs::SatObs) -> bool {
    let r_cp1 = r.get_observable_phase(1).is_some();
    let r_cp2 = r.get_observable_phase(2).is_some()
        || r.get_observable_phase(6).is_some()
        || r.get_observable_phase(7).is_some();
    let b_cp1 = b.get_observable_phase(1).is_some();
    let b_cp2 = b.get_observable_phase(2).is_some()
        || b.get_observable_phase(6).is_some()
        || b.get_observable_phase(7).is_some();
    let no_slip = r.get_lli(1).unwrap_or(0) & 1 == 0;
    let good_snr = r.get_snr(1).unwrap_or(0) >= 30;
    r_cp1 && r_cp2 && b_cp1 && b_cp2 && no_slip && good_snr
}

fn sat_tier2_eligible(r: &gneiss_core::obs::SatObs, b: &gneiss_core::obs::SatObs) -> bool {
    let r_cp = [1, 2, 6, 7].iter().any(|&band| r.get_observable_phase(band).is_some());
    let b_cp = [1, 2, 6, 7].iter().any(|&band| b.get_observable_phase(band).is_some());
    let max_snr = [1, 2, 6, 7].iter().filter_map(|&band| r.get_snr(band)).max().unwrap_or(0);
    r_cp && b_cp && max_snr >= 25
}

fn classify_candidate(sat: SatelliteId, r: &gneiss_core::obs::SatObs, b: &gneiss_core::obs::SatObs) -> usize {
    if sat_tier1_eligible(r, b) {
        if is_beidou_geo(sat) {
            1
        } else {
            0
        }
    } else if sat_tier2_eligible(r, b) {
        1
    } else {
        2
    }
}

/// Filter candidate reference satellites: prioritize satellites tracked on both
/// rover and base with valid dual-frequency carrier phase and SNR >= 30 dB-Hz.
pub fn filter_reference_candidates(
    const_sats: &[(SatelliteId, Vector3<f64>)],
    rover: &gneiss_core::obs::EpochObs,
    base: &gneiss_core::obs::EpochObs,
    const_id: u8,
) -> Vec<(SatelliteId, Vector3<f64>)> {
    let mut tiers: [Vec<(SatelliteId, Vector3<f64>)>; 3] = Default::default();
    for &(sat, pos) in const_sats {
        let prn = sat_to_prn_u16(sat);
        let rov = rover.satellites.iter().find(|s| sat_matches_id(s.sat, const_id, prn));
        let bas = base.satellites.iter().find(|s| sat_matches_id(s.sat, const_id, prn));
        if let (Some(r), Some(b)) = (rov, bas) {
            let tier = classify_candidate(sat, r, b);
            if tier < 2 {
                tiers[tier].push((sat, pos));
            }
            tiers[2].push((sat, pos));
        }
    }
    for t in tiers {
        if !t.is_empty() {
            return t;
        }
    }
    const_sats.to_vec()
}

/// Select reference satellite for a constellation with elevation hysteresis.
pub fn select_ref_sat_with_hysteresis(
    const_id: u8,
    sats: &[(SatelliteId, Vector3<f64>)],
    rx_pos: Vector3<f64>,
    ref_sats: &mut HashMap<u8, u16>,
) -> u16 {
    let rx_llh = ecef_to_llh(rx_pos);
    let best = sats.iter().max_by(|a, b| {
        let (_, ea) = az_el(rx_llh, rx_pos, a.1);
        let (_, eb) = az_el(rx_llh, rx_pos, b.1);
        ea.partial_cmp(&eb).unwrap_or(std::cmp::Ordering::Equal)
    });
    let Some((best_sv, best_pos)) = best else { return 0 };
    let (_, best_el) = az_el(rx_llh, rx_pos, *best_pos);
    if let Some(&prev) = ref_sats.get(&const_id) {
        if let Some((_, prev_p)) = sats.iter().find(|(s, _)| sat_to_prn_u16(*s) == prev) {
            let (_, prev_el) = az_el(rx_llh, rx_pos, *prev_p);
            if prev_el >= 0.35 && prev_el >= best_el - 0.26 {
                return prev;
            }
        }
    }
    let chosen = sat_to_prn_u16(*best_sv);
    ref_sats.insert(const_id, chosen);
    chosen
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
    use gneiss_core::time::GpsTime;

    #[test]
    fn test_sat_to_prn_u16_gps_and_qzss() {
        let gps1 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let gps32 = SatelliteId { constellation: Constellation::Gps, prn: 32 };
        let qzss2 = SatelliteId { constellation: Constellation::Qzss, prn: 2 };
        let gal5 = SatelliteId { constellation: Constellation::Galileo, prn: 5 };

        assert_eq!(sat_to_prn_u16(gps1), 1);
        assert_eq!(sat_to_prn_u16(gps32), 32);
        assert_eq!(sat_to_prn_u16(qzss2), 194);
        assert_eq!(sat_to_prn_u16(gal5), 5);
    }

    #[test]
    fn test_sat_matches_id_and_roundtrip() {
        let gps2 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let qzss2 = SatelliteId { constellation: Constellation::Qzss, prn: 2 };

        assert!(sat_matches_id(gps2, 0, 2));
        assert!(!sat_matches_id(gps2, 0, 194));
        assert!(sat_matches_id(qzss2, 0, 194));
        assert!(!sat_matches_id(qzss2, 0, 2));

        assert_eq!(prn_u16_to_sat(0, 2), gps2);
        assert_eq!(prn_u16_to_sat(0, 194), qzss2);
    }

    #[test]
    fn test_select_constellations_merges_qzss_into_gps() {
        let sats = vec![
            (SatelliteId { constellation: Constellation::Qzss, prn: 2 }, Vector3::zeros()),
            (SatelliteId { constellation: Constellation::Galileo, prn: 1 }, Vector3::zeros()),
        ];
        let c = select_constellations(&sats, false);
        assert_eq!(c, vec![0, 2]);
    }

    // -----------------------------------------------------------------
    // Helpers: build satellites at a chosen elevation from a real station
    // -----------------------------------------------------------------

    const RX: Vector3<f64> = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);

    fn obs(band: u8, kind: ObsType, value: f64) -> Observation {
        Observation {
            code: ObsCode { obs_type: kind, signal: SignalCode { freq_band: band, attribute: 'C' } },
            value,
            lock_time: None,
            lli: Some(0),
        }
    }

    fn phase(band: u8, lli: u8) -> Observation {
        let mut o = obs(band, ObsType::CarrierPhase, 1.0e8);
        o.lli = Some(lli);
        o
    }

    fn snr(band: u8, dbhz: u8) -> Observation {
        obs(band, ObsType::Snr, f64::from(dbhz))
    }

    /// Satellite observed on rover and base with the requested phase/SNR mix.
    fn sat_obs(sat: SatelliteId, bands: &[u8], lli1: u8, snr1: u8) -> SatObs {
        let mut o: Vec<Observation> = bands.iter().map(|b| phase(*b, lli1)).collect();
        o.push(snr(1, snr1));
        SatObs { sat, observations: o }
    }

    fn epoch(sats: Vec<SatObs>) -> EpochObs {
        EpochObs { time: GpsTime::new(2200, 300_000.0), satellites: sats }
    }

    /// Satellite position placed at `el_rad` elevation in a given azimuth
    /// (radians) as seen from RX.
    fn sat_at(el_rad: f64, az_rad: f64) -> Vector3<f64> {
        let (east, north, up) = enu_basis(RX);
        let horiz = north * el_rad.cos() * az_rad.cos() + east * el_rad.cos() * az_rad.sin();
        RX + (up * el_rad.sin() + horiz).normalize() * 2.2e7
    }

    /// Local east / north / up unit basis at `pos` (ECEF), from the geodetic
    /// latitude and longitude (radians) returned by `ecef_to_llh`.
    fn enu_basis(pos: Vector3<f64>) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>) {
        let llh = ecef_to_llh(pos);
        let (sp_lat, cp_lat) = llh.x.sin_cos();
        let (sp_lon, cp_lon) = llh.y.sin_cos();
        let east = Vector3::new(-sp_lon, cp_lon, 0.0);
        let north = Vector3::new(-sp_lat * cp_lon, -sp_lat * sp_lon, cp_lat);
        let up = Vector3::new(cp_lat * cp_lon, cp_lat * sp_lon, sp_lat);
        (east, north, up)
    }

    fn gps(prn: u8) -> SatelliteId {
        SatelliteId { constellation: Constellation::Gps, prn }
    }

    fn pair(sat: SatelliteId, pos: Vector3<f64>) -> (SatelliteId, Vector3<f64>) {
        (sat, pos)
    }

    // -----------------------------------------------------------------
    // is_beidou_geo
    // -----------------------------------------------------------------

    #[test]
    fn beidou_geo_covers_prn_1_to_5_and_59_to_63_only() {
        for prn in [1u8, 3, 5, 59, 61, 63] {
            assert!(is_beidou_geo(SatelliteId { constellation: Constellation::Beidou, prn }), "PRN {prn}");
        }
        for prn in [0u8, 6, 30, 58, 64] {
            assert!(!is_beidou_geo(SatelliteId { constellation: Constellation::Beidou, prn }), "PRN {prn}");
        }
        // Same PRN in another constellation is never GEO BeiDou.
        assert!(!is_beidou_geo(gps(1)));
    }

    // -----------------------------------------------------------------
    // filter_reference_candidates
    // -----------------------------------------------------------------

    #[test]
    fn tier1_needs_dual_frequency_phase_on_both_receivers_without_slip() {
        let s = gps(6);
        let clean = sat_obs(s, &[1, 2], 0, 45);
        let pos = sat_at(1.0, 0.3);
        let cands = vec![pair(s, pos)];
        // Tier 1 requires rover AND base both dual-frequency with phase and
        // no L1 slip, so a clean pair must be returned.
        let got = filter_reference_candidates(&cands, &epoch(vec![clean.clone()]), &epoch(vec![clean.clone()]), 0);
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].0, s);
    }

    #[test]
    fn tier2_catches_a_single_frequency_arc_but_never_a_snr_poor_one() {
        let s = gps(7);
        let pos = sat_at(1.0, 0.0);
        let cands = vec![pair(s, pos)];
        // L1 phase only -> not tier 1, but L1 SNR 30 >= 25 -> tier 2 kept.
        let single = sat_obs(s, &[1], 0, 30);
        let got = filter_reference_candidates(&cands, &epoch(vec![single.clone()]), &epoch(vec![single]), 0);
        assert_eq!(got.len(), 1, "single-frequency arc is still usable");
        // L1 phase with SNR 20 < 25 -> no tier at all -> tier-2 fallback keeps
        // every satellite observed on both sides, which is still this one.
        let weak = sat_obs(s, &[1], 0, 20);
        let fallback = filter_reference_candidates(&cands, &epoch(vec![weak.clone()]), &epoch(vec![weak]), 0);
        assert_eq!(fallback.len(), 1);
    }

    #[test]
    fn nothing_qualifying_falls_back_to_the_whole_candidate_list() {
        // The base has no observation for this satellite, so no tier can be
        // filled. The documented fallback is to return the candidate list
        // unchanged rather than to yield an empty pivot.
        let s = gps(9);
        let cands = vec![pair(s, sat_at(1.0, 0.0))];
        let rover = epoch(vec![sat_obs(s, &[1, 2], 0, 45)]);
        let base = EpochObs { time: rover.time, satellites: Vec::new() };
        assert_eq!(filter_reference_candidates(&cands, &rover, &base, 0), cands);
    }

    #[test]
    fn an_empty_candidate_list_falls_through_to_a_copy_of_the_input() {
        assert!(filter_reference_candidates(&[], &epoch(vec![]), &epoch(vec![]), 0).is_empty());
    }

    #[test]
    fn beidou_geo_is_demoted_below_an_ordinary_igso_of_equal_quality() {
        // Both are BeiDou (constellation id 3) with identical clean dual-
        // frequency arcs, so only the GEO demotion separates them.
        let geo = SatelliteId { constellation: Constellation::Beidou, prn: 3 };
        let igso = SatelliteId { constellation: Constellation::Beidou, prn: 12 };
        let cands = vec![pair(geo, sat_at(1.0, 0.0)), pair(igso, sat_at(1.0, 0.5))];
        let g = sat_obs(geo, &[1, 2], 0, 45);
        let i = sat_obs(igso, &[1, 2], 0, 45);
        let got = filter_reference_candidates(
            &cands,
            &epoch(vec![g.clone(), i.clone()]),
            &epoch(vec![g, i]),
            3,
        );
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].0, igso, "GEO must never be the first choice");
    }

    // -----------------------------------------------------------------
    // select_ref_sat_with_hysteresis
    // -----------------------------------------------------------------

    #[test]
    fn an_empty_geometry_selects_no_reference() {
        let mut refs = HashMap::new();
        assert_eq!(select_ref_sat_with_hysteresis(0, &[], RX, &mut refs), 0);
        assert!(refs.is_empty());
    }

    #[test]
    fn the_highest_elevation_satellite_is_chosen_and_remembered() {
        let high = pair(gps(6), sat_at(1.2, 0.0));
        let low = pair(gps(7), sat_at(0.3, 1.0));
        let mut refs = HashMap::new();
        let got = select_ref_sat_with_hysteresis(0, &[low, high], RX, &mut refs);
        assert_eq!(got, 6);
        assert_eq!(refs.get(&0), Some(&6));
    }

    #[test]
    fn a_still_usable_previous_reference_is_kept_over_a_marginal_gain() {
        // prev at 0.40 rad; best at 0.50 rad. Hysteresis requires
        // prev_el >= 0.35 AND prev_el >= best_el - 0.26 = 0.24, so 0.40
        // satisfies both and the reference must NOT move.
        let prev = pair(gps(6), sat_at(0.40, 0.2));
        let best = pair(gps(7), sat_at(0.50, 1.0));
        let mut refs = HashMap::new();
        refs.insert(0u8, 6u16);
        let sats = vec![prev, best];
        assert_eq!(select_ref_sat_with_hysteresis(0, &sats, RX, &mut refs), 6);
        assert_eq!(refs.get(&0), Some(&6), "a kept reference is not rewritten");
    }

    #[test]
    fn a_reference_that_has_set_below_the_floor_is_replaced() {
        // prev at 0.20 rad < the 0.35 rad floor: switch, whatever the gain.
        let prev = pair(gps(6), sat_at(0.20, 0.2));
        let best = pair(gps(7), sat_at(0.22, 1.0));
        let mut refs = HashMap::new();
        refs.insert(0u8, 6u16);
        let sats = vec![prev, best];
        assert_eq!(select_ref_sat_with_hysteresis(0, &sats, RX, &mut refs), 7);
        assert_eq!(refs.get(&0), Some(&7));
    }

    #[test]
    fn a_reference_that_has_set_disappeared_is_replaced() {
        let mut refs = HashMap::new();
        refs.insert(0u8, 6u16);
        let sats = vec![pair(gps(7), sat_at(1.0, 0.0))];
        assert_eq!(select_ref_sat_with_hysteresis(0, &sats, RX, &mut refs), 7);
    }

    #[test]
    fn handover_releases_the_incumbent_only_past_both_margins() {
        // best at 0.61 rad. The incumbent survives while
        //   prev_el >= 0.35  AND  prev_el >= 0.61 - 0.26 = 0.35
        // so 0.40 rad is held and 0.30 rad is released.
        let best = pair(gps(7), sat_at(0.61, 1.0));
        let held = pair(gps(6), sat_at(0.40, 0.2));
        let mut refs = HashMap::new();
        refs.insert(0u8, 6u16);
        assert_eq!(select_ref_sat_with_hysteresis(0, &[held, best], RX, &mut refs), 6);
        let dropped = pair(gps(6), sat_at(0.30, 0.2));
        let mut refs2 = HashMap::new();
        refs2.insert(0u8, 6u16);
        assert_eq!(select_ref_sat_with_hysteresis(0, &[dropped, best], RX, &mut refs2), 7);
    }

    #[test]
    fn reference_selection_is_per_constellation() {
        let mut refs = HashMap::new();
        let g = vec![pair(gps(6), sat_at(1.0, 0.0))];
        let e = vec![pair(SatelliteId { constellation: Constellation::Galileo, prn: 11 }, sat_at(1.0, 0.0))];
        assert_eq!(select_ref_sat_with_hysteresis(0, &g, RX, &mut refs), 6);
        assert_eq!(select_ref_sat_with_hysteresis(2, &e, RX, &mut refs), 11);
        assert_eq!(refs.len(), 2, "each constellation keeps its own pivot");
    }
}
