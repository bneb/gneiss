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

fn sat_tier1_eligible(r: &gneiss_core::obs::SatObs, b: &gneiss_core::obs::SatObs) -> bool {
    let r_cp1 = r.get_observable_phase(1).is_some();
    let r_cp2 = r.get_observable_phase(2).is_some() || r.get_observable_phase(7).is_some();
    let b_cp1 = b.get_observable_phase(1).is_some();
    let b_cp2 = b.get_observable_phase(2).is_some() || b.get_observable_phase(7).is_some();
    let no_slip = r.get_lli(1).unwrap_or(0) & 1 == 0;
    let good_snr = r.get_snr(1).unwrap_or(0) >= 30;
    r_cp1 && r_cp2 && b_cp1 && b_cp2 && no_slip && good_snr
}

fn sat_tier2_eligible(r: &gneiss_core::obs::SatObs, b: &gneiss_core::obs::SatObs) -> bool {
    let r_cp = r.get_observable_phase(1).is_some() || r.get_observable_phase(2).is_some();
    let b_cp = b.get_observable_phase(1).is_some() || b.get_observable_phase(2).is_some();
    let max_snr = [1, 2, 7].iter().filter_map(|&band| r.get_snr(band)).max().unwrap_or(0);
    r_cp && b_cp && max_snr >= 25
}

fn classify_candidate(r: &gneiss_core::obs::SatObs, b: &gneiss_core::obs::SatObs) -> usize {
    if sat_tier1_eligible(r, b) {
        0
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
            let tier = classify_candidate(r, b);
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
}
