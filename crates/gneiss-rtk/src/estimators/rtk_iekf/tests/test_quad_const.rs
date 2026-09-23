#![allow(clippy::unwrap_used)]

use super::*;
use gneiss_core::obs::{Observation, ObsCode, SatObs};
use gneiss_core::sat::{Constellation, SatelliteId};

#[test]
fn test_canonical_constellation_bands() {
    let gps_bands = formation::canonical_bands_for_constellation(Constellation::Gps);
    assert_eq!(gps_bands, &[1, 2, 5]);

    let gal_bands = formation::canonical_bands_for_constellation(Constellation::Galileo);
    assert_eq!(gal_bands, &[1, 5, 7]);

    let bds_bands = formation::canonical_bands_for_constellation(Constellation::Beidou);
    assert_eq!(bds_bands, &[1, 6, 7]);

    let glo_bands = formation::canonical_bands_for_constellation(Constellation::Glonass);
    assert_eq!(glo_bands, &[1, 2]);
}

#[test]
fn test_beidou_geo_identification() {
    for prn in 1..=5 {
        let sat = SatelliteId { constellation: Constellation::Beidou, prn };
        assert!(ref_sat::is_beidou_geo(sat), "C0{prn} must be identified as GEO");
    }
    for prn in 59..=63 {
        let sat = SatelliteId { constellation: Constellation::Beidou, prn };
        assert!(ref_sat::is_beidou_geo(sat), "C{prn} must be identified as GEO");
    }
    for prn in [6, 7, 8, 9, 10, 13, 16, 19, 21, 37] {
        let sat = SatelliteId { constellation: Constellation::Beidou, prn };
        assert!(!ref_sat::is_beidou_geo(sat), "C{prn} is IGSO/MEO, not GEO");
    }
}

#[test]
fn test_beidou_geo_demotion_in_reference_selection() {
    let geo = SatelliteId { constellation: Constellation::Beidou, prn: 1 };
    let igso = SatelliteId { constellation: Constellation::Beidou, prn: 6 };
    let pos_geo = Vector3::new(3.8e7, 1.8e7, 0.0);
    let pos_igso = Vector3::new(3.5e7, 1.5e7, 2.0e7);

    use std::str::FromStr;
    let mk_obs = |sat: SatelliteId| SatObs {
        sat,
        observations: vec![
            Observation { code: ObsCode::from_str("C1I").unwrap(), value: 3.8e7, lock_time: None, lli: None },
            Observation { code: ObsCode::from_str("L1I").unwrap(), value: 2.0e8, lock_time: None, lli: Some(0) },
            Observation { code: ObsCode::from_str("C7I").unwrap(), value: 3.8e7, lock_time: None, lli: None },
            Observation { code: ObsCode::from_str("L7I").unwrap(), value: 1.5e8, lock_time: None, lli: Some(0) },
            Observation { code: ObsCode::from_str("S1I").unwrap(), value: 45.0, lock_time: None, lli: None },
        ],
    };

    let rov = EpochObs {
        time: GpsTime::new(2200, 100.0),
        satellites: vec![mk_obs(geo), mk_obs(igso)],
    };
    let base = rov.clone();
    let const_sats = vec![(geo, pos_geo), (igso, pos_igso)];
    let bds_id = Constellation::Beidou as u8;

    let cands = ref_sat::filter_reference_candidates(&const_sats, &rov, &base, bds_id);
    assert_eq!(cands.len(), 1, "Only IGSO should qualify for Tier 0");
    assert_eq!(cands[0].0, igso, "Dynamic IGSO satellite must be selected over GEO");
}

#[test]
fn test_reference_transfer_quad_constellation_bands() {
    let mut state = RtkState::new(Vector3::new(1.0, 2.0, 3.0), GpsTime::new(2200, 100.0));
    let const_id = Constellation::Beidou as u8;
    for &band in &[1, 6, 7] {
        let k = DoubleDiffKey { constellation_id: const_id, sat: 8, ref_sat: 6, freq_band: band };
        state.ensure_ambiguity(k, 10.0 * band as f64, 1.0);
    }
    assert_eq!(state.ambiguities.len(), 3);

    for &freq_band in &formation::ALL_RTK_BANDS {
        state.transfer_reference_satellite(const_id, freq_band, 6, 8);
    }

    for &band in &[1, 6, 7] {
        let old_k = DoubleDiffKey { constellation_id: const_id, sat: 8, ref_sat: 6, freq_band: band };
        let new_k = DoubleDiffKey { constellation_id: const_id, sat: 6, ref_sat: 8, freq_band: band };
        assert!(state.get_amb_idx(&old_k).is_none(), "Old key must be removed");
        assert!(state.get_amb_idx(&new_k).is_some(), "New inverted key must be present for band {band}");
        let idx = state.get_amb_idx(&new_k).unwrap() - state.amb_offset();
        assert_eq!(state.ambiguities[idx].1, -10.0 * band as f64, "Inverted ambiguity sign must hold");
    }
}
