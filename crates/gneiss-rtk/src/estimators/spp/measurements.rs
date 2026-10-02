//! SPP raw observation extraction and satellite state calculation.

use gneiss_core::coords::{Coordinate, Datum, Frame};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_core::time::GpsTime;
use super::{SppConfig, SppMeasurement, LIGHT_SPEED};

pub(crate) fn build_single_measurement(
    sat_obs: &gneiss_core::obs::SatObs,
    ephemerides: &[Ephemeris],
    epoch_time: GpsTime,
) -> Option<SppMeasurement> {
    let eph = ephemerides
        .iter()
        .filter(|e| e.sat() == sat_obs.sat)
        .min_by(|a, b| {
            (a.toe().tow - epoch_time.tow)
                .abs()
                .partial_cmp(&(b.toe().tow - epoch_time.tow).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;

    let (f1, mut f2) = gneiss_core::signal::satellite_frequencies(sat_obs.sat, eph.freq_num());
    if f2 == 0.0 {
        f2 = f1;
    }

    let mut freq_band = 1;
    let p1_opt = sat_obs.get_observable(1);
    let p2_opt = match sat_obs.sat.constellation {
        gneiss_core::sat::Constellation::Galileo => {
            if let Some(obs) = sat_obs.get_observable(7) {
                freq_band = 7;
                Some(obs)
            } else if let Some(obs) = sat_obs.get_observable(5) {
                freq_band = 5;
                Some(obs)
            } else {
                None
            }
        }
        gneiss_core::sat::Constellation::Beidou => {
            if let Some(obs) = sat_obs.get_observable(7) {
                freq_band = 7;
                Some(obs)
            } else if let Some(obs) = sat_obs.get_observable(6) {
                freq_band = 6;
                Some(obs)
            } else {
                None
            }
        }
        _ => {
            if let Some(obs) = sat_obs.get_observable(2) {
                freq_band = 2;
                Some(obs)
            } else {
                None
            }
        }
    };

    let (raw_pr, is_iono_free) = if let (Some(p1), Some(p2)) = (p1_opt, p2_opt) {
        let f1_sq = f1 * f1;
        let f2_sq = f2 * f2;
        ((f1_sq * p1 - f2_sq * p2) / (f1_sq - f2_sq), true)
    } else {
        freq_band = 1;
        (p1_opt?, false)
    };

    Some(SppMeasurement {
        constellation: sat_obs.sat.constellation,
        raw_pr,
        snr: sat_obs.get_snr(1).unwrap_or(45) as f64,
        doppler: sat_obs.get_doppler(1).unwrap_or(0.0),
        time: epoch_time,
        eph: eph.clone(),
        is_iono_free,
        freq_band,
    })
}

pub fn build_measurements(
    epoch: &EpochObs,
    ephemerides: &[Ephemeris],
    _config: &SppConfig,
) -> Vec<SppMeasurement> {
    epoch
        .satellites
        .iter()
        .filter_map(|s| build_single_measurement(s, ephemerides, epoch.time))
        .collect()
}

pub(crate) fn compute_sat_state(m: &SppMeasurement, receiver_cdt: f64) -> (Coordinate, f64) {
    let pr_time = m.raw_pr / LIGHT_SPEED;
    let t_rcv = m.time.tow - (receiver_cdt / LIGHT_SPEED);
    let t_tx_sat = t_rcv - pr_time;
    let t_tx_sat_gps = GpsTime::new(m.time.week, t_tx_sat);

    let (_, _, sat_clk_err_rough, _) = if m.is_iono_free {
        m.eph.position_iono_free(t_tx_sat_gps)
    } else if m.freq_band == 7 {
        m.eph.position_e5b(t_tx_sat_gps)
    } else {
        m.eph.position(t_tx_sat_gps)
    };

    let t_tx_true = t_tx_sat - sat_clk_err_rough;
    let t_tx_true_gps = GpsTime::new(m.time.week, t_tx_true);

    let (sat_pos, _, sat_clk_err, _) = if m.is_iono_free {
        m.eph.position_iono_free(t_tx_true_gps)
    } else if m.freq_band == 7 {
        m.eph.position_e5b(t_tx_true_gps)
    } else {
        m.eph.position(t_tx_true_gps)
    };
    let corrected_pr = m.raw_pr + (sat_clk_err * LIGHT_SPEED);

    (
        Coordinate::new(sat_pos, Datum::WGS84, Frame::ECEF, m.time),
        corrected_pr,
    )
}

#[cfg(test)]
mod tests {
    //! Golden vectors for the dual-frequency observation combination.
    //!
    //! The ionosphere-free (IF) combination used here is the textbook
    //!     P_IF = (f1^2 P1 - f2^2 P2) / (f1^2 - f2^2)
    //! which, for a first-order ionospheric delay I(f) = I1 * f1^2/f^2 applied to
    //! both frequencies, cancels exactly:
    //!     P1 = R + I,  P2 = R + I f1^2/f2^2
    //!   => P_IF = [f1^2 (R + I) - f2^2 (R + I f1^2/f2^2)] / (f1^2 - f2^2)
    //!           = R (f1^2 - f2^2) / (f1^2 - f2^2) = R.
    //! So feeding a known R and I must return exactly R; a swapped numerator
    //! (f1^2 P2 - f2^2 P1) or a flipped sign returns R + 2 I f1^2/(f1^2-f2^2),
    //! which is metres off.

    use super::*;
    use crate::estimators::spp::fixture::*;
    use gneiss_core::obs::{EpochObs, Observation, SatObs};
    use gneiss_core::sat::{Constellation, SatelliteId};

    const R_TRUE: f64 = 20_000_000.0;
    const I_L1: f64 = 10.0;

    fn eph_for(constellation: Constellation, prn: u8, t: GpsTime) -> Ephemeris {
        match constellation {
            Constellation::Galileo => galileo_ephemeris(prn, t),
            Constellation::Beidou => beidou_ephemeris(prn, t),
            _ => gps_ephemeris(prn, t),
        }
    }

    fn measure(constellation: Constellation, prn: u8, obs: Vec<Observation>, t: GpsTime) -> Option<SppMeasurement> {
        let eph = eph_for(constellation, prn, t);
        build_single_measurement(&sat_obs(constellation, prn, obs), &[eph], t)
    }

    fn obs_on(band: u8, attr: char, value: f64) -> Observation {
        pseudorange(band, attr, value)
    }

    /// IF combination over two frequencies f1, f2 must cancel a first-order
    /// ionospheric delay and return the geometric range exactly.
    fn assert_if_cancels_iono(f1: f64, f2: f64, band2: u8, constellation: Constellation) {
        let t = epoch_time();
        let p1 = R_TRUE + I_L1;
        let p2 = R_TRUE + I_L1 * (f1 * f1) / (f2 * f2);
        let m = measure(
            constellation,
            1,
            vec![obs_on(1, 'C', p1), obs_on(band2, 'X', p2)],
            t,
        )
        .expect("dual-frequency observation must produce a measurement");
        assert!(m.is_iono_free, "two frequencies must form an ionosphere-free combination");
        assert!(
            (m.raw_pr - R_TRUE).abs() < 1e-6,
            "IF combination gave {} instead of {R_TRUE}",
            m.raw_pr
        );
        assert_eq!(m.freq_band, band2);
    }

    #[test]
    fn gps_l1_l2_iono_free_combination_cancels_a_first_order_delay() {
        assert_if_cancels_iono(L1_HZ, L2_HZ, 2, Constellation::Gps);
    }

    #[test]
    fn galileo_e1_e5b_iono_free_combination_cancels_a_first_order_delay() {
        assert_if_cancels_iono(L1_HZ, E5B_HZ, 7, Constellation::Galileo);
    }

    #[test]
    fn galileo_e1_e5a_falls_back_to_band_5() {
        // Galileo prefers E5b (RINEX band 7); with only E5a present the
        // constellation branch must select band 5 instead of giving up.
        let t = epoch_time();
        let p1 = R_TRUE + I_L1;
        let p2 = R_TRUE + I_L1 * (L1_HZ * L1_HZ) / (E5B_HZ * E5B_HZ);
        let m = measure(
            Constellation::Galileo,
            1,
            vec![obs_on(1, 'C', p1), obs_on(5, 'X', p2)],
            t,
        )
        .expect("E1 + E5a must produce a measurement");
        assert!(m.is_iono_free);
        assert_eq!(m.freq_band, 5);
        assert!((m.raw_pr - R_TRUE).abs() < 1e-6, "got {}", m.raw_pr);
    }

    #[test]
    fn beidou_b1i_b2i_selects_band_6_via_the_constellation_mapping() {
        // BeiDou band 6 is B2I. `get_observable(6)` matches only an observation
        // already reported on band 6, so this exercises that arm directly.
        let t = epoch_time();
        let p1 = R_TRUE + I_L1;
        let p2 = R_TRUE + I_L1 * (B1I_HZ * B1I_HZ) / (B2I_HZ * B2I_HZ);
        let m = measure(
            Constellation::Beidou,
            1,
            vec![obs_on(1, 'I', p1), obs_on(6, 'I', p2)],
            t,
        )
        .expect("B1I + B2I must produce a measurement");
        assert!(m.is_iono_free);
        assert_eq!(m.freq_band, 6);
        assert!((m.raw_pr - R_TRUE).abs() < 1e-6, "got {}", m.raw_pr);
    }

    #[test]
    fn single_frequency_gps_uses_p1_verbatim_and_is_not_iono_free() {
        let t = epoch_time();
        let p1 = R_TRUE + I_L1;
        let m = measure(Constellation::Gps, 1, vec![obs_on(1, 'C', p1)], t)
            .expect("single-frequency observation must produce a measurement");
        assert!(!m.is_iono_free);
        assert_eq!(m.freq_band, 1, "band resets to 1 when only P1 is available");
        assert!((m.raw_pr - p1).abs() < 1e-12, "raw_pr = {}", m.raw_pr);
        // No carrier observation and no SNR -> documented defaults.
        assert!((m.snr - 45.0).abs() < 1e-12);
        assert!((m.doppler).abs() < 1e-15);
    }

    #[test]
    fn galileo_without_a_second_frequency_degrades_to_e1() {
        // Galileo with only E1 must fall through the E5b/E5a chain to the
        // single-frequency path rather than returning None.
        let t = epoch_time();
        let m = measure(
            Constellation::Galileo,
            3,
            vec![obs_on(1, 'C', R_TRUE)],
            t,
        )
        .expect("E1-only Galileo must produce a measurement");
        assert!(!m.is_iono_free);
        assert_eq!(m.freq_band, 1);
        assert!((m.raw_pr - R_TRUE).abs() < 1e-12);
    }

    #[test]
    fn beidou_b1i_b2_selects_band_7_when_reported_on_band_7() {
        // BeiDou checks band 7 (B2, modern) before band 6 (B2I). A
        // pseudorange already reported on band 7 must therefore win, and the
        // ionosphere-free combination must be formed from B1I + B2.
        let t = epoch_time();
        let p1 = R_TRUE + I_L1;
        let p2 = R_TRUE + I_L1 * (B1I_HZ * B1I_HZ) / (E5B_HZ * E5B_HZ);
        let m = measure(
            Constellation::Beidou,
            6,
            vec![obs_on(1, 'I', p1), obs_on(7, 'X', p2)],
            t,
        )
        .expect("B1I + B2 must produce a measurement");
        assert!(m.is_iono_free);
        assert_eq!(m.freq_band, 7);
        assert!((m.raw_pr - R_TRUE).abs() < 1e-6, "got {}", m.raw_pr);
    }

    #[test]
    fn beidou_without_a_second_frequency_degrades_to_b1i() {
        let t = epoch_time();
        let m = measure(
            Constellation::Beidou,
            4,
            vec![obs_on(1, 'I', R_TRUE)],
            t,
        )
        .expect("B1I-only BeiDou must produce a measurement");
        assert!(!m.is_iono_free);
        assert_eq!(m.freq_band, 1);
    }

    #[test]
    fn observations_without_a_matching_ephemeris_or_without_p1_are_dropped() {
        let t = epoch_time();
        let no_eph = SatObs {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 99 },
            observations: vec![pseudorange(1, 'C', R_TRUE)],
        };
        assert!(build_single_measurement(&no_eph, &[gps_ephemeris(1, t)], t).is_none());

        // Right satellite, but the only observation is carrier phase, so there
        // is no usable pseudorange at all.
        let carrier_only = SatObs {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            observations: vec![Observation {
                code: gneiss_core::obs::ObsCode {
                    obs_type: gneiss_core::obs::ObsType::CarrierPhase,
                    signal: gneiss_core::obs::SignalCode { freq_band: 1, attribute: 'L' },
                },
                value: R_TRUE,
                lock_time: None,
                lli: None,
            }],
        };
        assert!(build_single_measurement(&carrier_only, &[gps_ephemeris(1, t)], t).is_none());
    }

    #[test]
    fn build_measurements_keeps_only_satellites_with_ephemeris_and_pseudorange() {
        let t = epoch_time();
        let eph = gps_ephemeris(1, t);
        let epoch = EpochObs {
            time: t,
            satellites: vec![
                sat_obs(Constellation::Gps, 1, vec![pseudorange(1, 'C', R_TRUE)]),
                // No ephemeris for PRN 2.
                sat_obs(Constellation::Gps, 2, vec![pseudorange(1, 'C', R_TRUE)]),
            ],
        };
        let m = build_measurements(&epoch, std::slice::from_ref(&eph), &open_config());
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].eph.sat(), eph.sat());
        assert!((m[0].raw_pr - R_TRUE).abs() < 1e-12);
    }

    #[test]
    fn satellite_state_uses_the_transmission_epoch_and_adds_the_clock_bias() {
        // compute_sat_state corrects the clock by ADDING sat_clk * c to the raw
        // pseudorange (the broadcast clock correction shortens the measured
        // range), and evaluates the ephemeris at the signal transmit time, one
        // light time before reception.
        let t = epoch_time();
        let eph = gps_ephemeris(1, t);
        let raw_pr = 20_100_000.0;
        let m = SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr,
            snr: 45.0,
            doppler: 0.0,
            time: t,
            eph: eph.clone(),
            is_iono_free: false,
            freq_band: 1,
        };
        let (coord, corrected) = compute_sat_state(&m, 0.0);
        let t_tx = GpsTime::new(t.week, t.tow - raw_pr / LIGHT_SPEED);
        let (sat_pos, _, clk, _) = eph.position(t_tx);
        // The two-pass clock refinement moves the epoch by the rough satellite
        // clock bias (up to ~1e-5 s, a few cm of satellite motion at e = 0.01),
        // so 0.1 m of slack still pins the TRANSMIT epoch.
        assert!((coord.vector - sat_pos).norm() < 0.1, "satellite moved to the wrong epoch");
        assert!((corrected - (raw_pr + clk * LIGHT_SPEED)).abs() < 1.0);
        // A receiver clock offset shifts the transmit epoch earlier by
        // cdt / c = 300 / 2.998e8 = 1.0e-6 s, so the satellite must move by
        // v * 1 us ~ 4 mm -- small, but strictly non-zero.
        let (coord_b, _) = compute_sat_state(&m, 300.0);
        let shift = (coord_b.vector - coord.vector).norm();
        assert!(shift > 1e-3 && shift < 1e-1, "shift = {shift} m (expect ~v * 300/c ~ 4 mm)");
    }
}
