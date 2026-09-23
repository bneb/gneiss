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
