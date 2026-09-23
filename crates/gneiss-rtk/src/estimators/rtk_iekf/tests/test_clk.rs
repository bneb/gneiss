#![allow(clippy::unwrap_used)]

use std::sync::Arc;
use super::*;
use gneiss_parsers::clk_centering::CenteredClock;
use gneiss_parsers::rinex_clk::{ClockRecord, RinexClock};

fn centered(bias_us: f64, spread_us: f64) -> Option<CenteredClock> {
    Some(CenteredClock {
        bias_s: bias_us * 1e-6,
        spread_s: spread_us * 1e-6,
    })
}

#[test]
fn centered_pair_correction_healthy_pair_applies_centered_delta() {
    let (corr, tripped) = centered_pair_correction(centered(10.0, 5.0), centered(-14.0, 6.0));
    assert!(!tripped);
    assert!((corr - SPEED_OF_LIGHT_M_S * 24e-6).abs() < 1e-9);
}

#[test]
fn centered_pair_correction_gates_when_either_side_spread_trips() {
    for (a, b) in [
        (centered(1.0, 150.0), centered(2.0, 5.0)),
        (centered(1.0, 5.0), centered(2.0, 150.0)),
    ] {
        let (corr, tripped) = centered_pair_correction(a, b);
        assert!(tripped, "pathological side must trip the gate");
        assert_eq!(corr, 0.0, "gated correction must be suppressed");
    }
}

#[test]
fn centered_pair_correction_threshold_is_strictly_greater() {
    let ok = centered_pair_correction(centered(1.0, 100.0), centered(2.0, 99.999));
    assert!(!ok.1);
    assert!(ok.0.abs() > 0.0);
    let bad = centered_pair_correction(centered(1.0, 100.001), centered(2.0, 5.0));
    assert!(bad.1 && bad.0 == 0.0);
}

#[test]
fn centered_pair_correction_missing_side_stays_silent_zero() {
    assert_eq!(centered_pair_correction(None, centered(2.0, 5.0)), (0.0, false));
    assert_eq!(centered_pair_correction(centered(1.0, 5.0), None), (0.0, false));
    assert_eq!(centered_pair_correction(None, None), (0.0, false));
}

fn clk_product(biases_us: &[(u8, f64)]) -> Arc<RinexClock> {
    let mut rc = RinexClock::default();
    for (prn, us) in biases_us {
        rc.satellites.insert(
            gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: *prn,
            },
            vec![ClockRecord {
                time: GpsTime::new(2200, 100.0),
                bias: us * 1e-6,
            }],
        );
    }
    Arc::new(rc)
}

fn dd_probe(eng: &GnssRtkIekf, sat: u8, reference: u8) -> f64 {
    let rx = Vector3::zeros();
    eng.formation_clock_corr_m(
        gneiss_core::sat::SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: sat,
        },
        u16::from(reference),
        rx,
        Vector3::new(2.0e7, 0.0, 0.0),
        Vector3::new(2.4e7, 0.0, 0.0),
    )
}

#[test]
fn precise_clock_dd_m_without_product_is_zero() {
    let eng = test_engine(GpsTime::new(2200, 100.0));
    assert!(eng.precise_clocks.is_none());
    assert_eq!(dd_probe(&eng, 5, 9), 0.0);
}

#[test]
fn precise_clock_dd_m_healthy_product_preserves_raw_delta() {
    let mut eng = test_engine(GpsTime::new(2200, 100.0));
    eng.precise_clocks = Some(clk_product(&[
        (27, 500.0),
        (28, 520.0),
        (5, 480.0),
        (10, 510.0),
    ]));
    let corr = dd_probe(&eng, 28, 27);
    let raw = SPEED_OF_LIGHT_M_S * 20e-6;
    assert!((corr - raw).abs() < 1e-9, "corr {corr} vs raw {raw}");
    assert!(!eng.clk_gate_warned.load(std::sync::atomic::Ordering::Relaxed));
}

#[test]
fn precise_clock_dd_m_pathological_product_returns_zero_and_latches() {
    let mut eng = test_engine(GpsTime::new(2200, 100.0));
    eng.precise_clocks = Some(clk_product(&[
        (1, 600.0),
        (2, -600.0),
        (3, 590.0),
        (4, -590.0),
    ]));
    assert_eq!(dd_probe(&eng, 1, 3), 0.0);
    assert!(
        eng.clk_gate_warned.load(std::sync::atomic::Ordering::Relaxed),
        "gate trip must latch the warning flag"
    );

    let mut eng = test_engine(GpsTime::new(2200, 100.0));
    eng.precise_clocks = Some(clk_product(&[(1, 600.0), (2, -600.0)]));
    assert_eq!(dd_probe(&eng, 1, 2), 0.0);
    assert!(!eng.clk_gate_warned.load(std::sync::atomic::Ordering::Relaxed));
}

mod obs_side_clk_tests {
    use super::*;

    #[test]
    fn formation_subtracts_clock_delta_from_measurements() {
        let t = GpsTime::new(2370, 43_200.0);
        let mut eng = GnssRtkIekf::new(Vector3::zeros(), t, 1.0);
        eng.state.iono_enabled = false;
        eng.precise_clocks = Some(std::sync::Arc::new(
            gneiss_parsers::rinex_clk::RinexClock::parse(&synthetic_clk_content()),
        ));
        let corr = eng.formation_clock_corr_m(
            gneiss_core::sat::SatelliteId {
                constellation: gneiss_core::sat::Constellation::Gps,
                prn: 1,
            },
            2,
            Vector3::zeros(),
            Vector3::new(2.0e7, 0.0, 0.0),
            Vector3::new(2.1e7, 0.0, 0.0),
        );
        assert!(
            (corr.abs() - 5_995.8).abs() < 10.0,
            "correction {corr} m != expected c·(−20 µs) = −5995.8 m"
        );
    }

    fn synthetic_clk_content() -> String {
        "     3.00           C                                       RINEX VERSION / TYPE\n\
         2    AS    AR                                          # / TYPES OF DATA\n\
         AS G01  2025  6  8 12  0  0.000000  1   -0.000110000000E+00\n\
         AS G02  2025  6  8 12  0  0.000000  1   -0.000090000000E+00\n\
         AS G03  2025  6  8 12  0  0.000000  1   -0.000090000000E+00\n".to_string()
    }
}
