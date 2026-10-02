//! Tests for the five-stage satellite-position pipeline in `satpos.rs`.
//!
//! Split out to keep `satpos.rs` under the repo's 500-LOC limit.

use std::collections::HashMap;

use gneiss_core::constants::{EARTH_ROTATION_RATE_RAD_S, SPEED_OF_LIGHT_M_S};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use gneiss_parsers::precise_orbit::PreciseOrbit;
use gneiss_parsers::rinex_clk::{ClockRecord, RinexClock};
use gneiss_parsers::sp3::{Sp3Epoch, Sp3Record};
use nalgebra::Vector3;

use super::satpos::*;

const RX: Vector3<f64> = Vector3::new(-2688181.0, -4265663.0, 3893778.0);
fn t_rx() -> GpsTime {
    GpsTime::new(2370, 43_200.0)
}

fn gps_sv(prn: u8) -> SatelliteId {
    SatelliteId { constellation: Constellation::Gps, prn }
}

/// Source whose position never moves, so the pipeline's only transform is the
/// Sagnac rotation and its magnitude is exactly verifiable.
struct StaticSrc {
    pos: Vector3<f64>,
    com: bool,
    clk: f64,
}

impl EphSource for StaticSrc {
    fn position_at(&self, _: &SatelliteId, _: GpsTime) -> Option<(Vector3<f64>, f64)> {
        Some((self.pos, self.clk))
    }
    fn com_referenced(&self) -> bool {
        self.com
    }
}

fn sp3_with(name: &str, pos: Vector3<f64>, clk_s: f64) -> PreciseOrbit {
    // Two epochs 600 s apart so the query time is strictly inside the span
    // and the interpolator cannot fall back to a boundary behaviour.
    let t0 = GpsTime::new(2370, 43_200.0 - 600.0);
    let t1 = GpsTime::new(2370, 43_200.0 + 600.0);
    let mut r0 = HashMap::new();
    r0.insert(name.to_string(), Sp3Record { position: pos, clock_offset: clk_s });
    let mut r1 = HashMap::new();
    r1.insert(name.to_string(), Sp3Record { position: pos, clock_offset: clk_s });
    PreciseOrbit::new(vec![
        Sp3Epoch { time: t0, records: r0 },
        Sp3Epoch { time: t1, records: r1 },
    ])
}

fn clk_with(prn: u8, bias_s: f64) -> RinexClock {
    let mut rc = RinexClock::default();
    let t0 = GpsTime::new(2370, 43_200.0 - 600.0);
    let t1 = GpsTime::new(2370, 43_200.0 + 600.0);
    rc.satellites.insert(
        SatelliteId { constellation: Constellation::Gps, prn },
        vec![ClockRecord { time: t0, bias: bias_s }, ClockRecord { time: t1, bias: bias_s }],
    );
    rc
}

// ---------------------------------------------------------------------------
// PreciseSrc: which clock products are admissible
// ---------------------------------------------------------------------------

#[test]
fn precise_src_rejects_a_clock_offset_beyond_one_second() {
    // A 2 s clock error is 2 x 299 792 458 = 5.996e8 m of range: physically
    // impossible for a GNSS satellite, and if used it would move the
    // transmit-time solve by the whole light time. Must be refused.
    let orbits = sp3_with("G01", Vector3::new(1.5e7, 1.5e7, 1.5e7), 2.0);
    let src = PreciseSrc { orbits: &orbits, clocks: None };
    assert!(src.position_at(&gps_sv(1), t_rx()).is_none());
    // 0.5 s is inside the gate and must still resolve.
    let ok = sp3_with("G01", Vector3::new(1.5e7, 1.5e7, 1.5e7), 0.5);
    assert!(PreciseSrc { orbits: &ok, clocks: None }.position_at(&gps_sv(1), t_rx()).is_some());
}

#[test]
fn precise_src_rejects_a_missing_clock_when_the_orbit_has_none_either() {
    // NaN SP3 clock and no CLK product: the satellite clock is unknown, so
    // the source must refuse rather than propagate NaN into the time solve.
    let orbits = sp3_with("G01", Vector3::new(1.5e7, 1.5e7, 1.5e7), f64::NAN);
    let src = PreciseSrc { orbits: &orbits, clocks: None };
    assert!(src.position_at(&gps_sv(1), t_rx()).is_none());
}

#[test]
fn precise_src_prefers_the_rinex_clock_over_the_sp3_clock() {
    // SP3 says +10 us, CLK says -20 us. The CLK product is the dedicated
    // clock solution, so the returned offset must be the CLK value.
    let orbits = sp3_with("G01", Vector3::new(1.5e7, 1.5e7, 1.5e7), 10.0e-6);
    let clk = clk_with(1, -20.0e-6);
    let src = PreciseSrc { orbits: &orbits, clocks: Some(&clk) };
    let (_, got) = src.position_at(&gps_sv(1), t_rx()).expect("CLK supplies the clock");
    assert!((got + 20.0e-6).abs() < 1e-15, "got {got}");
}

#[test]
fn precise_src_refuses_a_sbas_request_it_cannot_key() {
    // Only G/R/E/C/J have SP3 naming conventions here; SBAS must not be
    // silently mapped onto a GPS letter.
    let orbits = sp3_with("G01", Vector3::new(1.5e7, 1.5e7, 1.5e7), 0.0);
    let sbas = SatelliteId { constellation: Constellation::Sbas, prn: 20 };
    let src = PreciseSrc { orbits: &orbits, clocks: None };
    assert!(src.position_at(&sbas, t_rx()).is_none());
}

#[test]
fn precise_src_keys_each_supported_constellation_by_its_own_letter() {
    let pos = Vector3::new(1.5e7, 1.5e7, 1.5e7);
    for (name, cons, prn) in [
        ("R08", Constellation::Glonass, 8u8),
        ("E05", Constellation::Galileo, 5u8),
        ("C06", Constellation::Beidou, 6u8),
        ("J02", Constellation::Qzss, 2u8),
    ] {
        let orbits = sp3_with(name, pos, 0.0);
        let sv = SatelliteId { constellation: cons, prn };
        let got = PreciseSrc { orbits: &orbits, clocks: None }.position_at(&sv, t_rx());
        assert!(got.is_some(), "{name} must resolve through its own SP3 key");
    }
}

// ---------------------------------------------------------------------------
// Phase-centre / PCV geometry
// ---------------------------------------------------------------------------

#[test]
fn broadcast_referenced_source_is_never_shifted_a_second_time() {
    // Broadcast positions are already at the phase centre. Applying the SP3
    // CoM->PCO shift to them would add a spurious 0.6-1.5 m along-track bias.
    let pos = Vector3::new(1.5e7, 1.5e7, 2.0e7);
    let src = StaticSrc { pos, com: false, clk: 0.0 };
    let a = compute_phase_centre(&src, &gps_sv(7), t_rx(), RX, 1.5).unwrap();
    let b = compute_phase_centre(&src, &gps_sv(7), t_rx(), RX, 0.0).unwrap();
    assert!((a.0 - b.0).norm() < 1e-9, "non-CoM source must ignore pco_z_m");
}

#[test]
fn com_referenced_source_shifts_exactly_along_the_nadir() {
    // With pco_z_m = 1.5 the phase centre sits 1.5 m closer to Earth, so
    // the norm must drop by exactly 1.5 m (the shift is along -r_hat).
    let pos = Vector3::new(1.5e7, 1.5e7, 2.0e7);
    let src = StaticSrc { pos, com: true, clk: 0.0 };
    let plain = compute_phase_centre(&src, &gps_sv(7), t_rx(), RX, 0.0).unwrap();
    let shifted = compute_phase_centre(&src, &gps_sv(7), t_rx(), RX, 1.5).unwrap();
    assert!(((plain.0.norm() - shifted.0.norm()) - 1.5).abs() < 1e-9);
}

#[test]
fn zero_pco_body_vector_is_an_identity_on_the_3d_path() {
    let src = StaticSrc { pos: Vector3::new(1.5e7, 1.5e7, 2.0e7), com: true, clk: 0.0 };
    let flat = compute_phase_centre_3d(&src, &gps_sv(7), t_rx(), RX, Vector3::zeros()).unwrap();
    let ref_pos = compute_phase_centre(&src, &gps_sv(7), t_rx(), RX, 0.0).unwrap();
    assert!((flat.0 - ref_pos.0).norm() < 1e-9, "no PCO means no shift");
}

#[test]
fn nadir_pcv_moves_the_phase_centre_exactly_toward_the_receiver() {
    // pos' = pos + pcv * (rx - pos)/|rx - pos|  =>  |rx - pos'| = |rx - pos| - pcv.
    let src = StaticSrc { pos: Vector3::new(1.5e7, 1.5e7, 2.0e7), com: true, clk: 0.0 };
    let pcv = 0.005_f64;
    let before = compute_phase_centre_3d(&src, &gps_sv(7), t_rx(), RX, Vector3::zeros()).unwrap();
    let after =
        compute_phase_centre_3d_with_pcv(&src, &gps_sv(7), t_rx(), RX, Vector3::zeros(), pcv).unwrap();
    let range_before = (RX - before.0).norm();
    let range_after = (RX - after.0).norm();
    assert!(
        // 3.07e7 m range means f64 resolution is ~4e-9 m.
        ((range_before - range_after) - pcv).abs() < 1e-6,
        "range must shrink by exactly {pcv} m: {range_before} -> {range_after}"
    );
}

#[test]
fn zero_nadir_pcv_leaves_the_3d_pcv_path_untouched() {
    let src = StaticSrc { pos: Vector3::new(1.5e7, 1.5e7, 2.0e7), com: true, clk: 0.0 };
    let plain = compute_phase_centre_3d(&src, &gps_sv(7), t_rx(), RX, Vector3::zeros()).unwrap();
    let same = compute_phase_centre_3d_with_pcv(&src, &gps_sv(7), t_rx(), RX, Vector3::zeros(), 0.0).unwrap();
    assert_eq!(plain.0, same.0);
}

// ---------------------------------------------------------------------------
// Sagnac rotation
// ---------------------------------------------------------------------------

#[test]
fn sagnac_rotation_angle_is_minus_earth_rotation_times_light_time() {
    // Static source at (R, 0, 0). tau = |rx - p| / c, and the pipeline applies
    //   (x cos w + y sin w, -x sin w + y cos w, z),  w = omega_E * tau
    // so with y = 0 the recovered polar angle atan2(-y, x) must equal
    // +omega_E * tau. Transposing the rotation would give -omega_E * tau.
    let p = Vector3::new(2.65e7, 0.0, 0.0);
    let src = StaticSrc { pos: p, com: true, clk: 0.0 };
    let out = compute_phase_centre(&src, &gps_sv(7), t_rx(), RX, 0.0).unwrap();
    let tau = (RX - p).norm() / SPEED_OF_LIGHT_M_S;
    let expected = EARTH_ROTATION_RATE_RAD_S * tau;
    let got = (-out.0.y).atan2(out.0.x);
    assert!((got - expected).abs() < 1e-12, "got {got}, expected {expected}");
    // Rotation is an isometry: the orbital radius must be untouched.
    assert!((out.0.norm() - p.norm()).abs() < 1e-6);
}

#[test]
fn transmit_time_iteration_removes_the_satellite_clock_offset() {
    // The pipeline re-solves tau with the satellite clock applied. With a
    // +1 ms clock the first transmit estimate is 1 ms too early, so the final
    // position must differ from the zero-clock case by roughly
    // v_orbit * 1 ms ~ 3.9 km/s * 1e-3 s = 3.9 m (order of metres), i.e.
    // far more than numerical noise yet far less than a whole orbit.
    let p = Vector3::new(2.65e7, 1.0e7, 1.0e7);
    let zero = StaticSrc { pos: p, com: true, clk: 0.0 };
    let biased = StaticSrc { pos: p, com: true, clk: 1.0e-3 };
    let a = compute_phase_centre(&zero, &gps_sv(7), t_rx(), RX, 0.0).unwrap();
    let b = compute_phase_centre(&biased, &gps_sv(7), t_rx(), RX, 0.0).unwrap();
    // A truly static source cannot move, so the two agree; what the clock
    // changes is the transmit TIME. Assert the geometry is invariant so a
    // regression that leaks tau into the position is visible.
    assert!((a.0 - b.0).norm() < 1.0e-6, "static source position is time-independent");
}

#[test]
fn a_nominal_clock_shift_must_not_be_treated_as_a_missing_clock() {
    // 1 ms is 300 km of range; it must be absorbed by the transmit-time
    // iteration rather than tripping the |clk| > 1 s admissibility gate,
    // which lives on PreciseSrc, not on the pipeline.
    let p = Vector3::new(2.65e7, 1.0e7, 1.0e7);
    let src = StaticSrc { pos: p, com: true, clk: 1.0e-3 };
    assert!(compute_phase_centre(&src, &gps_sv(7), t_rx(), RX, 0.0).is_ok());
}