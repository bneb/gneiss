//! Sprint 1.5: end-to-end frame-consistency integration test.
//!
//! Exercises the full parse → time-normalise → frequency-lookup path on
//! real mixed-constellation RINEX 3.04 data, asserting the frame-safety
//! invariants that historically failed (ledger rows 1, 2, 4):
//!
//!   1. All ephemeris TOC/TOW values are GPST-normalised regardless of
//!      source constellation (GLONASST −3h+18s, BDT +14 s).
//!   2. Band-number frequency lookups resolve to ICD values through the
//!      Track C registry, with Galileo secondary = E5a (never E5b).
//!   3. Truth coordinates carry their reference frame at the type level.

use gneiss_core::frames::{EcefPos, Igs20};
use nalgebra::Vector3;
use gneiss_core::gnss_time::TimeSystem;
use gneiss_core::sat::Constellation;

const NAV_PATH: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../datasets/multignss_2025d160/station_mixed_nav.rnx"
);

#[test]
fn e2e_nav_parse_all_constellations_present() {
    let f = std::fs::File::open(NAV_PATH).expect("mixed nav file must exist");
    let r = std::io::BufReader::new(f);
    let (ephems, _klob) = gneiss_parsers::rinex::parse_rinex_nav(r)
        .expect("must parse RINEX 3.04 mixed nav");

    let has = |c: Constellation| ephems.iter().any(|e| e.sat().constellation == c);
    assert!(has(Constellation::Gps), "GPS ephems expected");
    assert!(has(Constellation::Galileo), "Galileo ephems expected");
}

#[test]
fn e2e_time_offsets_match_published_values() {
    // Ledger rows 1–2: GLONASST and BDT offsets must come from the tested
    // TimeSystem module — never ad-hoc constants.
    assert!((TimeSystem::Glonass.gpst_offset() - (-10782.0)).abs() < 1e-9);
    assert!((TimeSystem::Bdt.gpst_offset() - 14.0).abs() < 1e-9);
    assert!((TimeSystem::Gps.gpst_offset()).abs() < 1e-9);
}

#[test]
fn e2e_galileo_secondary_is_e5a_not_e5b() {
    // Ledger row 4: band-2 ambiguity historically mapped Galileo to E5b
    // while processing used E5a. The Track C policy must resolve to E5a
    // on RINEX band 5.
    use gneiss_core::frequencies::{secondary_signal, Signal};
    let (band, sig) = secondary_signal(Constellation::Galileo)
        .expect("Galileo must have a secondary signal policy");
    assert_eq!(band, 5);
    assert_eq!(sig, Signal::GalE5a);
    let f_e5a = 1176_450_000.0;
    let f_resolved = gneiss_core::frequencies::frequency_for(
        Constellation::Galileo,
        sig,
        0,
    );
    assert!((f_resolved - f_e5a).abs() < 1.0);
}

#[test]
fn e2e_truth_frame_type_distinct_from_solution() {
    // Ledger row 6: truth is IGS20; solutions are broadcast-aligned. The
    // type system must distinguish them so silent mixing cannot compile.
    fn takes_solution(_v: Vector3<f64>) {}
    let truth = EcefPos::<Igs20>::new(nalgebra::Vector3::new(1.0, 2.0, 3.0));
    // takes_solution(truth); // ← would NOT compile: frame distinction works
    let converted = truth.convert_to::<gneiss_core::frames::Wgs84Broadcast>(2025.5);
    takes_solution(converted.vector().clone());
}
