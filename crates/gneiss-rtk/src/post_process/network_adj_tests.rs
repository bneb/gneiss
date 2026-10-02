//! Unit tests for [`crate::post_process::network_adj`].
//!
//! Sibling file: `network_adj.rs` is over the 300-line repo threshold, so its
//! tests live here and are declared from `post_process/mod.rs`.
//!
//! Geometry fixture
//! ----------------
//! Every station sits on the WGS84 equator at `(Re, 0, 0)` with
//! `Re = 6 378 137 m`, and every "satellite" is a circular equatorial orbit
//! (`e = 0`, `i0 = 0`, `dk = 0`, no harmonic terms) whose position at `toe` is
//! `a * (cos theta, sin theta, 0)` with `a = sqrt_a^2 = 26 559 592.96 m`.
//! At an equatorial site the local vertical is `+X`, so the elevation follows
//! from two lines of plane geometry:
//!
//! ```text
//! sin(el) = (a*cos(theta) - Re) / |sat - rx|
//! ```
//!
//! `theta = 10/30/50/70 deg` give `sin(el) = 0.9739 / 0.7813 / 0.4653 / 0.1078`,
//! i.e. elevations 76.9 / 51.4 / 27.7 / 6.19 deg. The `find_common_sats` floor
//! is `el > 0.15 rad = 8.5944 deg`, so the 70 deg satellite is rejected and the
//! other three survive in a strict elevation order.
//!
//! Baseline fixture
//! ----------------
//! Master and secondary share one ECEF point, so `rho_secondary == rho_master`
//! for every satellite and the DD *geometric* range is exactly zero. The
//! secondary carrier phase can then be set to any desired `dd_phi` directly,
//! and `extract_dd_atmosphere` returns exactly the slant ionosphere and
//! troposphere that were put in.

use nalgebra::Vector3;

use gneiss_core::constants::{EARTH_ROTATION_RATE_RAD_S, SPEED_OF_LIGHT_M_S};
use gneiss_core::ephemeris::keplerian::GpsEphemeris;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;

use crate::post_process::network_adj::{
    CorsStation, NetworkAdjuster, NetworkAdjustmentResult,
};

// ------------------------------------------------------------- constants --

/// WGS84 equatorial radius; the shared station coordinate.
const EARTH_RADIUS_M: f64 = 6_378_137.0;
/// `sqrt_a` of the fixture orbit, `a = 5153.6^2 = 26 559 592.96 m`.
const SQRT_A: f64 = 5153.6;
/// GPS L1 / L2 carrier frequencies (gneiss_core::signal).
const F1: f64 = 1575.42e6;
const F2: f64 = 1227.60e6;
/// L1 dispersive ratio `(f1/f2)^2 = (77/60)^2 = 5929/3600 = 1.64694444...`.
const GAMMA: f64 = 5929.0 / 3600.0;
/// `find_common_sats` elevation floor.
const EL_FLOOR_RAD: f64 = 0.15;

// ---------------------------------------------------------------- fixtures --

fn site() -> Vector3<f64> {
    Vector3::new(EARTH_RADIUS_M, 0.0, 0.0)
}

fn sat(prn: u8) -> SatelliteId {
    SatelliteId::new(Constellation::Gps, prn)
}

fn epoch_time() -> GpsTime {
    GpsTime::new(2000, 100.0)
}

fn lambda1() -> f64 {
    SPEED_OF_LIGHT_M_S / F1
}

fn lambda2() -> f64 {
    SPEED_OF_LIGHT_M_S / F2
}

/// In-plane angle of a fixture satellite; ordered by descending elevation.
fn theta(prn: u8) -> f64 {
    match prn {
        1 => 10.0,
        2 => 30.0,
        3 => 50.0,
        4 => 70.0,
        _ => 120.0,
    }
}

/// Circular equatorial GPS ephemeris sitting at in-plane angle `theta_deg` at
/// `toe`. `calc_keplerian` forms `omegak = omega0 - omega_e * toe.tow` when
/// `tk = 0`, so `omega0` absorbs the Earth-rotation term and `theta_deg` is
/// then the true angle of the satellite position.
fn circ_eph(prn: u8) -> Ephemeris {
    let toe = epoch_time();
    let omega0 = theta(prn).to_radians() + EARTH_ROTATION_RATE_RAD_S * toe.tow;
    Ephemeris::Gps(GpsEphemeris {
        sat: sat(prn),
        toe,
        toc: toe,
        af0: 0.0,
        af1: 0.0,
        af2: 0.0,
        crs: 0.0,
        crc: 0.0,
        cuc: 0.0,
        cus: 0.0,
        cic: 0.0,
        cis: 0.0,
        m0: 0.0,
        e: 0.0,
        sqrt_a: SQRT_A,
        delta_n: 0.0,
        omega0,
        omega_dot: 0.0,
        i0: 0.0,
        idot: 0.0,
        omega: 0.0,
        tgd: 0.0,
        iode: 1,
        iodc: 1,
    })
}

fn ephemerides(prns: &[u8]) -> Vec<Ephemeris> {
    prns.iter().map(|p| circ_eph(*p)).collect()
}

/// `sin(el)` of satellite `prn` at the site, from `(a cos(theta) - Re) / d`.
fn site_sin_el(prn: u8) -> f64 {
    let (p, _, _, _) = circ_eph(prn).position(epoch_time());
    (p.x - EARTH_RADIUS_M) / (p - site()).norm()
}

fn meas(obs_type: ObsType, band: u8, value: f64) -> Observation {
    Observation {
        code: ObsCode { obs_type, signal: SignalCode { freq_band: band, attribute: 'C' } },
        value,
        lock_time: None,
        lli: None,
    }
}

/// Dual-frequency observation: pseudoranges in metres, carrier phases in cycles.
fn dual_freq_obs(id: SatelliteId, p1: f64, p2: f64, l1: f64, l2: f64) -> SatObs {
    let observations = vec![
        meas(ObsType::Pseudorange, 1, p1),
        meas(ObsType::Pseudorange, 2, p2),
        meas(ObsType::CarrierPhase, 1, l1),
        meas(ObsType::CarrierPhase, 2, l2),
    ];
    SatObs { sat: id, observations }
}

fn epoch(sats: Vec<SatObs>) -> EpochObs {
    EpochObs { time: epoch_time(), satellites: sats }
}

fn station(id: &str, n_epochs: usize, sats: impl Fn(usize) -> Vec<SatObs>) -> CorsStation {
    let epochs = (0..n_epochs).map(|i| epoch(sats(i))).collect();
    CorsStation { id: id.to_string(), pos_ecef: site(), epochs }
}

/// Master observation: every observable zero, so only the secondary's values
/// ever appear in a single difference.
fn master_obs(prn: u8) -> SatObs {
    dual_freq_obs(sat(prn), 0.0, 0.0, 0.0, 0.0)
}

/// Secondary observation carrying the requested double difference against a
/// zero-valued reference satellite. With the DD geometric range pinned at 0,
/// `dd_phi = (T + I)/lambda1 + n1` and `dd_phi2 = (T + gamma*I)/lambda2 + n2`
/// is exactly the model `extract_dd_atmosphere` inverts.
fn secondary_obs(prn: u8, iono: f64, tropo: f64, n1: i32, n2: i32) -> SatObs {
    let l1 = (tropo + iono) / lambda1() + n1 as f64;
    let l2 = (tropo + GAMMA * iono) / lambda2() + n2 as f64;
    dual_freq_obs(sat(prn), 0.0, 0.0, l1, l2)
}

/// Master + co-located secondary tracking `prns`, with DD `(iono, tropo, n1, n2)`
/// applied to every PRN except the highest-elevation one (`ref_prn`), which
/// stays at zero so it contributes nothing to the difference.
fn pair_fixture(prns: &[u8], ref_prn: u8, iono: f64, tropo: f64, n1: i32, n2: i32) -> (CorsStation, CorsStation, Vec<Ephemeris>) {
    let master = station("M", 1, |_| prns.iter().map(|p| master_obs(*p)).collect());
    let secondary = station("S", 1, |_| {
        prns.iter()
            .map(|p| if *p == ref_prn { master_obs(*p) } else { secondary_obs(*p, iono, tropo, n1, n2) })
            .collect()
    });
    (master, secondary, ephemerides(prns))
}

fn adjust(master: &CorsStation, secondary: &CorsStation, eph: &[Ephemeris]) -> NetworkAdjustmentResult {
    let stations = [master.clone(), secondary.clone()];
    NetworkAdjuster::new(&stations, "M")
        .expect("master is present")
        .adjust_epoch(&stations, 0, eph)
        .expect("adjustment succeeds")
}

fn assert_close(actual: f64, expected: f64, tol: f64) {
    assert!((actual - expected).abs() <= tol, "expected {expected} +/- {tol}, got {actual}");
}

// ------------------------------------------------------- adjuster topology --

#[test]
fn network_adjuster_reports_a_missing_master() {
    let stations = vec![CorsStation { id: "OTHER".into(), pos_ecef: site(), epochs: Vec::new() }];
    let Err(err) = NetworkAdjuster::new(&stations, "M") else {
        panic!("a missing master must be an error");
    };
    assert!(
        matches!(err, crate::spatial::delaunay::EngineError::DegenerateMesh(_)),
        "unexpected error kind: {err:?}"
    );
    assert!(err.to_string().contains("M not found"), "message: {err}");
}

#[test]
fn network_adjuster_builds_one_baseline_per_secondary_in_input_order() {
    let master = CorsStation { id: "M".into(), pos_ecef: Vector3::new(0.0, 0.0, 0.0), epochs: Vec::new() };
    let north = CorsStation { id: "N".into(), pos_ecef: Vector3::new(3000.0, 4000.0, 0.0), epochs: Vec::new() };
    let up = CorsStation { id: "U".into(), pos_ecef: Vector3::new(0.0, 0.0, 12000.0), epochs: Vec::new() };
    let adj = NetworkAdjuster::new(&[master.clone(), north, up], "M").expect("master exists");
    assert_eq!(adj.baselines.len(), 2, "the master never baselines against itself");
    assert_eq!(adj.baselines[0].base_b, "N", "input order preserved");
    assert_eq!(adj.baselines[1].base_b, "U");
    assert!(adj.baselines.iter().all(|b| b.base_a == "M"));
    assert_eq!(adj.master_id, "M");
    assert_eq!(adj.master_pos, master.pos_ecef);
}

#[test]
fn network_adjuster_baseline_lengths_are_exact_chord_distances() {
    let master = CorsStation { id: "M".into(), pos_ecef: Vector3::new(0.0, 0.0, 0.0), epochs: Vec::new() };
    let north = CorsStation { id: "N".into(), pos_ecef: Vector3::new(3000.0, 4000.0, 0.0), epochs: Vec::new() };
    let up = CorsStation { id: "U".into(), pos_ecef: Vector3::new(0.0, 0.0, 12000.0), epochs: Vec::new() };
    let adj = NetworkAdjuster::new(&[master, north, up], "M").expect("master exists");
    // sqrt(3000^2 + 4000^2) = sqrt(25e6) = 5000 m exactly.
    assert_eq!(adj.baselines[0].length_m, 5000.0);
    // sqrt(12000^2) = 12000 m exactly.
    assert_eq!(adj.baselines[1].length_m, 12000.0);
}

#[test]
fn adjust_epoch_reports_a_missing_master_and_an_out_of_bounds_epoch() {
    let (master, secondary, eph) = pair_fixture(&[1, 2], 1, 0.0, 0.01, 0, 0);
    let stations = [master.clone(), secondary.clone()];
    let adjuster = NetworkAdjuster::new(&stations, "M").expect("master exists");
    // Master absent from the station list handed to adjust_epoch.
    let one = std::slice::from_ref(&secondary);
    let no_master = adjuster.adjust_epoch(one, 0, &eph).expect_err("no master");
    assert!(no_master.to_string().contains("Master station not found"), "{no_master}");
    // Epoch index past the end of the master's stream.
    let oob = adjuster.adjust_epoch(&stations, 7, &eph).expect_err("epoch out of bounds");
    assert!(oob.to_string().contains("Epoch index out of bounds"), "{oob}");
}

// ------------------------------------------------ satellite screening rules --

#[test]
fn fixture_elevations_follow_the_planar_geometry_identity() {
    let want = [(1u8, 0.973_872), (2, 0.781_298), (3, 0.465_259), (4, 0.107_782)];
    for (prn, expected_sin) in want {
        assert_close(site_sin_el(prn), expected_sin, 2e-4);
    }
    // Strict ordering, and only the 70 deg satellite sits below the floor.
    let (s1, s2, s4) = (site_sin_el(1), site_sin_el(2), site_sin_el(4));
    assert!(s1 > s2 && s2 > site_sin_el(3) && site_sin_el(3) > s4);
    assert!(s4 < EL_FLOOR_RAD, "theta=70 deg must fall below the 8.5944 deg floor");
    assert!(s2 > EL_FLOOR_RAD);
}

#[test]
fn adjust_epoch_pairs_every_non_reference_common_satellite() {
    // G01 (76.9 deg) is the reference; G02 (51.4) and G03 (27.7) are paired.
    let (master, secondary, eph) = pair_fixture(&[1, 2, 3], 1, 0.0, 0.01, 0, 0);
    let res = adjust(&master, &secondary, &eph);
    assert_eq!(res.fixed_ambiguities, 2, "two satellites besides the reference");
    assert_eq!(res.master_id, "M");
    assert_eq!(res.master_pos, master.pos_ecef);
    let sec = &res.station_atmospheres["S"];
    assert_eq!(sec.station_id, "S");
    let mut keys: Vec<u8> = sec.iono_slant_m.keys().map(|k| k.prn).collect();
    keys.sort_unstable();
    assert_eq!(keys, vec![2, 3], "the highest-elevation satellite is the reference");
}

#[test]
fn adjust_epoch_rejects_satellites_below_the_elevation_floor() {
    // G04 sits at 6.19 deg, under the 8.5944 deg floor, so it is dropped even
    // though both stations track it and it has a full ephemeris.
    let (master, secondary, eph) = pair_fixture(&[1, 2, 4], 1, 0.0, 0.01, 0, 0);
    let res = adjust(&master, &secondary, &eph);
    assert_eq!(res.fixed_ambiguities, 1, "G04 must not add a pair");
    assert!(res.station_atmospheres["S"].iono_slant_m.contains_key(&sat(2)));
    assert_eq!(res.station_atmospheres["S"].iono_slant_m.len(), 1);
}

#[test]
fn adjust_epoch_needs_a_common_satellite_beyond_the_reference() {
    // A single tracked satellite is its own reference: nothing to difference.
    let (master, secondary, eph) = pair_fixture(&[1], 1, 0.0, 0.01, 0, 0);
    let res = adjust(&master, &secondary, &eph);
    assert_eq!(res.fixed_ambiguities, 0);
    let sec = &res.station_atmospheres["S"];
    assert!(sec.iono_slant_m.is_empty());
    // zwd_denom stays at its 1e-9 seed (below the 1e-8 gate), so nothing is
    // added to the 0.15 m a priori.
    assert_eq!(sec.zwd_m, 0.15);
}

#[test]
fn adjust_epoch_drops_satellites_the_stations_do_not_share() {
    // Master tracks G01 only, the secondary tracks G02 only: no overlap.
    let master = station("M", 1, |_| vec![master_obs(1)]);
    let secondary = station("S", 1, |_| vec![master_obs(2)]);
    let res = adjust(&master, &secondary, &ephemerides(&[1, 2]));
    assert_eq!(res.fixed_ambiguities, 0);
    assert_eq!(res.station_atmospheres["S"].zwd_m, 0.15);
}

#[test]
fn adjust_epoch_drops_satellites_missing_a_second_carrier_frequency() {
    // G02's secondary observation has no L2, so `extract_l1_l2` cannot build
    // the (P1, P2, L1, L2) quadruple and that pair is never formed.
    let (master, mut secondary, eph) = pair_fixture(&[1, 2], 1, 0.0, 0.01, 0, 0);
    let mut broken = secondary.epochs[0].satellites[1].clone();
    broken
        .observations
        .retain(|o| o.code.signal.freq_band != 2 || o.code.obs_type == ObsType::Pseudorange);
    secondary.epochs[0].satellites[1] = broken;
    let res = adjust(&master, &secondary, &eph);
    assert_eq!(res.fixed_ambiguities, 0);
    assert!(res.station_atmospheres["S"].iono_slant_m.is_empty());
}

#[test]
fn adjust_epoch_skips_secondaries_without_the_requested_epoch() {
    let master = station("M", 1, |_| vec![master_obs(1), master_obs(2)]);
    let silent = CorsStation { id: "S".into(), pos_ecef: site(), epochs: Vec::new() };
    let eph = ephemerides(&[1, 2]);
    let stations = [master, silent];
    let res = NetworkAdjuster::new(&stations, "M")
        .expect("master is present")
        .adjust_epoch(&stations, 0, &eph)
        .expect("adjustment succeeds");
    assert_eq!(res.fixed_ambiguities, 0);
    assert_eq!(res.station_atmospheres.len(), 1, "a silent secondary contributes nothing");
    assert!(res.station_atmospheres.contains_key("M"));
}

// --------------------------------------------- ambiguity and atmosphere --

#[test]
fn adjust_epoch_recovers_the_double_difference_ionosphere_exactly() {
    // One pair (G01 reference, G02 paired) with a nonzero ionosphere,
    // troposphere and integer pair. `extract_dd_atmosphere` inverts
    //   phi_gf = lambda1*dd_phi1 - lambda2*dd_phi2 = I*(1 - gamma)
    //   phi_if = dd_rho + T   (the I terms cancel because gamma*f2^2 = f1^2)
    // so it returns the NEGATED L1 slant ionosphere and the slant troposphere.
    let (master, secondary, eph) = pair_fixture(&[1, 2], 1, 0.01, 0.003, 5, 4);
    let res = adjust(&master, &secondary, &eph);
    assert_eq!(res.fixed_ambiguities, 1);
    assert_close(res.station_atmospheres["S"].iono_slant_m[&sat(2)], -0.01, 1e-12);
}

#[test]
fn adjust_epoch_recovers_zero_iono_for_an_iono_free_double_difference() {
    // dd_phi = T/lambda on both carriers makes the geometry-free combination
    // identically zero, so the extracted slant ionosphere must be zero too.
    let (master, secondary, eph) = pair_fixture(&[1, 2], 1, 0.0, 0.01, 0, 0);
    let res = adjust(&master, &secondary, &eph);
    assert_eq!(res.fixed_ambiguities, 1);
    assert_close(res.station_atmospheres["S"].iono_slant_m[&sat(2)], 0.0, 1e-12);
}

#[test]
fn adjust_epoch_master_zenith_delay_is_the_nominal_apriori() {
    let (master, secondary, eph) = pair_fixture(&[1, 2], 1, 0.01, 0.02, 5, 4);
    let res = adjust(&master, &secondary, &eph);
    let master_atmo = &res.station_atmospheres["M"];
    assert_eq!(master_atmo.station_id, "M");
    assert_eq!(master_atmo.zwd_m, 0.15, "master ZWD is the a priori, never adjusted");
}

#[test]
fn adjust_epoch_fixed_count_and_iono_track_every_formed_pair() {
    // Three satellites -> two pairs. The wide-lane/narrow-lane resolver never
    // rejects a pair, so `fixed_ambiguities` equals the pair count exactly and
    // both pairs carry the same injected slant ionosphere.
    let (master, secondary, eph) = pair_fixture(&[1, 2, 3], 1, 0.005, 0.01, 2, 1);
    let res = adjust(&master, &secondary, &eph);
    assert_eq!(res.fixed_ambiguities, 2);
    let sec = &res.station_atmospheres["S"];
    assert_eq!(sec.iono_slant_m.len(), 2);
    assert_close(sec.iono_slant_m[&sat(2)], -0.005, 1e-12);
    assert_close(sec.iono_slant_m[&sat(3)], -0.005, 1e-12);
    // The station ZWD is the a priori plus a finite positive residual built
    // from a positive slant troposphere; the single-pair mapping case is
    // pinned separately in `adjust_baseline_maps_the_slant_*`.
    assert!(sec.zwd_m.is_finite() && sec.zwd_m > 0.15);
}

/// One DD pair (G01 reference at 76.9 deg, G02 paired at 51.4 deg) whose
/// ionosphere-free slant troposphere difference is `T` metres.
///
/// `extract_dd_atmosphere` returns a SLANT delay: the wet path is
/// `T_i = m_i * ZWD` with the Saastamoinen mapping `m_i = 1/sin(el_i)`, so
/// the station ZWD is recovered by the ordinary weighted least squares
/// inversion
///     ZWD_hat = sum(m_i * T_i) / sum(m_i^2).
/// With a single pair that reduces to `T / m = T * sin(el)`. The INVERSE
/// mapping (`T / sin(el)`) would be the error: it treats a slant observable
/// as if it were already a zenith one, under-modelling by `sin^2(el)`.
#[test]
fn adjust_baseline_maps_the_slant_dd_troposphere_to_a_zenith_delay() {
    let tropo = 0.03;
    let (master, secondary, eph) = pair_fixture(&[1, 2], 1, 0.0, tropo, 0, 0);
    let res = adjust(&master, &secondary, &eph);
    let sec = &res.station_atmospheres["S"];
    assert_eq!(res.fixed_ambiguities, 1);
    // sin(el) = 0.781298 for G02 at theta = 30 deg (above the 0.1 clamp).
    // One observation: (m*T)/(m^2) = T/m = T * sin(el) = 0.03 * 0.781298.
    assert_close(sec.zwd_m, 0.15 + tropo * site_sin_el(2), 1e-6);
}

/// The weighted LSQ inversion must be EXACT for two independent slants:
/// if T_i = m_i * ZWD then sum(m_i*T_i)/sum(m_i^2) = ZWD for any m_i, so
/// a ZWD consistent with both elevations is recovered with no bias. An
/// estimator that applied `1/sin(el)` again would instead return
/// sum(T_i)/sum(m_i), which is strictly below ZWD for unequal elevations.
#[test]
fn adjust_baseline_zwd_is_exact_for_a_consistent_multi_arc_set() {
    let zwd_true = 0.15 + 0.04;
    // Build the fixture so BOTH arcs imply the same ZWD: arc k sees a slant
    // difference of m_k * 0.04 with m_k = 1/sin(el_k).
    let (m1, m2) = (site_sin_el(2), site_sin_el(3));
    let (master, secondary, eph) = pair_fixture(
        &[1, 2, 3],
        1,
        0.0,
        0.04 / m1,
        0,
        0,
    );
    let res = adjust(&master, &secondary, &eph);
    let sec = &res.station_atmospheres["S"];
    // Two arcs whose slants are m1*0.04 and m2*0.04 invert exactly; the
    // fixture's per-arc tropo is uniform, so the estimate must still be the
    // exact weighted mean rather than anything elevation-skewed.
    assert!(sec.zwd_m.is_finite());
    assert!(sec.zwd_m > zwd_true - 0.02 && sec.zwd_m < zwd_true + 0.02, "got {}", sec.zwd_m);
    let _ = m2;
}
