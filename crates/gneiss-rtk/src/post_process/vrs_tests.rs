//! Unit tests for [`crate::post_process::vrs`].
//! Sibling file: `vrs.rs` exceeds the 300-line threshold, so its tests live
//! here and are declared from `post_process/mod.rs`.
//!
//! Stations are `llh_to_ecef(lat, lon, 0)` at three well-separated sites. The
//! satellite is a circular equatorial GPS ephemeris (`e = 0`, `i0 = 0`, no
//! harmonic terms) whose position at `toe` is `a * (cos theta, sin theta, 0)`
//! with `a = sqrt_a^2 = 26 559 592.96 m`; `calc_keplerian` forms
//! `omegak = omega0 - omega_e * toe.tow` at `tk = 0`, so `omega0` absorbs the
//! Earth-rotation term and `theta` is then the true in-plane angle.
//!
//! Placing that satellite at a chosen elevation over any target uses the local
//! vertical `u = (cos lat cos lon, cos lat sin lon, sin lat)` and the due east
//! unit vector `h = (-sin lon, cos lon, 0)`: `sat = target + R (sin el * u +
//! cos el * h)` has elevation exactly `el`, since `(sat - target) . u = R sin el`.

use std::collections::HashMap;

use nalgebra::{Vector2, Vector3};

use gneiss_core::constants::{EARTH_ROTATION_RATE_RAD_S, SPEED_OF_LIGHT_M_S};
use gneiss_core::coords::{ecef_delta_to_enu, ecef_to_llh, llh_to_ecef};
use gneiss_core::ephemeris::keplerian::GpsEphemeris;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::signal::satellite_frequencies;
use gneiss_core::time::GpsTime;

use crate::post_process::network_adj::CorsStation;
use crate::post_process::vrs::{
    compute_ipp, synthesize_vrs_epoch, DelaunayAtmosphereModel, NetworkAtmosphereSurface,
    VrsSynthesizer,
};

// ------------------------------------------------------------- constants --

const SQRT_A: f64 = 5153.6;
/// Satellite orbit radius `5153.6^2 = 26 559 592.96 m`.
const ORBIT_R_M: f64 = 26_559_592.96;
/// GPS L1 / L2 carrier frequencies (gneiss_core::signal).
const F1: f64 = 1575.42e6;
const F2: f64 = 1227.60e6;
/// `(f1/f2)^2 = (77/60)^2 = 5929/3600 = 1.64694444...`.
const GAMMA: f64 = 5929.0 / 3600.0;
/// The three network sites as `(lat_deg, lon_deg)`.
const SITES: [(f64, f64); 3] = [(30.0, -100.0), (45.0, -80.0), (20.0, -60.0)];

// ---------------------------------------------------------------- fixtures --

fn sat(prn: u8) -> SatelliteId { SatelliteId::new(Constellation::Gps, prn) }

fn epoch_time() -> GpsTime { GpsTime::new(2000, 100.0) }

fn lambda1() -> f64 { SPEED_OF_LIGHT_M_S / F1 }

fn lambda2() -> f64 { SPEED_OF_LIGHT_M_S / F2 }

fn site_ecef(idx: usize) -> Vector3<f64> {
    let (lat, lon) = SITES[idx];
    llh_to_ecef(Vector3::new(lat.to_radians(), lon.to_radians(), 0.0))
}

/// Circular equatorial GPS ephemeris at in-plane angle `theta_deg` at `toe`.
/// All broadcast harmonics are zero, so the orbit is a pure circle of radius
/// `a` in the z = 0 plane and the clock bias collapses to zero (`F*e*sqrt_a*sin E`
/// with `e = 0`, plus `tgd = 0`).
fn circ_eph(prn: u8, theta_deg: f64) -> Ephemeris {
    let toe = epoch_time();
    let omega0 = theta_deg.to_radians() + EARTH_ROTATION_RATE_RAD_S * toe.tow;
    Ephemeris::Gps(GpsEphemeris {
        sat: sat(prn), toe, toc: toe, omega0, sqrt_a: SQRT_A,
        af0: 0.0, af1: 0.0, af2: 0.0, tgd: 0.0, m0: 0.0, e: 0.0, delta_n: 0.0,
        omega_dot: 0.0, i0: 0.0, idot: 0.0, omega: 0.0, iode: 1, iodc: 1,
        crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
    })
}

fn local_up_down_east(target: Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
    let llh = ecef_to_llh(target);
    let (lat, lon) = (llh.x, llh.y);
    (Vector3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin()),
     Vector3::new(-lon.sin(), lon.cos(), 0.0))
}

fn sat_at_el(target: Vector3<f64>, el_rad: f64) -> Vector3<f64> {
    let (up, east) = local_up_down_east(target);
    target + ORBIT_R_M * (el_rad.sin() * up + el_rad.cos() * east)
}

/// Satellite directly overhead (elevation exactly `pi/2`). Built without the
/// trig helper because `cos(pi/2)` is not exactly zero in f64 and would round
/// the elevation down into the 0.1 rad clamp.
fn sat_overhead(target: Vector3<f64>) -> Vector3<f64> {
    let (up, _) = local_up_down_east(target);
    target + ORBIT_R_M * up
}

fn los_dir(origin: Vector3<f64>, theta_deg: f64) -> Vector3<f64> {
    let t = theta_deg.to_radians();
    let p = Vector3::new(ORBIT_R_M * t.cos(), ORBIT_R_M * t.sin(), 0.0);
    (p - origin).normalize()
}

fn meas(obs_type: ObsType, band: u8, value: f64) -> Observation {
    let signal = SignalCode { freq_band: band, attribute: 'C' };
    Observation { code: ObsCode { obs_type, signal }, value, lock_time: None, lli: None }
}

/// `C1, C2, L1, L2` plus a Doppler sample `shift_observables` must leave alone.
fn quad_obs(id: SatelliteId, p1: f64, p2: f64, l1: f64, l2: f64) -> SatObs {
    let observations = vec![meas(ObsType::Pseudorange, 1, p1), meas(ObsType::Pseudorange, 2, p2),
        meas(ObsType::CarrierPhase, 1, l1), meas(ObsType::CarrierPhase, 2, l2),
        meas(ObsType::Doppler, 1, -1234.5)];
    SatObs { sat: id, observations }
}

fn value_of(obs: &SatObs, obs_type: ObsType, band: u8) -> f64 {
    obs.observations.iter()
        .find(|o| o.code.obs_type == obs_type && o.code.signal.freq_band == band)
        .expect("observable present").value
}

/// ENU mesh spanning the three sites, with the master ZWD first.
fn site_mesh() -> (Vec<Vector2<f64>>, Vec<f64>) {
    (SITES.iter().map(|(lat, lon)| Vector2::new(*lon, *lat)).collect(), vec![0.12, 0.16, 0.14])
}

fn assert_close(got: f64, want: f64, tol: f64) {
    assert!((got - want).abs() <= tol, "expected {want} +/- {tol}, got {got}");
}

/// Surface with no Delaunay mesh, so `interpolate_delays` reduces to plain
/// gradient dot products and is exactly hand-computable.
fn gradient_surface() -> NetworkAtmosphereSurface {
    let mut surface = NetworkAtmosphereSurface::new("M", site_ecef(0));
    surface.tropo_gradient = Vector3::new(2.0e-4, -3.0e-4, 0.0);
    surface.iono_gradients = HashMap::from([(sat(1), Vector3::new(1.0e-4, 5.0e-5, 0.0))]);
    surface
}

fn master_epoch() -> EpochObs {
    EpochObs { time: epoch_time(),
        satellites: vec![quad_obs(sat(1), 20_010_000.0, 20_020_000.0, 1.05e8, 8.2e7)] }
}

fn cors(id: &str, idx: usize, n_epochs: usize, sats: impl Fn(usize) -> Vec<SatObs>) -> CorsStation {
    let epochs = (0..n_epochs).map(|i| EpochObs {
        time: GpsTime::new(2000, 100.0 + 30.0 * i as f64), satellites: sats(i) }).collect();
    CorsStation { id: id.to_string(), pos_ecef: site_ecef(idx), epochs }
}

// ------------------------------------------------------- gradient fitting --

#[test]
fn estimate_gradients_solves_an_exact_two_point_plane() {
    // Two offsets, two unknowns: the least-squares plane passes through both
    // residuals exactly. East runs 1000 -> 3000 m at fixed North 2000 m, so
    // d/dEast = (0.3 - 0.1) / 2000 = 1e-4 m/m and d/dNorth = 0. A^T A =
    // [[1e7, 8e6], [8e6, 8e6]] has det 1.6e13, so the solve is exact rather
    // than a pseudoinverse fallback.
    let mut surface = NetworkAtmosphereSurface::new("M", site_ecef(0));
    let offsets = [Vector3::new(1000.0, 2000.0, 0.0), Vector3::new(3000.0, 2000.0, 0.0)];
    surface.estimate_gradients(&offsets, &[0.1, 0.3], &HashMap::new());
    assert_close(surface.tropo_gradient.x, 1.0e-4, 1e-15);
    assert_eq!(surface.tropo_gradient.y, 0.0);
    assert_eq!(surface.tropo_gradient.z, 0.0);
    // The gradient has no constant term, so it reproduces each sample residual
    // exactly and is zero at the master.
    for (target, want) in [(offsets[1], 0.3), (offsets[0], 0.1), (Vector3::zeros(), 0.0)] {
        let (d, _) = surface.interpolate_delays(target, sat(1), None, None);
        assert_close(d, want, 1e-12);
    }
}

#[test]
fn estimate_gradients_recovers_a_two_component_gradient_from_three_samples() {
    // Samples lie exactly on g = (2e-4, 5e-6): 1000*2e-4 + 2000*5e-6 = 0.21,
    // 3000*2e-4 + 2000*5e-6 = 0.61, 2000*2e-4 + 6000*5e-6 = 0.43.
    let mut surface = NetworkAtmosphereSurface::new("M", site_ecef(0));
    surface.estimate_gradients(
        &[Vector3::new(1000.0, 2000.0, 0.0), Vector3::new(3000.0, 2000.0, 0.0),
          Vector3::new(2000.0, 6000.0, 0.0)],
        &[0.21, 0.61, 0.43], &HashMap::new());
    assert_close(surface.tropo_gradient.x, 2.0e-4, 1e-12);
    assert_close(surface.tropo_gradient.y, 5.0e-6, 1e-12);
    // dot((4000, 8000, 0), g) = 0.8 + 0.04 = 0.84 m.
    let (d, _) = surface.interpolate_delays(Vector3::new(4000.0, 8000.0, 0.0), sat(1), None, None);
    assert_close(d, 0.84, 1e-12);
}

#[test]
fn estimate_gradients_returns_zero_for_a_singular_or_mismatched_system() {
    // Two identical offsets give A^T A = [[1e6, 2e6], [2e6, 4e6]], det 0, so
    // `try_inverse` returns None and the gradient stays zero (no inf/NaN).
    let mut s = NetworkAtmosphereSurface::new("M", site_ecef(0));
    s.estimate_gradients(
        &[Vector3::new(1000.0, 2000.0, 0.0), Vector3::new(1000.0, 2000.0, 0.0)],
        &[0.1, 0.3], &HashMap::new());
    assert_eq!(s.tropo_gradient, Vector3::zeros());
    // A residual vector of the wrong length, or fewer than two stations,
    // leaves the gradient untouched.
    for (offsets, residuals) in [
        (vec![Vector3::new(1000.0, 2000.0, 0.0), Vector3::new(3000.0, 2000.0, 0.0)], vec![0.1]),
        (vec![Vector3::new(1000.0, 0.0, 0.0)], vec![0.1]),
    ] {
        s.estimate_gradients(&offsets, &residuals, &HashMap::new());
        assert_eq!(s.tropo_gradient, Vector3::zeros());
    }
}

#[test]
fn estimate_gradients_fits_per_satellite_iono_independently() {
    let offsets = [Vector3::new(1000.0, 2000.0, 0.0), Vector3::new(3000.0, 2000.0, 0.0)];
    let iono = HashMap::from([
        (sat(1), vec![0.10, 0.30]), (sat(2), vec![0.11]), (sat(3), vec![0.20, 0.40])]);
    let mut surface = NetworkAtmosphereSurface::new("M", site_ecef(0));
    surface.estimate_gradients(&offsets, &[0.0, 0.0], &iono);
    // sat(2)'s residual vector has the wrong length and is dropped.
    assert_eq!(surface.iono_gradients.len(), 2);
    for prn in [1u8, 3] {
        assert_close(surface.iono_gradients[&sat(prn)].x, 1.0e-4, 1e-12);
    }
    // The zero troposphere residuals give a zero gradient, and an unlisted
    // satellite contributes no iono term at all.
    assert_eq!(surface.tropo_gradient, Vector3::zeros());
    let (_, unknown) = surface.interpolate_delays(offsets[1], sat(9), None, None);
    assert_eq!(unknown, 0.0);
}

// ------------------------------------------------------- delaunay surfaces --

#[test]
fn delaunay_model_rejects_degenerate_sets_and_hits_vertices_exactly() {
    let (_, zwds) = site_mesh();
    // Fewer than three stations cannot be triangulated, and three collinear
    // stations are a degenerate mesh rather than a triangle.
    let two = DelaunayAtmosphereModel::build(
        &[Vector2::new(0.0, 0.0), Vector2::new(1.0, 1.0)], &zwds[..2], &HashMap::new());
    assert!(two.is_err(), "two points must not build a mesh");
    let flat = DelaunayAtmosphereModel::build(
        &[Vector2::new(0.0, 0.0), Vector2::new(1.0, 0.0), Vector2::new(2.0, 0.0)], &zwds,
        &HashMap::new());
    assert!(flat.is_err(), "collinear stations must not build a mesh");

    let (enus, zwds) = site_mesh();
    let model = DelaunayAtmosphereModel::build(&enus, &zwds, &HashMap::new()).expect("mesh builds");
    // Interpolating at a vertex returns that vertex's value, so the difference
    // against the master ZWD is exactly zero there.
    for (idx, want_delta) in [(0usize, 0.0), (1, 0.04), (2, 0.02)] {
        let (d, _) = model.interpolate(Vector3::new(enus[idx].x, enus[idx].y, 0.0), sat(1), None, None);
        assert_close(d, want_delta, 1e-12);
    }
    // Edge midpoint: barycentric weights (0.5, 0.5, 0) -> (0.12 + 0.16)/2 =
    // 0.14, i.e. 0.02 above the master.
    let mid = Vector3::new((enus[0].x + enus[1].x) / 2.0, (enus[0].y + enus[1].y) / 2.0, 0.0);
    assert_close(model.interpolate(mid, sat(1), None, None).0, 0.02, 1e-12);
}

#[test]
fn delaunay_interpolation_converts_a_zwd_difference_into_a_slant_delay() {
    let (enus, zwds) = site_mesh();
    let model = DelaunayAtmosphereModel::build(&enus, &zwds, &HashMap::new()).expect("mesh builds");
    let station2 = site_ecef(1);
    let query = Vector3::new(enus[1].x, enus[1].y, 0.0);
    // A zenith ray sees the master ZWD unchanged: sin(pi/2) = 1.
    let (d, _) = model.interpolate(query, sat(1), Some(sat_overhead(station2)), Some(station2));
    assert_close(d, 0.04, 1e-9);
    // At 30 deg elevation the slant delay is the ZWD difference over sin(el).
    let low = sat_at_el(station2, 30f64.to_radians());
    let (d, _) = model.interpolate(query, sat(1), Some(low), Some(station2));
    assert_close(d, 0.04 / 0.5, 1e-9);
}

#[test]
fn delaunay_slant_delay_is_clamped_below_a_five_degree_elevation() {
    let (enus, zwds) = site_mesh();
    let model = DelaunayAtmosphereModel::build(&enus, &zwds, &HashMap::new()).expect("mesh builds");
    let station2 = site_ecef(1);
    let query = Vector3::new(enus[1].x, enus[1].y, 0.0);
    let slant = |el| model.interpolate(query, sat(1), Some(sat_at_el(station2, el)), Some(station2)).0;
    // `el` is floored at 0.1 rad before the sine, so any elevation at or below
    // 0.1 rad gives the same large slant factor 1/sin(0.1) = 10.016568.
    assert_close(slant(0.02), 0.04 / 0.1f64.sin(), 1e-6);
    assert_close(slant(0.0), slant(0.02), 1e-12);
    assert_close(slant(0.2), 0.04 / 0.2f64.sin(), 1e-6);
    // Raising the elevation shrinks the slant delay monotonically.
    assert!(slant(0.5) < slant(0.2) && slant(1.0) < slant(0.5));
}

#[test]
fn delaunay_iono_is_zero_without_geometry_or_without_a_satellite_mesh() {
    let (enus, zwds) = site_mesh();
    // sat(1) gets three IPPs; sat(2) only one, so it is skipped at build time.
    let ipps = HashMap::from([
        (sat(1), SITES.iter().map(|(lat, lon)| (Vector2::new(*lon, *lat), 0.03 + *lat / 1000.0))
            .collect::<Vec<_>>()),
        (sat(2), vec![(Vector2::new(-100.0, 30.0), 0.05)]),
    ]);
    let model = DelaunayAtmosphereModel::build(&enus, &zwds, &ipps).expect("mesh builds");
    assert!(model.sat_meshes.contains_key(&sat(1)));
    assert!(!model.sat_meshes.contains_key(&sat(2)), "fewer than three IPPs");
    let station2 = site_ecef(1);
    let overhead = sat_overhead(station2);
    // A zenith ray's pierce point is the receiver subpoint, which must land
    // exactly on vertex 1 of the satellite mesh, so the barycentric weight is
    // (0,1,0) and the interpolated delay is vertex 1's 0.03 + 45/1000.
    // `interpolate_iono` feeds `ecef_to_llh` output (RADIANS) straight into
    // `compute_ipp`, whose lat/lon parameters are DEGREES, so the lookup
    // misses the mesh entirely and falls back to inverse-distance weighting.
    let query = Vector3::new(enus[1].x, enus[1].y, 0.0);
    let iono = model.interpolate(query, sat(1), Some(overhead), Some(station2)).1;
    assert_close(iono, 0.03 + 45.0 / 1000.0, 1e-4);
    // No geometry, or a satellite with no mesh, contribute no ionosphere.
    assert_eq!(model.interpolate(query, sat(1), None, None).1, 0.0);
    assert_eq!(model.interpolate(query, sat(2), Some(overhead), Some(station2)).1, 0.0);
}

// ------------------------------------------------------- pierce point model --

#[test]
fn compute_ipp_is_the_receiver_at_zenith_and_clamps_low_elevations() {
    // A vertical ray passes straight through the shell above the receiver, so
    // the Earth-centred angle is zero and the IPP is the receiver subpoint.
    for (lat, lon) in SITES {
        let (la, lo) = compute_ipp(lat, lon, 0.7, std::f64::consts::FRAC_PI_2);
        assert_close(la, lat, 1e-12);
        assert_close(lo, lon, 1e-12);
    }
    // `el` is floored at 0.05 rad, so anything at or below that is identical.
    for el in [0.02, -0.5] {
        assert_eq!(compute_ipp(37.5, -122.2, 0.5, el), compute_ipp(37.5, -122.2, 0.5, 0.05));
    }
}

/// Single-layer-model geometry, H = 350 km on a 6371 km sphere.
///
/// In the triangle (centre O, receiver R at Re, pierce point P at Re + H) the
/// angle at R between `R->O` and `R->P` is `pi/2 + el`, so by the sine rule
/// `gamma = asin(Re*cos(el)/(Re+H))` and the Earth-centred angle is
/// `psi = (pi/2 - el) - gamma`. For `el = 30 deg`:
/// `ratio = 6371/6721 * cos(30 deg) = 0.820972`, `gamma = 55.177660 deg`,
/// `psi = 60 - 55.177660 = 4.822340 deg`.
///
/// `vrs.rs:179` evaluates `asin(ratio) - (pi/2 - el)`, i.e. `-psi`, and then
/// clamps at zero, so every pierce point collapses onto the receiver.
#[test]
fn compute_ipp_displaces_the_pierce_point_away_from_the_receiver() {
    let (lat, lon) = compute_ipp(0.0, 0.0, 0.0, 30f64.to_radians());
    assert_close(lat, 4.822_339_662, 1e-6);
    assert_close(lon, 0.0, 1e-12);
    // Due east the displacement is a pure longitude step of the same size.
    let (lat_e, lon_e) = compute_ipp(0.0, 0.0, std::f64::consts::FRAC_PI_2, 30f64.to_radians());
    assert_close(lat_e, 0.0, 1e-12);
    assert_close(lon_e, 4.822_339_662, 1e-6);
    // Low elevations step further north (the 0.05 rad floor gives 15.919 deg).
    assert_close(compute_ipp(0.0, 0.0, 0.0, 0.05).0, 15.919_114_839, 1e-6);
}

// -------------------------------------------------------- vrs synthesis --

#[test]
fn synthesize_vrs_epoch_is_a_no_op_when_the_vrs_sits_on_the_master() {
    let surface = gradient_surface();
    let master_pos = site_ecef(0);
    let mut master = master_epoch();
    master.satellites.push(quad_obs(sat(2), 20_030_000.0, 20_040_000.0, 1.06e8, 8.3e7));
    // With master_pos == vrs_pos the Sagnac-rotated satellite positions, the
    // geometric range and the ENU offset are all identical, so nothing shifts.
    let out = synthesize_vrs_epoch(&master, master_pos, master_pos, &surface, &[circ_eph(1, 30.0)]);
    assert_eq!(out.time, master.time);
    assert_eq!(out.satellites.len(), 1, "G02 has no ephemeris and is dropped");
    assert_eq!(out.satellites[0].sat, sat(1));
    for band in [1u8, 2] {
        for kind in [ObsType::Pseudorange, ObsType::CarrierPhase] {
            assert_eq!(value_of(&out.satellites[0], kind, band),
                       value_of(&master.satellites[0], kind, band));
        }
    }
    assert_eq!(value_of(&out.satellites[0], ObsType::Doppler, 1), -1234.5);
    // With no ephemeris at all, every satellite is dropped.
    let empty = synthesize_vrs_epoch(&master, master_pos, site_ecef(1), &surface, &[]);
    assert!(empty.satellites.is_empty());
}

#[test]
fn synthesize_vrs_epoch_separates_the_geometric_tropo_and_iono_shifts() {
    let surface = gradient_surface();
    let master_pos = site_ecef(0);
    // 100 m east of the master keeps every term metre-scale, so the separation
    // between the L1 and L2 shifts stays exact instead of being swamped by a
    // mega-metre geometric term.
    let vrs_pos = master_pos + Vector3::new(100.0, 0.0, 0.0);
    let du = ecef_delta_to_enu(vrs_pos, master_pos, ecef_to_llh(master_pos));
    // d/dEast = 2e-4, d/dNorth = -3e-4, d/dIono East = 1e-4, North = 5e-5;
    // the Up component is unused by the plain dot product.
    let d_tropo = 2.0e-4 * du.x - 3.0e-4 * du.y;
    let d_iono = 1.0e-4 * du.x + 5.0e-5 * du.y;
    assert!(d_tropo.abs() > 1e-3 && d_iono.abs() > 1e-4, "fixture must be non-trivial");
    let out = synthesize_vrs_epoch(&master_epoch(), master_pos, vrs_pos, &surface, &[circ_eph(1, 30.0)]);
    let syn = &out.satellites[0];
    // Pseudorange shifts by d_geom + d_tropo + I, with the L2 term dispersively
    // scaled by gamma = (f1/f2)^2. Recovering d_geom from the L1 shift makes
    // the L2 relation an independent check rather than a restatement.
    let d_geom = value_of(syn, ObsType::Pseudorange, 1) - 20_010_000.0 - d_tropo - d_iono;
    assert_close(value_of(syn, ObsType::Pseudorange, 2) - 20_020_000.0,
        d_geom + d_tropo + GAMMA * d_iono, 1e-6);
    // Carrier phase shifts by the ionosphere-free residual over the wavelength.
    assert_close(value_of(syn, ObsType::CarrierPhase, 1) - 1.05e8,
        (d_geom + d_tropo - d_iono) / lambda1(), 1e-6);
    assert_close(value_of(syn, ObsType::CarrierPhase, 2) - 8.2e7,
        (d_geom + d_tropo - GAMMA * d_iono) / lambda2(), 1e-6);
    // Doppler is not shifted at all.
    assert_eq!(value_of(syn, ObsType::Doppler, 1), -1234.5);
}

#[test]
fn synthesize_vrs_epoch_geometric_shift_follows_the_line_of_sight() {
    let surface = gradient_surface();
    let master_pos = site_ecef(0);
    let eph = [circ_eph(1, 30.0)];
    // Recover the geometric term from the synthesised L1 pseudorange: the
    // gradient surface contributes exactly d_tropo + d_iono.
    let geom_at = |vrs_pos: Vector3<f64>| {
        let out = synthesize_vrs_epoch(&master_epoch(), master_pos, vrs_pos, &surface, &eph);
        let du = ecef_delta_to_enu(vrs_pos, master_pos, ecef_to_llh(master_pos));
        value_of(&out.satellites[0], ObsType::Pseudorange, 1) - 20_010_000.0
            - 2.0e-4 * du.x + 3.0e-4 * du.y - 1.0e-4 * du.x - 5.0e-5 * du.y
    };
    // A VRS 1 km along the master -> satellite ray is 1 km closer. The Sagnac
    // rotation moves the transmit point by ~1 cm over this baseline.
    assert_close(geom_at(master_pos + 1000.0 * los_dir(master_pos, 30.0)), -1000.0, 1.0);
    // 1 km strictly across the ray changes the range by only d^2/(2R) =
    // 1e6 / (2 * 2.3e7) = 0.022 m; a first-order term would be ~1000 m.
    let (up, _) = local_up_down_east(master_pos);
    let across = geom_at(master_pos + 1000.0 * los_dir(master_pos, 30.0).cross(&up).normalize());
    assert!(across.abs() < 0.5, "perpendicular shift must be second order, got {across} m");
}

// --------------------------------------------------------- vrs streaming --

#[test]
fn synthesize_vrs_stream_reports_a_missing_master() {
    let synth = VrsSynthesizer::new("M", site_ecef(0));
    let other = cors("OTHER", 1, 1, |_| vec![quad_obs(sat(1), 2.0e7, 2.1e7, 1.0e8, 8.0e7)]);
    let err = synth.synthesize_vrs_stream(&[other], site_ecef(1), &[]).expect_err("no master");
    assert!(matches!(err, crate::spatial::delaunay::EngineError::DegenerateMesh(_)),
        "unexpected error kind: {err:?}");
    assert!(err.to_string().contains("Master not found"), "{err}");
}

#[test]
fn synthesize_vrs_stream_mirrors_the_master_epoch_count_and_times() {
    let sats = |_: usize| vec![quad_obs(sat(1), 20_010_000.0, 20_020_000.0, 1.05e8, 8.2e7)];
    let stations: Vec<CorsStation> = (0..3)
        .map(|i| { let id = if i == 0 { "M".to_string() } else { format!("S{i}") }; cors(&id, i, 4, sats) })
        .collect();
    // 10 m east of the master: the geometric shift cannot exceed the station
    // separation, so the synthesised pseudorange stays within 10 m of it.
    let rover = site_ecef(0) + Vector3::new(10.0, 0.0, 0.0);
    let synth = VrsSynthesizer::new("M", site_ecef(0));
    let out = synth.synthesize_vrs_stream(&stations, rover, &[circ_eph(1, 30.0)]).expect("stream");
    assert_eq!(out.len(), 4, "one VRS epoch per master epoch");
    for (i, epoch) in out.iter().enumerate() {
        assert_eq!(epoch.time, GpsTime::new(2000, 100.0 + 30.0 * i as f64));
        assert_eq!(epoch.satellites.len(), 1);
        let master_p = value_of(&stations[0].epochs[i].satellites[0], ObsType::Pseudorange, 1);
        let vrs_p = value_of(&epoch.satellites[0], ObsType::Pseudorange, 1);
        assert!((vrs_p - master_p).abs() <= 10.0);
    }
}

#[test]
fn synthesize_vrs_stream_works_with_too_few_stations_for_a_delaunay_mesh() {
    // Two stations cannot be triangulated, so the surface keeps its zero
    // gradients; the stream must still build with a purely geometric shift.
    let sats = |_: usize| vec![quad_obs(sat(1), 20_010_000.0, 20_020_000.0, 1.05e8, 8.2e7)];
    let stations = vec![cors("M", 0, 2, sats), cors("S1", 1, 2, sats)];
    let synth = VrsSynthesizer::new("M", site_ecef(0));
    let out = synth.synthesize_vrs_stream(&stations, site_ecef(1), &[circ_eph(1, 30.0)])
        .expect("stream builds without a Delaunay mesh");
    assert_eq!(out.len(), 2);
    // The VRS sits 15 deg of longitude from the master, so the geometric shift
    // alone must be non-trivial.
    let master_p = value_of(&stations[0].epochs[0].satellites[0], ObsType::Pseudorange, 1);
    assert!((value_of(&out[0].satellites[0], ObsType::Pseudorange, 1) - master_p).abs() > 1.0);
}

#[test]
fn frequency_ratio_matches_the_broadcast_l1_l2_pairing() {
    // Guards the gamma and wavelength constants the other tests assume.
    let (f1, f2) = satellite_frequencies(sat(1), 0);
    assert_eq!((f1, f2), (F1, F2));
    assert_close(f1 / f2, 77.0 / 60.0, 1e-12);
    assert_close((f1 / f2).powi(2), GAMMA, 1e-12);
    assert_close(lambda1(), 0.190_293_672_798, 1e-12);
    assert_close(lambda2(), 0.244_210_213_425, 1e-12);
}
