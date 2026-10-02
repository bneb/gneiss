//! Continuous Phase Windup Modeling for Circularly Polarized GNSS Transmissions.
//!
//! Implements the Wu et al. (1993) phase windup formulation for RHCP satellite
//! transmissions with continuous cycle wrap tracking across passes.
//!
//! # Frame convention (the part that is easy to get backwards)
//!
//! With `s` the unit line of sight **from the receiver toward the satellite**
//! and `P(v) = v - (s.v)s` the projection onto the plane normal to `s`,
//!
//! ```text
//!   D' = P(X) + s x Y          transmitter effective dipole   (PLUS)
//!   D  = P(x) - s x y          receiver    effective dipole   (MINUS)
//!   phi_u = atan2( s . (D x D'), D' . D )
//! ```
//!
//! Skeens, York, Petrov, Herrity & Ji-Cathriner (2025), *J. Geodesy* 99:71,
//! doi:10.1007/s00190-025-01993-z, Eqs. (2),(4),(7),(8),(26),(28). The same
//! construction appears in ESA Navipedia "Carrier Phase Wind-up Effect" Eq. (4)
//! and in RTKLIB `windupcorr`; the latter two use the opposite line-of-sight
//! direction, which negates *both* cross terms and so describes the very same
//! physical vectors.
//!
//! The opposite signs are forced by geometry, not a convention: both boresights
//! point "toward the Earth" (satellite nadir, receiver zenith) and both triads
//! are right-handed **about their own boresight**. Proof that this matters:
//! when the line of sight coincides with the receiver boresight, `s x y = -x`,
//! so the MINUS form gives `|D| = 2` while a PLUS form gives `|D| = 0` — a
//! spurious degeneracy.
//!
//! Consequence to respect at every call site: the receiver pair must be
//! **right-handed about the antenna boresight**, i.e. `(East, North)` — never
//! `(North, East)`, whose cross product is `-Up`. Getting this backwards does
//! not fail loudly: it negates `d(phi)/d(theta)` and moves the answer by up to
//! half a cycle.

use core::f64::consts::PI;
use nalgebra::Vector3;

/// Degeneracy floor for any vector that must define a direction.
const EPS: f64 = 1e-6;

#[inline]
fn vec_norm(v: &Vector3<f64>) -> f64 {
    libm::sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
}

/// Normalise, or `None` when the vector is too short to define a direction.
#[inline]
fn unit_or_none(v: &Vector3<f64>) -> Option<Vector3<f64>> {
    let n = vec_norm(v);
    if n < EPS {
        None
    } else {
        Some(v / n)
    }
}

/// Projection of `v` onto the plane normal to the line of sight `s`.
#[inline]
fn project_normal(s: &Vector3<f64>, v: &Vector3<f64>) -> Vector3<f64> {
    v - s * s.dot(v)
}

/// Transmitter effective dipole: `Z` is satellite nadir, `Y` is transverse to
/// the satellite-to-Sun direction, `X = Y x Z` completes the right-handed triad.
fn transmitter_dipole(
    sat_pos: &Vector3<f64>,
    sun_pos: &Vector3<f64>,
    s: &Vector3<f64>,
) -> Option<Vector3<f64>> {
    let z = unit_or_none(&(-*sat_pos))?;
    // Satellite-to-Sun, not geocentric-to-Sun.
    let e_sun = unit_or_none(&(sun_pos - sat_pos))?;
    let y = unit_or_none(&z.cross(&e_sun))?;
    let x = y.cross(&z);
    Some(project_normal(s, &x) + s.cross(&y))
}

/// Receiver effective dipole from a right-handed `(boresight, transverse)` pair.
fn receiver_dipole(
    s: &Vector3<f64>,
    boresight: &Vector3<f64>,
    transverse: &Vector3<f64>,
) -> Vector3<f64> {
    project_normal(s, boresight) - s.cross(transverse)
}

/// Signed windup angle in `(-pi, pi]`, or `None` at a dipole degeneracy.
fn windup_angle(
    s: &Vector3<f64>,
    d_prime: &Vector3<f64>,
    d_rx: &Vector3<f64>,
) -> Option<f64> {
    let (lp, lr) = (vec_norm(d_prime), vec_norm(d_rx));
    if lp < EPS || lr < EPS {
        return None;
    }
    let dot = (d_prime.dot(d_rx) / (lp * lr)).clamp(-1.0, 1.0);
    let turn = s.dot(&d_rx.cross(d_prime)); // >= 0 selects +acos
    let sign = if turn >= 0.0 { 1.0 } else { -1.0 };
    Some(sign * libm::acos(dot))
}

/// Continuous phase windup tracker across epochs.
#[derive(Debug, Clone, Copy, Default)]
pub struct PhaseWindupTracker {
    pub prev_windup_rad: f64,
    pub is_initialized: bool,
}

impl PhaseWindupTracker {
    pub fn new() -> Self {
        Self {
            prev_windup_rad: 0.0,
            is_initialized: false,
        }
    }

    /// Computes continuous phase windup correction in radians with optional receiver attitude.
    ///
    /// `rx_east` is the *boresight-side* dipole axis and `rx_north` the
    /// transverse one: `(East, North)` is right-handed about the zenith
    /// boresight, `(North, East)` is not. See the module docs.
    pub fn update_with_attitude(
        &mut self,
        sat_pos: &Vector3<f64>,
        sun_pos: &Vector3<f64>,
        rx_pos: &Vector3<f64>,
        rx_north: &Vector3<f64>,
        rx_east: &Vector3<f64>,
        r_body_to_ecef: Option<&nalgebra::Matrix3<f64>>,
    ) -> f64 {
        let (x_ax, y_ax) = match r_body_to_ecef {
            Some(r) => (r * rx_east, r * rx_north),
            None => (*rx_east, *rx_north),
        };
        self.step(sat_pos, sun_pos, rx_pos, &x_ax, &y_ax)
    }

    /// Shared core: build both dipoles, measure the angle, unwrap to continuity.
    fn step(
        &mut self,
        sat_pos: &Vector3<f64>,
        sun_pos: &Vector3<f64>,
        rx_pos: &Vector3<f64>,
        rx_x: &Vector3<f64>,
        rx_y: &Vector3<f64>,
    ) -> f64 {
        // s points receiver -> satellite.
        let Some(s) = unit_or_none(&(sat_pos - rx_pos)) else {
            return self.prev_windup_rad;
        };
        let Some(d_prime) = transmitter_dipole(sat_pos, sun_pos, &s) else {
            return self.prev_windup_rad;
        };
        let d_rx = receiver_dipole(&s, rx_x, rx_y);
        let Some(phi) = windup_angle(&s, &d_prime, &d_rx) else {
            return self.prev_windup_rad;
        };
        if !self.is_initialized {
            self.prev_windup_rad = phi;
            self.is_initialized = true;
        } else {
            self.prev_windup_rad = self.unwrap(phi);
        }
        self.prev_windup_rad
    }

    /// Lift the wrapped `phi` by whole turns so it stays nearest the previous value.
    fn unwrap(&self, phi: f64) -> f64 {
        let turns = libm::round((self.prev_windup_rad - phi) / (2.0 * PI));
        phi + 2.0 * PI * turns
    }

    /// Computes continuous phase windup correction in radians for a static/topocentric receiver.
    pub fn update(
        &mut self,
        sat_pos: &Vector3<f64>,
        sun_pos: &Vector3<f64>,
        rx_pos: &Vector3<f64>,
        _rx_up: &Vector3<f64>,
        rx_north: &Vector3<f64>,
        rx_east: &Vector3<f64>,
    ) -> f64 {
        self.update_with_attitude(sat_pos, sun_pos, rx_pos, rx_north, rx_east, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Matrix3;

    const AU: f64 = 1.495_978_707e11;
    const RE: f64 = 6_378_137.0;

    fn enu(lat_deg: f64, lon_deg: f64) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>) {
        let (la, lo) = (lat_deg.to_radians(), lon_deg.to_radians());
        let up = Vector3::new(la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin());
        let east = Vector3::new(-lo.sin(), lo.cos(), 0.0);
        let north = Vector3::new(-la.sin() * lo.cos(), -la.sin() * lo.sin(), la.cos());
        (up, east, north)
    }

    fn geodetic(lat_deg: f64, lon_deg: f64, h: f64) -> Vector3<f64> {
        let (la, lo) = (lat_deg.to_radians(), lon_deg.to_radians());
        let e2 = 6.694_379_990_141_316e-3;
        let n = RE / libm::sqrt(1.0 - e2 * la.sin() * la.sin());
        Vector3::new(
            (n + h) * la.cos() * lo.cos(),
            (n + h) * la.cos() * lo.sin(),
            (n * (1.0 - e2) + h) * la.sin(),
        )
    }

    /// Three independent scenes: mid-latitude, equatorial and southern-hemisphere.
    type Scene = (&'static str, f64, f64, f64, f64, f64);

    fn scenes() -> [Scene; 3] {
        [
            ("mid-lat", 45.0, 9.0, 45.0, 95.0, 20.0),
            ("equator", 0.0, 0.0, 10.0, 60.0, 80.0),
            ("south", -30.0, -70.0, -20.0, 120.0, 150.0),
        ]
    }

    /// Proper right-handed rotation of `theta` about `axis` (Rodrigues), built
    /// row-major so the element order is unambiguous.
    fn rot(axis: Vector3<f64>, theta: f64) -> Matrix3<f64> {
        let k = unit_or_none(&axis).unwrap();
        let (c, s) = (theta.cos(), theta.sin());
        let t = 1.0 - c;
        Matrix3::from_row_slice(&[
            t * k.x * k.x + c,
            t * k.x * k.y - s * k.z,
            t * k.x * k.z + s * k.y,
            t * k.x * k.y + s * k.z,
            t * k.y * k.y + c,
            t * k.y * k.z - s * k.x,
            t * k.x * k.z - s * k.y,
            t * k.y * k.z + s * k.x,
            t * k.z * k.z + c,
        ])
    }

    fn sun_at(lambda_deg: f64) -> Vector3<f64> {
        let l = lambda_deg.to_radians();
        Vector3::new(AU * l.cos(), AU * l.sin(), 0.3e11)
    }

    #[allow(clippy::type_complexity)] // test fixture: (rx, sat, sun, east, north)
    fn scene(
        lat: f64,
        lon: f64,
        slat: f64,
        slon: f64,
        slambda: f64,
    ) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>, Vector3<f64>, Vector3<f64>) {
        let (_, east, north) = enu(lat, lon);
        (
            geodetic(lat, lon, 100.0),
            geodetic(slat, slon, 20_200e3),
            sun_at(slambda),
            east,
            north,
        )
    }

    /// Independent reference: a literal transcription of the published
    /// equations, written from the paper rather than from the code under test.
    fn reference_windup(
        sat: &Vector3<f64>,
        sun: &Vector3<f64>,
        rx: &Vector3<f64>,
        east: &Vector3<f64>,
        north: &Vector3<f64>,
    ) -> Option<f64> {
        let u = |v: &Vector3<f64>| unit_or_none(v);
        let proj = |s: &Vector3<f64>, v: &Vector3<f64>| v - s * s.dot(v);
        let s = u(&(sat - rx))?; // Eq. (2) receiver -> satellite
        let z = u(&(-*sat))?; // Eq. (26) satellite nadir
        let e_sun = u(&(sun - sat))?; // Eq. (26) satellite -> Sun
        let y = u(&z.cross(&e_sun))?; // Eq. (26) transverse
        let x = y.cross(&z); // Eq. (26) aligned
        let t = proj(&s, &x) + s.cross(&y); // Eq. (8) transmitter, PLUS
        let r = proj(&s, east) - s.cross(north); // Eq. (8) receiver, MINUS
        Some(libm::atan2(s.dot(&r.cross(&t)), t.dot(&r))) // Eq. (7)
    }

    #[test]
    fn matches_the_published_formulation_on_every_test_scene() {
        for (name, lat, lon, slat, slon, sl) in scenes() {
            let (rx, sat, sun, east, north) = scene(lat, lon, slat, slon, sl);
            let want = reference_windup(&sat, &sun, &rx, &east, &north).unwrap();
            let mut t = PhaseWindupTracker::new();
            let got = t.update(&sat, &sun, &rx, &east, &north, &east); // up unused
            assert!(
                (got - want).abs() < 1e-12,
                "{name}: windup {got} != reference {want} (delta {:e} rad)",
                got - want
            );
        }
    }

    #[test]
    fn matches_the_published_formulation_through_the_attitude_path() {
        for (name, lat, lon, slat, slon, sl) in scenes() {
            let (rx, sat, sun, east, north) = scene(lat, lon, slat, slon, sl);
            let want = reference_windup(&sat, &sun, &rx, &east, &north).unwrap();
            let yaw = rot(Vector3::new(0.0, 0.0, 1.0), 0.4);
            let mut t = PhaseWindupTracker::new();
            let got = t.update_with_attitude(&sat, &sun, &rx, &north, &east, None);
            let yawed = t.update_with_attitude(&sat, &sun, &rx, &north, &east, Some(&yaw));
            let want_att = reference_windup(&sat, &sun, &rx, &(yaw * east), &(yaw * north)).unwrap();
            assert!((got - want).abs() < 1e-12, "{name}: plain path {got} != {want}");
            assert!(
                (yawed - want_att).abs() < 1e-12,
                "{name}: attitude path {yawed} != reference {want_att}"
            );
        }
    }

    /// The handedness identity that fixes the receiver triad: East x North = +Up
    /// but North x East = -Up. A left-handed receiver pair silently negates the
    /// yaw response, which is exactly what this test pins down.
    #[test]
    fn receiver_dipole_pair_must_be_right_handed_about_the_boresight() {
        for (name, lat, lon, slat, slon, sl) in scenes() {
            let (up, east, north) = enu(lat, lon);
            assert!(vec_norm(&(east.cross(&north) - up)) < 1e-12, "{name}: E x N != Up");
            assert!(vec_norm(&(north.cross(&east) + up)) < 1e-12, "{name}: N x E != -Up");
            let (rx, sat, sun, _, _) = scene(lat, lon, slat, slon, sl);
            let mut t = PhaseWindupTracker::new();
            let good = t.update(&sat, &sun, &rx, &up, &north, &east);
            let mut t2 = PhaseWindupTracker::new();
            let swapped = t2.update(&sat, &sun, &rx, &up, &east, &north);
            assert!(
                (good - swapped).abs() > 1e-3,
                "{name}: swapping (East, North) must change the answer"
            );
        }
    }

    /// A rigid +theta yaw of the antenna about its boresight must change the
    /// windup by exactly -theta radians: the published Phi1 convention has
    /// d(phi)/d(theta) = -1. Asserted *signed*; the magnitude-only version of
    /// this test passes for +1 as well and cannot catch a left-handed receiver
    /// frame.
    #[test]
    fn rigid_yaw_about_the_boresight_changes_windup_by_minus_theta() {
        let (up, east, north) = enu(45.0, 9.0);
        let (rx, sat, sun, _, _) = scene(45.0, 9.0, 45.0, 95.0, 20.0);
        let theta = 0.3_f64;
        let yaw = rot(up, theta);
        let mut a = PhaseWindupTracker::new();
        let base = a.update(&sat, &sun, &rx, &up, &north, &east);
        let mut b = PhaseWindupTracker::new();
        let yawed = b.update_with_attitude(&sat, &sun, &rx, &north, &east, Some(&yaw));
        let delta = yawed - base;
        assert!(
            (delta + theta).abs() < 1e-9,
            "d(phi)/d(theta) must be -1, got {:e}",
            delta / theta
        );
    }

    /// Windup is a purely geometric angle: rigidly rotating the whole scene
    /// (satellite, Sun, receiver and antenna triad) must leave it invariant.
    #[test]
    fn invariant_under_a_rigid_rotation_of_the_entire_scene() {
        let (up, east, north) = enu(-30.0, -70.0);
        let (rx, sat, sun, _, _) = scene(-30.0, -70.0, -20.0, 120.0, 150.0);
        let mut a = PhaseWindupTracker::new();
        let base = a.update(&sat, &sun, &rx, &up, &north, &east);
        for (axis, th) in [
            (Vector3::new(0.3, 0.5, 0.9), 0.7),
            (Vector3::new(-1.1, 0.2, 0.7), 1.3),
            (Vector3::new(2.0, -0.4, -0.6), -2.1),
        ] {
            let rot = rot(axis, th);
            let mut t = PhaseWindupTracker::new();
            let got = t.update(
                &(rot * sat),
                &(rot * sun),
                &(rot * rx),
                &(rot * up),
                &(rot * north),
                &(rot * east),
            );
            assert!(
                (got - base).abs() < 1e-9,
                "rotation by {th} about {axis:?} changed windup by {:e}",
                got - base
            );
        }
    }

    /// Integer-turn continuity: sweeping the attitude through a full turn must
    /// track continuously and accumulate exactly -2*pi.
    #[test]
    fn unwrapping_tracks_across_a_full_turn_of_attitude() {
        let (up, east, north) = enu(45.0, 9.0);
        let (rx, sat, sun, _, _) = scene(45.0, 9.0, 45.0, 95.0, 20.0);
        let mut tracker = PhaseWindupTracker::new();
        let base = tracker.update(&sat, &sun, &rx, &up, &north, &east);
        let mut prev = base;
        let mut max_step = 0.0_f64;
        let n = 512;
        for i in 1..=n {
            let th = 2.0 * PI * (i as f64) / (n as f64);
            let w = tracker.update_with_attitude(&sat, &sun, &rx, &north, &east, Some(&rot(up, th)));
            max_step = max_step.max((w - prev).abs());
            prev = w;
        }
        assert!(max_step < 0.1, "unwrapped windup stepped by {max_step} rad");
        // d(phi)/d(theta) = -1 exactly, so a full turn accumulates exactly -2*pi
        // ON TOP OF the starting angle (which itself lies in (-pi, pi]).
        assert!(
            ((prev - base) + 2.0 * PI).abs() < 1e-6,
            "one 360-deg yaw must accumulate exactly -2*pi, got {:e}",
            prev - base
        );
    }

    /// A satellite at zenith on the geocentric ray puts the receiver boresight
    /// on the line of sight. The published receiver dipole is then `D = 2x`
    /// (well conditioned); a left-handed frame or a PLUS sign would collapse
    /// `|D|` to 0 and report a spurious degeneracy.
    #[test]
    fn zenith_satellite_is_ill_conditioned_not_degenerate_for_the_receiver() {
        let (up, east, north) = enu(0.0, 0.0);
        let rx = geodetic(0.0, 0.0, 0.0);
        let sat = Vector3::new(RE + 20_200e3, 0.0, 0.0);
        let sun = sun_at(20.0);
        let s = unit_or_none(&(sat - rx)).unwrap();
        let d = receiver_dipole(&s, &east, &north);
        assert!(
            (vec_norm(&d) - 2.0).abs() < 1e-9,
            "at zenith the receiver dipole must have |D| = 2, got {}",
            vec_norm(&d)
        );
        let mut t = PhaseWindupTracker::new();
        assert!(t.update(&sat, &sun, &rx, &up, &north, &east).is_finite());
    }

    /// Sun-synchronous geometry collapses the *transmitter* transverse axis.
    /// That is a genuine degeneracy: hold the previous value, never emit NaN.
    #[test]
    fn transmitter_degeneracy_holds_the_previous_value() {
        let (up, east, north) = enu(0.0, 0.0);
        let (rx, sat, sun, _, _) = scene(0.0, 0.0, 10.0, 40.0, 35.0);
        let mut t = PhaseWindupTracker::new();
        assert!(t.update(&sat, &sun, &rx, &up, &north, &east).is_finite());
        let nadir = -unit_or_none(&sat).unwrap() * AU;
        let held = t.update(&sat, &nadir, &rx, &up, &north, &east);
        assert!(held.is_finite(), "degeneracy must not emit NaN");
        assert_eq!(held, t.prev_windup_rad, "degeneracy must hold the last value");
    }

    /// Adversarial inputs must never produce NaN.
    #[test]
    fn degenerate_inputs_never_produce_nan() {
        let (up, east, north) = enu(45.0, 9.0);
        let sat = geodetic(45.0, 95.0, 20_200e3);
        let rx = geodetic(45.0, 9.0, 100.0);
        let sun = sun_at(20.0);
        for (tag, s, u, r) in [
            ("zero sun", sat, Vector3::zeros(), rx),
            ("zero sat", Vector3::zeros(), sun, rx),
            ("rx at sat", sat, sun, sat),
            ("all zero", Vector3::zeros(), Vector3::zeros(), Vector3::zeros()),
        ] {
            let mut t = PhaseWindupTracker::new();
            let w = t.update(&s, &u, &r, &up, &north, &east);
            assert!(w.is_finite(), "{tag}: windup = {w}");
        }
    }

    /// Small satellite motion must not step the unwrapped windup.
    #[test]
    fn small_satellite_motion_stays_continuous() {
        let (up, east, north) = enu(0.0, 0.0);
        let sun = Vector3::new(AU, 0.0, 0.0);
        let rx = Vector3::new(RE, 0.0, 0.0);
        let mut t = PhaseWindupTracker::new();
        let w1 = t.update(&Vector3::new(0.0, 26_560_000.0, 0.0), &sun, &rx, &up, &north, &east);
        let w2 = t.update(&Vector3::new(100.0, 26_560_000.0, 0.0), &sun, &rx, &up, &north, &east);
        assert!(w1.is_finite() && w2.is_finite());
        assert!((w2 - w1).abs() < 0.1, "windup stepped {w1:.4} -> {w2:.4}");
    }
}