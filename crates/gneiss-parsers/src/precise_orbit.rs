//! Precise orbit interpolation from SP3 ephemeris data.
//!
//! Wraps parsed SP3 epochs and provides position/clock lookup at arbitrary
//! epochs via degree-8 Lagrange polynomial interpolation. This replaces
//! broadcast Keplerian propagation with IGS final/rapid products for cm-level
//! orbit accuracy.
//!
//! Frame handling: SP3 positions are Earth-fixed at their own epoch tag, so
//! each node is rotated about Z by the Earth rotation angle into the target
//! epoch's frame before fitting (Sagnac pre-rotation). Clocks use two-point
//! linear interpolation between bracketing epochs; degree-8 fits ring on
//! noisy 15-min clock data.

use crate::sp3::{Sp3Epoch, Sp3Record};
use gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S as OMEGA_E_GPS;
use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::time::GpsTime;
use nalgebra::Vector3;
use std::collections::HashMap;

/// One sample prepared for fitting: time (s since GPS week start), ECEF
/// position rotated into the target-epoch frame, and raw clock offset (s).
type Node = (f64, Vector3<f64>, f64);

/// Interpolated precise orbit store.
#[derive(Debug, Clone)]
pub struct PreciseOrbit {
    /// Per-satellite time-ordered samples.
    tracks: HashMap<String, Vec<(GpsTime, Sp3Record)>>,
}

impl PreciseOrbit {
    /// Build from parsed SP3 epochs. Keys are SV identifiers ("G01" etc).
    pub fn new(epochs: Vec<Sp3Epoch>) -> Self {
        let mut tracks: HashMap<String, Vec<(GpsTime, Sp3Record)>> = HashMap::new();
        for epoch in epochs {
            for (sv, rec) in epoch.records {
                tracks.entry(sv).or_default().push((epoch.time, rec));
            }
        }
        for v in tracks.values_mut() {
            v.sort_by_key(|(t, _)| (t.week, (t.tow * 1000.0) as i64));
        }
        PreciseOrbit { tracks }
    }

    /// Interpolated position (m) in the ECEF frame AT `t`, and clock offset
    /// (s). Position is Lagrange-fitted through nodes pre-rotated by the
    /// Earth rotation angle into `t`'s frame; the clock is interpolated
    /// linearly between bracketing epochs. Returns None if satellite unknown
    /// or insufficient data.
    pub fn position_at(&self, sv: &str, t: GpsTime) -> Option<(Vector3<f64>, f64)> {
        self.interpolate(sv, t)
    }

    /// Interpolated position (m), velocity (m/s), and clock offset (s) including
    /// periodic relativistic eccentricity clock correction (-\frac{2 r \cdot v}{c^2}).
    pub fn position_and_velocity_at(
        &self,
        sv: &str,
        t: GpsTime,
    ) -> Option<(Vector3<f64>, Vector3<f64>, f64)> {
        let (pos, clk_base) = self.interpolate(sv, t)?;
        let dt = 0.5; // seconds
        let t_fwd = GpsTime::new(t.week, t.tow + dt);
        let t_bwd = GpsTime::new(t.week, t.tow - dt);
        let (pos_fwd, _) = self.interpolate(sv, t_fwd)?;
        let (pos_bwd, _) = self.interpolate(sv, t_bwd)?;
        let vel = (pos_fwd - pos_bwd) / (2.0 * dt);
        let rel_corr = -2.0 * pos.dot(&vel) / (SPEED_OF_LIGHT_M_S * SPEED_OF_LIGHT_M_S);
        Some((pos, vel, clk_base + rel_corr))
    }

    /// As [`PreciseOrbit::position_at`], but when a receiver position hint is
    /// supplied the satellite state is evaluated at TRANSMIT time rather than
    /// receiver time: a first pass at `t` yields a rough range to the hint,
    /// then the returned state is re-evaluated at `t - range/c`. At GPS
    /// orbital speed this accounts for the 270-350 m the satellite travels
    /// during the ~70-90 ms signal flight; the residual after one refinement
    /// is second order (~(v/c) * 90 ms ~ 4 mm, and shrinks geometrically on
    /// further passes). With `None`, `t` is used directly.
    pub fn position_at_with_hint(
        &self,
        sv: &str,
        t: GpsTime,
        rx_pos_hint: Option<Vector3<f64>>,
    ) -> Option<(Vector3<f64>, f64)> {
        let rough = self.interpolate(sv, t)?;
        match rx_pos_hint {
            None => Some(rough),
            Some(rx_pos) => {
                let travel_s = (rx_pos - rough.0).norm() / SPEED_OF_LIGHT_M_S;
                self.interpolate(sv, t - travel_s)
            }
        }
    }

    /// Number of satellites in this store.
    pub fn len(&self) -> usize {
        self.tracks.len()
    }

    pub fn is_empty(&self) -> bool {
        self.tracks.is_empty()
    }

    /// Windowed interpolation shared by every public lookup path.
    fn interpolate(&self, sv: &str, t: GpsTime) -> Option<(Vector3<f64>, f64)> {
        let track = self.tracks.get(sv)?;
        let window = select_window(track, t)?;
        let target_s = s_time(t);
        let nodes = rotate_nodes_to_frame(window, target_s);
        let ts: Vec<f64> = nodes.iter().map(|n| n.0).collect();
        let px = lagrange(&ts, &component_values(&nodes, 0), target_s);
        let py = lagrange(&ts, &component_values(&nodes, 1), target_s);
        let pz = lagrange(&ts, &component_values(&nodes, 2), target_s);
        let clk = linear_clock(&ts, &component_values(&nodes, 3), target_s);
        Some((Vector3::new(px, py, pz), clk))
    }
}

/// Select up to `degree + 1` samples centred on the insertion point of `t`.
fn select_window(
    track: &[(GpsTime, Sp3Record)],
    t: GpsTime,
) -> Option<&[(GpsTime, Sp3Record)]> {
    if track.is_empty() {
        return None;
    }
    // Find insertion point and select degree+1 bracketing samples.
    let degree = 8.min(track.len());
    let half = degree / 2;
    let mut start = 0;
    for (i, (st, _)) in track.iter().enumerate() {
        if st.week > t.week || (st.week == t.week && st.tow > t.tow) {
            start = i.saturating_sub(half);
            break;
        }
        start = i;
    }
    let end = (start + degree + 1).min(track.len());
    let start = end.saturating_sub(degree + 1);
    let window = &track[start..end];
    if window.is_empty() { None } else { Some(window) }
}

/// Express each node's ECEF position in the TARGET epoch's earth-fixed frame.
///
/// An SP3 record tagged t_i is valid in the ECEF frame of orientation
/// theta(t_i) = w_e * t_i. Given p_i = R_z(-theta_i) * p_inertial, the target-
/// frame coordinates are R_z(theta_i - theta*) * p_i, so each node is rotated
/// about Z by phi = w_e * (t_node - t_target):
///   q_x = cos(phi) * x - sin(phi) * y
///   q_y = sin(phi) * x + cos(phi) * y
/// Clocks are frame-independent and pass through untouched.
fn rotate_nodes_to_frame(window: &[(GpsTime, Sp3Record)], target_s: f64) -> Vec<Node> {
    window
        .iter()
        .map(|(t, rec)| {
            let phi = OMEGA_E_GPS * (s_time(*t) - target_s);
            let (s, c) = phi.sin_cos();
            let pos = Vector3::new(
                c * rec.position.x - s * rec.position.y,
                s * rec.position.x + c * rec.position.y,
                rec.position.z,
            );
            (s_time(*t), pos, rec.clock_offset)
        })
        .collect()
}

/// Flatten one component (0=x, 1=y, 2=z, otherwise clock) out of `nodes`.
fn component_values(nodes: &[Node], component: usize) -> Vec<f64> {
    nodes
        .iter()
        .map(|n| match component {
            0 => n.1.x,
            1 => n.1.y,
            2 => n.1.z,
            _ => n.2,
        })
        .collect()
}

/// Degree-(n-1) Lagrange interpolation of `vals` over abscissae `ts`.
fn lagrange(ts: &[f64], vals: &[f64], target_s: f64) -> f64 {
    let mut result = 0.0;
    for (i, ti) in ts.iter().enumerate() {
        let mut term = vals[i];
        for (j, tj) in ts.iter().enumerate() {
            if j != i && (ti - tj).abs() >= 1e-12 {
                term *= (target_s - tj) / (ti - tj);
            }
        }
        result += term;
    }
    result
}

/// Two-point linear interpolation over the segment bracketing `target_s`
/// (end segments extend linearly when `target_s` falls outside the window).
/// SP3 clocks are noisy at 15-min sampling; high-order polynomial fits ring
/// on that noise instead of smoothing it, so clocks stay linear (the choice
/// RTKLIB makes in preceph.c).
fn linear_clock(ts: &[f64], vals: &[f64], target_s: f64) -> f64 {
    if ts.len() < 2 {
        return vals.first().copied().unwrap_or(0.0);
    }
    let mut k = 0;
    for (i, ti) in ts.iter().enumerate() {
        if ti <= &target_s {
            k = i;
        } else {
            break;
        }
    }
    if k >= ts.len() - 1 {
        k = ts.len() - 2;
    }
    let span = ts[k + 1] - ts[k];
    if span.abs() < 1e-12 {
        return vals[k];
    }
    vals[k] + (vals[k + 1] - vals[k]) * ((target_s - ts[k]) / span)
}

fn s_time(t: GpsTime) -> f64 {
    t.week as f64 * 604800.0 + t.tow
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S as OMEGA_E_GPS;
    use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    use nalgebra::Vector3;

    fn circular_orbit(radius_m: f64, period_s: f64, n: usize, dt: f64) -> Vec<Sp3Epoch> {
        (0..n)
            .map(|i| {
                let t = i as f64 * dt;
                let th = 2.0 * std::f64::consts::PI * t / period_s;
                let pos = Vector3::new(radius_m * th.cos(), radius_m * th.sin(), 0.0);
                let clk = 1e-6 + 1e-9 * t;
                let mut recs = HashMap::new();
                recs.insert("G01".to_string(), Sp3Record { position: pos, clock_offset: clk });
                Sp3Epoch { time: GpsTime::new(2105, t), records: recs }
            })
            .collect()
    }

    #[test]
    fn test_node_interpolation_exact() {
        let dt = 900.0;
        let eps = circular_orbit(26_560_000.0, 43_200.0, 20, dt);
        let orb = PreciseOrbit::new(eps);
        let t = GpsTime::new(2105, 5.0 * dt);
        let (pos, _) = orb.position_at("G01", t).unwrap();
        let expected_theta = 2.0 * std::f64::consts::PI * 5.0 * dt / 43_200.0;
        let expected = Vector3::new(
            26_560_000.0 * expected_theta.cos(),
            26_560_000.0 * expected_theta.sin(),
            0.0,
        );
        assert!(
            (pos - expected).norm() < 1.0,
            "node error {}",
            (pos - expected).norm()
        );
    }

    #[test]
    fn test_between_nodes_sub_metre() {
        let dt = 900.0;
        let eps = circular_orbit(26_560_000.0, 43_200.0, 96, dt);
        let orb = PreciseOrbit::new(eps);
        let mut max_err = 0.0_f64;
        for i in 10..80 {
            let t_mid = i as f64 * dt + dt / 2.0;
            let t = GpsTime::new(2105, t_mid);
            if let Some((pos, _)) = orb.position_at("G01", t) {
                let th = 2.0 * std::f64::consts::PI * t_mid / 43_200.0;
                let exp = Vector3::new(26_560_000.0 * th.cos(), 26_560_000.0 * th.sin(), 0.0);
                max_err = max_err.max((pos - exp).norm());
            }
        }
        assert!(max_err < 1.0, "mid-node error {:.4} m", max_err);
    }

    #[test]
    fn test_clock_linear_exact() {
        let dt = 900.0;
        let eps = circular_orbit(26_560_000.0, 43_200.0, 20, dt);
        let orb = PreciseOrbit::new(eps);
        let t = GpsTime::new(2105, 5.5 * dt);
        let (_, clk) = orb.position_at("G01", t).unwrap();
        let expected = 1e-6 + 1e-9 * 5.5 * dt;
        assert!((clk - expected).abs() < 1e-14);
    }

    #[test]
    fn test_unknown_sat_returns_none() {
        let eps = circular_orbit(26_560_000.0, 43_200.0, 5, 900.0);
        let orb = PreciseOrbit::new(eps);
        assert!(orb.position_at("G99", GpsTime::new(2105, 100.0)).is_none());
    }

    #[test]
    fn test_outside_range_clamps() {
        let eps = circular_orbit(26_560_000.0, 43_200.0, 10, 900.0);
        let orb = PreciseOrbit::new(eps);
        assert!(orb.position_at("G01", GpsTime::new(2104, 100.0)).is_some());
        assert!(orb.position_at("G01", GpsTime::new(2106, 100.0)).is_some());
    }

    /// SP3 nodes of a body FIXED in inertial space: node `i` stores the
    /// ECEF-at-t_i coordinates `R_z(-w_e*t_i) * P0`, i.e. the point appears to
    /// rotate westward through the earth-fixed frame between epochs.
    fn inertial_fixed_orbit(radius_m: f64, n: usize, dt: f64) -> Vec<Sp3Epoch> {
        (0..n)
            .map(|i| {
                let t = i as f64 * dt;
                let th = OMEGA_E_GPS * t;
                let pos = Vector3::new(radius_m * th.cos(), -(radius_m) * th.sin(), 0.0);
                let mut recs = HashMap::new();
                recs.insert(
                    "G01".to_string(),
                    Sp3Record { position: pos, clock_offset: 1e-6 },
                );
                Sp3Epoch { time: GpsTime::new(2105, t), records: recs }
            })
            .collect()
    }

    #[test]
    fn test_sagnac_nodes_rotated_to_target_frame() {
        // Hourly sampling makes the per-node frame drift w_e*dt ~ 0.26 rad, so
        // fitting raw (un-rotated) ECEF nodes cannot recover the target-frame
        // position; rotated nodes must reproduce it exactly.
        let dt = 3600.0;
        let r = 26_560_000.0;
        let orb = PreciseOrbit::new(inertial_fixed_orbit(r, 20, dt));
        let frac = 5.5;
        let t_target = GpsTime::new(2105, frac * dt);
        let (pos, _) = orb.position_at("G01", t_target).unwrap();
        // Truth: ECEF-at-target-frame coordinates of the inertially fixed point.
        let th = OMEGA_E_GPS * frac * dt;
        let expected = Vector3::new(r * th.cos(), -(r) * th.sin(), 0.0);
        let err = (pos - expected).norm();
        assert!(err < 1e-3, "sagnac frame error {:.6} m at midpoint", err);
    }

    #[test]
    fn test_clock_interpolation_is_linear_not_lagrange() {
        // Alternating +/-100 ns clock dither (15-min sampling): degree-8
        // Lagrange rings violently on this Nyquist-frequency content, while
        // linear interpolation returns exactly the bracketing chord average.
        let dt = 900.0;
        let mut eps = circular_orbit(26_560_000.0, 43_200.0, 20, dt);
        let base = 1.0e-6_f64;
        let dither = 1.0e-7_f64;
        for (i, e) in eps.iter_mut().enumerate() {
            let rec = e.records.get_mut("G01").expect("sv record present");
            rec.clock_offset = base + if i % 2 == 0 { dither } else { -dither };
        }
        let orb = PreciseOrbit::new(eps);
        let t = GpsTime::new(2105, 5.5 * dt);
        let (_, clk) = orb.position_at("G01", t).unwrap();
        // Nodes 5 (+dither) and 6 (-dither) bracket the query: linear interp
        // must return their exact average.
        assert!(
            (clk - base).abs() < 1e-15,
            "clock {} deviates from chord average {}",
            clk,
            base
        );
    }

    #[test]
    fn test_transmit_time_refinement_with_rx_hint() {
        // Satellite FIXED in ECEF (identical nodes) with a linearly drifting
        // clock: clock interpolation is then exact, so any deviation isolates
        // the transmit-time shift. A receiver at the origin sees the satellite
        // at range r; light travel time is r/c ~ 88.6 ms, during which this
        // clock drifts by 1e-9 * 0.0886 ~ 89 ps.
        let dt = 900.0;
        let r = 26_560_000.0_f64;
        let eps: Vec<Sp3Epoch> = (0..30)
            .map(|i| {
                let t = i as f64 * dt;
                let mut recs = HashMap::new();
                recs.insert(
                    "G01".to_string(),
                    Sp3Record {
                        position: Vector3::new(r, 0.0, 0.0),
                        clock_offset: 1.0e-6 + 1.0e-9 * t,
                    },
                );
                Sp3Epoch { time: GpsTime::new(2105, t), records: recs }
            })
            .collect();
        let orb = PreciseOrbit::new(eps);

        let t_rx = GpsTime::new(2105, 4500.0);
        let travel_s = r / SPEED_OF_LIGHT_M_S;

        // With a receiver-position hint the state must refer to transmit time
        // t_rx - range/c.
        let (pos, clk) = orb
            .position_at_with_hint("G01", t_rx, Some(Vector3::zeros()))
            .unwrap();
        let expected_clk = 1.0e-6 + 1.0e-9 * (4500.0 - travel_s);
        assert!(
            (pos - Vector3::new(r, 0.0, 0.0)).norm() < 1e-6,
            "position must be invariant under transmit-time refinement"
        );
        assert!(
            (clk - expected_clk).abs() < 1e-15,
            "clock {} not evaluated at transmit time {}",
            clk,
            expected_clk
        );

        // Without a hint the receiver-time behaviour is preserved unchanged.
        let (_, clk_nohint) = orb.position_at("G01", t_rx).unwrap();
        let expected_nohint = 1.0e-6 + 1.0e-9 * 4500.0;
        assert!((clk_nohint - expected_nohint).abs() < 1e-15);
    }

    #[test]
    fn test_position_and_velocity_relativistic_correction() {
        let dt = 900.0;
        let r = 26_560_000.0_f64;
        let eps: Vec<Sp3Epoch> = (0..30)
            .map(|i| {
                let t = i as f64 * dt;
                let mut recs = HashMap::new();
                recs.insert(
                    "G01".to_string(),
                    Sp3Record {
                        position: Vector3::new(r, 0.0, 0.0),
                        clock_offset: 1.0e-6,
                    },
                );
                Sp3Epoch { time: GpsTime::new(2105, t), records: recs }
            })
            .collect();
        let orb = PreciseOrbit::new(eps);

        let t_query = GpsTime::new(2105, 4500.0);
        let (pos, _vel, clk) = orb.position_and_velocity_at("G01", t_query).unwrap();
        assert!((pos.x - r).abs() < 1e-3);
        assert!((clk - 1.0e-6).abs() < 1e-10);
    }
}

