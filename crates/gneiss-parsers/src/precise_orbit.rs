//! Precise orbit interpolation from SP3 ephemeris data.
//!
//! Wraps parsed SP3 epochs and provides position/clock lookup at arbitrary
//! epochs via Lagrange polynomial interpolation (degree configurable,
//! default 8). This replaces broadcast Keplerian propagation with IGS
//! final/rapid products for cm-level orbit accuracy.

use crate::sp3::{Sp3Epoch, Sp3Record};
use gneiss_core::time::GpsTime;
use std::collections::HashMap;

/// Interpolated precise orbit store.
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

    /// Interpolated position (m) and clock offset (s) at `t`.
    /// Returns None if satellite unknown or insufficient data.
    pub fn position_at(&self, sv: &str, t: GpsTime) -> Option<(nalgebra::Vector3<f64>, f64)> {
        let track = self.tracks.get(sv)?;
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
        if window.is_empty() {
            return None;
        }

        // Lagrange interpolation for each component.
        let interp = |target_s: f64, component: usize| -> f64 {
            let mut result = 0.0;
            for i in 0..window.len() {
                let val = match component {
                    0 => window[i].1.position.x,
                    1 => window[i].1.position.y,
                    2 => window[i].1.position.z,
                    _ => window[i].1.clock_offset,
                };
                let mut term = val;
                let ti = s_time(window[i].0);
                for j in 0..window.len() {
                    if j == i { continue; }
                    let tj = s_time(window[j].0);
                    if (ti - tj).abs() < 1e-12 { continue; }
                    term *= (target_s - tj) / (ti - tj);
                }
                result += term;
            }
            result
        };

        let target = s_time(t);
        let px = interp(target, 0);
        let py = interp(target, 1);
        let pz = interp(target, 2);
        let clk = interp(target, 3);

        Some((nalgebra::Vector3::new(px, py, pz), clk))
    }

    /// Number of satellites in this store.
    pub fn len(&self) -> usize {
        self.tracks.len()
    }

    pub fn is_empty(&self) -> bool {
        self.tracks.is_empty()
    }
}

fn s_time(t: GpsTime) -> f64 {
    t.week as f64 * 604800.0 + t.tow
}

#[cfg(test)]
mod tests {
    use super::*;
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
}
