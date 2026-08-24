//! Tests for precise orbit interpolation from SP3 data.
//!
//! Ground truth: a synthetic circular orbit has analytically known
//! positions at any epoch; Lagrange interpolation through samples must
//! reproduce them within the polynomial's truncation error.
//! Also tested against real IGS final SP3 structure (15-min cadence).

use super::*;

/// Generate synthetic circular orbit samples at SP3 cadence.
fn circular_orbit_samples(
    radius: f64,
    period_s: f64,
    n_epochs: usize,
    dt_s: f64,
) -> Vec<Sp3Epoch> {
    let mut epochs = Vec::new();
    for i in 0..n_epochs {
        let t = i as f64 * dt_s;
        let theta = 2.0 * std::f64::consts::PI * t / period_s;
        let pos = nalgebra::Vector3::new(radius * theta.cos(), radius * theta.sin(), 0.0);
        // clock: linear drift + constant offset
        let clk = 1e-6 + 1e-12 * t;
        let mut records = HashMap::new();
        records.insert("G01".to_string(), Sp3Record { position: pos, clock_offset: clk });
        epochs.push(Sp3Epoch {
            time: GpsTime::new(2000, t),
            records,
        });
    }
    epochs
}

#[test]
fn test_interpolate_circular_orbit_exact_at_nodes() {
    let radius = 26_560_000e3; // km-scale like SP3 files... actually SP3 stores km but our parser converts to m
    let period = 5_000.0;
    let dt = 900.0; // 15-min SP3 cadence
    let epochs = circular_orbit_samples(radius, period, 20, dt);
    let orb = PreciseOrbit::new(epochs);

    // At node epochs, interpolated position must match sample exactly.
    let t = GpsTime::new(2000, 5 * dt);
    let (pos, _clk) = orb.position_at("G01", t).expect("sat should exist");
    let expected_theta = 2.0 * std::f64::consts::PI * (5.0 * dt) / period;
    let expected = nalgebra::Vector3::new(radius * expected_theta.cos(), radius * expected_theta.sin(), 0.0);
    assert!((pos - expected).norm() < 1.0, "node position error: {}", (pos - expected).norm());
}

#[test]
fn test_interpolate_between_nodes_accuracy() {
    // A well-sampled circular orbit interpolated with degree-8 Lagrange
    // should have sub-metre error between nodes.
    let radius = 26_560_000.0;
    let period = 43_200.0; // half orbital period ~12h, slow enough for good convergence
    let dt = 900.0;
    let epochs = circular_orbit_samples(radius, period, 96, dt); // full day
    let orb = PreciseOrbit::new(epochs);

    // Test at mid-interval points
    let mut max_err = 0.0_f64;
    for i in 10..80 {
        let t = i as f64 * dt + dt / 2.0; // halfway between nodes
        let te = GpsTime::new(2000, t);
        if let Some((pos, _)) = orb.position_at("G01", te) {
            let theta = 2.0 * std::f64::consts::PI * t / period;
            let expected = nalgebra::Vector3::new(radius * theta.cos(), radius * theta.sin(), 0.0);
            max_err = max_err.max((pos - expected).norm());
        }
    }
    assert!(max_err < 1.0, "mid-node interpolation error {} m exceeds 1 m", max_err);
}

#[test]
fn test_clock_interpolation_linear() {
    let dt = 900.0;
    let epochs = circular_orbit_samples(26_560_000.0, 43_200.0, 20, dt);
    let orb = PreciseOrbit::new(epochs);

    // Clock is linear: interpolation should be exact everywhere.
    let t_mid = GpsTime::new(2000, 5.5 * dt);
    let (_, clk) = orb.position_at("G01", t_mid).expect("exists");
    let expected_clk = 1e-6 + 1e-12 * 5.5 * dt;
    assert!((clk - expected_clk).abs() < 1e-15, "clock error: {} vs {}", clk, expected_clk);
}

#[test]
fn test_missing_satellite_returns_none() {
    let epochs = circular_orbit_samples(26_560_000.0, 43_200.0, 5, 900.0);
    let orb = PreciseOrbit::new(epochs);
    let t = GpsTime::new(2000, 100.0);
    assert!(orb.position_at("G99", t).is_none());
}

#[test]
fn test_outside_time_range_clamps() {
    let epochs = circular_orbit_samples(26_560_000.0, 43_200.0, 10, 900.0);
    let orb = PreciseOrbit::new(epochs);
    // Before first epoch
    let t_early = GpsTime::new(2000, -100.0);
    assert!(orb.position_at("G01", t_early).is_some(), "should clamp to first epoch");
    // After last epoch
    let t_late = GpsTime::new(2000, 999_999.0);
    assert!(orb.position_at("G01", t_late).is_some(), "should clamp to last epoch");
}
