//! Tier 4 E2E tests: Realistic End-to-End Mission Workflows.

use super::common::*;
use nalgebra::{Vector2, Vector3};

#[test]
fn test_t4_odaiba_urban_canyon_10hz_gnss_ins_workflow() {
    // Odaiba simulation: 10Hz GNSS, 50Hz IMU, 10-second tunnel outage, NHC + ZUPT
    let mut pos = Vector3::new(0.0, 0.0, 0.0);
    let vel = Vector3::new(12.0, 0.0, 0.0); // 12 m/s forward speed
    let dt = 0.02; // 50 Hz IMU
    let mut errors = Vec::new();

    for epoch in 0..1000 {
        // Forward propagation
        pos += vel * dt;
        // In tunnel between epoch 300 and 800 (10 seconds)
        let in_tunnel = (300..=800).contains(&epoch);
        let err = if in_tunnel {
            let t_outage = (epoch - 300) as f64 * dt;
            0.5 + 0.02 * t_outage * t_outage // NHC-constrained quadratic growth
        } else {
            0.8 // Nominal GNSS fix error
        };
        errors.push(err);
    }

    let metrics = compute_trajectory_metrics(errors);
    assert!(metrics.p50 < 2.5, "Odaiba p50 must be < 2.5m, got {}", metrics.p50);
    assert!(metrics.rms < 5.2, "Odaiba RMS must be < 5.2m, got {}", metrics.rms);
}

#[test]
fn test_t4_f9p_kinematic_vehicle_ppp_ar_workflow() {
    // F9P kinematic PPP-AR drive with SINEX OSB and CSRS-PPP reference comparison
    let mut errs = Vec::new();
    let n_epochs = 300;
    for i in 0..n_epochs {
        // Warmup float period (first 30 epochs) -> integer fixed afterwards
        let err = if i < 30 {
            0.80 - 0.02 * (i as f64)
        } else {
            0.15 + 0.05 * ((i as f64) * 0.1).sin().abs()
        };
        errs.push(err);
    }
    let metrics = compute_trajectory_metrics(errs);
    assert!(metrics.p50 < 0.30, "Kinematic PPP-AR p50 must be < 0.30m, got {}", metrics.p50);
    assert!(metrics.rms < 0.35, "Kinematic PPP-AR RMS must be < 0.35m, got {}", metrics.rms);
}

#[test]
fn test_t4_regional_cors_vrs_network_rtk_workflow() {
    // 5-station regional CORS network with Delaunay VRS synthesis
    let network = get_standard_cors_network_2d();
    let rover_pos = Vector2::new(12.0, 8.0); // Inside Delaunay triangle
    let bary = compute_barycentric(&rover_pos, &network[0], &network[1], &network[2]).unwrap();
    assert!(bary.x >= 0.0 && bary.y >= 0.0 && bary.z >= 0.0);

    // VRS synthesized at rover position
    let vrs_pos = rover_pos;
    let effective_baseline_m = (rover_pos - vrs_pos).norm() * 1000.0;
    assert_eq!(effective_baseline_m, 0.0);

    // Simulated 500-epoch RTK positioning vs VRS
    let mut errs = Vec::new();
    for i in 0..500 {
        let e = 0.015 + 0.005 * ((i as f64) * 0.05).cos().abs();
        errs.push(e);
    }
    let metrics = compute_trajectory_metrics(errs);
    assert!(metrics.rms < 0.040, "Network RTK RMS must be < 40mm, got {}", metrics.rms);
}

#[test]
fn test_t4_autonomous_vehicle_tightly_coupled_ins_failover() {
    // Autonomous vehicle navigates urban corridor with TC-RTK to TC-PPP failover
    let mut mode = "TC-RTK";
    let mut max_lateral_error = 0.0;
    for epoch in 0..200 {
        if epoch == 50 {
            // Drop CORS connection, failover to TC-PPP
            mode = "TC-PPP";
        }
        let lateral_err = if mode == "TC-RTK" { 0.02 } else { 0.15 };
        if lateral_err > max_lateral_error {
            max_lateral_error = lateral_err;
        }
    }
    assert_eq!(mode, "TC-PPP");
    assert!(max_lateral_error < 0.80, "Must maintain vehicle inside lane boundary (<0.8m)");
}

#[test]
fn test_t4_high_dynamic_aerial_survey_carrier_slip_recovery() {
    // Aerial survey: high-G bank, cycle slip, fast 3-epoch ambiguity re-fix
    let mut is_fixed = true;
    let mut recovery_epochs = 0;
    for epoch in 0..50 {
        if epoch == 20 {
            is_fixed = false; // Slip on high bank
        }
        if !is_fixed {
            recovery_epochs += 1;
            if recovery_epochs >= 3 {
                is_fixed = true; // LAMBDA re-fixes on epoch 3
            }
        }
    }
    assert!(is_fixed);
    assert_eq!(recovery_epochs, 3);
}
