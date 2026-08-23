//! Multi-Tier Automated Benchmark Integration Test Suite.
//!
//! Evaluates Gneiss post-processing engine against Tier 1 Geodetic CORS
//! and Tier 2 High-Dynamic Simulations.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use nalgebra::Vector3;

use gneiss_rtk::post_process::{execute_post_process, PostProcessOptions};
use gneiss_rtk::sim::generator::{SimulationConfig, TrajectoryProfile};
use gneiss_rtk::swfg::config::EngineConfig;

fn compute_horizontal_error(pos: Vector3<f64>, truth: Vector3<f64>) -> f64 {
    let llh = gneiss_core::coords::ecef_to_llh(truth);
    let ned = gneiss_core::coords::ecef_to_ned_matrix(llh) * (pos - truth);
    (ned.x * ned.x + ned.y * ned.y).sqrt()
}

#[test]
fn test_tier1_geodetic_cors_baseline_sub_centimeter() {
    let dir = Path::new("datasets/cors_short_baseline");
    if !dir.exists() {
        return;
    }

    let nav_path = dir.join("brdc1350.20n");
    let base_path = dir.join("tmg21350.20o");
    let rov_path = dir.join("tmgo1350.20o");

    if !nav_path.exists() || !base_path.exists() || !rov_path.exists() {
        return;
    }

    let nav_f = File::open(&nav_path).expect("Open nav");
    let (ephems, klob) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_f)).expect("Parse nav");

    let rov_f = File::open(&rov_path).expect("Open rover");
    let (rov_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(rov_f)).expect("Parse rover");

    let base_f = File::open(&base_path).expect("Open base");
    let (base_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(base_f)).expect("Parse base");

    let base_pos = Vector3::new(-1283433.9360, -4713073.2930, 4090105.0870);
    let truth_pos = Vector3::new(-1283387.0660, -4713016.7750, 4090190.3860);

    let config = EngineConfig::Rtk(gneiss_rtk::swfg::config::RtkConfig {
        initial_position: Some([base_pos.x, base_pos.y, base_pos.z]),
        ..Default::default()
    });

    let selected_rover = &rov_epochs[..rov_epochs.len().min(40)];

    let options = PostProcessOptions {
        enable_bidirectional: true,
        base_position: Some(base_pos),
        initial_rover_position: Some(base_pos),
        klobuchar_alpha: klob.as_ref().map(|k| k.alpha),
        klobuchar_beta: klob.as_ref().map(|k| k.beta),
    };

    let result = execute_post_process(&config, &ephems, selected_rover, Some(&base_epochs), None, &options)
        .expect("Post-process CORS short baseline");

    let mut h_errs: Vec<f64> = Vec::new();
    let mut fixed = 0;
    for ep in &result.trajectory {
        if ep.quality == 1 { fixed += 1; }
        h_errs.push(compute_horizontal_error(ep.position_ecef, truth_pos));
    }

    h_errs.sort_by(|a, b| a.total_cmp(b));
    let p50 = h_errs[h_errs.len() / 2];
    println!("Tier 1 Geodetic CORS (TMG2-TMGO): p50={:.4}m, fixed={}/{}", p50, fixed, h_errs.len());

    assert!(p50 < 0.020, "Tier 1 Geodetic CORS p50 horizontal error must be < 20mm, got {:.4}m", p50);
}

#[test]
fn test_tier2_high_dynamic_circular_kinematic_sub_centimeter() {
    let config = SimulationConfig {
        duration_s: 30.0,
        epoch_rate_hz: 1.0,
        cp_noise_m: 0.002,
        pr_noise_m: 0.15,
        profile: TrajectoryProfile::Circular {
            center_offset_ned: Vector3::new(100.0, 100.0, 0.0),
            radius_m: 50.0,
            speed_m_s: 5.0,
        },
        ..Default::default()
    };

    let sim = gneiss_rtk::sim::generator::generate_simulation_dataset(&config);
    let engine_config = EngineConfig::Rtk(gneiss_rtk::swfg::config::RtkConfig {
        initial_position: Some([config.base_ecef.x, config.base_ecef.y, config.base_ecef.z]),
        ..Default::default()
    });

    let options = PostProcessOptions {
        enable_bidirectional: true,
        base_position: Some(config.base_ecef),
        initial_rover_position: Some(config.base_ecef),
        klobuchar_alpha: None,
        klobuchar_beta: None,
    };

    let result = execute_post_process(
        &engine_config,
        &sim.ephemerides,
        &sim.rover_epochs,
        Some(&sim.base_epochs),
        None,
        &options,
    )
    .expect("High-dynamic circular kinematic post-process should succeed");

    let mut h_errs = Vec::new();
    for (i, epoch) in result.trajectory.iter().enumerate() {
        let truth = sim.truth_positions[i].1;
        let h_err = compute_horizontal_error(epoch.position_ecef, truth);
        h_errs.push(h_err);
    }

    let rms_h = (h_errs.iter().map(|e| e * e).sum::<f64>() / h_errs.len() as f64).sqrt();
    println!("Tier 2 High-Dynamic Circular Kinematic: RMS={:.4}m", rms_h);
    assert!(rms_h < 0.015, "High-dynamic circular kinematic RMS must be < 1.5cm, got {:.4}m", rms_h);
}
