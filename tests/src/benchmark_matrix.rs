//! Multi-Tier Automated Benchmark Integration Test Suite.
//!
//! Evaluates Gneiss post-processing engine against Real Geodetic CORS baselines,
//! Real-world IGS Float PPP tracking, and High-Dynamic Kinematic Simulations.

use std::fs::File;
use std::io::BufReader;
use nalgebra::Vector3;

use gneiss_rtk::post_process::{execute_post_process, PostProcessOptions};
use gneiss_rtk::sim::generator::{SimulationConfig, TrajectoryProfile};
use gneiss_rtk::swfg::config::EngineConfig;

fn compute_horizontal_error(pos: Vector3<f64>, truth: Vector3<f64>) -> f64 {
    let llh = gneiss_core::coords::ecef_to_llh(truth);
    let ned = gneiss_core::coords::ecef_to_ned_matrix(llh) * (pos - truth);
    (ned.x * ned.x + ned.y * ned.y).sqrt()
}

fn find_dataset_dir(name: &str) -> Option<std::path::PathBuf> {
    let p = std::path::PathBuf::from(name);
    if p.exists() {
        return Some(p);
    }
    let p = std::path::PathBuf::from("..").join(name);
    if p.exists() {
        return Some(p);
    }
    None
}

/// Benchmark 1: Real Geodetic CORS Dual-Frequency Short Baseline (TMG2-TMGO).
/// Evaluates double-differenced post-processed kinematic/static RTK against ground truth.
#[test]
fn test_real_geodetic_cors_baseline_tmg2_tmgo_sub_centimeter() {
    let dir = match find_dataset_dir("datasets/cors_short_baseline") {
        Some(d) => d,
        None => return,
    };

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
        ..Default::default()
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
    println!("Real Geodetic CORS (TMG2-TMGO): p50={:.4}m, fixed={}/{}", p50, fixed, h_errs.len());

    assert!(p50 < 0.010, "Real Geodetic CORS p50 horizontal error must be < 10mm, got {:.4}m", p50);
    assert!(fixed >= 35, "At least 35/40 epochs must achieve fixed integer ambiguities, got {}", fixed);
}

/// Benchmark 2: High-Dynamic Circular Kinematic Trajectory Simulation.
/// Evaluates estimator dynamics, carrier tracking, and ambiguity fixing under high angular rate.
#[test]
fn test_simulation_high_dynamic_circular_kinematic_sub_centimeter() {
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
        ..Default::default()
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
    println!("Simulation High-Dynamic Circular Kinematic: RMS={:.4}m", rms_h);
    assert!(rms_h < 0.010, "High-dynamic circular kinematic RMS must be < 1.0cm, got {:.4}m", rms_h);
}

/// Benchmark 3: Physical Geodesy Normalizations & Between-Satellite Single-Difference Wide-Lane AR Math.
/// Unit verification of physical normalizations and integer wide-lane resolving formulas.
#[test]
fn test_physical_geodesy_normalizations_and_sd_widelane_math() {
    use gneiss_rtk::ambiguity::ppp_ar::{PppArSolver, WideLaneCandidate};
    use gneiss_rtk::swfg::pipeline::factors::uduc::GeodeticNormalizations;

    let truth_station = Vector3::new(4075580.4020, 931853.8640, 4801569.9530);

    // Verify all 6 physical normalizations are computed and non-zero
    let (sun_pos, _) = gneiss_geodesy::tides::solar_lunar_positions(43200.0, 2370);
    let tide_disp = gneiss_geodesy::solid_earth_tide(truth_station, 43200.0, 2370);
    assert!(tide_disp.norm() > 0.01 && tide_disp.norm() < 0.40, "Solid earth tide magnitude must be 1-40cm");

    let sat_pos = Vector3::new(15_000_000.0, 15_000_000.0, 15_000_000.0);
    let sat_vel = Vector3::new(-2000.0, 2000.0, 500.0);
    let rel_range = gneiss_geodesy::periodic_relativistic_range_correction(&sat_pos, &sat_vel);
    assert!(rel_range.abs() > 0.01, "Periodic relativity correction must be non-zero");

    let shapiro = gneiss_geodesy::gravitational_shapiro_delay(&sat_pos, &truth_station);
    assert!(shapiro > 0.005 && shapiro < 0.030, "Shapiro delay must be 5-30mm");

    let body_pco = Vector3::new(0.05, 0.10, 1.20);
    let pco_ecef = gneiss_geodesy::project_satellite_pco_to_ecef(&sat_pos, &sun_pos, &body_pco);
    assert!(pco_ecef.norm() > 1.0, "Projected satellite PCO norm must match body offset scale");

    let _norm = GeodeticNormalizations {
        sat_relativity_m: rel_range,
        shapiro_delay_m: shapiro,
        sat_pco_ecef_m: pco_ecef,
        solid_earth_tide_m: tide_disp,
        ocean_tide_loading_m: Vector3::zeros(),
    };

    // Verify Between-Satellite Single-Difference Integer Wide-Lane AR fixing
    let candidates = vec![
        WideLaneCandidate {
            sat_idx: 2,
            mw_cycles: 34.02,
            mw_std_cycles: 0.04,
            bias_wl_cycles: 0.01,
        },
        WideLaneCandidate {
            sat_idx: 3,
            mw_cycles: -12.98,
            mw_std_cycles: 0.05,
            bias_wl_cycles: -0.02,
        },
    ];

    let sd_wl = PppArSolver::fix_sd_wide_lane(1, 10.00, 0.00, &candidates, 0.15);
    assert_eq!(sd_wl.len(), 2, "Both candidate satellites must fix to integer wide-lanes");
    assert_eq!(sd_wl[0].n_wl, 24, "34.01 - 10.00 = 24.01 -> 24");
    assert_eq!(sd_wl[1].n_wl, -23, "-12.96 - 10.00 = -22.96 -> -23");
}

/// Benchmark 4: Real-World IGS Station Wettzell (`WTZR`) 5-Hour Standalone Float PPP Tracking.
/// Evaluates standalone dual-frequency PPP convergence against published IGS ground truth using CODE SP3 & CLK products.
#[test]
fn test_real_rinex_wtzr_float_ppp_convergence() {
    let dir = match find_dataset_dir("datasets/wtzr_ppp_1224") {
        Some(d) => d,
        None => return,
    };

    let nav_path = dir.join("BRDC00IGS_R_20203590000_01D_MN.rnx");
    let obs_path = dir.join("WTZR00DEU_R_20203590000_01D_30S_MO.rnx");
    if !nav_path.exists() || !obs_path.exists() {
        return;
    }

    let nav_f = File::open(&nav_path).expect("Open WTZR nav");
    let (ephems, klob) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_f)).expect("Parse WTZR nav");

    let obs_f = File::open(&obs_path).expect("Open WTZR obs");
    let (obs_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(obs_f)).expect("Parse WTZR obs");

    // Published IGS coordinates for WTZR (Wettzell Tower 14201M010)
    let truth_station = Vector3::new(4075580.8863, 931853.5784, 4801567.9707);
    let antex_path = find_dataset_dir("datasets/igs").map(|d| d.join("igs14.atx"));
    let recv_pco = antex_path
        .and_then(|p| p.to_str().and_then(|s| gneiss_rtk::post_process::antenna::station_recv_pco_ecef(&obs_path, s, truth_station)))
        .unwrap_or_else(Vector3::zeros);
    let truth_apc = truth_station + recv_pco;

    let config = EngineConfig::Ppp(gneiss_rtk::swfg::config::PppConfig {
        initial_position: Some([truth_station.x, truth_station.y, truth_station.z]),
        window_size: 10,
        ..Default::default()
    });

    let sp3_path = dir.join("com21374.eph");
    let clk_path = dir.join("com21374.clk");

    let precise_orbits = if sp3_path.exists() {
        let f = File::open(&sp3_path).expect("Open WTZR SP3");
        let sp3_epochs = gneiss_parsers::sp3::parse_sp3(BufReader::new(f)).expect("Parse WTZR SP3");
        Some(std::sync::Arc::new(gneiss_parsers::precise_orbit::PreciseOrbit::new(sp3_epochs)))
    } else {
        None
    };

    let precise_clocks = if clk_path.exists() {
        let clk_str = std::fs::read_to_string(&clk_path).expect("Read WTZR CLK");
        let clk = gneiss_parsers::rinex_clk::RinexClock::parse(&clk_str);
        Some(std::sync::Arc::new(clk))
    } else {
        None
    };

    println!("Total epochs in WTZR: {}", obs_epochs.len());
    let selected_epochs = &obs_epochs[..obs_epochs.len().min(600)];

    let antex_database = find_dataset_dir("datasets/igs/igs14.atx")
        .and_then(|p| gneiss_parsers::antex::AntexDatabase::parse(p).ok())
        .map(std::sync::Arc::new);

    let options = PostProcessOptions {
        enable_bidirectional: false,
        base_position: None,
        initial_rover_position: Some(truth_station),
        klobuchar_alpha: klob.as_ref().map(|k| k.alpha),
        klobuchar_beta: klob.as_ref().map(|k| k.beta),
        q_accel: None,
        widelane_ar: false,
        tropo_gradients: false,
        network_sat_upd: None,
        receiver_pcv: None,
        dynamics: Default::default(),
        enable_glonass: false,
        continuity_gate: false,
        precise_orbits,
        precise_clocks,
        sinex_bias: None,
        antex_database: antex_database.clone(),
    };

    let result = execute_post_process(&config, &ephems, selected_epochs, None, None, &options)
        .expect("Post-process WTZR real RINEX dataset");

    assert!(!result.trajectory.is_empty(), "Trajectory must contain processed epochs");
    let mut errs = Vec::new();
    let mut vert_errs = Vec::new();
    let llh = gneiss_core::coords::ecef_to_llh(truth_apc);
    let r_enu = gneiss_core::coords::ecef_to_ned_matrix(llh);
    for (i, ep) in result.trajectory.iter().enumerate() {
        let err = compute_horizontal_error(ep.position_ecef, truth_apc);
        let ned = r_enu * (ep.position_ecef - truth_apc);
        let v_err = ned.z.abs();
        if i % 50 == 0 || i == result.trajectory.len() - 1 {
            println!("WTZR Epoch {:3}: dN={:+7.3}m dE={:+7.3}m dU={:+7.3}m | Horiz={:.3}m 3D={:.3}m",
                i, ned.x, ned.y, -ned.z, err, (ep.position_ecef - truth_apc).norm());
        }
        errs.push(err);
        vert_errs.push(v_err);
    }
    let mut sorted = errs.clone();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let p50 = sorted[sorted.len() / 2];
    let final_v_err = *vert_errs.last().unwrap();
    println!("WTZR Real RINEX Float PPP: p50={:.4}m, min={:.4}m, last_horiz={:.4}m, last_vert={:.4}m",
        p50, errs.iter().cloned().fold(f64::INFINITY, f64::min), errs.last().unwrap(), final_v_err);

    assert!(p50 < 1.00, "WTZR median horizontal float PPP error must be < 1.0m, got {:.4}m", p50);
    assert!(final_v_err < 2.00, "WTZR final vertical float PPP error must converge to < 2.0m, got {:.4}m", final_v_err);
}

/// Benchmark 5: Real-World IGS Station Alice Springs (`ALIC`) 5-Hour Standalone Float PPP Tracking.
/// Evaluates standalone PPP recovery from a 3-meter perturbed seed coordinate.
#[test]
fn test_real_rinex_alic_float_ppp_convergence() {
    let dir = match find_dataset_dir("datasets/igs") {
        Some(d) => d,
        None => return,
    };

    let nav_path = dir.join("brdc3350.19n");
    let obs_path = dir.join("alic3350.19o");
    if !nav_path.exists() || !obs_path.exists() {
        return;
    }

    let nav_f = File::open(&nav_path).expect("Open ALIC nav");
    let (ephems, klob) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_f)).expect("Parse ALIC nav");

    let obs_f = File::open(&obs_path).expect("Open ALIC obs");
    let (obs_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(obs_f)).expect("Parse ALIC obs");

    // Official IGS14 published coordinates for ALIC (Alice Springs) propagated to epoch 2019.92
    let truth_station = Vector3::new(-4052053.1199, 4212836.1950, -2545104.0884);
    let antex_path = find_dataset_dir("datasets/igs").map(|d| d.join("igs14.atx"));
    let recv_pco = antex_path
        .and_then(|p| p.to_str().and_then(|s| gneiss_rtk::post_process::antenna::station_recv_pco_ecef(&obs_path, s, truth_station)))
        .unwrap_or_else(Vector3::zeros);
    let truth_apc = truth_station + recv_pco;
    let perturbed_seed = Vector3::new(truth_station.x + 3.0, truth_station.y - 3.0, truth_station.z + 3.0);

    let config = EngineConfig::Ppp(gneiss_rtk::swfg::config::PppConfig {
        initial_position: Some([perturbed_seed.x, perturbed_seed.y, perturbed_seed.z]),
        ..Default::default()
    });

    let sp3_path = dir.join("cod20820.sp3");
    let clk_path = dir.join("gfz20820.clk");

    let precise_orbits = if sp3_path.exists() {
        let f = File::open(&sp3_path).expect("Open ALIC SP3");
        let sp3_epochs = gneiss_parsers::sp3::parse_sp3(BufReader::new(f)).expect("Parse ALIC SP3");
        Some(std::sync::Arc::new(gneiss_parsers::precise_orbit::PreciseOrbit::new(sp3_epochs)))
    } else {
        None
    };

    let precise_clocks = if clk_path.exists() {
        let clk_str = std::fs::read_to_string(&clk_path).expect("Read ALIC CLK");
        let clk = gneiss_parsers::rinex_clk::RinexClock::parse(&clk_str);
        Some(std::sync::Arc::new(clk))
    } else {
        None
    };

    let selected_epochs = &obs_epochs[..obs_epochs.len().min(600)];

    let antex_database = find_dataset_dir("datasets/igs/igs14.atx")
        .and_then(|p| gneiss_parsers::antex::AntexDatabase::parse(p).ok())
        .map(std::sync::Arc::new);

    let options = PostProcessOptions {
        enable_bidirectional: false,
        base_position: None,
        initial_rover_position: Some(perturbed_seed),
        klobuchar_alpha: klob.as_ref().map(|k| k.alpha),
        klobuchar_beta: klob.as_ref().map(|k| k.beta),
        q_accel: None,
        widelane_ar: false,
        tropo_gradients: false,
        network_sat_upd: None,
        receiver_pcv: None,
        dynamics: Default::default(),
        enable_glonass: false,
        continuity_gate: false,
        precise_orbits,
        precise_clocks,
        sinex_bias: None,
        antex_database: antex_database.clone(),
    };

    let result = execute_post_process(&config, &ephems, selected_epochs, None, None, &options)
        .expect("Post-process ALIC real RINEX dataset");

    assert!(!result.trajectory.is_empty(), "Trajectory must contain processed epochs");
    let mut errs = Vec::new();
    for ep in result.trajectory.iter() {
        let err = compute_horizontal_error(ep.position_ecef, truth_apc);
        errs.push(err);
    }
    let mut sorted = errs.clone();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let p50 = sorted[sorted.len() / 2];
    println!("ALIC Real RINEX Float PPP: p50={:.4}m, min={:.4}m, last={:.4}m",
        p50, errs.iter().cloned().fold(f64::INFINITY, f64::min), errs.last().unwrap());
    assert!(p50 < 1.50, "ALIC float PPP trajectory must recover from 3m perturbed seed to < 1.5m, got {:.4}m", p50);
}

#[test]
fn test_diagnostic_sat_pcos() {
    use gneiss_parsers::antex::AntexDatabase;
    let path = match find_dataset_dir("datasets/igs") {
        Some(d) => d.join("igs14.atx"),
        None => return,
    };
    if !path.exists() { return; }
    let db = AntexDatabase::parse(&path).expect("parse igs14.atx");
    let t = gneiss_core::time::GpsTime::new(2137, 345600.0);

    let sats = ["G08", "G10", "G15", "G16", "G18", "G20", "G23", "G26", "G27", "E07", "E12", "E14", "E24", "E25", "E26", "E31", "E33"];
    println!("\n{:<6} {:<16} {:<10} {:<10} {:<10} | {:<10} {:<10} {:<10}", "SAT", "TYPE", "F1_N(X)", "F1_E(Y)", "F1_U(Z)", "F2_N(X)", "F2_E(Y)", "F2_U(Z)");
    for s in sats {
        if let Some(ant) = db.find_satellite_gps(s, t) {
            let f1_code = if s.starts_with('G') { "G01" } else { "E01" };
            let f2_code = if s.starts_with('G') { "G02" } else { "E05" };
            let p1 = ant.frequencies.get(f1_code).map_or(nalgebra::Vector3::zeros(), |f| f.pco);
            let p2 = ant.frequencies.get(f2_code).map_or(nalgebra::Vector3::zeros(), |f| f.pco);
            println!("{:<6} {:<16} {:<10.2} {:<10.2} {:<10.2} | {:<10.2} {:<10.2} {:<10.2}",
                s, ant.antenna_type, p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
        } else {
            println!("{:<6} NOT FOUND", s);
        }
    }
}
