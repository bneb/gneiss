//! Benchmark evaluation of Gneiss Standalone PPP & Kinematic PPP engine.
//!
//! Evaluates static IGS geodetic observatories (WTZR, ALIC) and moving kinematic rovers (F9P)
//! against published millimeter ground truth using CODE precise orbit (SP3) and clock (CLK) products.
//!
//! Usage: cargo run --release --bin eval_ppp

use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;
use nalgebra::Vector3;

use gneiss_core::time::GpsTime;
use gneiss_parsers::precise_orbit::PreciseOrbit;
use gneiss_parsers::rinex_clk::RinexClock;
use gneiss_rtk::post_process::{execute_post_process, PostProcessOptions};
use gneiss_rtk::swfg::config::EngineConfig;

struct PppDatasetSpec {
    name: &'static str,
    dir: &'static str,
    obs_file: &'static str,
    nav_file: &'static str,
    sp3_path: &'static str,
    clk_path: &'static str,
    bia_path: Option<&'static str>,
    truth_file: Option<&'static str>,
    static_truth: Option<Vector3<f64>>,
    max_epochs: usize,
    is_kinematic: bool,
}

fn compute_horizontal_error(pos: Vector3<f64>, truth: Vector3<f64>) -> f64 {
    let llh = gneiss_core::coords::ecef_to_llh(truth);
    let r_enu = gneiss_core::coords::ecef_to_ned_matrix(llh);
    let diff = r_enu * (pos - truth);
    (diff.x * diff.x + diff.y * diff.y).sqrt()
}

fn compute_3d_error(pos: Vector3<f64>, truth: Vector3<f64>) -> f64 {
    (pos - truth).norm()
}

fn parse_pos_truth(path: &Path) -> BTreeMap<u32, Vector3<f64>> {
    let mut truth = BTreeMap::new();
    let content = match std::fs::read_to_string(path) {
        Ok(c) => c,
        Err(_) => return truth,
    };
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('%') { continue; }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 5 {
            let date_parts: Vec<&str> = parts[0].split('/').collect();
            let time_parts: Vec<&str> = parts[1].split(':').collect();
            if date_parts.len() == 3 && time_parts.len() == 3 {
                if let (Ok(y), Ok(m), Ok(d), Ok(hr), Ok(min), Ok(sec), Ok(px), Ok(py), Ok(pz)) = (
                    date_parts[0].parse::<i32>(), date_parts[1].parse::<i32>(), date_parts[2].parse::<i32>(),
                    time_parts[0].parse::<i32>(), time_parts[1].parse::<i32>(), time_parts[2].parse::<f64>(),
                    parts[2].parse::<f64>(), parts[3].parse::<f64>(), parts[4].parse::<f64>(),
                ) {
                    let gps_time = GpsTime::from_calendar(y, m, d, hr, min, sec);
                    truth.insert(gps_time.tow.round() as u32, Vector3::new(px, py, pz));
                }
            }
        }
    }
    truth
}

fn print_stats(name: &str, mut h_errs: Vec<f64>, mut d3_errs: Vec<f64>) {
    if h_errs.is_empty() { return; }
    h_errs.sort_by(|a, b| a.total_cmp(b));
    d3_errs.sort_by(|a, b| a.total_cmp(b));
    let n = h_errs.len();
    let p10 = h_errs[(n as f64 * 0.10) as usize];
    let p25 = h_errs[(n as f64 * 0.25) as usize];
    let p50 = h_errs[n / 2];
    let p68 = h_errs[(n as f64 * 0.68) as usize];
    let p75 = h_errs[(n as f64 * 0.75) as usize];
    let p90 = h_errs[(n as f64 * 0.90) as usize];
    let p95 = h_errs[(n as f64 * 0.95) as usize];
    let p99 = h_errs[(n as f64 * 0.99) as usize];
    let max = *h_errs.last().unwrap_or(&0.0);
    let mean = h_errs.iter().sum::<f64>() / n as f64;
    let rms = (h_errs.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();

    let p50_3d = d3_errs[n / 2];
    let p95_3d = d3_errs[(n as f64 * 0.95) as usize];
    let rms_3d = (d3_errs.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();

    println!("=== {} (N={}) ===", name, n);
    println!("Horizontal Error:  p50={:.3}m,  p68={:.3}m,  p95={:.3}m,  RMS={:.3}m", p50, p68, p95, rms);
    println!("  Horizontal CDF:  p10={:.3}m, p25={:.3}m, p50={:.3}m, p68={:.3}m, p75={:.3}m, p90={:.3}m, p95={:.3}m, p99={:.3}m, max={:.3}m, mean={:.3}m",
        p10, p25, p50, p68, p75, p90, p95, p99, max, mean);
    println!("3D Position Error: p50={:.3}m,  p95={:.3}m,  RMS={:.3}m", p50_3d, p95_3d, rms_3d);
}

type PreciseProducts = (
    Option<Arc<PreciseOrbit>>,
    Option<Arc<RinexClock>>,
    Option<Arc<gneiss_parsers::sinex_bia::SinexBias>>,
);

fn load_precise_products(
    sp3_path: &str,
    clk_path: &str,
    bia_path: Option<&str>,
) -> PreciseProducts {
    let precise_orbits = if Path::new(sp3_path).exists() {
        File::open(sp3_path).ok()
            .and_then(|f| gneiss_parsers::sp3::parse_sp3(BufReader::new(f)).ok())
            .map(|sp3_epochs| Arc::new(PreciseOrbit::new(sp3_epochs)))
    } else {
        None
    };

    let precise_clocks = if Path::new(clk_path).exists() {
        std::fs::read_to_string(clk_path).ok()
            .map(|clk_str| Arc::new(RinexClock::parse(&clk_str)))
    } else {
        None
    };

    let sinex_bias = bia_path.and_then(|p| {
        if Path::new(p).exists() {
            File::open(p).ok()
                .and_then(|f| gneiss_parsers::sinex_bia::SinexBias::parse(BufReader::new(f)).ok())
                .map(Arc::new)
        } else {
            None
        }
    });

    (precise_orbits, precise_clocks, sinex_bias)
}

fn evaluate_ppp_dataset(spec: &PppDatasetSpec) {
    println!("\n========================================================");
    println!("Evaluating PPP Dataset: {}", spec.name);
    println!("========================================================");

    let dir = Path::new(spec.dir);
    let nav_f = match File::open(dir.join(spec.nav_file)) {
        Ok(f) => f,
        Err(e) => { eprintln!("Failed to open nav {}: {}", spec.nav_file, e); return; }
    };
    let (ephemerides, klobuchar) = match gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_f)) {
        Ok(r) => r,
        Err(e) => { eprintln!("Failed to parse nav: {}", e); return; }
    };

    let obs_f = match File::open(dir.join(spec.obs_file)) {
        Ok(f) => f,
        Err(e) => { eprintln!("Failed to open obs {}: {}", spec.obs_file, e); return; }
    };
    let (obs_epochs, obs_header) = match gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(obs_f)) {
        Ok(r) => r,
        Err(e) => { eprintln!("Failed to parse obs: {}", e); return; }
    };
    let (precise_orbits, precise_clocks, sinex_bias) = load_precise_products(spec.sp3_path, spec.clk_path, spec.bia_path);
    let truth_map = spec.truth_file.map(|f| parse_pos_truth(&dir.join(f))).unwrap_or_default();

    let init_pos = obs_header.approx_position
        .map(|p| Vector3::new(p[0], p[1], p[2]))
        .or(spec.static_truth)
        .unwrap_or_else(|| Vector3::new(0.0, 0.0, 0.0));

    let is_kinematic = spec.is_kinematic;
    let config = EngineConfig::Ppp(gneiss_rtk::swfg::config::PppConfig {
        initial_position: Some([init_pos.x, init_pos.y, init_pos.z]),
        window_size: 10,
        is_kinematic,
        ..Default::default()
    });

    let max_epochs = std::env::var("MAX_EPOCHS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(spec.max_epochs);
    let selected_epochs = &obs_epochs[..obs_epochs.len().min(max_epochs)];
    let antex_path = std::path::Path::new("datasets/igs/igs14.atx");
    let antex_db = if antex_path.exists() {
        gneiss_parsers::antex::AntexDatabase::parse(antex_path).ok().map(std::sync::Arc::new)
    } else {
        None
    };
    let options = PostProcessOptions {
        enable_bidirectional: false,
        base_position: None,
        initial_rover_position: Some(init_pos),
        klobuchar_alpha: klobuchar.as_ref().map(|k| k.alpha),
        klobuchar_beta: klobuchar.as_ref().map(|k| k.beta),
        q_accel: None,
        widelane_ar: false,
        tropo_gradients: false,
        network_sat_upd: None,
        receiver_pcv: None,
        dynamics: if is_kinematic {
            gneiss_rtk::post_process::dynamics::ProcessingDynamics::Kinematic
        } else {
            gneiss_rtk::post_process::dynamics::ProcessingDynamics::Static
        },
        enable_glonass: false,
        continuity_gate: false,
        precise_orbits: precise_orbits.clone(),
        precise_clocks: precise_clocks.clone(),
        sinex_bias: sinex_bias.clone(),
        antex_database: antex_db.clone(),
    };

    let result = match execute_post_process(&config, &ephemerides, selected_epochs, None, None, &options) {
        Ok(r) => r,
        Err(e) => { eprintln!("PPP post-processing failed: {}", e); return; }
    };

    let tide0 = gneiss_geodesy::tides::solid_earth_tide(init_pos, 345600.0, 2137);
    let tide599 = gneiss_geodesy::tides::solid_earth_tide(init_pos, 363570.0, 2137);
    println!("SOLID EARTH TIDE at Epoch 0:   [{:.4}, {:.4}, {:.4}], norm={:.4}m", tide0.x, tide0.y, tide0.z, tide0.norm());
    println!("SOLID EARTH TIDE at Epoch 599: [{:.4}, {:.4}, {:.4}], norm={:.4}m", tide599.x, tide599.y, tide599.z, tide599.norm());
    let antex_path = std::path::Path::new("datasets/igs/igs14.atx");
    let recv_pco = if antex_path.exists() {
        gneiss_rtk::post_process::antenna::station_recv_pco_ecef(
            &dir.join(spec.obs_file),
            antex_path.to_str().expect("Valid UTF-8 ANTEX path"),
            init_pos,
        ).unwrap_or_else(Vector3::zeros)
    } else {
        Vector3::zeros()
    };
    if recv_pco.norm() > 0.0 {
        println!("Receiver PCO (ECEF): [{:.4}, {:.4}, {:.4}], norm={:.4}m", recv_pco.x, recv_pco.y, recv_pco.z, recv_pco.norm());
        if let Some(st) = spec.static_truth {
            let apc = st + recv_pco;
            println!("Truth Monument: [{:.4}, {:.4}, {:.4}]", st.x, st.y, st.z);
            println!("Truth APC:      [{:.4}, {:.4}, {:.4}]", apc.x, apc.y, apc.z);
        }
    }
    let n_total = result.trajectory.len();
    let mut h_errs = Vec::new();
    if let Some(truth) = spec.static_truth {
        for (label, target_ep) in [("Epoch 0", obs_epochs.first()), ("Last Epoch", selected_epochs.last())] {
            if let Some(ep) = target_ep {
                if let Some(ref orbits) = precise_orbits {
                    if let Ok(raw_sats) = gneiss_rtk::swfg::engine::epoch::extract_raw_observations_with_source(
                        ep,
                        &gneiss_rtk::estimators::rtk_iekf::satpos::PreciseSrc {
                            orbits,
                            clocks: precise_clocks.as_deref(),
                        },
                    Some(truth),
                    sinex_bias.as_deref(),
                    antex_db.as_deref(),
                ) {
                    println!("--- PREFIT RESIDUALS AT TRUTH ({}) ---", label);
                    for s in &raw_sats {
                        if let Some(pr2) = s.pr_l2 {
                            let gamma = (s.f1 / s.f2).powi(2);
                            let pr_if = (gamma * s.pr_l1 - pr2) / (gamma - 1.0);
                            let range = (s.sat_pos_ecef - truth).norm();
                            let modeled = range - s.sat_clock_m + s.tropo_dry_m;
                            let res = pr_if - modeled;
                            let cp_res_str = if let (Some(cp1), Some(cp2)) = (s.cp_l1, s.cp_l2) {
                                let l1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / s.f1;
                                let l2 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / s.f2;
                                let cp_if = (gamma * cp1 * l1 - cp2 * l2) / (gamma - 1.0);
                                format!("{:8.3}m", cp_if - modeled)
                            } else {
                                "     N/A".to_string()
                            };
                            let (el_deg, az_deg) = (s.elevation_rad.to_degrees(), s.azimuth_rad.to_degrees());
                            println!("Sat c={} p={:2}: el={:4.1} az={:5.1} | pr_res={:8.3}m  cp_res={}",
                                s.constellation_id, s.satellite, el_deg, az_deg, res, cp_res_str);
                        }
                    }
                }
            }
        }
    }
    }
    let mut d3_errs = Vec::new();
    for (i, ep) in result.trajectory.iter().enumerate() {
        let tow = ep.time.tow.round() as u32;
        let truth_opt = spec.static_truth.or_else(|| truth_map.get(&tow).copied());
        if let Some(truth) = truth_opt {
            let h = compute_horizontal_error(ep.position_ecef, truth);
            let d3 = compute_3d_error(ep.position_ecef, truth);
            if i < 3 || i + 3 >= n_total {
                println!("Epoch {} (TOW={}): sol=[{:.3}, {:.3}, {:.3}], truth=[{:.3}, {:.3}, {:.3}], diff=[{:.3}, {:.3}, {:.3}], h={:.3}m, 3D={:.3}m",
                    i, tow, ep.position_ecef.x, ep.position_ecef.y, ep.position_ecef.z,
                    truth.x, truth.y, truth.z,
                    ep.position_ecef.x - truth.x, ep.position_ecef.y - truth.y, ep.position_ecef.z - truth.z,
                    h, d3);
            }
            h_errs.push(h);
            d3_errs.push(d3);
        }
    }

    print_stats(&format!("PPP Solution - {}", spec.name), h_errs, d3_errs);
}

fn main() {
    tracing_subscriber::fmt()
        .with_env_filter("warn")
        .with_target(false)
        .without_time()
        .try_init()
        .ok();

    println!("========================================================");
    println!("Gneiss Standalone & Kinematic PPP Evaluation Harness");
    println!("========================================================");

    let only = std::env::var("PPP_ONLY").ok();

    // 1. Wettzell Geodetic Observatory (WTZR, Germany) - 5 Hours
    let wtzr_spec = PppDatasetSpec {
        name: "WTZR Geodetic Observatory (Germany, 30s, Float PPP)",
        dir: "datasets/wtzr_ppp_1224",
        obs_file: "WTZR00DEU_R_20203590000_01D_30S_MO.rnx",
        nav_file: "BRDC00IGS_R_20203590000_01D_MN.rnx",
        sp3_path: "datasets/wtzr_ppp_1224/com21374.eph",
        clk_path: "datasets/wtzr_ppp_1224/com21374.clk",
        bia_path: Some("datasets/wtzr_ppp_1224/com21374.bia"),
        truth_file: None,
        static_truth: Some(Vector3::new(4075580.8863, 931853.5784, 4801567.9707)),
        max_epochs: 600,
        is_kinematic: false,
    };
    if only.as_deref().is_none_or(|w| w == "wtzr") && Path::new(wtzr_spec.dir).exists() {
        evaluate_ppp_dataset(&wtzr_spec);
    }

    // 2. Alice Springs (ALIC, Australia) - 5 Hours
    let alic_spec = PppDatasetSpec {
        name: "ALIC Geodetic Observatory (Australia, 30s, Float PPP)",
        dir: "datasets/igs",
        obs_file: "alic3350.19o",
        nav_file: "brdc3350.19n",
        sp3_path: "datasets/igs/cod20820.sp3",
        clk_path: "datasets/igs/gfz20820.clk",
        bia_path: None,
        truth_file: None,
        static_truth: Some(Vector3::new(-4052052.7533, 4212835.9866, -2545104.6062)),
        max_epochs: 600,
        is_kinematic: false,
    };
    if only.as_deref().is_none_or(|w| w == "alic") && Path::new(alic_spec.dir).exists() {
        evaluate_ppp_dataset(&alic_spec);
    }

    // 3. RTK Explorer F9P Kinematic Vehicle Drive (1Hz, Float PPP)
    let f9p_spec = PppDatasetSpec {
        name: "RTK Explorer F9P Kinematic Vehicle Drive (1Hz, Float PPP)",
        dir: "datasets/rtkexplorer/sample_1/f9p_ppp_1224",
        obs_file: "rover.obs",
        nav_file: "BRDC00IGS_R_20203590000_01D_MN.rnx",
        sp3_path: "datasets/rtkexplorer/sample_1/f9p_ppp_1224/ESA0MGNFIN_20203590000_01D_05M_ORB.SP3",
        clk_path: "datasets/rtkexplorer/sample_1/f9p_ppp_1224/ESA0MGNFIN_20203590000_01D_30S_CLK.CLK",
        bia_path: Some("datasets/rtkexplorer/sample_1/f9p_ppp_1224/COD0MGXFIN_20203590000_01D_01D_OSB.BIA"),
        truth_file: Some("rover_ppk.pos"),
        static_truth: None,
        max_epochs: 600,
        is_kinematic: true,
    };
    if only.as_deref().is_none_or(|w| w == "f9p") && Path::new(f9p_spec.dir).exists() {
        evaluate_ppp_dataset(&f9p_spec);
    }
}
