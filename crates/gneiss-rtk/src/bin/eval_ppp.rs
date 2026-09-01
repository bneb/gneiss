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

    let init_pos = spec.static_truth.or_else(|| {
        obs_header.approx_position.map(|p| Vector3::new(p[0], p[1], p[2]))
    }).unwrap_or_else(|| Vector3::new(0.0, 0.0, 0.0));

    let config = EngineConfig::Ppp(gneiss_rtk::swfg::config::PppConfig {
        initial_position: Some([init_pos.x, init_pos.y, init_pos.z]),
        window_size: 10,
        ..Default::default()
    });

    let selected_epochs = &obs_epochs[..obs_epochs.len().min(spec.max_epochs)];
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
        dynamics: Default::default(),
        enable_glonass: false,
        continuity_gate: false,
        precise_orbits,
        precise_clocks,
        sinex_bias,
    };

    let result = match execute_post_process(&config, &ephemerides, selected_epochs, None, None, &options) {
        Ok(r) => r,
        Err(e) => { eprintln!("PPP post-processing failed: {}", e); return; }
    };

    let mut h_errs = Vec::new();
    let mut d3_errs = Vec::new();
    for ep in &result.trajectory {
        let tow = ep.time.tow.round() as u32;
        let truth_opt = spec.static_truth.or_else(|| truth_map.get(&tow).copied());
        if let Some(truth) = truth_opt {
            h_errs.push(compute_horizontal_error(ep.position_ecef, truth));
            d3_errs.push(compute_3d_error(ep.position_ecef, truth));
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
    };
    if Path::new(wtzr_spec.dir).exists() {
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
    };
    if Path::new(alic_spec.dir).exists() {
        evaluate_ppp_dataset(&alic_spec);
    }

    // 3. RTK Explorer F9P Kinematic Vehicle Drive (1Hz, Float PPP)
    let f9p_spec = PppDatasetSpec {
        name: "RTK Explorer F9P Kinematic Vehicle Drive (1Hz, Float PPP)",
        dir: "datasets/rtkexplorer/sample_1/f9p_ppp_1224",
        obs_file: "rover.obs",
        nav_file: "BRDC00IGS_R_20203590000_01D_MN.rnx",
        sp3_path: "datasets/rtkexplorer/sample_1/f9p_ppp_1224/COD0MGXFIN_20203590000_01D_05M_ORB.SP3",
        clk_path: "datasets/rtkexplorer/sample_1/f9p_ppp_1224/COD0MGXFIN_20203590000_01D_30S_CLK.CLK",
        bia_path: Some("datasets/rtkexplorer/sample_1/f9p_ppp_1224/COD0MGXFIN_20203590000_01D_01D_OSB.BIA"),
        truth_file: Some("rover_ppk.pos"),
        static_truth: None,
        max_epochs: 600,
    };
    if Path::new(f9p_spec.dir).exists() {
        evaluate_ppp_dataset(&f9p_spec);
    }
}
