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
    csrs_file: Option<&'static str>,
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

fn load_csrs_truth(path: &Path) -> BTreeMap<u32, Vector3<f64>> {
    let mut map = BTreeMap::new();
    if let Ok(f) = File::open(path) {
        if let Ok(epochs) = gneiss_parsers::csrs_pos::parse_csrs_pos(BufReader::new(f)) {
            for ep in epochs {
                map.insert(ep.time.tow.round() as u32, ep.pos_ecef);
            }
        }
    }
    map
}

#[derive(Default)]
struct CsrsComparisonStats {
    g_vs_csrs_h: Vec<f64>,
    g_vs_csrs_3d: Vec<f64>,
    csrs_vs_t_h: Vec<f64>,
    csrs_vs_t_3d: Vec<f64>,
    g_vs_t_h: Vec<f64>,
    g_vs_t_3d: Vec<f64>,
}

fn record_csrs_epoch(
    ep: &gneiss_rtk::post_process::SmoothedEpoch,
    c_pos: Vector3<f64>,
    truth_map: &BTreeMap<u32, Vector3<f64>>,
    stats: &mut CsrsComparisonStats,
) {
    let llh = gneiss_core::coords::ecef_to_llh(c_pos);
    let r_ned = gneiss_core::coords::ecef_to_ned_matrix(llh);
    let ned = r_ned * (ep.position_ecef - c_pos);
    let (north, east, up) = (ned.x, ned.y, -ned.z);
    let tow = ep.time.tow.round() as u32;
    let h = compute_horizontal_error(ep.position_ecef, c_pos);
    let d3 = compute_3d_error(ep.position_ecef, c_pos);
    println!("Epoch {} (TOW={}): E={:+.4}m, N={:+.4}m, U={:+.4}m, H={:.4}m, 3D={:.4}m",
        ep.time.tow, tow, east, north, up, h, d3);
    stats.g_vs_csrs_h.push(h);
    stats.g_vs_csrs_3d.push(d3);
    if let Some(&t_pos) = truth_map.get(&tow) {
        stats.csrs_vs_t_h.push(compute_horizontal_error(c_pos, t_pos));
        stats.csrs_vs_t_3d.push(compute_3d_error(c_pos, t_pos));
        stats.g_vs_t_h.push(compute_horizontal_error(ep.position_ecef, t_pos));
        stats.g_vs_t_3d.push(compute_3d_error(ep.position_ecef, t_pos));
    }
}

fn evaluate_csrs_comparison(
    traj: &[gneiss_rtk::post_process::SmoothedEpoch],
    csrs_map: &BTreeMap<u32, Vector3<f64>>,
    truth_map: &BTreeMap<u32, Vector3<f64>>,
) {
    if csrs_map.is_empty() { return; }
    let mut stats = CsrsComparisonStats::default();
    for ep in traj {
        let tow = ep.time.tow.round() as u32;
        if let Some(&c_pos) = csrs_map.get(&tow) {
            record_csrs_epoch(ep, c_pos, truth_map, &mut stats);
        }
    }
    println!("\n=== Commercial Tier-1 Benchmark vs. CSRS-PPP (Canada Geodetic Service) ===");
    print_stats("Discrepancy: Gneiss PPP vs CSRS-PPP", stats.g_vs_csrs_h, stats.g_vs_csrs_3d);
    if !stats.csrs_vs_t_h.is_empty() {
        print_stats("Commercial Baseline: CSRS-PPP vs RTK Truth", stats.csrs_vs_t_h, stats.csrs_vs_t_3d);
        print_stats("Gneiss PPP vs RTK Truth (Matched Subset)", stats.g_vs_t_h, stats.g_vs_t_3d);
    }
}

fn print_stats(name: &str, mut h_errs: Vec<f64>, mut d3_errs: Vec<f64>) {
    if h_errs.is_empty() { return; }
    h_errs.sort_by(|a, b| a.total_cmp(b));
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

    println!("=== {} (N={}) ===", name, n);
    println!("Horizontal Error:  p50={:.3}m,  p68={:.3}m,  p95={:.3}m,  RMS={:.3}m", p50, p68, p95, rms);
    println!("  Horizontal CDF:  p10={:.3}m, p25={:.3}m, p50={:.3}m, p68={:.3}m, p75={:.3}m, p90={:.3}m, p95={:.3}m, p99={:.3}m, max={:.3}m, mean={:.3}m",
        p10, p25, p50, p68, p75, p90, p95, p99, max, mean);

    if !d3_errs.is_empty() {
        d3_errs.sort_by(|a, b| a.total_cmp(b));
        let m = d3_errs.len();
        let p50_3d = d3_errs[m / 2];
        let p95_3d = d3_errs[(m as f64 * 0.95) as usize];
        let rms_3d = (d3_errs.iter().map(|e| e * e).sum::<f64>() / m as f64).sqrt();
        println!("3D Position Error: p50={:.3}m,  p95={:.3}m,  RMS={:.3}m", p50_3d, p95_3d, rms_3d);
    }
}

type PreciseProducts = (
    Option<Arc<PreciseOrbit>>,
    Option<Arc<RinexClock>>,
    Option<Arc<gneiss_parsers::sinex_bia::SinexBias>>,
);

fn load_sp3_orbit(path: &str) -> Option<Arc<PreciseOrbit>> {
    if !Path::new(path).exists() { return None; }
    File::open(path).ok()
        .and_then(|f| gneiss_parsers::sp3::parse_sp3(BufReader::new(f)).ok())
        .map(|epochs| Arc::new(PreciseOrbit::new(epochs)))
}

fn load_rinex_clock(path: &str) -> Option<Arc<RinexClock>> {
    if !Path::new(path).exists() { return None; }
    std::fs::read_to_string(path).ok().map(|s| Arc::new(RinexClock::parse(&s)))
}

fn load_sinex_bias(path: Option<&str>) -> Option<Arc<gneiss_parsers::sinex_bia::SinexBias>> {
    let p = path?;
    if !Path::new(p).exists() { return None; }
    File::open(p).ok()
        .and_then(|f| gneiss_parsers::sinex_bia::SinexBias::parse(BufReader::new(f)).ok())
        .map(Arc::new)
}

fn load_precise_products(
    sp3_path: &str,
    clk_path: &str,
    bia_path: Option<&str>,
) -> PreciseProducts {
    (load_sp3_orbit(sp3_path), load_rinex_clock(clk_path), load_sinex_bias(bia_path))
}

type RinexNavData = (
    Vec<gneiss_core::ephemeris::Ephemeris>,
    Option<gneiss_core::atmosphere::KlobucharParams>,
);

fn load_rinex_inputs(spec: &PppDatasetSpec) -> Option<(RinexNavData, Vec<gneiss_core::obs::EpochObs>, gneiss_parsers::rinex::RinexObsHeader)> {
    let dir = Path::new(spec.dir);
    let nav_f = File::open(dir.join(spec.nav_file)).ok()?;
    let nav_data = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_f)).ok()?;
    let obs_f = File::open(dir.join(spec.obs_file)).ok()?;
    let (obs_epochs, obs_header) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(obs_f)).ok()?;
    Some((nav_data, obs_epochs, obs_header))
}

fn determine_initial_position(
    spec: &PppDatasetSpec,
    obs_header: &gneiss_parsers::rinex::RinexObsHeader,
    obs_epochs: &[gneiss_core::obs::EpochObs],
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    truth_map: &BTreeMap<u32, Vector3<f64>>,
    csrs_map: &BTreeMap<u32, Vector3<f64>>,
) -> Vector3<f64> {
    spec.static_truth
        .or_else(|| truth_map.values().next().copied())
        .or_else(|| csrs_map.values().next().copied())
        .or_else(|| obs_epochs.iter().find_map(|e| gneiss_rtk::swfg::engine::epoch::compute_spp_seeding(e, ephemerides)))
        .or_else(|| obs_header.approx_position.map(|p| Vector3::new(p[0], p[1], p[2])))
        .unwrap_or_else(|| Vector3::new(0.0, 0.0, 0.0))
}

fn build_post_process_options(
    spec: &PppDatasetSpec,
    init_pos: Vector3<f64>,
    klobuchar: Option<&gneiss_core::atmosphere::KlobucharParams>,
    products: &PreciseProducts,
    antex_db: Option<Arc<gneiss_parsers::antex::AntexDatabase>>,
) -> PostProcessOptions {
    let dynamics = if spec.is_kinematic {
        gneiss_rtk::post_process::dynamics::ProcessingDynamics::Kinematic
    } else {
        gneiss_rtk::post_process::dynamics::ProcessingDynamics::Static
    };
    PostProcessOptions {
        enable_bidirectional: false,
        base_position: None,
        initial_rover_position: Some(init_pos),
        klobuchar_alpha: klobuchar.map(|k| k.alpha),
        klobuchar_beta: klobuchar.map(|k| k.beta),
        q_accel: None,
        widelane_ar: true,
        tropo_gradients: false,
        network_sat_upd: None,
        receiver_pcv: None,
        dynamics,
        enable_glonass: false,
        continuity_gate: false,
        precise_orbits: products.0.clone(),
        precise_clocks: products.1.clone(),
        sinex_bias: products.2.clone(),
        antex_database: antex_db,
        calibration: None,
    }
}

fn log_tides_and_pco(spec: &PppDatasetSpec, dir: &Path, init_pos: Vector3<f64>) {
    let tide0 = gneiss_geodesy::tides::solid_earth_tide(init_pos, 345600.0, 2137);
    let tide599 = gneiss_geodesy::tides::solid_earth_tide(init_pos, 363570.0, 2137);
    println!("SOLID EARTH TIDE at Epoch 0:   [{:.4}, {:.4}, {:.4}], norm={:.4}m", tide0.x, tide0.y, tide0.z, tide0.norm());
    println!("SOLID EARTH TIDE at Epoch 599: [{:.4}, {:.4}, {:.4}], norm={:.4}m", tide599.x, tide599.y, tide599.z, tide599.norm());
    let antex_path = Path::new("datasets/igs/igs14.atx");
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
            println!("Truth Monument: [{:.4}, {:.4}, {:.4}]", st.x, st.y, st.z);
            println!("Truth APC:      [{:.4}, {:.4}, {:.4}]", (st + recv_pco).x, (st + recv_pco).y, (st + recv_pco).z);
        }
    }
}

fn log_sat_prefit_residual(s: &gneiss_rtk::swfg::pipeline::RawObservation, truth: Vector3<f64>) {
    let Some(pr2) = s.pr_l2 else { return; };
    let gamma = (s.f1 / s.f2).powi(2);
    let pr_if = (gamma * s.pr_l1 - pr2) / (gamma - 1.0);
    let modeled = (s.sat_pos_ecef - truth).norm() - s.sat_clock_m + s.tropo_dry_m;
    let res = pr_if - modeled;
    let cp_res_str = if let (Some(cp1), Some(cp2)) = (s.cp_l1, s.cp_l2) {
        let l1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / s.f1;
        let l2 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / s.f2;
        let cp_if = (gamma * cp1 * l1 - cp2 * l2) / (gamma - 1.0);
        format!("{:8.3}m", cp_if - modeled)
    } else {
        "     N/A".to_string()
    };
    let mw_str = if let (Some(cp1), Some(cp2)) = (s.cp_l1, s.cp_l2) {
        let lambda_wl = gneiss_core::constants::SPEED_OF_LIGHT_M_S / (s.f1 - s.f2);
        let p_nl = (s.f1 * s.pr_l1 + s.f2 * pr2) / (s.f1 + s.f2);
        format!("mw={:10.3}cyc", (cp1 - cp2) - p_nl / lambda_wl)
    } else {
        "mw=       N/A".to_string()
    };
    println!("Sat c={} p={:2}: el={:4.1} az={:5.1} | pr_res={:8.3}m  cp_res={}  {}",
        s.constellation_id, s.satellite, s.elevation_rad.to_degrees(), s.azimuth_rad.to_degrees(), res, cp_res_str, mw_str);
}

fn log_epoch_prefit_residuals(
    label: &str,
    ep: &gneiss_core::obs::EpochObs,
    truth: Vector3<f64>,
    products: &PreciseProducts,
    antex_db: Option<&Arc<gneiss_parsers::antex::AntexDatabase>>,
) {
    let Some(ref orbits) = products.0 else { return; };
    let satpos_src = gneiss_rtk::estimators::rtk_iekf::satpos::PreciseSrc {
        orbits,
        clocks: products.1.as_deref(),
    };
    if let Ok(raw_sats) = gneiss_rtk::swfg::engine::epoch::extract_raw_observations_with_source(
        ep, &satpos_src, Some(truth), products.2.as_deref(), antex_db.map(|a| a.as_ref()),
    ) {
        println!("--- PREFIT RESIDUALS AT TRUTH ({}) ---", label);
        for s in &raw_sats {
            log_sat_prefit_residual(s, truth);
        }
    }
}

fn check_prefits_for_dataset(
    spec: &PppDatasetSpec,
    obs_epochs: &[gneiss_core::obs::EpochObs],
    selected_epochs: &[gneiss_core::obs::EpochObs],
    products: &PreciseProducts,
    antex_db: Option<&Arc<gneiss_parsers::antex::AntexDatabase>>,
    truth_map: &BTreeMap<u32, Vector3<f64>>,
    csrs_map: &BTreeMap<u32, Vector3<f64>>,
) {
    for (label, target_ep) in [("Epoch 0", obs_epochs.first()), ("Last Epoch", selected_epochs.last())] {
        if let Some(ep) = target_ep {
            let ep_tow = ep.time.tow.round() as u32;
            let ep_truth = spec.static_truth
                .or_else(|| truth_map.get(&ep_tow).copied())
                .or_else(|| csrs_map.get(&ep_tow).copied());
            if let Some(truth) = ep_truth {
                log_epoch_prefit_residuals(label, ep, truth, products, antex_db);
            }
        }
    }
}



fn evaluate_trajectory_errors(
    traj: &[gneiss_rtk::post_process::SmoothedEpoch],
    static_truth: Option<Vector3<f64>>,
    truth_map: &BTreeMap<u32, Vector3<f64>>,
) -> (Vec<f64>, Vec<f64>) {
    let mut h_errs = Vec::new();
    let mut d3_errs = Vec::new();
    let n_total = traj.len();
    for (i, ep) in traj.iter().enumerate() {
        let tow = ep.time.tow.round() as u32;
        if let Some(truth) = static_truth.or_else(|| truth_map.get(&tow).copied()) {
            let h = compute_horizontal_error(ep.position_ecef, truth);
            let d3 = compute_3d_error(ep.position_ecef, truth);
            if i < 3 || i + 3 >= n_total {
                let diff = ep.position_ecef - truth;
                println!("Epoch {} (TOW={}): sol=[{:.3}, {:.3}, {:.3}], truth=[{:.3}, {:.3}, {:.3}], diff=[{:.3}, {:.3}, {:.3}], h={:.3}m, 3D={:.3}m",
                    i, tow, ep.position_ecef.x, ep.position_ecef.y, ep.position_ecef.z,
                    truth.x, truth.y, truth.z, diff.x, diff.y, diff.z, h, d3);
            }
            h_errs.push(h);
            d3_errs.push(d3);
        }
    }
    (h_errs, d3_errs)
}

fn evaluate_ppp_dataset(spec: &PppDatasetSpec) {
    println!("\n========================================================");
    println!("Evaluating PPP Dataset: {}", spec.name);
    println!("========================================================");
    let Some(((ephemerides, klobuchar), obs_epochs, obs_header)) = load_rinex_inputs(spec) else { return; };
    let products = load_precise_products(spec.sp3_path, spec.clk_path, spec.bia_path);
    let dir = Path::new(spec.dir);
    let truth_map = spec.truth_file.map(|f| parse_pos_truth(&dir.join(f))).unwrap_or_default();
    let csrs_map = spec.csrs_file.map(|f| load_csrs_truth(&dir.join(f))).unwrap_or_default();
    let init_pos = determine_initial_position(spec, &obs_header, &obs_epochs, &ephemerides, &truth_map, &csrs_map);
    println!("Dataset: {}, init_pos: [{:.3}, {:.3}, {:.3}]", spec.name, init_pos.x, init_pos.y, init_pos.z);

    let max_epochs = std::env::var("MAX_EPOCHS").ok().and_then(|v| v.parse::<usize>().ok()).unwrap_or(spec.max_epochs);
    let selected_epochs = &obs_epochs[..obs_epochs.len().min(max_epochs)];
    let antex_path = Path::new("datasets/igs/igs14.atx");
    let antex_db = antex_path.exists().then(|| gneiss_parsers::antex::AntexDatabase::parse(antex_path).ok().map(Arc::new)).flatten();
    let options = build_post_process_options(spec, init_pos, klobuchar.as_ref(), &products, antex_db.clone());
    let config = EngineConfig::Ppp(gneiss_rtk::swfg::config::PppConfig {
        initial_position: Some([init_pos.x, init_pos.y, init_pos.z]),
        window_size: 10,
        is_kinematic: spec.is_kinematic,
        ..Default::default()
    });
    let Ok(result) = execute_post_process(&config, &ephemerides, selected_epochs, None, None, &options) else {
        eprintln!("PPP post-processing failed");
        return;
    };
    log_tides_and_pco(spec, dir, init_pos);
    check_prefits_for_dataset(spec, &obs_epochs, selected_epochs, &products, antex_db.as_ref(), &truth_map, &csrs_map);
    let (h_errs, d3_errs) = evaluate_trajectory_errors(&result.trajectory, spec.static_truth, &truth_map);
    print_stats(&format!("PPP Solution - {}", spec.name), h_errs, d3_errs);
    evaluate_csrs_comparison(&result.trajectory, &csrs_map, &truth_map);
}

fn wtzr_spec() -> PppDatasetSpec {
    PppDatasetSpec {
        name: "WTZR Geodetic Observatory (Germany, 30s, Float PPP)",
        dir: "datasets/wtzr_ppp_1224",
        obs_file: "WTZR00DEU_R_20203590000_01D_30S_MO.rnx",
        nav_file: "BRDC00IGS_R_20203590000_01D_MN.rnx",
        sp3_path: "datasets/wtzr_ppp_1224/com21374.eph",
        clk_path: "datasets/wtzr_ppp_1224/com21374.clk",
        bia_path: Some("datasets/wtzr_ppp_1224/com21374.bia"),
        truth_file: None,
        csrs_file: None,
        static_truth: Some(Vector3::new(4075580.8863, 931853.5784, 4801567.9707)),
        max_epochs: 600,
        is_kinematic: false,
    }
}

fn alic_spec() -> PppDatasetSpec {
    PppDatasetSpec {
        name: "ALIC Geodetic Observatory (Australia, 30s, Float PPP)",
        dir: "datasets/igs",
        obs_file: "alic3350.19o",
        nav_file: "brdc3350.19n",
        sp3_path: "datasets/igs/cod20820.sp3",
        clk_path: "datasets/igs/gfz20820.clk",
        bia_path: None,
        truth_file: None,
        csrs_file: None,
        static_truth: Some(Vector3::new(-4052052.7533, 4212835.9866, -2545104.6062)),
        max_epochs: 600,
        is_kinematic: false,
    }
}

fn f9p_spec() -> PppDatasetSpec {
    PppDatasetSpec {
        name: "RTK Explorer F9P Kinematic Vehicle Drive (1Hz, Integer PPP-AR)",
        dir: "datasets/rtkexplorer/sample_1/f9p_ppp_1224",
        obs_file: "rover.obs",
        nav_file: "BRDC00IGS_R_20203590000_01D_MN.rnx",
        sp3_path: "datasets/rtkexplorer/sample_1/f9p_ppp_1224/ESA0MGNFIN_20203590000_01D_05M_ORB.SP3",
        clk_path: "datasets/rtkexplorer/sample_1/f9p_ppp_1224/ESA0MGNFIN_20203590000_01D_30S_CLK.CLK",
        bia_path: None,
        truth_file: Some("rover_ppk.pos"),
        csrs_file: Some("rover_csrs.pos"),
        static_truth: None,
        max_epochs: 60,
        is_kinematic: true,
    }
}

fn main() {
    tracing_subscriber::fmt().with_env_filter("warn").with_target(false).without_time().try_init().ok();
    println!("========================================================");
    println!("Gneiss Standalone & Kinematic PPP Evaluation Harness");
    println!("========================================================");

    let only = std::env::var("PPP_ONLY").ok();
    let specs = [("wtzr", wtzr_spec()), ("alic", alic_spec()), ("f9p", f9p_spec())];
    for (key, spec) in specs {
        if only.as_deref().is_none_or(|w| w == key) && Path::new(spec.dir).exists() {
            evaluate_ppp_dataset(&spec);
        }
    }
}
