//! Benchmark evaluation of Gneiss Multi-GNSS Standalone and Kinematic PPP engine.
//!
//! Evaluates static observatories (WTZR, ALIC) and moving kinematic rovers (F9P)
//! against CSRS-PPP commercial truth, RTK ground truth, and calibrated LocalDatumTie.
//! Uses precise products (SP3, CLK, SINEX BIA, Bernese DCB) with multi-pass initialization.
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

#[derive(Default)]
struct PppDatasetSpec {
    name: &'static str,
    dir: &'static str,
    obs_file: &'static str,
    nav_file: &'static str,
    sp3_path: &'static str,
    clk_path: &'static str,
    bia_path: Option<&'static str>,
    dcb_path: Option<&'static str>,
    truth_file: Option<&'static str>,
    csrs_file: Option<&'static str>,
    static_truth: Option<Vector3<f64>>,
    max_epochs: usize,
    is_kinematic: bool,
    enable_glonass: bool,
    enable_galileo: bool,
    truth_is_nad83: bool,
}

fn compute_horizontal_error(pos: Vector3<f64>, truth: Vector3<f64>) -> f64 {
    let diff = gneiss_core::coords::ecef_to_ned_matrix(gneiss_core::coords::ecef_to_llh(truth)) * (pos - truth);
    (diff.x * diff.x + diff.y * diff.y).sqrt()
}

fn compute_3d_error(pos: Vector3<f64>, truth: Vector3<f64>) -> f64 {
    (pos - truth).norm()
}

fn parse_pos_line(line: &str) -> Option<(u32, Vector3<f64>)> {
    let p: Vec<&str> = line.split_whitespace().collect();
    if p.len() < 5 { return None; }
    let (dp, tp): (Vec<&str>, Vec<&str>) = (p[0].split('/').collect(), p[1].split(':').collect());
    if dp.len() != 3 || tp.len() != 3 { return None; }
    let (y, m, d) = (dp[0].parse().ok()?, dp[1].parse().ok()?, dp[2].parse().ok()?);
    let (hr, min, sec) = (tp[0].parse().ok()?, tp[1].parse().ok()?, tp[2].parse().ok()?);
    let pos = Vector3::new(p[2].parse().ok()?, p[3].parse().ok()?, p[4].parse().ok()?);
    Some((GpsTime::from_calendar(y, m, d, hr, min, sec).tow.round() as u32, pos))
}

fn parse_pos_truth(path: &Path, is_nad83: bool) -> BTreeMap<u32, Vector3<f64>> {
    std::fs::read_to_string(path).ok().map_or_else(BTreeMap::new, |c| {
        c.lines().filter(|l| !l.trim().is_empty() && !l.trim().starts_with('%')).filter_map(parse_pos_line).map(|(tow, pos)| {
            let p = if is_nad83 && std::env::var("NORMALIZE_DATUM").is_ok_and(|v| v == "1") {
                gneiss_core::frames::EcefPos::<gneiss_core::frames::Nad83_2011>::new(pos).convert_to::<gneiss_core::frames::Itrf2014>(2020.98).into_vector()
            } else { pos };
            (tow, p)
        }).collect()
    })
}

fn load_csrs_truth(path: &Path) -> BTreeMap<u32, Vector3<f64>> {
    let mut map = BTreeMap::new();
    if let Ok(f) = File::open(path) {
        if let Ok(epochs) = gneiss_parsers::csrs_pos::parse_csrs_pos(BufReader::new(f)) {
            for ep in epochs { map.insert(ep.time.tow.round() as u32, ep.pos_ecef); }
        }
    }
    map
}

#[derive(Default)]
struct CsrsComparisonStats {
    g_vs_csrs_h: Vec<f64>, g_vs_csrs_3d: Vec<f64>,
    csrs_vs_t_h: Vec<f64>, csrs_vs_t_3d: Vec<f64>,
    g_vs_t_h: Vec<f64>, g_vs_t_3d: Vec<f64>,
}

fn record_csrs_epoch(
    ep: &gneiss_rtk::post_process::SmoothedEpoch,
    c_pos: Vector3<f64>,
    truth_map: &BTreeMap<u32, Vector3<f64>>,
    stats: &mut CsrsComparisonStats,
    debug: bool,
) {
    let tow = ep.time.tow.round() as u32;
    let h = compute_horizontal_error(ep.position_ecef, c_pos);
    let d3 = compute_3d_error(ep.position_ecef, c_pos);
    if debug {
        let llh = gneiss_core::coords::ecef_to_llh(c_pos);
        let ned = gneiss_core::coords::ecef_to_ned_matrix(llh) * (ep.position_ecef - c_pos);
        println!("Epoch {} (TOW={}): E={:+.4}m, N={:+.4}m, U={:+.4}m, H={:.4}m, 3D={:.4}m",
            ep.time.tow, tow, ned.y, ned.x, -ned.z, h, d3);
    }
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
    let n = traj.len();
    for (i, ep) in traj.iter().enumerate() {
        let tow = ep.time.tow.round() as u32;
        if let Some(&c_pos) = csrs_map.get(&tow) {
            record_csrs_epoch(ep, c_pos, truth_map, &mut stats, i < 3 || i + 3 >= n);
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
    let pct = |p: f64| h_errs[(n as f64 * p) as usize];
    let (p10, p25, p50, p68, p75) = (pct(0.10), pct(0.25), h_errs[n / 2], pct(0.68), pct(0.75));
    let (p90, p95, p99, max) = (pct(0.90), pct(0.95), pct(0.99), *h_errs.last().unwrap_or(&0.0));
    let (mean, rms) = (h_errs.iter().sum::<f64>() / n as f64, (h_errs.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt());

    println!("=== {} (N={}) ===", name, n);
    println!("Horizontal Error:  p50={:.3}m,  p68={:.3}m,  p95={:.3}m,  RMS={:.3}m", p50, p68, p95, rms);
    println!("  Horizontal CDF:  p10={:.3}m, p25={:.3}m, p50={:.3}m, p68={:.3}m, p75={:.3}m, p90={:.3}m, p95={:.3}m, p99={:.3}m, max={:.3}m, mean={:.3}m",
        p10, p25, p50, p68, p75, p90, p95, p99, max, mean);

    if !d3_errs.is_empty() {
        d3_errs.sort_by(|a, b| a.total_cmp(b));
        let m = d3_errs.len();
        let (p50_3d, p95_3d) = (d3_errs[m / 2], d3_errs[(m as f64 * 0.95) as usize]);
        let rms_3d = (d3_errs.iter().map(|e| e * e).sum::<f64>() / m as f64).sqrt();
        println!("3D Position Error: p50={:.3}m,  p95={:.3}m,  RMS={:.3}m", p50_3d, p95_3d, rms_3d);
    }
}

type PreciseProducts = (Option<Arc<PreciseOrbit>>, Option<Arc<RinexClock>>, Option<Arc<gneiss_parsers::sinex_bia::SinexBias>>);

fn load_sp3_orbit(path: &str) -> Option<Arc<PreciseOrbit>> {
    Path::new(path).exists().then(|| File::open(path).ok().and_then(|f| gneiss_parsers::sp3::parse_sp3(BufReader::new(f)).ok()).map(|e| Arc::new(PreciseOrbit::new(e)))).flatten()
}

fn load_rinex_clock(path: &str) -> Option<Arc<RinexClock>> {
    Path::new(path).exists().then(|| std::fs::read_to_string(path).ok().map(|s| Arc::new(RinexClock::parse(&s)))).flatten()
}

fn load_sinex_bias(bia_path: Option<&str>, dcb_path: Option<&str>) -> Option<Arc<gneiss_parsers::sinex_bia::SinexBias>> {
    let mut bia = bia_path
        .filter(|p| Path::new(p).exists())
        .and_then(|p| File::open(p).ok())
        .and_then(|f| gneiss_parsers::sinex_bia::SinexBias::parse(BufReader::new(f)).ok());

    if let Some(f) = dcb_path.filter(|p| Path::new(p).exists()).and_then(|p| File::open(p).ok()) {
        let reader = BufReader::new(f);
        if let Some(b) = bia.as_mut() {
            let _ = b.load_bernese_dcb(reader);
        } else if let Ok(recs) = gneiss_parsers::bernese_dcb::parse_bernese_dcb(reader) {
            bia = Some(gneiss_parsers::sinex_bia::SinexBias::new(recs));
        }
    }
    bia.map(Arc::new)
}

fn load_precise_products(sp3: &str, clk: &str, bia: Option<&str>, dcb: Option<&str>) -> PreciseProducts {
    (load_sp3_orbit(sp3), load_rinex_clock(clk), load_sinex_bias(bia, dcb))
}

type RinexNavData = (Vec<gneiss_core::ephemeris::Ephemeris>, Option<gneiss_core::atmosphere::KlobucharParams>);

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
    let enable_bidirectional = std::env::var("BIDIRECTIONAL").map_or(spec.is_kinematic, |v| v == "1");
    PostProcessOptions {
        enable_bidirectional,
        initial_rover_position: Some(init_pos),
        klobuchar_alpha: klobuchar.map(|k| k.alpha),
        klobuchar_beta: klobuchar.map(|k| k.beta),
        widelane_ar: true,
        dynamics,
        enable_glonass: spec.enable_glonass,
        precise_orbits: products.0.clone(),
        precise_clocks: products.1.clone(),
        sinex_bias: products.2.clone(),
        antex_database: antex_db,
        init_passes: std::env::var("INIT_PASSES").ok().and_then(|v| v.parse().ok()).unwrap_or(1),
        ..Default::default()
    }
}

fn log_tides_and_pco(spec: &PppDatasetSpec, dir: &Path, init_pos: Vector3<f64>) {
    let (t0, t1) = (gneiss_geodesy::tides::solid_earth_tide(init_pos, 345600.0, 2137), gneiss_geodesy::tides::solid_earth_tide(init_pos, 363570.0, 2137));
    println!("SOLID EARTH TIDE at Epoch 0:   [{:.4}, {:.4}, {:.4}], norm={:.4}m", t0.x, t0.y, t0.z, t0.norm());
    println!("SOLID EARTH TIDE at Epoch 599: [{:.4}, {:.4}, {:.4}], norm={:.4}m", t1.x, t1.y, t1.z, t1.norm());
    let antex_path = Path::new("datasets/igs/igs14.atx");
    let recv_pco = antex_path.to_str().and_then(|ap| {
        antex_path.exists().then(|| gneiss_rtk::post_process::antenna::station_recv_pco_ecef(&dir.join(spec.obs_file), ap, init_pos)).flatten()
    }).unwrap_or_else(Vector3::zeros);
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
    let (cp_res, mw) = if let (Some(cp1), Some(cp2)) = (s.cp_l1, s.cp_l2) {
        let (l1, l2) = (gneiss_core::constants::SPEED_OF_LIGHT_M_S / s.f1, gneiss_core::constants::SPEED_OF_LIGHT_M_S / s.f2);
        let cp_if = (gamma * cp1 * l1 - cp2 * l2) / (gamma - 1.0);
        let l_wl = gneiss_core::constants::SPEED_OF_LIGHT_M_S / (s.f1 - s.f2);
        (format!("{:8.3}m", cp_if - modeled), format!("mw={:10.3}cyc", (cp1 - cp2) - (s.f1 * s.pr_l1 + s.f2 * pr2) / ((s.f1 + s.f2) * l_wl)))
    } else { ("     N/A".to_string(), "mw=       N/A".to_string()) };
    println!("Sat c={} p={:2}: el={:4.1} az={:5.1} | pr_res={:8.3}m  cp_res={}  {}",
        s.constellation_id, s.satellite, s.elevation_rad.to_degrees(), s.azimuth_rad.to_degrees(), pr_if - modeled, cp_res, mw);
}

fn log_epoch_prefit_residuals(
    label: &str,
    ep: &gneiss_core::obs::EpochObs,
    truth: Vector3<f64>,
    spec: &PppDatasetSpec,
    products: &PreciseProducts,
    antex_db: Option<&Arc<gneiss_parsers::antex::AntexDatabase>>,
    ephemerides: Option<&[gneiss_core::ephemeris::Ephemeris]>,
) {
    let Some(ref orbits) = products.0 else { return; };
    let satpos_src = gneiss_rtk::estimators::rtk_iekf::satpos::PreciseSrc {
        orbits,
        clocks: products.1.as_deref(),
    };
    if let Ok(raw_sats) = gneiss_rtk::swfg::engine::epoch::extract_raw_observations_with_source(
        ep, &satpos_src, Some(truth), products.2.as_deref(), antex_db.map(|a| a.as_ref()),
        ephemerides, spec.enable_glonass, spec.enable_galileo,
    ) {
        println!("--- PREFIT RESIDUALS AT TRUTH ({}) ---", label);
        for s in &raw_sats {
            log_sat_prefit_residual(s, truth);
        }
    }
}

fn check_prefits_for_dataset(
    spec: &PppDatasetSpec,
    epochs: (&[gneiss_core::obs::EpochObs], &[gneiss_core::obs::EpochObs]),
    products: &PreciseProducts,
    antex_db: Option<&Arc<gneiss_parsers::antex::AntexDatabase>>,
    truth_map: &BTreeMap<u32, Vector3<f64>>,
    csrs_map: &BTreeMap<u32, Vector3<f64>>,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
) {
    for (label, target_ep) in [("Epoch 0", epochs.0.first()), ("Last Epoch", epochs.1.last())] {
        if let Some(ep) = target_ep {
            let ep_tow = ep.time.tow.round() as u32;
            let ep_truth = spec.static_truth
                .or_else(|| truth_map.get(&ep_tow).copied())
                .or_else(|| csrs_map.get(&ep_tow).copied());
            if let Some(truth) = ep_truth {
                log_epoch_prefit_residuals(label, ep, truth, spec, products, antex_db, Some(ephemerides));
            }
        }
    }
}
fn evaluate_trajectory_errors(
    traj: &[gneiss_rtk::post_process::SmoothedEpoch],
    static_truth: Option<Vector3<f64>>,
    truth_map: &BTreeMap<u32, Vector3<f64>>,
    tie: Option<gneiss_geodesy::LocalDatumTie>,
) -> (Vec<f64>, Vec<f64>) {
    let mut h_errs = Vec::new();
    let mut d3_errs = Vec::new();
    let n_total = traj.len();
    for (i, ep) in traj.iter().enumerate() {
        let tow = ep.time.tow.round() as u32;
        if let Some(truth) = static_truth.or_else(|| truth_map.get(&tow).copied()) {
            let pos = tie.map_or(ep.position_ecef, |t| t.transform(ep.position_ecef));
            let h = compute_horizontal_error(pos, truth);
            let d3 = compute_3d_error(pos, truth);
            if (i < 3 || i + 3 >= n_total) && tie.is_none() {
                let llh = gneiss_core::coords::ecef_to_llh(truth);
                let r_enu = gneiss_core::coords::ecef_to_ned_matrix(llh);
                let enu = r_enu * (pos - truth);
                println!("Epoch {} (TOW={}): E={:+.3}m, N={:+.3}m, U={:+.3}m, h={:.3}m, 3D={:.3}m",
                    i, tow, enu.y, enu.x, -enu.z, h, d3);
            }
            h_errs.push(h);
            d3_errs.push(d3);
        }
    }
    (h_errs, d3_errs)
}

fn run_ppp_engine(
    spec: &PppDatasetSpec,
    init_pos: Vector3<f64>,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    epochs: &[gneiss_core::obs::EpochObs],
    options: &PostProcessOptions,
) -> Option<gneiss_rtk::post_process::PostProcessResult> {
    let config = EngineConfig::Ppp(gneiss_rtk::swfg::config::PppConfig {
        initial_position: Some([init_pos.x, init_pos.y, init_pos.z]),
        window_size: 10,
        is_kinematic: spec.is_kinematic,
        enable_glonass: spec.enable_glonass,
        enable_galileo: spec.enable_galileo,
        ..Default::default()
    });
    execute_post_process(&config, ephemerides, epochs, None, None, options).ok()
}

fn evaluate_calibrated_trajectory(
    spec: &PppDatasetSpec,
    traj: &[gneiss_rtk::post_process::SmoothedEpoch],
    truth_map: &BTreeMap<u32, Vector3<f64>>,
) {
    if !spec.is_kinematic || truth_map.is_empty() { return; }
    let cal_epochs = std::env::var("INIT_CAL_EPOCHS").ok().and_then(|v| v.parse().ok()).unwrap_or(30);
    let pairs: Vec<(Vector3<f64>, Vector3<f64>)> = traj.iter().take(cal_epochs)
        .filter_map(|ep| truth_map.get(&(ep.time.tow.round() as u32)).map(|&t| (ep.position_ecef, t))).collect();
    if let Some(tie) = gneiss_geodesy::LocalDatumTie::estimate(&pairs) {
        println!("Initial-Pass Local Datum Tie (ECEF): [{:+.4}, {:+.4}, {:+.4}], norm={:.4}m",
            tie.translation.x, tie.translation.y, tie.translation.z, tie.translation.norm());
        let (c_h, c_d3) = evaluate_trajectory_errors(traj, spec.static_truth, truth_map, Some(tie));
        print_stats(&format!("PPP Solution (Initial-Pass Calibrated) - {}", spec.name), c_h, c_d3);
    }
}

fn evaluate_ppp_dataset(spec: &PppDatasetSpec) {
    println!("\n========================================================");
    println!("Evaluating PPP Dataset: {}", spec.name);
    println!("========================================================");
    let Some(((ephemerides, klobuchar), obs_epochs, obs_header)) = load_rinex_inputs(spec) else { return; };
    let products = load_precise_products(spec.sp3_path, spec.clk_path, spec.bia_path, spec.dcb_path);
    let dir = Path::new(spec.dir);
    let truth_map = spec.truth_file.map(|f| parse_pos_truth(&dir.join(f), spec.truth_is_nad83)).unwrap_or_default();
    let csrs_map = spec.csrs_file.map(|f| load_csrs_truth(&dir.join(f))).unwrap_or_default();
    let init_pos = determine_initial_position(spec, &obs_header, &obs_epochs, &ephemerides, &truth_map, &csrs_map);
    println!("Dataset: {}, init_pos: [{:.3}, {:.3}, {:.3}]", spec.name, init_pos.x, init_pos.y, init_pos.z);

    let max_epochs = std::env::var("MAX_EPOCHS").ok().and_then(|v| v.parse::<usize>().ok()).unwrap_or(spec.max_epochs);
    let selected_epochs = &obs_epochs[..obs_epochs.len().min(max_epochs)];
    let antex_path = Path::new("datasets/igs/igs14.atx");
    let antex_db = antex_path.exists().then(|| gneiss_parsers::antex::AntexDatabase::parse(antex_path).ok().map(Arc::new)).flatten();
    let options = build_post_process_options(spec, init_pos, klobuchar.as_ref(), &products, antex_db.clone());
    let Some(result) = run_ppp_engine(spec, init_pos, &ephemerides, selected_epochs, &options) else {
        eprintln!("PPP post-processing failed");
        return;
    };
    log_tides_and_pco(spec, dir, init_pos);
    check_prefits_for_dataset(spec, (&obs_epochs, selected_epochs), &products, antex_db.as_ref(), &truth_map, &csrs_map, &ephemerides);
    let (h_errs, d3_errs) = evaluate_trajectory_errors(&result.trajectory, spec.static_truth, &truth_map, None);
    print_stats(&format!("PPP Solution - {}", spec.name), h_errs, d3_errs);
    evaluate_calibrated_trajectory(spec, &result.trajectory, &truth_map);
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
        static_truth: Some(Vector3::new(4075580.8863, 931853.5784, 4801567.9707)),
        max_epochs: 600,
        ..Default::default()
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
        static_truth: Some(Vector3::new(-4052052.7533, 4212835.9866, -2545104.6062)),
        max_epochs: 600,
        ..Default::default()
    }
}

fn f9p_spec() -> PppDatasetSpec {
    PppDatasetSpec {
        name: "RTK Explorer F9P Kinematic Vehicle Drive (1Hz, Integer PPP-AR)",
        dir: "datasets/rtkexplorer/sample_1/f9p_ppp_1224",
        obs_file: "rover.obs",
        nav_file: "rover.nav",
        sp3_path: "datasets/rtkexplorer/sample_1/f9p_ppp_1224/ESA0MGNFIN_20203590000_01D_05M_ORB.SP3",
        clk_path: "datasets/rtkexplorer/sample_1/f9p_ppp_1224/ESA0MGNFIN_20203590000_01D_30S_CLK.CLK",
        bia_path: Some("datasets/rtkexplorer/sample_1/f9p_ppp_1224/COD0MGXFIN_20203590000_01D_01D_OSB.BIA"),
        dcb_path: Some("datasets/rtkexplorer/sample_1/f9p_ppp_1224/P2C22011_RINEX.DCB"),
        truth_file: Some("rover_ppk.pos"),
        csrs_file: Some("rover_csrs.pos"),
        max_epochs: 600,
        is_kinematic: true,
        enable_glonass: true,
        enable_galileo: true,
        truth_is_nad83: true,
        ..Default::default()
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
