//! Multi-base network RTK benchmark on a NOAA CORS network (2020-05-14).
//!
//! Rover: P224 (Sibley Volcanic, Berkeley Hills) with NGS monument truth.
//! Bases: P181/OHLN/CAPO/P225/P222/SLAC at 15.0-49.7 km baselines.
//! Base positions: NOAA coord_14 ITRF2014 (ARP) propagated to the data epoch.
//! Run each base independently through forward + smoothed PPK to establish
//! the per-base error surface before multi-base fusion.
//!
//! Frame contract (all DOF):
//! - Base positions and rover truth: ITRF2014 @ 2020.367 (coord_14 + velocity).
//!   The RINEX header approx positions are NAD83(2011) and are ~1.5-1.8 m off;
//!   they must NOT be used as truth or base positions.
//! - Solution frame: broadcast ephemeris (WGS84), ~cm-consistent with ITRF2014.
//! - Time: GPS time in all obs/nav/truth files.
//! - Static monuments: zero velocity; no attitude DOF in this benchmark.
//! - Rover filter init: SPP in broadcast frame (None here); never seed from
//!   the NAD83 RINEX header.
//!
//! Usage: cargo run --release --bin eval_network_ppk
//!        MAX_EPOCHS=1440 cargo run --release --bin eval_network_ppk

use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use nalgebra::Vector3;

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_core::time::GpsTime;
use gneiss_rtk::post_process::{execute_post_process, network, PostProcessOptions, SmoothedEpoch};
use gneiss_rtk::swfg::config::EngineConfig;

const ROVER_FILE: &str = "p2241350.20o";
const TRUTH_FILE: &str = "p224_truth.pos";
const NAV_FILE: &str = "brdc1350.20n";

struct NetworkBase {
    id: &'static str,
    base_file: &'static str,
    base_pos: Vector3<f64>,
    baseline_km: f64,
}

// ITRF2014 ARP positions propagated from epoch 2010.0 to 2020.367 with
// the per-station velocities published in the NOAA coord_14 files.
const BASES: &[NetworkBase] = &[
    NetworkBase { id: "P181", base_file: "p1811350.20o", base_pos: Vector3::new(-2697941.2851, -4255089.1805, 3898009.7146), baseline_km: 14.97 },
    NetworkBase { id: "OHLN", base_file: "ohln1350.20o", base_pos: Vector3::new(-2686856.6743, -4254625.7745, 3905990.5995), baseline_km: 16.50 },
    NetworkBase { id: "CAPO", base_file: "capo1350.20o", base_pos: Vector3::new(-2693675.7831, -4273829.9413, 3880383.2888), baseline_km: 16.63 },
    NetworkBase { id: "P225", base_file: "p2251350.20o", base_pos: Vector3::new(-2681519.0648, -4281621.5811, 3880440.4336), baseline_km: 21.86 },
    NetworkBase { id: "P222", base_file: "p2221350.20o", base_pos: Vector3::new(-2689640.5315, -4290437.1233, 3865051.0346), baseline_km: 37.97 },
    NetworkBase { id: "SLAC", base_file: "slac1350.20o", base_pos: Vector3::new(-2703116.3008, -4291766.8064, 3854248.0730), baseline_km: 49.67 },
];

type Truth = BTreeMap<u32, Vector3<f64>>;

fn parse_truth(path: &Path) -> Truth {
    let mut truth = Truth::new();
    let content = match std::fs::read_to_string(path) {
        Ok(c) => c,
        Err(_) => return truth,
    };
    for line in content.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 5 || parts[0].starts_with('%') {
            continue;
        }
        let date: Vec<&str> = parts[0].split('/').collect();
        let time: Vec<&str> = parts[1].split(':').collect();
        let parsed = (
            date[0].parse::<i32>(), date[1].parse::<i32>(), date[2].parse::<i32>(),
            time[0].parse::<i32>(), time[1].parse::<i32>(), time[2].parse::<f64>(),
            parts[2].parse::<f64>(), parts[3].parse::<f64>(), parts[4].parse::<f64>(),
        );
        let (y, m, d, hr, min, sec, px, py, pz) = match (parsed.0, parsed.1, parsed.2, parsed.3, parsed.4, parsed.5, parsed.6, parsed.7, parsed.8) {
            (Ok(y), Ok(m), Ok(d), Ok(hr), Ok(min), Ok(sec), Ok(px), Ok(py), Ok(pz)) => (y, m, d, hr, min, sec, px, py, pz),
            _ => continue,
        };
        let gps_time = GpsTime::from_calendar(y, m, d, hr, min, sec);
        truth.insert(gps_time.tow.round() as u32, Vector3::new(px, py, pz));
    }
    truth
}

fn horizontal_error(pos: Vector3<f64>, truth: Vector3<f64>) -> f64 {
    let llh = gneiss_core::coords::ecef_to_llh(truth);
    let ned = gneiss_core::coords::ecef_to_ned_matrix(llh) * (pos - truth);
    (ned.x * ned.x + ned.y * ned.y).sqrt()
}

/// Signed up-axis error (m): positive = solution above truth. The vertical
/// axis is where residual troposphere shows up first in RTK.
fn vertical_error(pos: Vector3<f64>, truth: Vector3<f64>) -> f64 {
    let llh = gneiss_core::coords::ecef_to_llh(truth);
    let ned = gneiss_core::coords::ecef_to_ned_matrix(llh) * (pos - truth);
    -ned.z
}

fn print_stats(
    name: &str,
    mut h_errs: Vec<f64>,
    mut d3_errs: Vec<f64>,
    mut up_errs: Vec<f64>,
    fix_count: usize,
    total_count: usize,
) -> [f64; 4] {
    if h_errs.is_empty() {
        return [0.0; 4];
    }
    h_errs.sort_by(|a, b| a.total_cmp(b));
    d3_errs.sort_by(|a, b| a.total_cmp(b));
    up_errs.sort_by(|a, b| a.total_cmp(b));
    let n = h_errs.len();
    let p50 = h_errs[n / 2];
    let p68 = h_errs[(n as f64 * 0.68) as usize];
    let p95 = h_errs[(n as f64 * 0.95) as usize];
    let rms = (h_errs.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();
    let up_p50 = up_errs[n / 2];
    let up_rms = (up_errs.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();
    let fix_pct = (fix_count as f64 / total_count.max(1) as f64) * 100.0;
    println!("=== {} (N={}, Fixed={}/{} [{:.1}%]) ===", name, n, fix_count, total_count, fix_pct);
    println!("Horizontal Error:  p50={:.3}m,  p68={:.3}m,  p95={:.3}m,  RMS={:.3}m", p50, p68, p95, rms);
    println!("Vertical Error:    p50={:+.3}m,  RMS={:.3}m", up_p50, up_rms);
    println!("3D Position Error: p50={:.3}m,  p95={:.3}m", d3_errs[n / 2], d3_errs[(n as f64 * 0.95) as usize]);
    [p50, p68, p95, rms]
}

fn print_fixed_stats(name: &str, mut h_errs: Vec<f64>) -> f64 {
    if h_errs.is_empty() {
        return 0.0;
    }
    h_errs.sort_by(|a, b| a.total_cmp(b));
    let p50 = h_errs[h_errs.len() / 2];
    println!("  [{} fixed-only horizontal p50: {:.3}m, N={}]", name, p50, h_errs.len());
    p50
}

fn collect_errors(traj: &[SmoothedEpoch], truth: &Truth) -> (Vec<f64>, Vec<f64>, Vec<f64>, usize, Vec<f64>) {
    let mut h_errs = Vec::new();
    let mut d3_errs = Vec::new();
    let mut up_errs = Vec::new();
    let mut fixed_errs = Vec::new();
    let mut fix = 0;
    for ep in traj {
        let tow = ep.time.tow.round() as u32;
        if ep.quality == 1 {
            fix += 1;
        }
        if let Some(&t) = truth.get(&tow) {
            let h = horizontal_error(ep.position_ecef, t);
            if h < 100.0 {
                h_errs.push(h);
                d3_errs.push((ep.position_ecef - t).norm());
                up_errs.push(vertical_error(ep.position_ecef, t));
                if ep.quality == 1 {
                    fixed_errs.push(h);
                }
            }
        }
    }
    (h_errs, d3_errs, up_errs, fix, fixed_errs)
}

struct RunContext<'a> {
    ephemerides: &'a [Ephemeris],
    truth: &'a Truth,
    rover_init: Option<Vector3<f64>>,
    klob: Option<([f64; 4], [f64; 4])>,
}

fn run_pass(
    config: &EngineConfig,
    ctx: &RunContext,
    rover: &[EpochObs],
    base_epochs: &[EpochObs],
    base: &NetworkBase,
    label: &str,
    bidir: bool,
) -> ([f64; 4], Vec<SmoothedEpoch>) {
    let options = PostProcessOptions {
        enable_bidirectional: bidir,
        base_position: Some(base.base_pos),
        initial_rover_position: ctx.rover_init,
        klobuchar_alpha: ctx.klob.map(|k| k.0),
        klobuchar_beta: ctx.klob.map(|k| k.1),
        // Static monuments: q=1.0 re-randomizes position ~55 m per 30 s
        // epoch and keeps the float solution from converging (ambiguity
        // floats sit 0.4+ cycles off, blocking AR). 1e-6 allows slow drift.
        q_accel: Some(1e-6),
        // Long baselines need iono-immune fixing: MW wide-lane cascade AR
        // unlocks the iono-free stage beyond ~20 km.
        widelane_ar: std::env::var("WL_DISABLE").is_err(),
    };
    let res = match execute_post_process(config, ctx.ephemerides, rover, Some(base_epochs), None, &options) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("{} pass failed for {}: {}", label, base.id, e);
            return ([0.0; 4], Vec::new());
        }
    };
    // Honesty gating on the fused bidirectional product: an excursion that
    // jumps beyond what a static monument can do is reported as float, both
    // in these stats and downstream in network consensus.
    let mut traj = res.trajectory;
    if bidir && options.widelane_ar {
        traj = gneiss_rtk::post_process::network::apply_continuity_gate(
            traj,
            gneiss_rtk::post_process::network::CONTINUITY_JUMP_M,
            gneiss_rtk::post_process::network::CONTINUITY_MAX_DT_S,
        );
    }
    let (h, d3, up_errs, fix, fixed_errs) = collect_errors(&traj, ctx.truth);
    if std::env::var("WL_OUTLIERS").is_ok() {
        for ep in &traj {
            if let Some(&t) = ctx.truth.get(&(ep.time.tow.round() as u32)) {
                let herr = horizontal_error(ep.position_ecef, t);
                let verr = vertical_error(ep.position_ecef, t);
                if herr > 1.0 || verr.abs() > 0.30 {
                    eprintln!(
                        "OUTLIER {} {} tow={:.0} h={:.2} v={:+.3} q={} sep={:.1} nsat={}",
                        label, base.id, ep.time.tow, herr, verr, ep.quality,
                        ep.separation_3d, ep.n_satellites,
                    );
                }
            }
        }
    }
    let stats = print_stats(&format!("{} RTK [{}] ({:.1} km)", label, base.id, base.baseline_km), h, d3, up_errs, fix, traj.len());
    print_fixed_stats(label, fixed_errs);
    (stats, traj)
}

fn run_base(base: &NetworkBase, dir: &Path, ctx: &RunContext, rover: &[EpochObs]) -> ([f64; 8], Vec<SmoothedEpoch>) {
    let base_f = match File::open(dir.join(base.base_file)) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Failed to open {}: {}", base.base_file, e);
            return ([0.0; 8], Vec::new());
        }
    };
    let Ok((base_epochs, _)) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(base_f)) else {
        eprintln!("Failed to parse {}", base.base_file);
        return ([0.0; 8], Vec::new());
    };
    let config = EngineConfig::Rtk(gneiss_rtk::swfg::config::RtkConfig {
        initial_position: Some([base.base_pos.x, base.base_pos.y, base.base_pos.z]),
        ..Default::default()
    });
    let (fwd, _fwd_traj) = run_pass(&config, ctx, rover, &base_epochs, base, "Forward", false);
    let (smooth, smooth_traj) = run_pass(&config, ctx, rover, &base_epochs, base, "Smoothed", true);
    let mut out = [0.0; 8];
    out[..4].copy_from_slice(&fwd);
    out[4..].copy_from_slice(&smooth);
    (out, smooth_traj)
}

fn print_summary(results: &[(&NetworkBase, [f64; 8])], n_epochs: usize) {
    println!("========================================================");
    println!("Multi-Base Network RTK Benchmark (P224 rover, {:.1} h, {} epochs)", n_epochs as f64 * 30.0 / 3600.0, n_epochs);
    println!("========================================================");
    println!("\n{:>6} {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}", "base", "km", "fwd.p50", "fwd.p95", "fwd.rms", "sm.p50", "sm.p68", "sm.p95", "sm.rms");
    for (base, stats) in results {
        println!("{:>6} {:>7.1} {:>9.3} {:>9.3} {:>9.3} {:>9.3} {:>9.3} {:>9.3} {:>9.3}",
            base.id, base.baseline_km, stats[0], stats[2], stats[3], stats[4], stats[5], stats[6], stats[7]);
    }
}

fn select_rover_epochs(epochs: &[EpochObs]) -> &[EpochObs] {
    let max = std::env::var("MAX_EPOCHS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(2880);
    &epochs[..epochs.len().min(max)]
}

struct Dataset {
    ephemerides: Vec<Ephemeris>,
    klob: Option<([f64; 4], [f64; 4])>,
    truth: Truth,
    rover_epochs: Vec<EpochObs>,
}

fn load_dataset(dir: &Path) -> Result<Dataset, String> {
    if !dir.join(ROVER_FILE).exists() {
        return Err("Dataset missing: run scripts/fetch_high_fidelity_suite.py first".to_string());
    }
    let nav_f = File::open(dir.join(NAV_FILE)).map_err(|e| format!("Failed to open nav: {}", e))?;
    let (ephemerides, klobuchar) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_f))
        .map_err(|e| format!("Failed to parse nav: {}", e))?;
    let rov_f = File::open(dir.join(ROVER_FILE)).map_err(|e| format!("Failed to open rover obs: {}", e))?;
    let (rover_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(rov_f))
        .map_err(|e| format!("Failed to parse rover obs: {}", e))?;
    let truth = parse_truth(&dir.join(TRUTH_FILE));
    Ok(Dataset {
        ephemerides,
        klob: klobuchar.map(|k| (k.alpha, k.beta)),
        truth,
        rover_epochs,
    })
}

fn main() {
    tracing_subscriber::fmt()
        .with_env_filter(std::env::var("RUST_LOG").unwrap_or_else(|_| "warn".to_string()))
        .with_target(false)
        .without_time()
        .try_init()
        .ok();

    let dir = Path::new("datasets/cors_short_baseline");
    let data = match load_dataset(dir) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1);
        }
    };

    // Rover filter init is left to the SPP in broadcast frame (None): the
    // RINEX header approx is NAD83 and must not seed the filter.
    let ctx = RunContext { ephemerides: &data.ephemerides, truth: &data.truth, rover_init: None, klob: data.klob };
    let selected_rover = select_rover_epochs(&data.rover_epochs);
    let mut results = Vec::new();
    let mut base_trajs: Vec<Vec<SmoothedEpoch>> = Vec::new();
    let only = std::env::var("WL_ONLY_BASE").ok();
    for base in BASES {
        if let Some(want) = &only { if base.id != want.as_str() { continue; } }
        let (stats, traj) = run_base(base, dir, &ctx, selected_rover);
        if stats.iter().any(|s| *s > 0.0) {
            results.push((base, stats));
        }
        if !traj.is_empty() {
            base_trajs.push(traj);
        }
    }
    print_summary(&results, selected_rover.len());

    // Multi-base network consensus: combine all per-base smoothed
    // solutions into a single network product and report it.
    if base_trajs.len() >= 2 {
        let fused = network::fuse_network_solutions(
            &base_trajs, &network::NetworkConsensusConfig::default(),
        );
        let gated = network::apply_continuity_gate(fused, 0.20, 90.0);
        let (h, d3, up_errs, fix, fixed_errs) = collect_errors(&gated, ctx.truth);
        print_stats(&format!("NETWORK FUSED [{} bases]", base_trajs.len()), h, d3, up_errs, fix, gated.len());
        print_fixed_stats("Network", fixed_errs);
    }
}
