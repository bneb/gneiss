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

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use nalgebra::Vector3;

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_core::time::GpsTime;
use gneiss_rtk::estimators::rtk_iekf::DoubleDiffKey;
use gneiss_rtk::post_process::forward;
use gneiss_rtk::post_process::{execute_post_process, network, PostProcessOptions, ReceiverPcvPair, SmoothedEpoch};
use gneiss_rtk::swfg::config::EngineConfig;

const ROVER_FILE: &str = "p2241350.20o";
const TRUTH_FILE: &str = "p224_truth.pos";
const NAV_FILE: &str = "brdc1350.20n";

// --- 2025 multi-GNSS profile (GNEISS_DATASET=multi2025) -----------------
// Rover P224 + bases P181/P225/P222, all upgraded to multi-GNSS receivers;
// 20-observable mixed RINEX 2.11 exports (GPS+GLONASS+Galileo). Truth and
// coordinates from UNR IGS20 medians for 2025-06 (see
// scripts/gen_multignss_truth.py). Nav: RINEX 3.04 mixed broadcast from the
// BKG IGS mirror (decompress station_mixed_nav.rnx.gz first).
const M25_DIR: &str = "datasets/multignss_2025d160";
const M25_ROVER: &str = "p2241600.25o";
const M25_TRUTH: &str = "p224_truth.pos";
const M25_NAV: &str = "station_mixed_nav.rnx";
const M25_BASES: &[NetworkBase] = &[
    NetworkBase { id: "P181", base_file: "p1811600.25o", base_pos: Vector3::new(-2697941.1641, -4255089.2918, 3898009.6542), baseline_km: 14.97 },
    NetworkBase { id: "P225", base_file: "p2251600.25o", base_pos: Vector3::new(-2681518.9776, -4281621.6523, 3880440.4195), baseline_km: 21.86 },
    NetworkBase { id: "P222", base_file: "p2221600.25o", base_pos: Vector3::new(-2689640.4153, -4290437.2046, 3865050.9752), baseline_km: 37.97 },
];

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
    let up_p95 = up_errs[(n as f64 * 0.95) as usize];
    let up_rms = (up_errs.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();
    let fix_pct = (fix_count as f64 / total_count.max(1) as f64) * 100.0;
    println!("=== {} (N={}, Fixed={}/{} [{:.1}%]) ===", name, n, fix_count, total_count, fix_pct);
    println!("Horizontal Error:  p50={:.3}m,  p68={:.3}m,  p95={:.3}m,  RMS={:.3}m", p50, p68, p95, rms);
    // Guard scripts parse floats positionally off this line; new stats go
    // at the END only so existing indices stay valid.
    println!("Vertical Error:    p50={:+.3}m,  RMS={:.3}m,  p95={:.3}m", up_p50, up_rms, up_p95);
    println!("3D Position Error: p50={:.3}m,  p95={:.3}m", d3_errs[n / 2], d3_errs[(n as f64 * 0.95) as usize]);
    [p50, p68, p95, rms]
}

fn print_fixed_stats(name: &str, mut h_errs: Vec<f64>) -> f64 {
    if h_errs.is_empty() {
        return 0.0;
    }
    h_errs.sort_by(|a, b| a.total_cmp(b));
    let p50 = h_errs[h_errs.len() / 2];
    let n = h_errs.len();
    let p95 = h_errs[(n as f64 * 0.95) as usize];
    let rms = (h_errs.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();
    println!("  [{} fixed-only horizontal p50: {:.3}m, p95: {:.3}m, RMS: {:.3}m, N={}]", name, p50, p95, rms, n);
    p50
}

fn collect_errors(traj: &[SmoothedEpoch], truth: &Truth) -> (Vec<f64>, Vec<f64>, Vec<f64>, usize, Vec<f64>) {
    let mut h_errs = Vec::new();
    let mut d3_errs = Vec::new();
    let mut up_errs = Vec::new();
    let mut fixed_errs = Vec::new();
    let mut wrong_tows: Vec<u32> = Vec::new();
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
                    if h > 0.30 {
                        wrong_tows.push(tow);
                    }
                    fixed_errs.push(h);
                }
            }
        }
    }
    // Wrong-fix diagnostic: fixed epochs with gross position error.
    if !wrong_tows.is_empty() && std::env::var("WL_DIAG").is_ok() {
        let n_fixed = fixed_errs.len();
        eprintln!(
            "  [DIAG] wrong fixes (h>0.30m): {}/{} | TOW {}..{} (span {:.1} min)",
            wrong_tows.len(), n_fixed,
            wrong_tows[0], wrong_tows[wrong_tows.len()-1],
            (wrong_tows[wrong_tows.len()-1] - wrong_tows[0]) as f64 / 60.0
        );
    }
    (h_errs, d3_errs, up_errs, fix, fixed_errs)
}

struct RunContext<'a> {
    ephemerides: &'a [Ephemeris],
    truth: &'a Truth,
    rover_init: Option<Vector3<f64>>,
    klob: Option<([f64; 4], [f64; 4])>,
    tropo_gradients: bool,
}

#[allow(clippy::too_many_arguments)]

fn run_pass(
    config: &EngineConfig,
    ctx: &RunContext,
    rover: &[EpochObs],
    base_epochs: &[EpochObs],
    base: &NetworkBase,
    label: &str,
    bidir: bool,
    network_upd: Option<HashMap<u16, f64>>,
    base_pos_eff: Vector3<f64>,
    receiver_pcv: Option<std::sync::Arc<ReceiverPcvPair>>,
) -> ([f64; 4], Vec<SmoothedEpoch>) {
    // Receiver antenna PCV loading (gated by GNEISS_RECV_PCV=1).
    let options = PostProcessOptions {
        enable_bidirectional: bidir,
        base_position: Some(base_pos_eff),
        initial_rover_position: ctx.rover_init,
        klobuchar_alpha: ctx.klob.map(|k| k.0),
        klobuchar_beta: ctx.klob.map(|k| k.1),
        // Static monuments: q=1.0 re-randomizes position ~55 m per 30 s
        // epoch and keeps the float solution from converging. Loose phase
        // (1e-6) lets float ambiguities converge; the engine's two-phase
        // lock tightens to 1e-8 after 15 min so wrong fixes cannot hide.
        q_accel: Some(1e-6),
        network_sat_upd: network_upd.clone(),
        // Long baselines need iono-immune fixing: MW wide-lane cascade AR
        // unlocks the iono-free stage beyond ~20 km.
        widelane_ar: std::env::var("WL_DISABLE").is_err(),
        tropo_gradients: ctx.tropo_gradients,
        receiver_pcv,
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
    if std::env::var("WL_DUMP").is_ok() {
        // Full per-epoch error dump for cross-base correlation studies.
        let path = format!("/tmp/dump_{}_{}.csv", base.id, label);
        if let Ok(mut f) = std::fs::File::create(&path) {
            use std::io::Write;
            let _ = writeln!(f, "tow,h,v,q,sep,nsat");
            for ep in &traj {
                if let Some(&t) = ctx.truth.get(&(ep.time.tow.round() as u32)) {
                    let _ = writeln!(
                        f, "{:.0},{:.4},{:+.4},{},{:.3},{}",
                        ep.time.tow,
                        horizontal_error(ep.position_ecef, t),
                        vertical_error(ep.position_ecef, t),
                        ep.quality, ep.separation_3d, ep.n_satellites,
                    );
                }
            }
        }
    }
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

/// Antenna family + radome from a RINEX observation header (`ANT # / TYPE`),
/// radome defaulting to NONE when absent.
fn rinex_ant_type(rinex_path: &Path) -> Option<(String, String)> {
    use std::io::BufRead;
    let f = File::open(rinex_path).ok()?;
    for line in BufReader::new(f).lines().take(80).flatten() {
        if line.len() >= 60 && line[60..].trim() == "ANT # / TYPE" {
            let fields: Vec<&str> = line[20..40].split_whitespace().collect();
            let fam = fields.first()?.to_string();
            let rad = fields.get(1).copied().unwrap_or("NONE").to_string();
            return Some((fam, rad));
        }
    }
    None
}

#[allow(clippy::too_many_arguments)]
/// Receiver L1 PCO (ECEF, m) for a base station: RINEX header antenna
/// type/radome looked up in an ANTEX file. Opt-in via GNEISS_RECV_PCO=1
/// (file via GNEISS_ANTEX, default datasets/igs14.atx).
fn base_recv_pco_ecef(rinex_path: &Path, antex_path: &str, arp: Vector3<f64>) -> Option<Vector3<f64>> {
    use gneiss_parsers::antex::AntexDatabase;

    // 1. antenna type + radome from the RINEX2 header
    let (fam, rad) = rinex_ant_type(rinex_path)?;

    // 2. ANTEX receiver entry for that family/radome
    let db = match AntexDatabase::parse(antex_path) { Ok(d) => d, Err(e) => { eprintln!("DEBUG: ANTEX parse failed: {:?}", e); return None; } };
    let ant = gneiss_parsers::receiver_antenna::ReceiverAntenna::lookup(&db, &fam, &rad)?;
    // Parser reorders the ANTEX north/east/up columns to east/north/up.
    let [east_mm, north_mm, up_mm] = ant.pco_enu_mm;

    // 3. ENU -> ECEF at the ARP
    let llh = gneiss_core::coords::ecef_to_llh(arp);
    let (lat, lon) = (llh.x, llh.y);
    let (slat, clat) = (lat.sin(), lat.cos());
    let (slon, clon) = (lon.sin(), lon.cos());
    let east = Vector3::new(-slon, clon, 0.0);
    let north = Vector3::new(-slat * clon, -slat * slon, clat);
    let up = Vector3::new(clat * clon, clat * slon, slat);
    let m = 1e-3;
    Some(north * (north_mm * m) + east * (east_mm * m) + up * (up_mm * m))
}

/// Receiver antenna PCV models (rover, base) from the ANTEX database,
/// opt-in via GNEISS_PCV=1. Both headers must resolve to calibrations;
/// otherwise None keeps the legacy uncorrected path.
fn load_receiver_pcv(rover_path: &Path, base_path: &Path, antex_path: &str) -> Option<std::sync::Arc<ReceiverPcvPair>> {
    use gneiss_parsers::antex::AntexDatabase;
    use gneiss_parsers::receiver_antenna::ReceiverAntenna;

    let db = match AntexDatabase::parse(antex_path) { Ok(d) => d, Err(e) => { eprintln!("DEBUG: ANTEX parse failed: {:?}", e); return None; } };
    let (rfam, rrad) = match rinex_ant_type(rover_path) { Some(x) => x, None => { eprintln!("DEBUG: ant type not found in {}", rover_path.display()); return None; } };
    let rover = ReceiverAntenna::lookup(&db, &rfam, &rrad)?;
    let (bfam, brad) = rinex_ant_type(base_path)?;
    let base = ReceiverAntenna::lookup(&db, &bfam, &brad)?;
    println!(
        "RECV-PCV enabled: rover [{}] base [{}] ({})",
        rover.antenna_type(),
        base.antenna_type(),
        antex_path
    );
    Some(std::sync::Arc::new(ReceiverPcvPair {
        rover: std::sync::Arc::new(rover),
        base: std::sync::Arc::new(base),
    }))
}

fn run_base(
    base: &NetworkBase,
    dir: &Path,
    rover_file: &str,
    ctx: &RunContext,
    rover: &[EpochObs],
    network_upd: Option<HashMap<u16, f64>>,
) -> ([f64; 8], Vec<SmoothedEpoch>) {
    // Receiver-side PCO: shift the base reference point from ARP to its L1
    // phase centre so DD geometry references real antenna positions. The
    // rover stays ARP-referenced (truth datum), so this corrects the
    // cross-family differential without introducing an evaluation offset.
    let recv_pco_on = std::env::var("GNEISS_RECV_PCO").is_ok();
    let base_pos_eff = if recv_pco_on {
        // Differential datum correction: subtract the rover antenna's own
        // PCO so same-family baselines stay put and only cross-family
        // mismatches (e.g. CAPO's Leica +39.7 mm Up) move. Rover ARP
        // anchor = first truth epoch.
        let rover_arp = ctx.truth.values().next().copied();
        let antex = std::env::var("GNEISS_ANTEX")
            .unwrap_or_else(|_| "datasets/igs14.atx".into());
        match (
            base_recv_pco_ecef(&dir.join(base.base_file), &antex, base.base_pos),
            rover_arp.and_then(|arp| {
                base_recv_pco_ecef(Path::new("datasets/cors_short_baseline/p2241350.20o"), &antex, arp)
            }),
        ) {
            (Some(d_base), Some(d_rover)) => {
                let diff = d_base - d_rover;
                println!(
                    "RECV-PCO [{}]: differential |d|={:.1} mm",
                    base.id,
                    diff.norm() * 1000.0
                );
                base.base_pos + diff
            }
            _ => {
                eprintln!("RECV-PCO [{}]: lookup failed; using ARP", base.id);
                base.base_pos
            }
        }
    } else {
        base.base_pos
    };
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
    // Opt-in elevation-dependent receiver PCV correction (GNEISS_PCV=1):
    // strips the differential antenna signature from every DD phase.
    let receiver_pcv = if std::env::var("GNEISS_PCV").is_ok() {
        eprintln!("DEBUG: GNEISS_PCV detected");
        let antex = std::env::var("GNEISS_ANTEX")
            .unwrap_or_else(|_| "datasets/igs14.atx".into());
        load_receiver_pcv(&dir.join(rover_file), &dir.join(base.base_file), &antex)
    } else {
        None
    };
    let config = EngineConfig::Rtk(gneiss_rtk::swfg::config::RtkConfig {
        initial_position: Some([base_pos_eff.x, base_pos_eff.y, base_pos_eff.z]),
        ..Default::default()
    });
    let (fwd, _fwd_traj) = run_pass(&config, ctx, rover, &base_epochs, base, "Forward", false, network_upd.clone(), base_pos_eff, receiver_pcv.clone());
    let (smooth, smooth_traj) = run_pass(&config, ctx, rover, &base_epochs, base, "Smoothed", true, network_upd.clone(), base_pos_eff, receiver_pcv);
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

fn load_dataset(dir: &Path, rover_file: &str, truth_file: &str, nav_file: &str) -> Result<Dataset, String> {
    if !dir.join(rover_file).exists() {
        return Err(format!(
            "Dataset missing {} in {} (fetch/decompress first)",
            rover_file,
            dir.display()
        ));
    }
    let nav_f = File::open(dir.join(nav_file)).map_err(|e| format!("Failed to open nav: {}", e))?;
    let (ephemerides, klobuchar) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_f))
        .map_err(|e| format!("Failed to parse nav: {}", e))?;
    let rov_f = File::open(dir.join(rover_file)).map_err(|e| format!("Failed to open rover obs: {}", e))?;
    let (rover_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(rov_f))
        .map_err(|e| format!("Failed to parse rover obs: {}", e))?;
    // System filter for controlled experiments: GNEISS_SYSTEMS is a set of
    // constellation letters (e.g. "G" or "GE"). Default G preserves the
    // historical GPS-only behaviour on every dataset.
    let allowed: Vec<char> = std::env::var("GNEISS_SYSTEMS")
        .unwrap_or_else(|_| "G".into())
        .chars()
        .collect();
    let keep = |c: gneiss_core::sat::Constellation| match c {
        gneiss_core::sat::Constellation::Gps => allowed.contains(&'G'),
        gneiss_core::sat::Constellation::Glonass => allowed.contains(&'R'),
        gneiss_core::sat::Constellation::Galileo => allowed.contains(&'E'),
        gneiss_core::sat::Constellation::Beidou => allowed.contains(&'C'),
        _ => false,
    };
    let n_before: usize = rover_epochs.iter().map(|e| e.satellites.len()).sum();
    let rover_epochs: Vec<EpochObs> = rover_epochs
        .into_iter()
        .map(|mut e| {
            e.satellites.retain(|s| keep(s.sat.constellation));
            e
        })
        .collect();
    let n_after: usize = rover_epochs.iter().map(|e| e.satellites.len()).sum();
    if n_before != n_after {
        println!("SYSTEM FILTER: {}/{} rover satellites kept ({})", n_after, n_before, allowed.iter().collect::<String>());
    }
    let truth = parse_truth(&dir.join(truth_file));
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

    let multi2025 = std::env::var("GNEISS_DATASET").as_deref() == Ok("multi2025");
    let (dir, rover_file, truth_file, nav_file, bases): (&Path, &str, &str, &str, &[NetworkBase]) =
        if multi2025 {
            (Path::new(M25_DIR), M25_ROVER, M25_TRUTH, M25_NAV, M25_BASES)
        } else {
            (Path::new("datasets/cors_short_baseline"), ROVER_FILE, TRUTH_FILE, NAV_FILE, BASES)
        };
    let data = match load_dataset(dir, rover_file, truth_file, nav_file) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1);
        }
    };

    // Rover filter init is left to the SPP in broadcast frame (None): the
    // RINEX header approx is NAD83 and must not seed the filter.
    let tropo_gradients = multi2025 && std::env::var("GNEISS_TROPO_GRAD").as_deref() != Ok("0");
    let ctx = RunContext { ephemerides: &data.ephemerides, truth: &data.truth, rover_init: None, klob: data.klob, tropo_gradients };
    let selected_rover = select_rover_epochs(&data.rover_epochs);
    let mut results = Vec::new();
    let mut base_trajs: Vec<Vec<SmoothedEpoch>> = Vec::new();
    // Phase A: per-base forward passes collecting phase-only wide-lane
    // arc means; the cross-base least-squares decomposition solves the
    // satellite wide-lane UPDs that make MW rounding trustworthy.
    let only = std::env::var("WL_ONLY_BASE").ok();

    // Phase A: per-base forward passes collecting phase-only wide-lane
    // arc means; cross-base least squares solves satellite wide-lane UPDs
    // that Phase B applies before wide-lane rounding.
    let mut network_upd: Option<HashMap<u16, f64>> = None;
    if std::env::var("WL_NO_UPD").is_err() {
        let mut per_base: Vec<HashMap<DoubleDiffKey, f64>> = Vec::new();
        for base in bases {
            if let Some(want) = &only { if base.id != want.as_str() { continue; } }
            let Ok(f) = File::open(dir.join(base.base_file)) else { continue };
            let Ok((base_epochs, _)) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(f)) else { continue };
            let (_, _wl, pw) = forward::run_forward_pass_collecting(
                ctx.ephemerides, selected_rover, &base_epochs,
                base.base_pos, 1e-6,
            );
            let means: HashMap<DoubleDiffKey, f64> = pw
                .arc_means()
                .into_iter()
                .map(|(k, (m, _))| (k, m))
                .collect();
            println!("UPD pre-pass [{}]: {} converged pairs", base.id, means.len());
            per_base.push(means);
        }
        let sol = gneiss_rtk::estimators::rtk_iekf::mw::solve_network_upd(&per_base);
        println!(
            "Network WL UPD solution (residual RMS {:.3} cyc, {} sats):",
            sol.residual_rms, sol.sat_upd.len()
        );
        let mut rows: Vec<_> = sol.sat_upd.iter().collect();
        rows.sort_by_key(|(s, _)| **s);
        for (s, u) in rows {
            println!("  G{:02}: {:+.3}", s, u);
        }
        network_upd = Some(sol.sat_upd);
    }

    for base in bases {
        if let Some(want) = &only { if base.id != want.as_str() { continue; } }
        let (stats, traj) = run_base(base, dir, rover_file, &ctx, selected_rover, network_upd.clone());
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
