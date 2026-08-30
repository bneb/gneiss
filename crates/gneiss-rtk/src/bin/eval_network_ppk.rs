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
use rayon::prelude::*;

use nalgebra::Vector3;

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::metrics::{horizontal_error, vertical_error};
use gneiss_core::obs::EpochObs;
use gneiss_core::time::GpsTime;
use gneiss_rtk::estimators::rtk_iekf::DoubleDiffKey;
use gneiss_rtk::post_process::dynamics::{KINEMATIC_Q_ACCEL, STATIC_Q_ACCEL, ProcessingDynamics};
use gneiss_rtk::post_process::forward;
use gneiss_rtk::post_process::{execute_post_process, network, sidereal, PostProcessOptions, ReceiverPcvPair, SmoothedEpoch};
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

/// Truth coordinates in the UNR/IGS20 reference frame. Type-tagged so
/// that comparison against broadcast-frame solutions requires an
/// explicit conversion (FRAME_SAFETY_PLAN ledger row 6).
type Truth = BTreeMap<u32, gneiss_core::frames::EcefPos<gneiss_core::frames::Igs20>>;

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
        truth.insert(
            gps_time.tow.round() as u32,
            gneiss_core::frames::EcefPos::new(Vector3::new(px, py, pz)),
        );
    }
    truth
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

/// Express an IGS20 truth coordinate in the solution (broadcast) frame.
///
/// Both frames are ITRF2020-aligned to first order, so the Helmert is
/// near-identity; the documented ~50 mm vertical residual between UNR
/// IGS20 medians and broadcast-frame solutions is NOT captured by these
/// parameters and remains as a systematic offset in reported errors.
fn truth_in_solution_frame(
    t: gneiss_core::frames::EcefPos<gneiss_core::frames::Igs20>,
    epoch_yr: f64,
) -> Vector3<f64> {
    *t.convert_to::<gneiss_core::frames::Wgs84Broadcast>(epoch_yr)
        .vector()
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
            let t = truth_in_solution_frame(t, 2025.5);
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
    /// Rover motion model (GNEISS_DYNAMICS=kinematic opts in).
    dynamics: ProcessingDynamics,
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
    // receiver_pcv is already Some/None per the caller's GNEISS_PCV=1 gate
    // (load_receiver_pcv); the engine applies it whenever it's Some.
    // Long baselines need iono-immune fixing: MW wide-lane cascade AR
    // unlocks the iono-free stage beyond ~20 km.
    let widelane_ar = std::env::var("WL_DISABLE").is_err();
    let options = PostProcessOptions {
        enable_bidirectional: bidir,
        base_position: Some(base_pos_eff),
        initial_rover_position: ctx.rover_init,
        klobuchar_alpha: ctx.klob.map(|k| k.0),
        klobuchar_beta: ctx.klob.map(|k| k.1),
        // Static monuments (legacy): q=1e-6 — 1.0 re-randomizes position
        // ~55 m per 30 s epoch and keeps the float solution from
        // converging; the engine's two-phase lock tightens to 1e-8 after
        // 15 min so wrong fixes cannot hide. Kinematic profile: q=1.0
        // (correct mobile prior) and no monument lock.
        q_accel: Some(if ctx.dynamics.is_kinematic() {
            KINEMATIC_Q_ACCEL
        } else {
            STATIC_Q_ACCEL
        }),
        network_sat_upd: network_upd.clone(),
        widelane_ar,
        tropo_gradients: ctx.tropo_gradients,
        receiver_pcv,
        dynamics: ctx.dynamics,
        // Honesty gating only makes sense once the bidirectional fuse has
        // run; widelane_ar is this profile's signal that long-baseline
        // iono-free products are in play.
        continuity_gate: bidir && widelane_ar,
        // Moved from inside run_forward_iekf/run_backward_iekf, where it
        // was nested in a baseline-length gate meant for ZWD state (so
        // >25km baselines never got GLONASS regardless of this env var)
        // and, in the backward pass, gated differently than forward.
        // Same env var, now applies uniformly to both passes and all
        // baseline lengths -- see PostProcessOptions::enable_glonass.
        enable_glonass: std::env::var("GNEISS_GLONASS").is_ok(),
    };
    let res = match execute_post_process(config, ctx.ephemerides, rover, Some(base_epochs), None, &options) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("{} pass failed for {}: {}", label, base.id, e);
            return ([0.0; 4], Vec::new());
        }
    };
    // Honesty gating (excursion jumps beyond what the motion model
    // allows downgraded to float) now happens inside execute_post_process
    // itself via options.continuity_gate, set above.
    let mut traj = res.trajectory;
    // Sidereal stacking (GNEISS_SIDEREAL=1): diagnostic fold + strictly
    // causal first-half mitigation on the smoothed product, before any
    // error collection. Default OFF keeps the legacy path bit-identical.
    if std::env::var("GNEISS_SIDEREAL").is_ok() {
        traj = apply_sidereal_option(traj, ctx, base.id, label);
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
                    let t = truth_in_solution_frame(t, 2025.5);
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
                let t = truth_in_solution_frame(t, 2025.5);
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

/// Sidereal stacking option (GNEISS_SIDEREAL=1): fold truth-referenced
/// E/N/U errors of the smoothed trajectory onto sidereal phase bins,
/// print the per-channel verdict, dump a per-bin CSV, and subtract
/// causally-fitted first-half mean corrections from second-half epochs.
/// Knobs: GNEISS_SIDEREAL_BINS (default 240), GNEISS_SIDEREAL_MINCNT
/// (default 2), GNEISS_SIDEREAL_DIAG_DIR (default /tmp).
fn apply_sidereal_option(
    traj: Vec<SmoothedEpoch>,
    ctx: &RunContext,
    base_id: &str,
    label: &str,
) -> Vec<SmoothedEpoch> {
    let n_bins = std::env::var("GNEISS_SIDEREAL_BINS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(sidereal::DEFAULT_BINS);
    let min_cnt = std::env::var("GNEISS_SIDEREAL_MINCNT")
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .unwrap_or(sidereal::MIN_BIN_COUNT_DEFAULT);
    let truth_at = |tow: u32| ctx.truth.get(&tow).map(|t| truth_in_solution_frame(*t, 2025.5));
    let target = format!("{label}/{base_id}");
    let (traj, report) = sidereal::apply_to_trajectory(traj, truth_at, n_bins, min_cnt);
    println!("{}", report.summary(&target));
    let dir = std::env::var("GNEISS_SIDEREAL_DIAG_DIR").unwrap_or_else(|_| "/tmp".into());
    let path = Path::new(&dir).join(format!("sidereal_diag_{base_id}_{label}.csv"));
    match sidereal::write_diag_csv(&path, &target, &report.channels, n_bins, report.split_index) {
        Ok(_) => println!("SIDEREAL CSV -> {}", path.display()),
        Err(e) => eprintln!("SIDEREAL CSV write failed {}: {}", path.display(), e),
    }
    traj
}

use gneiss_rtk::post_process::antenna::{load_receiver_pcv, station_recv_pco_ecef};

fn run_base(
    base: &NetworkBase,
    dir: &Path,
    rover_file: &str,
    ctx: &RunContext,
    rover: &[EpochObs],
    base_epochs: &[EpochObs],
    network_upd: Option<HashMap<u16, f64>>,
) -> ([f64; 8], Vec<SmoothedEpoch>) {
    let recv_pco_on = std::env::var("GNEISS_RECV_PCO_DISABLE").is_err();
    let base_pos_eff = if recv_pco_on {
        let rover_arp = ctx.truth.values().next().copied();
        let antex = std::env::var("GNEISS_ANTEX")
            .unwrap_or_else(|_| "datasets/igs14.atx".into());
        match (
            station_recv_pco_ecef(&dir.join(base.base_file), &antex, base.base_pos),
            rover_arp.and_then(|arp| {
                let v = truth_in_solution_frame(arp, 2025.5);
                station_recv_pco_ecef(Path::new("datasets/cors_short_baseline/p2241350.20o"), &antex, v)
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
    let (fwd, _fwd_traj) = run_pass(&config, ctx, rover, base_epochs, base, "Forward", false, network_upd.clone(), base_pos_eff, receiver_pcv.clone());
    let (smooth, smooth_traj) = run_pass(&config, ctx, rover, base_epochs, base, "Smoothed", true, network_upd.clone(), base_pos_eff, receiver_pcv);
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


/// Apply Hatch filter to rover L1 pseudoranges in-place.
fn apply_hatch_filter(rover_epochs: &mut [EpochObs], window: usize) {
    use gneiss_core::hatch::HatchFilter;
    let mut hf = HatchFilter::new(window);
    for epoch in rover_epochs.iter_mut() {
        for sat_obs in epoch.satellites.iter_mut() {
            let lambda = gneiss_core::frequencies::frequency_for(
                sat_obs.sat.constellation,
                gneiss_core::frequencies::Signal::GpsL1Ca,
                1,
            ) / 299_792_458.0;
            if let Some(pr_val) = sat_obs.get_observable(1) {
                if let Some(cp_val) = sat_obs.get_observable_phase(1) {
                    let smoothed = hf.apply(sat_obs.sat, pr_val, cp_val, lambda);
                    if let Some(o) = sat_obs.observations.iter_mut().find(|o| {
                        o.code.obs_type == gneiss_core::obs::ObsType::Pseudorange
                            && o.code.signal.freq_band == 1
                    }) {
                        o.value = smoothed;
                    }
                }
            }
        }
    }
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
    let allowed = std::env::var("GNEISS_SYSTEMS").unwrap_or_else(|_| "G".into());
    let n_before: usize = rover_epochs.iter().map(|e| e.satellites.len()).sum();
    let mut rover_epochs = rover_epochs;
    gneiss_core::obs::filter_constellations(&mut rover_epochs, &allowed);
    let n_after: usize = rover_epochs.iter().map(|e| e.satellites.len()).sum();
    if n_before != n_after {
        println!("SYSTEM FILTER: {}/{} rover satellites kept ({})", n_after, n_before, allowed);
    }
    // Carrier-smoothed-code (Hatch) filter, opt-in via GNEISS_HATCH=N.
    if let Ok(win) = std::env::var("GNEISS_HATCH") {
        if let Ok(n) = win.parse::<usize>() {
            if n > 0 {
                apply_hatch_filter(&mut rover_epochs, n);
                println!("HATCH: {}-epoch carrier-smoothed code applied", n);
            }
        }
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

    // Processing-dynamics profile: GNEISS_DYNAMICS=kinematic opts in;
    // unset keeps the legacy static-monument path byte-identical.
    let dynamics = ProcessingDynamics::from_env();
    if dynamics.is_kinematic() {
        println!(
            "DYNAMICS: kinematic (q_accel={KINEMATIC_Q_ACCEL}, no monument lock, \
             innov gate x{:.0}, sigma-scaled combiner)",
            dynamics.innovation_gate_scale(),
        );
    }

    let multi2025 = std::env::var("GNEISS_DATASET").as_deref() == Ok("multi2025");
    // Generic acquisition override: point at any directory laid out like
    // multignss_2025d160 (P224 rover + truth + mixed nav + per-base obs).
    // Base filenames are expected under their DOY160 names — acquisition
    // scripts provide day->canonical symlinks.
    let custom_dir = std::env::var("GNEISS_DATA_DIR").ok();
    let (dir, rover_file, truth_file, nav_file, bases): (&Path, &str, &str, &str, &[NetworkBase]) =
        if let Some(d) = &custom_dir {
            (Path::new(d), M25_ROVER, M25_TRUTH, M25_NAV, M25_BASES)
        } else if multi2025 {
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
    let ctx = RunContext { ephemerides: &data.ephemerides, truth: &data.truth, rover_init: None, klob: data.klob, tropo_gradients, dynamics };
    let selected_rover = select_rover_epochs(&data.rover_epochs);
    let mut results = Vec::new();
    let mut base_trajs: Vec<Vec<SmoothedEpoch>> = Vec::new();
    // Phase A: per-base forward passes collecting phase-only wide-lane
    // arc means; the cross-base least-squares decomposition solves the
    // satellite wide-lane UPDs that make MW rounding trustworthy.
    let only = std::env::var("WL_ONLY_BASE").ok();

    let active_bases: Vec<_> = bases
        .iter()
        .filter(|b| only.as_ref().is_none_or(|want| b.id == want.as_str()))
        .collect();

    // Preload all active base observation files once in parallel
    let loaded_base_obs: HashMap<&str, std::sync::Arc<Vec<EpochObs>>> = active_bases
        .par_iter()
        .filter_map(|base| {
            let f = File::open(dir.join(base.base_file)).ok()?;
            let (mut epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(f)).ok()?;
            let allowed = std::env::var("GNEISS_SYSTEMS").unwrap_or_else(|_| "G".into());
            gneiss_core::obs::filter_constellations(&mut epochs, &allowed);
            Some((base.id, std::sync::Arc::new(epochs)))
        })
        .collect();

    // Phase A: per-base forward passes collecting phase-only wide-lane
    // arc means; cross-base least squares solves satellite wide-lane UPDs
    // that Phase B applies before wide-lane rounding.
    let mut network_upd: Option<HashMap<u16, f64>> = None;
    if std::env::var("WL_NO_UPD").is_err() {
        let per_base_means: Vec<_> = active_bases.par_iter().filter_map(|base| {
            let base_epochs = loaded_base_obs.get(base.id)?;
            let (_, _wl, pw) = forward::run_forward_pass_collecting(
                ctx.ephemerides, selected_rover, base_epochs,
                base.base_pos, 1e-6,
            );
            let means: HashMap<DoubleDiffKey, f64> = pw
                .arc_means()
                .into_iter()
                .map(|(k, (m, _))| (k, m))
                .collect();
            Some((base.id, means))
        }).collect();

        let mut per_base: Vec<HashMap<DoubleDiffKey, f64>> = Vec::new();
        for (base_id, means) in &per_base_means {
            println!("UPD pre-pass [{}]: {} converged pairs", base_id, means.len());
            per_base.push(means.clone());
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

    let base_results: Vec<_> = active_bases
        .par_iter()
        .filter_map(|base| {
            let base_epochs = loaded_base_obs.get(base.id)?;
            let (stats, traj) = run_base(base, dir, rover_file, &ctx, selected_rover, base_epochs, network_upd.clone());
            Some((*base, stats, traj))
        })
        .collect();

    for (base, stats, traj) in base_results {
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
        // Network consensus continuity gate: profile-aware (static
        // 0.20 m rule vs kinematic velocity+sigma allowance).
        let fused_dyn = if dynamics.is_kinematic() { dynamics } else { ProcessingDynamics::Static };
        let gated = network::apply_continuity_gate_dynamics(fused, 90.0, fused_dyn);
        let (h, d3, up_errs, fix, fixed_errs) = collect_errors(&gated, ctx.truth);
        print_stats(&format!("NETWORK FUSED [{} bases]", base_trajs.len()), h, d3, up_errs, fix, gated.len());
        print_fixed_stats("Network", fixed_errs);
    }
}
