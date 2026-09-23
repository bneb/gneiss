use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use nalgebra::Vector3;

use gneiss_rtk::swfg::imu_preintegration::ImuSample;

#[derive(Clone, Copy)]
pub struct ImuRecord {
    pub sample: ImuSample,
    pub speed: f64,
}

pub fn parse_imu_line(line: &str) -> Option<ImuRecord> {
    let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
    if parts.len() < 8 {
        return None;
    }
    let tow: f64 = parts[0].parse().ok()?;
    let ax: f64 = parts[2].parse().ok()?;
    let ay: f64 = parts[3].parse().ok()?;
    let az: f64 = parts[4].parse().ok()?;
    let gx: f64 = parts[5].parse().ok()?;
    let gy: f64 = parts[6].parse().ok()?;
    let gz: f64 = parts[7].parse().ok()?;
    let speed = parts.get(8).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);

    let time_us = (tow * 1_000_000.0).round() as u64;
    Some(ImuRecord {
        sample: ImuSample {
            accel: Vector3::new(ax, ay, az),
            gyro: Vector3::new(gx, gy, gz),
            time_us,
        },
        speed,
    })
}

pub fn parse_imu_csv(path: &Path) -> Vec<ImuRecord> {
    let file = File::open(path).expect("failed to open imu.csv");
    BufReader::new(file)
        .lines()
        .filter_map(|l| parse_imu_line(&l.ok()?))
        .collect()
}

pub type GnssFixMap = BTreeMap<u32, (Vector3<f64>, usize, bool)>;

pub fn load_cached_gnss_fixes(path: &Path) -> Option<GnssFixMap> {
    let content = std::fs::read_to_string(path).ok()?;
    let mut map = BTreeMap::new();
    for line in content.lines().skip(1) {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() >= 6 {
            let k: u32 = parts[0].trim().parse().ok()?;
            let x: f64 = parts[1].trim().parse().ok()?;
            let y: f64 = parts[2].trim().parse().ok()?;
            let z: f64 = parts[3].trim().parse().ok()?;
            let ns: usize = parts[4].trim().parse().ok()?;
            let fixed: bool = parts[5].trim().parse().ok()?;
            map.insert(k, (Vector3::new(x, y, z), ns, fixed));
        }
    }
    if map.len() > 1000 { Some(map) } else { None }
}

pub fn save_cached_gnss_fixes(path: &Path, solutions: &GnssFixMap) {
    if let Ok(mut f) = File::create(path) {
        let _ = writeln!(f, "epoch_key,x,y,z,ns,fixed");
        for (k, (pos, ns, fixed)) in solutions {
            let _ = writeln!(f, "{},{},{},{},{},{}", k, pos.x, pos.y, pos.z, ns, fixed);
        }
    }
}

pub fn parse_truth_line(line: &str) -> Option<(u32, (f64, f64, f64))> {
    let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
    if parts.len() < 8 {
        return None;
    }
    let tow: f64 = parts[0].parse().ok()?;
    let x: f64 = parts[5].parse().ok()?;
    let y: f64 = parts[6].parse().ok()?;
    let z: f64 = parts[7].parse().ok()?;
    Some(((tow * 10.0).round() as u32, (x, y, z)))
}

pub fn load_truth(path: &Path) -> BTreeMap<u32, (f64, f64, f64)> {
    let ref_csv = std::fs::read_to_string(path).expect("read reference.csv");
    ref_csv.lines().skip(1).filter_map(parse_truth_line).collect()
}

pub fn compute_base_pos(approx: Option<[f64; 3]>) -> Vector3<f64> {
    let default_base_arp = Vector3::new(-3961904.4341, 3348994.266, 3698211.7067);
    let default_base_llh = gneiss_core::coords::ecef_to_llh(default_base_arp);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(default_base_llh).transpose();
    let computed_base = default_base_arp + ned_to_ecef * Vector3::new(0.0, 0.0, 0.0855);
    approx.map(|p| Vector3::new(p[0], p[1], p[2])).unwrap_or(computed_base)
}

pub struct RinexInputs {
    pub rover_epochs: Vec<gneiss_core::obs::EpochObs>,
    pub base_index: BTreeMap<u32, gneiss_core::obs::EpochObs>,
    pub base_pos: Vector3<f64>,
    pub ephems: Vec<gneiss_core::ephemeris::Ephemeris>,
    pub klob: Option<([f64; 4], [f64; 4])>,
}

pub fn load_rinex_inputs(dataset: &Path) -> RinexInputs {
    let nav_file = File::open(dataset.join("base.nav")).expect("open base.nav");
    let (ephems, klob) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_file)).expect("parse nav");
    let rover_file = File::open(dataset.join("rover_trimble.obs")).expect("open rover_trimble.obs");
    let (rover_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(rover_file)).expect("parse rover");
    let base_file = File::open(dataset.join("base_trimble.obs")).expect("open base_trimble.obs");
    let (base_epochs, base_approx) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(base_file)).expect("parse base");
    let base_pos = compute_base_pos(base_approx.approx_position);
    let base_index = base_epochs.into_iter().map(|e| ((e.time.tow * 1000.0).round() as u32, e)).collect();
    RinexInputs { rover_epochs, base_index, base_pos, ephems, klob: klob.map(|k| (k.alpha, k.beta)) }
}

pub fn compute_horizontal_errors(
    positions: &BTreeMap<u32, Vector3<f64>>,
    truth: &BTreeMap<u32, (f64, f64, f64)>,
) -> (Vec<f64>, Vec<f64>) {
    let mut sorted = Vec::new();
    let mut chrono = Vec::new();
    for (tow, p) in positions {
        if let Some(&(tx, ty, tz)) = truth.get(tow) {
            let t_pos = Vector3::new(tx, ty, tz);
            let enu = gneiss_core::coords::ecef_delta_to_enu(*p, t_pos, gneiss_core::coords::ecef_to_llh(t_pos));
            let e = (enu.x * enu.x + enu.y * enu.y).sqrt();
            sorted.push(e);
            chrono.push(e);
        }
    }
    sorted.sort_by(|a, b| a.total_cmp(b));
    (sorted, chrono)
}

pub fn print_quartile_stats(chrono_errs: &[f64]) {
    if chrono_errs.len() <= 1000 { return; }
    let q = chrono_errs.len() / 4;
    for i in 0..4 {
        let mut s = chrono_errs[i * q..(i + 1) * q].to_vec();
        s.sort_by(|a, b| a.total_cmp(b));
        let c_rms = (s.iter().map(|e| e * e).sum::<f64>() / s.len() as f64).sqrt();
        println!("  Q{}: p50={:.3}m, RMS={:.3}m", i + 1, s[s.len() / 2], c_rms);
    }
}

pub fn print_trajectory_stats(
    name: &str,
    positions: &BTreeMap<u32, Vector3<f64>>,
    truth: &BTreeMap<u32, (f64, f64, f64)>,
) {
    let (h_errs, chrono_errs) = compute_horizontal_errors(positions, truth);
    if h_errs.is_empty() { return; }
    let n = h_errs.len();
    let p50 = h_errs[n / 2];
    let p68 = h_errs[(n as f64 * 0.68) as usize];
    let p95 = h_errs[(n as f64 * 0.95) as usize];
    let rms = (h_errs.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();
    println!("=== {} (N={}) ===", name, n);
    println!("Horizontal error: p50={:.3}m, p68={:.3}m, p95={:.3}m, RMS={:.3}m", p50, p68, p95, rms);
    print_quartile_stats(&chrono_errs);
}

pub fn print_evaluation_summary(
    gnss_positions: &BTreeMap<u32, Vector3<f64>>,
    forward_map: &BTreeMap<u32, Vector3<f64>>,
    smoothed_map: &BTreeMap<u32, Vector3<f64>>,
    truth: &BTreeMap<u32, (f64, f64, f64)>,
) {
    println!("Solutions count: GNSS={}, Forward={}, Smoothed={}", gnss_positions.len(), forward_map.len(), smoothed_map.len());
    print_trajectory_stats("GNSS-Only RTK (Raw Fixes)", gnss_positions, truth);
    print_trajectory_stats("Forward Inertial Filter", forward_map, truth);
    print_trajectory_stats("RTS Smoothed GNSS/INS", smoothed_map, truth);
    let smoothed_at_gnss: BTreeMap<u32, Vector3<f64>> = gnss_positions.keys()
        .filter_map(|k| smoothed_map.get(k).map(|&p| (*k, p)))
        .collect();
    print_trajectory_stats("RTS Smoothed (at GNSS Epochs)", &smoothed_at_gnss, truth);
}
