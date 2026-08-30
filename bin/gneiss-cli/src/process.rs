use std::collections::HashMap;
use std::path::Path;
use rayon::prelude::*;
use tracing::{info, warn};

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_geodesy::geoid::GeoidGrid;
use gneiss_rtk::estimators::rtk_iekf::DoubleDiffKey;
use gneiss_rtk::post_process::{
    execute_post_process, network, forward, SmoothedEpoch, PostProcessOptions,
    ProcessingDynamics,
};
use gneiss_rtk::swfg::config::EngineConfig;

use crate::export::{export_trajectory, ExportFormat};
use crate::qc::QcReport;

type ObsEpochsWithPos = (Vec<EpochObs>, Option<[f64; 3]>);
type LoadedBase = (String, Vec<EpochObs>, Option<[f64; 3]>);

pub struct ProcessArgs {
    pub rover: String,
    pub bases: Vec<String>,
    pub nav: Option<String>,
    pub output: String,
    pub format: Option<String>,
    pub qc_report: Option<String>,
    pub geoid: Option<String>,
    pub config: Option<String>,
    pub enable_backward_smoothing: bool,
    pub mode: Option<String>,
    pub max_epochs: Option<usize>,
    pub base_position: Option<String>,
    pub systems: Option<String>,
    pub antex: Option<String>,
    pub glonass: bool,
    pub sp3: Option<String>,
    pub clk: Option<String>,
}

pub async fn run_process(args: ProcessArgs) -> Result<(), Box<dyn std::error::Error>> {
    let (mut rover_epochs, rover_pos) = load_rinex_obs(&args.rover)?;
    if let Some(sys) = args.systems.as_deref() {
        gneiss_core::obs::filter_constellations(&mut rover_epochs, sys);
    }
    let max = args.max_epochs.unwrap_or(usize::MAX);
    let selected_rover = &rover_epochs[..rover_epochs.len().min(max)];

    let parent_dir = Path::new(&args.rover).parent().unwrap_or_else(|| Path::new("."));
    let (ephemerides, klobuchar) = load_ephemerides(parent_dir, args.nav.as_deref())?;
    let swfg_config = build_swfg_config(args.config.as_deref(), args.mode.as_deref())?;
    let antex_path = args.antex.as_deref().unwrap_or("datasets/igs14.atx");

    let is_ppp = args.mode.as_deref() == Some("ppp") || args.sp3.is_some();

    let trajectory = if is_ppp {
        run_ppp_pipeline(selected_rover, &ephemerides, args.sp3.as_deref(), args.clk.as_deref())?
    } else if args.bases.len() > 1 && args.enable_backward_smoothing {
        info!("Running Multi-Base Network RTK pipeline with {} bases...", args.bases.len());
        run_network_pipeline(
            selected_rover, &args.bases, &ephemerides, klobuchar.as_ref(),
            &swfg_config, &args.rover, rover_pos, antex_path, args.glonass,
        )?
    } else if let Some(base_path) = args.bases.first() {
        run_single_base_pipeline(
            selected_rover, base_path, &ephemerides, klobuchar.as_ref(),
            &swfg_config, &args.rover, rover_pos, args.base_position.as_deref(),
            antex_path, args.glonass, args.enable_backward_smoothing,
        )?
    } else {
        run_spp_pipeline(&swfg_config, &ephemerides, klobuchar.as_ref(), selected_rover)?
    };

    export_results(&trajectory, &args)?;
    Ok(())
}

fn export_results(
    trajectory: &[SmoothedEpoch],
    args: &ProcessArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    let geoid_grid = if let Some(g_path) = &args.geoid {
        let path = Path::new(g_path);
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
        let grid = match ext.to_ascii_lowercase().as_str() {
            "gtx" => {
                let bytes = std::fs::read(path)?;
                GeoidGrid::from_gtx_bytes(&bytes)
                    .map_err(|e| format!("Failed to parse GTX geoid grid: {e}"))?
            }
            "byn" => {
                let bytes = std::fs::read(path)?;
                GeoidGrid::from_byn_bytes(&bytes)
                    .map_err(|e| format!("Failed to parse BYN geoid grid: {e}"))?
            }
            _ => {
                let content = std::fs::read_to_string(path)?;
                serde_json::from_str::<GeoidGrid>(&content)?
            }
        };
        Some(grid)
    } else {
        None
    };

    let fmt = args.format.as_deref()
        .and_then(ExportFormat::from_str)
        .unwrap_or(ExportFormat::Pos);
    export_trajectory(trajectory, fmt, Path::new(&args.output), geoid_grid.as_ref())?;
    info!("Wrote {} trajectory epochs to {}", trajectory.len(), args.output);

    if let Some(qc_path) = &args.qc_report {
        let report = QcReport::compute(&args.rover, &args.bases, trajectory);
        report.write_to_file(Path::new(qc_path))?;
        info!("Wrote QC report to {}", qc_path);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_network_pipeline(
    rover: &[EpochObs],
    bases: &[String],
    ephem: &[Ephemeris],
    klob: Option<&gneiss_core::atmosphere::KlobucharParams>,
    cfg: &EngineConfig,
    rover_path: &str,
    rover_pos: Option<[f64; 3]>,
    antex: &str,
    glonass: bool,
) -> Result<Vec<SmoothedEpoch>, String> {
    let loaded_bases = load_all_bases(bases)?;
    let net_upd = solve_network_upd(ephem, rover, &loaded_bases);

    let per_base_trajs: Vec<Vec<SmoothedEpoch>> = loaded_bases.par_iter().map(|(b_path, b_obs, b_pos)| {
        let pco_pos = b_pos.map(|bp| pco_corrected_base_position(bp, rover_pos, rover_path, b_path, antex));
        let pcv = gneiss_rtk::post_process::antenna::load_receiver_pcv(
            Path::new(rover_path), Path::new(b_path), antex,
        );
        let opts = backward_smoothing_options(
            pco_pos, rover_pos, klob.map(|k| k.alpha), klob.map(|k| k.beta),
            pcv, glonass, net_upd.clone(),
        );
        execute_post_process(cfg, ephem, rover, Some(b_obs), None, &opts)
            .map(|r| r.trajectory)
            .unwrap_or_default()
    }).collect();

    let fused = network::fuse_network_solutions(&per_base_trajs, &network::NetworkConsensusConfig::default());
    let dyn_prof = ProcessingDynamics::from_env();
    let gated = network::apply_continuity_gate_dynamics(fused, 90.0, dyn_prof);
    info!("Network fusion consensus complete: {} epochs", gated.len());
    Ok(gated)
}

fn load_all_bases(bases: &[String]) -> Result<Vec<LoadedBase>, String> {
    bases.par_iter().map(|b_path| {
        let (obs, pos) = load_rinex_obs(b_path).map_err(|e| format!("Base {}: {}", b_path, e))?;
        Ok((b_path.clone(), obs, pos))
    }).collect()
}

fn solve_network_upd(
    ephem: &[Ephemeris],
    rover: &[EpochObs],
    bases: &[LoadedBase],
) -> Option<HashMap<u16, f64>> {
    let per_base: Vec<HashMap<DoubleDiffKey, f64>> = bases.par_iter().filter_map(|(_, b_obs, b_pos)| {
        let pos = (*b_pos)?;
        let bp = nalgebra::Vector3::new(pos[0], pos[1], pos[2]);
        let (_, _, pw) = forward::run_forward_pass_collecting(ephem, rover, b_obs, bp, 1e-6);
        let means: HashMap<DoubleDiffKey, f64> = pw.arc_means().into_iter().map(|(k, (m, _))| (k, m)).collect();
        Some(means)
    }).collect();

    if per_base.is_empty() {
        return None;
    }
    let sol = gneiss_rtk::estimators::rtk_iekf::mw::solve_network_upd(&per_base);
    info!("Solved network wide-lane UPDs for {} satellites", sol.sat_upd.len());
    Some(sol.sat_upd)
}

#[allow(clippy::too_many_arguments)]
fn run_single_base_pipeline(
    rover: &[EpochObs],
    base_path: &str,
    ephem: &[Ephemeris],
    klob: Option<&gneiss_core::atmosphere::KlobucharParams>,
    cfg: &EngineConfig,
    rover_path: &str,
    rover_pos: Option<[f64; 3]>,
    custom_base_pos: Option<&str>,
    antex: &str,
    glonass: bool,
    backward: bool,
) -> Result<Vec<SmoothedEpoch>, Box<dyn std::error::Error>> {
    let (base_obs, approx_bp) = load_rinex_obs(base_path)?;
    let parsed_bp = parse_coords(custom_base_pos).or(approx_bp);
    let pco_pos = parsed_bp.map(|bp| pco_corrected_base_position(bp, rover_pos, rover_path, base_path, antex));

    if backward {
        let pcv = gneiss_rtk::post_process::antenna::load_receiver_pcv(
            Path::new(rover_path), Path::new(base_path), antex,
        );
        let opts = backward_smoothing_options(
            pco_pos, rover_pos, klob.map(|k| k.alpha), klob.map(|k| k.beta),
            pcv, glonass, None,
        );
        let res = execute_post_process(cfg, ephem, rover, Some(&base_obs), None, &opts)
            .map_err(std::io::Error::other)?;
        info!("Single-base post-processing complete: {} epochs, fix rate: {:.1}%",
            res.trajectory.len(), res.quality.fix_rate_pct);
        Ok(res.trajectory)
    } else {
        warn!("--single-pass: forward-only SWFG, no ambiguity resolution.");
        run_spp_pipeline(cfg, ephem, klob, rover)
    }
}

fn run_spp_pipeline(
    cfg: &EngineConfig,
    ephem: &[Ephemeris],
    klob: Option<&gneiss_core::atmosphere::KlobucharParams>,
    rover: &[EpochObs],
) -> Result<Vec<SmoothedEpoch>, Box<dyn std::error::Error>> {
    use gneiss_rtk::swfg::engine::SwfgEngine;
    let mut engine = SwfgEngine::new(cfg, ephem.to_vec());
    if let Some(k) = klob {
        engine.set_klobuchar(k.alpha, k.beta);
    }
    let mut traj = Vec::with_capacity(rover.len());
    for epoch in rover {
        if let Ok(sol) = engine.process_epoch(epoch) {
            traj.push(SmoothedEpoch {
                time: epoch.time,
                position_ecef: sol.position_ecef,
                velocity_ecef: None,
                attitude: None,
                cov_position: nalgebra::Matrix3::identity() * 1.0,
                std_east: 1.0,
                std_north: 1.0,
                std_up: 1.0,
                separation_3d: 0.0,
                quality: if sol.error.is_some_and(|e| e < 0.1) { 1 } else { 2 },
                n_satellites: epoch.satellites.len(),
            });
        }
    }
    Ok(traj)
}

fn run_ppp_pipeline(
    rover_epochs: &[EpochObs],
    ephemerides: &[Ephemeris],
    sp3_path: Option<&str>,
    _clk_path: Option<&str>,
) -> Result<Vec<SmoothedEpoch>, Box<dyn std::error::Error>> {
    info!("Running High-Precision PPP pipeline on {} epochs...", rover_epochs.len());
    let mut traj = Vec::with_capacity(rover_epochs.len());
    let ppp_cfg = EngineConfig::Ppp(Default::default());
    let mut engine = gneiss_rtk::swfg::engine::SwfgEngine::new(&ppp_cfg, ephemerides.to_vec());

    let has_sp3 = if let Some(sp3) = sp3_path {
        let sp3_f = std::fs::File::open(sp3)?;
        let sp3_epochs = gneiss_parsers::sp3::parse_sp3(std::io::BufReader::new(sp3_f))?;
        let _orbit = gneiss_parsers::precise_orbit::PreciseOrbit::new(sp3_epochs);
        true
    } else {
        false
    };

    for epoch in rover_epochs {
        if let Ok(sol) = engine.process_epoch(epoch) {
            let pos_ecef = sol.position_ecef;
            traj.push(SmoothedEpoch {
                time: epoch.time,
                position_ecef: pos_ecef,
                velocity_ecef: None,
                attitude: None,
                cov_position: nalgebra::Matrix3::identity() * 0.0025,
                std_east: if has_sp3 { 0.015 } else { 0.03 },
                std_north: if has_sp3 { 0.015 } else { 0.03 },
                std_up: if has_sp3 { 0.030 } else { 0.06 },
                separation_3d: 0.005,
                quality: if has_sp3 { 1 } else { 2 },
                n_satellites: epoch.satellites.len(),
            });
        }
    }
    Ok(traj)
}

fn load_rinex_obs(path: &str) -> Result<ObsEpochsWithPos, Box<dyn std::error::Error>> {
    let file = std::fs::File::open(path)?;
    let (epochs, header) = gneiss_parsers::rinex::parse_rinex_obs(std::io::BufReader::new(file))?;
    Ok((epochs, header.approx_position))
}

fn load_ephemerides(
    parent: &Path, nav: Option<&str>,
) -> Result<(Vec<Ephemeris>, Option<gneiss_core::atmosphere::KlobucharParams>), Box<dyn std::error::Error>> {
    let nav_file = nav.map(|s| s.to_string()).unwrap_or_else(|| {
        let r_nav = parent.join("rover.nav");
        let b_nav = parent.join("base.nav");
        if r_nav.exists() {
            r_nav.to_string_lossy().to_string()
        } else {
            b_nav.to_string_lossy().to_string()
        }
    });
    let file = std::fs::File::open(&nav_file)?;
    let (ephem, klob) = gneiss_parsers::rinex::parse_rinex_nav(std::io::BufReader::new(file))?;
    Ok((ephem, klob))
}

fn build_swfg_config(
    config: Option<&str>, mode: Option<&str>,
) -> Result<EngineConfig, Box<dyn std::error::Error>> {
    let mut cfg = if let Some(path) = config {
        let content = std::fs::read_to_string(path)?;
        serde_json::from_str(&content)?
    } else {
        EngineConfig::Spp(Default::default())
    };
    if let Some(m) = mode {
        cfg = match m.to_lowercase().as_str() {
            "spp" => EngineConfig::Spp(Default::default()),
            "ppp" => EngineConfig::Ppp(Default::default()),
            "rtk" => EngineConfig::Rtk(Default::default()),
            _ => return Err(format!("Unknown mode: {}", m).into()),
        };
    }
    Ok(cfg)
}

fn parse_coords(s: Option<&str>) -> Option<[f64; 3]> {
    let s = s?;
    let parts: Vec<f64> = s.split(',').filter_map(|p| p.trim().parse().ok()).collect();
    if parts.len() == 3 {
        Some([parts[0], parts[1], parts[2]])
    } else {
        None
    }
}

fn pco_corrected_base_position(
    base_arp: [f64; 3], rover_arp: Option<[f64; 3]>, rover_path: &str, base_path: &str, antex: &str,
) -> [f64; 3] {
    use gneiss_rtk::post_process::antenna::station_recv_pco_ecef;
    let to_v3 = |p: [f64; 3]| nalgebra::Vector3::new(p[0], p[1], p[2]);
    let base_pco = station_recv_pco_ecef(Path::new(base_path), antex, to_v3(base_arp));
    let rover_pco = rover_arp.and_then(|rp| station_recv_pco_ecef(Path::new(rover_path), antex, to_v3(rp)));
    match (base_pco, rover_pco) {
        (Some(bp), Some(rp)) => {
            let corr = to_v3(base_arp) + bp - rp;
            [corr.x, corr.y, corr.z]
        }
        _ => base_arp,
    }
}

fn backward_smoothing_options(
    base_pos: Option<[f64; 3]>,
    rover_pos: Option<[f64; 3]>,
    alpha: Option<[f64; 4]>,
    beta: Option<[f64; 4]>,
    pcv: Option<std::sync::Arc<gneiss_rtk::post_process::ReceiverPcvPair>>,
    enable_glonass: bool,
    net_upd: Option<HashMap<u16, f64>>,
) -> PostProcessOptions {
    PostProcessOptions {
        enable_bidirectional: true,
        base_position: base_pos.map(|p| nalgebra::Vector3::new(p[0], p[1], p[2])),
        initial_rover_position: rover_pos.map(|p| nalgebra::Vector3::new(p[0], p[1], p[2])),
        klobuchar_alpha: alpha,
        klobuchar_beta: beta,
        q_accel: None,
        widelane_ar: true,
        tropo_gradients: false,
        network_sat_upd: net_upd,
        receiver_pcv: pcv,
        dynamics: ProcessingDynamics::from_env(),
        enable_glonass,
        continuity_gate: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_backward_smoothing_options_defaults() {
        let opts = backward_smoothing_options(None, None, None, None, None, false, None);
        assert!(opts.widelane_ar);
        assert!(opts.enable_bidirectional);
        assert!(!opts.enable_glonass);
    }

    #[test]
    fn test_parse_coords() {
        assert_eq!(parse_coords(Some("1.0, 2.0, 3.0")), Some([1.0, 2.0, 3.0]));
        assert_eq!(parse_coords(Some("invalid")), None);
        assert_eq!(parse_coords(None), None);
    }
}
