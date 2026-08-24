use tracing::{info, error};
use gneiss_rtk::swfg::config::EngineConfig;
use gneiss_rtk::swfg::engine::SwfgEngine;
type ObsEpochsWithPos = (Option<Vec<gneiss_core::obs::EpochObs>>, Option<[f64; 3]>);

#[allow(clippy::too_many_arguments)]
pub async fn run_process(
    rover: String, base: Option<String>, nav: Option<String>, output: String, config: Option<String>,
    _enable_backward_smoothing: bool, mode: Option<String>,
    _lambda_ratio: Option<f64>, _lambda_subset: Option<usize>, max_epochs: Option<usize>,
    _lever_arm: String, _calibrate_imu: bool,
    _raim_outlier_m: Option<f64>, _chi_square_pr: Option<f64>, _chi_square_cp: Option<f64>, _nominal_snr: Option<f64>,
    _base_position: Option<String>, _systems: Option<String>, _sp3: Option<String>, _clk: Option<String>,
    _antex: Option<String>, _clock_jump_threshold: Option<f64>, _disable_doppler: bool, _bia: Option<String>
) -> Result<(), Box<dyn std::error::Error>> {
    info!("Starting SWFG Processing Pipeline...");

    let (rover_rinex_epochs, _rover_approx_pos) = load_rover_epochs(&rover)?;
    let (base_rinex_epochs, _approx_base_pos) = load_base_epochs(&base)?;

    let swfg_config = build_swfg_config(config, mode)?;

    let parent_dir = std::path::Path::new(&rover).parent().unwrap_or_else(|| std::path::Path::new("."));
    let (ephemerides, klobuchar) = load_ephemerides(parent_dir, nav)?;

    let mut engine = SwfgEngine::new(&swfg_config, ephemerides.clone());
    if let Some(ref k) = klobuchar {
        engine.set_klobuchar(k.alpha, k.beta);
    }

    let rover_epochs = rover_rinex_epochs.ok_or("no rover data")?;
    let max = max_epochs.unwrap_or(usize::MAX);
    let selected_rover = &rover_epochs[..rover_epochs.len().min(max)];

    let parsed_base_pos: Option<[f64; 3]> = _base_position.and_then(|bp_str| {
        let parts: Vec<f64> = bp_str.split(',').filter_map(|s| s.trim().parse().ok()).collect();
        if parts.len() == 3 {
            Some([parts[0], parts[1], parts[2]])
        } else {
            None
        }
    }).or(_approx_base_pos);

    let mut trajectory = Vec::new();
    if _enable_backward_smoothing {
        info!("Running Qinertia-grade 4-pass offline post-processing pipeline...");
        let post_options = gneiss_rtk::post_process::PostProcessOptions {
            enable_bidirectional: true,
            base_position: parsed_base_pos.map(|p| nalgebra::Vector3::new(p[0], p[1], p[2])),
            initial_rover_position: _rover_approx_pos.map(|p| nalgebra::Vector3::new(p[0], p[1], p[2])),
            klobuchar_alpha: klobuchar.as_ref().map(|k| k.alpha),
            klobuchar_beta: klobuchar.as_ref().map(|k| k.beta),
            q_accel: None,
            widelane_ar: false,
            tropo_gradients: false,
            network_sat_upd: None,
        };
        let post_res = gneiss_rtk::post_process::execute_post_process(
            &swfg_config,
            &ephemerides,
            selected_rover,
            base_rinex_epochs.as_deref(),
            None,
            &post_options,
        ).map_err(std::io::Error::other)?;

        info!("Post-processing complete: {} epochs, fix rate: {:.1}%, median sep: {:.3}m",
            post_res.trajectory.len(), post_res.quality.fix_rate_pct, post_res.quality.median_separation_m);

        for ep in &post_res.trajectory {
            trajectory.push((ep.time.week as u16, ep.time.tow, ep.position_ecef, ep.quality));
        }
    } else {
        for epoch in selected_rover {
            let base_opt = base_rinex_epochs.as_ref().and_then(|b_epochs| {
                b_epochs.iter().min_by(|a, b|
                    (a.time.tow - epoch.time.tow).abs()
                        .partial_cmp(&(b.time.tow - epoch.time.tow).abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                )
            });
            let res = if let Some(base_ep) = base_opt {
                let bp = parsed_base_pos.unwrap_or([0.0; 3]);
                engine.process_rtk_epoch(epoch, base_ep, nalgebra::Vector3::new(bp[0], bp[1], bp[2]))
            } else {
                engine.process_epoch(epoch)
            };
            match res {
                Ok(sol) => {
                    let q: u8 = if sol.error.is_some_and(|e| e < 0.1) { 1 } else { 2 };
                    trajectory.push((epoch.time.week as u16, epoch.time.tow, sol.position_ecef, q));
                }
                Err(e) => error!("Epoch {} failed: {}", epoch.time.tow, e),
            }
        }
    }

    let n_processed = trajectory.len();
    write_results(&trajectory, &output, n_processed).await?;
    Ok(())
}

fn load_rover_epochs(rover: &str) -> Result<ObsEpochsWithPos, Box<dyn std::error::Error>> {
    let file = std::fs::File::open(rover)?;
    let (epochs, header) = gneiss_parsers::rinex::parse_rinex_obs(std::io::BufReader::new(file))?;
    info!("Loaded {} RINEX rover epochs.", epochs.len());
    Ok((Some(epochs), header.approx_position))
}

fn load_base_epochs(base: &Option<String>) -> Result<ObsEpochsWithPos, Box<dyn std::error::Error>> {
    if let Some(base_file) = base {
        let file = std::fs::File::open(base_file)?;
        let (epochs, header) = gneiss_parsers::rinex::parse_rinex_obs(std::io::BufReader::new(file))?;
        info!("Loaded {} RINEX base epochs.", epochs.len());
        Ok((Some(epochs), header.approx_position))
    } else {
        Ok((None, None))
    }
}

fn build_swfg_config(
    config: Option<String>, mode: Option<String>,
) -> Result<EngineConfig, Box<dyn std::error::Error>> {
    let mut engine_config = if let Some(config_path) = config {
        let content = std::fs::read_to_string(&config_path)?;
        serde_json::from_str(&content)?
    } else {
        EngineConfig::Spp(Default::default())
    };

    if let Some(m) = mode {
        engine_config = match m.to_lowercase().as_str() {
            "spp" => EngineConfig::Spp(Default::default()),
            "ppp" => EngineConfig::Ppp(Default::default()),
            "rtk" => EngineConfig::Rtk(Default::default()),
            _ => return Err(format!("Unknown mode: {}", m).into()),
        };
    }

    Ok(engine_config)
}

fn load_ephemerides(
    parent_dir: &std::path::Path, nav: Option<String>,
) -> Result<(Vec<gneiss_core::ephemeris::Ephemeris>, Option<gneiss_core::atmosphere::KlobucharParams>), Box<dyn std::error::Error>> {
    let nav_file = nav.unwrap_or_else(|| {
        let r_nav = parent_dir.join("rover.nav");
        let b_nav = parent_dir.join("base.nav");
        if r_nav.exists() {
            r_nav.to_string_lossy().to_string()
        } else {
            b_nav.to_string_lossy().to_string()
        }
    });
    let file = std::fs::File::open(&nav_file)?;
    let (ephemerides, klobuchar) = gneiss_parsers::rinex::parse_rinex_nav(std::io::BufReader::new(file))?;
    info!("Loaded {} ephemerides", ephemerides.len());
    Ok((ephemerides, klobuchar))
}

async fn write_results(
    trajectory: &[(u16, f64, nalgebra::Vector3<f64>, u8)],
    output: &str,
    n_processed: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    use tokio::io::AsyncWriteExt;
    let mut file = tokio::fs::File::create(output).await?;
    file.write_all(b"% GPST-Week TOW(s) x-ecef(m) y-ecef(m) z-ecef(m) Q\n").await?;
    for (week, tow, pos, q) in trajectory {
        let line = format!("{} {:.3} {:.4} {:.4} {:.4} {}\n", week, tow, pos.x, pos.y, pos.z, q);
        file.write_all(line.as_bytes()).await?;
    }
    info!("Wrote SWFG result ({} epochs processed) to {}", n_processed, output);
    Ok(())
}
