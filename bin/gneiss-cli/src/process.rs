use tracing::{info, error};
use gneiss_rtk::swfg::config::EngineConfig;
use gneiss_rtk::swfg::engine::SwfgEngine;
type ObsEpochsWithPos = (Option<Vec<gneiss_core::obs::EpochObs>>, Option<[f64; 3]>);

#[allow(clippy::too_many_arguments)]
pub async fn run_process(
    rover: String, base: Option<String>, nav: Option<String>, output: String, config: Option<String>,
    _enable_backward_smoothing: bool, mode: Option<String>, max_epochs: Option<usize>,
    _base_position: Option<String>, systems: Option<String>, _sp3: Option<String>, _clk: Option<String>,
    antex: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    if _sp3.is_some() || _clk.is_some() {
        return Err("--sp3/--clk are not yet wired into gneiss-cli: precise ephemeris \
            requires satellite PCO and precise clock together or it measurably degrades \
            accuracy (see docs/NETWORK_RTK_NEXT_STEPS.md, \"Precise ephemeris\"); refusing \
            rather than silently running broadcast-only".into());
    }

    let (mut rover_rinex_epochs, _rover_approx_pos) = load_rover_epochs(&rover)?;
    let (mut base_rinex_epochs, _approx_base_pos) = load_base_epochs(&base)?;
    if let Some(systems) = systems.as_deref() {
        if let Some(epochs) = rover_rinex_epochs.as_mut() {
            gneiss_core::obs::filter_constellations(epochs, systems);
        }
        if let Some(epochs) = base_rinex_epochs.as_mut() {
            gneiss_core::obs::filter_constellations(epochs, systems);
        }
    }

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
        let antex_path = antex.unwrap_or_else(|| "datasets/igs14.atx".to_string());
        let receiver_pcv = base.as_deref().and_then(|base_path| {
            gneiss_rtk::post_process::antenna::load_receiver_pcv(
                std::path::Path::new(&rover),
                std::path::Path::new(base_path),
                &antex_path,
            )
        });
        let pco_base_pos = match (base.as_deref(), parsed_base_pos) {
            (Some(base_path), Some(bp)) => Some(pco_corrected_base_position(
                bp, _rover_approx_pos, &rover, base_path, &antex_path,
            )),
            _ => parsed_base_pos,
        };
        let post_options = backward_smoothing_options(
            pco_base_pos,
            _rover_approx_pos,
            klobuchar.as_ref().map(|k| k.alpha),
            klobuchar.as_ref().map(|k| k.beta),
            receiver_pcv,
        );
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
        tracing::warn!("--single-pass: forward-only SWFG, no ambiguity resolution. \
            Typically 5-10x less accurate than the default 4-pass pipeline. \
            Drop --single-pass unless you specifically need it.");
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
    info!("Wrote result ({} epochs processed) to {}", n_processed, output);
    Ok(())
}

/// Differential receiver PCO: shift the base ARP by `(base_pco - rover_pco)`
/// so a cross-family baseline (e.g. a Leica base against a Trimble rover)
/// references real antenna phase centres instead of ARPs. Same-family
/// pairs move by ~0mm (the two corrections cancel; this is physics, not a
/// special case). Falls back to the raw ARP when either antenna's
/// calibration can't be resolved (unlisted type, missing ANTEX file) —
/// this measurably fixed a -54mm vertical bias on a real cross-family
/// baseline in docs/NETWORK_RTK_NEXT_STEPS.md ("Receiver PCO").
fn pco_corrected_base_position(
    base_arp: [f64; 3],
    rover_arp: Option<[f64; 3]>,
    rover_path: &str,
    base_path: &str,
    antex_path: &str,
) -> [f64; 3] {
    use gneiss_rtk::post_process::antenna::station_recv_pco_ecef;
    let to_v3 = |p: [f64; 3]| nalgebra::Vector3::new(p[0], p[1], p[2]);

    let base_pco = station_recv_pco_ecef(std::path::Path::new(base_path), antex_path, to_v3(base_arp));
    let rover_pco = rover_arp.and_then(|rp| {
        station_recv_pco_ecef(std::path::Path::new(rover_path), antex_path, to_v3(rp))
    });
    match (base_pco, rover_pco) {
        (Some(bp), Some(rp)) => {
            let corrected = to_v3(base_arp) + bp - rp;
            [corrected.x, corrected.y, corrected.z]
        }
        _ => base_arp,
    }
}

/// Options for the `--enable-backward-smoothing` (4-pass DD-RTK) path.
///
/// `widelane_ar` is fixed on: cadence-aware slip gating, reverse-safe
/// process noise and innovation-gated slip re-seeding are measured
/// strictly positive with no regressions across every CORS baseline on
/// record (see docs/NETWORK_RTK_NEXT_STEPS.md). `eval_network_ppk`
/// already defaults this on; the CLI must match its validated behavior
/// rather than silently running the legacy path.
fn backward_smoothing_options(
    base_position: Option<[f64; 3]>,
    rover_approx_pos: Option<[f64; 3]>,
    klobuchar_alpha: Option<[f64; 4]>,
    klobuchar_beta: Option<[f64; 4]>,
    receiver_pcv: Option<std::sync::Arc<gneiss_rtk::post_process::ReceiverPcvPair>>,
) -> gneiss_rtk::post_process::PostProcessOptions {
    gneiss_rtk::post_process::PostProcessOptions {
        enable_bidirectional: true,
        base_position: base_position.map(|p| nalgebra::Vector3::new(p[0], p[1], p[2])),
        initial_rover_position: rover_approx_pos.map(|p| nalgebra::Vector3::new(p[0], p[1], p[2])),
        klobuchar_alpha,
        klobuchar_beta,
        q_accel: None,
        widelane_ar: true,
        tropo_gradients: false,
        network_sat_upd: None,
        receiver_pcv,
        dynamics: gneiss_rtk::post_process::ProcessingDynamics::from_env(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn backward_smoothing_enables_widelane_ar_by_default() {
        let opts = backward_smoothing_options(None, None, None, None, None);
        assert!(
            opts.widelane_ar,
            "gneiss-cli must use the validated widelane_ar default, not the legacy off-path"
        );
        assert!(opts.enable_bidirectional);
        assert!(opts.receiver_pcv.is_none());
    }

    /// Real CORS data end-to-end: a cross-family pair (P224 rover, Trimble;
    /// CAPO base, Leica) must resolve through the CLI's own antex wiring,
    /// gated on the dataset being checked out.
    #[test]
    fn receiver_pcv_resolves_for_real_cross_family_pair() {
        let antex = std::path::PathBuf::from("../../datasets/igs14.atx");
        let rover = std::path::PathBuf::from("../../datasets/cors_short_baseline/p2241350.20o");
        let base = std::path::PathBuf::from("../../datasets/cors_short_baseline/capo1350.20o");
        if !antex.exists() || !rover.exists() || !base.exists() {
            return;
        }
        let pair = gneiss_rtk::post_process::antenna::load_receiver_pcv(
            &rover, &base, antex.to_str().unwrap(),
        ).expect("CAPO (LEIAR20) vs the Trimble rover must resolve against igs14.atx");
        assert_eq!(pair.base.ant_type, "LEIAR20");

        let opts = backward_smoothing_options(None, None, None, None, Some(pair));
        assert!(opts.receiver_pcv.is_some());
    }

    /// The actual bias fix (docs/NETWORK_RTK_NEXT_STEPS.md: CAPO v_p50
    /// -54 -> -14mm) lives in this shift, not in `receiver_pcv` above.
    #[test]
    fn pco_corrected_base_position_moves_for_real_cross_family_pair() {
        let antex = std::path::PathBuf::from("../../datasets/igs14.atx");
        let rover_path = "../../datasets/cors_short_baseline/p2241350.20o";
        let base_path = "../../datasets/cors_short_baseline/capo1350.20o";
        if !antex.exists()
            || !std::path::Path::new(rover_path).exists()
            || !std::path::Path::new(base_path).exists()
        {
            return;
        }
        let base_arp = [-2693675.7831, -4273829.9413, 3880383.2888];
        let rover_arp = Some([-2688181.50, -4265663.45, 3893784.80]);
        let corrected = pco_corrected_base_position(
            base_arp, rover_arp, rover_path, base_path, antex.to_str().unwrap(),
        );
        let shift = ((corrected[0] - base_arp[0]).powi(2)
            + (corrected[1] - base_arp[1]).powi(2)
            + (corrected[2] - base_arp[2]).powi(2))
        .sqrt();
        assert!(shift > 1e-4, "cross-family (Leica/Trimble) shift should be non-trivial, got {shift} m");
        assert!(shift < 1.0, "shift should be centimetre-scale, not metres, got {shift} m");
    }

    #[test]
    fn pco_corrected_base_position_falls_back_when_antex_missing() {
        let base_arp = [-2693675.7831, -4273829.9413, 3880383.2888];
        let corrected = pco_corrected_base_position(
            base_arp, Some([0.0, 0.0, 0.0]), "rover.obs", "base.obs", "/nonexistent.atx",
        );
        assert_eq!(corrected, base_arp);
    }
}
