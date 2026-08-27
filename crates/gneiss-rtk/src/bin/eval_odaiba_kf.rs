//! EKF + RTS smoother evaluation on Odaiba dataset.
//!
//! Matches Qinertia's GNSS-only post-processing architecture:
//!   1. Forward EKF pass with constant-velocity process model
//!   2. RTS backward smoother
//!   3. Output: smoothed positions with formal covariances
//!
//! Usage: cargo run --release --bin eval_odaiba_kf

use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use nalgebra::Vector3;

use gneiss_rtk::swfg::kalman_smoother::{GnssKalmanSmoother, ProcessNoise};
use gneiss_rtk::swfg::pipeline::{MeasurementPipeline, ReceiverState};

fn main() {
    tracing_subscriber::fmt()
        .with_env_filter("warn")
        .with_target(false)
        .without_time()
        .try_init()
        .ok();

    let dataset = Path::new("datasets/urbannav/tokyo/Tokyo_Data/Odaiba");
    println!("=== Odaiba EKF + RTS Smoother Benchmark ===");

    // Load navigation data
    let nav_file = File::open(dataset.join("base.nav")).expect("open base.nav");
    let (ephemerides, klobuchar) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_file))
        .expect("parse nav");

    // Load reference truth
    let ref_csv = std::fs::read_to_string(dataset.join("reference.csv")).expect("read reference.csv");
    let mut truth: BTreeMap<u32, (f64, f64, f64)> = BTreeMap::new();
    for line in ref_csv.lines().skip(1) {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 8 { continue; }
        if let (Ok(tow), Ok(x), Ok(y), Ok(z)) = (
            parts[0].trim().parse::<f64>(),
            parts[5].trim().parse::<f64>(),
            parts[6].trim().parse::<f64>(),
            parts[7].trim().parse::<f64>(),
        ) {
            truth.insert(tow.floor() as u32, (x, y, z));
        }
    }

    // Load RINEX observations
    let rover_file = File::open(dataset.join("rover_trimble.obs")).expect("open rover_trimble.obs");
    let (rover_epochs, _) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(rover_file))
        .expect("parse rover RINEX");

    let base_file = File::open(dataset.join("base_trimble.obs")).expect("open base_trimble.obs");
    let (base_epochs, base_approx) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(base_file))
        .expect("parse base RINEX");

    // Base station position (with antenna delta)
    let default_base_arp = Vector3::new(-3961904.4341, 3348994.266, 3698211.7067);
    let default_base_llh = gneiss_core::coords::ecef_to_llh(default_base_arp);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(default_base_llh).transpose();
    let computed_base = default_base_arp + ned_to_ecef * Vector3::new(0.0, 0.0, 0.0855);
    let base_pos = base_approx.approx_position
        .map(|p| Vector3::new(p[0], p[1], p[2]))
        .unwrap_or(computed_base);

    // Measurement pipeline (SPP-mode: broadcast clock + tropo + iono + variance)
    let mut pipeline = MeasurementPipeline::spp_mode();
    if let Some(ref k) = klobuchar {
        pipeline.set_klobuchar(k.alpha, k.beta);
    }

    // Seed position from first epoch's SPP solution
    let seed_pos = gneiss_rtk::swfg::engine::epoch::compute_spp_seeding(
        &rover_epochs[0], &ephemerides,
    ).unwrap_or(base_pos);

    println!("Seed position: ({:.1}, {:.1}, {:.1})", seed_pos.x, seed_pos.y, seed_pos.z);
    println!("Processing {} rover epochs...", rover_epochs.len());

    let mut kf = GnssKalmanSmoother::new(seed_pos, ProcessNoise::default());

    // ── Forward EKF pass ──
    let start = std::time::Instant::now();
    for (i, epoch) in rover_epochs.iter().enumerate() {
        let time = epoch.time;

        // Predict
        kf.predict(time);

        // Extract raw observations and apply corrections
        let rx_pos = Vector3::new(kf.x[0], kf.x[1], kf.x[2]);
        let raw_obs = match gneiss_rtk::swfg::engine::epoch::extract_raw_observations(
            epoch, &ephemerides, Some(rx_pos),
        ) {
            Ok(obs) => obs,
            Err(_) => {
                kf.update(time, &[]);
                continue;
            }
        };

        let rx_llh = gneiss_core::coords::ecef_to_llh(rx_pos);
        let rx_state = ReceiverState {
            position_ecef: rx_pos,
            llh_rad: rx_llh,
            clock_bias_m: vec![kf.x[6]],
            zwd_m: 0.1,
            ifb_glo: 0.0,
            time,
        };

        let corrected = pipeline.process(&raw_obs, &rx_state);

        // Form DD observations against base
        let base_epoch = base_epochs.iter()
            .filter(|b| (b.time.tow - epoch.time.tow).abs() < 0.1)
            .min_by(|a, b| {
                (a.time.tow - epoch.time.tow).abs()
                    .partial_cmp(&(b.time.tow - epoch.time.tow).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

        // Build measurement vector: (sat_pos, corrected_pr, variance)
        let observations: Vec<(Vector3<f64>, f64, f64)> = if let Some(base) = base_epoch {
            let base_raw = gneiss_rtk::swfg::engine::epoch::extract_raw_observations(
                base, &ephemerides, Some(base_pos),
            ).unwrap_or_default();

            // SD pseudoranges (rover - base + base_range)
            corrected.iter().filter_map(|rov| {
                base_raw.iter().find(|b| b.satellite == rov.satellite).map(|bas| {
                    let base_range = (bas.sat_pos_ecef - base_pos).norm();
                    let sd_pr = (rov.pr_l1 - bas.pr_l1) + base_range;
                    (rov.sat_pos_ecef, sd_pr, rov.variance_m2 * 2.0) // SD doubles variance
                })
            }).collect()
        } else {
            // Standalone: use corrected pseudoranges directly
            corrected.iter().map(|obs| {
                let pr_corrected = obs.pr_l1 - obs.sat_clock_m + obs.tropo_dry_m + obs.iono_l1_m;
                (obs.sat_pos_ecef, pr_corrected, obs.variance_m2)
            }).collect()
        };

        kf.update(time, &observations);

        if (i + 1) % 1000 == 0 || i == rover_epochs.len() - 1 {
            let tow = epoch.time.tow.floor() as u32;
            println!("Forward: epoch {}/{} (TOW {}) — {} sats",
                i + 1, rover_epochs.len(), tow, observations.len());
        }
    }

    let fwd_elapsed = start.elapsed();
    println!("Forward pass: {:.1}s ({} epochs stored)", fwd_elapsed.as_secs_f64(), kf.history.len());

    // ── RTS backward smoother ──
    let rts_start = std::time::Instant::now();
    let smoothed = kf.smooth();
    let rts_elapsed = rts_start.elapsed();
    println!("RTS smoother: {:.1}s ({} smoothed states)", rts_elapsed.as_secs_f64(), smoothed.len());

    // ── Evaluate accuracy ──
    let mut h_errors_fwd = Vec::new();
    let mut h_errors_smooth = Vec::new();
    let mut errors_3d_smooth = Vec::new();
    let mut matched = 0usize;

    for (i, state) in smoothed.iter().enumerate() {
        if state.n_obs < 4 { continue; }

        let tow_sec = state.time.tow.floor() as u32;
        if let Some(&(tx, ty, tz)) = truth.get(&tow_sec) {
            matched += 1;
            let truth_pos = Vector3::new(tx, ty, tz);
            let truth_llh = gneiss_core::coords::ecef_to_llh(truth_pos);

            // Smoothed error
            let enu_s = gneiss_core::coords::ecef_delta_to_enu(state.position_ecef, truth_pos, truth_llh);
            let h_s = (enu_s.x * enu_s.x + enu_s.y * enu_s.y).sqrt();
            h_errors_smooth.push(h_s);
            errors_3d_smooth.push(enu_s.norm());

            // Forward-only error for comparison
            let fwd_pos = Vector3::new(
                kf.history[i].x_updated[0],
                kf.history[i].x_updated[1],
                kf.history[i].x_updated[2],
            );
            let enu_f = gneiss_core::coords::ecef_delta_to_enu(fwd_pos, truth_pos, truth_llh);
            let h_f = (enu_f.x * enu_f.x + enu_f.y * enu_f.y).sqrt();
            h_errors_fwd.push(h_f);

            if matched.is_multiple_of(100) {
                let pos_sigma = (state.p_smoothed[(0,0)] + state.p_smoothed[(1,1)]).sqrt();
                println!("TOW={} fwd={:.2}m smooth={:.2}m σ_hz={:.2}m",
                    tow_sec, h_f, h_s, pos_sigma);
            }
        }
    }

    // Print statistics
    if h_errors_smooth.is_empty() {
        eprintln!("No matched epochs!");
        return;
    }

    h_errors_fwd.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    h_errors_smooth.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    errors_3d_smooth.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let n = h_errors_smooth.len();
    let percentile = |v: &[f64], p: f64| v[((n as f64 * p) as usize).min(n - 1)];
    let rms = |v: &[f64]| (v.iter().map(|e| e * e).sum::<f64>() / v.len() as f64).sqrt();

    println!("\n=== Odaiba EKF + RTS Smoother Results ===");
    println!("Total time: {:.1}s (forward {:.1}s + RTS {:.1}s)",
        fwd_elapsed.as_secs_f64() + rts_elapsed.as_secs_f64(),
        fwd_elapsed.as_secs_f64(), rts_elapsed.as_secs_f64());
    println!("Matched epochs: {} (of {} with ≥4 sats)", n, smoothed.iter().filter(|s| s.n_obs >= 4).count());

    println!("\nForward EKF (horizontal):");
    println!("  p50:  {:.3}m", percentile(&h_errors_fwd, 0.50));
    println!("  p68:  {:.3}m", percentile(&h_errors_fwd, 0.68));
    println!("  p95:  {:.3}m", percentile(&h_errors_fwd, 0.95));
    println!("  RMS:  {:.3}m", rms(&h_errors_fwd));

    println!("\nSmoothed (horizontal):");
    println!("  p50:  {:.3}m", percentile(&h_errors_smooth, 0.50));
    println!("  p68:  {:.3}m", percentile(&h_errors_smooth, 0.68));
    println!("  p95:  {:.3}m", percentile(&h_errors_smooth, 0.95));
    println!("  RMS:  {:.3}m", rms(&h_errors_smooth));

    println!("\nSmoothed (3D):");
    println!("  p50:  {:.3}m", percentile(&errors_3d_smooth, 0.50));
    println!("  p68:  {:.3}m", percentile(&errors_3d_smooth, 0.68));
    println!("  p95:  {:.3}m", percentile(&errors_3d_smooth, 0.95));
    println!("  RMS:  {:.3}m", rms(&errors_3d_smooth));
}
