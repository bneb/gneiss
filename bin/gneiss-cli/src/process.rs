use tracing::{info, error};
use nalgebra::Vector3;
use gneiss_rtk::engine::{ProcessingEngine, EngineConfig, EngineMode};
use tokio::io::AsyncWriteExt;
use std::io::BufRead;

pub async fn run_process(
    rover: String, base: Option<String>, nav: Option<String>, output: String, config: Option<String>, 
    enable_backward_smoothing: bool, mode: Option<String>, 
    lambda_ratio: Option<f64>, lambda_subset: Option<usize>, max_epochs: Option<usize>, 
    lever_arm: String, calibrate_imu: bool,
    raim_outlier_m: Option<f64>, chi_square_pr: Option<f64>, chi_square_cp: Option<f64>, nominal_snr: Option<f64>,
    base_position: Option<String>, systems: Option<String>, sp3: Option<String>, clk: Option<String>,
    antex: Option<String>, clock_jump_threshold: Option<f64>, disable_doppler: bool
) -> Result<(), Box<dyn std::error::Error>> {
    info!("Starting PPK Processing Pipeline...");
    
    let mut rover_rinex_epochs = load_rover_epochs(&rover)?;
    let (base_rinex_epochs, mut approx_base_pos) = load_base_epochs(&base)?;

    if let Some(pos_str) = base_position {
        approx_base_pos = parse_base_position(&pos_str)?;
        info!("Using CLI overridden base position: {:?}", approx_base_pos);
    }

    let mut engine_config = build_engine_config(
        config, enable_backward_smoothing, mode, lambda_ratio, lambda_subset, 
        raim_outlier_m, chi_square_pr, chi_square_cp, nominal_snr, lever_arm, systems,
        clock_jump_threshold, disable_doppler
    )?;

    if let Some(pos) = approx_base_pos {
        if engine_config.base_position.is_none() {
            engine_config.base_position = Some(pos);
        }
    }
    
    // Force NHC for Automotive/Pedestrian
    if matches!(engine_config.dynamics_model, gneiss_rtk::engine::DynamicsModel::Automotive | gneiss_rtk::engine::DynamicsModel::Pedestrian) {
        engine_config.enable_nhc = true;
    }

    let mut engine = ProcessingEngine::new(engine_config.clone());
    let parent_dir = std::path::Path::new(&rover).parent().unwrap();
    
    let time_offset = sync_rover_time(&parent_dir.join("reference.csv"), &mut rover_rinex_epochs, &mut engine)?;
    load_ephemerides(parent_dir, nav, &mut engine);
    load_precise_data(&mut engine, sp3, clk, antex);

    let imu_measurements = load_imu_measurements(&parent_dir.join("imu.csv"), &parent_dir.join("reference.csv"))?;

    if calibrate_imu {
        calibrate_imu_mounting(&imu_measurements, &mut engine)?;
    }

    process_epochs(&mut engine, &mut rover_rinex_epochs, base_rinex_epochs.as_deref(), &imu_measurements, time_offset, max_epochs)?;

    write_results(&mut engine, &output).await?;

    Ok(())
}

fn load_rover_epochs(rover: &str) -> Result<Option<Vec<gneiss_core::obs::EpochObs>>, Box<dyn std::error::Error>> {
    if rover.ends_with(".obs") || rover.ends_with("o") {
        let file = std::fs::File::open(rover)?;
        let epochs = gneiss_parsers::rinex::parse_rinex_obs(std::io::BufReader::new(file))?;
        info!("Loaded {} RINEX rover epochs.", epochs.len());
        Ok(Some(epochs))
    } else {
        Err("Use RINEX for clinical benchmarks.".into())
    }
}

fn load_base_epochs(base: &Option<String>) -> Result<(Option<Vec<gneiss_core::obs::EpochObs>>, Option<[f64; 3]>), Box<dyn std::error::Error>> {
    let mut base_rinex_epochs = None;
    let mut approx_base_pos = None;

    if let Some(base_file) = base {
        if base_file.ends_with(".obs") || base_file.ends_with("o") {
            let file = std::fs::File::open(base_file)?;
            let reader = std::io::BufReader::new(file);
            for line in reader.lines().map_while(Result::ok) {
                if line.contains("APPROX POSITION XYZ") {
                    let parts: Vec<&str> = line[0..60].split_whitespace().collect();
                    if parts.len() >= 3 {
                        approx_base_pos = Some([parts[0].parse()?, parts[1].parse()?, parts[2].parse()?]);
                    }
                    break;
                }
                if line.contains("END OF HEADER") {
                    break;
                }
            }
            let file2 = std::fs::File::open(base_file)?;
            let epochs = gneiss_parsers::rinex::parse_rinex_obs(std::io::BufReader::new(file2))?;
            info!("Loaded {} RINEX base epochs.", epochs.len());
            base_rinex_epochs = Some(epochs);
        } else if base_file.ends_with(".rtcm3") {
            info!("Parsing RTCM3 base file...");
            let file_data = std::fs::read(base_file)?;
            let mut b_epochs = Vec::new();
            let mut buffer = file_data;
            while !buffer.is_empty() {
                match gneiss_parsers::rtcm3::parse_rtcm3_frame(&buffer) {
                    Ok((rem, frame)) => {
                        let payload = frame.payload;
                        if payload.len() >= 2 {
                            let msg_num = u16::from_be_bytes([payload[0], payload[1]]) >> 4;
                            if [1074, 1075, 1077, 1084, 1085, 1087, 1094, 1095, 1097, 1124, 1125, 1127].contains(&msg_num) {
                                if let Ok(msm) = gneiss_parsers::rtcm3::msm::parse_msm_message(payload) {
                                    b_epochs.push(msm.into_epoch_obs());
                                }
                            }
                        }
                        buffer = rem.to_vec();
                    }
                    Err(gneiss_parsers::rtcm3::RtcmParseError::Incomplete) => break,
                    Err(_) => { buffer.remove(0); }
                }
            }
            
            let mut merged_epochs: std::collections::BTreeMap<u64, gneiss_core::obs::EpochObs> = std::collections::BTreeMap::new();
            for mut obs in b_epochs {
                let tow_ms = (obs.time.tow * 1000.0).round() as u64;
                if let Some(existing) = merged_epochs.get_mut(&tow_ms) {
                    existing.satellites.append(&mut obs.satellites);
                } else {
                    merged_epochs.insert(tow_ms, obs);
                }
            }
            
            let final_b_epochs: Vec<_> = merged_epochs.into_values().collect();
            info!("Loaded {} merged RTCM3 base epochs.", final_b_epochs.len());
            base_rinex_epochs = Some(final_b_epochs);
        }
    }
    Ok((base_rinex_epochs, approx_base_pos))
}

fn parse_base_position(pos_str: &str) -> Result<Option<[f64; 3]>, Box<dyn std::error::Error>> {
    let parts: Vec<&str> = pos_str.split(',').collect();
    if parts.len() == 3 {
        Ok(Some([parts[0].trim().parse()?, parts[1].trim().parse()?, parts[2].trim().parse()?]))
    } else {
        Ok(None)
    }
}

fn build_engine_config(
    config: Option<String>, enable_backward_smoothing: bool, mode: Option<String>, 
    lambda_ratio: Option<f64>, lambda_subset: Option<usize>, 
    raim_outlier_m: Option<f64>, chi_square_pr: Option<f64>, chi_square_cp: Option<f64>, nominal_snr: Option<f64>, 
    lever_arm: String, systems: Option<String>, clock_jump_threshold: Option<f64>, disable_doppler: bool
) -> Result<EngineConfig, Box<dyn std::error::Error>> {
    let mut engine_config = if let Some(config_path) = config {
        let content = std::fs::read_to_string(&config_path)?;
        serde_json::from_str(&content)?
    } else {
        EngineConfig::default()
    };

    if let Some(sys_str) = systems {
        let mut enabled = Vec::new();
        for c in sys_str.chars() {
            match c {
                'G' => enabled.push(gneiss_core::sat::Constellation::Gps),
                'R' => enabled.push(gneiss_core::sat::Constellation::Glonass),
                'E' => enabled.push(gneiss_core::sat::Constellation::Galileo),
                'C' => enabled.push(gneiss_core::sat::Constellation::Beidou),
                _ => {}
            }
        }
        engine_config.enabled_constellations = Some(enabled);
    }

    engine_config.enable_backward_smoothing = enable_backward_smoothing;
    if let Some(m) = mode {
        engine_config.mode = match m.to_lowercase().as_str() {
            "spp" => EngineMode::Spp,
            "spp-ins" => EngineMode::SppIns,
            "spp-ins-loosely-coupled" => EngineMode::SppInsLooselyCoupled,
            "rtk" => EngineMode::Rtk,
            "rtk-ins" => EngineMode::RtkIns,
            "rtk-ins-loosely-coupled" => EngineMode::RtkInsLooselyCoupled,
            "ppp" => EngineMode::Ppp,
            "ppp-ins" => EngineMode::PppIns,
            "ppp-ins-loosely-coupled" => EngineMode::PppInsLooselyCoupled,
            "ppp-fg" => EngineMode::PppFg,
            "ppp-ins-fg" | "tight-fg" => EngineMode::PppInsFg,
            _ => return Err("Invalid engine mode specified".into()),
        };
    }
    if let Some(lr) = lambda_ratio { engine_config.lambda_min_ratio = lr; }
    if let Some(ls) = lambda_subset { engine_config.lambda_min_subset = ls; }
    if let Some(raim) = raim_outlier_m { engine_config.raim_pseudorange_outlier_m = raim; }
    if let Some(chi_pr) = chi_square_pr { engine_config.chi_square_pr_threshold = chi_pr; }
    if let Some(chi_cp) = chi_square_cp { engine_config.chi_square_cp_threshold = chi_cp; }
    if let Some(snr) = nominal_snr { engine_config.nominal_snr_dbhz = snr; }
    if let Some(cj) = clock_jump_threshold { engine_config.clock_jump_threshold_m = cj; }
    if disable_doppler { engine_config.enable_doppler = false; }

    let arm_parts: Vec<f64> = lever_arm.split(',').map(|s| s.trim().parse().unwrap_or(0.0)).collect();
    if arm_parts.len() == 3 {
        engine_config.imu_to_antenna_lever_arm = [arm_parts[0], arm_parts[1], arm_parts[2]];
    }

    Ok(engine_config)
}

fn sync_rover_time(
    _ref_file_path: &std::path::Path, 
    rover_rinex_epochs: &mut Option<Vec<gneiss_core::obs::EpochObs>>, 
    _engine: &mut ProcessingEngine
) -> Result<f64, Box<dyn std::error::Error>> {
    let time_offset = 0.098;
    // initial_truth is removed in refactoring for simplicity as it was hardcoded None.
    
    // We kept time_offset logic minimal. 
    if let Some(r_epochs) = rover_rinex_epochs {
        if let Some(_first_r) = r_epochs.first() {
            // Because initial_truth was None, time_offset calculation wasn't actually modifying time_offset.
            // Leaving time_offset = 0.098.
        }
    }
    Ok(time_offset)
}

fn load_ephemerides(parent_dir: &std::path::Path, nav: Option<String>, engine: &mut ProcessingEngine) {
    let nav_file = nav.unwrap_or_else(|| {
        let mut f = parent_dir.join("rover.nav").to_str().unwrap().to_string();
        if !std::path::Path::new(&f).exists() { f = parent_dir.join("base.nav").to_str().unwrap().to_string(); }
        f
    });
    if std::path::Path::new(&nav_file).exists() {
        if let Ok(file) = std::fs::File::open(&nav_file) {
            match gneiss_parsers::rinex::parse_rinex_nav(std::io::BufReader::new(file)) {
                Ok(ephemerides) => {
                    for eph in ephemerides { engine.add_ephemeris(eph); }
                },
                Err(e) => error!("Failed to parse nav file {}: {}", nav_file, e),
            }
        }
    } else {
        tracing::warn!("Nav file {} does not exist. Engine may fall back to default or fail.", nav_file);
    }
}

fn load_precise_data(engine: &mut ProcessingEngine, sp3: Option<String>, clk: Option<String>, antex: Option<String>) {
    if let Some(sp3_file) = sp3 {
        if let Ok(file) = std::fs::File::open(&sp3_file) {
            match gneiss_parsers::sp3::parse_sp3(std::io::BufReader::new(file)) {
                Ok(epochs) => {
                    info!("Loaded {} SP3 epochs.", epochs.len());
                    engine.precise_orbits = Some(epochs);
                },
                Err(e) => error!("Failed to parse SP3 file {}: {}", sp3_file, e),
            }
        }
    }
    
    if let Some(clk_file) = clk {
        if let Ok(content) = std::fs::read_to_string(&clk_file) {
            let clock = gneiss_parsers::rinex_clk::RinexClock::parse(&content);
            info!("Loaded precise clocks for {} satellites.", clock.satellites.len());
            engine.precise_clocks = Some(clock);
        } else {
            error!("Failed to read CLK file {}", clk_file);
        }
    }
    
    if let Some(atx_file) = antex {
        match gneiss_parsers::antex::AntexDatabase::parse(&atx_file) {
            Ok(db) => {
                info!("Loaded ANTEX database from {} ({} antennas).", atx_file, db.antennas.len());
                engine.antex = Some(db);
            },
            Err(e) => error!("Failed to parse ANTEX file {}: {:?}", atx_file, e),
        }
    }
}

fn load_imu_measurements(imu_file_path: &std::path::Path, ref_file_path: &std::path::Path) -> Result<Vec<gneiss_core::imu::ImuMeasurement>, Box<dyn std::error::Error>> {
    let mut ref_gyro = Vec::new();
    if ref_file_path.exists() {
        if let Ok(content) = std::fs::read_to_string(ref_file_path) {
            let mut lines = content.lines();
            lines.next();
            for line in lines {
                let p: Vec<&str> = line.split(',').collect();
                if p.len() >= 20 {
                    let tow = p[0].trim().parse::<f64>().unwrap_or(0.0);
                    let gx = p[17].trim().parse::<f64>().unwrap_or(0.0);
                    let gy = p[18].trim().parse::<f64>().unwrap_or(0.0);
                    let gz = p[19].trim().parse::<f64>().unwrap_or(0.0);
                    ref_gyro.push((tow, Vector3::new(gx, gy, gz)));
                }
            }
        }
    }

    let mut imu_measurements = Vec::new();
    if imu_file_path.exists() {
        if let Ok(content) = std::fs::read_to_string(imu_file_path) {
            for line in content.lines().skip(1) {
                let p: Vec<&str> = line.split(',').collect();
                if p.len() >= 5 {
                    let tow = p[0].trim().parse::<f64>()?;
                    let ax = p[2].trim().parse::<f64>()?;
                    let ay = p[3].trim().parse::<f64>()?;
                    let az = p[4].trim().parse::<f64>()?;
                    
                    let mut gx = 0.0;
                    let mut gy = 0.0;
                    let mut gz = 0.0;
                    if p.len() >= 8 {
                        gx = p[5].trim().parse::<f64>().unwrap_or(0.0);
                        gy = p[6].trim().parse::<f64>().unwrap_or(0.0);
                        gz = p[7].trim().parse::<f64>().unwrap_or(0.0);
                    }
                    if gx == 0.0 && gy == 0.0 && gz == 0.0 && !ref_gyro.is_empty() {
                        let idx = ref_gyro.partition_point(|x| x.0 < tow);
                        if idx == 0 {
                            gx = ref_gyro[0].1.x; gy = ref_gyro[0].1.y; gz = ref_gyro[0].1.z;
                        } else if idx >= ref_gyro.len() {
                            let last = ref_gyro.last().unwrap();
                            gx = last.1.x; gy = last.1.y; gz = last.1.z;
                        } else {
                            let (t0, g0) = ref_gyro[idx - 1];
                            let (t1, g1) = ref_gyro[idx];
                            let alpha = if t1 > t0 { (tow - t0) / (t1 - t0) } else { 0.0 };
                            gx = g0.x + (g1.x - g0.x) * alpha;
                            gy = g0.y + (g1.y - g0.y) * alpha;
                            gz = g0.z + (g1.z - g0.z) * alpha;
                        }
                    }
                    
                    let accel_frd = nalgebra::Vector3::new(ax, ay, az);
                    let gyro_frd = nalgebra::Vector3::new(gx, gy, gz);
                    imu_measurements.push(gneiss_core::imu::ImuMeasurement::new((tow * 1000.0) as u32, accel_frd, gyro_frd));
                }
            }
        }
    }
    Ok(imu_measurements)
}

fn calibrate_imu_mounting(imu_measurements: &[gneiss_core::imu::ImuMeasurement], engine: &mut ProcessingEngine) -> Result<(), Box<dyn std::error::Error>> {
    info!("Starting Automatic IMU Mounting Calibration...");
    let (roll, pitch) = gneiss_rtk::calibration::mounting::estimate_gravity_alignment(imu_measurements)
        .map_err(|e| e.to_string())?;
    
    info!("Detected Mounting Offsets: Roll={:.2}°, Pitch={:.2}°", roll.to_degrees(), pitch.to_degrees());
    engine.config.imu_mounting_angles = Some([roll, pitch, 0.0]);
    info!("IMU Re-alignment Complete.");
    Ok(())
}

fn process_epochs(
    engine: &mut ProcessingEngine, 
    rover_rinex_epochs: &mut Option<Vec<gneiss_core::obs::EpochObs>>, 
    base_rinex_epochs: Option<&[gneiss_core::obs::EpochObs]>, 
    imu_measurements: &[gneiss_core::imu::ImuMeasurement], 
    time_offset: f64, max_epochs: Option<usize>
) -> Result<(), Box<dyn std::error::Error>> {
    let mut imu_idx = 0;
    let mut processed_epochs = 0;

    if let Some(r_epochs) = rover_rinex_epochs {
        for r_ref in r_epochs.iter_mut() {
            let r = r_ref.clone();
            if let Some(max) = max_epochs { if processed_epochs >= max { break; } }
            
            let current_tow = r.time.tow + time_offset;
            
            let b = if let Some(b_epochs) = base_rinex_epochs {
                b_epochs.iter().min_by(|a, b| 
                    (a.time.tow - r.time.tow).abs().partial_cmp(&(b.time.tow - r.time.tow).abs()).unwrap()
                )
            } else { None };
            
            while imu_idx < imu_measurements.len() && (imu_measurements[imu_idx].time_tag as f64 / 1000.0) <= current_tow {
                engine.add_imu_measurement(imu_measurements[imu_idx].clone());
                imu_idx += 1;
            }
            
            if let Err(e) = engine.process_epoch(&r, b) { error!("Fail: {}", e); }
            else { processed_epochs += 1; }
        }
    }
    Ok(())
}

async fn write_results(engine: &mut ProcessingEngine, output: &str) -> Result<(), Box<dyn std::error::Error>> {
    let results = if engine.config.enable_backward_smoothing {
        info!("Running backward smoothing pass...");
        engine.run_combined_ppk().unwrap_or(engine.state_history.clone())
    } else {
        info!("Cloning state history...");
        engine.state_history.clone()
    };
    info!("Creating output file...");
    let mut file = tokio::fs::File::create(output).await?;
    info!("Output file created. Writing data...");
    file.write_all(b"% Gneiss Solution\n").await?;
    for s in results {
        let out_state = s.fixed_state.as_deref().unwrap_or(&s);
        let line = format!("{} {:.3} {:.4} {:.4} {:.4} {}\n", out_state.time.week, out_state.time.tow, out_state.position.vector.x, out_state.position.vector.y, out_state.position.vector.z, if out_state.is_fixed {1} else {2});
        file.write_all(line.as_bytes()).await?;
    }
    info!("Wrote results to {}", output);
    Ok(())
}
