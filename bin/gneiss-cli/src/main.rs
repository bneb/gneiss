use tracing::{info, error};
use nalgebra::Vector3;
use clap::{Parser, Subcommand};
use gneiss_core::coords::{Coordinate, Datum, Frame};
use gneiss_rtk::engine::{ProcessingEngine, EngineConfig};
use tokio::io::{AsyncWriteExt};

mod evaluator;
mod live;

#[derive(Parser, Debug)]
#[command(name = "gneiss", about = "Gneiss Navigation CLI")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Run the engine in real-time using a serial port and NTRIP caster
    Live {
        #[arg(short, long, help = "Serial port device (e.g., /dev/ttyACM0, COM3)")]
        port: String,
        #[arg(short, long, help = "Baud rate", default_value_t = 460800)]
        baud: u32,
        #[arg(long, help = "NTRIP Caster URL (e.g., rtk2go.com)")]
        ntrip_url: Option<String>,
        #[arg(long, help = "NTRIP Mountpoint")]
        ntrip_mount: Option<String>,
        #[arg(long, help = "NTRIP Username")]
        ntrip_user: Option<String>,
        #[arg(long, help = "NTRIP Password")]
        ntrip_pass: Option<String>,
        #[arg(long, help = "Engine configuration file (.json)")]
        config: Option<String>,
        #[arg(long, help = "Engine mode (spp, spp-ins, rtk, rtk-ins, ppp, ppp-ins)")]
        mode: Option<String>,
        #[arg(short, long, help = "Output trajectory stream to file (.pos)")]
        output: Option<String>,
    },
    /// Process GNSS/IMU raw data to produce a trajectory
    Process {
        #[arg(short, long, help = "Path to rover raw data (.ubx, .obs)")]
        rover: String,
        #[arg(short, long, help = "Path to base station raw data (.rtcm3, .obs)")]
        base: Option<String>,
        #[arg(short, long, help = "Path to ephemeris/nav file (.nav, .rnx)")]
        nav: Option<String>,
        #[arg(short, long, help = "Path to output trajectory file (.pos)")]
        output: String,
        #[arg(long, help = "Path to engine configuration file (.json)")]
        config: Option<String>,
        #[arg(long, help = "Enable multi-pass backward smoothing")]
        enable_backward_smoothing: bool,
        #[arg(long, help = "Engine mode (spp, spp-ins, rtk, rtk-ins, ppp, ppp-ins)")]
        mode: Option<String>,
        #[arg(long, help = "LAMBDA PAR min ratio threshold")]
        lambda_ratio: Option<f64>,
        #[arg(long, help = "LAMBDA PAR minimum subset size")]
        lambda_subset: Option<usize>,
        #[arg(long, help = "Maximum number of epochs to process")]
        max_epochs: Option<usize>,
        #[arg(long, help = "Lever arm from IMU to GNSS antenna (x,y,z in meters)", default_value = "0,0,0")]
        lever_arm: String,
        #[arg(long, help = "Automatically detect and calibrate IMU mounting offsets")]
        calibrate_imu: bool,
        #[arg(long, help = "SPP RAIM Outlier Rejection Threshold (m)")]
        raim_outlier_m: Option<f64>,
        #[arg(long, help = "Pseudorange Chi-Square Reject Threshold")]
        chi_square_pr: Option<f64>,
        #[arg(long, help = "Carrier Phase Chi-Square Reject Threshold")]
        chi_square_cp: Option<f64>,
        #[arg(long, help = "Nominal SNR reference (dB-Hz)")]
        nominal_snr: Option<f64>,
        #[arg(long, help = "Surveyed base station ECEF coordinate override (x,y,z in meters)")]
        base_position: Option<String>,
        #[arg(long, help = "Enabled constellations (e.g. G,R,E,C). Default: all")]
        systems: Option<String>,
    },
    /// Evaluate the error CDFs against a ground truth RTKLIB POS file
    Eval {
        #[arg(short, long, help = "Path to our solution file (.pos)")]
        solution: String,
        #[arg(short, long, help = "Path to ground truth reference file (.pos, .csv)")]
        truth: String,
    },
    /// Fetch global ephemeris and local base station data
    Fetch {
        #[arg(long, help = "Path to rover observation file to determine location and time")]
        rover_obs: String,
        #[arg(long, help = "Data source (e.g., noaa, cddis)")]
        source: String,
        #[arg(long, help = "Output directory for downloaded files", default_value = ".")]
        out_dir: String,
    },
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    dotenv::dotenv().ok();
    tracing_subscriber::fmt::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Live { port, baud, ntrip_url, ntrip_mount, ntrip_user, ntrip_pass, config, mode, output } => {
            let mut engine_config = if let Some(config_path) = config {
                let content = std::fs::read_to_string(&config_path)?;
                serde_json::from_str(&content)?
            } else {
                EngineConfig::default()
            };
            engine_config.enable_backward_smoothing = false;
            if let Some(m) = mode {
                engine_config.mode = match m.to_lowercase().as_str() {
                    "spp" => gneiss_rtk::engine::EngineMode::Spp,
                    "spp-ins" => gneiss_rtk::engine::EngineMode::SppIns,
                    "spp-ins-loosely-coupled" => gneiss_rtk::engine::EngineMode::SppInsLooselyCoupled,
                    "rtk" => gneiss_rtk::engine::EngineMode::Rtk,
                    "rtk-ins" => gneiss_rtk::engine::EngineMode::RtkIns,
                    "rtk-ins-loosely-coupled" => gneiss_rtk::engine::EngineMode::RtkInsLooselyCoupled,
                    "ppp" => gneiss_rtk::engine::EngineMode::Ppp,
                    "ppp-ins" => gneiss_rtk::engine::EngineMode::PppIns,
                    _ => return Err("Invalid engine mode specified".into()),
                };
            }
            
            live::run_live(port, baud, ntrip_url, ntrip_mount, ntrip_user, ntrip_pass, engine_config, output).await?;
            Ok(())
        },
        Commands::Process { 
            rover, base, nav, output, config, 
            enable_backward_smoothing, mode, 
            lambda_ratio, lambda_subset, max_epochs, 
            lever_arm, calibrate_imu,
            raim_outlier_m, chi_square_pr, chi_square_cp, nominal_snr,
            base_position,
            systems
        } => {
            info!("Starting PPK Processing Pipeline...");
            
            let mut rover_rinex_epochs = if rover.ends_with(".obs") || rover.ends_with("o") {
                let file = std::fs::File::open(&rover)?;
                let epochs = gneiss_parsers::rinex::parse_rinex_obs(std::io::BufReader::new(file))?;
                info!("Loaded {} RINEX rover epochs.", epochs.len());
                Some(epochs)
            } else {
                return Err("Use RINEX for clinical benchmarks.".into());
            };
            
            let mut base_rinex_epochs = None;
            let mut approx_base_pos = None;

            if let Some(base_file) = &base {
                if base_file.ends_with(".obs") || base_file.ends_with("o") {
                    let file = std::fs::File::open(base_file)?;
                    let reader = std::io::BufReader::new(file);
                    use std::io::BufRead;
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
                }
            }

            if let Some(pos_str) = base_position {
                let parts: Vec<&str> = pos_str.split(',').collect();
                if parts.len() == 3 {
                    approx_base_pos = Some([
                        parts[0].trim().parse()?, 
                        parts[1].trim().parse()?, 
                        parts[2].trim().parse()?
                    ]);
                    info!("Using CLI overridden base position: {:?}", approx_base_pos);
                }
            }

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
                    "spp" => gneiss_rtk::engine::EngineMode::Spp,
                    "spp-ins" => gneiss_rtk::engine::EngineMode::SppIns,
                    "spp-ins-loosely-coupled" => gneiss_rtk::engine::EngineMode::SppInsLooselyCoupled,
                    "rtk" => gneiss_rtk::engine::EngineMode::Rtk,
                    "rtk-ins" => gneiss_rtk::engine::EngineMode::RtkIns,
                    "rtk-ins-loosely-coupled" => gneiss_rtk::engine::EngineMode::RtkInsLooselyCoupled,
                    "ppp" => gneiss_rtk::engine::EngineMode::Ppp,
                    "ppp-ins" => gneiss_rtk::engine::EngineMode::PppIns,
                    _ => return Err("Invalid engine mode specified".into()),
                };
            }
            if let Some(lr) = lambda_ratio { engine_config.lambda_min_ratio = lr; }
            if let Some(ls) = lambda_subset { engine_config.lambda_min_subset = ls; }
            
            if let Some(raim) = raim_outlier_m { engine_config.raim_pseudorange_outlier_m = raim; }
            if let Some(chi_pr) = chi_square_pr { engine_config.chi_square_pr_threshold = chi_pr; }
            if let Some(chi_cp) = chi_square_cp { engine_config.chi_square_cp_threshold = chi_cp; }
            if let Some(snr) = nominal_snr { engine_config.nominal_snr_dbhz = snr; }

            let arm_parts: Vec<f64> = lever_arm.split(',').map(|s| s.trim().parse().unwrap_or(0.0)).collect();
            if arm_parts.len() == 3 {
                engine_config.imu_to_antenna_lever_arm = [arm_parts[0], arm_parts[1], arm_parts[2]];
            }

            if let Some(pos) = approx_base_pos {
                if engine_config.base_position.is_none() {
                    engine_config.base_position = Some(pos);
                }
            }
            
            let mut engine = ProcessingEngine::new(engine_config.clone());
            let parent_dir = std::path::Path::new(&rover).parent().unwrap();
            
            let mut time_offset = 0.098;
            let initial_truth: Option<(f64, gneiss_rtk::filter::RtkState)> = None;
            let mut ref_gyro: Vec<(f64, Vector3<f64>)> = Vec::new();

            let ref_file_path = parent_dir.join("reference.csv");
            if ref_file_path.exists() {
                if let Ok(content) = std::fs::read_to_string(&ref_file_path) {
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

            if let (Some(r_epochs), Some((truth_tow, mut state))) = (&mut rover_rinex_epochs, initial_truth) {
                if let Some(first_r) = r_epochs.first() {
                    time_offset = truth_tow - first_r.time.tow;
                    info!("Detected Time Offset: {:.3}s. Aligning Observations...", time_offset);
                    
                    state.time = first_r.time; // Prevent negative dt on first epoch
                    engine.current_state = Some(state);
                }
            }

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

            let mut imu_measurements: Vec<gneiss_core::imu::ImuMeasurement> = Vec::new();
            let imu_file_path = parent_dir.join("imu.csv");
            if imu_file_path.exists() {
                if let Ok(content) = std::fs::read_to_string(&imu_file_path) {
                    for line in content.lines().skip(1) {
                        let p: Vec<&str> = line.split(',').collect();
                        if p.len() >= 5 {
                            let tow = p[0].trim().parse::<f64>()?;
                            let ax = p[2].trim().parse::<f64>()?;
                            let ay = p[3].trim().parse::<f64>()?;
                            let az = p[4].trim().parse::<f64>()?;
                            
                            // Interpolate Gyro from reference.csv
                            let mut gx = 0.0;
                            let mut gy = 0.0;
                            let mut gz = 0.0;
                            if !ref_gyro.is_empty() {
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
                            
                            // The IMU and Gyro data in this specific dataset are already pre-aligned to the vehicle FRD frame.
                            let accel_frd = nalgebra::Vector3::new(ax, ay, az);
                            let gyro_frd = nalgebra::Vector3::new(gx, gy, gz);
                            imu_measurements.push(gneiss_core::imu::ImuMeasurement::new((tow * 1000.0) as u32, accel_frd, gyro_frd));
                        }
                    }
                }
            }

            if calibrate_imu {
                info!("Starting Automatic IMU Mounting Calibration...");
                let (roll, pitch) = gneiss_rtk::calibration::mounting::estimate_gravity_alignment(&imu_measurements)
                    .map_err(|e| e.to_string())?;
                
                let roll_deg: f64 = roll.to_degrees();
                let pitch_deg: f64 = pitch.to_degrees();
                info!("Detected Mounting Offsets: Roll={:.2}°, Pitch={:.2}°", roll_deg, pitch_deg);
                engine.config.imu_mounting_angles = Some([roll, pitch, 0.0]);
                info!("IMU Re-alignment Complete.");
            }

            let mut imu_idx = 0;
            let mut processed_epochs = 0;

            if let Some(r_epochs) = &mut rover_rinex_epochs {
                let b_epochs = base_rinex_epochs.as_ref();
                for r_ref in r_epochs.iter_mut() {
                    let r = r_ref.clone();
                    if let Some(max) = max_epochs { if processed_epochs >= max { break; } }
                    
                    let current_tow = r.time.tow + time_offset; // Extract IMU up to the true aligned time
                    
                    let b = if let Some(b_epochs) = b_epochs {
                        b_epochs.iter().min_by(|a, b| 
                            (a.time.tow - r.time.tow).abs().partial_cmp(&(b.time.tow - r.time.tow).abs()).unwrap()
                        ).filter(|be| {
                            let diff = (be.time.tow - r.time.tow).abs();
                            if processed_epochs < 5 {
                                tracing::info!("Matching rover {} with base {}: diff = {}", r.time.tow, be.time.tow, diff);
                            }
                            diff < 1.0
                        })
                    } else { None };
                    
                    while imu_idx < imu_measurements.len() && (imu_measurements[imu_idx].time_tag as f64 / 1000.0) <= current_tow {
                        engine.add_imu_measurement(imu_measurements[imu_idx].clone());
                        imu_idx += 1;
                    }
                    
                    if let Err(e) = engine.process_epoch(&r, b) { error!("Fail: {}", e); }
                    else { processed_epochs += 1; }
                }
            }

            let results = if engine.config.enable_backward_smoothing {
                info!("Running backward smoothing pass...");
                engine.run_combined_ppk().unwrap_or(engine.state_history.clone())
            } else {
                engine.state_history.clone()
            };
            let mut file = tokio::fs::File::create(&output).await?;
            file.write_all(b"% Gneiss Solution\n").await?;
            for s in results {
                let out_state = s.fixed_state.as_deref().unwrap_or(&s);
                let line = format!("{} {:.3} {:.4} {:.4} {:.4} {}\n", out_state.time.week, out_state.time.tow, out_state.position.vector.x, out_state.position.vector.y, out_state.position.vector.z, if out_state.is_fixed {1} else {2});
                file.write_all(line.as_bytes()).await?;
            }
            info!("Wrote {} epochs to {}", processed_epochs, output);
            Ok(())
        },
        Commands::Eval { solution, truth } => {
            evaluator::evaluate(&solution, &truth).map_err(|e| e.into())
        },
        Commands::Fetch { rover_obs, source, out_dir } => {
            info!("Fetching data for {} using source {}", rover_obs, source);
            
            use std::io::{BufRead, BufReader};
            let file = std::fs::File::open(&rover_obs)?;
            let reader = BufReader::new(file);
            
            let mut approx_pos = None;
            let mut header_and_first_epoch = String::new();
            let mut in_header = true;
            let mut post_header_lines = 0;

            for line_res in reader.lines() {
                let line = line_res?;
                header_and_first_epoch.push_str(&line);
                header_and_first_epoch.push('\n');

                if in_header {
                    if line.contains("APPROX POSITION XYZ") {
                        let parts: Vec<&str> = line[0..60].split_whitespace().collect();
                        if parts.len() >= 3 {
                            let ecef = nalgebra::Vector3::new(parts[0].parse()?, parts[1].parse()?, parts[2].parse()?);
                            approx_pos = Some(ecef);
                        }
                    } else if line.contains("END OF HEADER") {
                        in_header = false;
                    }
                } else {
                    post_header_lines += 1;
                    // 100 lines is safely enough to capture the first epoch block in both RINEX 2 and 3
                    if post_header_lines > 100 {
                        break;
                    }
                }
            }

            let epochs = gneiss_parsers::rinex::parse_rinex_obs(std::io::Cursor::new(header_and_first_epoch))
                .map_err(|e| format!("Failed to parse first epoch from RINEX: {}", e))?;
            let first_epoch = epochs.first().ok_or("No epochs found in rover obs file")?;
            let time = first_epoch.time;
            
            let coord = if let Some(ecef) = approx_pos {
                Coordinate::new(ecef, Datum::WGS84, Frame::ECEF, time)
            } else {
                return Err("Could not determine approximate position from RINEX header".into());
            };

            use gneiss_fetch::provider::DataSource;
            use gneiss_fetch::sources::cddis::CddisProvider;
            use gneiss_fetch::sources::noaa::NoaaCorsProvider;
            
            let out_path = std::path::Path::new(&out_dir);
            std::fs::create_dir_all(out_path)?;

            let token = std::env::var("EARTHDATA_TOKEN").ok();
            
            if source.to_lowercase() == "cddis" || source.to_lowercase() == "all" {
                let cddis = CddisProvider { auth_token: token.clone() };
                cddis.fetch_ephemeris(time, out_path).await?;
            }
            
            if source.to_lowercase() == "noaa" || source.to_lowercase() == "all" {
                let noaa = NoaaCorsProvider;
                noaa.fetch_base_obs(coord, time, out_path).await?;
            }

            info!("Fetch completed successfully.");
            Ok(())
        }
    }
}
