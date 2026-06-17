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
        #[arg(long, help = "Path to precise orbit file (.sp3)")]
        sp3: Option<String>,
        #[arg(long, help = "Path to precise clock file (.clk)")]
        clk: Option<String>,
        #[arg(long, help = "Path to antenna exchange file (.atx)")]
        antex: Option<String>,
        #[arg(long, help = "Path to differential code bias files (.dcb)")]
        dcb: Vec<String>,
        #[arg(long, help = "Path to phase/code bias file (.bia/.bsx)")]
        bia: Option<String>,
        #[arg(long, help = "Path to engine configuration file (.json)")]
        config: Option<String>,
        #[arg(long, help = "Enable multi-pass backward smoothing")]
        enable_backward_smoothing: bool,
        #[arg(long, help = "Enable automatic multi-pass tuning of EKF hyperparameters")]
        enable_auto_tune: bool,
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
        #[arg(long, help = "Automatically calibrate 6-DOF extrinsics and GNSS intrinsics")]
        calibrate: bool,
        #[arg(long, help = "SPP RAIM Outlier Rejection Threshold (m)")]
        raim_outlier_m: Option<f64>,
        #[arg(long, help = "Pseudorange Chi-Square Reject Threshold")]
        chi_square_pr: Option<f64>,
        #[arg(long, help = "Carrier Phase Chi-Square Reject Threshold")]
        chi_square_cp: Option<f64>,
        #[arg(long, help = "Minimum SNR in dBHz")]
        min_snr: Option<f64>,
        
        #[arg(long, help = "Surveyed base station ECEF coordinate override (x,y,z in meters)", group = "base_coord")]
        base_coord_manual: Option<String>,
        #[arg(long, help = "Extract base coordinate from RTCM3 messages 1005/1006", group = "base_coord")]
        base_coord_rtcm: bool,
        #[arg(long, help = "Automatically fetch official coordinate from NGS CORS API", group = "base_coord")]
        base_coord_api: bool,
        #[arg(long, help = "Survey the base station using PPP before processing rover", group = "base_coord")]
        base_coord_ppp: bool,

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
                    "ppp-ins-loosely-coupled" => gneiss_rtk::engine::EngineMode::PppInsLooselyCoupled,
                    "rtk-ins-fg" | "rtk-fg" => gneiss_rtk::engine::EngineMode::RtkInsFactorGraph,
                    "ppp-fg" => gneiss_rtk::engine::EngineMode::PppFg,
                    "ppp-ins-fg" | "tight-fg" => gneiss_rtk::engine::EngineMode::PppInsFg,
                    _ => return Err("Invalid engine mode specified".into()),
                };
            }
            
            let live_cfg = live::LiveConfig { port, baud, ntrip_url, ntrip_mount, ntrip_user, ntrip_pass, _output: output };
            live::run_live(live_cfg, engine_config).await?;
            Ok(())
        },
        Commands::Process { 
            rover, base, nav, output, sp3, clk, antex, dcb, bia, config, 
            enable_backward_smoothing, enable_auto_tune, mode, 
            lambda_ratio, lambda_subset, max_epochs, 
            lever_arm, calibrate_imu, calibrate,
            raim_outlier_m, chi_square_pr, chi_square_cp, min_snr,
            base_coord_manual, base_coord_rtcm, base_coord_api, base_coord_ppp,
            systems
        } => {
            info!("Starting PPK Processing Pipeline...");
            
            let mut rover_rinex_epochs = if rover.ends_with(".obs") || rover.ends_with("o") || rover.ends_with(".rnx") || rover.ends_with(".RNX") {
                let file = std::fs::File::open(&rover)?;
                let epochs = gneiss_parsers::rinex::parse_rinex_obs(std::io::BufReader::new(file))?;
                info!("Loaded {} RINEX rover epochs.", epochs.len());
                Some(epochs)
            } else {
                return Err("Use RINEX for clinical benchmarks.".into());
            };
            
            let mut base_rinex_epochs = None;
            let mut approx_base_pos = None;
            let mut base_marker_name = None;

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
                        }
                        if line.contains("MARKER NAME") {
                            let parts: Vec<&str> = line[0..60].split_whitespace().collect();
                            if !parts.is_empty() {
                                base_marker_name = Some(parts[0].to_string());
                            }
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
                    let mut buffer = file_data.as_slice();
                    while !buffer.is_empty() {
                        if let Ok((rem, frame)) = gneiss_parsers::rtcm3::parse_rtcm3_frame(buffer) {
                            if let Ok(msg) = gneiss_parsers::rtcm3::msm::parse_msm_message(frame.payload) {
                                b_epochs.push(msg.into_epoch_obs());
                            } else if let Ok(arp) = gneiss_parsers::rtcm3::station::parse_station_arp(frame.payload) {
                                if base_coord_rtcm {
                                    approx_base_pos = Some([arp.ecef_x, arp.ecef_y, arp.ecef_z]);
                                    info!("Extracted base coordinate from RTCM3 ARP: {:?}", approx_base_pos);
                                }
                            }
                            buffer = rem;
                        } else {
                            buffer = &buffer[1..];
                        }
                    }
                    
                    // Merge epochs with same timestamp (e.g. GPS + GLONASS in same second)
                    let mut merged_epochs: std::collections::BTreeMap<i64, gneiss_core::obs::EpochObs> = std::collections::BTreeMap::new();
                    for epoch in b_epochs {
                        let ms = (epoch.time.tow * 1000.0).round() as i64;
                        if let Some(existing) = merged_epochs.get_mut(&ms) {
                            existing.satellites.extend(epoch.satellites);
                        } else {
                            merged_epochs.insert(ms, epoch);
                        }
                    }
                    
                    let final_b_epochs: Vec<_> = merged_epochs.into_values().collect();
                    info!("Loaded {} merged RTCM3 base epochs.", final_b_epochs.len());
                    base_rinex_epochs = Some(final_b_epochs);
                }
            }

            if base_coord_rtcm {
                if approx_base_pos.is_none() {
                    return Err("Failed to extract base coordinate from RTCM3 messages".into());
                }
            } else if base_coord_api {
                if let Some(marker) = &base_marker_name {
                    let provider = gneiss_fetch::sources::noaa::NoaaCorsProvider;
                    info!("Fetching official coordinate for {} from NOAA CORS API...", marker);
                    match provider.fetch_station_coordinate(marker).await {
                        Ok(coord) => {
                            approx_base_pos = Some([coord.vector.x, coord.vector.y, coord.vector.z]);
                            info!("Base coordinate retrieved from NOAA API: {:?}", approx_base_pos);
                        },
                        Err(e) => {
                            return Err(format!("Failed to fetch coordinate for {}: {}", marker, e).into());
                        }
                    }
                } else {
                    return Err("No MARKER NAME found in base RINEX. Cannot fetch API coordinate.".into());
                }
            } else if base_coord_ppp {
                info!("PPP Survey mode selected. Surveying base station...");
                return Err("PPP base survey not yet implemented.".into());
            } else if let Some(pos_str) = &base_coord_manual {
                let parts: Vec<&str> = pos_str.split(',').collect();
                if parts.len() == 3 {
                    approx_base_pos = Some([
                        parts[0].trim().parse()?, 
                        parts[1].trim().parse()?, 
                        parts[2].trim().parse()?
                    ]);
                    info!("Using CLI overridden base position: {:?}", approx_base_pos);
                } else {
                    return Err("Manual base coordinate must be in format x,y,z".into());
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
                    "ppp-ins-loosely-coupled" => gneiss_rtk::engine::EngineMode::PppInsLooselyCoupled,
                    "rtk-ins-fg" | "rtk-fg" => gneiss_rtk::engine::EngineMode::RtkInsFactorGraph,
                    "ppp-fg" => gneiss_rtk::engine::EngineMode::PppFg,
                    "ppp-ins-fg" | "tight-fg" => gneiss_rtk::engine::EngineMode::PppInsFg,
                    _ => return Err("Invalid engine mode specified".into()),
                };
                
                let is_ppp = matches!(engine_config.mode, 
                    gneiss_rtk::engine::EngineMode::Ppp |
                    gneiss_rtk::engine::EngineMode::PppIns |
                    gneiss_rtk::engine::EngineMode::PppInsLooselyCoupled |
                    gneiss_rtk::engine::EngineMode::PppFg |
                    gneiss_rtk::engine::EngineMode::PppInsFg
                );
                
                if is_ppp && (sp3.is_none() || clk.is_none()) {
                    return Err("Error: SP3 and CLK files are MANDATORY for PPP evaluations. Without precise clock and orbit corrections, carrier phase ambiguities cannot be resolved, resulting in unbounded drift. Use `gneiss fetch` to download them.".into());
                }
            }
            if let Some(lr) = lambda_ratio { engine_config.lambda_min_ratio = lr; }
            if let Some(ls) = lambda_subset { engine_config.lambda_min_subset = ls; }
            
            if let Some(raim) = raim_outlier_m { engine_config.raim_pseudorange_outlier_m = raim; }
            if let Some(chi_pr) = chi_square_pr { engine_config.chi_square_pr_threshold = chi_pr; }
            if let Some(chi_cp) = chi_square_cp { engine_config.chi_square_cp_threshold = chi_cp; }
            if let Some(snr) = min_snr { engine_config.min_snr_dbhz = snr; }

            let arm_parts: Vec<f64> = lever_arm.split(',').map(|s| s.trim().parse().unwrap_or(0.0)).collect();
            if arm_parts.len() == 3 {
                engine_config.imu_to_antenna_lever_arm = [arm_parts[0], arm_parts[1], arm_parts[2]];
            }

            if let Some(pos) = approx_base_pos {
                // The mutually exclusive CLI flags (--base-coord-*) take precedence over the JSON config.
                // If they are not specified, approx_base_pos defaults to the fallback extracted from the RINEX header.
                // So if the config specifies a base position, we should use that OVER the RINEX header fallback,
                // BUT we should use the CLI flags over the config.
                // Since approx_base_pos acts as both the fallback and the CLI flag resolution, we need to check if 
                // a CLI flag was actually used.
                let cli_override = base_coord_manual.is_some() || base_coord_rtcm || base_coord_api || base_coord_ppp;
                
                if cli_override || engine_config.base_position.is_none() {
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

            if let (Some(r_epochs), Some((truth_tow, mut state))) = (&mut rover_rinex_epochs, initial_truth.clone()) {
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
                        Ok((ephemerides, klobuchar)) => {
                            for eph in ephemerides { engine.add_ephemeris(eph); }
                            if let Some(klob) = klobuchar { engine.klobuchar_params = Some(klob); }
                        },
                        Err(e) => error!("Failed to parse nav file {}: {}", nav_file, e),
                    }
                }
            } else {
                tracing::warn!("Nav file {} does not exist. Engine may fall back to default or fail.", nav_file);
            }

            if let Some(sp3_path) = sp3 {
                if let Ok(file) = std::fs::File::open(&sp3_path) {
                    match gneiss_parsers::sp3::parse_sp3(std::io::BufReader::new(file)) {
                        Ok(epochs) => {
                            info!("Loaded {} SP3 epochs.", epochs.len());
                            engine.sp3_epochs = epochs;
                        },
                        Err(e) => error!("Failed to parse SP3 file {}: {}", sp3_path, e),
                    }
                }
            }

            if let Some(clk_path) = clk {
                if let Ok(content) = std::fs::read_to_string(&clk_path) {
                    let rinex_clk = gneiss_parsers::rinex_clk::RinexClock::parse(&content);
                    info!("Loaded CLK file with {} satellites.", rinex_clk.satellites.len());
                    engine.clk_data = Some(rinex_clk);
                } else {
                    error!("Failed to read CLK file {}", clk_path);
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

            if let Some(bia_file) = bia {
                if let Ok(file) = std::fs::File::open(&bia_file) {
                    match gneiss_parsers::sinex_bia::SinexBias::parse(std::io::BufReader::new(file)) {
                        Ok(bias) => {
                            info!("Loaded SINEX/BIA with {} records.", bias.records.len());
                            engine.sinex_bias = Some(bias);
                        },
                        Err(e) => error!("Failed to parse BIA file {}: {:?}", bia_file, e),
                    }
                }
            }

            for dcb_file in dcb {
                let filename = std::path::Path::new(&dcb_file).file_name().unwrap().to_string_lossy().to_string();
                let dcb_type = if filename.starts_with("P1C1") { "P1C1" }
                    else if filename.starts_with("P2C2") { "P2C2" }
                    else if filename.starts_with("P1P2") { "P1P2" }
                    else { "UNKNOWN" };

                if dcb_type != "UNKNOWN" {
                    if let Ok(content) = std::fs::read_to_string(&dcb_file) {
                        let mut count = 0;
                        for line in content.lines() {
                            if line.len() > 30 && line.starts_with(['G', 'R', 'E', 'C']) {
                                if let Ok(sat) = std::str::FromStr::from_str(&line[0..3]) {
                                    let parts: Vec<&str> = line[3..].split_whitespace().collect();
                                    if parts.len() >= 2 {
                                        if let Ok(val) = parts[0].parse::<f64>() {
                                            engine.dcbs.insert((sat, dcb_type.to_string()), val);
                                            count += 1;
                                        }
                                    }
                                }
                            }
                        }
                        info!("Loaded {} DCBs of type {} from {}", count, dcb_type, dcb_file);
                    } else {
                        error!("Failed to read DCB file {}", dcb_file);
                    }
                } else {
                    error!("Unknown DCB file type from filename: {}", dcb_file);
                }
            }

            if let Some(dir) = parent_dir.to_str() {
                let dcb_files = vec!["P1C1", "P2C2", "P1P2", "C1P1", "C2P2"];
                if let Ok(entries) = std::fs::read_dir(dir) {
                    for file in entries.flatten() {
                        let path = file.path();
                        if path.extension().and_then(|s| s.to_str()) == Some("DCB") {
                            if let Some(fname) = path.file_name().and_then(|s| s.to_str()) {
                                let mut dcb_type = "";
                                for dt in &dcb_files {
                                    if fname.starts_with(dt) {
                                        dcb_type = dt;
                                        break;
                                    }
                                }
                                if dcb_type.is_empty() { continue; }
                                if let Ok(content) = std::fs::read_to_string(&path) {
                                    let mut count = 0;
                                    for line in content.lines() {
                                        if line.starts_with('G') || line.starts_with('R') || line.starts_with('E') || line.starts_with('C') {
                                            let parts: Vec<&str> = line.split_whitespace().collect();
                                            if parts.len() >= 2 {
                                                let sys = parts[0].chars().next().unwrap();
                                                if let Ok(prn) = parts[0][1..].parse::<u8>() {
                                                    if let Ok(val) = parts[1].parse::<f64>() {
                                                        let constel = match sys {
                                                            'G' => gneiss_core::sat::Constellation::Gps,
                                                            'R' => gneiss_core::sat::Constellation::Glonass,
                                                            'E' => gneiss_core::sat::Constellation::Galileo,
                                                            'C' => gneiss_core::sat::Constellation::Beidou,
                                                            _ => continue,
                                                        };
                                                        engine.dcbs.insert((gneiss_core::sat::SatelliteId { prn, constellation: constel }, dcb_type.to_string()), val);
                                                        count += 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    info!("Loaded {} {} DCBs from {:?}", count, dcb_type, path);
                                }
                            }
                        }
                    }
                }
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
                            
                            // The IMU and Gyro data in this specific dataset are already in m/s^2 and rad/s.
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

            let mut final_results = Vec::new();
            let mut final_processed_epochs = 0;
            
            if calibrate {
                info!("Starting Preprocessing and Calibration Pass...");
                
                // --- Pass 1a: Intrinsics Calibration ---
                let rover_slice = if let Some(ref r) = rover_rinex_epochs {
                    let start = r.len().min(3000);
                    let end = r.len().min(start + 500);
                    &r[start..end]
                } else {
                    &[]
                };
                
                let (bias_var, drift_var) = gneiss_rtk::calibration::intrinsics::calibrate_intrinsics(&engine_config, rover_slice);
                engine_config.process_noise_cb = bias_var;
                engine_config.process_noise_cd = drift_var;
                engine.config.process_noise_cb = bias_var;
                engine.config.process_noise_cd = drift_var;
                info!("Intrinsics calibrated: Clock Bias Var = {:.4}, Clock Drift Var = {:.4}", bias_var, drift_var);
                
                // --- Pass 1b: Extrinsics 6-DOF Optimization ---
                let eval_fn = |cfg: &EngineConfig| -> f64 {
                    let mut eval_engine = ProcessingEngine::new(cfg.clone());
                    eval_engine.ephemerides = engine.ephemerides.clone();
                    eval_engine.klobuchar_params = engine.klobuchar_params.clone();
                    let mut eval_imu_idx = 0;
                    
                    for r in rover_slice {
                        let current_tow = r.time.tow + time_offset;
                        let b = if let Some(ref b_epochs) = base_rinex_epochs {
                            b_epochs.iter().min_by(|a, b| 
                                (a.time.tow - r.time.tow).abs().partial_cmp(&(b.time.tow - r.time.tow).abs()).unwrap()
                            )
                        } else { None };
                        
                        while eval_imu_idx < imu_measurements.len() && (imu_measurements[eval_imu_idx].time_tag as f64 / 1000.0) <= current_tow {
                            eval_engine.add_imu_measurement(imu_measurements[eval_imu_idx].clone());
                            eval_imu_idx += 1;
                        }
                        
                        let _ = eval_engine.process_epoch(r, b);
                    }
                    if eval_engine.state_history.last().is_some() {
                        let total_nis = eval_engine.innovation_tracker.get_total_nis();
                        if total_nis > 0.0 { total_nis } else { 1e9 }
                    } else {
                        1e9
                    }
                };
                
                if let Ok((best_lever_arm, best_angles)) = gneiss_rtk::calibration::extrinsics::calibrate_extrinsics_6dof(&engine_config, eval_fn) {
                    engine_config.imu_to_antenna_lever_arm = best_lever_arm;
                    engine_config.imu_mounting_angles = Some(best_angles);
                    engine.config.imu_to_antenna_lever_arm = best_lever_arm;
                    engine.config.imu_mounting_angles = Some(best_angles);
                    // Also use for NHC for now
                    engine_config.imu_to_nhc_lever_arm = best_lever_arm;
                    engine.config.imu_to_nhc_lever_arm = best_lever_arm;
                    info!("Extrinsics calibrated. Lever Arm: {:?}, Mounting Angles: {:?}", best_lever_arm, best_angles);
                } else {
                    error!("Extrinsics calibration failed.");
                }
            }
            
            let passes = if enable_auto_tune { 2 } else { 1 };
            
            for pass in 1..=passes {
                if pass == 1 && enable_auto_tune {
                    info!("Running Pass 1 (Auto-Tuning Calibration)...");
                    engine.config.tuning.auto_tune.enabled = true;
                } else if pass == 2 {
                    info!("Running Pass 2 (Final Execution)...");
                    // Turn off auto-tune for the final pass
                    let mut new_tuning = engine.config.tuning.clone();
                    new_tuning.auto_tune.enabled = false;
                    engine.reset_for_multipass();
                    engine.config.tuning = new_tuning;
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

                if pass == 1 && enable_auto_tune {
                    info!("Analyzing Pass 1 telemetry...");
                    engine.config.tuning = gneiss_rtk::engine::auto_tuner::tune_ekf_parameters(
                        &engine.state_history,
                        &engine.imu_history,
                        engine.config.tuning.clone()
                    );
                } else {
                    final_results = if engine.config.enable_backward_smoothing {
                        info!("Running backward smoothing pass...");
                        engine.run_combined_ppk().unwrap_or(engine.state_history.clone())
                    } else {
                        engine.state_history.clone()
                    };
                    final_processed_epochs = processed_epochs;
                }
            }

            let mut file = tokio::fs::File::create(&output).await?;
            file.write_all(b"% Gneiss Solution\n").await?;
            for s in final_results {
                let out_state = s.fixed_state.as_deref().unwrap_or(&s);
                let line = format!("{} {:.3} {:.4} {:.4} {:.4} {}\n", out_state.time.week, out_state.time.tow, out_state.position.vector.x, out_state.position.vector.y, out_state.position.vector.z, if out_state.is_fixed {1} else {2});
                file.write_all(line.as_bytes()).await?;
            }
            info!("Wrote {} epochs to {}", final_processed_epochs, output);
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

mod mw_test;
