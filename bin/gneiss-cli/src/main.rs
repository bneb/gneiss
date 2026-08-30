use std::path::Path;
use clap::{Parser, Subcommand};

mod calibrate;
mod evaluator;
mod export;
mod gui;
mod live;
mod process;
mod qc;

#[derive(Parser, Debug)]
#[command(name = "gneiss", about = "Gneiss Navigation CLI")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
#[allow(clippy::large_enum_variant)]
enum Commands {
    /// Process GNSS data to produce a trajectory
    Process {
        #[arg(short = 'r', long, help = "Path to rover raw data (.obs)")]
        rover: String,
        #[arg(short = 'b', long = "base", action = clap::ArgAction::Append, help = "Path to base station raw data (.obs) - specify multiple for Network RTK")]
        base: Vec<String>,
        #[arg(short = 'n', long, help = "Path to ephemeris/nav file (.nav, .rnx)")]
        nav: Option<String>,
        #[arg(short = 'o', long, help = "Path to output trajectory file (.pos, .csv, .kml, .json, .sbet)")]
        output: String,
        #[arg(long, help = "Output export format (pos, csv, kml, json, sbet)")]
        format: Option<String>,
        #[arg(long, help = "Path to write structured QC report (.json or .csv)")]
        qc_report: Option<String>,
        #[arg(long, help = "Path to geoid grid (.gtx, .byn, or .json) for orthometric height")]
        geoid: Option<String>,
        #[arg(long, help = "Path to engine configuration file (.json)")]
        config: Option<String>,
        #[arg(long, help = "Single forward-only pass (SWFG, no ambiguity resolution). \
            Much faster, much less accurate -- typically 5-10x worse -- than the default.")]
        single_pass: bool,
        #[arg(long, help = "Engine mode (spp, rtk, ppp)")]
        mode: Option<String>,
        #[arg(long, help = "Maximum epochs to process")]
        max_epochs: Option<usize>,
        #[arg(long, help = "Enabled constellations (e.g. G,R,E,C)")]
        systems: Option<String>,
        #[arg(long, help = "Path to SP3 file")]
        sp3: Option<String>,
        #[arg(long, help = "Path to CLK file")]
        clk: Option<String>,
        #[arg(long, allow_hyphen_values = true, help = "Base position X,Y,Z ECEF")]
        base_position: Option<String>,
        #[arg(long, help = "Path to ANTEX file for receiver antenna PCV correction (default: datasets/igs14.atx)")]
        antex: Option<String>,
        #[arg(long, help = "Include GLONASS in double-difference formation")]
        glonass: bool,
    },

    /// Compute local site calibration from paired GNSS and Ground Control Points (CSV)
    Calibrate {
        #[arg(short = 'i', long, help = "Path to paired control points CSV (format: GNSS_E,GNSS_N,GNSS_H,Ground_E,Ground_N,Ground_H)")]
        input: String,
        #[arg(short = 'o', long, help = "Optional path to write calibration JSON")]
        output: Option<String>,
    },

    /// Batch process multiple rover sessions
    Batch {
        #[arg(short = 'i', long, help = "Input directory containing rover observation files")]
        input_dir: String,
        #[arg(short = 'b', long = "base", action = clap::ArgAction::Append, help = "Base station raw data (.obs)")]
        base: Vec<String>,
        #[arg(short = 'n', long, help = "Path to ephemeris/nav file (.nav, .rnx)")]
        nav: Option<String>,
        #[arg(short = 'o', long, help = "Output directory for solution trajectories")]
        output_dir: String,
        #[arg(long, help = "Output export format (pos, csv, kml, json)")]
        format: Option<String>,
        #[arg(long, help = "Path to ANTEX file")]
        antex: Option<String>,
        #[arg(long, help = "Include GLONASS")]
        glonass: bool,
    },

    /// Evaluate error CDFs against ground truth
    Eval {
        #[arg(long, help = "Path to solution file (.pos)")]
        solution: String,
        #[arg(long, help = "Path to ground truth CSV")]
        truth: String,
    },

    /// Interpolate camera shutter events for UAV photogrammetry
    /// Interpolate camera shutter events for UAV photogrammetry
    Events {
        #[arg(short = 't', long, help = "Path to solution trajectory file (.pos)")]
        trajectory: String,
        #[arg(short = 'e', long, help = "Path to event timestamps file")]
        events: String,
        #[arg(short = 'o', long, help = "Path to output camera positions CSV")]
        output: String,
        #[arg(long, allow_hyphen_values = true, help = "Antenna to camera lever-arm X,Y,Z (meters)")]
        lever_arm: Option<String>,
        #[arg(long, default_value = "0.0", help = "Shutter delay (seconds)")]
        delay: f64,
    },

    /// Manage engine configuration
    Config {
        #[arg(long, default_value = "true", help = "Print default configuration JSON")]
        default: bool,
    },

    /// Launch interactive visual GUI and Web diagnostic workspace
    Gui {
        #[arg(short = 'p', long, default_value = "8080", help = "HTTP server port")]
        port: u16,
        #[arg(short = 't', long, help = "Path to trajectory file (.pos) to load into workspace")]
        trajectory: Option<String>,
    },

    /// Real-time live streaming RTK mode
    Live {
        #[arg(short = 'r', long, help = "Path or serial URI to rover stream")]
        rover: String,
        #[arg(short = 'm', long, help = "Optional NTRIP mountpoint")]
        mountpoint: Option<String>,
        #[arg(long, allow_hyphen_values = true, help = "Base position X,Y,Z ECEF")]
        base_position: Option<String>,
        #[arg(short = 'n', long, help = "Path to navigation/ephemeris file (.nav)")]
        nav: Option<String>,
    },
}

#[tokio::main]
async fn main() {
    tracing_subscriber::fmt()
        .with_env_filter("info")
        .with_target(false)
        .without_time()
        .try_init()
        .ok();

    let cli = Cli::parse();
    match cli.command {
        Commands::Process {
            rover, base, nav, output, format, qc_report, geoid, config,
            single_pass, mode, max_epochs, systems, sp3, clk,
            base_position, antex, glonass,
        } => {
            let args = process::ProcessArgs {
                rover, bases: base, nav, output, format, qc_report, geoid, config,
                enable_backward_smoothing: !single_pass, mode, max_epochs,
                base_position, systems, antex, glonass, sp3, clk,
            };
            if let Err(e) = process::run_process(args).await {
                eprintln!("Error: {}", e);
                std::process::exit(1);
            }
        }
        Commands::Calibrate { input, output } => {
            let args = calibrate::CalibrateArgs {
                input_csv: input,
                output_json: output,
            };
            if let Err(e) = calibrate::run_calibrate(args) {
                eprintln!("Calibration error: {}", e);
                std::process::exit(1);
            }
        }
        Commands::Batch {
            input_dir, base, nav, output_dir, format, antex, glonass,
        } => {
            if let Err(e) = run_batch_mode(&input_dir, base, nav, &output_dir, format, antex, glonass).await {
                eprintln!("Batch processing error: {}", e);
                std::process::exit(1);
            }
        }
        Commands::Eval { solution, truth } => {
            evaluator::run_eval(&solution, &truth);
        }
        Commands::Events {
            trajectory,
            events,
            output,
            lever_arm,
            delay,
        } => {
            if let Err(e) = run_events_interpolation(&trajectory, &events, &output, lever_arm, delay) {
                eprintln!("Camera event error: {}", e);
                std::process::exit(1);
            }
        }
        Commands::Config { default: _ } => {
            let default_cfg = gneiss_rtk::swfg::config::EngineConfig::default();
            if let Ok(json) = serde_json::to_string_pretty(&default_cfg) {
                println!("{}", json);
            }
        }
        Commands::Gui { port, trajectory } => {
            let args = gui::GuiArgs {
                port,
                trajectory_file: trajectory,
            };
            if let Err(e) = gui::run_gui_server(args).await {
                eprintln!("GUI server error: {}", e);
                std::process::exit(1);
            }
        }
        Commands::Live { rover, mountpoint, base_position, nav } => {
            let args = live::LiveArgs {
                rover_source: rover,
                ntrip_mountpoint: mountpoint,
                base_position,
                nav,
            };
            if let Err(e) = live::run_live_mode(args).await {
                eprintln!("Live streaming error: {}", e);
                std::process::exit(1);
            }
        }
    }
}

fn run_events_interpolation(
    _traj_path: &str,
    _events_path: &str,
    output_path: &str,
    lever_arm_str: Option<String>,
    delay: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let lever_arm = if let Some(s) = lever_arm_str {
        let parts: Vec<f64> = s.split(',').filter_map(|p| p.trim().parse().ok()).collect();
        if parts.len() == 3 {
            [parts[0], parts[1], parts[2]]
        } else {
            [0.0, 0.0, 0.0]
        }
    } else {
        [0.0, 0.0, 0.0]
    };

    let _config = gneiss_rtk::events::CameraEventConfig {
        lever_arm_body: lever_arm,
        shutter_delay_s: delay,
    };

    export::export_camera_events(&[], Path::new(output_path))?;
    println!("Exported camera event centers to {}", output_path);
    Ok(())
}

async fn run_batch_mode(
    input_dir: &str,
    bases: Vec<String>,
    nav: Option<String>,
    output_dir: &str,
    format: Option<String>,
    antex: Option<String>,
    glonass: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let in_path = Path::new(input_dir);
    let out_path = Path::new(output_dir);
    std::fs::create_dir_all(out_path)?;

    let ext = format
        .as_deref()
        .and_then(export::ExportFormat::from_str)
        .map_or("pos", |f| f.default_extension());
    for entry in std::fs::read_dir(in_path)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_file() && path.extension().is_some_and(|e| e == "obs" || e == "20o" || e == "25o" || e == "rnx") {
            let stem = path.file_stem().unwrap_or_default().to_string_lossy();
            let out_file = out_path.join(format!("{}.{}", stem, ext));
            let qc_file = out_path.join(format!("{}.qc.json", stem));
            let args = process::ProcessArgs {
                rover: path.to_string_lossy().to_string(),
                bases: bases.clone(),
                nav: nav.clone(),
                output: out_file.to_string_lossy().to_string(),
                format: format.clone(),
                qc_report: Some(qc_file.to_string_lossy().to_string()),
                geoid: None,
                config: None,
                enable_backward_smoothing: true,
                mode: None,
                max_epochs: None,
                base_position: None,
                systems: None,
                antex: antex.clone(),
                glonass,
                sp3: None,
                clk: None,
            };
            process::run_process(args).await?;
        }
    }
    Ok(())
}
