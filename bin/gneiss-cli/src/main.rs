use clap::{Parser, Subcommand};

mod evaluator;
mod process;

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
        #[arg(short = 'b', long, help = "Path to base station raw data (.obs)")]
        base: Option<String>,
        #[arg(short = 'n', long, help = "Path to ephemeris/nav file (.nav, .rnx)")]
        nav: Option<String>,
        #[arg(short = 'o', long, help = "Path to output trajectory file (.pos)")]
        output: String,
        #[arg(long, help = "Path to engine configuration file (.json)")]
        config: Option<String>,
        #[arg(long, help = "Single forward-only pass (SWFG, no ambiguity resolution). \
            Much faster, much less accurate -- typically 5-10x worse -- than the default. \
            Use only when you specifically need a quick look or true single-pass behavior.")]
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
    },

    /// Evaluate error CDFs against ground truth
    Eval {
        #[arg(long, help = "Path to solution file (.pos)")]
        solution: String,
        #[arg(long, help = "Path to ground truth CSV")]
        truth: String,
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
            rover, base, nav, output, config, single_pass, mode,
            max_epochs, systems, sp3, clk, base_position, antex,
        } => {
            if let Err(e) = process::run_process(
                rover, base, nav, output, config,
                !single_pass, mode, max_epochs,
                base_position, systems, sp3, clk,
                antex,
            ).await {
                eprintln!("Error: {}", e);
                std::process::exit(1);
            }
        }
        Commands::Eval { solution, truth } => {
            evaluator::run_eval(&solution, &truth);
        }
    }
}
