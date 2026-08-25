// Ephemeral validation example: prove fetched precise products flow through
// gneiss-parsers end to end (SP3 -> PreciseOrbit interpolation, CLK lookup).
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use gneiss_parsers::precise_orbit::PreciseOrbit;
use gneiss_parsers::rinex_clk::RinexClock;
use gneiss_parsers::sp3::parse_sp3;

fn main() {
    let repo = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(2)
        .expect("crate inside repo")
        .to_path_buf();
    let sp3_path = repo.join("datasets/precise/2025/160/GFZ0MGXRAP_20250609.sp3");
    let clk_path = repo.join("datasets/precise/2025/160/GFZ0MGXRAP_20250609.clk");

    let sp3_text = std::fs::read_to_string(&sp3_path).expect("sp3 present");
    let epochs = parse_sp3(sp3_text.as_bytes()).expect("sp3 parses");
    let first: Vec<&String> = epochs[0].records.keys().collect();
    let gal = first.iter().filter(|s| s.starts_with('E')).count();
    let gps = first.iter().filter(|s| s.starts_with('G')).count();
    println!(
        "SP3: {} epochs, {} sats in first epoch (GPS={}, GAL={})",
        epochs.len(),
        first.len(),
        gps,
        gal
    );
    assert!(epochs.len() >= 280, "expected ~289 five-minute epochs");
    assert!(gal >= 20 && gps >= 28, "multi-GNSS coverage missing");

    let orbit = PreciseOrbit::new(epochs.clone());
    // 2025-06-09 is Monday of GPS week 2370 => local noon = TOW 129600.
    let t = GpsTime {
        week: 2370,
        tow: 129_600.0,
    };
    let pos = orbit.position_at("E02", t).expect("E02 interpolates");
    println!(
        "E02 @12h GPS: x={:.1} y={:.1} z={:.1} m, |r|={:.1} m",
        pos.0.x,
        pos.0.y,
        pos.0.z,
        pos.0.norm()
    );


    let clk_text = std::fs::read_to_string(&clk_path).expect("clk present");
    let clock = RinexClock::parse(&clk_text);
    let n_sats = clock.satellites.len();
    let e02 = SatelliteId {
        constellation: Constellation::Galileo,
        prn: 2,
    };
    let bias = clock.get_clock_bias(e02, t).expect("E02 clock record");
    println!("CLK: {} satellites tracked; E02 bias @12h = {:.4e} s", n_sats, bias);
    assert!(n_sats >= 80, "expected full multi-GNSS clock coverage");
    assert!(bias.abs() < 1e-3, "clock bias should be sub-millisecond");
    println!("PRECISE PRODUCT PIPELINE OK");
}
