//! Which frequency bands do Galileo observations actually populate?
use gneiss_parsers::rinex::parse_rinex_obs;
use std::io::BufReader;

fn main() {
    let f = std::fs::File::open("datasets/multignss_2025d160/p1811600.25o").unwrap();
    let (epochs, _) = parse_rinex_obs(BufReader::new(f)).unwrap();
    println!("epochs: {}", epochs.len());
    // Scan all epochs; count band-population for Galileo satellites.
    let mut band_count = [0usize; 9];
    let mut e_sats = 0usize;
    for ep in &epochs {
        for so in &ep.satellites {
            if so.sat.constellation != gneiss_core::sat::Constellation::Galileo { continue; }
            e_sats += 1;
            for o in &so.observations {
                if o.code.obs_type == gneiss_core::obs::ObsType::CarrierPhase {
                    let b = o.code.signal.freq_band as usize;
                    if b < 9 { band_count[b] += 1; }
                }
            }
        }
    }
    println!("Galileo obs-records: {e_sats}");
    for (b, n) in band_count.iter().enumerate() {
        if *n > 0 { println!("  band {b}: {n} phases"); }
    }
}
