//! GLONASS FDMA diagnostic: k-channels, λ, and implied-DD-N sanity.
use gneiss_parsers::rinex::{parse_rinex_nav, parse_rinex_obs};
use std::io::BufReader;

fn main() {
    let nav = std::fs::File::open("datasets/multignss_2025d160/station_mixed_nav.rnx")
        .expect("nav file must exist");
    let (ephems, _) = parse_rinex_nav(BufReader::new(nav))
        .expect("valid RINEX nav file required");
    println!("total ephems: {}", ephems.len());
    let glo: Vec<_> = ephems.iter().filter_map(|e| match e {
        gneiss_core::ephemeris::Ephemeris::Glonass(g) => Some((g.sat, g.freq_num)),
        _ => None,
    }).collect();
    println!("GLONASS ephems: {}", glo.len());
    let mut seen = std::collections::BTreeSet::new();
    for (sat, k) in &glo {
        if seen.insert(sat.prn) {
            let lam1 = 299_792_458.0 / (1_602.000e6 + *k as f64 * 562.5e3);
            println!("R{:02} k={:+3} λ1={:.6} m", sat.prn, k, lam1);
        }
    }
    // Observation side: which R-sats appear, with phase on band1?
    let obs = std::fs::File::open("datasets/multignss_2025d160/p2241600.25o")
        .expect("obs file must exist");
    let (epochs, _) = parse_rinex_obs(BufReader::new(obs))
        .expect("valid RINEX obs file required");
    let mut r_phases = std::collections::BTreeMap::new();
    for ep in epochs.iter().take(50) {
        for so in &ep.satellites {
            if so.sat.constellation == gneiss_core::sat::Constellation::Glonass {
                *r_phases.entry(so.sat.prn).or_insert(0usize) +=
                    so.get_observable_phase(1).is_some() as usize;
            }
        }
    }
    println!("R-sat band-1 phase counts (first 50 epochs): {:?}", r_phases);
}
