use std::fs::File;
use std::io::BufReader;

fn main() {
    let file = File::open("datasets/gsdc/Pixel4_GnssLog.20o").unwrap();
    let reader = BufReader::new(file);
    let epochs = gneiss_parsers::rinex::parse_rinex_obs(reader).unwrap();
    for epoch in epochs.iter().take(2) {
        println!("Epoch {}", epoch.time);
        for sat in &epoch.satellites {
            let cp1 = sat.get_observable_phase(sat.primary_band());
            println!("Sat {:?}, cp1={:?}", sat.sat, cp1);
        }
    }
}
