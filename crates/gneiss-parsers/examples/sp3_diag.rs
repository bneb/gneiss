//! Diagnostic: broadcast vs SP3 satellite position agreement.
use gneiss_parsers::{precise_orbit::PreciseOrbit, rinex::parse_rinex_nav, sp3::parse_sp3};

fn main() {
    let nav = std::fs::File::open("datasets/multignss_2025d160/station_mixed_nav.rnx")
        .expect("nav file must exist");
    let (ephems, _) = parse_rinex_nav(std::io::BufReader::new(nav))
        .expect("valid RINEX nav file required");

    let sp3f = std::fs::File::open("datasets/precise/2025/160/GFZ0MGXRAP_20250609.sp3")
        .expect("SP3 file must exist");
    let epochs = parse_sp3(std::io::BufReader::new(sp3f))
        .expect("valid SP3 file required");
    // Evaluate at the first GPS TOE so broadcast fits are inside their
    // validity window; compare SP3 interpolation at that same instant.
    let t_mid = ephems.iter()
        .find_map(|e| match e {
            gneiss_core::ephemeris::Ephemeris::Gps(g) => Some(g.toe),
            _ => None,
        })
        .expect("GPS ephemeris required");
    println!("Evaluating at first GPS TOE: week={} tow={:.0}", t_mid.week, t_mid.tow);
    let precise = PreciseOrbit::new(epochs);

    let mut shown = 0;
    for e in ephems.iter() {
        if shown >= 6 { break; }
        let sat = e.sat();
        if sat.constellation != gneiss_core::sat::Constellation::Gps { continue; }
        // Broadcast position via the per-variant position() method
        // Select the GPS ephemeris with TOE nearest to the evaluation time
        // (broadcast fits are only valid ~±2 h around TOE).
        let mut best: Option<(&gneiss_core::ephemeris::Ephemeris, f64)> = None;
        for cand in ephems.iter() {
            if cand.sat().prn != sat.prn { continue; }
            if !matches!(cand, gneiss_core::ephemeris::Ephemeris::Gps(_)) { continue; }
            let toe = match cand {
                gneiss_core::ephemeris::Ephemeris::Gps(g) => g.toe,
                _ => continue,
            };
            let dt = if toe.week == t_mid.week { (toe.tow - t_mid.tow).abs() } else { f64::INFINITY };
            if best.is_none_or(|(_, d)| dt < d) {
                best = Some((cand, dt));
            }
        }
        let (cand_best, dt_best) = match best { Some(x) => x, None => continue };
        let (bc_pos, _vel, _clk_bias, _clk_drift) = match cand_best {
            gneiss_core::ephemeris::Ephemeris::Gps(g) => g.position(t_mid),
            _ => unreachable!(),
        };
        let _ = dt_best;
        let name = format!("{}{:02}", "G", sat.prn);
        if let Some((sp3_pos, clk)) = precise.position_at(&name, t_mid) {
            let d = (bc_pos - sp3_pos).norm();
            println!("{}: |bcast-sp3| = {:9.3} m  |dt_toe|={:6.0}s  sp3_clk={:+9.1} ns", name, d, dt_best, clk * 1e9);
            shown += 1;
        }
    }
}
