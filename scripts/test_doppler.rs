use std::io::BufRead;

fn main() {
    let rover = "/Users/kevin/projects/gneiss/tests/datasets/gsdc/2021-04-29-US-SJC-2/Pixel4/supplemental/Pixel4_GnssLog.obs";
    let file = std::fs::File::open(rover).unwrap();
    let epochs = gneiss_parsers::rinex::parse_rinex_obs(std::io::BufReader::new(file)).unwrap();
    
    let mut last_cp = std::collections::HashMap::new();
    let mut last_time = std::collections::HashMap::new();

    for obs in epochs.iter().take(200) {
        let t = obs.time.tow;
        for sat in &obs.satellites {
            if let Some(cp) = sat.get_observable_phase(1) {
                let dop = sat.observations.iter().find(|o| o.code.obs_type == gneiss_core::obs::ObsType::Doppler && o.code.signal.freq_band == 1).map(|o| o.value).unwrap_or(0.0);
                if dop != 0.0 {
                    if let Some(&prev_cp) = last_cp.get(&sat.sat) {
                        if let Some(&prev_time) = last_time.get(&sat.sat) {
                            let dt = t - prev_time;
                            if dt > 0.0 && dt < 10.0 {
                                let pred_cp = prev_cp - dop * dt;
                                let diff = (cp - pred_cp).abs();
                                if diff > 7.0 {
                                    println!("SLIP DETECTED! diff={:.2} sat={:?} dt={:.2} dop={:.2} prev_cp={:.2} cp={:.2}", diff, sat.sat, dt, dop, prev_cp, cp);
                                }
                            }
                        }
                    }
                }
                last_cp.insert(sat.sat, cp);
                last_time.insert(sat.sat, t);
            }
        }
    }
}
