use gneiss_core::obs::ObsType;
use gneiss_parsers::rinex::parse_rinex_obs;
use gneiss_parsers::sinex_bia::SinexBias;
use std::fs::File;
use std::io::BufReader;
use std::str::FromStr;

pub fn test_mw() {
    let obs_file = File::open("datasets/wtzr_ppp_1224/WTZR00DEU_R_20203590000_01D_30S_MO.rnx").unwrap();
    let epochs = parse_rinex_obs(BufReader::new(obs_file)).unwrap();
    let bia_file = File::open("datasets/wtzr_ppp_1224/com21374.bia").unwrap();
    let bias = SinexBias::parse(BufReader::new(bia_file)).unwrap();

    let c = 299792458.0;
    let f1 = 1575.42e6;
    let f2 = 1227.60e6;
    let lam1 = c / f1;
    let lam2 = c / f2;
    let lam_wl = c / (f1 - f2);

    for ep in epochs.iter().take(50) {
        for sat_obs in &ep.satellites {
            if sat_obs.sat.constellation != gneiss_core::sat::Constellation::Gps { continue; }
            let prn = sat_obs.sat.prn;
            
            let l1c = sat_obs.observations.iter().find(|o| o.code.to_string() == "L1C");
            let l2s = sat_obs.observations.iter().find(|o| o.code.to_string() == "L2S");
            let c1c = sat_obs.observations.iter().find(|o| o.code.to_string() == "C1C");
            let c2s = sat_obs.observations.iter().find(|o| o.code.to_string() == "C2S");
            let c2w = sat_obs.observations.iter().find(|o| o.code.to_string() == "C2W");
            let l2w = sat_obs.observations.iter().find(|o| o.code.to_string() == "L2W");

            if l1c.is_none() || c1c.is_none() { continue; }
            let l2_obs = l2w.or(l2s);
            let c2_obs = c2w.or(c2s);
            if l2_obs.is_none() || c2_obs.is_none() { continue; }

            let l1 = l1c.unwrap().value;
            let l2 = l2_obs.unwrap().value;
            let p1 = c1c.unwrap().value;
            let p2 = c2_obs.unwrap().value;

            let b_l1 = bias.get_exact_bias(sat_obs.sat, l1c.unwrap().code, ep.time).unwrap_or(0.0) * 1e-9 * c;
            let b_l2 = bias.get_exact_bias(sat_obs.sat, l2_obs.unwrap().code, ep.time).unwrap_or(0.0) * 1e-9 * c;
            let b_p1 = bias.get_exact_bias(sat_obs.sat, c1c.unwrap().code, ep.time).unwrap_or(0.0) * 1e-9 * c;
            let b_p2 = bias.get_exact_bias(sat_obs.sat, c2_obs.unwrap().code, ep.time).unwrap_or(0.0) * 1e-9 * c;

            let l1_corr = l1 - b_l1 / lam1;
            let l2_corr = l2 - b_l2 / lam2;
            let p1_corr = p1 - b_p1;
            let p2_corr = p2 - b_p2;

            let mw = (l1_corr - l2_corr) - (f1 * p1_corr + f2 * p2_corr) / (f1 + f2) / lam_wl;
            println!("Epoch {} Sat G{:02} MW = {:.4} (L2: {}, C2: {})", ep.time.tow, prn, mw, l2_obs.unwrap().code.to_string(), c2_obs.unwrap().code.to_string());
        }
    }
}
