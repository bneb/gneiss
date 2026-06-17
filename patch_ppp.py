import sys

with open("crates/gneiss-rtk/src/engine/ppp.rs", "r") as f:
    content = f.read()

target = """        let mut osb_p1 = 0.0;
        let mut osb_p2 = 0.0;
        let mut osb_cp1 = 0.0;
        let mut osb_cp2 = 0.0;
        let mut is_iono_free = false;
        let mut is_cp2_fallback = false;

        if let Some(sinex) = &engine.sinex_bias {
            let convert = |bias: f64| bias * 1e-9 * LIGHT_SPEED;
            if let Some((_, code)) = p1_obs { if let Some(b) = sinex.get_bias(sat_obs.sat, code, rover_obs.time) { osb_p1 = convert(b); } }
            if let Some((_, code)) = p2_obs { if let Some(b) = sinex.get_bias(sat_obs.sat, code, rover_obs.time) { osb_p2 = convert(b); } }
            if let Some((_, code)) = cp1_obs { if let Some(b) = sinex.get_bias(sat_obs.sat, code, rover_obs.time) { osb_cp1 = convert(b); } }
            if let Some((_, code)) = cp2_obs { 
                if let Some(b) = sinex.get_bias(sat_obs.sat, code, rover_obs.time) { 
                    osb_cp2 = convert(b); 
                    if sinex.get_exact_bias(sat_obs.sat, code, rover_obs.time).is_none() {
                        is_cp2_fallback = true;
                    }
                } 
            }
            
            if let Some(p1_val) = p1.as_mut() { *p1_val -= osb_p1; }
            if let Some(p2_val) = p2.as_mut() { *p2_val -= osb_p2; }
            if let Some(cp1_val) = cp1.as_mut() { *cp1_val -= osb_cp1 / (LIGHT_SPEED / f1); }
            if let Some(cp2_val) = cp2.as_mut() { 
                *cp2_val -= osb_cp2 / (LIGHT_SPEED / f2); 
                
                if is_cp2_fallback && sat_obs.sat.constellation == Constellation::Gps {
                    if let Some((_, code)) = cp2_obs {
                        let cstr = code.to_string();
                        if cstr == "L2L" || cstr == "L2S" || cstr == "L2X" {
                            // L2C tracks L2P(Y) by +0.25 cycles.
                            // Since we applied L2W bias, we must apply the algorithmic shift.
                            *cp2_val -= 0.25;
                            tracing::debug!("Sat {} L2C phase shift (+0.25c) applied", sat_obs.sat);
                        }
                    }
                }
            }
            tracing::debug!("Sat {:?} OSBs: p1={:.3} p2={:.3} cp1={:.3} cp2={:.3}", sat_obs.sat, osb_p1, osb_p2, osb_cp1, osb_cp2);
            
            if (!engine.sp3_epochs.is_empty() || engine.clk_data.is_some()) && !engine.config.uduc_ar {
                if let (Some(pv1), Some(pv2)) = (p1, p2) {
                    let gamma = (f1 * f1) / (f2 * f2);
                    p1 = Some((gamma * pv1 - pv2) / (gamma - 1.0));
                    is_iono_free = true;
                }
                if let (Some(l1), Some(l2)) = (cp1, cp2) {
                    let gamma = (f1 * f1) / (f2 * f2);
                    let l1_m = l1 * (LIGHT_SPEED / f1);
                    let l2_m = l2 * (LIGHT_SPEED / f2);
                    cp1 = Some((gamma * l1_m - l2_m) / (gamma - 1.0) / (LIGHT_SPEED / f1));
                }
            }
        }"""

replacement = """        let mut is_iono_free = false;

        if let Some(sinex) = &engine.sinex_bias {
            let (np1, np2, ncp1, ncp2) = apply_osb_corrections(
                sinex, sat_obs.sat, rover_obs.time, f1, f2,
                p1_obs, p2_obs, cp1_obs, cp2_obs
            );
            p1 = np1; p2 = np2; cp1 = ncp1; cp2 = ncp2;
            
            if (!engine.sp3_epochs.is_empty() || engine.clk_data.is_some()) && !engine.config.uduc_ar {
                if let (Some(pv1), Some(pv2)) = (p1, p2) {
                    let gamma = (f1 * f1) / (f2 * f2);
                    p1 = Some((gamma * pv1 - pv2) / (gamma - 1.0));
                    is_iono_free = true;
                }
                if let (Some(l1), Some(l2)) = (cp1, cp2) {
                    let gamma = (f1 * f1) / (f2 * f2);
                    let l1_m = l1 * (LIGHT_SPEED / f1);
                    let l2_m = l2 * (LIGHT_SPEED / f2);
                    cp1 = Some((gamma * l1_m - l2_m) / (gamma - 1.0) / (LIGHT_SPEED / f1));
                }
            }
        }"""

if target in content:
    content = content.replace(target, replacement)
else:
    print("TARGET NOT FOUND!")
    sys.exit(1)

with open("crates/gneiss-rtk/src/engine/ppp.rs", "w") as f:
    f.write(content)

with open("crates/gneiss-rtk/src/engine/ppp.rs", "a") as f:
    f.write("""

/// Applies OSB Phase and Code biases from SINEX files, and handles algorithmic fallbacks (e.g. +0.25 L2C shift)
pub fn apply_osb_corrections(
    sinex: &gneiss_parsers::sinex_bia::SinexBias,
    sat: gneiss_core::sat::SatelliteId,
    time: gneiss_core::time::GpsTime,
    f1: f64, f2: f64,
    p1_obs: Option<(f64, gneiss_core::obs::ObsCode)>,
    p2_obs: Option<(f64, gneiss_core::obs::ObsCode)>,
    cp1_obs: Option<(f64, gneiss_core::obs::ObsCode)>,
    cp2_obs: Option<(f64, gneiss_core::obs::ObsCode)>
) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>) {
    let mut p1 = p1_obs.map(|o| o.0);
    let mut p2 = p2_obs.map(|o| o.0);
    let mut cp1 = cp1_obs.map(|o| o.0);
    let mut cp2 = cp2_obs.map(|o| o.0);

    let convert = |bias: f64| bias * 1e-9 * LIGHT_SPEED;
    
    let mut osb_p1 = 0.0;
    let mut osb_p2 = 0.0;
    let mut osb_cp1 = 0.0;
    let mut osb_cp2 = 0.0;
    let mut is_cp2_fallback = false;

    if let Some((_, code)) = p1_obs { if let Some(b) = sinex.get_bias(sat, code, time) { osb_p1 = convert(b); } }
    if let Some((_, code)) = p2_obs { if let Some(b) = sinex.get_bias(sat, code, time) { osb_p2 = convert(b); } }
    if let Some((_, code)) = cp1_obs { if let Some(b) = sinex.get_bias(sat, code, time) { osb_cp1 = convert(b); } }
    if let Some((_, code)) = cp2_obs { 
        if let Some(b) = sinex.get_bias(sat, code, time) { 
            osb_cp2 = convert(b); 
            if sinex.get_exact_bias(sat, code, time).is_none() {
                is_cp2_fallback = true;
            }
        } 
    }
    
    if let Some(p1_val) = p1.as_mut() { *p1_val -= osb_p1; }
    if let Some(p2_val) = p2.as_mut() { *p2_val -= osb_p2; }
    if let Some(cp1_val) = cp1.as_mut() { *cp1_val -= osb_cp1 / (LIGHT_SPEED / f1); }
    if let Some(cp2_val) = cp2.as_mut() { 
        *cp2_val -= osb_cp2 / (LIGHT_SPEED / f2); 
        
        if is_cp2_fallback && sat.constellation == Constellation::Gps {
            if let Some((_, code)) = cp2_obs {
                let cstr = code.to_string();
                if cstr == "L2L" || cstr == "L2S" || cstr == "L2X" {
                    // L2C tracks L2P(Y) by +0.25 cycles.
                    // Since we applied L2W bias, we must apply the algorithmic shift.
                    *cp2_val -= 0.25;
                    tracing::debug!("Sat {} L2C phase shift (+0.25c) applied", sat);
                }
            }
        }
    }
    tracing::debug!("Sat {:?} OSBs: p1={:.3} p2={:.3} cp1={:.3} cp2={:.3}", sat, osb_p1, osb_p2, osb_cp1, osb_cp2);

    (p1, p2, cp1, cp2)
}

#[cfg(test)]
mod osb_tests {
    use super::*;
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;
    use gneiss_core::obs::ObsCode;
    use gneiss_parsers::sinex_bia::SinexBias;
    use std::io::Cursor;
    use std::str::FromStr;

    #[test]
    fn test_apply_osb_shift() {
        let content = r#"%=BIA 1.00
+BIAS/SOLUTION
*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___
 OSB  G002 G02           C1W       2021:123:00000 2021:123:86400 ns            1.0000000000    0.000000
 OSB  G002 G02           L1W       2021:123:00000 2021:123:86400 ns            2.0000000000    0.000000
 OSB  G002 G02           C2W       2021:123:00000 2021:123:86400 ns            3.0000000000    0.000000
 OSB  G002 G02           L2W       2021:123:00000 2021:123:86400 ns            4.0000000000    0.000000
-BIAS/SOLUTION
"#;
        let cursor = Cursor::new(content);
        let bias = SinexBias::parse(cursor).unwrap();
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let time = GpsTime::new(2156, 129600.0);
        let f1 = 1575.42e6;
        let f2 = 1227.60e6;
        let wl1 = LIGHT_SPEED / f1;
        let wl2 = LIGHT_SPEED / f2;

        let p1_obs = Some((10.0, ObsCode::from_str("C1C").unwrap()));
        let p2_obs = Some((20.0, ObsCode::from_str("C2L").unwrap()));
        let cp1_obs = Some((30.0, ObsCode::from_str("L1C").unwrap()));
        let cp2_obs = Some((40.0, ObsCode::from_str("L2L").unwrap()));

        let (np1, np2, ncp1, ncp2) = apply_osb_corrections(
            &bias, sat, time, f1, f2, p1_obs, p2_obs, cp1_obs, cp2_obs
        );

        // p1 -> C1C falls back to C1W (1.0 ns). 1 ns = 0.299792458 m.
        let bias_m = 1.0 * 1e-9 * LIGHT_SPEED;
        assert!((np1.unwrap() - (10.0 - bias_m)).abs() < 1e-6);

        // p2 -> C2L falls back to C2W (3.0 ns)
        let bias_m2 = 3.0 * 1e-9 * LIGHT_SPEED;
        assert!((np2.unwrap() - (20.0 - bias_m2)).abs() < 1e-6);

        // cp1 -> L1C falls back to L1W (2.0 ns)
        let bias_cp1_m = 2.0 * 1e-9 * LIGHT_SPEED;
        assert!((ncp1.unwrap() - (30.0 - bias_cp1_m / wl1)).abs() < 1e-6);

        // cp2 -> L2L falls back to L2W (4.0 ns)
        let bias_cp2_m = 4.0 * 1e-9 * LIGHT_SPEED;
        // MUST shift by +0.25 cycles because L2L tracked L2W!
        // We subtract the phase bias, and then subtract 0.25.
        assert!((ncp2.unwrap() - (40.0 - bias_cp2_m / wl2 - 0.25)).abs() < 1e-6);
    }
}
""")
