import sys

with open("crates/gneiss-rtk/src/engine/ppp.rs", "r") as f:
    content = f.read()

target1 = """            let (np1, np2, ncp1, ncp2) = apply_osb_corrections(
                sinex, sat_obs.sat, rover_obs.time, f1, f2,
                p1_obs, p2_obs, cp1_obs, cp2_obs
            );
            p1 = np1; p2 = np2; cp1 = ncp1; cp2 = ncp2;"""

replacement1 = """            let osb_results = apply_osb_corrections(
                sinex, sat_obs.sat, rover_obs.time, f1, f2,
                p1_obs, p2_obs, cp1_obs, cp2_obs
            );
            p1 = osb_results.p1; p2 = osb_results.p2; cp1 = osb_results.cp1; cp2 = osb_results.cp2;
            let osb_p1 = osb_results.osb_p1; let osb_p2 = osb_results.osb_p2; let osb_cp1 = osb_results.osb_cp1; let osb_cp2 = osb_results.osb_cp2;"""

if target1 in content:
    content = content.replace(target1, replacement1)

# Add fallback variable definitions in else branch
target2 = """            // Legacy DCB IF Combination fallback
            if let Some(p1_val) = p1.as_mut() {"""

replacement2 = """            let osb_p1 = 0.0; let osb_p2 = 0.0; let osb_cp1 = 0.0; let osb_cp2 = 0.0;
            // Legacy DCB IF Combination fallback
            if let Some(p1_val) = p1.as_mut() {"""

if target2 in content:
    content = content.replace(target2, replacement2)

target3 = """pub fn apply_osb_corrections(
    sinex: &gneiss_parsers::sinex_bia::SinexBias,
    sat: gneiss_core::sat::SatelliteId,
    time: gneiss_core::time::GpsTime,
    f1: f64, f2: f64,
    p1_obs: Option<(f64, gneiss_core::obs::ObsCode)>,
    p2_obs: Option<(f64, gneiss_core::obs::ObsCode)>,
    cp1_obs: Option<(f64, gneiss_core::obs::ObsCode)>,
    cp2_obs: Option<(f64, gneiss_core::obs::ObsCode)>
) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>) {"""

replacement3 = """pub struct OsbCorrections {
    pub p1: Option<f64>, pub p2: Option<f64>, pub cp1: Option<f64>, pub cp2: Option<f64>,
    pub osb_p1: f64, pub osb_p2: f64, pub osb_cp1: f64, pub osb_cp2: f64,
}

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
) -> OsbCorrections {"""

if target3 in content:
    content = content.replace(target3, replacement3)

target4 = """    tracing::debug!("Sat {:?} OSBs: p1={:.3} p2={:.3} cp1={:.3} cp2={:.3}", sat, osb_p1, osb_p2, osb_cp1, osb_cp2);

    (p1, p2, cp1, cp2)"""

replacement4 = """    tracing::debug!("Sat {:?} OSBs: p1={:.3} p2={:.3} cp1={:.3} cp2={:.3}", sat, osb_p1, osb_p2, osb_cp1, osb_cp2);

    OsbCorrections { p1, p2, cp1, cp2, osb_p1, osb_p2, osb_cp1, osb_cp2 }"""

if target4 in content:
    content = content.replace(target4, replacement4)

target5 = """        let (np1, np2, ncp1, ncp2) = apply_osb_corrections(
            &bias, sat, time, f1, f2, p1_obs, p2_obs, cp1_obs, cp2_obs
        );"""

replacement5 = """        let osb_results = apply_osb_corrections(
            &bias, sat, time, f1, f2, p1_obs, p2_obs, cp1_obs, cp2_obs
        );
        let np1 = osb_results.p1; let np2 = osb_results.p2; let ncp1 = osb_results.cp1; let ncp2 = osb_results.cp2;"""

if target5 in content:
    content = content.replace(target5, replacement5)

target6 = """        let mut is_iono_free = false;

        if let Some(sinex) = &engine.sinex_bias {"""

replacement6 = """        let mut is_iono_free = false;
        let mut osb_p1 = 0.0;
        let mut osb_p2 = 0.0;
        let mut osb_cp1 = 0.0;
        let mut osb_cp2 = 0.0;

        if let Some(sinex) = &engine.sinex_bias {"""

if target6 in content:
    content = content.replace(target6, replacement6)

target7 = """            p1 = osb_results.p1; p2 = osb_results.p2; cp1 = osb_results.cp1; cp2 = osb_results.cp2;
            let osb_p1 = osb_results.osb_p1; let osb_p2 = osb_results.osb_p2; let osb_cp1 = osb_results.osb_cp1; let osb_cp2 = osb_results.osb_cp2;"""

replacement7 = """            p1 = osb_results.p1; p2 = osb_results.p2; cp1 = osb_results.cp1; cp2 = osb_results.cp2;
            osb_p1 = osb_results.osb_p1; osb_p2 = osb_results.osb_p2; osb_cp1 = osb_results.osb_cp1; osb_cp2 = osb_results.osb_cp2;"""

if target7 in content:
    content = content.replace(target7, replacement7)

target8 = """            let osb_p1 = 0.0; let osb_p2 = 0.0; let osb_cp1 = 0.0; let osb_cp2 = 0.0;
            // Legacy DCB IF Combination fallback"""

replacement8 = """            // Legacy DCB IF Combination fallback"""

if target8 in content:
    content = content.replace(target8, replacement8)

with open("crates/gneiss-rtk/src/engine/ppp.rs", "w") as f:
    f.write(content)

