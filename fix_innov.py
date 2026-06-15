import re

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'r') as f:
    content = f.read()

old_in = """pub fn compute_innovations(
    state: &mut RtkState,
    group: &[(DdObservation, DdObservation)],
    ref_rover_orig: &DdObservation,
    ref_base_orig: &DdObservation,
    env: &MeasurementEnvironment,
) -> Option<(
    Vec<f64>,
    Vec<Vec<f64>>,
    Vec<f64>,
    Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
)> {
    let mut z_vals = Vec::new();
    let mut h_rows = Vec::new();
    let mut r_vals = Vec::new();
    let mut meas_type = Vec::new();

    let geom = EkfGeometryContext::new(state, env);
    let time_tow = state.time.tow;

    let find_eph = |sat| {
        env.ephemerides
            .iter()
            .filter(|e| e.sat() == sat)
            .min_by(|a, b| {
                let da = (a.toe().tow - time_tow).abs();
                let db = (b.toe().tow - time_tow).abs();
                da.partial_cmp(&db).unwrap()
            })
    };

    let ref_eph = find_eph(ref_rover_orig.sat)?;
    let ref_state = compute_sat_state(ref_eph, ref_rover_orig, ref_base_orig, state.time, &geom, env.base_time);

    let ref_idx_l1 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 1);
    let ref_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 2);

    for (rover_sat_orig, base_sat_orig) in group {
        if let Some(sat_eph) = find_eph(rover_sat_orig.sat) {
            process_single_satellite_pair(
                state, rover_sat_orig, base_sat_orig, ref_rover_orig, ref_base_orig,
                sat_eph, &ref_state, &geom, env, ref_idx_l1, ref_idx_l2,
                &mut z_vals, &mut h_rows, &mut r_vals, &mut meas_type
            );
        }
    }

    Some((z_vals, h_rows, r_vals, meas_type))
}"""

new_in = """fn find_ephemeris<'a>(
    ephemerides: &'a [gneiss_core::ephemeris::Ephemeris],
    sat: gneiss_core::sat::SatelliteId,
    time_tow: f64,
) -> Option<&'a gneiss_core::ephemeris::Ephemeris> {
    ephemerides.iter().filter(|e| e.sat() == sat).min_by(|a, b| {
        let da = (a.toe().tow - time_tow).abs();
        let db = (b.toe().tow - time_tow).abs();
        da.partial_cmp(&db).unwrap()
    })
}

pub fn compute_innovations(
    state: &mut RtkState, group: &[(DdObservation, DdObservation)],
    ref_rover_orig: &DdObservation, ref_base_orig: &DdObservation,
    env: &MeasurementEnvironment,
) -> Option<(Vec<f64>, Vec<Vec<f64>>, Vec<f64>, Vec<(gneiss_core::sat::SatelliteId, u8, f64)>)> {
    let mut z_vals = Vec::new(); let mut h_rows = Vec::new();
    let mut r_vals = Vec::new(); let mut meas_type = Vec::new();

    let geom = EkfGeometryContext::new(state, env);
    let ref_eph = find_ephemeris(env.ephemerides, ref_rover_orig.sat, state.time.tow)?;
    let ref_state = compute_sat_state(ref_eph, ref_rover_orig, ref_base_orig, state.time, &geom, env.base_time);

    let ref_idx_l1 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 1);
    let ref_idx_l2 = state.ambiguity_keys.iter().position(|&(s, f)| s == ref_rover_orig.sat && f == 2);

    for (rover_sat_orig, base_sat_orig) in group {
        if let Some(sat_eph) = find_ephemeris(env.ephemerides, rover_sat_orig.sat, state.time.tow) {
            process_single_satellite_pair(
                state, rover_sat_orig, base_sat_orig, ref_rover_orig, ref_base_orig,
                sat_eph, &ref_state, &geom, env, ref_idx_l1, ref_idx_l2,
                &mut z_vals, &mut h_rows, &mut r_vals, &mut meas_type
            );
        }
    }
    Some((z_vals, h_rows, r_vals, meas_type))
}"""

content = content.replace(old_in, new_in)

with open('crates/gneiss-rtk/src/engine/measurement.rs', 'w') as f:
    f.write(content)
