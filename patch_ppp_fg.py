import re

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "r") as f:
    content = f.read()

# build_h_row signature
content = re.sub(r'fn build_h_row\(los: &Vector3<f64>, map_wet: f64, amb_idx: Option<usize>, size: usize\) -> DVector<f64> \{',
                 r'fn build_h_row(los: &Vector3<f64>, map_wet: f64, amb_idx: Option<usize>, size: usize, constel: gneiss_core::ephemeris::Constellation) -> DVector<f64> {', content)

# build_h_row body
body1_old = """    h[15] = 1.0;
    if size > 17 { h[17] = map_wet; }
    if let Some(idx) = amb_idx { h[idx] = 1.0; }"""
body1_new = """    h[15] = 1.0;
    if size > 18 {
        match constel {
            gneiss_core::ephemeris::Constellation::Glonass => h[16] = 1.0,
            gneiss_core::ephemeris::Constellation::Galileo => h[17] = 1.0,
            gneiss_core::ephemeris::Constellation::BeiDou => h[18] = 1.0,
            _ => {}
        }
    }
    if size > 20 { h[20] = map_wet; }
    if let Some(idx) = amb_idx { h[idx] = 1.0; }"""
content = content.replace(body1_old, body1_new)

# build_measurements
meas_old = """        let ztd = if x_i.len() > 17 { x_i[17] } else { 0.0 };

        for sat in sats {
            let dist = (sat.sat_pos_rot - rcv_pos).norm();
            let los = (sat.sat_pos_rot - rcv_pos) / dist;
            let expected_base = dist + x_i[15] - sat.dt_sat_m + sat.tropo_dry + ztd * sat.map_wet;"""
meas_new = """        let ztd = if x_i.len() > 20 { x_i[20] } else { 0.0 };

        for sat in sats {
            let dist = (sat.sat_pos_rot - rcv_pos).norm();
            let los = (sat.sat_pos_rot - rcv_pos) / dist;
            let isb = if x_i.len() > 18 {
                match sat.sat_obs.sat.constellation {
                    gneiss_core::ephemeris::Constellation::Glonass => x_i[16],
                    gneiss_core::ephemeris::Constellation::Galileo => x_i[17],
                    gneiss_core::ephemeris::Constellation::BeiDou => x_i[18],
                    _ => 0.0,
                }
            } else { 0.0 };
            let expected_base = dist + x_i[15] + isb - sat.dt_sat_m + sat.tropo_dry + ztd * sat.map_wet;"""
content = content.replace(meas_old, meas_new)

# h_row call PR
content = content.replace("h_row: build_h_row(&los, sat.map_wet, None, x_i.len()),", 
                          "h_row: build_h_row(&los, sat.map_wet, None, x_i.len(), sat.sat_obs.sat.constellation),")
# h_row call CP
content = content.replace("h_row: build_h_row(&los, sat.map_wet, Some(CORE_STATE_SIZE + amb_idx), x_i.len()),",
                          "h_row: build_h_row(&los, sat.map_wet, Some(CORE_STATE_SIZE + amb_idx), x_i.len(), sat.sat_obs.sat.constellation),")

# build_h_row_doppler signature
content = content.replace("fn build_h_row_doppler(los: &Vector3<f64>, size: usize) -> DVector<f64> {",
                          "fn build_h_row_doppler(los: &Vector3<f64>, size: usize) -> DVector<f64> {") # no ISB in doppler?
# Wait! Does Doppler have ISB? Clock drift is just 1. It shouldn't differ between constellations because all receiver channels run off the same oscillator! So Doppler is fine.

# fix build_h_row_doppler body for size > 19
dop_old = """    if size > 16 {
        h[3] = -los.x; h[4] = -los.y; h[5] = -los.z;
        h[16] = 1.0;
    }"""
dop_new = """    if size > 19 {
        h[3] = -los.x; h[4] = -los.y; h[5] = -los.z;
        h[19] = 1.0;
    }"""
content = content.replace(dop_old, dop_new)

with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "w") as f:
    f.write(content)

