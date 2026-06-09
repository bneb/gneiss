import sys

with open("crates/gneiss-rtk/src/spp.rs", "r") as f:
    content = f.read()

# 1. Update struct
content = content.replace(
"""    pub cdt_bds: f64,
}

impl SppState {
    pub fn new(position: Coordinate, cdt: f64, cdt_gal: f64, cdt_bds: f64) -> Self {
        Self { position, cdt, cdt_gal, cdt_bds }
    }""",
"""    pub cdt_bds: f64,
    /// Receiver clock bias in meters for GLONASS.
    pub cdt_glo: f64,
}

impl SppState {
    pub fn new(position: Coordinate, cdt: f64, cdt_gal: f64, cdt_bds: f64, cdt_glo: f64) -> Self {
        Self { position, cdt, cdt_gal, cdt_bds, cdt_glo }
    }"""
)

# 2. Update seed_initial_state
content = content.replace(
"""    // Estimate initial clock bias from the first satellite
    let m0 = &measurements[0];
    let dx0 = seed_x - m0.sat_coord.vector.x;
    let dy0 = seed_y - m0.sat_coord.vector.y;
    let dz0 = seed_z - m0.sat_coord.vector.z;
    let geom_r0 = f64::sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
    let seed_cdt = m0.pseudorange - geom_r0;

    SppState::new(
        Coordinate::new(Vector3::new(seed_x, seed_y, seed_z), Datum::WGS84, Frame::ECEF, measurements[0].time),
        seed_cdt,
        seed_cdt,
        seed_cdt
    )""",
"""    let mut cdt_gps = None;
    let mut cdt_gal = None;
    let mut cdt_bds = None;
    let mut cdt_glo = None;

    for m in measurements {
        let dx = seed_x - m.sat_coord.vector.x;
        let dy = seed_y - m.sat_coord.vector.y;
        let dz = seed_z - m.sat_coord.vector.z;
        let geom_r = f64::sqrt(dx * dx + dy * dy + dz * dz);
        let cdt = m.pseudorange - geom_r;
        
        match m.constellation {
            gneiss_core::sat::Constellation::Gps => if cdt_gps.is_none() { cdt_gps = Some(cdt); },
            gneiss_core::sat::Constellation::Galileo => if cdt_gal.is_none() { cdt_gal = Some(cdt); },
            gneiss_core::sat::Constellation::Beidou => if cdt_bds.is_none() { cdt_bds = Some(cdt); },
            gneiss_core::sat::Constellation::Glonass => if cdt_glo.is_none() { cdt_glo = Some(cdt); },
            _ => {},
        }
    }

    let default_cdt = cdt_gps.or(cdt_gal).or(cdt_bds).or(cdt_glo).unwrap_or(0.0);

    SppState::new(
        Coordinate::new(Vector3::new(seed_x, seed_y, seed_z), Datum::WGS84, Frame::ECEF, measurements[0].time),
        cdt_gps.unwrap_or(default_cdt),
        cdt_gal.unwrap_or(default_cdt),
        cdt_bds.unwrap_or(default_cdt),
        cdt_glo.unwrap_or(default_cdt)
    )"""
)

# 3. Update compute_spp RAIM initialization
content = content.replace(
"""                let mut clean_state = SppState::new(
                    Coordinate::new(Vector3::new(seed_x, seed_y, seed_z), Datum::WGS84, Frame::ECEF, measurements[0].time),
                    seed_cdt,
                    seed_cdt_gal,
                    seed_cdt_bds
                );""",
"""                let mut clean_state = SppState::new(
                    Coordinate::new(Vector3::new(seed_x, seed_y, seed_z), Datum::WGS84, Frame::ECEF, measurements[0].time),
                    seed_cdt,
                    seed_cdt_gal,
                    seed_cdt_bds,
                    seed_cdt_glo
                );"""
)

# Also need to add seed_cdt_glo extraction before the RAIM loop:
content = content.replace(
"""    let seed_cdt_gal = state.cdt_gal;
    let seed_cdt_bds = state.cdt_bds;""",
"""    let seed_cdt_gal = state.cdt_gal;
    let seed_cdt_bds = state.cdt_bds;
    let seed_cdt_glo = state.cdt_glo;"""
)

# 4. Update spp_wnlls_step end
content = content.replace(
"""    let next_cdt = current_state.cdt + gps_col.map(|c| dx_vec[c]).unwrap_or(0.0);
    let next_cdt_gal = current_state.cdt_gal + gal_col.map(|c| dx_vec[c]).unwrap_or(0.0);
    let next_cdt_bds = current_state.cdt_bds + bds_col.map(|c| dx_vec[c]).unwrap_or(0.0);

    Ok(SppState::new(
        Coordinate::new(
            Vector3::new(
                current_state.position.vector.x + dx_vec[0],
                current_state.position.vector.y + dx_vec[1],
                current_state.position.vector.z + dx_vec[2],
            ),
            Datum::WGS84,
            Frame::ECEF,
            measurements[0].time
        ),
        next_cdt,
        next_cdt_gal,
        next_cdt_bds
    ))""",
"""    let next_cdt = current_state.cdt + gps_col.map(|c| dx_vec[c]).unwrap_or(0.0);
    let next_cdt_gal = current_state.cdt_gal + gal_col.map(|c| dx_vec[c]).unwrap_or(0.0);
    let next_cdt_bds = current_state.cdt_bds + bds_col.map(|c| dx_vec[c]).unwrap_or(0.0);
    let next_cdt_glo = current_state.cdt_glo + glo_col.map(|c| dx_vec[c]).unwrap_or(0.0);

    Ok(SppState::new(
        Coordinate::new(
            Vector3::new(
                current_state.position.vector.x + dx_vec[0],
                current_state.position.vector.y + dx_vec[1],
                current_state.position.vector.z + dx_vec[2],
            ),
            Datum::WGS84,
            Frame::ECEF,
            measurements[0].time
        ),
        next_cdt,
        next_cdt_gal,
        next_cdt_bds,
        next_cdt_glo
    ))"""
)

# 5. Update test case in spp.rs
content = content.replace(
"""        let mut state = SppState::new(
            Coordinate::new(Vector3::new(0.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time),
            0.0, 0.0, 0.0
        );""",
"""        let mut state = SppState::new(
            Coordinate::new(Vector3::new(0.0, 0.0, 0.0), Datum::WGS84, Frame::ECEF, time),
            0.0, 0.0, 0.0, 0.0
        );"""
)

with open("crates/gneiss-rtk/src/spp.rs", "w") as f:
    f.write(content)

