//! Observation parsing and SPP seeding for the SWFG engine.

use nalgebra::Vector3;

use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_core::sat::SatelliteId;
use gneiss_core::time::GpsTime;

use crate::swfg::pipeline::RawObservation;

pub fn select_best_ephemeris(
    ephemerides: &[Ephemeris],
    sat: SatelliteId,
    time: GpsTime,
) -> Option<&Ephemeris> {
    ephemerides
        .iter()
        .filter(|e| e.sat() == sat)
        .min_by(|a, b| {
            (a.toe().tow - time.tow)
                .abs()
                .partial_cmp(&(b.toe().tow - time.tow).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
}

pub fn extract_raw_observations(
    rover: &EpochObs,
    ephemerides: &[Ephemeris],
    approx_rx_pos: Option<Vector3<f64>>,
) -> Result<Vec<RawObservation>, String> {
    use gneiss_core::signal::satellite_frequencies;

    let mut raw_obs = Vec::new();
    let rx_pos = approx_rx_pos.unwrap_or_else(|| Vector3::new(
        gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0,
    ));

    for sat_obs in &rover.satellites {
        let eph = match select_best_ephemeris(ephemerides, sat_obs.sat, rover.time) {
            Some(e) => e,
            None => continue,
        };

        let (f1, f2) = satellite_frequencies(sat_obs.sat, eph.freq_num());

        // Read L1 pseudorange. BeiDou uses RINEX band 2 for B1I.
        let p1_band = match sat_obs.sat.constellation {
            gneiss_core::sat::Constellation::Beidou => 2,
            _ => 1,
        };
        let pr_m = match sat_obs.get_observable(p1_band) {
            Some(pr) if pr > 1e6 => pr,
            _ => continue,
        };

        // Two-pass light-time correction (follows spp.rs:206-236).
        // Pass 1: approximate transmit time from pseudorange.
        let t_tx = rover.time - (pr_m / gneiss_core::constants::SPEED_OF_LIGHT_M_S);

        // Get rough satellite clock to refine transmit time.
        let (_, _, sat_clk_err_rough, _) = eph.position(t_tx);

        // Pass 2: refine transmit time with satellite clock correction.
        let t_tx_true = t_tx - sat_clk_err_rough;
        let (mut sat_pos, sat_vel, clk_err, _clk_drift) = eph.position(t_tx_true);
        let sat_clock_m = clk_err * gneiss_core::constants::SPEED_OF_LIGHT_M_S;

        // Apply Earth rotation (Sagnac) correction to satellite position
        let flight_time = pr_m / gneiss_core::constants::SPEED_OF_LIGHT_M_S;
        let omega_tau = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * flight_time;
        let cos_wt = omega_tau.cos();
        let sin_wt = omega_tau.sin();
        sat_pos = Vector3::new(
            sat_pos.x * cos_wt + sat_pos.y * sin_wt,
            -sat_pos.x * sin_wt + sat_pos.y * cos_wt,
            sat_pos.z,
        );

        let pr_l2 = sat_obs.get_observable(2);
        let (cp_l1, cp_l1_lli) = match sat_obs.get_observable_phase_lli(1) {
            Some((cp, lli)) => (Some(cp), lli),
            None => (None, None),
        };
        let cp_l2 = sat_obs.get_observable_phase(2);
        let snr = sat_obs.get_snr(1).map_or(40.0, |s| s as f64);

        let (el_rad, az_rad) = if sat_pos.norm() < 1e6 {
            (std::f64::consts::FRAC_PI_2, 0.0)
        } else {
            let ref_llh = gneiss_core::coords::ecef_to_llh(rx_pos);
            let dx = sat_pos.x - rx_pos.x;
            let dy = sat_pos.y - rx_pos.y;
            let dz = sat_pos.z - rx_pos.z;
            let sin_lat = ref_llh.x.sin();
            let cos_lat = ref_llh.x.cos();
            let sin_lon = ref_llh.y.sin();
            let cos_lon = ref_llh.y.cos();

            let e = -sin_lon * dx + cos_lon * dy;
            let n = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz;
            let u = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz;
            let dist = (e * e + n * n + u * u).sqrt().max(1.0);
            let el = (u / dist).clamp(-1.0, 1.0).asin();
            let az = e.atan2(n);
            (el, az)
        };

        // Filter satellites below elevation mask
        if sat_pos.norm() >= 1e6 && el_rad < 10.0_f64.to_radians() {
            continue;
        }

        let sin_el = el_rad.sin().max(0.17);
        let snr_penalty = if snr < 40.0 {
            10.0_f64.powf((40.0 - snr) / 10.0)
        } else {
            1.0
        };
        let variance_m2 = (1.0 / (sin_el * sin_el) * snr_penalty).clamp(0.5, 100.0);
        let cp_variance_m2 = (9e-6 / (sin_el * sin_el) * snr_penalty).clamp(9e-6, 1e-3);

        raw_obs.push(RawObservation {
            satellite: sat_obs.sat.prn as u16,
            constellation_id: sat_obs.sat.constellation as u8,
            pr_l1: pr_m,
            pr_l2,
            cp_l1,
            cp_l1_lli,
            cp_l2,
            doppler: sat_obs.get_doppler(1).unwrap_or(0.0),
            snr_dbhz: snr,
            sat_pos_ecef: sat_pos,
            sat_vel_ecef: sat_vel,
            sat_clock_m,
            f1,
            f2,
            freq_num: eph.freq_num(),
            elevation_rad: el_rad,
            azimuth_rad: az_rad,
            tropo_dry_m: 0.0,
            tropo_map_wet: 1.0 / sin_el,
            iono_l1_m: 0.0,
            variance_m2,
            cp_variance_m2,
        });
    }

    Ok(raw_obs)
}

/// Extract raw observations using a generic ephemeris source (e.g. PreciseSrc with SP3/CLK).
fn lookup_satellite_bias(
    bias: Option<&gneiss_parsers::sinex_bia::SinexBias>,
    sat: SatelliteId,
    code: gneiss_core::obs::ObsCode,
    time: GpsTime,
) -> f64 {
    bias.and_then(|b| b.get_bias(sat, code, time))
        .unwrap_or(0.0) * gneiss_core::constants::SPEED_OF_LIGHT_M_S * 1e-9
}

fn compute_sat_pco_body(sat: &SatelliteId) -> Vector3<f64> {
    match sat.constellation {
        gneiss_core::sat::Constellation::Gps => match sat.prn {
            1 | 3 | 6 | 8 | 9 | 10 | 24 | 25 | 26 | 27 | 30 | 32 => Vector3::new(0.394, 0.0, 1.500),
            4 | 5 | 14 | 18 | 23 => Vector3::new(0.0, 0.0, 1.500),
            _ => Vector3::new(0.0, 0.0, 0.952),
        },
        gneiss_core::sat::Constellation::Galileo => Vector3::new(0.0, 0.0, 1.000),
        _ => Vector3::zeros(),
    }
}

fn sat_antex_freq_codes(
    c: gneiss_core::sat::Constellation,
    ant: &gneiss_parsers::antex::AntennaPcv,
) -> (&'static str, &'static str) {
    match c {
        gneiss_core::sat::Constellation::Galileo => {
            let f2 = if ant.frequencies.contains_key("E05") { "E05" }
            else if ant.frequencies.contains_key("E07") { "E07" }
            else { "E08" };
            ("E01", f2)
        }
        gneiss_core::sat::Constellation::Glonass => ("R01", "R02"),
        gneiss_core::sat::Constellation::Beidou => {
            let f1 = if ant.frequencies.contains_key("C02") { "C02" } else { "C01" };
            let f2 = if ant.frequencies.contains_key("C06") { "C06" }
            else if ant.frequencies.contains_key("C07") { "C07" }
            else { "C02" };
            (f1, f2)
        }
        _ => ("G01", "G02"),
    }
}

fn compute_sat_pco_from_antex(
    antex: &gneiss_parsers::antex::AntexDatabase,
    sat: &SatelliteId,
    rover_time: GpsTime,
    f1: f64,
    f2: f64,
) -> Vector3<f64> {
    let sys_char = match sat.constellation {
        gneiss_core::sat::Constellation::Gps => 'G',
        gneiss_core::sat::Constellation::Glonass => 'R',
        gneiss_core::sat::Constellation::Galileo => 'E',
        gneiss_core::sat::Constellation::Beidou => 'C',
        _ => return compute_sat_pco_body(sat),
    };
    let ant = match antex.find_satellite_gps(&format!("{}{:02}", sys_char, sat.prn), rover_time) {
        Some(a) => a,
        None => return compute_sat_pco_body(sat),
    };
    let (c1, c2) = sat_antex_freq_codes(sat.constellation, ant);
    let p1 = ant.frequencies.get(c1).map(|f| f.pco);
    let p2 = ant.frequencies.get(c2).map(|f| f.pco);
    let pco_mm = match (p1, p2) {
        (Some(p1), Some(p2)) if f1 > 0.0 && f2 > 0.0 && (f1 - f2).abs() > 1.0 => {
            let gamma = (f1 / f2).powi(2);
            (p1 * gamma - p2) / (gamma - 1.0)
        }
        (Some(p1), _) => p1,
        (_, Some(p2)) => p2,
        _ => return compute_sat_pco_body(sat),
    };
    pco_mm / 1000.0
}

fn compute_sat_clock_delay(
    src: &dyn crate::estimators::rtk_iekf::satpos::EphSource,
    sat: &SatelliteId,
    rx_eff: Vector3<f64>,
    sat_pos: Vector3<f64>,
    time: GpsTime,
) -> f64 {
    let tau = (rx_eff - sat_pos).norm() / gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    let t_tx = time - tau;
    let sat_clock_s = src.position_at(sat, t_tx).map_or(0.0, |(_, c)| c);
    let rel_corr = if src.com_referenced() {
        let t_plus = t_tx + 0.5;
        let t_minus = t_tx - 0.5;
        let (sat_p_plus, _) = src.position_at(sat, t_plus).unwrap_or((sat_pos, 0.0));
        let (sat_p_minus, _) = src.position_at(sat, t_minus).unwrap_or((sat_pos, 0.0));
        let sat_vel = sat_p_plus - sat_p_minus;
        gneiss_geodesy::relativity::periodic_relativistic_range_correction(&sat_pos, &sat_vel)
    } else {
        // Broadcast ephemeris clock fit already includes F * e * sqrt(a) * sin(E)
        0.0
    };
    let shapiro = gneiss_geodesy::relativity::gravitational_shapiro_delay(&sat_pos, &rx_eff);
    sat_clock_s * gneiss_core::constants::SPEED_OF_LIGHT_M_S + rel_corr - shapiro
}

fn compute_el_az(rx_pos: Vector3<f64>, sat_pos: Vector3<f64>) -> (f64, f64) {
    if sat_pos.norm() < 1e6 {
        return (std::f64::consts::FRAC_PI_2, 0.0);
    }
    let ref_llh = gneiss_core::coords::ecef_to_llh(rx_pos);
    let (dx, dy, dz) = (sat_pos.x - rx_pos.x, sat_pos.y - rx_pos.y, sat_pos.z - rx_pos.z);
    let (sin_lat, cos_lat) = (ref_llh.x.sin(), ref_llh.x.cos());
    let (sin_lon, cos_lon) = (ref_llh.y.sin(), ref_llh.y.cos());
    let e = -sin_lon * dx + cos_lon * dy;
    let n = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz;
    let u = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz;
    let dist = (e * e + n * n + u * u).sqrt().max(1.0);
    let el = (u / dist).clamp(-1.0, 1.0).asin();
    let az = e.atan2(n);
    (el, az)
}

fn compute_tropo_delay(rx_pos: Vector3<f64>, sin_el: f64) -> (f64, f64) {
    let ref_llh = gneiss_core::coords::ecef_to_llh(rx_pos);
    let h_m = ref_llh.z.clamp(-1000.0, 10000.0);
    let p_hpa = 1013.25 * (1.0 - 2.2557e-5 * h_m).powi(5);
    let zdry_m = 0.002277 * p_hpa;
    let m_dry = 1.001 / (0.002001 + sin_el * sin_el).sqrt();
    let m_wet = 1.001 / (0.002001 + sin_el * sin_el).sqrt();
    (zdry_m * m_dry, m_wet)
}

fn select_sat_bands(sat_obs: &gneiss_core::obs::SatObs) -> (u8, u8, f64, f64) {
    let p1_band = match sat_obs.sat.constellation {
        gneiss_core::sat::Constellation::Beidou => 2,
        _ => 1,
    };
    let p2_band = match sat_obs.sat.constellation {
        gneiss_core::sat::Constellation::Galileo => {
            if sat_obs.observations.iter().any(|o| o.code.signal.freq_band == 7) { 7 }
            else if sat_obs.observations.iter().any(|o| o.code.signal.freq_band == 5) { 5 }
            else { 2 }
        }
        gneiss_core::sat::Constellation::Beidou => {
            if sat_obs.observations.iter().any(|o| o.code.signal.freq_band == 7) { 7 }
            else if sat_obs.observations.iter().any(|o| o.code.signal.freq_band == 6) { 6 }
            else { 2 }
        }
        _ => {
            if sat_obs.observations.iter().any(|o| o.code.signal.freq_band == 2) { 2 }
            else if sat_obs.observations.iter().any(|o| o.code.signal.freq_band == 5) { 5 }
            else { 2 }
        }
    };
    let f1 = gneiss_core::signal::get_frequency(sat_obs.sat, p1_band, 0);
    let f2 = gneiss_core::signal::get_frequency(sat_obs.sat, p2_band, 0);
    (p1_band, p2_band, f1, f2)
}

fn extract_code_obs(
    sat_obs: &gneiss_core::obs::SatObs,
    p1_b: u8,
    p2_b: u8,
    bias: Option<&gneiss_parsers::sinex_bia::SinexBias>,
    t: GpsTime,
) -> Option<(f64, Option<f64>)> {
    let obs_p1 = sat_obs.observations.iter().find(|o| {
        o.code.obs_type == gneiss_core::obs::ObsType::Pseudorange
            && o.code.signal.freq_band == p1_b
            && o.value > 1e6
    })?;
    let pr_m = obs_p1.value - lookup_satellite_bias(bias, sat_obs.sat, obs_p1.code, t);
    let obs_p2 = sat_obs.observations.iter().find(|o| {
        o.code.obs_type == gneiss_core::obs::ObsType::Pseudorange && o.code.signal.freq_band == p2_b
    });
    let pr_l2 = obs_p2.map(|o| o.value - lookup_satellite_bias(bias, sat_obs.sat, o.code, t));
    Some((pr_m, pr_l2))
}

fn extract_phase_obs(
    sat_obs: &gneiss_core::obs::SatObs,
    p1_b: u8,
    p2_b: u8,
    f1: f64,
    f2: f64,
    bias: Option<&gneiss_parsers::sinex_bia::SinexBias>,
    t: GpsTime,
) -> (Option<f64>, Option<u8>, Option<f64>) {
    let wl_1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / f1;
    let wl_2 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / f2;
    let obs_l1 = sat_obs.observations.iter().find(|o| {
        o.code.obs_type == gneiss_core::obs::ObsType::CarrierPhase && o.code.signal.freq_band == p1_b
    });
    let (cp_l1, cp_l1_lli) = match obs_l1 {
        Some(o) => (Some(o.value - lookup_satellite_bias(bias, sat_obs.sat, o.code, t) / wl_1), o.lli),
        None => (None, None),
    };
    let obs_l2 = sat_obs.observations.iter().find(|o| {
        o.code.obs_type == gneiss_core::obs::ObsType::CarrierPhase && o.code.signal.freq_band == p2_b
    });
    let cp_l2 = obs_l2.map(|o| o.value - lookup_satellite_bias(bias, sat_obs.sat, o.code, t) / wl_2);
    (cp_l1, cp_l1_lli, cp_l2)
}

fn extract_single_sat_with_source(
    sat_obs: &gneiss_core::obs::SatObs,
    src: &dyn crate::estimators::rtk_iekf::satpos::EphSource,
    rx_pos: Vector3<f64>,
    rover_time: GpsTime,
    sinex_bias: Option<&gneiss_parsers::sinex_bia::SinexBias>,
    antex: Option<&gneiss_parsers::antex::AntexDatabase>,
) -> Option<RawObservation> {
    let is_supp = matches!(sat_obs.sat.constellation, gneiss_core::sat::Constellation::Gps | gneiss_core::sat::Constellation::Galileo);
    if !is_supp { return None; }
    let (p1_b, p2_b, f1, f2) = select_sat_bands(sat_obs);
    let (pr_l1, pr_l2) = extract_code_obs(sat_obs, p1_b, p2_b, sinex_bias, rover_time)?;
    let (cp_l1, cp_l1_lli, cp_l2) = extract_phase_obs(sat_obs, p1_b, p2_b, f1, f2, sinex_bias, rover_time);

    let tide_m = if rx_pos.norm() > 1e6 {
        gneiss_geodesy::tides::solid_earth_tide(rx_pos, rover_time.tow, rover_time.week)
    } else { Vector3::zeros() };
    let rx_eff = rx_pos + tide_m;
    let pco_body = match antex {
        Some(db) => compute_sat_pco_from_antex(db, &sat_obs.sat, rover_time, f1, f2),
        None => compute_sat_pco_body(&sat_obs.sat),
    };
    let pc = crate::estimators::rtk_iekf::satpos::compute_phase_centre_3d(src, &sat_obs.sat, rover_time, rx_eff, pco_body).ok()?;
    let sat_pos = pc.0;
    let sat_clock_m = compute_sat_clock_delay(src, &sat_obs.sat, rx_eff, sat_pos, rover_time);

    let (el_rad, az_rad) = compute_el_az(rx_pos, sat_pos);
    if sat_pos.norm() >= 1e6 && el_rad < 10.0_f64.to_radians() { return None; }
    let snr = sat_obs.get_snr(1).map_or(40.0, |s| s as f64);
    let sin_el = el_rad.sin().max(0.17);
    let snr_penalty = if snr < 40.0 { 10.0_f64.powf((40.0 - snr) / 10.0) } else { 1.0 };
    let (tropo_dry_m, tropo_map_wet) = compute_tropo_delay(rx_pos, sin_el);
    let variance_m2 = (1.0 / (sin_el * sin_el) * snr_penalty).clamp(0.5, 100.0);
    let cp_variance_m2 = (9e-6 / (sin_el * sin_el) * snr_penalty).clamp(9e-6, 1e-3);

    Some(RawObservation {
        satellite: sat_obs.sat.prn as u16, constellation_id: sat_obs.sat.constellation as u8,
        pr_l1, pr_l2, cp_l1, cp_l1_lli, cp_l2, doppler: sat_obs.get_doppler(1).unwrap_or(0.0),
        snr_dbhz: snr, sat_pos_ecef: sat_pos, sat_vel_ecef: Vector3::zeros(),
        sat_clock_m, f1, f2, freq_num: 0, elevation_rad: el_rad, azimuth_rad: az_rad,
        tropo_dry_m, tropo_map_wet, iono_l1_m: 0.0, variance_m2, cp_variance_m2,
    })
}

pub fn extract_raw_observations_with_source(
    rover: &EpochObs,
    src: &dyn crate::estimators::rtk_iekf::satpos::EphSource,
    approx_rx_pos: Option<Vector3<f64>>,
    sinex_bias: Option<&gneiss_parsers::sinex_bia::SinexBias>,
    antex: Option<&gneiss_parsers::antex::AntexDatabase>,
) -> Result<Vec<RawObservation>, String> {
    let rx_pos = approx_rx_pos.unwrap_or_else(|| Vector3::new(
        gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0,
    ));
    let raw_obs = rover.satellites.iter()
        .filter_map(|sat_obs| extract_single_sat_with_source(sat_obs, src, rx_pos, rover.time, sinex_bias, antex))
        .collect();
    Ok(raw_obs)
}

pub fn compute_spp_seeding(
    rover: &EpochObs,
    ephemerides: &[Ephemeris],
) -> Option<Vector3<f64>> {
    let config = crate::estimators::spp::SppConfig::default();
    crate::estimators::spp::compute_spp(rover, ephemerides, None, &config, None)
        .ok()
        .map(|state| state.position.vector)
}
