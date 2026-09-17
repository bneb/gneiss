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

fn apply_sagnac(sat_pos: Vector3<f64>, pr_m: f64) -> Vector3<f64> {
    let flight_time = pr_m / gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    let omega_tau = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * flight_time;
    let (cos_wt, sin_wt) = (omega_tau.cos(), omega_tau.sin());
    Vector3::new(sat_pos.x * cos_wt + sat_pos.y * sin_wt, -sat_pos.x * sin_wt + sat_pos.y * cos_wt, sat_pos.z)
}

fn compute_obs_noise(el_rad: f64, snr: f64) -> (f64, f64, f64) {
    let sin_el = el_rad.sin().max(0.17);
    let snr_pen = if snr < 40.0 { 10.0_f64.powf((40.0 - snr) / 10.0) } else { 1.0 };
    let v_pr = (1.0 / (sin_el * sin_el) * snr_pen).clamp(0.5, 100.0);
    let v_cp = (9e-6 / (sin_el * sin_el) * snr_pen).clamp(9e-6, 1e-3);
    (v_pr, v_cp, sin_el)
}

fn extract_broadcast_sat(
    sat_obs: &gneiss_core::obs::SatObs,
    ephemerides: &[Ephemeris],
    rover: &EpochObs,
    rx_pos: Vector3<f64>,
) -> Option<RawObservation> {
    let eph = select_best_ephemeris(ephemerides, sat_obs.sat, rover.time)?;
    let (f1, f2) = gneiss_core::signal::satellite_frequencies(sat_obs.sat, eph.freq_num());
    let p1_band = match sat_obs.sat.constellation {
        gneiss_core::sat::Constellation::Beidou => 2,
        _ => 1,
    };
    let pr_m = match sat_obs.get_observable(p1_band) {
        Some(pr) if pr > 1e6 => pr,
        _ => return None,
    };
    let t_tx = rover.time - (pr_m / gneiss_core::constants::SPEED_OF_LIGHT_M_S);
    let (_, _, sat_clk_rough, _) = eph.position(t_tx);
    let (sat_pos_raw, sat_vel, clk_err, _) = eph.position(t_tx - sat_clk_rough);
    let sat_clock_m = clk_err * gneiss_core::constants::SPEED_OF_LIGHT_M_S;
    let sat_pos = apply_sagnac(sat_pos_raw, pr_m);
    let (el_rad, az_rad) = compute_el_az(rx_pos, sat_pos);
    if sat_pos.norm() >= 1e6 && el_rad < 10.0_f64.to_radians() { return None; }
    let snr = sat_obs.get_snr(1).map_or(40.0, |s| s as f64);
    let (variance_m2, cp_variance_m2, sin_el) = compute_obs_noise(el_rad, snr);
    Some(RawObservation {
        satellite: sat_obs.sat.prn as u16, constellation_id: sat_obs.sat.constellation as u8,
        pr_l1: pr_m, pr_l2: sat_obs.get_observable(2),
        cp_l1: sat_obs.get_observable_phase_lli(1).map(|(cp, _)| cp),
        cp_l1_lli: sat_obs.get_observable_phase_lli(1).and_then(|(_, lli)| lli),
        cp_l2: sat_obs.get_observable_phase(2),
        doppler: sat_obs.get_doppler(1).unwrap_or(0.0), snr_dbhz: snr,
        sat_pos_ecef: sat_pos, sat_vel_ecef: sat_vel, sat_clock_m, f1, f2,
        freq_num: eph.freq_num(), elevation_rad: el_rad, azimuth_rad: az_rad,
        tropo_dry_m: 0.0, tropo_map_wet: 1.0 / sin_el, iono_l1_m: 0.0,
        variance_m2, cp_variance_m2,
    })
}

pub fn extract_raw_observations(
    rover: &EpochObs,
    ephemerides: &[Ephemeris],
    approx_rx_pos: Option<Vector3<f64>>,
) -> Result<Vec<RawObservation>, String> {
    let rx_pos = approx_rx_pos.unwrap_or_else(|| Vector3::new(
        gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0,
    ));
    let raw_obs = rover.satellites.iter()
        .filter_map(|sat_obs| extract_broadcast_sat(sat_obs, ephemerides, rover, rx_pos))
        .collect();
    Ok(raw_obs)
}

/// Extract raw observations using a generic ephemeris source (e.g. PreciseSrc with SP3/CLK).
fn lookup_satellite_bias(
    bias: Option<&gneiss_parsers::sinex_bia::SinexBias>,
    sat: SatelliteId,
    code: gneiss_core::obs::ObsCode,
    time: GpsTime,
) -> f64 {
    bias.and_then(|b| b.lookup_bias_m(sat, code, time)).unwrap_or(0.0)
}

fn compute_sat_pco_body(sat: &SatelliteId) -> Vector3<f64> {
    match sat.constellation {
        gneiss_core::sat::Constellation::Gps => match sat.prn {
            1 | 3 | 6 | 8 | 9 | 10 | 24 | 25 | 26 | 27 | 30 | 32 => Vector3::new(0.394, 0.0, 1.500),
            4 | 5 | 14 | 18 | 23 => Vector3::new(0.0, 0.0, 1.500),
            _ => Vector3::new(0.0, 0.0, 0.952),
        },
        gneiss_core::sat::Constellation::Galileo => Vector3::new(0.0, 0.0, 1.000),
        gneiss_core::sat::Constellation::Qzss => Vector3::new(0.0, 0.0, 1.000),
        _ => Vector3::zeros(),
    }
}

fn sat_antex_freq_codes(
    c: gneiss_core::sat::Constellation, ant: &gneiss_parsers::antex::AntennaPcv, f2: f64,
) -> (&'static str, &'static str) {
    match c {
        gneiss_core::sat::Constellation::Galileo => ("E01", if (f2 - 1207.140e6).abs() < 1e5 && ant.frequencies.contains_key("E07") { "E07" } else if ant.frequencies.contains_key("E05") { "E05" } else if ant.frequencies.contains_key("E07") { "E07" } else { "E08" }),
        gneiss_core::sat::Constellation::Glonass => ("R01", "R02"),
        gneiss_core::sat::Constellation::Beidou => (if ant.frequencies.contains_key("C02") { "C02" } else { "C01" },
            if (f2 - 1207.140e6).abs() < 1e5 && ant.frequencies.contains_key("C07") { "C07" } else if ant.frequencies.contains_key("C06") { "C06" } else if ant.frequencies.contains_key("C07") { "C07" } else { "C02" }),
        gneiss_core::sat::Constellation::Qzss => (if ant.frequencies.contains_key("J01") { "J01" } else { "G01" },
            if ant.frequencies.contains_key("J02") { "J02" } else { "G02" }),
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
        gneiss_core::sat::Constellation::Qzss => 'J',
        _ => return compute_sat_pco_body(sat),
    };
    let ant = match antex.find_satellite_gps(&format!("{}{:02}", sys_char, sat.prn), rover_time) {
        Some(a) => a,
        None => return compute_sat_pco_body(sat),
    };
    let (c1, c2) = sat_antex_freq_codes(sat.constellation, ant, f2);
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
        let (p_p, _) = src.position_at(sat, t_tx + 0.5).unwrap_or((sat_pos, 0.0));
        let (p_m, _) = src.position_at(sat, t_tx - 0.5).unwrap_or((sat_pos, 0.0));
        gneiss_geodesy::relativity::periodic_relativistic_range_correction(&sat_pos, &(p_p - p_m))
    } else { 0.0 };
    let shapiro = gneiss_geodesy::relativity::gravitational_shapiro_delay(&sat_pos, &rx_eff);
    sat_clock_s * gneiss_core::constants::SPEED_OF_LIGHT_M_S + rel_corr - shapiro
}

fn compute_el_az(rx_pos: Vector3<f64>, sat_pos: Vector3<f64>) -> (f64, f64) {
    if sat_pos.norm() < 1e6 { return (std::f64::consts::FRAC_PI_2, 0.0); }
    let ref_llh = gneiss_core::coords::ecef_to_llh(rx_pos);
    let d = sat_pos - rx_pos;
    let (slat, clat, slon, clon) = (ref_llh.x.sin(), ref_llh.x.cos(), ref_llh.y.sin(), ref_llh.y.cos());
    let (e, n, u) = (-slon * d.x + clon * d.y, -slat * clon * d.x - slat * slon * d.y + clat * d.z, clat * clon * d.x + clat * slon * d.y + slat * d.z);
    let dist = (e * e + n * n + u * u).sqrt().max(1.0);
    ((u / dist).clamp(-1.0, 1.0).asin(), e.atan2(n))
}

fn compute_tropo_delay(rx_pos: Vector3<f64>, sin_el: f64) -> (f64, f64) {
    let ref_llh = gneiss_core::coords::ecef_to_llh(rx_pos);
    let p_hpa = 1013.25 * (1.0 - 2.2557e-5 * ref_llh.z.clamp(-1000.0, 10000.0)).powi(5);
    let m = 1.001 / (0.002001 + sin_el * sin_el).sqrt();
    (0.002277 * p_hpa * m, m)
}

fn select_sat_bands(sat_obs: &gneiss_core::obs::SatObs, freq_num: i8) -> (u8, Option<u8>, f64, f64) {
    let obs = &sat_obs.observations;
    match sat_obs.sat.constellation {
        gneiss_core::sat::Constellation::Beidou => {
            let p1 = if obs.iter().any(|o| o.code.signal.freq_band == 2) { 2 } else { 1 };
            let (p2, f2) = if obs.iter().any(|o| o.code.signal.freq_band == 7) { (Some(7), 1207.140e6) }
            else if obs.iter().any(|o| o.code.signal.freq_band == 6) { (Some(6), 1268.520e6) }
            else if obs.iter().any(|o| o.code.signal.freq_band == 5) { (Some(5), 1176.450e6) }
            else { (None, 0.0) };
            (p1, p2, 1561.098e6, f2)
        }
        gneiss_core::sat::Constellation::Galileo => {
            let (p2, f2) = if obs.iter().any(|o| o.code.signal.freq_band == 5) { (Some(5), 1176.450e6) }
            else if obs.iter().any(|o| o.code.signal.freq_band == 7) { (Some(7), 1207.140e6) }
            else if obs.iter().any(|o| o.code.signal.freq_band == 6) { (Some(6), 1278.750e6) }
            else { (None, 0.0) };
            (1, p2, 1575.42e6, f2)
        }
        gneiss_core::sat::Constellation::Glonass => {
            let k = freq_num as f64;
            let f1 = 1602.0e6 + k * 0.5625e6;
            let f2 = 1246.0e6 + k * 0.4375e6;
            let p2 = obs.iter().find(|o| o.code.signal.freq_band == 2).map(|_| 2);
            (1, p2, f1, f2)
        }
        _ => {
            let (p2, f2) = if obs.iter().any(|o| o.code.signal.freq_band == 2) { (Some(2), 1227.60e6) }
            else if obs.iter().any(|o| o.code.signal.freq_band == 5) { (Some(5), 1176.450e6) }
            else { (None, 0.0) };
            (1, p2, 1575.42e6, f2)
        }
    }
}

fn extract_code_obs(
    sat_obs: &gneiss_core::obs::SatObs,
    p1_b: u8,
    p2_b: Option<u8>,
    bias: Option<&gneiss_parsers::sinex_bia::SinexBias>,
    t: GpsTime,
) -> Option<(f64, Option<f64>)> {
    let obs_p1 = sat_obs.observations.iter().find(|o| {
        o.code.obs_type == gneiss_core::obs::ObsType::Pseudorange
            && o.code.signal.freq_band == p1_b
            && o.value > 1e6
    })?;
    let pr_m = obs_p1.value - lookup_satellite_bias(bias, sat_obs.sat, obs_p1.code, t);
    let pr_l2 = p2_b.and_then(|b2| {
        sat_obs.observations.iter().find(|o| {
            o.code.obs_type == gneiss_core::obs::ObsType::Pseudorange && o.code.signal.freq_band == b2
        })
    }).map(|o| o.value - lookup_satellite_bias(bias, sat_obs.sat, o.code, t));
    Some((pr_m, pr_l2))
}

fn extract_phase_obs(
    sat_obs: &gneiss_core::obs::SatObs,
    p1_b: u8,
    p2_b: Option<u8>,
    f1: f64,
    f2: f64,
    bias: Option<&gneiss_parsers::sinex_bia::SinexBias>,
    t: GpsTime,
) -> (Option<f64>, Option<u8>, Option<f64>) {
    let wl_1 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / f1;
    let obs_l1 = sat_obs.observations.iter().find(|o| {
        o.code.obs_type == gneiss_core::obs::ObsType::CarrierPhase && o.code.signal.freq_band == p1_b
    });
    let (cp_l1, cp_l1_lli) = match obs_l1 {
        Some(o) => (Some(o.value - lookup_satellite_bias(bias, sat_obs.sat, o.code, t) / wl_1), o.lli),
        None => (None, None),
    };
    let cp_l2 = p2_b.filter(|_| f2 > 0.0).and_then(|b2| {
        let wl_2 = gneiss_core::constants::SPEED_OF_LIGHT_M_S / f2;
        sat_obs.observations.iter().find(|o| {
            o.code.obs_type == gneiss_core::obs::ObsType::CarrierPhase && o.code.signal.freq_band == b2
        }).map(|o| o.value - lookup_satellite_bias(bias, sat_obs.sat, o.code, t) / wl_2)
    });
    (cp_l1, cp_l1_lli, cp_l2)
}

fn compute_sat_nadir_pcv(
    antex: Option<&gneiss_parsers::antex::AntexDatabase>,
    sat: &SatelliteId,
    rover_time: GpsTime,
    rx_pos: Vector3<f64>,
    sat_pos: Vector3<f64>,
) -> f64 {
    let db = match antex {
        Some(d) => d,
        None => return 0.0,
    };
    let sys_char = match sat.constellation {
        gneiss_core::sat::Constellation::Gps => 'G',
        gneiss_core::sat::Constellation::Glonass => 'R',
        gneiss_core::sat::Constellation::Galileo => 'E',
        gneiss_core::sat::Constellation::Beidou => 'C',
        gneiss_core::sat::Constellation::Qzss => 'J',
        _ => return 0.0,
    };
    let prn_str = format!("{}{:02}", sys_char, sat.prn);
    let ant = match db.find_satellite_gps(&prn_str, rover_time) {
        Some(a) => a,
        None => return 0.0,
    };
    let (c1, _) = sat_antex_freq_codes(sat.constellation, ant, 0.0);
    let (el_rad, _) = compute_el_az(rx_pos, sat_pos);
    let ratio = if sat_pos.norm() > 1e6 { (rx_pos.norm() / sat_pos.norm()) * el_rad.cos() } else { 0.0 };
    let nadir_deg = ratio.clamp(-1.0, 1.0).asin().to_degrees();
    db.satellite_pcv_nadir_m(&prn_str, rover_time, nadir_deg, c1)
}

fn compute_sat_pos_and_clock(
    src: &dyn crate::estimators::rtk_iekf::satpos::EphSource,
    sat: &SatelliteId,
    rx_pos: Vector3<f64>,
    rover_time: GpsTime,
    f1: f64,
    f2: f64,
    antex: Option<&gneiss_parsers::antex::AntexDatabase>,
) -> Option<(Vector3<f64>, f64)> {
    let tide_m = if rx_pos.norm() > 1e6 {
        gneiss_geodesy::tides::solid_earth_tide(rx_pos, rover_time.tow, rover_time.week)
    } else {
        Vector3::zeros()
    };
    let rx_eff = rx_pos + tide_m;
    let pco_body = antex.map_or_else(|| compute_sat_pco_body(sat), |db| compute_sat_pco_from_antex(db, sat, rover_time, f1, f2));
    let pc_init = crate::estimators::rtk_iekf::satpos::compute_phase_centre_3d(src, sat, rover_time, rx_eff, pco_body).ok()?;
    let pcv_nadir_m = compute_sat_nadir_pcv(antex, sat, rover_time, rx_pos, pc_init.0);
    let pc = crate::estimators::rtk_iekf::satpos::compute_phase_centre_3d_with_pcv(src, sat, rover_time, rx_eff, pco_body, pcv_nadir_m).ok()?;
    let sat_pos = pc.0;
    let sat_clock_m = compute_sat_clock_delay(src, sat, rx_eff, sat_pos, rover_time);
    Some((sat_pos, sat_clock_m))
}

fn galileo_bgd_offset_m(ephs: Option<&[Ephemeris]>, sat: SatelliteId, t: GpsTime, f2: f64) -> f64 {
    if (f2 - 1207.140e6).abs() > 1e5 { return 0.0; }
    let Some(e) = ephs.and_then(|e| select_best_ephemeris(e, sat, t)) else { return 0.0; };
    let k_5a = 1176.45_f64.powi(2) / (1575.42_f64.powi(2) - 1176.45_f64.powi(2));
    let k_5b = 1207.14_f64.powi(2) / (1575.42_f64.powi(2) - 1207.14_f64.powi(2));
    -gneiss_core::constants::SPEED_OF_LIGHT_M_S * (k_5a * e.tgd() - k_5b * e.bgd_e5b())
}

#[allow(clippy::too_many_arguments)]
fn extract_single_sat_with_source(
    sat_obs: &gneiss_core::obs::SatObs,
    src: &dyn crate::estimators::rtk_iekf::satpos::EphSource,
    rx_pos: Vector3<f64>,
    rover_time: GpsTime,
    sinex_bias: Option<&gneiss_parsers::sinex_bia::SinexBias>,
    antex: Option<&gneiss_parsers::antex::AntexDatabase>,
    ephemerides: Option<&[Ephemeris]>,
    enable_glonass: bool,
    enable_galileo: bool,
    ep: u32,
) -> Option<RawObservation> {
    let c = sat_obs.sat.constellation;
    let is_supp = match c {
        gneiss_core::sat::Constellation::Gps => true,
        gneiss_core::sat::Constellation::Glonass => enable_glonass,
        gneiss_core::sat::Constellation::Galileo => enable_galileo,
        _ => false,
    };
    if !is_supp { return None; }
    let freq_num = if c == gneiss_core::sat::Constellation::Glonass {
        ephemerides.and_then(|ephs| select_best_ephemeris(ephs, sat_obs.sat, rover_time))
            .map(|e| e.freq_num()).unwrap_or(0)
    } else { 0 };
    let (p1_b, p2_b, f1, f2) = select_sat_bands(sat_obs, freq_num);
    let (pr_l1, pr_l2) = extract_code_obs(sat_obs, p1_b, p2_b, sinex_bias, rover_time)?;
    let (cp_l1, cp_l1_lli, cp_l2) = extract_phase_obs(sat_obs, p1_b, p2_b, f1, f2, sinex_bias, rover_time);
    if pr_l2.is_none() || cp_l2.is_none() { return None; }
    record_mw_if_dual_freq(sat_obs.sat, cp_l1, cp_l2, pr_l1, pr_l2, f1, f2, cp_l1_lli, ep);

    let (sat_pos, base_clk_m) = compute_sat_pos_and_clock(src, &sat_obs.sat, rx_pos, rover_time, f1, f2, antex)?;
    let sat_clock_m = base_clk_m + if c == gneiss_core::sat::Constellation::Galileo {
        galileo_bgd_offset_m(ephemerides, sat_obs.sat, rover_time, f2)
    } else { 0.0 };
    let (el_rad, az_rad) = compute_el_az(rx_pos, sat_pos);
    if sat_pos.norm() >= 1e6 && el_rad < 15.0_f64.to_radians() { return None; }
    let snr = sat_obs.get_snr(1).map_or(40.0, |s| s as f64);
    let (variance_m2, cp_variance_m2, sin_el) = compute_obs_noise(el_rad, snr);
    let (tropo_dry_m, tropo_map_wet) = compute_tropo_delay(rx_pos, sin_el);

    Some(RawObservation {
        satellite: sat_obs.sat.prn as u16, constellation_id: sat_obs.sat.constellation as u8,
        pr_l1, pr_l2, cp_l1, cp_l1_lli, cp_l2, doppler: sat_obs.get_doppler(1).unwrap_or(0.0),
        snr_dbhz: snr, sat_pos_ecef: sat_pos, sat_vel_ecef: Vector3::zeros(),
        sat_clock_m, f1, f2, freq_num, elevation_rad: el_rad, azimuth_rad: az_rad,
        tropo_dry_m, tropo_map_wet, iono_l1_m: 0.0, variance_m2, cp_variance_m2,
    })
}

thread_local! {
    static EPOCH_TRACKER: std::cell::RefCell<(f64, u32)> = const { std::cell::RefCell::new((-1.0, 0)) };
}

pub fn reset_epoch_tracker() {
    EPOCH_TRACKER.with(|c| *c.borrow_mut() = (-1.0, 0));
}

fn advance_epoch(tow: f64) -> u32 {
    EPOCH_TRACKER.with(|c| {
        let mut g = c.borrow_mut();
        if (tow - g.0).abs() > 0.001 {
            g.0 = tow;
            g.1 += 1;
        }
        g.1
    })
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn record_mw_if_dual_freq(
    sat: SatelliteId, cp1: Option<f64>, cp2: Option<f64>, pr1: f64, pr2: Option<f64>,
    f1: f64, f2: f64, lli: Option<u8>, ep: u32,
) {
    if let (Some(c1), Some(c2), Some(p2)) = (cp1, cp2, pr2) {
        if f1 > f2 && f2 > 0.0 {
            let lambda_wl = gneiss_core::constants::SPEED_OF_LIGHT_M_S / (f1 - f2);
            let p_nl = (f1 * pr1 + f2 * p2) / (f1 + f2);
            let mw = (c1 - c2) - p_nl / lambda_wl;
            super::ar_handler::record_mw_sample(
                sat.constellation as u8, sat.prn as u16, mw, ep, lli.unwrap_or(0) & 1 != 0,
            );
        }
    }
}

fn estimate_rx_clock_bias_s(
    rover: &EpochObs,
    src: &dyn crate::estimators::rtk_iekf::satpos::EphSource,
    rx_pos: Vector3<f64>,
) -> f64 {
    if rx_pos.norm() < 1e6 { return 0.0; }
    let mut diffs = Vec::new();
    for sat in &rover.satellites {
        let pr = match sat.get_observable(1) {
            Some(p) if p > 1e6 => p,
            _ => continue,
        };
        if let Some((p0, _)) = src.position_at(&sat.sat, rover.time) {
            diffs.push(pr - (p0 - rx_pos).norm());
        }
    }
    if diffs.is_empty() { return 0.0; }
    diffs.sort_by(|a, b| a.total_cmp(b));
    diffs[diffs.len() / 2] / gneiss_core::constants::SPEED_OF_LIGHT_M_S
}

#[allow(clippy::too_many_arguments)]
pub fn extract_raw_observations_with_source(
    rover: &EpochObs,
    src: &dyn crate::estimators::rtk_iekf::satpos::EphSource,
    approx_rx_pos: Option<Vector3<f64>>,
    sinex_bias: Option<&gneiss_parsers::sinex_bia::SinexBias>,
    antex: Option<&gneiss_parsers::antex::AntexDatabase>,
    ephemerides: Option<&[Ephemeris]>,
    enable_glonass: bool,
    enable_galileo: bool,
) -> Result<Vec<RawObservation>, String> {
    let rx_pos = approx_rx_pos.unwrap_or_else(|| Vector3::new(
        gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0,
    ));
    let rx_clk_s = estimate_rx_clock_bias_s(rover, src, rx_pos);
    let true_time = rover.time - rx_clk_s;
    let ep = advance_epoch(rover.time.tow);
    let raw_obs = rover.satellites.iter()
        .filter_map(|sat_obs| extract_single_sat_with_source(
            sat_obs, src, rx_pos, true_time, sinex_bias, antex, ephemerides, enable_glonass, enable_galileo, ep,
        ))
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
        .and_then(|s| ((s.position.vector.norm() - 6_371_000.0).abs() < 25_000.0).then_some(s.position.vector))
}
