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
        let t_tx = gneiss_core::time::GpsTime::new(
            rover.time.week,
            rover.time.tow - pr_m / gneiss_core::constants::SPEED_OF_LIGHT_M_S,
        );

        // Get rough satellite clock to refine transmit time.
        let (_, _, sat_clk_err_rough, _) = eph.position(t_tx);

        // Pass 2: refine transmit time with satellite clock correction.
        let t_tx_true = gneiss_core::time::GpsTime::new(
            rover.time.week,
            t_tx.tow - sat_clk_err_rough,
        );
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

pub fn compute_spp_seeding(
    rover: &EpochObs,
    ephemerides: &[Ephemeris],
) -> Option<Vector3<f64>> {
    let config = crate::estimators::spp::SppConfig::default();
    crate::estimators::spp::compute_spp(rover, ephemerides, None, &config, None)
        .ok()
        .map(|state| state.position.vector)
}
