use crate::engine::processed_sat::ProcessedSat;
use crate::engine::ProcessingEngine;
use chrono::TimeZone;
use gneiss_core::obs::EpochObs;
use gneiss_core::sat::Constellation;
use nalgebra::Vector3;

const LIGHT_SPEED: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

// ---------------------------------------------------------------------------
// Receiver antenna functions
// ---------------------------------------------------------------------------

/// Compute receiver antenna phase center offset in ECEF.
/// Returns zero vector if ANTEX not loaded or antenna type not found.
pub(crate) fn compute_receiver_pco(
    antex: Option<&gneiss_parsers::antex::AntexDatabase>,
    antenna_type: Option<&str>,
    freq_code: &str,
    rcv_llh: Vector3<f64>,
) -> Vector3<f64> {
    let db = match antex {
        Some(d) => d,
        None => return Vector3::zeros(),
    };
    let ant_type = match antenna_type {
        Some(t) => t,
        None => return Vector3::zeros(),
    };
    let antenna = match db.get_antenna(ant_type) {
        Some(a) => a,
        None => return Vector3::zeros(),
    };
    let freq = match antenna.frequencies.get(freq_code) {
        Some(f) => f,
        None => return Vector3::zeros(),
    };
    let pco_neu = &freq.pco;
    let (lat, lon) = (rcv_llh.x, rcv_llh.y);
    let (clat, slat) = (lat.cos(), lat.sin());
    let (clon, slon) = (lon.cos(), lon.sin());
    Vector3::new(
        -slat * clon * pco_neu.x - slon * pco_neu.y + clat * clon * pco_neu.z,
        -slat * slon * pco_neu.x + clon * pco_neu.y + clat * slon * pco_neu.z,
        clat * pco_neu.x + slat * pco_neu.z,
    )
}

/// Compute receiver Phase Center Variation (PCV) in meters based on satellite
/// elevation angle.  ANTEX `noazi` values are sampled at zenith angles from
/// `zen1` to `zen2` with step `dzen` (all in degrees).  The correction is
/// returned in meters (ANTEX stores in mm).
///
/// Returns 0.0 if ANTEX is not loaded, antenna type not found, or no noazi data.
pub(crate) fn compute_receiver_pcv(
    antex: Option<&gneiss_parsers::antex::AntexDatabase>,
    antenna_type: Option<&str>,
    freq_code: &str,
    el_rad: f64,
) -> f64 {
    let db = match antex {
        Some(d) => d,
        None => return 0.0,
    };
    let ant_type = match antenna_type {
        Some(t) => t,
        None => return 0.0,
    };
    let antenna = match db.get_antenna(ant_type) {
        Some(a) => a,
        None => return 0.0,
    };
    let freq = match antenna.frequencies.get(freq_code) {
        Some(f) => f,
        None => return 0.0,
    };
    if freq.noazi.is_empty() || antenna.dzen <= 0.0 {
        return 0.0;
    }

    // Convert elevation to zenith angle in degrees
    let zenith_deg = 90.0 - el_rad.to_degrees();
    let zenith_deg = zenith_deg.clamp(antenna.zen1, antenna.zen2);

    let idx_f = (zenith_deg - antenna.zen1) / antenna.dzen;
    let idx0 = idx_f.floor() as usize;
    let idx1 = (idx0 + 1).min(freq.noazi.len() - 1);

    if idx0 >= freq.noazi.len() {
        return 0.0;
    }

    let w1 = idx_f - idx0 as f64;
    let w0 = 1.0 - w1;
    let pcv_mm = if idx0 == idx1 {
        freq.noazi[idx0]
    } else {
        w0 * freq.noazi[idx0] + w1 * freq.noazi[idx1]
    };

    pcv_mm / 1000.0 // mm → m
}

// ---------------------------------------------------------------------------
// Satellite processing pipeline
// ---------------------------------------------------------------------------

pub(crate) fn build_sats<'a>(
    engine: &ProcessingEngine,
    r_obs: &'a EpochObs,
) -> Vec<ProcessedSat<'a>> {
    let state = engine.current_state.as_ref().unwrap();
    let base_rcv_pos = Vector3::new(
        state.position.vector.x,
        state.position.vector.y,
        state.position.vector.z,
    );
    let rcv_pos =
        base_rcv_pos + gneiss_core::tides::solid_earth_tides_ecef(r_obs.time, base_rcv_pos);
    let rcv_llh = gneiss_core::coords::ecef_to_llh(rcv_pos);
    // Apply receiver antenna PCO if available, using GPS L1 code as reference
    let rcv_pco = compute_receiver_pco(
        engine.antex.as_ref(),
        engine.config.receiver_antenna_type.as_deref(),
        "G01",
        rcv_llh,
    );
    let rcv_pos = rcv_pos + rcv_pco;

    r_obs
        .satellites
        .iter()
        .filter_map(|sat_obs| process_single_sat(engine, r_obs, sat_obs, rcv_pos, rcv_llh))
        .collect()
}

pub(crate) fn process_single_sat<'a>(
    engine: &ProcessingEngine,
    r_obs: &'a EpochObs,
    sat_obs: &'a gneiss_core::obs::SatObs,
    rcv_pos: Vector3<f64>,
    rcv_llh: Vector3<f64>,
) -> Option<ProcessedSat<'a>> {
    let eph = engine.ephemerides.iter().find(|e| e.sat() == sat_obs.sat)?;
    let (f1, mut f2) = gneiss_core::signal::satellite_frequencies(sat_obs.sat, eph.freq_num());
    if f2 == 0.0 {
        f2 = f1;
    }

    let (p1, mut p2, cp1, mut cp2, osb, is_if, actual_f2) =
        get_obs_and_corrections(engine, sat_obs, r_obs.time, f1, f2);
    f2 = actual_f2;
    let tau_pr = p1.unwrap_or(0.0) / LIGHT_SPEED;
    let t_nom = gneiss_core::time::GpsTime::new(r_obs.time.week, r_obs.time.tow - tau_pr);

    let (t_tx, dt_s, raw_pos, raw_vel) = crate::engine::ssr::compute_sat_state(
        &engine.sp3_epochs,
        engine.clk_data.as_ref(),
        eph,
        sat_obs.sat,
        t_nom,
    )?;
    let (sat_pos, sat_vel) =
        crate::engine::ppp_math::apply_earth_rotation(raw_pos, raw_vel, rcv_pos);

    let dist = (sat_pos - rcv_pos).norm();
    let (az, el) = gneiss_core::coords::az_el(rcv_llh, rcv_pos, sat_pos);
    // 5° elevation mask — industry standard for urban PPP (RTKLIB, NovAtel).
    // At 15°, most satellites are rejected in dense urban canyons.
    if el < 5.0_f64.to_radians() {
        return None;
    }

    let (tropo_dry, map_wet) = crate::engine::ppp_math::compute_tropo_dry(
        rcv_llh,
        el,
        r_obs.time,
        engine.tropo_mapper.as_ref(),
    );
    let iono_delay = if is_if {
        0.0
    } else {
        match engine.config.iono_model {
            crate::engine::types::IonosphereModel::Klobuchar => {
                let klobuchar = engine.klobuchar_params.unwrap_or_default();
                gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(
                    &klobuchar, rcv_llh, az, el, r_obs.time,
                )
            }
            crate::engine::types::IonosphereModel::Ionex => {
                if let Some(ref grid) = engine.ionex_grid {
                    if !grid.tec_maps.is_empty() {
                        let first_map = &grid.tec_maps[0].tec;
                        if !first_map.is_empty() && !first_map[0].is_empty() {
                            let maps_ref: Vec<_> = engine
                                .ionex_maps
                                .iter()
                                .map(|(t, m)| (*t, m))
                                .collect();
                            gneiss_core::atmosphere::AtmosphereModel::iono_ionex(
                                &maps_ref,
                                grid.lat1,
                                grid.lat2,
                                grid.dlat,
                                grid.lon1,
                                grid.lon2,
                                grid.dlon,
                                grid.height_km,
                                rcv_llh,
                                az,
                                el,
                                r_obs.time,
                            )
                        } else {
                            let klobuchar = engine.klobuchar_params.unwrap_or_default();
                            gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(
                                &klobuchar, rcv_llh, az, el, r_obs.time,
                            )
                        }
                    } else {
                        let klobuchar = engine.klobuchar_params.unwrap_or_default();
                        gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(
                            &klobuchar, rcv_llh, az, el, r_obs.time,
                        )
                    }
                } else {
                    let klobuchar = engine.klobuchar_params.unwrap_or_default();
                    gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(
                        &klobuchar, rcv_llh, az, el, r_obs.time,
                    )
                }
            }
        }
    };

    let pcv = compute_pcv(
        engine, sat_obs, r_obs.time, f1, f2, is_if, rcv_pos, sat_pos, &mut p2, &mut cp2,
    );
    // Bug 12 fix: apply receiver elevation-dependent phase center variation.
    let rcv_pcv = compute_receiver_pcv(
        engine.antex.as_ref(),
        engine.config.receiver_antenna_type.as_deref(),
        "G01",
        el,
    );
    tracing::trace!(
        "PPP sat={}, is_if={}, dist={:.1}, clk_m={:.3}, pcv={:.4}, rcv_pcv={:.4}",
        sat_obs.sat,
        is_if,
        dist,
        dt_s * LIGHT_SPEED,
        pcv,
        rcv_pcv
    );

    let snr = sat_obs.get_snr(1).unwrap_or(45) as f64;
    let doppler = sat_obs.get_doppler(1).unwrap_or(0.0);

    Some(ProcessedSat {
        sat_obs,
        dt_sat_m: dt_s * LIGHT_SPEED,
        p1: p1?,
        p2,
        is_iono_free: is_if,
        cp1,
        cp2,
        osb_p1: osb.osb_p1,
        osb_p2: osb.osb_p2,
        osb_cp1: osb.osb_cp1,
        osb_cp2: osb.osb_cp2,
        los: (sat_pos - rcv_pos) / dist,
        dist: dist - pcv,
        el,
        snr,
        doppler,
        lam1: LIGHT_SPEED / f1,
        lam2: LIGHT_SPEED / f2,
        tropo_dry,
        map_wet,
        iono_delay,
        f1,
        f2,
        sat_pos_rot: sat_pos,
        sat_vel,
        sat_clock_drift: eph.position(t_tx).3,
        rcv_pos_ecef: rcv_pos,
        pcv_correction: pcv,
    })
}

#[allow(clippy::type_complexity)]
pub(crate) fn get_obs_and_corrections(
    engine: &ProcessingEngine,
    sat_obs: &gneiss_core::obs::SatObs,
    time: gneiss_core::time::GpsTime,
    _f1: f64,
    _f2: f64,
) -> (
    Option<f64>,
    Option<f64>,
    Option<f64>,
    Option<f64>,
    crate::engine::ppp_math::OsbCorrections,
    bool,
    f64, // actual f2 (may differ from satellite_frequencies if L5 fallback)
) {
    let f1_b = if sat_obs.sat.constellation == Constellation::Beidou {
        2
    } else {
        1
    };
    let f2_b = match sat_obs.sat.constellation {
        Constellation::Galileo | Constellation::Beidou => 7,
        _ => 2,
    };
    // Try standard bands first
    let mut osb = crate::engine::ppp_math::apply_osb_corrections(
        engine.sinex_bias.as_ref(),
        sat_obs,
        time,
        _f1,
        _f2,
        f1_b,
        f2_b,
    );
    // L5 fallback: if no L2 observation (smartphones track L5, not L2),
    // retry with band 5 (GPS L5 / Galileo E5a / BDS B2a at 1176.45 MHz)
    let mut actual_f2 = _f2;
    if osb.p2.is_none() && osb.cp2.is_none() {
        let l5_band = 5u8;
        let l5_freq = gneiss_core::signal::get_frequency(sat_obs.sat, l5_band, 0);
        if l5_freq != _f2 {
            osb = crate::engine::ppp_math::apply_osb_corrections(
                engine.sinex_bias.as_ref(),
                sat_obs,
                time,
                _f1,
                l5_freq,
                f1_b,
                l5_band,
            );
            if osb.p2.is_some() || osb.cp2.is_some() {
                actual_f2 = l5_freq;
            }
        }
    }

    let mut p1 = osb.p1;
    let mut p2 = osb.p2;
    let mut cp1 = osb.cp1;
    let cp2 = osb.cp2;
    if engine.sinex_bias.is_none() {
        if let Some(v) = p1.as_mut() {
            if let Some(&d) = engine.dcbs.get(&(sat_obs.sat, "P1C1".to_string())) {
                *v -= d * 1e-9 * LIGHT_SPEED;
            }
        }
        if let Some(v) = p2.as_mut() {
            if let Some(&d) = engine.dcbs.get(&(sat_obs.sat, "P2C2".to_string())) {
                *v -= d * 1e-9 * LIGHT_SPEED;
            }
        }
    }

    let mut is_if = false;
    let precise = !engine.sp3_epochs.is_empty() || engine.clk_data.is_some();
    // When AR is enabled, skip iono-free combination: use raw single-frequency
    // observables + estimated iono state (UDUC mode). The IF combination destroys
    // integer ambiguity information needed for WL/NL cascade AR.
    if precise && !engine.config.uduc_ar && !engine.config.enable_ar {
        if let (Some(v1), Some(v2)) = (p1, p2) {
            p1 = Some(crate::engine::ppp_math::compute_iono_free(
                _f1, actual_f2, v1, v2,
            ));
            is_if = true;
        }
        if let (Some(l1), Some(l2)) = (cp1, cp2) {
            cp1 = Some(
                crate::engine::ppp_math::compute_iono_free(
                    _f1,
                    actual_f2,
                    l1 * LIGHT_SPEED / _f1,
                    l2 * LIGHT_SPEED / actual_f2,
                ) / (LIGHT_SPEED / _f1),
            );
        }
    }
    (p1, p2, cp1, cp2, osb, is_if, actual_f2)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_pcv(
    engine: &ProcessingEngine,
    sat_obs: &gneiss_core::obs::SatObs,
    t: gneiss_core::time::GpsTime,
    _f1: f64,
    f2: f64,
    is_if: bool,
    rcv_pos: Vector3<f64>,
    sat_pos: Vector3<f64>,
    p2: &mut Option<f64>,
    cp2: &mut Option<f64>,
) -> f64 {
    let precise = !engine.sp3_epochs.is_empty() || engine.clk_data.is_some();
    if !precise || engine.antex.is_none() {
        return 0.0;
    }

    let db = engine.antex.as_ref().unwrap();
    let utc_time = chrono::Utc.with_ymd_and_hms(1980, 1, 6, 0, 0, 0).unwrap()
        + chrono::Duration::seconds(t.week as i64 * 604800 + t.tow as i64);
    let ant = match db.find_satellite(&sat_obs.sat.to_string(), utc_time) {
        Some(a) => a,
        None => return 0.0,
    };

    let code1 = match sat_obs.sat.constellation {
        Constellation::Glonass => "R01",
        Constellation::Galileo => "E01",
        Constellation::Beidou => "C02",
        _ => "G01",
    };
    let code2 = match sat_obs.sat.constellation {
        Constellation::Glonass => "R02",
        Constellation::Galileo => "E05",
        Constellation::Beidou => "C07",
        _ => "G02",
    };

    let get_pco = |c| {
        ant.frequencies
            .get(c)
            .map(|p| p.pco / 1000.0)
            .unwrap_or_else(|| {
                ant.frequencies
                    .values()
                    .next()
                    .map(|p| p.pco / 1000.0)
                    .unwrap_or(Vector3::zeros())
            })
    };
    let pco1 = get_pco(code1);

    let k = (rcv_pos - sat_pos).normalize();
    let sat_z = -sat_pos.normalize();
    let e_sun = (gneiss_core::sun::sun_position_ecef(t) - sat_pos).normalize();
    let sat_y = sat_z.cross(&e_sun).normalize();
    let sat_x = sat_y.cross(&sat_z).normalize();

    let pco_ecef = sat_x * pco1.x + sat_y * pco1.y + sat_z * pco1.z;
    let mut pcv1 = -k.dot(&pco_ecef);

    // Bug 11 fix: apply nadir-angle-dependent satellite PCV from ANTEX noazi table.
    // The ANTEX "zenith" angle for satellite antennas is the nadir angle measured
    // from the satellite's nadir direction (toward Earth center).
    let nadir_deg = sat_z.dot(&k).clamp(-1.0, 1.0).acos().to_degrees();
    if let Some(freq1) = ant.frequencies.get(code1) {
        if ant.dzen > 0.0 && !freq1.noazi.is_empty() {
            let idx_f = ((nadir_deg - ant.zen1) / ant.dzen).max(0.0);
            let idx0 = libm::floor(idx_f) as usize;
            let idx1 = idx0 + 1;
            let pcv_mm = if idx0 >= freq1.noazi.len() {
                *freq1.noazi.last().unwrap_or(&0.0)
            } else if idx1 >= freq1.noazi.len() {
                freq1.noazi[idx0]
            } else {
                let w = idx_f - idx0 as f64;
                freq1.noazi[idx0] * (1.0 - w) + freq1.noazi[idx1] * w
            };
            pcv1 += pcv_mm / 1000.0; // mm → m
        }
    }

    if !is_if && sat_obs.get_observable_phase(2).is_some() {
        let pco2 = get_pco(code2);
        let pco2_ecef = sat_x * pco2.x + sat_y * pco2.y + sat_z * pco2.z;
        let mut diff = -k.dot(&pco2_ecef) - pcv1;
        // Apply nadir-dependent PCV for L2 as well
        if let Some(freq2) = ant.frequencies.get(code2) {
            if ant.dzen > 0.0 && !freq2.noazi.is_empty() {
                let idx_f = ((nadir_deg - ant.zen1) / ant.dzen).max(0.0);
                let idx0 = libm::floor(idx_f) as usize;
                let idx1 = idx0 + 1;
                let pcv_mm = if idx0 >= freq2.noazi.len() {
                    *freq2.noazi.last().unwrap_or(&0.0)
                } else if idx1 >= freq2.noazi.len() {
                    freq2.noazi[idx0]
                } else {
                    let w = idx_f - idx0 as f64;
                    freq2.noazi[idx0] * (1.0 - w) + freq2.noazi[idx1] * w
                };
                diff += pcv_mm / 1000.0;
            }
        }
        if let Some(p) = p2.as_mut() {
            *p -= diff;
        }
        if let Some(cp) = cp2.as_mut() {
            *cp -= diff / (LIGHT_SPEED / f2);
        }
    }
    pcv1
}
