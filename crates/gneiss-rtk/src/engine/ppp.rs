use crate::engine::ppp_fg::PppIteratedEkf;
use crate::engine::processed_sat::ProcessedSat;
use crate::engine::{EngineError, ProcessingEngine};
use crate::filter::RtkState;
use chrono::TimeZone;
use gneiss_core::obs::EpochObs;
use gneiss_core::sat::Constellation;
use nalgebra::Vector3;

const LIGHT_SPEED: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

pub fn process_ppp<'a>(
    engine: &'a mut ProcessingEngine,
    rover_obs: &'a EpochObs,
) -> Result<&'a RtkState, EngineError> {
    if !valid_pos(engine) {
        return engine.process_spp(rover_obs);
    }
    let dt = (rover_obs.time.tow - engine.current_state.as_ref().unwrap().time.tow).max(0.0);
    engine.predict_state(dt);
    let state = engine.current_state.as_mut().unwrap();
    state.time = rover_obs.time;
    state.position.epoch = rover_obs.time;
    let sats = build_sats(engine, rover_obs);
    if sats.is_empty() {
        return Err(EngineError::InsufficientSatellites);
    }
    let state = engine.current_state.as_mut().unwrap();
    update_phase_ambiguities(state, &sats, rover_obs.time);
    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
    PppIteratedEkf::new().solve(state, &sats)?;
    state.epoch_count = state.epoch_count.saturating_add(1);
    let final_state = engine.current_state.as_ref().unwrap().clone();
    engine.state_history.push(final_state);
    engine.obs_history.push((rover_obs.clone(), None));
    Ok(engine.current_state.as_ref().unwrap())
}

pub(crate) fn valid_pos(engine: &ProcessingEngine) -> bool {
    if let Some(state) = &engine.current_state {
        state.position.vector.norm().is_normal() && state.position.vector.norm() >= 1000.0
    } else {
        false
    }
}

pub(crate) fn build_sats<'a>(engine: &ProcessingEngine, r_obs: &'a EpochObs) -> Vec<ProcessedSat<'a>> {
    let state = engine.current_state.as_ref().unwrap();
    let base_rcv_pos = Vector3::new(
        state.position.vector.x,
        state.position.vector.y,
        state.position.vector.z,
    );
    let rcv_pos =
        base_rcv_pos + gneiss_core::tides::solid_earth_tides_ecef(r_obs.time, base_rcv_pos);
    let rcv_llh = gneiss_core::coords::ecef_to_llh(rcv_pos);

    r_obs
        .satellites
        .iter()
        .filter_map(|sat_obs| process_single_sat(engine, r_obs, sat_obs, rcv_pos, rcv_llh))
        .collect()
}

fn process_single_sat<'a>(
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

    let (p1, mut p2, cp1, mut cp2, osb, is_if) =
        get_obs_and_corrections(engine, sat_obs, r_obs.time, f1, f2);
    let tau_pr = p1.unwrap_or(0.0) / LIGHT_SPEED;
    let t_nom = gneiss_core::time::GpsTime::new(r_obs.time.week, r_obs.time.tow - tau_pr);

    let (t_tx, dt_s, raw_pos, raw_vel) = compute_sat_state(engine, eph, sat_obs.sat, t_nom)?;
    let (sat_pos, sat_vel) =
        crate::engine::ppp_math::apply_earth_rotation(raw_pos, raw_vel, rcv_pos);

    let dist = (sat_pos - rcv_pos).norm();
    let (az, el) = gneiss_core::coords::az_el(rcv_llh, rcv_pos, sat_pos);
    if el < 0.261799 {
        return None;
    }

    let (tropo_dry, map_wet) = crate::engine::ppp_math::compute_tropo_dry(rcv_llh, el, r_obs.time);
    let klobuchar = engine.klobuchar_params.unwrap_or_default();
    let iono_delay = if is_if {
        0.0
    } else {
        gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(
            &klobuchar, rcv_llh, az, el, r_obs.time,
        )
    };

    let pcv = compute_pcv(
        engine, sat_obs, r_obs.time, f1, f2, is_if, rcv_pos, sat_pos, &mut p2, &mut cp2,
    );
    tracing::trace!(
        "PPP sat={}, is_if={}, dist={:.1}, clk_m={:.3}, pcv={:.4}",
        sat_obs.sat,
        is_if,
        dist,
        dt_s * LIGHT_SPEED,
        pcv
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

fn get_obs_and_corrections(
    engine: &ProcessingEngine,
    sat_obs: &gneiss_core::obs::SatObs,
    time: gneiss_core::time::GpsTime,
    f1: f64,
    f2: f64,
) -> (
    Option<f64>,
    Option<f64>,
    Option<f64>,
    Option<f64>,
    crate::engine::ppp_math::OsbCorrections,
    bool,
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

    let osb = crate::engine::ppp_math::apply_osb_corrections(
        engine.sinex_bias.as_ref(),
        sat_obs,
        time,
        f1,
        f2,
        f1_b,
        f2_b,
    );

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
    if precise && !engine.config.uduc_ar {
        if let (Some(v1), Some(v2)) = (p1, p2) {
            p1 = Some(crate::engine::ppp_math::compute_iono_free(f1, f2, v1, v2));
            is_if = true;
        }
        if let (Some(l1), Some(l2)) = (cp1, cp2) {
            cp1 = Some(
                crate::engine::ppp_math::compute_iono_free(
                    f1,
                    f2,
                    l1 * LIGHT_SPEED / f1,
                    l2 * LIGHT_SPEED / f2,
                ) / (LIGHT_SPEED / f1),
            );
        }
    }
    (p1, p2, cp1, cp2, osb, is_if)
}

fn compute_sat_state(
    engine: &ProcessingEngine,
    eph: &gneiss_core::ephemeris::Ephemeris,
    sat: gneiss_core::sat::SatelliteId,
    t_nom: gneiss_core::time::GpsTime,
) -> Option<(gneiss_core::time::GpsTime, f64, Vector3<f64>, Vector3<f64>)> {
    let precise = !engine.sp3_epochs.is_empty() || engine.clk_data.is_some();
    let mut dt_s = engine
        .clk_data
        .as_ref()
        .and_then(|c| c.get_clock_bias(sat, t_nom))
        .unwrap_or(0.0);
    let brdc_clk = eph.position(t_nom).2;

    let mut clk_found = dt_s != 0.0;
    if !precise {
        dt_s = brdc_clk;
        clk_found = true;
    }

    if precise && !clk_found {
        if let Some((_, _, sp3_clk)) =
            crate::engine::ssr::get_precise_orbit(&engine.sp3_epochs, sat, t_nom, 10)
        {
            if !sp3_clk.is_nan() && sp3_clk != 0.0 {
                dt_s = sp3_clk;
                clk_found = true;
            }
        }
    }
    if precise && !clk_found {
        return None;
    }

    let t_tx = gneiss_core::time::GpsTime::new(t_nom.week, t_nom.tow - dt_s);
    let (brdc_pos, brdc_vel, _, _) = eph.position(t_tx);
    let mut sat_pos: Vector3<f64> = brdc_pos;
    let mut sat_vel: Vector3<f64> = brdc_vel;

    if precise {
        let dt_rel = -2.0 * brdc_pos.dot(&brdc_vel) / (LIGHT_SPEED * LIGHT_SPEED);
        dt_s += dt_rel;
        if let Some((sp3_p, sp3_v, _)) =
            crate::engine::ssr::get_precise_orbit(&engine.sp3_epochs, sat, t_tx, 10)
        {
            sat_pos = sp3_p;
            sat_vel = sp3_v;
        } else {
            return None;
        }
    }
    Some((t_tx, dt_s, sat_pos, sat_vel))
}

fn compute_pcv(
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
    let pcv1 = -k.dot(&pco_ecef);

    if !is_if && sat_obs.get_observable_phase(2).is_some() {
        let pco2 = get_pco(code2);
        let pco2_ecef = sat_x * pco2.x + sat_y * pco2.y + sat_z * pco2.z;
        let diff = -k.dot(&pco2_ecef) - pcv1;
        if let Some(p) = p2.as_mut() {
            *p -= diff;
        }
        if let Some(cp) = cp2.as_mut() {
            *cp -= diff / (LIGHT_SPEED / f2);
        }
    }
    pcv1
}

pub(crate) fn update_phase_ambiguities(
    state: &mut RtkState,
    sats: &[ProcessedSat],
    t: gneiss_core::time::GpsTime,
) {
    for sat in sats.iter().filter(|s| s.cp1.unwrap_or(0.0) != 0.0) {
        let cp1 = sat.cp1.unwrap();
        let wup = gneiss_core::windup::phase_windup(
            sat.sat_pos_rot,
            gneiss_core::sun::sun_position_ecef(t),
            sat.rcv_pos_ecef,
            *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0),
        );
        state.windup.insert(sat.sat_obs.sat, wup);
        let l_meas = if sat.is_iono_free && sat.cp2.is_some() {
            crate::engine::ppp_math::compute_iono_free(
                (cp1 + wup) * sat.lam1,
                (sat.cp2.unwrap() + wup) * sat.lam2,
                sat.f1,
                sat.f2,
            )
        } else {
            (cp1 + wup) * sat.lam1
        };
        let prev = *state.locktimes.get(&(sat.sat_obs.sat, 1)).unwrap_or(&0);
        let (slip, new_lk) = crate::engine::ppp_math::detect_cycle_slip(sat.sat_obs, prev as u32);
        state.locktimes.insert((sat.sat_obs.sat, 1), new_lk as u16);
        if slip {
            for i in 0..4 {
                state.remove_ambiguity(sat.sat_obs.sat, i);
            }
        }
        let isb = match sat.sat_obs.sat.constellation {
            Constellation::Glonass => state.isb_glo,
            Constellation::Galileo => state.isb_gal,
            Constellation::Beidou => state.isb_bds,
            _ => 0.0,
        };
        let expected_base = sat.dist + state.rcv_clk_bias + isb - sat.dt_sat_m
            + sat.tropo_dry
            + state.zwd * sat.map_wet;
        if !sat.is_iono_free && sat.cp2.is_some() && sat.p2.is_some() {
            add_uduc_ambiguities(state, sat, cp1, wup, expected_base);
        } else {
            let exp = if sat.is_iono_free && sat.cp2.is_some() {
                expected_base
            } else {
                expected_base - sat.iono_delay
            };
            if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 0)) {
                state.add_ambiguity(sat.sat_obs.sat, 0, l_meas - exp, 10000.0);
            }
            state
                .last_observed
                .insert((sat.sat_obs.sat, 0), state.epoch_count as u32);
        }
    }
}

fn add_uduc_ambiguities(
    state: &mut RtkState,
    sat: &ProcessedSat,
    cp1: f64,
    wup: f64,
    expected_base: f64,
) {
    let p2 = sat.p2.unwrap();
    let p1 = sat.p1;
    let gamma = (sat.f1 * sat.f1) / (sat.f2 * sat.f2);
    let mut i1_est = (p2 - p1) / (gamma - 1.0);
    if i1_est.is_nan() || i1_est.abs() > 100.0 {
        i1_est = 0.0;
    }

    let l1_meas = (cp1 + wup) * sat.lam1;
    let l2_meas = (sat.cp2.unwrap() + wup) * sat.lam2;

    if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 3)) {
        state.add_ambiguity(sat.sat_obs.sat, 3, i1_est, 100.0);
    }
    if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 1)) {
        state.add_ambiguity(
            sat.sat_obs.sat,
            1,
            l1_meas - (expected_base - i1_est),
            10000.0,
        );
    }
    if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 2)) {
        state.add_ambiguity(
            sat.sat_obs.sat,
            2,
            l2_meas - (expected_base - gamma * i1_est),
            10000.0,
        );
    }

    for i in 1..4 {
        state
            .last_observed
            .insert((sat.sat_obs.sat, i), state.epoch_count as u32);
    }
}

#[cfg(test)]
mod osb_tests {
    use super::*;
    use gneiss_core::obs::ObsCode;
    use gneiss_core::obs::SatObs;
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;
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
        let bias = SinexBias::parse(Cursor::new(content)).unwrap();
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 2,
        };
        let f1 = 1575.42e6;
        let wl1 = LIGHT_SPEED / f1;
        let f2 = 1227.60e6;
        let wl2 = LIGHT_SPEED / f2;

        let mut obs = SatObs {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            observations: vec![],
        };
        obs.sat = sat;
        obs.observations.push(gneiss_core::obs::Observation {
            code: ObsCode::from_str("C1C").unwrap(),
            value: 10.0,
            lli: None,
            lock_time: None,
        });
        obs.observations.push(gneiss_core::obs::Observation {
            code: ObsCode::from_str("C2L").unwrap(),
            value: 20.0,
            lli: None,
            lock_time: None,
        });
        obs.observations.push(gneiss_core::obs::Observation {
            code: ObsCode::from_str("L1C").unwrap(),
            value: 30.0,
            lli: None,
            lock_time: None,
        });
        obs.observations.push(gneiss_core::obs::Observation {
            code: ObsCode::from_str("L2L").unwrap(),
            value: 40.0,
            lli: None,
            lock_time: None,
        });

        let res = crate::engine::ppp_math::apply_osb_corrections(
            Some(&bias),
            &obs,
            GpsTime::new(2156, 129600.0),
            f1,
            f2,
            1,
            2,
        );

        assert!((res.p1.unwrap() - (10.0 - 1.0 * 1e-9 * LIGHT_SPEED)).abs() < 1e-6);
        assert!((res.p2.unwrap() - (20.0 - 3.0 * 1e-9 * LIGHT_SPEED)).abs() < 1e-6);
        assert!((res.cp1.unwrap() - (30.0 - (2.0 * 1e-9 * LIGHT_SPEED) / wl1)).abs() < 1e-6);
        assert!((res.cp2.unwrap() - (40.0 - (4.0 * 1e-9 * LIGHT_SPEED) / wl2 - 0.25)).abs() < 1e-6);
    }
}

#[cfg(test)]
mod ppp_tests {
    use super::*;
    use crate::engine::{EngineConfig, ProcessingEngine};
    use gneiss_core::coords::{Coordinate, Datum, Frame};

    #[test]
    fn test_valid_pos() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        assert!(!valid_pos(&engine));

        let state = RtkState::new(
            gneiss_core::time::GpsTime::new(0, 0.0),
            Coordinate::new(
                Vector3::new(0.0, 0.0, f64::NAN),
                Datum::WGS84,
                Frame::ECEF,
                gneiss_core::time::GpsTime::new(0, 0.0),
            ),
            1.0,
        );
        engine.current_state = Some(state.clone());
        assert!(!valid_pos(&engine));

        engine.current_state.as_mut().unwrap().position.vector = Vector3::new(500.0, 500.0, 0.0);
        assert!(!valid_pos(&engine)); // norm = 707.1 < 1000.0

        engine.current_state.as_mut().unwrap().position.vector = Vector3::new(1000.0, 0.0, 0.0);
        assert!(valid_pos(&engine)); // norm = 1000.0 >= 1000.0
    }

    #[test]
    fn test_process_ppp_empty_sats() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let mut state = RtkState::new(
            gneiss_core::time::GpsTime::new(2156, 129000.0),
            Coordinate::new(
                Vector3::new(2000.0, 2000.0, 2000.0),
                Datum::WGS84,
                Frame::ECEF,
                gneiss_core::time::GpsTime::new(2156, 129000.0),
            ),
            1.0,
        );
        state.covariance[(0, 0)] = 10.0;
        engine.current_state = Some(state);

        let time = gneiss_core::time::GpsTime::new(2156, 129600.0);
        let obs = EpochObs {
            time,
            satellites: vec![],
        };

        // Should return InsufficientSatellites
        let res = process_ppp(&mut engine, &obs);
        assert!(matches!(res, Err(EngineError::InsufficientSatellites)));

        // Check state was predicted and time updated
        let state = engine.current_state.as_ref().unwrap();
        assert_eq!(state.time, time);

        // dt = 129600 - 129000 = 600.0
        assert_eq!(state.covariance[(0, 0)], 756000010.0);
    }

    #[test]
    fn test_process_single_sat() {
        use gneiss_core::ephemeris::Ephemeris;
        use gneiss_core::obs::{Observation, SatObs};
        use gneiss_core::sat::{Constellation, SatelliteId};
        

        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let t = gneiss_core::time::GpsTime::new(2156, 0.0);

        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat,
            toe: t,
            toc: t,
            af0: 1e-5,
            af1: 0.0,
            af2: 0.0,
            crs: 0.0,
            crc: 0.0,
            cuc: 0.0,
            cus: 0.0,
            cic: 0.0,
            cis: 0.0,
            m0: 0.0,
            e: 0.0,
            sqrt_a: 5500.0,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: std::f64::consts::PI / 4.0,
            idot: 0.0,
            omega: 0.0,
            tgd: 0.0,
            iode: 0,
            iodc: 0,
        });
        engine.ephemerides.push(eph);

        let mut sat_obs = SatObs {
            sat,
            observations: vec![],
        };
        sat_obs.observations.push(Observation {
            code: "C1C".parse().unwrap(),
            value: crate::engine::ppp_math::LIGHT_SPEED * 0.07,
            lli: None,
            lock_time: None,
        });
        sat_obs.observations.push(Observation {
            code: "L1C".parse().unwrap(),
            value: 1000.0,
            lli: None,
            lock_time: None,
        });

        let obs = EpochObs {
            time: t,
            satellites: vec![sat_obs.clone()],
        };

        let rcv_pos = Vector3::new(6000000.0, 0.0, 0.0);
        let rcv_llh = gneiss_core::coords::ecef_to_llh(rcv_pos);

        let res = process_single_sat(&engine, &obs, &sat_obs, rcv_pos, rcv_llh);
        assert!(res.is_some());
        let psat = res.unwrap();
        assert_eq!(psat.p1, crate::engine::ppp_math::LIGHT_SPEED * 0.07);
        assert_eq!(
            psat.f1,
            gneiss_core::signal::satellite_frequencies(sat, 0).0
        );
        assert_eq!(
            psat.f2,
            gneiss_core::signal::satellite_frequencies(sat, 0).1
        );

        println!("dist = {}", psat.dist);
        println!("dt_sat_m = {}", psat.dt_sat_m);
        println!("los.x = {}", psat.los.x);

        assert_eq!(psat.dist, 24250000.000301756);
        assert_eq!(psat.dt_sat_m, 2997.9245800000003);
        assert_eq!(psat.los.x, 0.9999999999372633);
        assert_eq!(
            psat.lam1,
            crate::engine::ppp_math::LIGHT_SPEED
                / gneiss_core::signal::satellite_frequencies(sat, 0).0
        );
        assert_eq!(
            psat.lam2,
            crate::engine::ppp_math::LIGHT_SPEED
                / gneiss_core::signal::satellite_frequencies(sat, 0).1
        );
    }

    #[test]
    fn test_process_single_sat_low_el() {
        use gneiss_core::ephemeris::Ephemeris;
        use gneiss_core::obs::{Observation, SatObs};
        use gneiss_core::sat::{Constellation, SatelliteId};

        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let t = gneiss_core::time::GpsTime::new(2156, 0.0);
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat,
            toe: t,
            toc: t,
            af0: 0.0,
            af1: 0.0,
            af2: 0.0,
            crs: 0.0,
            crc: 0.0,
            cuc: 0.0,
            cus: 0.0,
            cic: 0.0,
            cis: 0.0,
            m0: 0.0,
            e: 0.0,
            sqrt_a: 5500.0,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 0.0,
            idot: 0.0,
            omega: 0.0,
            tgd: 0.0,
            iode: 0,
            iodc: 0,
        });
        engine.ephemerides.push(eph);

        let mut sat_obs = SatObs {
            sat,
            observations: vec![],
        };
        sat_obs.observations.push(Observation {
            code: "C1C".parse().unwrap(),
            value: crate::engine::ppp_math::LIGHT_SPEED * 0.07,
            lli: None,
            lock_time: None,
        });

        let obs = EpochObs {
            time: t,
            satellites: vec![sat_obs.clone()],
        };
        // Place receiver on the opposite side of the earth, or just exactly under it but rotate so elevation is very low or negative
        // The sat is at x ~ 30000000, y = 0, z = 0.
        // If receiver is at y = 6378000, x = 0, elevation will be low.
        let rcv_pos = Vector3::new(0.0, 6378000.0, 0.0);
        let rcv_llh = gneiss_core::coords::ecef_to_llh(rcv_pos);

        let res = process_single_sat(&engine, &obs, &sat_obs, rcv_pos, rcv_llh);
        assert!(res.is_none());
    }

    #[test]
    fn test_get_obs_and_corrections() {
        use gneiss_core::obs::{Observation, SatObs};
        use gneiss_core::sat::{Constellation, SatelliteId};
        let engine = ProcessingEngine::new(EngineConfig::default());
        let t = gneiss_core::time::GpsTime::new(2156, 0.0);
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let mut sat_obs = SatObs {
            sat,
            observations: vec![],
        };
        sat_obs.observations.push(Observation {
            code: "C1C".parse().unwrap(),
            value: 1000.0,
            lli: None,
            lock_time: None,
        });
        sat_obs.observations.push(Observation {
            code: "C2W".parse().unwrap(),
            value: 2000.0,
            lli: None,
            lock_time: None,
        });
        sat_obs.observations.push(Observation {
            code: "L1C".parse().unwrap(),
            value: 3000.0,
            lli: None,
            lock_time: None,
        });
        sat_obs.observations.push(Observation {
            code: "L2W".parse().unwrap(),
            value: 4000.0,
            lli: None,
            lock_time: None,
        });

        let (p1, p2, cp1, cp2, _, is_if) =
            get_obs_and_corrections(&engine, &sat_obs, t, 1.5e9, 1.2e9);
        assert_eq!(p1, Some(1000.0));
        assert_eq!(p2, Some(2000.0));
        assert_eq!(cp1, Some(3000.0));
        assert_eq!(cp2, Some(4000.0));
        assert!(!is_if);
    }
}
