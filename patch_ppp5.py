import os

with open("crates/gneiss-rtk/src/engine/ppp.rs", "r") as f:
    orig = f.read()

start_str = "fn build_sats<'a>(engine: &ProcessingEngine, rover_obs: &'a EpochObs) -> Vec<ProcessedSat<'a>> {"
end_str = "pub struct OsbCorrections {"

build_sats_new = """fn build_sats<'a>(engine: &ProcessingEngine, r_obs: &'a EpochObs) -> Vec<ProcessedSat<'a>> {
    let state = engine.current_state.as_ref().unwrap();
    let mut rcv_pos = Vector3::new(state.position.vector.x, state.position.vector.y, state.position.vector.z);
    rcv_pos += gneiss_core::tides::solid_earth_tides_ecef(r_obs.time, rcv_pos);
    let rcv_llh = gneiss_core::coords::ecef_to_llh(rcv_pos);
    
    r_obs.satellites.iter()
        .filter_map(|sat_obs| process_single_sat(engine, r_obs, sat_obs, rcv_pos, rcv_llh))
        .collect()
}

fn process_single_sat<'a>(
    engine: &ProcessingEngine, r_obs: &'a EpochObs, sat_obs: &'a gneiss_core::obs::SatObs,
    rcv_pos: Vector3<f64>, rcv_llh: Vector3<f64>
) -> Option<ProcessedSat<'a>> {
    let eph = engine.ephemerides.iter().find(|e| e.sat() == sat_obs.sat)?;
    let (f1, mut f2) = gneiss_core::signal::satellite_frequencies(sat_obs.sat, eph.freq_num());
    if f2 == 0.0 { f2 = f1; }

    let (mut p1, mut p2, mut cp1, mut cp2, osb, is_if) = get_obs_and_corrections(engine, sat_obs, r_obs.time, f1, f2);
    let tau_pr = p1.unwrap_or(0.0) / LIGHT_SPEED;
    let t_nom = gneiss_core::time::GpsTime::new(r_obs.time.week, r_obs.time.tow - tau_pr);

    let (t_tx, dt_s, raw_pos, raw_vel) = compute_sat_state(engine, eph, sat_obs.sat, t_nom)?;
    let (sat_pos, sat_vel) = crate::engine::ppp_math::apply_earth_rotation(raw_pos, raw_vel, rcv_pos);
    
    let dist = (sat_pos - rcv_pos).norm();
    let (az, el) = gneiss_core::coords::az_el(rcv_llh, rcv_pos, sat_pos);
    if el < 0.261799 { return None; }

    let (tropo_dry, map_wet) = crate::engine::ppp_math::compute_tropo_dry(rcv_llh, el, r_obs.time);
    let klobuchar = engine.klobuchar_params.unwrap_or_default();
    let iono_delay = if is_if { 0.0 } else { gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar(&klobuchar, rcv_llh, az, el, r_obs.time) };
    
    let pcv = compute_pcv(engine, sat_obs, r_obs.time, f1, f2, is_if, rcv_pos, sat_pos, &mut p2, &mut cp2);
    tracing::trace!("PPP sat={}, is_if={}, dist={:.1}, clk_m={:.3}, pcv={:.4}", sat_obs.sat, is_if, dist, dt_s * LIGHT_SPEED, pcv);

    let snr = sat_obs.get_snr(1).unwrap_or(45) as f64;
    let doppler = sat_obs.get_doppler(1).unwrap_or(0.0);
    
    Some(ProcessedSat {
        sat_obs, dt_sat_m: dt_s * LIGHT_SPEED, p1: p1?, p2, is_iono_free: is_if,
        cp1, cp2, osb_p1: osb.osb_p1, osb_p2: osb.osb_p2, osb_cp1: osb.osb_cp1, osb_cp2: osb.osb_cp2,
        los: (sat_pos - rcv_pos) / dist, dist: dist - pcv, el, snr, doppler,
        lam1: LIGHT_SPEED / f1, lam2: LIGHT_SPEED / f2, tropo_dry, map_wet, iono_delay,
        f1, f2, sat_pos_rot: sat_pos, sat_vel, sat_clock_drift: eph.position(t_tx).3,
        rcv_pos_ecef: rcv_pos, pcv_correction: pcv,
    })
}

fn get_obs_and_corrections(
    engine: &ProcessingEngine, sat_obs: &gneiss_core::obs::SatObs, time: gneiss_core::time::GpsTime, f1: f64, f2: f64
) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>, crate::engine::ppp_math::OsbCorrections, bool) {
    let get = |t, fb| sat_obs.observations.iter().find(|o| o.code.obs_type == t && o.code.signal.freq_band == fb).map(|o| (o.value, o.code));
    let f1_b = if sat_obs.sat.constellation == Constellation::Beidou { 2 } else { 1 };
    let f2_b = match sat_obs.sat.constellation { Constellation::Galileo | Constellation::Beidou => 7, _ => 2 };

    let osb = if let Some(sinex) = &engine.sinex_bias {
        crate::engine::ppp_math::apply_osb_corrections(sinex, sat_obs.sat, time, f1, f2, get(gneiss_core::obs::ObsType::Pseudorange, f1_b), get(gneiss_core::obs::ObsType::Pseudorange, f2_b), get(gneiss_core::obs::ObsType::CarrierPhase, f1_b), get(gneiss_core::obs::ObsType::CarrierPhase, f2_b))
    } else {
        crate::engine::ppp_math::OsbCorrections { p1: get(gneiss_core::obs::ObsType::Pseudorange, f1_b).map(|o| o.0), p2: get(gneiss_core::obs::ObsType::Pseudorange, f2_b).map(|o| o.0), cp1: get(gneiss_core::obs::ObsType::CarrierPhase, f1_b).map(|o| o.0), cp2: get(gneiss_core::obs::ObsType::CarrierPhase, f2_b).map(|o| o.0), osb_p1: 0.0, osb_p2: 0.0, osb_cp1: 0.0, osb_cp2: 0.0 }
    };

    let mut p1 = osb.p1; let mut p2 = osb.p2; let mut cp1 = osb.cp1; let mut cp2 = osb.cp2;
    if engine.sinex_bias.is_none() {
        if let Some(v) = p1.as_mut() { if let Some(&d) = engine.dcbs.get(&(sat_obs.sat, "P1C1".to_string())) { *v -= d * 1e-9 * LIGHT_SPEED; } }
        if let Some(v) = p2.as_mut() { if let Some(&d) = engine.dcbs.get(&(sat_obs.sat, "P2C2".to_string())) { *v -= d * 1e-9 * LIGHT_SPEED; } }
    }

    let mut is_if = false;
    let precise = !engine.sp3_epochs.is_empty() || engine.clk_data.is_some();
    if precise && !engine.config.uduc_ar {
        if let (Some(v1), Some(v2)) = (p1, p2) { p1 = Some(crate::engine::ppp_math::compute_iono_free(f1, f2, v1, v2)); is_if = true; }
        if let (Some(l1), Some(l2)) = (cp1, cp2) { cp1 = Some(crate::engine::ppp_math::compute_iono_free(f1, f2, l1 * LIGHT_SPEED / f1, l2 * LIGHT_SPEED / f2) / (LIGHT_SPEED / f1)); }
    }
    (p1, p2, cp1, cp2, osb, is_if)
}

fn compute_sat_state(
    engine: &ProcessingEngine, eph: &gneiss_core::ephemeris::EphemerisData, sat: gneiss_core::sat::SatelliteId, t_nom: gneiss_core::time::GpsTime
) -> Option<(gneiss_core::time::GpsTime, f64, Vector3<f64>, Vector3<f64>)> {
    let precise = !engine.sp3_epochs.is_empty() || engine.clk_data.is_some();
    let mut dt_s = engine.clk_data.as_ref().and_then(|c| c.get_clock_bias(sat, t_nom)).unwrap_or(0.0);
    let brdc_clk = eph.position(t_nom).2;
    crate::engine::ppp_math::check_clock_diff(sat, dt_s, brdc_clk);

    let mut clk_found = dt_s != 0.0;
    if !precise { dt_s = brdc_clk; clk_found = true; }
    
    if precise && !clk_found {
        if let Some((_, _, sp3_clk)) = crate::engine::ssr::get_precise_orbit(&engine.sp3_epochs, sat, t_nom, 10) {
            if !sp3_clk.is_nan() && sp3_clk != 0.0 { dt_s = sp3_clk; clk_found = true; }
        }
    }
    if precise && !clk_found { return None; }

    let t_tx = gneiss_core::time::GpsTime::new(t_nom.week, t_nom.tow - dt_s);
    let (brdc_pos, brdc_vel, _, _) = eph.position(t_tx);
    let mut sat_pos = brdc_pos; let mut sat_vel = brdc_vel;
    
    if precise {
        let dt_rel = -2.0 * brdc_pos.dot(&brdc_vel) / (LIGHT_SPEED * LIGHT_SPEED);
        dt_s += dt_rel;
        if let Some((sp3_p, sp3_v, _)) = crate::engine::ssr::get_precise_orbit(&engine.sp3_epochs, sat, t_tx, 10) {
            sat_pos = sp3_p; sat_vel = sp3_v;
            crate::engine::ppp_math::check_pos_diff(sat, sp3_p, brdc_pos);
        } else {
            return None;
        }
    }
    Some((t_tx, dt_s, sat_pos, sat_vel))
}

fn compute_pcv(
    engine: &ProcessingEngine, sat_obs: &gneiss_core::obs::SatObs, t: gneiss_core::time::GpsTime,
    f1: f64, f2: f64, is_if: bool, rcv_pos: Vector3<f64>, sat_pos: Vector3<f64>, p2: &mut Option<f64>, cp2: &mut Option<f64>
) -> f64 {
    let precise = !engine.sp3_epochs.is_empty() || engine.clk_data.is_some();
    if !precise || engine.antex.is_none() { return 0.0; }
    
    let db = engine.antex.as_ref().unwrap();
    let utc_time = chrono::Utc.with_ymd_and_hms(1980, 1, 6, 0, 0, 0).unwrap() + chrono::Duration::seconds(t.week as i64 * 604800 + t.tow as i64);
    let ant = match db.find_satellite(&sat_obs.sat.to_string(), utc_time) { Some(a) => a, None => return 0.0 };

    let code1 = match sat_obs.sat.constellation { Constellation::Glonass => "R01", Constellation::Galileo => "E01", Constellation::Beidou => "C02", _ => "G01" };
    let code2 = match sat_obs.sat.constellation { Constellation::Glonass => "R02", Constellation::Galileo => "E05", Constellation::Beidou => "C07", _ => "G02" };
    
    let get_pco = |c| ant.frequencies.get(c).map(|p| p.pco / 1000.0).unwrap_or_else(|| ant.frequencies.values().next().map(|p| p.pco / 1000.0).unwrap_or(Vector3::zeros()));
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
        if let Some(p) = p2.as_mut() { *p -= diff; }
        if let Some(cp) = cp2.as_mut() { *cp -= diff / (LIGHT_SPEED / f2); }
    }
    pcv1
}

fn update_phase_ambiguities(state: &mut RtkState, sats: &[ProcessedSat], t: gneiss_core::time::GpsTime) {
    for sat in sats {
        if let Some(cp1) = sat.cp1 {
            if cp1 == 0.0 { continue; }
            let wup = gneiss_core::windup::phase_windup(sat.sat_pos_rot, gneiss_core::sun::sun_position_ecef(t), sat.rcv_pos_ecef, *state.windup.get(&sat.sat_obs.sat).unwrap_or(&0.0));
            state.windup.insert(sat.sat_obs.sat, wup);
            
            let l_meas = if sat.is_iono_free && sat.cp2.is_some() { crate::combinations::iono_free((cp1 + wup) * sat.lam1, (sat.cp2.unwrap() + wup) * sat.lam2, sat.f1, sat.f2) } else { (cp1 + wup) * sat.lam1 };
            let prev = *state.locktimes.get(&(sat.sat_obs.sat, 1)).unwrap_or(&0);
            let (slip, new_lk) = crate::engine::ppp_math::detect_cycle_slip(sat.sat_obs, prev);
            state.locktimes.insert((sat.sat_obs.sat, 1), new_lk);
            
            if slip {
                for i in 0..4 { state.remove_ambiguity(sat.sat_obs.sat, i); }
            }
            
            let isb = match sat.sat_obs.sat.constellation { Constellation::Glonass => state.isb_glo, Constellation::Galileo => state.isb_gal, Constellation::Beidou => state.isb_bds, _ => 0.0 };
            let expected_base = sat.dist + state.rcv_clk_bias + isb - sat.dt_sat_m + sat.tropo_dry + state.zwd * sat.map_wet;
            
            if !sat.is_iono_free && sat.cp2.is_some() && sat.p2.is_some() {
                add_uduc_ambiguities(state, sat, cp1, wup, expected_base);
            } else {
                let exp = if sat.is_iono_free && sat.cp2.is_some() { expected_base } else { expected_base - sat.iono_delay };
                if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 0)) { state.add_ambiguity(sat.sat_obs.sat, 0, l_meas - exp, 10000.0); }
                state.last_observed.insert((sat.sat_obs.sat, 0), state.epoch_count as u32);
            }
        }
    }
}

fn add_uduc_ambiguities(state: &mut RtkState, sat: &ProcessedSat, cp1: f64, wup: f64, expected_base: f64) {
    let p2 = sat.p2.unwrap();
    let p1 = sat.p1;
    let gamma = (sat.f1 * sat.f1) / (sat.f2 * sat.f2);
    let mut i1_est = (p2 - p1) / (gamma - 1.0);
    if i1_est.is_nan() || i1_est.abs() > 100.0 { i1_est = 0.0; }
    
    let l1_meas = (cp1 + wup) * sat.lam1;
    let l2_meas = (sat.cp2.unwrap() + wup) * sat.lam2;

    if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 3)) { state.add_ambiguity(sat.sat_obs.sat, 3, i1_est, 100.0); }
    if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 1)) { state.add_ambiguity(sat.sat_obs.sat, 1, l1_meas - (expected_base - i1_est), 10000.0); }
    if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 2)) { state.add_ambiguity(sat.sat_obs.sat, 2, l2_meas - (expected_base - gamma * i1_est), 10000.0); }
    
    for i in 1..4 { state.last_observed.insert((sat.sat_obs.sat, i), state.epoch_count as u32); }
}
"""

start_idx = orig.find(start_str)
end_idx = orig.find(end_str)

if start_idx != -1 and end_idx != -1:
    new_orig = orig[:start_idx] + build_sats_new + orig[end_idx:]
    new_orig = new_orig.replace("crate::combinations::iono_free", "crate::engine::ppp_math::compute_iono_free")
    
    # We must remove OsbCorrections and apply_osb_corrections entirely since they're in ppp_math.rs
    rm_start = "pub struct OsbCorrections {"
    rm_end = "OsbCorrections { p1, p2, cp1, cp2, osb_p1, osb_p2, osb_cp1, osb_cp2 }\n}\n"
    r1 = new_orig.find(rm_start)
    r2 = new_orig.find(rm_end) + len(rm_end)
    if r1 != -1 and r2 != -1:
        new_orig = new_orig[:r1] + new_orig[r2:]
        
    new_orig = new_orig.replace("apply_osb_corrections(", "crate::engine::ppp_math::apply_osb_corrections(")
    new_orig = new_orig.replace("OsbCorrections {", "crate::engine::ppp_math::OsbCorrections {")
    with open("crates/gneiss-rtk/src/engine/ppp.rs", "w") as f:
        f.write(new_orig)
