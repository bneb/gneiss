use gneiss_core::obs::{EpochObs, SatObs};
use gneiss_core::coords::Coordinate;
use crate::filter::RtkState;
use crate::engine::{EngineError, ProcessingEngine};
use nalgebra::{Vector3, DMatrix, DVector};
use chrono::TimeZone;
use crate::engine::processed_sat::ProcessedSat;
use crate::factor_graph::{FactorGraphOptimizer, gnss_factors::{ErrorStatePseudorangeFactor, ErrorStateCarrierPhaseFactor}};

use gneiss_core::constants::{SPEED_OF_LIGHT_M_S, EARTH_ROTATION_RATE_RAD_S};

pub fn snr_scale(snr: i32) -> f64 {
    (10.0f64).powf((45.0 - snr as f64) / 10.0)
}

#[cfg_attr(test, mutants::skip)]
pub fn process_ppp_fg<'a>(engine: &'a mut ProcessingEngine, rover_obs: &EpochObs) -> Result<&'a RtkState, EngineError> {
    if !ensure_valid_state(engine, rover_obs) {
        return engine.process_spp(rover_obs);
    }
    let mut spp_cdt = None;
    let mut spp_pos = None;
    let epoch_count = engine.current_state.as_ref().unwrap().epoch_count;
    if epoch_count == 0 {
        if let Ok(s) = crate::spp::compute_spp(&rover_obs, &engine.ephemerides, engine.klobuchar.as_ref(), &crate::spp::SppConfig::default(), None) {
            tracing::info!("Initializing PPP state from SPP at {:?}", rover_obs.time);
            spp_cdt = Some(s.cdt);
            spp_pos = Some(s.position);
            let st = engine.current_state.as_mut().unwrap();
            st.position = s.position;
            st.rcv_clk_bias = s.cdt;
        } else {
            return Err(EngineEr
        println!("Clk correction: gps={:.2}, gal={:.2}, bds={:.2}, glo={:.2}", opt_delta[15], 
            if crate::filter::CORE_STATE_SIZE >= 21 { opt_delta[18] } else { 0.0 }, 
            if crate::filter::CORE_STATE_SIZE >= 21 { opt_delta[19] } else { 0.0 }, 
            if crate::filter::CORE_STATE_SIZE >= 21 { opt_delta[20] } else { 0.0 });
    }
    
    if state.epoch_count == 100 { println!("P_pos: {:.3}, P_amb: {:.3}", cov[(0,0)], cov[(crate::filter::CORE_STATE_SIZE, crate::filter::CORE_STATE_SIZE)]); }
    
    state.position.vector.x += opt_delta[0];
    state.position.vector.y += opt_delta[1];
    state.position.vector.z += opt_delta[2];
    state.velocity.x += opt_delta[3];
    state.velocity.y += opt_delta[4];
    state.velocity.z += opt_delta[5];
    state.attitude *= nalgebra::UnitQuaternion::from_scaled_axis(Vector3::new(opt_delta[6], opt_delta[7], opt_delta[8]));
    state.accel_bias.x += opt_delta[9];
    state.accel_bias.y += opt_delta[10];
    state.accel_bias.z += opt_delta[11];
    state.gyro_bias.x += opt_delta[12];
    state.gyro_bias.y += opt_delta[13];
    state.gyro_bias.z += opt_delta[14];
    state.rcv_clk_bias += opt_delta[15];
    state.rcv_clk_drift += opt_delta[16];
    state.zwd += opt_delta[17];
        if crate::filter::CORE_STATE_SIZE >= 21 {
        state.cdt_gal += opt_delta[18];
        state.cdt_bds += opt_delta[19];
        state.cdt_glo += opt_delta[20];
    }
    for i in 0..state.ambiguities.len() { state.ambiguities[i] += opt_delta[crate::filter::CORE_STATE_SIZE + i]; }
    
    state.covariance = cov;
    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
    state.time = rover_obs.time;
    state.position.epoch = rover_obs.time;
}
