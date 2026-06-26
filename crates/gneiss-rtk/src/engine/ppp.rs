use crate::engine::ppp_iekf::PppIteratedEkf;
use crate::engine::processed_sat::ProcessedSat;
use crate::engine::{EngineError, EngineMode, ProcessingEngine};
use crate::filter::RtkState;
use gneiss_core::obs::EpochObs;
use gneiss_core::sat::Constellation;
use nalgebra::Vector3;

use super::ppp_antenna::*;

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
    // SPP seed: compute position and use as soft prior in the IEKF.
    // Cold start (epoch 0): hard-reset position to SPP.
    // Subsequent epochs: inject SPP as a prior measurement with
    // variance that decreases as the filter converges. This allows
    // multi-epoch carrier-phase convergence while staying anchored.
    //
    // NOTE: Removing the SPP anchor after convergence was attempted
    // (Phase 2) but caused severe divergence (9m→44m Hz, confirming
    // POST_MORTEM hypothesis #8). Multi-epoch convergence requires a
    // sliding-window factor graph — not just disabling the prior.
    let mut position_prior: Option<(Vector3<f64>, f64)> = None;
    let mut spp_pos_for_recovery: Option<gneiss_core::coords::Coordinate> = None;
    if let Ok(spp) = crate::spp::compute_spp(
        rover_obs,
        &engine.ephemerides,
        engine.klobuchar_params.as_ref(),
        &crate::spp::SppConfig::default(),
        None,
    ) {
        spp_pos_for_recovery = Some(spp.position);
        engine.last_spp_position = Some(spp.position);
        let is_cold_start = state.epoch_count < 2;
        if is_cold_start {
            state.position = spp.position;
            state.rcv_clk_bias = spp.cdt;
            // Do NOT reset position covariance here. The initial predict_state
            // call spans a huge dt (GPS epoch 0 → rover time), inflating the
            // position variance to the clamp ceiling. This large variance is
            // beneficial: it tells the IEKF to trust measurements over the
            // position prior, allowing convergence from a potentially poor SPP
            // seed. The covariance shrinks naturally as measurements are
            // assimilated over subsequent epochs.
        } else {
            // Prior variance clamped to [9, 100] m².
            // 9 m² floor (σ=3m): prevents the prior from dominating
            // carrier-phase after rapid convergence, which would lock
            // the filter to a potentially wrong SPP seed (critical at
            // equatorial stations where SPP can be off by 30-70m).
            // 100 m² cap (σ=10m): prevents SPP outliers from
            // destabilising the filter during early convergence.
            let pos_cov = state.covariance[(0, 0)]
                .min(state.covariance[(1, 1)])
                .min(state.covariance[(2, 2)]);
            let prior_var = pos_cov.clamp(9.0, 100.0);
            position_prior = Some((spp.position.vector, prior_var));
        }
    }

    let _has_precise = !engine.sp3_epochs.is_empty() || engine.clk_data.is_some();

    let sats = build_sats(engine, rover_obs);
    if sats.is_empty() {
        return Err(EngineError::InsufficientSatellites);
    }
    let state = engine.current_state.as_mut().unwrap();
    update_phase_ambiguities(state, &sats, rover_obs.time);
    state.prune_stale_ambiguities(state.epoch_count as u32, 10);

    // Dispatch to appropriate solver based on engine mode
    let solve_result = if engine.config.mode == EngineMode::PppMultiEpoch {
        let opt = engine.ppp_multi_epoch_opt.take();
        let mut solver = opt.unwrap_or_else(|| {
            crate::engine::ppp_multi_epoch::MultiEpochOptimizer::new(2)
        });
        let result = solver.solve(state, &sats, position_prior);
        engine.ppp_multi_epoch_opt = Some(solver);
        result
    } else {
        let mut opt = engine.ppp_factor_opt.take();
        let result = if let Some(ref mut solver) = opt {
            solver.solve(state, &sats, position_prior)
        } else {
            PppIteratedEkf::new()
                .with_iono_model(engine.config.iono_model)
                .with_lambda_min_ratio(engine.config.lambda_min_ratio)
                .solve(state, &sats, position_prior)
        };
        engine.ppp_factor_opt = opt;
        result
    };
    state.epoch_count = state.epoch_count.saturating_add(1);

    // Innovation gate: if the IEKF solution has diverged catastrophically
    // from the predicted state, reject it and keep the prediction. Uses a
    // high threshold (500 m or 50× sqrt(trace_P), whichever is larger) to
    // catch true divergence without rejecting legitimate position updates
    // from satellite geometry changes.
    if let Some(ref pred_pos) = state.predicted_position {
        let innov = state.position.vector - pred_pos.vector;
        let p00 = state.covariance[(0, 0)];
        let p11 = state.covariance[(1, 1)];
        let p22 = state.covariance[(2, 2)];
        let trace_p = p00 + p11 + p22;
        if trace_p > 0.0 {
            let threshold = (50.0_f64 * trace_p.sqrt()).max(500.0);
            if innov.norm() > threshold {
                tracing::warn!(
                    "IEKF position divergence: |innov|={:.0} m > threshold={:.0} m -- keeping predicted state",
                    innov.norm(),
                    threshold
                );
                state.position = *pred_pos;
                if let Some(pred_vel) = state.predicted_velocity {
                    state.velocity = pred_vel;
                } else {
                    // Velocity estimate is likely corrupted — zero it
                    // and force re-estimation from measurements.
                    state.velocity = nalgebra::Vector3::zeros();
                }
                // Inflate position AND velocity covariance to force
                // the next IEKF solve to re-estimate both from
                // measurements rather than trusting the prediction.
                for i in 0..3 {
                    state.covariance[(i, i)] = state.covariance[(i, i)].max(100.0);
                }
                for i in 3..6 {
                    state.covariance[(i, i)] = state.covariance[(i, i)].max(100.0);
                }
            }
        }
    }

    // Always push to history — predict_state() was already called at the
    // top of process_ppp, so the propagated state is valid even when solve
    // returns InsufficientSatellites.  The smoother bridges the gap.
    let final_state = engine.current_state.as_ref().unwrap().clone();
    engine.state_history.push(final_state);
    engine.obs_history.push((rover_obs.clone(), None));

    // Check for covariance explosion after the IEKF solve. At equatorial
    // stations (e.g. NKLG), frequent cycle slips from ionospheric
    // scintillation multiply position/velocity covariance by 4x per slip.
    // Without this guard, covariance grows unchecked past 1e13, causing
    // total measurement rejection and permanent divergence.
    // The RTK/SPP pipelines already call check_covariance_divergence;
    // PPP was missing it — this was the root cause of NKLG's 1600
    // coasting events.
    if let Some(ref mut state) = engine.current_state {
        ProcessingEngine::check_covariance_divergence(
            state,
            spp_pos_for_recovery,
            None,
            false,
        );
    }

    // Attempt INS alignment for tightly-coupled PPP-INS modes.
    // This must run after the IEKF solve so the state velocity is
    // available for kinematic or static alignment.
    if engine.config.mode == EngineMode::PppIns {
        engine.attempt_kinematic_alignment();
    }

    solve_result?;
    Ok(engine.current_state.as_ref().unwrap())
}

pub(crate) fn valid_pos(engine: &ProcessingEngine) -> bool {
    if let Some(state) = &engine.current_state {
        state.position.vector.norm().is_normal() && state.position.vector.norm() >= 1000.0
    } else {
        false
    }
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
        let prev = *state.locktimes.get(&(sat.sat_obs.sat, 1)).unwrap_or(&0);
        let mut gf_prev = state.gf_prev.get(&sat.sat_obs.sat).copied();
        let mut mw_prev = state.mw_prev.get(&sat.sat_obs.sat).copied();
        let (slip, new_lk) = crate::engine::ppp_math::detect_slip_combined(
            sat.sat_obs,
            prev as u32,
            sat.cp1,
            sat.lam1,
            sat.cp2,
            sat.lam2,
            Some(sat.p1),
            sat.p2,
            &mut gf_prev,
            &mut mw_prev,
        );
        state.locktimes.insert((sat.sat_obs.sat, 1), new_lk as u16);
        if let Some(v) = gf_prev {
            state.gf_prev.insert(sat.sat_obs.sat, v);
        }
        if let Some(v) = mw_prev {
            state.mw_prev.insert(sat.sat_obs.sat, v);
        }
        if slip {
            for i in 0..4 {
                state.remove_ambiguity(sat.sat_obs.sat, i);
            }
            // Inflate position and velocity covariance after losing phase
            // constraints.  Over-confidence in the current coordinate estimate
            // prevents re-convergence on new phase observations.
            // Multiply position (0..3) and velocity (3..6) diagonal elements
            // by 4, but clamp at 1e4 m² (σ=100m) to prevent exponential growth
            // from frequent ionospheric scintillation at equatorial stations.
            // Without the cap, N slips produce 4^N growth (e.g., 10 slips
            // = 4^10 = 1,048,576×), causing total measurement rejection.
            for i in 0..6 {
                state.covariance[(i, i)] = (state.covariance[(i, i)] * 4.0).min(5e4);
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
        // Compute Melbourne-Wübbena widelane for AR seeding.
        // Uses RAW observables from the RINEX file (not iono-free combined)
        // because MW requires single-frequency measurements.
        // The geometric range cancels in the MW combination, so any
        // common-mode errors (clock, tropo) are eliminated.
        let raw_l1 = sat.sat_obs.get_observable_phase(1);
        let raw_l2 = sat.sat_obs.get_observable_phase(2);
        let raw_p1 = sat.sat_obs.get_observable(1);
        let raw_p2 = sat.sat_obs.get_observable(2);
        // BUGFIX: In IF mode, push_cp_measurement() looks up band-0 ambiguity
        // via find_ambiguity_index().  If we create UDUC-style bands 1/2/3
        // (which the old code always did when raw L1/L2 exist), the CP
        // measurement is silently dropped — the IEKF degrades to PR-only.
        // Check is_iono_free FIRST to ensure band-0 exists in IF mode.
        if sat.is_iono_free {
            let l_meas = (cp1 - wup) * sat.lam1;
            let exp = expected_base;
            if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 0)) {
                state.add_ambiguity(sat.sat_obs.sat, 0, l_meas - exp, 10000.0);
            }
            state
                .last_observed
                .insert((sat.sat_obs.sat, 0), state.epoch_count as u32);
        } else if let (Some(l1), Some(l2), Some(p1), Some(p2)) =
            (raw_l1, raw_l2, raw_p1, raw_p2)
        {
            let l1_m = (l1 - wup) * sat.lam1;
            let l2_m = (l2 - wup) * sat.lam2;
            let geo = sat.dist; // geometric range from ProcessedSat
            let l1_res = l1_m - geo;
            let l2_res = l2_m - geo;
            let p1_res = p1 - geo;
            let p2_res = p2 - geo;
            let mw_m = (sat.f1 * l1_res - sat.f2 * l2_res) / (sat.f1 - sat.f2)
                - (sat.f1 * p1_res + sat.f2 * p2_res) / (sat.f1 + sat.f2);
            let mw_cycles = mw_m * (sat.f1 - sat.f2) / LIGHT_SPEED;
            state.update_mw(sat.sat_obs.sat, mw_cycles);
            add_uduc_ambiguities(state, sat, cp1, wup, expected_base);
        } else {
            // Fallback: non-IF, no raw L2/P2 — use band-0
            // (e.g., single-frequency receivers)
            let l_meas = (cp1 - wup) * sat.lam1;
            let exp = expected_base - sat.iono_delay;
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

    let l1_meas = (cp1 - wup) * sat.lam1;
    let l2_meas = (sat.cp2.unwrap() - wup) * sat.lam2;

    // Use MW widelane to reduce initial ambiguity variance when available.
    // Require 50+ samples: the EMA first-sample weight drops to ~4% at N=50,
    // giving ~0.06 cycle WL precision — tight enough for safe LAMBDA.
    let mw_confident = state
        .mw_sd_counts
        .get(&sat.sat_obs.sat)
        .copied()
        .unwrap_or(0)
        > 50;
    let init_var = if mw_confident { 0.04 } else { 10000.0 }; // 0.2 cycle or 100m std
    if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 3)) {
        state.add_ambiguity(sat.sat_obs.sat, 3, i1_est, 100.0);
    }
    if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 1)) {
        state.add_ambiguity(
            sat.sat_obs.sat,
            1,
            l1_meas - (expected_base - i1_est),
            init_var,
        );
    }
    if !state.ambiguity_keys.contains(&(sat.sat_obs.sat, 2)) {
        state.add_ambiguity(
            sat.sat_obs.sat,
            2,
            l2_meas - (expected_base - gamma * i1_est),
            init_var,
        );
    }

    for i in 1..4 {
        state
            .last_observed
            .insert((sat.sat_obs.sat, i), state.epoch_count as u32);
    }
}

#[cfg(test)]
#[path = "ppp_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "ppp_antenna_tests.rs"]
mod antenna_tests;
