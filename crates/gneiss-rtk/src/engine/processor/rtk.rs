use super::ProcessingEngine;
use crate::engine::matcher::match_observations;
use crate::engine::updater_math::{CouplingStrategy, LooseCoupling, TightCoupling};
use crate::engine::{EngineConfig, EngineError, EngineMode};
use crate::filter::RtkState;
use gneiss_core::coords::{Coordinate, Datum, Frame};
use gneiss_core::obs::EpochObs;
use gneiss_core::sat::SatelliteId;
use nalgebra::{DMatrix, DVector, Vector3};

/// Apply adaptive R scaling based on innovation history.
/// Updates the tracker with current innovations and inflates R diagonals
/// for measurements with historically large normalized innovations.
fn apply_adaptive_r_scaling(
    tracker: &mut crate::engine::adaptive::InnovationTracker,
    state: &RtkState,
    z: &DVector<f64>,
    h: &DMatrix<f64>,
    r: &mut DMatrix<f64>,
    meas_types: &[(SatelliteId, u8, f64)],
    matched_obs: &[(crate::filter::DdObservation, crate::filter::DdObservation)],
) {
    let state_size = state.covariance.nrows();
    for i in 0..z.nrows() {
        let h_row_vec: Vec<f64> = (0..state_size).map(|c| h[(i, c)]).collect();
        let h_row_mat = DMatrix::from_row_slice(1, state_size, &h_row_vec);
        let s_ii = (&h_row_mat * &state.covariance * h_row_mat.transpose())[(0, 0)] + r[(i, i)];
        let (sat, mtype, _) = meas_types[i];
        // Determine freq band from measurement type: 0=PR_L1, 1=CP_L1, 2=CP_L2, 3=Dop
        let freq = match mtype {
            2 => 2,
            _ => 1,
        };

        let mut snr = None;
        if let Some((rov_obs, _)) = matched_obs.iter().find(|(rov, _)| rov.sat == sat) {
            snr = Some(rov_obs.snr);
        }

        let scale = tracker.update_and_scale(sat, freq, z[i], s_ii, snr);
        r[(i, i)] *= scale;
    }
}

impl ProcessingEngine {
    fn init_spp_state(
        &mut self,
        rover_obs: &EpochObs,
    ) -> Result<Option<crate::spp::SppState>, EngineError> {
        let spp_res = crate::spp::compute_spp(
            rover_obs,
            &self.ephemerides,
            self.klobuchar_params.as_ref(),
            &crate::spp::SppConfig::default(),
            None,
        )
        .ok();

        if spp_res.is_none() {
            tracing::warn!("compute_spp failed in process_rtk!");
        }

        if self.current_state.is_none() {
            if let Some(spp) = &spp_res {
                self.current_state = Some(RtkState::new(rover_obs.time, spp.position, 100.0));
            } else {
                return Err(EngineError::InitialSppFailed);
            }
        }
        Ok(spp_res)
    }

    fn evaluate_gnss_only_coasting(
        config: &EngineConfig,
        state: &mut RtkState,
        spp_pos: Option<Coordinate>,
        spp_state_ref: Option<&crate::spp::SppState>,
        had_imu_data: bool,
    ) {
        let use_gnss_only_seed = matches!(config.mode, EngineMode::Rtk | EngineMode::Ppp)
            || (matches!(
                config.mode,
                EngineMode::RtkIns
                    | EngineMode::PppIns
                    | EngineMode::SppIns
                    | EngineMode::RtkInsLooselyCoupled
                    | EngineMode::SppInsLooselyCoupled
                    | EngineMode::PppInsLooselyCoupled
            ) && !had_imu_data);

        if use_gnss_only_seed {
            tracing::debug!("epoch_count = {}", state.epoch_count);
            let need_spp_reset = if spp_pos.is_some() {
                state.epoch_count < 3
            } else {
                false
            };

            if need_spp_reset {
                if let Some(pos) = spp_pos {
                    tracing::warn!("Resetting EKF state to SPP due to divergence/startup.");
                    state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp());
                }
            }
        }
    }

    fn filter_valid_base<'a>(
        config: &EngineConfig,
        rover_obs: &EpochObs,
        base_obs: Option<&'a EpochObs>,
    ) -> Option<&'a EpochObs> {
        base_obs.filter(|b| {
            let age = (rover_obs.time.tow - b.time.tow).abs();
            if age > config.max_base_age_s {
                tracing::trace!(
                    "Base observation rejected due to Age of Differential ({:.1}s > {:.1}s)",
                    age,
                    config.max_base_age_s
                );
                false
            } else {
                true
            }
        })
    }

    fn perform_spp_fallback_update(config: &EngineConfig, state: &mut RtkState, pos: Coordinate) {
        tracing::debug!("No valid base data. Falling back to SPP.");
        tracing::debug!("RTK base missing or stale. Falling back to SPP update.");
        let z_diff = pos.vector - state.position.vector;
        let z_vec = nalgebra::DVector::from_column_slice(z_diff.as_slice());

        let mut h_mat = nalgebra::DMatrix::zeros(3, state.covariance.ncols());
        h_mat.view_mut((0, 0), (3, 3)).fill_diagonal(1.0);

        let mut r_mat = nalgebra::DMatrix::zeros(3, 3);
        r_mat.fill_diagonal(900.0);

        if let Err(e) = crate::engine::updater::update::<LooseCoupling>(
            state,
            &z_vec,
            &h_mat,
            &r_mat,
            config.spp_consistency_threshold_m,
            None,
            &config.tuning,
        ) {
            tracing::debug!("SPP Fallback update failed: {:?}", e);
        }
    }

    fn manage_imu_buffer_on_start(&mut self) {
        // Do not clear the buffer, so that it can be pushed to imu_history and used for alignment
    }

    fn update_state_time(&mut self, time: gneiss_core::time::GpsTime) -> Result<(), EngineError> {
        let state = self
            .current_state
            .as_mut()
            .ok_or(EngineError::StateDisappeared)?;
        state.is_reset = false;
        state.time = time;
        state.position.epoch = time;
        Ok(())
    }

    fn get_base_coord(
        config: &EngineConfig,
        rover_obs: &EpochObs,
    ) -> Result<Coordinate, EngineError> {
        let mut base_coord = if let Some(base_pos_arr) = config.base_position {
            Coordinate::new(
                Vector3::new(base_pos_arr[0], base_pos_arr[1], base_pos_arr[2]),
                Datum::WGS84,
                Frame::ECEF,
                rover_obs.time,
            )
        } else {
            return Err(EngineError::MissingBasePosition);
        };

        if let Some(helmert) = &config.base_datum_transform {
            let obs_epoch = rover_obs.time.to_fractional_year();
            let transformed_vec = helmert.transform(base_coord.vector, obs_epoch);
            base_coord =
                Coordinate::new(transformed_vec, Datum::WGS84, Frame::ECEF, rover_obs.time);
        }
        Ok(base_coord)
    }

    fn evaluate_gnn(
        gnn_raim: &Option<crate::engine::ml::gnn_raim::GnnRaimModel>,
        ephemerides: &[gneiss_core::ephemeris::Ephemeris],
        rover_obs: &EpochObs,
        matched_obs: &[(crate::filter::DdObservation, crate::filter::DdObservation)],
        state: &RtkState,
    ) -> std::collections::HashMap<SatelliteId, f64> {
        if let Some(gnn) = gnn_raim {
            let pos_apc = state.position.vector;
            let rov_llh = gneiss_core::coords::ecef_to_llh(pos_apc);
            crate::engine::ml::gnn_raim::evaluate_gnn_raim(
                gnn,
                matched_obs,
                rov_llh,
                pos_apc,
                ephemerides,
                rover_obs.time,
            )
        } else {
            std::collections::HashMap::new()
        }
    }

    fn apply_observations(
        &mut self,
        rover_obs: &EpochObs,
        base_obs: Option<&EpochObs>,
        spp_pos: Option<Coordinate>,
        spp_state_ref: Option<&crate::spp::SppState>,
    ) -> Result<(), EngineError> {
        let valid_base = Self::filter_valid_base(&self.config, rover_obs, base_obs);
        let base_coord_res = Self::get_base_coord(&self.config, rover_obs);
        let state = self
            .current_state
            .as_mut()
            .ok_or(EngineError::StateDisappeared)?;

        if let Some(base) = valid_base {
            state.epoch_count += 1;
            let base_coord = base_coord_res?;
            let matched_obs = match_observations(rover_obs, base, &self.ephemerides);

            if matched_obs.len() >= 5 {
                let gnn_variances = Self::evaluate_gnn(
                    &self.gnn_raim,
                    &self.ephemerides,
                    rover_obs,
                    &matched_obs,
                    state,
                );
                let ctx = RtkUpdateContext {
                    config: &self.config,
                    ephemerides: &self.ephemerides,
                    imu_history: &self.imu_history,
                    rover_obs,
                    base_obs: base,
                    matched_obs: &matched_obs,
                    base_coord: &base_coord,
                    spp_pos,
                    spp_state_ref,
                    gnn_variances,
                };
                process_rtk_update::<TightCoupling>(state, &mut self.innovation_tracker, &ctx);
            } else {
                state.consecutive_rejections += 1;
            }
        } else if let Some(pos) = spp_pos {
            Self::perform_spp_fallback_update(&self.config, state, pos);
        }
        Ok(())
    }

    pub fn process_rtk(
        &mut self,
        rover_obs: &EpochObs,
        base_obs: Option<&EpochObs>,
    ) -> Result<&RtkState, EngineError> {
        let spp_res = self.init_spp_state(rover_obs)?;
        let spp_pos = spp_res.as_ref().map(|s| s.position);
        let spp_state_ref = spp_res.as_ref();

        let mut rover_smoothed = rover_obs.clone();
        self.hatch_filter.smooth_epoch(&mut rover_smoothed);

        let dt = rover_obs.time.tow
            - self
                .current_state
                .as_ref()
                .ok_or(EngineError::StateDisappeared)?
                .time
                .tow;
        self.manage_imu_buffer_on_start();
        let had_imu_data = !self.imu_buffer.is_empty();

        self.predict_state(dt);
        self.update_state_time(rover_obs.time)?;

        let state = self.current_state.as_mut().expect("current_state is Some after ok_or early return");
        Self::check_covariance_divergence(
            state,
            spp_pos,
            spp_state_ref,
            !self.config.mode.is_ppp(),
        );
        Self::evaluate_gnss_only_coasting(
            &self.config,
            state,
            spp_pos,
            spp_state_ref,
            had_imu_data,
        );

        self.apply_observations(&rover_smoothed, base_obs, spp_pos, spp_state_ref)?;

        let state = self.current_state.as_mut().expect("current_state is Some after apply_observations");
        Self::apply_nhc_updates(&self.config, &self.imu_history, state);

        self.attempt_kinematic_alignment();

        if let Some(state) = &self.current_state {
            self.state_history.push(RtkState::clone(state));
        }
        self.obs_history
            .push((rover_obs.clone(), base_obs.cloned()));
        self.current_state
            .as_ref()
            .ok_or(EngineError::StateDisappeared)
    }
}

/// RTK measurement update — extracted as a free function to avoid borrow conflicts.
fn build_measurement_environment<'a>(
    config: &'a EngineConfig,
    imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
    state: &RtkState,
    ephemerides: &'a [gneiss_core::ephemeris::Ephemeris],
    base_coord: &'a Coordinate,
    base_time: gneiss_core::time::GpsTime,
) -> crate::engine::measurement::MeasurementEnvironment<'a> {
    let omega_ib_b = if let Some(imu_buf) = imu_history.last() {
        if let Some(last_imu) = imu_buf.last() {
            last_imu.gyro - state.gyro_bias
        } else {
            nalgebra::Vector3::zeros()
        }
    } else {
        nalgebra::Vector3::zeros()
    };
    let omega_ie_e =
        nalgebra::Vector3::new(0.0, 0.0, gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S);
    let r_e_b = state.attitude.to_rotation_matrix().transpose();
    let lever_arm = if state.ins_aligned {
        Vector3::from_column_slice(&config.imu_to_antenna_lever_arm)
    } else {
        Vector3::zeros()
    };

    crate::engine::measurement::MeasurementEnvironment {
        ephemerides,
        base_coord,
        base_time,
        lever_arm,
        omega_b: omega_ib_b - r_e_b * omega_ie_e,
        tuning: &config.tuning,
        gnn_variances: std::collections::HashMap::new(),
    }
}

fn handle_ekf_rejection(
    state: &mut RtkState,
    config: &EngineConfig,
    spp_pos: Option<Coordinate>,
    spp_state_ref: Option<&crate::spp::SppState>,
    reason: &str,
) {
    state.consecutive_rejections += 1;
    tracing::warn!(
        "GNSS EKF rejected for {} epochs ({}).",
        state.consecutive_rejections,
        reason
    );
    if state.ins_aligned {
        let inflate_factor = 1.0 + (state.consecutive_rejections as f64 * 0.02).min(0.5);
        for i in 0..6 {
            state.covariance[(i, i)] *= inflate_factor;
        }
        for i in 0..3 {
            state.covariance[(i, i)] += 1.0;
        }
        for i in 3..6 {
            state.covariance[(i, i)] += 0.1;
        }
    }
    let pos_var = state.covariance[(0, 0)] + state.covariance[(1, 1)] + state.covariance[(2, 2)];
    if state.ins_aligned && pos_var > 10000.0 {
        tracing::warn!(
            "Extreme divergence detected (pos_var {:.2}): resetting EKF to SPP fallback.",
            pos_var
        );
        if let Some(pos) = spp_pos {
            state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp());
            state.consecutive_rejections = 0;
        }
    } else if !state.ins_aligned && state.consecutive_rejections >= 3 {
        tracing::warn!(
            "Loosely coupled GNSS EKF rejected for 3 epochs: resetting to SPP fallback."
        );
        if let Some(pos) = spp_pos {
            state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp());
            state.consecutive_rejections = 0;
        }
    }
}

fn handle_ekf_acceptance(
    state: &mut RtkState,
    config: &EngineConfig,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    spp_pos: Option<Coordinate>,
    spp_state_ref: Option<&crate::spp::SppState>,
) {
    let pos_var = state.covariance[(0, 0)] + state.covariance[(1, 1)] + state.covariance[(2, 2)];
    if pos_var > 10000.0 {
        tracing::warn!(
            "Position variance extremely large ({:.2}): resetting EKF to SPP fallback.",
            pos_var
        );
        if let Some(pos) = spp_pos {
            state.reset_to_spp(pos, spp_state_ref, !config.mode.is_ppp());
        }
    }
    state.consecutive_rejections = 0;
    if let Ok(res) = state.resolve_ambiguities(ephemerides, config) {
        let fixed_state = res.fixed_state;
        tracing::debug!("Integer ambiguities resolved: {} sats", fixed_state.ambiguities.len());
        // Enable tight ambiguity process noise for fixed ambiguities.
        // Don't replace position/covariance — the AR position correction
        // may not be more accurate than the float solution.
        state.is_fixed = true;
        state.fixed_state = Some(Box::new(fixed_state));
    } else {
        state.is_fixed = false;
        state.fixed_state = None;
    }
}

pub struct RtkUpdateContext<'a> {
    pub config: &'a EngineConfig,
    pub ephemerides: &'a [gneiss_core::ephemeris::Ephemeris],
    pub imu_history: &'a [Vec<gneiss_core::imu::ImuMeasurement>],
    pub rover_obs: &'a EpochObs,
    pub base_obs: &'a EpochObs,
    pub matched_obs: &'a [(crate::filter::DdObservation, crate::filter::DdObservation)],
    pub base_coord: &'a Coordinate,
    pub spp_pos: Option<Coordinate>,
    pub spp_state_ref: Option<&'a crate::spp::SppState>,
    pub gnn_variances: std::collections::HashMap<gneiss_core::sat::SatelliteId, f64>,
}

fn process_rtk_update<C: CouplingStrategy>(
    state: &mut RtkState,
    tracker: &mut crate::engine::adaptive::InnovationTracker,
    ctx: &RtkUpdateContext,
) {
    crate::engine::ambiguity::manage_ambiguities_and_slips(
        state,
        ctx.config,
        ctx.matched_obs,
        ctx.ephemerides,
        ctx.base_coord,
        ctx.rover_obs.time,
        ctx.base_obs.time,
    );
    update_last_observed(state, ctx.matched_obs);

    let mut env = build_measurement_environment(
        ctx.config,
        ctx.imu_history,
        state,
        ctx.ephemerides,
        ctx.base_coord,
        ctx.base_obs.time,
    );
    env.gnn_variances = ctx.gnn_variances.clone();

    if let Some(mut m) = crate::engine::measurement::build_measurement_model(
        state,
        ctx.matched_obs,
        &env,
        ctx.config.chi_square_pr_threshold,
        ctx.config.chi_square_cp_threshold,
    ) {
        apply_adaptive_r_scaling(tracker, state, &m.z, &m.h, &mut m.r, &m.mt, ctx.matched_obs);
        execute_ekf_update::<C>(state, ctx, &m);
    } else {
        handle_ekf_rejection(
            state,
            ctx.config,
            ctx.spp_pos,
            ctx.spp_state_ref,
            "not enough measurements",
        );
    }
    state.prune_stale_ambiguities(state.epoch_count as u32, 10);
}

fn update_last_observed(
    state: &mut RtkState,
    matched_obs: &[(crate::filter::DdObservation, crate::filter::DdObservation)],
) {
    let current_epoch = state.epoch_count as u32;
    for (r_obs, _) in matched_obs {
        if r_obs.cp_l1.is_some() {
            state.last_observed.insert((r_obs.sat, 1), current_epoch);
        }
        if r_obs.cp_l2.is_some() {
            state.last_observed.insert((r_obs.sat, 2), current_epoch);
        }
    }
}

fn execute_ekf_update<C: CouplingStrategy>(
    state: &mut RtkState,
    ctx: &RtkUpdateContext,
    m: &crate::engine::measurement::EkfMeasurementMatrices,
) {
    let pr_thresh = ctx.config.chi_square_pr_threshold;
    let type_stripped: Vec<_> = m.mt.iter().map(|&(s, t, _)| (s, t)).collect();

    match crate::engine::updater::update::<C>(
        state,
        &m.z,
        &m.h,
        &m.r,
        pr_thresh,
        Some(&type_stripped),
        &ctx.config.tuning,
    ) {
        Err(_) => handle_ekf_rejection(
            state,
            ctx.config,
            ctx.spp_pos,
            ctx.spp_state_ref,
            "failed chi-square",
        ),
        Ok((valid_indices, dx)) => {
            export_gnn_dataset(state, ctx, m, &valid_indices, &dx);
            handle_ekf_acceptance(
                state,
                ctx.config,
                ctx.ephemerides,
                ctx.spp_pos,
                ctx.spp_state_ref,
            );
        }
    }
}

fn export_gnn_dataset(
    state: &RtkState,
    ctx: &RtkUpdateContext,
    m: &crate::engine::measurement::EkfMeasurementMatrices,
    valid_indices: &[usize],
    dx: &DVector<f64>,
) {
    if let Some(path) = &ctx.config.export_gnn_dataset_path {
        let rov_llh = gneiss_core::coords::ecef_to_llh(state.position.vector);
        let pos_apc = state.position.vector
            + state.attitude * nalgebra::Vector3::from(ctx.config.imu_to_antenna_lever_arm);

        crate::engine::ml::dataset::export_epoch_to_csv(
            path,
            state.epoch_count as u32,
            ctx.matched_obs,
            rov_llh,
            pos_apc,
            ctx.ephemerides,
            ctx.rover_obs.time,
            &m.z,
            &m.h,
            valid_indices,
            &m.mt,
            dx,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{EngineConfig, EngineError, EngineMode};
    use crate::filter::{CORE_STATE_SIZE, DdObservation};
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::obs::EpochObs;
    use gneiss_core::sat::Constellation;
    use gneiss_core::time::GpsTime;
    use gneiss_geodesy::helmert::HelmertParams;
    use nalgebra::{DMatrix, DVector, Vector3};

    fn make_test_sat(prn: u8) -> SatelliteId {
        SatelliteId {
            constellation: Constellation::Gps,
            prn,
        }
    }

    fn make_test_epoch(time: GpsTime) -> EpochObs {
        EpochObs {
            time,
            satellites: Vec::new(),
        }
    }

    fn make_dd_obs(sat: SatelliteId, cp_l1: Option<f64>, cp_l2: Option<f64>) -> DdObservation {
        DdObservation {
            sat,
            pr_l1: 0.0,
            pr_l2: None,
            cp_l1,
            cp_l2,
            doppler: 0.0,
            snr: 45.0,
            locktime: Some(100),
        }
    }

    #[test]
    fn test_filter_valid_base_accepts_recent() {
        let mut config = EngineConfig::default();
        config.max_base_age_s = 5.0;
        let rover_obs = make_test_epoch(GpsTime::new(0, 100.0));
        let base_obs = make_test_epoch(GpsTime::new(0, 101.0));
        let result =
            ProcessingEngine::filter_valid_base(&config, &rover_obs, Some(&base_obs));
        assert!(result.is_some());
    }

    #[test]
    fn test_filter_valid_base_rejects_old() {
        let mut config = EngineConfig::default();
        config.max_base_age_s = 5.0;
        let rover_obs = make_test_epoch(GpsTime::new(0, 100.0));
        let base_obs = make_test_epoch(GpsTime::new(0, 120.0));
        let result =
            ProcessingEngine::filter_valid_base(&config, &rover_obs, Some(&base_obs));
        assert!(result.is_none());
    }

    #[test]
    fn test_filter_valid_base_no_base() {
        let config = EngineConfig::default();
        let rover_obs = make_test_epoch(GpsTime::new(0, 100.0));
        let result = ProcessingEngine::filter_valid_base(&config, &rover_obs, None);
        assert!(result.is_none());
    }

    #[test]
    fn test_get_base_coord_ok() {
        let mut config = EngineConfig::default();
        config.base_position = Some([100.0, 200.0, 300.0]);
        let rover_obs = make_test_epoch(GpsTime::new(0, 100.0));
        let coord = ProcessingEngine::get_base_coord(&config, &rover_obs).unwrap();
        assert!((coord.vector.x - 100.0).abs() < 1e-6);
        assert!((coord.vector.y - 200.0).abs() < 1e-6);
        assert!((coord.vector.z - 300.0).abs() < 1e-6);
        assert_eq!(coord.datum, Datum::WGS84);
    }

    #[test]
    fn test_get_base_coord_missing() {
        let config = EngineConfig::default();
        let rover_obs = make_test_epoch(GpsTime::new(0, 100.0));
        let err = ProcessingEngine::get_base_coord(&config, &rover_obs).unwrap_err();
        assert!(matches!(err, EngineError::MissingBasePosition));
    }

    #[test]
    fn test_update_last_observed_inserts_all() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 5;

        let sat1 = make_test_sat(1);
        let sat2 = make_test_sat(2);
        let matched = vec![
            (make_dd_obs(sat1, Some(100.0), Some(200.0)), make_dd_obs(sat1, Some(100.0), Some(200.0))),
            (make_dd_obs(sat2, Some(300.0), None), make_dd_obs(sat2, Some(300.0), None)),
        ];

        update_last_observed(&mut state, &matched);

        assert_eq!(state.last_observed[&(sat1, 1)], 5);
        assert_eq!(state.last_observed[&(sat1, 2)], 5);
        assert_eq!(state.last_observed[&(sat2, 1)], 5);
        assert!(!state.last_observed.contains_key(&(sat2, 2)));
    }

    #[test]
    fn test_update_last_observed_skips_none() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 3;

        let sat1 = make_test_sat(1);
        let matched = vec![
            (make_dd_obs(sat1, None, None), make_dd_obs(sat1, None, None)),
        ];

        update_last_observed(&mut state, &matched);

        assert!(!state.last_observed.contains_key(&(sat1, 1)));
        assert!(!state.last_observed.contains_key(&(sat1, 2)));
    }

    #[test]
    fn test_evaluate_gnn_no_model() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);
        let gnn_variances = ProcessingEngine::evaluate_gnn(
            &None,
            &[],
            &make_test_epoch(time),
            &[],
            &state,
        );
        assert!(gnn_variances.is_empty());
    }

    #[test]
    fn test_evaluate_gnss_only_coasting_rtk_resets_when_young() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::Rtk;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 0;

        let spp_pos = Some(Coordinate::new(
            Vector3::new(10.0, 10.0, 10.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ));

        ProcessingEngine::evaluate_gnss_only_coasting(&config, &mut state, spp_pos, None, false);

        // Should have been reset to SPP because epoch_count < 3
        assert!((state.position.vector.x - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_evaluate_gnss_only_coasting_rtk_no_reset_when_mature() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::Rtk;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(5.0, 5.0, 5.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 5;

        let spp_pos = Some(Coordinate::new(
            Vector3::new(10.0, 10.0, 10.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ));

        ProcessingEngine::evaluate_gnss_only_coasting(&config, &mut state, spp_pos, None, false);

        // Should NOT reset because epoch_count >= 3
        assert!((state.position.vector.x - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_evaluate_gnss_only_coasting_ins_with_imu_data() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::RtkIns;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(5.0, 5.0, 5.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 0;

        let spp_pos = Some(Coordinate::new(
            Vector3::new(10.0, 10.0, 10.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ));

        // With IMU data, INS modes should NOT trigger GNSS-only coasting
        ProcessingEngine::evaluate_gnss_only_coasting(&config, &mut state, spp_pos, None, true);
        assert!((state.position.vector.x - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_evaluate_gnss_only_coasting_ins_without_imu_data() {
        let mut config = EngineConfig::default();
        config.mode = EngineMode::RtkIns;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(5.0, 5.0, 5.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 0;

        let spp_pos = Some(Coordinate::new(
            Vector3::new(10.0, 10.0, 10.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ));

        // Without IMU data, INS modes should trigger GNSS-only coasting
        ProcessingEngine::evaluate_gnss_only_coasting(&config, &mut state, spp_pos, None, false);
        assert!((state.position.vector.x - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_update_state_time_sets_time_and_clears_reset() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.is_reset = true;

        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.current_state = Some(state);

        let new_time = GpsTime::new(0, 10.0);
        engine.update_state_time(new_time).unwrap();

        let s = engine.current_state.as_ref().unwrap();
        assert!(!s.is_reset);
        assert_eq!(s.time.tow, 10.0);
    }

    #[test]
    fn test_update_state_time_no_state_errors() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let err = engine.update_state_time(GpsTime::new(0, 10.0)).unwrap_err();
        assert!(matches!(err, EngineError::StateDisappeared));
    }

    #[test]
    fn test_handle_ekf_rejection_increments() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let config = EngineConfig::default();

        handle_ekf_rejection(&mut state, &config, None, None, "test");
        assert_eq!(state.consecutive_rejections, 1);

        handle_ekf_rejection(&mut state, &config, None, None, "test");
        assert_eq!(state.consecutive_rejections, 2);
    }

    #[test]
    fn test_handle_ekf_rejection_ins_aligned_inflates_covariance() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;

        let cov_0_0_before = state.covariance[(0, 0)];
        let cov_3_3_before = state.covariance[(3, 3)];

        handle_ekf_rejection(&mut state, &EngineConfig::default(), None, None, "test");

        // Position covariances should be inflated
        assert!(state.covariance[(0, 0)] > cov_0_0_before);
        assert!(state.covariance[(3, 3)] > cov_3_3_before);
    }

    #[test]
    fn test_handle_ekf_rejection_not_aligned_no_reset_below_3() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(5.0, 5.0, 5.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = false;
        state.consecutive_rejections = 2;

        handle_ekf_rejection(&mut state, &EngineConfig::default(), None, None, "test");
        assert_eq!(state.consecutive_rejections, 3);
        // Position should NOT be reset because we didn't provide spp_pos
        assert!((state.position.vector.x - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_handle_ekf_rejection_not_aligned_resets_at_3_with_spp() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(5.0, 5.0, 5.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = false;
        state.consecutive_rejections = 2;

        let spp_pos = Coordinate::new(
            Vector3::new(50.0, 50.0, 50.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        handle_ekf_rejection(&mut state, &EngineConfig::default(), Some(spp_pos), None, "test");
        // Should have been reset
        assert!((state.position.vector.x - 50.0).abs() < 1e-6);
        assert_eq!(state.consecutive_rejections, 0);
    }

    #[test]
    fn test_handle_ekf_acceptance_clears_rejections() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.consecutive_rejections = 5;

        handle_ekf_acceptance(&mut state, &EngineConfig::default(), &[], None, None);
        assert_eq!(state.consecutive_rejections, 0);
    }

    #[test]
    fn test_handle_ekf_acceptance_resets_on_extreme_variance() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(5.0, 5.0, 5.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.covariance[(0, 0)] = 20000.0;
        state.covariance[(1, 1)] = 0.0;
        state.covariance[(2, 2)] = 0.0;
        state.consecutive_rejections = 5;

        let spp_pos = Coordinate::new(
            Vector3::new(50.0, 50.0, 50.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        handle_ekf_acceptance(&mut state, &EngineConfig::default(), &[], Some(spp_pos), None);
        // Should be reset
        assert!((state.position.vector.x - 50.0).abs() < 1e-6);
        assert_eq!(state.consecutive_rejections, 0);
    }

    #[test]
    fn test_build_measurement_environment_constructs() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);
        let base_coord = Coordinate::new(
            Vector3::new(110.0, 210.0, 310.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        let env = build_measurement_environment(
            &config,
            &[],
            &state,
            &[],
            &base_coord,
            time,
        );

        assert_eq!(env.base_coord.vector.x, 110.0);
        assert!(env.gnn_variances.is_empty());
    }

    #[test]
    fn test_perform_spp_fallback_update_does_not_panic() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::new(100.0, 200.0, 300.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.covariance = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);

        let new_pos = Coordinate::new(
            Vector3::new(101.0, 201.0, 301.0),
            Datum::WGS84,
            Frame::ECEF,
            GpsTime::new(0, 1.0),
        );

        // Should not panic
        ProcessingEngine::perform_spp_fallback_update(&config, &mut state, new_pos);
    }

    #[test]
    fn test_apply_adaptive_r_scaling_basic() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);

        let mut r = DMatrix::identity(2, 2);
        let z = DVector::from_vec(vec![1.0; 2]);
        let mut h = DMatrix::zeros(2, CORE_STATE_SIZE);
        h[(0, 0)] = 1.0;
        h[(1, 0)] = 1.0;

        let sat = make_test_sat(1);
        let meas_types = vec![(sat, 0, 0.0), (sat, 0, 0.0)];
        let matched_obs = vec![
            (make_dd_obs(sat, Some(100.0), None), make_dd_obs(sat, Some(100.0), None)),
        ];

        let mut tracker = crate::engine::adaptive::InnovationTracker::new();

        // Just ensure it doesn't panic and modifies r
        apply_adaptive_r_scaling(&mut tracker, &state, &z, &h, &mut r, &meas_types, &matched_obs);
        // R should be >= 1.0 after scaling (since min scale is 1.0)
        assert!(r[(0, 0)] >= 1.0);
        // Since innovation is 1.0 and predicted_var = h*P*h' + r = 1.0*1.0*1.0 + 1.0 = 2.0
        // nis = 1.0/2.0 = 0.5, so scale should be close to 1.0
        assert!((r[(0, 0)] - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_evaluate_gnss_only_coasting_no_spp_pos_does_nothing() {
        // When spp_pos is None, the function should not reset regardless of epoch_count.
        let mut config = EngineConfig::default();
        config.mode = EngineMode::Rtk;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(5.0, 5.0, 5.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.epoch_count = 0; // Would trigger reset IF spp_pos were Some

        ProcessingEngine::evaluate_gnss_only_coasting(&config, &mut state, None, None, false);

        // Position should remain unchanged
        assert!((state.position.vector.x - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_handle_ekf_rejection_ins_aligned_extreme_divergence_resets() {
        // When ins_aligned and position variance exceeds 10000,
        // handle_ekf_rejection should reset to SPP.
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(5.0, 5.0, 5.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;
        // Make position variance extreme
        state.covariance[(0, 0)] = 5000.0;
        state.covariance[(1, 1)] = 5000.0;
        state.covariance[(2, 2)] = 5000.0;

        let spp_pos = Coordinate::new(
            Vector3::new(50.0, 50.0, 50.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        handle_ekf_rejection(&mut state, &EngineConfig::default(), Some(spp_pos), None, "extreme");
        // Position should have been reset to SPP
        assert!((state.position.vector.x - 50.0).abs() < 1e-6);
        assert_eq!(state.consecutive_rejections, 0);
    }

    #[test]
    fn test_handle_ekf_rejection_ins_aligned_extreme_no_spp_does_not_panic() {
        // When ins_aligned and variance is extreme but no spp_pos provided,
        // the reset branch is skipped; consecutive_rejections still increments.
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(5.0, 5.0, 5.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;
        state.covariance[(0, 0)] = 20000.0;

        handle_ekf_rejection(&mut state, &EngineConfig::default(), None, None, "no_rescue");
        // Without spp_pos, the reset is skipped but consecutive_rejections is incremented
        assert_eq!(state.consecutive_rejections, 1);
        // Original position should be preserved
        assert!((state.position.vector.x - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_init_spp_state_when_state_already_exists() {
        // When current_state already exists, init_spp_state should not fail
        // even if SPP compute fails (no ephemerides).
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));

        let rover = EpochObs { time, satellites: vec![] };
        let result = engine.init_spp_state(&rover);
        // Should not error because state already exists
        assert!(result.is_ok());
        // SPP result should be None (no ephemerides to compute)
        assert!(result.unwrap().is_none());
    }

    #[test]
    fn test_apply_observations_with_fewer_than_5_matched_obs() {
        // When fewer than 5 observations match, consecutive_rejections should increase.
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));
        engine.config.base_position = Some([100.0, 200.0, 300.0]);

        let state = engine.current_state.as_mut().unwrap();
        state.consecutive_rejections = 0;

        let rover = EpochObs { time, satellites: vec![] };
        let base = EpochObs { time, satellites: vec![] };

        // This will call apply_observations with empty observations, resulting in
        // matched_obs being empty (< 5), which increments rejections.
        let result = engine.apply_observations(&rover, Some(&base), None, None);
        assert!(result.is_ok());
        assert_eq!(
            engine.current_state.as_ref().unwrap().consecutive_rejections,
            1
        );
    }

    #[test]
    fn test_apply_observations_with_valid_base_and_spp_fallback() {
        // When no valid base (age too large) but spp_pos is available,
        // apply_observations should call perform_spp_fallback_update (does not panic).
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        let time = GpsTime::new(0, 100.0);
        let pos = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));
        let state = engine.current_state.as_mut().unwrap();
        state.covariance = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);

        // Base obs with large time difference (age > max_base_age_s = 5.0)
        let old_time = GpsTime::new(0, 1.0);
        let rover = EpochObs { time, satellites: vec![] };
        let base = EpochObs { time: old_time, satellites: vec![] };

        let spp_pos = Some(Coordinate::new(
            Vector3::new(101.0, 201.0, 301.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        ));

        // Should not panic — SPP fallback update is applied
        let result = engine.apply_observations(&rover, Some(&base), spp_pos, None);
        assert!(result.is_ok());
    }

    #[test]
    fn test_process_rtk_basic_flow() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.base_position = Some([100.0, 200.0, 300.0]);

        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));

        let rover = EpochObs { time, satellites: vec![] };
        let base = EpochObs { time, satellites: vec![] };

        let result = engine.process_rtk(&rover, Some(&base));
        assert!(result.is_ok(), "process_rtk should succeed: {:?}", result.err());
        assert_eq!(engine.state_history.len(), 1);
        assert_eq!(engine.obs_history.len(), 1);
    }

    #[test]
    fn test_process_rtk_update_rejection() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let base_coord = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        let ctx = RtkUpdateContext {
            config: &config,
            ephemerides: &[],
            imu_history: &[],
            rover_obs: &EpochObs { time, satellites: vec![] },
            base_obs: &EpochObs { time, satellites: vec![] },
            matched_obs: &[],
            base_coord: &base_coord,
            spp_pos: None,
            spp_state_ref: None,
            gnn_variances: std::collections::HashMap::new(),
        };

        let mut tracker = crate::engine::adaptive::InnovationTracker::new();
        let rejections_before = state.consecutive_rejections;
        process_rtk_update::<TightCoupling>(&mut state, &mut tracker, &ctx);

        assert!(
            state.consecutive_rejections > rejections_before,
            "Rejections should increase after empty update"
        );
    }

    #[test]
    fn test_execute_ekf_update_rejection_nan() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        let base_coord = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        let ctx = RtkUpdateContext {
            config: &config,
            ephemerides: &[],
            imu_history: &[],
            rover_obs: &EpochObs { time, satellites: vec![] },
            base_obs: &EpochObs { time, satellites: vec![] },
            matched_obs: &[],
            base_coord: &base_coord,
            spp_pos: None,
            spp_state_ref: None,
            gnn_variances: std::collections::HashMap::new(),
        };

        let core_size = CORE_STATE_SIZE;
        let m = crate::engine::measurement::EkfMeasurementMatrices {
            z: DVector::from_vec(vec![f64::NAN]),
            h: DMatrix::zeros(1, core_size),
            r: DMatrix::identity(1, 1),
            mt: vec![(make_test_sat(1), 0, 1575.42e6)],
        };

        let rejections_before = state.consecutive_rejections;
        execute_ekf_update::<TightCoupling>(&mut state, &ctx, &m);

        assert!(
            state.consecutive_rejections > rejections_before,
            "NaN innovation should trigger rejection"
        );
    }

    #[test]
    fn test_execute_ekf_update_acceptance() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, pos, 1.0);
        state.covariance = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        state.consecutive_rejections = 3;

        let base_coord = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        let ctx = RtkUpdateContext {
            config: &config,
            ephemerides: &[],
            imu_history: &[],
            rover_obs: &EpochObs { time, satellites: vec![] },
            base_obs: &EpochObs { time, satellites: vec![] },
            matched_obs: &[],
            base_coord: &base_coord,
            spp_pos: None,
            spp_state_ref: None,
            gnn_variances: std::collections::HashMap::new(),
        };

        let core_size = CORE_STATE_SIZE;
        let mut h = DMatrix::zeros(1, core_size);
        h[(0, 0)] = 1.0;
        let m = crate::engine::measurement::EkfMeasurementMatrices {
            z: DVector::from_vec(vec![0.0]),
            h,
            r: DMatrix::identity(1, 1),
            mt: vec![(make_test_sat(1), 0, 1575.42e6)],
        };

        execute_ekf_update::<TightCoupling>(&mut state, &ctx, &m);

        assert_eq!(state.consecutive_rejections, 0, "Acceptance should clear rejections");
    }

    #[test]
    fn test_build_measurement_environment_with_imu_and_lever_arm() {
        let mut config = EngineConfig::default();
        config.imu_to_antenna_lever_arm = [0.5, 0.0, 1.0];

        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.ins_aligned = true;
        state.gyro_bias = Vector3::new(0.01, 0.02, 0.03);

        let base_coord = Coordinate::new(
            Vector3::new(110.0, 210.0, 310.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );

        let imu_history = vec![vec![gneiss_core::imu::ImuMeasurement {
            accel: Vector3::new(0.0, 0.0, 9.8),
            gyro: Vector3::new(0.1, 0.0, 0.0),
            time_tag: 0,
            temperature: None,
        }]];

        let env = build_measurement_environment(
            &config,
            &imu_history,
            &state,
            &[],
            &base_coord,
            time,
        );

        assert!(env.lever_arm.norm() > 0.0, "Lever arm should be non-zero when ins_aligned");
        assert!(env.omega_b.norm() > 0.0, "Omega should be non-zero with gyro data and bias");
    }

    #[test]
    fn test_export_gnn_dataset_no_path_does_nothing() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);
        let config = EngineConfig::default();
        let base_coord = Coordinate::new(
            Vector3::new(110.0, 210.0, 310.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let ctx = RtkUpdateContext {
            config: &config,
            ephemerides: &[],
            imu_history: &[],
            rover_obs: &EpochObs { time, satellites: vec![] },
            base_obs: &EpochObs { time, satellites: vec![] },
            matched_obs: &[],
            base_coord: &base_coord,
            spp_pos: None,
            spp_state_ref: None,
            gnn_variances: std::collections::HashMap::new(),
        };
        export_gnn_dataset(
            &state,
            &ctx,
            &crate::engine::measurement::EkfMeasurementMatrices {
                z: DVector::from_vec(vec![0.0]),
                h: DMatrix::zeros(CORE_STATE_SIZE, CORE_STATE_SIZE),
                r: DMatrix::identity(1, 1),
                mt: vec![(make_test_sat(1), 0, 1575.42e6)],
            },
            &[],
            &DVector::from_vec(vec![0.0]),
        );
    }

    #[test]
    fn test_init_spp_state_fails_when_no_state_and_spp_fails() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.current_state = None;
        let time = GpsTime::new(0, 0.0);
        let rover = EpochObs { time, satellites: vec![] };
        let result = engine.init_spp_state(&rover);
        assert!(matches!(result, Err(EngineError::InitialSppFailed)));
    }

    #[test]
    fn test_get_base_coord_with_helmert_transform() {
        let mut config = EngineConfig::default();
        config.base_position = Some([100.0, 200.0, 300.0]);
        config.base_datum_transform = Some(HelmertParams {
            tx: 1.0, ty: 2.0, tz: 3.0,
            rx: 0.0, ry: 0.0, rz: 0.0,
            s: 0.0,
            dtx: 0.0, dty: 0.0, dtz: 0.0,
            drx: 0.0, dry: 0.0, drz: 0.0,
            ds: 0.0,
            ref_epoch: 2000.0,
        });
        let rover_obs = EpochObs {
            time: GpsTime::new(0, 0.0),
            satellites: vec![],
        };
        let coord = ProcessingEngine::get_base_coord(&config, &rover_obs).unwrap();
        // Helmert applies [1,2,3]m translation with zero rates
        assert!((coord.vector.x - 101.0).abs() < 1.0);
        assert!((coord.vector.y - 202.0).abs() < 1.0);
        assert!((coord.vector.z - 303.0).abs() < 1.0);
    }

    #[test]
    fn test_process_rtk_fails_without_base_position() {
        let mut engine = ProcessingEngine::new(EngineConfig::default());
        engine.config.base_position = None;
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(1.0, 2.0, 3.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        engine.current_state = Some(RtkState::new(time, pos, 1.0));

        let rover = EpochObs { time, satellites: vec![] };
        let base = EpochObs { time, satellites: vec![] };
        let err = engine.process_rtk(&rover, Some(&base)).unwrap_err();
        assert!(matches!(err, EngineError::MissingBasePosition));
    }

    #[test]
    fn test_perform_spp_fallback_update_large_diff_does_not_panic() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(100.0, 200.0, 300.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.covariance = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);

        let new_pos = Coordinate::new(
            Vector3::new(1000000.0, 2000000.0, 3000000.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        ProcessingEngine::perform_spp_fallback_update(&config, &mut state, new_pos);
    }

    #[test]
    fn test_handle_ekf_acceptance_extreme_variance_no_spp_does_not_reset() {
        let config = EngineConfig::default();
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(
            Vector3::new(5.0, 5.0, 5.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.covariance[(0, 0)] = 20000.0;
        state.consecutive_rejections = 5;

        handle_ekf_acceptance(&mut state, &config, &[], None, None);
        assert_eq!(state.consecutive_rejections, 0);
        assert!((state.position.vector.x - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_apply_adaptive_r_scaling_with_freq2() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);

        let mut r = DMatrix::identity(1, 1);
        let z = DVector::from_vec(vec![1.0]);
        let mut h = DMatrix::zeros(1, CORE_STATE_SIZE);
        h[(0, 0)] = 1.0;

        let sat = make_test_sat(1);
        let meas_types = vec![(sat, 2, 1227.6e6)];
        let matched_obs = vec![(
            make_dd_obs(sat, Some(100.0), None),
            make_dd_obs(sat, Some(100.0), None),
        )];

        let mut tracker = crate::engine::adaptive::InnovationTracker::new();
        apply_adaptive_r_scaling(
            &mut tracker, &state, &z, &h, &mut r, &meas_types, &matched_obs,
        );
        assert!(r[(0, 0)] >= 1.0);
    }

    #[test]
    fn test_apply_adaptive_r_scaling_no_matching_snr() {
        let time = GpsTime::new(0, 0.0);
        let pos = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, time);
        let state = RtkState::new(time, pos, 1.0);

        let mut r = DMatrix::identity(1, 1);
        let z = DVector::from_vec(vec![1.0]);
        let mut h = DMatrix::zeros(1, CORE_STATE_SIZE);
        h[(0, 0)] = 1.0;

        let sat1 = make_test_sat(1);
        let sat2 = make_test_sat(2);
        let meas_types = vec![(sat2, 0, 1575.42e6)];
        let matched_obs = vec![(
            make_dd_obs(sat1, Some(100.0), None),
            make_dd_obs(sat1, Some(100.0), None),
        )];

        let mut tracker = crate::engine::adaptive::InnovationTracker::new();
        apply_adaptive_r_scaling(
            &mut tracker, &state, &z, &h, &mut r, &meas_types, &matched_obs,
        );
        assert!(r[(0, 0)] >= 1.0);
    }
}
