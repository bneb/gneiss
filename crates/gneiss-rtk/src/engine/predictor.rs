use crate::engine::fgo::factors::skew_symmetric;
use crate::filter::RtkState;
use nalgebra::{DMatrix, DVector, Matrix3, UnitQuaternion, Vector3};

use crate::engine::{DynamicsModel, EngineConfig};

pub fn integrate_imu_mechanization(
    state: &mut RtkState,
    dt: f64,
    imu_buffer: &[gneiss_core::imu::ImuMeasurement],
) {
    let imu_dt = dt / (imu_buffer.len() as f64);
    let omega_ie = Vector3::new(0.0, 0.0, gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S);

    for meas in imu_buffer {
        let f_b = meas.accel - state.accel_bias;
        let omega_b = meas.gyro - state.gyro_bias;

        let zeta = (omega_b - state.attitude.inverse() * omega_ie) * imu_dt;
        let angle = zeta.norm();
        let dq = if angle > 1e-12 {
            UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_unchecked(zeta / angle), angle)
        } else {
            UnitQuaternion::identity()
        };

        let dq_mid = if angle > 1e-12 {
            UnitQuaternion::from_axis_angle(
                &nalgebra::Unit::new_unchecked(zeta / angle),
                angle * 0.5,
            )
        } else {
            UnitQuaternion::identity()
        };
        let r_mid = state.attitude * dq_mid;

        state.attitude *= dq;
        state.attitude.renormalize();

        let f_e = r_mid * f_b;
        let gravity = gravity_wgs84(state.position.vector);
        let coriolis = 2.0 * omega_ie.cross(&state.velocity);
        let centrifugal = omega_ie.cross(&(omega_ie.cross(&state.position.vector)));

        let v_dot = f_e + gravity - coriolis - centrifugal;

        if state.epoch_count == 60 {
            tracing::info!(
                "IMU MECHANIZATION: f_b={:.2?} f_e={:.2?} gravity={:.2?} v_dot={:.2?} att={:.3?}",
                f_b.as_slice(),
                f_e.as_slice(),
                gravity.as_slice(),
                v_dot.as_slice(),
                state.attitude.coords.as_slice()
            );
        }

        let v_mid = state.velocity + v_dot * (imu_dt * 0.5);
        state.velocity += v_dot * imu_dt;
        state.position.vector += v_mid * imu_dt;
    }
}

pub fn compute_transition_matrix(
    state: &RtkState,
    dt: f64,
    imu_buffer: &[gneiss_core::imu::ImuMeasurement],
) -> DMatrix<f64> {
    let n = state.covariance.nrows();
    let mut phi = DMatrix::<f64>::identity(n, n);
    let omega_ie = Vector3::new(0.0, 0.0, gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S);

    if imu_buffer.is_empty() {
        phi[(0, 3)] = dt;
        phi[(1, 4)] = dt;
        phi[(2, 5)] = dt;
    } else {
        let r_b_e = state.attitude.to_rotation_matrix();
        let f_e = state.attitude * (imu_buffer.last().expect("imu_buffer is non-empty").accel - state.accel_bias);
        let f_e_skew = skew_symmetric(&f_e);
        let omega_ie_skew = skew_symmetric(&omega_ie);

        for i in 0..3 {
            phi[(i, 3 + i)] = dt;
        }

        // Sign convention: the error state stored at indices 6-8 is d_theta
        // (the correction axis-angle) such that R_true = (I + [d_theta×]) R_est.
        // This d_theta = -ψ where ψ is the conventional attitude error vector.
        // Under this convention the velocity-attitude coupling IS
        // -[f_e×]*dt = -f_e_skew * dt (since δv = -[f_e×] dt d_theta derives
        // from δf_e = -[d_theta×] f_e = +[f_e×] d_theta integrated over dt).
        let vel_att = -f_e_skew * dt;
        for r in 0..3 {
            for c in 0..3 {
                phi[(3 + r, 6 + c)] = vel_att[(r, c)];
            }
        }

        let vel_vel = Matrix3::identity() - 2.0 * omega_ie_skew * dt;
        for r in 0..3 {
            for c in 0..3 {
                phi[(3 + r, 3 + c)] = vel_vel[(r, c)];
            }
        }

        let vel_abias = -r_b_e.matrix() * dt;
        for r in 0..3 {
            for c in 0..3 {
                phi[(3 + r, 9 + c)] = vel_abias[(r, c)];
            }
        }

        let att_att = Matrix3::identity() - omega_ie_skew * dt;
        for r in 0..3 {
            for c in 0..3 {
                phi[(6 + r, 6 + c)] = att_att[(r, c)];
            }
        }

        let att_gbias = -r_b_e.matrix() * dt;
        for r in 0..3 {
            for c in 0..3 {
                phi[(6 + r, 12 + c)] = att_gbias[(r, c)];
            }
        }
    }

    if crate::filter::CORE_STATE_SIZE > 15 {
        // Random-walk clock model (φ[15,15]=1.0) preserves temporal
        // correlation of clock bias across epochs.  With process_noise_cb
        // tuned to the receiver oscillator (~1 m²/s for TCXO), the clock
        // covariance stays small after convergence, anchoring absolute
        // position through clock-position measurement coupling.
        // RTKLIB uses white-noise clock PER EPOCH, but that works because
        // their SPP re-anchors position each epoch — Gneiss carries
        // position forward, so the clock must carry forward too.
        phi[(15, 19)] = dt;
    }

    phi
}

pub fn compute_process_noise(
    dt: f64,
    config: &EngineConfig,
    is_imu_active: bool,
    is_fixed: bool,
    ambiguity_keys: &[(gneiss_core::sat::SatelliteId, u8)],
) -> DMatrix<f64> {
    let n = crate::filter::CORE_STATE_SIZE + ambiguity_keys.len();
    let mut q = DMatrix::<f64>::zeros(n, n);
    let dt_abs = dt.abs();

    if !is_imu_active {
        if config.dynamics_model == DynamicsModel::Static {
            // Static: position process noise is a small constant (σ≈1mm/√s),
            // matching the RTKLIB port. The velocity-integration model
            // (q_pos ∝ dt³) produces 9 m²/epoch at 30s, destroying any
            // tight position prior from RINEX header.
            for i in 0..3 {
                q[(i, i)] = 1e-6 * dt_abs;
            }
            for i in 3..6 {
                q[(i, i)] = 1e-6 * dt_abs; // velocity: small noise (static receiver)
            }
        } else {
            let q_acc = match config.dynamics_model {
                DynamicsModel::Static => 0.001, // unreachable
                DynamicsModel::Pedestrian => 1.0,
                DynamicsModel::Marine => 2.0,
                DynamicsModel::Automotive => 10.0,
                DynamicsModel::Airborne => 50.0,
            };
            let q_pos = q_acc * dt_abs.powi(3) / 3.0;
            let q_vel = q_acc * dt_abs;
            let q_pos_vel = q_acc * dt_abs.powi(2) / 2.0;
            for i in 0..3 {
                q[(i, i)] = q_pos;
                q[(i + 3, i + 3)] = q_vel;
                q[(i, i + 3)] = q_pos_vel;
                q[(i + 3, i)] = q_pos_vel;
            }
        }
        for i in 6..9 {
            q[(i, i)] = 1e-7 * dt_abs;
        }
    } else {
        let q_vel = config.tuning.sigma_v * config.tuning.sigma_v * dt_abs;
        let q_att = config.tuning.sigma_phi * config.tuning.sigma_phi * dt_abs;
        let q_ab = config.tuning.sigma_ab * config.tuning.sigma_ab * dt_abs;
        let q_gb = config.tuning.sigma_gb * config.tuning.sigma_gb * dt_abs;
        let q_pos = q_vel * dt_abs * dt_abs / 3.0; // position uncertainty from velocity noise integration
        for i in 0..3 {
            q[(i, i)] = q_pos;
            q[(3 + i, 3 + i)] = q_vel;
            q[(6 + i, 6 + i)] = q_att;
            q[(9 + i, 9 + i)] = q_ab;
            q[(12 + i, 12 + i)] = q_gb;
        }
    }

    if crate::filter::CORE_STATE_SIZE > 15 {
        q[(15, 15)] = config.process_noise_cb * dt_abs;
        q[(16, 16)] = config.process_noise_isb * dt_abs;
        q[(17, 17)] = config.process_noise_isb * dt_abs;
        q[(18, 18)] = config.process_noise_isb * dt_abs;
        q[(19, 19)] = config.process_noise_cd * dt_abs;
        q[(20, 20)] = config.process_noise_zwd * dt_abs;
    }

    for (i, key) in ambiguity_keys.iter().enumerate() {
        let idx = crate::filter::CORE_STATE_SIZE + i;
        if key.1 == 3 {
            q[(idx, idx)] = config.process_noise_iono * dt_abs;
        } else {
            q[(idx, idx)] = if is_fixed {
                config.process_noise_amb_fixed * dt_abs
            } else {
                config.process_noise_amb_float * dt_abs
            };
        }
    }

    q
}

pub fn predict(
    state: &mut RtkState,
    dt: f64,
    config: &EngineConfig,
    imu_buffer: &[gneiss_core::imu::ImuMeasurement],
) {
    if imu_buffer.is_empty() {
        state.position.vector += state.velocity * dt;
    } else {
        integrate_imu_mechanization(state, dt, imu_buffer);
    }

    if crate::filter::CORE_STATE_SIZE > 15 {
        state.rcv_clk_bias += state.rcv_clk_drift * dt;
    }

    let phi = compute_transition_matrix(state, dt, imu_buffer);
    let q = compute_process_noise(
        dt,
        config,
        !imu_buffer.is_empty(),
        state.is_fixed,
        &state.ambiguity_keys,
    );

    state.core_phi = Some(
        phi.view(
            (0, 0),
            (
                crate::filter::CORE_STATE_SIZE,
                crate::filter::CORE_STATE_SIZE,
            ),
        )
        .into_owned(),
    );

    let mut phi_full = DMatrix::identity(state.covariance.nrows(), state.covariance.ncols());
    phi_full
        .view_mut(
            (0, 0),
            (
                crate::filter::CORE_STATE_SIZE,
                crate::filter::CORE_STATE_SIZE,
            ),
        )
        .copy_from(state.core_phi.as_ref().expect("core_phi should be initialized before use"));

    // Clamp extreme covariance values to prevent overflow in Phi*P*Phi^T.
    // Cycle slip inflation (×4 per slip) can push position variance toward
    // f64 overflow.  Cap at 1e10 m² (σ=100 km) — any real filter would
    // have been reset long before reaching this.
    for i in 0..state.covariance.nrows() {
        for j in 0..state.covariance.ncols() {
            let v = &mut state.covariance[(i, j)];
            if v.abs() > 1e10 { *v = v.signum() * 1e10; }
        }
    }

    state.covariance = &phi_full * &state.covariance * phi_full.transpose() + q;

    // NaN can enter through numerical overflow in Phi*P*Phi^T.  Once
    // present, it cascades into StateDisappeared errors that cause 30m
    // position jumps.  Replace the entire matrix with a conservative
    // diagonal (fast Cholesky path, no SVD hang).  This discards
    // off-diagonal information but preserves filter continuity — far
    // better than a full state reset.
    if state.covariance.iter().any(|x| x.is_nan()) {
        let n = state.covariance.nrows();
        let core = crate::filter::CORE_STATE_SIZE.min(n);
        let nan_count = state.covariance.iter().filter(|x| x.is_nan()).count();
        tracing::warn!(
            "NaN in predicted covariance ({}×{}, {} entries). Reinit diagonal.",
            n, n, nan_count
        );
        state.covariance = DMatrix::zeros(n, n);
        // Position: 25 m² (σ=5m), Velocity: 1 m²/s², Attitude: 0.01 rad²
        for i in 0..3  { state.covariance[(i, i)] = 25.0; }
        for i in 3..6  { state.covariance[(i, i)] = 1.0; }
        for i in 6..9  { state.covariance[(i, i)] = 0.01; }
        // Accel bias, gyro bias
        for i in 9..15 { state.covariance[(i, i)] = 0.01; }
        // Clock bias: 10000 m², ISBs: 100, clock drift: 100, ZWD: 0.01
        if n > 15 { state.covariance[(15, 15)] = 10000.0; }
        if n > 16 { state.covariance[(16, 16)] = 100.0; }
        if n > 17 { state.covariance[(17, 17)] = 100.0; }
        if n > 18 { state.covariance[(18, 18)] = 100.0; }
        if n > 19 { state.covariance[(19, 19)] = 100.0; }
        if n > 20 { state.covariance[(20, 20)] = 0.01; }
        // Ambiguities: 10000 m²
        for i in core..n { state.covariance[(i, i)] = 10000.0; }
    }
    state.full_p_predict = Some(state.covariance.clone());

    let mut x_pred = DVector::zeros(state.covariance.nrows());
    x_pred.rows_mut(0, 3).copy_from(&state.position.vector);
    x_pred.rows_mut(3, 3).copy_from(&state.velocity);
    if state.covariance.nrows() > 6 {
        // Attitude indices 6-8 are intentionally omitted (left as zero).
        // After injection the error-state attitude is expected to be zero,
        // and the smoother handles attitude separately via the quaternion
        // predicted_attitude field rather than through x_pred.
        x_pred.rows_mut(9, 3).copy_from(&state.accel_bias);
        x_pred.rows_mut(12, 3).copy_from(&state.gyro_bias);
    }
    if crate::filter::CORE_STATE_SIZE > 15 {
        x_pred[15] = state.rcv_clk_bias;
        x_pred[16] = state.isb_glo;
        x_pred[17] = state.isb_gal;
        x_pred[18] = state.isb_bds;
        x_pred[19] = state.rcv_clk_drift;
        x_pred[20] = state.zwd;
    }
    for i in 0..state.ambiguities.len() {
        x_pred[crate::filter::CORE_STATE_SIZE + i] = state.ambiguities[i];
    }
    state.full_x_predict = Some(x_pred);

    state.predicted_position = Some(state.position);
    state.predicted_velocity = Some(state.velocity);
    if state.covariance.nrows() > 6 {
        state.predicted_attitude = Some(state.attitude);
        state.predicted_accel_bias = Some(state.accel_bias);
        state.predicted_gyro_bias = Some(state.gyro_bias);
    }
}

pub fn gravity_wgs84(pos_ecef: Vector3<f64>) -> Vector3<f64> {
    let x = pos_ecef.x;
    let y = pos_ecef.y;
    let z = pos_ecef.z;
    let r = pos_ecef.norm();
    if r < 1.0 {
        return Vector3::zeros();
    }

    let r2 = r * r;
    let r3 = r2 * r;
    let a = gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M;
    let mu = 3.986005e14;
    let j2 = 1.082627e-3;

    let a_r_2 = (a / r) * (a / r);
    let z_r_2 = (z / r) * (z / r);

    let g_base = -mu / r3;
    let g_j2_common = 1.5 * j2 * a_r_2;

    let gx = g_base * x * (1.0 - g_j2_common * (5.0 * z_r_2 - 1.0));
    let gy = g_base * y * (1.0 - g_j2_common * (5.0 * z_r_2 - 1.0));
    let gz = g_base * z * (1.0 - g_j2_common * (5.0 * z_r_2 - 3.0));

    Vector3::new(gx, gy, gz)
}
#[cfg(test)]
mod tests {
    use crate::engine::predictor::{
        compute_process_noise, compute_transition_matrix, gravity_wgs84,
        integrate_imu_mechanization, predict,
    };
    use crate::engine::{DynamicsModel, EngineConfig};
    use crate::filter::RtkState;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::imu::ImuMeasurement;
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;
    use nalgebra::{DMatrix, UnitQuaternion, Vector3};

    #[test]
    fn test_predictor_indices() {
        // Test that process noise and state transition matrix use the correct indices for ISBs and clock drift
        let mut state = RtkState::new(
            gneiss_core::time::GpsTime::new(2000, 0.0),
            gneiss_core::coords::Coordinate::new(
                nalgebra::Vector3::zeros(),
                gneiss_core::coords::Datum::WGS84,
                gneiss_core::coords::Frame::ECEF,
                gneiss_core::time::GpsTime::new(2000, 0.0),
            ),
            10.0,
        );
        state.covariance = DMatrix::zeros(21, 21);
        state.isb_glo = 10.0;
        state.isb_gal = 20.0;
        state.isb_bds = 30.0;
        state.rcv_clk_bias = 100.0;
        state.rcv_clk_drift = 2.0;
        state.zwd = 0.5;

        let config = EngineConfig {
            mode: crate::engine::EngineMode::Ppp,
            initial_position: None,
            base_position: None,
            base_datum_transform: None,
            receiver_antenna_type: None,
            imu_to_antenna_lever_arm: [0.0, 0.0, 0.0],
            imu_mounting_angles: None,
            imu_to_nhc_lever_arm: [0.0, 0.0, 0.0],
            enable_nhc: false,
            enable_backward_smoothing: false,
            lambda_min_ratio: 3.0,
            lambda_min_subset: 4,
            enabled_constellations: None,
            raim_pseudorange_outlier_m: 10.0,
            chi_square_pr_threshold: 15.0,
            chi_square_cp_threshold: 15.0,
            phase_windup_enabled: true,
            min_snr_dbhz: 0.0,
            dynamics_model: DynamicsModel::Static,
            doppler_slip_threshold_cycles: 5.0,
            max_reject_count: 3,
            max_base_age_s: 30.0,
            spp_consistency_threshold_m: 10.0,
            initial_ambiguity_variance: 100.0,
            ar_min_epoch_count: 10,
            ar_min_lock: 3,
            ar_ffrt_prob: 0.001,
            iono_model: crate::engine::types::IonosphereModel::Klobuchar,
            enable_tropo_gradients: false,
            enable_ar: false,
            elevation_mask_deg: 5.0,
            auto_detect_dynamics: false,
            process_noise_cb: 100.0,
            process_noise_cd: 10.0,
            process_noise_isb: 0.1,
            process_noise_zwd: 1e-8,
            process_noise_iono: 1e-6,
            enable_multi_base_rtk: false,
            enable_gnn_raim: false,
            export_gnn_dataset_path: None,
            process_noise_amb_float: 1e-4,
            process_noise_amb_fixed: 1e-7,
            uduc_ar: false,
            tropo_mapping: gneiss_core::atmosphere::TropoMapping::default(),
            tuning: crate::engine::config::EkfTuningConfig::default(),
        };

        crate::engine::predictor::predict(&mut state, 1.0, &config, &[]);

        // Check that clock bias is updated by drift (nominal state update unchanged)
        assert_eq!(state.rcv_clk_bias, 102.0); // 100.0 + 2.0 * 1.0

        let x_pred = state.full_x_predict.as_ref().unwrap();

        // Ensure ISBs and zwd are unchanged by predict
        assert_eq!(x_pred[16], 10.0);
        assert_eq!(x_pred[17], 20.0);
        assert_eq!(x_pred[18], 30.0);
        assert_eq!(x_pred[19], 2.0);
        assert_eq!(x_pred[20], 0.5);

        // Ensure covariance is updated correctly
        let p_pred = state.full_p_predict.as_ref().unwrap();
        // Clock drift noise goes to 19
        assert_eq!(p_pred[(19, 19)], 10.0);
        // ZWD noise goes to 20
        assert_eq!(p_pred[(20, 20)], 1e-8);
        // Clock bias: white-noise prediction (φ=0 so P_pred ≈ q + dt²·P_drift)
        assert!(p_pred[(15, 15)] >= 100.0);
        // ISB noises: process_noise_isb * dt = 0.1
        assert_eq!(p_pred[(16, 16)], 0.1);
        assert_eq!(p_pred[(17, 17)], 0.1);
        assert_eq!(p_pred[(18, 18)], 0.1);
    }

    fn default_config() -> EngineConfig {
        EngineConfig {
            mode: crate::engine::EngineMode::Ppp,
            initial_position: None,
            base_position: None,
            base_datum_transform: None,
            receiver_antenna_type: None,
            imu_to_antenna_lever_arm: [0.0, 0.0, 0.0],
            imu_mounting_angles: None,
            imu_to_nhc_lever_arm: [0.0, 0.0, 0.0],
            enable_nhc: false,
            enable_backward_smoothing: false,
            lambda_min_ratio: 3.0,
            lambda_min_subset: 4,
            enabled_constellations: None,
            raim_pseudorange_outlier_m: 10.0,
            chi_square_pr_threshold: 15.0,
            chi_square_cp_threshold: 15.0,
            phase_windup_enabled: true,
            min_snr_dbhz: 0.0,
            dynamics_model: DynamicsModel::Static,
            doppler_slip_threshold_cycles: 5.0,
            max_reject_count: 3,
            max_base_age_s: 30.0,
            spp_consistency_threshold_m: 10.0,
            initial_ambiguity_variance: 100.0,
            ar_min_epoch_count: 10,
            ar_min_lock: 3,
            ar_ffrt_prob: 0.001,
            iono_model: crate::engine::types::IonosphereModel::Klobuchar,
            enable_tropo_gradients: false,
            enable_ar: false,
            elevation_mask_deg: 5.0,
            auto_detect_dynamics: false,
            process_noise_cb: 100.0,
            process_noise_cd: 10.0,
            process_noise_isb: 0.1,
            process_noise_zwd: 1e-8,
            process_noise_iono: 1e-6,
            enable_multi_base_rtk: false,
            enable_gnn_raim: false,
            export_gnn_dataset_path: None,
            process_noise_amb_float: 1e-4,
            process_noise_amb_fixed: 1e-7,
            uduc_ar: false,
            tropo_mapping: gneiss_core::atmosphere::TropoMapping::default(),
            tuning: crate::engine::config::EkfTuningConfig::default(),
        }
    }

    #[test]
    fn test_integrate_imu_mechanization_nonzero_angle() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(
            Vector3::new(6378137.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 10.0);
        state.attitude = UnitQuaternion::identity();
        state.velocity = Vector3::zeros();
        state.gyro_bias = Vector3::new(0.05, 0.0, 0.0);

        let imu_buffer = [ImuMeasurement {
            time_tag: 0,
            gyro: Vector3::new(0.1, 0.0, 0.0),
            accel: Vector3::new(0.0, 0.0, 9.8),
            temperature: None,
        }];

        integrate_imu_mechanization(&mut state, 1.0, &imu_buffer);

        // Attitude changed from identity (rotation due to non-zero omega_b)
        assert!(
            (state.attitude.to_rotation_matrix().matrix() - nalgebra::Matrix3::identity())
                .norm()
                > 1e-6
        );

        // Velocity changed due to gravity + coriolis + centrifugal
        assert!(state.velocity.norm() > 1e-6);

        // Position changed due to velocity integration
        assert!(
            (state.position.vector - Vector3::new(6378137.0, 0.0, 0.0)).norm()
                > 1e-6
        );
    }

    #[test]
    fn test_integrate_imu_mechanization_zero_angle() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(
            Vector3::new(6378137.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 10.0);
        state.attitude = UnitQuaternion::identity();
        state.velocity = Vector3::zeros();
        // Set gyro_bias so that omega_b == omega_ie, making zeta = [0,0,0]
        let omega_ie_z: f64 = 7.2921151467e-5;
        state.gyro_bias = Vector3::new(0.1, 0.0, -omega_ie_z);

        let imu_buffer = [ImuMeasurement {
            time_tag: 0,
            gyro: Vector3::new(0.1, 0.0, 0.0),
            accel: Vector3::new(0.0, 0.0, 9.8),
            temperature: None,
        }];

        integrate_imu_mechanization(&mut state, 1.0, &imu_buffer);

        // Attitude remains identity (angle = 0, dq stays identity)
        assert!(
            (state.attitude.to_rotation_matrix().matrix() - nalgebra::Matrix3::identity())
                .norm()
                < 1e-10
        );

        // Velocity changed due to gravity
        assert!(state.velocity.norm() > 1e-6);
    }

    #[test]
    fn test_compute_transition_matrix_with_imu_buffer() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(
            Vector3::new(6378137.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 10.0);
        state.velocity = Vector3::new(10.0, 5.0, 1.0);
        state.accel_bias = Vector3::new(0.1, 0.2, 0.3);
        state.gyro_bias = Vector3::new(0.01, 0.02, 0.03);
        state.attitude = UnitQuaternion::from_axis_angle(&Vector3::z_axis(), 0.3);

        let imu_buffer = [ImuMeasurement {
            time_tag: 0,
            gyro: Vector3::new(0.1, 0.0, 0.0),
            accel: Vector3::new(0.0, 0.0, 9.8),
            temperature: None,
        }];

        let phi = compute_transition_matrix(&state, 1.0, &imu_buffer);

        assert_eq!(phi.nrows(), 21);
        assert_eq!(phi.ncols(), 21);

        // Position-velocity coupling (phi[i, 3+i] = dt)
        assert!((phi[(0, 3)] - 1.0).abs() < 1e-12);
        assert!((phi[(1, 4)] - 1.0).abs() < 1e-12);
        assert!((phi[(2, 5)] - 1.0).abs() < 1e-12);

        // Velocity-attitude coupling is populated (f_e is non-zero)
        // Diagonal of skew-symmetric is zero, so check off-diagonal
        assert!(phi[(3, 7)].abs() > 1e-10);
        assert!(phi[(4, 6)].abs() > 1e-10);

        // Velocity-accel_bias coupling is populated
        assert!(phi[(3, 9)].abs() > 1e-10);

        // Clock bias-drift coupling (CORE_STATE_SIZE=21 > 15)
        assert!((phi[(15, 19)] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_compute_transition_matrix_empty_imu() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(
            Vector3::new(6378137.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 10.0);
        state.velocity = Vector3::new(10.0, 0.0, 0.0);

        let phi = compute_transition_matrix(&state, 1.0, &[]);

        assert_eq!(phi.nrows(), 21);
        assert_eq!(phi.ncols(), 21);

        // Position-velocity coupling
        assert!((phi[(0, 3)] - 1.0).abs() < 1e-12);
        assert!((phi[(1, 4)] - 1.0).abs() < 1e-12);
        assert!((phi[(2, 5)] - 1.0).abs() < 1e-12);

        // Velocity-attitude coupling is zero (IMU false path)
        assert_eq!(phi[(3, 6)], 0.0);
        assert_eq!(phi[(4, 7)], 0.0);
        assert_eq!(phi[(5, 8)], 0.0);

        // Clock bias-drift coupling (CORE_STATE_SIZE=21 > 15)
        assert!((phi[(15, 19)] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_compute_process_noise_various_dynamics_models() {
        let models = [
            (DynamicsModel::Static, 0.001),
            (DynamicsModel::Pedestrian, 1.0),
            (DynamicsModel::Marine, 2.0),
            (DynamicsModel::Automotive, 10.0),
            (DynamicsModel::Airborne, 50.0),
        ];

        for (model, q_acc) in models {
            let config = EngineConfig {
                dynamics_model: model,
                ..default_config()
            };
            let q = compute_process_noise(1.0, &config, false, false, &[]);

            // Static: constant 1e-6×dt, not velocity-integration model
            let q_pos_expected = if model == DynamicsModel::Static {
                1e-6
            } else {
                q_acc / 3.0
            };
            assert!(
                (q[(0, 0)] - q_pos_expected).abs() < 1e-12,
                "Mismatch for {:?}: expected q[0,0] = {}, got {}",
                model,
                q_pos_expected,
                q[(0, 0)]
            );

            if model == DynamicsModel::Static {
                // Static: no velocity/acceleration coupling, large velocity noise
                assert!((q[(0, 3)] - 0.0).abs() < 1e-12);
                assert!((q[(3, 3)] - 100.0).abs() < 1e-12);
                assert!((q[(6, 6)] - 1e-7).abs() < 1e-15);
            } else {
                // Position-velocity cross term: q_pos_vel = q_acc * dt^2 / 2
                let q_pos_vel_expected = q_acc * 0.5;
                assert!((q[(0, 3)] - q_pos_vel_expected).abs() < 1e-12);
                // Velocity variance: q_vel = q_acc * dt
                assert!((q[(3, 3)] - q_acc).abs() < 1e-12);
                // Attitude variance set to 1e-7 * dt
                assert!((q[(6, 6)] - 1e-7).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn test_compute_process_noise_imu_active() {
        let config = default_config();
        let q = compute_process_noise(1.0, &config, true, false, &[]);

        let sigma_v = config.tuning.sigma_v;
        let sigma_phi = config.tuning.sigma_phi;
        let sigma_ab = config.tuning.sigma_ab;
        let sigma_gb = config.tuning.sigma_gb;

        // q_pos = sigma_v^2 * dt^3 / 3
        let q_pos_expected = sigma_v * sigma_v / 3.0;
        assert!((q[(0, 0)] - q_pos_expected).abs() < 1e-16);

        // q_vel = sigma_v^2 * dt
        let q_vel_expected = sigma_v * sigma_v;
        assert!((q[(3, 3)] - q_vel_expected).abs() < 1e-16);

        // q_att = sigma_phi^2 * dt
        let q_att_expected = sigma_phi * sigma_phi;
        assert!((q[(6, 6)] - q_att_expected).abs() < 1e-16);

        // q_ab = sigma_ab^2 * dt
        let q_ab_expected = sigma_ab * sigma_ab;
        assert!((q[(9, 9)] - q_ab_expected).abs() < 1e-20);

        // q_gb = sigma_gb^2 * dt
        let q_gb_expected = sigma_gb * sigma_gb;
        assert!((q[(12, 12)] - q_gb_expected).abs() < 1e-20);

        // Clock and ISB entries (CORE_STATE_SIZE=21 > 15)
        assert!((q[(15, 15)] - 100.0).abs() < 1e-10);
        assert!((q[(16, 16)] - 0.1).abs() < 1e-15);
        assert!((q[(17, 17)] - 0.1).abs() < 1e-15);
        assert!((q[(18, 18)] - 0.1).abs() < 1e-15);
        assert!((q[(19, 19)] - 10.0).abs() < 1e-10);
        assert!((q[(20, 20)] - 1e-8).abs() < 1e-15);
    }

    #[test]
    fn test_compute_process_noise_with_ambiguity_keys() {
        let sat1 = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let keys = vec![(sat1, 0), (sat1, 3), (sat1, 1)];
        let config = default_config();
        let cs = crate::filter::CORE_STATE_SIZE;

        // Test with is_fixed: false
        let q_float = compute_process_noise(1.0, &config, false, false, &keys);
        // key with .1 != 3 uses process_noise_amb_float
        assert!((q_float[(cs, cs)] - config.process_noise_amb_float).abs() < 1e-16);
        // key with .1 == 3 uses process_noise_iono
        assert!((q_float[(cs + 1, cs + 1)] - config.process_noise_iono).abs() < 1e-16);
        // key with .1 != 3 uses process_noise_amb_float
        assert!((q_float[(cs + 2, cs + 2)] - config.process_noise_amb_float).abs() < 1e-16);

        // Test with is_fixed: true
        let q_fixed = compute_process_noise(1.0, &config, false, true, &keys);
        // keys with .1 != 3 use process_noise_amb_fixed
        assert!((q_fixed[(cs, cs)] - config.process_noise_amb_fixed).abs() < 1e-16);
        // key with .1 == 3 still uses process_noise_iono
        assert!((q_fixed[(cs + 1, cs + 1)] - config.process_noise_iono).abs() < 1e-16);
        assert!((q_fixed[(cs + 2, cs + 2)] - config.process_noise_amb_fixed).abs() < 1e-16);
    }

    #[test]
    fn test_gravity_wgs84_at_origin() {
        // Vector with norm < 1.0 returns zeros
        let g0 = gravity_wgs84(Vector3::new(0.5, 0.5, 0.5));
        assert_eq!(g0, Vector3::zeros());

        // At Earth's surface on equator
        let g_eq = gravity_wgs84(Vector3::new(6378137.0, 0.0, 0.0));
        assert!(g_eq.x < 0.0); // gravity points toward center of Earth
        assert!(g_eq.y.abs() < 1e-10);
        assert!(g_eq.z.abs() < 1e-10);
        let mag = g_eq.norm();
        assert!((mag - 9.8).abs() < 0.1, "gravity magnitude {} not near 9.8", mag);
    }

    #[test]
    fn test_predict_with_imu_data() {
        let time = GpsTime::new(2000, 0.0);
        let pos = Coordinate::new(
            Vector3::new(6378137.0, 0.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            time,
        );
        let mut state = RtkState::new(time, pos, 10.0);
        state.velocity = Vector3::new(10.0, 0.0, 0.0);
        state.attitude = UnitQuaternion::identity();

        let imu_buffer = [
            ImuMeasurement {
                time_tag: 0,
                gyro: Vector3::new(0.1, 0.0, 0.0),
                accel: Vector3::new(0.0, 0.0, 9.8),
                temperature: None,
            },
            ImuMeasurement {
                time_tag: 1,
                gyro: Vector3::new(0.1, 0.0, 0.0),
                accel: Vector3::new(0.0, 0.0, 9.8),
                temperature: None,
            },
        ];

        let config = default_config();
        let cov_before = state.covariance.clone();

        predict(&mut state, 1.0, &config, &imu_buffer);

        // Position changed (via integrate_imu_mechanization)
        assert!(
            (state.position.vector - Vector3::new(6378137.0, 0.0, 0.0)).norm()
                > 1e-6
        );

        // Covariance updated (phi * P * phi^T + Q)
        let cov_diff = (&state.covariance - &cov_before).norm();
        assert!(cov_diff > 1e-10);

        // full_x_predict was set
        assert!(state.full_x_predict.is_some());

        // full_p_predict was set
        assert!(state.full_p_predict.is_some());

        // predicted_position was set
        assert!(state.predicted_position.is_some());

        // predicted_velocity was set
        assert!(state.predicted_velocity.is_some());

        // predicted_attitude was set
        assert!(state.predicted_attitude.is_some());
    }
}
