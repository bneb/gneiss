//! IMU preintegration for the sliding-window factor graph.
//!
//! Implements the Forster et al. (2015) preintegration method:
//!   - Accumulate Δp, Δv, Δq between GNSS epochs from raw IMU samples
//!   - Maintain Jacobians of deltas w.r.t. IMU biases
//!   - During LM iterations, correct deltas via first-order Taylor expansion
//!     (no reintegration needed)
//!
//! The ImuPreintegrationFactor connects:
//!   Pose(tₖ₋₁), Velocity(tₖ₋₁), ImuBias
//!   →
//!   Pose(tₖ), Velocity(tₖ)

use nalgebra::{DMatrix, DVector, Matrix3, UnitQuaternion, Vector3};
use std::fmt;

use crate::swfg::factor::Factor;
use crate::swfg::variables::{VariableId, VariableValues};

pub mod smoother;
pub mod stationary;

/// Preintegrated IMU measurements between two GNSS epochs.
///
/// The preintegration is performed once when IMU data arrives (not during
/// optimization).  During LM iterations, the bias-dependent corrections
/// use the stored Jacobians via first-order Taylor expansion:
///
///   Δp_new ≈ Δp_old + J_{p,ba} · δba + J_{p,bg} · δbg
///   Δv_new ≈ Δv_old + J_{v,ba} · δba + J_{v,bg} · δbg
///   Δq_new ≈ Δq_old ⊗ exp(J_{q,bg} · δbg)
#[derive(Debug, Clone)]
pub struct ImuPreintegration {
    /// Position change (ECEF, meters).
    pub dp: Vector3<f64>,
    /// Velocity change (ECEF, m/s).
    pub dv: Vector3<f64>,
    /// Attitude change (body-frame rotation).
    pub dq: UnitQuaternion<f64>,
    /// Total integration time (seconds).
    pub dt: f64,

    // Jacobians of deltas w.r.t. IMU biases (evaluated at the linearization point).
    pub dp_dba: Matrix3<f64>,
    pub dp_dbg: Matrix3<f64>,
    pub dv_dba: Matrix3<f64>,
    pub dv_dbg: Matrix3<f64>,
    pub dq_dbg: Matrix3<f64>,

    /// Preintegration covariance (15×15: [dp, dv, dq_axis_angle]).
    pub covariance: DMatrix<f64>,
    /// Whether this epoch interval was identified as stationary (for ZUPT).
    pub is_stationary: bool,
}

impl ImuPreintegration {
    pub fn new() -> Self {
        Self {
            dp: Vector3::zeros(),
            dv: Vector3::zeros(),
            dq: UnitQuaternion::identity(),
            dt: 0.0,
            dp_dba: Matrix3::zeros(),
            dp_dbg: Matrix3::zeros(),
            dv_dba: Matrix3::zeros(),
            dv_dbg: Matrix3::zeros(),
            dq_dbg: Matrix3::zeros(),
            covariance: DMatrix::identity(15, 15) * 1e-4,
            is_stationary: false,
        }
    }

    /// Integrate raw IMU samples between two GNSS epochs.
    ///
    /// `ba_i` and `bg_i` are the current best estimates of accel and gyro
    /// biases (the linearization point).  The Jacobians are accumulated
    /// so that later bias corrections can be applied without reintegration.
    pub fn integrate(&mut self, imu_data: &[ImuSample], ba_i: &Vector3<f64>, bg_i: &Vector3<f64>) {
        if imu_data.len() < 2 {
            return;
        }
        self.is_stationary = stationary::detect_stationary(imu_data);
        let mut prev_time = imu_data[0].time_us;

        for m in imu_data.iter().skip(1) {
            let dt = time_diff_us(prev_time, m.time_us);
            prev_time = m.time_us;
            let dt = dt.clamp(1e-4, 0.1); // 0.1ms to 100ms
            self.dt += dt;

            // Corrected acceleration in world frame at current attitude
            let accel_corrected = m.accel - ba_i;
            let a_world = self.dq * accel_corrected;

            // Standard mid-point integration
            self.dp += self.dv * dt + 0.5 * a_world * dt * dt;
            self.dv += a_world * dt;

            // Attitude update: q ← q ⊗ exp((ω - bg) * dt)
            let gyro_corrected = m.gyro - bg_i;
            let dq_inc = UnitQuaternion::from_scaled_axis(gyro_corrected * dt);
            let r_new = (self.dq * dq_inc).to_rotation_matrix().into_inner();
            self.dq *= dq_inc;

            // Accumulate bias Jacobians (Forster eq. 36-40)
            self.dp_dba += self.dv_dba * dt - 0.5 * r_new * dt * dt;
            self.dp_dbg += self.dv_dbg * dt;
            self.dv_dba -= r_new * dt;
            // dq_dbg is not accumulated here — it's computed from the final dq
        }

        // dq_dbg: Jacobian of quaternion w.r.t. gyro bias.
        // Approximated as -I * dt_total (small-angle approximation).
        // For proper implementation, see Forster eq. 40.
        self.dq_dbg = -Matrix3::identity() * self.dt;
    }

    /// Apply bias corrections to the preintegrated deltas using first-order
    /// Taylor expansion.  Called during LM iterations when biases change.
    pub fn correct(&self, dba: &Vector3<f64>, dbg: &Vector3<f64>) -> CorrectedPreintegration {
        let dp_corr = self.dp + self.dp_dba * dba + self.dp_dbg * dbg;
        let dv_corr = self.dv + self.dv_dba * dba + self.dv_dbg * dbg;
        let dq_corr = self.dq * UnitQuaternion::from_scaled_axis(self.dq_dbg * dbg);
        CorrectedPreintegration {
            dp: dp_corr,
            dv: dv_corr,
            dq: dq_corr,
            dt: self.dt,
            covariance: self.covariance.clone(),
            is_stationary: self.is_stationary,
        }
    }
}

impl Default for ImuPreintegration {
    fn default() -> Self {
        Self::new()
    }
}

/// Bias-corrected preintegration result (computed once per LM iteration).
#[derive(Debug, Clone)]
pub struct CorrectedPreintegration {
    pub dp: Vector3<f64>,
    pub dv: Vector3<f64>,
    pub dq: UnitQuaternion<f64>,
    pub dt: f64,
    pub covariance: DMatrix<f64>,
    pub is_stationary: bool,
}

/// A single IMU sample at a known time.
#[derive(Debug, Clone, Copy)]
pub struct ImuSample {
    /// Accelerometer measurement (m/s², body frame).
    pub accel: Vector3<f64>,
    /// Gyroscope measurement (rad/s, body frame).
    pub gyro: Vector3<f64>,
    /// Time tag in microseconds.
    pub time_us: u64,
}

/// Compute time difference between two microsecond time tags, handling rollover or wrap.
fn time_diff_us(prev: u64, curr: u64) -> f64 {
    const WEEK_US: u64 = 604_800_000_000;
    if curr >= prev {
        (curr - prev) as f64 / 1_000_000.0
    } else if prev <= WEEK_US && curr < 30_000_000 {
        ((WEEK_US - prev) + curr) as f64 / 1_000_000.0
    } else {
        ((u64::MAX - prev) + curr + 1) as f64 / 1_000_000.0
    }
}

// ---- Factor ---------------------------------------------------------------

/// Factor connecting IMU state at epoch `i` to epoch `j`.
///
/// Residual (15-DOF): difference between the predicted motion (from
/// integrating the nominal dynamics) and the preintegrated IMU deltas
/// (corrected for current bias estimates).
///
/// Connected variables (in order):
///   [0] Pose(epoch_i)    — 6-DOF [x, y, z, qx, qy, qz]
///   [1] Velocity(epoch_i) — 3-DOF
///   [2] Pose(epoch_j)    — 6-DOF
///   [3] Velocity(epoch_j) — 3-DOF
///   [4] ImuBias           — 6-DOF [ba_x, ba_y, ba_z, bg_x, bg_y, bg_z]
pub struct ImuPreintegrationFactor {
    /// Preintegrated deltas at the linearization point.
    pub preint: ImuPreintegration,
    /// Gravity vector in ECEF (m/s²).
    pub gravity: Vector3<f64>,

    // Nominal (linearization-point) values for each variable.
    pub nominal_p_i: Vector3<f64>,
    pub nominal_v_i: Vector3<f64>,
    pub nominal_q_i: UnitQuaternion<f64>,
    pub nominal_p_j: Vector3<f64>,
    pub nominal_v_j: Vector3<f64>,
    pub nominal_q_j: UnitQuaternion<f64>,
    pub nominal_ba: Vector3<f64>,
    pub nominal_bg: Vector3<f64>,

    // Variable IDs.
    var_pose_i: VariableId,
    var_vel_i: VariableId,
    var_pose_j: VariableId,
    var_vel_j: VariableId,
    var_bias: VariableId,

    // Cached variable list for the Factor trait.
    variables: Vec<VariableId>,
}

impl ImuPreintegrationFactor {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        preint: ImuPreintegration,
        gravity: Vector3<f64>,
        nominal_p_i: Vector3<f64>,
        nominal_v_i: Vector3<f64>,
        nominal_q_i: UnitQuaternion<f64>,
        nominal_p_j: Vector3<f64>,
        nominal_v_j: Vector3<f64>,
        nominal_q_j: UnitQuaternion<f64>,
        nominal_ba: Vector3<f64>,
        nominal_bg: Vector3<f64>,
        var_pose_i: VariableId,
        var_vel_i: VariableId,
        var_pose_j: VariableId,
        var_vel_j: VariableId,
        var_bias: VariableId,
    ) -> Self {
        let variables = vec![
            var_pose_i, var_vel_i, var_pose_j, var_vel_j, var_bias,
        ];
        Self {
            preint,
            gravity,
            nominal_p_i,
            nominal_v_i,
            nominal_q_i,
            nominal_p_j,
            nominal_v_j,
            nominal_q_j,
            nominal_ba,
            nominal_bg,
            var_pose_i,
            var_vel_i,
            var_pose_j,
            var_vel_j,
            var_bias,
            variables,
        }
    }
}

impl fmt::Debug for ImuPreintegrationFactor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ImuPreintegrationFactor")
            .field("dt", &self.preint.dt)
            .field("dp", &self.preint.dp)
            .finish()
    }
}

impl Factor for ImuPreintegrationFactor {
    fn variables(&self) -> &[VariableId] {
        &self.variables
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        // Extract current estimates from the variable graph
        let pose_i = values.get(self.var_pose_i).expect("pose_i in graph");
        let vel_i = values.get(self.var_vel_i).expect("vel_i in graph");
        let pose_j = values.get(self.var_pose_j).expect("pose_j in graph");
        let vel_j = values.get(self.var_vel_j).expect("vel_j in graph");
        let bias = values.get(self.var_bias).expect("bias in graph");

        let p_i = self.nominal_p_i + Vector3::new(pose_i[0], pose_i[1], pose_i[2]);
        let q_i = self.nominal_q_i
            * UnitQuaternion::from_scaled_axis(Vector3::new(pose_i[3], pose_i[4], pose_i[5]));
        let v_i = self.nominal_v_i + Vector3::new(vel_i[0], vel_i[1], vel_i[2]);

        let p_j = self.nominal_p_j + Vector3::new(pose_j[0], pose_j[1], pose_j[2]);
        let q_j = self.nominal_q_j
            * UnitQuaternion::from_scaled_axis(Vector3::new(pose_j[3], pose_j[4], pose_j[5]));
        let v_j = self.nominal_v_j + Vector3::new(vel_j[0], vel_j[1], vel_j[2]);

        let ba = self.nominal_ba + Vector3::new(bias[0], bias[1], bias[2]);
        let bg = self.nominal_bg + Vector3::new(bias[3], bias[4], bias[5]);

        let dba = ba - self.nominal_ba;
        let dbg = bg - self.nominal_bg;

        // Correct preintegration for current bias estimates
        let corr = self.preint.correct(&dba, &dbg);
        let dt = corr.dt;

        // Predicted motion from nominal dynamics:
        //   p_j = p_i + v_i*dt + 0.5*g*dt^2 + dp_imu  →  dp_pred = p_j - p_i - v_i*dt - 0.5*g*dt^2
        //   v_j = v_i + g*dt + dv_imu                  →  dv_pred = v_j - v_i - g*dt
        let g_dt = self.gravity * dt;
        let dp_pred_world = p_j - p_i - v_i * dt - 0.5 * g_dt * dt;
        let dv_pred_world = v_j - v_i - g_dt;
        let dq_pred = q_i.inverse() * q_j;

        let dp_pred_body = q_i.inverse() * dp_pred_world;
        let dv_pred_body = q_i.inverse() * dv_pred_world;

        // Residual = predicted - measured (preintegrated)
        let mut r = DVector::zeros(15);
        r.rows_mut(0, 3).copy_from(&(dp_pred_body - corr.dp));
        r.rows_mut(3, 3).copy_from(&(dv_pred_body - corr.dv));
        let dq_err = corr.dq.inverse() * dq_pred;
        r.rows_mut(6, 3).copy_from(&dq_err.scaled_axis());
        // Bias random walk (small penalty to keep estimates stable)
        r.rows_mut(9, 3).copy_from(&dba);
        r.rows_mut(12, 3).copy_from(&dbg);

        r
    }

    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64> {
        let total_dim = values.total_dim();
        let mut j = DMatrix::zeros(15, total_dim);

        // Get start indices for each variable
        let (start_pi, _) = values.index_of(self.var_pose_i).expect("pose_i in graph");
        let (start_vi, _) = values.index_of(self.var_vel_i).expect("vel_i in graph");
        let (start_pj, _) = values.index_of(self.var_pose_j).expect("pose_j in graph");
        let (start_vj, _) = values.index_of(self.var_vel_j).expect("vel_j in graph");
        let (start_bias, _) = values.index_of(self.var_bias).expect("bias in graph");

        let pose_i = values.get(self.var_pose_i).expect("pose_i in graph");
        let vel_i = values.get(self.var_vel_i).expect("vel_i in graph");
        let pose_j = values.get(self.var_pose_j).expect("pose_j in graph");
        let vel_j = values.get(self.var_vel_j).expect("vel_j in graph");

        let p_i = self.nominal_p_i + Vector3::new(pose_i[0], pose_i[1], pose_i[2]);
        let q_i = self.nominal_q_i
            * UnitQuaternion::from_scaled_axis(Vector3::new(pose_i[3], pose_i[4], pose_i[5]));
        let v_i = self.nominal_v_i + Vector3::new(vel_i[0], vel_i[1], vel_i[2]);

        let p_j = self.nominal_p_j + Vector3::new(pose_j[0], pose_j[1], pose_j[2]);
        let v_j = self.nominal_v_j + Vector3::new(vel_j[0], vel_j[1], vel_j[2]);

        let dt = self.preint.dt;
        let g_dt = self.gravity * dt;
        let dp_pred_world = p_j - p_i - v_i * dt - 0.5 * g_dt * dt;
        let dv_pred_world = v_j - v_i - g_dt;
        let dp_pred_body = q_i.inverse() * dp_pred_world;
        let dv_pred_body = q_i.inverse() * dv_pred_world;

        let r_i_t = q_i.inverse().to_rotation_matrix().into_inner();
        let r_corr_dq_inv = self.preint.dq.inverse().to_rotation_matrix().into_inner();

        fn skew(v: &Vector3<f64>) -> Matrix3<f64> {
            Matrix3::new(
                0.0, -v.z, v.y,
                v.z, 0.0, -v.x,
                -v.y, v.x, 0.0,
            )
        }

        // ∂r_p / ∂p_i = -R_i^T
        j.view_mut((0, start_pi), (3, 3)).copy_from(&-r_i_t);
        // ∂r_p / ∂q_i = [dp_pred_body]_x
        j.view_mut((0, start_pi + 3), (3, 3)).copy_from(&skew(&dp_pred_body));
        // ∂r_p / ∂v_i = -R_i^T * dt
        j.view_mut((0, start_vi), (3, 3)).copy_from(&(-r_i_t * dt));
        // ∂r_p / ∂p_j = R_i^T
        j.view_mut((0, start_pj), (3, 3)).copy_from(&r_i_t);

        // ∂r_v / ∂q_i = [dv_pred_body]_x
        j.view_mut((3, start_pi + 3), (3, 3)).copy_from(&skew(&dv_pred_body));
        // ∂r_v / ∂v_i = -R_i^T
        j.view_mut((3, start_vi), (3, 3)).copy_from(&-r_i_t);
        // ∂r_v / ∂v_j = R_i^T
        j.view_mut((3, start_vj), (3, 3)).copy_from(&r_i_t);

        // ∂r_q / ∂q_i = -R(corr.dq^{-1})
        j.view_mut((6, start_pi + 3), (3, 3)).copy_from(&-r_corr_dq_inv);
        // ∂r_q / ∂q_j = I
        j.view_mut((6, start_pj + 3), (3, 3)).copy_from(&Matrix3::identity());

        // ∂r/∂bias: from preintegration Jacobians
        for r_idx in 0..3 {
            for c_idx in 0..3 {
                // Position residual derivatives
                j[(r_idx, start_bias + c_idx)] = -self.preint.dp_dba[(r_idx, c_idx)];
                j[(r_idx, start_bias + 3 + c_idx)] = -self.preint.dp_dbg[(r_idx, c_idx)];
                // Velocity residual derivatives
                j[(3 + r_idx, start_bias + c_idx)] = -self.preint.dv_dba[(r_idx, c_idx)];
                j[(3 + r_idx, start_bias + 3 + c_idx)] = -self.preint.dv_dbg[(r_idx, c_idx)];
                // Attitude residual derivatives w.r.t. gyro bias
                j[(6 + r_idx, start_bias + 3 + c_idx)] = -self.preint.dq_dbg[(r_idx, c_idx)];
            }
        }
        // Bias random walk penalty
        for k in 0..6 {
            j[(9 + k, start_bias + k)] = 1.0;
        }

        j
    }

    fn information(&self) -> DMatrix<f64> {
        // Preintegration covariance inverted, plus bias prior
        let mut info = DMatrix::zeros(15, 15);
        // Use the preintegration covariance for the first 9 DOF (dp, dv, dq)
        let cov_9 = self.preint.covariance.view((0, 0), (9, 9)).into_owned();
        if let Some(inv) = cov_9.try_inverse() {
            info.view_mut((0, 0), (9, 9)).copy_from(&inv);
        }
        // Bias prior: 1e-2 rad/s and 1e-2 m/s² per sqrt(s) — loosely informative
        for k in 0..6 {
            info[(9 + k, 9 + k)] = 100.0; // σ = 0.1
        }
        info
    }

    fn robust_threshold(&self) -> Option<f64> {
        Some(3.0) // Huber with k=3 for IMU outliers
    }
}

// ---- Tests ----------------------------------------------------------------

#[cfg(test)]
mod tests;
