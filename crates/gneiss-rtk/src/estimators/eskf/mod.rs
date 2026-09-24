//! 15-State Error-State Kalman Filter (ESKF/MEKF) for GNSS/INS Integration.
//!
//! State formulation:
//! - delta_p^e (0..3): ECEF position error
//! - delta_v^e (3..6): ECEF velocity error
//! - delta_theta (6..9): attitude error (left-multiplied global frame)
//! - delta_b_a (9..12): body accelerometer bias
//! - delta_b_g (12..15): body gyroscope bias
//!
//! Features:
//! - Multiplicative error-quaternion reset: q <- delta_q * q
//! - Closed-loop online bias estimation driven by GNSS innovations
//! - Coupled Non-Holonomic Constraints (NHC) with attitude Jacobian
//! - Zero-Velocity Updates (ZUPT)
//! - Full 15-state backward Rauch-Tung-Striebel (RTS) smoother

pub mod alignment;
pub mod condition;
pub mod constraints;
pub mod dd_update;
pub mod predict;
pub mod smoother;
pub mod types;
pub mod update;

pub use alignment::{
    compute_gyro_bias, compute_initial_attitude, compute_leveling_angles, init_eskf_filter,
};
pub use condition::{apply_integer_conditioning, ConditionSummary};
pub use constraints::{update_body_velocity, update_nhc, update_zupt};
pub use dd_update::{
    compute_dd_jacobian_15, update_dd_scalar, DdMeasurementKind, DdSatGeometry,
    DdScalarUpdateResult, RowVector15,
};
pub use predict::{
    compute_process_noise, compute_transition_matrix, predict, predict_preintegrated,
    predict_with_phi, propagate_nominal_state,
};
pub use smoother::EskfSmoother;
pub use types::{
    clamp_vector, earth_rotation_rate_ecef, normal_gravity_ecef, skew_symmetric, EngineError,
    EskfSnapshot, EskfState, Matrix15, Vector15,
};
pub use update::{
    apply_error_injection, build_doppler_velocity_system, build_gnss_pos_system,
    joseph_form_update, update_doppler_velocity, update_gnss_pos_vel, update_gnss_position,
    Matrix3x15,
};

