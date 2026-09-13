use std::fmt;
use nalgebra::{Matrix3, SMatrix, SVector, UnitQuaternion, Vector3};

pub const WGS84_EARTH_ROTATION_RATE: f64 = 7.292_115_146_7e-5;

pub type Vector15<T = f64> = SVector<T, 15>;
pub type Matrix15<T = f64> = SMatrix<T, 15, 15>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EngineError {
    InversionError,
    SingularState,
    InvalidMeasurement(String),
    EmptyHistory,
    Internal(String),
}

impl fmt::Display for EngineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InversionError => write!(f, "Matrix inversion failed"),
            Self::SingularState => write!(f, "Singular covariance or state"),
            Self::InvalidMeasurement(msg) => write!(f, "Invalid measurement: {msg}"),
            Self::EmptyHistory => write!(f, "Empty history for smoother"),
            Self::Internal(msg) => write!(f, "Internal engine error: {msg}"),
        }
    }
}

impl std::error::Error for EngineError {}

#[derive(Debug, Clone, PartialEq)]
pub struct EskfState {
    pub pos_ecef: Vector3<f64>,
    pub vel_ecef: Vector3<f64>,
    pub attitude: UnitQuaternion<f64>,
    pub accel_bias: Vector3<f64>,
    pub gyro_bias: Vector3<f64>,
    pub cov: Matrix15<f64>,
}

impl EskfState {
    pub fn new(
        pos_ecef: Vector3<f64>,
        vel_ecef: Vector3<f64>,
        attitude: UnitQuaternion<f64>,
    ) -> Self {
        Self::with_cov(
            pos_ecef,
            vel_ecef,
            attitude,
            Vector3::zeros(),
            Vector3::zeros(),
            default_covariance(),
        )
    }

    pub fn with_cov(
        pos_ecef: Vector3<f64>,
        vel_ecef: Vector3<f64>,
        attitude: UnitQuaternion<f64>,
        accel_bias: Vector3<f64>,
        gyro_bias: Vector3<f64>,
        cov: Matrix15<f64>,
    ) -> Self {
        Self {
            pos_ecef,
            vel_ecef,
            attitude,
            accel_bias,
            gyro_bias,
            cov,
        }
    }
}

pub fn default_covariance() -> Matrix15<f64> {
    let mut cov = Matrix15::identity();
    for i in 0..3 {
        cov[(i, i)] = 1.0; // 1 m pos
        cov[(i + 3, i + 3)] = 0.1; // 0.3 m/s vel
        cov[(i + 6, i + 6)] = 0.01; // ~5.7 deg attitude
        cov[(i + 9, i + 9)] = 0.04; // 0.2 m/s^2 accel bias
        cov[(i + 12, i + 12)] = 1e-4; // 0.01 rad/s gyro bias
    }
    cov
}

#[derive(Debug, Clone)]
pub struct EskfSnapshot {
    pub time: gneiss_core::time::GpsTime,
    pub state_pred: EskfState,
    pub state_post: EskfState,
    pub phi: Matrix15<f64>,
    pub is_gnss_available: bool,
}

pub fn skew_symmetric(v: &Vector3<f64>) -> Matrix3<f64> {
    Matrix3::new(
        0.0, -v.z, v.y,
        v.z, 0.0, -v.x,
        -v.y, v.x, 0.0,
    )
}

pub fn normal_gravity_ecef(pos: &Vector3<f64>) -> Vector3<f64> {
    let norm = pos.norm();
    if norm < 1e-3 {
        return Vector3::new(0.0, 0.0, -9.780_325);
    }
    let sin_lat = (pos.z / norm).clamp(-1.0, 1.0);
    let sin2_lat = sin_lat * sin_lat;
    let g0 = 9.780_325_335_9 * (1.0 + 0.001_931_852_652_41 * sin2_lat)
        / (1.0 - 0.006_694_379_990_14 * sin2_lat).sqrt();
    -g0 * (pos / norm)
}

pub fn earth_rotation_rate_ecef() -> Vector3<f64> {
    Vector3::new(0.0, 0.0, WGS84_EARTH_ROTATION_RATE)
}

pub fn clamp_vector(v: &mut Vector3<f64>, max_abs: f64) {
    v.x = v.x.clamp(-max_abs, max_abs);
    v.y = v.y.clamp(-max_abs, max_abs);
    v.z = v.z.clamp(-max_abs, max_abs);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_skew_symmetric_properties() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        let s = skew_symmetric(&v);
        assert!((s + s.transpose()).norm() < 1e-12);
        let prod = s * v;
        assert!(prod.norm() < 1e-12);
    }

    #[test]
    fn test_normal_gravity_ecef_direction() {
        let pos = Vector3::new(6378137.0, 0.0, 0.0);
        let g = normal_gravity_ecef(&pos);
        assert!((g.norm() - 9.780_325).abs() < 1e-4);
        assert!(g.x < 0.0);
    }

    #[test]
    fn test_normal_gravity_near_zero() {
        let g0 = normal_gravity_ecef(&Vector3::zeros());
        assert!((g0.z + 9.780_325).abs() < 1e-4);
    }

    #[test]
    fn test_clamp_vector() {
        let mut v = Vector3::new(5.0, -10.0, 0.5);
        clamp_vector(&mut v, 2.0);
        assert_eq!(v, Vector3::new(2.0, -2.0, 0.5));
    }

    #[test]
    fn test_eskf_state_initialization() {
        let p = Vector3::new(1.0, 2.0, 3.0);
        let v = Vector3::new(0.1, 0.2, 0.3);
        let q = UnitQuaternion::identity();
        let state = EskfState::new(p, v, q);
        assert_eq!(state.pos_ecef, p);
        assert_eq!(state.vel_ecef, v);
        assert_eq!(state.accel_bias, Vector3::zeros());
        assert_eq!(state.gyro_bias, Vector3::zeros());
        assert_eq!(state.cov[(0, 0)], 1.0);
    }

    #[test]
    fn test_engine_error_formatting() {
        let err = EngineError::InvalidMeasurement("bad dim".into());
        let s = format!("{err}");
        assert!(s.contains("bad dim"));
    }
}
