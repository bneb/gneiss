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

    #[test]
    fn every_error_variant_renders_a_distinct_diagnostic() {
        // These strings are the operator-facing diagnostics, and the payloads
        // must actually reach the message: swapping two Display arms, or
        // dropping the payload, fails here.
        let cases = [
            (EngineError::InversionError, "Matrix inversion failed"),
            (EngineError::SingularState, "Singular covariance or state"),
            (
                EngineError::InvalidMeasurement("3x5 vs 5x3".into()),
                "Invalid measurement: 3x5 vs 5x3",
            ),
            (EngineError::EmptyHistory, "Empty history for smoother"),
            (
                EngineError::Internal("phi not square".into()),
                "Internal engine error: phi not square",
            ),
        ];
        for (err, expect) in cases {
            assert_eq!(format!("{err}"), expect);
            // Every variant is a std::error::Error, so ?-propagation composes.
            let dyn_err: &dyn std::error::Error = &err;
            assert_eq!(dyn_err.to_string(), expect);
        }
        // Same variant, different payload -> different message.
        assert_ne!(
            EngineError::InvalidMeasurement("a".into()),
            EngineError::InvalidMeasurement("b".into())
        );
        assert_eq!(
            EngineError::InvalidMeasurement("a".into()),
            EngineError::InvalidMeasurement("a".into())
        );
    }

    #[test]
    fn default_covariance_is_the_documented_diagonal_in_state_order() {
        // 1 m position, 0.1 m^2/s^2 velocity (sigma = 0.3162 m/s),
        // 0.01 rad^2 attitude (sigma = 0.1 rad = 5.730 deg),
        // 0.04 (m/s^2)^2 accel bias (sigma = 0.2), 1e-4 (rad/s)^2 gyro bias
        // (sigma = 0.01). All cross terms are exactly zero.
        let p = default_covariance();
        let want = [1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01, 0.04, 0.04, 0.04, 1e-4, 1e-4, 1e-4];
        for i in 0..15 {
            assert!((p[(i, i)] - want[i]).abs() < 1e-18, "P[{i},{i}] = {}", p[(i, i)]);
            for j in 0..15 {
                if i != j {
                    assert_eq!(p[(i, j)], 0.0, "P[{i},{j}] must be zero");
                }
            }
        }
        assert!((p[(6, 6)].sqrt() - 5.7295779513_f64.to_radians()).abs() < 1e-12);
        assert!((p[(3, 3)].sqrt() - 0.1_f64.sqrt()).abs() < 1e-15);
        assert!((p[(9, 9)].sqrt() - 0.2).abs() < 1e-15);
        assert!((p[(12, 12)].sqrt() - 0.01).abs() < 1e-15);
    }

    #[test]
    fn with_cov_preserves_every_field_including_biases() {
        let cov = default_covariance() * 7.0;
        let s = EskfState::with_cov(
            Vector3::new(1.0, 2.0, 3.0),
            Vector3::new(4.0, 5.0, 6.0),
            UnitQuaternion::from_axis_angle(&Vector3::z_axis(), 0.3),
            Vector3::new(0.01, 0.02, 0.03),
            Vector3::new(0.001, 0.002, 0.003),
            cov,
        );
        assert_eq!(s.pos_ecef, Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(s.vel_ecef, Vector3::new(4.0, 5.0, 6.0));
        assert_eq!(s.accel_bias, Vector3::new(0.01, 0.02, 0.03));
        assert_eq!(s.gyro_bias, Vector3::new(0.001, 0.002, 0.003));
        assert_eq!(s.cov, cov);
        assert_eq!(s.attitude.angle(), 0.3);
        // new() must be the zero-bias, default-covariance special case.
        let n = EskfState::new(Vector3::new(1.0, 2.0, 3.0), Vector3::new(4.0, 5.0, 6.0), s.attitude);
        assert_eq!(n.accel_bias, Vector3::zeros());
        assert_eq!(n.gyro_bias, Vector3::zeros());
        assert_eq!(n.cov, default_covariance());
    }

    #[test]
    fn skew_symmetric_reproduces_the_cross_product() {
        // skew(v) * w must equal v x w exactly for every v, w:
        //   (1,2,3) x (4,5,6) = (2*6-3*5, 3*4-1*6, 1*5-2*4) = (-3, 6, -3)
        let v = Vector3::new(1.0, 2.0, 3.0);
        let w = Vector3::new(4.0, 5.0, 6.0);
        let got = skew_symmetric(&v) * w;
        assert!((got - Vector3::new(-3.0, 6.0, -3.0)).norm() < 1e-15, "got {got:?}");
        assert!((skew_symmetric(&v) * w - v.cross(&w)).norm() < 1e-15);
        // Antisymmetry and the standard basis images: skew(e_z) * e_x = e_y.
        let ez = skew_symmetric(&Vector3::z());
        assert!((ez * Vector3::x() - Vector3::y()).norm() < 1e-15);
        assert!((ez * Vector3::y() + Vector3::x()).norm() < 1e-15);
    }

    #[test]
    fn normal_gravity_matches_somigliana_and_points_to_the_centre() {
        // Somigliana (closed form, no iteration needed):
        //   g0(lat) = gE * (1 + k sin^2 lat) / sqrt(1 - e^2 sin^2 lat)
        // with gE = 9.7803253359, k = 0.00193185265241, e^2 = 0.00669437999014.
        const G_E: f64 = 9.780_325_335_9;
        const K: f64 = 0.001_931_852_652_41;
        const E2: f64 = 0.006_694_379_990_14;

        // Equator (sin^2 = 0): g0 = gE exactly, directed along -ECEF x here.
        let eq = normal_gravity_ecef(&Vector3::new(6_378_137.0, 0.0, 0.0));
        assert!((eq - Vector3::new(-G_E, 0.0, 0.0)).norm() < 1e-12, "got {eq:?}");

        // Pole (sin^2 = 1): g0 = gE (1 + k) / sqrt(1 - e^2) = 9.83218... m/s^2.
        let b = 6_378_137.0 * (1.0 - 1.0 / 298.257223563);
        let gp = normal_gravity_ecef(&Vector3::new(0.0, 0.0, b));
        let expect_pole = G_E * (1.0 + K) / (1.0 - E2).sqrt();
        assert!((gp.norm() - expect_pole).abs() < 1e-9, "polar g = {}", gp.norm());
        assert!((gp.norm() - 9.832_185).abs() < 1e-6, "polar g = {}", gp.norm());
        assert!((gp - Vector3::new(0.0, 0.0, -gp.norm())).norm() < 1e-9, "must point at the centre");

        // 45 deg: sin_lat = 1/sqrt(2) is built in exactly by this position
        // (x^2 + y^2 + z^2 = 0.25 + 0.25 + 0.5 = 1), giving
        //   g0 = gE (1 + k/2) / sqrt(1 - e^2/2) = 9.806197... m/s^2.
        let r = 6_378_137.0;
        let p45 = Vector3::new(0.5, 0.5, 0.5_f64.sqrt()) * r;
        let g45 = normal_gravity_ecef(&p45);
        let expect_45 = G_E * (1.0 + K / 2.0) / (1.0 - E2 / 2.0).sqrt();
        assert!((g45.norm() - expect_45).abs() < 1e-9, "45 deg g = {}", g45.norm());
        assert!((g45.norm() - 9.806_198).abs() < 1e-6, "45 deg g = {}", g45.norm());
        // Gravity strictly increases from equator to pole.
        assert!(g45.norm() > eq.norm() && gp.norm() > g45.norm());

        // Gravity must always be antiparallel to the radius: a non-radial
        // component would be a spurious horizontal pull.
        let p = Vector3::new(4_027_893.0, 307_041.0, 4_919_475.0);
        let g = normal_gravity_ecef(&p);
        let sin_angle = g.cross(&p).norm() / (g.norm() * p.norm());
        assert!(sin_angle < 1e-15, "gravity is not radial at {p:?}: sin = {sin_angle}");
        assert!(g.dot(&p) < 0.0, "gravity must point inward");
        // (1,1,2) has geocentric latitude asin(2/sqrt(6)) = 54.7 deg > 45 deg,
        // so its gravity must exceed the 45 deg value. (Using (1,1,1) here would
        // be wrong: that is only 35.3 deg, below 45 deg.)
        assert!(normal_gravity_ecef(&Vector3::new(1.0, 1.0, 2.0)).norm() > g45.norm());
    }

    #[test]
    fn earth_rotation_rate_is_the_iers_nominal_value_along_ecef_z() {
        // IERS nominal Earth rotation rate, 7.2921151467e-5 rad/s, directed
        // along the ECEF z axis (the spin axis).
        let w = earth_rotation_rate_ecef();
        assert!((w.x).abs() < 1e-18 && (w.y).abs() < 1e-18);
        assert!((w.z - 7.292_115_146_7e-5).abs() < 1e-19, "w_z = {}", w.z);
        assert!((w.z - WGS84_EARTH_ROTATION_RATE).abs() < 1e-19);
        // In one hour the Earth turns 3600 * 7.2921151467e-5 = 0.2625161453 rad
        // = 15.0410670 deg (the familiar "15 deg per hour" to within 0.05 %).
        let rad_per_hour = 3600.0 * w.z;
        assert!((rad_per_hour - 0.262_516_145_281_2).abs() < 1e-13, "rad/hr = {rad_per_hour}");
        assert!((rad_per_hour.to_degrees() - 15.041_067_03).abs() < 1e-6,
            "deg/hr = {}", rad_per_hour.to_degrees());
    }

    #[test]
    fn clamp_vector_bounds_each_axis_independently() {
        let mut v = Vector3::new(0.05, -2.0, -0.05);
        clamp_vector(&mut v, 0.1);
        assert_eq!(v, Vector3::new(0.05, -0.1, -0.05));
        // Exactly on the limit is unchanged.
        let mut v = Vector3::new(0.1, -0.1, 0.1);
        clamp_vector(&mut v, 0.1);
        assert_eq!(v, Vector3::new(0.1, -0.1, 0.1));
        // Non-finite input is left to the comparison semantics, not silently
        // turned into zero.
        let mut v = Vector3::new(f64::NAN, 5.0, -5.0);
        clamp_vector(&mut v, 1.0);
        assert!(v.x.is_nan());
        assert_eq!(v.y, 1.0);
        assert_eq!(v.z, -1.0);
    }
}
