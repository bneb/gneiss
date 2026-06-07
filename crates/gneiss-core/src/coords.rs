use nalgebra::Vector3;
use crate::time::GpsTime;

pub trait GeodeticModel {
    fn a(&self) -> f64;
    fn f(&self) -> f64;
    fn b(&self) -> f64 { self.a() * (1.0 - self.f()) }
    fn e2(&self) -> f64 { 
        let a = self.a();
        let b = self.b();
        (a * a - b * b) / (a * a)
    }
}

pub struct Wgs84Model;
impl GeodeticModel for Wgs84Model {
    fn a(&self) -> f64 { 6378137.0 }
    fn f(&self) -> f64 { 1.0 / 298.257223563 }
}

/// Represents a geodetic reference datum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Datum {
    WGS84,
    ITRF2014,
    ITRF2020,
    JGD2011,
    PZ90,
    GTRF,
    CGCS2000,
}

/// Represents the coordinate frame for the vector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Frame {
    /// Earth-Centered, Earth-Fixed
    ECEF,
    /// Local Tangent Plane: East, North, Up (requires origin ECEF)
    ENU,
    /// Geodetic Latitude (rad), Longitude (rad), Height (m)
    LLH,
    /// Vehicle Body Frame: Forward, Right, Down
    Body,
}

/// A type-safe Coordinate wrapper preventing datum and frame mismatches.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Coordinate {
    pub vector: Vector3<f64>,
    pub datum: Datum,
    pub frame: Frame,
    pub epoch: GpsTime,
}

impl Coordinate {
    pub fn new(vector: Vector3<f64>, datum: Datum, frame: Frame, epoch: GpsTime) -> Self {
        Self { vector, datum, frame, epoch }
    }

    /// Ensures another coordinate matches this one exactly before mathematical operations.
    pub fn ensure_aligned(&self, other: &Coordinate) -> Result<(), &'static str> {
        if self.datum != other.datum {
            return Err("Datum mismatch");
        }
        if self.frame != other.frame {
            return Err("Frame mismatch");
        }
        // Strict epoch alignment requires them to be within 1 microsecond.
        if (self.epoch - other.epoch).abs() > 1e-6 {
            return Err("Epoch mismatch");
        }
        Ok(())
    }
}

/// Converts ECEF (meters) to Geodetic LLH (Latitude, Longitude, Height) in radians and meters.
pub fn ecef_to_llh_with_model<M: GeodeticModel>(model: &M, ecef: Vector3<f64>) -> Vector3<f64> {
    let a = model.a();
    let b = model.b();
    let e2 = model.e2();
    let ep2 = (a * a - b * b) / (b * b);
    
    let p = libm::sqrt(ecef.x * ecef.x + ecef.y * ecef.y);
    
    if p < 1e-10 {
        let lat = if ecef.z > 0.0 { core::f64::consts::FRAC_PI_2 } else { -core::f64::consts::FRAC_PI_2 };
        let height = libm::fabs(ecef.z) - b;
        return Vector3::new(lat, 0.0, height);
    }

    let theta = libm::atan2(ecef.z * a, p * b);
    let sin_theta = libm::sin(theta);
    let cos_theta = libm::cos(theta);

    let mut lat = libm::atan2(
        ecef.z + ep2 * b * sin_theta * sin_theta * sin_theta,
        p - e2 * a * cos_theta * cos_theta * cos_theta,
    );

    let lon = libm::atan2(ecef.y, ecef.x);

    // Iterative refinement for high-altitude precision
    let mut sin_lat = libm::sin(lat);
    let mut n = a / libm::sqrt(1.0 - e2 * sin_lat * sin_lat);
    let mut i = 0;
    while i < 10 {
        let lat_prev = lat;
        let height = if libm::fabs(lat) < core::f64::consts::FRAC_PI_4 {
            p / libm::cos(lat) - n
        } else {
            ecef.z / sin_lat - n * (1.0 - e2)
        };
        lat = libm::atan2(ecef.z, p * (1.0 - e2 * n / (n + height)));
        sin_lat = libm::sin(lat);
        n = a / libm::sqrt(1.0 - e2 * sin_lat * sin_lat);
        if libm::fabs(lat - lat_prev) < 1e-14 {
            break;
        }
        i += 1;
    }

    let height = if libm::fabs(lat) < core::f64::consts::FRAC_PI_4 {
        p / libm::cos(lat) - n
    } else {
        ecef.z / sin_lat - n * (1.0 - e2)
    };

    Vector3::new(lat, lon, height)
}

pub fn ecef_to_llh(ecef: Vector3<f64>) -> Vector3<f64> {
    ecef_to_llh_with_model(&Wgs84Model, ecef)
}

/// Converts Geodetic LLH (Latitude, Longitude, Height) in radians and meters to ECEF (meters).
pub fn llh_to_ecef_with_model<M: GeodeticModel>(model: &M, llh: Vector3<f64>) -> Vector3<f64> {
    let a = model.a();
    let e2 = model.e2();
    let lat = llh.x;
    let lon = llh.y;
    let height = llh.z;

    let sin_lat = libm::sin(lat);
    let cos_lat = libm::cos(lat);
    let sin_lon = libm::sin(lon);
    let cos_lon = libm::cos(lon);

    let n = a / libm::sqrt(1.0 - e2 * sin_lat * sin_lat);

    let x = (n + height) * cos_lat * cos_lon;
    let y = (n + height) * cos_lat * sin_lon;
    let z = (n * (1.0 - e2) + height) * sin_lat;

    Vector3::new(x, y, z)
}

pub fn llh_to_ecef(llh: Vector3<f64>) -> Vector3<f64> {
    llh_to_ecef_with_model(&Wgs84Model, llh)
}

/// Calculates Azimuth (radians) and Elevation (radians) of a satellite from a receiver's LLH and ECEF.
pub fn az_el(pos_llh: Vector3<f64>, pos_ecef: Vector3<f64>, sat_ecef: Vector3<f64>) -> (f64, f64) {
    let lat = pos_llh.x;
    let lon = pos_llh.y;

    let sin_lat = libm::sin(lat);
    let cos_lat = libm::cos(lat);
    let sin_lon = libm::sin(lon);
    let cos_lon = libm::cos(lon);

    // Vector from receiver to satellite in ECEF
    let dx = sat_ecef.x - pos_ecef.x;
    let dy = sat_ecef.y - pos_ecef.y;
    let dz = sat_ecef.z - pos_ecef.z;

    // Transform ECEF vector to ENU (East, North, Up) local tangent plane
    let e = -sin_lon * dx + cos_lon * dy;
    let n = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz;
    let u = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz;

    let horizontal_dist = libm::sqrt(e * e + n * n);
    let az = libm::fmod(libm::atan2(e, n) + 2.0 * core::f64::consts::PI, 2.0 * core::f64::consts::PI);
    let el = libm::atan2(u, horizontal_dist);

    (az, el)
}

/// Returns the 3x3 rotation matrix from ECEF to the local NED (North-East-Down) tangent plane.
pub fn ecef_to_ned_matrix(llh: Vector3<f64>) -> nalgebra::Matrix3<f64> {
    let lat = llh.x;
    let lon = llh.y;

    let sin_lat = libm::sin(lat);
    let cos_lat = libm::cos(lat);
    let sin_lon = libm::sin(lon);
    let cos_lon = libm::cos(lon);

    nalgebra::Matrix3::new(
        -sin_lat * cos_lon, -sin_lat * sin_lon,  cos_lat,
        -sin_lon,            cos_lon,           0.0,
        -cos_lat * cos_lon, -cos_lat * sin_lon, -sin_lat,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coordinate_alignment_validation() {
        let t1 = GpsTime::new(100, 100.0);
        let t2 = GpsTime::new(100, 100.0);
        let t3 = GpsTime::new(100, 101.0);

        let c1 = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, t1);
        
        // Exact match
        let c2 = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, t2);
        assert!(c1.ensure_aligned(&c2).is_ok());

        // Datum mismatch
        let c3 = Coordinate::new(Vector3::zeros(), Datum::PZ90, Frame::ECEF, t1);
        assert_eq!(c1.ensure_aligned(&c3).unwrap_err(), "Datum mismatch");

        // Frame mismatch
        let c4 = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ENU, t1);
        assert_eq!(c1.ensure_aligned(&c4).unwrap_err(), "Frame mismatch");

        // Epoch mismatch
        let c5 = Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, t3);
        assert_eq!(c1.ensure_aligned(&c5).unwrap_err(), "Epoch mismatch");
    }
}
