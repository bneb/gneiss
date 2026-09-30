//! Relational local tangent plane derived from a single reference position.

use super::markers::{Enu, Ned};
use super::primitives::{EcefCovariance, EcefPos, EnuCovariance, NedCovariance, SpatialVector};
use super::realizations::ReferenceFrame;
use crate::coords::{ecef_to_llh, ecef_to_ned_matrix};
use nalgebra::{Matrix3, Vector3};

/// Local tangent plane constructed from a single reference ECEF position.
///
/// Disagreement between origin ECEF and LLH is unrepresentable because all
/// transformation matrices are derived directly from the single typed origin.
///
/// Attempting to project an ECEF position from a different datum fails at compile time:
/// ```compile_fail
/// use gneiss_core::frames::{EcefPos, Itrf2014, LocalTangentPlane, Nad83_2011};
/// use nalgebra::Vector3;
/// let origin = EcefPos::<Itrf2014>::new(Vector3::new(4e6, 3e5, 5e6));
/// let plane = LocalTangentPlane::from_origin(origin);
/// let target_nad = EcefPos::<Nad83_2011>::new(Vector3::new(4e6, 3e5, 5e6));
/// let _ = plane.to_enu(&target_nad); // Compile error: mismatched datums
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocalTangentPlane<F: ReferenceFrame> {
    origin: EcefPos<F>,
    origin_llh: Vector3<f64>,
    r_ned: Matrix3<f64>,
}

impl<F: ReferenceFrame> LocalTangentPlane<F> {
    /// Constructs a local tangent plane anchored at the specified origin.
    #[must_use]
    pub fn from_origin(origin: EcefPos<F>) -> Self {
        let origin_llh = ecef_to_llh(*origin.coords());
        let r_ned = ecef_to_ned_matrix(origin_llh);
        Self {
            origin,
            origin_llh,
            r_ned,
        }
    }

    /// Reference origin position in ECEF.
    #[inline]
    #[must_use]
    pub const fn origin(&self) -> &EcefPos<F> {
        &self.origin
    }

    /// Origin geodetic coordinates (latitude [rad], longitude [rad], height [m]).
    #[inline]
    #[must_use]
    pub const fn origin_llh(&self) -> &Vector3<f64> {
        &self.origin_llh
    }

    /// Rotation matrix from ECEF to NED local tangent plane.
    #[inline]
    #[must_use]
    pub const fn r_ned(&self) -> &Matrix3<f64> {
        &self.r_ned
    }

    /// Converts target ECEF position to an ENU (East, North, Up) offset vector.
    #[must_use]
    pub fn to_enu(&self, target: &EcefPos<F>) -> SpatialVector<Enu> {
        let d_ecef = target.coords() - self.origin.coords();
        let d_ned = self.r_ned * d_ecef;
        SpatialVector::new(d_ned.y, d_ned.x, -d_ned.z)
    }

    /// Converts target ECEF position to a NED (North, East, Down) offset vector.
    #[must_use]
    pub fn to_ned(&self, target: &EcefPos<F>) -> SpatialVector<Ned> {
        let d_ecef = target.coords() - self.origin.coords();
        let d_ned = self.r_ned * d_ecef;
        SpatialVector::from_vector(d_ned)
    }

    /// Converts an ENU offset vector back to an absolute ECEF position.
    #[must_use]
    pub fn from_enu(&self, enu: &SpatialVector<Enu>) -> EcefPos<F> {
        let v = enu.vector();
        let ned = Vector3::new(v.y, v.x, -v.z);
        let d_ecef = self.r_ned.transpose() * ned;
        EcefPos::new(self.origin.coords() + d_ecef)
    }

    /// Converts a NED offset vector back to an absolute ECEF position.
    #[must_use]
    pub fn from_ned(&self, ned: &SpatialVector<Ned>) -> EcefPos<F> {
        let d_ecef = self.r_ned.transpose() * ned.vector();
        EcefPos::new(self.origin.coords() + d_ecef)
    }

    /// Projects an ECEF error covariance into the local ENU frame.
    #[must_use]
    pub fn project_cov_to_enu(&self, cov: &EcefCovariance<F>) -> EnuCovariance {
        let c = self.r_ned * cov.matrix() * self.r_ned.transpose();
        let mut m = Matrix3::zeros();
        m[(0, 0)] = c[(1, 1)]; // East
        m[(1, 1)] = c[(0, 0)]; // North
        m[(2, 2)] = c[(2, 2)]; // Up
        m[(0, 1)] = c[(1, 0)];
        m[(1, 0)] = c[(0, 1)];
        m[(0, 2)] = -c[(1, 2)];
        m[(2, 0)] = -c[(2, 1)];
        m[(1, 2)] = -c[(0, 2)];
        m[(2, 1)] = -c[(2, 0)];
        EnuCovariance::new(m)
    }

    /// Projects an ECEF error covariance into the local NED frame.
    #[must_use]
    pub fn project_cov_to_ned(&self, cov: &EcefCovariance<F>) -> NedCovariance {
        let c = self.r_ned * cov.matrix() * self.r_ned.transpose();
        NedCovariance::new(c)
    }

    /// Computes Azimuth and Elevation (radians) of a satellite from this origin.
    #[must_use]
    pub fn az_el(&self, sat: &EcefPos<F>) -> (f64, f64) {
        let enu = self.to_enu(sat);
        let e = enu.x();
        let n = enu.y();
        let u = enu.z();
        let h = libm::sqrt(e * e + n * n);
        let two_pi = 2.0 * core::f64::consts::PI;
        let az = libm::fmod(libm::atan2(e, n) + two_pi, two_pi);
        let el = libm::atan2(u, h);
        (az, el)
    }
}
