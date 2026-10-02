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

#[cfg(test)]
mod tests {
    //! Golden-vector tests for the relational local tangent plane.
    //!
    //! The reference case is the equatorial prime-meridian origin
    //! (lat = 0, lon = 0, h = 0) whose ECEF position is exactly
    //! (a, 0, 0) = (6378137, 0, 0). There the ECEF->NED rotation collapses to the
    //! exact permutation
    //!     R_ned = [[ 0, 0, 1],     (North = ECEF z)
    //!             [ 0, 1, 0],     (East  = ECEF y)
    //!             [-1, 0, 0]]     (Down  = -ECEF x)
    //! because sin(0) = 0 and cos(0) = 1 exactly, so every expected NED/ENU
    //! component below is an exact integer permutation of the ECEF offset and
    //! catches any axis swap, sign flip, or NED<->ENU confusion.

    use super::*;
    use crate::frames::realizations::Itrf2014;
    use alloc::vec::Vec;
    use core::f64::consts::{FRAC_PI_2, PI};
    use nalgebra::Vector3;

    const A: f64 = 6_378_137.0;

    fn equator_origin() -> LocalTangentPlane<Itrf2014> {
        LocalTangentPlane::from_origin(EcefPos::<Itrf2014>::new(Vector3::new(A, 0.0, 0.0)))
    }

    #[test]
    fn origin_llh_is_exactly_lat0_lon0_h0() {
        let plane = equator_origin();
        assert!(plane.origin_llh().x.abs() < 1e-12);
        assert!(plane.origin_llh().y.abs() < 1e-12);
        assert!(plane.origin_llh().z.abs() < 1e-6);
        assert_eq!(*plane.origin().coords(), Vector3::new(A, 0.0, 0.0));
    }

    #[test]
    fn ned_matrix_is_the_exact_permutation_at_the_equator() {
        let plane = equator_origin();
        let expect = Matrix3::new(0.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0);
        assert!((plane.r_ned() - expect).norm() < 1e-12, "got {:?}", plane.r_ned());
        // R_ned is orthogonal with det = +1 (a change of basis, not a mirror).
        assert!((plane.r_ned() * plane.r_ned().transpose() - Matrix3::identity()).norm() < 1e-12);
        assert!((plane.r_ned().determinant() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ecef_offset_maps_to_the_exact_ned_and_enu_components() {
        // d_ecef = (10, 20, 30):
        //   N =  30 (ECEF z), E = 20 (ECEF y), D = -10 (ECEF x)
        //   ENU = (E, N, U) = (20, 30, +10)   because U = -D
        let plane = equator_origin();
        let target = EcefPos::<Itrf2014>::new(Vector3::new(A + 10.0, 20.0, 30.0));
        let ned = plane.to_ned(&target);
        assert!((ned.x() - 30.0).abs() < 1e-9, "N = {}", ned.x());
        assert!((ned.y() - 20.0).abs() < 1e-9, "E = {}", ned.y());
        assert!((ned.z() + 10.0).abs() < 1e-9, "D = {}", ned.z());
        let enu = plane.to_enu(&target);
        assert!((enu.x() - 20.0).abs() < 1e-9, "E = {}", enu.x());
        assert!((enu.y() - 30.0).abs() < 1e-9, "N = {}", enu.y());
        assert!((enu.z() - 10.0).abs() < 1e-9, "U = {}", enu.z());
    }

    #[test]
    fn origin_maps_to_the_zero_offset() {
        let plane = equator_origin();
        let o = *plane.origin();
        assert!(plane.to_ned(&o).norm() < 1e-9);
        assert!(plane.to_enu(&o).norm() < 1e-9);
    }

    #[test]
    fn ned_and_enu_round_trips_return_the_input_position_to_nanometre() {
        // For any origin, R_ned is orthogonal, so R_ned^T * (R_ned * d) = d exactly
        // up to rounding; the round trip must return the absolute position.
        let origins = [
            Vector3::new(A, 0.0, 0.0),
            Vector3::new(4_027_893.0, 307_041.0, 4_919_475.0),
            Vector3::new(-1_200_000.0, -4_300_000.0, 3_900_000.0),
        ];
        let mut deltas: Vec<Vector3<f64>> = Vec::new();
        for k in 0..3 {
            deltas.push(Vector3::new(10.0 * (k as f64 + 1.0), -250.0, 33.0));
        }
        for o in origins {
            let plane = LocalTangentPlane::<Itrf2014>::from_origin(EcefPos::new(o));
            for d in &deltas {
                let target = EcefPos::<Itrf2014>::new(o + d);
                let via_ned = plane.from_ned(&plane.to_ned(&target));
                assert!((via_ned.into_vector() - (o + d)).norm() < 1e-9);
                let via_enu = plane.from_enu(&plane.to_enu(&target));
                assert!((via_enu.into_vector() - (o + d)).norm() < 1e-9);
            }
        }
    }

    #[test]
    fn project_cov_to_ned_permutes_the_ecef_axes_at_the_equator() {
        // C_ned = R P R^T with R the permutation above and P = diag(p1,p2,p3)
        // gives C_ned = diag(p3, p2, p1): North<-ECEF z, East<-ECEF y,
        // Down<-ECEF x (the sign of the Down axis does not change a variance).
        let plane = equator_origin();
        let p = EcefCovariance::<Itrf2014>::from_diagonal(Vector3::new(1.0, 2.0, 4.0));
        let c = plane.project_cov_to_ned(&p);
        assert!((c.matrix() - Matrix3::from_diagonal(&Vector3::new(4.0, 2.0, 1.0))).norm() < 1e-12,
            "got {:?}", c.matrix());
        assert!((c.std_north() - 2.0).abs() < 1e-12);
        assert!((c.std_east() - 2.0_f64.sqrt()).abs() < 1e-12);
        assert!((c.std_down() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn project_cov_to_enu_is_the_ned_projection_reindexed_with_flipped_cross_terms() {
        // With C_ned = [[c,q,-r],[q,b,-p],[-r,-p,a]] (derived from a dense ECEF P
        // at the equator) the implementation's ENU remap yields
        //   [[b, q, p],
        //    [q, c, r],
        //    [p, r, a]]
        // i.e. a pure (x,y,z)_ECEF -> (N,E,U) permutation of the entries, with
        // both cross-term signs flipped because Up = -Down.
        let plane = equator_origin();
        let (a, b, c, p, q, r) = (1.0, 2.0, 3.0, 0.1, 0.2, 0.3);
        let ecef = EcefCovariance::<Itrf2014>::new(Matrix3::new(
            a, p, r,
            p, b, q,
            r, q, c,
        ));
        let got = plane.project_cov_to_enu(&ecef);
        let expect = Matrix3::new(b, q, p, q, c, r, p, r, a);
        assert!((got.matrix() - expect).norm() < 1e-12, "got {:?}", got.matrix());
        assert!((got.matrix() - got.matrix().transpose()).norm() < 1e-12);
        assert!((got.std_east() - b.sqrt()).abs() < 1e-12);
        assert!((got.std_north() - c.sqrt()).abs() < 1e-12);
        assert!((got.std_up() - a.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn az_el_cardinal_directions_at_the_equator() {
        // ENU axis targets, all at 1000 m from the equatorial origin:
        //   +ECEF x is Up    -> az = atan2(0,0) = 0, el = atan2(1000,0) = pi/2
        //   +ECEF y is East  -> az = atan2(1000,0) = pi/2, el = atan2(0,1000) = 0
        //   +ECEF z is North -> az = atan2(0,1000) = 0, el = atan2(0,1000) = 0
        let plane = equator_origin();
        let cases = [
            (Vector3::new(1000.0, 0.0, 0.0), 0.0, FRAC_PI_2),
            (Vector3::new(0.0, 1000.0, 0.0), FRAC_PI_2, 0.0),
            (Vector3::new(0.0, 0.0, 1000.0), 0.0, 0.0),
        ];
        for (d, az_exp, el_exp) in cases {
            let sat = EcefPos::<Itrf2014>::new(Vector3::new(A, 0.0, 0.0) + d);
            let (az, el) = plane.az_el(&sat);
            assert!((az - az_exp).abs() < 1e-12, "az for {d:?} = {az}, want {az_exp}");
            assert!((el - el_exp).abs() < 1e-12, "el for {d:?} = {el}, want {el_exp}");
        }
    }

    #[test]
    fn az_el_forty_five_degree_diagonal_north_east() {
        // ENU = (1000, 1000, 1000): horizontal = sqrt(1e6+1e6) = 1000*sqrt(2).
        //   az = atan2(1000, 1000) = pi/4
        //   el = atan2(1000, 1000*sqrt(2)) = atan(1/sqrt(2)) = 0.61547970867 rad (35.264 deg)
        let plane = equator_origin();
        let sat = EcefPos::<Itrf2014>::new(Vector3::new(A + 1000.0, 1000.0, 1000.0));
        let (az, el) = plane.az_el(&sat);
        assert!((az - PI / 4.0).abs() < 1e-12, "az = {az}");
        let expect_el = libm::atan(1.0 / 2.0_f64.sqrt());
        assert!((el - expect_el).abs() < 1e-12, "el = {el}, want {expect_el}");
    }

    #[test]
    fn negative_azimuth_is_wrapped_into_zero_to_two_pi() {
        // A satellite due West gives atan2(-1, 0) = -pi/2, which must wrap to
        // 3*pi/2 rather than stay negative.
        let plane = equator_origin();
        let sat = EcefPos::<Itrf2014>::new(Vector3::new(A, -1000.0, 0.0));
        let (az, el) = plane.az_el(&sat);
        assert!((az - 1.5 * PI).abs() < 1e-12, "az = {az}");
        assert!(el.abs() < 1e-12);
    }
}
