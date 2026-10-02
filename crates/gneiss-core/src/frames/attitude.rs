//! Attitude direction cosine typestates and vehicle body antenna lever arms.

use super::markers::{BodyFrd, CoordinateFrame, Ecef, Ned};
use super::primitives::{SpatialCovariance, SpatialVector, SpatialVelocity};
use core::marker::PhantomData;
use nalgebra::{Matrix3, Rotation3, UnitQuaternion};

/// Antenna lever arm vector strictly typed in vehicle Body FRD (Forward, Right, Down) frame.
///
/// Unrotated addition to an ECEF position fails to compile:
/// ```compile_fail
/// use gneiss_core::frames::{AntennaLeverArm, EcefPos, Itrf2014};
/// use nalgebra::Vector3;
/// let pos = EcefPos::<Itrf2014>::new(Vector3::zeros());
/// let arm = AntennaLeverArm::new(0.5, 0.0, 0.0);
/// let _ = pos + arm; // Compile error: mismatched frames
/// ```
#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AntennaLeverArm(pub SpatialVector<BodyFrd>);

impl AntennaLeverArm {
    #[inline]
    #[must_use]
    pub const fn new(forward_m: f64, right_m: f64, down_m: f64) -> Self {
        Self(SpatialVector::new(forward_m, right_m, down_m))
    }

    #[inline]
    #[must_use]
    pub fn zero() -> Self {
        Self(SpatialVector::zero())
    }

    #[inline]
    #[must_use]
    pub const fn from_body_vector(v: SpatialVector<BodyFrd>) -> Self {
        Self(v)
    }

    #[inline]
    #[must_use]
    pub const fn as_body_vector(&self) -> &SpatialVector<BodyFrd> {
        &self.0
    }

    #[inline]
    #[must_use]
    pub const fn into_body_vector(self) -> SpatialVector<BodyFrd> {
        self.0
    }

    #[inline]
    #[must_use]
    pub fn forward(&self) -> f64 { self.0.x() }

    #[inline]
    #[must_use]
    pub fn right(&self) -> f64 { self.0.y() }

    #[inline]
    #[must_use]
    pub fn down(&self) -> f64 { self.0.z() }
}

impl From<SpatialVector<BodyFrd>> for AntennaLeverArm {
    #[inline]
    fn from(v: SpatialVector<BodyFrd>) -> Self { Self(v) }
}

impl From<AntennaLeverArm> for SpatialVector<BodyFrd> {
    #[inline]
    fn from(arm: AntennaLeverArm) -> Self { arm.0 }
}

/// Direction cosine / rotation attitude transforming vectors from `From` frame to `To` frame.
#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Attitude<From: CoordinateFrame, To: CoordinateFrame> {
    q: UnitQuaternion<f64>,
    _frames: PhantomData<(From, To)>,
}

impl<From: CoordinateFrame, To: CoordinateFrame> Attitude<From, To> {
    #[inline]
    #[must_use]
    pub const fn from_unit_quaternion(q: UnitQuaternion<f64>) -> Self {
        Self { q, _frames: PhantomData }
    }

    #[inline]
    #[must_use]
    pub fn from_rotation_matrix(r: &Matrix3<f64>) -> Self {
        let rot = Rotation3::from_matrix_unchecked(*r);
        Self::from_unit_quaternion(UnitQuaternion::from_rotation_matrix(&rot))
    }

    #[inline]
    #[must_use]
    pub fn identity() -> Self {
        Self::from_unit_quaternion(UnitQuaternion::identity())
    }

    #[inline]
    #[must_use]
    pub fn rotate(&self, v: &SpatialVector<From>) -> SpatialVector<To> {
        SpatialVector::from_vector(self.q * v.vector())
    }

    #[inline]
    #[must_use]
    pub fn rotate_vector(&self, v: SpatialVector<From>) -> SpatialVector<To> {
        self.rotate(&v)
    }

    #[inline]
    #[must_use]
    pub fn rotate_velocity(&self, v: SpatialVelocity<From>) -> SpatialVelocity<To> {
        SpatialVelocity::from_vector(self.q * v.vector())
    }

    #[inline]
    #[must_use]
    pub fn rotate_cov(&self, cov: &SpatialCovariance<From>) -> SpatialCovariance<To> {
        let r = self.q.to_rotation_matrix().into_inner();
        SpatialCovariance::new(r * cov.matrix() * r.transpose())
    }

    #[inline]
    #[must_use]
    pub fn inverse(&self) -> Attitude<To, From> {
        Attitude::from_unit_quaternion(self.q.inverse())
    }

    #[inline]
    #[must_use]
    pub const fn quaternion(&self) -> &UnitQuaternion<f64> {
        &self.q
    }

    #[inline]
    #[must_use]
    pub fn rotation_matrix(&self) -> Matrix3<f64> {
        self.q.to_rotation_matrix().into_inner()
    }
}

pub type BodyToEcef<R> = Attitude<BodyFrd, Ecef<R>>;
pub type EcefToBody<R> = Attitude<Ecef<R>, BodyFrd>;
pub type BodyToNed = Attitude<BodyFrd, Ned>;
pub type NedToBody = Attitude<Ned, BodyFrd>;
pub type NedToEcef<R> = Attitude<Ned, Ecef<R>>;
pub type EcefToNed<R> = Attitude<Ecef<R>, Ned>;

#[cfg(test)]
mod tests {
    //! Golden-vector tests for the attitude typestate.
    //!
    //! Every rotation expectation below is the closed-form right-handed rotation
    //! about a principal axis,
    //!   R_x(t): (x, y cos t - z sin t, y sin t + z cos t)
    //!   R_y(t): (x cos t + z sin t, y, -x sin t + z cos t)
    //!   R_z(t): (x cos t - y sin t, x sin t + y cos t, z)
    //! so a flipped sign, a transposed matrix, or an index swap fails loudly.

    use super::*;
    use crate::frames::realizations::Itrf2014;
    use alloc::format;
    use core::f64::consts::{FRAC_PI_2, PI};
    use nalgebra::{Unit, Vector3};

    type B2E = Attitude<BodyFrd, Ecef<Itrf2014>>;

    fn body(x: f64, y: f64, z: f64) -> SpatialVector<BodyFrd> {
        SpatialVector::new(x, y, z)
    }

    fn axis_rot(axis: Vector3<f64>, angle: f64) -> B2E {
        let unit = Unit::new_normalize(axis);
        B2E::from_unit_quaternion(UnitQuaternion::from_axis_angle(&unit, angle))
    }

    #[test]
    fn lever_arm_accessors_keep_frd_axis_order() {
        // BodyFrd is (Forward, Right, Down): x=forward, y=right, z=down.
        // Any permutation of these accessors is a real mounting error, so pin
        // each component against a distinct value.
        let arm = AntennaLeverArm::new(0.5, -0.2, 1.25);
        assert_eq!(arm.forward(), 0.5);
        assert_eq!(arm.right(), -0.2);
        assert_eq!(arm.down(), 1.25);
        let want = SpatialVector::new(0.5, -0.2, 1.25);
        assert_eq!(*arm.as_body_vector(), want);
        assert_eq!(arm.into_body_vector(), want);
    }

    #[test]
    fn zero_lever_arm_is_the_frd_origin() {
        let arm = AntennaLeverArm::zero();
        assert_eq!(*arm.as_body_vector(), SpatialVector::zero());
        assert_eq!(arm.forward(), 0.0);
        assert_eq!(arm.right(), 0.0);
        assert_eq!(arm.down(), 0.0);
    }

    #[test]
    fn lever_arm_from_and_into_body_vector_are_inverse() {
        let v = body(1.0, 2.0, 3.0);
        assert_eq!(AntennaLeverArm::from_body_vector(v).as_body_vector(), &v);
        assert_eq!(SpatialVector::from(AntennaLeverArm::from_body_vector(v)), v);
        assert_eq!(AntennaLeverArm::from(v).into_body_vector(), v);
    }

    #[test]
    fn identity_rotation_matrix_is_exactly_the_identity() {
        let r = Attitude::<BodyFrd, Ned>::identity().rotation_matrix();
        assert_eq!(r, Matrix3::<f64>::identity());
        // identity() rotates nothing.
        let att = Attitude::<BodyFrd, Ned>::identity();
        assert_eq!(att.rotate(&body(1.5, -2.5, 3.5)), SpatialVector::<Ned>::new(1.5, -2.5, 3.5));
    }

    #[test]
    fn rotation_by_pi_over_two_about_z_maps_x_to_y() {
        // R_z(pi/2): (1,0,0) -> (cos, sin, 0) = (0, 1, 0);
        //              (0,1,0) -> (-sin, cos, 0) = (-1, 0, 0);
        //              (0,0,1) -> (0, 0, 1)  (z is invariant about z).
        let att = axis_rot(Vector3::z(), FRAC_PI_2);
        let ex = att.rotate(&body(1.0, 0.0, 0.0));
        assert!(ex.x().abs() < 1e-15 && (ex.y() - 1.0).abs() < 1e-15 && ex.z().abs() < 1e-15,
            "R_z(90) x = ({}, {}, {})", ex.x(), ex.y(), ex.z());
        let ey = att.rotate(&body(0.0, 1.0, 0.0));
        assert!((ey.x() + 1.0).abs() < 1e-15 && ey.y().abs() < 1e-15 && ey.z().abs() < 1e-15,
            "R_z(90) y = ({}, {}, {})", ey.x(), ey.y(), ey.z());
        let ez = att.rotate(&body(0.0, 0.0, 1.0));
        assert!((ez.x()).abs() < 1e-15 && ez.y().abs() < 1e-15 && (ez.z() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn rotation_by_pi_about_z_negates_the_horizontal_plane() {
        // R_z(pi): cos(pi) = -1, sin(pi) = 0, so (x,y,z) -> (-x,-y,z).
        // NOTE: a half turn about z sends ECEF +X to ECEF -X, not to +Y;
        // mapping +X onto +Y takes a quarter turn (pi/2), not pi.
        let att = axis_rot(Vector3::z(), PI);
        let ex = att.rotate(&body(1.0, 0.0, 0.0));
        assert!((ex.x() + 1.0).abs() < 1e-15 && ex.y().abs() < 1e-15 && ex.z().abs() < 1e-15,
            "R_z(180) x = ({}, {}, {})", ex.x(), ex.y(), ex.z());
        let ey = att.rotate(&body(0.0, 1.0, 0.0));
        assert!(ey.x().abs() < 1e-15 && (ey.y() + 1.0).abs() < 1e-15 && ey.z().abs() < 1e-15,
            "R_z(180) y = ({}, {}, {})", ey.x(), ey.y(), ey.z());
        let ez = att.rotate(&body(0.0, 0.0, 1.0));
        assert!(ez.x().abs() < 1e-15 && ez.y().abs() < 1e-15 && (ez.z() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn rotation_about_x_leaves_x_untouched() {
        // R_x(t): x unchanged; (0,1,0) -> (0, cos t, sin t).
        let att = axis_rot(Vector3::x(), FRAC_PI_2);
        let v = att.rotate(&body(3.0, 0.0, 0.0));
        assert!((v.x() - 3.0).abs() < 1e-15 && v.y().abs() < 1e-15 && v.z().abs() < 1e-15);
        let v = att.rotate(&body(0.0, 1.0, 0.0));
        assert!(v.x().abs() < 1e-15 && v.y().abs() < 1e-15 && (v.z() - 1.0).abs() < 1e-15,
            "R_x(90) y = ({}, {}, {})", v.x(), v.y(), v.z());
    }

    #[test]
    fn rotation_about_y_leaves_y_untouched() {
        // R_y(t): y unchanged; (1,0,0) -> (cos t, 0, -sin t).
        let att = axis_rot(Vector3::y(), FRAC_PI_2);
        let v = att.rotate(&body(0.0, 7.0, 0.0));
        assert!(v.x().abs() < 1e-15 && (v.y() - 7.0).abs() < 1e-15 && v.z().abs() < 1e-15);
        let v = att.rotate(&body(1.0, 0.0, 0.0));
        assert!(v.x().abs() < 1e-15 && v.y().abs() < 1e-15 && (v.z() + 1.0).abs() < 1e-15,
            "R_y(90) x = ({}, {}, {})", v.x(), v.y(), v.z());
    }

    #[test]
    fn quaternion_is_unit_norm_and_orthonormal() {
        let att = axis_rot(Vector3::new(1.0, 2.0, 3.0).normalize(), 0.7);
        assert!((att.quaternion().norm() - 1.0).abs() < 1e-15);
        let r = att.rotation_matrix();
        assert!((r * r.transpose() - Matrix3::<f64>::identity()).norm() < 1e-12,
            "R^T R must equal I");
        // A proper rotation has determinant +1 (never -1: that would mirror).
        assert!((r.determinant() - 1.0).abs() < 1e-12, "det = {}", r.determinant());
    }

    #[test]
    fn inverse_attitude_undoes_the_rotation() {
        let att = axis_rot(Vector3::new(0.3, -0.5, 0.8).normalize(), 1.1);
        let inv: Attitude<Ecef<Itrf2014>, BodyFrd> = att.inverse();
        let v = body(4.0, -5.0, 6.0);
        assert!((inv.rotate(&att.rotate(&v)) - v).norm() < 1e-12);
        // Conjugate quaternion relation: q^-1 = (w, -x, -y, -z).
        let inv_att = att.inverse();
        let qi = inv_att.quaternion();
        assert!((qi.coords.x + att.quaternion().coords.x).abs() < 1e-15);
        assert!((qi.coords.y + att.quaternion().coords.y).abs() < 1e-15);
        assert!((qi.coords.z + att.quaternion().coords.z).abs() < 1e-15);
        assert!((qi.coords.w - att.quaternion().coords.w).abs() < 1e-15);
    }

    #[test]
    fn from_rotation_matrix_matches_the_axis_rotation() {
        let att = axis_rot(Vector3::new(1.0, 1.0, 1.0).normalize(), 0.9);
        let rebuilt = B2E::from_rotation_matrix(&att.rotation_matrix());
        let v = body(2.0, 3.0, 4.0);
        assert!((rebuilt.rotate(&v) - att.rotate(&v)).norm() < 1e-12);
    }

    #[test]
    fn rotate_covariance_is_r_p_r_transpose() {
        // R = R_z(90) = [[0,-1,0],[1,0,0],[0,0,1]], P = [[2,1,0],[1,3,1],[0,1,4]].
        //   R P   = [[-1,-3,-1],[2,1,0],[0,1,4]]
        //   R P R^T = [[3,-1,-1],[-1,2,0],[-1,0,4]]   (trace preserved: 9 = 9)
        // The transpose order R^T P R would give [[-3,1,1],[1,2,0],[1,0,4]], so
        // this pins the covariance rotation direction, not just symmetry.
        let att = axis_rot(Vector3::z(), FRAC_PI_2);
        let p = SpatialCovariance::<BodyFrd>::new(Matrix3::new(
            2.0, 1.0, 0.0,
            1.0, 3.0, 1.0,
            0.0, 1.0, 4.0,
        ));
        let got = att.rotate_cov(&p).into_matrix();
        let expect = Matrix3::new(
            3.0, -1.0, -1.0,
            -1.0, 2.0, 0.0,
            -1.0, 0.0, 4.0,
        );
        assert!((got - expect).norm() < 1e-12, "R P R^T = {got:?}");
        assert!((got - got.transpose()).norm() < 1e-12);
    }

    #[test]
    fn rotate_covariance_permutes_axes_without_changing_variance_magnitudes() {
        // Diagonal ECEF covariance diag(1,2,3) under R_z(90) must become
        // diag(2,1,3): ECEF x maps to ECEF y and vice versa, ECEF z is fixed.
        let att = axis_rot(Vector3::z(), FRAC_PI_2);
        let p = SpatialCovariance::<BodyFrd>::from_diagonal(Vector3::new(1.0, 2.0, 3.0));
        let got = att.rotate_cov(&p).into_matrix();
        let expect = Matrix3::from_diagonal(&Vector3::new(2.0, 1.0, 3.0));
        assert!((got - expect).norm() < 1e-12, "got {got:?}");
    }

    #[test]
    fn rotate_velocity_matches_rotate_for_the_same_quaternion() {
        // A velocity is a vector in the same frame, so the same rotation must
        // apply: R * v is the same triple whether carried as a vector or a velocity.
        let att = axis_rot(Vector3::new(-1.0, 0.5, 0.25).normalize(), 2.0);
        let v = SpatialVelocity::<BodyFrd>::new(11.0, -7.0, 3.5);
        let as_vector = att.rotate(&SpatialVector::<BodyFrd>::from(v.into_vector()));
        let as_velocity: SpatialVelocity<Ecef<Itrf2014>> = SpatialVelocity::from(as_vector.into_vector());
        assert!((att.rotate_velocity(v) - as_velocity).norm() < 1e-12);
    }

    #[test]
    fn rotate_vector_consumes_by_value_and_agrees_with_borrowed_form() {
        let att = axis_rot(Vector3::z(), 0.3);
        let v = body(1.0, 2.0, 3.0);
        assert_eq!(att.rotate_vector(v), att.rotate(&v));
    }

    #[test]
    fn ecef_body_round_trip_is_lossless() {
        let att = axis_rot(Vector3::new(0.0, 0.0, 1.0), 1.2);
        let forward_ecef = att.rotate(&body(1.0, 0.0, 0.0));
        let back = att.inverse().rotate(&forward_ecef);
        assert!((back - body(1.0, 0.0, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn exported_type_aliases_name_the_same_instantiations() {
        let _: EcefToNed<Itrf2014> = Attitude::identity();
        let _: NedToEcef<Itrf2014> = Attitude::identity();
        let _: BodyToNed = Attitude::identity();
        let _: NedToBody = Attitude::identity();
        let _: EcefToBody<Itrf2014> = Attitude::identity();
        let _: BodyToEcef<Itrf2014> = Attitude::identity();
        let _ = format!("{:?}", AntennaLeverArm::zero());
    }
}
