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
