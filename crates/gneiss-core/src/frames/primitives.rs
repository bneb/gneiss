//! Spatial points, vectors, velocities, covariances, lever arms, and attitudes.

use super::helmert::HelmertParams;
use super::markers::{BodyFrd, CoordinateFrame, Ecef, Enu, Ned};
use super::realizations::ReferenceFrame;
use core::marker::PhantomData;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use nalgebra::Vector3;

/// Affine 3D position point in coordinate frame `F`.
///
/// Addition of two `Point3` instances is structurally prohibited at compile time.
///
/// ```compile_fail
/// use gneiss_core::frames::{EcefPos, Itrf2014};
/// use nalgebra::Vector3;
/// let p1 = EcefPos::<Itrf2014>::new(Vector3::zeros());
/// let p2 = EcefPos::<Itrf2014>::new(Vector3::zeros());
/// let _ = p1 + p2; // Compile error: cannot add two points
/// ```
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq)]
pub struct Point3<F: CoordinateFrame> {
    coords: Vector3<f64>,
    _frame: PhantomData<F>,
}

impl<F: CoordinateFrame> Point3<F> {
    #[inline]
    #[must_use]
    pub const fn new(coords: Vector3<f64>) -> Self {
        Self { coords, _frame: PhantomData }
    }

    #[inline]
    #[must_use]
    pub const fn from_coords(x: f64, y: f64, z: f64) -> Self {
        Self::new(Vector3::new(x, y, z))
    }

    #[inline]
    #[must_use]
    pub const fn coords(&self) -> &Vector3<f64> {
        &self.coords
    }

    #[inline]
    #[must_use]
    pub const fn vector(&self) -> &Vector3<f64> {
        &self.coords
    }

    #[inline]
    #[must_use]
    pub const fn into_coords(self) -> Vector3<f64> {
        self.coords
    }

    #[inline]
    #[must_use]
    pub const fn into_vector(self) -> Vector3<f64> {
        self.coords
    }

    #[inline]
    #[must_use]
    pub fn x(&self) -> f64 { self.coords.x }

    #[inline]
    #[must_use]
    pub fn y(&self) -> f64 { self.coords.y }

    #[inline]
    #[must_use]
    pub fn z(&self) -> f64 { self.coords.z }

    #[inline]
    #[must_use]
    pub fn norm(&self) -> f64 { self.coords.norm() }

    #[inline]
    #[must_use]
    pub fn distance_to(&self, other: &Self) -> f64 {
        (self.coords - other.coords).norm()
    }
}

impl<R: ReferenceFrame> Point3<Ecef<R>> {
    #[must_use]
    pub fn convert_to<R2: ReferenceFrame>(&self, t_epoch_yr: f64) -> Point3<Ecef<R2>> {
        let to_hub = params_at(R::HELMERT_TO_ITRF2014, t_epoch_yr);
        let from_hub = params_at(R2::HELMERT_TO_ITRF2014, t_epoch_yr);
        Point3::new(from_hub.apply_inverse(to_hub.apply(self.coords)))
    }
}

fn params_at(p: Option<HelmertParams>, t_yr: f64) -> HelmertParams {
    p.map_or_else(|| HelmertParams::identity_at(t_yr), |params| params.at(t_yr))
}

impl<F: CoordinateFrame> core::fmt::Debug for Point3<F> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "Point3<{}>([{}, {}, {}])",
            F::NAME,
            self.coords.x,
            self.coords.y,
            self.coords.z
        )
    }
}

impl<F: CoordinateFrame> From<Vector3<f64>> for Point3<F> {
    #[inline]
    fn from(v: Vector3<f64>) -> Self { Self::new(v) }
}

impl<F: CoordinateFrame> From<Point3<F>> for Vector3<f64> {
    #[inline]
    fn from(p: Point3<F>) -> Self { p.coords }
}

impl<F: CoordinateFrame> Sub for Point3<F> {
    type Output = SpatialVector<F>;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        SpatialVector::from_vector(self.coords - rhs.coords)
    }
}

impl<F: CoordinateFrame> Add<SpatialVector<F>> for Point3<F> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: SpatialVector<F>) -> Self {
        Self::new(self.coords + rhs.vector)
    }
}

impl<F: CoordinateFrame> Sub<SpatialVector<F>> for Point3<F> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: SpatialVector<F>) -> Self {
        Self::new(self.coords - rhs.vector)
    }
}

impl<F: CoordinateFrame> AddAssign<SpatialVector<F>> for Point3<F> {
    #[inline]
    fn add_assign(&mut self, rhs: SpatialVector<F>) {
        self.coords += rhs.vector;
    }
}

impl<F: CoordinateFrame> SubAssign<SpatialVector<F>> for Point3<F> {
    #[inline]
    fn sub_assign(&mut self, rhs: SpatialVector<F>) {
        self.coords -= rhs.vector;
    }
}

/// 3D spatial displacement or direction vector in coordinate frame `F`.
///
/// Cross-frame vector operations are strictly prohibited at compile time.
///
/// ```compile_fail
/// use gneiss_core::frames::{BodyFrd, Ecef, Itrf2014, SpatialVector};
/// let ecef_v = SpatialVector::<Ecef<Itrf2014>>::zero();
/// let body_v = SpatialVector::<BodyFrd>::zero();
/// let _ = ecef_v + body_v; // Compile error: mismatched frames
/// ```
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq)]
pub struct SpatialVector<F: CoordinateFrame> {
    vector: Vector3<f64>,
    _frame: PhantomData<F>,
}

impl<F: CoordinateFrame> SpatialVector<F> {
    #[inline]
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { vector: Vector3::new(x, y, z), _frame: PhantomData }
    }

    #[inline]
    #[must_use]
    pub const fn from_vector(vector: Vector3<f64>) -> Self {
        Self { vector, _frame: PhantomData }
    }

    #[inline]
    #[must_use]
    pub fn zero() -> Self {
        Self::from_vector(Vector3::zeros())
    }

    #[inline]
    #[must_use]
    pub const fn vector(&self) -> &Vector3<f64> { &self.vector }

    #[inline]
    #[must_use]
    pub const fn coords(&self) -> &Vector3<f64> { &self.vector }

    #[inline]
    #[must_use]
    pub const fn into_vector(self) -> Vector3<f64> { self.vector }

    #[inline]
    #[must_use]
    pub const fn into_coords(self) -> Vector3<f64> { self.vector }

    #[inline]
    #[must_use]
    pub fn x(&self) -> f64 { self.vector.x }

    #[inline]
    #[must_use]
    pub fn y(&self) -> f64 { self.vector.y }

    #[inline]
    #[must_use]
    pub fn z(&self) -> f64 { self.vector.z }

    #[inline]
    #[must_use]
    pub fn norm(&self) -> f64 { self.vector.norm() }

    #[inline]
    #[must_use]
    pub fn dot(&self, other: &Self) -> f64 { self.vector.dot(&other.vector) }

    #[inline]
    #[must_use]
    pub fn cross(&self, other: &Self) -> Self {
        Self::from_vector(self.vector.cross(&other.vector))
    }
}

impl<F: CoordinateFrame> core::fmt::Debug for SpatialVector<F> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "SpatialVector<{}>([{}, {}, {}])", F::NAME, self.vector.x, self.vector.y, self.vector.z)
    }
}

impl<F: CoordinateFrame> From<Vector3<f64>> for SpatialVector<F> {
    #[inline]
    fn from(v: Vector3<f64>) -> Self { Self::from_vector(v) }
}

impl<F: CoordinateFrame> From<SpatialVector<F>> for Vector3<f64> {
    #[inline]
    fn from(v: SpatialVector<F>) -> Self { v.vector }
}

impl<F: CoordinateFrame> Add for SpatialVector<F> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self { Self::from_vector(self.vector + rhs.vector) }
}

impl<F: CoordinateFrame> Sub for SpatialVector<F> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self { Self::from_vector(self.vector - rhs.vector) }
}

impl<F: CoordinateFrame> Neg for SpatialVector<F> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self { Self::from_vector(-self.vector) }
}

impl<F: CoordinateFrame> Mul<f64> for SpatialVector<F> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f64) -> Self { Self::from_vector(self.vector * rhs) }
}

impl<F: CoordinateFrame> Mul<SpatialVector<F>> for f64 {
    type Output = SpatialVector<F>;
    #[inline]
    fn mul(self, rhs: SpatialVector<F>) -> SpatialVector<F> { rhs * self }
}

impl<F: CoordinateFrame> Div<f64> for SpatialVector<F> {
    type Output = Self;
    #[inline]
    fn div(self, rhs: f64) -> Self { Self::from_vector(self.vector / rhs) }
}

impl<F: CoordinateFrame> AddAssign for SpatialVector<F> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) { self.vector += rhs.vector; }
}

impl<F: CoordinateFrame> SubAssign for SpatialVector<F> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) { self.vector -= rhs.vector; }
}

impl<F: CoordinateFrame> MulAssign<f64> for SpatialVector<F> {
    #[inline]
    fn mul_assign(&mut self, rhs: f64) { self.vector *= rhs; }
}

impl<F: CoordinateFrame> DivAssign<f64> for SpatialVector<F> {
    #[inline]
    fn div_assign(&mut self, rhs: f64) { self.vector /= rhs; }
}

/// 3D spatial velocity in coordinate frame `F` (m/s).
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq)]
pub struct SpatialVelocity<F: CoordinateFrame> {
    vector: Vector3<f64>,
    _frame: PhantomData<F>,
}

impl<F: CoordinateFrame> SpatialVelocity<F> {
    #[inline]
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { vector: Vector3::new(x, y, z), _frame: PhantomData }
    }

    #[inline]
    #[must_use]
    pub const fn from_vector(vector: Vector3<f64>) -> Self {
        Self { vector, _frame: PhantomData }
    }

    #[inline]
    #[must_use]
    pub fn zero() -> Self { Self::from_vector(Vector3::zeros()) }

    #[inline]
    #[must_use]
    pub const fn vector(&self) -> &Vector3<f64> { &self.vector }

    #[inline]
    #[must_use]
    pub const fn into_vector(self) -> Vector3<f64> { self.vector }

    #[inline]
    #[must_use]
    pub fn x(&self) -> f64 { self.vector.x }

    #[inline]
    #[must_use]
    pub fn y(&self) -> f64 { self.vector.y }

    #[inline]
    #[must_use]
    pub fn z(&self) -> f64 { self.vector.z }

    #[inline]
    #[must_use]
    pub fn norm(&self) -> f64 { self.vector.norm() }

    #[inline]
    #[must_use]
    pub fn displacement(&self, dt_seconds: f64) -> SpatialVector<F> {
        SpatialVector::from_vector(self.vector * dt_seconds)
    }
}

impl<F: CoordinateFrame> core::fmt::Debug for SpatialVelocity<F> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "SpatialVelocity<{}>([{}, {}, {}])", F::NAME, self.vector.x, self.vector.y, self.vector.z)
    }
}

impl<F: CoordinateFrame> From<Vector3<f64>> for SpatialVelocity<F> {
    #[inline]
    fn from(v: Vector3<f64>) -> Self { Self::from_vector(v) }
}

impl<F: CoordinateFrame> From<SpatialVelocity<F>> for Vector3<f64> {
    #[inline]
    fn from(v: SpatialVelocity<F>) -> Self { v.vector }
}

impl<F: CoordinateFrame> Add for SpatialVelocity<F> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self { Self::from_vector(self.vector + rhs.vector) }
}

impl<F: CoordinateFrame> Sub for SpatialVelocity<F> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self { Self::from_vector(self.vector - rhs.vector) }
}

impl<F: CoordinateFrame> Neg for SpatialVelocity<F> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self { Self::from_vector(-self.vector) }
}

impl<F: CoordinateFrame> Mul<f64> for SpatialVelocity<F> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f64) -> Self { Self::from_vector(self.vector * rhs) }
}

impl<F: CoordinateFrame> Div<f64> for SpatialVelocity<F> {
    type Output = Self;
    #[inline]
    fn div(self, rhs: f64) -> Self { Self::from_vector(self.vector / rhs) }
}

pub use super::covariances::SpatialCovariance;

pub use super::attitude::{
    AntennaLeverArm, Attitude, BodyToEcef, BodyToNed, EcefToBody, EcefToNed, NedToBody, NedToEcef,
};

// Convenient type aliases
pub type EcefPos<R> = Point3<Ecef<R>>;
pub type NedPos = Point3<Ned>;
pub type EnuPos = Point3<Enu>;

pub type EcefVector<R> = SpatialVector<Ecef<R>>;
pub type NedVector = SpatialVector<Ned>;
pub type EnuVector = SpatialVector<Enu>;
pub type BodyVector = SpatialVector<BodyFrd>;

pub type EcefVelocity<R> = SpatialVelocity<Ecef<R>>;
pub type NedVelocity = SpatialVelocity<Ned>;
pub type EnuVelocity = SpatialVelocity<Enu>;
pub type BodyVelocity = SpatialVelocity<BodyFrd>;

pub type EcefCovariance<R> = SpatialCovariance<Ecef<R>>;
pub type NedCovariance = SpatialCovariance<Ned>;
pub type EnuCovariance = SpatialCovariance<Enu>;
pub type BodyCovariance = SpatialCovariance<BodyFrd>;

// Tests live in a sibling file because this module is already at the 500-LOC
// ceiling; `mod tests` is still compiled only under cfg(test).
#[cfg(test)]
mod tests;

