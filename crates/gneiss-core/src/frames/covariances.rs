//! Spatial error covariance matrices with frame tagging.

use super::markers::{CoordinateFrame, Enu, Ned};
use core::marker::PhantomData;
use nalgebra::{Matrix3, Vector3};

/// 3x3 spatial error covariance matrix in coordinate frame `F`.
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq)]
pub struct SpatialCovariance<F: CoordinateFrame> {
    matrix: Matrix3<f64>,
    _frame: PhantomData<F>,
}

impl<F: CoordinateFrame> SpatialCovariance<F> {
    #[inline]
    #[must_use]
    pub const fn new(matrix: Matrix3<f64>) -> Self {
        Self { matrix, _frame: PhantomData }
    }

    #[inline]
    #[must_use]
    pub fn from_diagonal(diag: Vector3<f64>) -> Self {
        Self::new(Matrix3::from_diagonal(&diag))
    }

    #[inline]
    #[must_use]
    pub fn from_stds(sx: f64, sy: f64, sz: f64) -> Self {
        Self::from_diagonal(Vector3::new(sx * sx, sy * sy, sz * sz))
    }

    #[inline]
    #[must_use]
    pub const fn matrix(&self) -> &Matrix3<f64> { &self.matrix }

    #[inline]
    #[must_use]
    pub const fn into_matrix(self) -> Matrix3<f64> { self.matrix }

    #[inline]
    #[must_use]
    pub fn std_axis(&self, axis: usize) -> f64 {
        libm::sqrt(self.matrix[(axis, axis)].max(0.0))
    }
}

impl SpatialCovariance<Ned> {
    #[inline]
    #[must_use]
    pub fn std_north(&self) -> f64 { self.std_axis(0) }

    #[inline]
    #[must_use]
    pub fn std_east(&self) -> f64 { self.std_axis(1) }

    #[inline]
    #[must_use]
    pub fn std_down(&self) -> f64 { self.std_axis(2) }
}

impl SpatialCovariance<Enu> {
    #[inline]
    #[must_use]
    pub fn std_east(&self) -> f64 { self.std_axis(0) }

    #[inline]
    #[must_use]
    pub fn std_north(&self) -> f64 { self.std_axis(1) }

    #[inline]
    #[must_use]
    pub fn std_up(&self) -> f64 { self.std_axis(2) }
}

impl<F: CoordinateFrame> core::fmt::Debug for SpatialCovariance<F> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "SpatialCovariance<{}>({:?})", F::NAME, self.matrix)
    }
}
