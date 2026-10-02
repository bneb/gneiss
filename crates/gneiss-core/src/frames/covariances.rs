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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frames::realizations::Itrf2014;
    use crate::frames::Ecef;
    use alloc::format;

    #[test]
    fn from_stds_squares_each_standard_deviation() {
        // std_axis returns sqrt(variance): 2^2 = 4, 3^2 = 9, 4^2 = 16.
        let c = SpatialCovariance::<Ned>::from_stds(2.0, 3.0, 4.0);
        assert_eq!(*c.matrix(), Matrix3::from_diagonal(&Vector3::new(4.0, 9.0, 16.0)));
        assert!((c.std_north() - 2.0).abs() < 1e-15);
        assert!((c.std_east() - 3.0).abs() < 1e-15);
        assert!((c.std_down() - 4.0).abs() < 1e-15);
    }

    #[test]
    fn ned_and_enu_named_axes_are_permutations_of_the_same_diagonal() {
        // The three per-frame accessors must read axis 0/1/2 in the frame's own
        // documented order: NED = (N,E,D), ENU = (E,N,U).
        let c = SpatialCovariance::<Enu>::from_stds(10.0, 20.0, 30.0);
        assert!((c.std_east() - 10.0).abs() < 1e-15);
        assert!((c.std_north() - 20.0).abs() < 1e-15);
        assert!((c.std_up() - 30.0).abs() < 1e-15);
        // Swapping N/E or N/U in either frame would break this:
        let n = SpatialCovariance::<Ned>::from_diagonal(Vector3::new(100.0, 400.0, 900.0));
        assert!((n.std_north() - 10.0).abs() < 1e-15);
        assert!((n.std_east() - 20.0).abs() < 1e-15);
        assert!((n.std_down() - 30.0).abs() < 1e-15);
    }

    #[test]
    fn std_axis_clamps_negative_variance_to_zero() {
        // A round-off-negated diagonal must not produce NaN: sqrt(max(s2, 0)) = 0.
        let c = SpatialCovariance::<Ecef<Itrf2014>>::new(Matrix3::from_diagonal(&Vector3::new(
            -1e-18, 0.25, 4.0,
        )));
        assert_eq!(c.std_axis(0), 0.0);
        assert!((c.std_axis(1) - 0.5).abs() < 1e-15);
        assert!((c.std_axis(2) - 2.0).abs() < 1e-15);
    }

    #[test]
    fn matrix_accessors_agree_and_into_matrix_is_a_move() {
        let m = Matrix3::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
        let c = SpatialCovariance::<Ned>::new(m);
        assert_eq!(*c.matrix(), m);
        assert_eq!(c.into_matrix(), m);
    }

    #[test]
    fn debug_names_the_coordinate_frame() {
        let c = SpatialCovariance::<Enu>::from_stds(1.0, 1.0, 1.0);
        let s = format!("{c:?}");
        assert!(s.contains("SpatialCovariance<ENU>"), "got {s}");
        let e = SpatialCovariance::<Ecef<Itrf2014>>::from_stds(1.0, 1.0, 1.0);
        assert!(format!("{e:?}").contains("SpatialCovariance<ITRF2014>"));
    }
}
