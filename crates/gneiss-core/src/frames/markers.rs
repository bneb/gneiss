//! Coordinate frame marker types and traits.

use super::realizations::ReferenceFrame;
use core::marker::PhantomData;

/// Marker trait for spatial coordinate reference frames.
pub trait CoordinateFrame: 'static + Send + Sync + Copy + PartialEq + Eq {
    /// Human-readable identifier for the coordinate frame.
    const NAME: &'static str;
    /// Indicates whether the coordinate frame is fixed to the rotating Earth.
    const IS_EARTH_FIXED: bool;
}

/// Earth-Centered, Earth-Fixed (ECEF) Cartesian coordinate frame tagged by reference realization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ecef<R: ReferenceFrame>(PhantomData<R>);

impl<R: ReferenceFrame> CoordinateFrame for Ecef<R> {
    const NAME: &'static str = R::NAME;
    const IS_EARTH_FIXED: bool = true;
}

/// Local North-East-Down (NED) tangent plane coordinate frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ned;

impl CoordinateFrame for Ned {
    const NAME: &'static str = "NED";
    const IS_EARTH_FIXED: bool = false;
}

/// Local East-North-Up (ENU) tangent plane coordinate frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Enu;

impl CoordinateFrame for Enu {
    const NAME: &'static str = "ENU";
    const IS_EARTH_FIXED: bool = false;
}

/// Vehicle or sensor Body frame: Forward (+X), Right (+Y), Down (+Z).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BodyFrd;

impl CoordinateFrame for BodyFrd {
    const NAME: &'static str = "Body-FRD";
    const IS_EARTH_FIXED: bool = false;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frames::realizations::Itrf2014;

    #[test]
    fn test_coordinate_frame_markers() {
        assert_eq!(<Ecef<Itrf2014> as CoordinateFrame>::NAME, "ITRF2014");
        let ecef_fixed = <Ecef<Itrf2014> as CoordinateFrame>::IS_EARTH_FIXED;
        assert!(ecef_fixed);

        assert_eq!(<Ned as CoordinateFrame>::NAME, "NED");
        let ned_fixed = <Ned as CoordinateFrame>::IS_EARTH_FIXED;
        assert!(!ned_fixed);

        assert_eq!(<Enu as CoordinateFrame>::NAME, "ENU");
        let enu_fixed = <Enu as CoordinateFrame>::IS_EARTH_FIXED;
        assert!(!enu_fixed);

        assert_eq!(<BodyFrd as CoordinateFrame>::NAME, "Body-FRD");
        let body_fixed = <BodyFrd as CoordinateFrame>::IS_EARTH_FIXED;
        assert!(!body_fixed);
    }
}
