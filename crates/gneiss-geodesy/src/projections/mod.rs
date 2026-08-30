//! Geodetic map projections for converting ellipsoidal coordinates (lat, lon) to planar grid (Easting, Northing).

pub mod lambert_conformal;
pub mod transverse_mercator;

pub use lambert_conformal::LambertConformalConic;
pub use transverse_mercator::{Ellipsoid, TransverseMercator};
