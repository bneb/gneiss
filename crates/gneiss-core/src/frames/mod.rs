//! Frame-tagged coordinates, spatial vectors, attitudes, and local tangent planes.
//!
//! Enforces zero-cost compile-time frame safety:
//! - Spatial coordinate markers: [`Ecef`], [`Ned`], [`Enu`], [`BodyFrd`].
//! - Reference frame realizations: [`Itrf2014`], [`Itrf2020`], [`Igs20`], [`Wgs84`], [`Nad83`], [`Jgd2011`], [`Pz90`].
//! - Spatial primitives: [`Point3`], [`SpatialVector`], [`SpatialVelocity`], [`SpatialCovariance`], [`AntennaLeverArm`], [`Attitude`].
//! - Relational local tangent plane: [`LocalTangentPlane`].

pub mod attitude;
pub mod covariances;
pub mod helmert;
pub mod markers;
pub mod positions;
pub mod primitives;
pub mod realizations;
pub mod tangent;

#[cfg(test)]
mod tests;

pub use helmert::HelmertParams;
pub use markers::{BodyFrd, CoordinateFrame, Ecef, Enu, Ned};
pub use positions::{AntennaReference, Apc, Arp, EpochPosition, GroundMonument};
pub use primitives::{
    AntennaLeverArm, Attitude, BodyCovariance, BodyToEcef, BodyToNed, BodyVector, BodyVelocity,
    EcefCovariance, EcefPos, EcefToBody, EcefToNed, EcefVector, EcefVelocity, EnuCovariance,
    EnuPos, EnuVector, EnuVelocity, NedCovariance, NedPos, NedToBody, NedToEcef, NedVector,
    NedVelocity, Point3, SpatialCovariance, SpatialVector, SpatialVelocity,
};
pub use realizations::{
    Etrs89, Gda2020, Igs20, Itrf2014, Itrf2020, Jgd2011, Nad83, Nad83_2011, Pz90, ReferenceFrame,
    Wgs84, Wgs84Broadcast,
};
pub use tangent::LocalTangentPlane;

#[cfg(test)]
pub(crate) use realizations::ITRF2020_TO_ITRF2014;
