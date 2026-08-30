//! Frame-tagged ECEF positions: compile-time prevention of datum mismatches.
//!
//! Every [`EcefPos`] carries its reference realization (`ITRF2014`, `IGS20`,
//! broadcast `WGS84`, ...) as a type parameter. Mixing frames without an
//! explicit [`EcefPos::convert_to`] call is a **compile error**.

pub mod helmert;
pub mod realizations;
pub mod positions;
#[cfg(test)]
mod tests;

pub use helmert::{HelbertParams, HelmertParams};
pub use realizations::{
    Etrs89, Gda2020, Igs20, Itrf2014, Itrf2020, Jgd2011, Nad83_2011, ReferenceFrame,
    Wgs84Broadcast,
};
#[cfg(test)]
pub(crate) use realizations::ITRF2020_TO_ITRF2014;
pub use positions::{AntennaReference, Apc, Arp, EcefPos, EpochPosition, GroundMonument};
