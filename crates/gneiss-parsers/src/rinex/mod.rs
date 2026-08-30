//! RINEX observation and navigation file parsing, split by RINEX file type:
//! `obs` handles OBS (measurement) files, `nav` handles NAV (broadcast
//! ephemeris) files. The two share no parsing state.

mod nav;
mod obs;

pub use nav::{parse_rinex_f64, parse_rinex_nav};
pub use obs::{parse_rinex_obs, parse_rinex_obs_epochs, RinexObsHeader};
