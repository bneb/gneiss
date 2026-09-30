//! Strongly-typed antenna point references and epoch-qualified geodetic positions.

use core::marker::PhantomData;
use nalgebra::Vector3;
use super::primitives::EcefPos;
use super::realizations::ReferenceFrame;

// Re-export EcefPos for callers importing from positions
pub use super::primitives::EcefPos as EcefPosition;

/// Marker trait for physical points of reference on a station.
pub trait AntennaReference: 'static {
    const NAME: &'static str;
}

/// Physical Antenna Reference Point (bottom of antenna mount / ARP).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Arp;
impl AntennaReference for Arp {
    const NAME: &'static str = "ARP";
}

/// Electromagnetic Antenna Phase Center for a specific frequency band.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Apc<const BAND: u8>;
impl<const BAND: u8> AntennaReference for Apc<BAND> {
    const NAME: &'static str = "APC";
}

/// Ground monument / survey nail (below tripod).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GroundMonument;
impl AntennaReference for GroundMonument {
    const NAME: &'static str = "Monument";
}

/// Fully-qualified geodetic position combining reference frame realization,
/// physical reference marker, epoch timestamp, and tectonic velocity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EpochPosition<F: ReferenceFrame, R: AntennaReference = Arp> {
    pub pos: EcefPos<F>,
    pub epoch_yr: f64,
    pub velocity_m_yr: Option<Vector3<f64>>,
    _marker: PhantomData<R>,
}

impl<F: ReferenceFrame, R: AntennaReference> EpochPosition<F, R> {
    pub const fn new(pos: EcefPos<F>, epoch_yr: f64, velocity_m_yr: Option<Vector3<f64>>) -> Self {
        Self {
            pos,
            epoch_yr,
            velocity_m_yr,
            _marker: PhantomData,
        }
    }

    #[must_use]
    pub fn at_epoch(&self, target_epoch_yr: f64) -> Self {
        let dt = target_epoch_yr - self.epoch_yr;
        let new_pos = if let Some(v) = self.velocity_m_yr {
            *self.pos.coords() + v * dt
        } else {
            *self.pos.coords()
        };
        Self {
            pos: EcefPos::new(new_pos),
            epoch_yr: target_epoch_yr,
            velocity_m_yr: self.velocity_m_yr,
            _marker: PhantomData,
        }
    }

    #[must_use]
    pub fn convert_frame<TargetFrame: ReferenceFrame>(&self) -> EpochPosition<TargetFrame, R> {
        let new_pos = self.pos.convert_to::<TargetFrame>(self.epoch_yr);
        EpochPosition {
            pos: new_pos,
            epoch_yr: self.epoch_yr,
            velocity_m_yr: self.velocity_m_yr,
            _marker: PhantomData,
        }
    }
}

impl<F: ReferenceFrame> EpochPosition<F, Arp> {
    pub fn to_apc<const BAND: u8>(
        &self,
        pco_neu_mm: Vector3<f64>,
        rx_llh_rad: Vector3<f64>,
    ) -> EpochPosition<F, Apc<BAND>> {
        let pco_ecef = neu_to_ecef(pco_neu_mm * 1e-3, rx_llh_rad);
        EpochPosition {
            pos: EcefPos::new(*self.pos.coords() + pco_ecef),
            epoch_yr: self.epoch_yr,
            velocity_m_yr: self.velocity_m_yr,
            _marker: PhantomData,
        }
    }
}

impl<F: ReferenceFrame, const BAND: u8> EpochPosition<F, Apc<BAND>> {
    pub fn to_arp(
        &self,
        pco_neu_mm: Vector3<f64>,
        rx_llh_rad: Vector3<f64>,
    ) -> EpochPosition<F, Arp> {
        let pco_ecef = neu_to_ecef(pco_neu_mm * 1e-3, rx_llh_rad);
        EpochPosition {
            pos: EcefPos::new(*self.pos.coords() - pco_ecef),
            epoch_yr: self.epoch_yr,
            velocity_m_yr: self.velocity_m_yr,
            _marker: PhantomData,
        }
    }
}

pub(crate) fn neu_to_ecef(neu: Vector3<f64>, llh_rad: Vector3<f64>) -> Vector3<f64> {
    let lat = llh_rad[0];
    let lon = llh_rad[1];
    let s_lat = libm::sin(lat);
    let c_lat = libm::cos(lat);
    let s_lon = libm::sin(lon);
    let c_lon = libm::cos(lon);

    let n = neu[0];
    let e = neu[1];
    let u = neu[2];

    let dx = -s_lat * c_lon * n - s_lon * e + c_lat * c_lon * u;
    let dy = -s_lat * s_lon * n + c_lon * e + c_lat * s_lon * u;
    let dz = c_lat * n + s_lat * u;

    Vector3::new(dx, dy, dz)
}
