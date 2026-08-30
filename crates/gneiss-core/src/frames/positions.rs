//! Strongly-typed ECEF coordinate containers and antenna point references.

use core::marker::PhantomData;
use nalgebra::Vector3;
use super::helmert::HelmertParams;
use super::realizations::ReferenceFrame;

/// Frame-tagged ECEF coordinates preventing accidental cross-datum mixing.
pub struct EcefPos<F: ReferenceFrame>(pub Vector3<f64>, pub PhantomData<F>);

impl<F: ReferenceFrame> EcefPos<F> {
    #[must_use]
    pub fn new(v: Vector3<f64>) -> Self {
        Self(v, PhantomData)
    }

    #[must_use]
    pub fn norm(&self) -> f64 {
        self.0.norm()
    }

    #[must_use]
    pub const fn vector(&self) -> &Vector3<f64> {
        &self.0
    }

    #[must_use]
    pub const fn into_vector(self) -> Vector3<f64> {
        self.0
    }

    pub fn convert_to<F2: ReferenceFrame>(&self, t_epoch_yr: f64) -> EcefPos<F2> {
        let to_hub = params_at(F::HELMERT_TO_ITRF2014, t_epoch_yr);
        let from_hub = params_at(F2::HELMERT_TO_ITRF2014, t_epoch_yr);
        EcefPos::new(from_hub.apply_inverse(to_hub.apply(self.0)))
    }
}

fn params_at(p: Option<HelmertParams>, t_yr: f64) -> HelmertParams {
    p.map_or_else(|| HelmertParams::identity_at(t_yr), |params| params.at(t_yr))
}

impl<F: ReferenceFrame> core::ops::Deref for EcefPos<F> {
    type Target = Vector3<f64>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<F: ReferenceFrame> From<Vector3<f64>> for EcefPos<F> {
    fn from(v: Vector3<f64>) -> Self {
        Self::new(v)
    }
}

impl<F: ReferenceFrame> From<EcefPos<F>> for Vector3<f64> {
    fn from(p: EcefPos<F>) -> Self {
        p.0
    }
}

impl<F: ReferenceFrame> Clone for EcefPos<F> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<F: ReferenceFrame> core::fmt::Debug for EcefPos<F> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "EcefPos<{}>({})", F::NAME, self.0)
    }
}

impl<F: ReferenceFrame> Copy for EcefPos<F> {}

impl<F: ReferenceFrame> PartialEq for EcefPos<F> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

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
            self.pos.0 + v * dt
        } else {
            self.pos.0
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
            pos: EcefPos::new(self.pos.0 + pco_ecef),
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
            pos: EcefPos::new(self.pos.0 - pco_ecef),
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
