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

#[cfg(test)]
mod tests {
    //! Golden vectors for the antenna-point typestate.
    //!
    //! `neu_to_ecef` is the transpose of the ECEF->NED basis, so at two
    //! analytically exact sites it collapses to a pure axis permutation:
    //!   lat = 0, lon = 0   (sin = 0, cos = 1): NED (n,e,u) -> ECEF (u, e, n)
    //!   lat = 0, lon = 90  (s_lat 0, c_lat 1, s_lon 1, c_lon 0):
    //!                        NED (n,e,u) -> ECEF (-e, u, n)
    //! Both match the physical local basis: at (0,0) North=+x, East=+y,
    //! Up=+z; at (0,90) North=+z, East=-x, Up=+y.

    use super::*;
    use crate::frames::realizations::Nad83_2011;
    use core::f64::consts::FRAC_PI_2;

    const LLH_EQ_PRIME: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
    const LLH_EQ_QUARTER: Vector3<f64> = Vector3::new(0.0, FRAC_PI_2, 0.0);

    #[test]
    fn station_without_velocity_does_not_move_between_epochs() {
        // A monument with no published tectonic velocity is stationary: the
        // `None` branch must return the position unchanged, not zero it.
        let pos = EcefPos::<Nad83_2011>::new(Vector3::new(1_000.0, -2_000.0, 3_000.0));
        let p = EpochPosition::<Nad83_2011, Arp>::new(pos, 2010.0, None);
        let later = p.at_epoch(2025.0);
        assert_eq!(*later.pos.coords(), Vector3::new(1_000.0, -2_000.0, 3_000.0));
        assert!((later.epoch_yr - 2025.0).abs() < 1e-15);
        assert_eq!(later.velocity_m_yr, None);
        // Same epoch must be a bit-exact identity too.
        assert_eq!(*p.at_epoch(2010.0).pos.coords(), Vector3::new(1_000.0, -2_000.0, 3_000.0));
    }

    #[test]
    fn station_with_velocity_extrapolates_linearly() {
        // p(2020) = p(2010) + v * 10 with v = (0.01, -0.02, 0.03) m/yr
        //          = (1000 + 0.1, -2000 - 0.2, 3000 + 0.3)
        let pos = EcefPos::<Nad83_2011>::new(Vector3::new(1_000.0, -2_000.0, 3_000.0));
        let v = Vector3::new(0.01, -0.02, 0.03);
        let p = EpochPosition::<Nad83_2011, Arp>::new(pos, 2010.0, Some(v));
        let got = *p.at_epoch(2020.0).pos.coords();
        let expect = Vector3::new(1_000.1, -2_000.2, 3_000.3);
        assert!((got - expect).norm() < 1e-12, "got {got:?}");
        // Going back in time is the exact negative increment.
        let back = *p.at_epoch(2000.0).pos.coords();
        assert!((back - Vector3::new(999.9, -1_999.8, 2_999.7)).norm() < 1e-12);
    }

    #[test]
    fn neu_to_ecef_at_the_equator_prime_meridian_is_a_permutation() {
        let got = neu_to_ecef(Vector3::new(11.0, 22.0, 33.0), LLH_EQ_PRIME);
        assert!((got - Vector3::new(33.0, 22.0, 11.0)).norm() < 1e-12, "got {got:?}");
    }

    #[test]
    fn neu_to_ecef_at_the_equator_quarter_meridian_is_a_permutation() {
        let got = neu_to_ecef(Vector3::new(11.0, 22.0, 33.0), LLH_EQ_QUARTER);
        assert!((got - Vector3::new(-22.0, 33.0, 11.0)).norm() < 1e-12, "got {got:?}");
    }

    #[test]
    fn arp_to_apc_offset_follows_the_local_ned_basis() {
        // PCO of (1, 2, 3) m in NED at the equator / prime meridian must add
        // exactly ECEF (3, 2, 1) m to the ARP.
        let arp = EpochPosition::<Nad83_2011, Arp>::new(
            EcefPos::new(Vector3::new(0.0, 0.0, 0.0)),
            2010.0,
            None,
        );
        let apc = arp.to_apc::<5>(Vector3::new(1_000.0, 2_000.0, 3_000.0), LLH_EQ_PRIME);
        assert!((apc.pos.into_vector() - Vector3::new(3.0, 2.0, 1.0)).norm() < 1e-9,
            "got {:?}", apc.pos.into_vector());
        // Converting back must undo it exactly.
        let back = apc.to_arp(Vector3::new(1_000.0, 2_000.0, 3_000.0), LLH_EQ_PRIME);
        assert!((back.pos.into_vector().norm()) < 1e-9, "got {:?}", back.pos.into_vector());
    }

    #[test]
    fn antenna_reference_names_are_distinct_per_marker() {
        use crate::frames::AntennaReference;
        assert_eq!(<Arp as AntennaReference>::NAME, "ARP");
        assert_eq!(<Apc<1> as AntennaReference>::NAME, "APC");
        assert_eq!(<GroundMonument as AntennaReference>::NAME, "Monument");
    }
}
