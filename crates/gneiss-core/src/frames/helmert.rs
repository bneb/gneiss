//! 14-parameter Helmert transformation between terrestrial reference frames.

use crate::constants::MILLIARCSEC_TO_RAD;
use nalgebra::{Matrix3, Vector3};

const MM_TO_M: f64 = 1.0e-3;
const PPB: f64 = 1.0e-9;

/// 14-parameter Helmert set mapping one frame's coordinates into another.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelmertParams {
    pub tx_mm: f64,
    pub ty_mm: f64,
    pub tz_mm: f64,
    pub scale_ppb: f64,
    pub rx_mas: f64,
    pub ry_mas: f64,
    pub rz_mas: f64,
    pub ref_epoch_yr: f64,
    pub tx_rate: f64,
    pub ty_rate: f64,
    pub tz_rate: f64,
    pub scale_rate: f64,
    pub rx_rate: f64,
    pub ry_rate: f64,
    pub rz_rate: f64,
}

impl HelmertParams {
    #[must_use]
    pub const fn identity_at(t_yr: f64) -> Self {
        Self {
            tx_mm: 0.0,
            ty_mm: 0.0,
            tz_mm: 0.0,
            scale_ppb: 0.0,
            rx_mas: 0.0,
            ry_mas: 0.0,
            rz_mas: 0.0,
            ref_epoch_yr: t_yr,
            tx_rate: 0.0,
            ty_rate: 0.0,
            tz_rate: 0.0,
            scale_rate: 0.0,
            rx_rate: 0.0,
            ry_rate: 0.0,
            rz_rate: 0.0,
        }
    }

    #[must_use]
    pub fn at(self, t_yr: f64) -> Self {
        let dt = t_yr - self.ref_epoch_yr;
        Self {
            tx_mm: self.tx_mm + self.tx_rate * dt,
            ty_mm: self.ty_mm + self.ty_rate * dt,
            tz_mm: self.tz_mm + self.tz_rate * dt,
            scale_ppb: self.scale_ppb + self.scale_rate * dt,
            rx_mas: self.rx_mas + self.rx_rate * dt,
            ry_mas: self.ry_mas + self.ry_rate * dt,
            rz_mas: self.rz_mas + self.rz_rate * dt,
            ref_epoch_yr: t_yr,
            ..self
        }
    }

    #[must_use]
    pub fn apply(self, v: Vector3<f64>) -> Vector3<f64> {
        self.rotation_rad().cross(&v) * self.scale_factor()
            + v * self.scale_factor()
            + self.translation_m()
    }

    #[must_use]
    pub fn apply_inverse(self, v: Vector3<f64>) -> Vector3<f64> {
        let w = self.rotation_rad();
        let s = self.scale_factor();
        let m = Matrix3::new(
            s, -s * w[2], s * w[1],
            s * w[2], s, -s * w[0],
            -s * w[1], s * w[0], s,
        );
        let rhs = v - self.translation_m();
        match m.try_inverse() {
            Some(inv) => inv * rhs,
            None => {
                let rel = rhs * (1.0 / s);
                rel + (-w).cross(&rel)
            }
        }
    }

    fn translation_m(self) -> Vector3<f64> {
        Vector3::new(self.tx_mm, self.ty_mm, self.tz_mm) * MM_TO_M
    }

    fn rotation_rad(self) -> Vector3<f64> {
        Vector3::new(self.rx_mas, self.ry_mas, self.rz_mas) * MILLIARCSEC_TO_RAD
    }

    fn scale_factor(self) -> f64 {
        1.0 + self.scale_ppb * PPB
    }
}
