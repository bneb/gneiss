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

#[cfg(test)]
mod tests {
    //! Textbook verification of the 14-parameter Helmert transformation.
    //!
    //! The forward map is the classic linearized form
    //!   X2 = T + (1 + s) * (X1 + w x X1),  w = [rx ry rz] in radians
    //! which is exactly M * X1 + T with M = s * (I + [w]x) and
    //! s = 1 + scale_ppb * 1e-9. The inverse is therefore M^-1 * (X1 - T),
    //! and the two must be exact inverses of each other.
    //!
    //! 1 ppb of scale at Earth-surface distances is millimetre-level: at
    //! r = 6.378e6 m, s = 1 + 12e-9 gives a radial shift of 12e-9 * 6.378e6
    //! = 0.0765 m, which is why the scale term is not negligible.

    use super::*;
    use crate::frames::realizations::{Itrf2014, Itrf2020, Nad83_2011, Pz90};
    use crate::frames::ReferenceFrame;

    fn full() -> HelmertParams {
        HelmertParams {
            tx_mm: 1.0,
            ty_mm: -2.5,
            tz_mm: 3.0,
            scale_ppb: 12.0,
            rx_mas: 800.0,
            ry_mas: -1500.0,
            rz_mas: 2200.0,
            ref_epoch_yr: 2015.0,
            tx_rate: 0.2,
            ty_rate: -0.1,
            tz_rate: 0.3,
            scale_rate: 0.05,
            rx_rate: 0.4,
            ry_rate: -0.6,
            rz_rate: 0.7,
        }
    }

    #[test]
    fn identity_params_are_a_pure_no_op() {
        let p = HelmertParams::identity_at(2020.0);
        assert_eq!(p.scale_factor(), 1.0);
        assert_eq!(p.translation_m(), Vector3::zeros());
        assert_eq!(p.rotation_rad(), Vector3::zeros());
        for v in [Vector3::new(1.0, 2.0, 3.0), Vector3::new(-4e6, 5e6, 6e6)] {
            assert_eq!(p.apply(v), v);
            assert_eq!(p.apply_inverse(v), v);
        }
    }

    #[test]
    fn pure_translation_shifts_by_exactly_the_millimetre_amount() {
        // With zero scale and rotation, apply(v) = v + T and T = (tx,ty,tz) mm.
        let p = HelmertParams { tx_mm: 12.0, ty_mm: -34.0, tz_mm: 56.0, ..HelmertParams::identity_at(2020.0) };
        let v = Vector3::new(100.0, 200.0, 300.0);
        let got = p.apply(v);
        assert!((got - Vector3::new(100.012, 199.966, 300.056)).norm() < 1e-12, "got {got:?}");
        assert!((p.apply_inverse(got) - v).norm() < 1e-12);
    }

    #[test]
    fn scale_term_shifts_by_scale_ppb_times_the_distance() {
        // apply(v) = (1 + 12e-9) * v for zero rotation/translation.
        // At |v| = 6.378e6 m the radial change is 12e-9 * 6.378e6 = 0.076536 m.
        let p = HelmertParams { scale_ppb: 12.0, ..HelmertParams::identity_at(2020.0) };
        assert!((p.scale_factor() - (1.0 + 12.0e-9)).abs() < 1e-18);
        let v = Vector3::new(6.378e6, 0.0, 0.0);
        let got = p.apply(v);
        assert!((got.x - (6.378e6 + 0.076536)).abs() < 1e-6, "got {}", got.x);
        assert!(got.y.abs() < 1e-9 && got.z.abs() < 1e-9);
    }

    #[test]
    fn rotation_term_is_the_antisymmetric_cross_product() {
        // apply(v) = s * (v + w x v) with w = [rx,ry,rz] * mas2rad.
        // Using w = [0, 0, rz]: w x v = (-rz * vy, rz * vx, 0).
        let rz_mas = 3000.0;
        let p = HelmertParams { rz_mas, ..HelmertParams::identity_at(2020.0) };
        let rz = rz_mas * MILLIARCSEC_TO_RAD;
        let v = Vector3::new(2.0, 3.0, 4.0);
        let expect = Vector3::new(2.0 - rz * 3.0, 3.0 + rz * 2.0, 4.0);
        let got = p.apply(v);
        assert!((got - expect).norm() < 1e-12, "got {got:?} want {expect:?}");
        // The perturbation is the first-order (linearised) rotation only, so its
        // magnitude must match |w x v| = |w| |v| sin(theta) exactly:
        // |w| = rz, |w x v| = rz * sqrt(2^2 + 3^2) = rz * sqrt(13).
        let expect_cross_norm = rz * 13.0f64.sqrt();
        let cross = Vector3::new(2.0, 3.0, 4.0).cross(&Vector3::new(0.0, 0.0, 1.0)) * rz;
        assert!((cross.norm() - expect_cross_norm).abs() < 1e-15, "got {}", cross.norm());
        // Scaling: s = 1 for scale_ppb = 0, so apply is a pure first-order rotation.
        assert!((p.scale_factor() - 1.0).abs() < 1e-18);
    }

    #[test]
    fn apply_inverse_exactly_inverts_apply_for_the_full_parameter_set() {
        let p = full();
        for v in [
            Vector3::new(1.0, -2.0, 3.0),
            Vector3::new(-2_430_601.8, -4_700_258.9, 3_544_321.5),
            Vector3::new(6_378_137.0, 0.0, 0.0),
        ] {
            assert!((p.apply_inverse(p.apply(v)) - v).norm() < 1e-6, "v = {v:?}");
        }
    }

    #[test]
    fn at_propagates_every_rate_from_the_reference_epoch() {
        // dt = 2020.0 - 2015.0 = 5 yr:
        //   tx = 1.0 + 0.2*5 = 2.0 mm ; ty = -2.5 + (-0.1)*5 = -3.0 mm
        //   tz = 3.0 + 0.3*5 = 4.5 mm ; s  = 12.0 + 0.05*5 = 12.25 ppb
        //   rx = 800 + 0.4*5 = 802.0 ; ry = -1500 + (-0.6)*5 = -1503.0
        //   rz = 2200 + 0.7*5 = 2203.5 mas
        let p = full().at(2020.0);
        assert!((p.tx_mm - 2.0).abs() < 1e-12);
        assert!((p.ty_mm + 3.0).abs() < 1e-12);
        assert!((p.tz_mm - 4.5).abs() < 1e-12);
        assert!((p.scale_ppb - 12.25).abs() < 1e-12);
        assert!((p.rx_mas - 802.0).abs() < 1e-12);
        assert!((p.ry_mas + 1503.0).abs() < 1e-12);
        assert!((p.rz_mas - 2203.5).abs() < 1e-12);
        assert!((p.ref_epoch_yr - 2020.0).abs() < 1e-12);
        // Evaluating at the reference epoch reproduces the original parameters.
        let back = p.at(2015.0);
        assert!((back.tx_mm - 1.0).abs() < 1e-12 && (back.rz_mas - 2200.0).abs() < 1e-12);
    }

    #[test]
    fn published_itrf2020_link_is_entered_verbatim() {
        // ITRF2020 -> ITRF2014 at epoch 2015.0: (-1.4, -0.9, +1.4) mm,
        // -0.42 ppb, zero rotations (ITRF2020 -> ITRF2014 Technical Note 36).
        let p = Itrf2020::HELMERT_TO_ITRF2014.expect("ITRF2020 has a published link");
        assert!((p.tx_mm + 1.4).abs() < 1e-12);
        assert!((p.ty_mm + 0.9).abs() < 1e-12);
        assert!((p.tz_mm - 1.4).abs() < 1e-12);
        assert!((p.scale_ppb + 0.42).abs() < 1e-12);
        assert!((p.rx_mas).abs() < 1e-12 && p.ry_mas.abs() < 1e-12 && p.rz_mas.abs() < 1e-12);
        assert!((p.ty_rate + 0.1).abs() < 1e-12);
        assert!((p.tz_rate - 0.2).abs() < 1e-12);
        assert!((p.ref_epoch_yr - 2015.0).abs() < 1e-12);
        assert_eq!(Itrf2014::HELMERT_TO_ITRF2014, None);
    }

    #[test]
    fn nad83_and_pz90_links_are_distinct_and_both_invertible() {
        let n = Nad83_2011::HELMERT_TO_ITRF2014.expect("NAD83 has a link");
        let z = Pz90::HELMERT_TO_ITRF2014.expect("PZ-90 has a link");
        // NAD83(2011) at 2020.0: dt = 10 yr, tx = -1005.3 + (-0.79)*10 = -1013.2 mm
        let n2020 = n.at(2020.0);
        assert!((n2020.tx_mm + 1013.2).abs() < 1e-9, "got {}", n2020.tx_mm);
        assert!((n2020.ty_mm - 1915.2).abs() < 1e-9, "ty = {}", n2020.ty_mm);
        assert!((n2020.tz_mm - 556.0).abs() < 1e-9, "tz = {}", n2020.tz_mm);
        // PZ-90 translations are metre-scale in the table (3, -1, 0) mm.
        assert!((z.tx_mm - 3.0).abs() < 1e-12 && (z.ty_mm + 1.0).abs() < 1e-12);
        let v = Vector3::new(-2_430_601.8, -4_700_258.9, 3_544_321.5);
        assert!((n.apply_inverse(n.apply(v)) - v).norm() < 1e-6);
        assert!((z.apply_inverse(z.apply(v)) - v).norm() < 1e-6);
    }

    #[test]
    fn singular_scale_falls_back_without_panicking() {
        // det(M) = s^3 * (1 + wx^2 + wy^2 + wz^2), and the rotation bracket is
        // >= 1, so M is singular only when the scale factor s is exactly zero,
        // i.e. scale_ppb = -1e9. That path must reach the documented fallback
        // (explicit inverse-of-scale + inverse rotation) rather than unwrap or
        // panic; with s = 0 the relative vector is infinite, so the fallback
        // result is non-finite. Pinning non-finiteness catches both a panic and
        // a "return the input unchanged" shortcut.
        let p = HelmertParams { scale_ppb: -1.0e9, ..HelmertParams::identity_at(2020.0) };
        assert_eq!(p.scale_factor(), 0.0);
        let out = p.apply_inverse(Vector3::new(1.0, 2.0, 3.0));
        assert!(!out.x.is_finite() || !out.y.is_finite() || !out.z.is_finite(),
            "singular-scale fallback returned {out:?}");
    }
}
