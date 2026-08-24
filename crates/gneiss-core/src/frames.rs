//! Frame-tagged ECEF positions: compile-time prevention of datum mismatches.
//!
//! Every [`EcefPos`] carries its reference realization (`ITRF2014`, `IGS20`,
//! broadcast `WGS84`, ...) as a type parameter. Mixing frames without an
//! explicit [`EcefPos::convert_to`] call is a **compile error** — the fix for
//! the CAPO 39.7 mm vertical bias (phase centers referenced to different
//! realizations).
//!
//! # Transformation model
//!
//! Conversions route through an ITRF2014 hub: each frame stores the
//! 14-parameter Helmert set mapping *its own* coordinates **to ITRF2014**
//! ([`ReferenceFrame::HELMERT_TO_ITRF2014`]); `A -> B` composes as
//! `inv(H_B) ∘ H_A`. Parameters follow the Altamimi/ITRF position convention
//! (same as `gneiss-geodesy::helmert`):
//!
//! ```text
//! X_target = T + (1 + s) · (X_source + ω × X_source)
//! ```
//!
//! `T` translation (mm), `s` scale (ppb), `ω = (rx, ry, rz)` rotation (mas);
//! positive `rz` rotates `+X` toward `+Y`. Epochs are fractional years;
//! parameters propagate linearly from their reference epoch (IERS Conventions
//! 2010 ch. 4, 14-parameter model). The inverse inverts the full
//! `(1+s)·(I + [ω]×)` matrix exactly, so round trips stay exact to rounding.
//!
//! # Type safety: what the compiler rejects
//!
//! ```rust,ignore
//! use gneiss_core::frames::{EcefPos, Igs20, Itrf2014};
//! fuse_base_rover(base, rover_igs20);
//! // ^^^ error[E0308]: expected `EcefPos<Itrf2014>`, found `EcefPos<Igs20>`
//! let tagged: EcefPos<Itrf2014> = rover_igs20;
//! // ^^^ error[E0308]: no implicit coercion between frame tags either
//! ```
//!
//! # Provenance of numeric constants
//!
//! * `ITRF2020 -> ITRF2014`: ITRF center table (Altamimi et al., epoch 2015.0):
//!   T = (-1.4, -0.9, +1.4) mm, D = -0.40 ppb, rates pending verification
//!   (`gneiss-geodesy/helmert.rs` carries a conflicting variant; reconcile).
//! * `IGS20`, `WGS84(G2296)`: aligned to ITRF2020 by construction / NGA STP
//!   (< 1 mm and cm-level respectively); both links equal the ITRF2020 row.

use core::marker::PhantomData;
use crate::constants::MILLIARCSEC_TO_RAD;
use nalgebra::{Matrix3, Vector3};

/// Millimetres to metres.
const MM_TO_M: f64 = 1.0e-3;
/// Parts per billion to dimensionless.
const PPB: f64 = 1.0e-9;

/// 14-parameter Helmert set mapping one frame's coordinates into another;
/// units mirror the published ITRF/IGS tables and `*_rate` fields are annual
/// rates for propagation from [`Self::ref_epoch_yr`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelmertParams {
    /// Translation X at `ref_epoch_yr` (mm).
    pub tx_mm: f64,
    /// Translation Y at `ref_epoch_yr` (mm).
    pub ty_mm: f64,
    /// Translation Z at `ref_epoch_yr` (mm).
    pub tz_mm: f64,
    /// Scale at `ref_epoch_yr` (ppb).
    pub scale_ppb: f64,
    /// Rotation about X at `ref_epoch_yr` (mas).
    pub rx_mas: f64,
    /// Rotation about Y at `ref_epoch_yr` (mas).
    pub ry_mas: f64,
    /// Rotation about Z at `ref_epoch_yr` (mas).
    pub rz_mas: f64,
    /// Reference epoch of the published parameters (fractional years).
    pub ref_epoch_yr: f64,
    /// Translation X rate (mm/yr).
    pub tx_rate: f64,
    /// Translation Y rate (mm/yr).
    pub ty_rate: f64,
    /// Translation Z rate (mm/yr).
    pub tz_rate: f64,
    /// Scale rate (ppb/yr).
    pub scale_rate: f64,
    /// Rotation X rate (mas/yr).
    pub rx_rate: f64,
    /// Rotation Y rate (mas/yr).
    pub ry_rate: f64,
    /// Rotation Z rate (mas/yr).
    pub rz_rate: f64,
}

/// Spec-spelled alias for [`HelmertParams`] (draft-API compatibility).
pub type HelbertParams = HelmertParams;

impl HelmertParams {
    /// All-zero parameter set (exact identity transform) anchored at `t_yr`.
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

    /// Propagates the parameters linearly to observation epoch `t_yr`.
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

    /// Forward transform `X_target = T + (1+s)(X + ω×X)`.
    #[must_use]
    pub fn apply(self, v: Vector3<f64>) -> Vector3<f64> {
        self.rotation_rad().cross(&v) * self.scale_factor()
            + v * self.scale_factor()
            + self.translation_m()
    }

    /// Inverse transform; inverts the full `(1+s)·(I + [ω]×)` matrix exactly
    /// (first-order shortcuts accumulate O(ω²·r) error). Falls back to a
    /// first-order formula only if numerically singular — impossible for
    /// physical parameters (det ≈ (1+s)³(1+‖ω‖²) > 0).
    #[must_use]
    pub fn apply_inverse(self, v: Vector3<f64>) -> Vector3<f64> {
        let w = self.rotation_rad();
        let s = self.scale_factor();
        let m = Matrix3::new(
            s, -s * w[2], s * w[1], //
            s * w[2], s, -s * w[0], //
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

/// Marker trait binding a coordinate realization to its published Helmert
/// link into the ITRF2014 hub.
pub trait ReferenceFrame {
    /// Human-readable realization name (logs, `Debug` output).
    const NAME: &'static str;

    /// Parameters transforming coordinates FROM this frame TO ITRF2014.
    /// `None` means "coincident with ITRF2014" (reserved for ITRF2014 itself);
    /// new frames must publish measured parameters, or an offline-precomposed
    /// constant if their published parameters target a different hub frame.
    const HELMERT_TO_ITRF2014: Option<HelmertParams>;
}

/// Published ITRF2020 -> ITRF2014 link (epoch 2015.0; rotations zero).
/// WARNING: rates zeroed pending confirmation at itrf.ign.fr (some sources
/// report Tz drift ≈ −0.2 mm/yr here ⇒ ~2 mm Z error by 2025 if true).
const ITRF2020_TO_ITRF2014: HelmertParams = HelmertParams {
    tx_mm: -1.4,
    ty_mm: -0.9,
    tz_mm: 1.4,
    scale_ppb: -0.40,
    rx_mas: 0.0,
    ry_mas: 0.0,
    rz_mas: 0.0,
    ref_epoch_yr: 2015.0,
    tx_rate: 0.0,
    ty_rate: 0.0,
    tz_rate: 0.0,
    scale_rate: 0.0,
    rx_rate: 0.0,
    ry_rate: 0.0,
    rz_rate: 0.0,
};

/// ITRF2014 itself: the hub frame, exact identity (`None` above).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Itrf2014;

impl ReferenceFrame for Itrf2014 {
    const NAME: &'static str = "ITRF2014";
    const HELMERT_TO_ITRF2014: Option<HelmertParams> = None;
}

/// ITRF2020 (Altamimi et al., 2023).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Itrf2020;

impl ReferenceFrame for Itrf2020 {
    const NAME: &'static str = "ITRF2020";
    const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(ITRF2020_TO_ITRF2014);
}

/// IGS cumulative frame IGS20 (Rebischung et al., 2022), ITRF2020-aligned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Igs20;

impl ReferenceFrame for Igs20 {
    const NAME: &'static str = "IGS20";
    const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(ITRF2020_TO_ITRF2014);
}

/// Broadcast WGS84 realization G2296, treated as ITRF2020-aligned (cm level).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wgs84Broadcast;

impl ReferenceFrame for Wgs84Broadcast {
    const NAME: &'static str = "WGS84(Broadcast)";
    const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(ITRF2020_TO_ITRF2014);
}

/// ECEF position vector tagged with its reference frame at the type level.
/// Manual trait impls let any [`ReferenceFrame`] marker qualify.
pub struct EcefPos<F: ReferenceFrame>(pub Vector3<f64>, pub PhantomData<F>);

impl<F: ReferenceFrame> EcefPos<F> {
    /// Tags `v` as being expressed in frame `F`.
    #[must_use]
    pub fn new(v: Vector3<f64>) -> Self {
        Self(v, PhantomData)
    }

    /// Magnitude of the position vector (metres).
    #[must_use]
    pub fn norm(&self) -> f64 {
        self.0.norm()
    }

    /// Read-only access to the underlying coordinates.
    #[must_use]
    pub const fn vector(&self) -> &Vector3<f64> {
        &self.0
    }

    /// Converts to frame `F2` via the ITRF2014 hub; both Helmert sets are
    /// propagated to observation epoch `t_epoch_yr` (fractional years).
    pub fn convert_to<F2: ReferenceFrame>(&self, t_epoch_yr: f64) -> EcefPos<F2> {
        let to_hub = params_at(F::HELMERT_TO_ITRF2014, t_epoch_yr);
        let from_hub = params_at(F2::HELMERT_TO_ITRF2014, t_epoch_yr);
        EcefPos::new(from_hub.apply_inverse(to_hub.apply(self.0)))
    }
}

fn params_at(p: Option<HelmertParams>, t_yr: f64) -> HelmertParams {
    p.map_or_else(|| HelmertParams::identity_at(t_yr), |params| params.at(t_yr))
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
        // Marker types carry no data; equal tags are guaranteed by the types.
        self.0 == other.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::format;

    const SITE: Vector3<f64> = Vector3::new(4_027_893.0, 307_041.0, 4_919_475.0);

    struct TestRotZ;
    impl ReferenceFrame for TestRotZ {
        const NAME: &'static str = "TEST-ROTZ";
        /// Distinct rotations about all three axes pin the sign convention.
        const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(HelmertParams {
            rx_mas: 1000.0,
            ry_mas: -2000.0,
            rz_mas: 3000.0,
            ..HelmertParams::identity_at(2015.0)
        });
    }

    struct TestRates;
    impl ReferenceFrame for TestRates {
        const NAME: &'static str = "TEST-RATES";
        const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(HelmertParams {
            tz_rate: 1.0,
            ..HelmertParams::identity_at(2015.0)
        });
    }

    fn full_params() -> HelmertParams {
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
    fn identity_conversion_is_bit_exact() {
        // Hub frame (None params) converts exactly at any epoch; ITRF2014
        // must stay the unique None frame (None => identity semantics).
        let p = EcefPos::<Itrf2014>::new(SITE);
        assert_eq!(Itrf2014::HELMERT_TO_ITRF2014, None);
        for t in [2015.0_f64, 2020.0, 2040.0] {
            assert_eq!(p.convert_to::<Itrf2014>(t).0, SITE);
        }
    }

    #[test]
    fn round_trip_preserves_position_below_spec_bound() {
        let cases = [
            (SITE, 2015.0_f64),
            (Vector3::new(500_000.0, -5_800_000.0, 2_100_000.0), 2040.0),
        ];
        for (v, t) in cases {
            assert_round_trip::<Igs20, Wgs84Broadcast>(v, t);
            assert_round_trip::<Itrf2020, Igs20>(v, t);
            assert_round_trip::<Wgs84Broadcast, Itrf2014>(v, t);
            assert_round_trip::<Itrf2020, Itrf2014>(v, t);
        }
    }

    fn assert_round_trip<A: ReferenceFrame, B: ReferenceFrame>(v: Vector3<f64>, t: f64) {
        let orig = EcefPos::<A>::new(v);
        let back = orig.convert_to::<B>(t).convert_to::<A>(t);
        let err = (back.vector() - orig.vector()).norm();
        assert!(err < 1.0e-6, "{} -> {} drifted {err:e} m", A::NAME, B::NAME);
    }

    #[test]
    fn itrf2020_to_itrf2014_shift_matches_published_parameters() {
        // Independent scalar derivation of X_14 = T + (1+s)·X_20, T in mm.
        let expected_shift = Vector3::new(
            -1.4e-3 + -0.40e-9 * SITE[0],
            -0.9e-3 + -0.40e-9 * SITE[1],
            1.4e-3 + -0.40e-9 * SITE[2],
        );
        let got = EcefPos::<Itrf2020>::new(SITE)
            .convert_to::<Itrf2014>(2015.0)
            .0
            - SITE;
        assert!(
            (got - expected_shift).abs().max() < 1.0e-8,
            "shift {got:?} != published {expected_shift:?}"
        );
        // IGS20 / broadcast WGS84 are ITRF2020-aligned: same published row.
        assert_eq!(Igs20::HELMERT_TO_ITRF2014, Some(ITRF2020_TO_ITRF2014));
        assert_eq!(
            Wgs84Broadcast::HELMERT_TO_ITRF2014,
            Some(ITRF2020_TO_ITRF2014)
        );
    }

    #[test]
    fn rotation_axes_follow_altamimi_sign_convention() {
        // X' = X + ω×X with ω = (1000, −2000, 3000) mas. Expected values are
        // hardcoded (rad per mas literal) so a corrupted unit conversion or
        // flipped sign on ANY axis cannot self-verify.
        const K: f64 = 4.848_136_811_095_36e-9; // rad per mas
        let cases = [
            (
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(1.0, 3000.0 * K, 2000.0 * K),
            ),
            (
                Vector3::new(0.0, 1.0, 0.0),
                Vector3::new(-3000.0 * K, 1.0, 1000.0 * K),
            ),
            (
                Vector3::new(0.0, 0.0, 1.0),
                Vector3::new(-2000.0 * K, -1000.0 * K, 1.0),
            ),
        ];
        for (v, expected) in cases {
            let out = EcefPos::<TestRotZ>::new(v).convert_to::<Itrf2014>(2015.0).0;
            assert!(
                (out - expected).abs().max() < 1.0e-12,
                "rotation {v:?} -> {out:?}, expected {expected:?}"
            );
        }
    }

    #[test]
    fn epoch_propagation_is_linear_in_years_from_reference_epoch() {
        let shift_at = |t| EcefPos::<TestRates>::new(SITE).convert_to::<Itrf2014>(t).0 - SITE;
        let at_ref = shift_at(2015.0);
        let plus10 = shift_at(2025.0);
        let minus10 = shift_at(2005.0);
        // Rate-only frame: z-shift = 1 mm/yr · dt, other axes untouched. The
        // 1e-8 tolerance is the nm-scale cancellation floor at |r| ≈ 5e6 m.
        assert!(at_ref.abs().max() < 1.0e-15);
        assert!((plus10[2] - 0.010).abs() < 1.0e-8 && (minus10[2] + 0.010).abs() < 1.0e-8);
        assert!(
            plus10[0].abs() + plus10[1].abs() + minus10[0].abs() + minus10[1].abs()
                < 1.0e-7
        );
        // `at` must re-anchor the reference epoch and propagate ALL seven
        // quantities with their signed rates (literals kill sign mutants).
        let q = full_params().at(2025.0);
        assert_eq!(q.ref_epoch_yr, 2025.0);
        let drift = (q.tx_mm - 3.0).abs()
            + (q.ty_mm + 3.5).abs()
            + (q.tz_mm - 6.0).abs()
            + (q.scale_ppb - 12.5).abs()
            + (q.rx_mas - 804.0).abs()
            + (q.ry_mas + 1506.0).abs()
            + (q.rz_mas - 2207.0).abs();
        assert!(drift < 1.0e-9, "propagation drift {drift:e}");
    }

    #[test]
    fn apply_inverse_inverts_apply_for_full_parameter_set() {
        for t in [1990.0_f64, 2015.0, 2035.0] {
            let p = full_params().at(t);
            for v in [
                SITE,
                Vector3::new(-6.378e6, 0.0, 0.0),
                Vector3::new(1.0e5, -2.0e5, 3.0e5),
            ] {
                let back = p.apply_inverse(p.apply(v));
                // f64 noise floor at |r| ≈ 5e6 m is ~1 nm/op; spec bound 0.1 mm.
                assert!(
                    (back - v).abs().max() < 1.0e-8,
                    "inverse failed at t={t}: {back:?} vs {v:?}"
                );
            }
        }
    }

    #[test]
    fn accessors_and_frame_metadata_work() {
        let p = EcefPos::<Igs20>::new(Vector3::new(3.0, 4.0, 0.0));
        assert_eq!(p.norm(), 5.0);
        assert_eq!(p.vector(), &Vector3::new(3.0, 4.0, 0.0));
        assert_eq!(<Igs20 as ReferenceFrame>::NAME, "IGS20");
        assert!(format!("{p:?}").contains("IGS20"), "Debug must show frame");
        let copied = p;
        let different = EcefPos::<Igs20>::new(Vector3::new(3.0, 4.0, 1.0));
        assert!(copied == p && different != p, "Copy/PartialEq broken");
    }
}
