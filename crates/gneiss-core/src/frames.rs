//! Frame-tagged ECEF positions: compile-time prevention of datum mismatches.
//!
//! Every [`EcefPos`] carries its reference realization (`ITRF2014`, `IGS20`,
//! broadcast `WGS84`, ...) as a type parameter. Mixing frames without an
//! explicit [`EcefPos::convert_to`] call is a **compile error**, closing the
//! bug class behind the CAPO 39.7 mm vertical bias (LEIAR20 vs TRM phase
//! centers referenced to different realizations).
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
//! positive `rz` rotates `+X` toward `+Y`. Epochs are fractional years (cf.
//! `GpsTime::to_fractional_year`); parameters propagate linearly from their
//! reference epoch (IERS Conventions 2010 ch. 4, 14-parameter model). The
//! inverse inverts the full `(1+s)·(I + [ω]×)` matrix, so round trips stay
//! exact to floating-point rounding regardless of rotation magnitude.
//!
//! # Type safety: what the compiler rejects
//!
//! ```rust,ignore
//! use gneiss_core::frames::{EcefPos, Igs20, Itrf2014};
//!
//! fuse_base_rover(base, rover_igs20);
//! // ^^^ error[E0308]: expected `EcefPos<Itrf2014>`, found `EcefPos<Igs20>`
//! let tagged: EcefPos<Itrf2014> = rover_igs20;
//! // ^^^ error[E0308]: no implicit coercion between frame tags exists either
//! ```
//!
//! # Provenance of numeric constants
//!
//! * `ITRF2020 -> ITRF2014`: ITRF center table (Altamimi et al., epoch 2015.0):
//!   T = (-1.4, -0.9, +1.4) mm, D = -0.40 ppb, R = 0, rates 0. Verify against
//!   itrf.ign.fr before relying on sub-mm fidelity (`gneiss-geodesy/helmert.rs`
//!   carries a variant set (-1.4, -1.2, +1.2 mm plus rates); reconcile).
//! * `IGS20 -> ITRF2020`: aligned by construction (Rebischung et al. 2022);
//!   residuals < 1 mm neglected, so its link equals the ITRF2020 link.
//! * `WGS84(G2296)`: NGA STP aligns it to ITRF2020 at cm level; treated as
//!   identical to ITRF2020 here.

use core::marker::PhantomData;
use nalgebra::{Matrix3, Vector3};

/// Millimetres to metres.
const MM_TO_M: f64 = 1.0e-3;
/// Parts per billion to dimensionless.
const PPB: f64 = 1.0e-9;
/// Milliarcseconds to radians.
const MAS_TO_RAD: f64 = core::f64::consts::PI / 648_000_000.0;

/// 14-parameter Helmert set mapping one frame's coordinates into another.
///
/// Field units mirror the published ITRF/IGS tables; the `*_rate` fields are
/// annual rates used to propagate parameters from [`Self::ref_epoch_yr`] to
/// the observation epoch.
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

/// Spec-spelled alias for [`HelmertParams`], kept so code written against the
/// original API draft compiles unchanged.
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
            tx_rate: self.tx_rate,
            ty_rate: self.ty_rate,
            tz_rate: self.tz_rate,
            scale_rate: self.scale_rate,
            rx_rate: self.rx_rate,
            ry_rate: self.ry_rate,
            rz_rate: self.rz_rate,
        }
    }

    /// Forward transform `X_target = T + (1+s)(X + ω×X)`.
    #[must_use]
    pub fn apply(self, v: Vector3<f64>) -> Vector3<f64> {
        self.rotation_rad().cross(&v) * self.scale_factor()
            + v * self.scale_factor()
            + self.translation_m()
    }

    /// Inverse transform, computed exactly by inverting the full
    /// `(1+s)·(I + [ω]×)` matrix (the closed-form first-order shortcut would
    /// accumulate O(ω²·r) error, i.e. hundreds of µm for mas-scale rotations).
    /// Falls back to the first-order formula only if the matrix were
    /// numerically singular, which cannot happen for physical parameters
    /// (det ≈ (1+s)³(1+‖ω‖²) > 0).
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
        Vector3::new(self.rx_mas, self.ry_mas, self.rz_mas) * MAS_TO_RAD
    }

    fn scale_factor(self) -> f64 {
        1.0 + self.scale_ppb * PPB
    }
}

/// Composes links `src -> mid` then `mid -> dst` into one `src -> dst` set,
/// both evaluated at a common epoch. Translation and scale combine exactly;
/// rotations add component-wise (exact for mas-level angles). Rate fields are
/// taken from `mid_to_dst`; composing two rated parameter sets is not
/// supported (no such pair is published today).
///
/// Used to precompute the `TO_ITRF2014` constants of frames whose published
/// parameters target a different hub epoch/frame.
#[must_use]
pub fn chain(src_to_mid: HelmertParams, mid_to_dst: HelmertParams) -> HelmertParams {
    let combined_t_mm = mid_to_dst.apply(src_to_mid.translation_m()) * (1.0 / MM_TO_M);
    let mut out = mid_to_dst;
    out.tx_mm = combined_t_mm[0];
    out.ty_mm = combined_t_mm[1];
    out.tz_mm = combined_t_mm[2];
    out.scale_ppb =
        src_to_mid.scale_ppb + mid_to_dst.scale_ppb * src_to_mid.scale_factor();
    out.rx_mas += src_to_mid.rx_mas;
    out.ry_mas += src_to_mid.ry_mas;
    out.rz_mas += src_to_mid.rz_mas;
    out
}

/// Marker trait binding a coordinate realization to its published Helmert
/// link into the ITRF2014 hub.
pub trait ReferenceFrame {
    /// Human-readable realization name (logs, `Debug` output).
    const NAME: &'static str;

    /// Parameters transforming coordinates FROM this frame TO ITRF2014.
    /// `None` means "coincident with ITRF2014" and is reserved for ITRF2014
    /// itself; new frames must publish measured parameters instead.
    const HELMERT_TO_ITRF2014: Option<HelmertParams>;
}

/// Published ITRF2020 -> ITRF2014 link (epoch 2015.0; rotations and rates 0).
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
    /// Equals `chain(identity_at(2020.0), ITRF2020_TO_ITRF2014)`; the test
    /// `composed_constants_equal_runtime_chaining` pins that equivalence.
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
///
/// `Clone`/`Copy`/`PartialEq` are implemented manually so any marker type
/// implementing [`ReferenceFrame`] qualifies, even without deriving traits.
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

    /// Converts to frame `F2` via the ITRF2014 hub, propagating both Helmert
    /// sets to observation epoch `t_epoch_yr` (fractional years).
    pub fn convert_to<F2: ReferenceFrame>(&self, t_epoch_yr: f64) -> EcefPos<F2> {
        let to_hub = params_at(F::HELMERT_TO_ITRF2014, t_epoch_yr);
        let from_hub = params_at(F2::HELMERT_TO_ITRF2014, t_epoch_yr);
        EcefPos::new(from_hub.apply_inverse(to_hub.apply(self.0)))
    }
}

fn params_at(p: Option<HelmertParams>, t_yr: f64) -> HelmertParams {
    match p {
        Some(params) => params.at(t_yr),
        None => HelmertParams::identity_at(t_yr),
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
        /// +1000 mas about Z pins down the rotation sign convention.
        const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(HelmertParams {
            rz_mas: 1000.0,
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

    struct TestUnknown;
    impl ReferenceFrame for TestUnknown {
        const NAME: &'static str = "TEST-UNKNOWN";
        const HELMERT_TO_ITRF2014: Option<HelmertParams> = None;
    }

    fn full_params(ref_epoch_yr: f64) -> HelmertParams {
        HelmertParams {
            tx_mm: 1.0,
            ty_mm: -2.5,
            tz_mm: 3.0,
            scale_ppb: 12.0,
            rx_mas: 800.0,
            ry_mas: -1500.0,
            rz_mas: 2200.0,
            ref_epoch_yr,
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
        let p = EcefPos::<Itrf2014>::new(SITE);
        for t in [2015.0_f64, 2020.0, 2040.0] {
            assert_eq!(p.convert_to::<Itrf2014>(t).0, SITE);
        }
    }

    #[test]
    fn unknown_frame_without_params_is_identity() {
        let p = EcefPos::<TestUnknown>::new(SITE);
        assert_eq!(p.convert_to::<Itrf2014>(2030.0).0, SITE);
    }

    #[test]
    fn round_trip_preserves_position_below_spec_bound() {
        let cases = [
            (SITE, 2015.0_f64),
            (Vector3::new(-2_900_000.0, 1_300_000.0, 5_500_000.0), 2000.0),
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
        assert!(
            err < 1.0e-6,
            "{} -> {} round trip drifted {err:e} m",
            A::NAME,
            B::NAME
        );
    }

    #[test]
    fn itrf2020_to_itrf2014_shift_matches_published_parameters() {
        // Independent scalar derivation of X_14 = T + (1+s)·X_20 with
        // T = (-1.4, -0.9, +1.4) mm, s = -0.40 ppb, rotations zero.
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
        // Total effect is millimetre-level: big enough to corrupt PPK output.
        assert!(got.norm() > 1.0e-3 && got.norm() < 5.0e-3);
    }

    #[test]
    fn positive_rotation_moves_x_toward_y() {
        // Altamimi convention: X' = T + (1+s)(X + ω×X), so +rz sends +X toward +Y.
        let out = EcefPos::<TestRotZ>::new(Vector3::new(1.0, 0.0, 0.0))
            .convert_to::<Itrf2014>(2015.0)
            .0;
        let rz_rad = 1000.0 * MAS_TO_RAD;
        assert!((out[0] - 1.0).abs() < 1.0e-15);
        assert!((out[1] - rz_rad).abs() < 1.0e-18, "rotation sign flipped?");
        assert!(out[2].abs() < 1.0e-24);
    }

    #[test]
    fn epoch_propagation_is_linear_in_years_from_reference_epoch() {
        let shift_at = |t| {
            EcefPos::<TestRates>::new(SITE).convert_to::<Itrf2014>(t).0 - SITE
        };
        let at_ref = shift_at(2015.0);
        let plus10 = shift_at(2025.0);
        let minus10 = shift_at(2005.0);
        // Rate-only frame: z-shift = 1 mm/yr · dt, other axes untouched.
        // Tolerance is the ~1 nm cancellation floor of subtracting
        // coordinates ~5e6 m apart; a dropped or wrong-sign rate shifts z by
        // 10 mm and fails loudly.
        assert!(at_ref.abs().max() < 1.0e-15);
        assert!((plus10[2] - 0.010).abs() < 1.0e-8);
        assert!((minus10[2] + 0.010).abs() < 1.0e-8);
        assert!(plus10[0].abs() < 1.0e-8 && plus10[1].abs() < 1.0e-8);
        assert!(minus10[0].abs() < 1.0e-8 && minus10[1].abs() < 1.0e-8);
        // `at` must re-anchor the reference epoch and carry rates through.
        let propagated = full_params(2015.0).at(2025.0);
        assert_eq!(propagated.ref_epoch_yr, 2025.0);
        assert!((propagated.tx_rate - 0.2).abs() < 1.0e-15);
    }

    #[test]
    fn chaining_matches_sequential_application_and_shipped_constants() {
        // Property: the composed set transforms identically to applying the
        // two links in order (sequential application is ground truth).
        let h1 = HelmertParams {
            tx_mm: 2.0,
            ty_mm: -1.0,
            tz_mm: 4.0,
            scale_ppb: 5.0,
            rx_mas: 100.0,
            ry_mas: -60.0,
            rz_mas: 40.0,
            ref_epoch_yr: 2000.0,
            tx_rate: 0.1,
            ty_rate: 0.05,
            tz_rate: -0.02,
            scale_rate: 0.01,
            rx_rate: 0.3,
            ry_rate: 0.2,
            rz_rate: -0.1,
        };
        let v = Vector3::new(-6.378e6, 1.0e6, -2.0e6);
        let composed = chain(h1, ITRF2020_TO_ITRF2014);
        for t in [1995.0_f64, 2000.0, 2005.0] {
            let seq = ITRF2020_TO_ITRF2014.at(t).apply(h1.at(t).apply(v));
            assert!(
                (composed.at(t).apply(v) - seq).abs().max() < 1.0e-7,
                "chain diverges from sequential application at t={t}"
            );
        }
        // Shipped constants equal chaining an ITRF2020-aligned identity link.
        for frame_params in [
            Itrf2020::HELMERT_TO_ITRF2014,
            Igs20::HELMERT_TO_ITRF2014,
            Wgs84Broadcast::HELMERT_TO_ITRF2014,
        ] {
            let stored = frame_params.expect("shipped frames carry parameters");
            let chained =
                chain(HelmertParams::identity_at(stored.ref_epoch_yr), ITRF2020_TO_ITRF2014);
            for t in [stored.ref_epoch_yr, stored.ref_epoch_yr + 10.0] {
                let (p, q) = (chained.at(t), stored.at(t));
                assert!(
                    (p.apply(SITE) - q.apply(SITE)).abs().max() < 1.0e-9,
                    "chained {p:?} != stored {q:?}"
                );
            }
            assert_eq!(chained.ref_epoch_yr, stored.ref_epoch_yr);
        }
    }

    #[test]
    fn apply_inverse_inverts_apply_for_full_parameter_set() {
        for t in [1990.0_f64, 2015.0, 2035.0] {
            let p = full_params(2015.0).at(t);
            for v in [
                SITE,
                Vector3::new(-6.378e6, 0.0, 0.0),
                Vector3::new(1.0e5, -2.0e5, 3.0e5),
            ] {
                let back = p.apply_inverse(p.apply(v));
                // Tolerance = accumulated f64 rounding at Earth-radius
                // magnitudes (ulp ≈ 1 nm per op); far below the 0.1 mm spec.
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
        assert_eq!(Itrf2014::NAME, "ITRF2014");
        assert!(format!("{p:?}").contains("IGS20"));
        let copied = p;
        assert!(copied == p, "Copy/PartialEq must work for any tag type");
    }
}
