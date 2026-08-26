//! Processing-dynamics profile: static-monument vs kinematic rover motion.
//!
//! The DD-IEKF pipeline historically hard-coded static-monument
//! assumptions: tight acceleration process noise (`q_accel = 1e-6`),
//! a two-phase "monument lock", and absolute combiner disagreement
//! limits justified by a receiver bolted to bedrock. Those assumptions
//! are wrong for a moving rover and this module is the single source
//! of truth for both profiles:
//!
//! - [`ProcessingDynamics::Static`] reproduces the legacy behaviour
//!   bit-for-bit (default; selected whenever `GNEISS_DYNAMICS` is unset
//!   or not exactly `kinematic`).
//! - [`ProcessingDynamics::Kinematic`] re-randomizes position/velocity
//!   every epoch (`q_accel ≈ 1`, the value the legacy doc comments call
//!   "~55 m of position re-randomization per 30 s epoch"), suppresses
//!   the monument lock, widens the innovation gate before robust
//!   weighting engages, and scales combiner disagreement limits by the
//!   per-epoch formal sigma instead of fixed metres.
//!
//! Everything here is pure: env access happens only in
//! [`ProcessingDynamics::from_env`], so every decision is unit-testable
//! and no test mutates process state.

/// Environment variable selecting the kinematic profile (`kinematic`,
/// exact match after trimming). Any other value — including absence —
/// keeps the legacy static path byte-identical.
pub const ENV_VAR: &str = "GNEISS_DYNAMICS";

/// Legacy static-monument acceleration PSD (m/s^2): the value the eval
/// binaries have always passed for CORS monuments.
pub const STATIC_Q_ACCEL: f64 = 1e-6;

/// Kinematic acceleration PSD (m/s^2). At 30 s epochs the position
/// process noise scales as dt^3/3 * q ≈ 55 m of re-randomization, which
/// is the correct prior for a moving rover and catastrophic for a
/// monument — hence the profile split. Also >= 100x the static value,
/// keeping velocity states live (Q_vel = dt * q).
pub const KINEMATIC_Q_ACCEL: f64 = 1.0;

/// Legacy robust-innovation gate multiplier (no change).
pub const STATIC_INNOV_GATE_SCALE: f64 = 1.0;

/// Kinematic gate multiplier applied to
/// [`ROBUST_INNOVATION_THRESHOLD`] (nominal 9 ≈ 3-sigma per axis).
/// During coordinated acceleration the constant-velocity model
/// mismatch biases every DD innovation coherently; without the wider
/// gate the robust weighting deweights genuine measurements exactly
/// when dynamics are largest. 3x moves the knee to ~5-sigma while
/// still bounding gross outliers (slips enter orders of magnitude
/// above it and are handled by the cycle-slip machinery).
pub const KINEMATIC_INNOV_GATE_SCALE: f64 = 3.0;

/// Rover motion model assumed by the processing pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ProcessingDynamics {
    /// Static monument: legacy behaviour, byte-identical defaults.
    #[default]
    Static,
    /// Kinematic rover: mobile Q, no monument lock, scaled gates.
    Kinematic,
}

impl ProcessingDynamics {
    /// Pure mapping of an env value to the profile. Only the exact
    /// string `kinematic` opts in; everything else stays Static.
    pub fn from_env_value(value: Option<&str>) -> Self {
        match value.map(str::trim) {
            Some("kinematic") => Self::Kinematic,
            _ => Self::Static,
        }
    }

    /// Read the profile from the process environment (once per call;
    /// binaries should call this in `main` and pass the result down).
    pub fn from_env() -> Self {
        Self::from_env_value(std::env::var(ENV_VAR).ok().as_deref())
    }

    /// True when the rover is assumed mobile.
    pub fn is_kinematic(self) -> bool {
        matches!(self, Self::Kinematic)
    }

    /// Default acceleration process noise (m/s^2) for the profile.
    /// An explicit `q_accel` option overrides this via [`resolve_q_accel`].
    pub fn q_accel_default(self) -> f64 {
        match self {
            Self::Static => STATIC_Q_ACCEL,
            Self::Kinematic => KINEMATIC_Q_ACCEL,
        }
    }

    /// Multiplier on the robust-innovation threshold before Huber-style
    /// variance inflation engages.
    pub fn innovation_gate_scale(self) -> f64 {
        match self {
            Self::Static => STATIC_INNOV_GATE_SCALE,
            Self::Kinematic => KINEMATIC_INNOV_GATE_SCALE,
        }
    }
}

/// Pipeline fallback when no explicit `q_accel` option is set. The
/// legacy passes used `q_accel.unwrap_or(1.0)`; preserving that value
/// for BOTH profiles keeps every unset-option caller byte-identical.
/// Profiles select their Q explicitly at the options level
/// ([`STATIC_Q_ACCEL`] / [`KINEMATIC_Q_ACCEL`]).
pub const Q_ACCEL_UNSET_FALLBACK: f64 = 1.0;

/// Honesty-gate steepness for the kinematic combiner: forward/backward
/// passes estimate the SAME trajectory, so disagreement beyond k formal
/// sigmas means one of them is wrong and fixed claims must be dropped.
pub const KIN_DISAGREE_K_SIGMA: f64 = 6.0;

/// Fusion-window steepness when only one pass claims a fix: fuse while
/// the passes agree within this many combined sigmas, else keep the
/// fixed side alone.
pub const KIN_FUSE_CROSS_K_SIGMA: f64 = 6.0;

/// Fusion-window steepness when both passes claim a fix (tighter than
/// the honesty gate: independent integer sets that differ by even ~2
/// sigmas should not be averaged).
pub const KIN_FUSE_BOTH_FIXED_K_SIGMA: f64 = 2.0;

/// Floor (m) for the sigma-scaled honesty and one-sided-fusion limits,
/// set EQUAL to the audited static bound (0.50 m, justified for
/// monuments and measured against CORS day sets where honest fwd/bwd
/// separations reach p99 ~ 0.41 m). Measured on the multi2025 replay:
/// formal sigmas are >10x optimistic versus realized pass-to-pass
/// disagreement, so a lower floor rejected ~60% of honest epochs. The
/// kinematic rule therefore only ever WIDENS the audited tolerance;
/// tightening below it awaits covariance recalibration on true moving
/// data.
pub const KIN_THRESHOLD_FLOOR_M: f64 = 0.50;

/// Cap (m) for the sigma-scaled honesty and one-sided-fusion limits: a
/// diverged covariance must not excuse arbitrary disagreement.
pub const KIN_THRESHOLD_CAP_M: f64 = 10.0;

/// Floor (m) for the both-sides-fixed fusion window: averaging two
/// integer sets farther apart than the audited 0.20 m bound
/// manufactures a position that neither pass claims.
pub const KIN_BOTH_FIXED_FLOOR_M: f64 = 0.20;
/// Cap (m) for the both-sides-fixed fusion window.
pub const KIN_BOTH_FIXED_CAP_M: f64 = 2.0;

/// Sigma-scaled separation limit (m):
/// `clamp(k * max(sigma_a, sigma_b), floor_m, cap_m)`.
pub fn kinematic_sep_limit_m(
    k_sigma: f64,
    sigma_a: f64,
    sigma_b: f64,
    floor_m: f64,
    cap_m: f64,
) -> f64 {
    (k_sigma * sigma_a.max(sigma_b)).clamp(floor_m, cap_m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::estimators::rtk_iekf::update::ROBUST_INNOVATION_THRESHOLD;

    #[test]
    fn env_unset_and_garbage_select_static() {
        assert_eq!(ProcessingDynamics::from_env_value(None), ProcessingDynamics::Static);
        assert_eq!(ProcessingDynamics::from_env_value(Some("")), ProcessingDynamics::Static);
        assert_eq!(ProcessingDynamics::from_env_value(Some("static")), ProcessingDynamics::Static);
        assert_eq!(
            ProcessingDynamics::from_env_value(Some("KINEMATIC")),
            ProcessingDynamics::Static,
            "only exact lowercase opt-in keeps the default conservative"
        );
        assert_eq!(
            ProcessingDynamics::from_env_value(Some(" kinematic ")),
            ProcessingDynamics::Kinematic,
            "surrounding whitespace must not defeat the explicit opt-in"
        );
    }

    #[test]
    fn default_trait_is_static() {
        assert_eq!(ProcessingDynamics::default(), ProcessingDynamics::Static);
    }

    #[test]
    fn static_profile_reproduces_legacy_constants_bitwise() {
        // The unset-option pipeline fallback must stay the legacy 1.0.
        assert_eq!(Q_ACCEL_UNSET_FALLBACK.to_bits(), 1.0_f64.to_bits());
        assert_eq!(STATIC_Q_ACCEL.to_bits(), 1e-6_f64.to_bits(), "static default must stay the legacy 1e-6");
        assert_eq!(
            ProcessingDynamics::Static.innovation_gate_scale().to_bits(),
            1.0_f64.to_bits(),
            "static gate scale must be exactly 1.0"
        );
    }

    #[test]
    fn kinematic_q_is_at_least_100x_static_and_gate_wider() {
        let qs = ProcessingDynamics::Static.q_accel_default();
        let qk = ProcessingDynamics::Kinematic.q_accel_default();
        assert!(qk >= 100.0 * qs, "kinematic Q must exceed static by >=100x");
        assert!(qk == 1.0, "kinematic q documented as ~1.0");
        assert!(ProcessingDynamics::Kinematic.innovation_gate_scale() > 1.0);
    }

    #[test]
    fn kinematic_sep_limit_is_monotone_with_sigma_and_bounded() {
        // Monotone non-decreasing in each sigma until the cap.
        let mut prev = 0.0;
        for i in 0..20 {
            let s = 0.01 * (1 << i) as f64;
            let lim = kinematic_sep_limit_m(KIN_DISAGREE_K_SIGMA, s, s, KIN_THRESHOLD_FLOOR_M, KIN_THRESHOLD_CAP_M);
            assert!(lim >= prev, "limit must grow with sigma: {s} -> {lim}");
            assert!(
                (KIN_THRESHOLD_FLOOR_M..=KIN_THRESHOLD_CAP_M).contains(&lim),
                "limit must stay bounded at sigma {s}"
            );
            prev = lim;
        }
        // Exact values at representative points.
        assert!(
            (kinematic_sep_limit_m(6.0, 0.2, 0.1, KIN_THRESHOLD_FLOOR_M, KIN_THRESHOLD_CAP_M) - 1.2).abs() < 1e-12,
            "6 sigma of the larger sigma"
        );
        // Floor equals the audited static bound: tiny sigmas never make
        // the kinematic rule stricter than the validated tolerance.
        assert_eq!(
            kinematic_sep_limit_m(6.0, 1e-6, 1e-6, KIN_THRESHOLD_FLOOR_M, KIN_THRESHOLD_CAP_M)
                .to_bits(),
            KIN_THRESHOLD_FLOOR_M.to_bits()
        );
        assert_eq!(
            kinematic_sep_limit_m(6.0, 1e3, 1.0, KIN_THRESHOLD_FLOOR_M, KIN_THRESHOLD_CAP_M)
                .to_bits(),
            KIN_THRESHOLD_CAP_M.to_bits()
        );
        // Larger of the two sigmas drives the limit.
        assert_eq!(
            kinematic_sep_limit_m(6.0, 0.05, 0.30, 0.02, 10.0)
                .to_bits(),
            kinematic_sep_limit_m(6.0, 0.30, 0.05, 0.02, 10.0).to_bits()
        );
    }

    #[test]
    fn is_kinomatic_flag_tracks_variant() {
        assert!(!ProcessingDynamics::Static.is_kinematic());
        assert!(ProcessingDynamics::Kinematic.is_kinematic());
    }

    #[test]
    fn gate_scale_multiplies_documented_threshold() {
        let eff = ROBUST_INNOVATION_THRESHOLD * ProcessingDynamics::Kinematic.innovation_gate_scale();
        assert_eq!(eff.to_bits(), 27.0_f64.to_bits());
        assert_eq!(
            (ROBUST_INNOVATION_THRESHOLD * ProcessingDynamics::Static.innovation_gate_scale())
                .to_bits(),
            ROBUST_INNOVATION_THRESHOLD.to_bits()
        );
    }
}
