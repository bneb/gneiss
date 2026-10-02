//! Adversarial and property tests for the observation-variance models.
//!
//! Split out of `variance.rs` to keep that file under the 500-line budget. A
//! variance feeds a Cholesky factorisation, so the non-negativity and
//! finiteness invariants asserted here are load-bearing.

#[cfg(test)]
mod tests {
    use super::super::*;
    use core::f64::consts::FRAC_PI_2;

    // -----------------------------------------------------------------------
    // Property tests. A variance enters a Cholesky factorisation: a negative
    // value corrupts the factorisation *silently*, so the non-negativity and
    // finiteness invariants below are load-bearing, not cosmetic.
    // -----------------------------------------------------------------------

    const L1: f64 = 0.190_293_672_798_364_87;
    const L2: f64 = 0.244_210_213_424_568_25;

    /// Deterministic 64-bit LCG producing a repeatable sweep of [lo, hi].
    fn sweep(lo: f64, hi: f64, n: usize, seed: u64) -> impl Iterator<Item = f64> {
        let span = hi - lo;
        let mut state = seed | 1;
        (0..n).map(move |_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            lo + span * (((state >> 11) as f64) / ((1u64 << 53) as f64))
        })
    }

    /// No admissible combination may produce a negative, zero-but-finite-by-
    /// accident, NaN or infinite variance.
    #[test]
    fn variance_is_always_finite_and_non_negative() {
        let (a, b) = (0.1414213562373095_f64, 0.1414213562373095_f64);
        for el in sweep(-core::f64::consts::PI, core::f64::consts::PI, 512, 7) {
            for snr in [-200.0_f64, -20.0, 0.0, 10.0, 40.0, 80.0, 200.0] {
                let v = sigma_snr_variance(snr, el, false);
                assert!(v.is_finite() && v >= 0.0, "code snr={snr} el={el} -> {v}");
                let v = sigma_snr_variance(snr, el, true);
                assert!(v.is_finite() && v >= 0.0, "phase snr={snr} el={el} -> {v}");
                let v = sigma_snr_variance_with_coeffs(el, snr, a, b);
                assert!(v.is_finite() && v >= 0.0, "coeffs snr={snr} el={el} -> {v}");
            }
            let v = snr_variance_scale(30.0, a, b) * elevation_variance_scale(el);
            assert!(v.is_finite() && v >= 0.0, "legacy el={el} -> {v}");
        }
        // Phase variance scales with 1/lambda; any positive wavelength is fine.
        for lam in sweep(1e-3, 1.0, 128, 11) {
            let v = sigma_snr_phase_variance(35.0, 0.8, lam);
            assert!(v.is_finite() && v >= 0.0, "lambda={lam} -> {v}");
        }
    }

    /// A negative wavelength cannot occur physically, but must not produce a
    /// negative variance if one is ever passed: a and b enter squared.
    #[test]
    fn negative_or_degenerate_wavelength_never_yields_a_negative_variance() {
        for lam in [-L1, -0.0, 0.0] {
            let v = sigma_snr_phase_variance(35.0, FRAC_PI_2, lam);
            assert!(v >= 0.0, "lambda={lam} -> {v}");
            assert!(!v.is_nan(), "lambda={lam} produced NaN");
        }
    }

    /// Monotonicity: raising elevation or SNR must never increase the variance.
    /// Sampled on a dense grid so a small non-monotone excursion cannot hide
    /// between samples.
    #[test]
    fn variance_is_monotone_non_increasing_in_elevation_and_snr() {
        let mut prev = f64::INFINITY;
        for deg in 0..=90 {
            let v = sigma_snr_variance(38.0, (deg as f64).to_radians(), false);
            assert!(v <= prev + 1e-18, "code rose at {deg} deg: {prev} -> {v}");
            prev = v;
        }
        prev = f64::INFINITY;
        for db in (-200..=200).step_by(1) {
            let v = sigma_snr_variance(db as f64, 0.9, false);
            assert!(v <= prev + 1e-18, "code rose at {db} dB-Hz: {prev} -> {v}");
            prev = v;
        }
        prev = f64::INFINITY;
        for db in (-200..=200).step_by(1) {
            let v = sigma_snr_variance(db as f64, 0.9, true);
            assert!(v <= prev + 1e-18, "phase rose at {db} dB-Hz: {prev} -> {v}");
            prev = v;
        }
        // The horizon-regularised factor is flat below the horizon by design.
        let (a, b) = (0.1414213562373095_f64, 0.1414213562373095_f64);
        assert_eq!(elevation_factor(-0.5, a, b), elevation_factor(-0.2, a, b));
        assert!(elevation_factor(0.0, a, b) >= elevation_factor(FRAC_PI_2, a, b));
    }

    /// The minimum must sit at the zenith and be strictly positive: a variance
    /// of exactly zero at the zenith would make the weighting matrix singular.
    #[test]
    fn the_minimum_variance_is_at_the_zenith_and_is_strictly_positive() {
        for snr in [-50.0_f64, 0.0, 20.0, 40.0, 60.0, 120.0] {
            let zenith = sigma_snr_variance(snr, FRAC_PI_2, false);
            assert!(zenith > 0.0, "snr={snr}: zenith variance {zenith} is singular");
            for deg in [-90, -10, 0, 5, 15, 45, 75, 89] {
                let v = sigma_snr_variance(snr, (deg as f64).to_radians(), false);
                assert!(
                    v > zenith,
                    "snr={snr}: el={deg} deg {v} is not above the zenith {zenith}"
                );
            }
        }
    }

    /// Homogeneity: `elevation_factor` is quadratic in (a, b), so scaling both
    /// coefficients by `s` must scale the factor by exactly `s^2`.
    #[test]
    fn elevation_factor_is_homogeneous_of_degree_two_in_the_coefficients() {
        let (a, b) = (0.0021213203435596424_f64, 0.008);
        for el in [0.0_f64, 0.3, 1.0, FRAC_PI_2] {
            let base = elevation_factor(el, a, b);
            for s in [0.5_f64, 2.0, 10.0] {
                let got = elevation_factor(el, a * s, b * s);
                let want = base * s * s;
                assert!(
                    (got - want).abs() <= 8.0 * f64::EPSILON * want,
                    "el={el} s={s}: {got} != {want}"
                );
            }
        }
    }

    /// Exact golden vector at the zenith, with the arithmetic written out.
    /// sin(5 deg) = 0.08715574274765817, squared = 0.0075961234938959691.
    /// At el = pi/2 the denominator is 1 + sin^2 5 deg = 1.00759612349389597.
    #[test]
    fn zenith_elevation_factor_matches_hand_arithmetic() {
        const SIN_5_DEG: f64 = 0.087_155_742_747_658_17;
        let sin2_5 = SIN_5_DEG * SIN_5_DEG;
        assert!((sin2_5 - 0.007_596_123_493_895_995).abs() < 1e-15);
        let (a, b) = (0.1_f64 * core::f64::consts::SQRT_2, 0.1 * core::f64::consts::SQRT_2);
        // a = b = 0.1*sqrt(2), so a^2 = b^2 = 0.01 * 2 = 0.02 (to 1 ulp).
        assert!((a * a - 0.02).abs() <= 2.0 * f64::EPSILON * 0.02, "a^2 = {}", a * a);
        assert!((b * b - 0.02).abs() <= 2.0 * f64::EPSILON * 0.02, "b^2 = {}", b * b);
        let want = 0.02 + 0.02 / (1.0 + sin2_5);
        let got = elevation_factor(FRAC_PI_2, a, b);
        assert!(
            (got - want).abs() <= 4.0 * f64::EPSILON * want,
            "got {got} want {want}"
        );
    }

    /// `snr_factor` must stay inside [1, 1000] and reach 1 exactly at the top.
    #[test]
    fn snr_factor_is_bounded_and_reaches_the_floor() {
        for db in (-300..=300).step_by(1) {
            let f = snr_factor(db as f64);
            assert!(f.is_finite(), "snr={db} -> {f}");
            assert!((1.0..1000.0).contains(&f), "snr={db} -> {f} outside [1, 1000)");
        }
        // Above the cutoff it is exactly the floor, bit for bit.
        assert_eq!(snr_factor(200.0), 1.0);
        assert_eq!(snr_factor(1000.0), 1.0);
        // Monotone non-increasing across the whole range.
        let mut prev = f64::INFINITY;
        for db in (-300..=300).step_by(1) {
            let f = snr_factor(db as f64);
            assert!(f <= prev, "snr_factor rose at {db} dB-Hz: {prev} -> {f}");
            prev = f;
        }
    }

    /// Phase variance must scale as 1/lambda^2 exactly, since
    /// `sigma_cycles = 0.003/lambda` and the model is quadratic in (a, b).
    #[test]
    fn phase_variance_scales_as_one_over_lambda_squared() {
        for snr in [10.0_f64, 25.0, 45.0, 70.0] {
            let v1 = sigma_snr_phase_variance(snr, 0.6, L1);
            let v2 = sigma_snr_phase_variance(snr, 0.6, L2);
            let ratio = (L2 / L1) * (L2 / L1);
            assert!(
                (v1 / v2 - ratio).abs() < 1e-9 * ratio,
                "snr={snr}: ratio {} want {ratio}",
                v1 / v2
            );
            // Longer wavelength => fewer cycles per metre => smaller variance.
            assert!(v1 > v2, "L1 must have the larger cycles^2 variance");
        }
    }

    /// Adversarial inputs must not produce a *negative* variance. NaN is
    /// tolerated (and documented) but must never silently become negative.
    #[test]
    fn adversarial_inputs_never_produce_a_negative_variance() {
        let (a, b) = (0.1414213562373095_f64, 0.1414213562373095_f64);
        for el in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1e9, 1e9] {
            let v = elevation_factor(el, a, b);
            assert!(v >= 0.0 || v.is_nan(), "elevation_factor({el}) -> {v}");
            let v = elevation_variance_scale(el);
            assert!(v >= 0.0 || v.is_nan(), "elevation_variance_scale({el}) -> {v}");
        }
        for snr in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1e9, 1e9] {
            let v = snr_variance_scale(snr, a, b);
            assert!(v >= 0.0 || v.is_nan(), "snr_variance_scale({snr}) -> {v}");
        }
        for (x, y) in [(f64::NAN, 1.0), (1.0, f64::NAN), (f64::INFINITY, f64::NAN)] {
            let v = sigma_snr_variance_with_coeffs(x, y, a, b);
            assert!(v >= 0.0 || v.is_nan(), "sigma_snr_variance_with_coeffs({x},{y}) -> {v}");
        }
        // Negative coefficients must not flip the sign either (they are squared).
        let v = sigma_snr_variance_with_coeffs(0.5, 30.0, -a, -b);
        assert!(v > 0.0, "negative coefficients gave {v}");
        assert_eq!(
            v,
            sigma_snr_variance_with_coeffs(0.5, 30.0, a, b),
            "coefficient sign must not matter"
        );
    }

    /// `observation_variance` is the product of the SNR and elevation terms and
    /// must therefore be non-negative and non-increasing in both.
    #[test]
    fn observation_variance_is_non_negative_and_monotone() {
        let (a, b) = (1.0_f64, 150.0);
        let mut prev = f64::INFINITY;
        for deg in (0..=90).step_by(5) {
            let v = observation_variance(35.0, (deg as f64).to_radians(), a, b);
            assert!(v >= 0.0 && v.is_finite(), "el={deg} -> {v}");
            assert!(v <= prev + 1e-15, "rose at {deg} deg: {prev} -> {v}");
            prev = v;
        }
        prev = f64::INFINITY;
        for db in (0..=60).step_by(5) {
            let v = observation_variance(db as f64, 0.7, a, b);
            assert!(v <= prev + 1e-15, "rose at {db} dB-Hz: {prev} -> {v}");
            prev = v;
        }
    }
}
