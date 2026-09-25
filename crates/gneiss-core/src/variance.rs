// This module is #![no_std] compatible - use libm

/// Computes SNR-based variance scaling factor using the a^2 + b^2 / 10^(SNR/10) model.
/// Returns a multiplicative variance factor.
/// snr_dbhz: Signal-to-noise ratio in dB-Hz
pub fn snr_variance_scale(snr_dbhz: f64, snr_a: f64, snr_b: f64) -> f64 {
    let snr_safe = if snr_dbhz < 10.0 { 10.0 } else { snr_dbhz };
    // var = a^2 + b^2 / 10^(SNR/10)
    snr_a * snr_a + (snr_b * snr_b) / libm::pow(10.0, snr_safe / 10.0)
}

/// Computes elevation-based variance scaling factor.
/// Returns 1/sin²(el), clamped to prevent singularity at horizon.
pub fn elevation_variance_scale(el_rad: f64) -> f64 {
    let sin_el = libm::sin(el_rad);
    let sin_el_safe = if sin_el < 0.1 { 0.1 } else { sin_el };
    1.0 / (sin_el_safe * sin_el_safe)
}

/// Computes combined observation variance factor.
pub fn observation_variance(snr_dbhz: f64, el_rad: f64, snr_a: f64, snr_b: f64) -> f64 {
    snr_variance_scale(snr_dbhz, snr_a, snr_b) * elevation_variance_scale(el_rad)
}

/// Elevation mapping factor with smooth regularisation at horizon (theta_0 = 5 deg).
/// Evaluates a^2 + b^2 / (sin^2(el) + sin^2(theta_0)).
pub fn elevation_factor(el_rad: f64, a: f64, b: f64) -> f64 {
    const SIN2_THETA_0: f64 = 0.007596123493895995;
    let s = libm::sin(el_rad);
    let s_pos = if s < 0.0 { 0.0 } else { s };
    let sin2_eff = s_pos * s_pos + SIN2_THETA_0;
    a * a + (b * b) / sin2_eff
}

/// Smooth SIGMA-SNR noise factor with logistic activation and algebraic saturation ceiling.
/// f_SNR(S) scales noise smoothly for S < 40.0 dB-Hz, saturating asymptotically at f_max = 1000.0.
pub fn snr_factor(snr_dbhz: f64) -> f64 {
    const S_NOM: f64 = 40.0;
    const TAU: f64 = 1.5;
    const F_MAX: f64 = 1000.0;
    if snr_dbhz.is_nan() {
        return 1.0;
    }
    let delta_s = S_NOM - snr_dbhz;
    let u = delta_s / TAU;
    let sigmoid = if u < -40.0 {
        0.0
    } else if u > 40.0 {
        1.0
    } else {
        1.0 / (1.0 + libm::exp(-u))
    };
    let power_term = if delta_s < -40.0 {
        0.0
    } else if delta_s > 60.0 {
        1.0e6
    } else {
        libm::pow(10.0, delta_s / 10.0)
    };
    let x = power_term * sigmoid;
    let m = F_MAX - 1.0;
    1.0 + (m * x) / (m + x)
}

/// Computes unified C1-smooth SIGMA-SNR variance given arbitrary channel noise parameters a and b.
pub fn sigma_snr_variance_with_coeffs(el_rad: f64, snr_dbhz: f64, a: f64, b: f64) -> f64 {
    elevation_factor(el_rad, a, b) * snr_factor(snr_dbhz)
}

/// Unified SIGMA-SNR observation variance (m^2) for pseudorange (is_phase=false) or carrier (is_phase=true).
pub fn sigma_snr_variance(snr_dbhz: f64, el_rad: f64, is_phase: bool) -> f64 {
    let (a, b) = if is_phase {
        (0.0021213203435596424, 0.0021213203435596424)
    } else {
        (0.1414213562373095, 0.1414213562373095)
    };
    sigma_snr_variance_with_coeffs(el_rad, snr_dbhz, a, b)
}

/// Unified SIGMA-SNR carrier phase variance in cycles^2 for wavelength lambda (meters).
pub fn sigma_snr_phase_variance(snr_dbhz: f64, el_rad: f64, lambda: f64) -> f64 {
    let sigma_cycles = 0.003 / lambda;
    let a = sigma_cycles * core::f64::consts::FRAC_1_SQRT_2;
    let b = a;
    sigma_snr_variance_with_coeffs(el_rad, snr_dbhz, a, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::FRAC_PI_2;

    #[test]
    fn test_tier1_analytical_golden_vectors() {
        // At zenith and nominal high SNR (50 dB-Hz), f_SNR ~ 1.0
        let v_zenith_code = sigma_snr_variance(50.0, FRAC_PI_2, false);
        let expected_zenith = 0.02 + 0.02 / (1.0 + 0.007596123493895995);
        assert!((v_zenith_code - expected_zenith).abs() < 1e-4);

        // At 40 dB-Hz, f_SNR ~ 1.49975
        let f40 = snr_factor(40.0);
        assert!((f40 - 1.49975).abs() < 1e-3);

        // At 30 dB-Hz, f_SNR ~ 10.888
        let f30 = snr_factor(30.0);
        assert!((f30 - 10.888).abs() < 0.05);

        // At 30 deg elevation and 35 dB-Hz
        let el_30deg = 30.0f64.to_radians();
        let v_mid = sigma_snr_variance(35.0, el_30deg, false);
        assert!(v_mid > v_zenith_code);
    }

    #[test]
    fn test_tier2_finite_difference_elevation_derivative() {
        let delta = 1e-4;
        let mut deg: f64 = 5.0;
        while deg <= 85.0 {
            let el = deg.to_radians();
            let v_plus = sigma_snr_variance(38.0, el + delta, false);
            let v_minus = sigma_snr_variance(38.0, el - delta, false);
            let d_sigma_d_theta = (v_plus - v_minus) / (2.0 * delta);
            assert!(
                d_sigma_d_theta <= 0.0,
                "d_sigma/d_theta must be <= 0 at {deg} deg: got {d_sigma_d_theta}"
            );
            deg += 5.0;
        }
    }

    #[test]
    fn test_tier2_finite_difference_snr_derivative() {
        let delta = 1e-4;
        let el = 45.0f64.to_radians();
        let mut snr = 10.0;
        while snr <= 50.0 {
            let v_plus = sigma_snr_variance(snr + delta, el, false);
            let v_minus = sigma_snr_variance(snr - delta, el, false);
            let d_sigma_d_snr = (v_plus - v_minus) / (2.0 * delta);
            assert!(
                d_sigma_d_snr <= 0.0,
                "d_sigma/d_snr must be <= 0 at {snr} dB-Hz: got {d_sigma_d_snr}"
            );
            snr += 2.5;
        }
    }

    #[test]
    fn test_tier3_stability_bounds_and_asymptotics() {
        // Horizon elevation (0 deg) and negative SNR (-20 dB-Hz)
        let v_worst = sigma_snr_variance(-20.0, 0.0, false);
        assert!(v_worst.is_finite(), "Worst-case variance must be finite");
        assert!(v_worst > 0.0, "Worst-case variance must be strictly positive");
        assert!(v_worst < 3000.0, "Worst-case variance must saturate below 3000.0");

        // Extreme positive SNR should approach nominal floor
        let v_high = sigma_snr_variance(80.0, FRAC_PI_2, false);
        assert!(v_high > 0.038 && v_high < 0.042);

        // NaN handling
        let v_nan = snr_factor(f64::NAN);
        assert_eq!(v_nan, 1.0);
    }

    #[test]
    fn test_phase_variance_scaling_with_wavelength() {
        let l1 = 0.19029367;
        let l2 = 0.24421021;
        let v_l1 = sigma_snr_phase_variance(45.0, FRAC_PI_2, l1);
        let v_l2 = sigma_snr_phase_variance(45.0, FRAC_PI_2, l2);
        // Larger wavelength -> fewer cycles per meter -> smaller variance in cycles^2
        assert!(v_l1 > v_l2);
    }

    #[test]
    fn test_variance_monotonic_with_snr() {
        // Higher SNR should give lower variance
        let a = 1.0;
        let b = 150.0;
        let v1 = snr_variance_scale(30.0, a, b);
        let v2 = snr_variance_scale(40.0, a, b);
        let v3 = snr_variance_scale(45.0, a, b);
        assert!(v1 > v2, "30 dBHz should have higher variance than 40 dBHz");
        assert!(v2 > v3, "40 dBHz should have higher variance than 45 dBHz");
    }

    #[test]
    fn test_variance_monotonic_with_elevation() {
        use core::f64::consts::FRAC_PI_2;
        // Higher elevation should give lower variance
        let v_low = elevation_variance_scale(0.2); // ~11.5 degrees
        let v_mid = elevation_variance_scale(0.5); // ~28.6 degrees
        let v_high = elevation_variance_scale(FRAC_PI_2); // 90 degrees (zenith)
        assert!(v_low > v_mid);
        assert!(v_mid > v_high);
        // At zenith, should be exactly 1.0
        assert!((v_high - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_variance_boundary_values() {
        let a = 1.0;
        let b = 150.0;
        // SNR below 10 clamps to 10
        let v_low = snr_variance_scale(5.0, a, b);
        let v_10 = snr_variance_scale(10.0, a, b);
        assert!((v_low - v_10).abs() < 1e-6, "SNR below 10 should clamp");
    }
}
