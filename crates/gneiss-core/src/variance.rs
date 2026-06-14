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

#[cfg(test)]
mod tests {
    use super::*;

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
