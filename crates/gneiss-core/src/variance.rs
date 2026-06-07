// This module is #![no_std] compatible - use libm

/// Computes SNR-based variance scaling factor.
/// Returns a multiplicative factor: higher SNR = lower variance.
/// snr_dbhz: Signal-to-noise ratio in dB-Hz
/// nominal_snr: Reference SNR (typically 45.0 dBHz)
pub fn snr_variance_scale(snr_dbhz: f64, nominal_snr: f64) -> f64 {
    let snr_clamped = snr_dbhz.clamp(25.0, 50.0);
    let scale = libm::pow(10.0, (nominal_snr - snr_clamped) / 10.0);
    scale.min(100.0)
}

/// Computes elevation-based variance scaling factor.
/// Returns 1/sin²(el), clamped to prevent singularity at horizon.
pub fn elevation_variance_scale(el_rad: f64) -> f64 {
    let sin_el = libm::sin(el_rad);
    let sin_el_safe = if sin_el < 0.1 { 0.1 } else { sin_el };
    1.0 / (sin_el_safe * sin_el_safe)
}

/// Computes combined observation variance factor.
pub fn observation_variance(snr_dbhz: f64, el_rad: f64, nominal_snr: f64) -> f64 {
    snr_variance_scale(snr_dbhz, nominal_snr) * elevation_variance_scale(el_rad)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_variance_monotonic_with_snr() {
        // Higher SNR should give lower variance
        let v1 = snr_variance_scale(30.0, 45.0);
        let v2 = snr_variance_scale(40.0, 45.0);
        let v3 = snr_variance_scale(45.0, 45.0);
        assert!(v1 > v2, "30 dBHz should have higher variance than 40 dBHz");
        assert!(v2 > v3, "40 dBHz should have higher variance than 45 dBHz");
        // At nominal, scale should be 1.0
        assert!((v3 - 1.0).abs() < 1e-6, "At nominal SNR, scale should be 1.0");
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
        // SNR below 25 clamps to 25
        let v_low = snr_variance_scale(10.0, 45.0);
        let v_25 = snr_variance_scale(25.0, 45.0);
        assert!((v_low - v_25).abs() < 1e-6, "SNR below 25 should clamp");
        // SNR above 50 clamps to 50
        let v_high = snr_variance_scale(60.0, 45.0);
        let v_50 = snr_variance_scale(50.0, 45.0);
        assert!((v_high - v_50).abs() < 1e-6, "SNR above 50 should clamp");
        // Scale should never exceed 100
        assert!(snr_variance_scale(25.0, 45.0) <= 100.0);
    }
}
