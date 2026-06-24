/// Computes the dynamic measurement variance for a GNSS observation based on SNR and Elevation.
/// Uses the SIGMA-Epsilon formulation.
///
/// `snr_dbhz`: Signal-to-Noise Ratio in dB-Hz
/// `elevation_rad`: Elevation angle of the satellite in radians
/// `base_variance`: The theoretical minimum variance of the measurement (e.g., 0.0001 for Carrier Phase, 9.0 for Pseudorange)
pub fn dynamic_variance(snr_dbhz: f64, elevation_rad: f64, base_variance: f64) -> f64 {
    let snr_clamped = snr_dbhz.clamp(25.0, 50.0);
    let snr_scale = libm::pow(10.0, (45.0 - snr_clamped) / 10.0).clamp(1.0, 100.0);

    let sin_el = elevation_rad.sin().max(0.1);
    let el_scale = 1.0 / (sin_el * sin_el);

    base_variance * snr_scale * el_scale
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dynamic_variance_high_quality() {
        let snr = 50.0; // High SNR
        let el = core::f64::consts::FRAC_PI_2; // 90 degrees (Zenith)
        let base_var = 0.0001; // 1 cm^2 for Carrier Phase

        let var = dynamic_variance(snr, el, base_var);

        assert!(
            (var - base_var).abs() < 0.00001,
            "Expected roughly {}, got {}",
            base_var,
            var
        );
    }

    #[test]
    fn test_dynamic_variance_low_quality() {
        let snr = 25.0; // Terrible SNR
        let el = 0.1745; // 10 degrees (very low elevation)
        let base_var = 0.0001;

        let var = dynamic_variance(snr, el, base_var);

        assert!(
            var > base_var * 100.0,
            "Expected variance inflation > {}, got {}",
            base_var * 100.0,
            var
        );
    }

    #[test]
    fn test_dynamic_variance_boundary_snr_low_clamp() {
        // SNR below 25 should clamp to 25
        let el = core::f64::consts::FRAC_PI_2;
        let base_var = 1.0;
        let v_low = dynamic_variance(10.0, el, base_var);
        let v_clamp = dynamic_variance(25.0, el, base_var);
        assert!((v_low - v_clamp).abs() < 1e-12, "SNR below 25 should clamp to 25");
    }

    #[test]
    fn test_dynamic_variance_boundary_snr_high_clamp() {
        // SNR above 50 should clamp to 50
        let el = core::f64::consts::FRAC_PI_2;
        let base_var = 1.0;
        let v_high = dynamic_variance(60.0, el, base_var);
        let v_clamp = dynamic_variance(50.0, el, base_var);
        assert!((v_high - v_clamp).abs() < 1e-12, "SNR above 50 should clamp to 50");
    }

    #[test]
    fn test_dynamic_variance_snr_scale_min_clamp() {
        // When snr >= 45, snr_scale = 10^((45-45)/10) = 10^0 = 1.0 (clamped to min 1.0)
        let el = core::f64::consts::FRAC_PI_2;
        let base_var = 1.0;
        let v = dynamic_variance(45.0, el, base_var);
        let expected = base_var * 1.0 * 1.0; // snr_scale=1, el_scale=1/sin^2(pi/2)=1
        assert!((v - expected).abs() < 1e-12);
    }

    #[test]
    fn test_dynamic_variance_snr_scale_max_clamp() {
        // When snr = 25, snr_scale = 10^((45-25)/10) = 10^2 = 100, clamped to max 100
        let el = core::f64::consts::FRAC_PI_2;
        let base_var = 1.0;
        let v = dynamic_variance(25.0, el, base_var);
        let expected = base_var * 100.0 * 1.0;
        assert!((v - expected).abs() < 1e-12);
    }

    #[test]
    fn test_dynamic_variance_low_elevation_clamp() {
        // Elevation below 0.1 rad in sine should clamp sin_el to 0.1
        let snr = 50.0;
        let base_var = 1.0;
        let v_very_low = dynamic_variance(snr, 0.01, base_var);
        let v_at_clamp = dynamic_variance(snr, 0.1001, base_var);
        // At very low elevation (0.01 rad), sin(0.01) ≈ 0.01, which clamps to 0.1
        // So we get the same result as sin(0.1) ≈ 0.0998, actually 0.1001 rad sin ≈ 0.1
        // The clamp happens at 0.1, so sin(0.01) is clamped to 0.1
        // sin(0.1001) ≈ 0.1, so they should be close
        assert!((v_very_low - v_at_clamp).abs() < 1e-6);
    }

    #[test]
    fn test_dynamic_variance_pseudorange_base() {
        let snr = 40.0;
        let el = core::f64::consts::FRAC_PI_2;
        let base_var = 9.0; // Pseudorange base
        let var = dynamic_variance(snr, el, base_var);
        // SNR = 45 - 40 = 5, snr_scale = 10^(5/10) ≈ 3.162
        let expected_snr_scale = 10.0_f64.powf((45.0 - 40.0) / 10.0);
        let expected = base_var * expected_snr_scale;
        assert!((var - expected).abs() < 1e-10);
    }

    #[test]
    fn test_dynamic_variance_mid_elevation() {
        let snr = 50.0;
        let el = 0.5; // ~28.6 degrees, sin ≈ 0.479
        let base_var = 1.0;
        let var = dynamic_variance(snr, el, base_var);
        let sin_el = el.sin();
        let expected = base_var * 1.0 * (1.0 / (sin_el * sin_el));
        assert!((var - expected).abs() < 1e-10);
    }

    #[test]
    fn test_dynamic_variance_zero_base_variance() {
        let var = dynamic_variance(40.0, 1.0, 0.0);
        assert!((var - 0.0).abs() < 1e-12, "zero base variance gives zero");
    }

    #[test]
    fn test_dynamic_variance_extreme_elevation() {
        // Very low elevation just above the sin clamp (sin(0.1005) > 0.1)
        let el: f64 = 0.1005;
        let sin_el = el.sin().max(0.1); // match function's internal clamp
        let var = dynamic_variance(50.0, el, 1.0);
        let expected = 1.0 / (sin_el * sin_el);
        assert!((var - expected).abs() < 1e-10);

        // Effect of sin clamp at very low elevation - should match the clamped value
        let var_clamped = dynamic_variance(50.0, 0.001, 1.0);
        let clamp_expected = 1.0 / (0.1 * 0.1); // sin clamped to 0.1
        assert!((var_clamped - clamp_expected).abs() < 1e-6);
    }

    #[test]
    fn test_dynamic_variance_mid_snr_scale() {
        // SNR = 35: snr_scale = 10^((45-35)/10) = 10^1 = 10
        let var = dynamic_variance(35.0, core::f64::consts::FRAC_PI_2, 1.0);
        assert!((var - 10.0).abs() < 1e-10, "mid SNR gives scale of 10");
    }
}
