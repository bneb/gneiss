/// Computes the dynamic measurement variance for a GNSS observation based on SNR and Elevation.
/// Uses the SIGMA-Epsilon formulation.
/// 
/// `snr_dbhz`: Signal-to-Noise Ratio in dB-Hz
/// `elevation_rad`: Elevation angle of the satellite in radians
/// `base_variance`: The theoretical minimum variance of the measurement (e.g., 0.0001 for Carrier Phase, 9.0 for Pseudorange)
pub fn dynamic_variance(snr_dbhz: f64, elevation_rad: f64, base_variance: f64) -> f64 {
    let snr_clamped = snr_dbhz.max(25.0).min(50.0);
    let snr_scale = libm::pow(10.0, (45.0 - snr_clamped) / 10.0).max(1.0).min(100.0);
    
    let sin_el = elevation_rad.sin().max(0.1);
    let el_scale = 1.0 / (sin_el * sin_el);
    
    base_variance * snr_scale * el_scale
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dynamic_variance_high_quality() {
        let snr = 50.0; // Excellent SNR
        let el = core::f64::consts::FRAC_PI_2; // 90 degrees (Zenith)
        let base_var = 0.0001; // 1 cm^2 for Carrier Phase

        let var = dynamic_variance(snr, el, base_var);
        
        assert!((var - base_var).abs() < 0.00001, "Expected roughly {}, got {}", base_var, var);
    }

    #[test]
    fn test_dynamic_variance_low_quality() {
        let snr = 25.0; // Terrible SNR
        let el = 0.1745; // 10 degrees (very low elevation)
        let base_var = 0.0001;

        let var = dynamic_variance(snr, el, base_var);
        
        assert!(var > base_var * 100.0, "Expected variance inflation > {}, got {}", base_var * 100.0, var);
    }
}
