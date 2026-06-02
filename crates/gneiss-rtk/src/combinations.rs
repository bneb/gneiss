const LIGHT_SPEED: f64 = 299792458.0;

/// Calculates the Wide-Lane (WL) carrier phase combination in meters.
/// L_wl = (f1 * L1 - f2 * L2) / (f1 - f2)
/// Assumes L1 and L2 are in meters.
pub fn wide_lane_phase(l1: f64, l2: f64, f1: f64, f2: f64) -> f64 {
    (f1 * l1 - f2 * l2) / (f1 - f2)
}

/// Calculates the Narrow-Lane (NL) pseudorange combination in meters.
/// P_nl = (f1 * P1 + f2 * P2) / (f1 + f2)
pub fn narrow_lane_pseudorange(p1: f64, p2: f64, f1: f64, f2: f64) -> f64 {
    (f1 * p1 + f2 * p2) / (f1 + f2)
}

/// Calculates the Melbourne-Wübbena (MW) combination in meters.
/// MW = L_wl - P_nl
/// This combination is geometry-free, ionosphere-free, and clock-free.
/// It isolates the Wide-Lane ambiguity: MW = lambda_wl * N_wl + biases + noise.
pub fn melbourne_wubbena(l1: f64, l2: f64, p1: f64, p2: f64, f1: f64, f2: f64) -> f64 {
    wide_lane_phase(l1, l2, f1, f2) - narrow_lane_pseudorange(p1, p2, f1, f2)
}

/// Computes the Ionosphere-Free (IF) linear combination.
/// IF = (f1^2 * v1 - f2^2 * v2) / (f1^2 - f2^2)
/// Can be applied to either pseudoranges or carrier phases (in meters).
pub fn iono_free(v1: f64, v2: f64, f1: f64, f2: f64) -> f64 {
    let f1_2 = f1 * f1;
    let f2_2 = f2 * f2;
    (f1_2 * v1 - f2_2 * v2) / (f1_2 - f2_2)
}

/// Wide-Lane Wavelength
pub fn lambda_wl(f1: f64, f2: f64) -> f64 {
    LIGHT_SPEED / (f1 - f2)
}

/// Narrow-Lane Wavelength
pub fn lambda_nl(f1: f64, f2: f64) -> f64 {
    LIGHT_SPEED / (f1 + f2)
}

#[cfg(test)]
mod tests {
    use super::*;

    const F1_GPS: f64 = 1575.42e6;
    const F2_GPS: f64 = 1227.60e6;

    #[test]
    fn test_melbourne_wubbena_isolates_wl_ambiguity() {
        // TDD: Prove that the MW combination isolates the Wide-Lane ambiguity,
        // stripping out true geometry, clock errors, and ionosphere.
        
        let true_range = 20_000_000.0;
        let cdt = 300_000.0; // 1 ms clock error
        let iono_l1 = 5.0; // 5 meters iono delay on L1
        let f_ratio = (F1_GPS / F2_GPS).powi(2);
        let iono_l2 = iono_l1 * f_ratio; // Iono is dispersive (1/f^2)
        
        // True ambiguities (in cycles)
        let n1_cycles = 1000.0;
        let n2_cycles = 1000.0 + 5.0; // N_wl = N1 - N2 = -5.0
        let n_wl_cycles = n1_cycles - n2_cycles;

        let lambda_1 = LIGHT_SPEED / F1_GPS;
        let lambda_2 = LIGHT_SPEED / F2_GPS;

        // Construct raw measurements
        // P = rho + cdt + I
        let p1 = true_range + cdt + iono_l1;
        let p2 = true_range + cdt + iono_l2;

        // L = rho + cdt - I + lambda * N
        let l1 = true_range + cdt - iono_l1 + lambda_1 * n1_cycles;
        let l2 = true_range + cdt - iono_l2 + lambda_2 * n2_cycles;

        let mw = melbourne_wubbena(l1, l2, p1, p2, F1_GPS, F2_GPS);
        let expected_mw = lambda_wl(F1_GPS, F2_GPS) * n_wl_cycles;

        assert!((mw - expected_mw).abs() < 1e-6, "MW combination failed to isolate WL ambiguity. Got: {}, Expected: {}", mw, expected_mw);
    }
    
    #[test]
    fn test_iono_free_combination() {
        let true_range = 20_000_000.0;
        let cdt = 300_000.0; 
        let iono_l1 = 5.0; 
        let f_ratio = (F1_GPS / F2_GPS).powi(2);
        let iono_l2 = iono_l1 * f_ratio;
        
        let p1 = true_range + cdt + iono_l1;
        let p2 = true_range + cdt + iono_l2;
        
        let p_if = iono_free(p1, p2, F1_GPS, F2_GPS);
        
        // IF combination deletes iono completely
        assert!((p_if - (true_range + cdt)).abs() < 1e-6);
    }

    #[test]
    fn test_narrow_lane_cascade_math() {
        // TDD: Prove that if we know N_wl, we can substitute it into the IF phase combination
        // to isolate N1 multiplied by lambda_nl.
        
        let true_range = 20_000_000.0;
        let cdt = 300_000.0; 
        let iono_l1 = 5.0; 
        let f_ratio = (F1_GPS / F2_GPS).powi(2);
        let iono_l2 = iono_l1 * f_ratio;
        
        let n1_cycles = 1000.0;
        let n2_cycles = 1005.0; 
        let n_wl_cycles = n1_cycles - n2_cycles; // -5.0

        let lambda_1 = LIGHT_SPEED / F1_GPS;
        let lambda_2 = LIGHT_SPEED / F2_GPS;

        let l1 = true_range + cdt - iono_l1 + lambda_1 * n1_cycles;
        let l2 = true_range + cdt - iono_l2 + lambda_2 * n2_cycles;

        let l_if = iono_free(l1, l2, F1_GPS, F2_GPS);

        // Derivation: L_if = rho + cdt + lambda_nl * N1 + (f2 / (f1 - f2)) * lambda_nl * N_wl
        let lambda_nl = lambda_nl(F1_GPS, F2_GPS);
        let wl_correction = (F2_GPS / (F1_GPS - F2_GPS)) * lambda_nl * n_wl_cycles;
        
        let l_if_corrected = l_if - wl_correction;
        let expected_l_if_corrected = true_range + cdt + lambda_nl * n1_cycles;

        assert!((l_if_corrected - expected_l_if_corrected).abs() < 1e-6, 
            "NL cascade failed. Corrected: {}, Expected: {}", l_if_corrected, expected_l_if_corrected);
    }
}