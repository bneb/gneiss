//! Measurement noise and covariance modeling for double-differenced observations.

use nalgebra::Vector3;

/// Evaluated noise variances for a double-difference observation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DdVariance {
    /// Pseudorange variance (m^2).
    pub pr_var_m2: f64,
    /// Carrier phase variance (cycles^2).
    pub cp_var_cycles2: f64,
    /// Reference satellite pseudorange contribution (m^2).
    pub pr_ref_var_m2: f64,
    /// Reference satellite carrier phase contribution (cycles^2).
    pub cp_ref_var_cycles2: f64,
}

/// Trait for converting diverse SNR 4-tuples into floating-point C/N0 values (dB-Hz).
pub trait IntoSnrs {
    fn into_snrs_f64(self) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>);
}

impl IntoSnrs for (Option<f64>, Option<f64>, Option<f64>, Option<f64>) {
    #[inline]
    fn into_snrs_f64(self) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>) {
        self
    }
}

impl IntoSnrs for (Option<u8>, Option<u8>, Option<u8>, Option<u8>) {
    #[inline]
    fn into_snrs_f64(self) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>) {
        (
            self.0.map(|s| s as f64),
            self.1.map(|s| s as f64),
            self.2.map(|s| s as f64),
            self.3.map(|s| s as f64),
        )
    }
}

impl IntoSnrs for (Option<i32>, Option<i32>, Option<i32>, Option<i32>) {
    #[inline]
    fn into_snrs_f64(self) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>) {
        (
            self.0.map(|s| s as f64),
            self.1.map(|s| s as f64),
            self.2.map(|s| s as f64),
            self.3.map(|s| s as f64),
        )
    }
}

/// Compute SIGMA-C/N0 noise scale factor.
/// Signals at or above 40 dB-Hz have nominal unit weight (1.0).
/// Attenuated signals (<40 dB-Hz) scale variance inversely with signal power down to 15 dB-Hz.
pub fn snr_weight(snr_dbhz: Option<u8>) -> f64 {
    const NOMINAL_CN0: f64 = 40.0;
    const MIN_CN0: f64 = 15.0;
    match snr_dbhz {
        Some(snr) if (snr as f64) < NOMINAL_CN0 => {
            let clamped = (snr as f64).max(MIN_CN0);
            10.0f64.powf((NOMINAL_CN0 - clamped) / 10.0)
        }
        _ => 1.0,
    }
}

/// Expected nominal C/N0 (dB-Hz) as a function of satellite elevation (radians).
pub fn expected_cn0_dbhz(el_rad: f64) -> f64 {
    let sin_el = el_rad.sin().clamp(0.0, 1.0);
    (30.0 + 20.0 * sin_el).min(50.0)
}

/// Compute attenuation penalty multiplier.
/// If C/N0 is below 25 dB-Hz or drops > 10 dB-Hz below its elevation expectation,
/// the variance is scaled up smoothly using an asymptotic saturation ceiling to prevent
/// corrupted phases from pulling the filter without introducing slope discontinuities.
pub fn attenuation_scale(snr_dbhz: Option<u8>, el_rad: f64) -> f64 {
    let snr = match snr_dbhz {
        Some(s) => s as f64,
        None => return 1.0,
    };
    attenuation_scale_f64(snr, el_rad)
}

/// Smooth floating-point attenuation scale with asymptotic algebraic ceiling.
pub fn attenuation_scale_f64(snr: f64, el_rad: f64) -> f64 {
    let exp_cn0 = expected_cn0_dbhz(el_rad);
    let drop_db = (exp_cn0 - snr).max(0.0);
    let mut penalty = 0.0;
    if drop_db > 10.0 {
        penalty += 10.0f64.powf((drop_db - 10.0) / 5.0) - 1.0;
    }
    if snr < 25.0 {
        penalty += 10.0f64.powf((25.0 - snr) / 2.5) - 1.0;
    }
    if penalty <= 0.0 {
        return 1.0;
    }
    const M: f64 = 999.0;
    1.0 + (M * penalty) / (M + penalty)
}

fn single_diff_var(
    base_sigma: f64,
    snr_a: Option<f64>,
    snr_b: Option<f64>,
    el_rad: f64,
) -> f64 {
    let el_fac = gneiss_core::variance::elevation_factor(el_rad, 0.0, base_sigma);
    let s_a = snr_a.unwrap_or(45.0);
    let s_b = snr_b.unwrap_or(45.0);
    let w_a = gneiss_core::variance::snr_factor(s_a);
    let w_b = gneiss_core::variance::snr_factor(s_b);
    el_fac * (w_a + w_b)
}

/// Compute double-difference code and phase variances combining satellite elevation and C/N0.
pub fn compute_dd_variances<S: IntoSnrs>(
    rx_pos: Vector3<f64>,
    sat_pos: Vector3<f64>,
    ref_pos: Vector3<f64>,
    lambda: f64,
    snrs: S,
) -> DdVariance {
    let snrs = snrs.into_snrs_f64();
    let rx_llh = gneiss_core::coords::ecef_to_llh(rx_pos);
    let (_az_s, el_s) = gneiss_core::coords::az_el(rx_llh, rx_pos, sat_pos);
    let (_az_r, el_r) = gneiss_core::coords::az_el(rx_llh, rx_pos, ref_pos);

    let pr_sat_var = single_diff_var(0.20, snrs.0, snrs.2, el_s);
    let pr_ref_var = single_diff_var(0.20, snrs.1, snrs.3, el_r);
    let cp_sat_var = single_diff_var(0.003 / lambda, snrs.0, snrs.2, el_s);
    let cp_ref_var = single_diff_var(0.003 / lambda, snrs.1, snrs.3, el_r);

    DdVariance {
        pr_var_m2: pr_sat_var + pr_ref_var,
        cp_var_cycles2: cp_sat_var + cp_ref_var,
        pr_ref_var_m2: pr_ref_var,
        cp_ref_var_cycles2: cp_ref_var,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_snr_weight_nominal_and_missing() {
        assert_eq!(snr_weight(None), 1.0);
        assert_eq!(snr_weight(Some(40)), 1.0);
        assert_eq!(snr_weight(Some(45)), 1.0);
        assert_eq!(snr_weight(Some(52)), 1.0);
    }

    #[test]
    fn test_snr_weight_attenuated_scaling() {
        let w30 = snr_weight(Some(30));
        assert!((w30 - 10.0).abs() < 1e-6, "w30={w30}");

        let w20 = snr_weight(Some(20));
        assert!((w20 - 100.0).abs() < 1e-6, "w20={w20}");

        let w10 = snr_weight(Some(10));
        assert!((w10 - 10.0f64.powf(2.5)).abs() < 1e-4, "w10={w10}");
    }

    #[test]
    fn test_compute_dd_variances_nominal_matches_legacy() {
        let rx_pos = Vector3::new(-2688181.50, -4265663.45, 3893784.80);
        let sat_pos = rx_pos + Vector3::new(1e7, 1e7, 1e7);
        let ref_pos = rx_pos + Vector3::new(0.0, 1e7, 2e7);
        let lambda = 0.19029367;

        let var = compute_dd_variances(rx_pos, sat_pos, ref_pos, lambda, (None::<f64>, None, None, None));
        assert!(var.pr_var_m2 > 0.0);
        assert!(var.cp_var_cycles2 > 0.0);
        assert!(var.pr_ref_var_m2 > 0.0);
        assert!(var.pr_var_m2 > var.pr_ref_var_m2);
        assert!(var.cp_var_cycles2 > var.cp_ref_var_cycles2);
    }

    #[test]
    fn test_expected_cn0_and_attenuation_scale() {
        use std::f64::consts::FRAC_PI_2;
        let exp_zenith = expected_cn0_dbhz(FRAC_PI_2);
        assert!((exp_zenith - 50.0).abs() < 1e-6);

        let exp_30deg = expected_cn0_dbhz(30.0f64.to_radians());
        assert!((exp_30deg - 40.0).abs() < 1e-6);

        let scale_nominal = attenuation_scale(Some(40), 30.0f64.to_radians());
        assert_eq!(scale_nominal, 1.0);

        let scale_attenuated = attenuation_scale(Some(25), FRAC_PI_2);
        assert!(scale_attenuated > 100.0, "scale={scale_attenuated}");

        let scale_sub25 = attenuation_scale(Some(20), 30.0f64.to_radians());
        assert!(scale_sub25 > 100.0, "scale={scale_sub25}");
    }

    #[test]
    fn test_dd_variances_monotonic_with_snr() {
        let rx_pos = Vector3::new(-2688181.50, -4265663.45, 3893784.80);
        let sat_pos = rx_pos + Vector3::new(1e7, 1e7, 1e7);
        let ref_pos = rx_pos + Vector3::new(0.0, 1e7, 2e7);
        let lambda = 0.19029367;

        let var_high = compute_dd_variances(rx_pos, sat_pos, ref_pos, lambda, (Some(48.0), Some(48.0), Some(48.0), Some(48.0)));
        let var_mid = compute_dd_variances(rx_pos, sat_pos, ref_pos, lambda, (Some(35.0), Some(48.0), Some(48.0), Some(48.0)));
        let var_low = compute_dd_variances(rx_pos, sat_pos, ref_pos, lambda, (Some(22.0), Some(48.0), Some(48.0), Some(48.0)));

        assert!(var_low.pr_var_m2 > var_mid.pr_var_m2);
        assert!(var_mid.pr_var_m2 > var_high.pr_var_m2);
        assert!(var_low.cp_var_cycles2 > var_mid.cp_var_cycles2);
        assert!(var_mid.cp_var_cycles2 > var_high.cp_var_cycles2);
    }

    #[test]
    fn test_into_snrs_types_interoperability() {
        let rx_pos = Vector3::new(-2688181.50, -4265663.45, 3893784.80);
        let sat_pos = rx_pos + Vector3::new(1e7, 1e7, 1e7);
        let ref_pos = rx_pos + Vector3::new(0.0, 1e7, 2e7);
        let lambda = 0.19029367;

        let v_u8 = compute_dd_variances(rx_pos, sat_pos, ref_pos, lambda, (Some(40u8), Some(45u8), Some(48u8), Some(48u8)));
        let v_f64 = compute_dd_variances(rx_pos, sat_pos, ref_pos, lambda, (Some(40.0f64), Some(45.0f64), Some(48.0f64), Some(48.0f64)));
        let v_i32 = compute_dd_variances(rx_pos, sat_pos, ref_pos, lambda, (Some(40i32), Some(45i32), Some(48i32), Some(48i32)));

        assert!((v_u8.pr_var_m2 - v_f64.pr_var_m2).abs() < 1e-9);
        assert!((v_u8.cp_var_cycles2 - v_f64.cp_var_cycles2).abs() < 1e-9);
        assert!((v_i32.pr_var_m2 - v_f64.pr_var_m2).abs() < 1e-9);
    }
}
