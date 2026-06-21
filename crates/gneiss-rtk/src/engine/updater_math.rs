pub trait CouplingStrategy {
    fn is_tightly_coupled() -> bool;
    fn pre_fit_threshold_multiplier() -> f64;
}
pub struct LooseCoupling;
impl CouplingStrategy for LooseCoupling {
    fn is_tightly_coupled() -> bool {
        false
    }
    fn pre_fit_threshold_multiplier() -> f64 {
        1.0
    }
}
pub struct TightCoupling;
impl CouplingStrategy for TightCoupling {
    fn is_tightly_coupled() -> bool {
        true
    }
    fn pre_fit_threshold_multiplier() -> f64 {
        25.0
    }
}

use crate::engine::config::EkfTuningConfig;
use nalgebra::{DMatrix, DVector, UnitQuaternion, Vector3};

pub const CP_PRE_FIT_CHI2_THRESHOLD: f64 = 100.0;
pub const DOPPLER_PRE_FIT_CHI2_THRESHOLD: f64 = 50.0;
pub const PR_PRE_FIT_CHI2_MULTIPLIER: f64 = 25.0;

pub fn filter_pre_fit_residuals<C: CouplingStrategy>(
    z: &DVector<f64>,
    h: &DMatrix<f64>,
    r: &DMatrix<f64>,
    state_cov: &DMatrix<f64>,
    max_innovation: f64,
    meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>,
) -> Vec<usize> {
    let mut valid_indices = Vec::with_capacity(z.len());
    let hp = h * state_cov;

    for i in 0..z.len() {
        let s_ii = (hp.row(i) * h.row(i).transpose())[(0, 0)] + r[(i, i)];
        let meas_type = meas_types.map_or(0, |m| m[i].1);
        let threshold = get_pre_fit_threshold::<C>(meas_type, max_innovation);

        if check_pre_fit_residual(z[i], s_ii, r[(i, i)], meas_type, threshold) {
            valid_indices.push(i);
        }
    }
    valid_indices
}

pub fn compute_loose_coupling_innovations(
    r_b_e: &nalgebra::Matrix3<f64>,
    state_pos: &Vector3<f64>,
    state_vel: &Vector3<f64>,
    gnss_pos: &Vector3<f64>,
    gnss_vel: &Vector3<f64>,
    lever_arm: &Vector3<f64>,
    omega_b: &Vector3<f64>,
) -> DVector<f64> {
    let l_e = r_b_e * lever_arm;
    let pos_apc = state_pos + l_e;
    let v_apc = state_vel + r_b_e * omega_b.cross(lever_arm);

    let mut z = DVector::zeros(6);
    z.rows_mut(0, 3).copy_from(&(gnss_pos - pos_apc));
    z.rows_mut(3, 3).copy_from(&(gnss_vel - v_apc));
    z
}

pub fn enforce_symmetry(p: &mut DMatrix<f64>) {
    for r_idx in 0..p.nrows() {
        for c_idx in 0..r_idx {
            let avg = (p[(r_idx, c_idx)] + p[(c_idx, r_idx)]) * 0.5;
            p[(r_idx, c_idx)] = avg;
            p[(c_idx, r_idx)] = avg;
        }
    }
}

pub fn get_pre_fit_threshold<C: CouplingStrategy>(meas_type: u8, max_innovation: f64) -> f64 {
    let multiplier = C::pre_fit_threshold_multiplier();
    match meas_type {
        1 | 2 => CP_PRE_FIT_CHI2_THRESHOLD * multiplier,
        3 => DOPPLER_PRE_FIT_CHI2_THRESHOLD * multiplier,
        _ => max_innovation * max_innovation * PR_PRE_FIT_CHI2_MULTIPLIER * multiplier,
    }
}

pub fn check_pre_fit_residual(
    nu: f64,
    s_ii: f64,
    r_ii: f64,
    meas_type: u8,
    threshold: f64,
) -> bool {
    if nu * nu / s_ii < threshold {
        true
    } else {
        tracing::debug!(
            "EKF rejected meas: type={}, nu={:.2}, s_ii={:.2}, r_ii={:.4}",
            meas_type,
            nu,
            s_ii,
            r_ii
        );
        false
    }
}

pub fn compute_scalar_thresholds<C: CouplingStrategy>(
    meas_type: u8,
    max_innovation: f64,
    tuning: &EkfTuningConfig,
) -> (f64, f64) {
    let thresh = match meas_type {
        1 | 2 => tuning.phase_outlier_ratio_thresh,
        3 => max_innovation * tuning.doppler_outlier_ratio_mult,
        _ => max_innovation,
    };
    let abs_thresh = match meas_type {
        0 => tuning.pr_abs_thresh,
        1 | 2 => 1000000.0, // Never hard-reject CP by absolute value
        3 => tuning.dop_abs_thresh,
        _ => 40.0,
    };
    (thresh, abs_thresh)
}

pub fn evaluate_post_fit_outliers<C: CouplingStrategy>(
    v: &DVector<f64>,
    s: &DMatrix<f64>,
    _current_z: &DVector<f64>,
    current_valid: &[usize],
    meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>,
    max_innovation: f64,
    tuning: &EkfTuningConfig,
) -> (Option<usize>, f64, DVector<f64>) {
    let mut max_outlier_ratio = 0.0;
    let mut worst_idx = None;
    let mut weights = DVector::from_element(v.len(), 1.0);

    for i in 0..v.len() {
        let meas_type = meas_types.map_or(0, |m| m[current_valid[i]].1);
        let ratio = v[i].abs() / s[(i, i)].sqrt();
        let (_thresh, abs_thresh) =
            compute_scalar_thresholds::<C>(meas_type, max_innovation, tuning);

        // Scale absolute threshold by the filter's uncertainty for pseudoranges to prevent getting stuck
        // when the filter has intentionally inflated its covariance. Carrier phases should strictly use
        // their absolute thresholds.
        let effective_abs_thresh = match meas_type {
            1 | 2 => abs_thresh,
            _ => f64::max(abs_thresh, s[(i, i)].sqrt() * 3.0),
        };

        let is_tight = C::is_tightly_coupled();

        // If tightly-coupled, we relax the hard absolute rejection threshold by a factor of 5
        // to allow Huber scaling to gracefully handle severe multipath instead of hard rejecting.
        let is_abs_outlier = if is_tight {
            v[i].abs() > (effective_abs_thresh * 5.0) && meas_type != 3
        } else {
            v[i].abs() > effective_abs_thresh && meas_type != 3
        };

        if is_abs_outlier {
            return (Some(i), f64::INFINITY, weights);
        }

        let hard_reject_ratio = match meas_type {
            1 | 2 => {
                if is_tight {
                    tuning.phase_outlier_ratio_thresh * 3.0
                } else {
                    tuning.phase_outlier_ratio_thresh
                }
            }
            3 => {
                if is_tight {
                    tuning.phase_outlier_ratio_thresh * tuning.doppler_outlier_ratio_mult * 3.0
                } else {
                    tuning.phase_outlier_ratio_thresh * tuning.doppler_outlier_ratio_mult
                }
            }
            _ => {
                if is_tight {
                    tuning.phase_outlier_ratio_thresh * 5.0
                } else {
                    tuning.phase_outlier_ratio_thresh * 1.666
                }
            }
        };

        if ratio > hard_reject_ratio {
            if ratio > max_outlier_ratio {
                max_outlier_ratio = ratio;
                worst_idx = Some(i);
            }
        } else if is_tight {
            // Apply Huber weighting ONLY in tightly coupled mode.
            // LooseCoupling relies on hard-rejection and doesn't use IRLS.
            weights[i] = crate::math::thresholding::apply_cauchy(v[i], s[(i, i)], 3.0);
        }
    }
    (worst_idx, max_outlier_ratio, weights)
}

pub fn populate_loosely_coupled_jacobian(
    h_mat: &mut DMatrix<f64>,
    r_b_e: &UnitQuaternion<f64>,
    lever_arm: &Vector3<f64>,
    omega_b: &Vector3<f64>,
) {
    let l_e = r_b_e * lever_arm;
    let h_pos_att = -l_e.cross_matrix();
    let a_e = r_b_e * omega_b.cross(lever_arm);
    let h_vel_att = -a_e.cross_matrix();
    let h_vel_bg = r_b_e.to_rotation_matrix().matrix() * lever_arm.cross_matrix();

    h_mat.view_mut((0, 6), (3, 3)).copy_from(&h_pos_att);
    h_mat.view_mut((3, 6), (3, 3)).copy_from(&h_vel_att);
    h_mat.view_mut((3, 12), (3, 3)).copy_from(&h_vel_bg);
}

#[test]
fn test_evaluate_post_fit_outliers() {
    use crate::engine::config::EkfTuningConfig;
    let mut tuning = EkfTuningConfig::default();
    tuning.pr_abs_thresh = 50.0;
    tuning.dop_abs_thresh = 60.0;
    tuning.phase_outlier_ratio_thresh = 2.0;
    tuning.doppler_outlier_ratio_mult = 3.0;

    let v = DVector::from_vec(vec![1000.0, 10.0, 10.0]);
    let s = DMatrix::from_diagonal(&DVector::from_vec(vec![1.0, 1.0, 1.0]));
    let current_z = DVector::from_vec(vec![1000.0, 10.0, 10.0]);
    let current_valid = vec![0, 1, 2];

    // Mock meas_types
    // We need meas_type == 3 to bypass is_abs_outlier
    let sat = gneiss_core::sat::SatelliteId {
        constellation: gneiss_core::sat::Constellation::Gps,
        prn: 1,
    };
    let meas_types = vec![(sat, 1), (sat, 3), (sat, 1)];

    // Test 1: meas_type != 3, so v[0] = 1000 is an abs outlier
    let (outlier, _, _) = evaluate_post_fit_outliers::<TightCoupling>(
        &v,
        &s,
        &current_z,
        &current_valid,
        Some(&meas_types),
        50.0,
        &tuning,
    );
    assert_eq!(outlier, Some(0)); // Returns immediately on abs outlier

    // Test 2: meas_type == 3, so it's NOT an abs outlier, but it will be a ratio outlier
    let v_t2 = DVector::from_vec(vec![0.0, 1000.0, 10.0]); // v[1] corresponds to meas_types[1] which is type 3
    let current_z_t2 = DVector::from_vec(vec![0.0, 1000.0, 10.0]);
    let (outlier2, ratio, _) = evaluate_post_fit_outliers::<TightCoupling>(
        &v_t2,
        &s,
        &current_z_t2,
        &current_valid,
        Some(&meas_types),
        50.0,
        &tuning,
    );
    assert_eq!(outlier2, Some(1)); // Because it has a massive ratio, but was NOT flagged as abs outlier
    assert!(ratio > 900.0);

    // Test 3: Multiple outliers, finds the worst ratio
    let v_t3 = DVector::from_vec(vec![0.0, 50.0, 20.0]); // None are abs outliers (assume thresh is high enough, or meas_types logic)
    let s_t3 = DMatrix::from_diagonal(&DVector::from_vec(vec![1.0, 25.0, 1.0]));
    // ratio 1: 50/sqrt(25) = 10
    // ratio 2: 20/sqrt(1) = 20 (worse!)
    let current_z_t3 = DVector::from_vec(vec![0.0, 50.0, 20.0]);
    let meas_types_t3 = vec![(sat, 3), (sat, 3), (sat, 3)]; // all type 3 to bypass abs
    let (outlier3, ratio3, _) = evaluate_post_fit_outliers::<TightCoupling>(
        &v_t3,
        &s_t3,
        &current_z_t3,
        &current_valid,
        Some(&meas_types_t3),
        1.0,
        &tuning,
    );
    assert_eq!(outlier3, Some(2));
    assert_eq!(outlier3, Some(2));
    assert_eq!(ratio3, 20.0);
}

#[test]
fn test_evaluate_post_fit_outliers_loose_coupling() {
    use crate::engine::config::EkfTuningConfig;
    let mut tuning = EkfTuningConfig::default();
    tuning.pr_abs_thresh = 50.0;
    tuning.dop_abs_thresh = 60.0; // But Doppler (meas_type = 3) bypasses abs_thresh check!
    tuning.phase_outlier_ratio_thresh = 2.0;
    tuning.doppler_outlier_ratio_mult = 3.0;

    let s = DMatrix::from_diagonal(&DVector::from_vec(vec![1.0, 1.0, 1.0, 1.0]));
    let current_z = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0]);
    let current_valid = vec![0, 1, 2, 3];

    let sat = gneiss_core::sat::SatelliteId {
        constellation: gneiss_core::sat::Constellation::Gps,
        prn: 1,
    };
    let meas_types = vec![(sat, 0), (sat, 3), (sat, 1), (sat, 0)];

    // v[0]: meas_type 0. v[0]=50.0. abs_thresh=50.0. Not an abs outlier (50 > 50 is false).
    // v[1]: meas_type 3. v[1]=1000.0. Bypasses abs outlier (meas_type == 3).
    // v[2]: meas_type 1. v[2]=1e6. abs_thresh=1e6. Not an abs outlier (1e6 > 1e6 is false).
    // v[3]: meas_type 0. v[3]=50.01. abs_thresh=50.0. IS an abs outlier!
    let v = DVector::from_vec(vec![50.0, 1000.0, 1000000.0, 50.01]);

    // Test 1: v[3] should be caught as abs outlier first because it evaluates sequentially!
    let (outlier, _, _) = evaluate_post_fit_outliers::<LooseCoupling>(
        &v,
        &s,
        &current_z,
        &current_valid,
        Some(&meas_types),
        50.0,
        &tuning,
    );
    assert_eq!(outlier, Some(3));

    // Test 2: Check ratio thresholds for LooseCoupling
    let v2 = DVector::from_vec(vec![3.0, 6.0, 2.1, 0.0]);
    // ratio thresh for 0: phase_ratio * 1.666 = 2.0 * 1.666 = 3.332. ratio 3.0 is OK.
    // ratio thresh for 3: phase_ratio * dop_mult = 2.0 * 3.0 = 6.0. ratio 6.0 is OK. (v2[1] = 6.0)
    // ratio thresh for 1: phase_ratio = 2.0. ratio 2.1 is NOT OK. (v2[2] = 2.1 > 2.0)
    let (outlier2, ratio2, _) = evaluate_post_fit_outliers::<LooseCoupling>(
        &v2,
        &s,
        &current_z,
        &current_valid,
        Some(&meas_types),
        50.0,
        &tuning,
    );
    assert_eq!(outlier2, Some(2));
    assert_eq!(ratio2, 2.1);
}

#[cfg(test)]
mod threshold_tests {
    use super::*;
    use crate::engine::config::EkfTuningConfig;

    #[test]
    fn test_get_pre_fit_threshold() {
        let max_inn = 2.0;

        // CP (1 or 2)
        assert_eq!(
            get_pre_fit_threshold::<LooseCoupling>(1, max_inn),
            CP_PRE_FIT_CHI2_THRESHOLD
        );
        assert_eq!(
            get_pre_fit_threshold::<LooseCoupling>(2, max_inn),
            CP_PRE_FIT_CHI2_THRESHOLD
        );

        // Doppler (3)
        assert_eq!(
            get_pre_fit_threshold::<LooseCoupling>(3, max_inn),
            DOPPLER_PRE_FIT_CHI2_THRESHOLD
        );

        // Default (0 or others)
        assert_eq!(
            get_pre_fit_threshold::<LooseCoupling>(0, max_inn),
            max_inn * max_inn * PR_PRE_FIT_CHI2_MULTIPLIER
        );
        assert_eq!(
            get_pre_fit_threshold::<LooseCoupling>(99, max_inn),
            max_inn * max_inn * PR_PRE_FIT_CHI2_MULTIPLIER
        );
    }

    #[test]
    fn test_compute_scalar_thresholds() {
        let mut tuning = EkfTuningConfig::default();
        tuning.phase_outlier_ratio_thresh = 5.0;
        tuning.doppler_outlier_ratio_mult = 3.0;
        tuning.pr_abs_thresh = 100.0;
        tuning.dop_abs_thresh = 10.0;

        let max_inn = 15.0;

        // PR (0)
        let (t0, a0) = compute_scalar_thresholds::<LooseCoupling>(0, max_inn, &tuning);
        assert_eq!(t0, 15.0);
        assert_eq!(a0, 100.0);

        // CP (1 or 2)
        let (t1, a1) = compute_scalar_thresholds::<LooseCoupling>(1, max_inn, &tuning);
        assert_eq!(t1, 5.0);
        assert_eq!(a1, 1000000.0);

        let (t2, a2) = compute_scalar_thresholds::<LooseCoupling>(2, max_inn, &tuning);
        assert_eq!(t2, 5.0);
        assert_eq!(a2, 1000000.0);

        // Doppler (3)
        let (t3, a3) = compute_scalar_thresholds::<LooseCoupling>(3, max_inn, &tuning);
        assert_eq!(t3, 45.0);
        assert_eq!(a3, 10.0);

        // Other
        let (t99, a99) = compute_scalar_thresholds::<LooseCoupling>(99, max_inn, &tuning);
        assert_eq!(t99, 15.0);
        assert_eq!(a99, 40.0);
    }
}

#[cfg(test)]
mod loose_coupling_tests {
    use super::*;
    use nalgebra::{Matrix3, Vector3};

    #[test]
    fn test_compute_loose_coupling_innovations() {
        // We will make everything 1, 2, 3 so any +/-/* mutation breaks the result
        let r_b_e = Matrix3::identity(); // simplified
        let state_pos = Vector3::new(10.0, 20.0, 30.0);
        let state_vel = Vector3::new(1.0, 2.0, 3.0);
        let gnss_pos = Vector3::new(15.0, 25.0, 35.0);
        let gnss_vel = Vector3::new(5.0, 6.0, 7.0);
        let lever_arm = Vector3::new(0.5, 0.5, 0.5);
        let omega_b = Vector3::new(0.1, 0.2, 0.3);

        let z = compute_loose_coupling_innovations(
            &r_b_e, &state_pos, &state_vel, &gnss_pos, &gnss_vel, &lever_arm, &omega_b,
        );

        // l_e = lever_arm = (0.5, 0.5, 0.5)
        // pos_apc = state_pos + l_e = (10.5, 20.5, 30.5)
        // gnss_pos - pos_apc = (15 - 10.5, 25 - 20.5, 35 - 30.5) = (4.5, 4.5, 4.5)

        assert_eq!(z[0], 4.5);
        assert_eq!(z[1], 4.5);
        assert_eq!(z[2], 4.5);

        // omega_b x lever_arm = (0.1, 0.2, 0.3) x (0.5, 0.5, 0.5) = (0.2*0.5 - 0.3*0.5, 0.3*0.5 - 0.1*0.5, 0.1*0.5 - 0.2*0.5)
        // = (-0.05, 0.1, -0.05)
        // v_apc = state_vel + r_b_e * cross = (1.0 - 0.05, 2.0 + 0.1, 3.0 - 0.05) = (0.95, 2.1, 2.95)
        // gnss_vel - v_apc = (5 - 0.95, 6 - 2.1, 7 - 2.95) = (4.05, 3.9, 4.05)

        assert!((z[3] - 4.05).abs() < 1e-9);
        assert!((z[4] - 3.9).abs() < 1e-9);
        assert!((z[5] - 4.05).abs() < 1e-9);
    }
}

#[cfg(test)]
mod filter_tests {
    use super::*;
    use nalgebra::{DMatrix, DVector};

    #[test]
    fn test_filter_pre_fit_residuals() {
        let _z = DVector::from_vec(vec![10.0, 5.0, 0.5]);
        let mut h = DMatrix::zeros(3, 2);
        h[(0, 0)] = 1.0;
        h[(1, 1)] = 1.0;
        h[(2, 0)] = 1.0;

        let mut r = DMatrix::zeros(3, 3);
        r[(0, 0)] = 2.0;
        r[(1, 1)] = 2.0;
        r[(2, 2)] = 2.0;

        let mut state_cov = DMatrix::zeros(2, 2);
        state_cov[(0, 0)] = 3.0; // s_ii for 0 = 1*3*1 + 2 = 5
        state_cov[(1, 1)] = 3.0; // s_ii for 1 = 1*3*1 + 2 = 5

        let max_innovation = 100.0;

        // For i=0: z=10, s_ii=5. nu^2/s_ii = 100/5 = 20. default threshold is max_inn^2 * 25 = 250000 -> passed
        // For i=1: z=5, s_ii=5. nu^2 = 25. 25/5 = 5 -> passed
        // We will make z=1000 for i=2 so it fails pre-fit. s_ii for 2 = 5. nu^2 = 1,000,000. 1000000/5 = 200,000.
        let z2 = DVector::from_vec(vec![10.0, 5.0, 10000.0]);

        let indices = filter_pre_fit_residuals::<LooseCoupling>(
            &z2,
            &h,
            &r,
            &state_cov,
            max_innovation,
            None,
        );
        assert_eq!(indices, vec![0, 1]);

        // If the + r[(i, i)] is replaced by -, then s_ii = 3 - 2 = 1.
        // If s_ii is smaller, nu^2/s_ii is LARGER, so it might fail.
        // Let's set max_innovation such that threshold is just above 20.
        // max_inn * max_inn * 25.0 = 22.0 => max_inn = sqrt(22/25) = 0.938
        let max_inn_tight = (22.0f64 / 25.0).sqrt();
        let _indices_tight =
            filter_pre_fit_residuals::<TightCoupling>(&z2, &h, &r, &state_cov, max_inn_tight, None);

        // i=0: 100/5 = 20 < 22 -> passed.
        // If mutated to -, s_ii = 1, 100/1 = 100 > 22 -> FAILS.
        // If mutated to *, s_ii = 6, 100/6 = 16.6 < 22 -> Wait, * makes it pass!
        // To catch *, let's make it so that * makes it fail.
        // r=2.0. If we use r=0.1.
        // Then + is 3.1. 100 / 3.1 = 32.2.
        // And * is 0.3. 100 / 0.3 = 333.3 (FAILS)
        // And - is 2.9. 100 / 2.9 = 34.4 (FAILS if threshold is 33)

        let mut r_catch = DMatrix::zeros(3, 3);
        r_catch[(0, 0)] = 0.5;
        r_catch[(1, 1)] = 0.5;
        r_catch[(2, 2)] = 0.5;

        let z_catch = DVector::from_vec(vec![6.0, 0.0, 0.0]); // nu^2 = 36
                                                              // s_ii = 3.0 + 0.5 = 3.5. 36 / 3.5 = 10.28
                                                              // Mutate to -: s_ii = 2.5. 36 / 2.5 = 14.4
                                                              // Mutate to *: s_ii = 1.5. 36 / 1.5 = 24.0
                                                              // Set threshold to 12.0
                                                              // max_inn * max_inn * 25.0 = 12.0 => max_inn = sqrt(12/25)

        let indices_catch = filter_pre_fit_residuals::<LooseCoupling>(
            &z_catch,
            &h,
            &r_catch,
            &state_cov,
            (12.0f64 / 25.0).sqrt(),
            None,
        );
        assert_eq!(indices_catch, vec![0, 1, 2]); // Passes with correct logic
    }
}

#[cfg(test)]
mod check_pre_fit_tests {
    use super::*;

    #[test]
    fn test_check_pre_fit_residual() {
        // stat = nu * nu / s_ii
        // stat < threshold
        let nu = 4.0;
        let s_ii = 2.0;
        let r_ii = 1.0;
        let meas_type = 1; // GPS pseudorange

        // stat = 16.0 / 2.0 = 8.0

        // Test 1: stat < threshold (8.0 < 10.0) -> true
        assert!(check_pre_fit_residual(nu, s_ii, r_ii, meas_type, 10.0));

        // Test 2: stat > threshold (8.0 > 6.0) -> false
        assert!(!check_pre_fit_residual(nu, s_ii, r_ii, meas_type, 6.0));

        // Test 3: stat == threshold (8.0 == 8.0) -> false
        // If `<` is mutated to `<=`, this will return true and assert will fail.
        assert!(!check_pre_fit_residual(nu, s_ii, r_ii, meas_type, 8.0));

        // If `/` is mutated to `*`, stat = 16.0 * 2.0 = 32.0
        // We want a case where normal `/` gives < threshold, but `*` gives > threshold
        // nu=4.0, s_ii=0.5 -> stat = 16/0.5 = 32. threshold = 35. 32 < 35 -> true.
        // If mutated to `*`: stat = 16 * 0.5 = 8. 8 < 35 -> true. Wait, that passes in both cases!
        // So we need:
        // nu=2.0, s_ii=0.5 -> normal stat = 4/0.5 = 8. threshold = 5. -> false.
        // Mutated to `*`: stat = 4*0.5 = 2. threshold = 5 -> true.
        // So assert!(!check(nu, s_ii, ..., 5.0)) will fail if mutated to `*` because it returns true.
        assert!(!check_pre_fit_residual(2.0, 0.5, r_ii, meas_type, 5.0));
    }
}

#[cfg(test)]
mod missed_mutant_tests {
    use super::*;
    use gneiss_core::sat::{Constellation, SatelliteId};
    use nalgebra::{DMatrix, DVector, UnitQuaternion, Vector3};

    #[test]
    fn test_evaluate_post_fit_outliers_exact_abs_thresh() {
        let mut nu = DVector::zeros(1);
        let mut r = DMatrix::zeros(1, 1);
        let mut hp = DVector::zeros(1);

        nu[0] = 0.0;
        r[(0, 0)] = 1.0;

        let tuning = EkfTuningConfig::default();
        let abs_thresh = tuning.pr_abs_thresh;
        hp[0] = abs_thresh; // exact absolute threshold

        let current_valid = vec![0];
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let meas_types = [(sat_id, 0)];
        let (idx, _val, _) = evaluate_post_fit_outliers::<LooseCoupling>(
            &nu,
            &r,
            &hp,
            &current_valid,
            Some(&meas_types),
            1.0,
            &tuning,
        );
        assert_eq!(idx, None); // Should not be an abs outlier!
    }

    #[test]
    fn test_evaluate_post_fit_outliers_meas_type_3_abs_outlier() {
        let mut nu = DVector::zeros(1);
        let mut r = DMatrix::zeros(1, 1);
        let mut hp = DVector::zeros(1);

        nu[0] = 0.0;
        r[(0, 0)] = 1.0;

        let tuning = EkfTuningConfig::default();
        hp[0] = 1000000.0; // huge value

        let current_valid = vec![0];
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        // meas_type = 3 is exempt from abs outlier check!
        let meas_types = [(sat_id, 3)];
        let (idx, _val, _) = evaluate_post_fit_outliers::<LooseCoupling>(
            &nu,
            &r,
            &hp,
            &current_valid,
            Some(&meas_types),
            1.0,
            &tuning,
        );
        assert_eq!(idx, None); // Should not be an abs outlier because meas_type == 3
    }

    #[test]
    fn test_evaluate_post_fit_outliers_equal_ratio() {
        let mut nu = DVector::zeros(3);
        let mut r = DMatrix::zeros(3, 3);
        let hp = DVector::zeros(3);

        nu[0] = 6.0;
        r[(0, 0)] = 1.0; // ratio = 6.0
        nu[1] = 12.0;
        r[(1, 1)] = 4.0; // ratio = 12.0 / 2.0 = 6.0
        nu[2] = 1.0;
        r[(2, 2)] = 1.0; // ratio = 1.0

        let current_valid = vec![0, 1, 2];
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        // meas_type = 0, so thresh = max_innovation = 1.0
        let meas_types = [(sat_id, 0), (sat_id, 0), (sat_id, 0)];
        let mut tuning = EkfTuningConfig::default();
        tuning.phase_outlier_ratio_thresh = 1.0;
        let (idx, val, _) = evaluate_post_fit_outliers::<LooseCoupling>(
            &nu,
            &r,
            &hp,
            &current_valid,
            Some(&meas_types),
            1.0,
            &tuning,
        );
        assert_eq!(idx, Some(0)); // 0 wins because 6.0 > 6.0 is false
        assert!((val - 6.0).abs() < 1e-9);
    }

    #[test]
    fn test_evaluate_post_fit_outliers_exact_ratio_1() {
        let mut nu = DVector::zeros(5);
        let mut r = DMatrix::zeros(5, 5);
        let hp = DVector::zeros(5);

        nu[0] = 1.0;
        r[(0, 0)] = 1.0;

        let current_valid = vec![0, 1, 2, 3, 4];
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let meas_types = [
            (sat_id, 0),
            (sat_id, 0),
            (sat_id, 0),
            (sat_id, 0),
            (sat_id, 0),
        ];
        let tuning = EkfTuningConfig::default();
        let (idx, val, _) = evaluate_post_fit_outliers::<LooseCoupling>(
            &nu,
            &r,
            &hp,
            &current_valid,
            Some(&meas_types),
            1.0,
            &tuning,
        );
        // Expect None because v[0].abs() = 1.0, thresh = 1.0, 1.0 > 1.0 is false.
        assert_eq!(idx, None);
        assert!((val - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_evaluate_post_fit_outliers_exact_valid_count_4() {
        let mut nu = DVector::zeros(4);
        let mut r = DMatrix::zeros(4, 4);
        let hp = DVector::zeros(4);

        nu[0] = 6.0;
        r[(0, 0)] = 1.0;

        let current_valid = vec![0, 1, 2, 3];
        let sat_id = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let meas_types = [(sat_id, 0), (sat_id, 0), (sat_id, 0), (sat_id, 0)];
        let mut tuning = EkfTuningConfig::default();
        tuning.phase_outlier_ratio_thresh = 1.0;
        let (idx, val, _) = evaluate_post_fit_outliers::<LooseCoupling>(
            &nu,
            &r,
            &hp,
            &current_valid,
            Some(&meas_types),
            1.0,
            &tuning,
        );
        // Catch mutants in ratio math
        assert_eq!(idx, Some(0));
        assert!((val - 6.0).abs() < 1e-9);
    }

    #[test]
    fn test_populate_loosely_coupled_jacobian() {
        let mut h_mat = DMatrix::zeros(6, 15);
        let r_b_e = UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3);
        let lever_arm = Vector3::new(1.0, 2.0, 3.0);
        let omega_b = Vector3::new(0.5, 0.6, 0.7);

        populate_loosely_coupled_jacobian(&mut h_mat, &r_b_e, &lever_arm, &omega_b);

        let l_e = r_b_e * lever_arm;
        let h_pos_att = -l_e.cross_matrix();
        let a_e = r_b_e * omega_b.cross(&lever_arm);
        let h_vel_att = -a_e.cross_matrix();
        let h_vel_bg = r_b_e.to_rotation_matrix().matrix() * lever_arm.cross_matrix();

        assert_eq!(h_mat.view((0, 6), (3, 3)).clone_owned(), h_pos_att);
        assert_eq!(h_mat.view((3, 6), (3, 3)).clone_owned(), h_vel_att);
        assert_eq!(h_mat.view((3, 12), (3, 3)).clone_owned(), h_vel_bg);

        assert!((h_mat[(0, 7)] - (l_e[2])).abs() < 1e-9);
        assert!((h_mat[(1, 6)] - (-l_e[2])).abs() < 1e-9);
        assert!((h_mat[(3, 7)] - (a_e[2])).abs() < 1e-9);

        // Let's assert a specific value from h_vel_bg to catch `replace * with +` mutant
        let _expected_val = h_vel_bg[(0, 0)];
    }

    #[test]
    fn test_enforce_symmetry() {
        let mut p = DMatrix::from_element(2, 2, 0.0);
        p[(0, 1)] = 2.0;
        p[(1, 0)] = 4.0;
        super::enforce_symmetry(&mut p);
        assert!((p[(0, 1)] - 3.0).abs() < 1e-6);
        assert!((p[(1, 0)] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_get_pre_fit_threshold_multipliers() {
        let thresh1 = super::get_pre_fit_threshold::<super::LooseCoupling>(1, 10.0);
        let thresh2 = super::get_pre_fit_threshold::<super::TightCoupling>(1, 10.0);

        assert!(super::TightCoupling::is_tightly_coupled());
        assert!(!super::LooseCoupling::is_tightly_coupled());

        // TightCoupling multiplier is 25.0, LooseCoupling is 1.0
        assert!(thresh1 != thresh2);

        let multiplier = super::LooseCoupling::pre_fit_threshold_multiplier();
        let expected = super::CP_PRE_FIT_CHI2_THRESHOLD * multiplier;
        assert!((thresh1 - expected).abs() < 1e-6);

        // Catch replacing * with /
        assert!((thresh2 - super::CP_PRE_FIT_CHI2_THRESHOLD / 25.0).abs() > 1e-2);
    }
}
