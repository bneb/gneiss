//! Tests for the PAR candidate selection and subset machinery.
//!
//! A child module of `ar_subsets` so the private ranking helpers are reachable.

use gneiss_core::time::GpsTime;
use nalgebra::Vector3;

use super::*;

const RX: Vector3<f64> = Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067);

fn key(sat: u16, ref_sat: u16, cons: u8, band: u8) -> DoubleDiffKey {
    DoubleDiffKey { constellation_id: cons, sat, ref_sat, freq_band: band }
}

fn state_with(keys: &[DoubleDiffKey], vals: &[f64], vars: &[f64]) -> RtkState {
    let mut st = RtkState::new(RX, GpsTime::new(2200, 300_000.0));
    for ((k, v), var) in keys.iter().zip(vals).zip(vars) {
        st.ensure_ambiguity(*k, *v, *var);
    }
    st
}

/// DD measurement whose satellite sits at `el_rad` elevation and whose noise
/// fields are set directly, so the metadata mapping can be checked by hand.
fn dd_for(k: DoubleDiffKey, dir: Vector3<f64>, pr_var: f64, cp_var: f64) -> DoubleDiffMeasurement {
    let up = RX.normalize();
    let east = Vector3::new(-up.z, 0.0, up.x).normalize();
    DoubleDiffMeasurement {
        key: k,
        dd_pr_m: 0.0,
        dd_cp_cycles: None,
        sat_pos: RX + (up * dir.x + east * dir.y).normalize() * 2.4e7,
        ref_pos: RX + up * 2.4e7,
        base_pos: RX,
        lambda: 0.19,
        pr_var_m2: pr_var,
        cp_var_cycles2: cp_var,
        pr_ref_var_m2: 0.02,
        cp_ref_var_cycles2: 5.0e-5,
        dm_wet_rov: 0.0,
        dgrad_n_rov: 0.0,
        dgrad_e_rov: 0.0,
        tide_dd_m: 0.0,
        dd_pcv_m: 0.0,
    }
}

// ---------------------------------------------------------------------------
// compute_metadata_from_dd / is_severe_nlos
// ---------------------------------------------------------------------------

#[test]
fn metadata_maps_declared_noise_onto_the_documented_snr_estimate() {
    // snr_est = clamp(45 - 10*log10(cp_var/5e-4), 20, 50)
    //   cp_var = 5e-4  -> log10(1) = 0      -> 45.0
    //   cp_var = 5e-5  -> log10(0.1) = -1  -> 55 -> clamped to 50
    //   cp_var = 5e-3  -> log10(10) = 1    -> 35.0
    let k = key(2, 1, 0, 1);
    let st = state_with(&[k], &[1.0], &[0.04]);
    let dir = Vector3::new(0.8, 0.6, 0.0);
    for (cp_var, want) in [(5.0e-4, 45.0), (5.0e-5, 50.0), (5.0e-3, 35.0)] {
        let dd = vec![dd_for(k, dir, 0.04, cp_var)];
        let meta = compute_metadata_from_dd(&st, &dd);
        assert_eq!(meta.len(), 1);
        assert!((meta[0].snr_dbhz - want).abs() < 1e-9, "cp_var {cp_var}: got {}", meta[0].snr_dbhz);
    }
}

#[test]
fn metadata_clamps_the_cmc_sigma_to_the_documented_band() {
    // cmc_sigma = min(sqrt(max(pr_var, 1e-4)), 10)
    //   pr_var = 0.04 -> 0.2 ;  pr_var = 0 -> sqrt(1e-4) = 0.01 ;
    //   pr_var = 1000 -> 31.6 -> clamped to 10
    let k = key(2, 1, 0, 1);
    let st = state_with(&[k], &[1.0], &[0.04]);
    let dd = vec![dd_for(k, Vector3::new(1.0, 0.0, 0.0), 0.0, 5.0e-4)];
    assert!((compute_metadata_from_dd(&st, &dd)[0].cmc_sigma - 0.01).abs() < 1e-12);
    let dd2 = vec![dd_for(k, Vector3::new(1.0, 0.0, 0.0), 1000.0, 5.0e-4)];
    assert!((compute_metadata_from_dd(&st, &dd2)[0].cmc_sigma - 10.0).abs() < 1e-12);
}

#[test]
fn a_pair_without_a_measurement_gets_the_neutral_default() {
    let k = key(2, 1, 0, 1);
    let st = state_with(&[k], &[1.0], &[0.04]);
    let meta = compute_metadata_from_dd(&st, &[]);
    assert_eq!(meta.len(), 1);
    assert_eq!(meta[0], AmbiguityMetadata::default());
}

#[test]
fn severe_nlos_needs_both_a_low_elevation_and_a_degraded_signal() {
    let el_ok = 30.0_f64.to_radians();
    let el_bad = 10.0_f64.to_radians();
    // High elevation is never "severe", however bad the SNR.
    assert!(!is_severe_nlos(Some(&AmbiguityMetadata::new(el_ok, 5.0, 1.0, 9.0))));
    // Low elevation with a weak signal, or with heavy code multipath, is.
    assert!(is_severe_nlos(Some(&AmbiguityMetadata::new(el_bad, 27.9, 30.0, 0.1))));
    assert!(is_severe_nlos(Some(&AmbiguityMetadata::new(el_bad, 45.0, 30.0, 3.1))));
    // Low elevation but a healthy signal is not.
    assert!(!is_severe_nlos(Some(&AmbiguityMetadata::new(el_bad, 28.0, 30.0, 3.0))));
    // No metadata at all can never be flagged.
    assert!(!is_severe_nlos(None));
}

// ---------------------------------------------------------------------------
// kinematic candidate pool ladder
// ---------------------------------------------------------------------------

#[test]
fn kinematic_pool_keeps_only_variance_bounded_non_geo_non_nlos_pairs() {
    let keys = [
        key(2, 1, 0, 1), key(3, 1, 0, 1), key(6, 1, 0, 1), key(7, 1, 0, 1), // clean
        key(9, 1, 0, 1),    // var 4.0 -> rejected (q_ii > 1.0)
        key(4, 1, 3, 1),    // BeiDou GEO PRN 4 -> rejected
        key(11, 1, 0, 1),   // severe NLOS -> rejected
    ];
    let vals: Vec<f64> = (1..=7).map(|i| i as f64).collect();
    let vars = [0.25, 0.25, 0.25, 0.25, 4.0, 0.25, 0.25];
    let st = state_with(&keys, &vals, &vars);
    let q = DMatrix::from_diagonal(&DVector::from_vec(vars.to_vec()));
    let meta: Vec<AmbiguityMetadata> = (0..7).map(|i| {
        let bad = i == 6;
        AmbiguityMetadata::new(
            (if bad { 5.0_f64 } else { 30.0_f64 }).to_radians(),
            if bad { 20.0 } else { 45.0 },
            if bad { 1.0 } else { 30.0 },
            if bad { 8.0 } else { 0.2 },
        )
    }).collect();
    let pool = filter_kinematic_pool(&st, &q, 7, 4, Some(&meta));
    assert_eq!(pool, vec![0, 1, 2, 3], "only the clean, bounded, non-GEO pairs survive");
}

#[test]
fn kinematic_pool_relaxes_the_nlos_filter_before_it_relaxes_the_variance_filter() {
    // Every candidate is NLOS. With min_k = 2 the pool can still be filled
    // once the NLOS filter is dropped, and the variance filter is never
    // reached -- so the high-variance pair must come back.
    let keys = [key(2, 1, 0, 1), key(3, 1, 0, 1)];
    let st = state_with(&keys, &[1.0, 2.0], &[0.25, 4.0]);
    let a = DVector::from_vec(vec![1.0, 2.0]);
    let q = DMatrix::from_diagonal(&DVector::from_vec(vec![0.25, 4.0]));
    let meta = vec![
        AmbiguityMetadata::new(1.0_f64.to_radians(), 5.0, 1.0, 9.0),
        AmbiguityMetadata::new(1.0_f64.to_radians(), 5.0, 1.0, 9.0),
    ];
    let pool = filter_kinematic_pool(&st, &q, 2, 2, Some(&meta));
    assert_eq!(pool, vec![0, 1], "NLOS filter is the first rung to be dropped");
    let _ = a;
}

#[test]
fn kinematic_pool_falls_all_the_way_back_to_every_candidate() {
    // Only two candidates, both BeiDou GEO with huge variance. The GEO filter
    // is dropped next, and with nothing left to filter the pool is all of them.
    let keys = [key(3, 1, 3, 1), key(60, 1, 3, 1)];
    let st = state_with(&keys, &[1.0, 2.0], &[9.0, 9.0]);
    let q = DMatrix::from_diagonal(&DVector::from_vec(vec![9.0, 9.0]));
    let pool = filter_kinematic_pool(&st, &q, 2, 2, None);
    assert_eq!(pool, vec![0, 1]);
}

// ---------------------------------------------------------------------------
// select_par_candidates_*
// ---------------------------------------------------------------------------

#[test]
fn static_ranking_is_by_ascending_variance_and_caps_the_pool() {
    let n = 20;
    let keys: Vec<DoubleDiffKey> = (0..n).map(|i| key(2 + i as u16, 1, 0, 1)).collect();
    let vals: Vec<f64> = (0..n).map(|i| i as f64).collect();
    let vars: Vec<f64> = (0..n).map(|i| 1.0 + i as f64).collect();
    let st = state_with(&keys, &vals, &vars);
    let a = DVector::from_vec(vals.clone());
    let q = DMatrix::from_diagonal(&DVector::from_vec(vars));
    let (idx, m) = select_par_candidates_with_metadata(&st, &a, &q, 3, false, None);
    // The full ordering is returned; `m` is the caller's pool size, capped at
    // MAX_SUBSET_SIZE and always one short of the candidate count.
    assert_eq!(m, (n - 1).min(MAX_SUBSET_SIZE));
    assert_eq!(idx.len(), n);
    // Tightest covariance first: variances 1.0, 1.1, ... in index order.
    assert_eq!(idx[0], 0);
    assert_eq!(idx[1], 1);
    assert!(idx.windows(2).all(|w| w[0] < w[1]));
}

#[test]
fn dd_supplied_metadata_is_used_when_no_cqm_slice_is_given() {
    let keys = [key(2, 1, 0, 1), key(3, 1, 0, 1)];
    let st = state_with(&keys, &[10.001, 20.001], &[0.04, 0.04]);
    let a = DVector::from_vec(vec![10.001, 20.001]);
    let q = DMatrix::from_diagonal(&DVector::from_vec(vec![0.04, 0.04]));
    let dd = vec![
        dd_for(keys[0], Vector3::new(1.0, 0.0, 0.0), 0.04, 5.0e-4),
        dd_for(keys[1], Vector3::new(1.0, 0.0, 0.0), 9.0, 5.0e-3),
    ];
    let with_dd = select_par_candidates_with_dd(&st, &a, &q, 2, true, Some(&dd));
    let without = select_par_candidates_with_dd(&st, &a, &q, 2, true, None);
    assert_eq!(with_dd.0.len(), 2);
    // The noisier pair has the lower CQM and must be ranked second.
    assert_eq!(with_dd.0[0], 0);
    assert_eq!(without.0.len(), 2);
}

#[test]
fn the_plain_selector_matches_the_metadata_free_path() {
    let keys = [key(2, 1, 0, 1), key(3, 1, 0, 1)];
    let st = state_with(&keys, &[10.001, 20.001], &[0.04, 0.09]);
    let a = DVector::from_vec(vec![10.001, 20.001]);
    let q = DMatrix::from_diagonal(&DVector::from_vec(vec![0.04, 0.09]));
    assert_eq!(
        select_par_candidates(&st, &a, &q, 2, false),
        select_par_candidates_with_metadata(&st, &a, &q, 2, false, None)
    );
}

// ---------------------------------------------------------------------------
// subset generation
// ---------------------------------------------------------------------------

#[test]
fn constellation_partitions_keep_only_the_declared_clusters() {
    let keys = [
        key(2, 1, 0, 1),   // GPS
        key(3, 1, 0, 1),   // GPS
        key(5, 1, 2, 1),   // Galileo
        key(6, 1, 3, 1),   // BeiDou
        key(7, 1, 3, 1),   // BeiDou
    ];
    let st = state_with(&keys, &[1.0, 2.0, 3.0, 4.0, 5.0], &[0.04; 5]);
    let ranked: Vec<usize> = (0..5).collect();
    let parts = partition_constellation_subsets(&st, &ranked, 3);
    // cluster [0,2] = GPS+Galileo -> 3 ; [0,3] = GPS+BeiDou -> 4 ;
    // [2,3] = Galileo+BeiDou -> 3. All three clear min_k = 3.
    assert_eq!(parts.len(), 3);
    assert_eq!(parts[0].as_slice(), &[0, 1, 2]);
    assert_eq!(parts[1].as_slice(), &[0, 1, 3, 4]);
    assert_eq!(parts[2].as_slice(), &[2, 3, 4]);
}

#[test]
fn constellation_partitions_drop_clusters_below_the_floor() {
    let keys = [key(2, 1, 0, 1), key(5, 1, 2, 1)];
    let st = state_with(&keys, &[1.0, 2.0], &[0.04, 0.04]);
    // Only GPS+Galileo has 2 members; requiring 3 leaves nothing.
    assert!(partition_constellation_subsets(&st, &[0, 1], 3).is_empty());
    assert_eq!(partition_constellation_subsets(&st, &[0, 1], 2).len(), 1);
}

#[test]
fn two_omission_subsets_need_room_to_remove_two_and_keep_min_k() {
    let idx: Vec<usize> = (1..=6).collect();
    // pool 6, min_k 5 -> each subset keeps 4 < 5, so nothing survives.
    assert!(generate_two_omission_subsets(&idx, 6, 5).is_empty());
    // min_k 4 -> every pair of omissions leaves exactly 4 members.
    let got = generate_two_omission_subsets(&idx, 6, 4);
    assert_eq!(got.len(), 6 * 5 / 2);
    assert!(got.iter().all(|s| s.len == 4));
    // The pool is additionally capped at 10 entries.
    let long: Vec<usize> = (1..=20).collect();
    assert_eq!(generate_two_omission_subsets(&long, 20, 4).len(), 10 * 9 / 2);
}

#[test]
fn one_omission_subsets_are_capped_by_the_stack_limit() {
    let long: Vec<usize> = (0..30).collect();
    let got = generate_omission_subsets(&long, 30);
    assert_eq!(got.len(), MAX_SUBSET_SIZE);
    assert!(got.iter().all(|s| s.len == MAX_SUBSET_SIZE - 1));
}

#[test]
fn dd_subset_geometry_needs_three_pairs_of_real_spread() {
    let up = RX.normalize();
    let e1 = Vector3::new(-up.z, 0.0, up.x).normalize();
    let e2 = up.cross(&e1);
    let keys: Vec<DoubleDiffKey> = (2..6).map(|i| key(i, 1, 0, 1)).collect();
    let st = state_with(&keys, &[1.0, 2.0, 3.0, 4.0], &[0.04; 4]);
    // Four satellites on a 30 deg cone at four azimuths, each sharing one
    // overhead reference: five unique positions, PDOP ~2.5.
    let dd: Vec<DoubleDiffMeasurement> = keys.iter().enumerate().map(|(i, k)| {
        let az = i as f64 * std::f64::consts::FRAC_PI_2;
        let el = 30.0_f64.to_radians();
        let d = (up * el.sin() + (e1 * az.cos() + e2 * az.sin()) * el.cos()).normalize();
        let mut m = dd_for(*k, up, 0.04, 5.0e-4);
        m.sat_pos = RX + d * 2.2e7;
        m
    }).collect();
    // Sanity: the same five positions clear the DOP guard directly.
    let uniq: Vec<Vector3<f64>> = dd.iter().flat_map(|m| [m.sat_pos, m.ref_pos]).collect();
    assert!(validate_subset_geometry(RX, &uniq, MAX_ACCEPTABLE_PDOP), "fixture geometry too poor");
    assert!(validate_dd_subset_geometry(&st, &dd, &[0, 1, 2, 3]));
    assert!(!validate_dd_subset_geometry(&st, &dd, &[0, 1]));
    // Indices past the ambiguity list are skipped, which can starve the guard.
    assert!(!validate_dd_subset_geometry(&st, &dd, &[0, 1, 99]));
}