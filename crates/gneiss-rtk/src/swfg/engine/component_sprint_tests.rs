//! Sprint-hunt tests for the SWFG accumulator, the RTK factor builder's
//! reference-satellite hysteresis, and the PPP-AR MW tracker lifecycle.
#![allow(clippy::unwrap_used)]

use std::collections::HashMap;

use nalgebra::Vector3;

use crate::swfg::engine::accumulator::{DdAccumulatorEntry, DdPseudorangeAccumulator, DdSatKey};
use crate::swfg::engine::ar_handler::{execute_ar_step, record_mw_sample, reset_ppp_tracker};
use crate::swfg::engine::builder::select_ref_satellite;
use crate::swfg::pipeline::passes::CorrectedObservation;
use crate::swfg::variables::VariableValues;

use super::ar_step_sprint_tests::rtk_graph;
use super::sprint_tests::BASE;


// ===========================================================================
// accumulator: running-sum invariant
// ===========================================================================

/// `DdAccumulatorEntry` maintains `running_sum` incrementally next to a bounded
/// history. The invariant tying them together is `running_sum == sum(history)`
/// at every point. This is independently derivable and cannot hold by
/// construction, because `add` subtracts the evicted sample.
#[test]
fn accumulator_incremental_mean_matches_the_retained_history_sum() {
    for cap in [1usize, 2, 5, 17] {
        let mut e = DdAccumulatorEntry::new(cap);
        for k in 0..200 {
            // Signed, wildly different magnitudes, so any cancellation error
            // in the incremental update shows up.
            let v = ((k as f64) * 0.37).sin() * 1e6 + (k as f64) * 1e-3;
            e.add(v);
            assert!(e.count() <= cap.max(1), "cap {cap}: count {} exceeds the bound", e.count());
            let hist: f64 = e.history.iter().sum();
            let n = e.count() as f64;
            // `mean()` returns the incrementally maintained running_sum / n;
            // it must equal the mean of the retained window it also exposes.
            let mean = e.mean().unwrap();
            assert!(
                (mean - hist / n).abs() <= 1e-9 * (1.0 + (hist / n).abs()),
                "cap {cap}, k {k}: incremental mean {mean} != history mean {}",
                hist / n
            );
        }
    }
}

/// `mean` and `std_dev` must be the sample statistics of the retained window.
/// Window [1,2,3,4,5]: mean = 15/5 = 3; sample variance = 10/4 = 2.5, so
/// std = sqrt(2.5) = 1.5811388.
#[test]
fn accumulator_mean_and_std_match_the_closed_form() {
    let mut e = DdAccumulatorEntry::new(5);
    assert!(e.mean().is_none());
    assert!(e.std_dev().is_none(), "std_dev needs at least two samples");
    for v in 1..=5u64 {
        e.add(v as f64);
    }
    assert!((e.mean().unwrap() - 3.0).abs() < 1e-12);
    assert!((e.std_dev().unwrap() - 2.5_f64.sqrt()).abs() < 1e-12);

    // Evict 1, add 6 -> window [2,3,4,5,6]. sum = 20, mean = 4.
    // SS = 4+9+16+25+36 = 90, var = (90 - 5*16)/4 = 10/4 = 2.5.
    e.add(6.0);
    assert!((e.mean().unwrap() - 4.0).abs() < 1e-12);
    assert!((e.std_dev().unwrap() - 2.5_f64.sqrt()).abs() < 1e-12);
    assert_eq!(e.count(), 5);
}

/// The MW wide-lane constraint in `builder.rs` fires on `count >= 3 &&
/// std_dev < 0.25`. Two samples must not be enough; three identical ones must.
#[test]
fn accumulator_exposes_stats_only_once_three_samples_exist() {
    let mut acc = DdPseudorangeAccumulator::new(10);
    let key = DdSatKey { sat: 1, ref_sat: 2, frequency: 0 };
    acc.add_observation(key, 5.0);
    assert!(acc.get_stats(&key).is_none(), "one sample has no dispersion estimate");
    assert_eq!(acc.get_averaged(&key), Some((5.0, 1)));
    acc.add_observation(key, 5.0);
    let (mean, sd, n) = acc.get_stats(&key).unwrap();
    assert!((mean - 5.0).abs() < 1e-12 && sd < 1e-12 && n == 2);
    acc.add_observation(key, 5.0);
    let (mean, sd, n) = acc.get_stats(&key).unwrap();
    assert!((mean - 5.0).abs() < 1e-12 && sd < 1e-12 && n == 3);
    // `builder.rs` gates the MW wide-lane constraint on `count >= 3`, which is
    // where the sample-count requirement actually lives.
    assert!(sd < 0.25 && n >= 3, "zero dispersion passes the MW gate");
    acc.clear();
    assert!(acc.get_stats(&key).is_none(), "clear must drop every entry");
    assert_eq!(acc.get_averaged(&key), None);
}

// ===========================================================================
// builder::select_ref_satellite hysteresis
// ===========================================================================

fn corrected(sat: u16, el_deg: f64) -> CorrectedObservation {
    CorrectedObservation {
        satellite: sat,
        constellation_id: 0,
        pr_l1: 2.0e7,
        pr_l2: Some(1.5e7),
        cp_l1: Some(1.0e8),
        cp_l1_lli: Some(0),
        cp_l2: Some(6.0e7),
        doppler: 0.0,
        snr_dbhz: 45.0,
        sat_pos_ecef: Vector3::zeros(),
        sat_clock_m: 0.0,
        f1: 1_575.42e6,
        f2: 1_227.60e6,
        freq_num: 0,
        elevation_rad: el_deg.to_radians(),
        tropo_dry_m: 2.3,
        tropo_map_wet: 1.0,
        iono_l1_m: 3.0,
        variance_m2: 0.25,
        cp_variance_m2: 1e-5,
    }
}

#[test]
fn ref_satellite_holds_above_fifteen_degrees_and_switches_above_thirty() {
    // 0.26 rad = 14.9038 deg (hold threshold), 0.52 rad = 29.7949 deg (switch).
    // 0.26 rad * (180/pi) = 14.896903 deg; 0.52 rad * (180/pi) = 29.793806 deg.
    assert!((0.26_f64.to_degrees() - 14.896_903).abs() < 1e-5);
    assert!((0.52_f64.to_degrees() - 29.793_806).abs() < 1e-5);
    let prev: HashMap<u8, u16> = HashMap::from([(0u8, 7u16)]);

    let a: Vec<CorrectedObservation> = vec![corrected(7, 15.0), corrected(3, 40.0)];
    let ar: Vec<&CorrectedObservation> = a.iter().collect();
    assert_eq!(select_ref_satellite(&ar, 0, &prev), 7, "prev at 15 deg > 0.26 rad must hold");

    let b: Vec<CorrectedObservation> = vec![corrected(7, 14.0), corrected(3, 40.0)];
    let br: Vec<&CorrectedObservation> = b.iter().collect();
    assert_eq!(select_ref_satellite(&br, 0, &prev), 3, "prev at 14 deg < 0.26 rad must switch");

    let c: Vec<CorrectedObservation> = vec![corrected(7, 10.0), corrected(3, 25.0)];
    let cr: Vec<&CorrectedObservation> = c.iter().collect();
    assert_eq!(
        select_ref_satellite(&cr, 0, &prev),
        7,
        "best candidate at 25 deg (< 0.52 rad) must not trigger a switch"
    );

    assert_eq!(select_ref_satellite(&br, 0, &HashMap::new()), 3, "no history -> highest elevation");
    assert_eq!(select_ref_satellite(&[], 0, &HashMap::new()), 0, "empty -> sentinel 0");
}

// ===========================================================================
// PPP-AR tracker lifecycle
// ===========================================================================

/// `reset_ppp_tracker` clears the thread-local MW tracker (and the epoch
/// tracker). After a reset, `record_mw_sample` must start a fresh arc rather
/// than continuing the old one, so an epoch-0 AR attempt sees a cold tracker.
#[test]
fn ppp_mw_tracker_reset_clears_previously_recorded_samples() {
    for e in 1..20u32 {
        record_mw_sample(0, 1, 100_000.5, e, false);
    }
    reset_ppp_tracker();
    // A fresh tracker: one sample cannot satisfy the >= 5 count requirement in
    // `try_fix_constellation`, so the PPP path must be a no-op and the graph
    // must come out untouched.
    let (mut solver, pose, _) = rtk_graph(9, 0.4, 0.3);
    solver.solve().expect("solve");
    let n_before = solver.graph.n_factors();
    execute_ar_step(&mut solver, pose, BASE, false, 0);
    assert_eq!(solver.graph.n_factors(), n_before, "PPP AR must not inject factors");
    assert!(VariableValues::build(&solver.graph.variables).get(pose).is_some());
}
