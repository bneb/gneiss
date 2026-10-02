//! Unit tests for `post_process::backward` (pass 3).
//!
//! Companion file: `backward.rs` is over the 300-line threshold, so its tests
//! live here and are pulled in with `#[cfg(test)] #[path = ...] mod tests`.
//!
//! The preintegration expectations below use the closed-form midpoint
//! trapezoid that `ImuPreintegration::integrate` implements for
//! `dq = identity`, `ba = bg = 0`:
//!   dv_{k+1} = dv_k + a * dt
//!   dp_{k+1} = dp_k + dv_{k+1} * dt + 0.5 * a * dt^2

use super::*;
use crate::post_process::dynamics::ProcessingDynamics;
use crate::swfg::config::{ImuConfig, PppConfig, PppInsConfig, RtkConfig, RtkInsConfig, SppConfig};

fn rtk_ins(rtk: RtkConfig) -> EngineConfig {
    EngineConfig::RtkIns(RtkInsConfig {
        rtk,
        imu: ImuConfig::default(),
    })
}

fn ppp_ins(ppp: PppConfig) -> EngineConfig {
    EngineConfig::PppIns(PppInsConfig {
        ppp,
        imu: ImuConfig::default(),
    })
}

pub(super) fn epoch(tow: f64) -> EpochObs {
    EpochObs {
        time: gneiss_core::time::GpsTime::new(2000, tow),
        satellites: Vec::new(),
    }
}

pub(super) fn imu(t_us: u64, accel: Vector3<f64>, gyro: Vector3<f64>) -> ImuSample {
    ImuSample {
        accel,
        gyro,
        time_us: t_us,
    }
}

pub(super) fn seed() -> Vector3<f64> {
    Vector3::new(-3_961_907.0, 3_351_057.0, 3_694_313.0)
}

// ---------------------------------------------------------------------------
// configure_backward_engine
// ---------------------------------------------------------------------------

/// The backward pass is seeded from the end of the forward trajectory, so the
/// seed must actually reach the engine config for every profile -- not just
/// for `Rtk`, which is where it used to be wired. 10/20/30 m offsets from the
/// pre-existing value make a dropped assignment obvious.
#[test]
fn backward_seed_overrides_initial_position_for_every_profile() {
    let s = seed();
    let old = [1.0, 2.0, 3.0];

    let rtk = configure_backward_engine(
        &EngineConfig::Rtk(RtkConfig {
            initial_position: Some(old),
            ..Default::default()
        }),
        Some(s),
        ProcessingDynamics::Static,
    );
    match rtk {
        EngineConfig::Rtk(c) => assert_eq!(c.initial_position, Some([s.x, s.y, s.z])),
        _ => panic!("profile must be preserved"),
    }

    let ins = configure_backward_engine(
        &rtk_ins(RtkConfig {
            initial_position: Some(old),
            ..Default::default()
        }),
        Some(s),
        ProcessingDynamics::Static,
    );
    match ins {
        EngineConfig::RtkIns(c) => assert_eq!(c.rtk.initial_position, Some([s.x, s.y, s.z])),
        _ => panic!("profile must be preserved"),
    }

    let ppk = configure_backward_engine(
        &EngineConfig::Ppp(PppConfig {
            initial_position: Some(old),
            ..Default::default()
        }),
        Some(s),
        ProcessingDynamics::Static,
    );
    match ppk {
        EngineConfig::Ppp(c) => assert_eq!(c.initial_position, Some([s.x, s.y, s.z])),
        _ => panic!("profile must be preserved"),
    }

    let pins = configure_backward_engine(
        &ppp_ins(PppConfig {
            initial_position: Some(old),
            ..Default::default()
        }),
        Some(s),
        ProcessingDynamics::Static,
    );
    match pins {
        EngineConfig::PppIns(c) => assert_eq!(c.ppp.initial_position, Some([s.x, s.y, s.z])),
        _ => panic!("profile must be preserved"),
    }

    let spp = configure_backward_engine(
        &EngineConfig::Spp(SppConfig {
            initial_position: Some(old),
            ..Default::default()
        }),
        Some(s),
        ProcessingDynamics::Static,
    );
    match spp {
        EngineConfig::Spp(c) => assert_eq!(c.initial_position, Some([s.x, s.y, s.z])),
        _ => panic!("profile must be preserved"),
    }
}

/// `pos_arr.or(c.initial_position)` means a missing seed must leave the
/// configured value alone rather than overwrite it with `None`.
#[test]
fn backward_config_keeps_the_existing_seed_when_none_is_supplied() {
    let old = [7.0, 8.0, 9.0];
    let cfg = EngineConfig::Rtk(RtkConfig {
        initial_position: Some(old),
        ..Default::default()
    });
    let out = configure_backward_engine(&cfg, None, ProcessingDynamics::Static);
    match out {
        EngineConfig::Rtk(c) => assert_eq!(c.initial_position, Some(old)),
        _ => panic!("profile must be preserved"),
    }
}

/// The PPP profiles mirror the forward pass's kinematics: `is_kinematic` is
/// taken from the motion model, and a supplied seed pins the initial position
/// sigma to 0.15 m. The `Rtk`/`Spp` profiles have no sigma field at all, so
/// the PPP-only assignment must not be applied to them.
#[test]
fn backward_config_sets_ppp_kinematics_and_seed_sigma() {
    let s = seed();
    let dyn_ppp = configure_backward_engine(
        &EngineConfig::Ppp(Default::default()),
        Some(s),
        ProcessingDynamics::Kinematic,
    );
    match dyn_ppp {
        EngineConfig::Ppp(c) => {
            assert!(c.is_kinematic, "Kinematic dynamics must set is_kinematic");
            assert_eq!(
                c.initial_pos_sigma_m,
                Some(0.15),
                "a seed pins sigma to 0.15 m"
            );
        }
        _ => panic!("profile must be preserved"),
    }
    let stat_ppp = configure_backward_engine(
        &EngineConfig::Ppp(Default::default()),
        Some(s),
        ProcessingDynamics::Static,
    );
    match stat_ppp {
        EngineConfig::Ppp(c) => assert!(!c.is_kinematic, "Static dynamics must clear is_kinematic"),
        _ => panic!("profile must be preserved"),
    }
    // PppIns pins the same seed sigma, one level down. `is_kinematic` is
    // deliberately NOT asserted here: the `PppIns` arm of
    // `configure_backward_engine` only assigns `initial_position` and
    // `initial_pos_sigma_m`, and `run_forward_pass` likewise only touches
    // `is_kinematic` for the plain-`Ppp` variant -- so the two passes agree
    // with each other and that agreement is what a regression test should
    // pin. Whether PPP-INS *should* honour the motion model is a separate
    // design question, not something to be frozen here.
    let pins = configure_backward_engine(
        &ppp_ins(PppConfig::default()),
        Some(s),
        ProcessingDynamics::Kinematic,
    );
    match pins {
        EngineConfig::PppIns(c) => assert_eq!(c.ppp.initial_pos_sigma_m, Some(0.15)),
        _ => panic!("profile must be preserved"),
    }
}

/// Without a seed the sigma is left exactly as configured -- the
/// `if has_seed` guard must not blanket-assign 0.15.
#[test]
fn backward_config_does_not_touch_sigma_without_a_seed() {
    let cfg = EngineConfig::Ppp(PppConfig {
        initial_pos_sigma_m: Some(3.0),
        ..Default::default()
    });
    let out = configure_backward_engine(&cfg, None, ProcessingDynamics::Kinematic);
    match out {
        EngineConfig::Ppp(c) => assert_eq!(c.initial_pos_sigma_m, Some(3.0)),
        _ => panic!("profile must be preserved"),
    }
    // An explicitly-absent sigma must stay absent.
    let cfg = EngineConfig::Ppp(PppConfig {
        initial_pos_sigma_m: None,
        ..Default::default()
    });
    match configure_backward_engine(&cfg, None, ProcessingDynamics::Static) {
        EngineConfig::Ppp(c) => assert!(c.initial_pos_sigma_m.is_none()),
        _ => panic!("profile must be preserved"),
    }
}

// ---------------------------------------------------------------------------
// group_imu_by_epoch
// ---------------------------------------------------------------------------

/// Without IMU data the map is empty and every lookup misses.
#[test]
fn imu_grouping_without_data_is_empty() {
    assert!(group_imu_by_epoch(&[epoch(100.0)], None).is_empty());
}

/// Samples are consumed in order with `time_us <= cur_us`, so the 0.5 s and
/// 1.5 s samples land in the tow-1 and tow-2 buckets respectively. Epoch keys
/// are milliseconds: 1000 and 2000.
#[test]
fn imu_grouping_buckets_samples_by_epoch() {
    let samples = vec![
        imu(
            500_000,
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.1, 0.0, 0.0),
        ),
        imu(
            1_500_000,
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(0.2, 0.0, 0.0),
        ),
    ];
    let map = group_imu_by_epoch(&[epoch(1.0), epoch(2.0)], Some(&samples));
    assert_eq!(map.len(), 2, "both epochs received samples");
    assert_eq!(map[&1000].len(), 1);
    assert_eq!(map[&2000].len(), 1);
    assert_eq!(map[&1000][0].accel.x, 1.0);
    assert_eq!(map[&2000][0].accel.x, 2.0);
}

/// Epochs that receive no samples are not inserted at all, so a later lookup
/// misses rather than resolving to an empty slice.
#[test]
fn imu_grouping_skips_epochs_without_samples() {
    let samples = vec![
        imu(
            1_500_000,
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(0.2, 0.0, 0.0),
        ),
        imu(
            2_500_000,
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.3, 0.0, 0.0),
        ),
    ];
    let map = group_imu_by_epoch(&[epoch(1.0), epoch(2.0), epoch(3.0)], Some(&samples));
    assert!(!map.contains_key(&1000), "tow 1.0 s precedes every sample");
    assert_eq!(map[&2000].len(), 1);
    assert_eq!(map[&3000].len(), 1);
    assert_eq!(map[&2000][0].accel.x, 2.0);
    assert_eq!(map[&3000][0].accel.x, 3.0);
    // Samples strictly after the last epoch are never consumed.
    let late = vec![imu(9_000_000, Vector3::zeros(), Vector3::zeros())];
    assert!(group_imu_by_epoch(&[epoch(1.0)], Some(&late)).is_empty());
}

// ---------------------------------------------------------------------------
// extract_backward_imu_slice
// ---------------------------------------------------------------------------

/// A key that was never grouped, a single-sample bucket, and an all-zero-rate
/// bucket are all rejected -- the last one because a motionless window carries
/// no rotational information.
#[test]
fn backward_imu_slice_rejects_unusable_buckets() {
    let mut map = BTreeMap::new();
    assert!(
        extract_backward_imu_slice(100_000, &map).is_none(),
        "missing key"
    );
    map.insert(
        100_000,
        vec![imu(
            100_000,
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.1, 0.0, 0.0),
        )],
    );
    assert!(
        extract_backward_imu_slice(100_000, &map).is_none(),
        "one sample"
    );
    map.insert(
        200_000,
        vec![
            imu(200_000, Vector3::new(1.0, 0.0, 0.0), Vector3::zeros()),
            imu(250_000, Vector3::new(1.0, 0.0, 0.0), Vector3::zeros()),
        ],
    );
    assert!(
        extract_backward_imu_slice(200_000, &map).is_none(),
        "zero gyro rate"
    );
}

/// Backward integration negates the specific force and reverses the sample
/// order. Two samples 50 ms apart with `a = (1,0,0) m/s^2` become
/// `a = (-1,0,0)` in reverse order; the first step's `dt` saturates the
/// `clamp(1e-4, 0.1)` bound (the reversed timestamps read as a GPS week
/// wrap), so:
///   dt = 0.1 s
///   dv = (-1,0,0) * 0.1           = (-0.1, 0, 0)
///   dp = 0.5 * (-1,0,0) * 0.1^2  = (-0.005, 0, 0)
/// Both are the exact negation of the forward integral over the same 0.05 s
/// step would give scaled by the saturated interval.
#[test]
fn backward_imu_slice_negates_specific_force() {
    let mut map = BTreeMap::new();
    map.insert(
        100_000,
        vec![
            imu(
                100_000,
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.1, 0.0, 0.0),
            ),
            imu(
                150_000,
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.1, 0.0, 0.0),
            ),
        ],
    );
    let p = extract_backward_imu_slice(100_000, &map).expect("preintegration runs");
    assert!(
        p.dv.x < 0.0,
        "reversing time must negate the integrated force"
    );
    assert!(
        (p.dv - Vector3::new(-0.1, 0.0, 0.0)).norm() < 1e-12,
        "dv = {:?}",
        p.dv.as_slice()
    );
    assert!(
        (p.dp - Vector3::new(-0.005, 0.0, 0.0)).norm() < 1e-12,
        "dp = {:?}",
        p.dp.as_slice()
    );
    assert!(
        p.dt > 0.0,
        "a backward interval is still a positive duration"
    );
}

/// The long-window guard behaves exactly as in the forward pass: 30 samples
/// 100 ms apart with `a = 0.01 m/s^2` give dt = 2.9 s > 2.0 and an implied
/// speed well under 0.5 m/s, so both deltas are zeroed while `dt` survives.
#[test]
fn backward_imu_slice_zeroes_a_long_motionless_window() {
    let samples: Vec<ImuSample> = (0..30)
        .map(|i| {
            imu(
                i as u64 * 100_000,
                Vector3::new(0.01, 0.0, 0.0),
                Vector3::new(0.01, 0.0, 0.0),
            )
        })
        .collect();
    let mut map = BTreeMap::new();
    map.insert(0u64, samples);
    let p = extract_backward_imu_slice(0, &map).expect("preintegration runs");
    assert!(p.dt > 2.0, "dt = {}", p.dt);
    assert_eq!(p.dp, Vector3::zeros());
    assert_eq!(p.dv, Vector3::zeros());
}

/// A long window that really does accelerate keeps its deltas.
#[test]
fn backward_imu_slice_keeps_a_long_accelerating_window() {
    let samples: Vec<ImuSample> = (0..30)
        .map(|i| {
            imu(
                i as u64 * 100_000,
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.01, 0.0, 0.0),
            )
        })
        .collect();
    let mut map = BTreeMap::new();
    map.insert(0u64, samples);
    let p = extract_backward_imu_slice(0, &map).expect("preintegration runs");
    assert!(p.dt > 2.0);
    assert!(
        p.dv.norm() > 1.0,
        "a 1 m/s^2 ramp must leave a real dv, got {:?}",
        p.dv.as_slice()
    );
}

