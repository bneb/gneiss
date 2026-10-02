//! Unit tests for [`crate::post_process::iekf_pass`]. Sibling file:
//! `iekf_pass.rs` exceeds the 300-line threshold, so its tests live here and
//! are declared from `post_process/mod.rs`. Every expected number is derived by
//! hand in the comment above it.

use std::collections::HashMap;
use std::sync::{Mutex, MutexGuard};

use nalgebra::{DMatrix, DVector, Vector3};

use gneiss_core::obs::EpochObs;
use gneiss_core::time::GpsTime;

use crate::estimators::rtk_iekf::{DoubleDiffKey, GnssRtkIekf, IekfSnapshot};
use crate::post_process::dynamics::{ProcessingDynamics, KINEMATIC_INNOV_GATE_SCALE,
    STATIC_INNOV_GATE_SCALE};
use crate::post_process::iekf_pass::{configure_iekf, dump_amb_history,
    estimate_epoch_covariance, estimate_swfg_epoch_covariance, find_matched_base,
    ZWD_BASELINE_GATE_M};
use crate::post_process::ReceiverPcvPair;

// ---------------------------------------------------------------- helpers --

/// Every env var `configure_iekf` reads, cleared before and after each
/// guarded test so an opt-in can never leak into a sibling test.
const TRACKED_ENV: &[&str] = &["GNEISS_ELEV_DEG", "GNEISS_IONO_STATES", "GNEISS_SAT_IONO",
    "GNEISS_SP3", "GNEISS_CLK", "GNEISS_AMB_DUMP", "GNEISS_AMB_DUMP_DIR"];

static ENV_LOCK: Mutex<()> = Mutex::new(());

/// Serialises env mutation between test threads (the environment is
/// process-global) and scrubs [`TRACKED_ENV`] on entry and exit. `ENV_LOCK` is
/// not reentrant: never hold two at once.
struct EnvGuard(#[allow(dead_code)] MutexGuard<'static, ()>);

impl EnvGuard {
    fn acquire() -> Self {
        let guard = ENV_LOCK.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
        for key in TRACKED_ENV { std::env::remove_var(key); }
        Self(guard)
    }

    fn set(&self, key: &'static str, value: &str) {
        assert!(TRACKED_ENV.contains(&key), "untracked env var {key}");
        std::env::set_var(key, value);
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        for key in TRACKED_ENV { std::env::remove_var(key); }
    }
}

fn assert_close(actual: f64, expected: f64, tol: f64) {
    assert!((actual - expected).abs() <= tol, "expected {expected} +/- {tol}, got {actual}");
}

fn epoch(tow: f64) -> EpochObs {
    EpochObs { time: GpsTime::new(2000, tow), satellites: Vec::new() }
}

/// `n` epochs at `dt` spacing starting at `tow0`.
fn stream(tow0: f64, n: usize, dt: f64) -> Vec<EpochObs> {
    (0..n).map(|i| epoch(tow0 + i as f64 * dt)).collect()
}

/// Temp directory plus its path as a string (env vars are strings).
fn tmpdir() -> (tempfile::TempDir, String) {
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir.path().to_str().expect("utf8 path").to_string();
    (dir, path)
}

/// Rover/base geometry for one `configure_iekf` call. The baseline runs purely
/// along ECEF +X, so `baseline_m` is the exact `|rover - base|` norm.
struct Cfg {
    baseline_m: f64,
    widelane_ar: bool,
    tropo_grad: bool,
    enable_glonass: bool,
    rover_epochs: Vec<EpochObs>,
    sat_upd: Option<HashMap<u16, f64>>,
    pcv: Option<std::sync::Arc<ReceiverPcvPair>>,
}

impl Default for Cfg {
    fn default() -> Self {
        Self { baseline_m: 0.0, widelane_ar: false, tropo_grad: false, enable_glonass: false,
            rover_epochs: Vec::new(), sat_upd: None, pcv: None }
    }
}

fn build(c: Cfg) -> GnssRtkIekf {
    build_with(c, ProcessingDynamics::Static)
}

fn build_with(c: Cfg, dynamics: ProcessingDynamics) -> GnssRtkIekf {
    let init_pos = Vector3::new(-3_961_904.43, 3_348_994.27, 3_698_211.71);
    let base_pos = init_pos + Vector3::new(c.baseline_m, 0.0, 0.0);
    configure_iekf(init_pos, GpsTime::new(2000, 100.0), 1e-6, base_pos, dynamics, c.widelane_ar,
        c.tropo_grad, c.enable_glonass, c.pcv, &c.rover_epochs, c.sat_upd)
}

fn key(sat: u16, ref_sat: u16, band: u8) -> DoubleDiffKey {
    DoubleDiffKey { constellation_id: 0, sat, ref_sat, freq_band: band }
}

fn snapshot(tow: f64, keys: Vec<DoubleDiffKey>, x_post: Vec<f64>) -> IekfSnapshot {
    let n = x_post.len();
    let zero = DMatrix::zeros(n, n);
    IekfSnapshot {
        time: GpsTime::new(2000, tow),
        x_pred: DVector::zeros(n), p_pred: zero.clone(),
        x_post: DVector::from_vec(x_post), p_post: zero,
        f_mat: DMatrix::zeros(6, 6),
        is_fixed: false, n_sats: 8, quality: 1, amb_keys: keys,
    }
}

/// Minimal SP3: one epoch, one satellite, one 60-column record. `PG01` plus
/// four 14-column fields is 60 bytes, matching the parser's X/Y/Z slices.
fn sp3_text() -> String {
    let c = |v: f64| format!("{v:>14}");
    format!("#cP 2024 1 1 0 0 0.0 1\n*  2024  1  1  0  0  0.00000000\nPG01{}{}{}{}\nEOF\n",
        c(15_000_000.0), c(10_000_000.0), c(20_000_000.0), c(100.0))
}

/// Minimal RINEX CLK: two `AS` records, bias field at byte 40, 1e-4 s each.
fn clk_text() -> String {
    let rec = |prn: u8| {
        let head = format!("AS G{prn:02}  2024 01 01 00 00  0.000000  2");
        format!("{head:<40}{:<19}", "1.000000000000E-04")
    };
    format!("{}\n{}\n", rec(1), rec(5))
}

// ------------------------------------------------ covariance estimation --

/// Horizontal (or vertical) variance of a representative epoch covariance:
/// geom = max(8/n_sats.max(4), 0.5) times one of sigma fixed 0.01 / rtk 0.50 /
/// ppp 0.15 + 2.35/(1 + epochs/60) / else 2.5, then squared.
fn cov(n: usize, rtk: bool, fixed: bool, ppp: bool, epochs: usize) -> f64 {
    estimate_swfg_epoch_covariance(n, rtk, fixed, ppp, epochs)[(0, 0)]
}

fn cov_z(n: usize, rtk: bool, fixed: bool, ppp: bool, epochs: usize) -> f64 {
    estimate_swfg_epoch_covariance(n, rtk, fixed, ppp, epochs)[(2, 2)]
}

#[test]
fn covariance_branch_sigmas_match_hand_arithmetic() {
    assert_close(cov(8, true, true, false, 0), 1.0e-4, 1e-15);
    assert_eq!(cov(8, true, false, false, 0), 0.25);
    assert_eq!(cov(8, false, false, false, 0), 6.25);
    assert_close(cov(8, true, true, true, 60), 1.0e-4, 1e-15);
    assert_eq!(cov(8, true, false, true, 60), 0.25);
    assert_close(cov(8, false, false, true, 60), 1.325 * 1.325, 1e-14);
    // The test-only wrapper is the (is_ppp=false, epochs=0) specialisation;
    // six satellites give geom 4/3, sigma 2/3, var 4/9.
    // Precedence: is_fixed (0.01) > is_rtk (0.50) > is_ppp(60) (1.325).
    assert_close(estimate_epoch_covariance(8, true, true)[(0, 0)], 1.0e-4, 1e-14);
    assert_eq!(estimate_epoch_covariance(8, false, false)[(0, 0)], 6.25);
    assert_close(estimate_epoch_covariance(6, true, false)[(0, 0)], 4.0 / 9.0, 1e-14);
    // Diagonal, positive, finite, equal horizontals, and a vertical sigma
    // carrying an extra 1.5x (so its variance is 1.5^2 = 2.25x).
    for n in [4usize, 8, 12, 20, 31] {
        let m = estimate_swfg_epoch_covariance(n, false, false, true, 120);
        for (i, j) in [(0usize, 1usize), (0, 2), (1, 2)] {
            assert_eq!(m[(i, j)], 0.0, "covariance must stay diagonal");
        }
        let var = m[(0, 0)];
        assert!(var.is_finite() && var > 0.0);
        assert_eq!(m[(1, 1)], var, "horizontal variances must be equal");
        assert_close(m[(2, 2)] / var, 2.25, 1e-13);
    }
}

#[test]
fn ppp_covariance_decays_monotonically_toward_the_converged_sigma() {
    // decay = 1/(1 + n/60), sigma = 0.15 + 2.35*decay: n=0 -> 2.50 (var 6.25),
    // n=60 -> 1.325 (1.755625), n=120 -> 14/15 (196/225, zz*2.25 = 1.96),
    // n=59940 -> decay 1/1000 -> sigma 0.15235.
    assert_close(cov(8, false, false, true, 0), 6.25, 1e-14);
    assert_close(cov(8, false, false, true, 60), 1.755625, 1e-14);
    assert_close(cov(8, false, false, true, 120), 196.0 / 225.0, 1e-14);
    assert_close(cov_z(8, false, false, true, 120), 1.96, 1e-13);
    assert_close(cov(8, false, false, true, 59_940), 0.152_35 * 0.152_35, 1e-14);
    // Strictly decreasing, never reaching the converged sigma^2 = 0.0225.
    let mut prev = f64::INFINITY;
    for n in [0usize, 10, 60, 120, 600, 6000, 59_940] {
        let v = cov(8, false, false, true, n);
        assert!(v < prev, "PPP variance must fall with epochs (n={n})");
        assert!(v > 0.0225, "PPP variance must stay above the converged 0.0225");
        prev = v;
    }
}

#[test]
fn geometry_factor_floors_below_four_satellites_and_ceil_above_sixteen() {
    // n_sats.max(4) makes 0, 3 and 4 identical -> geom 2 -> fixed var 4e-4;
    // 8/16 = 0.5 exactly and .max(0.5) freezes 17 and 32 there too.
    let three = cov(3, true, true, false, 0);
    assert_close(three, 4.0e-4, 1e-14);
    assert_eq!(cov(0, true, true, false, 0), three);
    assert_eq!(cov(4, true, true, false, 0), three);
    let sixteen = cov(16, true, true, false, 0);
    assert_close(sixteen, 2.5e-5, 1e-14);
    assert_eq!(cov(17, true, true, false, 0), sixteen);
    assert_eq!(cov(32, true, true, false, 0), sixteen);
    let mut prev = f64::INFINITY;
    for n in 4..=40usize {
        let v = cov(n, true, true, false, 0);
        assert!(v <= prev, "variance must not grow with satellite count (n={n})");
        prev = v;
    }
}

// ------------------------------------------------- base epoch matching --

#[test]
fn find_matched_base_requires_base_data_within_the_window() {
    assert!(find_matched_base(100.0, None).is_none());
    assert!(find_matched_base(100.0, Some(&[])).is_none());
    // Strict |dtow| < 0.1 s; probing from tow = 0 keeps the subtraction exact.
    assert!(find_matched_base(0.1, Some(&[epoch(0.0)])).is_none());
    assert!(find_matched_base(0.0, Some(&[epoch(0.1)])).is_none());
    // 1 ns inside the window is accepted from either side.
    assert!(find_matched_base(0.0, Some(&[epoch(0.099_999_999)])).is_some());
    assert!(find_matched_base(0.1, Some(&[epoch(0.000_000_000_1)])).is_some());
    // Candidates at -0.05, +0.02, +0.09 s: abs() before comparing. Leading
    // epochs count too, and reordering must not change the winner.
    let three = vec![epoch(99.95), epoch(100.02), epoch(100.09)];
    assert_eq!(find_matched_base(100.0, Some(&three)).expect("match").time.tow, 100.02);
    for cand in [vec![epoch(99.95), epoch(99.98)], vec![epoch(100.02), epoch(99.95)]] {
        let want = if cand[0].time.tow < 99.99 { 99.98 } else { 100.02 };
        assert_eq!(find_matched_base(100.0, Some(&cand)).expect("match").time.tow, want);
    }
}

// -------------------------------------------------- filter configuration --

#[test]
fn configure_iekf_defaults_leave_every_opt_in_off() {
    let iekf = build(Cfg::default());
    let mobile = build_with(Cfg::default(), ProcessingDynamics::Kinematic);
    assert_eq!(mobile.robust_innov_scale, KINEMATIC_INNOV_GATE_SCALE);
    assert_eq!(KINEMATIC_INNOV_GATE_SCALE, 3.0);
    assert!(mobile.is_kinematic);
    // A kinematic rover has no monument to lock and no static ZWD state.
    assert_eq!(mobile.static_lock_after_s, None);
    assert!(!mobile.state.zwd_enabled);
    assert_eq!(iekf.q_accel, 1e-6);
    assert_eq!(iekf.start_tow, 100.0);
    assert_eq!(iekf.min_elevation_rad, 0.1745, "10 degree default mask");
    assert_eq!(iekf.robust_innov_scale, STATIC_INNOV_GATE_SCALE);
    assert!(!iekf.is_kinematic);
    assert!(!iekf.widelane_ar && !iekf.enable_glonass && !iekf.track_ambiguity_keys);
    assert!(!iekf.state.zwd_enabled && !iekf.state.grad_enabled);
    assert!(iekf.precise_orbits.is_none() && iekf.precise_clocks.is_none());
    assert!(iekf.receiver_pcv.is_none());
    assert_eq!(iekf.slip_detector.cadence_hint_s, None);
    assert!(iekf.wl_tracker.sat_upd.is_none());
}

#[test]
fn configure_iekf_static_lock_applies_only_to_static_sessions() {
    let locked = build(Cfg { widelane_ar: true, ..Default::default() });
    assert_eq!(locked.static_lock_after_s, Some(900.0));
    assert_eq!(locked.static_lock_q_accel, 1e-8);
    // Without widelane_ar the loose Q stays for the whole session.
    assert_eq!(build(Cfg::default()).static_lock_after_s, None);
    // Nor for a kinematic rover on the same configuration.
    let mobile = build_with(Cfg { widelane_ar: true, ..Default::default() }, ProcessingDynamics::Kinematic);
    assert_eq!(mobile.static_lock_after_s, None);
}

#[test]
fn configure_iekf_gates_the_zwd_state_on_baseline_length() {
    let inside = build(Cfg { widelane_ar: true, baseline_m: ZWD_BASELINE_GATE_M - 1.0, ..Default::default() });
    assert!(inside.state.zwd_enabled);
    // enable_zwd(0.0225) lands the ~15 cm zenith wet init variance at column 6.
    assert_eq!(inside.state.cov[(6, 6)], 0.0225);
    // At the gate the strict `<` excludes it, as does one metre beyond.
    for m in [ZWD_BASELINE_GATE_M, ZWD_BASELINE_GATE_M + 1.0] {
        assert!(!build(Cfg { widelane_ar: true, baseline_m: m, ..Default::default() }).state.zwd_enabled);
    }
    // Short baseline without widelane_ar stays legacy: no ZWD state.
    assert!(!build(Cfg { baseline_m: 1000.0, ..Default::default() }).state.zwd_enabled);
    // Gradients are opt-in and can only exist alongside the ZWD state, which
    // reserves columns 6 (ZWD), 7-8 (gradients), so ambiguities start at 9.
    assert!(!build(Cfg { widelane_ar: true, ..Default::default() }).state.grad_enabled);
    let grad = build(Cfg { widelane_ar: true, tropo_grad: true, ..Default::default() });
    assert!(grad.state.grad_enabled);
    assert_eq!(grad.state.grad_idx().expect("gradients enabled"), (7, 8));
    assert_eq!(grad.state.amb_offset(), 9, "6 pose/vel + 1 ZWD + 2 gradients");
    let far = Cfg { widelane_ar: true, tropo_grad: true, baseline_m: ZWD_BASELINE_GATE_M + 1.0, ..Default::default() };
    assert!(!build(far).state.grad_enabled, "no gradients without the ZWD state");
}

#[test]
fn configure_iekf_tropo_gradients_are_opt_in_inside_the_zwd_gate() {
    assert!(!build(Cfg { widelane_ar: true, ..Default::default() }).state.grad_enabled);
    let on = build(Cfg { widelane_ar: true, tropo_grad: true, ..Default::default() });
    assert!(on.state.grad_enabled);
    // Gradients cannot exist without the ZWD state they attach to: (7, 8).
    assert_eq!(on.state.grad_idx().expect("gradients enabled"), (7, 8));
    assert_eq!(on.state.amb_offset(), 9, "6 pose/vel + 1 ZWD + 2 gradients");
    // Outside the gate, requesting gradients still yields nothing.
    let gated = build(Cfg {
        widelane_ar: true, tropo_grad: true,
        baseline_m: ZWD_BASELINE_GATE_M + 1.0,
        ..Default::default()
    });
    assert!(!gated.state.grad_enabled);
}

#[test]
fn configure_iekf_glonass_is_independent_of_baseline_and_ambiguity_resolution() {
    // GLONASS used to be nested in the ZWD baseline gate; both escapes hold,
    // and the caller's opt-out still wins.
    assert!(build(Cfg { enable_glonass: true, baseline_m: 40_000.0, ..Default::default() }).enable_glonass);
    assert!(build(Cfg { enable_glonass: true, widelane_ar: true, ..Default::default() }).enable_glonass);
    assert!(!build(Cfg { widelane_ar: true, ..Default::default() }).enable_glonass);
}

#[test]
fn configure_iekf_cadence_hint_reflects_the_rover_stream_sampling_rate() {
    // Six epochs 30 s apart: five diffs, median index 5/2 = 2 -> 30.0 s.
    let slow = build(Cfg { rover_epochs: stream(0.0, 6, 30.0), ..Default::default() });
    assert_eq!(slow.slip_detector.cadence_hint_s, Some(30.0));
    assert_eq!(slow.base_slip_detector.cadence_hint_s, Some(30.0));
    // 1 s spacing stays below the 2 s threshold: no hint (legacy behaviour).
    let fast = build(Cfg { rover_epochs: stream(0.0, 6, 1.0), ..Default::default() });
    assert_eq!(fast.slip_detector.cadence_hint_s, None);
    // Four epochs give three diffs, under the five-difference minimum.
    let short = build(Cfg { rover_epochs: stream(0.0, 4, 30.0), ..Default::default() });
    assert_eq!(short.slip_detector.cadence_hint_s, None);
    // Set regardless of widelane_ar: it describes the data, not the AR mode.
    let fixed = build(Cfg { widelane_ar: true, rover_epochs: stream(0.0, 6, 30.0), ..Default::default() });
    assert_eq!(fixed.slip_detector.cadence_hint_s, Some(30.0));
}

#[test]
fn configure_iekf_installs_the_network_solved_satellite_upd_map() {
    let mut upd = HashMap::new();
    upd.insert(1u16, 0.5);
    upd.insert(7u16, -0.25);
    for ar in [false, true] { // installed verbatim, never truncated
        let iekf = build(Cfg { sat_upd: Some(upd.clone()), widelane_ar: ar, ..Default::default() });
        assert_eq!(iekf.wl_tracker.sat_upd.as_ref(), Some(&upd));
    }
    assert_eq!(upd.len(), 2);
}

#[test]
fn configure_iekf_installs_a_receiver_pcv_pair_when_supplied() {
    let Ok(db) = gneiss_parsers::antex::AntexDatabase::parse("../../datasets/igs14.atx") else {
        return; // dataset absent in this checkout
    };
    use gneiss_parsers::receiver_antenna::ReceiverAntenna;
    let rover = std::sync::Arc::new(
        ReceiverAntenna::lookup(&db, "TRM59800.00", "SCIT").expect("igs14 rover"));
    let base = std::sync::Arc::new(
        ReceiverAntenna::lookup(&db, "ASH701945B_M", "SCIT").expect("igs14 base"));
    let iekf = build(Cfg {
        pcv: Some(std::sync::Arc::new(ReceiverPcvPair { rover: rover.clone(), base: base.clone() })),
        ..Default::default()
    });
    let (got_rover, got_base) = iekf.receiver_pcv.as_ref().expect("pair installed");
    assert!(std::sync::Arc::ptr_eq(got_rover, &rover), "rover antenna identity");
    assert!(std::sync::Arc::ptr_eq(got_base, &base), "base antenna identity");
}

// ---------------------------------------------------- environment opt-ins --

#[test]
fn env_elevation_and_iono_state_overrides_default_off() {
    let env = EnvGuard::acquire();
    let plain = build(Cfg::default());
    assert!(!plain.state.iono_enabled && !plain.state.sat_iono_enabled);
    env.set("GNEISS_IONO_STATES", "1");
    env.set("GNEISS_SAT_IONO", "1");
    assert!(build(Cfg::default()).state.iono_enabled);
    env.set("GNEISS_IONO_STATES", "true"); // only the exact string "1" opts in
    env.set("GNEISS_SAT_IONO", "0");
    let partial = build(Cfg::default());
    assert!(!partial.state.iono_enabled && !partial.state.sat_iono_enabled);
    env.set("GNEISS_ELEV_DEG", "15");
    // 15 deg = 15 * pi/180 = 0.2617993877991494 rad.
    assert_close(build(Cfg::default()).min_elevation_rad, 0.261_799_387_799_149_4, 1e-15);
    env.set("GNEISS_ELEV_DEG", "not-a-number");
    // A parse failure leaves the 10 degree default in place.
    assert_eq!(build(Cfg::default()).min_elevation_rad, 0.1745);
    // GNEISS_AMB_DUMP is gated on widelane_ar; it alone is the switch then.
    env.set("GNEISS_AMB_DUMP", "1");
    assert!(!build(Cfg::default()).track_ambiguity_keys);
    assert!(build(Cfg { widelane_ar: true, ..Default::default() }).track_ambiguity_keys);
}

// -------------------------------------------------------- precise products --

#[test]
fn precise_orbits_load_only_when_gated_in_and_the_sp3_parses() {
    let (_dir, base) = tmpdir();
    let path = format!("{base}/test.sp3");
    std::fs::write(&path, sp3_text()).expect("write sp3");
    let env = EnvGuard::acquire();
    // Gate closed (no widelane_ar): the env var is never read at all.
    env.set("GNEISS_SP3", &path);
    assert!(build(Cfg::default()).precise_orbits.is_none());
    let loaded = build(Cfg { widelane_ar: true, ..Default::default() });
    assert_eq!(loaded.precise_orbits.as_ref().expect("SP3 loaded").len(), 1, "PG01 only");
    env.set("GNEISS_SP3", &format!("{base}/missing.sp3")); // error branch
    assert!(build(Cfg { widelane_ar: true, ..Default::default() }).precise_orbits.is_none());
}

#[test]
fn precise_clocks_load_from_a_rinex_clk_file() {
    let (_dir, base) = tmpdir();
    let path = format!("{base}/test.clk");
    std::fs::write(&path, clk_text()).expect("write clk");
    let env = EnvGuard::acquire();
    assert!(build(Cfg { widelane_ar: true, ..Default::default() }).precise_clocks.is_none());
    env.set("GNEISS_CLK", &path);
    let with_clk = build(Cfg { widelane_ar: true, ..Default::default() });
    assert_eq!(with_clk.precise_clocks.as_ref().expect("CLK loaded").satellites.len(), 2, "G01, G05");
    env.set("GNEISS_CLK", &format!("{base}/missing.clk"));
    assert!(build(Cfg { widelane_ar: true, ..Default::default() }).precise_clocks.is_none());
}

// ------------------------------------------------------ ambiguity history --

/// Union key order is first-seen across history (K1, then K2). Headers use the
/// raw key fields, so GPS cid 0 renders `0_1_5_b1`. offset = 6 + zwd +
/// 2*grad = 6 here, and each snapshot's own `amb_keys` index selects the
/// column, not the union index: t=100 [K1] -> local 0 -> x_post[6] = 1.5;
/// t=130 [] -> two empty fields; t=160 [K2, K1] -> K1 local 1 -> 4.5 and K2
/// local 0 -> 2.5; t=190 [K1] with x_post len 6 -> index out of range -> empty.
#[test]
fn dump_amb_history_writes_the_key_union_at_the_ambiguity_offset() {
    let (dir, base) = tmpdir();
    let env = EnvGuard::acquire();
    env.set("GNEISS_AMB_DUMP_DIR", &base);
    let (k1, k2) = (key(1, 5, 1), key(7, 5, 2));
    let mut iekf = build(Cfg::default()); // no widelane_ar => no ZWD/gradients
    assert_eq!(iekf.state.amb_offset(), 6, "no atmospheric states to skip");
    iekf.history = vec![
        snapshot(100.0, vec![k1], vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5]),
        snapshot(130.0, vec![], vec![0.0; 9]),
        snapshot(160.0, vec![k2, k1], vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.5, 4.5]),
        snapshot(190.0, vec![k1], vec![0.0; 6]),
    ];
    dump_amb_history(&iekf, "fwd");
    let csv = std::fs::read_to_string(dir.path().join("amb_fwd.csv")).expect("dump written");
    let want = concat!("tow,0_1_5_b1,0_7_5_b2\n", "100,1.500000,\n", "130,,\n",
        "160,4.500000,2.500000\n", "190,,\n");
    assert_eq!(csv, want);
}

#[test]
fn dump_amb_history_shifts_the_offset_by_the_atmospheric_state_columns() {
    let (dir, base) = tmpdir();
    let env = EnvGuard::acquire();
    env.set("GNEISS_AMB_DUMP_DIR", &base);
    let k1 = key(1, 5, 1);
    // Short-baseline static session with widelane_ar reserves column 6 for
    // ZWD, so ambiguities start at 7.
    let mut zwd = build(Cfg { widelane_ar: true, ..Default::default() });
    assert!(zwd.state.zwd_enabled && !zwd.state.grad_enabled);
    assert_eq!(zwd.state.amb_offset(), 7);
    zwd.history = vec![snapshot(100.0, vec![k1], vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 8.25])];
    dump_amb_history(&zwd, "zwd");
    let csv = std::fs::read_to_string(dir.path().join("amb_zwd.csv")).expect("dump written");
    assert_eq!(csv, "tow,0_1_5_b1\n100,8.250000\n");
    // Adding the two gradient columns moves the ambiguities to 9.
    let mut grad = build(Cfg { widelane_ar: true, tropo_grad: true, ..Default::default() });
    assert!(grad.state.zwd_enabled && grad.state.grad_enabled);
    assert_eq!(grad.state.amb_offset(), 9);
    grad.history = vec![snapshot(100.0, vec![k1], vec![9.0; 10])];
    dump_amb_history(&grad, "grad");
    let csv = std::fs::read_to_string(dir.path().join("amb_grad.csv")).expect("dump written");
    assert_eq!(csv, "tow,0_1_5_b1\n100,9.000000\n");
}

#[test]
fn dump_amb_history_writes_nothing_without_keys_or_without_a_directory() {
    let dir = tempfile::tempdir().expect("tempdir");
    let mut iekf = build(Cfg::default());
    // No tracked keys: early return before the env var is consulted.
    iekf.history = vec![snapshot(100.0, vec![], vec![0.0; 9])];
    dump_amb_history(&iekf, "nokeys");
    assert!(!dir.path().join("amb_nokeys.csv").exists());
    // Keys present but no GNEISS_AMB_DUMP_DIR: nowhere to write.
    iekf.history = vec![snapshot(100.0, vec![key(1, 5, 1)], vec![0.0; 9])];
    drop(EnvGuard::acquire());
    dump_amb_history(&iekf, "nodir");
    assert!(std::fs::read_dir(dir.path()).expect("readable").next().is_none());
}
