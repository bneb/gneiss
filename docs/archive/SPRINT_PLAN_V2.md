> **Superseded.** This document describes an earlier architecture (pre network-RTK/PPK pivot) and is kept for historical record only. For current status and roadmap, see `docs/PROJECT_STATUS.md` and `docs/NETWORK_RTK_NEXT_STEPS.md`.

# Sprint Plan v10 — The Conference Edition

**Produced by**: Dellaert, Teunissen, Humphreys, Molteno, Dampf, Bisnath
**Synthesized for**: Gneiss coding agent
**Date**: 2026-07-24

---

## Opening Remarks — Where We Actually Stand

**Dampf opens with a reality check.** "I've read your codebase. You have a single-epoch IEKF with per-epoch SPP reset. That architecture has a ~5m accuracy floor no matter what you bolt onto it. NovAtel killed that approach in 2004. The fact that you beat RTKLIB on some modes is impressive — RTKLIB has the same architecture. But Qinertia doesn't. Qinertia carries position and ambiguity states across epochs without resetting. Until you do that, nothing else matters."

**Teunissen nods.** "I also read your codebase. Your LAMBDA is correct. Your FFRT thresholds are correct. Your fix rate is 18.9% because your float ambiguity covariance is inflated by the per-epoch SPP reset. Every epoch, position jumps by ~5m, which feeds into the ambiguity update as process noise. The ambiguities never converge. Fix rate is a symptom, not a cause."

**Dellaert draws on the whiteboard.** "What you have is a graphical model where the prior on position is independent at each epoch. What you need is a graphical model where position at epoch k is connected to position at epoch k+1 through a motion model. That's a factor graph. You already have the optimizer — your `estimators/factor_graph` module is correct. You already have the factors — `gnss_factors.rs` has pseudorange, carrier phase, and doppler. What's missing is sliding-window state management and marginalization. That's ~2000 lines of code, not a rewrite."

**Bisnath adds:** "And once you have that, PPP-AR with precise products will get you to 4cm. The IGS final orbit/clock products are free. The CDDIS Earthdata registration is free. The problem isn't data availability — it's that your current architecture can't use the data effectively."

**Humphreys pushes back gently.** "4cm is great for surveyors. But the market for this software is autonomous vehicles and smartphones — low-cost receivers in urban environments. For that market, you need to handle multipath, not just orbit error. Your PR validation gate is a step in the right direction. But the real win is IMU tight coupling — when GNSS drops out under a bridge, the IMU keeps the solution alive. Qinertia's primary value proposition isn't PPP accuracy — it's robustness in challenging environments."

**Molteno, quietly:** "And none of this matters if the configuration is wrong. I've seen production GNSS systems fail because someone set `max_base_age_s` to 0.5 when the base station was 30 seconds away. Your `EngineConfig` has 60 flat fields with no validation. If you're going to redesign this, put the type system to work."

---

## Cocktail Hour — Raw, Unsolicited Individual Inputs

### Dampf (production firmware, 3 drinks in)

"Here's what actually matters when you ship GNSS software:"

```rust
// This is the only function signature your users care about:
fn process_epoch(
    rover: &EpochObs,
    base: Option<&EpochObs>,
    imu: Option<&[ImuSample]>,
) -> Result<PositionSolution, EngineError>;

// And this is all they want in the output:
struct PositionSolution {
    time: GpsTime,
    ecef: Vector3<f64>,
    covariance: Matrix3<f64>,  // 3x3 ECEF
    mode: SolutionMode,        // SPP, Float, Fixed, INS
    hpl: f64,
    vpl: f64,
    n_sats: usize,
    fix_ratio: Option<f64>,
    // That's it. Nothing else. No internal state leakage.
}
```

"Your current API exposes `RtkState` with 40 public fields. That's an abstraction violation. Every field you expose is a field you can never change. Hide the state. Expose the solution. Everything else is internal."

"Also, your error handling is wrong. GNSS is fundamentally unreliable — satellites go behind buildings, base stations drop out, IMUs drift. Your processing loop must never crash. Not on NaN. Not on empty observations. Not on singular matrices. Every error path must degrade gracefully to the next-best solution mode. If RTK fails, fall back to float. If float fails, fall back to SPP. If SPP fails, coast on IMU. If IMU diverges, output nothing and wait. This is the hierarchy. Every GNSS receiver implements it. Yours doesn't."

### Teunissen (AR, nursing a single scotch)

"The 18.9% fix rate. Here's why."

"Your float ambiguity covariance after the EKF update has two components: the measurement contribution (which shrinks with more epochs) and the process noise contribution (which grows with time). With per-epoch SPP reset, the process noise on position is ~100 m² per epoch. That propagates into the ambiguity states through the measurement geometry. After 10 epochs, your ambiguity variance is still dominated by process noise, not by measurement averaging."

"The fix: remove the position reset. After epoch 1, the position uncertainty from SPP is ~100 m². After epoch 10 with continuous filtering and no reset, it should be ~0.1 m². That's a 1000× reduction. The ambiguity variance follows the same trajectory. With σ_position = 0.3m and σ_ambiguity = 0.5 cycles, LAMBDA will fix correctly >99% of the time."

"But there's a subtlety. You're estimating ionosphere as random walk states in UDUC mode. Those states ALSO get reset by the SPP anchor. You need to carry them forward too. And the ZWD troposphere state. Basically anything in the state vector that's not white noise needs to persist across epochs."

"The multi-epoch factor graph handles this naturally — every variable persists for the duration of the window. But you can also fix it in the IEKF by simply not calling `reset_to_spp` after epoch 1. The SPP prior becomes a one-time initialization, not a per-epoch anchor."

### Humphreys (urban/low-cost, animated)

"Forget survey-grade. The market is u-blox F9P receivers in cars and drones. These receivers have 2-5m code multipath in urban environments. Single-epoch RTK is limited by code multipath, not by your estimation algorithm. Your own 2026-07-03 memo says this: 'AR is correct; float p50=1.5m is measurement-limited.'"

"So why are you spending 80% of your effort on features that only help when measurement quality is already good? IONEX, precise products, GLONASS IFB — these reduce the 5-20cm of residual ionosphere/orbit error. They don't touch the 2-5m of code multipath. For urban positioning, the highest-ROI improvements are:"

"1. **C/N0-based variance scaling.** Your measurement variance model uses `pr_base_var = 0.5 m²` for all satellites regardless of signal strength. A satellite at 25 dB-Hz has 10× the noise of one at 45 dB-Hz. Scale variance by C/N0. The SNR-based variance model already exists in `EkfTuningConfig` — use it."

"2. **Elevation-dependent variance.** Low-elevation satellites have more multipath and more atmospheric error. Scale variance by `1/sin(el)`. Every GNSS receiver does this. Yours doesn't."

"3. **IMU tight coupling.** In urban canyons, you have 2-4 visible satellites on each side of the street. That's not enough for a position fix. The IMU bridges outages and provides the motion constraint that makes the factor graph well-conditioned with partial satellite visibility. This is the single biggest accuracy improvement for urban positioning."

"4. **Robust estimation in the measurement domain, not the position domain.** Your chi-square gating rejects measurements that don't fit the current estimate. But in urban environments, the current estimate might be wrong because it's based on contaminated measurements. Use Huber/Cauchy M-estimation in the factor graph so that outliers are downweighted rather than rejected. You already have this infrastructure in `estimators/factor_graph` — use it for every measurement, not just the EKF fallback."

### Dellaert (factor graphs, at the whiteboard at 11pm)

"Let me show you the variable set. This is the whole design:"

```rust
/// A variable in the factor graph is identified by (type, satellite, time).
/// Variables are created once and persist across epochs until marginalized.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum VariableKey {
    /// Rover position at a specific epoch (0, 1, ..., window_size-1).
    Position(u32),
    /// Rover velocity at a specific epoch.
    Velocity(u32),
    /// Rover attitude (quaternion) at a specific epoch.
    Attitude(u32),
    /// Receiver clock bias per constellation per epoch.
    ClockBias { epoch: u32, constellation: Constellation },
    /// Troposphere zenith wet delay at a specific epoch.
    TropoZwd(u32),
    /// Ionosphere slant delay for a satellite-epoch pair.
    IonoSlant { epoch: u32, satellite: SatelliteId },
    /// Carrier-phase ambiguity for a satellite-frequency pair.
    /// Persists across ALL epochs — this is the key difference from an EKF.
    Ambiguity { satellite: SatelliteId, frequency: u8 },
    /// GLONASS IFB slope (one per receiver session).
    IfbGlonass,
    /// IMU accelerometer bias (slowly varying, shared across epochs).
    AccelBias,
    /// IMU gyroscope bias (slowly varying, shared across epochs).
    GyroBias,
}
```

"Ambiguities are NOT per-epoch. They're the same variable across the entire window. That's the fundamental insight. When you do EKF per-epoch, you're creating a new ambiguity variable every epoch and throwing away the previous one. The factor graph treats ambiguity as a single variable connected to carrier-phase factors at every epoch it appears. This is physically correct — the integer ambiguity doesn't change unless there's a cycle slip."

"Here's what a factor looks like:"

```rust
trait Factor: Debug {
    /// Which variables this factor connects to (in order).
    fn variables(&self) -> &[VariableKey];

    /// Residual r(x) = h(x) - z.  Dimension: measurement dim.
    fn residual(&self, values: &VariableValues) -> DVector<f64>;

    /// Jacobian J = ∂r/∂x evaluated at `values`.  Dimension: meas_dim × total_dim.
    /// The total_dim is sum of dimensions of all connected variables.
    fn jacobian(&self, values: &VariableValues) -> DMatrix<f64>;

    /// Information matrix (inverse measurement covariance).  meas_dim × meas_dim.
    fn information(&self) -> DMatrix<f64>;

    /// Optional robust loss function (Huber, Cauchy, none).
    fn robust_threshold(&self) -> Option<f64> { None }
}

/// A concrete factor: pseudorange on satellite G01 at epoch 3.
/// Connects to variables: Position(3), ClockBias(3, GPS), TropoZwd(3).
struct PseudorangeFactor {
    epoch: u32,
    satellite: SatelliteId,
    observed_pr_m: f64,
    variance_m2: f64,
    sat_position_ecef: Vector3<f64>,
    sat_clock_m: f64,
    elevation_rad: f64,
}
```

"The optimizer is standard Levenberg-Marquardt over the variable set:"

```rust
struct SlidingWindowOptimizer {
    window: VecDeque<EpochData>,
    variables: BTreeMap<VariableKey, VariableState>,
    max_window_size: usize,
}

struct VariableState {
    /// Current estimate (the value at which we linearize).
    value: DVector<f64>,
    /// Dimension of this variable (3 for position, 1 for clock, etc.).
    dim: usize,
    /// Whether this variable should be marginalized out.
    marginalized: bool,
}

impl SlidingWindowOptimizer {
    fn solve(&mut self) -> Result<DVector<f64>, SolveError> {
        let total_dim: usize = self.variables.values().map(|v| v.dim).sum();

        // Build linear system J^T W J Δx = -J^T W r
        for iteration in 0..self.max_iterations {
            let jtj = self.build_hessian(total_dim);
            let jtr = self.build_gradient(total_dim);

            // Solve with Cholesky, fall back to SVD if near-singular
            let delta = solve_linear_system(&jtj, &jtr)?;

            // Update variable values
            self.apply_delta(&delta);

            // Check convergence
            if delta.norm() < self.convergence_tol {
                break;
            }
        }

        // Marginalize old variables via Schur complement
        // (only when window is full and we're pushing a new epoch)
        self.marginalize_oldest_epoch()?;

        Ok(self.extract_current_state())
    }
}
```

"The Schur complement marginalization is the key computational trick. When the window has `N` epochs and you add epoch `N+1`, you don't want to keep all `N+1` epochs in memory. You want to marginalize out epoch 0 (the oldest), keeping only the information it provides about the remaining variables. The Schur complement does this:"

```rust
fn marginalize_oldest_epoch(&mut self) {
    // Identify variables that appear ONLY in the oldest epoch
    // (position, clock, tropo at epoch 0 — but NOT ambiguities
    //  since they connect to multiple epochs)
    let (marginalized, remaining): (Vec<_>, Vec<_>) =
        self.variables.keys().partition(|k| is_epoch_specific(*k, oldest_epoch));

    // Partition the Hessian: [A  B; B^T  C]
    // A = marginalized-marginalized, C = remaining-remaining
    // Schur complement of A: S = C - B^T A^{-1} B
    // This encodes the information the marginalized variables
    // provided about the remaining variables.
    let s = compute_schur_complement(&hessian, &marginalized, &remaining);

    // Remove marginalized variables from the active set.
    // Add S as a prior on the remaining variables.
    // The window now has N-1 explicit epochs + 1 implicit (marginalized).
}
```

"This is the same technique that iSAM2, GTSAM, and every modern SLAM system uses. For a GNSS-only window with 10 epochs and 20 satellites, the total dimension is ~250 variables × ~3 dim each = ~750 parameters. Marginalizing one epoch takes ~50×50 Cholesky ≈ microseconds. The whole solve-marginalize cycle takes < 10ms on a modern CPU — well within real-time constraints at 10 Hz."

### Molteno (Rust API design, over coffee the next morning)

"The biggest source of bugs in your current codebase isn't the math. It's the config. Here's what type-safe configuration looks like:"

```rust
// Config that can't express nonsense.
// Each mode has its own config struct with only the relevant fields.

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "mode")]
enum EngineMode {
    Spp(SppConfig),
    Ppp(PppConfig),
    Rtk(RtkConfig),
    RtkIns(RtkInsConfig),
    PppIns(PppInsConfig),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PppConfig {
    /// Path to precise products (SP3, CLK, BIA, ANTEX).
    #[serde(default)]
    precise_products: Option<PreciseProductConfig>,

    /// Ionosphere model.
    #[serde(default)]
    ionosphere: IonosphereConfig,

    /// Ambiguity resolution.
    #[serde(default)]
    ar: Option<ArConfig>,

    /// Window size for multi-epoch estimation.
    #[serde(default = "default_window")]
    window_size: usize,

    /// Initial position (from RINEX header or known survey mark).
    initial_position: Option<[f64; 3]>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct RtkConfig {
    /// Base station position (ECEF, meters). REQUIRED.
    base_position: [f64; 3],

    /// Base station observation source.
    base_source: BaseSource,

    /// Maximum base observation age (seconds).
    #[serde(default = "default_max_base_age")]
    max_base_age_s: f64,

    /// Ambiguity resolution.
    #[serde(default)]
    ar: Option<ArConfig>,

    /// IMU configuration (only for RtkIns mode).
    #[serde(default)]
    imu: Option<ImuConfig>,

    /// Initial position.
    initial_position: Option<[f64; 3]>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct RtkInsConfig {
    rtk: RtkConfig,
    imu: ImuConfig,  // REQUIRED — compile-time guarantee
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ImuConfig {
    /// IMU-to-antenna lever arm in body frame (meters). REQUIRED.
    lever_arm: [f64; 3],

    /// IMU mounting angles [roll, pitch, yaw] in radians.
    #[serde(default)]
    mounting_angles: [f64; 3],

    /// Enable non-holonomic constraints.
    #[serde(default)]
    enable_nhc: bool,

    /// NHC lever arm (IMU to rear axle).
    #[serde(default)]
    nhc_lever_arm: [f64; 3],

    /// Process noise tuning.
    #[serde(default)]
    tuning: ImuTuningConfig,
}

// Compile-time validation via constructor:
impl RtkInsConfig {
    pub fn new(rtk: RtkConfig, imu: ImuConfig) -> Self {
        // Validation happens here, at construction time.
        // If it constructs, it's valid.
        assert!(imu.lever_arm.iter().any(|&x| x != 0.0),
            "IMU lever arm must be non-zero");
        Self { rtk, imu }
    }
}
```

"The key pattern: impossible states are unrepresentable. You can't create an `RtkInsConfig` without an `ImuConfig`. You can't create a `PppConfig` with AR enabled but no window size. The type system enforces consistency. No runtime `EngineError::InvalidConfig` — if the config deserializes and constructs, it's valid."

"Second big API issue: your current code leaks internal state everywhere. `RtkState` has 40+ public fields. `ProcessingEngine` has 30+ public methods. This is a testing and maintenance nightmare. Hide the state:"

```rust
/// The ONLY public type users interact with.
pub struct Engine {
    // Private! Not accessible to users.
    state: Option<Box<dyn EstimatorState>>,
    config: EngineConfig,
    corrections: CorrectionStack,
}

impl Engine {
    /// Create a new engine from a validated configuration.
    pub fn new(config: EngineConfig) -> Result<Self, EngineError> { ... }

    /// Process one epoch of observations.  This is the MAIN entry point.
    pub fn process_epoch(
        &mut self,
        rover: EpochObs,
        base: Option<EpochObs>,
        imu: Option<&[ImuSample]>,
    ) -> Result<Solution, EngineError> { ... }

    /// Get the current solution without processing new data.
    pub fn current_solution(&self) -> Option<&Solution> { ... }

    /// Reset the estimator (e.g., after a cycle slip or large outage).
    pub fn reset(&mut self) { ... }
}

/// The ONLY public output type.
#[derive(Debug, Clone, Serialize)]
pub struct Solution {
    pub time: GpsTime,
    pub position: Position,
    pub velocity: Option<Velocity>,
    pub attitude: Option<Attitude>,
    pub mode: SolutionMode,
    pub protection_levels: Option<ProtectionLevels>,
    pub diagnostics: Diagnostics,
}

#[derive(Debug, Clone, Serialize)]
pub struct Diagnostics {
    pub n_satellites: usize,
    pub n_ambiguities_fixed: usize,
    pub fix_ratio: Option<f64>,
    pub rms_residual_pr_m: f64,
    pub rms_residual_cp_m: f64,
    pub pdop: f64,
    pub solver_iterations: usize,
    pub processing_time_us: u64,
}
```

"This API has exactly two public methods users care about: `new` and `process_epoch`. Everything else is an implementation detail. You can rewrite the entire internals without breaking any caller."

### Bisnath (PPP-AR, over dinner)

"The performance target is wrong. You're comparing yourself to Qinertia, which is a PPK engine for surveyors. The real competition is CSRS-PPP, magicGNSS, and the IGS combination products. These are the PPP services that researchers use. They achieve 2-4cm horizontal in static mode with daily observations, and 5-10cm in kinematic mode. That's the benchmark."

"For PPP-AR to work, you need four things:

1. **Precise satellite orbits and clocks** (SP3 + CLK files from IGS or CODE). You have these parsers. Use CODE MGEX products — they're multi-GNSS and include GLONASS, Galileo, BeiDou.

2. **Satellite phase biases** (BIA/SINEX bias files). These are the fractional-cycle biases that prevent PPP ambiguity resolution. Without them, you can only fix widelane, not narrowlane. The CODE bias products include these. Your `gneiss-parsers` already has SINEX bias parsing.

3. **Satellite and receiver antenna calibrations** (ANTEX files). The phase center variation is 1-3cm at L1 and needs to be corrected for centimeter PPP. You have ANTEX parsing. Are you applying PCV corrections to every measurement? I couldn't tell from the code.

4. **At least 30 minutes of continuous data.** PPP convergence takes time. The ionosphere, troposphere, and ambiguity states need enough geometry change to decorrelate from position. With multi-GNSS (GPS+GLO+GAL+BDS), 30 minutes gives ~80 satellites worth of geometry diversity. Single-GNSS takes 60+ minutes."

"The GICI-LIB dataset is perfect for validation — it includes CODE products, reference trajectories, and multiple environments. Run Gneiss on GICI-LIB with PPP-AR enabled and compare against the GICI reference solution. Until you've done that, you don't know whether your PPP implementation works."

---

## Roundtable — Heated Debates and Resolutions

### Debate 1: Factor Graph vs EKF

**Dellaert**: "Use the factor graph. It's the right architecture. Re-linearization matters."

**Dampf**: "The EKF works fine for 99% of use cases. I shipped it for 15 years at NovAtel. Don't over-engineer."

**Resolution — Hybrid strategy**:

The factor graph is the offline/post-processing path. The EKF is the real-time path. They share the same measurement pipeline and the same variable definitions. The difference is how state propagates across epochs.

```
                    ┌─────────────┐
                    │ Measurement │
                    │  Pipeline   │  ← Shared by both paths
                    └──────┬──────┘
                           │
              ┌────────────┴────────────┐
              │                         │
     ┌────────▼────────┐      ┌────────▼────────┐
     │  Real-time EKF  │      │  Batch Factor   │
     │  (low latency)  │      │  Graph          │
     │                 │      │  (high accuracy) │
     │  No SPP reset   │      │  Full window     │
     │  after epoch 1  │      │  optimization    │
     │  Cov carries     │      │  Re-linearize    │
     │  forward         │      │  Marginalize     │
     └────────┬────────┘      └────────┬────────┘
              │                         │
              └────────────┬────────────┘
                           │
                    ┌──────▼──────┐
                    │   Solution  │
                    └─────────────┘
```

For real-time: initialize position from SPP at epoch 1 only. After that, carry position forward via IMU mechanization (or constant-velocity if no IMU). The EKF predicts, updates, and passes covariance to the next epoch. No reset. This single change — removing the per-epoch SPP reset — is the highest-ROI improvement in the entire codebase.

For post-processing: run the sliding-window factor graph over the full trajectory. This gives the best possible accuracy by re-linearizing within a 10-30 epoch window. The EKF solution provides the initial guess.

### Debate 2: GNSS-only vs INS-first

**Humphreys**: "INS is the primary sensor. GNSS corrects INS drift. That's how commercial systems work."

**Dampf**: "Disagree. GNSS is the primary sensor for positioning. INS provides smoothing and outage bridging. Most users don't have a tactical-grade IMU — they have a $10 MEMS chip that drifts 100m in 10 seconds without GNSS."

**Resolution — The architecture supports both, but GNSS is the baseline.**

If `imu_config` is `None`, the motion model between epochs is constant-velocity with process noise proportional to the expected dynamics. This works for static receivers and slow-moving platforms.

If `imu_config` is `Some`, the motion model is IMU preintegration between epochs. This provides the tightest possible constraint and enables INS coasting during GNSS outages.

The INS path requires IMU data. The GNSS-only path doesn't. This is reflected in the type system — `RtkInsConfig` requires `ImuConfig`, `RtkConfig` doesn't.

### Debate 3: Precision vs Robustness

**Bisnath**: "The market for centimeter PPP is surveyors who will pay $5K for a software license. That's where the revenue is."

**Humphreys**: "The market for robust urban positioning is 100× larger. Every autonomous vehicle, every drone, every smartphone needs this. Revenue per user is lower but volume is enormous."

**Resolution — Both, via the same architecture.**

The measurement pipeline is identical regardless of target accuracy. The differences are in configuration:
- Survey mode: precise products enabled, AR enabled, static dynamics, long convergence
- Automotive mode: broadcast products, AR disabled (or partial), automotive dynamics, IMU enabled
- Smartphone mode: single-frequency, elevation mask raised, robust estimation, C/N0 weighting

The solver doesn't care. It minimizes the same cost function with different weights and different corrections.

---

## Bespoke Novel Technology

Three inventions emerged from the whiteboard sessions that don't exist in any current GNSS engine:

### 1. Retroactive Cycle Slip Repair

**Dellaert's insight**: "In a sliding-window factor graph, you can detect and repair cycle slips after the fact. If epoch 5 shows a sudden 1-cycle jump in carrier phase on G01, the factor graph can test both hypotheses — with and without the slip — and pick the one with lower residual. The EKF can't do this because it commits to a decision at the moment of the slip."

```rust
/// A cycle-slip hypothesis factor.  Tests whether adding an integer offset
/// to the ambiguity at epoch k produces a better fit.
struct CycleSlipRepairFactor {
    satellite: SatelliteId,
    frequency: u8,
    slip_epoch: u32,     // epoch where slip may have occurred
    test_offset: i32,    // +1, -1, +2, -2 cycles
}

impl Factor for CycleSlipRepairFactor {
    fn variables(&self) -> &[VariableKey] {
        // Connects to the ambiguity variable and all CP measurements
        // after the slip epoch.  The offset is applied as a constant
        // in the residual function.
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        // residual = (CP_meas - geometric_range) / lambda - (N + offset)
        // If offset is correct, residual ≈ 0.
        // If offset is wrong, residual ≈ 1 cycle.
    }
}
```

The optimizer tests each hypothesis (offset = 0, ±1, ±2). The one with the lowest chi-squared residual wins. This is possible because the factor graph can evaluate counterfactuals without committing to them. The EKF, by contrast, must decide at detection time whether a slip occurred.

### 2. Learned Measurement Covariance from C/N0 and Elevation

**Humphreys and Dellaert jointly**: "The variance model `σ² = a² + b² / 10^(C/N0/10)` is a parametric approximation. You can do better by learning the actual variance from residuals."

```rust
/// A measurement variance model learned from historical residuals.
/// Maps (C/N0, elevation, constellation) → expected variance.
struct LearnedVarianceModel {
    // Binned by C/N0 (5 dB-Hz bins), elevation (10° bins), constellation.
    // Each bin stores the empirical variance of historical residuals.
    bins: HashMap<(u8, u8, Constellation), RunningVariance>,
}

impl LearnedVarianceModel {
    fn update(&mut self, cn0: f64, el: f64, constl: Constellation, residual: f64) {
        let key = (cn0_bin(cn0), el_bin(el), constl);
        self.bins.entry(key).or_default().push(residual);
    }

    fn variance(&self, cn0: f64, el: f64, constl: Constellation) -> f64 {
        let key = (cn0_bin(cn0), el_bin(el), constl);
        self.bins.get(&key)
            .map(|v| v.variance())
            .unwrap_or(0.5) // fallback
    }
}
```

This adapts to the actual receiver and environment. A u-blox F9P in urban canyon will have different variance characteristics than a Septentrio in open sky. The model learns these differences from data. No tuning required.

### 3. Topology-Aware Reference Satellite Selection

**Dampf and Teunissen jointly**: "When the reference satellite changes, you have to reset all DD ambiguities. Current practice is to pick the highest-elevation satellite as reference. But in a sliding-window factor graph, you can evaluate the impact of reference changes and minimize ambiguity resets."

The insight: un-differenced measurements are the native representation. DD is a linear transformation applied for convenience. In a factor graph, you can keep un-differenced measurements internally and apply the DD transformation as part of the factor construction. When the reference satellite changes, you don't need to reset anything — you just apply a different linear transformation to the same underlying un-differenced ambiguities.

```rust
// Instead of storing DD ambiguities (fragile, reference-dependent):
//   N_dd = N_rov_sat - N_rov_ref
//
// Store un-differenced ambiguities with a gauge constraint:
//   N_ud[sat] = integer for each satellite
//   Σ N_ud[sat] = 0  (gauge constraint — fixes the reference)

struct UdAmbiguityFactor {
    satellite: SatelliteId,
    frequency: u8,
}

impl Factor for UdAmbiguityFactor {
    fn variables(&self) -> &[VariableKey] {
        // Connects to N_ud[self.satellite] and clock variables.
        // No reference satellite needed — the gauge constraint
        // handles the datum.
    }
}
```

When the reference satellite changes, nothing resets. The gauge constraint shifts from one satellite to another, but the individual satellite ambiguities maintain their integer nature.

---

## The Sprint Plan

### Phase 1: Remove the Ceiling (2 sessions)

**Goal**: Eliminate the per-epoch SPP reset. Single change, maximum impact.

1. Add `initialized: bool` to the EKF state. After epoch 1's SPP initialization, set `initialized = true`.

2. In `predict_state`, when `initialized && !state.is_reset`, skip the SPP position prior. Position carries forward via the state transition matrix (IMU mechanization or constant-velocity).

3. In `check_covariance_divergence`, keep the safety net: if position variance exceeds 10000 m², reset. But the reset threshold should degrade gracefully — first try removing the SPP prior, then fall back to SPP if the filter truly diverges.

4. Wire `process_noise_amb_float = 1e-7` (already the default) through the EKF without the per-epoch ambiguity inflation from the SPP reset.

5. Benchmark: Odaiba RTK with and without SPP reset. Expected: float p50 drops from 1.5m to <0.5m within 30 epochs.

```rust
// In processor/mod.rs or predictor.rs — the key change:
fn predict_state(&mut self, dt: f64) {
    let state = self.current_state.as_mut().unwrap();

    if state.initialized && !state.is_reset {
        // Carry state forward via dynamics model.  Position,
        // velocity, attitude, biases all persist.  Only clock
        // bias is white-noise (reset each epoch).
        self.propagate_dynamics(state, dt);
    } else {
        // First epoch or after reset: initialize from SPP.
        // This path runs exactly once per session (or after
        // a detected divergence).
        self.initialize_from_spp(state, dt);
        state.initialized = true;
    }
}
```

### Phase 2: The Variable-Based Estimator (3 sessions)

**Goal**: Replace the fixed index-based state vector with a variable-based factor graph. Reuse the existing `estimators/factor_graph` optimizer.

1. Define `VariableKey`, `VariableValues`, and `Factor` trait as specified in Dellaert's design above.

2. Implement `SlidingWindowOptimizer` with:
   - `push_epoch()` — adds new GNSS + IMU factors for the current epoch
   - `solve()` — Levenberg-Marquardt over the active variable set
   - `marginalize_oldest()` — Schur complement to remove epoch 0
   - `extract_current_state()` — returns `Solution` for the newest epoch

3. Port existing factors (`gnss_factors.rs`, `imu_factors.rs`) to the new `Factor` trait. The math is already implemented — the change is the interface.

4. Implement the un-differenced ambiguity representation (topology-aware reference satellite selection). Port the LAMBDA decorrelation and search to work with the new ambiguity variables.

```rust
struct SlidingWindowOptimizer {
    window: VecDeque<EpochData>,
    variables: BTreeMap<VariableKey, VariableState>,
    factors: Vec<Box<dyn Factor>>,
    max_window_size: usize,
    marginalized_prior: Option<MarginalizedPrior>,
}

struct MarginalizedPrior {
    /// Linear prior: J^T J * x = J^T r from marginalized epochs.
    /// Encodes all information from epochs no longer in the window.
    hessian: DMatrix<f64>,
    gradient: DVector<f64>,
}
```

### Phase 3: Measurement Pipeline Unification (2 sessions)

**Goal**: One measurement pipeline for all modes. Stacked corrections instead of branched code.

1. Define the correction stack as a pipeline of functions:

```rust
type CorrectionFn = Box<dyn Fn(&EpochObs, &SatState, &RxState) -> f64>;

struct CorrectionStack {
    layers: Vec<(String, CorrectionFn)>,
}

impl CorrectionStack {
    fn ppp_mode(precise: &PreciseProductConfig) -> Self {
        let mut stack = Self::new();
        stack.push("sat_clock",   precise_clock_correction);
        stack.push("sat_orbit",   precise_orbit_correction);
        stack.push("sat_antenna", ant ex_pcv_correction);
        stack.push("rcv_antenna", ant ex_pcv_correction);
        stack.push("tide",        solid_earth_tide);
        stack.push("phase_windup", phase_windup_correction);
        stack
    }

    fn rtk_mode(base_obs: &EpochObs, base_pos: Vector3<f64>) -> Self {
        let mut stack = Self::new();
        stack.push("broadcast_clock", broadcast_clock_correction);
        stack.push("broadcast_orbit", broadcast_orbit_correction);
        // DD handles the rest — troposphere and ionosphere cancel
        // on short baselines.
        stack
    }
}
```

2. The measurement builder is mode-agnostic. It takes a `CorrectionStack` and produces factors:

```rust
fn build_gnss_factors(
    obs: &EpochObs,
    corrections: &CorrectionStack,
    state: &VariableValues,
) -> Vec<Box<dyn Factor>> {
    let mut factors = Vec::new();
    for sat_obs in &obs.satellites {
        let sat_state = compute_sat_state(sat_obs, corrections);
        if let Some(pr_factor) = PseudorangeFactor::new(sat_obs, &sat_state) {
            factors.push(Box::new(pr_factor));
        }
        if let Some(cp_factor) = CarrierPhaseFactor::new(sat_obs, &sat_state) {
            factors.push(Box::new(cp_factor));
        }
        if let Some(dop_factor) = DopplerFactor::new(sat_obs, &sat_state) {
            factors.push(Box::new(dop_factor));
        }
    }
    factors
}
```

3. Delete the separate PPP and RTK measurement model code. Replace with the unified pipeline. This removes ~8000 lines of duplicate code.

### Phase 4: INS Integration (2 sessions)

**Goal**: First-class IMU support with preintegration factors.

1. Implement `ImuPreintegrationFactor` using the existing IMU mechanization from `predictor.rs`:

```rust
struct ImuPreintegrationFactor {
    epoch_from: u32,
    epoch_to: u32,
    delta_p: Vector3<f64>,  // position change (m)
    delta_v: Vector3<f64>,  // velocity change (m/s)
    delta_q: UnitQuaternion<f64>,  // attitude change
    cov: Matrix9<f64>,      // 9×9 preintegration covariance
    jacobian_wrt_bias: Matrix9x6<f64>,  // bias Jacobians
}

impl Factor for ImuPreintegrationFactor {
    fn variables(&self) -> &[VariableKey] {
        &[
            VariableKey::Position(self.epoch_from),
            VariableKey::Velocity(self.epoch_from),
            VariableKey::Attitude(self.epoch_from),
            VariableKey::Position(self.epoch_to),
            VariableKey::Velocity(self.epoch_to),
            VariableKey::Attitude(self.epoch_to),
            VariableKey::AccelBias,
            VariableKey::GyroBias,
        ]
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        // residual = predicted_motion - preintegrated_delta
        // Where predicted_motion = f(pos_k, vel_k, att_k, pos_{k+1}, vel_{k+1}, att_{k+1})
        // Accounting for bias corrections via the Jacobian.
    }
}
```

2. Implement NHC as a factor (not a separate update):

```rust
struct NhcFactor {
    epoch: u32,
    sigma_lateral: f64,
    sigma_vertical: f64,
}

impl Factor for NhcFactor {
    fn variables(&self) -> &[VariableKey] {
        &[VariableKey::Velocity(self.epoch), VariableKey::Attitude(self.epoch)]
    }

    fn residual(&self, values: &VariableValues) -> DVector<f64> {
        // residual = R_body_to_ecef * [v_x, v_y, v_z] - [0, 0, 0]
        // Penalizes non-zero lateral and vertical velocity in body frame.
    }
}
```

3. Benchmark on drive_imu dataset. Expected: RTK-INS with NHC should maintain <1m accuracy through 5-second GNSS outages.

### Phase 5: Benchmark Suite (1 session)

**Goal**: Continuous accuracy validation on reference datasets.

1. Create `crates/gneiss-bench/` with benchmark harness:

```rust
#[derive(Debug, Serialize)]
struct BenchmarkResult {
    name: String,
    mode: String,
    n_epochs: usize,
    horizontal_rms: f64,
    vertical_rms: f64,
    horizontal_p50: f64,
    horizontal_p95: f64,
    fix_rate: f64,
    time_to_first_fix_s: Option<f64>,
}

fn run_benchmark(dataset: &Dataset, config: &EngineConfig) -> BenchmarkResult {
    let mut engine = Engine::new(config.clone())?;
    let mut errors = Vec::new();

    for epoch in &dataset.epochs {
        let solution = engine.process_epoch(
            epoch.rover.clone(),
            epoch.base.clone(),
            epoch.imu.as_deref(),
        )?;

        if let Some(ref truth) = epoch.ground_truth {
            let err = (solution.position.ecef - truth.position).norm();
            errors.push(err);
        }
    }

    BenchmarkResult {
        horizontal_rms: rms(&errors),
        // ... compute all metrics
    }
}
```

2. Target datasets:
   - **Odaiba 4km** — RTK, short baseline, urban
   - **F9P 8km** — RTK, medium baseline, suburban
   - **GICI-LIB** — PPP-AR, open sky + urban, multi-GNSS
   - **drive_imu** — RTK-INS, highway
   - **WHU-Smartphone** — smartphone-grade, urban

3. CI integration: `cargo bench --workspace` runs all benchmarks. Results compared against stored baselines. Accuracy regression >5% fails CI.

### Phase 6: Hardening (1 session)

**Goal**: Production-grade error recovery and diagnostics.

1. Graceful degradation hierarchy:

```rust
fn process_epoch_fallible(&mut self, ...) -> Result<Solution, EngineError> {
    // Try RTK fixed
    if let Ok(sol) = self.try_rtk_fixed(rover, base, imu) {
        return Ok(sol);
    }
    // Try RTK float
    if let Ok(sol) = self.try_rtk_float(rover, base, imu) {
        return Ok(sol);
    }
    // Try SPP
    if let Ok(sol) = self.try_spp(rover) {
        return Ok(sol);
    }
    // Coast on IMU
    if let Some(imu_data) = imu {
        if let Ok(sol) = self.try_imu_coast(imu_data) {
            return Ok(sol);
        }
    }
    // Nothing works — return last known position with inflated covariance
    Ok(self.last_known_solution().with_inflated_covariance())
}
```

2. Every `unwrap()` removed. Every `expect()` has a documented invariant. Every matrix inverse has a fallback. Every NaN is caught and handled.

---

## Total Scope

| Phase | Sessions | LOC Impact |
|-------|----------|------------|
| 1. Remove SPP reset | 2 | +200, -50 |
| 2. Variable-based estimator | 3 | +3000, -2000 |
| 3. Unified measurement pipeline | 2 | +1500, -8000 |
| 4. INS integration | 2 | +1500, -500 |
| 5. Benchmark suite | 1 | +1000, -0 |
| 6. Hardening | 1 | +500, -200 |
| **Total** | **11 sessions** | **+7700, -10750** |

Net result: a ~3000-line reduction in codebase size with dramatically more capability. The variable-based estimator replaces the divergent PPP/RTK/SPP paths with a single code path. Every feature (IONEX, IFB, PR validation, robust estimation) benefits all modes automatically.
