# Project: Gneiss GNSS/INS Tier-1 Commercial Frontiers

## Architecture
The Gneiss navigation engine is being expanded to commercial Tier-1 parity across three frontiers and composite modes:
1. **Frontier R1 (Inertial Navigation)**: A 15-state Error-State Kalman Filter (ESKF/MEKF) in `crates/gneiss-rtk/src/estimators/eskf/` with states ($\delta \mathbf{p}^e, \delta \mathbf{v}^e, \delta \boldsymbol{\theta}, \delta \mathbf{b}_a, \delta \mathbf{b}_g$), error-quaternion feedback ($\mathbf{q} \leftarrow \mathbf{q} \otimes \delta\mathbf{q}$), closed-loop online accelerometer and gyro bias estimation, coupled Non-Holonomic Constraints (NHC) and ZUPT, and a full 15-state backward Rauch-Tung-Striebel (RTS) smoother over forward trajectory history.
2. **Frontier R2 (Precise Point Positioning)**: Integer PPP-AR engine in `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs` and `ar_handler.rs`, powered by fast indexed SINEX OSB/BIA phase/code bias ingestion in `crates/gneiss-parsers/src/sinex_bia.rs`, satellite nadir PCV and standalone receiver PCO/PCV, multi-constellation support (GPS, Galileo, BeiDou, QZSS), and single-differenced Wide-Lane (MW) and Narrow-Lane LAMBDA integer search.
3. **Frontier R3 (Network RTK / VRS)**: Multi-station regional CORS atmospheric engine in `crates/gneiss-rtk/src/post_process/` (`network_adj.rs`, `vrs.rs`) and `crates/gneiss-rtk/src/spatial/delaunay.rs`, performing multi-baseline double-difference integer ambiguity network adjustment, spatial 2D Delaunay triangulation for station troposphere ZWD and satellite ionospheric pierce points (IPP), and synthesizing localized Virtual Reference Station (VRS) observations (< 1 km effective baseline).
4. **Frontier R4 (Unified Composite Integration)**: Modular composite interfaces in `crates/gneiss-rtk/src/composite/` (`tc_ppp.rs` for Tightly-Coupled PPP/INS, `tc_rtk.rs` for Tightly-Coupled Network RTK/INS) coupling un-differenced/double-difference carrier-phase innovations with the 15-state ESKF.
5. **E2E Testing Track**: Parallel requirement-driven, opaque-box test suite across Tiers 1-4 testing all features independently of internal implementation.

## Feature Inventory
| # | Feature | Description | Milestone | Source |
|---|---------|-------------|-----------|--------|
| 1 | 15-State ESKF Formulation | States ($\delta\mathbf{p}^e, \delta\mathbf{v}^e, \delta\boldsymbol{\theta}, \delta\mathbf{b}_a, \delta\mathbf{b}_g$), transition $\boldsymbol{\Phi}_{15\times 15}$ with strictly positive `vel_att = +f_e_skew * dt`, spectral $\mathbf{Q}_{15\times 15}$ | M1 | Survey R1 |
| 2 | Error-Quaternion Feedback | Multiplicative error quaternion reset $\mathbf{q} \leftarrow \mathbf{q} \otimes \delta\mathbf{q}$ after measurement updates | M1 | Survey R1 |
| 3 | Online Bias Estimation | Closed-loop accelerometer ($\delta\mathbf{b}_a$) and gyroscope ($\delta\mathbf{b}_g$) bias updates driven by GNSS position & velocity innovations | M1 | Survey R1 |
| 4 | 15-State RTS Smoother | Backward Rauch-Tung-Striebel smoothing over forward filter history propagating corrections to attitude and biases | M1 | Survey R1 |
| 5 | Coupled NHC & ZUPT | Non-holonomic lateral/vertical constraints with attitude Jacobian $\mathbf{E}_{23}(\mathbf{R}_b^e)^T[\mathbf{v}^e\times]$ and zero-velocity updates | M1 | Survey R1 |
| 6 | Odaiba Benchmark Target | `eval_odaiba_ins` achieves $p_{50} < 2.5\text{ m}$ and $\text{RMS} < 5.2\text{ m}$ across 12,398 epochs | M1 | Survey R1 |
| 7 | Fast SINEX OSB Ingestion | Indexed `SinexBias` with $O(1)$ lookups and Wide-Lane/Narrow-Lane satellite bias helpers in `sinex_bia.rs` | M2 | Survey R2 |
| 8 | Antenna PCO/PCV Corrections | Satellite nadir PCV interpolation and standalone receiver PCO/PCV for un-differenced PPP | M2 | Survey R2 |
| 9 | Multi-Constellation PPP | Support GPS, Galileo, BeiDou, and QZSS across observation processing, ephemeris, and band selection | M2 | Survey R2 |
| 10 | Single-Differenced LAMBDA AR | Single-differenced Wide-Lane (MW) rounding and Narrow-Lane LAMBDA integer search decoupled from receiver phase biases | M2 | Survey R2 |
| 11 | Kinematic PPP-AR Benchmark | `eval_ppp` resolves integer ambiguities on F9P drive to sub-meter kinematic accuracy vs CSRS-PPP ($0.296\text{ m}$ RMS) | M2 | Survey R2 |
| 12 | Multi-Station CORS Ingestion | Ingest 5–10 regional CORS base station observation streams simultaneously from `datasets/cors_sf_bay_network` | M3 | Survey R3 |
| 13 | Network Baseline Adjustment | Formulate multi-baseline double-difference network adjustment to solve for integer ambiguities across CORS network | M3 | Survey R3 |
| 14 | Delaunay Atmospheric Models | Spatial 2D Delaunay triangulation for iono pierce points (IPP) and troposphere ZWD gradients | M3 | Survey R3 |
| 15 | Localized VRS Synthesis | Synthesize localized Virtual Reference Station (VRS) observations at rover's position (< 1 km effective baseline) | M3 | Survey R3 |
| 16 | Network RTK Benchmark | `eval_network_ppk` reduces baseline ppm error across CORS baselines (P181, P222, P225) toward Leica $8\text{ mm} + 1\text{ ppm}$ | M3 | Survey R3 |
| 17 | Tightly-Coupled PPP/INS | Modular composite interface `tc_ppp.rs` coupling un-differenced PPP carrier phase/pseudorange with 15-state ESKF | M4 | Survey R4 |
| 18 | Tightly-Coupled Network RTK/INS | Modular composite interface `tc_rtk.rs` coupling VRS double-difference carrier phase with 15-state ESKF | M4 | Survey R4 |
| 19 | Composite Execution Tests | Integration verification and benchmarks confirming successful execution of TC-PPP and TC-RTK pipelines | M4 | Survey R4 |
| 20 | E2E Test Suite (Tiers 1-4) | Opaque-box requirement-driven tests covering all 19 features with $\ge 5$ cases per feature, BVA, pairwise, and application workloads | E2E Track | ORIGINAL_REQUEST |
| 21 | Adversarial Hardening (Tier 5) | White-box stress-testing, edge case generation, and zero surviving mutant verification | M5 | Project Pattern |

## Milestones
| # | Name | Scope | Dependencies | Status |
|---|------|-------|-------------|--------|
| E2E | E2E Testing Suite | Requirements-driven opaque-box test suite (Tiers 1-4) + TEST_READY.md | none | DONE |
| M1 | 15-State ESKF GNSS/INS | 15-state ESKF, error-quaternion reset, online bias estimation, 15-state RTS smoother, coupled NHC/ZUPT, eval_odaiba_ins benchmark | none | DONE |
| M2 | Integer PPP-AR Engine | SINEX OSB parser indexing, nadir satellite PCV, standalone receiver PCO/PCV, multi-constellation, SD WL/NL LAMBDA AR, eval_ppp benchmark | none | IN_PROGRESS |
| M3 | Network RTK VRS Engine | Multi-CORS ingestion, multi-baseline DD network adjustment, Delaunay iono/tropo models, VRS synthesis, eval_network_ppk benchmark | none | DONE |
| M4 | Unified Composite Integration | Tightly-Coupled PPP/INS (tc_ppp.rs) and Tightly-Coupled Network RTK/INS (tc_rtk.rs) composite pipelines | M1, M2, M3 | DONE |
| M5 | Final Milestone: E2E & Tier 5 | Pass 100% of E2E test suite (Tiers 1-4) and complete Adversarial Coverage Hardening (Tier 5) | M4, E2E | PLANNED |

## Interface Contracts

### 1. 15-State ESKF Module (`crates/gneiss-rtk/src/estimators/eskf/`)
- `EskfState`:
  ```rust
  pub struct EskfState {
      pub pos_ecef: Vector3<f64>,
      pub vel_ecef: Vector3<f64>,
      pub attitude: UnitQuaternion<f64>,
      pub accel_bias: Vector3<f64>,
      pub gyro_bias: Vector3<f64>,
      pub cov: Matrix15<f64>,
  }
  ```
- Predict step:
  ```rust
  pub fn predict(state: &mut EskfState, imu: &ImuMeasurement, dt: f64, q_diag: &Vector15<f64>) -> Result<(), EngineError>;
  ```
- Measurement updates:
  ```rust
  pub fn update_gnss_pos_vel(state: &mut EskfState, pos_meas: &Vector3<f64>, vel_meas: &Vector3<f64>, lever_arm: &Vector3<f64>, r_pos: &Matrix3<f64>, r_vel: &Matrix3<f64>) -> Result<(), EngineError>;
  pub fn update_nhc(state: &mut EskfState, lever_arm: &Vector3<f64>, r_nhc: &Matrix2<f64>) -> Result<(), EngineError>;
  pub fn update_zupt(state: &mut EskfState, r_zupt: &Matrix3<f64>) -> Result<(), EngineError>;
  ```
- Backward RTS smoother:
  ```rust
  pub struct EskfSmoother { ... }
  impl EskfSmoother {
      pub fn smooth(&self) -> Result<Vec<EskfState>, EngineError>;
  }
  ```

### 2. SINEX OSB & PPP-AR (`crates/gneiss-parsers/src/sinex_bia.rs`, `crates/gneiss-rtk/src/ambiguity/ppp_ar.rs`)
- `SinexBias`:
  ```rust
  impl SinexBias {
      pub fn lookup_bias(&self, sat: &Satellite, obs_code: &str, time: Epoch) -> Option<f64>;
      pub fn wide_lane_satellite_bias(&self, sat: &Satellite, time: Epoch) -> Option<f64>;
      pub fn narrow_lane_satellite_bias(&self, sat: &Satellite, time: Epoch) -> Option<f64>;
  }
  ```
- `PppArSolver`:
  ```rust
  impl PppArSolver {
      pub fn fix_single_diff_ambiguities(
          &self,
          float_ambiguities: &[f64],
          cov: &DMatrix<f64>,
          wavelengths: &[f64],
      ) -> Result<Vec<i32>, EngineError>;
  }
  ```

### 3. Delaunay Triangulation & VRS (`crates/gneiss-rtk/src/spatial/delaunay.rs`, `crates/gneiss-rtk/src/post_process/vrs.rs`)
- `DelaunayTriangulation2D`:
  ```rust
  pub struct Delaunay2D { ... }
  impl Delaunay2D {
      pub fn new(points: &[Vector2<f64>]) -> Result<Self, EngineError>;
      pub fn interpolate_barycentric(&self, query: &Vector2<f64>, values: &[f64]) -> Option<f64>;
  }
  ```
- `VrsSynthesizer`:
  ```rust
  pub struct VrsSynthesizer { ... }
  impl VrsSynthesizer {
      pub fn synthesize_epoch(&self, rover_approx_ecef: &Vector3<f64>, base_epochs: &[StationEpoch]) -> Result<EpochObservation, EngineError>;
  }
  ```

### 4. Unified Composite Modules (`crates/gneiss-rtk/src/composite/`)
- `TightlyCoupledPppIns`:
  ```rust
  pub struct TightlyCoupledPppIns {
      eskf: EskfState,
      ppp_ar: PppArSolver,
  }
  impl TightlyCoupledPppIns {
      pub fn process_epoch(&mut self, imu_samples: &[ImuMeasurement], ppp_obs: &EpochObservation) -> Result<NavSolution, EngineError>;
  }
  ```
- `TightlyCoupledNetworkRtkIns`:
  ```rust
  pub struct TightlyCoupledNetworkRtkIns {
      eskf: EskfState,
      vrs_synth: VrsSynthesizer,
  }
  impl TightlyCoupledNetworkRtkIns {
      pub fn process_epoch(&mut self, imu_samples: &[ImuMeasurement], rover_obs: &EpochObservation, cors_obs: &[StationEpoch]) -> Result<NavSolution, EngineError>;
  }
  ```

## Code Layout
Strict adherence to `AGENTS.md` (< 500 LOC/file, < 32 LOC/function, < 3 nesting depth, 0 unwrap, 0 warnings):
```
crates/gneiss-rtk/src/
├── estimators/
│   └── eskf/
│       ├── mod.rs          # Module declarations and re-exports (< 100 LOC)
│       ├── types.rs        # EskfState, error covariance, ImuMeasurement (< 250 LOC)
│       ├── predict.rs      # Mechanization and covariance propagation Phi_15x15, Q_15x15 (< 300 LOC)
│       ├── update.rs       # GNSS position & velocity updates, error-quaternion reset (< 300 LOC)
│       ├── constraints.rs   # Coupled NHC with attitude Jacobian and ZUPT (< 250 LOC)
│       └── smoother.rs     # 15-state backward RTS smoother (< 350 LOC)
├── spatial/
│   ├── mod.rs              # Spatial algorithms module (< 50 LOC)
│   └── delaunay.rs         # 2D Delaunay triangulation & barycentric interpolation (< 350 LOC)
├── post_process/
│   ├── network_adj.rs      # Multi-baseline DD CORS network adjustment (< 400 LOC)
│   └── vrs.rs              # Delaunay-enhanced VRS observation synthesis (< 400 LOC)
├── composite/
│   ├── mod.rs              # Composite module declarations (< 50 LOC)
│   ├── tc_ppp.rs           # Tightly-coupled PPP/INS pipeline (< 350 LOC)
│   └── tc_rtk.rs           # Tightly-coupled Network RTK/INS pipeline (< 350 LOC)
└── ambiguity/
    └── ppp_ar.rs           # Enhanced Wide-Lane / Narrow-Lane LAMBDA AR (< 400 LOC)
crates/gneiss-parsers/src/
├── sinex_bia.rs            # Fast indexed SINEX OSB/BIA lookup & helpers (< 350 LOC)
└── receiver_pcv/           # Standalone receiver PCO/PCV evaluation (< 300 LOC)
```
