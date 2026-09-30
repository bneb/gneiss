# Project: Gneiss Frame Safety & Epoch Alignment Refactoring

## Architecture
Zero-cost compile-time typestates in `gneiss-core` enforcing spatial frame safety (ECEF, NED, ENU, Body FRD), geodetic datum consistency (WGS84, ITRF2014, ITRF2020, NAD83, JGD2011, PZ-90), temporal epoch alignment (`Epoch<Scale>`, `TimeDelta`, `EpochKey`), relational coupling on shared geometry (`DoubleDiffGeometry<F>`), and refactored estimator states (`EskfState`, `SwfgEngine`, `PostProcessOptions`, `eval_odaiba_ins`).

## Feature Inventory
| # | Feature | Description | Milestone | Source |
|---|---------|-------------|-----------|--------|
| 1 | Spatial typestate markers & wrappers | `Point3<F>`, `SpatialVector<F>`, `SpatialVelocity<F>`, `SpatialCovariance<F>`, `NedCovariance` | M1 | Survey 1 |
| 2 | Reference frame realizations | Add `Pz90` to realizations (Helmert aligned to ITRF2014) alongside existing datums | M1 | Survey 1 |
| 3 | Antenna lever arm & attitude typestates | `AntennaLeverArm(SpatialVector<BodyFrd>)`, `Attitude<From, To>`, prohibiting unrotated addition | M1 | Survey 1 |
| 4 | Non-leaky typestate deref & operators | Remove leaky `Deref<Target=Vector3<f64>>`; provide `.vector()` / `.coords()` and compile-fail cross-frame ops | M1 | Survey 1 |
| 5 | Relational local tangent plane | `LocalTangentPlane<F>` constructed from single `EcefPos<F>` replacing parallel bare float calls | M1 | Survey 1 |
| 6 | Temporal typestate scale markers | `Epoch<Scale: TimeScale>` with `GpsScale`, `BdtScale`, `GstScale`, `GlonassScale`, `UtcScale` | M2 | Survey 2 |
| 7 | Integer nanosecond time arithmetic | `(week: u32, tow_nanos: u64)`, `TimeDelta { nanos: i64 }`, week-rollover safe math | M2 | Survey 2 |
| 8 | BeiDou broadcast ephemeris scale fix | Align `toe` and `toc` time systems to eliminate 14-second satellite clock bias error | M2 | Survey 2 |
| 9 | Explicit leap second conversions | Explicit leap second handling in `UtcScale` / `GpsScale`, fixing `antex.rs` and `gneiss-fetch` | M2 | Survey 2 |
| 10 | Structured epoch alignment & key | `EpochKey` and tolerance-based matching `is_within(tolerance)` replacing ad-hoc float rounding | M2 | Survey 2 |
| 11 | Relational double-difference geometry | `DoubleDiffGeometry<F>` derived from single ephemeris and station positions | M3 | Survey 3 |
| 12 | Coordinate fix in receiver PCV | Fix `receiver_pcv::dd_correction_from_geometry` passing `rov_llh` as `pos_ecef` | M3 | Survey 3 |
| 13 | Relational receiver antenna coupling | Compute base antenna PCV from base elevation instead of borrowing rover elevation | M3 | Survey 3 |
| 14 | Typed ESKF estimator state | `EskfState` using `EcefPos<F>`, `SpatialVelocity<Ecef<F>>`, `BodyFrdVector` for biases | M4 | Survey 3 |
| 15 | Lever arm rotation enforcement | Enforce $C_b^e$ rotation in `build_gnss_pos_system`, `build_doppler_velocity_system`, `dd_update` | M4 | Survey 3 |
| 16 | SwfgEngine & PostProcessOptions refactoring | Migrate initial/prev positions to `EcefPos<F>` and lever arms to `AntennaLeverArm` | M4 | Survey 3 |
| 17 | eval_odaiba_ins refactoring | Migrate Odaiba INS pipeline to typed primitives while remaining strictly < 500 LOC | M4 | Survey 3 |
| 18 | Workspace test suite verification | All 789 existing workspace tests pass with 0 failures | M5 | Survey 1/2/3 |
| 19 | Dual CI smoke benchmark verification | `check_network_benchmark.py --smoke` & `check_multignss_benchmark.py --smoke` pass | M5 | Survey 3 |
| 20 | eval_odaiba_ins benchmark verification | RTS smoothed horizontal $p_{50} \le 1.80$m, $\text{RMS} \le 3.50$m, 0 false fixes | M5 | Survey 3 |
| 21 | AGENTS.md code standards audit | < 500 LOC/file, <= 32 LOC/fn, nesting < 3, 0 warnings under `-D warnings`, 0 unwrap() | M5 | Survey 1/2/3 |

## Milestones
| # | Name | Scope | Dependencies | Status |
|---|------|-------|-------------|--------|
| 1 | M1: Spatial Primitives & Frames | Features 1–5 in `gneiss-core::frames` & `coords` | none | IN_PROGRESS |
| 2 | M2: Temporal Frames & Epoch Safety | Features 6–10 in `gneiss-core::time`, `gneiss-parsers`, `gneiss-fetch` | none | IN_PROGRESS |
| 3 | M3: Relational Coupling & Geometry | Features 11–13 in `gneiss-core::coords`, `gneiss-parsers::receiver_pcv`, `receiver_antenna`, `gneiss-rtk::formation` | M1 | PLANNED |
| 4 | M4: Estimator State & Pipelines | Features 14–17 in `gneiss-rtk::estimators::eskf`, `SwfgEngine`, `PostProcessOptions`, `eval_odaiba_ins` | M1, M2, M3 | PLANNED |
| 5 | M5: Verification & Benchmarks | Features 18–21 workspace tests, CI smoke guards, Odaiba INS benchmark, AGENTS.md audit | M1, M2, M3, M4 | PLANNED |

## Interface Contracts
### `gneiss-core::frames`
- `pub struct SpatialVector<F: CoordinateFrame>(Vector3<f64>, PhantomData<F>);`
- `pub struct SpatialVelocity<F: CoordinateFrame>(Vector3<f64>, PhantomData<F>);`
- `pub struct SpatialCovariance<F: CoordinateFrame>(Matrix3<f64>, PhantomData<F>);`
- `pub struct Point3<F: CoordinateFrame>(Vector3<f64>, PhantomData<F>);`
- `pub type EcefPos<R> = Point3<Ecef<R>>;`
- `pub struct AntennaLeverArm(SpatialVector<BodyFrd>);`
- `pub struct Attitude<From: CoordinateFrame, To: CoordinateFrame>(UnitQuaternion<f64>, PhantomData<(From, To)>);`
- `impl<From: CoordinateFrame, To: CoordinateFrame> Attitude<From, To> { pub fn rotate(&self, v: &SpatialVector<From>) -> SpatialVector<To>; }`
- Addition between `Point3<F>` and `SpatialVector<F>` yields `Point3<F>`. Addition between mismatched frames does not compile.

### `gneiss-core::time`
- `pub trait TimeScale: Copy + Clone + 'static {}`
- `pub struct GpsScale; pub struct BdtScale; pub struct GstScale; pub struct GlonassScale; pub struct UtcScale;`
- `pub struct Epoch<Scale: TimeScale> { week: u32, tow_nanos: u64, _scale: PhantomData<Scale> }`
- `pub struct TimeDelta { nanos: i64 }`
- `Epoch<S> - Epoch<S> -> TimeDelta`. Cross-scale subtraction fails to compile without explicit conversion.

### `gneiss-core::coords` / `gneiss-rtk::formation`
- `pub struct DoubleDiffGeometry<F: ReferenceFrame>`
  Constructed from `(rov_pos: EcefPos<F>, bas_pos: EcefPos<F>, sat_pos: EcefPos<F>, ref_sat_pos: EcefPos<F>)`.
  Exposes: `rov_los_unit: Vector3<f64>`, `bas_los_unit: Vector3<f64>`, `rov_el_rad: f64`, `bas_el_rad: f64`, `rov_az_rad: f64`, `bas_az_rad: f64`, `dd_los_unit: Vector3<f64>`.

## Code Layout
- `crates/gneiss-core/src/frames/`:
  - `mod.rs`: re-exports
  - `markers.rs`: `CoordinateFrame`, `Ecef<R>`, `Ned`, `Enu`, `BodyFrd`
  - `realizations.rs`: reference frame realizations including `Pz90`
  - `primitives.rs`: `SpatialVector`, `SpatialVelocity`, `SpatialCovariance`, `Point3`, `AntennaLeverArm`, `Attitude`
  - `tangent.rs`: `LocalTangentPlane<F>`
- `crates/gneiss-core/src/time/`:
  - `epoch.rs`: `Epoch<Scale>`, `TimeDelta`, `EpochKey`
  - `scales.rs`: `GpsScale`, `BdtScale`, `GstScale`, `GlonassScale`, `UtcScale`
- `crates/gneiss-core/src/coords/`:
  - `geometry.rs`: `DoubleDiffGeometry<F>`
- `crates/gneiss-rtk/src/estimators/eskf/`:
  - `types.rs`: typed `EskfState`
  - `update.rs`, `dd_update.rs`: typed measurement updates
- `crates/gneiss-rtk/src/bin/eval_odaiba_ins/`:
  - `main.rs`: refactored evaluation pipeline (< 500 LOC)
