# Survey Report: R2 — Temporal Frame & Epoch Alignment Safety

**Author**: Survey Explorer 2  
**Date**: 2026-09-25  
**Working Directory**: `.agents/teamwork/teamwork_preview_explorer_survey_2`  
**Milestone**: Frame Safety & Epoch Alignment Refactoring (R2 Focus)  

---

## 1. Executive Summary

Gneiss currently relies almost universally on `gneiss_core::time::GpsTime` as an untyped container holding `(week: u32, tow: f64)`. While `GpsTime` was intended to represent GPS System Time (GPST), it is routinely used across the codebase to store calendar dates, GLONASS timestamps (UTC+3h), BeiDou timestamps (BDT), and raw millisecond tags from RTCM and IMU streams without compile-time distinction.

The investigation uncovered critical systemic risks across four major domains:
1. **Silent Cross-Constellation Scale Pollution**: In `gneiss-parsers/src/rinex/nav/builder.rs:144`, BeiDou broadcast ephemeris sets `toe` using raw BDT seconds of week, whereas `toc` has had a 14-second GPS offset applied. In `gneiss-core/src/ephemeris/keplerian.rs:178-181`, `BeidouEphemeris` computes `tk = t_bdt - toe` and `tc = t_bdt - toc`, introducing a 14-second clock-error discrepancy.
2. **Week-Rollover Catastrophic Failures (`tow - tow`)**: Callers frequently bypass `GpsTime::sub(other)` and instead compute `(t1.tow - t2.tow)` directly (e.g. in `rtk_iekf/predict.rs:28`, `tc_rtk.rs:165`, `swfg/engine/epoch.rs:21`). At GPS week rollover (Sunday 00:00:00 GPST), `dt` explodes to `-604799` seconds, which corrupts Kalman filter state prediction and ephemeris selection.
3. **Leap-Second Omissions in External Services**: `GpsTime::from_calendar` calculates Julian Days without applying leap seconds. In `gneiss-parsers/src/antex.rs:194` and `gneiss-fetch/src/sources/{bkg, noaa, cddis}.rs`, GPS time is converted to Unix/UTC timestamps without subtracting the 18 leap seconds ($\Delta t_{LS} = 18\text{ s}$), leading to corrupted PCV lookup windows and incorrect DOY/Year URL queries near midnight.
4. **Fragile Ad-Hoc Epoch Matching**: Across 16 identified sites, epoch matching is implemented through ad-hoc float rounding (`(tow * 1000.0).round() as u64`, `(tow * 10.0).round() as u32`, `tow.round() as u32`, `tow.floor() as u32`). In `eval_odaiba_ins/main.rs:48-49`, base station matching requires an exact integer millisecond match (`base_index.get(&exact_ms)`); any sub-millisecond clock drift or jitter drops the base observation entirely.

This report establishes the complete evidence base and proposes a zero-cost typestate architecture (`Epoch<Scale>`) and nanosecond-precision `TimeDelta` that structurally prevents these bugs at compile time.

---

## 2. Current Time Representations Across Workspace

### 2.1 `gneiss_core::time::GpsTime`
Defined in `crates/gneiss-core/src/time.rs:7-10`:
```rust
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct GpsTime {
    pub week: u32,
    pub tow: f64,
}
```

#### Identified Architectural Flaws:
1. **Public Mutable Fields**: Callers can directly mutate `week` and `tow`, bypassing normalization invariants (`tow ∈ [0.0, 604800.0)`).
2. **Naive Derived Ordering**: `PartialOrd` compares `week` first, then `tow`. If an instance is unnormalized (e.g. `week: 2200, tow: 604805.0` vs `week: 2201, tow: 5.0`), `PartialOrd` reports inequality and incorrect order for identical physical instants.
3. **Calendar Epoch Misconception**: `GpsTime::from_calendar(year, month, day, hour, minute, sec)` computes JD and subtracts Jan 6 1980 00:00:00 JD (`2444244.5`). It assumes the input calendar date is in GPS time, but callers consistently feed it civil UTC calendar dates without applying the 18-second leap-second correction.
4. **Float Jitter in Equality**: `PartialEq` performs exact `f64` comparison on `tow`, causing equality checks (`t1 == t2`) to fail when timestamps originate from different arithmetic paths.

### 2.2 `gneiss_core::gnss_time::{TimeSystem, GnssTime}`
Defined in `crates/gneiss-core/src/gnss_time.rs:44-53, 99-107`:
```rust
pub enum TimeSystem {
    Gps,
    Glonass,
    Bdt,
    Gst,
}

pub struct GnssTime {
    pub sys: TimeSystem,
    pub week: u32,
    pub tow: f64,
}
```

#### Analysis:
- `GnssTime` provides pinned ICD conversion offsets (`TimeSystem::gpst_offset()`), but it is a **runtime-checked enum**, not a compile-time typestate.
- Cross-system arithmetic is not prevented: a user can still freely access `tow` and add it to another timestamp.
- Outside of `crates/gneiss-rtk/src/events.rs` and its own unit tests, `GnssTime` is **almost completely unused** across the codebase. Parsers and estimators bypass it and use `GpsTime` directly.

### 2.3 Completely Missing Time Representations
The following essential time systems have **no dedicated types**:
| System | Required Representation | Current Workaround |
|---|---|---|
| **UTC** | `UtcTime` (civil time with leap seconds) | Bare calendar parsing, chrono `NaiveDateTime`, or untyped `GpsTime` |
| **GLONASS** | `GlonassTime` (UTC(SU)+3h, leap-second stepped) | `GpsTime` with raw `+ (18.0 - 10800.0)` offset |
| **BeiDou** | `BdtTime` (epoch 2006-01-01, GPST - BDT = 14s) | `GpsTime` with raw `- 14.0` offset |
| **Galileo** | `GstTime` (GST, nominal 0s offset vs GPST) | Treated identically as `GpsTime` |
| **Duration** | `TimeDelta` / `Duration` (signed nanoseconds) | Bare float `dt: f64` or `u64` microseconds |

---

## 3. Census of Ad-Hoc Float Rounding & Integer Truncation Sites

The codebase exhibits widespread ad-hoc quantization for matching epochs across rover, base, truth, and IMU streams.

### 3.1 Detailed Inventory of Matching Sites

| File & Line | Code Expression | Matching Target | Failure Mode / Vulnerability |
|---|---|---|---|
| `crates/gneiss-rtk/src/post_process/network.rs:92-94` | `(time.tow.round() as i64).rem_euclid(604800) as u32` | Base trajectory fusion (`by_tow`) | 1 Hz integer rounding; collapses 5Hz/10Hz epochs into key collisions; ignores week rollover |
| `crates/gneiss-rtk/src/post_process/iekf_pass.rs:248-251` | `(b.time.tow - tow).abs() < 0.1` | Base epoch matching | Ignores `week`; fails across GPS week boundaries |
| `crates/gneiss-rtk/src/post_process/combiner.rs:47` | `(fwd.time.tow * 1000.0).round() as u64` | Forward/backward fusion lookup | Exact integer millisecond lookup; fails if passes have sub-ms jitter; ignores week |
| `crates/gneiss-rtk/src/post_process/backward.rs:255-256` | `(epoch.time.tow * 1000.0).round() as u64`<br>`(epoch.time.tow * 1_000_000.0).round() as u64` | IMU grouping by epoch | Float multiplication; ignores week; sensitive to float roundoff |
| `crates/gneiss-rtk/src/streaming.rs:86, 101` | `(epoch.time.tow * 1000.0).round() as u64` | RTK streaming buffer | Rollover at week boundary; assumes integer millisecond alignment |
| `crates/gneiss-rtk/src/bin/eval_odaiba_ins/main.rs:48` | `(epoch.time.tow * 1000.0).round() as u32` | Base station lookup (`base_index.get`) | Exact millisecond match; drops base if timestamp drifts by 1 ms |
| `crates/gneiss-rtk/src/bin/eval_odaiba_ins/main.rs:50, 328` | `(epoch.time.tow * 10.0).round() as u32` | GNSS RTK fix cache (`gnss_map.get`) | 100 ms (deci-second) rounding; collision if sampling rate > 10 Hz |
| `crates/gneiss-rtk/src/bin/eval_odaiba_ins/main.rs:411` | `(epoch.time.tow * 1_000_000.0).round() as u64` | IMU sample advancement | Microsecond rounding; assumes IMU TOW in microseconds within week |
| `crates/gneiss-rtk/src/bin/eval_odaiba_ins/main.rs:471` | `((time.tow * 10.0).round() as u32)` | Smoothed trajectory output | Decisecond rounding; ignores week |
| `crates/gneiss-rtk/src/bin/eval_odaiba_ins/odaiba_helpers.rs:29` | `(tow * 1_000_000.0).round() as u64` | IMU CSV parsing (`time_us`) | Quantization of float TOW from CSV |
| `crates/gneiss-rtk/src/bin/eval_odaiba_ins/odaiba_helpers.rs:86` | `((tow * 10.0).round() as u32)` | Truth reference parsing | Discretizes truth to 10 Hz |
| `crates/gneiss-rtk/src/bin/eval_odaiba_ins/odaiba_helpers.rs:118` | `((e.time.tow * 1000.0).round() as u32)` | Base RINEX indexing | Integer millisecond map key |
| `crates/gneiss-rtk/src/bin/eval_odaiba.rs:72-73` | `tow_sec = epoch.time.tow.floor() as u32`<br>`exact_ms = (epoch.time.tow * 1000.0).round() as u32` | Base vs truth matching | Inconsistent: floor() for truth, round() for base |
| `crates/gneiss-rtk/src/bin/eval_f9p_rover.rs:219, 293` | `(tow_s * 1000.0).round() as i64` | Truth window matching (`±200ms`) | Ignores week |
| `crates/gneiss-rtk/src/bin/eval_f9p_rover.rs:398, 402` | `(b.time.tow * 10.0).round() as i64` | Rover/base synchronous filtering | Decisecond rounding set intersection; ignores week |
| `crates/gneiss-rtk/src/bin/eval_qinertia_ppk.rs:64, 87, 178` | `tow.round() as u32` | Truth matching | 1 Hz integer second discretization |
| `crates/gneiss-rtk/src/bin/eval_network_ppk.rs:110, 184` | `gps_time.tow.round() as u32` | Truth matching | 1 Hz integer second discretization |
| `crates/gneiss-rtk/src/bin/eval_ppp.rs:59, 77, 382` | `ep.time.tow.round() as u32` | POS / CSRS truth matching | 1 Hz integer second discretization |

### 3.2 Direct `tow - tow` Subtractions (Week-Rollover Vulnerability)
In addition to rounding, 16 critical sites subtract raw `tow` values directly rather than subtracting full `GpsTime` structs:
- `crates/gneiss-rtk/src/estimators/rtk_iekf/predict.rs:28`: `let dt = target_time.tow - state.time.tow;`
- `crates/gneiss-rtk/src/swfg/engine/epoch.rs:21`: `(a.toe().tow - time.tow).abs()` in `select_best_ephemeris`
- `crates/gneiss-rtk/src/composite/tc_rtk.rs:165`: `let dt = self.last_time.map_or(0.1, |t| (epoch_time.tow - t.tow).clamp(0.0, 10.0));`
- `crates/gneiss-rtk/src/composite/tc_ppp.rs:186`: `let dt = self.last_time.map_or(0.1, |t| (epoch_time.tow - t.tow).clamp(0.0, 10.0));`
- `crates/gneiss-rtk/src/swfg/kalman_smoother.rs:124, 221`: `let dt = self.history[k + 1].time.tow - self.history[k].time.tow;`
- `crates/gneiss-rtk/src/swfg/engine/mod.rs:162, 319`: `let dt_sec = self.prev_time.map_or(1.0, |t| (rover.time.tow - t.tow).abs());`
- `crates/gneiss-rtk/src/post_process/screening.rs:85, 230`: `let dt = (epoch.time.tow - prev_tow).abs();`

**Consequence**: When an epoch crosses the week boundary (Sunday 00:00:00 GPST, TOW = 0.0), `dt = 0.0 - 604799.0 = -604799.0` seconds. In the Kalman filter, multiplying velocity by `-604799.0` causes the state covariance and position to instantly diverge. In `select_best_ephemeris`, a 7-day-old ephemeris is selected over the freshly broadcast ephemeris.

---

## 4. Leap-Second Handling & Cross-System Offsets

### 4.1 GPS vs. UTC 18-Second Offset
- **Standard**: GPS Time was synchronous with UTC at epoch 1980-01-06 00:00:00 UTC. Since 2017-01-01, 18 leap seconds have been inserted into UTC by the IERS ($\text{TAI} - \text{UTC} = 37\text{ s}$, $\text{TAI} - \text{GPST} = 19\text{ s} \implies \text{GPST} - \text{UTC} = 18\text{ s}$).
- **Vulnerabilities in Gneiss**:
  1. `crates/gneiss-parsers/src/antex.rs:194`:
     ```rust
     let unix_s = 315_964_800 + (time.week as i64) * 604_800 + (time.tow as i64);
     let dt = DateTime::from_timestamp(unix_s, 0)?;
     self.find_satellite(prn, dt)
     ```
     `315_964_800` is the Unix timestamp of 1980-01-06 00:00:00 UTC. But Unix time does not count leap seconds; it tracks UTC. Adding elapsed GPS seconds without subtracting 18 leap seconds produces a `dt` that is **18 seconds in the future**.
  2. `crates/gneiss-fetch/src/sources/{bkg, noaa, cddis}.rs`:
     ```rust
     let gps_epoch = chrono::NaiveDate::from_ymd_opt(1980, 1, 6)...;
     let seconds = (time.week as i64 * 604800) + time.tow as i64;
     let utc_time = gps_epoch + chrono::Duration::seconds(seconds);
     let year = utc_time.format("%Y").to_string();
     let doy = utc_time.format("%j").to_string();
     ```
     Called `utc_time`, but it is actually GPS calendar time. For observations within the first 18 seconds of a UTC day or year, this queries the wrong Day-of-Year (DOY) or Year from the IGS / NOAA / CDDIS servers.
  3. `crates/gneiss-parsers/src/rinex/obs/header.rs`:
     RINEX headers often contain a `LEAP SECONDS` record (e.g. `18    LEAP SECONDS`). The parser ignores this line entirely.

### 4.2 GLONASS Time & Leap Seconds
- **Standard**: GLONASST is referenced to Moscow Time: $\text{GLONASST} = \text{UTC(SU)} + 3\text{ hours}$. Because UTC(SU) contains leap seconds, GLONASST periodically steps by leap seconds.
  $$\text{GPST} = \text{GLONASST} - 10800\text{ s} + \Delta t_{LS}$$
  Under $\Delta t_{LS} = 18\text{ s}$, $\text{GPST} = \text{GLONASST} - 10782\text{ s}$.
- **Vulnerabilities in Gneiss**:
  1. `crates/gneiss-parsers/src/rtcm3/msm/decoder.rs:12`:
     ```rust
     let time = GpsTime::new(0, self.header.epoch_time as f64 / 1000.0);
     ```
     For GLONASS MSM (messages 1081–1087), `epoch_time` is GLONASS Time of Day (or Day of Week + TOD in UTC(SU)+3h). Storing this in `GpsTime` without applying $-10782\text{ s}$ causes a 3-hour error in GLONASS pseudo-range alignment.
  2. `crates/gneiss-parsers/src/rinex/nav/mod.rs:95`:
     ```rust
     Constellation::Glonass => toc_gpst = toc_gpst + TimeSystem::Glonass.gpst_offset(),
     ```
     Hardcodes `GPS_LEAP_SECONDS = 18`. If past archival data (e.g. 2012, when $\Delta t_{LS} = 15$ or $16$) or future data is processed, the offset is incorrect.

### 4.3 BeiDou 14-Second Offset & Critical Ephemeris Bug
- **Standard**: BDT commenced at 2006-01-01 00:00:00 UTC. At that epoch, GPS had accumulated 14 leap seconds relative to UTC ($\text{GPST} - \text{UTC} = 14\text{ s}$). Since neither GPST nor BDT introduces leap seconds:
  $$\text{GPST} - \text{BDT} = 14.000\text{ s (exact and permanent)}$$
- **The Ephemeris Bug Discovered in Gneiss**:
  In `crates/gneiss-parsers/src/rinex/nav/mod.rs:96`:
  ```rust
  Constellation::Beidou => toc_gpst = toc_gpst + TimeSystem::Bdt.gpst_offset(), // +14.0s
  ```
  `toc` is converted to GPST.
  However, in `crates/gneiss-parsers/src/rinex/nav/builder.rs:144`:
  ```rust
  fn build_beidou_ephemeris(..., toc: GpsTime, vals: &[f64; 32]) -> Option<Ephemeris> {
      ...
      toe: GpsTime::new(toc.week, vals[8]), // vals[8] is BDT seconds of week!
  ```
  `vals[8]` in RINEX 3 NAV is BDT TOE (0 to 604800). It is assigned directly to `toe` without adding 14 seconds.
  Then in `crates/gneiss-core/src/ephemeris/keplerian.rs:178-186`:
  ```rust
  pub fn position(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
      let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
      calc_keplerian(t_bdt, self.toe, self.toc, ...)
  }
  ```
  And inside `calc_keplerian` (`crates/gneiss-core/src/keplerian.rs:143, 179`):
  ```rust
  let tk = t - toe; // t_bdt - self.toe -> (BDT - BDT) = correct
  let tc = t - toc; // t_bdt - self.toc -> (BDT - GPST) = -14.0s ERROR!
  let clk_err = af0 + af1 * tc + af2 * tc * tc ...;
  ```
  Because `toc` was shifted to GPST (+14s) while `t_bdt` was shifted to BDT (-14s), `tc` contains a **28-second (or 14-second) timing error**, directly corrupting satellite clock bias calculation for all BeiDou satellites!

---

## 5. IMU and GNSS Observable Synchronization

In `crates/gneiss-rtk/src/bin/eval_odaiba_ins/`:
- **Sensors**:
  - GNSS Rover: 10 Hz (`EpochObs` with `GpsTime`)
  - GNSS Base: 10 Hz (`EpochObs` with `GpsTime`)
  - IMU: 50 Hz / 100 Hz (`ImuRecord` with `time_us: u64`)
- **Pipeline Synchronization Mechanism**:
  1. `imu.csv` is parsed into `ImuRecord` where `time_us = (tow * 1_000_000.0).round() as u64`.
  2. Outer loop steps through `rover_epochs` at 10 Hz.
  3. `cur_us = (epoch.time.tow * 1_000_000.0).round() as u64;`
  4. `advance_imu_records` pushes all IMU samples whose `sample.time_us <= cur_us` into `acc_imu`.
  5. The last IMU sample is preserved across epochs (`acc_imu.push(last_sample)`) to ensure continuous integration boundaries.
  6. `preint.integrate` uses `time_diff_us(prev, curr)` with a heuristic wrap check at `604_800_000_000` us (`WEEK_US`).
- **Critical Fragilities Identified**:
  1. **Brittle Base Station Key Lookup**:
     ```rust
     let exact_ms = (epoch.time.tow * 1000.0).round() as u32;
     let base = base_index.get(&exact_ms).copied()?;
     ```
     This requires exact integer millisecond match. If rover and base receivers have an unsynchronized initial offset (e.g. rover at `.002s`, base at `.000s`) or clock drift jitter, `base_index.get(&exact_ms)` returns `None`.
  2. **Decisecond Quantization for Truth & Fix Cache**:
     `epoch_key = (epoch.time.tow * 10.0).round() as u32;`
     Used as key in `GnssFixMap` and truth reference map. This hardcodes the pipeline to $\le 10\text{ Hz}$ GNSS data; any 20 Hz or 50 Hz GNSS receiver causes collisions.
  3. **No Week Tracking in IMU**:
     `ImuSample.time_us` stores microseconds within the week. If an IMU data collection spans midnight Saturday / Sunday, `time_diff_us` executes its rollover branch, but if timestamps are outside the arbitrary `30_000_000` us window, it wraps mod `u64::MAX`.

---

## 6. Proposed Strictly Typed Temporal Architecture

To eliminate all ad-hoc rounding, truncation, and cross-system timing errors, we propose a zero-cost typestate temporal architecture in `crates/gneiss-core`.

### 6.1 Type Architecture

```rust
pub trait TimeScale: Copy + Clone + PartialEq + Eq + core::fmt::Debug + 'static {
    const NAME: &'static str;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GpsScale;
impl TimeScale for GpsScale { const NAME: &'static str = "GPST"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BdtScale;
impl TimeScale for BdtScale { const NAME: &'static str = "BDT"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GstScale;
impl TimeScale for GstScale { const NAME: &'static str = "GST"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GlonassScale;
impl TimeScale for GlonassScale { const NAME: &'static str = "GLONASST"; }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UtcScale;
impl TimeScale for UtcScale { const NAME: &'static str = "UTC"; }

/// Strictly typed, normalized epoch in a specific GNSS or civil time scale.
#[derive(Debug, Clone, Copy)]
pub struct Epoch<S: TimeScale> {
    week: u32,
    tow_nanos: u64, // [0, 604_800_000_000_000)
    _scale: core::marker::PhantomData<S>,
}

pub type GpsTime = Epoch<GpsScale>;
pub type BdtTime = Epoch<BdtScale>;
pub type GstTime = Epoch<GstScale>;
pub type GlonassTime = Epoch<GlonassScale>;
pub type UtcTime = Epoch<UtcScale>;
```

### 6.2 Nanosecond-Precision Integer Representation
By storing `week: u32` and `tow_nanos: u64` (integer nanoseconds in week, $0 \le \text{tow\_nanos} < 604{,}800 \times 10^9$):
- **Exact millisecond, microsecond, and nanosecond comparisons**: 0 float rounding errors.
- **Full range**: 1 picosecond or nanosecond resolution without f64 precision degradation.
- **Normalization Invariant**: Enforced by private fields and validated constructors (`Epoch::from_week_nanos(week, nanos)`).
- **Exact Ord and Eq**: `Ord` compares `(week, tow_nanos)` lexicographically, which is provably monotonic and consistent.

### 6.3 Strongly Typed Duration (`TimeDelta`)
```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct TimeDelta {
    nanos: i64,
}

impl TimeDelta {
    pub const fn from_nanos(nanos: i64) -> Self { Self { nanos } }
    pub const fn from_micros(micros: i64) -> Self { Self { nanos: micros * 1_000 } }
    pub const fn from_millis(millis: i64) -> Self { Self { nanos: millis * 1_000_000 } }
    pub fn from_seconds(secs: f64) -> Self { Self { nanos: (secs * 1e9).round() as i64 } }
    
    pub fn as_seconds(&self) -> f64 { self.nanos as f64 * 1e-9 }
    pub fn as_millis(&self) -> i64 { self.nanos / 1_000_000 }
    pub fn as_micros(&self) -> i64 { self.nanos / 1_000 }
    pub const fn as_nanos(&self) -> i64 { self.nanos }
}
```

### 6.4 Compile-Time Arithmetic Safety
- `Epoch<S> - Epoch<S> -> TimeDelta` (Same scale: allowed and accounts for week differences).
- `Epoch<GpsScale> - Epoch<BdtScale>` -> **COMPILE-TIME ERROR**.
- `Epoch<GpsScale> - Epoch<UtcScale>` -> **COMPILE-TIME ERROR**.
- `Epoch<S> + TimeDelta -> Epoch<S>`.
- `Epoch<S> - TimeDelta -> Epoch<S>`.

### 6.5 Conversion Semantics
1. **Infallible Scale Conversions**:
   - `BdtTime::to_gpst(&self) -> GpsTime`: adds exact 14.0 seconds ($14 \times 10^9\text{ ns}$).
   - `GpsTime::to_bdt(&self) -> BdtTime`: subtracts exact 14.0 seconds.
   - `GstTime::to_gpst(&self) -> GpsTime`: nominal 0.0 second offset.
2. **Leap-Second-Explicit Conversions**:
   - `GlonassTime::to_gpst(&self, leap_seconds: i32) -> GpsTime`:
     $$\text{GPST} = \text{GLONASST} - 10800\text{ s} + \text{leap\_seconds}$$
   - `UtcTime::to_gpst(&self, leap_seconds: i32) -> GpsTime`:
     $$\text{GPST} = \text{UTC} + \text{leap\_seconds}$$
   - `GpsTime::to_utc(&self, leap_seconds: i32) -> UtcTime`:
     $$\text{UTC} = \text{GPST} - \text{leap\_seconds}$$

### 6.6 Collision-Free Alignment Keys (`EpochKey`)
Instead of ad-hoc float rounding:
```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EpochKey<S: TimeScale> {
    continuous_ms: u64,
    _scale: core::marker::PhantomData<S>,
}

impl<S: TimeScale> Epoch<S> {
    /// Continuous millisecond timestamp since scale origin, immune to week rollover.
    pub fn continuous_ms(&self) -> u64 {
        self.week as u64 * 604_800_000 + (self.tow_nanos / 1_000_000)
    }

    pub fn to_key(&self) -> EpochKey<S> {
        EpochKey { continuous_ms: self.continuous_ms(), _scale: core::marker::PhantomData }
    }

    /// Check if another epoch is within a given tolerance window.
    pub fn is_within(&self, other: Self, tolerance: TimeDelta) -> bool {
        let diff = if *self >= other { *self - other } else { other - *self };
        diff <= tolerance
    }
}
```

---

## 7. Comprehensive Refactoring Roadmap & Module Impact

### 7.1 Detailed Module Impact Analysis

| Crate | Module / File | Current Issue | Required Refactoring |
|---|---|---|---|
| `gneiss-core` | `src/time.rs` | `GpsTime` bare floats, unnormalized ordering, incorrect calendar math | Implement `Epoch<Scale>`, `TimeScale` traits, `TimeDelta`, and integer nanosecond representation |
| `gneiss-core` | `src/gnss_time.rs` | Runtime enum `GnssTime`, public fields | Migrate to typestate `Epoch<Scale>`; preserve `TimeSystem` for dynamic serialization |
| `gneiss-core` | `src/ephemeris/keplerian.rs` | BeiDou `t_bdt` using `GpsTime::new(t.week, t.tow - 14.0)` | Replace with typed `BdtTime`; ensure `toe` and `toc` share consistent scale in `calc_keplerian` |
| `gneiss-core` | `src/keplerian.rs` | `calc_keplerian` takes raw `GpsTime` for all constellations | Enforce scale consistency; update Keplerian propagation signatures |
| `gneiss-core` | `src/ephemeris/glonass.rs` | `toe: GpsTime` | Type as `GlonassTime` or explicitly documented GPST |
| `gneiss-core` | `src/obs.rs` | `EpochObs.time: GpsTime` | Standardize on normalized `GpsTime` (`Epoch<GpsScale>`) |
| `gneiss-parsers` | `src/rinex/nav/mod.rs` | Ad-hoc `+ gpst_offset()` on parsed `GpsTime` | Parse directly into native `GlonassTime` / `BdtTime`, then convert explicitly via `.to_gpst()` |
| `gneiss-parsers` | `src/rinex/nav/builder.rs` | `toe` in BDT while `toc` in GPST in `build_beidou_ephemeris` | Fix BeiDou ephemeris: apply +14s to `toe` when converting to GPST, or store both in native BDT |
| `gneiss-parsers` | `src/rinex/obs/header.rs` | `TIME SYSTEM` and `LEAP SECONDS` ignored | Parse and store header time system and leap seconds |
| `gneiss-parsers` | `src/rtcm3/msm/decoder.rs` | Hardcoded `week: 0`, GLONASS TOD treated as GPST | Apply constellation-specific time decoding (GLONASS TOD vs BDT TOW vs GPST TOW) |
| `gneiss-parsers` | `src/antex.rs` | Missing 18 leap seconds in Unix timestamp | Subtract 18 leap seconds in `find_satellite_gps` |
| `gneiss-parsers` | `src/sp3.rs` | Header time system ignored | Parse SP3 time scale and convert to GPST |
| `gneiss-parsers` | `src/precise_orbit.rs` | `s_time` continuous float conversion | Replace with `TimeDelta` / `continuous_seconds()` |
| `gneiss-fetch` | `src/sources/{bkg, noaa, cddis}.rs` | GPS calendar dates used for UTC DOY / Year queries | Subtract 18 leap seconds before computing UTC DOY and Year |
| `gneiss-rtk` | `src/swfg/engine/epoch.rs` | `(a.toe().tow - time.tow).abs()` in `select_best_ephemeris` | Use `(a.toe() - time).abs()` to prevent week-rollover misselection |
| `gneiss-rtk` | `src/post_process/network.rs` | `tow_key` 1 Hz rounding with `rem_euclid(604800)` | Replace with continuous `EpochKey` |
| `gneiss-rtk` | `src/post_process/iekf_pass.rs` | `find_matched_base` ignores week | Use `(b.time - target).abs() < TimeDelta::from_millis(100)` |
| `gneiss-rtk` | `src/post_process/combiner.rs` | `(fwd.time.tow * 1000.0).round() as u64` | Replace with `EpochKey` |
| `gneiss-rtk` | `src/post_process/backward.rs` | Float rounding in `group_imu_by_epoch` | Replace with continuous `EpochKey` / `TimeDelta` |
| `gneiss-rtk` | `src/streaming.rs` | Float rounding in buffer keys | Replace with `EpochKey` |
| `gneiss-rtk` | `src/estimators/rtk_iekf/predict.rs` | `target_time.tow - state.time.tow` | Replace with `target_time - state.time` (`TimeDelta`) |
| `gneiss-rtk` | `src/swfg/imu_preintegration/mod.rs` | `time_diff_us` with `WEEK_US` heuristic | Use typed `TimeDelta` from continuous timestamps |
| `gneiss-rtk` | `src/bin/eval_odaiba_ins/main.rs` | `exact_ms`, `epoch_key`, `cur_us` float roundings | Replace with `EpochKey` and `tolerance` matching |
| `gneiss-rtk` | `src/bin/eval_f9p_rover.rs` | `* 10.0` decisecond rounding | Replace with `EpochKey` |
| `gneiss-rtk` | `src/bin/eval_network_ppk.rs` | `tow.round() as u32` truth matching | Replace with `EpochKey` |
| `gneiss-rtk` | `src/bin/eval_ppp.rs` | `tow.round() as u32` truth matching | Replace with `EpochKey` |

### 7.2 Implementation Phases (TDD Approach)

```
Phase 1: Core Primitives (gneiss-core)
├── Implement TimeScale traits (GpsScale, BdtScale, GstScale, GlonassScale, UtcScale)
├── Implement Epoch<S> with (week: u32, tow_nanos: u64)
├── Implement TimeDelta (nanoseconds i64)
├── Write unit tests pinning:
│   ├── ICD offsets (BDT +14s, GLONASS -10782s)
│   ├── Week carries at 0.0 and 604800.0
│   ├── Compile-time failure on cross-system subtraction
│   └── Nanosecond-exact arithmetic round trips
└── Provide backward-compatible type alias pub type GpsTime = Epoch<GpsScale>

Phase 2: Parsers & External Interfaces (gneiss-parsers & gneiss-fetch)
├── Fix BeiDou ephemeris builder: align toe and toc scales
├── Fix antex.rs PCV lookup: apply 18 leap seconds
├── Fix gneiss-fetch DOY / Year queries: apply 18 leap seconds
├── Update rinex/nav to construct native GlonassTime and BdtTime before explicit conversion
└── Update SP3 and RTCM decoders to respect constellation time scales

Phase 3: Estimators & Engine (gneiss-rtk)
├── Fix select_best_ephemeris to use (toe - time).abs() across week boundaries
├── Replace raw tow - tow subtractions in rtk_iekf, kalman_smoother, tc_rtk, tc_ppp
├── Replace tow_key in network.rs with EpochKey
├── Update find_matched_base in iekf_pass.rs to use TimeDelta tolerance
└── Update imu_preintegration to consume TimeDelta

Phase 4: Benchmarks & Pipelines
├── Refactor eval_odaiba_ins: replace ad-hoc float rounding with EpochKey
├── Refactor eval_f9p_rover, eval_network_ppk, eval_ppp truth matching
├── Verify all 789+ workspace tests pass (cargo test --workspace)
├── Run python3 scripts/check_network_benchmark.py --smoke
├── Run python3 scripts/check_multignss_benchmark.py --smoke
└── Verify eval_odaiba_ins reproduces benchmark metrics (p50 <= 1.80m, RMS <= 3.50m)
```

---

## 8. Conclusion

The Gneiss codebase currently relies on dangerous ad-hoc time handling that mixes different constellation time scales inside `GpsTime`, bypasses week arithmetic via direct `tow - tow` subtractions, drops leap seconds when interacting with civil time and PCV files, and relies on brittle float rounding for epoch matching.

By introducing the proposed `Epoch<Scale>` zero-cost typestate and `TimeDelta` architecture, the compiler will structurally reject cross-time-system arithmetic, week-rollover failures will be eradicated, and epoch matching across GNSS and IMU streams will be exact, robust, and mathematically sound.
