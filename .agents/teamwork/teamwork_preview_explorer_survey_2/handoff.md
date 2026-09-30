# Handoff Report: Survey 2 — Temporal Frame & Epoch Alignment Safety (R2)

**Author**: Survey Explorer 2  
**Date**: 2026-09-25  
**Working Directory**: `/Users/kevin/projects/gneiss/.agents/teamwork/teamwork_preview_explorer_survey_2`  
**Recipient**: Orchestrator (`db66ae0c-b21b-4e14-ac97-93509c51c4b0`)  
**Type**: Hard Handoff (Investigation Complete)  

---

## 1. Observation

Direct code observations across `gneiss-core`, `gneiss-parsers`, `gneiss-fetch`, and `gneiss-rtk`:

1. **Bare Float Time Representation**:
   In `crates/gneiss-core/src/time.rs:7-10`:
   ```rust
   #[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
   pub struct GpsTime {
       pub week: u32,
       pub tow: f64,
   }
   ```
   `week` and `tow` are public and mutable; `PartialOrd` compares `(week, tow)` lexicographically without normalization.

2. **Critical BeiDou Broadcast Ephemeris Scale Discrepancy**:
   In `crates/gneiss-parsers/src/rinex/nav/mod.rs:96`:
   ```rust
   Constellation::Beidou => toc_gpst = toc_gpst + gneiss_core::gnss_time::TimeSystem::Bdt.gpst_offset(),
   ```
   `toc` is shifted to GPST by adding +14.0s. But in `crates/gneiss-parsers/src/rinex/nav/builder.rs:144`:
   ```rust
   toe: GpsTime::new(toc.week, vals[8]),
   ```
   `vals[8]` is raw BDT seconds of week from RINEX, leaving `toe` in BDT while `toc` is in GPST.
   In `crates/gneiss-core/src/ephemeris/keplerian.rs:178`:
   ```rust
   let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
   calc_keplerian(t_bdt, self.toe, self.toc, ...)
   ```
   In `crates/gneiss-core/src/keplerian.rs:179`:
   ```rust
   let tc = t - toc;
   ```
   Here `t` is `t_bdt` (BDT), while `toc` is GPST (BDT + 14s). Thus `tc` has a 14-second offset error in satellite clock correction.

3. **Direct `tow - tow` Subtractions Bypassing Week Arithmetic**:
   At least 16 sites bypass `GpsTime` subtraction and subtract raw floats:
   - `crates/gneiss-rtk/src/estimators/rtk_iekf/predict.rs:28`:
     `let dt = target_time.tow - state.time.tow;`
   - `crates/gneiss-rtk/src/swfg/engine/epoch.rs:21-23`:
     `(a.toe().tow - time.tow).abs().partial_cmp(&(b.toe().tow - time.tow).abs())`
   - `crates/gneiss-rtk/src/composite/tc_rtk.rs:165`:
     `let dt = self.last_time.map_or(0.1, |t| (epoch_time.tow - t.tow).clamp(0.0, 10.0));`
   - `crates/gneiss-rtk/src/swfg/kalman_smoother.rs:124, 221`:
     `let dt = self.history[k + 1].time.tow - self.history[k].time.tow;`

4. **Leap-Second Omission in External Time Conversions**:
   - `crates/gneiss-parsers/src/antex.rs:194`:
     `let unix_s = 315_964_800 + (time.week as i64) * 604_800 + (time.tow as i64);`
     Converts GPS elapsed time directly to Unix epoch seconds without subtracting 18 leap seconds.
   - `crates/gneiss-fetch/src/sources/bkg.rs:25-26`:
     `let seconds = (time.week as i64 * 604800) + time.tow as i64;`
     `let utc_time = gps_epoch + chrono::Duration::seconds(seconds);`
     Labelled `utc_time`, but it is GPS time; for the first 18 seconds of a day/year, this queries the wrong Day-of-Year or Year URL.

5. **Ad-Hoc Float Rounding in Pipeline Alignment**:
   In `crates/gneiss-rtk/src/bin/eval_odaiba_ins/main.rs`:
   - Line 48: `let exact_ms = (epoch.time.tow * 1000.0).round() as u32;`
   - Line 49: `let base = base_index.get(&exact_ms).copied()?;`
   - Line 50, 328: `let epoch_key = (epoch.time.tow * 10.0).round() as u32;`
   - Line 411: `let cur_us = (epoch.time.tow * 1_000_000.0).round() as u64;`
   - Line 471: `(((time.tow * 10.0).round() as u32), s.pos_ecef)`
   In `crates/gneiss-rtk/src/post_process/network.rs:92-94`:
   - `fn tow_key(time: &GpsTime) -> u32 { (time.tow.round() as i64).rem_euclid(604800) as u32 }`

---

## 2. Logic Chain

1. **From Observation 1**: Because `GpsTime` has public `week: u32` and `tow: f64` fields with no time scale tag, developers treat `GpsTime` as a generic tuple. There is no type-level difference between GPST, BDT, GLONASST, or UTC.
2. **From Observation 2**: Because BDT and GPST share the same struct, BeiDou broadcast ephemeris construction in `rinex/nav/builder.rs` parsed `toe` in BDT and `toc` in GPST without compiler rejection. When propagated in `calc_keplerian`, `t_bdt - toc` subtracts GPST from BDT, introducing a silent 14-second clock-error offset.
3. **From Observation 3**: Because `tow` is an exposed float, developers intuitively write `t1.tow - t2.tow` rather than `t1 - t2`. At week boundary rollover (Sunday 00:00:00 GPST), `tow` resets from 604799 to 0. The expression evaluates to `-604799.0` seconds. In Kalman prediction (`rtk_iekf/predict.rs`), velocity propagation explodes. In ephemeris selection (`swfg/engine/epoch.rs`), a 7-day-old ephemeris is selected over the active one.
4. **From Observation 4**: Because `from_calendar` claims to parse UTC calendar dates but internally computes GPS week/TOW without subtracting leap seconds, callers in `antex.rs` and `gneiss-fetch` assume `GPST == UTC` (or Unix), introducing an uncorrected 18-second shift into antenna calibrations and archive downloads.
5. **From Observation 5**: Because there is no standardized discrete epoch key or alignment tolerance helper, every binary and post-processing pass invents its own float rounding (`* 1000.0`, `* 10.0`, `round()`, `floor()`). Exact millisecond map lookups (`base_index.get(&exact_ms)`) drop base observations upon tiny sub-millisecond receiver clock drift.

---

## 3. Caveats

- **RTCM MSM Message Timing**: GLONASS MSM time (DF416) uses Day of Week and Time of Day (in Moscow time UTC+3h). A complete fix for RTCM MSM requires tracking the full GPS week from an external reference or receiver navigation state, as RTCM does not carry full GPS week numbers in MSM headers.
- **Historic Leap Second Tables**: While the current offset between GPST and UTC is fixed at 18 seconds (and $\text{GPST} - \text{BDT} = 14\text{ s}$ is permanent), processing data prior to 2017 requires a leap second table rather than a single static constant `18`.
- **Scope Limit**: As an explorer subagent in read-only mode, no production files were modified. All proposals are preserved for the implementation agents.

---

## 4. Conclusion

The temporal architecture of Gneiss requires refactoring under Requirement R2:
1. **Implement `Epoch<Scale: TimeScale>` Typestate**: Formalize `GpsScale`, `BdtScale`, `GstScale`, `GlonassScale`, `UtcScale` in `gneiss-core`. Prohibit cross-scale subtraction at compile time (`epoch_gps - epoch_bdt` fails to compile).
2. **Nanosecond Integer Representation**: Represent epochs internally as `(week: u32, tow_nanos: u64)` and duration as `TimeDelta { nanos: i64 }`. This provides exact integer precision, monotonic ordering, and eliminates float roundoff jitter.
3. **Explicit Leap Second Conversions**: Make leap seconds explicit in all UTC and GLONASS conversions (`to_gpst_with_leap_seconds`). Fix `antex.rs` and `gneiss-fetch` to apply the 18-second leap offset.
4. **Fix BeiDou Ephemeris Scale Alignment**: Align `toe` and `toc` in `BeidouEphemeris` to eliminate the 14-second clock calculation error.
5. **Standardized `EpochKey` & Tolerance Matching**: Replace all 16 ad-hoc float rounding sites with continuous `EpochKey` and tolerance-based matching (`is_within(tolerance)`), making base matching and IMU synchronization robust against week rollovers and sub-millisecond drift.

---

## 5. Verification Method

To verify these findings and validate the refactoring when implemented:

1. **Verify Existing Tests**:
   ```bash
   cargo test --workspace
   ```
   (Must pass with 0 failures and 0 warnings).
2. **Verify Regression Guards**:
   ```bash
   cargo build --release --bin eval_network_ppk
   python3 scripts/check_network_benchmark.py --smoke
   python3 scripts/check_multignss_benchmark.py --smoke
   ```
3. **Verify `eval_odaiba_ins` Metric Compliance**:
   ```bash
   cargo run --release --bin eval_odaiba_ins
   ```
   Assert $p_{50} \le 1.80\text{ m}$, $\text{RMS} \le 3.50\text{ m}$, 0 false fixes.
4. **New Temporal Compile-Fail & Unit Invariant Tests**:
   - Write unit tests in `crates/gneiss-core/src/time.rs` asserting:
     - Cross-system subtraction (`Epoch<GpsScale> - Epoch<BdtScale>`) fails to compile (`compile_fail` test).
     - Subtraction across week boundaries yields exact elapsed seconds:
       `Epoch::new(2201, 10.0) - Epoch::new(2200, 604790.0) == TimeDelta::from_seconds(20.0)`.
     - `BeidouEphemeris` satellite position and clock bias are identical whether evaluated in GPST or BDT.
     - `EpochKey` remains monotonic across GPS week boundaries.
