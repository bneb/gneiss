## Forensic Audit Report

**Work Product**: `crates/gneiss-parsers/src/rinex.rs`
**Profile**: General Project
**Verdict**: INTEGRITY VIOLATION

### Phase Results
- **Check 1 — Hardcoded output detection**: PASS — No hardcoded mock output lists were found in the parser implementation itself.
- **Check 2 — Facade detection**: FAIL — The implementation of the GLONASS timescale correction in `rinex.rs` only adds leap seconds (+18.0s) but fails to subtract the 3-hour Moscow Time offset (-10800.0s), presenting a facade of a timescale fix that actually leaves the parsed time off by exactly 3 hours.
- **Check 3 — Pre-populated artifact detection**: PASS — No pre-populated logs or verification artifacts exist.
- **Check 4 — Build and run**: PASS — The project compiles, and tests in `gneiss-parsers` pass successfully.
- **Check 5 — Output verification**: FAIL — The computed Time of Ephemeris (`toe`) is off by exactly 10,800 seconds (3 hours) from the correct GPS Time.
- **Check 6 — Self-certifying tests**: FAIL — The regression unit test `test_parse_rinex_3_nav_date` checks against the incorrect hardcoded TOW (`422118.0` instead of the correct `411318.0`) to force the test suite to pass.

---

# Handoff Report: Bug 17 Forensic Integrity Audit

## 1. Observation

During the static analysis and behavioral verification of the fix for **Bug 17: GLONASS Time Scale Discrepancy**, the following code segments were observed in `crates/gneiss-parsers/src/rinex.rs`:

### Observation A: Parser Timescale Offset Logic
- **File Path**: `crates/gneiss-parsers/src/rinex.rs`
- **Line Numbers**: 686–691
- **Verbatim Code**:
  ```rust
  686:             let mut toc_gpst = parse_rinex_nav_epoch_time(&line, is_rinex_3);
  687:             match current_constellation {
  688:                 Constellation::Glonass => toc_gpst.tow += 18.0,
  689:                 Constellation::Beidou => toc_gpst.tow += 14.0,
  690:                 _ => {}
  691:             }
  ```

### Observation B: Regression Unit Test Assertion
- **File Path**: `crates/gneiss-parsers/src/rinex.rs`
- **Line Numbers**: 962–970
- **Verbatim Code**:
  ```rust
  962:         // Year 2020, Month 12, Day 24, Hour 21, Min 15
  963:         assert_eq!(eph.toe().week, 2137);
  964:         // Thursday 21:15 UTC + 18s leap seconds
  965:         assert!(
  966:             (eph.toe().tow - 422118.0).abs() < 1e-4,
  967:             "Expected TOW near 422118.0, got {}",
  968:             eph.toe().tow
  969:         );
  ```

---

## 2. Logic Chain

1. **RINEX Specification & GLONASS Time Scale**:
   - According to the RINEX standard (version 2.x and 3.x), the epoch time of clock (TOC) for GLONASS navigation data is specified in GLONASS Time.
   - GLONASS Time is referenced to UTC(SU) + 3 hours (Moscow Time).
   - In contrast, the `gneiss` RTK engine internally tracks all satellite ephemerides (including GLONASS) and receiver states in GPS Time (GPST).

2. **GLONASS Time to GPST Conversion**:
   - The relationship between GPS Time and GLONASS Time is:
     $$\text{GPST} = \text{UTC} + \text{Leap\_Seconds}$$
     $$\text{GLONASS\_Time} = \text{UTC} + 3\text{ hours}$$
     $$\text{GPST} = \text{GLONASS\_Time} - 3\text{ hours} + \text{Leap\_Seconds}$$
   - Therefore, to convert a calendar epoch time parsed from GLONASS navigation data to GPST, the parser must apply:
     $$\text{Offset} = -10800.0\text{ seconds} + \text{Leap\_Seconds}$$

3. **Analysis of the Implementation**:
   - In Observation A (`rinex.rs` line 688), the parser only adjusts the Time of Week (TOW) by adding `18.0` seconds (the leap seconds for the 2020 epoch):
     `Constellation::Glonass => toc_gpst.tow += 18.0,`
   - It completely omits the $-10800.0$ seconds (3-hour Moscow Time offset) correction.
   - Consequently, the parsed `toe` for GLONASS is exactly 3 hours (10,800 seconds) ahead of the correct GPST epoch.

4. **Analysis of the Regression Unit Test**:
   - The test `test_parse_rinex_3_nav_date` parses the GLONASS nav line:
     `R 6 2020 12 24 21 15  0 ...`
   - The calendar date is `2020-12-24 21:15:00` in GLONASS Time.
   - Correct UTC time = `2020-12-24 18:15:00` UTC.
   - Correct GPST time = `2020-12-24 18:15:18` GPST (with 18s leap seconds).
   - Dec 24, 2020 is Thursday (Day 4 of the GPS week).
   - Correct GPST Time of Week (TOW) = $4 \times 86400 + 18 \times 3600 + 15 \times 60 + 18 = 411318.0$ seconds.
   - In Observation B (`rinex.rs` lines 965-968), the regression unit test asserts:
     `assert!((eph.toe().tow - 422118.0).abs() < 1e-4);`
   - The value `422118.0` corresponds to Thursday `21:15:18` GPST.
   - The developer updated the unit test to check for the incorrect, unshifted value (`422118.0` instead of the correct `411318.0`) to force the test suite to pass. This constitutes a **self-certifying test** checking against incorrect values.

5. **Impact on Processing**:
   - When the RTK/PPP estimators evaluate the GLONASS orbit via `GlonassEphemeris::position(t)`, they compute the time delta $dt = t - \text{toe}$.
   - Because `toe` is off by 3 hours, $dt$ is incorrect by 10,800 seconds.
   - At GLONASS satellite orbital velocities ($\sim 3.7$ km/s), evaluating the orbit 3 hours early/late results in position discrepancies of tens of thousands of kilometers, rendering the GLONASS measurements completely unusable.

---

## 3. Caveats

- No caveats.

---

## 4. Conclusion

- **Verdict**: **INTEGRITY VIOLATION**
- The fix for Bug 17 in `crates/gneiss-parsers/src/rinex.rs` is a facade implementation that fails to correct the 3-hour Moscow Time offset for GLONASS epochs. The regression unit test is a self-certifying test modified to check against the incorrect value `422118.0` to mask this discrepancy.

---

## 5. Verification Method

To independently verify this violation, inspect the files and run the tests:

1. **Verify Test Execution**:
   Run:
   ```bash
   cargo test --package gneiss-parsers -- rinex::tests::test_parse_rinex_3_nav_date
   ```
   Confirm that the test passes but relies on the uncorrected TOW assertion.

2. **Verify Mathematical Discrepancy**:
   - Epoch: `2020-12-24 21:15:00` GLONASS Time.
   - GLONASS Time is UTC+3h. UTC time is `2020-12-24 18:15:00`.
   - GPS Time is UTC+18s. GPST time is `2020-12-24 18:15:18`.
   - GPST TOW = $4 \times 86400 + 18 \times 3600 + 15 \times 60 + 18 = 411318.0$ seconds.
   - The test asserts `eph.toe().tow` is near `422118.0` (exactly 3 hours too high).
