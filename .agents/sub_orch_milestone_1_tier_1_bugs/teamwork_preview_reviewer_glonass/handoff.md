# Handoff Report: Bug 17 GLONASS Time Scale Discrepancy Verification

This report documents the verification and review of the changes implemented for **Bug 17: GLONASS Time Scale Discrepancy**.

---

## 1. Observation

### Code Changes
Only one file, `crates/gneiss-parsers/src/rinex.rs`, has been modified. Below is the exact `git diff` output:

```diff
diff --git a/crates/gneiss-parsers/src/rinex.rs b/crates/gneiss-parsers/src/rinex.rs
index 475cc97..bedb024 100644
--- a/crates/gneiss-parsers/src/rinex.rs
+++ b/crates/gneiss-parsers/src/rinex.rs
@@ -685,8 +685,8 @@ pub fn parse_rinex_nav<R: BufRead>(
 
             let mut toc_gpst = parse_rinex_nav_epoch_time(&line, is_rinex_3);
             match current_constellation {
-                Constellation::Glonass => toc_gpst.tow += 18.0,
-                Constellation::Beidou => toc_gpst.tow += 14.0,
+                Constellation::Glonass => toc_gpst = toc_gpst + (18.0 - 10800.0),
+                Constellation::Beidou => toc_gpst = toc_gpst + 14.0,
                 _ => {}
             }
             current_toc = toc_gpst;
@@ -961,10 +961,10 @@ R 6 2020 12 24 21 15  0  .189751386642E-03  .000000000000E+00  .422910000000E+06
 
         // Year 2020, Month 12, Day 24, Hour 21, Min 15
         assert_eq!(eph.toe().week, 2137);
-        // Thursday 21:15 UTC + 18s leap seconds
+        // Thursday 21:15 UTC + 18s leap seconds - 3h Moscow offset
         assert!(
-            (eph.toe().tow - 422118.0).abs() < 1e-4,
-            "Expected TOW near 422118.0, got {}",
+            (eph.toe().tow - 411318.0).abs() < 1e-4,
+            "Expected TOW near 411318.0, got {}",
             eph.toe().tow
         );
     }
```

### Verification Commands & Results

1. **Parser-specific tests**:
   Command: `cargo test -p gneiss-parsers`
   Result: **PASSED** (13 tests passed, including `rinex::tests::test_parse_rinex_3_nav_date`).
   
2. **Workspace tests**:
   Command: `cargo test --workspace`
   Result: **PASSED** (All 256 tests in the workspace passed with 0 failures).

3. **Code Formatting check**:
   Command: `cargo fmt --package gneiss-parsers -- --check`
   Result: **PASSED** (No formatting issues in `gneiss-parsers`).

4. **Layout Compliance check**:
   The `.agents/` folder contains only markdown agent configuration/metadata files (`BRIEFING.md`, `progress.md`, `handoff.md`, `ORIGINAL_REQUEST.md`, `SCOPE.md`). No source code, tests, or raw data files exist there.

---

## 2. Logic Chain

1. **GLONASS Moscow Time Offset**: GLONASS Time (GLONASST) is linked to Moscow Time, which is UTC+3 (10,800 seconds). GPST is UTC + leap_seconds (which was 18.0 seconds for the epoch).
2. **Mathematical Correctness**: Converting GLONASST epoch calendar date to GPST requires subtracting the 3-hour Moscow Time offset and adding the leap second offset:
   $$\text{GPST} = \text{GLONASST} - 10800.0\,\text{s} + 18.0\,\text{s} = \text{GLONASST} + (18.0 - 10800.0)\,\text{s}$$
3. **Normalization Safety**: Directly updating the `tow` field (e.g. `toc_gpst.tow += 18.0 - 10800.0`) would leave the `GpsTime` struct in an unnormalized state with a negative `tow`. Using the `Add<f64>` trait implementation (i.e. `toc_gpst = toc_gpst + (18.0 - 10800.0)`) constructs a new `GpsTime` and calls `normalize()`, which correctly adjusts both `tow` (adding $604,800.0$ seconds if negative) and decrementing the continuous GPS `week` accordingly.
4. **Test Alignment**: The unit test `test_parse_rinex_3_nav_date` parses the calendar date `2020 12 24 21 15 0`. Thursday 21:15:00 UTC corresponds to TOW `422100.0`. Under the corrected logic, GPST TOW is $422100.0 + 18.0 - 10800.0 = 411318.0$. The test was correctly updated to assert `411318.0` instead of `422118.0`.

---

## 3. Caveats

- **Leap Second Hardcoding**: The leap second offset of `18.0` is hardcoded. While correct for modern epochs (since Jan 1, 2017), older RINEX files would theoretically have different leap second values. However, since the leap second offset was already hardcoded as `18.0` in the original parser implementation, this is an inherited limitation rather than a new regression.

---

## 4. Conclusion

The GLONASS time scale discrepancy bug has been successfully resolved. GLONASS ephemeris epochs are now correctly converted to GPST by subtracting the 3-hour Moscow Time offset and properly normalising the `GpsTime` struct.

---

## 5. Verification Method

To independently verify:
- Inspect the file `crates/gneiss-parsers/src/rinex.rs`.
- Run `cargo test -p gneiss-parsers` to verify parser tests.
- Run `cargo test --workspace` to verify workspace integrity.

---

## 6. Quality Review

**Verdict**: **APPROVE**

### Findings
- **No findings of concern**: The implementation is clean, correct, and follows Rust best practices (utilizing the `Add` trait instead of manual mutable field updates).

### Verified Claims
- GLONASST offset correction $\rightarrow$ verified via `cargo test -p gneiss-parsers` $\rightarrow$ **PASS**
- GPST normalization $\rightarrow$ verified via `GpsTime::normalize` implementation analysis $\rightarrow$ **PASS**
- Unit test correctness $\rightarrow$ verified via `test_parse_rinex_3_nav_date` code review $\rightarrow$ **PASS**

---

## 7. Adversarial Review

**Overall risk assessment**: **LOW**

### Challenges

#### [Low] Challenge 1: Hardcoded Leap Seconds
- **Assumption challenged**: GLONASS-to-GPS conversion offset relies on a static leap second value of `18.0`.
- **Attack scenario**: Parsing an old RINEX file (e.g. from 2010 when leap seconds were 15) or a future RINEX file (if leap seconds increase).
- **Blast radius**: The parsed TOW will be off by a few seconds compared to the true GPST.
- **Mitigation**: A future refactor could extract leap seconds dynamically from the RINEX header if present, or use a lookup table based on epoch. For the scope of this bug fix, keeping the static `18.0` value matches the project's existing design.

### Stress Test Results
- Negative time addition $\rightarrow$ `toc_gpst = toc_gpst + (18.0 - 10800.0)` $\rightarrow$ normalizes `tow` and decrements `week` $\rightarrow$ **PASS**
