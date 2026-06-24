# Quality Review Report — Bug 15 Verification

## Review Summary

**Verdict**: REQUEST_CHANGES (due to EKF syntax compile failure in `filter.rs` and mathematical omission in Beidou combined group delay).

*Note: While the parser mappings and Gps/Galileo/Qzss iono-free clock corrections are correctly implemented and verified, the `gneiss-rtk` package fails to compile on the main branch due to a syntax brace mismatch, and the Beidou combined group delay calculation for dual-frequency B1I/B2I iono-free combination is omitted/incorrect.*

## Findings

### [Critical] Finding 1: Syntax Mismatch causing Compilation Failure in `gneiss-rtk`
- **What**: The EKF estimator code `crates/gneiss-rtk/src/estimators/ekf/filter.rs` fails to compile due to an unexpected closing delimiter.
- **Where**: `crates/gneiss-rtk/src/estimators/ekf/filter.rs:1767:1`
- **Why**: The `mod tests` block at line 523 is closed at line 1413. However, additional tests were added at the end of the file (lines 1676-1766) and closed by a trailing `}` on line 1767. This extra closing brace causes a compiler error since it is outside any open block.
- **Suggestion**: Re-open the `mod tests` block or move the new tests before the closing brace at line 1413.

### [Major] Finding 2: Omission of Combined Group Delay for Beidou B1I/B2I Iono-Free Mode
- **What**: The combined group delay for Beidou B1I/B2I is not calculated or applied in the dual-frequency iono-free satellite clock correction.
- **Where**: `crates/gneiss-core/src/ephemeris.rs` line 688 (`BeidouEphemeris::position_iono_free`)
- **Why**: BeiDou broadcast clock parameters ($a_{f0}, a_{f1}, a_{f2}$) are referenced to the B3I frequency. When performing dual-frequency iono-free processing on B1I and B2I, the combined group delay $TGD_{IF} = \frac{f_1^2 \cdot TGD_1 - f_2^2 \cdot TGD_2}{f_1^2 - f_2^2}$ must be subtracted from the broadcast clock to obtain the correct B1I/B2I combination clock correction. Currently, `position_iono_free` passes `0.0` as the group delay, leaving a systematic bias of several meters (corresponding to the combined group delay of B1I/B2I).
- **Suggestion**: Compute $TGD_{IF}$ in `position_iono_free` for Beidou using B1I and B2I frequencies, and pass it as the group delay parameter to `calc_keplerian` instead of `0.0`.

### [Minor] Finding 3: 14-Second Discrepancy in `tc` Calculation for Beidou Ephemeris
- **What**: The clock term evaluation time `tc` is calculated with a 14-second offset for Beidou.
- **Where**: `crates/gneiss-core/src/ephemeris.rs` line 482 (`let tc = t - toc;` inside `calc_keplerian`)
- **Why**: `toc` for Beidou is stored in GPS Time (GPST) inside the parser, but `t` passed to `calc_keplerian` for Beidou is `t_bdt` (BDT, which is 14 seconds behind GPST). The subtraction `t_bdt - toc_gpst` results in a 14-second error in `tc`, causing minor sub-millimeter to centimeter errors in the clock drift correction.
- **Suggestion**: Convert `toc` to BDT before calculating `tc` in `calc_keplerian`, or pass BDT-referenced `toc` for Beidou.

## Verified Claims

- **Parser mapping for Beidou fields (`tgd1`, `tgd2`, `aode`, `aodc`)** → verified via code inspection and `test_rinex_3_nav_beidou` → **PASS**
- **Galileo E5b clock correction logic** → verified via code inspection and `test_galileo_bgd_e5b` and `test_galileo_position_e5b_regression` → **PASS**
- **GPS / QZSS iono-free clock correction logic** → verified via code inspection and `test_broadcast_clock_tgd_correct` → **PASS**

## Coverage Gaps

- **No test case for Beidou dual-frequency combined group delay** — risk level: High — recommendation: Add a regression test verifying that `position_iono_free` for Beidou applies the correct $TGD_{IF(B1I/B2I)}$ combined delay instead of `0.0`.

## Unverified Items

- **Verification of full workspace test execution (`cargo test --workspace`)** — reason not verified: `gneiss-rtk` compile error prevents running workspace-wide tests. Only `gneiss-core` (102 tests) and `gneiss-parsers` (104 tests) were successfully executed.
