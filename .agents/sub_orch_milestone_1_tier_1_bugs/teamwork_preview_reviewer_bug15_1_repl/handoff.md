# Handoff Report — Bug 15 Verification (Incorrect Broadcast Clock TGD Correction)

## 1. Observation

Direct observations and file inspections on `/Users/kevin/projects/gneiss` at commit `420ed89`:

1. **`crates/gneiss-core/src/ephemeris.rs`**:
   - Added `position_iono_free` to `Ephemeris` (lines 47-55) and to variants (e.g. `BeidouEphemeris` at line 688).
   - In `BeidouEphemeris::position_iono_free`, `0.0` is passed as the group delay parameter to `calc_keplerian`:
     ```rust
     pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
         let t_bdt = GpsTime::new(t.week, t.tow - 14.0);
         let is_bds_geo = self.sat.prn <= 5 || self.sat.prn >= 59;
         calc_keplerian(
             t_bdt, self.toe, self.toc, self.af0, self.af1, self.af2,
             self.crs, self.crc, self.cuc, self.cus, self.cic, self.cis,
             self.m0, self.e, self.sqrt_a, self.delta_n, self.omega0, self.omega_dot,
             self.i0, self.idot, self.omega,
             0.0, // Pass 0.0 to avoid subtracting TGD
             MU_BDS, OMEGA_E_BDS, is_bds_geo,
         )
     }
     ```
   - In `calc_keplerian` (line 482):
     ```rust
     let tc = t - toc;
     ```
     For BeiDou, `t` is `t_bdt` (in BDT scale) and `toc` is `self.toc` (in GPST scale), introducing a 14-second offset discrepancy.

2. **`crates/gneiss-parsers/src/rinex.rs`**:
   - Parses `tgd1`, `tgd2`, `aode`, and `aodc` from the BeiDou navigation records (lines 553-579):
     ```rust
     aode: vals[0] as u32,
     ...
     tgd1: vals[22],
     tgd2: vals[23],
     aodc: vals[25] as u32,
     ```

3. **`crates/gneiss-rtk/src/estimators/spp.rs`**:
   - Maps `tgd2: 0.0` inside tests (line 1109) and initializes it in struct construction.

4. **Verbatim Compilation Error in `crates/gneiss-rtk/src/estimators/ekf/filter.rs`**:
   - Proposing `cargo test --workspace` outputs:
     ```
     error: unexpected closing delimiter: `}`
        --> crates/gneiss-rtk/src/estimators/ekf/filter.rs:1767:1
         |
     1752 |     fn test_compute_candidate_variance_high_var_filtered() {
          |                                                            - this opening brace...
     ...
     1766 |     }
          |     - ...matches this closing brace
     1767 | }
          | ^ unexpected closing delimiter
     ```

---

## 2. Logic Chain

1. **Syntax Mismatch**: The EKF test code in `filter.rs` has an extra closing brace `}` at line 1767. Since the `mod tests` block was already closed at line 1413, the braces at the end of the file do not match, causing `gneiss-rtk` compilation to fail.
2. **BeiDou Clock Reference**: The broadcast clock parameters ($a_{f0}, a_{f1}, a_{f2}$) for BeiDou are referenced to the B3I frequency.
3. **BeiDou Combined Group Delay Correction**: In dual-frequency iono-free mode using B1I and B2I, the combined group delay $TGD_{IF(B1I/B2I)} = \frac{f_1^2 \cdot TGD_1 - f_2^2 \cdot TGD_2}{f_1^2 - f_2^2}$ must be subtracted from the clock correction to yield mathematically correct clock error values. 
4. **Omission of Correction**: Because `BeidouEphemeris::position_iono_free` passes `0.0` for `tgd`, it does not subtract the combined group delay, leaving a systematic bias of several meters in Beidou-only dual-frequency SPP/PPP calculations.
5. **14-Second Discrepancy**: Since `self.toc` is stored in GPST, subtracting it from `t_bdt` (which is in BDT, 14 seconds behind) introduces a 14-second clock offset error in Keplerian clock error evaluation.

---

## 3. Caveats

- We assumed standard BeiDou ICD rules where B3I is the reference frequency for clock parameters.
- We did not attempt to fix the compile error in `filter.rs` as it is out-of-scope for the Reviewer role (Review-only constraint).

---

## 4. Conclusion

- The implementation of Bug 15 has **critical compilation and major mathematical issues** and cannot be approved as-is.
- **Verdict**: FAIL / REQUEST_CHANGES

---

## 5. Verification Method

To verify the compilation failure, run:
```bash
cargo build --workspace
```
To verify that `gneiss-core` and `gneiss-parsers` pass their tests individually (excluding EKF compiler issues), run:
```bash
cargo test -p gneiss-core
cargo test -p gneiss-parsers
```
To inspect the code structure, examine `crates/gneiss-rtk/src/estimators/ekf/filter.rs` around lines 1410-1420 and lines 1760-1768.
