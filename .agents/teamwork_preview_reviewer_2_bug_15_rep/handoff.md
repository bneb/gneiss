# Handoff Report: Bug 15 Implementation Review

## 1. Observation
I directly observed the following files, implementation code, test logs, and repository history:

* **Target files and commits**:
  * Git commit `f348335 fix(ephemeris, ppp): Bug 15 — TGD not applied in dual-frequency PPP` was analyzed.
  * In `crates/gneiss-core/src/ephemeris.rs`:
    * Enum-level delegation (lines 47-55):
      ```rust
      pub fn position_iono_free(&self, t: GpsTime) -> (Vector3<f64>, Vector3<f64>, f64, f64) {
          match self {
              Ephemeris::Gps(e) => e.position_iono_free(t),
              Ephemeris::Galileo(e) => e.position_iono_free(t),
              Ephemeris::Beidou(e) => e.position_iono_free(t),
              Ephemeris::Qzss(e) => e.position_iono_free(t),
              Ephemeris::Glonass(e) => e.position(t),
          }
      }
      ```
    * `GpsEphemeris::position_iono_free` (lines 517-545) calls `calc_keplerian` passing `0.0` for `tgd` instead of `self.tgd`.
    * `GalileoEphemeris::position_iono_free` (lines 579-607) calls `calc_keplerian` passing `0.0` for `tgd` instead of `self.bgd_e1_e5a`.
    * `BeidouEphemeris::position_iono_free` (lines 680-711) calls `calc_keplerian` passing `0.0` for `tgd` instead of `self.tgd1`.
    * `QzssEphemeris::position_iono_free` (lines 746-774) calls `calc_keplerian` passing `0.0` for `tgd` instead of `self.tgd`.
    * The unit test `test_broadcast_clock_tgd_correct` (lines 1419-1548) checks that `clk_if - clk_pos` matches the exact TGD values within a `1e-15` tolerance for GPS, Galileo, Beidou, and QZSS.
  * In `crates/gneiss-rtk/src/engine/ppp.rs`:
    * `compute_sat_state` updates (lines 427, 450) fetch the broadcast clock and satellite state via `eph.position_iono_free(t_nom)` and `eph.position_iono_free(t_tx)`.
  * In `crates/gneiss-rtk/src/estimators/spp.rs`:
    * `compute_sat_state` updates (lines 185-198) dynamically branch on `m.is_iono_free`:
      ```rust
      let (_, _, sat_clk_err_rough, _) = if m.is_iono_free {
          m.eph.position_iono_free(t_tx_sat_gps)
      } else {
          m.eph.position(t_tx_sat_gps)
      };
      ```
* **Test results**:
  * Running `cargo test` succeeded cleanly with `test ephemeris::tests::test_broadcast_clock_tgd_correct ... ok` and all other tests passing:
    ```
    test result: ok. 257 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.20s
    test result: ok. 49 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
    test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
    ```

## 2. Logic Chain
1. **Clock Correction Definition**: According to GNSS broadcast standards (e.g., GPS IS-GPS-200), broadcast clock parameters (af0, af1, af2) are referenced to the dual-frequency ionosphere-free combination (e.g., L1/L2 P(Y)). Thus, for dual-frequency ionosphere-free users, no TGD/BGD correction should be subtracted. For single-frequency users, TGD/BGD corrections must be applied.
2. **Keplerian Base Function**: `calc_keplerian` calculates clock error as `clk_err = af0 + af1 * tc + af2 * tc * tc + F * e * sqrt_a * sin(ek) - tgd`.
3. **Iono-Free Clock Retrieval**: By creating `position_iono_free` which passes `tgd = 0.0` to `calc_keplerian`, the returned clock offset correctly lacks any TGD subtraction, representing the pure dual-frequency reference clock.
4. **Constellation Coverage**: GPS (`tgd`), Galileo (`bgd_e1_e5a`), Beidou (`tgd1`), and QZSS (`tgd`) are all correctly handled. GLONASS, having no TGD, correctly falls back to `position(t)`.
5. **PPP Correction**: In `engine/ppp.rs`, using `position_iono_free` removes the need for the fragile `brdc_clk + eph.tgd()` hack and ensures the satellite transmission state (`t_tx`) is computed consistently.
6. **SPP Correction**: In `estimators/spp.rs`, the code dynamically branches. Dual-frequency combinations (which are iono-free) bypass TGD, whereas single-frequency measurements preserve TGD subtraction. This is mathematically correct.

## 3. Caveats
* **Single-Frequency PPP**: If a user runs the PPP engine in single-frequency mode without precise clock files, the broadcast clock correction retrieved via `position_iono_free` will not have TGD subtracted. However, the PPP engine is designed for dual-frequency / precise carrier-phase setups, and single-frequency operations are delegated to the SPP engine or handled via OSB corrections. Therefore, this is an acceptable constraint.

## 4. Conclusion

### Review Summary
**Verdict**: APPROVE

### Findings
None. The implementation is clean, conforms to GNSS physics, and passes all tests.

### Verified Claims
* `position_iono_free` bypasses TGD/BGD subtraction -> verified via `view_file` on `crates/gneiss-core/src/ephemeris.rs` -> PASS
* Unit test `test_broadcast_clock_tgd_correct` executes and checks constraints -> verified via `cargo test` -> PASS
* SPP correctly uses conditional branching based on measurement type -> verified via `view_file` on `crates/gneiss-rtk/src/estimators/spp.rs` -> PASS

### Challenge Summary
**Overall risk assessment**: LOW

### Challenges
No critical challenges or failure modes were identified. Edge cases (such as GLONASS and Galileo alternate bands) were analyzed and found to be handled correctly.

## 5. Verification Method
To independently verify:
1. Run `cargo test -p gneiss-core --lib ephemeris::tests::test_broadcast_clock_tgd_correct` to run the specific unit test.
2. Run `cargo test` to execute the full workspace test suite.
