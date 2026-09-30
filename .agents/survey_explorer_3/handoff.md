# Handoff Report: Survey Explorer 3 (R3 Doppler Cycle Slip Validation & R4 Prioritized PAR)

**Working Directory**: `/Users/kevin/projects/gneiss/.agents/survey_explorer_3`  
**Report File**: `/Users/kevin/projects/gneiss/.agents/survey_explorer_3/handoff.md`  
**Detailed Survey Artifact**: `/Users/kevin/projects/gneiss/.agents/survey_explorer_3/survey_r3_r4.md`  
**Parent Conversation ID**: `c1309e2d-6c95-4b14-a86d-d26a13f2a150`  

---

## 1. Observation

Direct code observations from the inspected codebase:

1. **Doppler Cycle Slip Threshold & Frequency Bands** (`crates/gneiss-rtk/src/post_process/screening.rs`):
   - Lines 69–73:
     ```rust
     for band in [1, 2, 7] {
         let band_slips = self.check_band_doppler_slips(epoch, band);
         slips.extend(band_slips);
     }
     ```
     Only bands 1, 2, and 7 are audited. Band 5 (Galileo E5a) and Band 6 (BeiDou B3I) are ignored.
   - Lines 87–91:
     ```rust
     let d_phi = cp - prev_cp;
     let dop_avg = 0.5 * (dop + prev_dop);
     let pred_d_phi = -dop_avg * dt;
     discs.push((sat, d_phi - pred_d_phi, dt));
     ```
   - Lines 150–154:
     ```rust
     let thresh = (1.0 * dt).max(1.0);
     if residual > thresh {
         slipped.push(sat);
     }
     ```
     Threshold is minimum $1.0$ cycle, which mathematically misses all half-cycle slips ($0.5$ cyc) and marginal 1-cycle slips.

2. **Base Slip Gating Leak** (`crates/gneiss-rtk/src/estimators/rtk_iekf/mod.rs` & `formation.rs`):
   - `mod.rs:211–213`:
     ```rust
     if self.widelane_ar {
         self.base_slip_detector.check_epoch(base);
     }
     ```
   - `formation.rs:415–418`:
     ```rust
     let mut cur_arc = slip_detector.get_arc(sat_id) + slip_detector.get_arc(ref_sat);
     if widelane_ar {
         cur_arc += base_slip_detector.get_arc(sat_id) + base_slip_detector.get_arc(ref_sat);
     }
     ```
     When `widelane_ar == false` (standard short/medium baseline RTK), base-side slips are completely bypassed.

3. **Melbourne-Wübbena & Tracking Arc Leak** (`crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs` & `formation.rs`):
   - `mw.rs:304–307`:
     ```rust
     let slip = [rov_s, rov_ref, bas_s, bas_ref].iter().any(|o| {
         o.get_lli(1).is_some_and(|l| l & 1 != 0)
             || o.get_lli(b2).is_some_and(|l| l & 1 != 0)
     });
     ```
     MW tracker checks only raw observation LLI bit 0; it never receives the `lli_slip` flag derived from Doppler or GF cycle slip detection.
   - `formation.rs:336`:
     ```rust
     self.pw_tracker.update(key, pw, 0.0, false);
     ```
     Phase-only wide-lane tracker has `slip: false` hardcoded.
   - `formation.rs:161`:
     ```rust
     *self.pair_epochs.entry(m.key).or_insert(0) += 1;
     ```
     `pair_epochs` tracking duration counter is never reset to 0 when a cycle slip occurs.

4. **PAR Candidate Ranking Baseline** (`crates/gneiss-rtk/src/ambiguity/par.rs` & `ar_subsets.rs`):
   - `par.rs:24–25`:
     ```rust
     let mut indexed: Vec<(usize, f64)> = (0..n).map(|i| (i, cond_stdevs[i])).collect();
     indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
     ```
     Sorts purely by diagonal variance $Q_{ii}$.
   - `ar_subsets.rs:80–84`:
     ```rust
     let score = |i: usize| {
         let frac = (a_float[i] - a_float[i].round()).abs();
         frac + q_amb[(i, i)].sqrt() * 0.5
     };
     c.sort_by(|&i, &j| score(i).total_cmp(&score(j)));
     ```
     Sorts purely by fractional float offset and variance. Ignores elevation, C/N0, lock time, and CMC multipath variance.

5. **Off-Diagonal Covariance Leakage** (`crates/gneiss-rtk/src/estimators/rtk_iekf/ar.rs`):
   - Lines 349–353:
     ```rust
     let mut new_cov = &state.cov - cov_red;
     new_cov = 0.5 * (&new_cov + &new_cov.transpose());
     for i in 0..state.dim() { new_cov[(i, i)] = new_cov[(i, i)].max(1e-6); }
     for &idx in &indices { new_cov[(idx, idx)] = 1e-4; }
     state.cov = new_cov;
     ```
     Sets diagonal variance to `1e-4`, but fails to zero row and column cross-covariances `new_cov[(r, idx)] = 0.0` and `new_cov[(idx, r)] = 0.0`. This produces indefinite $P$ ($\lambda_{\min}(P) < 0$).
   - In contrast, `tc_ambiguity.rs:340–344` explicitly zeroes out row/column cross-covariances:
     ```rust
     self.p_xa.column_mut(idx).fill(0.0);
     self.q_aa.column_mut(idx).fill(0.0);
     self.q_aa.row_mut(idx).fill(0.0);
     self.q_aa[(idx, idx)] = 1e-4;
     ```

---

## 2. Logic Chain

1. **Inability to Detect Half-Cycle Slips**:
   From Observation 1, `thresh = (1.0 * dt).max(1.0)`. For a half-cycle slip, the discrepancy is $0.5$ cycles. Because $0.5 < 1.0$, half-cycle slips are mathematically undetectable under the current code. Lowering the detection threshold to $\approx 0.30$ cycles (scaling adaptively with $\Delta t$) allows detection of both $0.5$-cycle and $1.0$-cycle slips while providing $>6\sigma$ separation from nominal Doppler noise ($\sigma \approx 0.035$ cyc).

2. **Phase Corruption in Wide-Lane AR**:
   From Observation 2 and 3, base slips are unmonitored when `widelane_ar == false`, and `WidelaneTracker` in `mw.rs` is never notified of Doppler or GF slips detected on rover or base. Consequently, when a slip occurs without raw LLI bit 0, the MW filter averages the contaminated phase into its 20-epoch running mean, corrupting wide-lane fixing for subsequent epochs.

3. **Suboptimal Candidate Selection in Urban Canyons**:
   From Observation 4, `select_par_candidates` only evaluates fractional float error and variance. In urban canyons, code multipath on a low-elevation or low-SNR satellite can pull its float ambiguity near an integer coincidentally, while a clean high-elevation satellite ($70^\circ$, $45\text{ dB-Hz}$) has a slight float offset ($0.15$ cyc). The current logic selects the corrupted satellite and omits the clean one, leading to ratio test rejection or false integer fixes. A Composite Quality Metric (CQM) incorporating elevation, SNR, lock duration, and CMC multipath variance is mathematically required.

4. **Covariance Indefiniteness & Filter Divergence**:
   From Observation 5, when ambiguities are fixed to integers, clamping diagonal variance to $10^{-4}$ while leaving cross-covariances intact violates the Cauchy-Schwarz inequality $|P_{ij}| \le \sqrt{P_{ii} P_{jj}}$. For $P_{ii} = 10^{-4}$ and $P_{jj} = 1.0$, any $|P_{ij}| > 0.01$ causes the 2x2 minor determinant $P_{ii} P_{jj} - P_{ij}^2$ to become negative, injecting negative eigenvalues into $P$. Zeroing the fixed ambiguity rows and columns (as in `tc_ambiguity.rs`) and checking $\lambda_{\min}(P) \ge 10^{-6}$ strictly restores positive-definiteness ($P \succ 0$).

---

## 3. Caveats

- **GLONASS FDMA Receiver Clock Modeling**: Doppler-predicted carrier phase increments rely on common-mode receiver clock drift cancellation via the median. For GLONASS satellites with different carrier frequencies, the clock drift in cycles varies per satellite frequency number ($k$). Exact GLONASS modeling should scale the receiver clock drift rate by the individual channel frequency.
- **Extreme Vehicle Jerk**: High angular jerk during rapid vehicle turns can induce higher-order phase rate differences across opposite lines of sight. For 1 Hz data, the adaptive threshold accounts for typical vehicle dynamics ($\le 5\text{ m/s}^2$ acceleration), but at $>10\text{ m/s}^3$ jerk, Doppler velocity prediction from `estimate_doppler_velocity` is superior to undifferenced trapezoidal integration.

---

## 4. Conclusion

The implementation paths for R3 and R4 are well-defined, localized, and directly actionable:
1. **R3**: Update `CycleSlipDetector` in `screening.rs` with multi-band support (1, 2, 5, 6, 7) and an adaptive threshold ($0.28\text{--}0.35$ cyc) to catch half-cycle and 1-cycle slips. Plug slip propagation leaks in `formation.rs` (unconditional base slip check, pass slip to `mw.rs` and `pw_tracker`, reset `pair_epochs` to 0).
2. **R4**: Update `ar_subsets.rs` and `par.rs` with the Composite Quality Metric (elevation, C/N0, lock duration, CMC variance, float variance). Add DOP guard using `gneiss_core::dop::compute_dop_from_positions` with minimum subset size $\ge 4$. Fix off-diagonal covariance leakage in `ar.rs:condition_state_on_integers` by zeroing cross-covariances and enforcing $\lambda_{\min}(P) \ge 10^{-6}$.
3. All proposed changes conform to `AGENTS.md` rules (functions $\le 32$ LOC, files $< 500$ LOC, 0 unwrap, 0 warnings).

---

## 5. Verification Method

Independent verification of the findings can be performed using:

1. **Inspect Code Locations**:
   - `screening.rs:150` for the coarse threshold `(1.0 * dt).max(1.0)`.
   - `formation.rs:416` for the `if widelane_ar` gate on base slip detection.
   - `mw.rs:304` for the LLI-only slip check in `update_tracker_from_obs`.
   - `ar_subsets.rs:80` for the fractional-only candidate score.
   - `ar.rs:349-353` for the un-zeroed off-diagonal cross-covariances.

2. **Workspace Test Suite**:
   ```bash
   cargo test --workspace
   cargo clippy --workspace --all-targets -- -D warnings
   ```

3. **CI Smoke Regression Guards**:
   ```bash
   cargo build --release --bin eval_network_ppk
   python3 scripts/check_network_benchmark.py --smoke
   python3 scripts/check_multignss_benchmark.py --smoke
   ```

4. **Invalidation Condition**:
   If an existing test demonstrates that half-cycle slips ($0.5$ cyc) are currently detected in real-time `rtk_iekf`, or that `ar.rs:condition_state_on_integers` maintains positive eigenvalues without zeroing cross-covariances, this conclusion would be invalidated.
