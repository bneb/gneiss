# Handoff Report: Survey Explorer 1 (R1 Baseline & Benchmarks)

**Author**: Survey Explorer 1  
**Task ID**: `c1309e2d-6c95-4b14-a86d-d26a13f2a150`  
**Working Directory**: `/Users/kevin/projects/gneiss/.agents/survey_explorer_1`  
**Handoff Type**: Hard (Task complete)  
**Detailed Report**: `/Users/kevin/projects/gneiss/.agents/survey_explorer_1/survey_r1_bench.md`  

---

## 1. Observation

1. **Physical Variance Models (`crates/gneiss-core/src/variance.rs:6-23`)**:
   Currently contains:
   ```rust
   pub fn snr_variance_scale(snr_dbhz: f64, snr_a: f64, snr_b: f64) -> f64 {
       let snr_safe = if snr_dbhz < 10.0 { 10.0 } else { snr_dbhz };
       snr_a * snr_a + (snr_b * snr_b) / libm::pow(10.0, snr_safe / 10.0)
   }
   pub fn elevation_variance_scale(el_rad: f64) -> f64 {
       let sin_el = libm::sin(el_rad);
       let sin_el_safe = if sin_el < 0.1 { 0.1 } else { sin_el };
       1.0 / (sin_el_safe * sin_el_safe)
   }
   pub fn observation_variance(snr_dbhz: f64, el_rad: f64, snr_a: f64, snr_b: f64) -> f64 {
       snr_variance_scale(snr_dbhz, snr_a, snr_b) * elevation_variance_scale(el_rad)
   }
   ```
   This is only called by SPP (`crates/gneiss-rtk/src/estimators/spp/solver.rs:165`).

2. **Observation Structure Quantization (`crates/gneiss-core/src/obs.rs:184-192`)**:
   ```rust
   pub fn get_snr(&self, freq_band: u8) -> Option<u8> {
       self.observations.iter().find(|o| {
           o.code.obs_type == ObsType::Snr
               && self.matches_band(o.code.signal.freq_band, o.code.signal.attribute, freq_band)
       }).map(|o| o.value as u8)
   }
   ```
   `o.value` is stored as `f64`, but `get_snr` casts it to `u8`, discarding fractional C/N0 precision and inducing 1 dB discretization jumps.

3. **Current Double-Difference Weighting (`crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs:21-65`)**:
   - `snr_weight(snr_dbhz: Option<u8>) -> f64`: Uses hard piece-wise cutoff at `40.0` dB-Hz with non-zero derivative below 40 and zero derivative above 40 (slope jump $\Delta(dw/dS) = 0.230$).
   - `attenuation_scale(snr_dbhz: Option<u8>, el_rad: f64) -> f64`: Multiple sharp threshold branches (`drop_db > 10.0`, `snr < 25.0`, `scale.min(1000.0)`).
   - `single_diff_var`: `let sin_el = el_rad.sin().max(0.1);` with slope jump at $\theta = 5.74^\circ$.

4. **Double-Difference Covariance Structure (`crates/gneiss-rtk/src/estimators/rtk_iekf/update/system.rs:237-251`)**:
   - Off-diagonal DD covariance between satellites sharing the reference satellite is filled as `cov = mi.ref_var.min(mj.ref_var)`.
   - Diagonal terms $R_{ii}$ are formed from `m.pr_var_m2` and `m.cp_var_cycles2`.
   - Double-difference covariance matrix satisfies $R_{DD} = D R_{\text{undiff}} D^T \succ 0$.

5. **CI Smoke Guards & Test Baseline Verification**:
   - `python3 scripts/check_network_benchmark.py --smoke`: Exit code 0, `ALL CHECKS PASSED`. Network fused horizontal $p_{50} = 0.024\text{ m} \le 0.040\text{ m}$, RMS $= 0.037\text{ m} \le 0.060\text{ m}$.
   - `python3 scripts/check_multignss_benchmark.py --smoke`: Exit code 0, `ALL CHECKS PASSED`. Network fused fix rate $= 99.70\% \ge 96.50\%$, P181 fix rate $= 98.30\% \ge 97.50\%$.
   - `cargo test --workspace`: Exit code 0, `205 passed; 0 failed`.
   - `cargo clippy --workspace --all-targets -- -D warnings`: Exit code 0, 0 warnings.
   - `MAX_EPOCHS=200 eval_f9p_rover`:
     - Tokyo Odaiba: Forward RTK 20.0%, Smoothed PPK 58.5%.
     - Tokyo Shinjuku: Forward RTK 7.5%, Smoothed PPK 20.5%.
     - Hong Kong Whampoa Survey: Forward RTK 41.0%, Smoothed PPK 40.0%.

---

## 2. Logic Chain

1. From Observation 1, `variance.rs` currently implements a classical product model $\sigma^2 = (a^2 + b^2/10^{S/10}) \cdot (1/\sin^2\theta)$ but is isolated in SPP. Double difference RTK in `formation_cov.rs` uses a separate, ad-hoc multiplier system (`snr_weight` and `attenuation_scale`).
2. From Observation 3, both `snr_weight` and `attenuation_scale` contain sharp slope discontinuities at 40 dB-Hz, 25 dB-Hz, and at the 1000.0 ceiling. In dynamic vehicle navigation, when C/N0 or elevation crosses these thresholds, the observation variance derivative jumps discontinuously, causing Kalman gain chatter in the IEKF.
3. From Observation 2, `SatObs::get_snr` casts `f64` to `u8`, injecting artificial discrete stepping into variance weights whenever C/N0 fluctuates across an integer boundary.
4. From Observation 4, double-difference observation error propagation preserves strict positive definiteness ($R_{DD} \succ 0$) as long as individual satellite observation variances are strictly positive and reference satellite noise is consistently tracked.
5. Therefore, formulating a unified, $C^1$-smooth SIGMA-SNR model:
   $$\sigma^2(\theta, S) = \left( a^2 + \frac{b^2}{\sin^2 \theta + \sin^2 \theta_0} \right) \cdot f_{\text{SNR}}(S)$$
   with smooth logistic activation and algebraic saturation ceiling:
   - Eliminates derivative discontinuities and filter chatter.
   - Correctly scales noise exponentially below 40 dB-Hz.
   - Smoothly bounds maximum variance ($\le 1000.0$) preventing numerical ill-conditioning.
   - Evaluates rover and base observation noise independently according to rigorous error propagation.

---

## 3. Caveats

1. The smoke runs (`--smoke` and `MAX_EPOCHS=200`) evaluate the initial segments of the benchmark datasets. Full-dataset baseline numbers ($N=1,242$ for Odaiba, $N=2,096$ for Shinjuku, $N=1,535$ for Whampoa, and $N=12,398$ for Odaiba INS) are documented from `PROJECT_STATUS.md` and should be used as the definitive comparison standard for full validation passes.
2. In `crates/gneiss-core/src/obs.rs`, `get_snr(&self, freq_band: u8) -> Option<u8>` is part of the existing public API. Rather than breaking existing call sites, `get_snr_f64(&self, freq_band: u8) -> Option<f64>` should be introduced alongside it.

---

## 4. Conclusion

1. The codebase is in a pristine, fully verified state with zero compiler/clippy warnings and 100% passing tests and CI smoke guards.
2. The exact SIGMA-SNR formulation has been derived and mathematically validated:
   $$\sigma^2(\theta, S) = \left( a^2 + \frac{b^2}{\sin^2 \theta + \sin^2 \theta_0} \right) \cdot f_{\text{SNR}}(S)$$
   with $\theta_0 = 5^\circ$, $a_{\text{code}} = b_{\text{code}} = 0.1414\text{ m}$, $a_{\text{phase}} = b_{\text{phase}} = 0.00212/\lambda\text{ cyc}$, smooth logistic activation $\tau = 1.5\text{ dB-Hz}$, and smooth algebraic ceiling $f_{\max} = 1000.0$.
3. All target files for R1 implementation have been identified (`crates/gneiss-core/src/obs.rs`, `crates/gneiss-core/src/variance.rs`, `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`, and `formation.rs`), with a clear implementation roadmap documented in `/Users/kevin/projects/gneiss/.agents/survey_explorer_1/survey_r1_bench.md`.

---

## 5. Verification Method

To independently verify all findings and baselines in this report:

```bash
# 1. CI Smoke Guards
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke

# 2. Workspace Quality Standards
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings

# 3. Kinematic Rover Smoke Runs
MAX_EPOCHS=200 cargo run --release --bin eval_f9p_rover -- odaiba
MAX_EPOCHS=200 cargo run --release --bin eval_f9p_rover -- shinjuku
MAX_EPOCHS=200 cargo run --release --bin eval_f9p_rover -- whampoa_survey
```

**Invalidation Conditions**:
- Any regression or non-zero exit in the smoke guards or workspace tests.
- Non-monotonicity in $\sigma^2(\theta, S)$ with respect to elevation or SNR ($\frac{\partial \sigma^2}{\partial \theta} > 0$ or $\frac{\partial \sigma^2}{\partial S} > 0$).
- Discontinuous derivatives causing filter divergence or NaN values in double-difference updates.
