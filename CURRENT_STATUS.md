# Current Status — gneiss

**Report Date:** 2026-10-01T03:38:27-07:00  
**Review Window:** Last 24 Hours (2026-09-30T03:38:27-07:00 to 2026-10-01T03:38:27-07:00)  
**Status:** Active Development / Scientific Verification  
**Current Branch:** `main` (clean working tree)  
**Latest Tag:** `checkpoint/sprint-61`

---

## 1. Executive Summary

Over the last 24 hours, `gneiss` (a GNSS post-processing engine in Rust with double-difference EKF and LAMBDA integer ambiguity resolution) completed 32 commits. The sprint resolved a critical GLONASS ephemeris propagation defect, introduced statistical rigor to RTK evaluations via paired bootstrap CDF bands and Weibull tail models, resolved the zero-integer-fix anomaly on u-blox F9P receivers, and uncovered that urban canyon positioning errors are caused by bimodal convergence rather than progressive drift.

---

## 2. Developments & Progress (Last 24 Hours)

### GLONASS Ephemeris & Residual Defect Resolution (`e985330`, `1251d7a`, `aedeffc`, `3afa3c2`)
- **Root Cause Identified:** Stopped propagating GLONASS broadcast ephemerides as Keplerian orbital elements (GLONASS uses PZ-90 coordinate systems requiring Runge-Kutta numerical integration).
- **Residual Collapse:** Eliminated catastrophic double-difference code residuals (which previously spiked up to 20 km), resolving the mystery of 0% integer fix rates on u-blox F9P hardware where code screens had been vetoing every fix.

### Statistical Rigor & Evaluation Tooling (`037e53d`, `8555eeb`, `d4e8f5c`, `e51f9e2`)
- **`eval_compare` Tool:** Built dedicated CLI utility to evaluate GNSS positioning solutions.
- **Paired Bootstrap CDF Bands:** Implemented `bootstrap_cdf_band_paired` in `gneiss_core::stats` using deterministic xorshift PRNG to calculate 95% confidence intervals on quantile differences across shared epochs.
- **Empirical AR Trade-off:** Proved that on u-blox F9P, Integer Ambiguity Resolution (AR) crosses over at $\sim p70$: AR provides statistically significant improvements in bulk accuracy ($p10$ to $p50$), but degrades the extreme tail ($p95$) compared to pure float solutions under multipath.

### Urban Canyon Dynamics & NLOS Mitigation (`47b6644`, `d822a27`, `028f6b3`, `31dee7e`)
- **Bimodal Convergence:** Analyzed Shinjuku canyon datasets; refuted progressive drift hypotheses and proved that tail degradation is caused by bimodal convergence (the estimator locks onto a wrong integer mode and holds it for $\sim 1$ minute with gross single-epoch jumps).
- **Signal Predictor:** Identified that ambiguity covariance states predict wrong-mode transitions.
- **Median-Relative NLOS Gate:** Added per-epoch median-relative pseudorange residual screening to collapse urban canyon tail errors.

---

## 3. Current Verification & Test Status

- **Working Tree:** Clean.
- **Clippy:** `cargo clippy --workspace --all-targets` passes with 0 warnings.
- **Compilation:** `cargo build --release --bin eval_compare` compiles cleanly.

---

## 4. Blockers & Next Steps

1. Integrate the ambiguity state signal into an active integer validation veto to prevent wrong-mode convergence in urban environments.
2. Advance Tier-1 PPK/PPP roadmap benchmarks (Sprint 62).
