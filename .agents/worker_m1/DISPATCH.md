# Dispatch: Worker M1 (R1 Adaptive C/N0 & Elevation Observation Covariance)

## 2026-09-24T13:53:08Z

## Objective
Implement Milestone 1 (R1: Adaptive C/N0 (SNR) and elevation observation covariance weighting) across the Gneiss positioning engine.

## Required Reading
- `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (timestamp 2026-09-24T13:30:49Z)
- `/Users/kevin/projects/gneiss/AGENTS.md`
- `/Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon/PROJECT.md`
- `/Users/kevin/projects/gneiss/.agents/survey_explorer_1/survey_r1_bench.md`

## Mandatory Integrity Warning
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

## Code Ownership & Targets
You exclusively own and may modify:
1. `crates/gneiss-core/src/obs.rs`:
   - Add `pub fn get_snr_f64(&self, freq_band: u8) -> Option<f64>` returning exact floating-point C/N0 in dB-Hz without integer quantization. Keep existing `get_snr` intact for backward compatibility.
2. `crates/gneiss-core/src/variance.rs`:
   - Implement the physically grounded, $C^1$-smooth SIGMA-SNR elevation and noise model:
     $$\sigma^2(\theta, S) = \left( a^2 + \frac{b^2}{\sin^2 \theta + \sin^2 \theta_0} \right) \cdot f_{\text{SNR}}(S)$$
     where $\theta_0 = 5^\circ$ ($0.087266$ rad), smooth logistic activation ($\tau = 1.5\text{ dB-Hz}$ around threshold 40.0 dB-Hz), and smooth algebraic saturation ceiling ($f_{\max} = 1000.0$).
   - Ensure continuous derivatives with respect to elevation and SNR ($\frac{\partial \sigma^2}{\partial \theta} \le 0$ and $\frac{\partial \sigma^2}{\partial S} \le 0$).
   - Implement unit tests covering the Three-Tier verification standard (analytical golden vectors, finite differences, stability bounds).
3. `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`:
   - Refactor observation variance scaling to use `get_snr_f64` and the smooth SIGMA-SNR model.
   - Replace discontinuous piecewise branches (`snr_weight`, `attenuation_scale`, and hard `min(1000.0)`) with smooth formulation.
   - Ensure double-difference covariance maintains strict positive definiteness ($R_{DD} \succ 0$).

## Standards & Constraints
- File size < 500 LOC
- Function size <= 32 LOC
- Nesting depth < 3 levels
- No `unwrap()` in production code (use `?`, `match`, `if let`, or descriptive `.expect()`)
- 0 compiler warnings, 0 clippy warnings (`cargo clippy --workspace --all-targets -- -D warnings`)
- All tests pass: `cargo test --workspace`
- Both CI smoke guard scripts pass:
  - `python3 scripts/check_network_benchmark.py --smoke`
  - `python3 scripts/check_multignss_benchmark.py --smoke`

## Deliverables
- Implementation changes and passing unit tests.
- Deliver `handoff.md` with:
  - Observation
  - Logic Chain
  - Caveats
  - Conclusion
  - Verification Method (with exact build and test command outputs)
- Send completion message to parent orchestrator.

## 2026-09-24T14:30:19Z
**Context**: Orchestrator liveness check for Milestone 1
**Content**: Checking on current status of Milestone 1 implementation. What is the current progress, test/benchmark status, and are you blocked on anything?
**Action**: Please reply with a brief status update and update your progress.md.
