# BRIEFING — 2026-09-24T13:33:14Z

## Mission
Investigate and design R2: Code-Minus-Carrier (CMC) Multipath Detection & Down-Weighting for Gneiss RTK positioning engine.

## 🔒 My Identity
- Archetype: explorer
- Roles: survey, investigation, synthesis
- Working directory: /Users/kevin/projects/gneiss/.agents/survey_explorer_2
- Original parent: c1309e2d-6c95-4b14-a86d-d26a13f2a150
- Milestone: Urban Canyon Fix Rate Expansion and Multipath Mitigation - R2 CMC Multipath Mitigation

## 🔒 Key Constraints
- Read-only investigation — do NOT implement
- File size < 500 LOC, function size < 32 LOC, nesting depth < 3
- 0 compiler warnings, 0 clippy warnings
- No unwrap() in production code
- Frame safety & relational coupling invariants

## Current Parent
- Conversation ID: c1309e2d-6c95-4b14-a86d-d26a13f2a150
- Updated: 2026-09-24T13:33:14Z

## Investigation State
- **Explored paths**:
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs`
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/update/robust.rs`
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/ar_gate.rs`
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/screen.rs`
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/update/system.rs`
  - `crates/gneiss-rtk/src/estimators/rtk_iekf/update/kalman.rs`
  - `crates/gneiss-rtk/src/measurements/multipath.rs`
  - `crates/gneiss-rtk/src/post_process/screening.rs`
- **Key findings**:
  - **DD Formation**: Forms independent code and phase rows in $H$ and $R$ (`system.rs`). Off-diagonal covariances exist strictly within the same measurement kind. Inflating code variance decouples completely from phase precision.
  - **Tracking Arcs**: Bounded by cycle slips (`screening.rs` GF 5cm, Doppler, LLI, and gaps). Along continuous arcs, carrier ambiguity $\lambda N$ is invariant.
  - **Critical Defect in `screen.rs`**: Prefit gross error check (`screen_gross_pr_errors`) drops the entire `DoubleDiffMeasurement` on $>15\text{m}$ PR residual, throwing away clean carrier phase.
  - **Vulnerability in `mw.rs`**: 10m PR multipath jump causes a 6.5-cycle step in $MW$, tripping `SLIP_INNOVATION_CYCLES = 1.0` in `MwTrack::absorb`, falsely resetting wide-lane integer tracks.
  - **CMC Formulation**: Dual-frequency $MP_1 = P_1 - (1+\beta)\Phi_1 + \beta\Phi_2$ eliminates geometry, clocks, troposphere, and ionosphere. Single-frequency $CMC_1 = P_1 - \lambda_1 \Phi_1$ drifts smoothly ($<5\text{mm/s}$) but detects 5m–20m code multipath steps with $10\sigma\text{--}50\sigma$ sensitivity.
  - **Carrier Preservation**: Multipath variance $\sigma_{\text{mp}}^2$ inflates only $pr\_var\_m2$; carrier $cp\_var\_cycles2$ stays at nominal mm-level. Code row is down-weighted or omitted; phase row retains full weight.
- **Unexplored areas**:
  - Full implementation details of R1 (SNR weighting), R3 (Doppler cycle slips), and R4 (PAR prioritization) being handled by peer explorers.

## Key Decisions Made
- Completed deep dive into RTK double-difference engine.
- Formulated the exact two-stage CMC multipath detection and down-weighting algorithm.
- Documented findings in `survey_r2_cmc.md`.

## Artifact Index
- `/Users/kevin/projects/gneiss/.agents/survey_explorer_2/survey_r2_cmc.md` — Comprehensive findings on R2 CMC
- `/Users/kevin/projects/gneiss/.agents/survey_explorer_2/handoff.md` — 5-component handoff report

