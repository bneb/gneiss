# Task Assignment — Worker M2: Code-Minus-Carrier (CMC) Multipath Detection & Down-Weighting

## Role & Mission
You are an implementation worker (`teamwork_preview_worker`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon/`.
Your parent orchestrator is: `5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0`.

## Authoritative Documents to Read Before Starting Work
1. `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (authoritative user requirements, especially Section R2)
2. `/Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon_gen2/PROJECT.md` (project architecture, interface contracts, and code layout)
3. `/Users/kevin/projects/gneiss/.agents/survey_explorer_2/survey_r2_cmc.md` (comprehensive CMC analysis, mathematical derivations, code inspection, and implementation specifications)
4. `/Users/kevin/projects/gneiss/.agents/survey_explorer_2/handoff.md` (summary of decoupling, CMC tracking, and MW shielding)
5. `/Users/kevin/projects/gneiss/.agents/worker_m1/handoff.md` (Milestone 1 completed covariance and SNR extraction)
6. `/Users/kevin/projects/gneiss/AGENTS.md` (strict code quality standards)

## MANDATORY INTEGRITY WARNING
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

## Code Quality Standards (AGENTS.md)
- File size strictly < 500 LOC
- Function size strictly <= 32 LOC
- Nesting depth strictly < 3 levels
- Exactly 0 `unwrap()` calls in production code (`match`, `if let`, `ok_or()?`, or descriptive `.expect()` only; `unwrap()` in `#[cfg(test)]` is acceptable)
- Zero compiler warnings, zero clippy warnings (`cargo clippy --workspace --all-targets -- -D warnings`)
- All tests pass (`cargo test --workspace`)
- Both CI smoke guard scripts pass:
  - `python3 scripts/check_network_benchmark.py --smoke`
  - `python3 scripts/check_multignss_benchmark.py --smoke`

## Exclusive Write Ownership
You own ONLY:
- `crates/gneiss-rtk/src/estimators/rtk_iekf/screen.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/update/robust.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/formation_cov.rs` (if auxiliary helpers are needed)
Do NOT edit any files outside this set without coordinating.

## Scope & Implementation Objectives (Milestone 2: R2)

1. **Decouple Code Gross Error Screen (`screen.rs`)**:
   - Currently, `screen_gross_pr_errors` in `screen.rs` calls `measurements.remove(idx)` when pseudorange residual exceeds `GROSS_PR_ERROR_THRESHOLD_M` (15m), which discards the entire `DoubleDiffMeasurement` including carrier phase `dd_cp_cycles`!
   - Modify the screening logic so that code gross errors do NOT drop the carrier phase measurement. Instead, mark the code measurement as suppressed / de-weighted (e.g., set `pr_suppressed = true` or inflate `pr_var_m2` to an extreme de-weighting value such as `1e6` or mark code weight to 0) while keeping `dd_cp_cycles` active with nominal phase variance.
   - Ensure the function returns keys of de-weighted / suppressed code measurements or appropriately tracks them.

2. **Code-Minus-Carrier (CMC) Multipath Detection & Adaptive Down-Weighting (`formation.rs`, `update/robust.rs`)**:
   - Track Code-Minus-Carrier residuals:
     $$\text{CMC} = P - \Phi - 2I$$
     or dual-frequency geometry-free combinations across continuous carrier tracking arcs.
   - When sudden pseudorange deviations occur without corresponding carrier phase movement along a continuous carrier arc, identify the satellite pair as affected by code multipath / NLOS reflection.
   - Adaptively inflate the pseudorange measurement variance:
     $$R_{PP} \leftarrow R_{PP} + \sigma_{\text{mp}}^2$$
     or down-weight the contaminated observation row in the filter update (`robust.rs`), without inflating carrier phase variance $R_{\Phi\Phi}$.
   - Maintain zero false fixes: ensure multipath de-weighting does not falsely admit corrupted carrier phase tracking into the integer search.

3. **Melbourne-Wübbena (MW) Multipath Shielding (`mw.rs`)**:
   - In `mw.rs`, code multipath steps (e.g. 5m–15m pseudorange jumps) produce wide-lane innovation jumps exceeding `SLIP_INNOVATION_CYCLES = 1.0`, triggering false cycle slip resets in `absorb()`.
   - Protect the running wide-lane filter `MwTrack` from resetting on code multipath jumps along a continuous carrier arc. Gate or shield the MW update when code multipath is detected so wide-lane statistics remain stable across urban canyon blocks.

4. **Unit Tests & Invariants**:
   - Add unit tests in each modified file verifying:
     - Code blunders do not discard carrier phase.
     - Sudden code steps inflate code variance while phase variance remains nominal.
     - MW tracker is shielded from false resets on code multipath jumps.
   - Verify all functions <= 32 LOC and files < 500 LOC.

## Verification Commands
Before delivering your handoff, execute and verify:
```bash
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p gneiss-rtk --lib
cargo test --test test_urban_canyon_e2e
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke
```
Write your full report to `/Users/kevin/projects/gneiss/.agents/worker_m2_urban_canyon/handoff.md` following the Handoff Protocol (Observation, Logic Chain, Caveats, Conclusion, Verification Method).
Then notify parent with a concise message.

## 2026-09-24T18:48:46Z
User prompt received: You are Worker M2 for Milestone 2: Code-Minus-Carrier (CMC) Multipath Detection & Down-Weighting.
Implement F4 (decouple PR gross error screen in screen.rs so carrier phase is retained), F5/F6 (CMC arc tracking & adaptive code de-weighting in formation.rs / robust.rs), and F7 (MW multipath shielding in mw.rs).
