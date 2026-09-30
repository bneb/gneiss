# Task Assignment — Worker M3: Doppler-Assisted Cycle Slip Detection & Phase Continuity Validation

## Role & Mission
You are an implementation worker (`teamwork_preview_worker`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/worker_m3_urban_canyon/`.
Your parent orchestrator is: `5fc6ee4b-4008-46c0-b1a5-eb5f599dc4d0`.

## Authoritative Documents to Read Before Starting Work
1. `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md` (authoritative user requirements, Section R3)
2. `/Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon_gen2/PROJECT.md` (project specifications, interface contracts, and code layout)
3. `/Users/kevin/projects/gneiss/.agents/survey_explorer_3/survey_r3_r4.md` (detailed R3 mathematical derivations, code inspection, and implementation specifications)
4. `/Users/kevin/projects/gneiss/.agents/survey_explorer_3/handoff.md` (summary of Doppler slip detection and slip propagation leaks)
5. `/Users/kevin/projects/gneiss/AGENTS.md` (strict code quality standards)

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
- `crates/gneiss-rtk/src/post_process/screening.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/mw.rs`
- `crates/gneiss-rtk/src/estimators/rtk_iekf/mod.rs` (minimal edit: unconditional base slip check)
- `crates/gneiss-rtk/src/estimators/doppler.rs` (if auxiliary Doppler helpers are needed)
Do NOT edit any files outside this set without coordinating.
Keep all files strictly < 500 LOC (note `mod.rs` is at 482 LOC, `mw.rs` is at 485 LOC, `formation.rs` is at 452 LOC, `screening.rs` is at 425 LOC).

## Scope & Implementation Objectives (Milestone 3: R3)

1. **Multi-Band Doppler Slip Detection & Adaptive Threshold (`screening.rs`)**:
   - In `CycleSlipDetector`:
     - Extend audited frequency bands from `[1, 2, 7]` to `[1, 2, 5, 6, 7]` in `check_doppler_slips` to cover GPS L1/L2/L5, Galileo E1/E5a/E5b, and BeiDou B1/B2a/B3I.
     - Replace the coarse threshold `(1.0 * dt).max(1.0)` with a tuned adaptive threshold:
       ```rust
       let thresh = (0.30 * dt).clamp(0.28, 1.0);
       ```
       This detects half-cycle ($0.5$ cyc) and 1.0-cycle slips while keeping $> 6\sigma$ margin against nominal Doppler measurement noise ($\sigma \approx 0.035$ cyc).

2. **Unconditional Base Slip Detection (`mod.rs` & `formation.rs`)**:
   - In `crates/gneiss-rtk/src/estimators/rtk_iekf/mod.rs`:
     - Currently `self.base_slip_detector.check_epoch(base)` only runs `if self.widelane_ar`.
     - Remove the `if self.widelane_ar` condition so base slip detection runs unconditionally for all RTK modes.
   - In `crates/gneiss-rtk/src/estimators/rtk_iekf/formation.rs`:
     - In `check_pair_slip` / arc tracking calculation, query `base_slip_detector.get_arc(sat_id) + base_slip_detector.get_arc(ref_sat)` unconditionally (not gated behind `widelane_ar`).

3. **Plug Slip Propagation Leaks (`formation.rs`, `mw.rs`)**:
   - In `formation.rs`:
     - When a slip occurs on a satellite pair (`m.slip == true` or `check_pair_slip` returns a slip):
       - Reset `pair_epochs` tracking duration counter to 0 (`self.pair_epochs.insert(m.key, 0)`).
       - Ensure `self.pw_tracker.update(key, pw, 0.0, is_slip)` receives the actual slip flag rather than hardcoded `false`.
   - In `mw.rs`:
     - In `update_tracker_from_obs`, ensure the wide-lane tracker checks both observation LLI and external cycle slip flags (Doppler or GF cycle slips) so that wide-lane tracking arcs are reset immediately when a cycle slip occurs.

4. **Unit Tests & Invariants**:
   - Add unit tests verifying:
     - 0.5-cycle (half-cycle) slip detection on bands 1, 2, 5, 6, 7.
     - Unconditional base slip detection flags base slips in standard RTK mode.
     - `pair_epochs` resets to 0 upon cycle slip.
     - All functions <= 32 LOC, files < 500 LOC, 0 unwrap in prod.

## Verification Commands
Before delivering your handoff, execute and verify:
```bash
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p gneiss-rtk --lib
cargo test --test test_urban_canyon_e2e
python3 scripts/check_network_benchmark.py --smoke
python3 scripts/check_multignss_benchmark.py --smoke
```
Write your full report to `/Users/kevin/projects/gneiss/.agents/worker_m3_urban_canyon/handoff.md` following the Handoff Protocol (Observation, Logic Chain, Caveats, Conclusion, Verification Method).
Then notify parent with a concise message.

## 2026-09-25T06:19:15Z
You are Worker M3 for Milestone 3: Doppler-Assisted Cycle Slip Detection & Phase Continuity Validation.
Your working directory is: /Users/kevin/projects/gneiss/.agents/worker_m3_urban_canyon/
Read your task assignment in /Users/kevin/projects/gneiss/.agents/worker_m3_urban_canyon/DISPATCH.md, the authoritative user request in /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md, the project specifications in /Users/kevin/projects/gneiss/.agents/orchestrator_urban_canyon_gen2/PROJECT.md, and the deep survey findings in /Users/kevin/projects/gneiss/.agents/survey_explorer_3/survey_r3_r4.md and handoff.md.

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Strictly adhere to all AGENTS.md rules:
- File size < 500 LOC (watch out for mod.rs at 482, mw.rs at 485, formation.rs at 452, screening.rs at 425)
- Function size <= 32 LOC
- Nesting depth < 3 levels
- Zero unwrap() in production code
- Zero clippy/compiler warnings (cargo clippy --workspace --all-targets -- -D warnings)
- All workspace tests pass (cargo test -p gneiss-rtk --lib, cargo test --test test_urban_canyon_e2e)
- Both CI smoke guard scripts pass:
  python3 scripts/check_network_benchmark.py --smoke
  python3 scripts/check_multignss_benchmark.py --smoke

Implement F8 (multi-band Doppler slip detection across bands 1, 2, 5, 6, 7 with adaptive threshold 0.28-1.0 cyc in screening.rs) and F9 (unconditional base slip check in mod.rs and formation.rs, reset pair_epochs on slips, propagate slip flags to pw_tracker and mw.rs).

When all implementation and verifications pass, write /Users/kevin/projects/gneiss/.agents/worker_m3_urban_canyon/handoff.md and send a completion message to your parent.
