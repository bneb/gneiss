## 2026-06-21T09:59:27Z
You are teamwork_preview_explorer.
Your working directory is /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_bug_15.
Your task is to analyze Bug 15: Incorrect Broadcast Clock TGD Correction.
- File: crates/gneiss-core/src/ephemeris.rs (around lines 403, etc.)
- Bug: tgd is subtracted for dual-frequency/iono-free clock corrections. Under these modes, TGD cancels, so it shouldn't be subtracted.

Investigate the code in crates/gneiss-core/src/ephemeris.rs and other related files (like crates/gneiss-rtk/src/engine/ppp.rs) to:
1. Locate where `tgd` is applied/subtracted in the satellite clock calculation (e.g., `calc_keplerian` or clock error correction functions).
2. Understand how the codebase handles/distinguishes between single-frequency and dual-frequency/ionosphere-free clock corrections.
3. Propose a clean, correct fix strategy that avoids subtracting TGD for dual-frequency/ionosphere-free modes, while retaining the correct subtraction for single-frequency modes if applicable.
4. Recommend a regression test design (named `test_broadcast_clock_tgd_correct`) that verifies the fix.

Write your findings to /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_bug_15/analysis.md and notify your parent.
