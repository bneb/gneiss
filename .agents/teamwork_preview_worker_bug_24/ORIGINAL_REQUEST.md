## 2026-06-21T15:03:42Z
You are teamwork_preview_worker.
Your working directory is /Users/kevin/projects/gneiss/.agents/teamwork_preview_worker_bug_24.
Your task is to implement the fix for Bug 24: Outlier Tolerance in Precise Clock Gaps.

Please implement the following changes:
1. In `crates/gneiss-parsers/src/rinex_clk.rs` in the `get_clock_bias` method:
   - Change the binary search block to return the exact record bias immediately if an exact match is found:
     ```rust
     let idx = match records.binary_search_by(|r| r.time.partial_cmp(&t).unwrap()) {
         Ok(i) => return Some(records[i].bias),
         Err(i) => i,
     };
     ```
   - Update the gap checking block after `r1` and `r2` are defined. If the gap is too large to interpolate (e.g. either endpoint distance is greater than 900.0 seconds), return `None`.
     ```rust
     let dt = r2.time - r1.time;
     if dt == 0.0 || (t - r1.time).abs() > 900.0 || (r2.time - t).abs() > 900.0 {
         return None;
     }
     ```
   - Add a regression unit test named `test_precise_clock_gap_tolerance` inside the `tests` module in `rinex_clk.rs`. The test must check that:
     * An exact match at a timestamp (which is preceded by a gap > 900s) returns the correct bias (and doesn't incorrectly trigger the gap check).
     * Non-exact matches within gaps > 900s return `None`.
     * Non-exact matches within gaps <= 900s correctly interpolate.
     * Extrapolations beyond bounds (> 900s) return `None`.
2. Verify that the codebase builds cleanly (`cargo build`) and tests pass cleanly (`cargo test`).
3. Write your implementation report to /Users/kevin/projects/gneiss/.agents/teamwork_preview_worker_bug_24/handoff.md and notify your parent.

MANDATORY WARNING: DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A Forensic Auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.
