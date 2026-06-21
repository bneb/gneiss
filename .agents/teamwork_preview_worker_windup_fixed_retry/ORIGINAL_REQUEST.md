## 2026-06-21T15:00:08Z

<USER_REQUEST>
Resume work at /Users/kevin/projects/gneiss/.agents/teamwork_preview_worker_windup_fixed_retry.
Your role is to implement the phase wind-up correction sign fix (Bug 18) and verify it.

Task Details:
1. In `crates/gneiss-rtk/src/engine/measurement.rs`:
   - In `apply_windup_to_obs` (around line 30-37), change `*cp += windup;` and `*cp2 += windup;` to subtraction (`-=`).
   - Add the following unit test to verify the sign correction (e.g. under `mod tests` block):
     ```rust
     #[test]
     fn test_phase_windup_correction_sign_rtk() {
         use gneiss_core::sat::{Constellation, SatelliteId};
         
         let mut obs = DdObservation {
             sat: SatelliteId {
                 constellation: Constellation::Gps,
                 prn: 1,
             },
             pr_l1: 0.0,
             pr_l2: None,
             cp_l1: Some(10.0),
             cp_l2: Some(20.0),
             doppler: 0.0,
             snr: 45.0,
             locktime: None,
         };

         let windup = 0.25; // 0.25 cycles of positive wind-up
         apply_windup_to_obs(&mut obs, windup);

         // Corrected carrier phase = raw_cp - windup
         assert_eq!(obs.cp_l1.unwrap(), 9.75);
         assert_eq!(obs.cp_l2.unwrap(), 19.75);
     }
     ```

2. In `crates/gneiss-rtk/src/engine/ppp.rs`:
   - In `update_phase_ambiguities` (around lines 495-504, 545-546), change all instances where `wup` is added to `cp1` or `cp2` (e.g., `cp1 + wup`, `sat.cp2.unwrap() + wup`) to subtraction (`- wup`).
   - In `add_uduc_ambiguities` (around lines 588-589), change additions of `wup` to subtraction (`- wup`).

3. In `crates/gneiss-rtk/src/engine/ppp_iekf.rs`:
   - In `predict_carrier_phase` or where carrier phase predictions are made, change `cp1 + windup` to `cp1 - windup`.
   - In UDUC phase measurement residual construction (around lines 1026 and 1043), change `sat.cp1.unwrap() + windup` and `sat.cp2.unwrap() + windup` to subtract `windup` instead.

4. Verify:
   - Run tests (`cargo test -p gneiss-rtk` and `cargo test --workspace`) and verify they pass cleanly.
   - Run formatting check: `cargo fmt --check` (run `cargo fmt` if there are any errors).
   - Ensure the regression test passes with the fix applied but fails if the fix is reverted.

5. Outputs:
   - Document your changes in `changes.md` in your directory.
   - Write a detailed handoff report in `handoff.md` in your directory, detailing the files modified, exact changes made, commands executed, and output of the test results.

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A Forensic Auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Send a completion message to the parent (conversation ID: 947de8ff-f313-48f8-be52-d7ba9185b0cc) when done.
</USER_REQUEST>
