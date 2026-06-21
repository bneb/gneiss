export const meta = {
  name: 'gneiss-fix-all',
  description: 'Fan out to fix all remaining P0/P1/P2 issues discovered in the gap analysis',
  phases: [
    { title: 'Fix Tests', detail: 'Strengthen weak test assertions' },
    { title: 'Fix Dead Code', detail: 'Register or remove dead modules and orphaned files' },
    { title: 'Fix Code Quality', detail: 'Fix cfg(test), dedup code, gate test-only structs' },
    { title: 'Fix Comments', detail: 'Improve smoother/math documentation' },
    { title: 'Verify', detail: 'Run tests and confirm clean build' },
  ],
};

// ============================================================
// Phase 1: Parallel fixes
// ============================================================
phase('Fix Tests');
phase('Fix Dead Code');
phase('Fix Code Quality');
phase('Fix Comments');

const fixResults = await parallel([
  // Agent A: Strengthen weak test assertions
  function() {
    return agent(
      "You are strengthening weak test assertions in the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss. Fix the following tests to have meaningful assertions. Make minimal, targeted edits — don't rewrite whole tests unless necessary.\n\n1. tests/src/lib.rs line 56: Change assert!(matches!(engine.config.mode, EngineMode::RtkIns | EngineMode::SppIns | EngineMode::PppIns)) to assert_eq!(engine.config.mode, gneiss_rtk::engine::EngineMode::RtkIns) — the config was explicitly set to RtkIns.\n\n2. crates/gneiss-rtk/src/engine/measurement.rs: test_compute_dd_doppler only has assert!(!u.z.is_nan()). Add assertions that the doppler innovation is finite and within a reasonable range (abs < 1000 Hz for these inputs), and that the variance u.r is positive.\n\n3. crates/gneiss-rtk/src/engine/tests_updater.rs line 371: test_loosely_coupled_jacobian only checks res.is_ok(). The test name says 'Jacobian' but nothing is verified. Add assertions checking the innovation vector has 6 elements (3 position + 3 velocity) and at minimum that the state was modified (position/velocity changed from initial).\n\n4. crates/gneiss-rtk/src/engine/tests_updater.rs line 406: test_update_loosely_coupled_huber has absurd bound 'state.position.vector.x > 0.0 && state.position.vector.x < 15.0' after a 15m initial error. The Huber scaling should bring the correction down significantly. Replace with a tighter bound like 'state.position.vector.x > 0.1 && state.position.vector.x < 5.0'.\n\n5. crates/gneiss-rtk/src/engine/tests_updater.rs line 632: test_compute_update_iteration_s_threshold only checks res.is_ok(). When the function succeeds, verify the result contains valid innovation/covariance values with assert!(result.is_ok() && result.unwrap().dx.len() == 3).\n\n6. crates/gneiss-rtk/src/engine/measurement.rs: test_compute_phase_windup uses assert!(x != 10.0) which is weak. Replace with an assertion that the windup correction amount w_sat is non-zero and that the corrected phase differs from the original by exactly w_sat.\n\nFor each fix, read the surrounding code first to understand the test, then make the minimal edit. Run 'cargo test -p gneiss-rtk --lib' after ALL edits to verify.",
      { label: 'fix-tests', phase: 'Fix Tests', model: 'sonnet' }
    );
  },

  // Agent B: Fix dead/orphaned code
  function() {
    return agent(
      "You are fixing dead and orphaned code in the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss. Make these fixes:\n\n1. kinematics.rs: Read crates/gneiss-rtk/src/engine/kinematics.rs. If it contains useful production code, add 'pub mod kinematics;' to crates/gneiss-rtk/src/engine/mod.rs. If it's just stubs/scratch, delete the file. Make the call based on content.\n\n2. Root-level tests_updater.rs: crates/gneiss-rtk/src/tests_updater.rs is an orphaned 90-line duplicate of engine/tests_updater.rs (26k). Check if it has any unique tests not in engine/tests_updater.rs. If unique, move them to engine/tests_updater.rs. Then remove the 'mod tests_updater;' declaration from crates/gneiss-rtk/src/lib.rs and delete the file.\n\n3. test_sp3.rs: Read crates/gneiss-rtk/src/engine/test_sp3.rs. It has production code (Sp3Record, Sp3Epoch). Check if this is a test file or production code. If it's production code that should be in the parsers crate, note that. If it's a test, add '#[cfg(test)] mod test_sp3;' to engine/mod.rs. Fix any compilation issues.\n\n4. tests_ppp_fg_mutants.rs: Read crates/gneiss-rtk/src/engine/tests_ppp_fg_mutants.rs. It has a test module with broken import paths. Either fix the imports and declare it in engine/mod.rs, or delete it if the tests are worthless.\n\n5. ArMock struct: In crates/gneiss-rtk/src/engine/ppp_iekf.rs, find the ArMock struct definition (around line 14). If it's only used in tests, wrap it with #[cfg(test)].\n\nAfter each fix, verify with 'cargo build -p gneiss-rtk 2>&1 | grep -E \"error|warning\"' and fix any issues.",
      { label: 'fix-dead-code', phase: 'Fix Dead Code', model: 'sonnet' }
    );
  },

  // Agent C: Fix code quality issues
  function() {
    return agent(
      "You are fixing code quality issues in the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss. Make these fixes:\n\n1. Missing #[cfg(test)] on test modules:\n   - crates/gneiss-core/src/lib.rs line 10: Change 'mod geodetic_tests;' to '#[cfg(test)] mod geodetic_tests;'\n   - crates/gneiss-rtk/src/lib.rs: Find 'mod tests_ekf;', 'mod tests_predictor;', 'mod tests_updater;' and prefix each with '#[cfg(test)]'\n\n2. Duplicate snr_scale tests: test_snr_scale() is duplicated in three places:\n   - crates/gneiss-rtk/src/engine/ppp_common.rs (the definition site — KEEP this one)\n   - crates/gneiss-rtk/src/engine/ppp_iekf.rs (find and REMOVE the duplicate test_snr_scale)\n   - crates/gneiss-rtk/src/engine/ppp_ins_iekf.rs (find and REMOVE the duplicate test_snr_scale)\n   Only remove the test functions, leave the production code.\n\n3. skew_symmetric deduplication: The function is duplicated in 4 locations:\n   - crates/gneiss-rtk/src/engine/fgo/factors/mod.rs\n   - crates/gneiss-rtk/src/engine/fgo/factors/imu.rs\n   - crates/gneiss-rtk/src/engine/predictor.rs (private fn)\n   - crates/gneiss-rtk/src/measurements/nhc.rs\n   The cleanest approach: make the one in fgo/factors/mod.rs pub, and have the others call it. If any is private and only used locally, add a comment noting the duplication is intentional for module isolation.\n\nAfter each fix, verify with 'cargo build -p gneiss-rtk 2>&1 | grep -E \"error\"' and fix any issues. Then run 'cargo test -p gneiss-rtk --lib 2>&1 | grep -E \"test result|FAILED\"'.",
      { label: 'fix-code-quality', phase: 'Fix Code Quality', model: 'sonnet' }
    );
  },

  // Agent D: Fix comments/documentation issues
  function() {
    return agent(
      "You are improving documentation comments in the gneiss GNSS smoother and predictor. Working directory: /Users/kevin/projects/gneiss. Make these fixes:\n\n1. crates/gneiss-rtk/src/engine/smoother.rs line 117 area: The code zeros clock bias (15) and ISBs (16-18) as white-noise states but does NOT zero clock drift (19). Add a comment explaining WHY clock drift is safe to smooth (it has temporal correlation via Allan variance, unlike clock bias which jumps by km-equivalents between epochs). The comment at lines 117-126 already explains the white-noise rationale — extend it to mention clock drift and ZWD explicitly.\n\n2. crates/gneiss-rtk/src/engine/smoother.rs line 255 area: invert_p_pred excludes states 15-18 from inversion. Add a comment clarifying why clock drift (19) and ZWD (20) are INCLUDED in the active set for inversion (they have physically meaningful temporal correlation).\n\n3. crates/gneiss-rtk/src/engine/predictor.rs lines 86-90: The vel_att sign convention comment is misleading. The error state stored at indices 6-8 is actually d_theta (the correction) such that R_true = (I + [d_theta×]) R_est. This d_theta = -ψ where ψ is the conventional attitude error. Under this convention, the velocity-attitude coupling IS -[f_e×]*dt. Update the comment to explain this clearly instead of saying 'empirically, the negative sign produces stable convergence.'\n\n4. crates/gneiss-rtk/src/engine/predictor.rs lines 255-273: full_x_predict omits attitude indices 6-8 (they remain 0). Add a brief comment at line 258-259 explaining this is correct because the error-state attitude is expected to be zero after injection, and the smoother handles attitude separately via predicted_attitude quaternion.\n\nMake targeted comment edits only. Do not change any logic. Verify with 'cargo build -p gneiss-rtk 2>&1 | grep -E \"error\"'.",
      { label: 'fix-comments', phase: 'Fix Comments', model: 'sonnet' }
    );
  },
]);

// ============================================================
// Phase 2: Verify
// ============================================================
phase('Verify');

// Run tests
const testResult = await agent(
  "Run 'cargo test --workspace 2>&1' and report the test results. Also run 'cargo build --workspace 2>&1 | grep -E \"warning|error\"' to check for any remaining warnings or errors. Report the exact output of both commands.",
  { label: 'verify-tests', phase: 'Verify', model: 'haiku' }
);

const summary = {
  test_fixes: fixResults[0] || 'Agent A result',
  dead_code_fixes: fixResults[1] || 'Agent B result',
  code_quality_fixes: fixResults[2] || 'Agent C result',
  comment_fixes: fixResults[3] || 'Agent D result',
  test_output: testResult || 'Test result',
};

return summary;