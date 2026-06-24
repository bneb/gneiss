export const meta = {
  name: 'gneiss-red-team',
  description: 'Adversarially verify every fix from the gap-analysis/fix-all sessions',
  phases: [
    { title: 'Review Commits', detail: 'Per-commit adversarial review — did each fix introduce regressions?' },
    { title: 'Sign Convention Audit', detail: 'Cross-check windup sign, vel_att sign, and measurement Jacobian consistency' },
    { title: 'Dead Code Safety', detail: 'Verify promoted/deleted code doesn\'t break anything' },
    { title: 'Skew Symmetric Audit', detail: 'Verify all call sites reference the correct canonical function' },
    { title: 'Test Quality Audit', detail: 'Verify strengthened tests test the right things with correct expected values' },
    { title: 'Missing Anything?', detail: 'Completeness critic — what did we miss?' },
    { title: 'Verdict', detail: 'Synthesize findings, render pass/fail verdict' },
  ],
};

phase('Review Commits');
phase('Sign Convention Audit');
phase('Dead Code Safety');
phase('Skew Symmetric Audit');
phase('Test Quality Audit');
phase('Missing Anything?');

const reviews = await parallel([
  // Reviewer 1: Per-commit adversarial review
  function() {
    return agent(
      "You are an adversarial code reviewer for the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss.\n\nRun: git log --oneline -5\n\nThen for each of the 4 recent commits (f3fcd20, 8108b58, 6ed8fb8, cf9376f), run: git show <hash> --stat && git show <hash>\n\nFor EACH commit, answer:\n1. Does the change do what it claims to do?\n2. Could it introduce a regression? (check: does the fix match the surrounding code conventions? does it change behavior in unexpected ways?)\n3. Are there any off-by-one errors, sign errors, or index errors?\n4. If it removes code — is there any caller that still references the removed code?\n5. If it adds an import — does the import path resolve correctly? Is there a circular dependency risk?\n\nBe skeptical. Assume every change is wrong until proven correct. Flag anything suspicious, even if you're not 100% sure.",
      { label: 'review-commits', phase: 'Review Commits', schema: {
        type: 'object',
        properties: {
          commit_reviews: { type: 'array', items: { type: 'object', properties: {
            commit: { type: 'string' },
            verdict: { type: 'string' },
            issues_found: { type: 'array', items: { type: 'object', properties: {
              severity: { type: 'string' },
              file: { type: 'string' },
              line: { type: 'integer' },
              description: { type: 'string' }
            }, required: ['severity', 'description'] } }
          }, required: ['commit', 'verdict'] } }
        },
        required: ['commit_reviews']
      }}
    );
  },

  // Reviewer 2: Sign convention cross-audit
  function() {
    return agent(
      "You are auditing sign conventions across the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss.\n\nThe engine uses an error-state formulation. There have been recent sign-related changes:\n1. Windup: ppp_ins_iekf.rs changed from +windup to -windup (matching ppp_iekf.rs)\n2. vel_att: predictor.rs uses -f_e_skew*dt (empirically stable)\n\nYour job: verify cross-file sign consistency.\n\nStep 1: Read the windup computation in gneiss-core (crates/gneiss-core/src/windup.rs). Understand what sign the windup value carries.\n\nStep 2: Verify ALL call sites use windup consistently:\n- grep -rn \"windup\" crates/gneiss-rtk/src/ --include=\"*.rs\" | grep -v test | grep -v \".orig\"\n- Check: ppp_iekf.rs, ppp_ins_iekf.rs, measurement.rs, ppp.rs\n- Is windup ADDED or SUBTRACTED at every site? Are they ALL consistent now?\n\nStep 3: Check the vel_att sign convention chain:\n- predictor.rs: vel_att = -f_e_skew * dt\n- updater.rs / updater_math.rs: How is attitude error projected into measurements? Do the measurement Jacobians use +skew or -skew?\n- smoother.rs: How is the attitude correction (delta_x[6:9]) computed? Does it use the same sign convention?\n- grep for 'skew_symmetric' usage in measurement Jacobians\n\nStep 4: Look for any OTHER sign-sensitive operations:\n- Coriolis term: 2*omega_ie.cross(&velocity) — is the sign correct?\n- Gravity: subtract or add?\n- Lever arm: body-to-ECEF rotation direction?\n\nReport every inconsistency found. If everything is consistent, say so explicitly for each category.",
      { label: 'audit-signs', phase: 'Sign Convention Audit', schema: {
        type: 'object',
        properties: {
          windup_consistency: { type: 'string' },
          vel_att_consistency: { type: 'string' },
          coriolis_consistency: { type: 'string' },
          issues: { type: 'array', items: { type: 'object', properties: {
            severity: { type: 'string' },
            file: { type: 'string' },
            line: { type: 'integer' },
            description: { type: 'string' }
          }, required: ['severity', 'description'] } }
        },
        required: ['windup_consistency', 'vel_att_consistency', 'issues']
      }}
    );
  },

  // Reviewer 3: Dead code safety
  function() {
    return agent(
      "You are verifying that dead code removal and promotion in the gneiss engine didn't break anything. Working directory: /Users/kevin/projects/gneiss.\n\nChanges made:\n1. DELETED: crates/gneiss-rtk/src/tests_updater.rs (root-level duplicate)\n2. DELETED: crates/gneiss-rtk/src/engine/tests_ppp_fg_mutants.rs (worthless stubs)\n3. PROMOTED: kinematics.rs — added 'pub mod kinematics;' to engine/mod.rs\n4. PROMOTED: test_sp3.rs — added 'pub mod test_sp3;' to engine/mod.rs\n5. REMOVED: 'mod tests_updater;' from crates/gneiss-rtk/src/lib.rs\n\nYour job:\n\nStep 1: Verify deletions didn't orphan callers:\n- grep -rn \"tests_updater\" crates/ --include=\"*.rs\" | grep -v target\n- grep -rn \"tests_ppp_fg_mutants\" crates/ --include=\"*.rs\" | grep -v target\n\nStep 2: For promoted modules (kinematics.rs, test_sp3.rs):\n- Read the promoted files. Do they compile? Are their dependencies all available in the engine module?\n- Check: does kinematics.rs call any functions that don't exist?\n- Check: does test_sp3.rs have imports that resolve correctly?\n- Run: cargo build -p gneiss-rtk 2>&1 to verify\n\nStep 3: Check for naming issues:\n- test_sp3.rs has 'test_' prefix but is NOT a test module (it's production code). Could this cause issues with cargo test?\n- kinematics.rs — is 'apply_kinematic_constraints' called anywhere? grep for it.\n\nStep 4: Check if we created any unused imports or dead code by promoting modules:\n- After promotion, does anything in the promoted modules generate warnings?\n\nReport any breakage, warnings, or naming concerns.",
      { label: 'audit-dead-code', phase: 'Dead Code Safety', schema: {
        type: 'object',
        properties: {
          deletions_safe: { type: 'string' },
          promotions_safe: { type: 'string' },
          naming_issues: { type: 'array', items: { type: 'string' } },
          issues: { type: 'array', items: { type: 'object', properties: {
            severity: { type: 'string' },
            description: { type: 'string' }
          }, required: ['severity', 'description'] } }
        },
        required: ['deletions_safe', 'promotions_safe', 'issues']
      }}
    );
  },

  // Reviewer 4: Skew symmetric dedup audit
  function() {
    return agent(
      "You are verifying the skew_symmetric deduplication in the gneiss engine. Working directory: /Users/kevin/projects/gneiss.\n\nThe canonical skew_symmetric is now at: crates/gneiss-rtk/src/engine/fgo/factors/mod.rs\n\nThree files now import it instead of defining their own:\n- crates/gneiss-rtk/src/engine/fgo/factors/imu.rs (uses super::skew_symmetric or super::super::skew_symmetric)\n- crates/gneiss-rtk/src/engine/predictor.rs (uses crate::engine::fgo::factors::skew_symmetric)\n- crates/gneiss-rtk/src/measurements/nhc.rs (uses crate::engine::fgo::factors::skew_symmetric)\n\nYour job:\n\nStep 1: Verify the canonical implementation is correct:\n- Read crates/gneiss-rtk/src/engine/fgo/factors/mod.rs — find the skew_symmetric function\n- Verify it constructs the standard skew-symmetric matrix: [[0, -z, y], [z, 0, -x], [-y, x, 0]]\n\nStep 2: Verify every call site compiles and calls the right function:\n- Check imu.rs: does it use super::skew_symmetric? Does it compile? (The imu.rs module is inside fgo/factors/, so super refers to fgo/factors/mod.rs)\n- Check predictor.rs: does the import path resolve? Is the function visible? (predictor.rs is inside engine/, fgo/factors is engine::fgo::factors)\n- Check nhc.rs: does the import path resolve? (nhc.rs is in measurements/, which is in gneiss-rtk)\n\nStep 3: Check test code:\n- grep for 'fn skew_symmetric' in test modules. These should remain as test-local helpers. Are they still there?\n- Check fgo/factors/pseudorange.rs and carrier_phase.rs test modules — do they still have their own skew_symmetric?\n\nStep 4: Verify no call site uses the old import paths:\n- grep -rn \"use.*skew_symmetric\" crates/ --include=\"*.rs\"\n- All imports should point to the canonical location OR be test-local definitions\n\nStep 5: Run 'cargo build -p gneiss-rtk 2>&1' and verify zero errors.\n\nReport any broken imports, missing visibility, or test compilation failures.",
      { label: 'audit-skew', phase: 'Skew Symmetric Audit', schema: {
        type: 'object',
        properties: {
          canonical_correct: { type: 'string' },
          all_imports_resolve: { type: 'string' },
          tests_intact: { type: 'string' },
          issues: { type: 'array', items: { type: 'object', properties: {
            severity: { type: 'string' },
            file: { type: 'string' },
            description: { type: 'string' }
          }, required: ['severity', 'description'] } }
        },
        required: ['canonical_correct', 'all_imports_resolve', 'issues']
      }}
    );
  },

  // Reviewer 5: Test quality — are the new assertions correct?
  function() {
    return agent(
      "You are auditing the strengthened test assertions in gneiss. Working directory: /Users/kevin/projects/gneiss.\n\nThe following tests were modified to have stronger assertions:\n1. tests/src/lib.rs: Changed matches! to assert_eq! for EngineMode::RtkIns\n2. test_compute_dd_doppler: Added finiteness, magnitude, and variance assertions\n3. test_loosely_coupled_jacobian: Added position/velocity change checks\n4. test_update_loosely_coupled_huber: Tightened bound to 0.1 < x < 5.0\n5. test_compute_update_iteration_s_threshold: Added dx.len() == 3 check\n6. test_compute_phase_windup: Changed from != checks to exact windup value assertions\n\nYour job: For EACH modified test, read the source code and answer:\n\n1. Are the new expected values CORRECT? Could they be wrong?\n   - For test #4 (Huber): with a 15m initial error and Huber k=3.0, what should the expected position correction be? Is 0.1 < x < 5.0 a valid bound or does it mask errors?\n   - For test #6 (windup): are the exact values (10.0 + w_sat) correct given the measurement model?\n\n2. Could any of the new assertions be flaky? (non-deterministic, depends on floating point rounding, etc.)\n\n3. Do the assertions actually test what they claim to test? Or are they testing something else entirely?\n\n4. For test #3 (loosely_coupled_jacobian): the test now checks state.position.vector != initial. But with zero initial position and a synthetic test, does the update actually change position? If lever_arm and omega_b produce zero correction, the assertion would fail.\n\n5. Run: cargo test -p gneiss-rtk --lib 2>&1 | grep -E \"test.*test_compute_dd_doppler|test.*loosely_coupled_jacobian|test.*loosely_coupled_huber|test.*iteration_s_threshold|test.*phase_windup\" — do all the modified tests pass?\n\nBe skeptical. Flag any assertion that could be wrong, flaky, or testing the wrong thing.",
      { label: 'audit-tests', phase: 'Test Quality Audit', schema: {
        type: 'object',
        properties: {
          test_results: { type: 'array', items: { type: 'object', properties: {
            test_name: { type: 'string' },
            assertions_correct: { type: 'string' },
            flaky_risk: { type: 'string' },
            issues: { type: 'array', items: { type: 'string' } }
          }, required: ['test_name', 'assertions_correct'] } }
        },
        required: ['test_results']
      }}
    );
  },

  // Reviewer 6: Completeness critic — what did we miss?
  function() {
    return agent(
      "You are a completeness critic for the gneiss GNSS engine fix session. Working directory: /Users/kevin/projects/gneiss.\n\nTwo Ultracode workflows ran: gap-discovery (6 explorers) and fix-all (4 fixers + verification). They found and fixed ~20 issues.\n\nYour job: Find what was MISSED.\n\nStep 1: Re-scan for remaining issues in categories we covered:\n- Are there any remaining .unwrap() calls that could panic? (grep -rn \".unwrap()\" crates/gneiss-rtk/src/ --include=\"*.rs\" | grep -v \"#\\[cfg(test)\\]\" | grep -v \"test\" | head -20)\n- Are there any remaining duplicate function definitions? (check for duplicated code patterns)\n- Are there any remaining test files or modules not compiled?\n- Are there any files with suspicious names (test_ prefix on production code)?\n\nStep 2: Check for issues we DIDN'T cover in our categories:\n- Error handling: are there any .unwrap() calls on Results in non-test code?\n- Integer overflow: any unchecked arithmetic?\n- Resource leaks: any open file handles not closed?\n- Race conditions: any unsafe shared state? (grep for unsafe, Mutex, RefCell)\n- Documentation: are there any doc comments that are now wrong after our changes?\n\nStep 3: Check git status:\n- Are there any .orig or .rej files still present?\n- Are there any untracked files that should be committed or .gitignored?\n\nStep 4: Check for any new warnings or issues:\n- Run: cargo build --workspace 2>&1 | grep -E \"warning\"\n- Run: cargo test --workspace 2>&1 | grep -E \"ignored\" \n- Are the remaining ignored tests justified?\n\nReport EVERYTHING you find, even minor issues. This is a completeness review — no issue is too small.",
      { label: 'completeness-critic', phase: 'Missing Anything?', schema: {
        type: 'object',
        properties: {
          missed_issues: { type: 'array', items: { type: 'object', properties: {
            severity: { type: 'string' },
            category: { type: 'string' },
            file: { type: 'string' },
            line: { type: 'integer' },
            description: { type: 'string' }
          }, required: ['severity', 'category', 'description'] } }
        },
        required: ['missed_issues']
      }}
    );
  },
]);

// ============================================================
// Phase: Verdict
// ============================================================
phase('Verdict');

const commitReview = reviews[0]?.commit_reviews || [];
const signAudit = reviews[1] || {};
const deadCodeAudit = reviews[2] || {};
const skewAudit = reviews[3] || {};
const testAudit = reviews[4]?.test_results || [];
const completeness = reviews[5]?.missed_issues || [];

// Gather all issues
const allIssues = [].concat(
  (commitReview || []).flatMap(function(r) { return (r.issues_found || []).map(function(i) { return Object.assign({}, i, {source: 'commit-review', commit: r.commit}); }); }),
  (signAudit.issues || []).map(function(i) { return Object.assign({}, i, {source: 'sign-audit'}); }),
  (deadCodeAudit.issues || []).map(function(i) { return Object.assign({}, i, {source: 'dead-code-audit'}); }),
  (skewAudit.issues || []).map(function(i) { return Object.assign({}, i, {source: 'skew-audit'}); }),
  (testAudit || []).flatMap(function(t) { return (t.issues || []).map(function(i) { return Object.assign({}, {description: i, source: 'test-audit', test: t.test_name, severity: 'warning'}); }); }),
  (completeness || []).map(function(i) { return Object.assign({}, i, {source: 'completeness'}); })
);

const criticals = allIssues.filter(function(i) { return i.severity === 'critical'; });
const majors = allIssues.filter(function(i) { return i.severity === 'major'; });
const warnings = allIssues.filter(function(i) { return i.severity === 'warning' || i.severity === 'minor'; });

log('Red Team Verdict:');
log('  Commit review: ' + (commitReview.length || 0) + ' commits reviewed');
log('  Sign audit: windup=' + (signAudit.windup_consistency || '?') + ', vel_att=' + (signAudit.vel_att_consistency || '?'));
log('  Dead code: deletions=' + (deadCodeAudit.deletions_safe || '?') + ', promotions=' + (deadCodeAudit.promotions_safe || '?'));
log('  Skew sym: canonical=' + (skewAudit.canonical_correct || '?') + ', imports=' + (skewAudit.all_imports_resolve || '?'));
log('  Tests: ' + (testAudit.length || 0) + ' tests audited');
log('  Missed: ' + (completeness.length || 0) + ' missed issues found');
log('');
log('  Issues: ' + criticals.length + ' critical, ' + majors.length + ' major, ' + warnings.length + ' minor');

var verdict = criticals.length > 0 ? 'FAIL' : majors.length > 0 ? 'CONDITIONAL PASS' : 'PASS';
log('');
log('VERDICT: ' + verdict);

return {
  verdict: verdict,
  criticals: criticals,
  majors: majors,
  warnings: warnings,
  sign_audit: { windup: signAudit.windup_consistency, vel_att: signAudit.vel_att_consistency },
  dead_code: { deletions: deadCodeAudit.deletions_safe, promotions: deadCodeAudit.promotions_safe },
  skew: { canonical: skewAudit.canonical_correct, imports: skewAudit.all_imports_resolve },
  completeness_missed: completeness.length,
};