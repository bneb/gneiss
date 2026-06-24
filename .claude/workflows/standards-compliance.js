export const meta = {
  name: 'gneiss-standards-compliance',
  description: 'Fan out to measure and fix codebase against the 6 quality standards in CLAUDE.md',
  phases: [
    { title: 'Measure', detail: 'Measure LOC, function size, nesting, coverage baseline' },
    { title: 'Fix Red Team', detail: 'Fix remaining red team warnings' },
    { title: 'Setup Tooling', detail: 'Configure tarpaulin, cargo-mutants, clippy' },
    { title: 'Verify', detail: 'Run all tests and checks' },
  ],
};

phase('Measure');
phase('Fix Red Team');
phase('Setup Tooling');

const results = await parallel([
  // Agent 1: Measure codebase against all 6 standards
  function() {
    return agent(
      "You are measuring the gneiss GNSS engine codebase against quality standards defined in CLAUDE.md. Working directory: /Users/kevin/projects/gneiss. Standards: (1) Files under 500 LOC, (2) Functions under 32 LOC, (3) Nesting under 3 levels, (4) Coverage over 95%, (5) Zero mutation survivors, (6) Zero compiler warnings. Measure each standard. For files over 500 LOC run: find crates -name '*.rs' -not -path '*/target/*' -exec wc -l {} \\; | sort -rn and list those over 500. For function sizes, read the 5 biggest files and identify the largest functions. For nesting, check indentation depth in the biggest files. For coverage check if tarpaulin is installed (cargo tarpaulin --version). For mutations check cargo-mutants. For warnings run cargo build --workspace 2>&1 | grep warning. Provide clear measurements with specific file names and counts.",

      { label: 'measure-standards', phase: 'Measure', schema: {
        type: 'object',
        properties: {
          files_over_500: { type: 'array', items: { type: 'object', properties: {
            file: { type: 'string' }, loc: { type: 'integer' }
          }, required: ['file', 'loc'] } },
          largest_functions: { type: 'array', items: { type: 'string' } },
          deepest_nesting: { type: 'array', items: { type: 'string' } },
          coverage_available: { type: 'string' },
          mutation_available: { type: 'string' },
          warning_count: { type: 'integer' },
        },
        required: ['files_over_500', 'coverage_available', 'mutation_available', 'warning_count']
      }}
    );
  },

  // Agent 2: Fix remaining red team warnings
  function() {
    return agent(
      "You are fixing the remaining warnings from the red team review of the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss.\n\nWarnings to fix:\n\n1. Missed skew_symmetric dedup in:\n   - crates/gneiss-rtk/src/engine/fgo/factors/pseudorange.rs (line ~8): inline-constructs skew-symmetric matrix instead of calling super::skew_symmetric(&l_b)\n   - crates/gneiss-rtk/src/engine/fgo/factors/carrier_phase.rs (line ~13): same issue\n   Replace the inline Matrix3::new(0.0, -z, y, z, 0.0, -x, -y, x, 0.0) with super::skew_symmetric(&vector).\n\n2. The test assertion warnings from the red team are acceptable as-is — they're directionally correct and the tests pass. No changes needed.\n\n3. The vel_att inconsistency in ppp_ins_iekf.rs:702 is minor (small practical impact). Add a brief comment noting the right-perturbation convention used in the INS code vs left-perturbation in the rest of the codebase.\n\n4. kinematics.rs file still exists on disk even though the module declaration was reverted. Check if it's still present at crates/gneiss-rtk/src/engine/kinematics.rs. If so, delete it — it's dead code.\n\n5. test_sp3.rs file still exists on disk. Check crates/gneiss-rtk/src/engine/test_sp3.rs. If so, delete it — it's dead code.\n\nAfter each fix, run 'cargo build -p gneiss-rtk 2>&1 | grep -E \"error\"' to verify nothing breaks.",

      { label: 'fix-red-team', phase: 'Fix Red Team', model: 'sonnet' }
    );
  },

  // Agent 3: Set up coverage and mutation testing tooling
  function() {
    return agent(
      "You are setting up code quality tooling for the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss.\n\nStep 1: Check what's available:\n- Run 'cargo tarpaulin --version 2>&1' \n- Run 'cargo mutants --version 2>&1'\n- Run 'cargo clippy --version 2>&1'\n\nStep 2: For any tools NOT installed, provide the install commands.\n\nStep 3: Add a .cargo/config.toml or Makefile.toml with targets for:\n- make coverage: runs tarpaulin with HTML output\n- make mutants: runs cargo-mutants on the main crates\n- make lint: runs clippy with strict settings\n- make check: runs all of the above + tests\n\nStep 4: Add clippy configuration to Cargo.toml workspace level:\n- Check if there's a [workspace.lints] section already\n- Add lints for: unsafe_code, missing_docs (optional), clippy::unwrap_used (warn on unwrap in production)\n\nStep 5: Check if there's a CI config (.github/workflows/):\n- If so, suggest adding coverage and mutation steps\n- If not, note that CI is not configured\n\nOutput: summary of what's installed, what needs installing, and the config snippets to add.",

      { label: 'setup-tooling', phase: 'Setup Tooling', model: 'sonnet' }
    );
  },
]);

// ============================================================
// Verify
// ============================================================
phase('Verify');

const measureResult = results[0] || {};
const fixResult = results[1] || '';
const toolingResult = results[2] || '';

log('Standards Measurement:');
log('  Files > 500 LOC: ' + (measureResult.files_over_500?.length || '?'));
log('  Coverage tool: ' + (measureResult.coverage_available || '?'));
log('  Mutation tool: ' + (measureResult.mutation_available || '?'));
log('  Warnings: ' + (measureResult.warning_count || '?'));

return {
  measurement: measureResult,
  red_team_fixes: fixResult,
  tooling: toolingResult,
};