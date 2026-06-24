export const meta = {
  name: 'gneiss-standards-fix',
  description: 'Fan out to measure baseline and fix the worst standards violations in parallel',
  phases: [
    { title: 'Measure + Fix All', detail: 'Coverage/mutant baseline, flatten nesting, split big files, refactor long functions' },
    { title: 'Verify', detail: 'Run tests and confirm standards improvements' },
  ],
};

phase('Measure + Fix All');

const results = await parallel([
  // Agent 1: Measure coverage and mutation baseline
  function() {
    return agent(
      "Run code coverage and mutation testing baselines for the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss. Step 1: Run 'cargo tarpaulin --out Html --output-dir target/tarpaulin -p gneiss-rtk 2>&1 | tail -10' and report the coverage percentage. Step 2: Run 'cargo mutants -p gneiss-rtk -j 4 --timeout 15 2>&1 | tail -15' to get a mutation score baseline. If it takes too long, try with fewer files. Report: (a) line coverage %, (b) mutation score (killed/total), (c) any files with 0% coverage or high mutant survival.",
      { label: 'measure-baseline', phase: 'Measure + Fix All' }
    );
  },

  // Agent 2: Flatten 8-level nesting in AR cascade logic
  function() {
    return agent(
      "You are fixing deep nesting in the gneiss GNSS engine to meet the CLAUDE.md standard of <3 levels. Working directory: /Users/kevin/projects/gneiss. The worst offender is crates/gneiss-rtk/src/engine/ppp_iekf.rs around the AR cascade logic — it has up to 8 levels of nesting. Read the resolve_cascade_ar function and related code. Fix strategy: use guard clauses with early returns/continues, extract nested blocks into small private helper functions, use let-else patterns. Goal: bring all nesting from 8 to max 3-4 levels. Make minimal, safe edits. Do NOT change any logic, math, or control flow semantics — only restructure with early returns, guard clauses, and extracted helpers. After editing, run 'cargo test -p gneiss-rtk --lib 2>&1 | grep -E \"test result|FAILED\"' to verify. Run 'cargo build -p gneiss-rtk 2>&1 | grep error' to check compilation.",
      { label: 'fix-nesting', phase: 'Measure + Fix All', model: 'sonnet' }
    );
  },

  // Agent 3: Split measurement.rs (2199 LOC) into submodules
  function() {
    return agent(
      "You are splitting oversized files in the gneiss GNSS engine to meet the CLAUDE.md standard of <500 LOC per file. Working directory: /Users/kevin/projects/gneiss. Target: crates/gneiss-rtk/src/engine/measurement.rs (2199 LOC). This file contains multiple logical groups. Strategy: create an engine/measurement/ directory with mod.rs re-exporting everything, and sub-files for each group. Each sub-file must be under 500 LOC. Move tests with their corresponding code. Update engine/mod.rs to declare the new module structure. IMPORTANT: This is a large refactor. Be methodical — do one group at a time, verify compilation after each. Run 'cargo build -p gneiss-rtk 2>&1 | grep error' and 'cargo test -p gneiss-rtk --lib 2>&1 | grep \"FAILED\"' after the full refactor to verify. Do NOT change any logic — only reorganize into submodules.",
      { label: 'split-measurement', phase: 'Measure + Fix All', model: 'sonnet' }
    );
  },

  // Agent 4: Refactor functions over 32 LOC
  function() {
    return agent(
      "You are refactoring the longest functions in the gneiss GNSS engine to meet the CLAUDE.md standard of <32 LOC per function. Working directory: /Users/kevin/projects/gneiss. Top targets: (1) crates/gneiss-rtk/src/estimators/ekf/filter.rs: RtkState::new at ~77 LOC — break the struct literal into logical initialization groups (extract init_position, init_velocity, init_covariance, init_ambiguities helpers). (2) compute_dd_pseudorange at ~42 LOC — extract helper sub-functions. (3) update_windup_state_and_obs at ~40 LOC — split into update_windup and apply_windup. Make minimal, safe edits. Do NOT change any logic or math. After each refactor, run 'cargo test -p gneiss-rtk --lib 2>&1 | grep -E \"FAILED|test result\"' to verify.",
      { label: 'refactor-functions', phase: 'Measure + Fix All', model: 'sonnet' }
    );
  },
]);

// Phase 2: Verify
phase('Verify');

const verifyResult = await agent(
  "Run 'cargo test --workspace 2>&1 | grep -E \"test result|FAILED\"' and 'cargo build --workspace 2>&1 | grep -E \"warning|error\"' to verify everything passes. Report any failures or warnings. Also check: how many files are still over 500 LOC? Run 'find crates -name \"*.rs\" -not -path \"*/target/*\" -exec wc -l {} \\; | sort -rn | awk \"$1 > 500\" | wc -l' to count.",
  { label: 'verify-all', phase: 'Verify', model: 'haiku' }
);

return {
  baseline: results[0] || 'not run',
  nesting_fix: results[1] || 'done',
  file_split: results[2] || 'done',
  function_refactor: results[3] || 'done',
  verify: verifyResult || 'done',
};