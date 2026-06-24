export const meta = {
  name: 'gneiss-coverage-drive',
  description: 'Fan out to write targeted tests for the lowest-coverage files in gneiss-rtk',
  phases: [
    { title: 'Write Tests', detail: '6 agents writing tests for different file groups' },
    { title: 'Measure', detail: 'Re-run tarpaulin and report improvement' },
  ],
};

phase('Write Tests');

const results = await parallel([
  // Agent 1: Math & utility functions
  function() {
    return agent(
      "You are writing tests for low-coverage math and utility files in the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss. Current coverage report shows these files need tests: crates/gneiss-rtk/src/math/covariance.rs (5/14 lines), crates/gneiss-rtk/src/math/inversion.rs (9/12 lines), crates/gneiss-rtk/src/engine/types.rs (0/14 lines). Read each file, understand the functions, and write targeted unit tests. Add tests inside existing #[cfg(test)] mod tests blocks, or create them if they don't exist. Each test should exercise the function with known inputs and verify outputs. Do NOT change production code. Run 'cargo test -p gneiss-rtk --lib 2>&1 | grep -E \"FAILED|test result\"' after adding tests.",
      { label: 'test-math', phase: 'Write Tests', model: 'sonnet' }
    );
  },

  // Agent 2: Processor pipeline
  function() {
    return agent(
      "You are writing tests for the processor pipeline in the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss. Current coverage: crates/gneiss-rtk/src/engine/processor/rtk.rs (0/290), rtk_iekf.rs (0/7), spp.rs (14/121), ins.rs (16/139), mod.rs (15/199). Read each file and write tests for the functions with lowest coverage. Focus on functions that can be tested with minimal setup — pure functions, configuration builders, error paths. For complex integration functions that need full engine state, write tests that at minimum exercise the error paths and boundary conditions. Add tests inside existing #[cfg(test)] blocks. Run 'cargo test -p gneiss-rtk --lib 2>&1 | grep -E \"FAILED|test result\"' after adding tests.",
      { label: 'test-processor', phase: 'Write Tests', model: 'sonnet' }
    );
  },

  // Agent 3: Standalone algorithms
  function() {
    return agent(
      "You are writing tests for standalone algorithm files in the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss. Current coverage: crates/gneiss-rtk/src/engine/tcar.rs (0/61), tight_iekf.rs (0/68), spp_tight.rs (3/224), ssr.rs (15/66). Read each file, identify testable functions, and write targeted unit tests. For tcar.rs: test the triple-carrier ambiguity resolution math with known frequencies. For tight_iekf.rs: test the tight coupling state machine transitions. For spp_tight.rs: test SPP computation with synthetic observations. For ssr.rs: test SSR correction parsing and application. Add tests inside existing #[cfg(test)] blocks or create them. Run 'cargo test -p gneiss-rtk --lib 2>&1 | grep -E \"FAILED|test result\"' after.",
      { label: 'test-algorithms', phase: 'Write Tests', model: 'sonnet' }
    );
  },

  // Agent 4: Factor graph estimators
  function() {
    return agent(
      "You are writing tests for factor graph estimator files in the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss. Current coverage: crates/gneiss-rtk/src/estimators/factor_graph/gnss_factors.rs (78/179), mod.rs (70/87), imu_factors.rs (93/99). Focus on gnss_factors.rs which has the largest gap (101 uncovered lines). Read the file and write tests for the GNSS factor functions with known synthetic inputs. Also check if there are integration test gaps in the factor graph optimizer. Run 'cargo test -p gneiss-rtk --lib 2>&1 | grep -E \"FAILED|test result\"' after adding tests.",
      { label: 'test-factor-graph', phase: 'Write Tests', model: 'sonnet' }
    );
  },

  // Agent 5: Smoother, updater, and calibration
  function() {
    return agent(
      "You are writing tests for smoother, updater, and calibration code in the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss. Current coverage: crates/gneiss-rtk/src/engine/smoother.rs (143/179), updater.rs (175/237), updater_math.rs (98/99). Also check calibration/ directory. Focus on edge cases: empty input, error paths, boundary conditions, numerical stability. For smoother.rs: test the epoch matching, phi submatrix building, and covariance regularization. For updater.rs: test the innovation computation with known state/measurement pairs. Run 'cargo test -p gneiss-rtk --lib 2>&1 | grep -E \"FAILED|test result\"' after.",
      { label: 'test-smoother-updater', phase: 'Write Tests', model: 'sonnet' }
    );
  },

  // Agent 6: PPP math and estimators
  function() {
    return agent(
      "You are writing tests for PPP math and estimator files in the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss. Current coverage: crates/gneiss-rtk/src/engine/ppp_math.rs (108/136), ppp_common.rs (needs check), estimators/spp.rs (219/305), estimators/araim.rs (55/59), estimators/ekf/filter.rs (361/428). Read each file and write tests for the uncovered functions. For ppp_math.rs: test OSB correction application, iono constraint building, SNr scaling. For filter.rs: test the remaining uncovered paths in state management. Run 'cargo test -p gneiss-rtk --lib 2>&1 | grep -E \"FAILED|test result\"' after adding tests.",
      { label: 'test-ppp-math', phase: 'Write Tests', model: 'sonnet' }
    );
  },
]);

// Phase 2: Re-measure
phase('Measure');

const coverageResult = await agent(
  "Run 'cargo tarpaulin -p gneiss-rtk --out Stdout --exclude-files scratch/* 2>&1 | tail -5' to get the new coverage percentage. Also run 'cargo test --workspace 2>&1 | grep -E \"test result|FAILED\"' to verify all tests pass. Report: (a) new coverage %, (b) whether all tests pass, (c) files still at 0% coverage.",
  { label: 'measure-coverage', phase: 'Measure', model: 'haiku' }
);

return {
  agent_results: '6 agents wrote tests',
  final_coverage: coverageResult || 'measurement pending',
};