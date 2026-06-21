export const meta = {
  name: 'gneiss-gap-discovery',
  description: 'Fan out to discover every remaining gap in the gneiss GNSS engine',
  phases: [
    { title: 'Discover', detail: '6 parallel explorers across test infra, test quality, math, features, code quality, smoother' },
    { title: 'Synthesize', detail: 'Merge, deduplicate, prioritize all findings' },
  ],
};

phase('Discover');

const TEST_INFRA_SCHEMA = {
  type: 'object',
  properties: {
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          file: { type: 'string' },
          line: { type: 'integer' },
          category: { type: 'string' },
          description: { type: 'string' },
          fix: { type: 'string' },
          priority: { type: 'string' }
        },
        required: ['file', 'line', 'category', 'description', 'fix', 'priority']
      }
    }
  },
  required: ['findings']
};

const TEST_QUALITY_SCHEMA = {
  type: 'object',
  properties: {
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          file: { type: 'string' },
          line: { type: 'integer' },
          category: { type: 'string' },
          description: { type: 'string' },
          fix: { type: 'string' },
          priority: { type: 'string' }
        },
        required: ['file', 'line', 'category', 'description', 'fix', 'priority']
      }
    }
  },
  required: ['findings']
};

const MATH_SCHEMA = {
  type: 'object',
  properties: {
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          file: { type: 'string' },
          line: { type: 'integer' },
          severity: { type: 'string' },
          description: { type: 'string' },
          suggested_fix: { type: 'string' }
        },
        required: ['file', 'line', 'severity', 'description', 'suggested_fix']
      }
    }
  },
  required: ['findings']
};

const FEATURES_SCHEMA = {
  type: 'object',
  properties: {
    gaps: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          category: { type: 'string' },
          feature: { type: 'string' },
          current_state: { type: 'string' },
          competitors_have: { type: 'string' },
          impact_on_accuracy: { type: 'string' },
          priority: { type: 'string' }
        },
        required: ['category', 'feature', 'current_state', 'impact_on_accuracy', 'priority']
      }
    }
  },
  required: ['gaps']
};

const CODE_QUALITY_SCHEMA = {
  type: 'object',
  properties: {
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          file: { type: 'string' },
          line: { type: 'integer' },
          category: { type: 'string' },
          description: { type: 'string' },
          fix: { type: 'string' },
          priority: { type: 'string' }
        },
        required: ['file', 'line', 'category', 'description', 'fix', 'priority']
      }
    }
  },
  required: ['findings']
};

const SMOOTHER_SCHEMA = {
  type: 'object',
  properties: {
    hypotheses: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          symptom: { type: 'string' },
          file: { type: 'string' },
          line_approx: { type: 'integer' },
          hypothesis: { type: 'string' },
          confidence: { type: 'string' },
          suggested_fix: { type: 'string' }
        },
        required: ['symptom', 'hypothesis', 'confidence', 'suggested_fix']
      }
    }
  },
  required: ['hypotheses']
};

const results = await parallel([
  function() {
    return agent(
      "You are auditing the gneiss GNSS navigation engine test infrastructure. Working directory: /Users/kevin/projects/gneiss. Find EVERY issue: missing #[test] attributes on functions inside #[cfg(test)] mod tests blocks, dead test files that are never compiled/run, unregistered integration tests (check tests/src/lib.rs), test duplication (verbatim copies), redundant proc macro attributes (duplicate #[test]), ignored tests with #[ignore]. For each finding provide exact file path + line number, what's wrong, the fix, and priority (P0=test never runs, P1=ignored without reason, P2=duplication). Be thorough - check every crate directory. The suspicious_tests_report.md at repo root has initial findings to verify and expand.",
      { label: 'explore-test-infra', phase: 'Discover', schema: TEST_INFRA_SCHEMA }
    );
  },
  function() {
    return agent(
      "You are auditing test assertion quality in the gneiss GNSS navigation engine. Working directory: /Users/kevin/projects/gneiss. Find EVERY issue: silent tests with zero assertions or that discard results with _, weak/trivial assertions like tautologies (assert!(x || !x)), loose tolerances where GNSS math tolerances are too loose to catch bugs (tropo delay at +-60%, wavelength checks spanning multiple constellations, coordinate transforms 1000x too loose), logic bugs in assertions (message contradicting check, wrong comparison operator, missing field checks like week but not TOW), and missing assertions on critical math paths. For each finding: exact file+line, what's wrong, fix, priority (P0=silent/tautology, P1=loose tolerance, P2=missing coverage). suspicious_tests_report.md has initial findings to verify and expand.",
      { label: 'explore-test-quality', phase: 'Discover', schema: TEST_QUALITY_SCHEMA }
    );
  },
  function() {
    return agent(
      "You are auditing mathematical correctness in the gneiss GNSS navigation engine. Working directory: /Users/kevin/projects/gneiss. The engine uses error-state EKF with left-multiplied attitude error R_true = (I - [ψ×]) R_est. The predictor.rs recently changed vel_att sign from +f_e_skew*dt to -f_e_skew*dt. Check: (1) Is measurement Jacobian sign in updater.rs/updater_math.rs consistent? (2) Are there other skew-symmetric or cross-product sign errors? (3) Verify quaternion multiplication order matches error-state convention. (4) Check IMU mechanization Coriolis and gravity signs. (5) Verify all phi blocks in compute_transition_matrix - especially vel_vel (I-2*omega_ie_skew*dt) and att_att (I-omega_ie_skew*dt). (6) Check smoother.rs backward recursion math, ISB state handling. For each finding: exact file+line, mathematical issue, suggested fix, severity (critical/major/minor).",
      { label: 'explore-math', phase: 'Discover', schema: MATH_SCHEMA }
    );
  },
  function() {
    return agent(
      "You are analyzing the gneiss GNSS engine for missing functionality vs industry leaders. Working directory: /Users/kevin/projects/gneiss. Read POST_MORTEM.md, SPRINT_PLAN.md, FUTURE_TECH.md, BENCHMARKS.md, ARCHITECTURE.md for context. Identify gaps in: (1) Correction models - VMF1/VMF3 grid troposphere, GIM/IONEX ionosphere, ocean tide loading, atmospheric pressure loading, pole tide, higher-order iono. (2) Constellation support - GLONASS IFB, BeiDou phase bias, QZSS/NavIC/SBAS. (3) Processing modes - PPP-RTK, Network RTK (MAC/FKP/VRS), multi-base RTK. (4) Algorithms - PAR status, TCAR, multi-epoch factor graph, robust estimation. (5) Data formats - BINEX, NVS, RTCM3 MSM coverage, SSR formats. (6) Performance - real-time capability, memory for long sessions, multi-session. Rate: Critical (blocks beating RTKLIB), Important (Leica/Novatel tier), Nice-to-have (Qinertia tier).",
      { label: 'explore-features', phase: 'Discover', schema: FEATURES_SCHEMA }
    );
  },
  function() {
    return agent(
      "You are auditing code quality in the gneiss GNSS engine. Working directory: /Users/kevin/projects/gneiss. Run 'cargo build --workspace 2>&1' and 'cargo test --workspace 2>&1' to capture compiler warnings. Find: (1) All compiler warnings with categorization (unused imports, unused variables, dead code, etc.). (2) Dead code paths - functions/structs defined but never called. (3) Unnecessary .clone() on large objects (matrices, vectors). (4) unwrap() calls that could panic in non-test code. (5) Code duplication across files, especially ppp_iekf vs ppp_ins_iekf, measurement vs measurement_math, processor/*.rs. (6) TODO/FIXME/HACK comments. (7) Unsafe blocks. (8) Orphaned .orig and .rej files to clean up. For each: file+line, description, fix, priority (P0=panic risk, P1=dead code/dup, P2=style).",
      { label: 'explore-code-quality', phase: 'Discover', schema: CODE_QUALITY_SCHEMA }
    );
  },
  function() {
    return agent(
      "You are investigating root causes of PPP/smoother accuracy bugs in gneiss. Working directory: /Users/kevin/projects/gneiss. Read smoother.rs, ppp_iekf.rs, predictor.rs, updater.rs, ppp_math.rs. Known symptoms: (1) Smoother degrades horizontal: fwd 5.29m Hz vs smooth 6.42m on Odaiba. (2) Shinjuku smoother catastrophe: 16.7m Hz vs 5.0m fwd. (3) GSDC PPP divergence: 178m on smartphone L5 data. Investigate: In smoother.rs - how does backward pass handle SPP prior? Are core_phi matrices correctly sized/indexed? How does C_k relate to forward covariance? Clock/position correlation mishandling? ISB state zeroing correct or workaround? In ppp_iekf.rs - how is SPP prior applied each epoch? Any state resets causing discontinuities? In predictor.rs - is measurement-side sign consistent with new vel_att sign? How is core_phi saved? For GSDC - smartphone vs survey-grade differences? Variance model issues? Numerical overflow paths? Provide specific code locations, hypotheses, confidence (high/medium/low), suggested fixes.",
      { label: 'explore-smoother', phase: 'Discover', schema: SMOOTHER_SCHEMA }
    );
  },
]);

phase('Synthesize');

const testInfra = results[0] ? (results[0].findings || []) : [];
const testQuality = results[1] ? (results[1].findings || []) : [];
const mathAudit = results[2] ? (results[2].findings || []) : [];
const missingFeatures = results[3] ? (results[3].gaps || []) : [];
const codeQuality = results[4] ? (results[4].findings || []) : [];
const smootherHypotheses = results[5] ? (results[5].hypotheses || []) : [];

// Count by priority
const allFindings = [].concat(
  testInfra.map(function(f) { return Object.assign({}, f, {source: 'test-infra'}); }),
  testQuality.map(function(f) { return Object.assign({}, f, {source: 'test-quality'}); }),
  mathAudit.map(function(f) { return Object.assign({}, f, {source: 'math-audit', priority: f.severity}); }),
  codeQuality.map(function(f) { return Object.assign({}, f, {source: 'code-quality'}); })
);

const p0 = allFindings.filter(function(f) { return f.priority === 'P0'; }).length;
const p1 = allFindings.filter(function(f) { return f.priority === 'P1'; }).length;
const p2 = allFindings.filter(function(f) { return f.priority === 'P2'; }).length;

log('Discovery complete:');
log('  Test Infra: ' + testInfra.length + ' findings');
log('  Test Quality: ' + testQuality.length + ' findings');
log('  Math Audit: ' + mathAudit.length + ' findings');
log('  Missing Features: ' + missingFeatures.length + ' gaps');
log('  Code Quality: ' + codeQuality.length + ' findings');
log('  Smoother/PPP: ' + smootherHypotheses.length + ' hypotheses');
log('  Priority: ' + p0 + ' P0, ' + p1 + ' P1, ' + p2 + ' P2');

return {
  test_infra: testInfra,
  test_quality: testQuality,
  math_audit: mathAudit,
  missing_features: missingFeatures,
  code_quality: codeQuality,
  smoother_hypotheses: smootherHypotheses,
  summary: {
    total_findings: allFindings.length + missingFeatures.length + smootherHypotheses.length,
    p0_count: p0,
    p1_count: p1,
    p2_count: p2,
  }
};