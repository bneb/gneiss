# RTK Verification and Testing Skill

This skill governs the development life cycle of GNSS/INS RTK systems. It ensures that critical EKF math (such as analytical Jacobians) and file parsers (such as RINEX/UBX/RTCM3) are thoroughly tested, preventing silent regressions and EKF divergence.

## 1. Numerical Jacobian Verification
Analytical Jacobians in EKFs (e.g., measurement models, motion constraints, and process transition matrices) are highly prone to:
* Sign flips (especially w.r.t. attitude error or coordinate transformations).
* Incorrect scaling (e.g., using degrees instead of radians, or neglecting frequency factors).
* Mismatched state indices.

### The Rule of Numerical Verification
Whenever you modify or add a Jacobian in the filter or constraints, **you MUST write a unit test that verifies the analytical Jacobian against a numerical Jacobian computed using central finite differences.**

## 2. Parser Edge-Case Testing
Every GNSS parser (RINEX, UBX, RTCM) must have dedicated test coverage demonstrating:
1. **Flag Parsing**: Verification of cycle-slip indicators (LLI bits), lock-time counters, and quality flags.
2. **Malformed Input Safety**: Graceful error handling or exclusion of corrupt data lines instead of panic/crash.
3. **Empty Field Handling**: Correct identification of empty or missing observation values (e.g., whitespace spacing in RINEX, or missing fields in UBX).

## 3. Red-to-Green Bug Isolation
Before applying a fix for EKF divergence or parser failures:
1. **Reproduce the bug in a unit test first** (e.g., a test that asserts the incorrect Jacobian sign or the ignored field).
2. Verify the test fails (Red).
3. Apply the fix.
4. Verify the test passes (Green).
