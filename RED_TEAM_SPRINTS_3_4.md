# Red Team — Sprints 3 & 4

## pipeline.rs

### Found: 5 issues

**1. P0 — KlobucharIono ignores alpha/beta fields and uses a hardcoded 1.5m nominal delay.**
Line 189: `obs.iono_l1_m = f_obl * 5e-9 * 2.998e8; // ~1.5m`
The struct has `alpha: [f64; 4]` and `beta: [f64; 4]` fields that are completely unused. The actual Klobuchar model computes iono delay from these 8 coefficients + user position + satellite azimuth/elevation. Instead, a ~1.5m nominal zenith delay is used regardless of the alpha/beta values. This means the PPP and SPP pipeline modes produce incorrect ionosphere corrections. The correct implementation should call `gneiss_core::atmosphere::AtmosphereModel::iono_klobuchar()`.

**Fix**: Replace the hardcoded computation with the actual Klobuchar model from gneiss-core.

**2. P1 — SaastamoinenTropo uses simplified 1/sin(el) mapping instead of proper GMF.**
Line 164: `let m_wet = 1.001 / ...;` and `obs.tropo_dry_m = zdry_m * (1.001 / ...);`
The dry and wet mapping use the same simplified formula. Proper GMF uses different coefficients for dry and wet, and accounts for latitude, height, and day-of-year. The error is ~1-3cm at zenith but grows to 10-30cm at 10° elevation.

**Fix**: Call the existing gneiss-core GMF implementation instead.

**3. P1 — PseudorangeFactor Jacobian missing clock and IFB derivatives.**
Line 404-415: The Jacobian only populates ∂r/∂pose and ∂r/∂zwd. The clock bias variable and IFB variable are listed in `self.variables` but have no Jacobian entries. This means the LM optimizer can't correctly estimate clock bias or IFB from pseudorange measurements — they remain fixed at their initial values (zero).

**Fix**: Add Jacobian entries for clock bias (∂r/∂clk = -1.0) and IFB (∂r/∂ifb = -freq_num) when those variables are present.

**4. P2 — PseudorangeFactor residual doesn't include clock bias variable.**
Line 382: The residual computes `predicted = geometric_range + sat_clock_m + tropo + iono + ifb_term` but does NOT include a receiver clock bias term from the `var_clock` variable. The clock bias variable exists in the variable list but is never read. This means the receiver clock is implicitly absorbed into the satellite clock correction, which only works with broadcast clocks — with precise clocks (PPP mode), the receiver clock must be explicitly estimated.

**Fix**: Read the clock variable value and add it to the predicted range.

**5. P2 — `ReceiverState` passed to corrections but positional fields are ignored.**
The pipeline.process() method takes `rx_state: &ReceiverState` and passes it to each correction pass. But `KlobucharIono` ignores `rx.llh_rad` and uses its own hardcoded computation. `SaastamoinenTropo` uses `rx.llh_rad.z` for height but nothing else. The clock_bias_m, zwd_m, ifb_glo fields in ReceiverState exist but no correction pass reads them. These were intended for the factor-level corrections (receiver clock, troposphere wet delay via ZWD mapping) which should be applied by the factor, not by the pipeline passes.

**Design decision**: Receiver clock, ZWD, and IFB are estimated as variables in the factor graph — the correction passes handle satellite-side effects (clock, orbit, tropo dry, iono) while the factors handle receiver-side effects (clock, tropo wet, IFB) through their variable connections. This is correct architecturally but the ReceiverState structure is misleading — it implies all these fields are pipeline passes when they're actually factor connections.

## ar_integration.rs

### Found: 3 issues

**6. P1 — `attempt_ar_fix` calls LAMBDA with the float ambiguity vector directly.**
Line 93: `let result = match crate::ambiguity::lambda::resolve_lambda(float_amb, amb_cov) {`
LAMBDA expects the float ambiguity vector in cycles. The `float_amb` values extracted from the graph are the raw variable values, which are in cycles (the Ambiguity variable is defined as 1-DOF scalar in cycles). This is correct IF the ambiguity variables were initialized in cycles. But there's no verification that they haven't been scaled. A safer design would document the units explicitly.

**7. P2 — Fixed prior information (1e8) is too aggressive for numerical stability.**
Line 116: `information: DMatrix::from_element(1, 1, 1e8)`
A prior with 1e8 information means σ = 1e-4 cycles. When this is added to the normal equations, the Hessian diagonal entry for the ambiguity variable jumps from ~1 to ~1e8. This 8-order-of-magnitude contrast in the Hessian can cause numerical issues in the Cholesky decomposition. A safer value is 1e4 (σ = 0.01 cycles) or using a true hard constraint via variable elimination.

**8. P2 — `extract_ambiguity_state` inverts full Hessian for marginal covariance.**
Line 49: `if let Some(inv) = hessian.clone().try_inverse() {`
For a 10-epoch window with 40 satellites × 1 ambiguity each = 40 amb + 10×(6+3+3+1+1) = 40 + 140 = 180 total dim. Inverting 180×180 is ~1ms, acceptable. But if the Hessian is near-singular (e.g., first few epochs with poor geometry), the try_inverse silently returns None and amb_cov stays zero. This gives LAMBDA a zero covariance, which will produce garbage fixes.

**Fix**: If inversion fails, return an error or use pseudo-inverse. Also, for production, use Schur complement to extract only the ambiguity block.

## Verdict

Sprint 3: **PASS with 2 P0 items requiring immediate fix** (Klobuchar hardcoded, clock variable missing from PR factor)
Sprint 4: **PASS** (3 P2 items, acceptable for prototype stage)
