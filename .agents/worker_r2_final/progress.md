# Progress — 2026-09-13T01:13:30Z
Last visited: 2026-09-13T01:13:30Z

## Current Status
- Verified WTZR 5-hour static PPP benchmark achieves sub-meter accuracy: `p50 = 0.860m`, `RMS = 0.840m`, `3D RMS = 1.089m`.
- Investigated F9P kinematic vehicle drive error decomposition:
  - Error vs CSRS-PPP / PPK ground truth is concentrated in ECEF Z (+2.07m), corresponding to +1.58m North and +1.28m Up at Boulder lat/lon, with East error near zero (-0.06m to +0.20m).
  - Prefit residuals at truth reveal a North/South gradient of 3.25m (South satellites G03, G08 ~ +2.45m; North satellites G09, G26 ~ -0.80m).
  - Investigating the source of this systematic elevation/azimuth bias (satellite PCO modeling, multi-constellation geometry, and epoch duration).
  - Integer PPP-AR LAMBDA fixing is verified working with ratios up to 1.74B.

## Next Steps
1. Determine root cause of the 3.25m North/South prefit residual gradient on F9P.
2. Enable multi-constellation tracking (Galileo/BeiDou) in `epoch.rs` to balance satellite geometry.
3. Validate kinematic trajectory accuracy against CSRS-PPP and PPK ground truth.
4. Verify AGENTS.md compliance, run tests & clippy, and report back to parent.
