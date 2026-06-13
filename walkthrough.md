# Phase 1 & 4 Bug Hunt: Subtle Math Bugs

During a deeper audit of the math modeling as requested, I identified and fixed two subtle but critical errors in the measurement equations that significantly impact velocity and atmospheric estimation:

## 1. Doppler Velocity & Receiver Clock Drift (`doppler.rs`)

**The Bug:** The Doppler range-rate observation model (`compute_doppler_update`) was correctly computing the satellite relative velocity but completely **ignored the receiver's clock drift state**. Since clock drift in a standard TCXO can be 1-2 ppm (equivalent to ~300-600 m/s of range rate), any unmodeled clock drift was directly bleeding into the EKF's velocity estimates. Furthermore, the `H` matrix row for Doppler was missing the `1.0` at index 16 (clock drift).

**The Fix:**
*   Added `state.rcv_clk_drift` to `predicted_range_rate`.
*   Added `h_row[16] = 1.0` if `state_size > 16`.
*   Updated the tests in `doppler.rs` to expect the clock drift index in the `H` matrix.

## 2. ZWD RTK Measurement Omission (`measurement.rs`)

**The Bug:** The RTK engine utilizes Double Differences (DD) to cancel clock errors, and estimates Zenith Wet Delay (ZWD) as a core state (index 17). However, `compute_dd_pseudorange` and `compute_dd_carrier_phase` were **completely ignoring the ZWD state**. The `H` matrix for index 17 was hardcoded to `0.0`, and the `comp_pr_dd` predicted pseudorange didn't apply the ZWD. This meant the EKF was *never* updating its ZWD estimate from measurements, essentially running without atmospheric refinement.

**The Fix:**
*   Computed the mapping functions `m_w_rov_sat` and `m_w_rov_ref`.
*   Propagated `zwd_dd = (m_w_rov_sat - m_w_rov_ref) * state.zwd` into the predicted range equations.
*   Updated the `H` matrix at index 17: `h[17] = m_w_rov_sat - m_w_rov_ref`.
*   Updated the RTKLIB golden data unit test, since the properly modeled `zwd = 0.1` default initialization now properly alters the predicted innovations.

---

The benchmarks (`task-505`) are still running in the background. With these measurement fixes and the earlier smoother fixes, the EKF should now have completely sound mathematical foundations for high-precision convergence.
