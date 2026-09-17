# Task Assignment — Worker Milestone M1: 15-State ESKF/MEKF GNSS/INS

## Role & Mission
You are an implementation worker (`teamwork_preview_worker`).
Your working directory is: `/Users/kevin/projects/gneiss/.agents/worker_m1/`.
The authoritative user request is: `/Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md`.
The project specification is: `/Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md`.
The survey report is: `/Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r1/survey_r1.md`.

## MANDATORY INTEGRITY WARNING
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

## Code Quality Standards (AGENTS.md)
- File size strictly < 500 LOC
- Function size strictly < 32 LOC
- Nesting depth strictly < 3 levels
- Exactly 0 `unwrap()` calls in production code (`match`, `if let`, `ok_or()?`, or descriptive `.expect()` only)
- Zero clippy warnings (`cargo clippy --workspace --all-targets -- -D warnings`)
- Passing test coverage with comprehensive unit tests in `#[cfg(test)] mod tests`

## Exclusive Write Ownership
You own ONLY:
- `crates/gneiss-rtk/src/estimators/eskf/` (`mod.rs`, `types.rs`, `predict.rs`, `update.rs`, `constraints.rs`, `smoother.rs`)
- `crates/gneiss-rtk/src/estimators/mod.rs` (to register `pub mod eskf;`)
- `crates/gneiss-rtk/src/bin/eval_odaiba_ins.rs`
Do NOT edit any files outside this set.

## Implementation Requirements
Follow the exact mathematical equations from `survey_r1.md`:
1. **Types (`types.rs`)**:
   - `EskfState`: `pos_ecef` ($3$), `vel_ecef` ($3$), `attitude: UnitQuaternion<f64>`, `accel_bias` ($3$), `gyro_bias` ($3$), `cov: Matrix15<f64>`.
   - Layout: $\delta\mathbf{p}^e (0..3)$, $\delta\mathbf{v}^e (3..6)$, $\delta\boldsymbol{\theta} (6..9)$, $\delta\mathbf{b}_a (9..12)$, $\delta\mathbf{b}_g (12..15)$.
2. **Prediction (`predict.rs`)**:
   - Nominal attitude propagation: $\mathbf{q}_{k} = \mathbf{q}_{k-1} \otimes \Delta\mathbf{q}(\boldsymbol{\omega}_m - \mathbf{b}_g)$.
   - Discrete state transition $\boldsymbol{\Phi}_{15\times 15}$:
     - Position block: $\mathbf{I} + \mathbf{I}\Delta t$
     - Velocity block: $\mathbf{I} - 2[\boldsymbol{\omega}_{ie}^e\times]\Delta t$
     - **Velocity-attitude coupling**: `vel_att = +f_e_skew * dt` (MANDATORY POSITIVE SIGN per `AGENTS.md` and `ORIGINAL_REQUEST.md:90-106`).
     - Velocity-accel-bias coupling: $-\mathbf{R}_b^e \Delta t$.
     - Attitude-gyro-bias coupling: $-\mathbf{R}_b^e \Delta t$.
   - Spectral process noise $\mathbf{Q}_{15\times 15}$ using calibrated IMU noise densities.
3. **Updates (`update.rs`)**:
   - GNSS position and velocity innovation: $\mathbf{z}_p = \mathbf{p}_{gnss} - (\mathbf{p}^e + \mathbf{R}_b^e \mathbf{l}^b)$, $\mathbf{z}_v = \mathbf{v}_{gnss} - (\mathbf{v}^e + \mathbf{R}_b^e (\boldsymbol{\omega}\times\mathbf{l}^b))$.
   - Measurement Jacobians with lever arm $-[(\mathbf{R}_b^e\mathbf{l}^b)\times]$.
   - Kalman gain, Joseph-form covariance update: $\mathbf{P} = (\mathbf{I} - \mathbf{K}\mathbf{H})\mathbf{P}(\mathbf{I} - \mathbf{K}\mathbf{H})^T + \mathbf{K}\mathbf{R}\mathbf{K}^T$.
   - **Error state reset & attitude injection**: $\mathbf{p} \leftarrow \mathbf{p} + \delta\mathbf{p}$, $\mathbf{v} \leftarrow \mathbf{v} + \delta\mathbf{v}$, $\mathbf{q} \leftarrow \mathbf{q} \otimes \Delta\mathbf{q}(\delta\boldsymbol{\theta})$, $\mathbf{b}_a \leftarrow \mathbf{b}_a + \delta\mathbf{b}_a$, $\mathbf{b}_g \leftarrow \mathbf{b}_g + \delta\mathbf{b}_g$.
4. **Constraints (`constraints.rs`)**:
   - Coupled NHC: Lateral and vertical body velocities $\mathbf{z}_{nhc} = - \mathbf{E}_{23} (\mathbf{R}_b^e)^T \mathbf{v}^e$.
   - **Attitude coupling Jacobian**: $\mathbf{H}_{nhc,\theta} = \mathbf{E}_{23} (\mathbf{R}_b^e)^T [\mathbf{v}^e\times]$ (crucial for heading observability).
   - Zero-Velocity Updates (ZUPT) when stationary.
5. **RTS Smoother (`smoother.rs`)**:
   - Full 15-state backward Rauch-Tung-Striebel smoother: $\mathbf{C}_k = \mathbf{P}_{k|k} \boldsymbol{\Phi}_{k+1}^T \mathbf{P}_{k+1|k}^{-1}$.
   - Backward state and covariance recursion propagating corrections across all 15 states.
6. **Benchmark Validation (`eval_odaiba_ins.rs`)**:
   - Integrate 15-state ESKF with forward filtering and backward RTS smoothing.
   - Run `cargo run --release --bin eval_odaiba_ins`.
   - Verify: $p_{50} < 2.5\text{ m}$ and $\text{RMS} < 5.2\text{ m}$ across the full 12,398 epochs.

## Verification Commands
Before reporting completion, run:
```bash
cargo build --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p gneiss-rtk --lib estimators::eskf
cargo run --release --bin eval_odaiba_ins
```
Record all test and benchmark outputs in your `handoff.md`.
Send a completion message when done.

## 2026-09-12T16:48:52Z
You are Worker M1.
Your working directory is: /Users/kevin/projects/gneiss/.agents/worker_m1/
Read your task instructions in /Users/kevin/projects/gneiss/.agents/worker_m1/DISPATCH.md, the project scope in /Users/kevin/projects/gneiss/.agents/orchestrator_frontiers/PROJECT.md, the survey report in /Users/kevin/projects/gneiss/.agents/teamwork_preview_explorer_survey_r1/survey_r1.md, and the authoritative request in /Users/kevin/projects/gneiss/.agents/ORIGINAL_REQUEST.md.
Implement Frontier R1 (15-State ESKF/MEKF GNSS/INS, RTS Smoother, coupled NHC/ZUPT, eval_odaiba_ins benchmark) adhering strictly to AGENTS.md.
Run your verification commands and record all output in your handoff.md. Send a completion message when done.
