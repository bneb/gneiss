import re

with open("crates/gneiss-rtk/src/engine/ppp_ins_fg.rs", "r") as f:
    content = f.read()

# 1. Update solver.solve signature to take lever_arm
content = content.replace(
    """    match solver.solve(
        state,
        &rover_smoothed,
        engine.imu_history.as_slice(),
        engine.state_history.last(),
    ) {""",
    """    let lever_arm = nalgebra::Vector3::from_column_slice(&engine.config.imu_to_antenna_lever_arm);
    match solver.solve(
        state,
        &rover_smoothed,
        engine.imu_history.as_slice(),
        engine.state_history.last(),
        &lever_arm,
    ) {"""
)

# 2. Update solve method signature
content = content.replace(
    """    pub fn solve(
        &self,
        state: &mut RtkState,
        sats: &[ProcessedSat],
        imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
        last_state: Option<&RtkState>,
    ) -> Result<bool, EngineError> {""",
    """    pub fn solve(
        &self,
        state: &mut RtkState,
        sats: &[ProcessedSat],
        imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
        last_state: Option<&RtkState>,
        lever_arm: &nalgebra::Vector3<f64>,
    ) -> Result<bool, EngineError> {"""
)

# 3. Add lever_arm to compute_iteration_dx calls
content = content.replace(
    "self.compute_iteration_dx(state, sats, &x_i, &x_pred, &p_inv, _iter, imu_history, last_state)?",
    "self.compute_iteration_dx(state, sats, &x_i, &x_pred, &p_inv, _iter, imu_history, last_state, lever_arm)?"
)

content = content.replace(
    """    fn compute_iteration_dx(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
        x_i: &DVector<f64>,
        x_pred: &DVector<f64>,
        p_inv: &DMatrix<f64>,
        iter: usize,
        imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
        last_state: Option<&RtkState>,
    ) -> Result<Option<DVector<f64>>, EngineError> {""",
    """    fn compute_iteration_dx(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
        x_i: &DVector<f64>,
        x_pred: &DVector<f64>,
        p_inv: &DMatrix<f64>,
        iter: usize,
        imu_history: &[Vec<gneiss_core::imu::ImuMeasurement>],
        last_state: Option<&RtkState>,
        lever_arm: &nalgebra::Vector3<f64>,
    ) -> Result<Option<DVector<f64>>, EngineError> {"""
)

content = content.replace(
    """        let meas = self.build_measurements(state, sats, x_i, iter, &omega_eb_b);""",
    """        let meas = self.build_measurements(state, sats, x_i, iter, &omega_eb_b, lever_arm);"""
)

# 4. Add lever_arm to build_measurements signature
content = content.replace(
    """    fn build_measurements(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
        x_i: &DVector<f64>,
        iter: usize,
        omega_eb_b: &nalgebra::Vector3<f64>,
    ) -> Vec<FgMeasurement> {""",
    """    fn build_measurements(
        &self,
        state: &RtkState,
        sats: &[ProcessedSat],
        x_i: &DVector<f64>,
        iter: usize,
        omega_eb_b: &nalgebra::Vector3<f64>,
        lever_arm: &nalgebra::Vector3<f64>,
    ) -> Vec<FgMeasurement> {"""
)

# 5. Fix lever_arm reading in build_measurements
content = content.replace(
    """        let lever_arm = nalgebra::Vector3::from_column_slice(&self.config.imu_to_antenna_lever_arm);""",
    ""
)

# 6. Fix `self.build_measurements` missing parameters in solve
content = content.replace(
    """        let final_meas = self.build_measurements(state, sats, x_i, self.max_iterations);""",
    """        let omega_eb_b = imu_history.last()
            .and_then(|buf| buf.last())
            .map(|m| m.gyro - nalgebra::Vector3::new(x_i[12], x_i[13], x_i[14]))
            .unwrap_or(nalgebra::Vector3::zeros());
        let final_meas = self.build_measurements(state, sats, x_i, self.max_iterations, &omega_eb_b, lever_arm);"""
)

content = content.replace(
    """        let last_meas = self.build_measurements(state, sats, x_i, self.max_iterations);""",
    """        let omega_eb_b = imu_history.last()
            .and_then(|buf| buf.last())
            .map(|m| m.gyro - nalgebra::Vector3::new(x_i[12], x_i[13], x_i[14]))
            .unwrap_or(nalgebra::Vector3::zeros());
        let last_meas = self.build_measurements(state, sats, x_i, self.max_iterations, &omega_eb_b, lever_arm);"""
)

content = content.replace(
    """    let meas = solver.build_measurements(state, sats, x_i, solver.max_iterations - 1);""",
    """    let lever_arm = nalgebra::Vector3::from_column_slice(&engine.config.imu_to_antenna_lever_arm);
    let omega_eb_b = engine.imu_history.last()
        .and_then(|buf| buf.last())
        .map(|m| m.gyro - nalgebra::Vector3::new(x_i[12], x_i[13], x_i[14]))
        .unwrap_or(nalgebra::Vector3::zeros());
    let meas = solver.build_measurements(state, sats, x_i, solver.max_iterations - 1, &omega_eb_b, &lever_arm);"""
)


# 7. Add `_h_pos_att` usage or remove warning in push_uduc_measurements
content = content.replace(
    """    fn push_uduc_measurements(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        _iter: usize,
        los: &Vector3<f64>,
        expected_base: f64,
        _dist: f64,
        _isb: f64,
        h_pos_att: &nalgebra::Matrix3<f64>,
    ) {""",
    """    fn push_uduc_measurements(
        &self,
        meas: &mut Vec<FgMeasurement>,
        state: &RtkState,
        sat: &ProcessedSat,
        x_i: &DVector<f64>,
        _iter: usize,
        los: &Vector3<f64>,
        expected_base: f64,
        _dist: f64,
        _isb: f64,
        h_pos_att: &nalgebra::Matrix3<f64>,
    ) {"""
)

# Fix h_pos_att not found in build_h_row_uduc usage in build_measurements? No, it was missing in push_uduc_measurements.
# Wait, `pr_h_att` calculation in `build_h_row_uduc` missing. 
# Ah, the error `E0425: cannot find value h_pos_att in this scope` was at `crates/gneiss-rtk/src/engine/ppp_ins_fg.rs:954:42`.
# `let pr_h_att = los.transpose() * h_pos_att;` inside `build_h_row_uduc`?
# I see.

with open("crates/gneiss-rtk/src/engine/ppp_ins_fg.rs", "w") as f:
    f.write(content)

