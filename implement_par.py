import re

content = open("crates/gneiss-rtk/src/engine/ppp_fg.rs").read()

old_func = """    fn resolve_widelane_ar(&self, state: &RtkState, subset: &[((gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64), (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64))], x: &DVector<f64>) -> Result<(DVector<f64>, DMatrix<f64>), &'static str> {
        let mut d_wl = DMatrix::zeros(subset.len(), state.covariance.nrows());
        for (i, (c, ref_sat)) in subset.iter().enumerate() {
            d_wl[(i, CORE_STATE_SIZE + c.1)] = 1.0 / c.4;
            d_wl[(i, CORE_STATE_SIZE + c.2)] = -1.0 / c.5;
            d_wl[(i, CORE_STATE_SIZE + ref_sat.1)] = -1.0 / ref_sat.4;
            d_wl[(i, CORE_STATE_SIZE + ref_sat.2)] = 1.0 / ref_sat.5;
        }

        let a_wl = &d_wl * x;
        let q_wl = &d_wl * &state.covariance * d_wl.transpose();
        
        let mut diag = Vec::new();
        for i in 0..q_wl.nrows() {
            diag.push(q_wl[(i, i)].sqrt());
        }
        tracing::info!("WL Float Ambiguities: {:?}", a_wl.as_slice());
        tracing::info!("WL StdDevs: {:?}", diag);

        let (fixed_z, ratio) = gneiss_lambda::lambda_reduction_and_search(&a_wl, &q_wl, 2);
        
        if ratio < 2.0 {
            tracing::info!("WL ratio test failed. ratio={:.2}, success={:.2}", ratio, gneiss_lambda::success_rate(&q_wl));
            return Err("WL ratio test failed");
        }"""

new_func = """    fn resolve_widelane_ar(&self, state: &RtkState, subset: &[((gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64), (gneiss_core::sat::SatelliteId, usize, usize, f64, f64, f64))], x: &DVector<f64>) -> Result<(DVector<f64>, DMatrix<f64>), &'static str> {
        let mut d_wl_full = DMatrix::zeros(subset.len(), state.covariance.nrows());
        for (i, (c, ref_sat)) in subset.iter().enumerate() {
            d_wl_full[(i, CORE_STATE_SIZE + c.1)] = 1.0 / c.4;
            d_wl_full[(i, CORE_STATE_SIZE + c.2)] = -1.0 / c.5;
            d_wl_full[(i, CORE_STATE_SIZE + ref_sat.1)] = -1.0 / ref_sat.4;
            d_wl_full[(i, CORE_STATE_SIZE + ref_sat.2)] = 1.0 / ref_sat.5;
        }

        let mut a_wl_full = &d_wl_full * x;
        let mut q_wl_full = &d_wl_full * &state.covariance * d_wl_full.transpose();
        
        let mut diag = Vec::new();
        for i in 0..q_wl_full.nrows() {
            diag.push(q_wl_full[(i, i)].sqrt());
        }
        tracing::info!("WL Float Ambiguities: {:?}", a_wl_full.as_slice());
        tracing::info!("WL StdDevs: {:?}", diag);

        // Simple Partial Ambiguity Resolution (PAR): drop ambiguities with stddev > 0.15 cycles
        let mut keep_indices = Vec::new();
        for i in 0..q_wl_full.nrows() {
            if q_wl_full[(i, i)].sqrt() < 0.15 {
                keep_indices.push(i);
            }
        }
        if keep_indices.len() < 4 {
            return Err("Insufficient well-converged Widelane ambiguities");
        }

        let mut d_wl = DMatrix::zeros(keep_indices.len(), state.covariance.nrows());
        for (i, &idx) in keep_indices.iter().enumerate() {
            for j in 0..state.covariance.nrows() {
                d_wl[(i, j)] = d_wl_full[(idx, j)];
            }
        }
        
        let a_wl = &d_wl * x;
        let q_wl = &d_wl * &state.covariance * d_wl.transpose();

        let (fixed_z, ratio) = gneiss_lambda::lambda_reduction_and_search(&a_wl, &q_wl, 2);
        
        if ratio < 2.0 {
            tracing::info!("WL ratio test failed. ratio={:.2}, success={:.2}", ratio, gneiss_lambda::success_rate(&q_wl));
            return Err("WL ratio test failed");
        }"""

if old_func in content:
    content = content.replace(old_func, new_func)
    open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "w").write(content)
    print("Fixed PAR!")
else:
    print("Could not find the block to replace!")
