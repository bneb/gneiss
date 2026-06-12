import os

replacements = {
    "bin/gneiss-cli/src/evaluator.rs": [
        ('println!("   AEROSPACE METRICS: ERROR CDFs (meters / deg)   ");', 'println!("   ERROR CDFs (meters / deg)   ");')
    ],
    "crates/gneiss-rtk/src/spp.rs": [
        ('// If the variance is huge (e.g. we deweighted satellites and have poor remaining geometry), reject it.', '// If the variance exceeds the threshold (e.g. due to deweighted satellites or poor geometry), reject the measurement.')
    ],
    "crates/gneiss-rtk/src/filter.rs": [
        ('tracing::info!("Multi-Const PAR Fixed! Ratio={:.2}, Ps={:.4}", res.ratio, res.success_rate);', 'tracing::info!("Multi-constellation PAR fixed. Ratio={:.2}, Ps={:.4}", res.ratio, res.success_rate);')
    ],
    "crates/gneiss-rtk/src/factor_graph/mod.rs": [
        ('tracing::debug!("Cholesky FAILED for final covariance!");', 'tracing::debug!("Cholesky decomposition failed for final covariance.");')
    ],
    "crates/gneiss-rtk/src/engine/smoother.rs": [
        ('tracing::warn!("Smoother huge pos correction: {:.1}m at k={}", pos_corr_norm, k_idx);', 'tracing::warn!("Smoother large position correction: {:.1}m at k={}", pos_corr_norm, k_idx);')
    ],
    "crates/gneiss-rtk/src/engine/measurement.rs": [
        ('tracing::debug!("DD Doppler Innov huge! sat={} innov={:.3} obs={:.3} pred={:.3}", rov_sat.sat.to_string(), innov, observed_dd_rr, predicted_dd_rr);', 'tracing::debug!("DD Doppler innovation exceeds threshold. sat={} innov={:.3} obs={:.3} pred={:.3}", rov_sat.sat.to_string(), innov, observed_dd_rr, predicted_dd_rr);'),
        ('tracing::error!("MASSIVE Z PASSED PRE-FILTER! type: {}, z: {:.1}, chi2: {:.1}", type_all[i].1, z_all[i], chi2);', 'tracing::error!("Large innovation passed pre-filter. type: {}, z: {:.1}, chi2: {:.1}", type_all[i].1, z_all[i], chi2);')
    ],
    "crates/gneiss-rtk/src/engine/updater.rs": [
        ('tracing::debug!("EKF rejected Doppler/PR measurement! type={}, nu={:.2}, s_ii={:.2}, r_ii={:.4}", meas_type, nu, s_ii, r_ii);', 'tracing::debug!("EKF rejected Doppler/PR measurement. type={}, nu={:.2}, s_ii={:.2}, r_ii={:.4}", meas_type, nu, s_ii, r_ii);'),
        ('tracing::error!("EKF rejected ALL {} pseudoranges! Force reset to SPP.", total_pr);', 'tracing::error!("EKF rejected all {} pseudoranges. Falling back to SPP.", total_pr);'),
        ('// Phase: Massive ratio. EKF should NOT reject Phase, rely on `handle_cycle_slips`.', '// Phase: High threshold ratio. EKF does not reject phase directly; relies on `handle_cycle_slips`.')
    ],
    "crates/gneiss-rtk/src/engine/ppp.rs": [
        ('tracing::info!("Detected clock jump! Snapping clock state from {:.2} to {:.2}", st.rcv_clk_bias, st.rcv_clk_bias + clock_diff);', 'tracing::info!("Detected clock jump. Adjusting clock state from {:.2} to {:.2}", st.rcv_clk_bias, st.rcv_clk_bias + clock_diff);')
    ],
    "crates/gneiss-rtk/src/engine/ppp_fg.rs": [
        ('tracing::info!("Detected clock jump! Snapping clock state from {:.2} to {:.2}", st.rcv_clk_bias, st.rcv_clk_bias + clock_diff);', 'tracing::info!("Detected clock jump. Adjusting clock state from {:.2} to {:.2}", st.rcv_clk_bias, st.rcv_clk_bias + clock_diff);'),
        ('// We can approximate the phase jump by checking if the carrier phase residual would be massive.', '// We can approximate the phase jump by checking if the carrier phase residual exceeds the nominal range.'),
        ('tracing::warn!("Rejecting {:.2}m jump! med_res={:.2}m, sats={}", jump, med_res, sats_to_process.len());', 'tracing::warn!("Rejecting {:.2}m clock jump. med_res={:.2}m, sats={}", jump, med_res, sats_to_process.len());')
    ],
    "crates/gneiss-rtk/src/engine/mod.rs": [
        ('tracing::warn!("compute_spp failed in process_rtk!");', 'tracing::warn!("compute_spp failed in process_rtk.");'),
        ('tracing::debug!("Integer ambiguities resolved!");', 'tracing::debug!("Integer ambiguities resolved.");'),
        ('tracing::warn!("Smoother huge pos correction: {:.1}m at k={}. x_k1_n pos: {:.1}, x_pred_k1 pos: {:.1}", pos_corr_norm, k, x_k1_n.fixed_rows::<3>(0).norm(), x_pred_k1_sub.fixed_rows::<3>(0).norm());', 'tracing::warn!("Smoother large position correction: {:.1}m at k={}. x_k1_n pos: {:.1}, x_pred_k1 pos: {:.1}", pos_corr_norm, k, x_k1_n.fixed_rows::<3>(0).norm(), x_pred_k1_sub.fixed_rows::<3>(0).norm());')
    ],
    "crates/gneiss-rtk/src/engine/tests_updater.rs": [
        ('// Massive clock jump', '// Large clock jump')
    ],
    "crates/gneiss-rtk/src/engine/tight_fg.rs": [
        ('tracing::info!("Detected clock jump! Snapping clock state from {:.2} to {:.2}", st.rcv_clk_bias, st.rcv_clk_bias + clock_diff);', 'tracing::info!("Detected clock jump. Adjusting clock state from {:.2} to {:.2}", st.rcv_clk_bias, st.rcv_clk_bias + clock_diff);'),
        ('// If Carrier Phase has a massive cycle slip, it will have a huge residual here and get Cauchy-rejected!', '// If Carrier Phase has a cycle slip, it will produce a large residual and get Cauchy-rejected.')
    ]
}

for filepath, reps in replacements.items():
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        continue
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()
    
    for old, new in reps:
        if old in content:
            content = content.replace(old, new)
        else:
            print(f"Warning: could not find target text in {filepath}:\n{old}")
            
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)

print("Slop text replacements complete.")
