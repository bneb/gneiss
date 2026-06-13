with open("crates/gneiss-rtk/src/engine/predictor.rs", "r") as f:
    content = f.read()

config_new = """        let config = EngineConfig {
            mode: crate::engine::EngineMode::Ppp,
            initial_position: None,
            base_position: None,
            base_datum_transform: None,
            imu_to_antenna_lever_arm: [0.0, 0.0, 0.0],
            imu_mounting_angles: None,
            imu_to_nhc_lever_arm: [0.0, 0.0, 0.0],
            enable_nhc: false,
            enable_backward_smoothing: false,
            lambda_min_ratio: 3.0,
            lambda_min_subset: 4,
            enabled_constellations: None,
            raim_pseudorange_outlier_m: 10.0,
            chi_square_pr_threshold: 15.0,
            chi_square_cp_threshold: 15.0,
            nominal_snr_dbhz: 30.0,
            dynamics_model: DynamicsModel::Static,
            doppler_slip_threshold_cycles: 5.0,
            max_reject_count: 3,
            max_base_age_s: 30.0,
            spp_consistency_threshold_m: 10.0,
            initial_ambiguity_variance: 100.0,
            ar_min_epoch_count: 10,
            ar_min_lock: 10,
            process_noise_cb: 100.0,
            process_noise_cd: 10.0,
            process_noise_zwd: 0.1,
            process_noise_amb_float: 1e-4,
            process_noise_amb_fixed: 1e-7,
            tuning: crate::engine::config::EkfTuningConfig::default(),
        };"""

import re
content = re.sub(r'let config_str = r#"\{.*?\}";\n        let config: crate::engine::EngineConfig = serde_json::from_str\(config_str\)\.unwrap\(\);', config_new, content, flags=re.DOTALL)

with open("crates/gneiss-rtk/src/engine/predictor.rs", "w") as f:
    f.write(content)
