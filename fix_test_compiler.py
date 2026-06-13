with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "r") as f:
    ppp_fg = f.read()
ppp_fg = ppp_fg.replace("apply_state_vector(&mut state2, &x, state2.covariance.clone());", "let cov2 = state2.covariance.clone();\n        apply_state_vector(&mut state2, &x, cov2);")
with open("crates/gneiss-rtk/src/engine/ppp_fg.rs", "w") as f:
    f.write(ppp_fg)

with open("crates/gneiss-rtk/src/engine/predictor.rs", "r") as f:
    predictor = f.read()

config_old = """        let config = EngineConfig {
            mode: crate::engine::EngineMode::Ppp,
            enabled_constellations: None,
            dynamics_model: DynamicsModel::Static,
            lambda_min_ratio: 3.0,
            lambda_min_subset: 4,
            ar_min_epoch_count: 10,
            ar_min_lock: 10,
            imu_mounting_angles: None,
            process_noise_cb: 100.0,
            process_noise_cd: 10.0,
            process_noise_zwd: 0.1,
            process_noise_amb_float: 1e-4,
            process_noise_amb_fixed: 1e-7,
        };"""

config_new = """        let config_str = r#"{
            "mode": "ppp",
            "process_noise_cb": 100.0,
            "process_noise_cd": 10.0,
            "process_noise_zwd": 0.1,
            "process_noise_amb_float": 1e-4,
            "process_noise_amb_fixed": 1e-7
        }"#;
        let config: crate::engine::EngineConfig = serde_json::from_str(config_str).unwrap();"""
predictor = predictor.replace(config_old, config_new)

with open("crates/gneiss-rtk/src/engine/predictor.rs", "w") as f:
    f.write(predictor)

