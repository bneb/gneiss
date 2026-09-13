//! Kinematic Simulation Tests for Post-Processing RTK/PPK.

#[cfg(test)]
mod tests {
    use nalgebra::Vector3;

    use gneiss_rtk::post_process::{execute_post_process, PostProcessOptions};
    use gneiss_rtk::sim::generator::{
        generate_simulation_dataset, SimulationConfig, TrajectoryProfile,
    };
    use gneiss_rtk::swfg::config::{EngineConfig, RtkConfig};

    fn compute_horizontal_error(pos: Vector3<f64>, truth: Vector3<f64>) -> f64 {
        let llh = gneiss_core::coords::ecef_to_llh(truth);
        let ned = gneiss_core::coords::ecef_to_ned_matrix(llh) * (pos - truth);
        (ned.x * ned.x + ned.y * ned.y).sqrt()
    }

    #[test]
    fn test_post_process_open_sky_kinematic_sub_centimeter() {
        let cfg = SimulationConfig {
            duration_s: 30.0,
            epoch_rate_hz: 1.0,
            profile: TrajectoryProfile::Circular {
                center_offset_ned: Vector3::new(50.0, 50.0, 0.0),
                radius_m: 30.0,
                speed_m_s: 4.0,
            },
            pr_noise_m: 0.15,
            cp_noise_m: 0.002,
            doppler_noise_m_s: 0.02,
            num_satellites: 24,
            outages: Vec::new(),
            cycle_slips: Vec::new(),
            ..Default::default()
        };

        let sim = generate_simulation_dataset(&cfg);
        let engine_cfg = EngineConfig::Rtk(RtkConfig {
            initial_position: Some([cfg.base_ecef.x, cfg.base_ecef.y, cfg.base_ecef.z]),
            ..Default::default()
        });

        let options = PostProcessOptions {
            enable_bidirectional: true,
            base_position: Some(cfg.base_ecef),
            ..Default::default()
        };

        let result = execute_post_process(
            &engine_cfg,
            &sim.ephemerides,
            &sim.rover_epochs,
            Some(&sim.base_epochs),
            None,
            &options,
        )
        .expect("Post-processing should succeed");

        assert_eq!(result.trajectory.len(), 30);

        let mut h_errs = Vec::new();
        let mut fixed_count = 0;

        for (i, epoch) in result.trajectory.iter().enumerate() {
            let truth = sim.truth_positions[i].1;
            let h_err = compute_horizontal_error(epoch.position_ecef, truth);
            h_errs.push(h_err);
            if epoch.quality == 1 {
                fixed_count += 1;
            }
        }

        let rms_h = (h_errs.iter().map(|e| e * e).sum::<f64>() / h_errs.len() as f64).sqrt();
        println!("Open-sky Kinematic Post-Processed RTK: RMS={:.4}m, Fixed={}/{}", rms_h, fixed_count, result.trajectory.len());

        assert!(fixed_count >= 25, "At least 25/30 epochs should be fixed, got {}", fixed_count);
        assert!(rms_h < 0.010, "Horizontal RMS error must be < 1.0cm, got {:.4}m", rms_h);
    }

    #[test]
    fn test_post_process_cycle_slip_recovery() {
        let slips = vec![(10.0, 2, 15), (20.0, 5, -8)];

        let cfg = SimulationConfig {
            duration_s: 30.0,
            epoch_rate_hz: 1.0,
            profile: TrajectoryProfile::Linear {
                start_offset_ned: Vector3::new(10.0, 0.0, 0.0),
                velocity_ned: Vector3::new(2.0, 1.0, 0.0),
            },
            pr_noise_m: 0.15,
            cp_noise_m: 0.002,
            doppler_noise_m_s: 0.02,
            num_satellites: 24,
            outages: Vec::new(),
            cycle_slips: slips,
            ..Default::default()
        };

        let sim = generate_simulation_dataset(&cfg);
        let engine_cfg = EngineConfig::Rtk(RtkConfig {
            initial_position: Some([cfg.base_ecef.x, cfg.base_ecef.y, cfg.base_ecef.z]),
            ..Default::default()
        });

        let options = PostProcessOptions {
            enable_bidirectional: true,
            base_position: Some(cfg.base_ecef),
            ..Default::default()
        };

        let result = execute_post_process(
            &engine_cfg,
            &sim.ephemerides,
            &sim.rover_epochs,
            Some(&sim.base_epochs),
            None,
            &options,
        )
        .expect("Post-processing with cycle slips should succeed");

        let mut h_errs = Vec::new();
        for (i, epoch) in result.trajectory.iter().enumerate() {
            let truth = sim.truth_positions[i].1;
            let h_err = compute_horizontal_error(epoch.position_ecef, truth);
            h_errs.push(h_err);
        }

        let rms_h = (h_errs.iter().map(|e| e * e).sum::<f64>() / h_errs.len() as f64).sqrt();
        println!("Cycle Slip Post-Processed RTK: RMS={:.4}m", rms_h);
        assert!(rms_h < 0.010, "Cycle slip post-processing RMS must be < 1.0cm, got {:.4}m", rms_h);
    }

    #[test]
    fn test_post_process_outage_continuity() {
        // 5s outage on satellites 1, 2, 3 during t in [10.0, 15.0]
        let outages = vec![(10.0, 15.0, vec![1, 2, 3])];

        let cfg = SimulationConfig {
            duration_s: 25.0,
            epoch_rate_hz: 1.0,
            profile: TrajectoryProfile::Linear {
                start_offset_ned: Vector3::new(0.0, 0.0, 0.0),
                velocity_ned: Vector3::new(1.0, 0.5, 0.0),
            },
            pr_noise_m: 0.15,
            cp_noise_m: 0.002,
            doppler_noise_m_s: 0.02,
            num_satellites: 24,
            outages,
            cycle_slips: Vec::new(),
            ..Default::default()
        };

        let sim = generate_simulation_dataset(&cfg);
        let engine_cfg = EngineConfig::Rtk(RtkConfig {
            initial_position: Some([cfg.base_ecef.x, cfg.base_ecef.y, cfg.base_ecef.z]),
            ..Default::default()
        });

        let options = PostProcessOptions {
            enable_bidirectional: true,
            base_position: Some(cfg.base_ecef),
            ..Default::default()
        };

        let result = execute_post_process(
            &engine_cfg,
            &sim.ephemerides,
            &sim.rover_epochs,
            Some(&sim.base_epochs),
            None,
            &options,
        )
        .expect("Post-processing with outage should succeed");

        let mut h_errs = Vec::new();
        for (i, epoch) in result.trajectory.iter().enumerate() {
            let truth = sim.truth_positions[i].1;
            let h_err = compute_horizontal_error(epoch.position_ecef, truth);
            h_errs.push(h_err);
        }

        let rms_h = (h_errs.iter().map(|e| e * e).sum::<f64>() / h_errs.len() as f64).sqrt();
        println!("Outage Post-Processed RTK: RMS={:.4}m", rms_h);
        assert!(rms_h < 0.015, "Outage post-processing RMS must be < 1.5cm, got {:.4}m", rms_h);
    }
}
