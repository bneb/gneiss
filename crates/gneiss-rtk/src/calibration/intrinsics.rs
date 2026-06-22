use crate::engine::config::EngineConfig;
use gneiss_core::obs::EpochObs;

/// Analyzes receiver clock bias/drift from raw observations to tune EKF parameters.
/// Returns (clock_bias_variance, clock_drift_variance).
pub fn calibrate_intrinsics(config: &EngineConfig, obs: &[EpochObs]) -> (f64, f64) {
    if obs.len() < 2 {
        return (config.process_noise_cb, config.process_noise_cd);
    }

    let mut drift_estimates = Vec::new();

    for i in 1..obs.len() {
        let prev_epoch = &obs[i - 1];
        let curr_epoch = &obs[i];

        let dt = curr_epoch.time - prev_epoch.time;
        if dt <= 0.0 || dt > 10.0 {
            continue;
        }

        let mut pr_diffs = Vec::new();

        for curr_sat in &curr_epoch.satellites {
            if let Some(prev_sat) = prev_epoch.satellites.iter().find(|s| s.sat == curr_sat.sat) {
                if let (Some(curr_pr), Some(prev_pr)) =
                    (curr_sat.get_observable(1), prev_sat.get_observable(1))
                {
                    // Simple rate of change of pseudorange (m/s)
                    let diff = (curr_pr - prev_pr) / dt;
                    pr_diffs.push(diff);
                }
            }
        }

        if !pr_diffs.is_empty() {
            // The average PR change across all satellites is a rough estimate of the clock drift for this epoch
            let avg_drift = pr_diffs.iter().sum::<f64>() / pr_diffs.len() as f64;
            drift_estimates.push(avg_drift);
        }
    }

    if drift_estimates.len() < 2 {
        return (config.process_noise_cb, config.process_noise_cd);
    }

    let mean_drift = drift_estimates.iter().sum::<f64>() / drift_estimates.len() as f64;
    let mut drift_var = drift_estimates
        .iter()
        .map(|d| {
            let diff = d - mean_drift;
            diff * diff
        })
        .sum::<f64>()
        / (drift_estimates.len() - 1) as f64;

    // Scale it slightly to account for unmodeled satellite motion
    drift_var *= 1.5;

    // Apply bounds based on typical TCXO/OCXO characteristics
    if drift_var < 1.0 {
        drift_var = 1.0;
    }
    if drift_var > config.process_noise_cd * 10.0 {
        drift_var = config.process_noise_cd * 10.0;
    }

    // Bias variance is loosely coupled to drift variance, but typically much larger
    let mut bias_var = drift_var * 100.0;

    if bias_var < config.process_noise_cb * 0.1 {
        bias_var = config.process_noise_cb * 0.1;
    }
    if bias_var > config.process_noise_cb * 10.0 {
        bias_var = config.process_noise_cb * 10.0;
    }

    (bias_var, drift_var)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::obs::{ObsCode, ObsType, Observation, SatObs, SignalCode};
    use gneiss_core::sat::{Constellation, SatelliteId};
    use gneiss_core::time::GpsTime;

    fn make_obs(time: f64, pr1: f64, pr2: f64) -> EpochObs {
        let code = ObsCode {
            obs_type: ObsType::Pseudorange,
            signal: SignalCode {
                freq_band: 1,
                attribute: 'C',
            },
        };

        EpochObs {
            time: GpsTime::new(2200, time),
            satellites: vec![
                SatObs {
                    sat: SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 1,
                    },
                    observations: vec![Observation {
                        code,
                        value: pr1,
                        lock_time: None,
                        lli: None,
                    }],
                },
                SatObs {
                    sat: SatelliteId {
                        constellation: Constellation::Gps,
                        prn: 2,
                    },
                    observations: vec![Observation {
                        code,
                        value: pr2,
                        lock_time: None,
                        lli: None,
                    }],
                },
            ],
        }
    }

    #[test]
    fn test_calibrate_intrinsics_insufficient_data() {
        let config = EngineConfig::default();
        let obs = vec![make_obs(0.0, 20000000.0, 21000000.0)];
        let (cb, cd) = calibrate_intrinsics(&config, &obs);
        assert_eq!(cb, config.process_noise_cb);
        assert_eq!(cd, config.process_noise_cd);
    }

    #[test]
    fn test_calibrate_intrinsics_basic_drift() {
        let config = EngineConfig::default();
        let mut obs = Vec::new();

        let mut pr1 = 20000000.0;
        let mut pr2 = 21000000.0;

        for i in 0..10 {
            obs.push(make_obs(i as f64, pr1, pr2));
            pr1 += 100.0 + (i as f64 * 1.5);
            pr2 += 100.0 + (i as f64 * 0.5);
        }

        let (cb, cd) = calibrate_intrinsics(&config, &obs);
        assert!(cd > 1.0);
        assert!(cb > cd);
    }

    #[test]
    fn test_calibrate_intrinsics_large_gap_skipped() {
        let config = EngineConfig::default();
        let mut obs = Vec::new();

        // Gap of 15 seconds between epochs (>10 -> skipped)
        obs.push(make_obs(0.0, 20000000.0, 21000000.0));
        obs.push(make_obs(15.0, 20000100.0, 21000100.0));
        obs.push(make_obs(30.0, 20000200.0, 21000200.0));

        let (cb, cd) = calibrate_intrinsics(&config, &obs);
        // All drift estimates skipped due to large dt -> fewer than 2 drift_estimates
        // So fallback to config defaults
        assert_eq!(cb, config.process_noise_cb);
        assert_eq!(cd, config.process_noise_cd);
    }

    #[test]
    fn test_calibrate_intrinsics_no_common_satellites() {
        use gneiss_core::obs::{ObsCode, ObsType, Observation, SatObs, SignalCode};
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;

        let config = EngineConfig::default();

        let code = ObsCode {
            obs_type: ObsType::Pseudorange,
            signal: SignalCode { freq_band: 1, attribute: 'C' },
        };

        // Two epochs with different satellites
        let obs = vec![
            EpochObs {
                time: GpsTime::new(2200, 0.0),
                satellites: vec![SatObs {
                    sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
                    observations: vec![Observation { code, value: 20000000.0, lock_time: None, lli: None }],
                }],
            },
            EpochObs {
                time: GpsTime::new(2200, 1.0),
                satellites: vec![SatObs {
                    sat: SatelliteId { constellation: Constellation::Gps, prn: 2 },
                    observations: vec![Observation { code: code.clone(), value: 20000100.0, lock_time: None, lli: None }],
                }],
            },
        ];

        let (cb, cd) = calibrate_intrinsics(&config, &obs);
        // No common satellites -> pr_diffs is empty -> fewer than 2 drift_estimates
        assert_eq!(cb, config.process_noise_cb);
        assert_eq!(cd, config.process_noise_cd);
    }

    #[test]
    fn test_calibrate_intrinsics_bias_variance_clamping() {
        let config = EngineConfig::default();
        let mut obs = Vec::new();

        // Create very consistent drift (all PR diffs identical) -> drift_var = 0
        // but then it's clamped to 1.0
        for i in 0..5 {
            obs.push(make_obs(i as f64, 20000000.0 + (i as f64) * 100.0, 21000000.0 + (i as f64) * 100.0));
        }

        let (_cb, cd) = calibrate_intrinsics(&config, &obs);
        // drift_var should be clamped to min 1.0
        assert!(cd >= 1.0);
    }

    #[test]
    fn test_calibrate_intrinsics_drift_var_upper_bound() {
        let config = EngineConfig::default();
        let mut obs = Vec::new();

        // Create very erratic drift to push drift_var above the upper bound
        // Upper bound = config.process_noise_cd * 10.0
        for i in 0..5 {
            let pr1 = 20000000.0 + (i as f64 * 5000.0); // huge jumps
            let pr2 = 21000000.0 + (i as f64 * 5000.0);
            obs.push(make_obs(i as f64, pr1, pr2));
        }

        let (_cb, cd) = calibrate_intrinsics(&config, &obs);
        let upper_bound = config.process_noise_cd * 10.0;
        assert!(cd <= upper_bound, "drift_var {:.2} should be capped at {:.2}", cd, upper_bound);
    }

    #[test]
    fn test_calibrate_intrinsics_bias_var_lower_bound() {
        let config = EngineConfig::default();
        let mut obs = Vec::new();

        // drift_var = 0.0 (all identical), after scaling drift_var = 0.0
        // clamped to 1.0. Then bias_var = 1.0 * 100 = 100.0
        // Lower bound = config.process_noise_cb * 0.1
        for i in 0..5 {
            obs.push(make_obs(i as f64, 20000000.0 + i as f64, 21000000.0 + i as f64));
        }

        let (cb, _cd) = calibrate_intrinsics(&config, &obs);
        let lower_bound = config.process_noise_cb * 0.1;
        assert!(cb >= lower_bound, "bias_var {:.2} should be at least {:.2}", cb, lower_bound);
    }
}
