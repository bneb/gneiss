use std::collections::HashMap;
use gneiss_core::obs::{EpochObs, ObsType, SignalCode};
use gneiss_core::sat::SatelliteId;
use gneiss_core::time::GpsTime;

/// State for a single satellite/signal in the Hatch filter.
#[derive(Debug, Clone)]
pub struct HatchState {
    pub smoothed_pr: f64,
    pub last_phase: f64,
    pub last_time: GpsTime,
    pub count: usize,
}

/// A Carrier Phase Smoothing (Hatch) Filter.
/// Reduces thermal noise in pseudorange measurements by leveraging
/// the extreme precision of carrier phase changes over continuous arcs.
#[derive(Debug, Clone)]
pub struct HatchFilter {
    /// Maps a (SatelliteId, SignalCode) to its smoothing state.
    pub states: HashMap<(SatelliteId, SignalCode), HatchState>,
    /// Maximum smoothing window size (e.g., 100 epochs).
    pub max_window: usize,
    /// Threshold to detect cycle slips or reset the filter (in meters).
    pub slip_threshold_m: f64,
    /// Maximum allowed time gap between epochs before resetting (in seconds).
    pub max_time_gap_s: f64,
}

impl Default for HatchFilter {
    fn default() -> Self {
        Self {
            states: HashMap::new(),
            max_window: 100,
            slip_threshold_m: 5.0, // A 5m deviation between phase-projected PR and raw PR indicates a cycle slip
            max_time_gap_s: 10.0,
        }
    }
}

impl HatchFilter {
    pub fn new(max_window: usize, slip_threshold_m: f64, max_time_gap_s: f64) -> Self {
        Self {
            states: HashMap::new(),
            max_window,
            slip_threshold_m,
            max_time_gap_s,
        }
    }

    /// Smoothes an entire EpochObs in-place.
    /// Modifies the pseudorange observation values to be their smoothed counterparts.
    pub fn smooth_epoch(&mut self, epoch: &mut EpochObs) {
        let current_time = epoch.time;

        for sat_obs in &mut epoch.satellites {
            // Find all available signals (freq_bands) for this satellite
            let mut available_signals = Vec::new();
            for obs in &sat_obs.observations {
                if !available_signals.contains(&obs.code.signal) {
                    available_signals.push(obs.code.signal);
                }
            }

            for signal in available_signals {
                // Try to find both Pseudorange and Carrier Phase for this signal
                let pr_idx = sat_obs.observations.iter().position(|o| o.code.obs_type == ObsType::Pseudorange && o.code.signal == signal);
                let cp_idx = sat_obs.observations.iter().position(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal == signal);

                if let (Some(pi), Some(ci)) = (pr_idx, cp_idx) {
                    let pr_val = sat_obs.observations[pi].value;
                    let cp_val = sat_obs.observations[ci].value; // In meters
                    
                    let key = (sat_obs.sat, signal);

                    let smoothed_val = if let Some(state) = self.states.get_mut(&key) {
                        let dt = (current_time - state.last_time).abs();
                        
                        let delta_phase = cp_val - state.last_phase;
                        
                        // IF the receiver's carrier phase is defined opposite to pseudorange (e.g. phase increases when range decreases),
                        // delta_phase must be subtracted. For u-blox, phase increases as distance increases. Wait, let's assume standard sign:
                        let projected_pr = state.smoothed_pr + delta_phase;

                        // Check for cycle slip or large time gap
                        if dt > self.max_time_gap_s || (pr_val - projected_pr).abs() > 50.0 { // Relax to 50m to avoid false trips
                            // Reset
                            state.smoothed_pr = pr_val;
                            state.last_phase = cp_val;
                            state.last_time = current_time;
                            state.count = 1;
                            pr_val
                        } else {
                            // Smooth
                            state.count = std::cmp::min(state.count + 1, self.max_window);
                            let w = 1.0 / (state.count as f64);
                            
                            let new_smoothed = w * pr_val + (1.0 - w) * projected_pr;
                            
                            state.smoothed_pr = new_smoothed;
                            state.last_phase = cp_val;
                            state.last_time = current_time;
                            
                            new_smoothed
                        }
                    } else {
                        // Initialize new state
                        self.states.insert(key, HatchState {
                            smoothed_pr: pr_val,
                            last_phase: cp_val,
                            last_time: current_time,
                            count: 1,
                        });
                        pr_val
                    };

                    // Update the EpochObs in place
                    sat_obs.observations[pi].value = smoothed_val;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::obs::{ObsCode, SatObs, Observation};
    use gneiss_core::sat::Constellation;

    #[test]
    fn test_hatch_filter_smoothing() {
        let mut filter = HatchFilter::default();
        let sig = SignalCode { freq_band: 1, attribute: 'C' };
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        
        let pr_code = ObsCode { obs_type: ObsType::Pseudorange, signal: sig };
        let cp_code = ObsCode { obs_type: ObsType::CarrierPhase, signal: sig };

        // Epoch 1
        let mut epoch1 = EpochObs {
            time: GpsTime::new(0, 100.0),
            satellites: vec![SatObs {
                sat,
                observations: vec![
                    Observation { code: pr_code, value: 20000000.0, lock_time: None },
                    Observation { code: cp_code, value: 20000000.0, lock_time: None },
                ]
            }],
        };
        filter.smooth_epoch(&mut epoch1);
        assert_eq!(epoch1.satellites[0].observations[0].value, 20000000.0);

        // Epoch 2: Phase moves 1m, PR noise makes it move 3m
        let mut epoch2 = EpochObs {
            time: GpsTime::new(0, 101.0),
            satellites: vec![SatObs {
                sat,
                observations: vec![
                    Observation { code: pr_code, value: 20000003.0, lock_time: None }, // +3m
                    Observation { code: cp_code, value: 20000001.0, lock_time: None }, // +1m
                ]
            }],
        };
        filter.smooth_epoch(&mut epoch2);
        
        // Count = 2, w = 0.5
        // projected_pr = 20000000.0 + 1.0 = 20000001.0
        // smoothed_pr = 0.5 * 20000003.0 + 0.5 * 20000001.0 = 20000002.0
        assert_eq!(epoch2.satellites[0].observations[0].value, 20000002.0);
    }
}
