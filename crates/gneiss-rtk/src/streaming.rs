//! Real-time incremental streaming RTK engine.
//!
//! Processes GNSS epochs one-by-one as they arrive from live feeds
//! (serial/USB UBX, NTRIP RTCM3 streams) with time synchronization,
//! base observation buffering, and robust double-difference state updates.

use std::collections::BTreeMap;
use nalgebra::Vector3;

use gneiss_core::coords::ecef_cov_to_enu_std;
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::obs::EpochObs;
use gneiss_core::time::GpsTime;

use crate::estimators::rtk_iekf::GnssRtkIekf;
use crate::post_process::dynamics::ProcessingDynamics;

/// Real-time epoch solution emitted by the streaming engine.
#[derive(Debug, Clone, PartialEq)]
pub struct StreamingEpochSolution {
    pub time: GpsTime,
    pub position_ecef: Vector3<f64>,
    pub velocity_ecef: Option<Vector3<f64>>,
    pub std_east: f64,
    pub std_north: f64,
    pub std_up: f64,
    pub quality: u8,
    pub n_satellites: usize,
}

/// Configuration for real-time streaming RTK.
#[derive(Debug, Clone)]
pub struct StreamingConfig {
    pub base_position: Vector3<f64>,
    pub max_base_age_s: f64,
    pub dynamics: ProcessingDynamics,
    pub enable_glonass: bool,
}

impl Default for StreamingConfig {
    fn default() -> Self {
        Self {
            base_position: Vector3::zeros(),
            max_base_age_s: 2.0,
            dynamics: ProcessingDynamics::Kinematic,
            enable_glonass: false,
        }
    }
}

/// Real-time incremental RTK processing engine.
pub struct StreamingRtkEngine {
    config: StreamingConfig,
    iekf: Option<GnssRtkIekf>,
    ephemerides: Vec<Ephemeris>,
    base_buffer: BTreeMap<u64, EpochObs>,
}

impl StreamingRtkEngine {
    /// Creates a new streaming RTK engine.
    pub fn new(config: StreamingConfig, ephemerides: Vec<Ephemeris>) -> Self {
        Self {
            config,
            iekf: None,
            ephemerides,
            base_buffer: BTreeMap::new(),
        }
    }

    /// Number of loaded ephemerides.
    pub fn ephemerides_len(&self) -> usize {
        self.ephemerides.len()
    }

    /// Update ephemeris catalogue.
    pub fn push_ephemeris(&mut self, ephem: Ephemeris) {
        if let Some(pos) = self.ephemerides.iter().position(|e| e.sat() == ephem.sat()) {
            self.ephemerides[pos] = ephem;
        } else {
            self.ephemerides.push(ephem);
        }
    }

    /// Buffer an incoming base station epoch observation.
    pub fn push_base_obs(&mut self, base_epoch: EpochObs) {
        let tow_ms = (base_epoch.time.tow * 1000.0).round() as u64;
        self.base_buffer.insert(tow_ms, base_epoch);

        // Prune old base epochs older than 30 seconds
        if self.base_buffer.len() > 100 {
            if let Some(&first_key) = self.base_buffer.keys().next() {
                if tow_ms.saturating_sub(first_key) > 30_000 {
                    self.base_buffer.remove(&first_key);
                }
            }
        }
    }

    /// Process a live rover observation epoch and emit a real-time solution.
    pub fn process_rover_obs(&mut self, rover_epoch: &EpochObs) -> Option<StreamingEpochSolution> {
        let tow_ms = (rover_epoch.time.tow * 1000.0).round() as u64;
        let base_obs = self.find_closest_base(tow_ms)?.clone();

        let ephem = self.ephemerides.clone();
        let base_pos = self.config.base_position;
        let iekf = self.get_or_init_iekf(rover_epoch.time);

        let filtered = iekf.process_epoch(rover_epoch, &base_obs, base_pos, &ephem).ok()?;
        let (std_e, std_n, std_u) = ecef_cov_to_enu_std(filtered.position_ecef, filtered.cov_position);

        Some(StreamingEpochSolution {
            time: rover_epoch.time,
            position_ecef: filtered.position_ecef,
            velocity_ecef: filtered.velocity_ecef,
            std_east: std_e,
            std_north: std_n,
            std_up: std_u,
            quality: filtered.quality,
            n_satellites: filtered.n_satellites,
        })
    }

    fn find_closest_base(&self, rover_tow_ms: u64) -> Option<&EpochObs> {
        let max_dt_ms = (self.config.max_base_age_s * 1000.0) as u64;
        let mut best: Option<(&u64, &EpochObs)> = None;
        let mut min_diff = u64::MAX;

        for (b_tow, b_epoch) in &self.base_buffer {
            let diff = (*b_tow as i64 - rover_tow_ms as i64).unsigned_abs();
            if diff <= max_dt_ms && diff < min_diff {
                min_diff = diff;
                best = Some((b_tow, b_epoch));
            }
        }
        best.map(|(_, epoch)| epoch)
    }

    fn get_or_init_iekf(&mut self, start_time: GpsTime) -> &mut GnssRtkIekf {
        if self.iekf.is_none() {
            let q_accel = if self.config.dynamics.is_kinematic() { 1.0 } else { 1e-6 };
            let mut iekf = GnssRtkIekf::new(self.config.base_position, start_time, q_accel);
            iekf.enable_glonass = self.config.enable_glonass;
            self.iekf = Some(iekf);
        }
        match self.iekf.as_mut() {
            Some(i) => i,
            None => unreachable!("initialized above"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::obs::SatObs;
    use gneiss_core::sat::{Constellation, SatelliteId};

    fn dummy_epoch(tow: f64) -> EpochObs {
        EpochObs {
            time: GpsTime { week: 2100, tow },
            satellites: vec![
                SatObs {
                    sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
                    observations: vec![],
                }
            ],
        }
    }

    #[test]
    fn test_streaming_engine_base_buffering_and_matching() {
        let cfg = StreamingConfig::default();
        let mut engine = StreamingRtkEngine::new(cfg, vec![]);

        engine.push_base_obs(dummy_epoch(100.0));
        engine.push_base_obs(dummy_epoch(101.0));

        let rover = dummy_epoch(101.05);
        let base_matched = engine.find_closest_base((rover.time.tow * 1000.0) as u64);
        assert!(base_matched.is_some());
        assert_eq!(base_matched.unwrap().time.tow, 101.0);
    }

    #[test]
    fn test_streaming_engine_out_of_sync_returns_none() {
        let cfg = StreamingConfig {
            max_base_age_s: 1.0,
            ..Default::default()
        };
        let mut engine = StreamingRtkEngine::new(cfg, vec![]);
        engine.push_base_obs(dummy_epoch(100.0));

        // Rover is 5 seconds ahead of base
        let rover = dummy_epoch(105.0);
        let base_matched = engine.find_closest_base((rover.time.tow * 1000.0) as u64);
        assert!(base_matched.is_none());
    }
}
