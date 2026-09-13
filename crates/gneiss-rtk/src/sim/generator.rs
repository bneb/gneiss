//! GNSS Physical Observation & Trajectory Simulator.
//!
//! Generates multi-constellation, multi-frequency base and rover observations
//! with realistic orbital geometry, true carrier phase integer ambiguities,
//! cycle slips, Doppler, and signal dropouts.

use std::f64::consts::PI;
use nalgebra::Vector3;

use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::ephemeris::{Ephemeris, GpsEphemeris};
use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;

/// Trajectory profile for rover simulation.
#[derive(Debug, Clone)]
pub enum TrajectoryProfile {
    /// Stationary rover at fixed offset from base.
    Static { offset_ned: Vector3<f64> },
    /// Kinematic circular motion around a center point.
    Circular {
        center_offset_ned: Vector3<f64>,
        radius_m: f64,
        speed_m_s: f64,
    },
    /// Linear kinematic trajectory with constant velocity.
    Linear {
        start_offset_ned: Vector3<f64>,
        velocity_ned: Vector3<f64>,
    },
}

/// Simulation scenario configuration.
#[derive(Debug, Clone)]
pub struct SimulationConfig {
    pub base_ecef: Vector3<f64>,
    pub start_time: GpsTime,
    pub duration_s: f64,
    pub epoch_rate_hz: f64,
    pub profile: TrajectoryProfile,
    pub pr_noise_m: f64,
    pub cp_noise_m: f64,
    pub doppler_noise_m_s: f64,
    pub num_satellites: usize,
    /// Satellite outage interval: (start_s, end_s, dropped_sat_ids)
    pub outages: Vec<(f64, f64, Vec<u8>)>,
    /// Injected cycle slips: (time_s, sat_prn, slip_cycles)
    pub cycle_slips: Vec<(f64, u8, i32)>,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            base_ecef: Vector3::new(-3961904.4341, 3348994.2660, 3698211.7067),
            start_time: GpsTime::new(2200, 300000.0),
            duration_s: 60.0,
            epoch_rate_hz: 1.0,
            profile: TrajectoryProfile::Circular {
                center_offset_ned: Vector3::new(100.0, 100.0, 0.0),
                radius_m: 50.0,
                speed_m_s: 5.0,
            },
            pr_noise_m: 0.20,
            cp_noise_m: 0.002,
            doppler_noise_m_s: 0.02,
            num_satellites: 24,
            outages: Vec::new(),
            cycle_slips: Vec::new(),
        }
    }
}

/// Result of running a simulation scenario.
pub struct SimulationDataset {
    pub ephemerides: Vec<Ephemeris>,
    pub base_epochs: Vec<EpochObs>,
    pub rover_epochs: Vec<EpochObs>,
    pub truth_positions: Vec<(GpsTime, Vector3<f64>)>,
    pub truth_velocities: Vec<(GpsTime, Vector3<f64>)>,
}

/// Deterministic pseudo-random number generator for reproducible noise.
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self { state: seed.max(1) }
    }

    fn next_f64(&mut self) -> f64 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let x = ((self.state >> 33) ^ self.state).wrapping_mul(0xff51afd7ed558ccd);
        let y = ((x >> 33) ^ x) as u32;
        (y as f64) / (u32::MAX as f64)
    }

    fn next_gaussian(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-12);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }
}

/// Compute rover position and velocity at time t in ECEF.
fn compute_rover_kinematics(
    cfg: &SimulationConfig,
    t_s: f64,
    ned_to_ecef: &nalgebra::Matrix3<f64>,
) -> (Vector3<f64>, Vector3<f64>) {
    let (pos_ned, vel_ned) = match &cfg.profile {
        TrajectoryProfile::Static { offset_ned } => (*offset_ned, Vector3::zeros()),
        TrajectoryProfile::Linear { start_offset_ned, velocity_ned } => {
            (*start_offset_ned + *velocity_ned * t_s, *velocity_ned)
        }
        TrajectoryProfile::Circular { center_offset_ned, radius_m, speed_m_s } => {
            let omega = speed_m_s / radius_m.max(1.0);
            let angle = omega * t_s;
            let p = *center_offset_ned + Vector3::new(radius_m * angle.cos(), radius_m * angle.sin(), 0.0);
            let v = Vector3::new(-speed_m_s * angle.sin(), speed_m_s * angle.cos(), 0.0);
            (p, v)
        }
    };
    (base_ecef_offset(cfg.base_ecef, ned_to_ecef, pos_ned), ned_to_ecef * vel_ned)
}

fn base_ecef_offset(base: Vector3<f64>, ned_to_ecef: &nalgebra::Matrix3<f64>, ned: Vector3<f64>) -> Vector3<f64> {
    base + ned_to_ecef * ned
}

/// Generate synthetic ephemerides for GPS satellites distributed in a 6-plane constellation.
pub fn generate_synthetic_ephemerides(num_sats: usize, start_time: GpsTime) -> Vec<Ephemeris> {
    let mut ephems = Vec::with_capacity(num_sats);
    let r_orbit: f64 = 26560000.0;
    let n_planes = 6;
    let sats_per_plane = num_sats.div_ceil(n_planes);

    for i in 0..num_sats {
        let plane_idx = i % n_planes;
        let sat_in_plane = i / n_planes;
        let prn = (i + 1) as u8;
        let raan0 = 2.0 * PI * (plane_idx as f64) / (n_planes as f64);
        let m0 = 2.0 * PI * (sat_in_plane as f64) / (sats_per_plane.max(1) as f64) + (plane_idx as f64) * 0.6;
        let eph = GpsEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn,
            },
            toe: start_time,
            toc: start_time,
            af0: 0.0,
            af1: 0.0,
            af2: 0.0,
            crs: 0.0,
            delta_n: 0.0,
            m0,
            cuc: 0.0,
            e: 0.001,
            cus: 0.0,
            sqrt_a: r_orbit.sqrt(),
            cic: 0.0,
            omega0: raan0,
            cis: 0.0,
            i0: 55.0_f64.to_radians(),
            crc: 0.0,
            omega: 0.0,
            omega_dot: -2.6e-9,
            idot: 0.0,
            tgd: 0.0,
            iodc: 1,
            iode: 1,
        };
        ephems.push(Ephemeris::Gps(eph));
    }
    ephems
}

/// Run simulation and return full synthetic dataset.
pub fn generate_simulation_dataset(cfg: &SimulationConfig) -> SimulationDataset {
    let ephemerides = generate_synthetic_ephemerides(cfg.num_satellites, cfg.start_time);
    let base_llh = gneiss_core::coords::ecef_to_llh(cfg.base_ecef);
    let ned_to_ecef = gneiss_core::coords::ecef_to_ned_matrix(base_llh).transpose();

    let mut rng = SimpleRng::new(42);
    let total_steps = (cfg.duration_s * cfg.epoch_rate_hz).round() as usize;
    let dt = 1.0 / cfg.epoch_rate_hz;

    let mut base_epochs = Vec::with_capacity(total_steps);
    let mut rover_epochs = Vec::with_capacity(total_steps);
    let mut truth_pos = Vec::with_capacity(total_steps);
    let mut truth_vel = Vec::with_capacity(total_steps);

    // Fixed true ambiguities for L1 (cycles) per satellite
    let mut true_amb_l1: Vec<i32> = (0..cfg.num_satellites)
        .map(|i| 100_000 + (i as i32) * 5432)
        .collect();

    for step in 0..total_steps {
        let t_s = step as f64 * dt;
        let ep_time = GpsTime::new(cfg.start_time.week, cfg.start_time.tow + t_s);
        let (r_pos, r_vel) = compute_rover_kinematics(cfg, t_s, &ned_to_ecef);

        truth_pos.push((ep_time, r_pos));
        truth_vel.push((ep_time, r_vel));

        apply_cycle_slips(t_s, &cfg.cycle_slips, &mut true_amb_l1);
        let dropped = get_dropped_sats(t_s, &cfg.outages);

        let (base_ep, rover_ep) = generate_epoch_pair(
            ep_time, cfg.base_ecef, r_pos, r_vel, &ephemerides,
            &true_amb_l1, &dropped, cfg, &mut rng,
        );

        base_epochs.push(base_ep);
        rover_epochs.push(rover_ep);
    }

    SimulationDataset {
        ephemerides,
        base_epochs,
        rover_epochs,
        truth_positions: truth_pos,
        truth_velocities: truth_vel,
    }
}

fn apply_cycle_slips(t_s: f64, slips: &[(f64, u8, i32)], ambs: &mut [i32]) {
    for (slip_t, prn, delta) in slips {
        if (t_s - slip_t).abs() < 0.01 {
            let idx = (*prn as usize).saturating_sub(1);
            if idx < ambs.len() {
                ambs[idx] += delta;
            }
        }
    }
}

fn get_dropped_sats(t_s: f64, outages: &[(f64, f64, Vec<u8>)]) -> Vec<u8> {
    for (start, end, sats) in outages {
        if t_s >= *start && t_s <= *end {
            return sats.clone();
        }
    }
    Vec::new()
}

fn compute_satellite_range(eph: &Ephemeris, time: GpsTime, rx_pos: Vector3<f64>) -> f64 {
    let (sat_p_rough, _, _, _) = eph.position(time);
    let tau = (sat_p_rough - rx_pos).norm() / SPEED_OF_LIGHT_M_S;
    let (sat_p, _, _, _) = eph.position(GpsTime::new(time.week, time.tow - tau));
    let om = gneiss_core::constants::EARTH_ROTATION_RATE_RAD_S * tau;
    let sat_p_rot = Vector3::new(
        sat_p.x * om.cos() + sat_p.y * om.sin(),
        -sat_p.x * om.sin() + sat_p.y * om.cos(),
        sat_p.z,
    );
    (sat_p_rot - rx_pos).norm()
}

#[allow(clippy::too_many_arguments)]
fn generate_epoch_pair(
    time: GpsTime,
    base_pos: Vector3<f64>,
    rover_pos: Vector3<f64>,
    rover_vel: Vector3<f64>,
    ephems: &[Ephemeris],
    ambs: &[i32],
    dropped: &[u8],
    cfg: &SimulationConfig,
    rng: &mut SimpleRng,
) -> (EpochObs, EpochObs) {
    let mut base_sats = Vec::new();
    let mut rover_sats = Vec::new();

    let lambda_l1 = SPEED_OF_LIGHT_M_S / 1575.42e6;
    let lambda_l2 = SPEED_OF_LIGHT_M_S / 1227.60e6;

    for (i, eph) in ephems.iter().enumerate() {
        let prn = eph.sat().prn;
        if dropped.contains(&prn) {
            continue;
        }

        let r_base = compute_satellite_range(eph, time, base_pos);
        let r_rov = compute_satellite_range(eph, time, rover_pos);

        let dt = 0.005;
        let t_p = GpsTime::new(time.week, time.tow - dt);
        let t_n = GpsTime::new(time.week, time.tow + dt);
        let rov_p_prev = rover_pos - rover_vel * dt;
        let rov_p_next = rover_pos + rover_vel * dt;

        let r_dot_base = (compute_satellite_range(eph, t_n, base_pos) - compute_satellite_range(eph, t_p, base_pos)) / (2.0 * dt);
        let r_dot_rov = (compute_satellite_range(eph, t_n, rov_p_next) - compute_satellite_range(eph, t_p, rov_p_prev)) / (2.0 * dt);

        let dop_base = -r_dot_base / lambda_l1;
        let dop_rov = -r_dot_rov / lambda_l1;

        let amb1 = ambs[i] as f64;
        let amb2 = (80_000 + (i as i32) * 4231) as f64;

        base_sats.push(build_sat_obs(eph.sat(), r_base, dop_base, 0.0, 0.0, lambda_l1, lambda_l2, cfg, rng));
        rover_sats.push(build_sat_obs(eph.sat(), r_rov, dop_rov, amb1, amb2, lambda_l1, lambda_l2, cfg, rng));
    }

    (
        EpochObs { time, satellites: base_sats },
        EpochObs { time, satellites: rover_sats },
    )
}

#[allow(clippy::too_many_arguments)]
fn build_sat_obs(
    sat: SatelliteId,
    range: f64,
    dop: f64,
    amb_l1: f64,
    amb_l2: f64,
    lambda1: f64,
    lambda2: f64,
    cfg: &SimulationConfig,
    rng: &mut SimpleRng,
) -> SatObs {
    let pr_err = rng.next_gaussian() * cfg.pr_noise_m;
    let cp_err = rng.next_gaussian() * cfg.cp_noise_m;
    let dop_err = rng.next_gaussian() * cfg.doppler_noise_m_s;

    let c1c = ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 1, attribute: 'C' } };
    let l1c = ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band: 1, attribute: 'C' } };
    let d1c = ObsCode { obs_type: ObsType::Doppler, signal: SignalCode { freq_band: 1, attribute: 'C' } };
    let s1c = ObsCode { obs_type: ObsType::Snr, signal: SignalCode { freq_band: 1, attribute: 'C' } };

    let c2w = ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 2, attribute: 'W' } };
    let l2w = ObsCode { obs_type: ObsType::CarrierPhase, signal: SignalCode { freq_band: 2, attribute: 'W' } };

    let obs = vec![
        Observation { code: c1c, value: range + pr_err, lock_time: Some(100), lli: None },
        Observation { code: l1c, value: (range + cp_err) / lambda1 + amb_l1, lock_time: Some(100), lli: None },
        Observation { code: d1c, value: dop + dop_err / lambda1, lock_time: None, lli: None },
        Observation { code: s1c, value: 45.0, lock_time: None, lli: None },
        Observation { code: c2w, value: range + pr_err * 1.2, lock_time: Some(100), lli: None },
        Observation { code: l2w, value: (range + cp_err * 1.2) / lambda2 + amb_l2, lock_time: Some(100), lli: None },
    ];

    SatObs { sat, observations: obs }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simulation_generator_produces_consistent_ranges() {
        let cfg = SimulationConfig {
            duration_s: 5.0,
            ..Default::default()
        };
        let dataset = generate_simulation_dataset(&cfg);
        assert_eq!(dataset.rover_epochs.len(), 5);
        assert_eq!(dataset.base_epochs.len(), 5);
        assert_eq!(dataset.truth_positions.len(), 5);

        let first_base = &dataset.base_epochs[0];
        assert_eq!(first_base.satellites.len(), 24);
        let pr = first_base.satellites[0].get_observable(1);
        assert!(pr.is_some());
        assert!(pr.unwrap() > 20_000_000.0);
    }
}
