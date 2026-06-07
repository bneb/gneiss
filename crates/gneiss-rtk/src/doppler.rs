/// Doppler-based velocity estimation for EKF update.
///
/// Doppler measurements provide instantaneous velocity information by measuring
/// the rate of change of the carrier phase. The observation model is:
///   ḋ = -ė · (v_rcv - v_sat) + ċdt_drift
///
/// This module computes the Doppler velocity EKF update matrices (z, H, R)
/// for use in the SPP and RTK engines.

use nalgebra::{DMatrix, DVector, Vector3};
use gneiss_core::ephemeris::Ephemeris;
use gneiss_core::time::GpsTime;
use crate::filter::RtkState;

const LIGHT_SPEED: f64 = 299792458.0;

/// Represents a Doppler measurement for one satellite.
#[derive(Debug, Clone)]
pub struct DopplerMeasurement {
    pub sat: gneiss_core::sat::SatelliteId,
    /// Doppler shift in Hz (positive = satellite approaching)
    pub doppler_hz: f64,
    /// L1 frequency in Hz
    pub frequency: f64,
    /// Elevation angle in radians
    pub elevation: f64,
    /// SNR in dB-Hz
    pub snr: f64,
}

/// Computes the Doppler velocity EKF update matrices.
///
/// Returns (z, H, R) where:
/// - z: innovation vector (observed - predicted Doppler range-rate)
/// - H: observation matrix mapping state to Doppler measurements
/// - R: measurement noise covariance
///
/// Returns None if fewer than 4 measurements are available.
pub fn compute_doppler_update(
    state: &RtkState,
    measurements: &[DopplerMeasurement],
    ephemerides: &[Ephemeris],
) -> Option<(DVector<f64>, DMatrix<f64>, DMatrix<f64>)> {
    if measurements.len() < 4 {
        return None;
    }

    let state_size = state.covariance.ncols();
    let mut z_vals = Vec::new();
    let mut h_rows = Vec::new();
    let mut r_vals = Vec::new();

    for meas in measurements {
        let eph = ephemerides.iter().find(|e| e.sat() == meas.sat)?;

        // Get satellite position and velocity
        let (sat_pos, sat_vel, _dt_s, _) = eph.position(state.time);

        // Line-of-sight unit vector from receiver to satellite
        let los = sat_pos - state.position.vector;
        let range = los.norm();
        if range < 1.0 {
            continue;
        }
        let e_los = los / range;

        // Observed range-rate from Doppler (negative sign: positive Doppler = approaching = negative range-rate)
        let lambda = LIGHT_SPEED / meas.frequency;
        let observed_range_rate = -meas.doppler_hz * lambda;

        // Predicted range-rate: ė · (v_sat - v_rcv)
        // Note: This ignores Earth rotation correction for simplicity
        let rel_vel = sat_vel - state.velocity;
        let predicted_range_rate = e_los.dot(&rel_vel);

        // Innovation
        let z_i = observed_range_rate - predicted_range_rate;

        // H matrix row: ∂(range_rate)/∂(state)
        // Range-rate = ė · (v_sat - v_rcv), so ∂(range_rate)/∂v_rcv = -ė
        let mut h_row = vec![0.0; state_size];
        // Velocity is at state indices 3, 4, 5
        h_row[3] = -e_los.x;
        h_row[4] = -e_los.y;
        h_row[5] = -e_los.z;

        // Variance: scaled by elevation and SNR
        let var = 0.04 * gneiss_core::variance::observation_variance(meas.snr, meas.elevation, 45.0);

        z_vals.push(z_i);
        h_rows.push(h_row);
        r_vals.push(var);
    }

    if z_vals.len() < 4 {
        return None;
    }

    let n_meas = z_vals.len();
    let z = DVector::from_vec(z_vals);
    let mut h = DMatrix::zeros(n_meas, state_size);
    for (i, row) in h_rows.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            h[(i, j)] = val;
        }
    }
    let mut r = DMatrix::zeros(n_meas, n_meas);
    for (i, &val) in r_vals.iter().enumerate() {
        r[(i, i)] = val;
    }

    Some((z, h, r))
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::sat::{SatelliteId, Constellation};
    use gneiss_core::ephemeris::{Ephemeris, GpsEphemeris};
    use gneiss_core::time::GpsTime;
    use nalgebra::Vector3;

    fn make_test_ephemeris(prn: u8, m0: f64, omega0: f64) -> Ephemeris {
        let time = GpsTime::new(2137, 422922.0);
        Ephemeris::Gps(GpsEphemeris {
            sat: SatelliteId { constellation: Constellation::Gps, prn },
            toe: time, toc: time,
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
            omega0, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
            iode: 0, iodc: 0,
        })
    }

    #[test]
    fn test_doppler_velocity_stationary() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(6378137.0, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, time,
        );
        let state = RtkState::new(time, pos, 1.0);

        let ephemerides = vec![
            make_test_ephemeris(1, 0.0, 0.0),
            make_test_ephemeris(2, 1.0, 0.5),
            make_test_ephemeris(3, 2.0, 1.0),
            make_test_ephemeris(4, 3.0, 1.5),
            make_test_ephemeris(5, 0.5, 2.0),
        ];

        // For a stationary receiver, Doppler should reflect only satellite velocity.
        // The predicted range-rate = ė · (v_sat - 0) and observed should match if
        // we construct Doppler from the satellite velocity.
        let mut measurements = Vec::new();
        for eph in &ephemerides {
            let (sat_pos, sat_vel, _, _) = eph.position(time);
            let los = sat_pos - state.position.vector;
            let range = los.norm();
            let e_los = los / range;
            let range_rate = e_los.dot(&sat_vel);
            let freq = 1575.42e6;
            let lambda = LIGHT_SPEED / freq;
            let doppler_hz = -range_rate / lambda;

            let llh = gneiss_core::coords::ecef_to_llh(state.position.vector);
            let (_, el) = gneiss_core::coords::az_el(llh, state.position.vector, sat_pos);

            measurements.push(DopplerMeasurement {
                sat: eph.sat(),
                doppler_hz,
                frequency: freq,
                elevation: el,
                snr: 45.0,
            });
        }

        let result = compute_doppler_update(&state, &measurements, &ephemerides);
        assert!(result.is_some(), "Should compute Doppler update with 5 sats");

        let (z, _h, _r) = result.unwrap();
        // For a perfectly consistent stationary scenario, innovations should be near zero
        for i in 0..z.len() {
            assert!(z[i].abs() < 0.1,
                "Doppler innovation[{}] should be near zero for stationary, got {}", i, z[i]);
        }
    }

    #[test]
    fn test_doppler_velocity_known_motion() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(6378137.0, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        // State thinks receiver is stationary
        state.velocity = Vector3::zeros();

        let true_velocity = Vector3::new(10.0, 0.0, 0.0); // 10 m/s in ECEF X

        let ephemerides = vec![
            make_test_ephemeris(1, 0.0, 0.0),
            make_test_ephemeris(2, 1.0, 0.5),
            make_test_ephemeris(3, 2.0, 1.0),
            make_test_ephemeris(4, 3.0, 1.5),
            make_test_ephemeris(5, 0.5, 2.0),
        ];

        // Construct Doppler from true velocity
        let mut measurements = Vec::new();
        for eph in &ephemerides {
            let (sat_pos, sat_vel, _, _) = eph.position(time);
            let los = sat_pos - state.position.vector;
            let range = los.norm();
            let e_los = los / range;
            let rel_vel = sat_vel - true_velocity;
            let range_rate = e_los.dot(&rel_vel);
            let freq = 1575.42e6;
            let lambda = LIGHT_SPEED / freq;
            let doppler_hz = -range_rate / lambda;

            let llh = gneiss_core::coords::ecef_to_llh(state.position.vector);
            let (_, el) = gneiss_core::coords::az_el(llh, state.position.vector, sat_pos);

            measurements.push(DopplerMeasurement {
                sat: eph.sat(),
                doppler_hz,
                frequency: freq,
                elevation: el,
                snr: 45.0,
            });
        }

        let result = compute_doppler_update(&state, &measurements, &ephemerides);
        assert!(result.is_some());

        let (z, h, r) = result.unwrap();
        // Apply the Kalman update manually to verify velocity correction
        let h_t = h.transpose();
        let s = &h * &state.covariance * &h_t + &r;
        if let Some(s_inv) = s.try_inverse() {
            let k = &state.covariance * &h_t * s_inv;
            let dx = &k * &z;

            // The velocity correction should point toward the true velocity
            let dv = Vector3::new(dx[3], dx[4], dx[5]);
            // After correction: state.velocity + dv should be closer to true_velocity
            let corrected_vel = state.velocity + dv;
            let error_before = (state.velocity - true_velocity).norm();
            let error_after = (corrected_vel - true_velocity).norm();
            assert!(error_after < error_before,
                "Doppler should improve velocity estimate: before={:.3}, after={:.3}",
                error_before, error_after);
        }
    }

    #[test]
    fn test_doppler_insufficient_sats() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(6378137.0, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, time,
        );
        let state = RtkState::new(time, pos, 1.0);

        // Only 2 measurements — should return None
        let measurements = vec![
            DopplerMeasurement {
                sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
                doppler_hz: 1000.0,
                frequency: 1575.42e6,
                elevation: 0.5,
                snr: 45.0,
            },
            DopplerMeasurement {
                sat: SatelliteId { constellation: Constellation::Gps, prn: 2 },
                doppler_hz: -500.0,
                frequency: 1575.42e6,
                elevation: 0.8,
                snr: 40.0,
            },
        ];

        let ephemerides = vec![
            make_test_ephemeris(1, 0.0, 0.0),
            make_test_ephemeris(2, 1.0, 0.5),
        ];

        let result = compute_doppler_update(&state, &measurements, &ephemerides);
        assert!(result.is_none(), "Should return None with < 4 sats");
    }

    #[test]
    fn test_doppler_h_matrix_structure() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(6378137.0, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, time,
        );
        let state = RtkState::new(time, pos, 1.0);

        let ephemerides = vec![
            make_test_ephemeris(1, 0.0, 0.0),
            make_test_ephemeris(2, 1.0, 0.5),
            make_test_ephemeris(3, 2.0, 1.0),
            make_test_ephemeris(4, 3.0, 1.5),
        ];

        let mut measurements = Vec::new();
        for eph in &ephemerides {
            let (sat_pos, _, _, _) = eph.position(time);
            let llh = gneiss_core::coords::ecef_to_llh(state.position.vector);
            let (_, el) = gneiss_core::coords::az_el(llh, state.position.vector, sat_pos);
            measurements.push(DopplerMeasurement {
                sat: eph.sat(),
                doppler_hz: 0.0,
                frequency: 1575.42e6,
                elevation: el,
                snr: 45.0,
            });
        }

        let (_, h, _) = compute_doppler_update(&state, &measurements, &ephemerides).unwrap();

        // H matrix should only have non-zero entries in velocity columns (3, 4, 5)
        for i in 0..h.nrows() {
            for j in 0..h.ncols() {
                if (3..=5).contains(&j) {
                    // Velocity columns can be non-zero
                } else {
                    assert!(h[(i, j)].abs() < 1e-12,
                        "H[{},{}] should be zero for Doppler velocity, got {}", i, j, h[(i, j)]);
                }
            }
            // Each row should be a unit vector in velocity space
            let row_norm = (h[(i, 3)].powi(2) + h[(i, 4)].powi(2) + h[(i, 5)].powi(2)).sqrt();
            assert!((row_norm - 1.0).abs() < 1e-6,
                "H row {} velocity components should form unit vector, norm={}", i, row_norm);
        }
    }
}
