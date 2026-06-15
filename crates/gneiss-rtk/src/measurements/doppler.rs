/// Doppler-based velocity estimation for EKF update.
///
/// Doppler measurements provide instantaneous velocity information by measuring
/// the rate of change of the carrier phase. The observation model is:
///   ḋ = -ė · (v_rcv - v_sat) + ċdt_drift
///
/// This module computes the Doppler velocity EKF update matrices (z, H, R)
/// for use in the SPP and RTK engines.

use nalgebra::{DMatrix, DVector};
use gneiss_core::ephemeris::Ephemeris;
use crate::filter::RtkState;

const LIGHT_SPEED: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

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
    dop_base_var: f64,
) -> Option<(DVector<f64>, DMatrix<f64>, DMatrix<f64>)> {

    let state_size = state.covariance.ncols();
    let mut z_vals = Vec::new();
    let mut h_rows = Vec::new();
    let mut r_vals = Vec::new();

    for meas in measurements {
        let eph = match ephemerides.iter().find(|e| e.sat() == meas.sat) {
            Some(e) => e,
            None => continue,
        };

        // Get satellite position and velocity
        let (sat_pos, sat_vel, _dt_s, sat_drift_sec_per_sec) = eph.position(state.time);
        let sat_drift_ms = sat_drift_sec_per_sec * LIGHT_SPEED;

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

        // Predicted range-rate: ė · (v_sat - v_rcv) + c_drift - sat_drift_ms
        // Note: This ignores Earth rotation correction for simplicity
        let rel_vel = sat_vel - state.velocity;
        let predicted_range_rate = e_los.dot(&rel_vel) + state.rcv_clk_drift - sat_drift_ms;

        // Innovation
        let z_i = observed_range_rate - predicted_range_rate;

        // H matrix row: ∂(range_rate)/∂(state)
        // Range-rate = ė · (v_sat - v_rcv) + c_drift, so ∂(range_rate)/∂v_rcv = -ė
        let mut h_row = vec![0.0; state_size];
        // Velocity is at state indices 3, 4, 5
        h_row[3] = -e_los.x;
        h_row[4] = -e_los.y;
        h_row[5] = -e_los.z;
        if state_size > 19 {
            h_row[19] = 1.0; // Receiver clock drift
        }

        // Variance: scaled by elevation and SNR
        let var = dop_base_var * gneiss_core::variance::observation_variance(meas.snr, meas.elevation, 45.0);

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
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
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

        let result = compute_doppler_update(&state, &measurements, &ephemerides, 1.0);
        assert!(result.is_some(), "Should compute Doppler update with 5 sats");
        let (z, h, r) = result.unwrap();
        assert_eq!(z.len(), measurements.len(), "Innovation vector size should match measurements");
        assert_eq!(r.nrows(), measurements.len(), "R matrix size should match measurements");
        assert_eq!(h.nrows(), measurements.len(), "H matrix rows should match measurements");
        assert_eq!(h.ncols(), state.covariance.ncols(), "H matrix cols should match state size");
        // For a consistent stationary scenario, innovations should be near zero
        for i in 0..z.len() {
            assert!(z[i].abs() < 0.1,
                "Doppler innovation[{}] should be near zero for stationary, got {}", i, z[i]);
            let expected_var = 0.04 * gneiss_core::variance::observation_variance(45.0, measurements[i].elevation, 45.0);
            assert!((r[(i, i)] - expected_var).abs() < 1e-6, "r matrix should have correct variance");
        }
    }

    #[test]
    fn test_doppler_velocity_known_motion() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        // State thinks receiver is moving! This catches the relative velocity mutant.
        let receiver_velocity = Vector3::new(5.0, -2.0, 1.0);
        state.velocity = receiver_velocity;

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

        let result = compute_doppler_update(&state, &measurements, &ephemerides, 1.0);
        assert!(result.is_some());

        let (z, h, r) = result.unwrap();
        assert_eq!(z.len(), measurements.len(), "Innovation vector size should match measurements");
        assert_eq!(r.nrows(), measurements.len(), "R matrix size should match measurements");
        assert_eq!(h.nrows(), measurements.len(), "H matrix rows should match measurements");
        assert_eq!(h.ncols(), state.covariance.ncols(), "H matrix cols should match state size");

        // Apply the Kalman update manually to verify velocity correction
        let h_t = h.transpose();
        let s = &h * &state.covariance * &h_t + &r;
        let s_inv = s.try_inverse().expect("S matrix should be invertible");
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
        assert!(error_after < 0.1, "Velocity error should be very small, got {}", error_after);
    }

    #[test]
    fn test_doppler_insufficient_sats() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
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

        let result = compute_doppler_update(&state, &measurements, &ephemerides, 1.0);
        assert!(result.is_none(), "Should return None with < 4 sats");
    }

    #[test]
    fn test_doppler_h_matrix_structure() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
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

        let (_, h, _) = compute_doppler_update(&state, &measurements, &ephemerides, 1.0).unwrap();

        for i in 0..h.nrows() {
            let eph = &ephemerides[i];
            let (sat_pos, _, _, _) = eph.position(time);
            let los = sat_pos - state.position.vector;
            let range = los.norm();
            let e_los = los / range;

            for j in 0..h.ncols() {
                if j == 3 {
                    assert!((h[(i, j)] - (-e_los.x)).abs() < 1e-6, "H[{},3] should be -e_los.x", i);
                } else if j == 4 {
                    assert!((h[(i, j)] - (-e_los.y)).abs() < 1e-6, "H[{},4] should be -e_los.y", i);
                } else if j == 5 {
                    assert!((h[(i, j)] - (-e_los.z)).abs() < 1e-6, "H[{},5] should be -e_los.z", i);
                } else {
                    assert!(h[(i, j)].abs() < 1e-12,
                        "H[{},{}] should be zero for Doppler velocity, got {}", i, j, h[(i, j)]);
                }
            }
        }
    }

    #[test]
    fn test_doppler_exact_range_mutant() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::zeros();

        let ephemerides = vec![
            make_test_ephemeris(1, 0.0, 0.0),
            make_test_ephemeris(2, 1.0, 0.5),
            make_test_ephemeris(3, 2.0, 1.0),
            make_test_ephemeris(4, 3.0, 1.5),
        ];

        let (sat_pos, _, _, _) = ephemerides[0].position(time);
        state.position.vector = sat_pos + Vector3::new(1.0, 0.0, 0.0); // range = exactly 1.0!

        let mut measurements = Vec::new();
        for prn in 1..=4 {
            let sat = SatelliteId { constellation: Constellation::Gps, prn };
            measurements.push(DopplerMeasurement {
                sat,
                doppler_hz: 0.0,
                frequency: 1575.42e6,
                elevation: 0.5,
                snr: 45.0,
            });
        }

        let result = compute_doppler_update(&state, &measurements, &ephemerides, 1.0);
        // Since range is exactly 1.0, it is NOT < 1.0, so it should NOT be skipped!
        // Therefore, we have 4 valid measurements, so it should succeed.
        assert!(result.is_some(), "Should NOT skip if range == 1.0");
    }

    #[test]
    fn test_doppler_short_range_continue_mutant() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::zeros();

        let ephemerides = vec![
            make_test_ephemeris(1, 0.0, 0.0),
            make_test_ephemeris(2, 1.0, 0.5),
            make_test_ephemeris(3, 2.0, 1.0),
            make_test_ephemeris(4, 3.0, 1.5),
            make_test_ephemeris(5, 4.0, 2.0),
        ];

        // Make the FIRST satellite short range
        let (sat_pos, _, _, _) = ephemerides[0].position(time);
        state.position.vector = sat_pos + Vector3::new(0.5, 0.0, 0.0);

        let mut measurements = Vec::new();
        for prn in 1..=5 {
            let sat = SatelliteId { constellation: Constellation::Gps, prn };
            measurements.push(DopplerMeasurement {
                sat,
                doppler_hz: 0.0,
                frequency: 1575.42e6,
                elevation: 0.5,
                snr: 45.0,
            });
        }

        let result = compute_doppler_update(&state, &measurements, &ephemerides, 1.0);
        // It should SKIP the 1st one (continue), but process the remaining 4.
        // So it should return Some with 4 measurements!
        assert!(result.is_some(), "Should not break; should process remaining 4 sats");
        assert_eq!(result.unwrap().0.len(), 4, "Should have exactly 4 measurements processed");
    }

    #[test]
    fn test_missing_ephemeris_first_mutant() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::zeros();

        let ephemerides = vec![
            make_test_ephemeris(2, 1.0, 0.5),
            make_test_ephemeris(3, 2.0, 1.0),
            make_test_ephemeris(4, 3.0, 1.5),
            make_test_ephemeris(5, 4.0, 2.0),
        ];

        let mut measurements = Vec::new();
        for prn in 1..=5 {
            let sat = SatelliteId { constellation: Constellation::Gps, prn };
            measurements.push(DopplerMeasurement {
                sat,
                doppler_hz: 0.0,
                frequency: 1575.42e6,
                elevation: 0.5,
                snr: 45.0,
            });
        }

        let result = compute_doppler_update(&state, &measurements, &ephemerides, 1.0);
        // PRN 1 is missing, but 2,3,4,5 are present. It should continue on 1.
        assert!(result.is_some(), "Should not break on first missing ephemeris");
        assert_eq!(result.unwrap().0.len(), 4, "Should process the 4 valid ones");
    }
}

#[cfg(test)]
mod missing_eph_tests {
    use super::*;
    use nalgebra::Vector3;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::sat::{SatelliteId, Constellation};
    use gneiss_core::ephemeris::{Ephemeris, GpsEphemeris};
    use gneiss_core::time::GpsTime;

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
    fn test_missing_ephemeris() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::zeros();

        let ephemerides = vec![
            make_test_ephemeris(1, 0.0, 0.0),
            make_test_ephemeris(2, 1.0, 0.5),
            make_test_ephemeris(3, 2.0, 1.0),
            make_test_ephemeris(4, 3.0, 1.5),
        ];

        let mut measurements = Vec::new();
        for prn in 1..=5 {
            let sat = SatelliteId { constellation: Constellation::Gps, prn };
            measurements.push(DopplerMeasurement {
                sat,
                doppler_hz: 0.0,
                frequency: 1575.42e6,
                elevation: 0.5,
                snr: 45.0,
            });
        }

        let result = compute_doppler_update(&state, &measurements, &ephemerides, 1.0);
        assert!(result.is_some(), "Should succeed with 4 valid ephemerides out of 5 measurements");
        assert_eq!(result.unwrap().0.len(), 4, "Should process exactly 4 valid measurements");
    }

    #[test]
    fn test_doppler_short_range() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::zeros();

        let ephemerides = vec![
            make_test_ephemeris(1, 0.0, 0.0),
            make_test_ephemeris(2, 1.0, 0.5),
            make_test_ephemeris(3, 2.0, 1.0),
            make_test_ephemeris(4, 3.0, 1.5),
        ];

        let (sat_pos, _, _, _) = ephemerides[0].position(time);
        state.position.vector = sat_pos + Vector3::new(0.5, 0.0, 0.0); // range = 0.5 < 1.0

        let mut measurements = Vec::new();
        for prn in 1..=4 {
            let sat = SatelliteId { constellation: Constellation::Gps, prn };
            measurements.push(DopplerMeasurement {
                sat,
                doppler_hz: 0.0,
                frequency: 1575.42e6,
                elevation: 0.5,
                snr: 45.0,
            });
        }

        let result = compute_doppler_update(&state, &measurements, &ephemerides, 1.0);
        assert!(result.is_none(), "Should fail because short range satellite is skipped");
    }

    #[test]
    fn test_doppler_exact_range_mutant() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::zeros();

        let ephemerides = vec![
            make_test_ephemeris(1, 0.0, 0.0),
            make_test_ephemeris(2, 1.0, 0.5),
            make_test_ephemeris(3, 2.0, 1.0),
            make_test_ephemeris(4, 3.0, 1.5),
        ];

        let (sat_pos, _, _, _) = ephemerides[0].position(time);
        state.position.vector = sat_pos + Vector3::new(1.0, 0.0, 0.0); // range = exactly 1.0!

        let mut measurements = Vec::new();
        for prn in 1..=4 {
            let sat = SatelliteId { constellation: Constellation::Gps, prn };
            measurements.push(DopplerMeasurement {
                sat,
                doppler_hz: 0.0,
                frequency: 1575.42e6,
                elevation: 0.5,
                snr: 45.0,
            });
        }

        let result = compute_doppler_update(&state, &measurements, &ephemerides, 1.0);
        // Since range is exactly 1.0, it is NOT < 1.0, so it should NOT be skipped!
        // Therefore, we have 4 valid measurements, so it should succeed.
        assert!(result.is_some(), "Should NOT skip if range == 1.0");
    }

    #[test]
    fn test_doppler_short_range_continue_mutant() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::zeros();

        let ephemerides = vec![
            make_test_ephemeris(1, 0.0, 0.0),
            make_test_ephemeris(2, 1.0, 0.5),
            make_test_ephemeris(3, 2.0, 1.0),
            make_test_ephemeris(4, 3.0, 1.5),
            make_test_ephemeris(5, 4.0, 2.0),
        ];

        // Make the FIRST satellite short range
        let (sat_pos, _, _, _) = ephemerides[0].position(time);
        state.position.vector = sat_pos + Vector3::new(0.5, 0.0, 0.0);

        let mut measurements = Vec::new();
        for prn in 1..=5 {
            let sat = SatelliteId { constellation: Constellation::Gps, prn };
            measurements.push(DopplerMeasurement {
                sat,
                doppler_hz: 0.0,
                frequency: 1575.42e6,
                elevation: 0.5,
                snr: 45.0,
            });
        }

        let result = compute_doppler_update(&state, &measurements, &ephemerides, 1.0);
        // It should SKIP the 1st one (continue), but process the remaining 4.
        // So it should return Some with 4 measurements!
        assert!(result.is_some(), "Should not break; should process remaining 4 sats");
        assert_eq!(result.unwrap().0.len(), 4, "Should have exactly 4 measurements processed");
    }

    #[test]
    fn test_missing_ephemeris_first_mutant() {
        let time = GpsTime::new(2137, 422922.0);
        let pos = Coordinate::new(
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
            Datum::WGS84, Frame::ECEF, time,
        );
        let mut state = RtkState::new(time, pos, 1.0);
        state.velocity = Vector3::zeros();

        let ephemerides = vec![
            make_test_ephemeris(2, 1.0, 0.5),
            make_test_ephemeris(3, 2.0, 1.0),
            make_test_ephemeris(4, 3.0, 1.5),
            make_test_ephemeris(5, 4.0, 2.0),
        ];

        let mut measurements = Vec::new();
        for prn in 1..=5 {
            let sat = SatelliteId { constellation: Constellation::Gps, prn };
            measurements.push(DopplerMeasurement {
                sat,
                doppler_hz: 0.0,
                frequency: 1575.42e6,
                elevation: 0.5,
                snr: 45.0,
            });
        }

        let result = compute_doppler_update(&state, &measurements, &ephemerides, 1.0);
        // PRN 1 is missing, but 2,3,4,5 are present. It should continue on 1.
        assert!(result.is_some(), "Should not break on first missing ephemeris");
        assert_eq!(result.unwrap().0.len(), 4, "Should process the 4 valid ones");
    }
}
