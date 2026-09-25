//! Common test fixtures, synthetic geometry builders, and verification helpers
//! for the Urban Canyon E2E test suite.
#![allow(dead_code)]

pub use nalgebra::{DMatrix, Vector3};
pub use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
pub use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
pub use gneiss_core::sat::{Constellation, SatelliteId};
pub use gneiss_core::time::GpsTime;
pub use gneiss_rtk::estimators::rtk_iekf::state::DoubleDiffKey;
pub use gneiss_rtk::estimators::rtk_iekf::update::DoubleDiffMeasurement;

/// GPS L1 carrier frequency (Hz).
pub const GPS_L1_FREQ_HZ: f64 = 1575.42e6;

/// GPS L1 carrier wavelength (metres).
pub const GPS_L1_WAVELENGTH_M: f64 = SPEED_OF_LIGHT_M_S / GPS_L1_FREQ_HZ;

/// Nominal Tokyo reference station position in ECEF (metres).
pub fn ref_ecef() -> Vector3<f64> {
    Vector3::new(-3959950.0, 3352850.0, 3699700.0)
}

/// Nominal Tokyo rover position (500 m baseline).
pub fn rover_ecef() -> Vector3<f64> {
    ref_ecef() + Vector3::new(300.0, 400.0, 0.0)
}

/// Generate satellite ECEF position from azimuth (deg) and elevation (deg).
pub fn sat_pos_az_el(rx_pos: Vector3<f64>, az_deg: f64, el_deg: f64, range_m: f64) -> Vector3<f64> {
    let az = az_deg.to_radians();
    let el = el_deg.to_radians();
    let e = range_m * el.cos() * az.sin();
    let n = range_m * el.cos() * az.cos();
    let u = range_m * el.sin();
    let llh = gneiss_core::coords::ecef_to_llh(rx_pos);
    let lat = llh.x;
    let lon = llh.y;
    let dx = -lon.sin() * e - lat.sin() * lon.cos() * n + lat.cos() * lon.cos() * u;
    let dy = lon.cos() * e - lat.sin() * lon.sin() * n + lat.cos() * lon.sin() * u;
    let dz = lat.cos() * n + lat.sin() * u;
    rx_pos + Vector3::new(dx, dy, dz)
}

/// Compute true geometric double difference (metres).
pub fn compute_geometric_dd(
    rx_pos: Vector3<f64>,
    base_pos: Vector3<f64>,
    sat_pos: Vector3<f64>,
    ref_pos: Vector3<f64>,
) -> f64 {
    let r_sat = (sat_pos - rx_pos).norm();
    let r_ref = (ref_pos - rx_pos).norm();
    let b_sat = (sat_pos - base_pos).norm();
    let b_ref = (ref_pos - base_pos).norm();
    (r_sat - r_ref) - (b_sat - b_ref)
}

/// Create a satellite observation with code, phase, doppler, and SNR on L1.
pub fn make_sat_obs_l1(
    prn: u8,
    constellation: Constellation,
    pr_m: f64,
    cp_cycles: f64,
    doppler_hz: f64,
    snr_dbhz: f64,
    lli: Option<u8>,
) -> SatObs {
    let code_c1c = ObsCode {
        obs_type: ObsType::Pseudorange,
        signal: SignalCode { freq_band: 1, attribute: 'C' },
    };
    let code_l1c = ObsCode {
        obs_type: ObsType::CarrierPhase,
        signal: SignalCode { freq_band: 1, attribute: 'C' },
    };
    let code_d1c = ObsCode {
        obs_type: ObsType::Doppler,
        signal: SignalCode { freq_band: 1, attribute: 'C' },
    };
    let code_s1c = ObsCode {
        obs_type: ObsType::Snr,
        signal: SignalCode { freq_band: 1, attribute: 'C' },
    };

    let observations = vec![
        Observation { code: code_c1c, value: pr_m, lock_time: Some(100), lli: None },
        Observation { code: code_l1c, value: cp_cycles, lock_time: Some(100), lli },
        Observation { code: code_d1c, value: doppler_hz, lock_time: None, lli: None },
        Observation { code: code_s1c, value: snr_dbhz, lock_time: None, lli: None },
    ];

    SatObs {
        sat: SatelliteId { constellation, prn },
        observations,
    }
}

/// Create an EpochObs at the given TOW.
pub fn make_epoch(tow: f64, sats: Vec<SatObs>) -> EpochObs {
    EpochObs {
        time: GpsTime::new(2300, tow),
        satellites: sats,
    }
}

/// Construct a DoubleDiffMeasurement with custom variances and geometry.
pub fn make_dd_meas(
    pair: (u8, u16, u16),
    dd_pr_m: f64,
    dd_cp_cycles: Option<f64>,
    vars: (f64, f64),
    sat_pos: Vector3<f64>,
    ref_pos: Vector3<f64>,
) -> DoubleDiffMeasurement {
    DoubleDiffMeasurement {
        key: DoubleDiffKey {
            constellation_id: pair.0,
            sat: pair.1,
            ref_sat: pair.2,
            freq_band: 1,
        },
        dd_pr_m,
        dd_cp_cycles,
        sat_pos,
        ref_pos,
        base_pos: ref_ecef(),
        lambda: GPS_L1_WAVELENGTH_M,
        pr_var_m2: vars.0,
        cp_var_cycles2: vars.1,
        pr_ref_var_m2: vars.0 * 0.5,
        cp_ref_var_cycles2: vars.1 * 0.5,
        dm_wet_rov: 0.0,
        dgrad_n_rov: 0.0,
        dgrad_e_rov: 0.0,
        tide_dd_m: 0.0,
        dd_pcv_m: 0.0,
    }
}

/// Checks if a matrix is strictly symmetric and positive definite (all eigenvalues > 0).
pub fn is_matrix_positive_definite(m: &DMatrix<f64>, tol: f64) -> bool {
    let sym = 0.5 * (m + m.transpose());
    let eig = sym.symmetric_eigen();
    eig.eigenvalues.iter().all(|&v| v > tol)
}

/// Minimum eigenvalue of a symmetric matrix.
pub fn matrix_min_eigenvalue(m: &DMatrix<f64>) -> f64 {
    let sym = 0.5 * (m + m.transpose());
    let eig = sym.symmetric_eigen();
    eig.eigenvalues.iter().copied().fold(f64::INFINITY, f64::min)
}
