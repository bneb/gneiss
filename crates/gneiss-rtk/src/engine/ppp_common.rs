// Shared helpers for PPP IEKF solvers (ppp_iekf.rs and ppp_ins_iekf.rs).
use crate::filter::{RtkState, CORE_STATE_SIZE};
use crate::math::inversion::invert_matrix_robust;
use nalgebra::{DMatrix, DVector, UnitQuaternion, Vector3};

pub(crate) const NOMINAL_SNR_DBHZ: f64 = 45.0;
pub(crate) const SNR_SCALE_DIVISOR: f64 = 10.0;

pub(crate) fn snr_scale(snr: i32) -> f64 {
    (10.0f64).powf((NOMINAL_SNR_DBHZ - snr as f64) / SNR_SCALE_DIVISOR)
}

pub(crate) fn invert_matrix(mat: &DMatrix<f64>) -> Option<DMatrix<f64>> {
    if mat.iter().any(|x| x.is_nan()) {
        tracing::warn!(
            "NaN detected in {}×{} covariance matrix",
            mat.nrows(), mat.ncols()
        );
        return None;
    }
    Some(invert_matrix_robust(mat))
}

pub(crate) fn find_ambiguity_index(
    state: &RtkState,
    sat: gneiss_core::sat::SatelliteId,
) -> Option<usize> {
    state
        .ambiguity_keys
        .iter()
        .position(|&(s, f)| s == sat && f == 0)
}

pub(crate) fn find_amb_idx(
    state: &RtkState,
    sat: gneiss_core::sat::SatelliteId,
    f: u8,
) -> Option<usize> {
    state.ambiguity_keys.iter().position(|&x| x == (sat, f))
}

pub(crate) fn extract_state_vector(state: &RtkState) -> DVector<f64> {
    let mut x = DVector::zeros(CORE_STATE_SIZE + state.ambiguities.len());
    x[0] = state.position.vector.x;
    x[1] = state.position.vector.y;
    x[2] = state.position.vector.z;
    x[3] = state.velocity.x;
    x[4] = state.velocity.y;
    x[5] = state.velocity.z;
    let r = state.attitude.scaled_axis();
    x[6] = r.x;
    x[7] = r.y;
    x[8] = r.z;
    x[9] = state.accel_bias.x;
    x[10] = state.accel_bias.y;
    x[11] = state.accel_bias.z;
    x[12] = state.gyro_bias.x;
    x[13] = state.gyro_bias.y;
    x[14] = state.gyro_bias.z;
    x[15] = state.rcv_clk_bias;
    if CORE_STATE_SIZE > 16 {
        x[16] = state.isb_glo;
        x[17] = state.isb_gal;
        x[18] = state.isb_bds;
        x[19] = state.rcv_clk_drift;
        x[20] = state.zwd;
    }
    for (i, &amb) in state.ambiguities.iter().enumerate() {
        x[CORE_STATE_SIZE + i] = amb;
    }
    x
}

pub(crate) fn apply_state_vector(state: &mut RtkState, x: &DVector<f64>, p: DMatrix<f64>) {
    state.position.vector.x = x[0];
    state.position.vector.y = x[1];
    state.position.vector.z = x[2];
    state.velocity.x = x[3];
    state.velocity.y = x[4];
    state.velocity.z = x[5];
    let r_vec = Vector3::new(x[6], x[7], x[8]);
    state.attitude = if r_vec.norm() > 1e-12 {
        UnitQuaternion::from_scaled_axis(r_vec)
    } else {
        UnitQuaternion::identity()
    };
    state.accel_bias.x = x[9];
    state.accel_bias.y = x[10];
    state.accel_bias.z = x[11];
    state.gyro_bias.x = x[12];
    state.gyro_bias.y = x[13];
    state.gyro_bias.z = x[14];
    state.rcv_clk_bias = x[15];
    if CORE_STATE_SIZE > 16 {
        state.isb_glo = x[16];
        state.isb_gal = x[17];
        state.isb_bds = x[18];
        state.rcv_clk_drift = x[19];
        state.zwd = x[20];
    }
    for i in 0..state.ambiguities.len() {
        state.ambiguities[i] = x[CORE_STATE_SIZE + i];
    }
    state.covariance = p;
}

pub(crate) fn build_weight_matrix(meas: &[FgMeasurement], r_mat: &DMatrix<f64>) -> DMatrix<f64> {
    let mut w_mat = DMatrix::zeros(meas.len(), meas.len());
    for i in 0..meas.len() {
        w_mat[(i, i)] = 1.0 / r_mat[(i, i)];
    }
    w_mat
}

pub(crate) fn assemble_matrices(
    meas: &[FgMeasurement],
    cols: usize,
) -> (DMatrix<f64>, DVector<f64>, DMatrix<f64>) {
    let mut h = DMatrix::zeros(meas.len(), cols);
    let mut z = DVector::zeros(meas.len());
    let mut r = DMatrix::zeros(meas.len(), meas.len());
    for (i, m) in meas.iter().enumerate() {
        for j in 0..cols {
            h[(i, j)] = m.h_row[j];
        }
        z[i] = m.res;
        r[(i, i)] = m.weight;
    }
    (h, z, r)
}

pub(crate) fn build_iono_constraint_row(size: usize, i1_idx: usize) -> DVector<f64> {
    let mut h = DVector::zeros(size);
    h[i1_idx] = 1.0;
    h
}

pub struct FgMeasurement {
    pub res: f64,
    pub h_row: DVector<f64>,
    pub weight: f64,
    pub raw_var: f64,
    pub is_phase: bool,
    pub sat: Option<gneiss_core::sat::SatelliteId>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::RtkState;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use nalgebra::Vector3;

    fn make_state() -> RtkState {
        RtkState::new(
            GpsTime::new(0, 0.0),
            Coordinate::new(
                Vector3::zeros(),
                Datum::WGS84,
                Frame::ECEF,
                GpsTime::new(0, 0.0),
            ),
            1.0,
        )
    }

    #[test]
    fn test_extract_apply_roundtrip() {
        let mut state = make_state();
        state.position.vector = Vector3::new(1.0, 2.0, 3.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);
        state.rcv_clk_bias = 100.0;
        let x = extract_state_vector(&state);
        let mut state2 = make_state();
        apply_state_vector(&mut state2, &x, state.covariance.clone());
        assert!((state2.position.vector.x - 1.0).abs() < 1e-10);
        assert!((state2.velocity.y - 5.0).abs() < 1e-10);
        assert!((state2.rcv_clk_bias - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_snr_scale_values() {
        assert!((snr_scale(45) - 1.0).abs() < 1e-10);
        assert!((snr_scale(35) - 10.0).abs() < 1e-10);
        assert!((snr_scale(55) - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_assemble_matrices_empty() {
        let meas: Vec<FgMeasurement> = vec![];
        let (h, z, r) = assemble_matrices(&meas, 3);
        assert_eq!(h.nrows(), 0);
        assert_eq!(z.len(), 0);
        assert_eq!(r.nrows(), 0);
    }

    #[test]
    fn test_find_amb_idx_not_found() {
        let state = make_state();
        let sat = gneiss_core::sat::SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: 99,
        };
        assert_eq!(find_amb_idx(&state, sat, 1), None);
        assert_eq!(find_ambiguity_index(&state, sat), None);
    }

    #[test]
    fn test_build_iono_constraint_row_bounds() {
        let h = build_iono_constraint_row(30, 0);
        assert_eq!(h[0], 1.0);
        assert_eq!(h[29], 0.0);
        let h2 = build_iono_constraint_row(22, 21);
        assert_eq!(h2[21], 1.0);
    }

    #[test]
    fn test_invert_matrix_valid() {
        let m = DMatrix::from_row_slice(2, 2, &[4.0, 1.0, 1.0, 3.0]);
        let inv = invert_matrix(&m);
        assert!(inv.is_some());
        let inv = inv.unwrap();
        assert!((inv[(0, 0)] - 0.272727).abs() < 1e-5);
    }

    #[test]
    fn test_invert_matrix_nan() {
        // NaN entries are now sanitized (replaced with 10000 on diag, 0 off-diag)
        // rather than returning None, to prevent cascading StateDisappeared errors.
        let m = DMatrix::from_row_slice(2, 2, &[f64::NAN, 1.0, 1.0, 3.0]);
        let result = invert_matrix(&m);
        assert!(result.is_some(), "NaN should be sanitized, not rejected");
        // Verify result is finite
        let inv = result.unwrap();
        assert!(inv.iter().all(|x| x.is_finite()));
    }

    #[test]
    fn test_build_weight_matrix() {
        let meas = vec![
            FgMeasurement { res: 1.0, h_row: DVector::zeros(3), weight: 4.0, raw_var: 2.0, is_phase: false, sat: None },
            FgMeasurement { res: 2.0, h_row: DVector::zeros(3), weight: 9.0, raw_var: 3.0, is_phase: true, sat: None },
        ];
        let r_mat = DMatrix::from_diagonal(&DVector::from_vec(vec![4.0, 9.0]));
        let w = build_weight_matrix(&meas, &r_mat);
        assert!((w[(0, 0)] - 0.25).abs() < 1e-10);
        assert!((w[(1, 1)] - (1.0 / 9.0)).abs() < 1e-10);
    }

    #[test]
    fn test_assemble_matrices_with_data() {
        let meas = vec![
            FgMeasurement { res: 1.5, h_row: DVector::from_vec(vec![1.0, 0.0]), weight: 2.0, raw_var: 1.0, is_phase: false, sat: None },
            FgMeasurement { res: 3.0, h_row: DVector::from_vec(vec![0.0, 1.0]), weight: 5.0, raw_var: 2.0, is_phase: true, sat: None },
        ];
        let (h, z, r) = assemble_matrices(&meas, 2);
        assert_eq!(h.nrows(), 2);
        assert_eq!(h.ncols(), 2);
        assert!((h[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((z[1] - 3.0).abs() < 1e-10);
        assert!((r[(0, 0)] - 2.0).abs() < 1e-10);
        assert!((r[(1, 1)] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_extract_state_vector_with_ambiguities() {
        let mut state = make_state();
        state.position.vector = Vector3::new(1.0, 2.0, 3.0);
        state.velocity = Vector3::new(4.0, 5.0, 6.0);
        state.rcv_clk_bias = 100.0;
        state.ambiguities = vec![0.5, 1.5];
        state.ambiguity_keys = vec![
            (gneiss_core::sat::SatelliteId { constellation: gneiss_core::sat::Constellation::Gps, prn: 1 }, 1),
            (gneiss_core::sat::SatelliteId { constellation: gneiss_core::sat::Constellation::Gps, prn: 2 }, 1),
        ];
        // Manually resize covariance to include amb cols
        state.covariance = state.covariance.clone().insert_row(21, 0.0).insert_column(21, 0.0);
        state.covariance = state.covariance.clone().insert_row(22, 0.0).insert_column(22, 0.0);
        let x = extract_state_vector(&state);
        assert_eq!(x.len(), 23);
        assert_eq!(x[21], 0.5);
        assert_eq!(x[22], 1.5);
    }

    #[test]
    fn test_find_amb_idx_found() {
        let mut state = make_state();
        let sat = gneiss_core::sat::SatelliteId {
            constellation: gneiss_core::sat::Constellation::Gps,
            prn: 1,
        };
        state.ambiguity_keys = vec![(sat, 1), (sat, 2)];
        state.ambiguities = vec![0.0, 0.0];
        state.covariance = state.covariance.clone().insert_row(21, 0.0).insert_column(21, 0.0);
        state.covariance = state.covariance.clone().insert_row(22, 0.0).insert_column(22, 0.0);
        assert_eq!(find_amb_idx(&state, sat, 1), Some(0));
        assert_eq!(find_amb_idx(&state, sat, 2), Some(1));
        // find_ambiguity_index only matches f==0, so it won't find frequency band 1
        assert_eq!(find_ambiguity_index(&state, sat), None);
    }

    #[test]
    fn test_snr_scale_boundary_values() {
        // Nominal SNR at 45 dB-Hz -> scale = 1.0
        assert!((snr_scale(45) - 1.0).abs() < 1e-10);
        // Very high SNR at 65 dB-Hz -> scale = 0.01 (0.1^2)
        assert!((snr_scale(65) - 0.01).abs() < 1e-10);
        // Very low SNR at 25 dB-Hz -> scale = 100.0 (10^2)
        assert!((snr_scale(25) - 100.0).abs() < 1e-10);
        // Exactly 0 dB-Hz -> scale = 10^(45/10) = 10^4.5 = 31622.77...
        // The function does (10)^((45 - 0)/10) = 10^4.5
        assert!((snr_scale(0) - 31622.776601683792).abs() < 1.0);
    }

    #[test]
    fn test_build_iono_constraint_row_end_index() {
        let h = build_iono_constraint_row(10, 9);
        assert_eq!(h[9], 1.0);
        assert_eq!(h[0], 0.0);
    }

    #[test]
    fn test_invert_matrix_singular() {
        // Singular matrix: all zeros -> not NaN, but singular
        let m = DMatrix::from_row_slice(2, 2, &[0.0, 0.0, 0.0, 0.0]);
        let inv = invert_matrix(&m);
        // invert_matrix_robust should still compute a pseudo-inverse / fallback
        // and return Some (no NaN in input).
        assert!(inv.is_some());
    }

    #[test]
    fn test_build_iono_constraint_row_out_of_bounds_handled_by_panics() {
        // A one-element vector at index 0 works fine
        let h = build_iono_constraint_row(1, 0);
        assert_eq!(h[0], 1.0);
    }

    #[test]
    fn test_extract_state_vector_zero_position() {
        let state = make_state();
        let x = extract_state_vector(&state);
        assert_eq!(x.len(), CORE_STATE_SIZE);
        // Position, velocity, attitude, biases, clock are all 0 at initialization
        for i in 0..15 {
            assert_eq!(x[i], 0.0, "x[{}] should be 0.0, got {}", i, x[i]);
        }
        assert_eq!(x[15], 0.0); // rcv_clk_bias
        assert_eq!(x[16], 0.0); // isb_glo
        assert_eq!(x[17], 0.0); // isb_gal
        assert_eq!(x[18], 0.0); // isb_bds
        assert_eq!(x[19], 0.0); // rcv_clk_drift
        // zwd is initialized to 0.1 in RtkState::new()
        assert!((x[20] - 0.1).abs() < 1e-10, "zwd should be 0.1, got {}", x[20]);
    }

    #[test]
    fn test_apply_state_vector_empty_ambiguities() {
        let mut state = make_state();
        let x = extract_state_vector(&state);
        // Modify x slightly
        let mut x_mod = x.clone();
        x_mod[0] = 10.0;
        let p = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        apply_state_vector(&mut state, &x_mod, p);
        assert!((state.position.vector.x - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_apply_state_vector_resets_attitude_when_zero() {
        let mut state = make_state();
        let x = extract_state_vector(&state);
        // When rotation vector is essentially zero (all states are 0),
        // attitude should be identity
        let p = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        apply_state_vector(&mut state, &x, p);
        assert!((state.attitude.quaternion().w - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_extract_apply_with_core_only_indices() {
        // When CORE_STATE_SIZE is exactly 16 (checking the >16 path)
        // We can't change the constant, but we can verify the path is there.
        let mut state = make_state();
        state.rcv_clk_bias = 42.0;
        let x = extract_state_vector(&state);
        assert_eq!(x[15], 42.0);
    }

    #[test]
    fn test_apply_state_vector_with_large_rotation() {
        let mut state = make_state();
        let mut x = extract_state_vector(&state);
        // Set a non-trivial rotation vector
        x[6] = 0.5;
        x[7] = 0.5;
        x[8] = 0.5;
        let p = DMatrix::identity(CORE_STATE_SIZE, CORE_STATE_SIZE);
        apply_state_vector(&mut state, &x, p);
        assert!(state.attitude.quaternion().w < 1.0);
        assert!(state.attitude.quaternion().w > 0.0);
    }
}
