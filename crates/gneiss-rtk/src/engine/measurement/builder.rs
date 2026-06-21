use crate::engine::measurement::types::{EkfUpdates, MeasurementEnvironment};
use crate::engine::measurement::compute_innovations;
use crate::engine::measurement_math;
use crate::filter::{DdObservation, RtkState};
use nalgebra::{DMatrix, DVector};

pub struct EkfMeasurementMatrices {
    pub z: DVector<f64>,
    pub h: DMatrix<f64>,
    pub r: DMatrix<f64>,
    pub mt: Vec<(gneiss_core::sat::SatelliteId, u8, f64)>,
}

fn group_measurements_by_constellation(
    matched_obs: &[(DdObservation, DdObservation)],
) -> std::collections::HashMap<gneiss_core::sat::Constellation, Vec<(DdObservation, DdObservation)>>
{
    let mut const_groups = std::collections::HashMap::<
        gneiss_core::sat::Constellation,
        Vec<(DdObservation, DdObservation)>,
    >::new();
    for obs in matched_obs {
        const_groups
            .entry(obs.0.sat.constellation)
            .or_default()
            .push(obs.clone());
    }
    const_groups
}

fn select_reference_satellite(
    group: &[(DdObservation, DdObservation)],
    state: &RtkState,
    env: &MeasurementEnvironment,
) -> usize {
    let mut best_score = -1.0;
    let mut ref_idx = 0;
    let rov_llh = gneiss_core::coords::ecef_to_llh(state.position.vector);

    for (i, (r, _)) in group.iter().enumerate() {
        if let Some(eph) = env
            .ephemerides
            .iter()
            .filter(|e| e.sat() == r.sat)
            .min_by(|a, b| {
                let da = (a.toe().tow - state.time.tow).abs();
                let db = (b.toe().tow - state.time.tow).abs();
                da.partial_cmp(&db).unwrap()
            })
        {
            let tau = r.pr_l1 / gneiss_core::constants::SPEED_OF_LIGHT_M_S;
            let t_tx = gneiss_core::time::GpsTime::new(state.time.week, state.time.tow - tau);
            let (sat_pos, _, _, _): (nalgebra::Vector3<f64>, _, _, _) = eph.position(t_tx);
            let (_, el) = gneiss_core::coords::az_el(rov_llh, state.position.vector, sat_pos);

            let score = if r.cp_l1.is_some() { el + 100.0 } else { el };
            if score > best_score {
                best_score = score;
                ref_idx = i;
            }
        }
    }
    ref_idx
}

fn update_reject_counts(
    state: &mut RtkState,
    h_row: &DMatrix<f64>,
    state_size: usize,
    passed: bool,
) {
    for c in crate::filter::CORE_STATE_SIZE..state_size {
        if h_row[(0, c)] > 0.5 {
            let key = state.ambiguity_keys[c - crate::filter::CORE_STATE_SIZE];
            if passed {
                state.reject_counts.insert(key, 0);
            } else {
                let count = *state.reject_counts.get(&key).unwrap_or(&0) + 1;
                state.reject_counts.insert(key, count);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn filter_innovations_chi_squared(
    state: &mut RtkState,
    state_size: usize,
    chi_pr: f64,
    chi_cp: f64,
    z_all: &[f64],
    h_all: &[Vec<f64>],
    r_all: &[f64],
    type_all: &[(gneiss_core::sat::SatelliteId, u8, f64)],
) -> Vec<usize> {
    let mut safe_indices = Vec::new();
    for i in 0..z_all.len() {
        let mut h_row = DMatrix::zeros(1, state_size);
        for c in 0..state_size {
            h_row[(0, c)] = h_all[i][c];
        }
        let s_ii = (&h_row * &state.covariance * h_row.transpose())[(0, 0)] + r_all[i];
        let chi2 = z_all[i] * z_all[i] / s_ii;

        let threshold = match type_all[i].1 {
            0 => chi_pr * chi_pr,
            1 | 2 => chi_cp * chi_cp,
            3 => chi_pr * 1000.0,
            _ => chi_pr * chi_pr,
        };

        let passed = chi2 <= threshold;
        if passed {
            safe_indices.push(i);
        } else {
            tracing::debug!(
                "Rejected meas type {} with inn: {:.3}, chi2: {:.1}",
                type_all[i].1,
                z_all[i],
                chi2
            );
        }

        if type_all[i].1 == 1 || type_all[i].1 == 2 {
            update_reject_counts(state, &h_row, state_size, passed);
        }
    }
    safe_indices
}

fn build_final_measurement_matrices(
    state_size: usize,
    safe_indices: Vec<usize>,
    z_all: &[f64],
    h_all: &[Vec<f64>],
    r_all: &[f64],
    type_all: &[(gneiss_core::sat::SatelliteId, u8, f64)],
) -> Option<EkfMeasurementMatrices> {
    if safe_indices.len() >= 4 {
        let mut z_vec = DVector::zeros(safe_indices.len());
        let mut h_mat = DMatrix::zeros(safe_indices.len(), state_size);
        let mut r_diagonals = Vec::new();
        let mut t_vec = Vec::new();

        for (new_i, &old_i) in safe_indices.iter().enumerate() {
            z_vec[new_i] = z_all[old_i];
            for c in 0..state_size {
                h_mat[(new_i, c)] = h_all[old_i][c];
            }
            r_diagonals.push(r_all[old_i]);
            t_vec.push(type_all[old_i]);
        }
        let r_mat = measurement_math::build_dense_covariance_matrix(&r_diagonals, &t_vec);
        Some(EkfMeasurementMatrices {
            z: z_vec,
            h: h_mat,
            r: r_mat,
            mt: t_vec,
        })
    } else {
        tracing::warn!(
            "measurement model empty! all_z={}, safe_indices={}",
            z_all.len(),
            safe_indices.len()
        );
        None
    }
}

pub fn build_measurement_model(
    state: &mut RtkState,
    matched_obs: &[(DdObservation, DdObservation)],
    env: &MeasurementEnvironment,
    chi_square_pr_threshold: f64,
    chi_square_cp_threshold: f64,
) -> Option<EkfMeasurementMatrices> {
    let mut all = EkfUpdates::new();
    let const_groups = group_measurements_by_constellation(matched_obs);

    for (_, group) in const_groups {
        if group.len() < 2 {
            continue;
        }
        let ref_idx = select_reference_satellite(&group, state, env);
        let mut group_clone = group.clone();
        let (ref_rover, ref_base) = group_clone.remove(ref_idx);

        if let Some(updates) = compute_innovations(state, &group_clone, &ref_rover, &ref_base, env)
        {
            all.extend(updates);
        }
    }

    let state_size = crate::filter::CORE_STATE_SIZE + state.ambiguities.len();
    tracing::trace!("Pre-filter z_all len: {}", all.z.len());

    let safe_indices = filter_innovations_chi_squared(
        state,
        state_size,
        chi_square_pr_threshold,
        chi_square_cp_threshold,
        &all.z,
        &all.h,
        &all.r,
        &all.mt,
    );

    build_final_measurement_matrices(state_size, safe_indices, &all.z, &all.h, &all.r, &all.mt)
}
