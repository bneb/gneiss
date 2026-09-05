//! Pre-fit gross pseudorange error screening (fault detection/exclusion).
//!
//! Split out of `update.rs` (CLAUDE.md's 500-line file standard) as a
//! self-contained unit: production code plus its own already-isolated test
//! module, no shared test helpers with anything left behind — see
//! `docs/PROJECT_STATUS.md` Sprint 13.

use nalgebra::Vector3;

use super::state::DoubleDiffKey;
use super::update::{compute_tropo_dd, DoubleDiffMeasurement};

/// Reject threshold (metres) for a prefit DD pseudorange residual. Real
/// code noise (elevation-weighted, multipath-inflated) rarely exceeds
/// ~10m even at low elevation; a single bad ephemeris selection, a
/// GLONASS time-system slip, or severe urban multipath can instead
/// produce a residual thousands of metres wide on exactly one satellite.
/// Set well above plausible noise so this only ever fires on genuine
/// blunders.
pub const GROSS_PR_ERROR_THRESHOLD_M: f64 = 15.0;

/// Cap on removals per epoch: if more than a few satellites look like
/// blunders, the reference satellite or the position prior itself is the
/// more likely culprit, and further blind exclusion would just mask that.
pub const MAX_GROSS_PR_REJECTIONS_PER_EPOCH: usize = 3;

/// Iteratively drops the single worst-residual DD pseudorange pair while
/// its prefit residual (evaluated at `pos_pred`, i.e. *before* this
/// epoch's update -- a bad satellite can't hide behind an update it
/// corrupted itself) exceeds [`GROSS_PR_ERROR_THRESHOLD_M`], up to
/// [`MAX_GROSS_PR_REJECTIONS_PER_EPOCH`] removals.
///
/// The Huber-style robust weighting in `iekf_update_gated` already
/// down-weights moderately-large innovations, but that's a soft,
/// residual-based mechanism: with the deliberately loose position prior
/// this filter runs (so it can reacquire after outages), a single
/// kilometre-scale blunder can drag the *linearization point* itself
/// before robust weighting ever gets a chance to see it as an outlier.
/// This screen runs first, at a fixed prior position, so it can't be
/// fooled that way -- classic one-at-a-time "data snooping" FDE, as used
/// by RTKLIB, NovAtel Waypoint, and Leica Infinity for exactly this
/// reason. Returns the rejected keys, in rejection order, for logging.
pub fn screen_gross_pr_errors(
    measurements: &mut Vec<DoubleDiffMeasurement>,
    pos_pred: Vector3<f64>,
) -> Vec<DoubleDiffKey> {
    let mut rejected = Vec::new();
    for _ in 0..MAX_GROSS_PR_REJECTIONS_PER_EPOCH {
        let Some((idx, residual)) = worst_pr_residual(measurements, pos_pred) else { break };
        if residual.abs() <= GROSS_PR_ERROR_THRESHOLD_M {
            break;
        }
        rejected.push(measurements.remove(idx).key);
    }
    rejected
}

/// Index and signed prefit DD pseudorange residual of the pair with the
/// largest `|obs - predicted|`, or `None` if `measurements` is empty.
fn worst_pr_residual(measurements: &[DoubleDiffMeasurement], pos_pred: Vector3<f64>) -> Option<(usize, f64)> {
    measurements.iter()
        .enumerate()
        .map(|(i, m)| (i, pr_residual(m, pos_pred)))
        .max_by(|(_, a), (_, b)| a.abs().total_cmp(&b.abs()))
}

/// Prefit DD pseudorange residual (observed minus geometric model,
/// including the tropospheric DD correction) at a given position.
fn pr_residual(m: &DoubleDiffMeasurement, pos_pred: Vector3<f64>) -> f64 {
    let r_sat = (m.sat_pos - pos_pred).norm();
    let r_ref = (m.ref_pos - pos_pred).norm();
    let base_dd = (m.sat_pos - m.base_pos).norm() - (m.ref_pos - m.base_pos).norm();
    let geom_dd = (r_sat - r_ref) - base_dd
        + compute_tropo_dd(m.sat_pos, m.ref_pos, m.base_pos, pos_pred);
    m.dd_pr_m - geom_dd
}

#[cfg(test)]
mod gross_pr_error_tests {
    use super::*;

    fn make_meas(sat: u16, dd_pr_m: f64, sat_pos: Vector3<f64>, ref_pos: Vector3<f64>) -> DoubleDiffMeasurement {
        DoubleDiffMeasurement {
            key: DoubleDiffKey { constellation_id: 0, sat, ref_sat: 1, freq_band: 1 },
            dd_pr_m,
            dd_cp_cycles: None,
            sat_pos,
            ref_pos,
            base_pos: Vector3::zeros(),
            lambda: 0.190,
            pr_var_m2: 0.04,
            cp_var_cycles2: 0.0001,
            pr_ref_var_m2: 0.02,
            cp_ref_var_cycles2: 0.00005,
            dm_wet_rov: 0.0,
            dgrad_n_rov: 0.0,
            dgrad_e_rov: 0.0,
            tide_dd_m: 0.0,
            dd_pcv_m: 0.0,
        }
    }

    fn true_dd(pos: Vector3<f64>, sat_pos: Vector3<f64>, ref_pos: Vector3<f64>, base_pos: Vector3<f64>) -> f64 {
        let base_dd = (sat_pos - base_pos).norm() - (ref_pos - base_pos).norm();
        (sat_pos - pos).norm() - (ref_pos - pos).norm() - base_dd
    }

    #[test]
    fn keeps_all_pairs_when_residuals_are_small() {
        let pos = Vector3::new(100.0, 200.0, 300.0);
        let mut meas = vec![
            make_meas(2, true_dd(pos, Vector3::new(10_000.0, 20_000.0, 20_000.0), Vector3::new(5_000.0, 25_000.0, 20_000.0), Vector3::zeros()), Vector3::new(10_000.0, 20_000.0, 20_000.0), Vector3::new(5_000.0, 25_000.0, 20_000.0)),
            make_meas(3, true_dd(pos, Vector3::new(-8_000.0, 15_000.0, 22_000.0), Vector3::new(5_000.0, 25_000.0, 20_000.0), Vector3::zeros()), Vector3::new(-8_000.0, 15_000.0, 22_000.0), Vector3::new(5_000.0, 25_000.0, 20_000.0)),
        ];
        let rejected = screen_gross_pr_errors(&mut meas, pos);
        assert!(rejected.is_empty());
        assert_eq!(meas.len(), 2);
    }

    #[test]
    fn rejects_single_gross_outlier() {
        let pos = Vector3::new(100.0, 200.0, 300.0);
        let sat_pos = Vector3::new(10_000.0, 20_000.0, 20_000.0);
        let ref_pos = Vector3::new(5_000.0, 25_000.0, 20_000.0);
        let clean = true_dd(pos, sat_pos, ref_pos, Vector3::zeros());
        let mut meas = vec![
            make_meas(2, clean, sat_pos, ref_pos),
            make_meas(3, clean + 5_000.0, sat_pos, ref_pos), // 5km blunder
            make_meas(4, clean, sat_pos, ref_pos),
        ];
        let rejected = screen_gross_pr_errors(&mut meas, pos);
        assert_eq!(rejected, vec![DoubleDiffKey { constellation_id: 0, sat: 3, ref_sat: 1, freq_band: 1 }]);
        assert_eq!(meas.len(), 2);
        assert!(meas.iter().all(|m| m.key.sat != 3));
    }

    #[test]
    fn stops_after_max_rejections_even_if_more_are_bad() {
        let pos = Vector3::new(100.0, 200.0, 300.0);
        let sat_pos = Vector3::new(10_000.0, 20_000.0, 20_000.0);
        let ref_pos = Vector3::new(5_000.0, 25_000.0, 20_000.0);
        let clean = true_dd(pos, sat_pos, ref_pos, Vector3::zeros());
        let mut meas: Vec<DoubleDiffMeasurement> = (2..7)
            .map(|sat| make_meas(sat, clean + 1_000.0 * sat as f64, sat_pos, ref_pos))
            .collect();
        let rejected = screen_gross_pr_errors(&mut meas, pos);
        assert_eq!(rejected.len(), MAX_GROSS_PR_REJECTIONS_PER_EPOCH);
        assert_eq!(meas.len(), 5 - MAX_GROSS_PR_REJECTIONS_PER_EPOCH);
    }
}
