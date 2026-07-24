use crate::engine::EngineConfig;
use crate::filter::{PrBufferEntry, RtkState};
use gneiss_core::sat::SatelliteId;
use nalgebra::Vector3;

type DdKey = (SatelliteId, SatelliteId);

const MULTIPATH_FLOOR_M: f64 = 1.5;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeometryValidationDecision {
    Accepted,
    Rejected,
    InsufficientEvidence,
    NotApplicable,
}

#[derive(Debug, Clone)]
pub struct GeometryValidationOutcome {
    pub decision: GeometryValidationDecision,
    pub sample_count: usize,
    pub effective_sample_count: f64,
    pub normalized_residual: Option<f64>,
    pub threshold: Option<f64>,
    pub failing_key: Option<DdKey>,
    pub reason: &'static str,
}

impl GeometryValidationOutcome {
    fn new(decision: GeometryValidationDecision, reason: &'static str) -> Self {
        Self {
            decision,
            sample_count: 0,
            effective_sample_count: 0.0,
            normalized_residual: None,
            threshold: None,
            failing_key: None,
            reason,
        }
    }

    fn rejected(key: DdKey, normalized_residual: f64, threshold: f64, reason: &'static str) -> Self {
        Self {
            decision: GeometryValidationDecision::Rejected,
            sample_count: 0,
            effective_sample_count: 0.0,
            normalized_residual: Some(normalized_residual),
            threshold: Some(threshold),
            failing_key: Some(key),
            reason,
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct PairStatistic {
    sample_count: usize,
    effective_sample_count: f64,
    normalized_residual: f64,
}

pub fn validate_candidate_geometry(
    state: &RtkState,
    fixed_state: &RtkState,
    config: &EngineConfig,
    base_age_s: f64,
) -> GeometryValidationOutcome {
    if !config.enable_pr_validation {
        return GeometryValidationOutcome::new(GeometryValidationDecision::NotApplicable, "validation disabled");
    }
    if !valid_config(config) {
        return GeometryValidationOutcome::new(GeometryValidationDecision::InsufficientEvidence, "invalid validation configuration");
    }
    if !valid_base_age(base_age_s, config) {
        return GeometryValidationOutcome::new(GeometryValidationDecision::InsufficientEvidence, "base age exceeds validation limit");
    }

    let position_delta = fixed_state.position.vector - state.position.vector;
    let mut total_samples = 0;
    let mut total_effective_samples = 0.0;
    let mut max_normalized = 0.0_f64;
    let mut valid_pairs = 0;

    for (key, window) in &state.pr_dd_window {
        let statistic = match evaluate_pair(*key, window.entries(), position_delta, config) {
            Ok(Some(statistic)) => statistic,
            Ok(None) => continue,
            Err(outcome) => return outcome,
        };
        total_samples += statistic.sample_count;
        total_effective_samples += statistic.effective_sample_count;
        max_normalized = max_normalized.max(statistic.normalized_residual);
        valid_pairs += 1;
    }

    if valid_pairs < config.pr_validation_min_pairs {
        return GeometryValidationOutcome {
            sample_count: total_samples,
            effective_sample_count: total_effective_samples,
            reason: "insufficient independently keyed DD pairs",
            ..GeometryValidationOutcome::new(
                GeometryValidationDecision::InsufficientEvidence,
                "insufficient independently keyed DD pairs",
            )
        };
    }

    GeometryValidationOutcome {
        decision: GeometryValidationDecision::Accepted,
        sample_count: total_samples,
        effective_sample_count: total_effective_samples,
        normalized_residual: Some(max_normalized),
        threshold: Some(config.pr_validation_normalized_residual_threshold),
        failing_key: None,
        reason: "all DD geometry residuals accepted",
    }
}

fn valid_config(config: &EngineConfig) -> bool {
    config.pr_validation_min_samples > 0
        && config.pr_validation_min_pairs > 0
        && config.pr_validation_correlation_penalty.is_finite()
        && config.pr_validation_correlation_penalty >= 1.0
        && config.pr_validation_normalized_residual_threshold.is_finite()
        && config.pr_validation_normalized_residual_threshold > 0.0
        && config.pr_validation_gross_outlier_sigma.is_finite()
        && config.pr_validation_gross_outlier_sigma > 0.0
        && config.pr_validation_max_epoch_gap_s.is_finite()
        && config.pr_validation_max_epoch_gap_s > 0.0
}

fn valid_base_age(base_age_s: f64, config: &EngineConfig) -> bool {
    base_age_s.is_finite() && base_age_s.abs() <= config.max_base_age_s.min(0.5)
}

fn evaluate_pair<'a>(
    key: DdKey,
    entries: impl Iterator<Item = &'a PrBufferEntry>,
    position_delta: Vector3<f64>,
    config: &EngineConfig,
) -> Result<Option<PairStatistic>, GeometryValidationOutcome> {
    let entries: Vec<&PrBufferEntry> = entries.collect();
    if entries.len() < config.pr_validation_min_samples {
        return Ok(None);
    }
    if !entries_are_contiguous(&entries, config.pr_validation_max_epoch_gap_s) {
        return Err(GeometryValidationOutcome::rejected(key, 0.0, 0.0, "validation window contains a time gap"));
    }
    evaluate_residuals(key, &entries, position_delta, config)
}

fn entries_are_contiguous(entries: &[&PrBufferEntry], max_gap_s: f64) -> bool {
    entries.windows(2).all(|pair| {
        let gap = pair[1].time.tow - pair[0].time.tow;
        gap.is_finite() && gap > 0.0 && gap <= max_gap_s
    })
}

fn evaluate_residuals(
    key: DdKey,
    entries: &[&PrBufferEntry],
    position_delta: Vector3<f64>,
    config: &EngineConfig,
) -> Result<Option<PairStatistic>, GeometryValidationOutcome> {
    let mut weighted_sum = 0.0;
    let mut total_weight = 0.0;

    for entry in entries {
        let (residual, variance) = residual_and_variance(entry, position_delta);
        if !residual.is_finite() || !variance.is_finite() || variance <= 0.0 {
            return Err(GeometryValidationOutcome::rejected(key, 0.0, 0.0, "non-finite measurement variance"));
        }
        let normalized = residual.abs() / variance.sqrt();
        if normalized > config.pr_validation_gross_outlier_sigma {
            return Err(GeometryValidationOutcome::rejected(key, normalized, config.pr_validation_gross_outlier_sigma, "single-pair outlier"));
        }
        let weight = variance.recip();
        weighted_sum += weight * residual;
        total_weight += weight;
    }

    if !total_weight.is_finite() || total_weight <= 0.0 {
        return Ok(None);
    }
    let mean = weighted_sum / total_weight;
    let variance = config.pr_validation_correlation_penalty / total_weight;
    let normalized = mean.abs() / variance.sqrt();
    if normalized > config.pr_validation_normalized_residual_threshold {
        return Err(GeometryValidationOutcome::rejected(key, normalized, config.pr_validation_normalized_residual_threshold, "weighted DD residual exceeds threshold"));
    }

    Ok(Some(PairStatistic {
        sample_count: entries.len(),
        effective_sample_count: entries.len() as f64 / config.pr_validation_correlation_penalty,
        normalized_residual: normalized,
    }))
}

fn residual_and_variance(entry: &PrBufferEntry, position_delta: Vector3<f64>) -> (f64, f64) {
    let float_geometry = dd_geometry(entry.ekf_pos, entry);
    let candidate_geometry = dd_geometry(entry.ekf_pos + position_delta, entry);
    let residual = entry.z - (candidate_geometry - float_geometry);
    let variance = entry.variance_m2 + entry.position_variance_m2 + MULTIPATH_FLOOR_M.powi(2);
    (residual, variance)
}

fn dd_geometry(rover_pos: Vector3<f64>, entry: &PrBufferEntry) -> f64 {
    let rover_sat = (entry.sat_pos - rover_pos).norm();
    let rover_ref = (entry.ref_sat_pos - rover_pos).norm();
    let base_sat = (entry.base_sat_pos - entry.base_pos).norm();
    let base_ref = (entry.base_ref_sat_pos - entry.base_pos).norm();
    (rover_sat - base_sat) - (rover_ref - base_ref)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::PrRingBuffer;
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::sat::Constellation;
    use gneiss_core::time::GpsTime;

    fn sat(prn: u8) -> SatelliteId {
        SatelliteId { constellation: Constellation::Gps, prn }
    }

    fn state(position: Vector3<f64>) -> RtkState {
        let time = GpsTime::new(2200, 0.0);
        RtkState::new(time, Coordinate::new(position, Datum::WGS84, Frame::ECEF, time), 1.0)
    }

    fn sample(time: f64, rover: Vector3<f64>, delta: Vector3<f64>) -> PrBufferEntry {
        let base = Vector3::new(0.0, 0.0, 0.0);
        let sat_pos = Vector3::new(20_200_000.0, 1_000_000.0, 2_000_000.0);
        let ref_pos = Vector3::new(19_000_000.0, -3_000_000.0, 4_000_000.0);
        let entry = PrBufferEntry {
            time: GpsTime::new(2200, time), z: 0.0, ref_sat: sat(1), ekf_pos: rover,
            tdcp_pos: None, sat_pos, ref_sat_pos: ref_pos, base_pos: base,
            base_sat_pos: sat_pos, base_ref_sat_pos: ref_pos, variance_m2: 0.25,
            position_variance_m2: 0.01,
        };
        let innovation = dd_geometry(rover + delta, &entry) - dd_geometry(rover, &entry);
        PrBufferEntry { z: innovation, ..entry }
    }

    fn config() -> EngineConfig {
        EngineConfig {
            enable_pr_validation: true,
            pr_validation_min_samples: 3,
            pr_validation_min_pairs: 2,
            pr_validation_correlation_penalty: 2.0,
            pr_validation_max_epoch_gap_s: 2.0,
            ..EngineConfig::default()
        }
    }

    fn add_pair(state: &mut RtkState, target: SatelliteId, delta: Vector3<f64>) {
        let mut window = PrRingBuffer::new(4);
        for time in 0..3 {
            let entry = sample(time as f64 + 1.0, state.position.vector, delta);
            window.push(
                entry.time, entry.z, entry.ref_sat, entry.ekf_pos, entry.sat_pos,
                entry.ref_sat_pos, entry.base_pos, entry.base_sat_pos, entry.base_ref_sat_pos,
                entry.variance_m2, entry.position_variance_m2,
            );
        }
        state.pr_dd_window.insert((target, sat(1)), window);
    }

    #[test]
    fn accepts_time_consistent_candidate_across_moving_samples() {
        let mut float = state(Vector3::new(6_378_000.0, 0.0, 0.0));
        let delta = Vector3::new(0.4, -0.2, 0.1);
        add_pair(&mut float, sat(2), delta);
        add_pair(&mut float, sat(3), delta);
        let fixed = state(float.position.vector + delta);

        let outcome = validate_candidate_geometry(&float, &fixed, &config(), 0.1);
        assert_eq!(outcome.decision, GeometryValidationDecision::Accepted);
        assert_eq!(outcome.sample_count, 6);
        assert!(outcome.effective_sample_count < outcome.sample_count as f64);
    }

    #[test]
    fn rejects_wrong_candidate_without_mutating_float_state() {
        let mut float = state(Vector3::new(6_378_000.0, 0.0, 0.0));
        let correct_delta = Vector3::new(0.2, 0.0, 0.0);
        add_pair(&mut float, sat(2), correct_delta);
        add_pair(&mut float, sat(3), correct_delta);
        let before = float.position.vector;
        let wrong = state(before + Vector3::new(0.0, 50.0, 0.0));

        let outcome = validate_candidate_geometry(&float, &wrong, &config(), 0.1);
        assert_eq!(outcome.decision, GeometryValidationDecision::Rejected);
        assert_eq!(float.position.vector, before);
        assert!(!float.is_fixed);
    }

    #[test]
    fn rejects_windows_with_epoch_gaps() {
        let mut float = state(Vector3::new(6_378_000.0, 0.0, 0.0));
        let delta = Vector3::new(0.2, 0.0, 0.0);
        let mut window = PrRingBuffer::new(3);
        for time in [1.0, 2.0, 10.0] {
            let entry = sample(time, float.position.vector, delta);
            window.push(entry.time, entry.z, entry.ref_sat, entry.ekf_pos, entry.sat_pos,
                entry.ref_sat_pos, entry.base_pos, entry.base_sat_pos, entry.base_ref_sat_pos,
                entry.variance_m2, entry.position_variance_m2);
        }
        float.pr_dd_window.insert((sat(2), sat(1)), window);
        let fixed = state(float.position.vector + delta);

        let outcome = validate_candidate_geometry(&float, &fixed, &config(), 0.1);
        assert_eq!(outcome.decision, GeometryValidationDecision::Rejected);
        assert_eq!(outcome.reason, "validation window contains a time gap");
    }

    #[test]
    fn returns_insufficient_evidence_for_stale_base_or_too_few_pairs() {
        let float = state(Vector3::new(6_378_000.0, 0.0, 0.0));
        let fixed = state(float.position.vector);
        let stale = validate_candidate_geometry(&float, &fixed, &config(), 0.6);
        assert_eq!(stale.decision, GeometryValidationDecision::InsufficientEvidence);

        let no_pairs = validate_candidate_geometry(&float, &fixed, &config(), 0.1);
        assert_eq!(no_pairs.decision, GeometryValidationDecision::InsufficientEvidence);
    }
}
