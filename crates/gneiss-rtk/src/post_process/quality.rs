//! Quality Control and Statistical Verification for Post-Processed Solutions.
//!
//! Evaluates forward-backward separation, 1σ/2σ/95% uncertainty bounds,
//! fix rates, and assigns standard GNSS Quality indicators (Q1..Q5).

use crate::post_process::combiner::SmoothedEpoch;
use crate::post_process::percentile;

/// Comprehensive quality report for a processed session.
#[derive(Debug, Clone)]
pub struct QualityReport {
    pub total_epochs: usize,
    pub fixed_epochs: usize,
    pub float_epochs: usize,
    pub dgps_epochs: usize,
    pub spp_epochs: usize,
    pub fix_rate_pct: f64,
    pub median_separation_m: f64,
    pub p95_separation_m: f64,
    pub median_std_horizontal_m: f64,
    pub median_std_3d_m: f64,
}

/// Generate statistical quality report from smoothed trajectory.
pub fn generate_quality_report(epochs: &[SmoothedEpoch]) -> QualityReport {
    if epochs.is_empty() {
        return empty_report();
    }

    let mut separations = Vec::with_capacity(epochs.len());
    let mut h_stds = Vec::with_capacity(epochs.len());
    let mut d3_stds = Vec::with_capacity(epochs.len());

    let mut fixed = 0usize;
    let mut float = 0usize;
    let mut dgps = 0usize;
    let mut spp = 0usize;

    for ep in epochs {
        match ep.quality {
            1 => fixed += 1,
            2 => float += 1,
            3 => dgps += 1,
            _ => spp += 1,
        }
        separations.push(ep.separation_3d);
        let h_std = (ep.std_east * ep.std_east + ep.std_north * ep.std_north).sqrt();
        let d3_std = (h_std * h_std + ep.std_up * ep.std_up).sqrt();
        h_stds.push(h_std);
        d3_stds.push(d3_std);
    }

    separations.sort_by(|a, b| a.total_cmp(b));
    h_stds.sort_by(|a, b| a.total_cmp(b));
    d3_stds.sort_by(|a, b| a.total_cmp(b));

    let n = epochs.len();
    QualityReport {
        total_epochs: n,
        fixed_epochs: fixed,
        float_epochs: float,
        dgps_epochs: dgps,
        spp_epochs: spp,
        fix_rate_pct: (fixed as f64) / (n as f64) * 100.0,
        median_separation_m: percentile(&separations, 0.50),
        p95_separation_m: percentile(&separations, 0.95),
        median_std_horizontal_m: percentile(&h_stds, 0.50),
        median_std_3d_m: percentile(&d3_stds, 0.50),
    }
}

fn empty_report() -> QualityReport {
    QualityReport {
        total_epochs: 0,
        fixed_epochs: 0,
        float_epochs: 0,
        dgps_epochs: 0,
        spp_epochs: 0,
        fix_rate_pct: 0.0,
        median_separation_m: 0.0,
        p95_separation_m: 0.0,
        median_std_horizontal_m: 0.0,
        median_std_3d_m: 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::time::GpsTime;
    use nalgebra::{Matrix3, Vector3};

    #[test]
    fn test_quality_report_stats() {
        let ep = SmoothedEpoch {
            time: GpsTime::new(2000, 100.0),
            position_ecef: Vector3::zeros(),
            velocity_ecef: None,
            attitude: None,
            cov_position: Matrix3::identity(),
            std_east: 0.02,
            std_north: 0.02,
            std_up: 0.05,
            separation_3d: 0.01,
            quality: 1,
            n_satellites: 8,
        };
        let report = generate_quality_report(&[ep]);
        assert_eq!(report.fixed_epochs, 1);
        assert_eq!(report.fix_rate_pct, 100.0);
        assert!(report.median_separation_m < 0.02);
    }
}
