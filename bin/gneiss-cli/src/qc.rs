//! Quality Control (QC) Report Generator.
//!
//! Evaluates post-processing statistical metrics, satellite coverage,
//! coordinate precision envelopes, forward/backward consistency,
//! and checks tolerances for surveyor deliverables.

use std::fs::File;
use std::io::Write;
use std::path::Path;

use serde::{Deserialize, Serialize};

use gneiss_rtk::post_process::SmoothedEpoch;

/// Single tolerance criteria check result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QcCheck {
    pub name: String,
    pub criteria: String,
    pub measured: String,
    pub passed: bool,
}

/// Structured quality control summary for a post-processed session.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QcReport {
    pub rover: String,
    pub bases: Vec<String>,
    pub total_epochs: usize,
    pub fix_rate_pct: f64,
    pub float_rate_pct: f64,
    pub n_sat_min: usize,
    pub n_sat_mean: f64,
    pub n_sat_max: usize,
    pub sigma_east_p50_m: f64,
    pub sigma_north_p50_m: f64,
    pub sigma_up_p50_m: f64,
    pub sigma_east_p95_m: f64,
    pub sigma_north_p95_m: f64,
    pub sigma_up_p95_m: f64,
    pub separation_3d_p50_m: f64,
    pub separation_3d_p95_m: f64,
    pub separation_3d_max_m: f64,
    pub checks: Vec<QcCheck>,
    pub overall_passed: bool,
}

impl QcReport {
    /// Generates a QC report from a trajectory and base/rover labels.
    pub fn compute(rover: &str, bases: &[String], traj: &[SmoothedEpoch]) -> Self {
        if traj.is_empty() {
            return Self::empty(rover, bases);
        }
        let total = traj.len();
        let fixed_count = traj.iter().filter(|e| e.quality == 1).count();
        let fix_rate_pct = (fixed_count as f64 / total as f64) * 100.0;
        let float_rate_pct = 100.0 - fix_rate_pct;

        let (n_sat_min, n_sat_mean, n_sat_max) = compute_sat_stats(traj);
        let (se_p50, se_p95) = compute_p50_p95(traj, |e| e.std_east);
        let (sn_p50, sn_p95) = compute_p50_p95(traj, |e| e.std_north);
        let (su_p50, su_p95) = compute_p50_p95(traj, |e| e.std_up);
        let (sep_p50, sep_p95, sep_max) = compute_sep_stats(traj);

        let mut checks = Vec::new();
        checks.push(QcCheck {
            name: "Fix Rate".into(),
            criteria: ">= 80.0%".into(),
            measured: format!("{:.1}%", fix_rate_pct),
            passed: fix_rate_pct >= 80.0,
        });
        checks.push(QcCheck {
            name: "Horizontal 95% Precision".into(),
            criteria: "<= 0.050 m".into(),
            measured: format!("{:.3} m", se_p95.max(sn_p95)),
            passed: se_p95.max(sn_p95) <= 0.050,
        });
        checks.push(QcCheck {
            name: "Vertical 95% Precision".into(),
            criteria: "<= 0.100 m".into(),
            measured: format!("{:.3} m", su_p95),
            passed: su_p95 <= 0.100,
        });
        checks.push(QcCheck {
            name: "3D Separation P50".into(),
            criteria: "<= 0.050 m".into(),
            measured: format!("{:.3} m", sep_p50),
            passed: sep_p50 <= 0.050,
        });

        let overall_passed = checks.iter().all(|c| c.passed);

        Self {
            rover: rover.to_string(),
            bases: bases.to_vec(),
            total_epochs: total,
            fix_rate_pct,
            float_rate_pct,
            n_sat_min,
            n_sat_mean,
            n_sat_max,
            sigma_east_p50_m: se_p50,
            sigma_north_p50_m: sn_p50,
            sigma_up_p50_m: su_p50,
            sigma_east_p95_m: se_p95,
            sigma_north_p95_m: sn_p95,
            sigma_up_p95_m: su_p95,
            separation_3d_p50_m: sep_p50,
            separation_3d_p95_m: sep_p95,
            separation_3d_max_m: sep_max,
            checks,
            overall_passed,
        }
    }

    fn empty(rover: &str, bases: &[String]) -> Self {
        Self {
            rover: rover.to_string(),
            bases: bases.to_vec(),
            total_epochs: 0,
            fix_rate_pct: 0.0,
            float_rate_pct: 0.0,
            n_sat_min: 0,
            n_sat_mean: 0.0,
            n_sat_max: 0,
            sigma_east_p50_m: 0.0,
            sigma_north_p50_m: 0.0,
            sigma_up_p50_m: 0.0,
            sigma_east_p95_m: 0.0,
            sigma_north_p95_m: 0.0,
            sigma_up_p95_m: 0.0,
            separation_3d_p50_m: 0.0,
            separation_3d_p95_m: 0.0,
            separation_3d_max_m: 0.0,
            checks: Vec::new(),
            overall_passed: false,
        }
    }

    /// Writes the QC report to disk as JSON, CSV, or HTML based on extension.
    pub fn write_to_file(&self, path: &Path) -> Result<(), std::io::Error> {
        let mut file = File::create(path)?;
        let ext = path.extension().and_then(|s| s.to_str()).unwrap_or("json");
        match ext {
            "csv" => self.write_csv(&mut file)?,
            "html" | "htm" => self.write_html(&mut file)?,
            _ => {
                let json = serde_json::to_string_pretty(self)
                    .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
                file.write_all(json.as_bytes())?;
            }
        }
        Ok(())
    }

    fn write_csv(&self, w: &mut dyn Write) -> Result<(), std::io::Error> {
        writeln!(w, "metric,value")?;
        writeln!(w, "rover,{}", self.rover)?;
        writeln!(w, "bases,{}", self.bases.join(";"))?;
        writeln!(w, "total_epochs,{}", self.total_epochs)?;
        writeln!(w, "fix_rate_pct,{:.1}", self.fix_rate_pct)?;
        writeln!(w, "sigma_east_p95_m,{:.4}", self.sigma_east_p95_m)?;
        writeln!(w, "sigma_north_p95_m,{:.4}", self.sigma_north_p95_m)?;
        writeln!(w, "sigma_up_p95_m,{:.4}", self.sigma_up_p95_m)?;
        writeln!(w, "separation_3d_p50_m,{:.4}", self.separation_3d_p50_m)?;
        writeln!(w, "overall_passed,{}", self.overall_passed)?;
        Ok(())
    }

    fn write_html(&self, w: &mut dyn Write) -> Result<(), std::io::Error> {
        let status_color = if self.overall_passed { "#10b981" } else { "#ef4444" };
        let status_text = if self.overall_passed { "PASSED" } else { "WARNING" };

        writeln!(w, "<!DOCTYPE html><html><head><meta charset='utf-8'>")?;
        writeln!(w, "<title>Gneiss PPK Executive QC Report - {}</title>", self.rover)?;
        writeln!(w, "<style>")?;
        writeln!(w, "body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; background: #f8fafc; color: #1e293b; margin: 0; padding: 32px; }}")?;
        writeln!(w, ".card {{ background: white; border-radius: 12px; box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1); padding: 24px; max-width: 900px; margin: 0 auto; }}")?;
        writeln!(w, ".badge {{ display: inline-block; padding: 6px 14px; border-radius: 9999px; font-weight: 700; color: white; background: {}; }}", status_color)?;
        writeln!(w, ".grid {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 16px; margin: 24px 0; }}")?;
        writeln!(w, ".kpi {{ background: #f1f5f9; border-radius: 8px; padding: 16px; text-align: center; }}")?;
        writeln!(w, ".kpi-val {{ font-size: 24px; font-weight: 700; color: #0f172a; }}")?;
        writeln!(w, ".kpi-label {{ font-size: 12px; color: #64748b; text-transform: uppercase; margin-top: 4px; }}")?;
        writeln!(w, "table {{ width: 100%; border-collapse: collapse; margin-top: 16px; }}")?;
        writeln!(w, "th, td {{ text-align: left; padding: 10px 12px; border-bottom: 1px solid #e2e8f0; }}")?;
        writeln!(w, "th {{ background: #f8fafc; font-size: 12px; text-transform: uppercase; color: #64748b; }}")?;
        writeln!(w, "</style></head><body><div class='card'>")?;

        writeln!(w, "<div style='display:flex; justify-content:space-between; align-items:center;'>")?;
        writeln!(w, "  <h2>Gneiss PPK Executive Quality Report</h2>")?;
        writeln!(w, "  <span class='badge'>{}</span>", status_text)?;
        writeln!(w, "</div>")?;

        writeln!(w, "<p><strong>Rover:</strong> {} &nbsp;|&nbsp; <strong>Bases:</strong> {}</p>", self.rover, self.bases.join(", "))?;

        writeln!(w, "<div class='grid'>")?;
        writeln!(w, "  <div class='kpi'><div class='kpi-val'>{:.1}%</div><div class='kpi-label'>Fix Rate</div></div>", self.fix_rate_pct)?;
        writeln!(w, "  <div class='kpi'><div class='kpi-val'>{:.1} mm</div><div class='kpi-label'>Horiz Precision (p95)</div></div>", self.sigma_east_p95_m.hypot(self.sigma_north_p95_m) * 1000.0)?;
        writeln!(w, "  <div class='kpi'><div class='kpi-val'>{:.1} mm</div><div class='kpi-label'>Vert Precision (p95)</div></div>", self.sigma_up_p95_m * 1000.0)?;
        writeln!(w, "</div>")?;

        writeln!(w, "<h3>Surveyor Tolerance Verification</h3>")?;
        writeln!(w, "<table><tr><th>Metric / Check</th><th>Tolerance</th><th>Measured</th><th>Status</th></tr>")?;
        for check in &self.checks {
            let pass_icon = if check.passed { "<span style='color:#10b981; font-weight:700;'>PASS</span>" } else { "<span style='color:#ef4444; font-weight:700;'>FAIL</span>" };
            writeln!(w, "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>", check.name, check.criteria, check.measured, pass_icon)?;
        }
        writeln!(w, "</table>")?;

        writeln!(w, "<p style='margin-top:24px; font-size:12px; color:#94a3b8;'>Generated by Gneiss RTK/PPK Suite &bull; Normalized ITRF2020 / WGS84</p>")?;
        writeln!(w, "</div></body></html>")?;
        Ok(())
    }
}

fn compute_sat_stats(traj: &[SmoothedEpoch]) -> (usize, f64, usize) {
    let mut min_sat = usize::MAX;
    let mut max_sat = 0;
    let mut sum_sat = 0;
    for e in traj {
        min_sat = min_sat.min(e.n_satellites);
        max_sat = max_sat.max(e.n_satellites);
        sum_sat += e.n_satellites;
    }
    let mean_sat = sum_sat as f64 / traj.len() as f64;
    (min_sat, mean_sat, max_sat)
}

fn compute_p50_p95<F: Fn(&SmoothedEpoch) -> f64>(traj: &[SmoothedEpoch], f: F) -> (f64, f64) {
    let mut vals: Vec<f64> = traj.iter().map(f).collect();
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let p50 = vals[(vals.len() as f64 * 0.50).floor() as usize];
    let p95 = vals[(vals.len() as f64 * 0.95).floor() as usize];
    (p50, p95)
}

fn compute_sep_stats(traj: &[SmoothedEpoch]) -> (f64, f64, f64) {
    let mut vals: Vec<f64> = traj.iter().map(|e| e.separation_3d).collect();
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let p50 = vals[(vals.len() as f64 * 0.50).floor() as usize];
    let p95 = vals[(vals.len() as f64 * 0.95).floor() as usize];
    let max = vals.last().copied().unwrap_or(0.0);
    (p50, p95, max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::time::GpsTime;
    use nalgebra::{Matrix3, Vector3};
    use tempfile::NamedTempFile;

    fn sample_epoch(quality: u8, sep: f64) -> SmoothedEpoch {
        SmoothedEpoch {
            time: GpsTime { week: 2105, tow: 345600.0 },
            position_ecef: Vector3::new(-2688181.50, -4265663.45, 3893784.80),
            velocity_ecef: None,
            attitude: None,
            cov_position: Matrix3::identity() * 0.0001,
            std_east: 0.010,
            std_north: 0.010,
            std_up: 0.020,
            separation_3d: sep,
            quality,
            n_satellites: 12,
        }
    }

    #[test]
    fn test_qc_report_computation() {
        let traj = vec![
            sample_epoch(1, 0.01),
            sample_epoch(1, 0.02),
            sample_epoch(2, 0.15),
        ];
        let report = QcReport::compute("rover.obs", &["base1.obs".into()], &traj);
        assert_eq!(report.total_epochs, 3);
        assert!((report.fix_rate_pct - 66.666).abs() < 0.1);
        assert_eq!(report.n_sat_min, 12);
        assert_eq!(report.n_sat_max, 12);
        assert!(!report.overall_passed); // Fix rate < 80%
    }

    #[test]
    fn test_qc_report_json_and_csv_export() {
        let traj = vec![sample_epoch(1, 0.01)];
        let report = QcReport::compute("rover.obs", &["base1.obs".into()], &traj);

        let tmp_json = NamedTempFile::new().expect("tempfile");
        assert!(report.write_to_file(tmp_json.path()).is_ok());

        let tmp_csv = tempfile::Builder::new().suffix(".csv").tempfile().expect("tempfile");
        assert!(report.write_to_file(tmp_csv.path()).is_ok());
    }
}
