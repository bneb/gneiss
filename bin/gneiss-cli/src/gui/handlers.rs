//! API route request handlers and JSON serializers.

use serde_json::json;
use gneiss_core::coords::ecef_to_llh;
use gneiss_rtk::post_process::SmoothedEpoch;
use crate::export::ExportFormat;

/// Serializes trajectory epochs into a lightweight JSON array for Canvas visualization.
pub fn format_trajectory_json(trajectory: &[SmoothedEpoch]) -> String {
    let list: Vec<serde_json::Value> = trajectory.iter().map(|ep| {
        let llh = ecef_to_llh(ep.position_ecef);
        json!({
            "week": ep.time.week,
            "tow": ep.time.tow,
            "lat": llh.x.to_degrees(),
            "lon": llh.y.to_degrees(),
            "alt": llh.z,
            "x": ep.position_ecef.x,
            "y": ep.position_ecef.y,
            "z": ep.position_ecef.z,
            "quality": ep.quality,
            "nsat": ep.n_satellites,
            "sd_e": ep.std_east,
            "sd_n": ep.std_north,
            "sd_u": ep.std_up,
            "sep_3d": ep.separation_3d,
        })
    }).collect();

    serde_json::to_string(&list).unwrap_or_else(|_| "[]".to_string())
}

/// Formats a quick JSON summary of QC metrics.
pub fn format_qc_json(trajectory: &[SmoothedEpoch]) -> String {
    if trajectory.is_empty() {
        return json!({ "total_epochs": 0, "fix_rate_pct": 0.0, "h_rms_m": 0.0, "v_rms_m": 0.0 }).to_string();
    }

    let total = trajectory.len();
    let fixed_count = trajectory.iter().filter(|ep| ep.quality == 1).count();
    let fix_rate = (fixed_count as f64 / total as f64) * 100.0;

    let mut sum_h2 = 0.0;
    let mut sum_v2 = 0.0;
    for ep in trajectory {
        sum_h2 += ep.std_east * ep.std_east + ep.std_north * ep.std_north;
        sum_v2 += ep.std_up * ep.std_up;
    }
    let h_rms = libm::sqrt(sum_h2 / total as f64);
    let v_rms = libm::sqrt(sum_v2 / total as f64);

    json!({
        "total_epochs": total,
        "fixed_epochs": fixed_count,
        "fix_rate_pct": fix_rate,
        "h_rms_m": h_rms,
        "v_rms_m": v_rms,
    }).to_string()
}

/// Generates export bytes for live browser download.
pub fn generate_export_bytes(
    trajectory: &[SmoothedEpoch],
    format: ExportFormat,
) -> Result<Vec<u8>, std::io::Error> {
    let tmp_path = std::env::temp_dir().join(format!("gneiss_export_{}.{}", std::process::id(), format.default_extension()));
    crate::export::export_trajectory(trajectory, format, &tmp_path, None)?;
    let bytes = std::fs::read(&tmp_path)?;
    let _ = std::fs::remove_file(&tmp_path);
    Ok(bytes)
}
