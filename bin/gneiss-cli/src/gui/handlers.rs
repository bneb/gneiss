use serde::Serialize;
use serde_json::json;
use gneiss_core::coords::ecef_to_llh;
use gneiss_rtk::post_process::SmoothedEpoch;
use crate::export::ExportFormat;

#[derive(Debug, Clone, Serialize)]
pub struct BaseInfo {
    pub name: String,
    pub lat: f64,
    pub lon: f64,
    pub alt: f64,
}

/// Serializes base stations into a JSON array for map visualization.
pub fn format_bases_json(bases: &[BaseInfo]) -> String {
    serde_json::to_string(bases).unwrap_or_else(|_| "[]".to_string())
}

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

#[derive(Debug, Clone, serde::Deserialize)]
pub struct SolveRequest {
    pub rover: String,
    pub bases: Vec<String>,
    pub nav: Option<String>,
    pub antex: Option<String>,
    pub max_epochs: Option<usize>,
    pub auto_cors: Option<usize>,
    pub auto_products: bool,
}

pub async fn handle_solve_request(
    req: SolveRequest,
) -> Result<(Vec<SmoothedEpoch>, Vec<BaseInfo>), String> {
    let tmp_out = std::env::temp_dir().join(format!("gneiss_solve_{}.pos", std::process::id()));
    let args = crate::process::ProcessArgs {
        rover: req.rover,
        bases: req.bases.clone(),
        nav: req.nav,
        output: tmp_out.to_string_lossy().to_string(),
        format: Some("pos".to_string()),
        qc_report: None,
        geoid: None,
        config: None,
        enable_backward_smoothing: true,
        mode: None,
        max_epochs: req.max_epochs,
        base_position: None,
        systems: None,
        antex: req.antex,
        glonass: true,
        sp3: None,
        clk: None,
        auto_cors: req.auto_cors,
        auto_products: req.auto_products,
        calibrate_passes: None,
    };
    crate::process::run_process(args).await.map_err(|e| e.to_string())?;
    let epochs = super::server::load_trajectory_file(&tmp_out).map_err(|e| e.to_string())?;
    let _ = std::fs::remove_file(&tmp_out);
    let bases = super::server::load_base_stations(&req.bases);
    Ok((epochs, bases))
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::time::GpsTime;
    use nalgebra::{Matrix3, Vector3};

    fn make_test_epoch(tow: f64, quality: u8, sd_e: f64, sd_n: f64, sd_u: f64) -> SmoothedEpoch {
        SmoothedEpoch {
            time: GpsTime::new(2137, tow),
            position_ecef: Vector3::new(-2688179.79, -4265663.79, 3893784.58),
            velocity_ecef: None,
            attitude: None,
            cov_position: Matrix3::identity() * 0.0001,
            std_east: sd_e,
            std_north: sd_n,
            std_up: sd_u,
            separation_3d: 0.005,
            quality,
            n_satellites: 14,
        }
    }

    #[test]
    fn test_format_bases_json() {
        let bases = vec![BaseInfo {
            name: "P181".to_string(),
            lat: 37.9145,
            lon: -122.3767,
            alt: 72.74,
        }];
        let json_str = format_bases_json(&bases);
        assert!(json_str.contains("P181"));
        assert!(json_str.contains("37.9145"));
    }

    #[test]
    fn test_format_qc_json_accuracy() {
        let epochs = vec![
            make_test_epoch(100.0, 1, 0.01, 0.01, 0.02),
            make_test_epoch(101.0, 1, 0.01, 0.01, 0.02),
            make_test_epoch(102.0, 2, 0.02, 0.02, 0.04),
            make_test_epoch(103.0, 1, 0.01, 0.01, 0.02),
        ];
        let qc_str = format_qc_json(&epochs);
        let qc: serde_json::Value = serde_json::from_str(&qc_str).unwrap();
        assert_eq!(qc["total_epochs"], 4);
        assert_eq!(qc["fixed_epochs"], 3);
        assert_eq!(qc["fix_rate_pct"], 75.0);
    }

    #[test]
    fn test_format_trajectory_json() {
        let epochs = vec![make_test_epoch(100.0, 1, 0.008, 0.012, 0.024)];
        let traj_str = format_trajectory_json(&epochs);
        let list: serde_json::Value = serde_json::from_str(&traj_str).unwrap();
        assert_eq!(list.as_array().unwrap().len(), 1);
        assert_eq!(list[0]["quality"], 1);
        assert_eq!(list[0]["nsat"], 14);
    }

    #[test]
    fn test_solve_request_deserialization() {
        let json_payload = r#"{
            "rover": "datasets/rover.ubx",
            "bases": ["datasets/base.obs"],
            "nav": null,
            "antex": null,
            "max_epochs": 100,
            "auto_cors": null,
            "auto_products": true
        }"#;
        let req: SolveRequest = serde_json::from_str(json_payload).unwrap();
        assert_eq!(req.rover, "datasets/rover.ubx");
        assert_eq!(req.bases.len(), 1);
        assert_eq!(req.max_epochs, Some(100));
        assert!(req.auto_products);
    }
}
