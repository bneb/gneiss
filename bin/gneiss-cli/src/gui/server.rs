//! Embedded HTTP server for Gneiss Navigation Diagnostic Workspace.

use std::net::SocketAddr;
use std::path::Path;
use std::sync::Arc;
use tokio::io::{AsyncReadExt, AsyncWriteExt};
use tokio::net::TcpListener;
use tracing::info;

use gneiss_rtk::post_process::SmoothedEpoch;
use crate::export::ExportFormat;
use super::assets::HTML_INDEX;
use super::handlers::{format_qc_json, format_trajectory_json, generate_export_bytes};

pub struct GuiArgs {
    pub port: u16,
    pub trajectory_file: Option<String>,
}

pub async fn run_gui_server(args: GuiArgs) -> Result<(), Box<dyn std::error::Error>> {
    let trajectory: Arc<Vec<SmoothedEpoch>> = if let Some(path) = &args.trajectory_file {
        info!("Loading trajectory from {}...", path);
        Arc::new(load_trajectory_file(Path::new(path))?)
    } else {
        Arc::new(Vec::new())
    };

    let addr = SocketAddr::from(([127, 0, 0, 1], args.port));
    let listener = TcpListener::bind(addr).await?;

    info!("================================================================================");
    info!("             GNEISS DIAGNOSTIC WORKSPACE & VISUAL GUI READY                     ");
    info!("================================================================================");
    info!("  Local UI URL       : http://127.0.0.1:{}", args.port);
    info!("  Trajectory Epochs  : {}", trajectory.len());
    info!("================================================================================");

    loop {
        let (mut socket, _) = match listener.accept().await {
            Ok(s) => s,
            Err(_) => continue,
        };

        let traj = Arc::clone(&trajectory);
        tokio::spawn(async move {
            let mut buf = [0u8; 4096];
            let n = match socket.read(&mut buf).await {
                Ok(n) if n > 0 => n,
                _ => return,
            };

            let req_str = String::from_utf8_lossy(&buf[..n]);
            let first_line = req_str.lines().next().unwrap_or("");
            let parts: Vec<&str> = first_line.split_whitespace().collect();

            if parts.len() < 2 {
                return;
            }

            let path = parts[1];
            if path == "/" || path == "/index.html" {
                let resp = format!(
                    "HTTP/1.1 200 OK\r\nContent-Type: text/html; charset=utf-8\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
                    HTML_INDEX.len(),
                    HTML_INDEX
                );
                let _ = socket.write_all(resp.as_bytes()).await;
            } else if path == "/api/trajectory" {
                let json = format_trajectory_json(&traj);
                let resp = format!(
                    "HTTP/1.1 200 OK\r\nContent-Type: application/json\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
                    json.len(),
                    json
                );
                let _ = socket.write_all(resp.as_bytes()).await;
            } else if path == "/api/qc" {
                let json = format_qc_json(&traj);
                let resp = format!(
                    "HTTP/1.1 200 OK\r\nContent-Type: application/json\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
                    json.len(),
                    json
                );
                let _ = socket.write_all(resp.as_bytes()).await;
            } else if path.starts_with("/api/export") {
                let fmt_str = path.split("format=").nth(1).unwrap_or("pos");
                let fmt = ExportFormat::from_str(fmt_str).unwrap_or(ExportFormat::Pos);
                if let Ok(bytes) = generate_export_bytes(&traj, fmt) {
                    let header = format!(
                        "HTTP/1.1 200 OK\r\nContent-Type: application/octet-stream\r\nContent-Disposition: attachment; filename=\"trajectory.{}\"\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                        fmt.default_extension(),
                        bytes.len()
                    );
                    let _ = socket.write_all(header.as_bytes()).await;
                    let _ = socket.write_all(&bytes).await;
                }
            } else {
                let resp = "HTTP/1.1 404 NOT FOUND\r\nContent-Length: 0\r\nConnection: close\r\n\r\n";
                let _ = socket.write_all(resp.as_bytes()).await;
            }
        });
    }
}

fn load_trajectory_file(path: &Path) -> Result<Vec<SmoothedEpoch>, Box<dyn std::error::Error>> {
    use gneiss_core::coords::llh_to_ecef;
    use gneiss_core::time::GpsTime;
    use nalgebra::{Matrix3, Vector3};

    let content = std::fs::read_to_string(path)?;
    let mut epochs = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('%') || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 11 {
            let week: u32 = parts[0].parse().unwrap_or(0);
            let tow: f64 = parts[1].parse().unwrap_or(0.0);
            let lat_deg: f64 = parts[2].parse().unwrap_or(0.0);
            let lon_deg: f64 = parts[3].parse().unwrap_or(0.0);
            let h_m: f64 = parts[4].parse().unwrap_or(0.0);
            let quality: u8 = parts[6].parse().unwrap_or(1);
            let n_sat: usize = parts[7].parse().unwrap_or(10);
            let sd_e: f64 = parts[8].parse().unwrap_or(0.01);
            let sd_n: f64 = parts[9].parse().unwrap_or(0.01);
            let sd_u: f64 = parts[10].parse().unwrap_or(0.02);

            let pos = llh_to_ecef(Vector3::new(lat_deg.to_radians(), lon_deg.to_radians(), h_m));

            epochs.push(SmoothedEpoch {
                time: GpsTime::new(week, tow),
                position_ecef: pos,
                velocity_ecef: None,
                attitude: None,
                cov_position: Matrix3::identity() * 0.0001,
                std_east: sd_e,
                std_north: sd_n,
                std_up: sd_u,
                separation_3d: 0.005,
                quality,
                n_satellites: n_sat,
            });
        }
    }

    Ok(epochs)
}
