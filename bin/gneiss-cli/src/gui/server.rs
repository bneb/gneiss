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
use super::handlers::{format_bases_json, format_qc_json, format_trajectory_json, generate_export_bytes, BaseInfo};

use tokio::sync::RwLock;

pub struct GuiArgs {
    pub port: u16,
    pub trajectory_file: Option<String>,
    pub base_files: Vec<String>,
}

pub async fn run_gui_server(args: GuiArgs) -> Result<(), Box<dyn std::error::Error>> {
    let initial_epochs = if let Some(path) = &args.trajectory_file {
        info!("Loading trajectory from {}...", path);
        load_trajectory_file(Path::new(path))?
    } else {
        Vec::new()
    };

    let initial_bases = load_base_stations(&args.base_files);

    let trajectory: Arc<RwLock<Vec<SmoothedEpoch>>> = Arc::new(RwLock::new(initial_epochs));
    let bases: Arc<RwLock<Vec<BaseInfo>>> = Arc::new(RwLock::new(initial_bases));
    let is_solving: Arc<std::sync::atomic::AtomicBool> = Arc::new(std::sync::atomic::AtomicBool::new(false));
    let last_error: Arc<RwLock<Option<String>>> = Arc::new(RwLock::new(None));

    let addr = SocketAddr::from(([127, 0, 0, 1], args.port));
    let listener = TcpListener::bind(addr).await?;

    info!("================================================================================");
    info!("             GNEISS DIAGNOSTIC WORKSPACE & VISUAL GUI READY                     ");
    info!("================================================================================");
    info!("  Local UI URL       : http://127.0.0.1:{}", args.port);
    info!("================================================================================");

    loop {
        let (mut socket, _) = match listener.accept().await {
            Ok(s) => s,
            Err(_) => continue,
        };

        let traj = Arc::clone(&trajectory);
        let base_list = Arc::clone(&bases);
        let solving = Arc::clone(&is_solving);
        let err_state = Arc::clone(&last_error);

        tokio::spawn(async move {
            let mut buf = [0u8; 8192];
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
                let current = traj.read().await;
                let json = format_trajectory_json(&current);
                let resp = format!(
                    "HTTP/1.1 200 OK\r\nContent-Type: application/json\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
                    json.len(),
                    json
                );
                let _ = socket.write_all(resp.as_bytes()).await;
            } else if path == "/api/bases" {
                let current = base_list.read().await;
                let json = format_bases_json(&current);
                let resp = format!(
                    "HTTP/1.1 200 OK\r\nContent-Type: application/json\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
                    json.len(),
                    json
                );
                let _ = socket.write_all(resp.as_bytes()).await;
            } else if path == "/api/qc" {
                let current = traj.read().await;
                let json = format_qc_json(&current);
                let resp = format!(
                    "HTTP/1.1 200 OK\r\nContent-Type: application/json\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
                    json.len(),
                    json
                );
                let _ = socket.write_all(resp.as_bytes()).await;
            } else if path == "/api/status" {
                let is_busy = solving.load(std::sync::atomic::Ordering::SeqCst);
                let err = err_state.read().await;
                let json = serde_json::json!({
                    "is_solving": is_busy,
                    "error": err.as_ref()
                }).to_string();
                let resp = format!(
                    "HTTP/1.1 200 OK\r\nContent-Type: application/json\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
                    json.len(),
                    json
                );
                let _ = socket.write_all(resp.as_bytes()).await;
            } else if path == "/api/solve" && first_line.starts_with("POST") {
                if solving.compare_exchange(false, true, std::sync::atomic::Ordering::SeqCst, std::sync::atomic::Ordering::SeqCst).is_err() {
                    let resp = "HTTP/1.1 409 CONFLICT\r\nContent-Type: application/json\r\nContent-Length: 35\r\nConnection: close\r\n\r\n{\"error\":\"A solve is already running\"}";
                    let _ = socket.write_all(resp.as_bytes()).await;
                    return;
                }
                let body = req_str.split("\r\n\r\n").nth(1).unwrap_or("");
                if let Ok(req) = serde_json::from_str::<super::handlers::SolveRequest>(body) {
                    let traj_write = Arc::clone(&traj);
                    let base_write = Arc::clone(&base_list);
                    let solving_flag = Arc::clone(&solving);
                    let err_write = Arc::clone(&err_state);
                    tokio::spawn(async move {
                        *err_write.write().await = None;
                        match super::handlers::handle_solve_request(req).await {
                            Ok((new_epochs, new_bases)) => {
                                *traj_write.write().await = new_epochs;
                                *base_write.write().await = new_bases;
                            }
                            Err(e) => {
                                *err_write.write().await = Some(e);
                            }
                        }
                        solving_flag.store(false, std::sync::atomic::Ordering::SeqCst);
                    });
                    let resp = "HTTP/1.1 200 OK\r\nContent-Type: application/json\r\nContent-Length: 28\r\nConnection: close\r\n\r\n{\"status\":\"solving_started\"}";
                    let _ = socket.write_all(resp.as_bytes()).await;
                } else {
                    solving.store(false, std::sync::atomic::Ordering::SeqCst);
                    let resp = "HTTP/1.1 400 BAD REQUEST\r\nContent-Length: 0\r\nConnection: close\r\n\r\n";
                    let _ = socket.write_all(resp.as_bytes()).await;
                }
            } else if path.starts_with("/api/export") {
                let current = traj.read().await;
                let fmt_str = path.split("format=").nth(1).unwrap_or("pos");
                let fmt = ExportFormat::from_str(fmt_str).unwrap_or(ExportFormat::Pos);
                if let Ok(bytes) = generate_export_bytes(&current, fmt) {
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

pub(crate) fn load_trajectory_file(path: &Path) -> Result<Vec<SmoothedEpoch>, Box<dyn std::error::Error>> {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    if ext.eq_ignore_ascii_case("sbet") || ext.eq_ignore_ascii_case("out") {
        load_sbet_trajectory(path)
    } else {
        load_ascii_pos_trajectory(path)
    }
}

fn load_sbet_trajectory(path: &Path) -> Result<Vec<SmoothedEpoch>, Box<dyn std::error::Error>> {
    use gneiss_core::coords::llh_to_ecef;
    use gneiss_core::time::GpsTime;
    use gneiss_parsers::sbet::SbetRecord;
    use nalgebra::{Matrix3, Vector3};

    let mut file = std::fs::File::open(path)?;
    let mut epochs = Vec::new();
    while let Ok(rec) = SbetRecord::read_from(&mut file) {
        let pos = llh_to_ecef(Vector3::new(rec.latitude, rec.longitude, rec.altitude));
        epochs.push(SmoothedEpoch {
            time: GpsTime::new(2370, rec.time),
            position_ecef: pos,
            velocity_ecef: Some(Vector3::new(rec.x_vel, rec.y_vel, rec.z_vel)),
            attitude: None,
            cov_position: Matrix3::identity() * 0.0001,
            std_east: 0.01,
            std_north: 0.01,
            std_up: 0.02,
            separation_3d: 0.005,
            quality: 1,
            n_satellites: 12,
        });
    }
    Ok(epochs)
}

fn load_ascii_pos_trajectory(path: &Path) -> Result<Vec<SmoothedEpoch>, Box<dyn std::error::Error>> {
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

pub(crate) fn load_base_stations(paths: &[String]) -> Vec<BaseInfo> {
    use gneiss_core::coords::ecef_to_llh;
    let mut bases = Vec::new();
    for p in paths {
        let path = Path::new(p);
        let name = path.file_stem().unwrap_or_default().to_string_lossy().to_string();
        if let Ok((_, Some(pos))) = crate::ingest::UniversalObsReader::read_file(path) {
            let v = nalgebra::Vector3::new(pos[0], pos[1], pos[2]);
            let llh = ecef_to_llh(v);
            bases.push(BaseInfo {
                name,
                lat: llh.x.to_degrees(),
                lon: llh.y.to_degrees(),
                alt: llh.z,
            });
        }
    }
    bases
}
