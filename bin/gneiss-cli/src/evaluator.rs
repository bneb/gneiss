use tracing::{info, error};
use gneiss_core::time::GpsTime;
use gneiss_core::coords::{ecef_to_llh, llh_to_ecef};

pub fn evaluate(solution: &str, truth: &str) -> Result<(), String> {
    info!("Evaluating solution against ground truth...");

    let truth_data = std::fs::read_to_string(truth).map_err(|e| format!("Failed to read truth file: {}", e))?;
    let mut truth_epochs = Vec::new();

    let is_llh = truth_data.lines().any(|l| l.contains("latitude(deg)"));

    if truth.ends_with(".csv") {
        let mut lines = truth_data.lines();
        if let Some(header) = lines.next() {
            let headers: Vec<&str> = header.split(',').map(|s| s.trim()).collect();
            
            let is_gsdc = headers.contains(&"millisSinceGpsEpoch");
            
            if is_gsdc {
                let time_idx = headers.iter().position(|h| *h == "millisSinceGpsEpoch").unwrap_or(2);
                let lat_idx = headers.iter().position(|h| *h == "latDeg").unwrap_or(3);
                let lon_idx = headers.iter().position(|h| *h == "lngDeg").unwrap_or(4);
                let hgt_idx = headers.iter().position(|h| *h == "heightAboveWgs84EllipsoidM").unwrap_or(5);
                
                for line in lines {
                    if line.trim().is_empty() { continue; }
                    let parts: Vec<&str> = line.split(',').collect();
                    if parts.len() > hgt_idx {
                        let millis: u64 = parts[time_idx].trim().parse().unwrap_or(0);
                        if millis > 0 {
                            let tow = (millis % 604_800_000) as f64 / 1000.0;
                            
                            let lat: f64 = parts[lat_idx].trim().parse().unwrap_or(0.0);
                            let lon: f64 = parts[lon_idx].trim().parse().unwrap_or(0.0);
                            let hgt: f64 = parts[hgt_idx].trim().parse().unwrap_or(0.0);
                            
                            let ecef = llh_to_ecef(nalgebra::Vector3::new(lat.to_radians(), lon.to_radians(), hgt));
                            truth_epochs.push((tow, ecef));
                        }
                    }
                }
            } else {
                for line in lines {
                    if line.trim().is_empty() { continue; }
                    let parts: Vec<&str> = line.split(',').collect();
                    if parts.len() >= 8 {
                        let tow: f64 = parts[0].trim().parse().unwrap_or(0.0);
                        let x: f64 = parts[5].trim().parse().unwrap_or(0.0);
                        let y: f64 = parts[6].trim().parse().unwrap_or(0.0);
                        let z: f64 = parts[7].trim().parse().unwrap_or(0.0);
                        truth_epochs.push((tow, nalgebra::Vector3::new(x, y, z)));
                    }
                }
            }
        }
    } else {
        for line in truth_data.lines() {
            if line.starts_with('%') || line.trim().is_empty() { continue; }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 5 { continue; }

            let date_parts: Vec<&str> = parts[0].split('/').collect();
            let time_parts: Vec<&str> = parts[1].split(':').collect();
            if date_parts.len() == 3 && time_parts.len() == 3 {
                let year = date_parts[0].parse().unwrap_or(0);
                let month = date_parts[1].parse().unwrap_or(0);
                let day = date_parts[2].parse().unwrap_or(0);
                let hour = time_parts[0].parse().unwrap_or(0);
                let min = time_parts[1].parse().unwrap_or(0);
                let sec = time_parts[2].parse().unwrap_or(0.0);

                let gps_time = GpsTime::from_calendar(year, month, day, hour, min, sec);
                let c1: f64 = parts[2].parse().unwrap_or(0.0);
                let c2: f64 = parts[3].parse().unwrap_or(0.0);
                let c3: f64 = parts[4].parse().unwrap_or(0.0);

                let ecef = if is_llh {
                    llh_to_ecef(nalgebra::Vector3::new(c1.to_radians(), c2.to_radians(), c3))
                } else {
                    nalgebra::Vector3::new(c1, c2, c3)
                };

                truth_epochs.push((gps_time.tow, ecef));
            }
        }
    }
    
    info!("Loaded {} epochs from ground truth.", truth_epochs.len());

    let sol_data = std::fs::read_to_string(solution).map_err(|e| format!("Failed to read solution file: {}", e))?;

    let mut horiz_errors = Vec::new();
    let mut vert_errors = Vec::new();
    let mut err_3d = Vec::new();
    let mut heading_errors = Vec::new();

    let sol_is_llh = sol_data.lines().any(|l| l.contains("latitude(deg) longitude(deg)"));

    let mut sol_epochs = Vec::new();
    for line in sol_data.lines() {
        if line.starts_with('%') || line.trim().is_empty() { continue; }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 5 { continue; }

        let mut tow: f64 = 0.0;
        let mut x: f64 = 0.0;
        let mut y: f64 = 0.0;
        let mut z: f64 = 0.0;

        if parts[0].contains('/') {
            let date_parts: Vec<&str> = parts[0].split('/').collect();
            let time_parts: Vec<&str> = parts[1].split(':').collect();
            if date_parts.len() == 3 && time_parts.len() == 3 {
                let year = date_parts[0].parse().unwrap_or(0);
                let month = date_parts[1].parse().unwrap_or(0);
                let day = date_parts[2].parse().unwrap_or(0);
                let hour = time_parts[0].parse().unwrap_or(0);
                let min = time_parts[1].parse().unwrap_or(0);
                let sec = time_parts[2].parse().unwrap_or(0.0);

                let gps_time = GpsTime::from_calendar(year, month, day, hour, min, sec);
                tow = gps_time.tow;
                x = parts[2].parse().unwrap_or(0.0);
                y = parts[3].parse().unwrap_or(0.0);
                z = parts[4].parse().unwrap_or(0.0);
            }
        } else {
            tow = parts[1].parse().unwrap_or(0.0);
            x = parts[2].parse().unwrap_or(0.0);
            y = parts[3].parse().unwrap_or(0.0);
            z = parts[4].parse().unwrap_or(0.0);
        }

        if tow > 0.0 {
            let ecef = if sol_is_llh {
                llh_to_ecef(nalgebra::Vector3::new(x.to_radians(), y.to_radians(), z))
            } else {
                nalgebra::Vector3::new(x, y, z)
            };
            sol_epochs.push((tow, ecef));
        }
    }

    let mut prev_true_ecef: Option<nalgebra::Vector3<f64>> = None;
    let mut prev_sol_ecef: Option<nalgebra::Vector3<f64>> = None;

    for (sol_tow, sol_ecef) in sol_epochs {
        let mut best_diff = f64::MAX;
        let mut best_true_ecef = None;

        for (true_tow, true_ecef) in &truth_epochs {
            let min_diff = (sol_tow - true_tow).abs();
            if min_diff < 0.15
                && min_diff < best_diff {
                    best_diff = min_diff;
                    best_true_ecef = Some(true_ecef);
                }
        }
        
        if let Some(true_ecef) = best_true_ecef {
            let diff_ecef = sol_ecef - true_ecef;

            let true_llh = ecef_to_llh(*true_ecef);
            let lat = true_llh.x;
            let lon = true_llh.y;
            let sin_lat = f64::sin(lat);
            let cos_lat = f64::cos(lat);
            let sin_lon = f64::sin(lon);
            let cos_lon = f64::cos(lon);

            let dx = diff_ecef.x;
            let dy = diff_ecef.y;
            let dz = diff_ecef.z;

            let e = -sin_lon * dx + cos_lon * dy;
            let n = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz;
            let u = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz;

            let h_err = f64::sqrt(e * e + n * n);
            let v_err = u.abs();
            let d3 = f64::sqrt(h_err * h_err + v_err * v_err);

            horiz_errors.push(h_err);
            vert_errors.push(v_err);
            err_3d.push(d3);
            
            if horiz_errors.len() % 100 == 0 {
                println!("Epoch {} (tow {:.1}): h_err={:.3}m v_err={:.3}m dx={:.3} dy={:.3} dz={:.3}", horiz_errors.len(), sol_tow, h_err, v_err, dx, dy, dz);
            }

            if let (Some(p_true_ecef), Some(p_sol_ecef)) = (prev_true_ecef, prev_sol_ecef) {
                let v_true_ecef = true_ecef - p_true_ecef;
                let v_sol_ecef = sol_ecef - p_sol_ecef;
                
                let v_true_e = -sin_lon * v_true_ecef.x + cos_lon * v_true_ecef.y;
                let v_true_n = -sin_lat * cos_lon * v_true_ecef.x - sin_lat * sin_lon * v_true_ecef.y + cos_lat * v_true_ecef.z;
                
                let v_sol_e = -sin_lon * v_sol_ecef.x + cos_lon * v_sol_ecef.y;
                let v_sol_n = -sin_lat * cos_lon * v_sol_ecef.x - sin_lat * sin_lon * v_sol_ecef.y + cos_lat * v_sol_ecef.z;

                let true_heading = f64::atan2(v_true_e, v_true_n).to_degrees();
                let sol_heading = f64::atan2(v_sol_e, v_sol_n).to_degrees();

                let mut h_err_deg = (sol_heading - true_heading).abs();
                if h_err_deg > 180.0 { h_err_deg = 360.0 - h_err_deg; }

                if f64::sqrt(v_true_e * v_true_e + v_true_n * v_true_n) > 0.5 {
                    heading_errors.push(h_err_deg);
                }
            }

            prev_true_ecef = Some(*true_ecef);
            prev_sol_ecef = Some(sol_ecef);
        }
    }

    if err_3d.is_empty() {
        error!("No matching epochs found between solution and truth.");
        return Ok(());
    }

    horiz_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    vert_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    err_3d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    heading_errors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    if heading_errors.is_empty() { heading_errors.push(0.0); }

    let count = err_3d.len() as f64;
    let h_count = heading_errors.len() as f64;

    let p25 = (count * 0.25) as usize;
    let p50 = (count * 0.50) as usize;
    let p75 = (count * 0.75) as usize;
    let p90 = (count * 0.90) as usize;
    let p95 = (count * 0.95) as usize;
    let p99 = (count * 0.99) as usize;

    let hp25 = (h_count * 0.25) as usize;
    let hp50 = (h_count * 0.50) as usize;
    let hp75 = (h_count * 0.75) as usize;
    let hp90 = (h_count * 0.90) as usize;
    let hp95 = (h_count * 0.95) as usize;
    let hp99 = (h_count * 0.99) as usize;

    println!("==========================================================================================");
    println!("                       AEROSPACE METRICS: ERROR CDFs (meters / deg)                       ");
    println!("==========================================================================================");
    println!("Evaluated {} epochs. Moving epochs: {}", count, h_count);
    println!("------------------------------------------------------------------------------------------");
    println!("| Metric  |   25th %%  |   50th %%  |   75th %%  |   90th %%  |   95th %%  |   99th %%  |");
    println!("|---------|------------|------------|------------|------------|------------|------------|");
    println!("| Horiz   | {:>8.3} m | {:>8.3} m | {:>8.3} m | {:>8.3} m | {:>8.3} m | {:>8.3} m |", horiz_errors[p25], horiz_errors[p50], horiz_errors[p75], horiz_errors[p90], horiz_errors[p95], horiz_errors[p99]);
    println!("| Vert    | {:>8.3} m | {:>8.3} m | {:>8.3} m | {:>8.3} m | {:>8.3} m | {:>8.3} m |", vert_errors[p25], vert_errors[p50], vert_errors[p75], vert_errors[p90], vert_errors[p95], vert_errors[p99]);
    println!("| 3D      | {:>8.3} m | {:>8.3} m | {:>8.3} m | {:>8.3} m | {:>8.3} m | {:>8.3} m |", err_3d[p25], err_3d[p50], err_3d[p75], err_3d[p90], err_3d[p95], err_3d[p99]);
    println!("| Heading | {:>8.3}°  | {:>8.3}°  | {:>8.3}°  | {:>8.3}°  | {:>8.3}°  | {:>8.3}°  |", heading_errors[hp25], heading_errors[hp50], heading_errors[hp75], heading_errors[hp90], heading_errors[hp95], heading_errors[hp99]);
    println!("==================================================");

    Ok(())
}