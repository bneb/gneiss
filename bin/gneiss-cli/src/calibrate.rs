//! Interactive Local Site Calibration CLI.
//!
//! Fits a 7-parameter (4-parameter horizontal + 3-parameter vertical inclined plane)
//! transformation from paired GNSS projection coordinates to local Ground Control Points.

use std::fs::File;
use std::io::Write;
use std::path::Path;

use nalgebra::Vector3;
use gneiss_geodesy::site_calibration::SiteCalibration;

pub struct CalibrateArgs {
    pub input_csv: String,
    pub output_json: Option<String>,
}

pub fn run_calibrate(args: CalibrateArgs) -> Result<(), Box<dyn std::error::Error>> {
    let content = std::fs::read_to_string(&args.input_csv)?;
    let mut pairs = Vec::new();
    let mut point_names = Vec::new();

    for (idx, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with('%') {
            continue;
        }

        let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
        if parts.len() >= 6 {
            // Check if first column is point name
            let (name, gnss_e, gnss_n, gnss_h, gnd_e, gnd_n, gnd_h) = if parts.len() >= 7 {
                let name = parts[0].to_string();
                let g_e: f64 = parts[1].parse().map_err(|_| "invalid GNSS Easting")?;
                let g_n: f64 = parts[2].parse().map_err(|_| "invalid GNSS Northing")?;
                let g_h: f64 = parts[3].parse().map_err(|_| "invalid GNSS Height")?;
                let t_e: f64 = parts[4].parse().map_err(|_| "invalid Ground Easting")?;
                let t_n: f64 = parts[5].parse().map_err(|_| "invalid Ground Northing")?;
                let t_h: f64 = parts[6].parse().map_err(|_| "invalid Ground Height")?;
                (name, g_e, g_n, g_h, t_e, t_n, t_h)
            } else {
                let name = format!("GCP_{:02}", idx + 1);
                let g_e: f64 = parts[0].parse().map_err(|_| "invalid GNSS Easting")?;
                let g_n: f64 = parts[1].parse().map_err(|_| "invalid GNSS Northing")?;
                let g_h: f64 = parts[2].parse().map_err(|_| "invalid GNSS Height")?;
                let t_e: f64 = parts[3].parse().map_err(|_| "invalid Ground Easting")?;
                let t_n: f64 = parts[4].parse().map_err(|_| "invalid Ground Northing")?;
                let t_h: f64 = parts[5].parse().map_err(|_| "invalid Ground Height")?;
                (name, g_e, g_n, g_h, t_e, t_n, t_h)
            };

            point_names.push(name);
            pairs.push((Vector3::new(gnss_e, gnss_n, gnss_h), Vector3::new(gnd_e, gnd_n, gnd_h)));
        }
    }

    if pairs.len() < 3 {
        return Err(format!("Site calibration requires at least 3 valid point pairs, found {}", pairs.len()).into());
    }

    let calib = SiteCalibration::fit(&pairs)
        .ok_or("Collinear or degenerate geometry: failed to compute site calibration")?;

    let residuals = calib.compute_residuals(&pairs);

    println!("\n================================================================================");
    println!("                      GNEISS LOCAL SITE CALIBRATION REPORT                      ");
    println!("================================================================================");
    println!("Control Points Ingested : {}", pairs.len());
    println!("Source Grid Centroid (E): {:.4} m", calib.e0);
    println!("Source Grid Centroid (N): {:.4} m", calib.n0);
    println!("--------------------------------------------------------------------------------");
    println!("CALIBRATION PARAMETERS:");
    println!("  Translation East (dX)  : {:+10.4} m", calib.dx);
    println!("  Translation North (dY) : {:+10.4} m", calib.dy);
    println!("  Rotation Angle         : {:+10.6} deg ({:+8.2} arcsec)", calib.rotation_rad.to_degrees(), calib.rotation_rad.to_degrees() * 3600.0);
    let ppm = (calib.scale - 1.0) * 1e6;
    println!("  Horizontal Scale Factor: {:10.8} (Scale Error: {:+7.2} ppm)", calib.scale, ppm);
    println!("  Vertical Constant (dZ0): {:+10.4} m", calib.dz0);
    println!("  Vertical Slope East    : {:+10.4} ppm", calib.slope_east * 1e6);
    println!("  Vertical Slope North   : {:+10.4} ppm", calib.slope_north * 1e6);
    println!("--------------------------------------------------------------------------------");
    println!("CONTROL POINT RESIDUALS:");
    println!("  {:<10} {:>10} {:>10} {:>10} {:>10} {:>10}", "Point", "GNSS_E", "GNSS_N", "dE (mm)", "dN (mm)", "dH (mm)");

    let mut sum_h2 = 0.0;
    let mut sum_v2 = 0.0;
    let mut max_h_err = 0.0;
    let mut max_v_err = 0.0;

    for (i, (h_err, v_err)) in residuals.iter().enumerate() {
        let (src, tgt) = pairs[i];
        let transformed = calib.transform(src);
        let de_mm = (transformed.x - tgt.x) * 1000.0;
        let dn_mm = (transformed.y - tgt.y) * 1000.0;
        let dh_mm = v_err * 1000.0;

        sum_h2 += h_err * h_err;
        sum_v2 += v_err * v_err;
        if *h_err > max_h_err { max_h_err = *h_err; }
        if v_err.abs() > max_v_err { max_v_err = v_err.abs(); }

        println!("  {:<10} {:10.2} {:10.2} {:+10.2} {:+10.2} {:+10.2}",
            point_names[i], src.x, src.y, de_mm, dn_mm, dh_mm);
    }

    let h_rms_mm = libm::sqrt(sum_h2 / (pairs.len() as f64)) * 1000.0;
    let v_rms_mm = libm::sqrt(sum_v2 / (pairs.len() as f64)) * 1000.0;

    println!("--------------------------------------------------------------------------------");
    println!("ACCURACY SUMMARY:");
    println!("  Horizontal RMS Error   : {:6.2} mm  (Max: {:6.2} mm)", h_rms_mm, max_h_err * 1000.0);
    println!("  Vertical RMS Error     : {:6.2} mm  (Max: {:6.2} mm)", v_rms_mm, max_v_err * 1000.0);
    println!("================================================================================\n");

    if let Some(out_path) = args.output_json {
        let json_str = serde_json::to_string_pretty(&calib_to_json(&calib))?;
        let mut out_file = File::create(Path::new(&out_path))?;
        out_file.write_all(json_str.as_bytes())?;
        println!("Calibration saved to: {}", out_path);
    }

    Ok(())
}

fn calib_to_json(c: &SiteCalibration) -> serde_json::Value {
    serde_json::json!({
        "dx_m": c.dx,
        "dy_m": c.dy,
        "rotation_rad": c.rotation_rad,
        "rotation_deg": c.rotation_rad.to_degrees(),
        "scale": c.scale,
        "scale_ppm": (c.scale - 1.0) * 1e6,
        "e0_m": c.e0,
        "n0_m": c.n0,
        "dz0_m": c.dz0,
        "slope_east_ppm": c.slope_east * 1e6,
        "slope_north_ppm": c.slope_north * 1e6,
    })
}
