//! Multi-format trajectory exporter.
//!
//! Supports standard and extended POS, surveyor geodetic CSV/LLH,
//! Google Earth KML tracks, and GeoJSON with per-epoch uncertainties
//! and optional orthometric height conversion.

use std::fs::File;
use std::io::Write;
use std::path::Path;

use nalgebra::Vector3;

use gneiss_core::coords::ecef_to_llh;
use gneiss_geodesy::geoid::GeoidGrid;
use gneiss_rtk::post_process::SmoothedEpoch;

/// Output export format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportFormat {
    Pos,
    Csv,
    Kml,
    Json,
    Sbet,
}

impl ExportFormat {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "pos" => Some(Self::Pos),
            "csv" | "llh" => Some(Self::Csv),
            "kml" => Some(Self::Kml),
            "json" | "geojson" => Some(Self::Json),
            "sbet" | "out" => Some(Self::Sbet),
            _ => None,
        }
    }

    pub fn default_extension(&self) -> &'static str {
        match self {
            Self::Pos => "pos",
            Self::Csv => "csv",
            Self::Kml => "kml",
            Self::Json => "json",
            Self::Sbet => "sbet",
        }
    }
}

/// Formats a smoothed trajectory into the specified export format.
pub fn export_trajectory(
    trajectory: &[SmoothedEpoch],
    format: ExportFormat,
    output_path: &Path,
    geoid: Option<&GeoidGrid>,
) -> Result<(), std::io::Error> {
    match format {
        ExportFormat::Pos => {
            let mut file = File::create(output_path)?;
            write_pos(&mut file, trajectory, geoid)
        }
        ExportFormat::Csv => {
            let mut file = File::create(output_path)?;
            write_csv(&mut file, trajectory, geoid)
        }
        ExportFormat::Kml => {
            let mut file = File::create(output_path)?;
            write_kml(&mut file, trajectory)
        }
        ExportFormat::Json => {
            let mut file = File::create(output_path)?;
            write_json(&mut file, trajectory, geoid)
        }
        ExportFormat::Sbet => {
            let rms_path = output_path.with_extension("sbet.rms");
            gneiss_rtk::post_process::export_sbet_trajectory(trajectory, output_path, Some(&rms_path))
        }
    }
}

fn write_pos(
    w: &mut dyn Write,
    trajectory: &[SmoothedEpoch],
    geoid: Option<&GeoidGrid>,
) -> Result<(), std::io::Error> {
    writeln!(w, "% Program   : Gneiss PPK Engine")?;
    writeln!(w, "% Version   : 0.1.0")?;
    writeln!(w, "% Coordinate: WGS84 / IGS20")?;
    writeln!(w, "% GPST-Week TOW(s) Latitude(deg) Longitude(deg) H_ellips(m) H_ortho(m) Q Nsat sd_e(m) sd_n(m) sd_u(m) sep_3d(m) x-ecef(m) y-ecef(m) z-ecef(m)")?;
    for ep in trajectory {
        let (lat_deg, lon_deg, h_ellips, h_ortho) = compute_heights(ep.position_ecef, geoid);
        writeln!(
            w,
            "{} {:.3} {:.9} {:.9} {:.4} {:.4} {} {:2} {:.4} {:.4} {:.4} {:.4} {:.4} {:.4} {:.4}",
            ep.time.week, ep.time.tow, lat_deg, lon_deg, h_ellips, h_ortho,
            ep.quality, ep.n_satellites, ep.std_east, ep.std_north, ep.std_up,
            ep.separation_3d, ep.position_ecef.x, ep.position_ecef.y, ep.position_ecef.z,
        )?;
    }
    Ok(())
}

fn write_csv(
    w: &mut dyn Write,
    trajectory: &[SmoothedEpoch],
    geoid: Option<&GeoidGrid>,
) -> Result<(), std::io::Error> {
    writeln!(w, "week,tow_s,lat_deg,lon_deg,h_ellips_m,h_ortho_m,undulation_m,q,n_sat,sd_e_m,sd_n_m,sd_u_m,sep_3d_m,vx_m_s,vy_m_s,vz_m_s")?;
    for ep in trajectory {
        let (lat_deg, lon_deg, h_ellips, h_ortho) = compute_heights(ep.position_ecef, geoid);
        let undulation = h_ellips - h_ortho;
        let vx = ep.velocity_ecef.map_or(0.0, |v| v.x);
        let vy = ep.velocity_ecef.map_or(0.0, |v| v.y);
        let vz = ep.velocity_ecef.map_or(0.0, |v| v.z);
        writeln!(
            w,
            "{},{:.3},{:.9},{:.9},{:.4},{:.4},{:.4},{},{},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}",
            ep.time.week, ep.time.tow, lat_deg, lon_deg, h_ellips, h_ortho, undulation,
            ep.quality, ep.n_satellites, ep.std_east, ep.std_north, ep.std_up,
            ep.separation_3d, vx, vy, vz,
        )?;
    }
    Ok(())
}

fn write_kml(w: &mut dyn Write, trajectory: &[SmoothedEpoch]) -> Result<(), std::io::Error> {
    writeln!(w, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")?;
    writeln!(w, "<kml xmlns=\"http://www.opengis.net/kml/2.2\">")?;
    writeln!(w, "  <Document>")?;
    writeln!(w, "    <name>Gneiss Trajectory</name>")?;
    write_kml_styles(w)?;
    write_kml_linestring(w, trajectory)?;
    writeln!(w, "  </Document>")?;
    writeln!(w, "</kml>")?;
    Ok(())
}

fn write_kml_styles(w: &mut dyn Write) -> Result<(), std::io::Error> {
    writeln!(w, "    <Style id=\"trackStyle\">")?;
    writeln!(w, "      <LineStyle><color>ff00ff00</color><width>3</width></LineStyle>")?;
    writeln!(w, "    </Style>")?;
    Ok(())
}

fn write_kml_linestring(w: &mut dyn Write, trajectory: &[SmoothedEpoch]) -> Result<(), std::io::Error> {
    writeln!(w, "    <Placemark>")?;
    writeln!(w, "      <name>Path</name>")?;
    writeln!(w, "      <styleUrl>#trackStyle</styleUrl>")?;
    writeln!(w, "      <LineString>")?;
    writeln!(w, "        <altitudeMode>absolute</altitudeMode>")?;
    writeln!(w, "        <coordinates>")?;
    for ep in trajectory {
        let llh = ecef_to_llh(ep.position_ecef);
        let lat_deg = llh.x.to_degrees();
        let lon_deg = llh.y.to_degrees();
        writeln!(w, "          {:.9},{:.9},{:.3}", lon_deg, lat_deg, llh.z)?;
    }
    writeln!(w, "        </coordinates>")?;
    writeln!(w, "      </LineString>")?;
    writeln!(w, "    </Placemark>")?;
    Ok(())
}

fn write_json(
    w: &mut dyn Write,
    trajectory: &[SmoothedEpoch],
    geoid: Option<&GeoidGrid>,
) -> Result<(), std::io::Error> {
    writeln!(w, "{{")?;
    writeln!(w, "  \"type\": \"FeatureCollection\",")?;
    writeln!(w, "  \"features\": [")?;
    for (i, ep) in trajectory.iter().enumerate() {
        let (lat_deg, lon_deg, h_ellips, h_ortho) = compute_heights(ep.position_ecef, geoid);
        let comma = if i + 1 < trajectory.len() { "," } else { "" };
        writeln!(w, "    {{")?;
        writeln!(w, "      \"type\": \"Feature\",")?;
        writeln!(w, "      \"geometry\": {{ \"type\": \"Point\", \"coordinates\": [{:.9}, {:.9}, {:.4}] }},", lon_deg, lat_deg, h_ortho)?;
        writeln!(
            w,
            "      \"properties\": {{ \"week\": {}, \"tow\": {:.3}, \"h_ellips\": {:.4}, \"q\": {}, \"nsat\": {}, \"sd_e\": {:.4}, \"sd_n\": {:.4}, \"sd_u\": {:.4}, \"sep_3d\": {:.4} }}",
            ep.time.week, ep.time.tow, h_ellips, ep.quality, ep.n_satellites, ep.std_east, ep.std_north, ep.std_up, ep.separation_3d
        )?;
        writeln!(w, "    }}{}", comma)?;
    }
    writeln!(w, "  ]")?;
    writeln!(w, "}}")?;
    Ok(())
}

fn compute_heights(pos: Vector3<f64>, geoid: Option<&GeoidGrid>) -> (f64, f64, f64, f64) {
    let llh = ecef_to_llh(pos);
    let lat_deg = llh.x.to_degrees();
    let lon_deg = llh.y.to_degrees();
    let h_ellips = llh.z;
    let h_ortho = if let Some(g) = geoid {
        g.ellipsoidal_to_orthometric(llh).map_or(h_ellips, |o| o.z)
    } else {
        h_ellips
    };
    (lat_deg, lon_deg, h_ellips, h_ortho)
}

/// Exports photogrammetry camera event records to CSV format.
#[cfg(test)]
pub fn export_camera_events(
    events: &[gneiss_rtk::events::CameraEventRecord],
    output_path: &Path,
) -> Result<(), std::io::Error> {
    let mut file = File::create(output_path)?;
    writeln!(
        file,
        "# Photogrammetry Camera Center Positions - Gneiss PPK"
    )?;
    writeln!(
        file,
        "event_id,time_gpst,lat_deg,lon_deg,height_m,ecef_x,ecef_y,ecef_z,sd_e,sd_n,sd_u"
    )?;
    for ev in events {
        writeln!(
            file,
            "{},{:.6},{:.9},{:.9},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}",
            ev.event_id,
            ev.time_gpst_s,
            ev.lat_deg,
            ev.lon_deg,
            ev.height_m,
            ev.pos_ecef[0],
            ev.pos_ecef[1],
            ev.pos_ecef[2],
            ev.std_enu[0],
            ev.std_enu[1],
            ev.std_enu[2]
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::time::GpsTime;
    use nalgebra::Matrix3;
    use tempfile::NamedTempFile;

    fn sample_epoch() -> SmoothedEpoch {
        SmoothedEpoch {
            time: GpsTime { week: 2105, tow: 345600.0 },
            position_ecef: Vector3::new(-2688181.50, -4265663.45, 3893784.80),
            velocity_ecef: Some(Vector3::new(0.01, -0.02, 0.005)),
            attitude: None,
            cov_position: Matrix3::identity() * 0.0004,
            std_east: 0.012,
            std_north: 0.015,
            std_up: 0.025,
            separation_3d: 0.018,
            quality: 1,
            n_satellites: 14,
        }
    }

    #[test]
    fn test_format_parsing() {
        assert_eq!(ExportFormat::from_str("pos"), Some(ExportFormat::Pos));
        assert_eq!(ExportFormat::from_str("CSV"), Some(ExportFormat::Csv));
        assert_eq!(ExportFormat::from_str("llh"), Some(ExportFormat::Csv));
        assert_eq!(ExportFormat::from_str("kml"), Some(ExportFormat::Kml));
        assert_eq!(ExportFormat::from_str("geojson"), Some(ExportFormat::Json));
        assert_eq!(ExportFormat::from_str("invalid"), None);
    }

    #[test]
    fn test_export_all_formats() {
        let traj = vec![sample_epoch()];
        for fmt in [ExportFormat::Pos, ExportFormat::Csv, ExportFormat::Kml, ExportFormat::Json] {
            let tmp = NamedTempFile::new().expect("tempfile");
            let res = export_trajectory(&traj, fmt, tmp.path(), None);
            assert!(res.is_ok());
            let content = std::fs::read_to_string(tmp.path()).expect("read");
            assert!(!content.is_empty());
        }
    }

    #[test]
    fn test_export_with_geoid() {
        let traj = vec![sample_epoch()];
        let grid = GeoidGrid::new(30.0, 40.0, -130.0, -110.0, 1.0, 1.0, 11, 21, vec![32.0; 231])
            .expect("valid grid");
        let tmp = NamedTempFile::new().expect("tempfile");
        let res = export_trajectory(&traj, ExportFormat::Csv, tmp.path(), Some(&grid));
        assert!(res.is_ok());
        let content = std::fs::read_to_string(tmp.path()).expect("read");
        assert!(content.contains("345600.000"));
    }

    #[test]
    fn test_export_camera_events() {
        let tmp = NamedTempFile::new().expect("tempfile");
        let record = gneiss_rtk::events::CameraEventRecord {
            event_id: 1,
            time_gpst_s: 345600.0,
            pos_ecef: [-1283433.0, -4713073.0, 4090105.0],
            lat_deg: 35.0,
            lon_deg: -115.0,
            height_m: 500.0,
            std_enu: [0.01, 0.01, 0.02],
        };
        let res = export_camera_events(&[record], tmp.path());
        assert!(res.is_ok());
        let content = std::fs::read_to_string(tmp.path()).expect("read");
        assert!(content.contains("Camera Center Positions"));
    }
}
