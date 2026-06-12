use tracing::info;
use gneiss_core::coords::{Coordinate, Datum, Frame};
use std::io::{BufRead, BufReader};

pub async fn run_fetch(rover_obs: String, source: String, out_dir: String) -> Result<(), Box<dyn std::error::Error>> {
    info!("Fetching data for {} using source {}", rover_obs, source);
    
    let (time, coord) = parse_rover_info(&rover_obs)?;

    let out_path = std::path::Path::new(&out_dir);
    std::fs::create_dir_all(out_path)?;

    let token = std::env::var("EARTHDATA_TOKEN").ok();
    let source_lower = source.to_lowercase();
    
    if source_lower == "cddis" || source_lower == "all" {
        use gneiss_fetch::sources::cddis::CddisProvider;
        use gneiss_fetch::provider::DataSource;
        let cddis = CddisProvider { auth_token: token.clone() };
        cddis.fetch_ephemeris(time, out_path).await?;
    }
    
    if source_lower == "noaa" || source_lower == "all" {
        use gneiss_fetch::sources::noaa::NoaaCorsProvider;
        use gneiss_fetch::provider::DataSource;
        let noaa = NoaaCorsProvider;
        noaa.fetch_base_obs(coord, time, out_path).await?;
    }

    info!("Fetch completed successfully.");
    Ok(())
}

fn parse_rover_info(rover_obs: &str) -> Result<(gneiss_core::time::GpsTime, Coordinate), Box<dyn std::error::Error>> {
    let file = std::fs::File::open(rover_obs)?;
    let reader = BufReader::new(file);
    
    let mut approx_pos = None;
    let mut header_and_first_epoch = String::new();
    let mut in_header = true;
    let mut post_header_lines = 0;

    for line_res in reader.lines() {
        let line = line_res?;
        header_and_first_epoch.push_str(&line);
        header_and_first_epoch.push('\n');

        if in_header {
            if line.contains("APPROX POSITION XYZ") {
                let parts: Vec<&str> = line[0..60].split_whitespace().collect();
                if parts.len() >= 3 {
                    let ecef = nalgebra::Vector3::new(parts[0].parse()?, parts[1].parse()?, parts[2].parse()?);
                    approx_pos = Some(ecef);
                }
            } else if line.contains("END OF HEADER") {
                in_header = false;
            }
        } else {
            post_header_lines += 1;
            if post_header_lines > 100 {
                break;
            }
        }
    }

    let epochs = gneiss_parsers::rinex::parse_rinex_obs(std::io::Cursor::new(header_and_first_epoch))
        .map_err(|e| format!("Failed to parse first epoch from RINEX: {}", e))?;
    let first_epoch = epochs.first().ok_or("No epochs found in rover obs file")?;
    let time = first_epoch.time;
    
    if let Some(ecef) = approx_pos {
        let coord = Coordinate::new(ecef, Datum::WGS84, Frame::ECEF, time);
        Ok((time, coord))
    } else {
        Err("Could not determine approximate position from RINEX header".into())
    }
}
