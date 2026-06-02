use std::io::BufRead;
use std::collections::HashMap;
use gneiss_core::time::GpsTime;

#[derive(Debug, Clone)]
pub struct Sp3Epoch {
    pub time: GpsTime,
    pub records: HashMap<String, Sp3Record>,
}

#[derive(Debug, Clone)]
pub struct Sp3Record {
    pub position: nalgebra::Vector3<f64>, // meters
    pub clock_offset: f64,      // seconds
}

pub fn parse_sp3<R: BufRead>(reader: R) -> Result<Vec<Sp3Epoch>, String> {
    let mut epochs = Vec::new();
    let mut current_epoch: Option<Sp3Epoch> = None;

    for line_result in reader.lines() {
        let line = line_result.map_err(|e| e.to_string())?;
        if line.is_empty() {
            continue;
        }

        if line.starts_with('*') {
            // Epoch header line: *  YYYY MM DD HH MM SS.sssssss
            let parts: Vec<&str> = line[1..].split_whitespace().collect();
            if parts.len() < 6 {
                continue;
            }
            if let (Ok(year), Ok(month), Ok(day), Ok(hour), Ok(minute), Ok(sec)) = (
                parts[0].parse::<i32>(),
                parts[1].parse::<u32>(),
                parts[2].parse::<u32>(),
                parts[3].parse::<u32>(),
                parts[4].parse::<u32>(),
                parts[5].parse::<f64>(),
            ) {
                if let Some(mut epoch) = current_epoch.take() {
                    epochs.push(epoch);
                }

                let time = GpsTime::from_calendar(year, month as i32, day as i32, hour as i32, minute as i32, sec);
                current_epoch = Some(Sp3Epoch {
                    time,
                    records: HashMap::new(),
                });
            }
        } else if line.starts_with('P') {
            // Position record: PG01  X Y Z Clock
            if let Some(epoch) = &mut current_epoch {
                if line.len() >= 60 {
                    let sat_id = line[1..4].trim().to_string();
                    let x_str = line[4..18].trim();
                    let y_str = line[18..32].trim();
                    let z_str = line[32..46].trim();
                    let clk_str = line[46..60].trim();

                    if let (Ok(x), Ok(y), Ok(z)) = (
                        x_str.parse::<f64>(),
                        y_str.parse::<f64>(),
                        z_str.parse::<f64>(),
                    ) {
                        let clock_offset = if !clk_str.is_empty() && clk_str != "999999.999999" {
                            clk_str.parse::<f64>().unwrap_or(core::f64::NAN) * 1e-6 // microseconds to seconds
                        } else {
                            core::f64::NAN
                        };

                        epoch.records.insert(sat_id, Sp3Record {
                            position: nalgebra::Vector3::new(x * 1000.0, y * 1000.0, z * 1000.0), // km to meters
                            clock_offset,
                        });
                    }
                }
            }
        } else if line.starts_with("EOF") {
            if let Some(epoch) = current_epoch.take() {
                epochs.push(epoch);
            }
            break;
        }
    }

    if let Some(epoch) = current_epoch {
        epochs.push(epoch);
    }

    Ok(epochs)
}
