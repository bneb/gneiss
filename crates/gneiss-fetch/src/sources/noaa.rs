use crate::provider::{DataSource, FetchError};
use async_trait::async_trait;
use gneiss_core::coords::{Coordinate, Datum, Frame};
use gneiss_core::time::GpsTime;
use nalgebra::Vector3;
use std::path::{Path, PathBuf};

pub struct NoaaCorsProvider;

use serde::Deserialize;

#[derive(Deserialize, Debug)]
struct CorsStation {
    #[serde(rename = "corsId")]
    cors_id: String,
}

#[async_trait]
impl DataSource for NoaaCorsProvider {
    fn name(&self) -> &str {
        "NOAA_CORS"
    }

    async fn fetch_base_obs(
        &self,
        location: Coordinate,
        time: GpsTime,
        out_dir: &Path,
    ) -> Result<PathBuf, FetchError> {
        // Find nearest station
        let api_url = format!(
            "https://geodesy.noaa.gov/api/nde/ncors?x={:.0}&y={:.0}&z={:.0}",
            location.vector.x, location.vector.y, location.vector.z
        );
        tracing::info!("Querying NOAA CORS for nearest station: {}", api_url);

        let client = reqwest::Client::new();
        let stations: Vec<CorsStation> = client.get(&api_url).send().await?.json().await?;
        let station_id = stations
            .first()
            .map(|s| s.cors_id.to_lowercase())
            .ok_or_else(|| FetchError::NotFound("No nearby NOAA CORS stations found".into()))?;

        tracing::info!("Found nearest NOAA CORS station: {}", station_id);

        let gps_epoch = chrono::NaiveDate::from_ymd_opt(1980, 1, 6)
            .unwrap()
            .and_hms_opt(0, 0, 0)
            .unwrap();
        let seconds = (time.week as i64 * 604800) + time.tow as i64;
        let utc_time = gps_epoch + chrono::Duration::seconds(seconds);

        let year = utc_time.format("%Y").to_string(); // 2023
        let doy = utc_time.format("%j").to_string(); // 100
        let yy = utc_time.format("%y").to_string(); // 23

        // Download Hatanaka compressed RINEX observation
        let remote_filename = format!("{}{}0.{}d.gz", station_id, doy, yy);
        let unzipped_filename = format!("{}{}0.{}d", station_id, doy, yy);
        let url = format!(
            "https://geodesy.noaa.gov/corsdata/rinex/{}/{}/{}/{}",
            year, doy, station_id, remote_filename
        );

        tracing::info!("Fetching NOAA CORS Base Data: {}", url);

        let response = client.get(&url).send().await?;

        if !response.status().is_success() {
            return Err(FetchError::Network(format!(
                "HTTP Error {}: {}",
                response.status(),
                url
            )));
        }

        let unzipped_file = out_dir.join(&unzipped_filename);
        let mut dest = std::fs::File::create(&unzipped_file)
            .map_err(|e| FetchError::Network(e.to_string()))?;

        // Extract .gz while saving (leaves us with a .d hatanaka file)
        use flate2::read::GzDecoder;
        let bytes = response.bytes().await?;
        let mut decoder = GzDecoder::new(&bytes[..]);
        std::io::copy(&mut decoder, &mut dest)
            .map_err(|e| FetchError::Decompression(e.to_string()))?;

        tracing::info!(
            "Saved Hatanaka compressed base data to {}",
            unzipped_file.display()
        );

        tracing::info!("Decompressing Hatanaka to RINEX...");
        let final_file = crate::hatanaka::decompress(&unzipped_file, out_dir)?;

        Ok(final_file)
    }

    async fn fetch_ephemeris(
        &self,
        _time: GpsTime,
        _out_dir: &Path,
    ) -> Result<PathBuf, FetchError> {
        Err(FetchError::NotFound("Use CDDIS for ephemeris".into()))
    }
}

impl NoaaCorsProvider {
    pub async fn fetch_station_coordinate(
        &self,
        station_id: &str,
    ) -> Result<Coordinate, FetchError> {
        let url = format!(
            "https://geodesy.noaa.gov/corsdata/coord/coord_14/{}_14.coord.txt",
            station_id.to_lowercase()
        );
        tracing::info!("Fetching NOAA CORS Coordinate: {}", url);

        let client = reqwest::Client::new();
        let response = client.get(&url).send().await?;

        if !response.status().is_success() {
            return Err(FetchError::Network(format!(
                "Failed to fetch coordinate for station {}",
                station_id
            )));
        }

        let text = response.text().await?;

        let mut x = None;
        let mut y = None;
        let mut z = None;

        let mut in_itrf_block = false;
        for line in text.lines() {
            if line.contains("ITRF2014 POSITION") {
                in_itrf_block = true;
            } else if line.contains("NAD_83")
                || line.contains("VELOCITY")
                || line.contains("L1 Phase Center")
            {
                in_itrf_block = false;
            }

            if in_itrf_block {
                if line.contains("X =") {
                    if let Some(val_str) =
                        line.split("X =").nth(1).and_then(|s| s.split('m').next())
                    {
                        x = val_str.trim().parse::<f64>().ok();
                    }
                }
                if line.contains("Y =") {
                    if let Some(val_str) =
                        line.split("Y =").nth(1).and_then(|s| s.split('m').next())
                    {
                        y = val_str.trim().parse::<f64>().ok();
                    }
                }
                if line.contains("Z =") {
                    if let Some(val_str) =
                        line.split("Z =").nth(1).and_then(|s| s.split('m').next())
                    {
                        z = val_str.trim().parse::<f64>().ok();
                    }
                }
                if x.is_some() && y.is_some() && z.is_some() {
                    break;
                }
            }
        }

        if let (Some(x), Some(y), Some(z)) = (x, y, z) {
            Ok(Coordinate {
                vector: Vector3::new(x, y, z),
                datum: Datum::ITRF2014,
                frame: Frame::ECEF,
                epoch: GpsTime::new(0, 0.0),
            })
        } else {
            Err(FetchError::NotFound(format!(
                "Could not parse ITRF2014 XYZ from {}",
                url
            )))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_fetch_station_coordinate() {
        let provider = NoaaCorsProvider;
        let coord = provider.fetch_station_coordinate("P222").await.unwrap();
        // The ITRF2014 ARP coordinate for P222:
        // X = -2689640.291
        // Y = -4290437.370
        // Z = 3865050.934
        assert!((coord.vector.x - (-2689640.291)).abs() < 1e-3);
        assert!((coord.vector.y - (-4290437.370)).abs() < 1e-3);
        assert!((coord.vector.z - 3865050.934).abs() < 1e-3);
    }
}
