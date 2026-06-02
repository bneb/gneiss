use crate::provider::{DataSource, FetchError};
use async_trait::async_trait;
use std::path::{Path, PathBuf};
use gneiss_core::coords::Coordinate;
use gneiss_core::time::GpsTime;

pub struct CddisProvider {
    pub auth_token: Option<String>,
}

#[async_trait]
impl DataSource for CddisProvider {
    fn name(&self) -> &str {
        "NASA_CDDIS"
    }

    async fn fetch_base_obs(&self, _location: Coordinate, _time: GpsTime, _out_dir: &Path) -> Result<PathBuf, FetchError> {
        Err(FetchError::NotFound("Base obs fetch via CDDIS not yet implemented".into()))
    }

    async fn fetch_ephemeris(&self, time: GpsTime, out_dir: &Path) -> Result<PathBuf, FetchError> {
        let token = self.auth_token.as_ref().ok_or_else(|| FetchError::Network("CDDIS requires an Earthdata token".into()))?;

        // Calculate UTC Date from GPS Time
        let gps_epoch = chrono::NaiveDate::from_ymd_opt(1980, 1, 6)
            .unwrap()
            .and_hms_opt(0, 0, 0)
            .unwrap();
        let seconds = (time.week as i64 * 604800) + time.tow as i64;
        let utc_time = gps_epoch + chrono::Duration::seconds(seconds);

        let year = utc_time.format("%Y").to_string(); // 2023
        let doy = utc_time.format("%j").to_string();  // 100
        let yy = utc_time.format("%y").to_string();   // 23

        let filename = format!("BRDC00IGS_R_{}{}0000_01D_MN.rnx", year, doy);
        let gz_filename = format!("{}.gz", filename);
        
        let url = format!("https://cddis.nasa.gov/archive/gnss/data/daily/{}/{}/{}p/{}", year, doy, yy, gz_filename);
        
        tracing::info!("Fetching CDDIS Ephemeris: {}", url);

        let client = reqwest::Client::new();
        let response = client
            .get(&url)
            .bearer_auth(token)
            .send()
            .await?;

        if !response.status().is_success() {
            return Err(FetchError::Network(format!("HTTP Error {}: {}", response.status(), url)));
        }

        let out_file = out_dir.join(&filename);
        let mut dest = std::fs::File::create(&out_file).map_err(|e| FetchError::Network(e.to_string()))?;

        // Extract .gz
        use flate2::read::GzDecoder;
        let bytes = response.bytes().await?;
        let mut decoder = GzDecoder::new(&bytes[..]);
        std::io::copy(&mut decoder, &mut dest).map_err(|e| FetchError::Decompression(e.to_string()))?;

        tracing::info!("Saved ephemeris to {}", out_file.display());
        Ok(out_file)
    }
}
