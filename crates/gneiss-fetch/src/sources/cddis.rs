use crate::provider::{DataSource, FetchError};
use async_trait::async_trait;
use gneiss_core::coords::Coordinate;
use gneiss_core::time::GpsTime;
use std::path::{Path, PathBuf};

pub struct CddisProvider {
    pub auth_token: Option<String>,
}

impl CddisProvider {
    fn gps_time_to_utc(&self, time: GpsTime) -> chrono::NaiveDateTime {
        let gps_epoch = chrono::NaiveDate::from_ymd_opt(1980, 1, 6)
            .expect("GPS epoch 1980-01-06 is valid")
            .and_hms_opt(0, 0, 0)
            .expect("GPS epoch 00:00:00 is valid");
        let seconds = (time.week as i64 * 604800) + time.tow as i64;
        gps_epoch + chrono::Duration::seconds(seconds)
    }

    async fn download_and_decompress(
        &self,
        url: &str,
        filename: &str,
        out_dir: &Path,
    ) -> Result<PathBuf, FetchError> {
        let token = self
            .auth_token
            .as_ref()
            .ok_or_else(|| FetchError::Network("CDDIS requires an Earthdata token".into()))?;

        tracing::info!("Fetching CDDIS Product: {}", url);
        let client = reqwest::Client::new();
        let response = client.get(url).bearer_auth(token).send().await?;

        if !response.status().is_success() {
            return Err(FetchError::Network(format!(
                "HTTP Error {}: {}",
                response.status(),
                url
            )));
        }

        let out_file = out_dir.join(filename);
        let mut dest =
            std::fs::File::create(&out_file).map_err(|e| FetchError::Network(e.to_string()))?;

        use flate2::read::GzDecoder;
        let bytes = response.bytes().await?;
        let mut decoder = GzDecoder::new(&bytes[..]);
        std::io::copy(&mut decoder, &mut dest)
            .map_err(|e| FetchError::Decompression(e.to_string()))?;

        tracing::info!("Saved product to {}", out_file.display());
        Ok(out_file)
    }
}

#[async_trait]
impl DataSource for CddisProvider {
    fn name(&self) -> &str {
        "NASA_CDDIS"
    }

    async fn fetch_base_obs(
        &self,
        _location: Coordinate,
        _time: GpsTime,
        _out_dir: &Path,
    ) -> Result<PathBuf, FetchError> {
        Err(FetchError::NotFound(
            "Base obs fetch via CDDIS not yet implemented".into(),
        ))
    }

    async fn fetch_ephemeris(&self, time: GpsTime, out_dir: &Path) -> Result<PathBuf, FetchError> {
        let utc_time = self.gps_time_to_utc(time);
        let year = utc_time.format("%Y").to_string();
        let doy = utc_time.format("%j").to_string();
        let yy = utc_time.format("%y").to_string();

        let filename = format!("BRDC00IGS_R_{}{}0000_01D_MN.rnx", year, doy);
        let gz_filename = format!("{}.gz", filename);
        let url = format!(
            "https://cddis.nasa.gov/archive/gnss/data/daily/{}/{}/{}p/{}",
            year, doy, yy, gz_filename
        );

        self.download_and_decompress(&url, &filename, out_dir).await
    }

    async fn fetch_sp3(&self, time: GpsTime, out_dir: &Path) -> Result<PathBuf, FetchError> {
        let utc_time = self.gps_time_to_utc(time);
        let year = utc_time.format("%Y").to_string();
        let doy = utc_time.format("%j").to_string();

        let filename = format!("GFZ0MGXRAP_{}{}0000_01D_05M_ORB.SP3", year, doy);
        let gz_filename = format!("{}.gz", filename);
        let url = format!(
            "https://cddis.nasa.gov/archive/gnss/products/{}/{}",
            time.week, gz_filename
        );

        self.download_and_decompress(&url, &filename, out_dir).await
    }

    async fn fetch_clk(&self, time: GpsTime, out_dir: &Path) -> Result<PathBuf, FetchError> {
        let utc_time = self.gps_time_to_utc(time);
        let year = utc_time.format("%Y").to_string();
        let doy = utc_time.format("%j").to_string();

        let filename = format!("GFZ0MGXRAP_{}{}0000_01D_30S_CLK.CLK", year, doy);
        let gz_filename = format!("{}.gz", filename);
        let url = format!(
            "https://cddis.nasa.gov/archive/gnss/products/{}/{}",
            time.week, gz_filename
        );

        self.download_and_decompress(&url, &filename, out_dir).await
    }

    async fn fetch_bias(&self, time: GpsTime, out_dir: &Path) -> Result<PathBuf, FetchError> {
        let utc_time = self.gps_time_to_utc(time);
        let year = utc_time.format("%Y").to_string();
        let doy = utc_time.format("%j").to_string();

        let filename = format!("CAS0MGXRAP_{}{}0000_01D_01D_DCB.BIA", year, doy);
        let gz_filename = format!("{}.gz", filename);
        let url = format!(
            "https://cddis.nasa.gov/archive/gnss/products/bias/{}/{}",
            year, gz_filename
        );

        self.download_and_decompress(&url, &filename, out_dir).await
    }
}
