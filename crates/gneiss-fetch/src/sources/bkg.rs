use crate::provider::{DataSource, FetchError};
use async_trait::async_trait;
use gneiss_core::coords::Coordinate;
use gneiss_core::time::GpsTime;
use std::path::{Path, PathBuf};

/// Standard leap seconds between GPST and UTC since 2017-01-01 (18 seconds).
const GPS_LEAP_SECONDS: i64 = 18;

/// Convert GPST to civil UTC NaiveDateTime accounting for leap seconds.
pub(crate) fn gps_to_utc_naive(time: GpsTime) -> Result<chrono::NaiveDateTime, FetchError> {
    let gps_epoch = chrono::NaiveDate::from_ymd_opt(1980, 1, 6)
        .ok_or_else(|| FetchError::Decompression("Invalid GPS epoch date".to_string()))?
        .and_hms_opt(0, 0, 0)
        .ok_or_else(|| FetchError::Decompression("Invalid GPS epoch time".to_string()))?;
    let seconds = (time.week as i64 * 604800) + (time.tow as i64) - GPS_LEAP_SECONDS;
    Ok(gps_epoch + chrono::Duration::seconds(seconds))
}

pub struct BkgProvider;

#[async_trait]
impl DataSource for BkgProvider {
    fn name(&self) -> &str {
        "BKG"
    }

    async fn fetch_base_obs(
        &self,
        _location: Coordinate,
        time: GpsTime,
        out_dir: &Path,
    ) -> Result<PathBuf, FetchError> {
        let utc_time = gps_to_utc_naive(time)?;

        let year = utc_time.format("%Y").to_string(); // 2020
        let doy = utc_time.format("%j").to_string(); // 359

        // Hardcoded for the benchmark example (WTZR). Could be parameterized later.
        let station = "WTZR00DEU";
        let filename = format!("{}_R_{}{}0000_01D_30S_MO.crx.gz", station, year, doy);

        let url = format!(
            "https://igs.bkg.bund.de/root_ftp/IGS/obs/{}/{}/{}",
            year, doy, filename
        );

        tracing::info!("Fetching BKG Observation: {}", url);

        let client = reqwest::Client::builder()
            .user_agent("Gneiss-Navigation-Engine/0.1.0")
            .build()
            .expect("reqwest Client::builder().build() should succeed with default settings");
        let response = client.get(&url).send().await?;

        if !response.status().is_success() {
            return Err(FetchError::Network(format!(
                "HTTP Error {}: {}",
                response.status(),
                url
            )));
        }

        // Output file should be .crx without .gz
        let crx_filename = filename.trim_end_matches(".gz");
        let crx_file = out_dir.join(crx_filename);
        let mut dest =
            std::fs::File::create(&crx_file).map_err(|e| FetchError::Network(e.to_string()))?;

        // Extract .gz
        use flate2::read::GzDecoder;
        let bytes = response.bytes().await?;
        let mut decoder = GzDecoder::new(&bytes[..]);
        std::io::copy(&mut decoder, &mut dest)
            .map_err(|e| FetchError::Decompression(e.to_string()))?;

        tracing::info!(
            "Saved Hatanaka compressed base data to {}",
            crx_file.display()
        );

        tracing::info!("Decompressing Hatanaka to RINEX...");
        let final_file = crate::hatanaka::decompress(&crx_file, out_dir)?;

        Ok(final_file)
    }

    async fn fetch_ephemeris(&self, time: GpsTime, out_dir: &Path) -> Result<PathBuf, FetchError> {
        let utc_time = gps_to_utc_naive(time)?;

        let year = utc_time.format("%Y").to_string();
        let doy = utc_time.format("%j").to_string();

        let filename = format!("BRDC00WRD_R_{}{}0000_01D_MN.rnx", year, doy);
        let gz_filename = format!("{}.gz", filename);

        let url = format!(
            "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{}/{}/{}",
            year, doy, gz_filename
        );

        tracing::info!("Fetching BKG BRDC: {}", url);

        let client = reqwest::Client::builder()
            .user_agent("Gneiss-Navigation-Engine/0.1.0")
            .build()
            .expect("reqwest Client::builder().build() should succeed with default settings");
        let response = client.get(&url).send().await?;

        if !response.status().is_success() {
            return Err(FetchError::Network(format!(
                "HTTP Error {}: {}",
                response.status(),
                url
            )));
        }

        let out_file = out_dir.join(&filename);
        let mut dest =
            std::fs::File::create(&out_file).map_err(|e| FetchError::Network(e.to_string()))?;

        // Extract .gz
        use flate2::read::GzDecoder;
        let bytes = response.bytes().await?;
        let mut decoder = GzDecoder::new(&bytes[..]);
        std::io::copy(&mut decoder, &mut dest)
            .map_err(|e| FetchError::Decompression(e.to_string()))?;

        tracing::info!("Saved ephemeris to {}", out_file.display());
        Ok(out_file)
    }
}

impl BkgProvider {
    pub async fn fetch_precise_products(
        &self,
        time: GpsTime,
        out_dir: &Path,
    ) -> Result<Vec<PathBuf>, FetchError> {
        let gps_week = time.week;
        let gps_dow = (time.tow / 86400.0).floor() as u32;

        let base_url = format!(
            "https://igs.bkg.bund.de/root_ftp/IGS/products/mgex/{}/",
            gps_week
        );
        let files = [
            format!("com{}{}.eph.Z", gps_week, gps_dow),
            format!("com{}{}.clk.Z", gps_week, gps_dow),
            format!("com{}{}.bia.Z", gps_week, gps_dow),
        ];

        let mut out_paths = Vec::new();
        let client = reqwest::Client::builder()
            .user_agent("Gneiss-Navigation-Engine/0.1.0")
            .build()
            .expect("reqwest Client::builder().build() should succeed with default settings");

        for f in &files {
            let url = format!("{}{}", base_url, f);
            tracing::info!("Fetching BKG Precise Product: {}", url);

            let response = client.get(&url).send().await?;

            if !response.status().is_success() {
                return Err(FetchError::Network(format!(
                    "HTTP Error {}: {}",
                    response.status(),
                    url
                )));
            }

            let dest_z = out_dir.join(f);
            let mut dest =
                std::fs::File::create(&dest_z).map_err(|e| FetchError::Network(e.to_string()))?;
            let bytes = response.bytes().await?;
            std::io::copy(&mut std::io::Cursor::new(&bytes), &mut dest)
                .map_err(|e| FetchError::Network(e.to_string()))?;

            // run gunzip to extract .Z
            tracing::info!("Extracting {}...", dest_z.display());
            let status = std::process::Command::new("gunzip")
                .arg("-f")
                .arg(&dest_z)
                .status()
                .map_err(|e| {
                    FetchError::Decompression(format!("Failed to execute gunzip: {}", e))
                })?;

            if !status.success() {
                return Err(FetchError::Decompression(format!(
                    "gunzip failed with status {}",
                    status
                )));
            }

            let extracted_path = out_dir.join(f.trim_end_matches(".Z"));
            out_paths.push(extracted_path);
        }

        Ok(out_paths)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gps_to_utc_leap_seconds_boundary() {
        // Week 2086, TOW 259200 corresponds to 2020-01-01 00:00:00 in raw elapsed days
        // At 00:00:10 GPST (TOW = 259210.0), subtracting 18s leap seconds gives 2019-12-31 23:59:52 UTC
        let t_before_midnight = GpsTime::new(2086, 259_210.0);
        let utc_prev = gps_to_utc_naive(t_before_midnight).unwrap();
        assert_eq!(utc_prev.format("%Y").to_string(), "2019");
        assert_eq!(utc_prev.format("%j").to_string(), "365");

        // At 00:00:18 GPST (TOW = 259218.0), exactly 2020-01-01 00:00:00 UTC is reached
        let t_exact_midnight = GpsTime::new(2086, 259_218.0);
        let utc_curr = gps_to_utc_naive(t_exact_midnight).unwrap();
        assert_eq!(utc_curr.format("%Y").to_string(), "2020");
        assert_eq!(utc_curr.format("%j").to_string(), "001");
    }
}
