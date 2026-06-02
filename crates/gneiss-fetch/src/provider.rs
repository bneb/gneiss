use async_trait::async_trait;
use std::path::PathBuf;
use gneiss_core::coords::Coordinate;
use gneiss_core::time::GpsTime;

#[derive(Debug)]
pub enum FetchError {
    Network(String),
    NotFound(String),
    Decompression(String),
}

impl From<reqwest::Error> for FetchError {
    fn from(e: reqwest::Error) -> Self {
        FetchError::Network(e.to_string())
    }
}

impl std::fmt::Display for FetchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FetchError::Network(s) => write!(f, "Network Error: {}", s),
            FetchError::NotFound(s) => write!(f, "Not Found: {}", s),
            FetchError::Decompression(s) => write!(f, "Decompression Error: {}", s),
        }
    }
}

impl std::error::Error for FetchError {}

#[async_trait]
pub trait DataSource {
    /// Returns the name of the provider (e.g., "NOAA_CORS", "Proprietary_API")
    fn name(&self) -> &str;
    
    /// Fetches the closest RINEX observation file for a given coordinate and time
    async fn fetch_base_obs(&self, location: Coordinate, time: GpsTime, out_dir: &std::path::Path) -> Result<PathBuf, FetchError>;
    
    /// Fetches the global broadcast ephemeris for a specific time
    async fn fetch_ephemeris(&self, time: GpsTime, out_dir: &std::path::Path) -> Result<PathBuf, FetchError>;
}
