use std::path::{Path, PathBuf};
use flate2::read::GzDecoder;
use std::io::Read;
use tracing::info;
use crate::provider::FetchError;

/// ResourceResolver resolves local and remote URIs (HTTP/HTTPS, CORS schemas, S3)
/// to local cached and decompressed file paths.
pub struct ResourceResolver {
    cache_dir: PathBuf,
    client: reqwest::Client,
}

impl Default for ResourceResolver {
    fn default() -> Self {
        Self::new()
    }
}

impl ResourceResolver {
    /// Creates a new ResourceResolver with default cache path (~/.cache/gneiss).
    pub fn new() -> Self {
        let cache_dir = dirs_cache_dir().unwrap_or_else(|| PathBuf::from(".gneiss_cache"));
        std::fs::create_dir_all(&cache_dir).ok();
        Self {
            cache_dir,
            client: reqwest::Client::builder()
                .user_agent("gneiss-fetch/0.1.0")
                .build()
                .unwrap_or_default(),
        }
    }

    /// Resolves a path, URL, or CORS URI to a local decompressed file path.
    pub async fn resolve(&self, uri: &str) -> Result<PathBuf, FetchError> {
        let trimmed = uri.trim();
        if trimmed.starts_with("http://") || trimmed.starts_with("https://") {
            self.fetch_url(trimmed).await
        } else if let Some(station) = trimmed.strip_prefix("cors://").or_else(|| trimmed.strip_prefix("noaa://")) {
            self.fetch_noaa_cors(station).await
        } else if let Some(station) = trimmed.strip_prefix("bkg://").or_else(|| trimmed.strip_prefix("euref://")) {
            self.fetch_bkg_cors(station).await
        } else {
            let path = PathBuf::from(trimmed);
            self.decompress_local(&path)
        }
    }

    async fn fetch_url(&self, url: &str) -> Result<PathBuf, FetchError> {
        let hash = compute_url_hash(url);
        let ext = extract_extension(url);
        let cached_file = self.cache_dir.join(format!("{}.{}", hash, ext));

        if !cached_file.exists() {
            info!("Fetching remote GNSS resource: {}", url);
            let resp = self.client.get(url).send().await
                .map_err(|e| FetchError::Network(e.to_string()))?;
            if !resp.status().is_success() {
                return Err(FetchError::Network(format!("HTTP {} from {}", resp.status(), url)));
            }
            let bytes = resp.bytes().await
                .map_err(|e| FetchError::Network(e.to_string()))?;
            std::fs::write(&cached_file, &bytes)
                .map_err(|e| FetchError::Network(e.to_string()))?;
        }

        self.decompress_local(&cached_file)
    }

    async fn fetch_noaa_cors(&self, station: &str) -> Result<PathBuf, FetchError> {
        let stn = station.to_ascii_lowercase();
        let s3_url = format!(
            "https://noaa-cors-pds.s3.amazonaws.com/rinex/latest/{}.obs.gz",
            stn
        );
        info!("Resolving NOAA CORS station {}: {}", stn, s3_url);
        self.fetch_url(&s3_url).await
    }

    async fn fetch_bkg_cors(&self, station: &str) -> Result<PathBuf, FetchError> {
        let stn = station.to_ascii_lowercase();
        let url = format!(
            "https://igs.bkg.bund.de/root_ftp/IGS/obs/latest/{}.rnx.gz",
            stn
        );
        info!("Resolving BKG EUREF station {}: {}", stn, url);
        self.fetch_url(&url).await
    }

    fn decompress_local(&self, path: &Path) -> Result<PathBuf, FetchError> {
        if !path.exists() {
            return Err(FetchError::NotFound(format!("File does not exist: {}", path.display())));
        }

        let path_str = path.to_string_lossy();
        if path_str.ends_with(".gz") {
            self.decompress_gzip(path)
        } else if path_str.ends_with(".crx") || path_str.ends_with(".d") {
            crate::hatanaka::decompress(path, &self.cache_dir)
        } else {
            Ok(path.to_path_buf())
        }
    }

    fn decompress_gzip(&self, gz_path: &Path) -> Result<PathBuf, FetchError> {
        let stem = gz_path.file_stem().unwrap_or_default().to_string_lossy();
        let target = self.cache_dir.join(format!("decomp_{}", stem));

        if !target.exists() {
            let file = std::fs::File::open(gz_path)
                .map_err(|e| FetchError::Decompression(e.to_string()))?;
            let mut decoder = GzDecoder::new(file);
            let mut buffer = Vec::new();
            decoder.read_to_end(&mut buffer)
                .map_err(|e| FetchError::Decompression(e.to_string()))?;
            std::fs::write(&target, &buffer)
                .map_err(|e| FetchError::Decompression(e.to_string()))?;
        }

        if target.to_string_lossy().ends_with(".crx") || target.to_string_lossy().ends_with(".d") {
            crate::hatanaka::decompress(&target, &self.cache_dir)
        } else {
            Ok(target)
        }
    }
}

fn compute_url_hash(url: &str) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for b in url.as_bytes() {
        h ^= *b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{:016x}", h)
}

fn extract_extension(url: &str) -> String {
    let clean = url.split('?').next().unwrap_or(url);
    let last_segment = clean.split('/').next_back().unwrap_or(clean);
    if let Some(idx) = last_segment.rfind('.') {
        last_segment[idx + 1..].to_string()
    } else {
        "dat".to_string()
    }
}

fn dirs_cache_dir() -> Option<PathBuf> {
    if let Ok(home) = std::env::var("HOME") {
        Some(PathBuf::from(home).join(".cache").join("gneiss"))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_url_hash_deterministic() {
        let h1 = compute_url_hash("https://example.com/test.obs.gz");
        let h2 = compute_url_hash("https://example.com/test.obs.gz");
        assert_eq!(h1, h2);
        assert_eq!(h1.len(), 16);
    }

    #[test]
    fn test_extract_extension() {
        assert_eq!(extract_extension("http://site.org/data.25o.gz?token=123"), "gz");
        assert_eq!(extract_extension("http://site.org/data.ubx"), "ubx");
        assert_eq!(extract_extension("http://site.org/plain"), "dat");
    }

    #[tokio::test]
    async fn test_resolve_local_file() {
        let resolver = ResourceResolver::new();
        let path = PathBuf::from("Cargo.toml");
        let res = resolver.resolve("Cargo.toml").await;
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), path);
    }
}
