use std::path::{Path, PathBuf};
use crate::provider::FetchError;

/// Decompresses a Hatanaka-compressed RINEX file (.crx) to standard RINEX (.obs).
/// This typically requires the `crx2rnx` utility to be available on the system PATH.
pub fn decompress(input_path: &Path, out_dir: &Path) -> Result<PathBuf, FetchError> {
    // Generate output file name. e.g. cnmr1000.23d -> cnmr1000.23o
    let file_name = input_path.file_name()
        .ok_or_else(|| FetchError::Decompression("Invalid input path".into()))?
        .to_string_lossy();
    
    let out_file_name = if file_name.ends_with('d') {
        file_name.replace('d', "o")
    } else {
        format!("{}.obs", file_name)
    };

    let output_path = out_dir.join(out_file_name);

    tracing::info!("Decompressing Hatanaka file {} to {}", input_path.display(), output_path.display());

    // crx2rnx replaces the 'd' extension with 'o' natively.
    let status = std::process::Command::new("crx2rnx")
        .arg(input_path.file_name().unwrap_or(input_path.as_os_str()))
        .arg("-f") // Force overwrite
        .arg("-d") // Delete input file if successful
        .current_dir(out_dir)
        .status()
        .map_err(|e| FetchError::Decompression(format!("Failed to execute crx2rnx: {}. Is it installed?", e)))?;

    if !status.success() {
        return Err(FetchError::Decompression(format!("crx2rnx failed with status {}", status)));
    }

    Ok(output_path)
}
