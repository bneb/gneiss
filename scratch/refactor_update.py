import os
import re

file_path = "crates/gneiss-rtk/src/engine/updater.rs"

with open(file_path, "r") as f:
    content = f.read()

# Add evaluate_post_fit_outliers function
new_func = """
pub fn evaluate_post_fit_outliers(
    v: &DVector<f64>,
    s: &DMatrix<f64>,
    current_z: &DVector<f64>,
    current_valid: &[usize],
    meas_types: Option<&[(gneiss_core::sat::SatelliteId, u8)]>,
    max_innovation: f64,
    is_tightly_coupled: bool,
) -> (Option<usize>, f64) {
    let mut max_outlier_ratio = 0.0;
    let mut worst_idx = None;
    
    for i in 0..v.len() {
        let orig_idx = current_valid[i];
        let meas_type = meas_types.map_or(0, |m| m[orig_idx].1);
        
        let s_ii = s[(i, i)];
        let ratio = v[i].abs() / s_ii.sqrt();
        
        let thresh = match meas_type {
            0 => max_innovation, 
            1 | 2 => 5.0, // 5.0m for Phase to allow tracking vehicle dynamics over outages
            3 => max_innovation * 2.0, // Doppler max innovation (typically 15-30m/s)
            _ => max_innovation,
        };
        
        let abs_thresh = match meas_type {
            0 => 40.0, // Pseudorange max 40m error
            1 | 2 => 1.0, // Phase max 1m error
            3 => 15.0, // Doppler max 15m/s error
            _ => 40.0,
        };
        
        if (v[i].abs() > thresh && ratio > max_outlier_ratio) || (is_tightly_coupled && current_z[i].abs() > abs_thresh) {
            if v[i].abs() > thresh && ratio > max_outlier_ratio {
                max_outlier_ratio = ratio;
                worst_idx = Some(i);
            } else if is_tightly_coupled && current_z[i].abs() > abs_thresh {
                // Force rejection if absolute value is too high (only for INS to protect attitude)
                worst_idx = Some(i);
                max_outlier_ratio = f64::INFINITY;
            }
        }
    }
    
    (worst_idx, max_outlier_ratio)
}
"""

if "evaluate_post_fit_outliers" not in content:
    # Insert before update()
    content = content.replace("pub fn update(", new_func + "\npub fn update(")

# Replace the loop in update()
old_loop = """        let mut max_outlier_ratio = 0.0;
        let mut worst_idx = None;
        
        for i in 0..v.len() {
            let orig_idx = current_valid[i];
            let meas_type = meas_types.map_or(0, |m| m[orig_idx].1);
            
            let s_ii = s[(i, i)];
            let ratio = v[i].abs() / s_ii.sqrt();
            
            let thresh = match meas_type {
                0 => max_innovation, 
                1 | 2 => 5.0, // 5.0m for Phase to allow tracking vehicle dynamics over outages
                3 => max_innovation * 2.0, // Doppler max innovation (typically 15-30m/s)
                _ => max_innovation,
            };
            
            let abs_thresh = match meas_type {
                0 => 40.0, // Pseudorange max 40m error
                1 | 2 => 1.0, // Phase max 1m error
                3 => 15.0, // Doppler max 15m/s error
                _ => 40.0,
            };
            
            if (v[i].abs() > thresh && ratio > max_outlier_ratio) || (is_tightly_coupled && current_z[i].abs() > abs_thresh) {
                if v[i].abs() > thresh && ratio > max_outlier_ratio {
                    max_outlier_ratio = ratio;
                    worst_idx = Some(i);
                } else if is_tightly_coupled && current_z[i].abs() > abs_thresh {
                    // Force rejection if absolute value is too high (only for INS to protect attitude)
                    worst_idx = Some(i);
                    max_outlier_ratio = f64::INFINITY;
                }
            }
        }"""

new_call = """        let (worst_idx, max_outlier_ratio) = evaluate_post_fit_outliers(&v, &s, &current_z, &current_valid, meas_types, max_innovation, is_tightly_coupled);"""

content = content.replace(old_loop, new_call)

with open(file_path, "w") as f:
    f.write(content)
