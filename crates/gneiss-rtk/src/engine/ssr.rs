use nalgebra::Vector3;
use gneiss_core::time::GpsTime;
use gneiss_core::sat::SatelliteId;

pub fn format_sp3_id(sat: SatelliteId) -> String {
    let c = match sat.constellation {
        gneiss_core::sat::Constellation::Gps => 'G',
        gneiss_core::sat::Constellation::Glonass => 'R',
        gneiss_core::sat::Constellation::Galileo => 'E',
        gneiss_core::sat::Constellation::Beidou => 'C',
        gneiss_core::sat::Constellation::Qzss => 'J',
        _ => '?',
    };
    format!("{}{:02}", c, sat.prn)
}

/// Interpolates precise orbit using an N-point Lagrange polynomial.
/// SP3 epochs should be provided as a slice of (time, Vector3) sorted by time.
pub fn interpolate_orbit_lagrange(points: &[(GpsTime, Vector3<f64>)], target: GpsTime) -> Option<Vector3<f64>> {
    let n = points.len();
    if n == 0 {
        return None;
    }
    if n == 1 {
        return Some(points[0].1);
    }

    let mut result = Vector3::zeros();
    for i in 0..n {
        let mut term = points[i].1;
        for j in 0..n {
            if i != j {
                let num = target - points[j].0;
                let den = points[i].0 - points[j].0;
                if den == 0.0 {
                    continue;
                }
                term *= num / den;
            }
        }
        result += term;
    }

    Some(result)
}

/// Finds the best N points around the target time and interpolates.
/// SP3 files are usually sorted by time.
pub fn get_precise_orbit(
    sp3_epochs: &[gneiss_parsers::sp3::Sp3Epoch],
    sat: SatelliteId,
    t: GpsTime,
    degree: usize,
) -> Option<(Vector3<f64>, f64)> {
    let sat_id = format_sp3_id(sat);

    // Extract valid points for this satellite
    let mut valid_points = Vec::new();
    let mut clock_bias = None;
    let mut clock_diff = f64::MAX;

    for epoch in sp3_epochs {
        if let Some(record) = epoch.records.get(&sat_id) {
            valid_points.push((epoch.time, record.position));
            
            // For clock bias from SP3 (if RINEX CLK is unavailable), just use nearest neighbor or linear.
            // But usually we just take the nearest if it's within a threshold.
            let dt = (epoch.time - t).abs();
            if dt < clock_diff && !record.clock_offset.is_nan() {
                clock_diff = dt;
                clock_bias = Some(record.clock_offset);
            }
        }
    }

    if valid_points.is_empty() {
        return None;
    }

    // Find the closest index
    let mut closest_idx: usize = 0;
    let mut min_dt = f64::MAX;
    for (i, (pt_t, _)) in valid_points.iter().enumerate() {
        let dt = (*pt_t - t).abs();
        if dt < min_dt {
            min_dt = dt;
            closest_idx = i;
        }
    }

    // Select `n_points` around `closest_idx`
    let n_points = degree + 1;
    let mut start_idx = closest_idx.saturating_sub(n_points / 2);
    let mut end_idx = start_idx + n_points;

    if end_idx > valid_points.len() {
        end_idx = valid_points.len();
        start_idx = end_idx.saturating_sub(n_points);
    }

    let slice = &valid_points[start_idx..end_idx];
    
    // Check if the points are too far away in time (e.g., > 2 hours)
    if min_dt > 7200.0 {
        return None;
    }

    let pos = interpolate_orbit_lagrange(slice, t)?;
    let clk = clock_bias.unwrap_or(0.0);

    Some((pos, clk))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lagrange_interpolation() {
        let mut points = Vec::new();
        // y = x^2
        for i in -5..=5 {
            let t = GpsTime { week: 0, tow: i as f64 };
            let v = Vector3::new((i * i) as f64, 0.0, 0.0);
            points.push((t, v));
        }

        // Interpolate at x = 1.5. y should be 2.25
        let t_target = GpsTime { week: 0, tow: 1.5 };
        let result = interpolate_orbit_lagrange(&points, t_target).unwrap();
        assert!((result.x - 2.25).abs() < 1e-9);
        assert!((result.y - 0.0).abs() < 1e-9);
        assert!((result.z - 0.0).abs() < 1e-9);
    }
}
