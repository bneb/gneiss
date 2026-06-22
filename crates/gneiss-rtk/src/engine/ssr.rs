use gneiss_core::sat::SatelliteId;
use gneiss_core::time::GpsTime;
use nalgebra::Vector3;

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
pub fn interpolate_orbit_lagrange(
    points: &[(GpsTime, Vector3<f64>)],
    target: GpsTime,
) -> Option<Vector3<f64>> {
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
) -> Option<(Vector3<f64>, Vector3<f64>, f64)> {
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

    // Check if the points are too far away in time (e.g., > 2 hours)
    if min_dt > 7200.0 {
        return None;
    }

    let subset = &valid_points[start_idx..end_idx];
    let pos = interpolate_orbit_lagrange(subset, t);

    // Compute velocity via central difference (dt = 1.0 second is small enough for orbit, large enough for float precision)
    for p in subset {
        tracing::debug!("SP3 point: time={}, pos={:?}", p.0.tow, p.1);
    }
    let pos_plus = interpolate_orbit_lagrange(subset, t + 0.5);
    let pos_minus = interpolate_orbit_lagrange(subset, t - 0.5);

    if let (Some(p), Some(p_plus), Some(p_minus)) = (pos, pos_plus, pos_minus) {
        let vel = (p_plus - p_minus) / 1.0;
        Some((p, vel, clock_bias.unwrap_or(0.0)))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::sat::Constellation;

    #[test]
    fn test_format_sp3_id_gps() {
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        assert_eq!(format_sp3_id(sat), "G01");
    }

    #[test]
    fn test_format_sp3_id_glonass() {
        let sat = SatelliteId {
            constellation: Constellation::Glonass,
            prn: 24,
        };
        assert_eq!(format_sp3_id(sat), "R24");
    }

    #[test]
    fn test_format_sp3_id_galileo() {
        let sat = SatelliteId {
            constellation: Constellation::Galileo,
            prn: 7,
        };
        assert_eq!(format_sp3_id(sat), "E07");
    }

    #[test]
    fn test_format_sp3_id_beidou() {
        let sat = SatelliteId {
            constellation: Constellation::Beidou,
            prn: 10,
        };
        assert_eq!(format_sp3_id(sat), "C10");
    }

    #[test]
    fn test_format_sp3_id_qzss() {
        let sat = SatelliteId {
            constellation: Constellation::Qzss,
            prn: 3,
        };
        assert_eq!(format_sp3_id(sat), "J03");
    }

    #[test]
    fn test_format_sp3_id_sbas_fallback() {
        // SBAS and Navic are not explicitly listed -> fallback to '?'
        let sat = SatelliteId {
            constellation: Constellation::Sbas,
            prn: 5,
        };
        assert_eq!(format_sp3_id(sat), "?05");
    }

    #[test]
    fn test_lagrange_empty_slice() {
        let points: Vec<(GpsTime, Vector3<f64>)> = vec![];
        let target = GpsTime::new(0, 0.0);
        assert!(interpolate_orbit_lagrange(&points, target).is_none());
    }

    #[test]
    fn test_lagrange_single_point() {
        let points = vec![(
            GpsTime {
                week: 0,
                tow: 10.0,
            },
            Vector3::new(42.0, -17.0, 99.0),
        )];
        let target = GpsTime {
            week: 0,
            tow: 20.0,
        };
        // Single point is returned regardless of target time
        let result = interpolate_orbit_lagrange(&points, target).unwrap();
        assert!((result.x - 42.0).abs() < 1e-10);
        assert!((result.y + 17.0).abs() < 1e-10);
        assert!((result.z - 99.0).abs() < 1e-10);
    }

    #[test]
    fn test_lagrange_linear_interpolation() {
        // y = 3x + 1 over x=0..10
        let mut points = Vec::new();
        for x in 0..=10 {
            let t = GpsTime {
                week: 0,
                tow: x as f64,
            };
            points.push((t, Vector3::new(3.0 * x as f64 + 1.0, 0.0, 0.0)));
        }
        let target = GpsTime {
            week: 0,
            tow: 4.5,
        };
        let result = interpolate_orbit_lagrange(&points, target).unwrap();
        // Expected: 3*4.5 + 1 = 14.5
        assert!((result.x - 14.5).abs() < 1e-9);
    }

    #[test]
    fn test_lagrange_interpolation() {
        let mut points = Vec::new();
        // y = x^2
        for i in -5..=5 {
            let t = GpsTime {
                week: 0,
                tow: i as f64,
            };
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

    #[test]
    fn test_lagrange_duplicate_time_does_not_crash() {
        // Two points at the same time: denominator is 0, the function should
        // continue rather than panic from a division by zero.
        let points = vec![
            (
                GpsTime { week: 0, tow: 10.0 },
                Vector3::new(42.0, 0.0, 0.0),
            ),
            (
                GpsTime { week: 0, tow: 10.0 },
                Vector3::new(99.0, 0.0, 0.0),
            ),
        ];
        let target = GpsTime { week: 0, tow: 12.0 };
        // Should not panic; result is implementation-defined due to degenerate input
        let _result = interpolate_orbit_lagrange(&points, target);
    }

    #[test]
    fn test_lagrange_fewer_points_than_degree() {
        // Using 3 points for degree 4 pattern: should still work (n=3, interpolate uses all 3)
        let mut points = Vec::new();
        for x in 0..3 {
            let t = GpsTime { week: 0, tow: x as f64 };
            points.push((t, Vector3::new(x as f64, 0.0, 0.0)));
        }
        let target = GpsTime { week: 0, tow: 1.5 };
        let result = interpolate_orbit_lagrange(&points, target).unwrap();
        // Linear interpolation through [0,1,2] at 1.5 should give 1.5
        assert!((result.x - 1.5).abs() < 1e-9);
    }

    #[test]
    fn test_lagrange_vector_interpolation() {
        // Interpolate a vector with x=t, y=2t, z=3t
        let mut points = Vec::new();
        for x in 0..=5 {
            let t = GpsTime { week: 0, tow: x as f64 };
            points.push((t, Vector3::new(x as f64, 2.0 * x as f64, 3.0 * x as f64)));
        }
        let target = GpsTime { week: 0, tow: 2.5 };
        let result = interpolate_orbit_lagrange(&points, target).unwrap();
        assert!((result.x - 2.5).abs() < 1e-9);
        assert!((result.y - 5.0).abs() < 1e-9);
        assert!((result.z - 7.5).abs() < 1e-9);
    }

    // --- get_precise_orbit tests ---

    use std::collections::HashMap;
    use gneiss_parsers::sp3::Sp3Record;

    fn epoch_at(week: u32, tow: f64, sat_id: &str, pos: Vector3<f64>, clock: f64) -> gneiss_parsers::sp3::Sp3Epoch {
        let mut records = HashMap::new();
        records.insert(
            sat_id.to_string(),
            Sp3Record {
                position: pos,
                clock_offset: clock,
            },
        );
        gneiss_parsers::sp3::Sp3Epoch {
            time: GpsTime { week, tow },
            records,
        }
    }

    #[test]
    fn test_get_precise_orbit_simple() {
        // Create SP3 epochs with GPS satellite G01 at known positions
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let mut epochs = Vec::new();
        for i in 0..=5 {
            let t = GpsTime { week: 2200, tow: i as f64 * 300.0 };
            let pos = Vector3::new(
                10000000.0 + i as f64 * 100.0,
                20000000.0 + i as f64 * 50.0,
                15000000.0 + i as f64 * 75.0,
            );
            epochs.push(epoch_at(2200, i as f64 * 300.0, "G01", pos, 0.001));
        }
        // Target time at 600s (between epoch 2 at 600s and epoch 3 at 900s)
        let target = GpsTime { week: 2200, tow: 600.0 };
        let result = get_precise_orbit(&epochs, sat, target, 4);
        assert!(result.is_some(), "get_precise_orbit should find a result");
        let (pos, vel, clk) = result.unwrap();
        // Position at t=600 should be epoch 2's position (since 600 is exactly at epoch 2)
        assert!((pos.x - 10000200.0).abs() < 100.0, "pos.x should be near 10000200, got {}", pos.x);
        assert!((pos.y - 20000100.0).abs() < 100.0, "pos.y should be near 20000100, got {}", pos.y);
        assert!((pos.z - 15000150.0).abs() < 100.0, "pos.z should be near 15000150, got {}", pos.z);
        // Velocity should be non-zero (central difference of positions)
        assert!(vel.norm() > 0.0, "Velocity should be non-zero");
        // Clock bias should be from nearest epoch
        assert!((clk - 0.001).abs() < 1e-12, "Clock bias should be 0.001, got {}", clk);
    }

    #[test]
    fn test_get_precise_orbit_empty_epochs() {
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let epochs: Vec<gneiss_parsers::sp3::Sp3Epoch> = vec![];
        let target = GpsTime { week: 2200, tow: 0.0 };
        let result = get_precise_orbit(&epochs, sat, target, 4);
        assert!(result.is_none(), "Should return None for empty epochs");
    }

    #[test]
    fn test_get_precise_orbit_sat_not_found() {
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        // Only include a different satellite
        let epochs = vec![epoch_at(2200, 0.0, "G02", Vector3::new(1.0, 0.0, 0.0), 0.0)];
        let target = GpsTime { week: 2200, tow: 0.0 };
        let result = get_precise_orbit(&epochs, sat, target, 4);
        assert!(result.is_none(), "Should return None when satellite is not found");
    }

    #[test]
    fn test_get_precise_orbit_too_far_away() {
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        // Epochs too far from target (> 2 hours)
        let mut epochs = Vec::new();
        for i in 0..5 {
            let t = GpsTime { week: 2200, tow: 10000.0 + i as f64 * 300.0 };
            epochs.push(epoch_at(2200, 10000.0 + i as f64 * 300.0, "G01", Vector3::new(1.0, 0.0, 0.0), 0.0));
        }
        // Target at 0.0, but earliest epoch is at 10000 (2.78 hours away)
        let target = GpsTime { week: 2200, tow: 0.0 };
        let result = get_precise_orbit(&epochs, sat, target, 4);
        assert!(result.is_none(), "Should return None when target is too far from epochs");
    }

    #[test]
    fn test_get_precise_orbit_clock_nan() {
        // When clock_offset is NaN, the nearest valid clock should be used from the next epoch
        let sat = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let mut records_nan = HashMap::new();
        records_nan.insert("G01".to_string(), Sp3Record {
            position: Vector3::new(10000000.0, 0.0, 0.0),
            clock_offset: f64::NAN,
        });
        let mut records_valid = HashMap::new();
        records_valid.insert("G01".to_string(), Sp3Record {
            position: Vector3::new(10000200.0, 0.0, 0.0),
            clock_offset: 0.002,
        });
        let epochs = vec![
            gneiss_parsers::sp3::Sp3Epoch {
                time: GpsTime { week: 2200, tow: 0.0 },
                records: records_nan,
            },
            gneiss_parsers::sp3::Sp3Epoch {
                time: GpsTime { week: 2200, tow: 300.0 },
                records: records_valid,
            },
        ];
        let target = GpsTime { week: 2200, tow: 0.0 };
        let result = get_precise_orbit(&epochs, sat, target, 2);
        assert!(result.is_some());
        let (_, _, clk) = result.unwrap();
        // Clock should be 0.002 (from the valid epoch, since the one at t=0 has NaN)
        assert!((clk - 0.002).abs() < 1e-12, "Expected clock 0.002 from valid record, got {}", clk);
    }
}
