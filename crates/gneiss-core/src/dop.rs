use nalgebra::{DMatrix, Vector3};

/// Dilution of Precision values computed from satellite geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DopValues {
    /// Geometric DOP (3D position + time)
    pub gdop: f64,
    /// Position DOP (3D position)
    pub pdop: f64,
    /// Horizontal DOP (2D horizontal)
    pub hdop: f64,
    /// Vertical DOP
    pub vdop: f64,
    /// Time DOP
    pub tdop: f64,
}

/// Computes DOP values from a geometry matrix H.
/// H must be an n×4 matrix where each row is [dx/r, dy/r, dz/r, 1.0]
/// representing the unit line-of-sight vectors from receiver to each satellite
/// plus the clock column.
///
/// Returns None if the geometry is degenerate (matrix not invertible).
pub fn compute_dop(h: &DMatrix<f64>) -> Option<DopValues> {
    if h.nrows() < 4 || h.ncols() != 4 {
        return None;
    }

    let ht = h.transpose();
    let hth = &ht * h;
    let q = hth.try_inverse()?;

    // Ensure diagonal elements are non-negative
    if q[(0, 0)] < 0.0 || q[(1, 1)] < 0.0 || q[(2, 2)] < 0.0 || q[(3, 3)] < 0.0 {
        return None;
    }

    let gdop = libm::sqrt(q[(0, 0)] + q[(1, 1)] + q[(2, 2)] + q[(3, 3)]);
    let pdop = libm::sqrt(q[(0, 0)] + q[(1, 1)] + q[(2, 2)]);
    let tdop = libm::sqrt(q[(3, 3)]);

    // For HDOP/VDOP we need the local ENU frame decomposition.
    // If operating in ECEF, HDOP ≈ sqrt(q11 + q22) and VDOP ≈ sqrt(q33) is only
    // approximate. For a proper decomposition, we'd need the receiver position.
    // For now, provide the ECEF-based approximation.
    let hdop = libm::sqrt(q[(0, 0)] + q[(1, 1)]);
    let vdop = libm::sqrt(q[(2, 2)]);

    Some(DopValues {
        gdop,
        pdop,
        hdop,
        vdop,
        tdop,
    })
}

/// Computes DOP values given receiver ECEF position and satellite ECEF positions.
/// Automatically constructs the geometry matrix from the positions.
pub fn compute_dop_from_positions(
    rcv_ecef: Vector3<f64>,
    sat_positions: &[Vector3<f64>],
) -> Option<DopValues> {
    if sat_positions.len() < 4 {
        return None;
    }

    let n = sat_positions.len();
    let mut h = DMatrix::zeros(n, 4);

    for (i, sat_pos) in sat_positions.iter().enumerate() {
        let diff = rcv_ecef - *sat_pos;
        let range = diff.norm();
        if range < 1.0 {
            return None;
        }
        h[(i, 0)] = diff.x / range;
        h[(i, 1)] = diff.y / range;
        h[(i, 2)] = diff.z / range;
        h[(i, 3)] = 1.0;
    }

    compute_dop(&h)
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;
    use nalgebra::Vector3;

    #[test]
    fn test_dop_well_distributed_geometry() {
        // 6 satellites well-distributed in sky (approximate tetrahedron + 2)
        let rcv = Vector3::new(crate::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let sats = vec![
            Vector3::new(20000000.0, 5000000.0, 5000000.0),
            Vector3::new(22000000.0, -5000000.0, 5000000.0),
            Vector3::new(19000000.0, 5000000.0, -5000000.0),
            Vector3::new(21000000.0, -5000000.0, -5000000.0),
            Vector3::new(25000000.0, 0.0, 0.0),
            Vector3::new(15000000.0, 0.0, 10000000.0),
        ];

        let dop = compute_dop_from_positions(rcv, &sats).expect("Should compute DOP");

        // With good geometry, PDOP should be reasonable (< 5)
        assert!(
            dop.pdop > 0.0 && dop.pdop < 10.0,
            "PDOP should be reasonable, got {}",
            dop.pdop
        );
        assert!(dop.gdop > dop.pdop, "GDOP must be >= PDOP");
        assert!(dop.pdop >= dop.hdop, "PDOP must be >= HDOP");
        // GDOP^2 = PDOP^2 + TDOP^2
        let gdop_check = libm::sqrt(dop.pdop * dop.pdop + dop.tdop * dop.tdop);
        assert!(
            (dop.gdop - gdop_check).abs() < 1e-6,
            "GDOP^2 != PDOP^2 + TDOP^2"
        );
    }

    #[test]
    fn test_dop_coplanar_degenerate() {
        // All satellites in the same plane → poor VDOP, may not invert
        let rcv = Vector3::new(crate::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let sats = vec![
            Vector3::new(20000000.0, 5000000.0, 0.0),
            Vector3::new(20000000.0, -5000000.0, 0.0),
            Vector3::new(25000000.0, 5000000.0, 0.0),
            Vector3::new(25000000.0, -5000000.0, 0.0),
        ];

        // Coplanar sats (all z=0) should either return None or very high DOP
        if let Some(dop) = compute_dop_from_positions(rcv, &sats) {
            assert!(
                dop.vdop > 100.0,
                "Coplanar sats should have very high VDOP, got {}",
                dop.vdop
            );
        }
    }

    #[test]
    fn test_dop_insufficient_sats() {
        let rcv = Vector3::new(crate::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let sats = vec![
            Vector3::new(20000000.0, 0.0, 0.0),
            Vector3::new(22000000.0, 0.0, 0.0),
        ];

        assert!(
            compute_dop_from_positions(rcv, &sats).is_none(),
            "Should return None with < 4 sats"
        );
    }

    #[test]
    fn test_dop_more_sats_improves() {
        let rcv = Vector3::new(crate::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let base_sats = vec![
            Vector3::new(20000000.0, 5000000.0, 5000000.0),
            Vector3::new(22000000.0, -5000000.0, 5000000.0),
            Vector3::new(19000000.0, 5000000.0, -5000000.0),
            Vector3::new(21000000.0, -5000000.0, -5000000.0),
        ];

        let dop4 = compute_dop_from_positions(rcv, &base_sats).unwrap();

        let mut more_sats = base_sats.clone();
        more_sats.push(Vector3::new(25000000.0, 0.0, 0.0));
        more_sats.push(Vector3::new(15000000.0, 0.0, 10000000.0));

        let dop6 = compute_dop_from_positions(rcv, &more_sats).unwrap();

        // More well-distributed sats should improve (lower) GDOP
        assert!(
            dop6.gdop < dop4.gdop,
            "6 sats should have lower GDOP ({}) than 4 sats ({})",
            dop6.gdop,
            dop4.gdop
        );
    }
}
