//! Satellite L1 phase-centre Z-offset correction for precise orbits.
//!
//! SP3 positions reference the satellite centre of mass; observations
//! reference the antenna phase centre along body +Z (toward Earth).
//! Offsets are 0.6-1.55 m depending on block type (igs14.atx).

use nalgebra::Vector3;

const SAT_PCO_Z_MM: [f64; 33] = [
    0.0, 1501.8, 728.8, 1550.6, 1232.4, 778.0, 1467.0, 822.4, 1501.4, 1522.6, 1515.1, 1232.4,
    767.8, 1348.3, 1232.4, 622.8, 1468.7, 770.9, 1232.4, 808.2, 1313.5, 1359.1, 1304.5, 1232.4,
    1407.1, 1517.4, 1503.5, 1522.3, 1232.4, 791.8, 1522.1, 912.5, 1534.8,
];

/// Move SP3 centre-of-mass position to L1 phase centre via nadir projection.
pub fn apply_sat_pco_z(pos_com: Vector3<f64>, prn: u16) -> Vector3<f64> {
    if prn == 0 || prn as usize >= SAT_PCO_Z_MM.len() {
        return pos_com;
    }
    let pco_z_m = SAT_PCO_Z_MM[prn as usize] * 1e-3;
    pos_com + (-pos_com.normalize()) * pco_z_m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_correction_shifts_toward_earth() {
        let pos = Vector3::new(26_560_000.0, 0.0, 0.0);
        let corrected = apply_sat_pco_z(pos, 1);
        assert!(corrected.norm() < pos.norm());
        let diff = pos.norm() - corrected.norm();
        assert!((diff - 1.5018).abs() < 0.001);
    }

    #[test]
    fn test_invalid_prn_passthrough() {
        let pos = Vector3::new(1e7, 2e7, 3e7);
        assert_eq!(apply_sat_pco_z(pos, 0), pos);
    }

    #[test]
    fn test_all_prns_realistic_range() {
        for prn in 1..=32u16 {
            let pco = SAT_PCO_Z_MM[prn as usize];
            assert!(pco > 500.0 && pco < 1600.0, "PRN {} = {}", prn, pco);
        }
    }
}
