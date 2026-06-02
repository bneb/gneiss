#[cfg(test)]
mod tests {
    use super::super::coords::ecef_to_llh;
    use nalgebra::Vector3;
    

    #[test]
    fn test_ecef_to_llh_precision() {
        // WGS84 ECEF: [6378137.0, 0, 0] (Equator, Lon 0, Height 0)
        let pos = Vector3::new(6378137.0, 0.0, 0.0);
        let llh = ecef_to_llh(pos);
        assert!(llh.x.abs() < 1e-12); // Lat 0
        assert!(llh.y.abs() < 1e-12); // Lon 0
        assert!(llh.z.abs() < 1e-4);  // Height 0
    }
}
