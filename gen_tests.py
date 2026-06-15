test_code = """
#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Vector3, Matrix3};
    use gneiss_core::time::GpsTime;
    use gneiss_core::ephemeris::Ephemeris;
    use gneiss_core::sat::SatelliteId;
    use gneiss_core::coords::{Datum, Frame, Coordinate};
    
    #[test]
    fn test_compute_zwd_mapping() {
        let el_sat = 45.0_f64.to_radians();
        let el_ref = 60.0_f64.to_radians();
        let zwd = 0.1;
        let (h_zwd, zwd_dd) = compute_zwd_mapping(el_sat, el_ref, zwd);
        
        let m_w_sat = 1.0 / el_sat.sin();
        let m_w_ref = 1.0 / el_ref.sin();
        assert!((h_zwd - (m_w_sat - m_w_ref)).abs() < 1e-9);
        assert!((zwd_dd - (h_zwd * zwd)).abs() < 1e-9);
        
        // Edge case: elevation < 0.001
        let (h_zwd_edge, _) = compute_zwd_mapping(0.0001, 0.0001, zwd);
        let max_val = 1.0 / 0.001_f64.sin();
        assert!((h_zwd_edge - (max_val - max_val)).abs() < 1e-9);
    }
    
    #[test]
    fn test_compute_geometric_dd() {
        let pos_apc = Vector3::new(10.0, 20.0, 30.0);
        let base_coord_vec = Vector3::new(10.0, 20.0, 30.0);
        let sat_vec_rov = Vector3::new(100.0, 200.0, 300.0);
        let ref_sat_vec_rov = Vector3::new(-100.0, -200.0, -300.0);
        
        let dd = compute_geometric_dd(
            pos_apc, base_coord_vec, 
            sat_vec_rov, ref_sat_vec_rov, 
            sat_vec_rov, ref_sat_vec_rov
        );
        // Because rover and base are at the same position, DD should exactly cancel out to 0.
        assert!(dd.abs() < 1e-9);
        
        let pos_apc2 = pos_apc + Vector3::new(1.0, 0.0, 0.0);
        let dd2 = compute_geometric_dd(
            pos_apc2, base_coord_vec, 
            sat_vec_rov, ref_sat_vec_rov, 
            sat_vec_rov, ref_sat_vec_rov
        );
        let d_rov_sat = (pos_apc2 - sat_vec_rov).norm();
        let d_rov_ref = (pos_apc2 - ref_sat_vec_rov).norm();
        let d_bas_sat = (base_coord_vec - sat_vec_rov).norm();
        let d_bas_ref = (base_coord_vec - ref_sat_vec_rov).norm();
        assert!((dd2 - ((d_rov_sat - d_rov_ref) - (d_bas_sat - d_bas_ref))).abs() < 1e-9);
    }
    
    #[test]
    fn test_compute_variance_factors() {
        let (var, ref_var) = compute_variance_factors(
            40.0, 45.0, 45.0_f64.to_radians(), 60.0_f64.to_radians(), 
            45.0_f64.to_radians(), 60.0_f64.to_radians(),
            100.0, 100.0
        );
        
        let expected_ref_var = gneiss_core::variance::observation_variance(45.0, 60.0_f64.to_radians(), 100.0, 100.0) +
            gneiss_core::variance::elevation_variance_scale(60.0_f64.to_radians());
        let expected_var = gneiss_core::variance::observation_variance(40.0, 45.0_f64.to_radians(), 100.0, 100.0) +
            gneiss_core::variance::elevation_variance_scale(45.0_f64.to_radians()) + expected_ref_var;
            
        assert!((ref_var - expected_ref_var).abs() < 1e-9);
        assert!((var - expected_var).abs() < 1e-9);
    }
    
    #[test]
    fn test_range_attitude_jacobian() {
        let lever = Vector3::new(1.0, 2.0, 3.0);
        let h_r = Vector3::new(4.0, 5.0, 6.0);
        let expected = lever.cross(&h_r);
        let result = range_attitude_jacobian(&lever, &h_r);
        assert!((result - expected).norm() < 1e-9);
    }
    
    #[test]
    fn test_doppler_attitude_jacobian() {
        let r_b_e = Matrix3::identity();
        let omega = Vector3::new(0.1, 0.2, 0.3);
        let lever = Vector3::new(1.0, 2.0, 3.0);
        let h_r = Vector3::new(4.0, 5.0, 6.0);
        
        let a = r_b_e * omega.cross(&lever);
        let expected = a.cross(&h_r);
        let result = doppler_attitude_jacobian(&r_b_e, &omega, &lever, &h_r);
        assert!((result - expected).norm() < 1e-9);
    }
    
    #[test]
    fn test_compute_phase_windup() {
        let updates = compute_phase_windup(
            Vector3::new(10.0, 20.0, 30.0),
            Vector3::new(10.0, 20.0, 30.0),
            Vector3::new(1e11, 0.0, 0.0),
            Vector3::new(100.0, 0.0, 0.0),
            Vector3::new(0.0, 100.0, 0.0),
            Vector3::new(100.0, 0.0, 0.0),
            Vector3::new(0.0, 100.0, 0.0),
            0.1, 0.2, 0.3, 0.4
        );
        // Phase windup function exists, we just need to ensure the wrapper passes values identically.
        // It's tested elsewhere, but we ensure our mapping doesn't mangle it.
        assert_eq!(updates.w_sat, gneiss_core::windup::phase_windup(Vector3::new(100.0, 0.0, 0.0), Vector3::new(1e11, 0.0, 0.0), Vector3::new(10.0, 20.0, 30.0), 0.1));
        assert_eq!(updates.w_ref, gneiss_core::windup::phase_windup(Vector3::new(0.0, 100.0, 0.0), Vector3::new(1e11, 0.0, 0.0), Vector3::new(10.0, 20.0, 30.0), 0.2));
        assert_eq!(updates.w_bas_sat, gneiss_core::windup::phase_windup(Vector3::new(100.0, 0.0, 0.0), Vector3::new(1e11, 0.0, 0.0), Vector3::new(10.0, 20.0, 30.0), 0.3));
        assert_eq!(updates.w_bas_ref, gneiss_core::windup::phase_windup(Vector3::new(0.0, 100.0, 0.0), Vector3::new(1e11, 0.0, 0.0), Vector3::new(10.0, 20.0, 30.0), 0.4));
    }
}
"""

with open("crates/gneiss-rtk/src/engine/measurement_math.rs", "r") as f:
    text = f.read()

# remove old #[cfg(test)] if any
if "#[cfg(test)]" in text:
    text = text[:text.find("#[cfg(test)]")]

text += test_code

with open("crates/gneiss-rtk/src/engine/measurement_math.rs", "w") as f:
    f.write(text)

test_code_2 = """
    #[test]
    fn test_get_sat_state() {
        let eph = Ephemeris::default();
        let pr = 20000000.0;
        let t_rx = GpsTime::new(2000, 100000.0);
        let rx_pos = Vector3::new(10.0, 20.0, 30.0);
        let (pos, vel) = get_sat_state(&eph, pr, t_rx, rx_pos);
        // Default ephemeris returns 0, so pos and vel should be zero for this dummy test
        assert!(pos.norm() < 1e-9);
        assert!(vel.norm() < 1e-9);
    }

    #[test]
    fn test_compute_atmospheric_delays() {
        let state_time = GpsTime::new(2000, 100000.0);
        let pos_apc = Vector3::new(10.0, 20.0, 30.0);
        let base_coord_vec = Vector3::new(10.0, 20.0, 30.0);
        let sat_vec_rov = Vector3::new(100.0, 200.0, 300.0);
        let ref_sat_vec_rov = Vector3::new(-100.0, -200.0, -300.0);
        let sat_vec_bas = Vector3::new(100.0, 200.0, 300.0);
        let ref_sat_vec_bas = Vector3::new(-100.0, -200.0, -300.0);
        
        let sat_f1 = 1575.42e6;
        let sat_f2 = 1227.60e6;
        let ref_f1 = 1575.42e6;
        let ref_f2 = 1227.60e6;
        
        let (tropo, iono1, iono2) = compute_atmospheric_delays(
            state_time, pos_apc, base_coord_vec, sat_vec_rov, ref_sat_vec_rov, sat_vec_bas, ref_sat_vec_bas,
            sat_f1, sat_f2, ref_f1, ref_f2
        );
        // Since rover and base are at same position and vectors are identical, double difference should be 0.
        assert!(tropo.abs() < 1e-9);
        assert!(iono1.abs() < 1e-9);
        assert!(iono2.abs() < 1e-9);
    }
"""

with open("crates/gneiss-rtk/src/engine/measurement_math.rs", "r") as f:
    text = f.read()

text = text.replace("}\n", test_code_2 + "\n}\n")

with open("crates/gneiss-rtk/src/engine/measurement_math.rs", "w") as f:
    f.write(text)

