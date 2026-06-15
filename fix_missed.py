import re

with open("crates/gneiss-rtk/src/engine/measurement_math.rs", "r") as f:
    text = f.read()

new_test = """
    #[test]
    fn test_compute_atmospheric_delays_catches_tropo_mutation() {
        let state_time = GpsTime::new(2137, 422922.0);
        
        let pos_apc = Vector3::new(6378137.0, 0.0, 0.0);
        let base_coord_vec = Vector3::new(6378137.0, 1000.0, 0.0);
        
        let sat_vec_rov = Vector3::new(26000000.0, 0.0, 0.0);
        let ref_sat_vec_rov = Vector3::new(26000000.0, 10000000.0, 0.0);
        
        let sat_vec_bas = Vector3::new(26000000.0, 0.0, 0.0);
        let ref_sat_vec_bas = Vector3::new(26000000.0, 10000000.0, 0.0);

        let (tropo_dd, _, _) = compute_atmospheric_delays(
            state_time, pos_apc, base_coord_vec, 
            sat_vec_rov, ref_sat_vec_rov, sat_vec_bas, ref_sat_vec_bas, 
            1.0, 1.0, 1.0, 1.0
        );
        
        assert!((tropo_dd - 0.0).abs() > 1e-6);
        // Specifically, let's just assert that it exactly matches the logic:
        let tropo_params = gneiss_core::atmosphere::TropoParams::default();
        let base_llh = gneiss_core::coords::ecef_to_llh(base_coord_vec);
        let rov_llh = gneiss_core::coords::ecef_to_llh(pos_apc);
        let (_, el_rov_sat) = gneiss_core::coords::az_el(rov_llh, pos_apc, sat_vec_rov);
        let (_, el_rov_ref) = gneiss_core::coords::az_el(rov_llh, pos_apc, ref_sat_vec_rov);
        let (_, el_bas_sat) = gneiss_core::coords::az_el(base_llh, base_coord_vec, sat_vec_bas);
        let (_, el_bas_ref) = gneiss_core::coords::az_el(base_llh, base_coord_vec, ref_sat_vec_bas);

        let t1 = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&tropo_params, rov_llh, el_rov_sat);
        let t2 = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&tropo_params, rov_llh, el_rov_ref);
        let t3 = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&tropo_params, base_llh, el_bas_sat);
        let t4 = gneiss_core::atmosphere::AtmosphereModel::tropo_rtklib_saastamoinen(&tropo_params, base_llh, el_bas_ref);
        
        let expected_dd = (t1 - t2) - (t3 - t4);
        assert!((tropo_dd - expected_dd).abs() < 1e-9);
    }
"""

text = re.sub(r'#\[test\]\s+fn test_compute_atmospheric_delays_catches_tropo_mutation\(\) \{.*?(?=\n    }\n)\n    }', new_test.strip(), text, flags=re.DOTALL)

with open("crates/gneiss-rtk/src/engine/measurement_math.rs", "w") as f:
    f.write(text)

