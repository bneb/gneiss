with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "r") as f:
    content = f.read()

# remove trailing }
content = content.replace("    }\n}\n", "    }\n")

test_str = """
    #[test]
    fn test_get_sat_state() {
        use crate::engine::measurement::get_sat_state;
        let time = GpsTime::new(2137, 422922.0);
        let rx_pos = Vector3::new(1000.0, 2000.0, 3000.0);
        let sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 1.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
            iode: 0, iodc: 0,
        });

        let pr = 20000000.0;
        let (pos, vel) = get_sat_state(&eph, pr, time, rx_pos);
        
        // Assert non-zero output
        assert!(pos.norm() > 10000.0);
        assert!(vel.norm() > 10.0);
        
        let (pos0, _vel0) = get_sat_state(&eph, 0.0, time, rx_pos);
        assert!((pos.x - pos0.x).abs() > 0.0);
    }

    #[test]
    fn test_compute_atmospheric_delays() {
        use crate::engine::measurement::compute_atmospheric_delays;
        let time = GpsTime::new(2137, 422922.0);
        let pos_apc = Vector3::new(1000.0, 2000.0, 3000.0);
        let base_coord = Vector3::new(1005.0, 2005.0, 3005.0);
        
        let sat_vec_rov = Vector3::new(20000000.0, 1000000.0, 0.0);
        let ref_sat_vec_rov = Vector3::new(0.0, 20000000.0, 1000000.0);
        
        let sat_vec_bas = Vector3::new(20000000.0, 1000000.0, 0.0);
        let ref_sat_vec_bas = Vector3::new(0.0, 20000000.0, 1000000.0);
        
        let f1 = 1575.42e6;
        let f2 = 1227.60e6;
        
        let (tropo, iono1, iono2) = compute_atmospheric_delays(time, pos_apc, base_coord, sat_vec_rov, ref_sat_vec_rov, sat_vec_bas, ref_sat_vec_bas, f1, f2, f1, f2);
        
        assert!((tropo - -0.0039775f64).abs() < 1e-4);
        assert!((iono1 - -0.015112f64).abs() < 1e-4);
        assert!((iono2 - -0.02488f64).abs() < 1e-4);
    }
}
"""

content += test_str

with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "w") as f:
    f.write(content)
