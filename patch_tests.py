with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "r") as f:
    content = f.read()

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
        
        let (pos0, vel0) = get_sat_state(&eph, 0.0, time, rx_pos);
        assert!((pos.x - pos0.x).abs() > 0.0);
    }
}
"""

content = content.rstrip().rsplit('}', 1)[0] + test_str
with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "w") as f:
    f.write(content)

