with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "r") as f:
    content = f.read()

content = content.replace("    }\n}\n", "    }\n")

test_str = """
    #[test]
    fn test_compute_dd() {
        use crate::engine::measurement::{compute_dd_pseudorange, compute_dd_carrier_phase, compute_dd_doppler, apply_phase_windup};
        let rov_sat = crate::filter::DdObservation {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            pr_l1: 20000000.0,
            pr_l2: Some(20000010.0),
            cp_l1: Some(100000000.0),
            cp_l2: Some(80000000.0),
            doppler: 100.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let base_sat = crate::filter::DdObservation {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            pr_l1: 20005000.0,
            pr_l2: Some(20005010.0),
            cp_l1: Some(100020000.0),
            cp_l2: Some(80016000.0),
            doppler: 50.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let rov_ref = crate::filter::DdObservation {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 2 },
            pr_l1: 21000000.0,
            pr_l2: Some(21000010.0),
            cp_l1: Some(105000000.0),
            cp_l2: Some(84000000.0),
            doppler: -100.0,
            snr: 45.0,
            locktime: Some(100),
        };
        let ref_base = crate::filter::DdObservation {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 2 },
            pr_l1: 21005000.0,
            pr_l2: Some(21005010.0),
            cp_l1: Some(105020000.0),
            cp_l2: Some(84016000.0),
            doppler: -50.0,
            snr: 45.0,
            locktime: Some(100),
        };

        let f1 = 1575.42e6;
        let f2 = 1227.60e6;

        let (pr1, pr2) = compute_dd_pseudorange(&rov_sat, &base_sat, &rov_ref, &ref_base);
        assert!((pr1.unwrap() - 0.0).abs() < 1e-6);
        assert!((pr2.unwrap() - 0.0).abs() < 1e-6);

        let (cp1, cp2) = compute_dd_carrier_phase(&rov_sat, &base_sat, &rov_ref, &ref_base, f1, f2);
        assert!((cp1.unwrap() - 0.0).abs() < 1e-6);
        assert!((cp2.unwrap() - 0.0).abs() < 1e-6);

        let dop = compute_dd_doppler(&rov_sat, &base_sat, &rov_ref, &ref_base);
        assert!((dop.unwrap() - 100.0).abs() < 1e-6);

        let time = GpsTime::new(2137, 422922.0);
        let rx_pos = Vector3::new(1000.0, 2000.0, 3000.0);
        let sat_pos = Vector3::new(20000000.0, 1000000.0, 0.0);
        let mut prev_windup = 0.0;
        let wu = apply_phase_windup(time, rx_pos, sat_pos, &mut prev_windup);
        assert!(wu.abs() > 0.0);
    }
}
"""

content += test_str

with open("crates/gneiss-rtk/src/engine/tests_measurement.rs", "w") as f:
    f.write(content)
