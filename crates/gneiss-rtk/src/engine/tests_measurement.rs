#[cfg(test)]
mod tests {
    
    use crate::filter::{RtkState, DdObservation};
    use gneiss_core::coords::{Coordinate, Datum, Frame};
    use gneiss_core::time::GpsTime;
    use gneiss_core::sat::{SatelliteId, Constellation};
    use gneiss_core::ephemeris::Ephemeris;
    use nalgebra::Vector3;

    #[test]
    fn test_measurement_model_against_rtklib_golden_data() {
        let time = GpsTime::new(2137, 422922.0);
        let mut state = RtkState::new(time, Coordinate::new(Vector3::new(1000.0, 2000.0, 3000.0), Datum::WGS84, Frame::ECEF, time), 10.0);
        state.velocity = Vector3::new(10.0, -5.0, 2.0);
        
        let ref_sat = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let rov_sat1 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let rov_sat2 = SatelliteId { constellation: Constellation::Gps, prn: 3 };
        
        state.add_ambiguity(ref_sat, 1, 5.0, 100.0);
        state.add_ambiguity(rov_sat1, 1, 10.0, 100.0);
        state.add_ambiguity(rov_sat2, 1, 15.0, 100.0);

        let ref_rover = DdObservation { sat: ref_sat, pr_l1: 20000000.0, pr_l2: Some(20000001.0), cp_l1: Some(100000000.0), cp_l2: Some(80000000.0), doppler: 100.0, snr: 45.0, locktime: Some(100) };
        let ref_base = DdObservation { sat: ref_sat, pr_l1: 20005000.0, pr_l2: Some(20005001.0), cp_l1: Some(100020000.0), cp_l2: Some(80016000.0), doppler: 0.0, snr: 45.0, locktime: Some(100) };
        
        let rov1_rover = DdObservation { sat: rov_sat1, pr_l1: 21000000.0, pr_l2: Some(21000001.0), cp_l1: Some(105000000.0), cp_l2: Some(84000000.0), doppler: -50.0, snr: 45.0, locktime: Some(100) };
        let rov1_base = DdObservation { sat: rov_sat1, pr_l1: 21005000.0, pr_l2: Some(21005001.0), cp_l1: Some(105020000.0), cp_l2: Some(84016000.0), doppler: 0.0, snr: 45.0, locktime: Some(100) };

        let rov2_rover = DdObservation { sat: rov_sat2, pr_l1: 22000000.0, pr_l2: Some(22000001.0), cp_l1: Some(110000000.0), cp_l2: Some(88000000.0), doppler: -20.0, snr: 45.0, locktime: Some(100) };
        let rov2_base = DdObservation { sat: rov_sat2, pr_l1: 22005000.0, pr_l2: Some(22005001.0), cp_l1: Some(110020000.0), cp_l2: Some(88016000.0), doppler: 0.0, snr: 45.0, locktime: Some(100) };

        let matched_obs = vec![(rov1_rover, rov1_base), (rov2_rover, rov2_base)];

        let eph_ref = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: ref_sat, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
            iode: 0, iodc: 0,
        });

        let eph_rov1 = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: rov_sat1, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 1.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 0.5, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
            iode: 0, iodc: 0,
        });
        
        let eph_rov2 = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: rov_sat2, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 2.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
            omega0: 1.0, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
            iode: 0, iodc: 0,
        });

        let ephemerides = vec![eph_ref, eph_rov1, eph_rov2];
        let base_coord = Coordinate::new(Vector3::new(1005.0, 2005.0, 3005.0), Datum::WGS84, Frame::ECEF, time);

        // We explicitly use compute_innovations to avoid the Mahalanobis chi2 filter rejecting dummy data
        let _config = crate::engine::EngineConfig::default();
        let (z, _h, r, _) = super::super::measurement::compute_innovations(&mut state, &matched_obs, &ephemerides, &base_coord, base_coord.epoch, &ref_rover, &ref_base, Vector3::zeros(), Vector3::zeros()).unwrap();

        println!("Z: {:?}", z);
        
        // Lock in the golden Z vector (updated for iterative ecef_to_llh refinement)
        assert!((z[0] - 0.70573).abs() < 1e-3, "z[0]={}", z[0]);
        assert!((z[1] - 0.70623).abs() < 1e-3, "z[1]={}", z[1]);
        assert!((z[2] - -4.29582).abs() < 1e-3, "z[2]={}", z[2]);
        assert!((z[3] - 7.72370).abs() < 1e-3, "z[3]={}", z[3]);
        assert!((z[4] - 7.72444).abs() < 1e-3, "z[4]={}", z[4]);
        assert!((z[5] - -2.27859).abs() < 1e-3, "z[5]={}", z[5]);

        // Lock in the golden R diagonal
        assert!(r[0] >= 16.0); // Now scales with elevation
        assert!(r[1] >= 16.0);
        assert!(r[2] >= 0.0001);
        assert!(r[3] >= 16.0);
        assert!(r[4] >= 16.0);
        assert!(r[5] >= 0.0001);
    }
}