#[cfg(test)]
mod tests {
    use super::*;
    use gneiss_core::sat::Constellation;
    use gneiss_core::coords::{Datum, Frame};

    #[test]
    fn test_resolve_ambiguities_multi_constellation() {
        let time = GpsTime::new(2137, 422922.0);
        let initial_pos = Coordinate::new(Vector3::new(1000.0, 2000.0, 3000.0), Datum::WGS84, Frame::ECEF, time);
        let mut state = RtkState::new(time, initial_pos, 10.0);
        state.epoch_count() = 100;

        let gps_ref = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let gps_rov1 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let gps_rov2 = SatelliteId { constellation: Constellation::Gps, prn: 3 };
        let gps_rov3 = SatelliteId { constellation: Constellation::Gps, prn: 4 };
        let gal_ref = SatelliteId { constellation: Constellation::Galileo, prn: 10 };
        let gal_rov = SatelliteId { constellation: Constellation::Galileo, prn: 11 };

        let lam = 0.19029367279836487;
        state.add_ambiguity(gps_ref, 1, 10.1 * lam, 0.0001);
        state.add_ambiguity(gps_rov1, 1, 15.2 * lam, 0.0001);
        state.add_ambiguity(gps_rov2, 1, 20.3 * lam, 0.0001);
        state.add_ambiguity(gps_rov3, 1, 25.4 * lam, 0.0001);
        state.add_ambiguity(gal_ref, 1, 30.5 * lam, 0.0001);
        state.add_ambiguity(gal_rov, 1, 36.6 * lam, 0.0001);

        for &(sat, freq) in &[(gps_ref, 1), (gps_rov1, 1), (gps_rov2, 1), (gps_rov3, 1), (gal_ref, 1), (gal_rov, 1)] {
            state.locktimes_mut().insert((sat, freq), 100);
        }

        use gneiss_core::ephemeris::{Ephemeris, GpsEphemeris, GalileoEphemeris};
        let ephemerides = vec![
            Ephemeris::Gps(GpsEphemeris { 
                sat: gps_ref, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
                omega0: 0.0, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
                iode: 0, iodc: 0,
            }),
            Ephemeris::Gps(GpsEphemeris { 
                sat: gps_rov1, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
                omega0: 0.1, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
                iode: 0, iodc: 0,
            }),
            Ephemeris::Gps(GpsEphemeris { 
                sat: gps_rov2, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
                omega0: 0.2, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
                iode: 0, iodc: 0,
            }),
            Ephemeris::Gps(GpsEphemeris { 
                sat: gps_rov3, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.01, sqrt_a: 5153.6, delta_n: 0.0,
                omega0: 0.3, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, tgd: 0.0,
                iode: 0, iodc: 0,
            }),
            Ephemeris::Galileo(GalileoEphemeris { 
                sat: gal_ref, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.01, sqrt_a: 5440.6, delta_n: 0.0,
                omega0: 0.0, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, bgd_e1_e5a: 0.0,
                iod_nav: 0,
            }),
            Ephemeris::Galileo(GalileoEphemeris { 
                sat: gal_rov, toe: time, toc: time, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
                m0: 0.0, e: 0.01, sqrt_a: 5440.6, delta_n: 0.0,
                omega0: 0.5, omega_dot: 0.0, i0: 1.0, idot: 0.0, omega: 0.0, bgd_e1_e5a: 0.0,
                iod_nav: 0,
            }),
        ];

        state.resolve_ambiguities(&ephemerides, 4, 5, 3, 3.0).expect("AR should run");

        assert!(state.is_fixed(), "Should achieve fix with multi-constellation support");
        
        let idx_ref = state.ambiguity_keys().iter().position(|&(s, f)| s == gps_ref && f == 1).unwrap();
        let idx_rov = state.ambiguity_keys().iter().position(|&(s, f)| s == gps_rov1 && f == 1).unwrap();
        let dd_gps = (state.ambiguities()[idx_rov] - state.ambiguities()[idx_ref]) / lam;
        assert!((dd_gps.round() - 5.0).abs() < 1e-6);

        let idx_ref_gal = state.ambiguity_keys().iter().position(|&(s, f)| s == gal_ref && f == 1).unwrap();
        let idx_rov_gal = state.ambiguity_keys().iter().position(|&(s, f)| s == gal_rov && f == 1).unwrap();
        let dd_gal = (state.ambiguities()[idx_rov_gal] - state.ambiguities()[idx_ref_gal]) / lam;
        assert!((dd_gal.round() - 6.0).abs() < 1e-6);
    }

    #[test]
    fn test_double_difference_eliminates_clocks() {
        let sat_ref = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let sat_a = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let true_r_rover_ref = 20_000_000.0;
        let true_r_rover_a   = 21_000_000.0;
        let true_r_base_ref  = 20_005_000.0;
        let true_r_base_a    = 21_004_000.0;
        let rover_clk = 300.0; let base_clk = -150.0;
        let sat_ref_clk = 1000.0; let sat_a_clk = -500.0;
        let rover_ref_obs = DdObservation { sat: sat_ref, pr_l1: true_r_rover_ref + rover_clk - sat_ref_clk, pr_l2: None, cp_l1: 0.0, cp_l2: None, doppler: 0.0, snr: 45.0, locktime: 1000 };
        let rover_a_obs = DdObservation { sat: sat_a, pr_l1: true_r_rover_a + rover_clk - sat_a_clk, pr_l2: None, cp_l1: 0.0, cp_l2: None, doppler: 0.0, snr: 45.0, locktime: 1000 };
        let base_ref_obs = DdObservation { sat: sat_ref, pr_l1: true_r_base_ref + base_clk - sat_ref_clk, pr_l2: None, cp_l1: 0.0, cp_l2: None, doppler: 0.0, snr: 45.0, locktime: 1000 };
        let base_a_obs = DdObservation { sat: sat_a, pr_l1: true_r_base_a + base_clk - sat_a_clk, pr_l2: None, cp_l1: 0.0, cp_l2: None, doppler: 0.0, snr: 45.0, locktime: 1000 };
        let f1 = 1575.42e6;
        let f2 = 1227.60e6;
        let dd = compute_double_difference(&rover_ref_obs, &rover_a_obs, &base_ref_obs, &base_a_obs, f1, f2, f1, f2);
        assert!((dd - 1000.0).abs() < 1e-6);
    }
}
