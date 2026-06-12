    #[test]
    fn test_spp_not_enough_measurements() {
        use gneiss_core::time::GpsTime;
        use gneiss_core::obs::EpochObs;

        let t = GpsTime::new(2000, 100000.0);
        let epoch = EpochObs {
            time: t,
            satellites: vec![],
        };
        let config = SppConfig::default();
        let res = compute_spp(&epoch, &[], None, &config, None);
        assert_eq!(res.unwrap_err(), SppError::NotEnoughMeasurements);
    }

    #[test]
    fn test_spp_wnlls_step_not_enough_measurements() {
        use gneiss_core::time::GpsTime;
        use gneiss_core::sat::{Constellation, SatelliteId};

        let t = GpsTime::new(2000, 100000.0);
        let m1 = SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20000000.0, snr: 45.0, doppler: 0.0, time: t,
            eph: Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: SatelliteId { constellation: Constellation::Gps, prn: 1 }, toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0, sqrt_a: 5153.6, delta_n: 0.0,
                omega0: 0.0, omega_dot: 0.0, i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
            }),
        };

        let config = SppConfig::default();
        let state = SppState::new(Coordinate::new(nalgebra::Vector3::zeros(), Datum::WGS84, Frame::ECEF, t), 0.0, 0.0, 0.0, 0.0);
        
        let res = spp_wnlls_step(&state, &[m1], None, &config);
        assert_eq!(res.unwrap_err(), SppError::NotEnoughMeasurements);
    }
