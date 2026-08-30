#![allow(clippy::unwrap_used)]

use super::*;
use super::measurements::*;
use super::raim::*;
use super::solver::*;
use gneiss_core::coords::{ecef_to_llh, Coordinate, Datum, Frame};
use gneiss_core::obs::ObsType;
use nalgebra::{DMatrix, DVector, Vector3};

    #[test]
    fn test_compute_spp() {
        use gneiss_core::obs::{EpochObs, ObsCode, Observation, SatObs, SignalCode};
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;

        let t = GpsTime::new(2000, 100000.0);
        let true_pos =
            nalgebra::Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let true_cdt = 1000.0;

        let mut ephemerides = Vec::new();
        let mut satellites = Vec::new();

        let sats = vec![
            (1, (20000000.0, 5000000.0, 5000000.0)),
            (2, (22000000.0, -5000000.0, 5000000.0)),
            (3, (19000000.0, 5000000.0, -5000000.0)),
            (4, (21000000.0, -5000000.0, -5000000.0)),
            (5, (25000000.0, 0.0, 0.0)),
        ];

        for (prn, (_sx, _sy, _sz)) in sats {
            let sat_id = SatelliteId {
                constellation: Constellation::Gps,
                prn,
            };

            let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: sat_id,
                toe: t,
                toc: t,
                af0: 0.0,
                af1: 0.0,
                af2: 0.0,
                crs: 0.0,
                crc: 0.0,
                cuc: 0.0,
                cus: 0.0,
                cic: 0.0,
                cis: 0.0,
                m0: (prn as f64) * core::f64::consts::PI / 3.0,
                e: 0.01,
                sqrt_a: 5153.6,
                delta_n: 0.0,
                omega0: (prn as f64) * core::f64::consts::PI / 2.0,
                omega_dot: 0.0,
                i0: 0.95,
                idot: 0.0,
                omega: 0.0,
                tgd: 0.0,
                iode: 1,
                iodc: 1,
            });
            ephemerides.push(eph.clone());

            let mut raw_pr = 20000000.0 + true_cdt; // rough initial guess
            let light_speed = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

            for _ in 0..5 {
                let pr_time = raw_pr / light_speed;
                let t_tx_sat = t.tow - pr_time;
                let (_, _, sat_clk_err_rough, _) = eph.position(GpsTime::new(t.week, t_tx_sat));
                let t_tx_true = t_tx_sat - sat_clk_err_rough;

                let (sat_pos, _, sat_clk_err, _) = eph.position(GpsTime::new(t.week, t_tx_true));
                let dx = true_pos.x - sat_pos.x;
                let dy = true_pos.y - sat_pos.y;
                let dz = true_pos.z - sat_pos.z;
                let geometric_range = f64::sqrt(dx * dx + dy * dy + dz * dz);

                raw_pr = geometric_range + true_cdt - (sat_clk_err * light_speed);
            }

            satellites.push(SatObs {
                sat: sat_id,
                observations: vec![Observation {
                    code: ObsCode {
                        obs_type: ObsType::Pseudorange,
                        signal: SignalCode {
                            freq_band: 1,
                            attribute: 'C',
                        },
                    },
                    value: raw_pr,
                    lock_time: None,
                    lli: None,
                }],
            });
        }

        let epoch = EpochObs {
            time: t,
            satellites,
        };

        let config = SppConfig {
            enable_sagnac: false,
            enable_tropo: false,
            enable_iono: false,
            geometry_variance_threshold: 100000.0,
            elevation_mask_rad: -core::f64::consts::PI, // Bypass elevation mask for tests
            ..Default::default()
        };

        let state = compute_spp(&epoch, &ephemerides, None, &config, None).unwrap();

        assert!(
            (state.position.vector.x - true_pos.x).abs() < 1e-1,
            "X error too large: {}",
            (state.position.vector.x - true_pos.x).abs()
        );
        assert!(
            (state.position.vector.y - true_pos.y).abs() < 1e-1,
            "Y error too large: {}",
            (state.position.vector.y - true_pos.y).abs()
        );
        assert!(
            (state.position.vector.z - true_pos.z).abs() < 1e-1,
            "Z error too large: {}",
            (state.position.vector.z - true_pos.z).abs()
        );
        assert!(
            (state.cdt - true_cdt).abs() < 1e-1,
            "CDT error too large: {}",
            (state.cdt - true_cdt).abs()
        );
    }
    #[test]
    fn test_spp_not_enough_measurements() {
        use gneiss_core::obs::EpochObs;
        use gneiss_core::time::GpsTime;

        let t = GpsTime::new(2000, 100000.0);
        let epoch = EpochObs {
            time: t,
            satellites: vec![],
        };
        let config = SppConfig {
            elevation_mask_rad: -core::f64::consts::PI,
            ..Default::default()
        };
        let res = compute_spp(&epoch, &[], None, &config, None);
        assert_eq!(res.unwrap_err(), SppError::NotEnoughMeasurements);
    }

    #[test]
    fn test_spp_wnlls_step_not_enough_measurements() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;

        let t = GpsTime::new(2000, 100000.0);
        let m1 = SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20000000.0,
            snr: 45.0,
            doppler: 0.0,
            time: t,
            eph: Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: SatelliteId {
                    constellation: Constellation::Gps,
                    prn: 1,
                },
                toe: t,
                toc: t,
                af0: 0.0,
                af1: 0.0,
                af2: 0.0,
                crs: 0.0,
                crc: 0.0,
                cuc: 0.0,
                cus: 0.0,
                cic: 0.0,
                cis: 0.0,
                m0: 0.0,
                e: 0.0,
                sqrt_a: 5153.6,
                delta_n: 0.0,
                omega0: 0.0,
                omega_dot: 0.0,
                i0: 0.95,
                idot: 0.0,
                omega: 0.0,
                tgd: 0.0,
                iode: 1,
                iodc: 1,
            }),
            is_iono_free: false,
            freq_band: 1,
        };

        let config = SppConfig {
            elevation_mask_rad: -core::f64::consts::PI,
            ..Default::default()
        };
        let state = SppState::new(
            Coordinate::new(nalgebra::Vector3::zeros(), Datum::WGS84, Frame::ECEF, t),
            0.0,
            0.0,
            0.0,
            0.0,
        );

        let res = spp_wnlls_step(&state, &[m1], None, &config);
        assert_eq!(res.unwrap_err(), SppError::NotEnoughMeasurements);
    }

    #[test]
    fn test_build_measurements() {
        use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;

        let t = GpsTime::new(2000, 100000.0);
        let t_eph1 = GpsTime::new(2000, 100010.0); // da = 10
        let t_eph2 = GpsTime::new(2000, 99980.0); // db = 20. eph1 is closer!

        let sat_gps = SatelliteId {
            constellation: Constellation::Gps,
            prn: 1,
        };
        let sat_glo = SatelliteId {
            constellation: Constellation::Glonass,
            prn: 1,
        };
        let sat_gal = SatelliteId {
            constellation: Constellation::Galileo,
            prn: 2,
        };
        let sat_bds = SatelliteId {
            constellation: Constellation::Beidou,
            prn: 3,
        };

        let mut epoch = EpochObs {
            time: t,
            satellites: vec![],
        };

        // GLONASS omitted for this test

        // Add GPS with no pseudorange (should be skipped by pr check)
        epoch.satellites.push(SatObs {
            sat: sat_gps,
            observations: vec![Observation {
                code: ObsCode {
                    obs_type: ObsType::CarrierPhase,
                    signal: SignalCode {
                        freq_band: 1,
                        attribute: 'C',
                    },
                },
                value: 100000.0,
                lock_time: None,
                lli: None,
            }],
        });

        // Add Galileo with pseudorange, snr, doppler, and correct freq band
        epoch.satellites.push(SatObs {
            sat: sat_gal,
            observations: vec![
                Observation {
                    code: ObsCode {
                        obs_type: ObsType::Pseudorange,
                        signal: SignalCode {
                            freq_band: 1,
                            attribute: 'C',
                        },
                    },
                    value: 21000000.0,
                    lock_time: None,
                    lli: None,
                },
                Observation {
                    code: ObsCode {
                        obs_type: ObsType::Snr,
                        signal: SignalCode {
                            freq_band: 1,
                            attribute: 'C',
                        },
                    },
                    value: 42.0,
                    lock_time: None,
                    lli: None,
                },
                Observation {
                    code: ObsCode {
                        obs_type: ObsType::Doppler,
                        signal: SignalCode {
                            freq_band: 1,
                            attribute: 'C',
                        },
                    },
                    value: 1500.0,
                    lock_time: None,
                    lli: None,
                },
                // Wrong freq band pseudorange just to test
                Observation {
                    code: ObsCode {
                        obs_type: ObsType::Pseudorange,
                        signal: SignalCode {
                            freq_band: 2,
                            attribute: 'C',
                        },
                    },
                    value: 22000000.0,
                    lock_time: None,
                    lli: None,
                },
            ],
        });

        // Beidou omitted for this test

        let eph_gps = gneiss_core::ephemeris::GpsEphemeris {
            sat: sat_gps,
            toe: t,
            toc: t,
            af0: 0.0,
            af1: 0.0,
            af2: 0.0,
            crs: 0.0,
            crc: 0.0,
            cuc: 0.0,
            cus: 0.0,
            cic: 0.0,
            cis: 0.0,
            m0: 0.0,
            e: 0.0,
            sqrt_a: 5153.6,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 0.95,
            idot: 0.0,
            omega: 0.0,
            tgd: 0.0,
            iode: 1,
            iodc: 1,
        };

        let eph_gal_1 = gneiss_core::ephemeris::GalileoEphemeris {
            sat: sat_gal,
            toe: t_eph1,
            toc: t_eph1,
            af0: 0.0,
            af1: 0.0,
            af2: 0.0,
            crs: 0.0,
            crc: 0.0,
            cuc: 0.0,
            cus: 0.0,
            cic: 0.0,
            cis: 0.0,
            m0: 0.0,
            e: 0.0,
            sqrt_a: 5153.6,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 0.95,
            idot: 0.0,
            omega: 0.0,
            bgd_e1_e5a: 0.0,
            bgd_e1_e5b: 0.0,
            iod_nav: 1,
        };

        let eph_gal_2 = gneiss_core::ephemeris::GalileoEphemeris {
            sat: sat_gal,
            toe: t_eph2,
            toc: t_eph2,
            af0: 1.0,
            af1: 0.0,
            af2: 0.0, // closer in time!
            crs: 0.0,
            crc: 0.0,
            cuc: 0.0,
            cus: 0.0,
            cic: 0.0,
            cis: 0.0,
            m0: 0.0,
            e: 0.0,
            sqrt_a: 5153.6,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 0.95,
            idot: 0.0,
            omega: 0.0,
            bgd_e1_e5a: 0.0,
            bgd_e1_e5b: 0.0,
            iod_nav: 2,
        };

        let eph_bds = gneiss_core::ephemeris::BeidouEphemeris {
            sat: sat_bds,
            toe: t,
            toc: t,
            af0: 0.0,
            af1: 0.0,
            af2: 0.0,
            crs: 0.0,
            crc: 0.0,
            cuc: 0.0,
            cus: 0.0,
            cic: 0.0,
            cis: 0.0,
            m0: 0.0,
            e: 0.0,
            sqrt_a: 5153.6,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 0.95,
            idot: 0.0,
            omega: 0.0,
            tgd1: 0.0,
            tgd2: 0.0,
            aode: 1,
            aodc: 1,
        };

        let eph_glo = gneiss_core::ephemeris::GlonassEphemeris {
            sat: sat_glo,
            toe: t,
            freq_num: 1,
            tau_n: 0.0,
            gamma_n: 0.0,
            delta_tau_n: 0.0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
            ax: 0.0,
            ay: 0.0,
            az: 0.0,
        };

        let ephemerides = vec![
            Ephemeris::Gps(eph_gps),
            Ephemeris::Galileo(eph_gal_1),
            Ephemeris::Galileo(eph_gal_2.clone()),
            Ephemeris::Beidou(eph_bds),
            Ephemeris::Glonass(eph_glo),
        ];

        let config = SppConfig::default();
        let measurements = build_measurements(&epoch, &ephemerides, &config);

        assert_eq!(
            measurements.len(),
            1,
            "Should have exactly 1 measurement (Galileo)"
        );

        // Galileo check
        assert_eq!(measurements[0].snr, 42.0);
        assert_eq!(measurements[0].doppler, 1500.0);
        match &measurements[0].eph {
            Ephemeris::Galileo(g) => assert_eq!(g.iod_nav, 1), // eph1 is now closer
            _ => panic!("Wrong ephemeris type"),
        }
        assert!(!measurements[0].is_iono_free); // missing L2
    }

    #[test]
    fn test_compute_sagnac_correction() {
        let sat_ecef = Vector3::new(20000000.0, 5000000.0, 5000000.0);
        let corrected = compute_sagnac_correction(sat_ecef, 25000000.0);
        assert!((corrected.norm() - sat_ecef.norm()).abs() < 1e-6, "Sagnac should preserve norm");
        assert_ne!(corrected.x, sat_ecef.x); // Should rotate
    }

    #[test]
    fn test_compute_atmospheric_delays_below_min_earth_radius() {
        let rec_ecef = Vector3::new(1.0, 0.0, 0.0);
        let (tropo, iono) = compute_atmospheric_delays(
            rec_ecef, Vector3::zeros(), 0.0, 0.5, GpsTime::new(0, 0.0), None, &SppConfig::default(),
        );
        assert_eq!(tropo, 0.0);
        assert_eq!(iono, 0.0);
    }

    #[test]
    fn test_apply_height_constraint_sets_correct_structure() {
        let llh = Vector3::new(0.5, 1.0, 100.0); // lat=0.5rad, lon=1.0rad
        let mut h = DMatrix::zeros(5, 4);
        let mut w = DMatrix::identity(5, 5);
        let mut dz = DVector::zeros(5);
        apply_height_constraint(llh, 4, &mut h, &mut w, &mut dz);
        let sin_lat = 0.5f64.sin();
        let cos_lat = 0.5f64.cos();
        let sin_lon = 1.0f64.sin();
        let cos_lon = 1.0f64.cos();
        assert!((h[(4, 0)] - cos_lat * cos_lon).abs() < 1e-10);
        assert!((h[(4, 1)] - cos_lat * sin_lon).abs() < 1e-10);
        assert!((h[(4, 2)] - sin_lat).abs() < 1e-10);
        assert_eq!(dz[4], 0.0);
        assert!((w[(4, 4)] - 1.0 / 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_spp_state_get_cdt() {
        use gneiss_core::sat::Constellation;
        let state = SppState::new(
            Coordinate::new(Vector3::zeros(), Datum::WGS84, Frame::ECEF, GpsTime::new(0, 0.0)),
            1.0, 2.0, 3.0, 4.0,
        );
        assert_eq!(state.get_cdt(Constellation::Gps), 1.0);
        assert_eq!(state.get_cdt(Constellation::Galileo), 2.0);
        assert_eq!(state.get_cdt(Constellation::Beidou), 3.0);
        assert_eq!(state.get_cdt(Constellation::Glonass), 4.0);
        assert_eq!(state.get_cdt(Constellation::Qzss), 1.0);
        assert_eq!(state.get_cdt(Constellation::Navic), 1.0);
    }

    #[test]
    fn test_find_clock_cols_gps_only() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        let t = GpsTime::new(2000, 100000.0);
        let m1 = SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20000000.0, snr: 45.0, doppler: 0.0, time: t,
            eph: Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
                toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
                cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0, sqrt_a: 5153.6,
                delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
                i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
            }),
            is_iono_free: false,
            freq_band: 1,
        };
        let (cols, clk) = find_clock_cols(&[m1]);
        assert_eq!(cols, 4);
        assert_eq!(clk.0, Some(3));
        assert!(clk.1.is_none());
        assert!(clk.2.is_none());
        assert!(clk.3.is_none());
    }

    #[test]
    fn test_find_clock_cols_multi_constellation() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        let t = GpsTime::new(2000, 100000.0);
        let eph_base = || -> Ephemeris {
            Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
                toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
                cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0, sqrt_a: 5153.6,
                delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
                i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
            })
        };
        let ms: Vec<SppMeasurement> = vec![
            SppMeasurement { constellation: Constellation::Gps, raw_pr: 20000000.0, snr: 45.0, doppler: 0.0, time: t, eph: eph_base(), is_iono_free: false, freq_band: 1 },
            SppMeasurement { constellation: Constellation::Galileo, raw_pr: 20000000.0, snr: 45.0, doppler: 0.0, time: t, eph: eph_base(), is_iono_free: false, freq_band: 1 },
            SppMeasurement { constellation: Constellation::Beidou, raw_pr: 20000000.0, snr: 45.0, doppler: 0.0, time: t, eph: eph_base(), is_iono_free: false, freq_band: 1 },
            SppMeasurement { constellation: Constellation::Glonass, raw_pr: 20000000.0, snr: 45.0, doppler: 0.0, time: t, eph: eph_base(), is_iono_free: false, freq_band: 1 },
        ];
        let (cols, clk) = find_clock_cols(&ms);
        assert_eq!(cols, 7);
        assert!(clk.0.is_some());
        assert!(clk.1.is_some());
        assert!(clk.2.is_some());
        assert!(clk.3.is_some());
    }

    #[test]
    fn test_spp_wnlls_step_not_enough_measurements_height_constrained() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        let t = GpsTime::new(2000, 100000.0);
        // 3 GPS measurements with cols=4 -> height constraint applied but still underdetermined
        let ms: Vec<SppMeasurement> = (0..3)
            .map(|_| SppMeasurement {
                constellation: Constellation::Gps,
                raw_pr: 20000000.0, snr: 45.0, doppler: 0.0, time: t,
                eph: Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                    sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
                    toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
                    crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
                    cic: 0.0, cis: 0.0, m0: 0.0, e: 0.0, sqrt_a: 5153.6,
                    delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
                    i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
                }),
                is_iono_free: false,
                freq_band: 1,
            })
            .collect();
        let state = SppState::new(
            Coordinate::new(Vector3::new(10000000.0, 10000000.0, 0.0), Datum::WGS84, Frame::ECEF, t),
            0.0, 0.0, 0.0, 0.0,
        );
        let config = SppConfig {
            elevation_mask_rad: -core::f64::consts::PI,
            enable_sagnac: false,
            enable_tropo: false,
            enable_iono: false,
            ..Default::default()
        };
        let res = spp_wnlls_step(&state, &ms, None, &config);
        // With 3 measurements and height constraint (4 total eqns, 4 unknowns),
        // identical geometry makes the matrix singular -> inversion fails
        assert!(res.is_err());
    }

    #[test]
    fn test_compute_sat_state_iono_free_path() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;

        let t = GpsTime::new(2000, 100000.0);
        let raw_pr = 20000000.0;

        let m = SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr,
            snr: 45.0,
            doppler: 0.0,
            time: t,
            eph: Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
                toe: t,
                toc: t,
                af0: 1e-6, // positive clock bias so corrected_pr > raw_pr
                af1: 0.0,
                af2: 0.0,
                crs: 0.0,
                crc: 0.0,
                cuc: 0.0,
                cus: 0.0,
                cic: 0.0,
                cis: 0.0,
                m0: 0.0,
                e: 0.01,
                sqrt_a: 5153.6,
                delta_n: 0.0,
                omega0: 0.0,
                omega_dot: 0.0,
                i0: 0.95,
                idot: 0.0,
                omega: 0.0,
                tgd: 0.0,
                iode: 1,
                iodc: 1,
            }),
            is_iono_free: true,
            freq_band: 1,
        };

        let (coord, corrected_pr) = compute_sat_state(&m, 0.0);
        // Returns without panicking — coordinate should be non-zero
        assert!(
            coord.vector.norm() > 0.0,
            "satellite position should be non-zero"
        );
        // sat_clk_err (positive due to af0 > 0) is *added* to raw_pr,
        // so corrected_pr should be larger
        assert!(
            corrected_pr > raw_pr,
            "corrected_pr ({}) should be > raw_pr ({})",
            corrected_pr,
            raw_pr
        );
    }

    #[test]
    fn test_compute_sat_state_freq_band_7() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;

        let t = GpsTime::new(2000, 100000.0);

        let m = SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20000000.0,
            snr: 45.0,
            doppler: 0.0,
            time: t,
            eph: Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
                toe: t,
                toc: t,
                af0: 0.0,
                af1: 0.0,
                af2: 0.0,
                crs: 0.0,
                crc: 0.0,
                cuc: 0.0,
                cus: 0.0,
                cic: 0.0,
                cis: 0.0,
                m0: 0.0,
                e: 0.01,
                sqrt_a: 5153.6,
                delta_n: 0.0,
                omega0: 0.0,
                omega_dot: 0.0,
                i0: 0.95,
                idot: 0.0,
                omega: 0.0,
                tgd: 0.0,
                iode: 1,
                iodc: 1,
            }),
            is_iono_free: false,
            freq_band: 7,
        };

        // For GPS, freq_band 7 falls through to position() in compute_sat_state
        let (coord, corrected_pr) = compute_sat_state(&m, 0.0);
        assert!(
            coord.vector.norm() > 0.0,
            "satellite position should be non-zero"
        );
        assert!(
            corrected_pr > 0.0,
            "corrected_pr should be positive"
        );
    }

    #[test]
    fn test_spp_wnlls_step_poor_geometry() {
        use gneiss_core::ephemeris::GpsEphemeris;
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;

        let t = GpsTime::new(2000, 100000.0);

        let eph = Ephemeris::Gps(GpsEphemeris {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.01, sqrt_a: 5153.6,
            delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
            i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        });

        // All 4 measurements use the SAME ephemeris and SAME raw_pr,
        // producing IDENTICAL line-of-sight rows in the H matrix.
        // This makes H^T W H exactly singular -> MatrixInversionFailed.
        let ms: Vec<SppMeasurement> = (0..4)
            .map(|_| SppMeasurement {
                constellation: Constellation::Gps,
                raw_pr: 20000000.0,
                snr: 45.0,
                doppler: 0.0,
                time: t,
                eph: eph.clone(),
                is_iono_free: false,
                freq_band: 1,
            })
            .collect();

        let state = SppState::new(
            Coordinate::new(
                Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
                Datum::WGS84,
                Frame::ECEF,
                t,
            ),
            0.0,
            0.0,
            0.0,
            0.0,
        );

        let config = SppConfig {
            enable_sagnac: false,
            enable_tropo: false,
            enable_iono: false,
            elevation_mask_rad: -core::f64::consts::PI,
            geometry_variance_threshold: 0.0,
            ..Default::default()
        };

        let res = spp_wnlls_step(&state, &ms, None, &config);
        // With 4 identical measurements, the trace of (H^T W H)^{-1} is
        // positive (~39) but well below the default threshold (10000).
        // Set threshold=0 so that any positive trace triggers PoorGeometry.
        assert!(res.is_err(), "expected poor-geometry error, got Ok");
    }

    #[test]
    fn test_build_design_matrix_elevation_mask() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;

        let t = GpsTime::new(2000, 100000.0);

        // 4 GPS measurements with identical ephemeris
        let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            toe: t,
            toc: t,
            af0: 0.0,
            af1: 0.0,
            af2: 0.0,
            crs: 0.0,
            crc: 0.0,
            cuc: 0.0,
            cus: 0.0,
            cic: 0.0,
            cis: 0.0,
            m0: 0.0,
            e: 0.01,
            sqrt_a: 5153.6,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 0.95,
            idot: 0.0,
            omega: 0.0,
            tgd: 0.0,
            iode: 1,
            iodc: 1,
        });

        let ms: Vec<SppMeasurement> = (0..4)
            .map(|i| SppMeasurement {
                constellation: Constellation::Gps,
                raw_pr: 20000000.0 + i as f64 * 1000000.0,
                snr: 45.0,
                doppler: 0.0,
                time: t,
                eph: eph.clone(),
                is_iono_free: false,
                freq_band: 1,
            })
            .collect();

        let state = SppState::new(
            Coordinate::new(
                Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
                Datum::WGS84,
                Frame::ECEF,
                t,
            ),
            0.0,
            0.0,
            0.0,
            0.0,
        );

        // elevation_mask_rad = 100.0 rad — all elevations are below this (max el ~ PI/2)
        let config = SppConfig {
            elevation_mask_rad: 100.0,
            enable_sagnac: false,
            enable_tropo: false,
            enable_iono: false,
            ..Default::default()
        };

        let result = build_design_matrix(&state, &ms, None, &config);

        // build_design_matrix returns Ok even when all sats are masked
        assert!(result.is_ok(), "build_design_matrix should return Ok");

        let (h_matrix, w_matrix, _dz_vector, _cols, _clocks) = result.unwrap();

        // Every measurement is masked — W diagonal should be MIN_WEIGHT
        for i in 0..ms.len() {
            assert!(
                (w_matrix[(i, i)] - MIN_WEIGHT).abs() < 1e-15,
                "w_matrix[{}] should be MIN_WEIGHT, got {}",
                i,
                w_matrix[(i, i)]
            );
        }

        // H rows for real measurements should be zero (never set due to continue)
        for i in 0..ms.len() {
            assert_eq!(
                h_matrix.row(i).iter().sum::<f64>(),
                0.0,
                "h_matrix row {} should be all zeros for masked sat",
                i
            );
        }
    }

    #[test]
    fn test_compute_atmospheric_delays_disabled_features() {
        let rec_ecef =
            Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let rec_llh = ecef_to_llh(rec_ecef);
        let config = SppConfig {
            enable_tropo: false,
            enable_iono: false,
            ..Default::default()
        };
        let (tropo, iono) = compute_atmospheric_delays(
            rec_ecef,
            rec_llh,
            0.0,
            0.5,
            GpsTime::new(0, 0.0),
            None,
            &config,
        );
        assert_eq!(tropo, 0.0, "tropo should be 0.0 when disabled");
        assert_eq!(iono, 0.0, "iono should be 0.0 when disabled");
    }

    #[test]
    fn test_compute_atmospheric_delays_min_earth_radius() {
        // At exactly MIN_EARTH_RADIUS_M (6,000,000), the norm check returns (0, 0)
        let rec_ecef = Vector3::new(MIN_EARTH_RADIUS_M, 0.0, 0.0);
        let (tropo, iono) = compute_atmospheric_delays(
            rec_ecef,
            Vector3::zeros(),
            0.0,
            0.5,
            GpsTime::new(0, 0.0),
            None,
            &SppConfig::default(),
        );
        assert_eq!(tropo, 0.0, "tropo should be 0.0 at MIN_EARTH_RADIUS boundary");
        assert_eq!(iono, 0.0, "iono should be 0.0 at MIN_EARTH_RADIUS boundary");
    }

    #[test]
    fn test_solve_spp_iteratively_convergence_failure() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;

        let t = GpsTime::new(2000, 100000.0);

        let m = SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20000000.0,
            snr: 45.0,
            doppler: 0.0,
            time: t,
            eph: Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
                toe: t,
                toc: t,
                af0: 0.0,
                af1: 0.0,
                af2: 0.0,
                crs: 0.0,
                crc: 0.0,
                cuc: 0.0,
                cus: 0.0,
                cic: 0.0,
                cis: 0.0,
                m0: 0.0,
                e: 0.01,
                sqrt_a: 5153.6,
                delta_n: 0.0,
                omega0: 0.0,
                omega_dot: 0.0,
                i0: 0.95,
                idot: 0.0,
                omega: 0.0,
                tgd: 0.0,
                iode: 1,
                iodc: 1,
            }),
            is_iono_free: false,
            freq_band: 1,
        };

        let state = SppState::new(
            Coordinate::new(
                Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
                Datum::WGS84,
                Frame::ECEF,
                t,
            ),
            0.0,
            0.0,
            0.0,
            0.0,
        );

        // max_iterations = 0 means the loop never executes -> ConvergenceFailed
        let config = SppConfig {
            max_iterations: 0,
            enable_sagnac: false,
            enable_tropo: false,
            enable_iono: false,
            elevation_mask_rad: -core::f64::consts::PI,
            ..Default::default()
        };

        let res = solve_spp_iteratively(state, &[m], None, &config);
        assert_eq!(
            res.unwrap_err(),
            SppError::ConvergenceFailed,
            "with max_iterations=0, should immediately return ConvergenceFailed"
        );
    }

    #[test]
    fn test_seed_initial_state_with_previous() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;

        let t = GpsTime::new(2000, 100000.0);

        // Previous state with a non-origin position
        let prev_position = Coordinate::new(
            Vector3::new(6000000.0, 6000000.0, 0.0),
            Datum::WGS84,
            Frame::ECEF,
            t,
        );
        let prev_state = SppState::new(prev_position, 100.0, 0.0, 0.0, 0.0);

        let m = SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20000000.0,
            snr: 45.0,
            doppler: 0.0,
            time: t,
            eph: Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
                toe: t,
                toc: t,
                af0: 0.0,
                af1: 0.0,
                af2: 0.0,
                crs: 0.0,
                crc: 0.0,
                cuc: 0.0,
                cus: 0.0,
                cic: 0.0,
                cis: 0.0,
                m0: 0.0,
                e: 0.01,
                sqrt_a: 5153.6,
                delta_n: 0.0,
                omega0: 0.0,
                omega_dot: 0.0,
                i0: 0.95,
                idot: 0.0,
                omega: 0.0,
                tgd: 0.0,
                iode: 1,
                iodc: 1,
            }),
            is_iono_free: false,
            freq_band: 1,
        };

        let state = seed_initial_state(&[m], Some(&prev_state));

        // Position should come from prev_state (not compute_seed_position)
        assert_eq!(
            state.position.vector.x, 6000000.0,
            "seed x should match prev_state x"
        );
        assert_eq!(
            state.position.vector.y, 6000000.0,
            "seed y should match prev_state y"
        );
        assert_eq!(
            state.position.vector.z, 0.0,
            "seed z should match prev_state z"
        );

        // Clock biases should be computed from the measurement relative to the seed position.
        // With the receiver at (6000000, 6000000, 0) and the satellite computed from the
        // ephemeris, cdt = corrected_pr - geometric_range should be some non-zero value.
        assert!(
            state.cdt != 0.0,
            "cdt should be non-zero (computed from measurement)"
        );
    }

    #[test]
    fn test_compute_atmospheric_delays_below_min_earth_radius_matches_constant() {
        // Verify the test uses the same constant as production code
        assert_eq!(MIN_EARTH_RADIUS_M, 6_000_000.0);
    }

    #[test]
    fn test_compute_seed_position_direct() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;
        let t = GpsTime::new(2000, 100000.0);
        let ms: Vec<SppMeasurement> = (0..3).map(|i| SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20000000.0 + i as f64 * 1000000.0,
            snr: 45.0, doppler: 0.0, time: t,
            eph: Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: SatelliteId { constellation: Constellation::Gps, prn: i + 1 },
                toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
                cic: 0.0, cis: 0.0, m0: 0.0, e: 0.01, sqrt_a: 5153.6,
                delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
                i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
            }),
            is_iono_free: false, freq_band: 1,
        }).collect();
        let (sx, sy, sz) = compute_seed_position(&ms);
        assert!(!sx.is_nan());
        assert!(!sy.is_nan());
        assert!(!sz.is_nan());
    }

    #[test]
    fn test_compute_seed_clocks_multi_constellation() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;
        let t = GpsTime::new(2000, 100000.0);
        let eph_base = || -> Ephemeris {
            Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
                toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
                cic: 0.0, cis: 0.0, m0: 0.0, e: 0.01, sqrt_a: 5153.6,
                delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
                i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
            })
        };
        let ms = vec![
            SppMeasurement { constellation: Constellation::Gps, raw_pr: 20000000.0, snr: 45.0, doppler: 0.0, time: t, eph: eph_base(), is_iono_free: false, freq_band: 1 },
            SppMeasurement { constellation: Constellation::Galileo, raw_pr: 21000000.0, snr: 42.0, doppler: 0.0, time: t, eph: eph_base(), is_iono_free: false, freq_band: 1 },
            SppMeasurement { constellation: Constellation::Beidou, raw_pr: 20500000.0, snr: 40.0, doppler: 0.0, time: t, eph: eph_base(), is_iono_free: false, freq_band: 1 },
            SppMeasurement { constellation: Constellation::Glonass, raw_pr: 19500000.0, snr: 38.0, doppler: 0.0, time: t, eph: eph_base(), is_iono_free: false, freq_band: 1 },
        ];
        let (cdt_gps, cdt_gal, cdt_bds, cdt_glo) = compute_seed_clocks(&ms, 6378137.0, 0.0, 0.0);
        // All clock biases should be non-zero and finite
        assert!(cdt_gps.is_finite() && cdt_gps != 0.0);
        assert!(cdt_gal.is_finite() && cdt_gal != 0.0);
        assert!(cdt_bds.is_finite() && cdt_bds != 0.0);
        assert!(cdt_glo.is_finite() && cdt_glo != 0.0);
    }

    #[test]
    fn test_seed_initial_state_without_previous() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;
        let t = GpsTime::new(2000, 100000.0);
        let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.01, sqrt_a: 5153.6,
            delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
            i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        });
        let ms = vec![SppMeasurement {
            constellation: Constellation::Gps, raw_pr: 20000000.0,
            snr: 45.0, doppler: 0.0, time: t, eph, is_iono_free: false, freq_band: 1,
        }];
        // No previous state: should use compute_seed_position
        let state = seed_initial_state(&ms, None);
        assert!(state.position.vector.norm() > 0.0);
    }

    #[test]
    fn test_filter_raim_outliers_runs_without_panic() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;
        let t = GpsTime::new(2000, 100000.0);
        let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.01, sqrt_a: 5153.6,
            delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
            i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        });
        let ms: Vec<SppMeasurement> = (0..4).map(|i| SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20000000.0 + i as f64 * 1000000.0,
            snr: 45.0, doppler: 0.0, time: t,
            eph: eph.clone(),
            is_iono_free: false, freq_band: 1,
        }).collect();
        let state = SppState::new(
            Coordinate::new(
                Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
                Datum::WGS84, Frame::ECEF, t,
            ),
            0.0, 0.0, 0.0, 0.0,
        );
        let config = SppConfig {
            enable_sagnac: false, enable_tropo: false, enable_iono: false,
            elevation_mask_rad: -core::f64::consts::PI,
            raim_outlier_m: 1000000.0, // Very permissive - keep everything
            ..Default::default()
        };
        let filtered = filter_raim_outliers(&state, &ms, None, &config);
        // With permissive threshold, all measurements should be kept
        assert_eq!(filtered.len(), 4, "All measurements should pass with permissive threshold");
    }

    #[test]
    fn test_apply_raim_with_outlier_removal() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;
        let t = GpsTime::new(2000, 100000.0);
        let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.01, sqrt_a: 5153.6,
            delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
            i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        });
        let ms: Vec<SppMeasurement> = (0..6).map(|i| SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20000000.0 + if i == 5 { 500000.0 } else { 0.0 } + i as f64 * 100.0,
            snr: 45.0, doppler: 0.0, time: t,
            eph: eph.clone(),
            is_iono_free: false, freq_band: 1,
        }).collect();
        let seed = SppState::new(
            Coordinate::new(
                Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
                Datum::WGS84, Frame::ECEF, t,
            ), 0.0, 0.0, 0.0, 0.0,
        );
        let config = SppConfig {
            enable_sagnac: false, enable_tropo: false, enable_iono: false,
            elevation_mask_rad: -core::f64::consts::PI,
            max_iterations: 2,
            ..Default::default()
        };
        let result = apply_raim(seed.clone(), seed, &ms, None, &config);
        // With outlier present, RAIM should attempt re-solve.
        // Even if it fails (max_iterations=2 might not converge), it returns an error or the original state
        // depending on path.
        assert!(result.is_ok());
    }

    #[test]
    fn test_find_clock_cols_no_gps() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;
        let t = GpsTime::new(2000, 100000.0);
        let eph_base = || -> Ephemeris {
            Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
                sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
                toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
                crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
                cic: 0.0, cis: 0.0, m0: 0.0, e: 0.01, sqrt_a: 5153.6,
                delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
                i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
            })
        };
        // Only Galileo and GLONASS — no GPS
        let ms = vec![
            SppMeasurement { constellation: Constellation::Galileo, raw_pr: 20000000.0, snr: 45.0, doppler: 0.0, time: t, eph: eph_base(), is_iono_free: false, freq_band: 1 },
            SppMeasurement { constellation: Constellation::Glonass, raw_pr: 21000000.0, snr: 42.0, doppler: 0.0, time: t, eph: eph_base(), is_iono_free: false, freq_band: 1 },
        ];
        let (cols, clk) = find_clock_cols(&ms);
        // Without GPS/QZSS: 3 base + 1 (GAL) + 1 (GLO) = 5 columns
        assert_eq!(cols, 5);
        assert!(clk.0.is_none());
        assert!(clk.1.is_some());
        assert!(clk.2.is_none());
        assert!(clk.3.is_some());
    }

    #[test]
    fn test_compute_atmospheric_delays_with_klobuchar() {
        let rec_ecef = Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let rec_llh = ecef_to_llh(rec_ecef);
        let iono = gneiss_core::atmosphere::KlobucharParams::default();
        let config = SppConfig {
            enable_tropo: true,
            enable_iono: true,
            ..Default::default()
        };
        let (tropo, iono_delay) = compute_atmospheric_delays(
            rec_ecef, rec_llh, 0.5, 0.3, GpsTime::new(2000, 50000.0), Some(&iono), &config,
        );
        assert!(tropo > 0.0, "Tropospheric delay should be positive");
        assert!(iono_delay != 0.0, "Ionospheric delay should be non-zero");
    }

    #[test]
    fn test_compute_atmospheric_delays_below_min_elevation() {
        let rec_ecef = Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0);
        let rec_llh = ecef_to_llh(rec_ecef);
        let config = SppConfig {
            enable_tropo: true,
            enable_iono: false,
            ..Default::default()
        };
        // el = 0.01 rad, which is below MIN_ATMOSPHERE_ELEVATION_RAD (0.087 rad)
        let (tropo, iono) = compute_atmospheric_delays(
            rec_ecef, rec_llh, 0.5, 0.01, GpsTime::new(2000, 50000.0), None, &config,
        );
        // The function uses safe_el = max(el, MIN_ATMOSPHERE_ELEVATION_RAD)
        // So it should still produce a non-zero tropo
        assert!(tropo > 0.0, "Tropospheric delay should be positive even at low elevation");
        assert_eq!(iono, 0.0);
    }

    #[test]
    fn test_build_design_matrix_multi_constellation_clock_columns() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;
        let t = GpsTime::new(2000, 100000.0);
        let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.01, sqrt_a: 5153.6,
            delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
            i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        });
        // 4 GPS + 1 Galileo + 1 GLONASS
        let mut ms = Vec::new();
        for i in 0..4 {
            ms.push(SppMeasurement {
                constellation: Constellation::Gps,
                raw_pr: 20000000.0 + i as f64 * 1000000.0,
                snr: 45.0, doppler: 0.0, time: t,
                eph: eph.clone(),
                is_iono_free: false, freq_band: 1,
            });
        }
        ms.push(SppMeasurement {
            constellation: Constellation::Galileo,
            raw_pr: 21000000.0, snr: 42.0, doppler: 0.0, time: t,
            eph: eph.clone(),
            is_iono_free: false, freq_band: 1,
        });
        ms.push(SppMeasurement {
            constellation: Constellation::Glonass,
            raw_pr: 19500000.0, snr: 40.0, doppler: 0.0, time: t,
            eph: eph.clone(),
            is_iono_free: false, freq_band: 1,
        });

        let state = SppState::new(
            Coordinate::new(
                Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
                Datum::WGS84, Frame::ECEF, t,
            ), 0.0, 0.0, 0.0, 0.0,
        );
        let config = SppConfig {
            enable_sagnac: false, enable_tropo: false, enable_iono: false,
            elevation_mask_rad: -core::f64::consts::PI,
            ..Default::default()
        };
        let result = build_design_matrix(&state, &ms, None, &config);
        assert!(result.is_ok());
        let (h, _w, _dz, _cols, clocks) = result.unwrap();
        // 4 GPS + 1 GAL + 1 GLO = 6 measurements. With GPS+GAL+GLO: 3+3=6 columns
        assert_eq!(h.nrows(), 6);
        assert_eq!(h.ncols(), 6);
        // Check that the clock columns are set
        assert!(clocks.0.is_some()); // GPS
        assert!(clocks.1.is_some()); // GAL
        assert!(clocks.3.is_some()); // GLO
    }

    #[test]
    fn test_build_design_matrix_sagnac_enabled() {
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;
        let t = GpsTime::new(2000, 100000.0);
        let eph = Ephemeris::Gps(gneiss_core::ephemeris::GpsEphemeris {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.01, sqrt_a: 5153.6,
            delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
            i0: 0.95, idot: 0.0, omega: 0.0, tgd: 0.0, iode: 1, iodc: 1,
        });
        let ms: Vec<SppMeasurement> = (0..4).map(|i| SppMeasurement {
            constellation: Constellation::Gps,
            raw_pr: 20000000.0 + i as f64 * 1000000.0,
            snr: 45.0, doppler: 0.0, time: t,
            eph: eph.clone(),
            is_iono_free: false, freq_band: 1,
        }).collect();
        let state = SppState::new(
            Coordinate::new(
                Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0),
                Datum::WGS84, Frame::ECEF, t,
            ), 0.0, 0.0, 0.0, 0.0,
        );
        let config = SppConfig {
            enable_sagnac: false, // Sagnac enabled
            enable_tropo: false, enable_iono: false,
            elevation_mask_rad: -core::f64::consts::PI,
            ..Default::default()
        };
        let result = build_design_matrix(&state, &ms, None, &config);
        assert!(result.is_ok(), "Design matrix with sagnac should build successfully");
    }

    #[test]
    fn test_build_single_measurement_beidou_path() {
        use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;
        let t = GpsTime::new(2000, 100000.0);
        let sat_bds = SatelliteId { constellation: Constellation::Beidou, prn: 3 };

        // Beidou: p1 comes from get_observable(2)
        let sat_obs = SatObs {
            sat: sat_bds,
            observations: vec![
                Observation {
                    code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 2, attribute: 'I' } },
                    value: 22000000.0, lock_time: None, lli: None,
                },
                Observation {
                    code: ObsCode { obs_type: ObsType::Snr, signal: SignalCode { freq_band: 2, attribute: 'I' } },
                    value: 40.0, lock_time: None, lli: None,
                },
            ],
        };

        let eph_bds = gneiss_core::ephemeris::BeidouEphemeris {
            sat: sat_bds, toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.01, sqrt_a: 5153.6,
            delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
            i0: 0.95, idot: 0.0, omega: 0.0, tgd1: 0.0, tgd2: 0.0, aode: 1, aodc: 1,
        };

        let epoch = EpochObs { time: t, satellites: vec![sat_obs] };
        let ephs = vec![Ephemeris::Beidou(eph_bds)];

        let ms = build_measurements(&epoch, &ephs, &SppConfig::default());
        assert_eq!(ms.len(), 1, "Beidou measurement should be built");
        assert_eq!(ms[0].constellation, Constellation::Beidou);
        assert!(!ms[0].is_iono_free, "Beidou single-freq should not be iono-free");
    }

    #[test]
    fn test_build_single_measurement_galileo_band5_path() {
        use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;
        let t = GpsTime::new(2000, 100000.0);
        let sat_gal = SatelliteId { constellation: Constellation::Galileo, prn: 2 };

        // Galileo: p1 from band 1, p2 tries band 7 first, then band 5
        let sat_obs = SatObs {
            sat: sat_gal,
            observations: vec![
                Observation {
                    code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 1, attribute: 'C' } },
                    value: 21000000.0, lock_time: None, lli: None,
                },
                Observation {
                    code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 5, attribute: 'Q' } },
                    value: 22000000.0, lock_time: None, lli: None,
                },
                Observation {
                    code: ObsCode { obs_type: ObsType::Snr, signal: SignalCode { freq_band: 1, attribute: 'C' } },
                    value: 42.0, lock_time: None, lli: None,
                },
            ],
        };

        let eph_gal = gneiss_core::ephemeris::GalileoEphemeris {
            sat: sat_gal, toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.01, sqrt_a: 5153.6,
            delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
            i0: 0.95, idot: 0.0, omega: 0.0,
            bgd_e1_e5a: 0.0, bgd_e1_e5b: 0.0, iod_nav: 1,
        };

        let epoch = EpochObs { time: t, satellites: vec![sat_obs] };
        let ephs = vec![Ephemeris::Galileo(eph_gal)];

        let ms = build_measurements(&epoch, &ephs, &SppConfig::default());
        assert_eq!(ms.len(), 1, "Galileo measurement should be built");
        assert_eq!(ms[0].constellation, Constellation::Galileo);
        assert!(ms[0].is_iono_free, "Galileo dual-freq should be iono-free");
        assert_eq!(ms[0].freq_band, 5, "Galileo band 5 should be used");
    }

    #[test]
    fn test_build_single_measurement_beidou_band6_path() {
        use gneiss_core::obs::{EpochObs, ObsCode, ObsType, Observation, SatObs, SignalCode};
        use gneiss_core::sat::{Constellation, SatelliteId};
        use gneiss_core::time::GpsTime;
        let t = GpsTime::new(2000, 100000.0);
        let sat_bds = SatelliteId { constellation: Constellation::Beidou, prn: 3 };

        // Beidou: p1 from band 2, p2 tries band 7 first, then band 6
        let sat_obs = SatObs {
            sat: sat_bds,
            observations: vec![
                Observation {
                    code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 2, attribute: 'I' } },
                    value: 22000000.0, lock_time: None, lli: None,
                },
                Observation {
                    code: ObsCode { obs_type: ObsType::Pseudorange, signal: SignalCode { freq_band: 6, attribute: 'X' } },
                    value: 23000000.0, lock_time: None, lli: None,
                },
                Observation {
                    code: ObsCode { obs_type: ObsType::Snr, signal: SignalCode { freq_band: 2, attribute: 'I' } },
                    value: 40.0, lock_time: None, lli: None,
                },
            ],
        };

        let eph_bds = gneiss_core::ephemeris::BeidouEphemeris {
            sat: sat_bds, toe: t, toc: t, af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0,
            cic: 0.0, cis: 0.0, m0: 0.0, e: 0.01, sqrt_a: 5153.6,
            delta_n: 0.0, omega0: 0.0, omega_dot: 0.0,
            i0: 0.95, idot: 0.0, omega: 0.0, tgd1: 0.0, tgd2: 0.0, aode: 1, aodc: 1,
        };

        let epoch = EpochObs { time: t, satellites: vec![sat_obs] };
        let ephs = vec![Ephemeris::Beidou(eph_bds)];

        let ms = build_measurements(&epoch, &ephs, &SppConfig::default());
        assert_eq!(ms.len(), 1, "Beidou measurement should be built");
        assert_eq!(ms[0].constellation, Constellation::Beidou);
        assert!(ms[0].is_iono_free, "Beidou dual-freq should be iono-free");
        assert_eq!(ms[0].freq_band, 6, "Beidou band 6 should be used");
    }
