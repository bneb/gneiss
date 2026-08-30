#![allow(clippy::unwrap_used)]

use super::*;
use crate::sat::{Constellation, SatelliteId};
use crate::time::GpsTime;
use nalgebra::Vector3;

    #[test]
    fn test_gps_position_calculation() {
        let eph = GpsEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100010.0),
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
            omega_dot: -2.0e-9,
            i0: 0.95,
            idot: 0.0,
            omega: 0.0,
            tgd: 0.0,
            iode: 1,
            iodc: 1,
        };

        let (pos, vel, _clk_err, _clk_drift) = eph.position(GpsTime::new(2000, 100000.0));
        assert!(pos.norm() > 20_000_000.0);
        assert!(pos.norm() < 30_000_000.0);
        assert!(vel.norm() > 1000.0); // GPS satellites move at ~3.9 km/s
    }

    #[test]
    fn test_glonass_rk4_physics() {
        // Test GLONASS Cartesian RK4 numerical integrator
        let eph = GlonassEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Glonass,
                prn: 1,
            },
            toe: GpsTime::new(2000, 100000.0),
            freq_num: 1,
            tau_n: 1e-5,
            gamma_n: 1e-9,
            delta_tau_n: 0.0,
            x: 10_000_000.0,
            y: 15_000_000.0,
            z: 20_000_000.0,
            vx: -2000.0,
            vy: 1500.0,
            vz: 1000.0,
            ax: 0.0,
            ay: 0.0,
            az: 0.0, // Solar/lunar accels
        };

        let (pos, vel, clk_err, _clk_drift) = eph.position(GpsTime::new(2000, 100060.0)); // Propagate 60 seconds

        // 60s at roughly 2.5km/s gives about 150km change
        let dist_moved = (pos - Vector3::new(10_000_000.0, 15_000_000.0, 20_000_000.0)).norm();
        assert!(dist_moved > 100_000.0 && dist_moved < 200_000.0);
        assert!(vel.norm() > 1000.0);
        assert!((clk_err - (1e-5 + 1e-9 * 60.0)).abs() < 1e-12);
    }

    #[test]
    fn test_ephemeris_enum_dispatch() {
        let gal_eph = GalileoEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Galileo,
                prn: 2,
            },
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100010.0),
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
            sqrt_a: 5440.6,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: -2.0e-9,
            i0: 0.95,
            idot: 0.0,
            omega: 0.0,
            bgd_e1_e5a: 0.0,
            bgd_e1_e5b: 0.0,
            iod_nav: 1,
        };

        let enum_eph = Ephemeris::Galileo(gal_eph.clone());
        assert_eq!(enum_eph.sat().constellation, Constellation::Galileo);

        let (pos1, _, _, _) = gal_eph.position(GpsTime::new(2000, 100000.0));
        let (pos2, _, _, _) = enum_eph.position(GpsTime::new(2000, 100000.0));
        assert_eq!(pos1, pos2);
    }

    #[test]
    fn test_ephemeris_exact() {
        let gps_eph = GpsEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100010.0),
            af0: 1.0,
            af1: 2.0,
            af2: 3.0,
            crs: 4.0,
            crc: 5.0,
            cuc: 6.0,
            cus: 7.0,
            cic: 8.0,
            cis: 9.0,
            m0: 0.1,
            e: 0.01,
            sqrt_a: 5153.6,
            delta_n: 0.001,
            omega0: 0.2,
            omega_dot: -2.0e-9,
            i0: 0.95,
            idot: 0.01,
            omega: 0.3,
            tgd: 0.02,
            iode: 1,
            iodc: 1,
        };
        let gal_eph = GalileoEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Galileo,
                prn: 2,
            },
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100010.0),
            af0: 1.0,
            af1: 2.0,
            af2: 3.0,
            crs: 4.0,
            crc: 5.0,
            cuc: 6.0,
            cus: 7.0,
            cic: 8.0,
            cis: 9.0,
            m0: 0.1,
            e: 0.01,
            sqrt_a: 5153.6,
            delta_n: 0.001,
            omega0: 0.2,
            omega_dot: -2.0e-9,
            i0: 0.95,
            idot: 0.01,
            omega: 0.3,
            bgd_e1_e5a: 0.02,
            bgd_e1_e5b: 0.01,
            iod_nav: 1,
        };
        let bds_geo_eph = BeidouEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Beidou,
                prn: 1,
            }, // Geo
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100010.0),
            af0: 1.0,
            af1: 2.0,
            af2: 3.0,
            crs: 4.0,
            crc: 5.0,
            cuc: 6.0,
            cus: 7.0,
            cic: 8.0,
            cis: 9.0,
            m0: 0.1,
            e: 0.01,
            sqrt_a: 5153.6,
            delta_n: 0.001,
            omega0: 0.2,
            omega_dot: -2.0e-9,
            i0: 0.95,
            idot: 0.01,
            omega: 0.3,
            tgd1: 0.02,
            tgd2: 0.0,
            aode: 1,
            aodc: 1,
        };
        let bds_igso_eph = BeidouEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Beidou,
                prn: 6,
            }, // Not Geo
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100010.0),
            af0: 1.0,
            af1: 2.0,
            af2: 3.0,
            crs: 4.0,
            crc: 5.0,
            cuc: 6.0,
            cus: 7.0,
            cic: 8.0,
            cis: 9.0,
            m0: 0.1,
            e: 0.01,
            sqrt_a: 5153.6,
            delta_n: 0.001,
            omega0: 0.2,
            omega_dot: -2.0e-9,
            i0: 0.95,
            idot: 0.01,
            omega: 0.3,
            tgd1: 0.02,
            tgd2: 0.0,
            aode: 1,
            aodc: 1,
        };
        let qzss_eph = QzssEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Qzss,
                prn: 4,
            },
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100010.0),
            af0: 1.0,
            af1: 2.0,
            af2: 3.0,
            crs: 4.0,
            crc: 5.0,
            cuc: 6.0,
            cus: 7.0,
            cic: 8.0,
            cis: 9.0,
            m0: 0.1,
            e: 0.01,
            sqrt_a: 5153.6,
            delta_n: 0.001,
            omega0: 0.2,
            omega_dot: -2.0e-9,
            i0: 0.95,
            idot: 0.01,
            omega: 0.3,
            tgd: 0.02,
            iode: 1,
            iodc: 1,
        };
        let glo_eph_fwd = GlonassEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Glonass,
                prn: 5,
            },
            toe: GpsTime::new(2000, 100000.0),
            freq_num: 7,
            tau_n: 1e-5,
            gamma_n: 1e-9,
            delta_tau_n: 0.0,
            x: 10_000_000.0,
            y: 15_000_000.0,
            z: 20_000_000.0,
            vx: -2000.0,
            vy: 1500.0,
            vz: 1000.0,
            ax: 0.1,
            ay: 0.2,
            az: 0.3,
        };
        let glo_eph_bwd = glo_eph_fwd.clone();

        let t = GpsTime::new(2000, 100060.0);
        let e_gps = Ephemeris::Gps(gps_eph.clone());
        let e_gal = Ephemeris::Galileo(gal_eph.clone());
        let e_bds = Ephemeris::Beidou(bds_igso_eph.clone());
        let e_qzss = Ephemeris::Qzss(qzss_eph.clone());
        let e_glo = Ephemeris::Glonass(glo_eph_fwd.clone());

        assert_eq!(
            e_gps.sat(),
            SatelliteId {
                constellation: Constellation::Gps,
                prn: 1
            }
        );
        assert_eq!(e_gps.toe(), GpsTime::new(2000, 100000.0));
        assert_eq!(e_gps.freq_num(), 0);

        assert_eq!(
            e_gal.sat(),
            SatelliteId {
                constellation: Constellation::Galileo,
                prn: 2
            }
        );
        assert_eq!(e_gal.toe(), GpsTime::new(2000, 100000.0));
        assert_eq!(e_gal.freq_num(), 0);

        assert_eq!(
            e_bds.sat(),
            SatelliteId {
                constellation: Constellation::Beidou,
                prn: 6
            }
        );
        assert_eq!(e_bds.toe(), GpsTime::new(2000, 100000.0));
        assert_eq!(e_bds.freq_num(), 0);

        assert_eq!(
            e_qzss.sat(),
            SatelliteId {
                constellation: Constellation::Qzss,
                prn: 4
            }
        );
        assert_eq!(e_qzss.toe(), GpsTime::new(2000, 100000.0));
        assert_eq!(e_qzss.freq_num(), 0);

        assert_eq!(
            e_glo.sat(),
            SatelliteId {
                constellation: Constellation::Glonass,
                prn: 5
            }
        );
        assert_eq!(e_glo.toe(), GpsTime::new(2000, 100000.0));
        assert_eq!(e_glo.freq_num(), 7);

        // Calculate positions
        let _p_gps = gps_eph.position(t);
        let _p_gal = gal_eph.position(t);
        let _p_bds_geo = bds_geo_eph.position(t);
        let _p_bds_igso = bds_igso_eph.position(t);
        let _p_qzss = qzss_eph.position(t);
        let _p_glo_fwd = glo_eph_fwd.position(t);
        let _p_glo_bwd = glo_eph_bwd.position(GpsTime::new(2000, 99940.0));

        // Assert exact values to kill arithmetic mutants
        use nalgebra::Vector3;

        let p_gps = gps_eph.position(t);
        let p_gal = gal_eph.position(t);
        let p_bds_geo = bds_geo_eph.position(t);
        let p_bds_igso = bds_igso_eph.position(t);
        let p_qzss = qzss_eph.position(t);
        let p_glo_fwd = glo_eph_fwd.position(t);
        let p_glo_bwd = glo_eph_bwd.position(GpsTime::new(2000, 99940.0));

        // macro to assert approx equal for positions to avoid precision issue on diff architectures, but strict enough to catch mutation
        macro_rules! assert_vec_eq {
            ($a:expr, $b:expr) => {
                assert!(
                    ($a - $b).norm() < 1e-6,
                    "Vectors not equal: {:?} and {:?}",
                    $a,
                    $b
                );
            };
        }

        assert_vec_eq!(
            p_gps.0,
            Vector3::new(-20111686.91105315, 16173621.585368361, -5050844.794698619)
        );
        assert_vec_eq!(
            p_gps.1,
            Vector3::new(-42915.9866782879, -49950.87879784308, -59327.05664687349)
        );
        extern crate std;
        std::println!("p_gps.2 = {:.15}", p_gps.2);
        assert!((p_gps.2 - 7600.979999996116).abs() < 1e-12);
        assert!((p_gps.3 - 302.0).abs() < 1e-12);

        assert_vec_eq!(
            p_gal.0,
            Vector3::new(-20111686.90739729, 16173621.585960612, -5050844.807207882)
        );
        assert_vec_eq!(
            p_gal.1,
            Vector3::new(-42915.98666471951, -49950.87870889442, -59327.056414597435)
        );
        assert!((p_gal.2 - 7600.979999996116).abs() < 1e-12);
        assert!((p_gal.3 - 302.0).abs() < 1e-12);

        assert_vec_eq!(
            p_bds_geo.0,
            Vector3::new(-20530859.400528494, 16091202.889367301, -3331282.3000197206)
        );
        assert_vec_eq!(
            p_bds_geo.1,
            Vector3::new(-42143.353117060804, -43745.76008367346, -67353.04008649733)
        );
        assert!((p_bds_geo.2 - 3960.9799999964825).abs() < 1e-12);
        assert!((p_bds_geo.3 - 218.0).abs() < 1e-12);

        assert_vec_eq!(
            p_bds_igso.0,
            Vector3::new(-20532037.822592802, 15739895.005709056, -4715036.370416497)
        );
        assert_vec_eq!(
            p_bds_igso.1,
            Vector3::new(-42188.08236021536, -49443.13151989598, -63140.92184141916)
        );
        assert!((p_bds_igso.2 - 3960.9799999964825).abs() < 1e-12);
        assert!((p_bds_igso.3 - 218.0).abs() < 1e-12);

        assert_vec_eq!(
            p_qzss.0,
            Vector3::new(-20111686.91105315, 16173621.585368361, -5050844.794698619)
        );
        assert_vec_eq!(
            p_qzss.1,
            Vector3::new(-42915.9866782879, -49950.87879784308, -59327.05664687349)
        );
        assert!((p_qzss.2 - 7600.979999996116).abs() < 1e-12);
        assert!((p_qzss.3 - 302.0).abs() < 1e-12);

        assert_vec_eq!(
            p_glo_fwd.0,
            Vector3::new(9880305.169245299, 15090476.713825395, 20059805.549791936)
        );
        assert_vec_eq!(
            p_glo_fwd.1,
            Vector3::new(-1989.7750855530473, 1515.879262912734, 993.5292177774875)
        );
        assert!((p_glo_fwd.2 - 1.006e-5).abs() < 1e-12);
        assert!((p_glo_fwd.3 - 1e-9).abs() < 1e-12);

        assert_vec_eq!(
            p_glo_bwd.0,
            Vector3::new(10120298.841760032, 14910478.053419847, 19939804.282132257)
        );
        assert_vec_eq!(
            p_glo_bwd.1,
            Vector3::new(-2009.9085434547917, 1484.0537579365496, 1006.5341630852408)
        );
        assert!((p_glo_bwd.2 - 9.94e-6).abs() < 1e-12);
        assert!((p_glo_bwd.3 - 1e-9).abs() < 1e-12);
    }

    #[test]
    fn test_bds_boundaries() {
        let mut eph = BeidouEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Beidou,
                prn: 5,
            },
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100010.0),
            af0: 1.0,
            af1: 2.0,
            af2: 3.0,
            crs: 4.0,
            crc: 5.0,
            cuc: 6.0,
            cus: 7.0,
            cic: 8.0,
            cis: 9.0,
            m0: 0.1,
            e: 0.01,
            sqrt_a: 5153.6,
            delta_n: 0.001,
            omega0: 0.2,
            omega_dot: -2.0e-9,
            i0: 0.95,
            idot: 0.01,
            omega: 0.3,
            tgd1: 0.02,
            tgd2: 0.0,
            aode: 1,
            aodc: 1,
        };
        let t = GpsTime::new(2000, 100060.0);
        let p5 = eph.position(t); // geo

        eph.sat.prn = 59;
        let p59 = eph.position(t); // geo

        eph.sat.prn = 6;
        let p6 = eph.position(t); // igso

        // Ensure 5 and 59 behave like geo, 6 like igso
        assert_eq!(p5.0, p59.0); // Wait, position will be identical because prn is only used for geo check
        assert_ne!(p5.0, p6.0);
    }

    #[test]
    fn test_glonass_partial_step() {
        let glo_eph = GlonassEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Glonass,
                prn: 5,
            },
            toe: GpsTime::new(2000, 100000.0),
            freq_num: 7,
            tau_n: 1e-5,
            gamma_n: 1e-9,
            delta_tau_n: 0.0,
            x: 10_000_000.0,
            y: 15_000_000.0,
            z: 20_000_000.0,
            vx: -2000.0,
            vy: 1500.0,
            vz: 1000.0,
            ax: 0.1,
            ay: 0.2,
            az: 0.3,
        };
        let t = GpsTime::new(2000, 100050.0);
        let p = glo_eph.position(t);
        // We just assert on exactly what it produces so it locks it in.

        let expected_x = 9900211.557692498;
        let expected_y = 15075331.129011473;
        let expected_z = 20049864.88968714;

        // Use a tiny delta to lock it in
        assert!((p.0.x - expected_x).abs() < 1e-4);
        assert!((p.0.y - expected_y).abs() < 1e-4);
        assert!((p.0.z - expected_z).abs() < 1e-4);
    }

    /// Bug 15 regression test: In dual-frequency iono-free PPP, TGD must NOT be
    /// applied (it cancels in the IF combination).  The Ephemeris::tgd() method
    /// must return the exact value that calc_keplerian subtracted so the PPP code
    /// can add it back.
    #[test]
    fn test_tgd_not_applied_dual_frequency() {
        let tgd_val = 1.5e-8_f64; // ~4.5 m in seconds
        let gps_eph = GpsEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 5,
            },
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100000.0),
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
            i0: 0.0,
            idot: 0.0,
            omega: 0.0,
            tgd: tgd_val,
            iode: 1,
            iodc: 1,
        };
        let t = GpsTime::new(2000, 100000.0);
        let eph = Ephemeris::Gps(gps_eph);
        let (_, _, clk_with_tgd, _) = eph.position(t);
        // calc_keplerian subtracts tgd: clk = af0 + ... - tgd
        // When af0=0 and relativistic correction=0: clk = -tgd_val
        assert!(
            (clk_with_tgd - (-tgd_val)).abs() < 1e-15,
            "calc_keplerian should produce clk = -tgd when af0=0, got {clk_with_tgd}"
        );
        // Dual-frequency PPP path must undo TGD by adding it back:
        let clk_if = clk_with_tgd + eph.tgd();
        // After adding back TGD, the IF-corrected clock should be ≈ 0 (af0=0)
        assert!(
            clk_if.abs() < 1e-15,
            "Dual-freq IF clock (after TGD undo) must be ≈ 0, got {clk_if}"
        );
    }

    /// Bug 16 regression test: Galileo should use BGD(E1,E5b) for E5b observations
    /// (band 7) and BGD(E1,E5a) for E5a observations (band 5).
    #[test]
    fn test_galileo_bgd_e5b() {
        let gal_eph = GalileoEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Galileo,
                prn: 3,
            },
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100000.0),
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
            sqrt_a: 5440.6,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 0.0,
            idot: 0.0,
            omega: 0.0,
            bgd_e1_e5a: 1.0e-9,
            bgd_e1_e5b: 2.0e-9,
            iod_nav: 1,
        };
        let eph = Ephemeris::Galileo(gal_eph);
        // tgd() returns bgd_e1_e5a (used by position())
        assert!(
            (eph.tgd() - 1.0e-9).abs() < 1e-18,
            "tgd() must return bgd_e1_e5a"
        );
        // bgd_e5b() returns bgd_e1_e5b for E5b tracking
        assert!(
            (eph.bgd_e5b() - 2.0e-9).abs() < 1e-18,
            "bgd_e5b() must return bgd_e1_e5b, got {}",
            eph.bgd_e5b()
        );
        // The two should differ
        assert!(
            (eph.tgd() - eph.bgd_e5b()).abs() > 0.5e-9,
            "bgd_e1_e5a and bgd_e1_e5b must differ in this test"
        );
    }

    #[test]
    fn test_broadcast_clock_tgd_correct() {
        // Test GPS
        let gps_eph = GpsEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Gps,
                prn: 1,
            },
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100000.0),
            af0: 1e-3,
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
            i0: 0.0,
            idot: 0.0,
            omega: 0.0,
            tgd: 1.5e-8,
            iode: 1,
            iodc: 1,
        };
        let gps = Ephemeris::Gps(gps_eph.clone());
        let t = GpsTime::new(2000, 100000.0);
        let (_, _, clk_pos, _) = gps.position(t);
        let (_, _, clk_if, _) = gps.position_iono_free(t);
        assert!((clk_if - clk_pos - gps_eph.tgd).abs() < 1e-15);

        // Test Galileo
        let gal_eph = GalileoEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Galileo,
                prn: 1,
            },
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100000.0),
            af0: 1e-3,
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
            sqrt_a: 5440.6,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 0.0,
            idot: 0.0,
            omega: 0.0,
            bgd_e1_e5a: 2.0e-9,
            bgd_e1_e5b: 3.0e-9,
            iod_nav: 1,
        };
        let gal = Ephemeris::Galileo(gal_eph.clone());
        let (_, _, clk_pos, _) = gal.position(t);
        let (_, _, clk_if, _) = gal.position_iono_free(t);
        assert!((clk_if - clk_pos - gal_eph.bgd_e1_e5a).abs() < 1e-15);

        // Test Beidou
        let bds_eph = BeidouEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Beidou,
                prn: 1,
            },
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100000.0),
            af0: 1e-3,
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
            sqrt_a: 5440.6,
            delta_n: 0.0,
            omega0: 0.0,
            omega_dot: 0.0,
            i0: 0.0,
            idot: 0.0,
            omega: 0.0,
            tgd1: 4.0e-9,
            tgd2: 0.0,
            aode: 1,
            aodc: 1,
        };
        let bds = Ephemeris::Beidou(bds_eph.clone());
        let (_, _, clk_pos, _) = bds.position(t);
        let (_, _, clk_if, _) = bds.position_iono_free(t);
        assert!((clk_if - clk_pos - bds_eph.tgd1).abs() < 1e-15);

        // Test Qzss
        let qzss_eph = QzssEphemeris {
            sat: SatelliteId {
                constellation: Constellation::Qzss,
                prn: 193,
            },
            toe: GpsTime::new(2000, 100000.0),
            toc: GpsTime::new(2000, 100000.0),
            af0: 1e-3,
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
            i0: 0.0,
            idot: 0.0,
            omega: 0.0,
            tgd: 1.5e-8,
            iode: 1,
            iodc: 1,
        };
        let qzss = Ephemeris::Qzss(qzss_eph.clone());
        let (_, _, clk_pos, _) = qzss.position(t);
        let (_, _, clk_if, _) = qzss.position_iono_free(t);
        assert!((clk_if - clk_pos - qzss_eph.tgd).abs() < 1e-15);
    }
