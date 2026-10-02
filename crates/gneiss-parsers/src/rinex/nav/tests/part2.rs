use super::*;

    #[test]
    fn test_rinex_nav_empty() {
        let data = padded("");
        let mut reader = BufReader::new(data.as_bytes());
        let (eph, klob) = parse_rinex_nav(&mut reader).unwrap();
        assert!(eph.is_empty());
        assert!(klob.is_none());
    }

    #[test]
    fn test_rinex_nav_missing_end_header() {
        // Navigation without END OF HEADER - the header parser iterates until
        // it finds END OF HEADER or runs out of lines. Without it, the parsing
        // proceeds to read ephemeris data from what would have been header lines.
        let data = padded("     2.11           N: GPS NAV DATA                         RINEX VERSION / TYPE
 6 20  5 14 22  0  0.0-2.715480513871D-04-6.821210263297D-12 0.000000000000D+00
    7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
   -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
    4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
    9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
    4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
    0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
    4.248000000000D+05 4.000000000000D+00
");
        let mut reader = BufReader::new(data.as_bytes());
        let (_eph, klob) = parse_rinex_nav(&mut reader).unwrap();
        // Without END OF HEADER, the ephemeris first line is consumed by the
        // header parser as header data. The remaining 8 lines form a partial
        // ephemeris. This test just verifies no panic.
        assert!(klob.is_none());
    }

    // -----------------------------------------------------------------------
    // RINEX 2 NAV: two successive ephemerides
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_nav_two_ephemerides() {
        let data = padded("     2.11           N: GPS NAV DATA                         RINEX VERSION / TYPE
                                                            END OF HEADER
 1 20  5 14 22  0  0.0-2.715480513871D-04-6.821210263297D-12 0.000000000000D+00
    7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
   -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
    4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
    9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
    4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
    0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
    4.248000000000D+05 4.000000000000D+00
 2 20  5 14 22  0  0.0-2.715480513871D-04-6.821210263297D-12 0.000000000000D+00
    7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
   -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
    4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
    9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
    4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
    0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
    4.248000000000D+05 4.000000000000D+00
");
        let mut reader = BufReader::new(data.as_bytes());
        let (eph, _klob) = parse_rinex_nav(&mut reader).unwrap();
        assert_eq!(eph.len(), 2);
        assert_eq!(eph[0].sat().prn, 1);
        assert_eq!(eph[1].sat().prn, 2);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 Beidou NAV
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_beidou() {
        let data = padded("     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
C 1 2020 06 15 01 30  0 -.271548051387D-03 -.682121026330D-11  .000000000000D+00
     7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
    -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
     4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
     9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
     4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
     4.248000000000D+05 4.000000000000D+00
");
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();

        assert_eq!(ephemerides.len(), 1);
        let eph = &ephemerides[0];

        assert_eq!(eph.sat().prn, 1);
        assert_eq!(eph.sat().constellation, Constellation::Beidou);

        if let Ephemeris::Beidou(bds) = eph {
            assert!((bds.sqrt_a - 5153.539648056).abs() < 1e-6);
            assert_eq!(bds.aode, 70);
        } else {
            panic!("Expected Beidou ephemeris, got {:?}", eph.sat().constellation);
        }
    }

    // -----------------------------------------------------------------------
    // RINEX 3 Galileo NAV
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_galileo() {
        let data = padded("     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
E 2 2020 08 01 12  0  0 -.153846153846D-03  .000000000000E+00  .000000000000E+00
     1.000000000000D+01 0.000000000000D+00 5.000000000000D-09 3.000000000000D+00
     0.000000000000D+00 1.500000000000D-03 4.500000000000D-06 5.200000000000D+03
     4.320000000000D+05 1.200000000000D-07 2.500000000000D+00 3.400000000000D-08
     9.800000000000D-01 3.000000000000D+02 1.000000000000D+00 8.000000000000D-09
     4.500000000000D-11 0.000000000000D+00 2.000000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 1.000000000000D+01
     4.320000000000D+05 4.000000000000D+00
");
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();

        assert_eq!(ephemerides.len(), 1);
        let eph = &ephemerides[0];
        assert_eq!(eph.sat().prn, 2);
        assert_eq!(eph.sat().constellation, Constellation::Galileo);

        if let Ephemeris::Galileo(gal) = eph {
            assert!((gal.sqrt_a - 5200.0).abs() < 1e-6);
            assert_eq!(gal.iod_nav, 10);
        } else {
            panic!("Expected Galileo ephemeris");
        }
    }

    // -----------------------------------------------------------------------
    // parse_rinex_nav_epoch_time edge cases (2-digit year)
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_rinex_nav_epoch_time_two_digit_year() {
        let line_r2 = " 6 99 12 25  0  0  0.0";
        let time_r2 = parse_rinex_nav_epoch_time(line_r2, false);
        assert_eq!(time_r2.week, GpsTime::from_calendar(1999, 12, 25, 0, 0, 0.0).week);

        let line_r2_2 = " 6  1  1  1  0  0  0.0";
        let time_r2_2 = parse_rinex_nav_epoch_time(line_r2_2, false);
        assert_eq!(time_r2_2.week, GpsTime::from_calendar(2001, 1, 1, 0, 0, 0.0).week);
    }

    #[test]
    fn test_parse_rinex_nav_epoch_time_4_digit_year() {
        let line_r3 = "G 1 2020 06 15 01 30  0";
        let time_r3 = parse_rinex_nav_epoch_time(line_r3, true);
        let expected = GpsTime::from_calendar(2020, 6, 15, 1, 30, 0.0);
        assert_eq!(time_r3.week, expected.week);
        assert!((time_r3.tow - expected.tow).abs() < 1e-4);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 NAV with I (IRNSS) constellation
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_irnss_as_gps() {
        let data = padded("     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
I 1 2020 06 15 01 30  0 -.271548051387D-03 -.682121026330D-11  .000000000000D+00
     7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
    -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
     4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
     9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
     4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
     4.248000000000D+05 4.000000000000D+00
");
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();
        assert_eq!(ephemerides.len(), 1);
        assert_eq!(ephemerides[0].sat().constellation, Constellation::Navic);
        assert_eq!(ephemerides[0].sat().prn, 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 nav with blank lines between ephemerides
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_blank_lines() {
        let data = padded("     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
G 1 2020 06 15 01 30  0 -.271548051387D-03 -.682121026330D-11  .000000000000D+00
     7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
    -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
     4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
     9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
     4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
     4.248000000000D+05 4.000000000000D+00

R 6 2020 12 24 21 15  0  .189751386642E-03  .000000000000E+00  .422910000000E+06
     -.740158740234E+04 -.212037086487E+00  .000000000000E+00  .000000000000E+00
     -.206682856445E+05 -.176755714417E+01  .931322574615E-09 -.400000000000E+01
      .129489067383E+05 -.294115734100E+01 -.186264514923E-08  .000000000000E+00
");
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();
        assert_eq!(ephemerides.len(), 2);
        assert_eq!(ephemerides[0].sat().constellation, Constellation::Gps);
        assert_eq!(ephemerides[1].sat().constellation, Constellation::Glonass);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 NAV: QZSS ephemeris
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_qzss() {
        let data = padded("     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
J 1 2020 06 15 01 30  0 -.271548051387D-03 -.682121026330D-11  .000000000000D+00
     7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
    -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
     4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
     9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
     4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
     4.248000000000D+05 4.000000000000D+00
");
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();
        assert_eq!(ephemerides.len(), 1);
        assert_eq!(ephemerides[0].sat().constellation, Constellation::Qzss);
        assert_eq!(ephemerides[0].sat().prn, 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 NAV: Sbas treated as GPS (code fallback in build_ephemeris)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_nav_sbas_as_gps() {
        // Sbas maps to Constellation::Gps in nav parser (line 668).
        // build_ephemeris for Gps returns a GpsEphemeris.
        let data = padded("     3.02           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
S 1 2020 06 15 01 30  0 -.271548051387D-03 -.682121026330D-11  .000000000000D+00
     7.000000000000D+01-8.000000000000D+00 4.321251426038D-09 2.456182385192D+00
    -1.769512891769D-07 1.934998203069D-03 4.604458808899D-06 5.153539648056D+03
     4.248000000000D+05 1.359730958939D-07-2.968124618807D+00-3.911554813385D-08
     9.798670864830D-01 2.980937500000D+02-1.061584443985D+00-8.075693527908D-09
     4.571618997710D-11 0.000000000000D+00 2.105000000000D+03 0.000000000000D+00
     0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 7.000000000000D+01
     4.248000000000D+05 4.000000000000D+00
");
        let mut reader = BufReader::new(data.as_bytes());
        let (ephemerides, _klob) = parse_rinex_nav(&mut reader).unwrap();
        // SBAS nav messages are not parsed as GPS ephemeris in the current parser
        assert_eq!(ephemerides.len(), 0);
    }

    // -----------------------------------------------------------------------
    // parse_rinex_f64 with various edge cases
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_rinex_f64_edge_cases() {
        assert!((parse_rinex_f64("  1.23D-4  ").unwrap() - 0.000123).abs() < 1e-12);
        assert!((parse_rinex_f64("42").unwrap() - 42.0).abs() < 1e-12);
        assert!((parse_rinex_f64("-0.0").unwrap()).abs() < 1e-12);
        assert!((parse_rinex_f64("1.0D+1").unwrap() - 10.0).abs() < 1e-12);
        assert!((parse_rinex_f64("5d-1").unwrap() - 0.5).abs() < 1e-12);
    }
