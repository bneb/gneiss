use super::*;

    #[test]
    fn test_rinex_3_obs_types_list_continuation_ok() {
        let count = 14;
        let mut types_str = "C1C L1C D1C S1C C2L L2L D2L S2L C5Q L5Q D5Q S5Q ".to_string();
        // Continuation line: types start at character position 7 (0-indexed) per code
        let mut lines = vec![
            "G  14  C1W L1W                                           SYS / # / OBS TYPES".to_string(),
        ].into_iter();
        let mut current_line = String::new();

        let types = parse_rinex_3_obs_types_list(
            count, &mut types_str, &mut lines, &mut current_line
        ).unwrap();
        assert_eq!(types.len(), 14);
        assert_eq!(types[0], "C1C");
        assert_eq!(types[12], "C1W");
        assert_eq!(types[13], "L1W");
    }

    #[test]
    fn test_rinex_3_obs_types_list_continuation_err() {
        let count = 14;
        let mut types_str = "C1C L1C D1C S1C C2L L2L D2L S2L C5Q L5Q D5Q S5Q ".to_string();
        let mut lines = vec![
            "G  14 C1W L1W                                            WRONG HEADER     ".to_string(),
        ].into_iter();
        let mut current_line = String::new();

        let result = parse_rinex_3_obs_types_list(
            count, &mut types_str, &mut lines, &mut current_line
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Expected continuation"));
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS with LLI and signal strength
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_with_lli() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     4    C1    L1    D1    S1                              # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1G01
  25140323.324   125140323.32419      2514.032            45.0
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let g01 = &epochs[0].satellites[0];

        let l1_obs = g01.observations.iter()
            .find(|o| o.code.obs_type == ObsType::CarrierPhase && o.code.signal.freq_band == 1)
            .unwrap();
        assert_eq!(l1_obs.lli, Some(1));
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS with many satellites (continuation line)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_many_sats() {
        // Hard-code the data with proper alignment to avoid string continuation issues.
        // 32 spaces before G13G14 for the RINEX 2 sat-list continuation at column 32.
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     2    C1    L1                                                # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  2G01G02
  25140323.324   125140323.324
  25140000.000   125140000.000
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 2);
        assert_eq!(epochs[0].satellites[0].sat.prn, 1);
        assert_eq!(epochs[0].satellites[1].sat.prn, 2);
        assert!((epochs[0].satellites[1].get_observable(1).unwrap() - 25140000.000).abs() < 1e-6);
        assert!((epochs[0].satellites[0].get_observable_phase(1).unwrap() - 125140323.324).abs() < 1e-6);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 OBS: LLI parsing on RINEX 3 data
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_line_lli() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
G    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  0  1
G01  25140323.324   125140323.32419      2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let g01 = &epochs[0].satellites[0];
        let l1c = g01.observations.iter()
            .find(|o| o.code.obs_type == ObsType::CarrierPhase)
            .unwrap();
        assert_eq!(l1c.lli, Some(1));
    }

    // -----------------------------------------------------------------------
    // RINEX 3 OBS epoch with flag > 1 (skip epoch)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_skip_header() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
G    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  2  1
> 2020 06 15 02 00 00.0000000  0  1
G01  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 obs with QZSS
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_qzss() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
J    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  0  1
J01  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1);
        assert_eq!(
            epochs[0].satellites[0].sat,
            SatelliteId { constellation: Constellation::Qzss, prn: 1 }
        );
    }

    // -----------------------------------------------------------------------
    // RINEX 2 obs: P-code (P1, P2) type parsing
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_p_code_types() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    P1    P2    L1                                        # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1G01
  25140323.324   25140323.324   125140323.324
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let g01 = &epochs[0].satellites[0];
        assert!(g01.get_observable(1).is_some());
        assert!(g01.get_observable(2).is_some());
        assert!(g01.get_observable_phase(1).is_some());
    }

    // -----------------------------------------------------------------------
    // RINEX 3 obs: Beidou constellation
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_beidou() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
C    3 C1I L1I D1I                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  0  1
C01  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1);
        assert_eq!(
            epochs[0].satellites[0].sat,
            SatelliteId { constellation: Constellation::Beidou, prn: 1 }
        );
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: Glonass satellite ('R' prefix)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_glonass_sat() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1R06
  22100000.000   121000000.000        2200.000
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1);
        let r06 = &epochs[0].satellites[0];
        assert_eq!(r06.sat.constellation, Constellation::Glonass);
        assert_eq!(r06.sat.prn, 6);
        assert!((r06.get_observable(1).unwrap() - 22100000.000).abs() < 1e-6);
        assert!((r06.get_observable_phase(1).unwrap() - 121000000.000).abs() < 1.0);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: space prefix treated as GPS
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_space_prefix_gps() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1 01
  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites[0].sat.constellation, Constellation::Gps);
        assert_eq!(epochs[0].satellites[0].sat.prn, 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: Sbas satellite ('S' prefix)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_sbas_sat() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1S13
  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let s13 = &epochs[0].satellites[0];
        assert_eq!(s13.sat.constellation, Constellation::Sbas);
        assert_eq!(s13.sat.prn, 13);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: Galileo satellite ('E' prefix)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_galileo_sat() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1E02
  27123456.789   127123456.789        2712.345
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let e02 = &epochs[0].satellites[0];
        assert_eq!(e02.sat.constellation, Constellation::Galileo);
        assert_eq!(e02.sat.prn, 2);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: Beidou satellite ('C' prefix)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_beidou_sat() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1C01
  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let c01 = &epochs[0].satellites[0];
        assert_eq!(c01.sat.constellation, Constellation::Beidou);
        assert_eq!(c01.sat.prn, 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: QZSS satellite ('J' prefix)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_qzss_sat() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1J01
  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let j01 = &epochs[0].satellites[0];
        assert_eq!(j01.sat.constellation, Constellation::Qzss);
        assert_eq!(j01.sat.prn, 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: 6 obs types (multi-line observation values)
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_6_types() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     6    C1    L1    D1    S1    P2    L2              # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1G01
  25140323.324   125140323.324        2514.032            45.0      25140323.324
  125140323.324
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites[0].observations.len(), 6);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: unknown constellation char returns empty sat list
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_unknown_constellation() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     1    C1                                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1X01
  25140323.324
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 0);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 OBS: Sbas constellation
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_sbas() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
S    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> 2020 06 15 01 30 00.0000000  0  1
S01  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1);
        assert_eq!(epochs[0].satellites[0].sat.constellation, Constellation::Sbas);
        assert_eq!(epochs[0].satellites[0].sat.prn, 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 OBS: short epoch header (below 35 chars) is skipped
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_short_epoch_line() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
G    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
> short
> 2020 06 15 01 30 00.0000000  0  1
G01  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        assert_eq!(epochs[0].satellites.len(), 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 3 OBS: epoch with no '>' prefix is skipped
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_3_obs_skips_non_epoch_lines() {
        let data = "     3.03           O: GNSS OBS DATA    M: MIXED            RINEX VERSION / TYPE
G    3 C1C L1C D1C                                             SYS / # / OBS TYPES
                                                            END OF HEADER
random junk line
> 2020 06 15 01 30 00.0000000  0  1
G01  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: year < 80 -> 2000-based
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_year_below_80() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1G01
  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let expected = GpsTime::from_calendar(2020, 5, 14, 22, 0, 0.0);
        assert_eq!(epochs[0].time.week, expected.week);
        assert!((epochs[0].time.tow - expected.tow).abs() < 1e-6);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 OBS: year >= 80 -> 1900-based
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_year_above_80() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     3    C1    L1    D1                                    # / TYPES OF OBSERV
                                                            END OF HEADER
 99 12 25  0  0  0.0000000  0  1G01
  25140323.324   125140323.324        2514.032
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        let expected = GpsTime::from_calendar(1999, 12, 25, 0, 0, 0.0);
        assert_eq!(epochs[0].time.week, expected.week);
    }

    // -----------------------------------------------------------------------
    // RINEX 2 obs: unrecognized obs type is skipped
    // -----------------------------------------------------------------------
    #[test]
    fn test_rinex_2_obs_unrecognized_obs_type() {
        let data = "     2.11           O: GPS OBS DATA    M: Mixed            RINEX VERSION / TYPE
     2    X1    C1                                              # / TYPES OF OBSERV
                                                            END OF HEADER
 20  5 14 22  0  0.0000000  0  1G01
  25140323.324   25140323.324
";
        let mut reader = BufReader::new(data.as_bytes());
        let (epochs, _header) = parse_rinex_obs(&mut reader).unwrap();
        assert_eq!(epochs.len(), 1);
        // X1 is skipped by map_rinex_type, only C1 remains
        assert_eq!(epochs[0].satellites[0].observations.len(), 1);
        let obs = &epochs[0].satellites[0].observations[0];
        assert_eq!(obs.code.obs_type, ObsType::Pseudorange);
        assert_eq!(obs.code.signal.freq_band, 1);
    }
