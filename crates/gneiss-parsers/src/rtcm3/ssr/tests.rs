//! Tests for the RTCM 3 SSR parser.

use super::*;

    /// Pack (num_bits, value) pairs into a byte buffer in Msb0 order.
    fn pack_bits(pairs: &[(usize, u64)]) -> Vec<u8> {
        let total_bits: usize = pairs.iter().map(|(b, _)| *b).sum();
        let mut bytes = vec![0u8; total_bits.div_ceil(8)];
        let mut pos = 0;
        for &(bits, val) in pairs {
            for i in 0..bits {
                if (val >> (bits - 1 - i)) & 1 != 0 {
                    bytes[pos / 8] |= 1 << (7 - (pos % 8));
                }
                pos += 1;
            }
        }
        bytes
    }

    // -----------------------------------------------------------------------
    // parse_ssr_header
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_ssr_header_incomplete() {
        let bits = [0u8; 8].view_bits::<Msb0>(); // 64 bits < 68
        let result = parse_ssr_header(bits);
        assert!(matches!(result, Err(RtcmParseError::Incomplete)));
    }

    #[test]
    fn test_parse_ssr_header_valid() {
        let payload = pack_bits(&[
            (12, 1057),   // message_number = SSR Orbit Correction (GPS)
            (20, 123456), // epoch_time
            (4, 5),       // update_interval = 5 sec
            (1, 1),       // multiple_message_indicator
            (1, 0),       // satellite_reference_datum
            (4, 8),       // iod_ssr
            (16, 42),     // provider_id
            (4, 3),       // solution_id
            (6, 12),      // num_satellites = 12
        ]);
        let bits = payload.view_bits::<Msb0>();
        let (remaining, header) = parse_ssr_header(bits).unwrap();
        assert_eq!(header.message_number, 1057);
        assert_eq!(header.epoch_time, 123456);
        assert_eq!(header.update_interval, 5);
        assert!(header.multiple_message_indicator);
        assert!(!header.satellite_reference_datum);
        assert_eq!(header.iod_ssr, 8);
        assert_eq!(header.provider_id, 42);
        assert_eq!(header.solution_id, 3);
        assert_eq!(header.num_satellites, 12);
        assert_eq!(remaining.len(), bits.len() - 68);
    }

    // -----------------------------------------------------------------------
    // parse_ssr_orbit
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_ssr_orbit_empty() {
        // Header with num_satellites = 0, no orbit data
        let payload = pack_bits(&[
            (12, 1057), (20, 0), (4, 0), (1, 0), (1, 0), (4, 0), (16, 0), (4, 0), (6, 0),
        ]);
        let (header, sats) = parse_ssr_orbit(&payload).unwrap();
        assert_eq!(header.num_satellites, 0);
        assert!(sats.is_empty());
    }

    #[test]
    fn test_parse_ssr_orbit_one_sat() {
        // Single pack_bits call: 68 header + 135 sat = 203 bits = 26 bytes = 208 virtual
        let payload = pack_bits(&[
            (12, 1057), (20, 100), (4, 1), (1, 0), (1, 0), (4, 2), (16, 1), (4, 1), (6, 1),
            (6, 5), (8, 17), (22, 100), (20, 200), (20, (-50i32 as u64) & 0xFFFFF),
            (21, 10), (19, 5), (19, (-3i32 as u64) & 0x7FFFF),
        ]);
        let (_header, sats) = parse_ssr_orbit(&payload).unwrap();
        assert_eq!(sats.len(), 1);
        assert_eq!(sats[0].sat_id, 5);
        assert_eq!(sats[0].iode, 17);
    }

    #[test]
    fn test_parse_ssr_orbit_two_sats() {
        // Single pack_bits: 68 header + 135 + 135 = 338 bits = 43 bytes
        let payload = pack_bits(&[
            (12, 1057), (20, 100), (4, 1), (1, 0), (1, 0), (4, 2), (16, 1), (4, 1), (6, 2),
            (6, 1), (8, 10), (22, 50), (20, 100), (20, 30), (21, 5), (19, 2), (19, 1),
            (6, 2), (8, 20), (22, 60), (20, 110), (20, 40), (21, 6), (19, 3), (19, 2),
        ]);
        let (_header, sats) = parse_ssr_orbit(&payload).unwrap();
        assert_eq!(sats.len(), 2);
        assert_eq!(sats[0].sat_id, 1);
        assert_eq!(sats[1].sat_id, 2);
        assert_eq!(sats[0].iode, 10);
        assert_eq!(sats[1].iode, 20);
    }

    #[test]
    fn test_parse_ssr_orbit_incomplete() {
        // Header (68 bits) with num_satellites=1, no orbit data.
        // Total = 68 bits = 9 bytes = 72 bits virtual -> remaining after header = 4 bits
        // parse_orbit_sat checks for 135 bits -> fails
        let payload = pack_bits(&[
            (12, 1057), (20, 0), (4, 0), (1, 0), (1, 0), (4, 0), (16, 0), (4, 0), (6, 1),
        ]);
        let result = parse_ssr_orbit(&payload);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // parse_ssr_clock
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_ssr_clock_empty() {
        let payload = pack_bits(&[
            (12, 1058), (20, 0), (4, 0), (1, 0), (1, 0), (4, 0), (16, 0), (4, 0), (6, 0),
        ]);
        let (header, sats) = parse_ssr_clock(&payload).unwrap();
        assert_eq!(header.message_number, 1058);
        assert!(sats.is_empty());
    }

    #[test]
    fn test_parse_ssr_clock_one_sat() {
        // Single pack_bits: 68 header + 76 clock = 144 bits = 18 bytes = 144 virtual
        let payload = pack_bits(&[
            (12, 1058), (20, 200), (4, 2), (1, 0), (1, 0), (4, 3), (16, 1), (4, 1), (6, 1),
            (6, 7), (22, 500), (21, 10), (27, 5),
        ]);
        let (_header, sats) = parse_ssr_clock(&payload).unwrap();
        assert_eq!(sats.len(), 1);
        assert_eq!(sats[0].sat_id, 7);
        let expected_c0 = sign_extend_i32(500, 22) as f64 * RES_CLOCK_C0;
        let expected_c1 = sign_extend_i32(10, 21) as f64 * RES_CLOCK_C1;
        let expected_c2 = sign_extend_i32(5, 27) as f64 * RES_CLOCK_C2;
        assert!((sats[0].delta_clock_c0 - expected_c0).abs() < 1e-12);
        assert!((sats[0].delta_clock_c1 - expected_c1).abs() < 1e-12);
        assert!((sats[0].delta_clock_c2 - expected_c2).abs() < 1e-12);
    }

    #[test]
    fn test_parse_ssr_clock_incomplete() {
        // Header (68 bits) with num_satellites=1, no sat data = 9 bytes = 72 bits virtual
        // parse_clock_sat checks for 76 bits -> fails
        let payload = pack_bits(&[
            (12, 1058), (20, 0), (4, 0), (1, 0), (1, 0), (4, 0), (16, 0), (4, 0), (6, 1),
        ]);
        let result = parse_ssr_clock(&payload);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // parse_ssr_code_bias
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_ssr_code_bias_empty() {
        let payload = pack_bits(&[
            (12, 1059), (20, 0), (4, 0), (1, 0), (1, 0), (4, 0), (16, 0), (4, 0), (6, 0),
        ]);
        let (header, sats) = parse_ssr_code_bias(&payload).unwrap();
        assert_eq!(header.message_number, 1059);
        assert!(sats.is_empty());
    }

    #[test]
    fn test_parse_ssr_code_bias_one_sat() {
        // Single pack_bits: 68 header + 11 sat_header + 19*2 biases = 117 bits = 15 bytes = 120 virtual
        let payload = pack_bits(&[
            (12, 1059), (20, 300), (4, 1), (1, 0), (1, 0), (4, 1), (16, 1), (4, 1), (6, 1),
            (6, 3), (5, 2),
            (5, 1), (14, 50),
            (5, 2), (14, (-30i32 as u64) & 0x3FFF),
        ]);
        let (_header, sats) = parse_ssr_code_bias(&payload).unwrap();
        assert_eq!(sats.len(), 1);
        assert_eq!(sats[0].sat_id, 3);
        assert_eq!(sats[0].num_biases, 2);
        assert_eq!(sats[0].biases.len(), 2);
        assert_eq!(sats[0].biases[0].signal_and_tracking_mode, 1);
        assert!((sats[0].biases[0].bias - 50.0 * 0.01).abs() < 1e-12);
        assert_eq!(sats[0].biases[1].signal_and_tracking_mode, 2);
        assert!((sats[0].biases[1].bias - (-30.0 * 0.01)).abs() < 1e-12);
    }

    #[test]
    fn test_parse_ssr_code_bias_incomplete_header() {
        // 68 header bits only -> 9 bytes (72 bits virtual).
        // After header: 4 bits remaining. parse_code_bias_sat needs 11 -> Err
        let payload = pack_bits(&[
            (12, 1059), (20, 0), (4, 0), (1, 0), (1, 0), (4, 0), (16, 0), (4, 0), (6, 1),
        ]);
        let result = parse_ssr_code_bias(&payload);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_ssr_code_bias_incomplete_bias() {
        // 68 header + 11 sat_header + 10 bias_start = 89 bits = 12 bytes (96 virtual)
        // After header (68): 28 bits. After sat_header (11): 17 bits.
        // Bias check: 17 < 19 = true -> Err
        let payload = pack_bits(&[
            (12, 1059), (20, 0), (4, 0), (1, 0), (1, 0), (4, 0), (16, 0), (4, 0), (6, 1),
            (6, 1), (5, 1),
            (10, 0),
        ]);
        let result = parse_ssr_code_bias(&payload);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // SSR message number constants
    // -----------------------------------------------------------------------
    #[test]
    fn test_ssr_parse_orbit_with_update_interval_0() {
        // update_interval = 0 is valid (single epoch)
        let payload = pack_bits(&[
            (12, 1057), (20, 0), (4, 0), (1, 0), (1, 0), (4, 0), (16, 0), (4, 0), (6, 1),
            (6, 1), (8, 0), (22, 0), (20, 0), (20, 0), (21, 0), (19, 0), (19, 0),
        ]);
        let (_, sats) = parse_ssr_orbit(&payload).unwrap();
        assert_eq!(sats.len(), 1);
    }


    fn pack(pairs: &[(usize, u64)]) -> Vec<u8> {
        let total: usize = pairs.iter().map(|(b, _)| *b).sum();
        let mut bytes = vec![0u8; total.div_ceil(8)];
        let mut pos = 0;
        for &(bits, val) in pairs {
            for i in 0..bits {
                if (val >> (bits - 1 - i)) & 1 != 0 {
                    bytes[pos / 8] |= 1 << (7 - (pos % 8));
                }
                pos += 1;
            }
        }
        bytes
    }

    fn header(msg: u64, n_sat: u64) -> Vec<(usize, u64)> {
        vec![
            (12, msg),
            (20, 123_456),
            (4, 5),
            (1, 1),
            (1, 0),
            (4, 8),
            (16, 42),
            (4, 3),
            (6, n_sat),
        ]
    }

    /// A GPS SSR orbit record is exactly 135 bits: 6 PRN, 8 IODE,
    /// 22 + 20 + 20 position bits and 21 + 19 + 19 rate bits. Two records must
    /// therefore consume 270 bits and the second must decode cleanly.
    #[test]
    fn orbit_record_is_135_bits() {
        let mut pairs = header(1057, 2);
        for prn in [3u64, 17] {
            pairs.extend_from_slice(&[
                (6, prn),
                (8, 9),
                (22, 1),
                (20, 2),
                (20, 3),
                (21, 4),
                (19, 5),
                (19, 6),
            ]);
        }
        let (hdr, sats) = parse_ssr_orbit(&pack(&pairs)).unwrap();
        assert_eq!(hdr.num_satellites, 2);
        assert_eq!(sats.len(), 2);
        assert_eq!(sats[0].sat_id, 3);
        assert_eq!(sats[1].sat_id, 17);
        // Resolutions: radial 1e-4 m, along/cross 4e-4 m, radial rate 1e-6,
        // along/cross rate 4e-6 m/s (RTKLIB decode_ssr1, rtcm3.c).
        assert!((sats[0].delta_radial - 1.0e-4).abs() < 1e-18);
        assert!((sats[0].delta_along_track - 8.0e-4).abs() < 1e-18);
        assert!((sats[0].delta_cross_track - 1.2e-3).abs() < 1e-18);
        assert!((sats[0].dot_delta_radial - 4.0e-6).abs() < 1e-18);
        assert!((sats[0].dot_delta_along_track - 2.0e-5).abs() < 1e-18);
        assert!((sats[0].dot_delta_cross_track - 2.4e-5).abs() < 1e-18);
    }

    /// The trailing satellite of the declared count must be complete: dropping
    /// its last bit must be an error, not a truncated record.
    #[test]
    fn truncated_last_orbit_record_is_rejected() {
        let mut pairs = header(1057, 1);
        pairs.extend_from_slice(&[
            (6, 3),
            (8, 9),
            (22, 1),
            (20, 2),
            (20, 3),
            (21, 4),
            (19, 5),
            (19, 6),
        ]);
        let full = pack(&pairs);
        // 68 header bits + 135 record bits = 203, so 25 bytes is one bit short.
        assert_eq!(full.len(), 26);
        for nbytes in [25usize, 24, 20] {
            assert_eq!(
                parse_ssr_orbit(&full[..nbytes]),
                Err(RtcmParseError::Incomplete),
                "{nbytes} bytes must be rejected for a 1-satellite orbit message"
            );
        }
        assert!(parse_ssr_orbit(&full).is_ok());
    }

    /// A satellite count larger than the buffer must fail rather than loop over
    /// an exhausted slice.
    #[test]
    fn inflated_satellite_count_is_rejected() {
        let pairs = header(1057, 63);
        assert_eq!(parse_ssr_orbit(&pack(&pairs)), Err(RtcmParseError::Incomplete));
    }

    /// Only the GPS message numbers may be decoded: applying GPS field widths
    /// to a Galileo, GLONASS, QZSS or BeiDou SSR message would silently return
    /// wrong corrections.
    #[test]
    fn only_gps_ssr_messages_are_decoded() {
        let other: Vec<u64> = (1050u64..1100)
            .filter(|m| *m != GPS_SSR_ORBIT as u64
                && *m != GPS_SSR_CLOCK as u64
                && *m != GPS_SSR_CODE_BIAS as u64)
            .collect();
        for msg in other {
            let mut pairs = header(msg, 1);
            // Give every case a body long enough that only the type check can fail.
            pairs.extend_from_slice(&[
                (6, 1), (8, 0), (22, 0), (20, 0), (20, 0), (21, 0), (19, 0), (19, 0),
            ]);
            assert_eq!(
                parse_ssr_orbit(&pack(&pairs)),
                Err(RtcmParseError::UnsupportedMsmType),
                "message {msg} must not be decoded as a GPS SSR orbit message"
            );

            let mut pairs = header(msg, 1);
            pairs.extend_from_slice(&[(6, 1), (22, 0), (21, 0), (27, 0)]);
            assert_eq!(
                parse_ssr_clock(&pack(&pairs)),
                Err(RtcmParseError::UnsupportedMsmType),
                "message {msg} must not be decoded as a GPS SSR clock message"
            );

            let mut pairs = header(msg, 1);
            pairs.extend_from_slice(&[(6, 1), (5, 0)]);
            assert_eq!(
                parse_ssr_code_bias(&pack(&pairs)),
                Err(RtcmParseError::UnsupportedMsmType),
                "message {msg} must not be decoded as a GPS SSR code-bias message"
            );
        }
    }

    #[test]
    fn code_bias_resolution_is_one_centimetre() {
        let mut pairs = header(1059, 1);
        pairs.extend_from_slice(&[(6, 7), (5, 2), (5, 3), (14, 100), (5, 11), (14, (-50i64) as u64)]);
        let (hdr, sats) = parse_ssr_code_bias(&pack(&pairs)).unwrap();
        assert_eq!(hdr.num_satellites, 1);
        assert_eq!(sats[0].sat_id, 7);
        assert_eq!(sats[0].num_biases, 2);
        assert_eq!(sats[0].biases[0].signal_and_tracking_mode, 3);
        assert!((sats[0].biases[0].bias - 1.0).abs() < 1e-15); // 100 * 0.01 m
        assert!((sats[0].biases[1].bias + 0.5).abs() < 1e-15); // -50 * 0.01 m
    }

    #[test]
    fn clock_resolutions_match_the_reference_decoder() {
        let mut pairs = header(1058, 1);
        pairs.extend_from_slice(&[(6, 5), (22, 10), (21, 20), (27, 30)]);
        let (_, sats) = parse_ssr_clock(&pack(&pairs)).unwrap();
        assert_eq!(sats[0].sat_id, 5);
        assert!((sats[0].delta_clock_c0 - 1.0e-3).abs() < 1e-15); // 10 * 1e-4 s
        assert!((sats[0].delta_clock_c1 - 2.0e-5).abs() < 1e-18); // 20 * 1e-6 s/s
        assert!((sats[0].delta_clock_c2 - 6.0e-7).abs() < 1e-20); // 30 * 2e-8 s/s^2
    }
