use super::*;

    // -----------------------------------------------------------------------
    // MsmType::from_message_number
    // -----------------------------------------------------------------------
    #[test]
    fn test_msm_type_from_message_number() {
        assert_eq!(MsmType::from_message_number(1074), Some(MsmType::Msm4));
        assert_eq!(MsmType::from_message_number(1075), Some(MsmType::Msm5));
        assert_eq!(MsmType::from_message_number(1076), Some(MsmType::Msm6));
        assert_eq!(MsmType::from_message_number(1077), Some(MsmType::Msm7));
        assert_eq!(MsmType::from_message_number(1084), Some(MsmType::Msm4)); // GLONASS
        assert_eq!(MsmType::from_message_number(1095), Some(MsmType::Msm5)); // Galileo
        assert_eq!(MsmType::from_message_number(1127), Some(MsmType::Msm7)); // Beidou
        assert_eq!(MsmType::from_message_number(1071), None);  // MSM1
        assert_eq!(MsmType::from_message_number(1072), None);  // MSM2
        assert_eq!(MsmType::from_message_number(1073), None);  // MSM3
        assert_eq!(MsmType::from_message_number(1078), None);  // Unsupported
    }

    // -----------------------------------------------------------------------
    // parse_msm_header
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_msm_header_incomplete() {
        let bits = [0u8; 9].view_bits::<Msb0>(); // 72 bits < 73
        let result = parse_msm_header(bits);
        assert!(matches!(result, Err(RtcmParseError::Incomplete)));
    }

    #[test]
    fn test_parse_msm_header_valid() {
        let payload = pack_bits(&[
            (12, 1074),
            (12, 42),
            (30, 123456),
            (1, 1),
            (3, 5),
            (7, 0),
            (2, 3),
            (2, 1),
            (1, 0),
            (3, 4),
        ]);
        let bits = payload.view_bits::<Msb0>();
        let (remaining, header) = parse_msm_header(bits).unwrap();
        assert_eq!(header.message_number, 1074);
        assert_eq!(header.station_id, 42);
        assert_eq!(header.epoch_time, 123456);
        assert!(header.multiple_message);
        assert_eq!(header.iods, 5);
        assert_eq!(header.clock_steering, 3);
        assert_eq!(header.external_clock, 1);
        assert!(!header.smoothing_indicator);
        assert_eq!(header.smoothing_interval, 4);
        assert_eq!(remaining.len(), bits.len() - 73);
    }

    // -----------------------------------------------------------------------
    // parse_msm_masks
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_msm_masks_incomplete() {
        let bits = [0u8; 11].view_bits::<Msb0>(); // 88 bits < 96 (minimum for sat+signal masks)
        let result = parse_msm_masks(bits);
        assert!(matches!(result, Err(RtcmParseError::Incomplete)));
    }

    #[test]
    fn test_parse_msm_masks_valid_no_cells() {
        // Satellite mask with zero sats -> 0 cell mask bits
        let payload = pack_bits(&[
            (64, 0),     // satellite_mask: no satellites
            (32, 0),     // signal_mask: no signals
        ]);
        let bits = payload.view_bits::<Msb0>();
        let (remaining, masks) = parse_msm_masks(bits).unwrap();
        assert_eq!(masks.satellite_mask, 0);
        assert_eq!(masks.signal_mask, 0);
        assert!(masks.cell_mask.is_empty()); // 0*0 = 0 cells
        assert_eq!(remaining.len(), bits.len() - 96);
    }

    #[test]
    fn test_parse_msm_masks_missing_cell_bits() {
        // 64 sat bits + 32 sig bits = 96 bits (exactly). With n_sat=1, n_sig=1,
        // need 1 more cell bit, but there are none. Use 31-bit signal field
        // with bit 30 set (MSB of 31-bit field) to ensure signal_mask has 1 bit.
        let payload = pack_bits(&[
            (64, 1u64 << 63),
            (31, 1u64 << 30), // 31 bits: set bit 30 (MSB of field)
            // Only 95 bits total -> 12 bytes (96 bits virtual), 1 bit padding
        ]);
        let bits = payload.view_bits::<Msb0>();
        let result = parse_msm_masks(bits);
        assert!(matches!(result, Err(RtcmParseError::Incomplete)));
    }

    // -----------------------------------------------------------------------
    // parse_satellite_data error paths
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_satellite_data_incomplete_rough_ranges() {
        // 8 bits: exactly enough for DF397 (rough range, integer ms) alone,
        // but MSM4 still needs DF398 (10 more bits) right after it.
        let bits = [0u8; 1].view_bits::<Msb0>();
        let result = parse_satellite_data(bits, 1, MsmType::Msm4);
        assert!(matches!(result, Err(RtcmParseError::Incomplete)));
    }

    #[test]
    fn test_parse_satellite_data_incomplete_ext_info() {
        // DF419 (extended sat info) is MSM5/7 only, so this must use MSM5 to
        // actually exercise that field -- MSM4 never reads it at all.
        // 12 bits: exactly enough for DF397 (8) but not DF419 (4 more).
        let raw = [0u8; 2]; // 16 bits
        let bits = raw.view_bits::<Msb0>();
        let short_bits = &bits[..12];
        let result = parse_satellite_data(short_bits, 1, MsmType::Msm5);
        assert!(matches!(result, Err(RtcmParseError::Incomplete)));
    }

    #[test]
    fn test_parse_satellite_data_incomplete_rough_rate() {
        // 22 bits: exactly enough for DF397 (8) + DF419 (4) + DF398 (10), but
        // not enough for 14-bit DF399 rough_phase_range_rates (MSM5 only).
        let raw = [0u8; 3]; // 24 bits
        let bits = raw.view_bits::<Msb0>();
        let short_bits = &bits[..22];
        let result = parse_satellite_data(short_bits, 1, MsmType::Msm5);
        assert!(matches!(result, Err(RtcmParseError::Incomplete)));
    }

    // -----------------------------------------------------------------------
    // parse_signal_data error paths
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_signal_data_incomplete_fine_pr() {
        let bits = [0u8; 1].view_bits::<Msb0>();
        let result = parse_signal_data(bits, 1, MsmType::Msm4);
        assert!(matches!(result, Err(RtcmParseError::Incomplete)));
    }

    // -----------------------------------------------------------------------
    // parse_msm_message - full MSM4 parse
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_msm4_message() {
        let payload = payload_for_msm4_1sat_1sig();
        let msg = parse_msm_message(&payload).unwrap();
        assert_eq!(msg.msm_type, MsmType::Msm4);
        assert_eq!(msg.header.message_number, 1074);
        assert_eq!(msg.header.station_id, 1);
        assert_eq!(msg.header.epoch_time, 5000);
        assert_eq!(msg.masks.satellite_mask, 1u64 << 63);
        assert_eq!(msg.masks.signal_mask, 1u32 << 31);
        assert_eq!(msg.masks.cell_mask, vec![true]);
        assert_eq!(msg.satellite_data.rough_range_int_ms, vec![2]);
        assert_eq!(msg.satellite_data.rough_ranges, vec![100]);
        assert!(
            msg.satellite_data.extended_sat_info.is_empty(),
            "MSM4 has no DF419 extended sat info (MSM5/7 only)"
        );
        assert!(msg.satellite_data.rough_phase_range_rates.is_empty());
        assert_eq!(msg.signal_data.fine_pseudoranges, vec![50]);
        assert_eq!(msg.signal_data.fine_phase_ranges, vec![500]);
        assert_eq!(msg.signal_data.lock_time_indicators, vec![3]);
        assert_eq!(msg.signal_data.half_cycle_ambiguities, vec![false]);
        assert_eq!(msg.signal_data.cnrs, vec![35]);
        assert!(msg.signal_data.fine_phase_range_rates.is_empty());
    }

    // -----------------------------------------------------------------------
    // parse_msm_message - full MSM5 parse (with range rates)
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_msm5_message() {
        let payload = payload_for_msm5_1sat_2sig();
        let msg = parse_msm_message(&payload).unwrap();
        assert_eq!(msg.msm_type, MsmType::Msm5);
        assert_eq!(msg.header.message_number, 1075);
        assert_eq!(msg.satellite_data.rough_range_int_ms, vec![2]);
        assert_eq!(msg.satellite_data.rough_ranges, vec![100]);
        assert_eq!(msg.satellite_data.extended_sat_info, vec![5]);
        assert_eq!(msg.satellite_data.rough_phase_range_rates, vec![25]);
        assert_eq!(msg.signal_data.fine_pseudoranges.len(), 2);
        assert_eq!(msg.signal_data.fine_phase_ranges.len(), 2);
        assert_eq!(msg.signal_data.lock_time_indicators, vec![3, 5]);
        assert_eq!(msg.signal_data.half_cycle_ambiguities, vec![false, true]);
        assert_eq!(msg.signal_data.cnrs, vec![35, 40]);
        assert_eq!(msg.signal_data.fine_phase_range_rates, vec![10, 20]);
    }

    // -----------------------------------------------------------------------
    // parse_msm_message - full MSM6 parse (high-res, no rates)
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_msm6_message() {
        let payload = payload_for_msm6_1sat_1sig();
        let msg = parse_msm_message(&payload).unwrap();
        assert_eq!(msg.msm_type, MsmType::Msm6);
        assert_eq!(msg.header.message_number, 1076);
        assert_eq!(msg.signal_data.fine_pseudoranges, vec![50]);
        assert_eq!(msg.signal_data.fine_phase_ranges, vec![500]);
        assert_eq!(msg.signal_data.lock_time_indicators, vec![7]);
        assert_eq!(msg.signal_data.cnrs, vec![35]);
        assert!(msg.signal_data.fine_phase_range_rates.is_empty());
        assert!(msg.satellite_data.rough_phase_range_rates.is_empty());
    }

    // -----------------------------------------------------------------------
    // parse_msm_message - full MSM7 parse (high-res with rates)
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_msm7_message() {
        let payload = payload_for_msm7_1sat_1sig();
        let msg = parse_msm_message(&payload).unwrap();
        assert_eq!(msg.msm_type, MsmType::Msm7);
        assert_eq!(msg.header.message_number, 1077);
        assert_eq!(msg.satellite_data.rough_phase_range_rates, vec![25]);
        assert_eq!(msg.signal_data.fine_phase_range_rates, vec![10]);
    }

    // -----------------------------------------------------------------------
    // parse_msm_message - unsupported MSM type (MSM1)
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_msm_message_unsupported_type() {
        let payload = pack_bits(&[
            (12, 1071), // MSM1 (unsupported)
            (12, 0), (30, 0), (1, 0), (3, 0), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
        ]);
        let result = parse_msm_message(&payload);
        assert!(matches!(result, Err(RtcmParseError::UnsupportedMsmType)));
    }

    // -----------------------------------------------------------------------
    // MsmMessage::into_epoch_obs
    // -----------------------------------------------------------------------
    #[test]
    fn test_msm_into_epoch_obs() {
        let payload = payload_for_msm4_1sat_1sig();
        let msg = parse_msm_message(&payload).unwrap();
        let epoch = msg.into_epoch_obs();
        assert_eq!(epoch.satellites.len(), 1);
        assert_eq!(epoch.satellites[0].sat.prn, 1);
        // Constellation derived from message_number 1074/10 = 107 -> GPS
        assert_eq!(epoch.satellites[0].sat.constellation, Constellation::Gps);
        // rough_range_int_ms=2, rough_ranges(modulo)=100, fine_pr=50 (MSM4:
        // P2_24 scale) -- same formula as into_epoch_obs's own doc comment,
        // checked here so a typo in the implementation shows up as a test
        // failure rather than a silently different number.
        let expected_pr = 2.0 * RANGE_MS + 100.0 * P2_10 * RANGE_MS + 50.0 * P2_24 * RANGE_MS;
        let obs = &epoch.satellites[0].observations;
        let pr = obs.iter().find(|o| o.code.obs_type == ObsType::Pseudorange)
            .expect("MSM4 1sat/1sig payload has a valid cell, must produce a pseudorange");
        assert!((pr.value - expected_pr).abs() < 1e-6, "pr={} expected={}", pr.value, expected_pr);
        assert_eq!(pr.code.signal, SignalCode { freq_band: 1, attribute: 'C' }, "GPS sig_id=1 -> RINEX 1C");

        let expected_range_m = 2.0 * RANGE_MS + 100.0 * P2_10 * RANGE_MS + 500.0 * P2_29 * RANGE_MS;
        let expected_cycles = expected_range_m * gneiss_core::frequencies::track_c_frequency(Constellation::Gps, 1, 0) / SPEED_OF_LIGHT_M_S;
        let ph = obs.iter().find(|o| o.code.obs_type == ObsType::CarrierPhase)
            .expect("MSM4 1sat/1sig payload has a valid cell, must produce a carrier phase");
        assert!((ph.value - expected_cycles).abs() < 1e-6, "ph={} expected={}", ph.value, expected_cycles);
    }

    // -----------------------------------------------------------------------
    // into_epoch_obs: hand-verified formula edge cases
    // -----------------------------------------------------------------------

    /// All-zero rough range + fine pseudorange must produce exactly 0.0 --
    /// the simplest possible check that isn't sensitive to a units/scale
    /// error that happens to cancel out at nonzero values.
    #[test]
    fn test_into_epoch_obs_zero_range_is_zero() {
        let payload = pack_bits(&[
            (12, 1074), (12, 1), (30, 0), (1, 0), (3, 0), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, 0), (10, 0),
            (15, 0), (22, 0), (4, 0), (1, 0), (6, 0),
        ]);
        let epoch = parse_msm_message(&payload).unwrap().into_epoch_obs();
        let pr = epoch.satellites[0].observations.iter()
            .find(|o| o.code.obs_type == ObsType::Pseudorange).unwrap();
        assert_eq!(pr.value, 0.0);
        let ph = epoch.satellites[0].observations.iter()
            .find(|o| o.code.obs_type == ObsType::CarrierPhase).unwrap();
        assert_eq!(ph.value, 0.0);
    }

    /// DF397 == 255 (the "not available" sentinel) must suppress BOTH
    /// pseudorange and carrier phase for that satellite -- not just leave
    /// the rough part out of the sum.
    #[test]
    fn test_into_epoch_obs_rough_range_sentinel_suppresses_observations() {
        let payload = pack_bits(&[
            (12, 1074), (12, 1), (30, 0), (1, 0), (3, 0), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, 255), (10, 100),
            (15, 50), (22, 500), (4, 0), (1, 0), (6, 0),
        ]);
        let epoch = parse_msm_message(&payload).unwrap().into_epoch_obs();
        assert!(epoch.satellites[0].observations.is_empty());
    }

    /// Fine-pseudorange sentinel (-16384 for MSM4/5's 15-bit field) must
    /// suppress only the pseudorange, not the carrier phase computed from
    /// a separately-valid fine phase-range in the same cell.
    #[test]
    fn test_into_epoch_obs_fine_pseudorange_sentinel_suppresses_only_pr() {
        // -16384 is the most-negative representable 15-bit two's-complement
        // value (-2^14): its bit pattern is the sign bit alone, 1<<14.
        let sentinel_15bit: u64 = 1 << 14;
        let payload = pack_bits(&[
            (12, 1074), (12, 1), (30, 0), (1, 0), (3, 0), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, 1), (10, 0),
            (15, sentinel_15bit), (22, 0), (4, 0), (1, 0), (6, 0),
        ]);
        let epoch = parse_msm_message(&payload).unwrap().into_epoch_obs();
        let obs = &epoch.satellites[0].observations;
        assert!(obs.iter().all(|o| o.code.obs_type != ObsType::Pseudorange), "PR sentinel must suppress the pseudorange");
        assert!(obs.iter().any(|o| o.code.obs_type == ObsType::CarrierPhase), "phase is independently valid, must still be emitted");
    }

    /// GLONASS never emits carrier phase (FDMA channel not reliably known
    /// from MSM) but pseudorange, which needs no frequency, still works.
    #[test]
    fn test_into_epoch_obs_glonass_phase_is_scoped_out() {
        let payload = pack_bits(&[
            (12, 1084), (12, 1), (30, 0), (1, 0), (3, 0), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, 1), (10, 0),
            (15, 0), (22, 0), (4, 0), (1, 0), (6, 0),
        ]);
        let epoch = parse_msm_message(&payload).unwrap().into_epoch_obs();
        let obs = &epoch.satellites[0].observations;
        assert!(obs.iter().any(|o| o.code.obs_type == ObsType::Pseudorange), "GLONASS pseudorange needs no frequency, must still work");
        assert!(obs.iter().all(|o| o.code.obs_type != ObsType::CarrierPhase), "GLONASS phase is deliberately scoped out, see doc comment");
    }
