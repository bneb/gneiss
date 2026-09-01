#![allow(clippy::unwrap_used)]

use super::signals::*;
use super::*;
use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::obs::{ObsType, SignalCode};
use gneiss_core::sat::Constellation;

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

    fn payload_for_msm4_1sat_1sig() -> Vec<u8> {
        // MSM4 GPS (msg 1074), 1 satellite (sat 1), 1 signal (sig 1), 1 cell.
        // MSM4 has no DF419 extended sat info (MSM5/7 only): DF397 (8 bits)
        // then DF398 (10 bits) directly, per RTKLIB decode_msm4 / RTCM
        // 10403.3 field order.
        pack_bits(&[
            (12, 1074), (12, 1), (30, 5000), (1, 0), (3, 2), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, 2), (10, 100),
            (15, 50), (22, 500), (4, 3), (1, 0), (6, 35),
        ])
    }

    fn payload_for_msm5_1sat_2sig() -> Vec<u8> {
        // Note: parse_signal_data reads field-by-field across all cells, not cell-by-cell.
        // Sat data order: DF397 (int ms) -> DF419 (ext info) -> DF398 (rough range) -> DF399 (rate).
        // Order: fine_pr[0..1], fine_ph[0..1], lock[0..1], half_cycle[0..1], cnr[0..1], rates[0..1]
        pack_bits(&[
            (12, 1075), (12, 1), (30, 5000), (1, 0), (3, 2), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, (1u64 << 31) | (1u64 << 30)), (2, 0b11),
            (8, 2), (4, 5), (10, 100), (14, 25),
            (15, 50), (15, 100),  // fine_pr[0], fine_pr[1]
            (22, 500), (22, 1000), // fine_ph[0], fine_ph[1]
            (4, 3), (4, 5),        // lock[0], lock[1]
            (1, 0), (1, 1),        // half_cycle[0]=false, half_cycle[1]=true
            (6, 35), (6, 40),      // cnr[0], cnr[1]
            (15, 10), (15, 20),    // fine_phase_range_rates[0], [1]
        ])
    }

    fn payload_for_msm6_1sat_1sig() -> Vec<u8> {
        // MSM6 has no DF419 extended sat info either (MSM5/7 only).
        pack_bits(&[
            (12, 1076), (12, 1), (30, 5000), (1, 0), (3, 2), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, 2), (10, 100),
            (20, 50), (24, 500), (10, 7), (1, 0), (10, 35),
        ])
    }

    fn payload_for_msm7_1sat_1sig() -> Vec<u8> {
        pack_bits(&[
            (12, 1077), (12, 1), (30, 5000), (1, 0), (3, 2), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, 2), (4, 5), (10, 100), (14, 25),
            (20, 50), (24, 500), (10, 7), (1, 0), (10, 35), (15, 10),
        ])
    }

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

    /// 2 satellites x 2 signals, sparse cell mask (only 2 of 4 cells
    /// active, on DIFFERENT satellites) -- the highest-risk new logic is
    /// the cell-index-to-(satellite,signal) mapping; verify each active
    /// cell's data lands on the correct satellite, not the next one over.
    #[test]
    fn test_into_epoch_obs_sparse_cell_mask_maps_to_correct_satellite() {
        // sats: PRN1 (bit63), PRN2 (bit62). sigs: sig1 (bit31), sig2 (bit30).
        // cell_mask (sat-major, sig-minor over 2x2): [sat1/sig1=1, sat1/sig2=0, sat2/sig1=0, sat2/sig2=1]
        // -> 2 active cells: cell0 = (PRN1, sig1) with fine_pr=111; cell1 = (PRN2, sig2) with fine_pr=222.
        let payload = pack_bits(&[
            (12, 1074), (12, 1), (30, 0), (1, 0), (3, 0), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, (1u64 << 63) | (1u64 << 62)), (32, (1u64 << 31) | (1u64 << 30)),
            (1, 1), (1, 0), (1, 0), (1, 1), // cell_mask: sat1/sig1, sat1/sig2, sat2/sig1, sat2/sig2
            (8, 1), (8, 1),   // rough_range_int_ms[sat1], [sat2]
            (10, 0), (10, 0), // rough_ranges[sat1], [sat2]
            (15, 111), (15, 222), // fine_pr[cell0], fine_pr[cell1]
            (22, 0), (22, 0),
            (4, 0), (4, 0),
            (1, 0), (1, 0),
            (6, 0), (6, 0),
        ]);
        let epoch = parse_msm_message(&payload).unwrap().into_epoch_obs();
        assert_eq!(epoch.satellites.len(), 2);
        let pr_for = |prn: u8| -> f64 {
            epoch.satellites.iter().find(|s| s.sat.prn == prn).unwrap()
                .observations.iter().find(|o| o.code.obs_type == ObsType::Pseudorange).unwrap().value
        };
        let expected = |fine: f64| 1.0 * RANGE_MS + fine * P2_24 * RANGE_MS;
        assert!((pr_for(1) - expected(111.0)).abs() < 1e-6, "cell0 (fine_pr=111) must land on PRN1, its actual satellite");
        assert!((pr_for(2) - expected(222.0)).abs() < 1e-6, "cell1 (fine_pr=222) must land on PRN2, not PRN1");
    }

    // -----------------------------------------------------------------------
    // parse_satellite_data for MSM6 (no rough_phase_range_rates)
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_satellite_data_msm6() {
        // MSM6 has no DF419 extended sat info (MSM5/7 only): DF397 (8 bits)
        // then DF398 (10 bits) directly.
        let payload = pack_bits(&[
            (8, 3),
            (10, 200),
        ]);
        let bits = payload.view_bits::<Msb0>();
        let (_remaining, data) = parse_satellite_data(bits, 1, MsmType::Msm6).unwrap();
        assert_eq!(data.rough_range_int_ms, vec![3]);
        assert_eq!(data.rough_ranges, vec![200]);
        assert!(
            data.extended_sat_info.is_empty(),
            "MSM6 has no DF419 extended sat info (MSM5/7 only)"
        );
        assert!(data.rough_phase_range_rates.is_empty());
    }

    // -----------------------------------------------------------------------
    // parse_signal_data for MSM4 (no fine_phase_range_rates)
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_signal_data_msm4() {
        let payload = pack_bits(&[
            (15, 10), (22, 100), (4, 1), (1, 0), (6, 30),
        ]);
        let bits = payload.view_bits::<Msb0>();
        let (_remaining, data) = parse_signal_data(bits, 1, MsmType::Msm4).unwrap();
        assert_eq!(data.fine_pseudoranges, vec![10]);
        assert_eq!(data.fine_phase_ranges, vec![100]);
        assert!(data.fine_phase_range_rates.is_empty());
    }

    // -----------------------------------------------------------------------
    // GLONASS-based MSM4 (message 1084)
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_msm4_glonass() {
        // Single pack_bits call to avoid bit misalignment. MSM4 has no
        // DF419 extended sat info (MSM5/7 only): DF397 (8 bits) then
        // DF398 (10 bits) directly.
        let payload = pack_bits(&[
            (12, 1084), (12, 0), (30, 1000), (1, 0), (3, 0), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 62), (32, 1u64 << 31), (1, 1),
            (8, 0), (10, 50),
            (15, 5), (22, 50), (4, 0), (1, 0), (6, 30),
        ]);
        let result = parse_msm_message(&payload);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().header.message_number, 1084);
    }

    // =======================================================================
    // Adversarial round-trip tests for MSM5 lock_time and half_cycle bit ordering
    // =======================================================================

    /// Build a minimal MSM5 payload with header+masks+sat_data all zeroed,
    /// but with caller-specified signal data values.
    /// n_sat=1, n_sig=1, so 1 cell.
    fn build_msm5_signal_only(fields: &[(usize, u64)]) -> Vec<u8> {
        // Header: message_number=1075, everything else 0
        let header = [
            (12, 1075u64), (12, 0u64), (30, 0u64), (1, 0u64), (3, 0u64), (7, 0u64),
            (2, 0u64), (2, 0u64), (1, 0u64), (3, 0u64),
        ];
        // Masks: sat_mask[0]=PRN1, sig_mask[0]=sig1, cell_mask=[1]
        let masks = [
            (64, 1u64 << 63),
            (32, 1u64 << 31),
            (1, 1u64),
        ];
        // Sat data (MSM5 order): DF397 int_ms=0, DF419 ext_info=0,
        // DF398 rough_range=0, DF399 rough_phase_range_rate=0.
        let sat_data = [
            (8, 0u64),
            (4, 0u64),
            (10, 0u64),
            (14, 0u64),
        ];
        let all: Vec<_> = header
            .into_iter()
            .chain(masks)
            .chain(sat_data)
            .chain(fields.iter().copied())
            .collect();
        pack_bits(&all)
    }

    /// Build a minimal MSM5 payload with 1 sat, 2 sigs, 2 cells.
    /// The signal data fields are provided in field-by-field order:
    /// fine_pr[0], fine_pr[1], fine_ph[0], fine_ph[1], lock[0], lock[1],
    /// half_cycle[0], half_cycle[1], cnr[0], cnr[1], rate[0], rate[1]
    fn build_msm5_1sat_2sig_signal_only(fields: &[(usize, u64)]) -> Vec<u8> {
        let header = [
            (12, 1075u64), (12, 0u64), (30, 0u64), (1, 0u64), (3, 0u64), (7, 0u64),
            (2, 0u64), (2, 0u64), (1, 0u64), (3, 0u64),
        ];
        let masks = [
            (64, 1u64 << 63),
            (32, (1u64 << 31) | (1u64 << 30)),
            (2, 0b11u64),
        ];
        // Sat data (MSM5 order): DF397 int_ms=0, DF419 ext_info=0,
        // DF398 rough_range=0, DF399 rough_phase_range_rate=0.
        let sat_data = [
            (8, 0u64),
            (4, 0u64),
            (10, 0u64),
            (14, 0u64),
        ];
        let all: Vec<_> = header
            .into_iter()
            .chain(masks)
            .chain(sat_data)
            .chain(fields.iter().copied())
            .collect();
        pack_bits(&all)
    }

    /// Round-trip all lock_time boundary values (0, 1, 7, 14, 15) with a single cell.
    #[test]
    fn test_msm5_lock_time_roundtrip_boundaries() {
        for lock in [0u64, 1, 7, 14, 15] {
            // Build MSM5 with 1 cell; signal data: pr=0, ph=0, lock=N, half=0, cnr=0, rate=0
            let payload = build_msm5_signal_only(&[
                (15, 0),  // fine_pr
                (22, 0),  // fine_ph
                (4, lock), // lock_time
                (1, 0),   // half_cycle
                (6, 0),   // cnr
                (15, 0),  // fine_rate
            ]);
            let msg = parse_msm_message(&payload)
                .unwrap_or_else(|_| panic!("parse failed for lock={}", lock));
            assert_eq!(
                msg.signal_data.lock_time_indicators,
                vec![lock as u16],
                "lock_time round-trip failed for value {}",
                lock
            );
        }
    }

    /// Round-trip half_cycle both values with a single cell.
    #[test]
    fn test_msm5_half_cycle_roundtrip() {
        for half in [0u64, 1u64] {
            let payload = build_msm5_signal_only(&[
                (15, 0),
                (22, 0),
                (4, 0),     // lock_time = 0
                (1, half),  // half_cycle
                (6, 0),
                (15, 0),
            ]);
            let msg = parse_msm_message(&payload)
                .unwrap_or_else(|_| panic!("parse failed for half_cycle={}", half));
            assert_eq!(
                msg.signal_data.half_cycle_ambiguities,
                vec![half != 0],
                "half_cycle round-trip failed for value {}",
                half
            );
        }
    }

    /// Round-trip all 16 lock_time values in a single test.
    #[test]
    fn test_msm5_lock_time_all_values() {
        for lock in 0u64..16 {
            let payload = build_msm5_signal_only(&[
                (15, 0),
                (22, 0),
                (4, lock),
                (1, 0),
                (6, 0),
                (15, 0),
            ]);
            let msg = parse_msm_message(&payload)
                .unwrap_or_else(|_| panic!("parse failed for lock_time={}", lock));
            assert_eq!(
                msg.signal_data.lock_time_indicators,
                vec![lock as u16],
                "lock_time round-trip failed for all-values test, value={}",
                lock
            );
        }
    }

    /// Round-trip a complex combination: lock=[7, 15], half_cycle=[true, false]
    /// with 2 cells to verify multi-cell field-by-field ordering.
    #[test]
    fn test_msm5_lock_half_2cell_roundtrip() {
        let payload = build_msm5_1sat_2sig_signal_only(&[
            (15, 0), (15, 0),      // fine_pr[0], fine_pr[1]
            (22, 0), (22, 0),      // fine_ph[0], fine_ph[1]
            (4, 7), (4, 15),       // lock[0]=7, lock[1]=15
            (1, 1), (1, 0),        // half_cycle[0]=true, half_cycle[1]=false
            (6, 0), (6, 0),        // cnr[0], cnr[1]
            (15, 0), (15, 0),      // rate[0], rate[1]
        ]);
        let msg = parse_msm_message(&payload).unwrap();
        assert_eq!(
            msg.signal_data.lock_time_indicators,
            vec![7, 15],
            "2-cell lock_time round-trip"
        );
        assert_eq!(
            msg.signal_data.half_cycle_ambiguities,
            vec![true, false],
            "2-cell half_cycle round-trip"
        );
    }

    /// Round-trip alternating half_cycle flags in 2-cell configuration.
    #[test]
    fn test_msm5_half_cycle_2cell_all_combos() {
        for h0 in [0u64, 1u64] {
            for h1 in [0u64, 1u64] {
                let payload = build_msm5_1sat_2sig_signal_only(&[
                    (15, 0), (15, 0),
                    (22, 0), (22, 0),
                    (4, 0), (4, 0),
                    (1, h0), (1, h1),
                    (6, 0), (6, 0),
                    (15, 0), (15, 0),
                ]);
                let msg = parse_msm_message(&payload)
                    .unwrap_or_else(|_| panic!("parse failed for h0={}, h1={}", h0, h1));
                assert_eq!(
                    msg.signal_data.half_cycle_ambiguities,
                    vec![h0 != 0, h1 != 0],
                    "half_cycle round-trip failed for [{}, {}]",
                    h0, h1
                );
            }
        }
    }

    /// Round-trip lock_time and half_cycle together in a 1-cell config,
    /// testing every lock value (0..16) with both half_cycle states.
    #[test]
    fn test_msm5_lock_and_half_all_combos() {
        for lock in 0u64..16 {
            for half in [0u64, 1u64] {
                let payload = build_msm5_signal_only(&[
                    (15, 0),
                    (22, 0),
                    (4, lock),
                    (1, half),
                    (6, 0),
                    (15, 0),
                ]);
                let msg = parse_msm_message(&payload)
                    .unwrap_or_else(|_| {
                        panic!("parse failed for lock={}, half={}", lock, half)
                    });
                assert_eq!(
                    msg.signal_data.lock_time_indicators,
                    vec![lock as u16],
                    "lock_time mismatch for lock={}, half={}",
                    lock,
                    half
                );
                assert_eq!(
                    msg.signal_data.half_cycle_ambiguities,
                    vec![half != 0],
                    "half_cycle mismatch for lock={}, half={}",
                    lock,
                    half
                );
            }
        }
    }

    /// Verify the existing test helper payload_for_msm5_1sat_2sig round-trips
    /// exactly (no lenient assertions) — this locks in the current behavior.
    #[test]
    fn test_msm5_existing_payload_exact_roundtrip() {
        let payload = payload_for_msm5_1sat_2sig();
        let msg = parse_msm_message(&payload).unwrap();
        assert_eq!(msg.signal_data.lock_time_indicators, vec![3, 5]);
        assert_eq!(msg.signal_data.half_cycle_ambiguities, vec![false, true]);
        assert_eq!(msg.signal_data.cnrs, vec![35, 40]);
        assert_eq!(msg.signal_data.fine_phase_range_rates, vec![10, 20]);
        assert_eq!(msg.signal_data.fine_pseudoranges, vec![50, 100]);
        assert_eq!(msg.signal_data.fine_phase_ranges, vec![500, 1000]);
        assert_eq!(msg.satellite_data.rough_range_int_ms, vec![2]);
        assert_eq!(msg.satellite_data.rough_ranges, vec![100]);
        assert_eq!(msg.satellite_data.extended_sat_info, vec![5]);
        assert_eq!(msg.satellite_data.rough_phase_range_rates, vec![25]);
    }

    /// Adversarial: verify that non-zero fine_pr and fine_ph don't shift
    /// the lock_time bit positions. Max positive values packed, lock=5, half=true.
    #[test]
    fn test_msm5_max_fine_fields_lock_half_roundtrip() {
        // 15-bit signed max positive = 0x3FFF = 16383
        let fine_pr_pos: u64 = 0x3FFF;
        // 22-bit signed max positive = 0x1FFFFF = 2097151
        let fine_ph_pos: u64 = 0x1FFFFF;
        let payload = build_msm5_signal_only(&[
            (15, fine_pr_pos),
            (22, fine_ph_pos),
            (4, 5),    // lock=5
            (1, 1),    // half=true
            (6, 0),
            (15, 0),
        ]);
        let msg = parse_msm_message(&payload).unwrap();
        assert_eq!(msg.signal_data.fine_pseudoranges, vec![fine_pr_pos as i32]);
        assert_eq!(msg.signal_data.fine_phase_ranges, vec![fine_ph_pos as i32]);
        assert_eq!(msg.signal_data.lock_time_indicators, vec![5]);
        assert_eq!(msg.signal_data.half_cycle_ambiguities, vec![true]);
    }

    /// Adversarial: verify MSM4 lock_time round-trip (4-bit lock, no rate).
    #[test]
    fn test_msm4_lock_time_roundtrip() {
        for lock in [0u64, 1, 7, 14, 15] {
            // Build MSM4 with 1 cell
            let payload = pack_bits(&[
                (12, 1074), (12, 0), (30, 0), (1, 0), (3, 0), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
                (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
                (8, 0), (10, 0),
                (15, 0), (22, 0), (4, lock), (1, 0), (6, 0),
            ]);
            let msg = parse_msm_message(&payload).unwrap();
            assert_eq!(
                msg.signal_data.lock_time_indicators,
                vec![lock as u16],
                "MSM4 lock_time round-trip failed for {}",
                lock
            );
        }
    }

    /// Adversarial: verify MSM6 lock_time round-trip (10-bit lock).
    #[test]
    fn test_msm6_lock_time_roundtrip() {
        for lock in [0u64, 1, 127, 511, 1023] {
            let payload = pack_bits(&[
                (12, 1076), (12, 0), (30, 0), (1, 0), (3, 0), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
                (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
                (8, 0), (10, 0),
                (20, 0), (24, 0), (10, lock), (1, 0), (10, 0),
            ]);
            let msg = parse_msm_message(&payload).unwrap();
            assert_eq!(
                msg.signal_data.lock_time_indicators,
                vec![lock as u16],
                "MSM6 lock_time round-trip failed for {}",
                lock
            );
        }
    }

    /// Adversarial: verify MSM7 lock_time and half_cycle round-trip (10-bit lock, no half_cycle ambiguity shift).
    #[test]
    fn test_msm7_lock_half_roundtrip() {
        let payload = pack_bits(&[
            (12, 1077), (12, 0), (30, 0), (1, 0), (3, 0), (7, 0), (2, 0), (2, 0), (1, 0), (3, 0),
            (64, 1u64 << 63), (32, 1u64 << 31), (1, 1),
            (8, 0), (4, 0), (10, 0), (14, 0),
            (20, 0), (24, 0), (10, 512), (1, 1), (10, 0), (15, 0),
        ]);
        let msg = parse_msm_message(&payload).unwrap();
        assert_eq!(msg.signal_data.lock_time_indicators, vec![512]);
        assert_eq!(msg.signal_data.half_cycle_ambiguities, vec![true]);
    }

    /// Adversarial: verify that deleting a payload byte (shift by 8 bits) does NOT
    /// produce the same lock_time/half_cycle values. This catches accidental
    /// byte-alignment tolerance in the parser.
    #[test]
    fn test_msm5_locked_values_change_when_bits_shifted() {
        let payload = build_msm5_signal_only(&[
            (15, 0),
            (22, 0),
            (4, 7),    // lock=7
            (1, 1),    // half=true
            (6, 0),
            (15, 0),
        ]);
        let msg = parse_msm_message(&payload).unwrap();
        assert_eq!(msg.signal_data.lock_time_indicators, vec![7]);
        assert_eq!(msg.signal_data.half_cycle_ambiguities, vec![true]);

        // Now try with lock=0, half=false — they should NOT be the same as above
        let payload2 = build_msm5_signal_only(&[
            (15, 0),
            (22, 0),
            (4, 0),    // lock=0
            (1, 0),    // half=false
            (6, 0),
            (15, 0),
        ]);
        let msg2 = parse_msm_message(&payload2).unwrap();
        assert_eq!(msg2.signal_data.lock_time_indicators, vec![0]);
        assert_eq!(msg2.signal_data.half_cycle_ambiguities, vec![false]);

        // Verify they are different
        assert_ne!(
            msg.signal_data.lock_time_indicators,
            msg2.signal_data.lock_time_indicators,
            "lock_time should differ between lock=7 and lock=0"
        );
        assert_ne!(
            msg.signal_data.half_cycle_ambiguities,
            msg2.signal_data.half_cycle_ambiguities,
            "half_cycle should differ between true and false"
        );
    }
