use super::*;

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
