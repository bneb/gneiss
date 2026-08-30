#![allow(clippy::unwrap_used)]

    use super::*;

    #[test]
    fn test_ubx_parsing() {
        // Example UBX-ACK-ACK frame
        let msg = [0xB5, 0x62, 0x05, 0x01, 0x02, 0x00, 0x06, 0x01, 0x0F, 0x38];
        let (rem, frame) = parse_ubx_frame(&msg).unwrap();
        assert_eq!(rem.len(), 0);
        assert_eq!(frame.class, 0x05);
        assert_eq!(frame.id, 0x01);
        assert_eq!(frame.payload, &[0x06, 0x01]);
    }

    #[test]
    fn test_parse_rxm_rawx() {
        let mut payload = Vec::new();
        // 16 byte header
        payload.extend_from_slice(&100000.5f64.to_le_bytes()); // rcv_tow
        payload.extend_from_slice(&2100u16.to_le_bytes()); // week
        payload.push(18); // leapS
        payload.push(1); // numMeas = 1
        payload.push(0x01); // recStat
        payload.push(0x01); // version
        payload.extend_from_slice(&[0, 0]); // reserved

        // 32 byte measurement block
        payload.extend_from_slice(&20000000.5f64.to_le_bytes()); // prMes
        payload.extend_from_slice(&1000000.5f64.to_le_bytes()); // cpMes
        payload.extend_from_slice(&123.4f32.to_le_bytes()); // doMes
        payload.push(0); // gnssId (GPS)
        payload.push(12); // svId (PRN 12)
        payload.push(0); // sigId (L1C/A)
        payload.push(0); // freqId
        payload.extend_from_slice(&5000u16.to_le_bytes()); // locktime
        payload.push(45); // cno
        payload.push(2); // prStdev (0.01 * 2^2 = 0.04)
        payload.push(3); // cpStdev (0.004 * 2^3 = 0.032)
        payload.push(1); // doStdev (0.002 * 2^1 = 0.004)
        payload.push(0x03); // trkStat (prValid | cpValid)
        payload.push(0); // reserved

        let rawx = parse_rxm_rawx(&payload).unwrap();
        assert_eq!(rawx.rcv_tow, 100000.5);
        assert_eq!(rawx.week, 2100);
        assert_eq!(rawx.num_meas, 1);

        let meas = &rawx.measurements[0];
        assert_eq!(meas.pr_mes, 20000000.5);
        assert_eq!(meas.sv_id, 12);
        assert_eq!(meas.cno, 45);
        assert!(meas.pr_valid);
        assert!(meas.cp_valid);
        assert!(!meas.half_cycle_valid);
        assert_eq!(meas.pr_stdev, 0.04);
        assert_eq!(meas.cp_stdev, 0.032);
    }

    #[test]
    fn test_parse_rxm_sfrbx() {
        let mut payload = vec![
            0,  // gnssId = GPS
            12, // svId = 12
            0,  // sigId = L1C/A
            0,  // freqId
            2,  // numWords = 2
            1,  // chn = 1
            2,  // version = 2
            0,  // reserved
        ];

        // Words
        payload.extend_from_slice(&0xDEADBEEFu32.to_le_bytes());
        payload.extend_from_slice(&0xCAFEBABEu32.to_le_bytes());

        let sfrbx = parse_rxm_sfrbx(&payload).unwrap();
        assert_eq!(sfrbx.gnss_id, 0);
        assert_eq!(sfrbx.sv_id, 12);
        assert_eq!(sfrbx.num_words, 2);
        assert_eq!(sfrbx.words.len(), 2);
        assert_eq!(sfrbx.words[0], 0xDEADBEEF);
        assert_eq!(sfrbx.words[1], 0xCAFEBABE);
    }

    #[test]
    fn test_parse_esf_meas() {
        let mut payload = Vec::new();
        // 8 byte header
        payload.extend_from_slice(&1234567u32.to_le_bytes()); // time_tag
        payload.extend_from_slice(&0x0001u16.to_le_bytes()); // flags (time mark sent)
        payload.extend_from_slice(&0x0000u16.to_le_bytes()); // id (provider)

        // 1 measurement (4 bytes)
        // Data: 0xABCDEF (24 bits)
        // DataType: 5 (z-axis gyro) (5 << 24 = 0x05000000)
        let meas_raw = 0xABCDEF | (5u32 << 24);
        payload.extend_from_slice(&meas_raw.to_le_bytes());

        let esf_meas = parse_esf_meas(&payload).unwrap();
        assert_eq!(esf_meas.time_tag, 1234567);
        assert_eq!(esf_meas.flags, 1);
        assert_eq!(esf_meas.id, 0);
        assert_eq!(esf_meas.measurements.len(), 1);

        let meas = &esf_meas.measurements[0];
        assert_eq!(meas.data_type, 5);
        assert_eq!(meas.data, 0xFFABCDEF_u32 as i32);
    }

    #[test]
    fn test_parse_esf_status() {
        let mut payload = Vec::new();
        // 16 byte header minimum
        payload.extend_from_slice(&9876543u32.to_le_bytes()); // time_tag
        payload.push(0); // version
        payload.push(1); // fusion_mode (1 = Fusion)
        payload.extend_from_slice(&[0; 5]); // reserved
        payload.push(2); // numSensors
        payload.extend_from_slice(&[0; 4]); // reserved

        let esf_status = parse_esf_status(&payload).unwrap();
        assert_eq!(esf_status.time_tag, 9876543);
        assert_eq!(esf_status.fusion_mode, 1);
        assert_eq!(esf_status.num_sensors, 2);
    }

    // -----------------------------------------------------------------------
    // parse_ubx_frame error paths
    // -----------------------------------------------------------------------
    #[test]
    fn test_ubx_frame_incomplete_short() {
        let msg = [0xB5];
        let result = parse_ubx_frame(&msg);
        assert!(matches!(result, Err(UbxParseError::Incomplete)));

        // Exactly 7 bytes (needs 8 minimum)
        let msg = [0xB5, 0x62, 0x05, 0x01, 0x02, 0x00, 0x06];
        let result = parse_ubx_frame(&msg);
        assert!(matches!(result, Err(UbxParseError::Incomplete)));
    }

    #[test]
    fn test_ubx_frame_invalid_sync() {
        let msg = [0x00, 0x00, 0x05, 0x01, 0x02, 0x00, 0x06, 0x01, 0x0F, 0x38];
        let result = parse_ubx_frame(&msg);
        assert!(matches!(result, Err(UbxParseError::InvalidSync)));
    }

    #[test]
    fn test_ubx_frame_checksum_mismatch() {
        // Good frame with bad checksum byte
        let msg = [0xB5, 0x62, 0x05, 0x01, 0x02, 0x00, 0x06, 0x01, 0xFF, 0xFF];
        let result = parse_ubx_frame(&msg);
        assert!(matches!(result, Err(UbxParseError::ChecksumMismatch)));
    }

    #[test]
    fn test_ubx_frame_incomplete_payload() {
        // Claims 4 bytes of payload but only has 1
        let msg = [0xB5, 0x62, 0x05, 0x01, 0x04, 0x00, 0x06];
        let result = parse_ubx_frame(&msg);
        assert!(matches!(result, Err(UbxParseError::Incomplete)));
    }

    // -----------------------------------------------------------------------
    // ubx_checksum verification
    // -----------------------------------------------------------------------
    #[test]
    fn test_ubx_checksum_known_values() {
        // From UBX-ACK-ACK frame (class=0x05, id=0x01, payload=[0x06, 0x01])
        let data = [0x05, 0x01, 0x02, 0x00, 0x06, 0x01];
        let (ck_a, ck_b) = ubx_checksum(&data);
        assert_eq!(ck_a, 0x0F);
        assert_eq!(ck_b, 0x38);
    }

    #[test]
    fn test_ubx_checksum_empty() {
        let (ck_a, ck_b) = ubx_checksum(&[]);
        assert_eq!(ck_a, 0);
        assert_eq!(ck_b, 0);
    }

    // -----------------------------------------------------------------------
    // parse_rxm_rawx error paths
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_rxm_rawx_invalid_length() {
        // Too short (< 16)
        let payload = vec![0u8; 10];
        let result = parse_rxm_rawx(&payload);
        assert!(matches!(result, Err(UbxParseError::InvalidLength)));

        // Payload length doesn't match expected (num_meas=1 but only 16+32-1 bytes)
        let mut payload = vec![0u8; 16 + 31];
        payload[11] = 1; // num_meas = 1
        let result = parse_rxm_rawx(&payload);
        assert!(matches!(result, Err(UbxParseError::InvalidLength)));
    }

    // -----------------------------------------------------------------------
    // parse_rxm_sfrbx error paths
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_rxm_sfrbx_invalid_length() {
        // Too short (< 8)
        let payload = vec![0u8; 5];
        let result = parse_rxm_sfrbx(&payload);
        assert!(matches!(result, Err(UbxParseError::InvalidLength)));

        // num_words=3 but payload doesn't match (8 + 3*4 = 20)
        let mut payload = vec![0u8; 19];
        payload[4] = 3; // num_words = 3
        let result = parse_rxm_sfrbx(&payload);
        assert!(matches!(result, Err(UbxParseError::InvalidLength)));
    }

    // -----------------------------------------------------------------------
    // parse_esf_meas error paths
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_esf_meas_invalid_length() {
        // Too short (< 8)
        let payload = vec![0u8; 5];
        let result = parse_esf_meas(&payload);
        assert!(matches!(result, Err(UbxParseError::InvalidLength)));

        // Non-multiple of 4 after header
        let payload = vec![0u8; 9];
        let result = parse_esf_meas(&payload);
        assert!(matches!(result, Err(UbxParseError::InvalidLength)));
    }

    // -----------------------------------------------------------------------
    // parse_esf_status error paths
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_esf_status_invalid_length() {
        // Too short (< 16)
        let payload = vec![0u8; 10];
        let result = parse_esf_status(&payload);
        assert!(matches!(result, Err(UbxParseError::InvalidLength)));
    }

    // -----------------------------------------------------------------------
    // EsfMeasData::scaled_value for all sensor types
    // -----------------------------------------------------------------------
    #[test]
    fn test_esf_meas_data_scaled_gyro() {
        // Type 5: z-axis gyro, val=1024 -> 1024 * 2^-12 * pi/180
        let m = EsfMeasData { data: 1024, data_type: 5 };
        let expected = 1024.0 * 2.0f64.powi(-12) * (std::f64::consts::PI / 180.0);
        assert!((m.scaled_value() - expected).abs() < 1e-12);

        // Type 13: also gyro
        let m = EsfMeasData { data: -512, data_type: 13 };
        assert!((m.scaled_value() - (-512.0 * 2.0f64.powi(-12) * (std::f64::consts::PI / 180.0))).abs() < 1e-12);

        // Type 14: also gyro
        let m = EsfMeasData { data: 2048, data_type: 14 };
        assert!((m.scaled_value() - (2048.0 * 2.0f64.powi(-12) * (std::f64::consts::PI / 180.0))).abs() < 1e-12);
    }

    #[test]
    fn test_esf_meas_data_scaled_accel() {
        // Type 16: accelerometer x-axis, val=512 -> 512 * 2^-10
        let m = EsfMeasData { data: 512, data_type: 16 };
        assert!((m.scaled_value() - (512.0 * 2.0f64.powi(-10))).abs() < 1e-12);

        // Type 17
        let m = EsfMeasData { data: -256, data_type: 17 };
        assert!((m.scaled_value() - (-256.0 * 2.0f64.powi(-10))).abs() < 1e-12);

        // Type 18
        let m = EsfMeasData { data: 1024, data_type: 18 };
        assert!((m.scaled_value() - (1024.0 * 2.0f64.powi(-10))).abs() < 1e-12);
    }

    #[test]
    fn test_esf_meas_data_scaled_other() {
        // Type 0: unknown returns raw value
        let m = EsfMeasData { data: 12345, data_type: 0 };
        assert!((m.scaled_value() - 12345.0).abs() < 1e-12);

        // Type 255: unknown returns raw value
        let m = EsfMeasData { data: -42, data_type: 255 };
        assert!((m.scaled_value() - (-42.0)).abs() < 1e-12);
    }

    // -----------------------------------------------------------------------
    // UbxRxmRawx::into_epoch_obs with various gnss_ids
    // -----------------------------------------------------------------------
    fn make_rawx_measurement(gnss_id: u8, sv_id: u8, sig_id: u8) -> RxmRawxMeas {
        RxmRawxMeas {
            pr_mes: 20000000.0,
            cp_mes: 100000000.0,
            do_mes: 1000.0,
            gnss_id,
            sv_id,
            sig_id,
            freq_id: 0,
            locktime: 100,
            cno: 45,
            pr_stdev: 0.01,
            cp_stdev: 0.004,
            do_stdev: 0.002,
            pr_valid: true,
            cp_valid: true,
            half_cycle_valid: true,
            sub_half_cycle: false,
        }
    }

    #[test]
    fn test_ubx_into_epoch_obs_gps() {
        let rawx = UbxRxmRawx {
            rcv_tow: 100000.0,
            week: 2100,
            leap_s: 18,
            num_meas: 1,
            rec_stat: 0,
            version: 1,
            measurements: vec![make_rawx_measurement(0, 12, 0)], // gnss_id=0=GPS
        };
        let epoch = rawx.into_epoch_obs();
        assert_eq!(epoch.satellites.len(), 1);
        assert_eq!(epoch.satellites[0].sat.constellation, Constellation::Gps);
        assert_eq!(epoch.satellites[0].sat.prn, 12);
        assert_eq!(epoch.satellites[0].observations.len(), 4); // pr + cp + do + snr
    }

    #[test]
    fn test_ubx_into_epoch_obs_glonass() {
        let rawx = UbxRxmRawx {
            rcv_tow: 100000.0,
            week: 2100,
            leap_s: 18,
            num_meas: 1,
            rec_stat: 0,
            version: 1,
            measurements: vec![make_rawx_measurement(6, 5, 0)], // gnss_id=6=GLONASS
        };
        let epoch = rawx.into_epoch_obs();
        assert_eq!(epoch.satellites[0].sat.constellation, Constellation::Glonass);
        assert_eq!(epoch.satellites[0].sat.prn, 5);
    }

    #[test]
    fn test_ubx_into_epoch_obs_galileo() {
        let rawx = UbxRxmRawx {
            rcv_tow: 100000.0,
            week: 2100,
            leap_s: 18,
            num_meas: 1,
            rec_stat: 0,
            version: 1,
            measurements: vec![make_rawx_measurement(2, 1, 0)], // gnss_id=2=Galileo
        };
        let epoch = rawx.into_epoch_obs();
        assert_eq!(epoch.satellites[0].sat.constellation, Constellation::Galileo);
    }

    #[test]
    fn test_ubx_into_epoch_obs_beidou() {
        let rawx = UbxRxmRawx {
            rcv_tow: 100000.0,
            week: 2100,
            leap_s: 18,
            num_meas: 1,
            rec_stat: 0,
            version: 1,
            measurements: vec![make_rawx_measurement(3, 10, 0)], // gnss_id=3=Beidou
        };
        let epoch = rawx.into_epoch_obs();
        assert_eq!(epoch.satellites[0].sat.constellation, Constellation::Beidou);
    }

    #[test]
    fn test_ubx_into_epoch_obs_qzss() {
        let rawx = UbxRxmRawx {
            rcv_tow: 100000.0,
            week: 2100,
            leap_s: 18,
            num_meas: 1,
            rec_stat: 0,
            version: 1,
            measurements: vec![make_rawx_measurement(5, 7, 0)], // gnss_id=5=QZSS
        };
        let epoch = rawx.into_epoch_obs();
        assert_eq!(epoch.satellites[0].sat.constellation, Constellation::Qzss);
    }

    #[test]
    fn test_ubx_into_epoch_obs_sbas() {
        let rawx = UbxRxmRawx {
            rcv_tow: 100000.0,
            week: 2100,
            leap_s: 18,
            num_meas: 1,
            rec_stat: 0,
            version: 1,
            measurements: vec![make_rawx_measurement(1, 20, 0)], // gnss_id=1=SBAS
        };
        let epoch = rawx.into_epoch_obs();
        assert_eq!(epoch.satellites[0].sat.constellation, Constellation::Sbas);
    }

    #[test]
    fn test_ubx_into_epoch_obs_unknown_gnss_skipped() {
        let rawx = UbxRxmRawx {
            rcv_tow: 100000.0,
            week: 2100,
            leap_s: 18,
            num_meas: 2,
            rec_stat: 0,
            version: 1,
            measurements: vec![
                make_rawx_measurement(0, 1, 0),   // GPS (valid)
                make_rawx_measurement(99, 2, 0),  // Unknown (skipped)
            ],
        };
        let epoch = rawx.into_epoch_obs();
        assert_eq!(epoch.satellites.len(), 1);
        assert_eq!(epoch.satellites[0].sat.prn, 1);
    }

    #[test]
    fn test_ubx_into_epoch_obs_pr_invalid() {
        let mut meas = make_rawx_measurement(0, 1, 0);
        meas.pr_valid = false;
        let rawx = UbxRxmRawx {
            rcv_tow: 100000.0,
            week: 2100,
            leap_s: 18,
            num_meas: 1,
            rec_stat: 0,
            version: 1,
            measurements: vec![meas],
        };
        let epoch = rawx.into_epoch_obs();
        // No pseudorange observation if pr_valid is false
        let has_pr = epoch.satellites[0].observations.iter()
            .any(|o| o.code.obs_type == ObsType::Pseudorange);
        assert!(!has_pr);
    }

    #[test]
    fn test_ubx_into_epoch_obs_cp_invalid() {
        let mut meas = make_rawx_measurement(0, 1, 0);
        meas.cp_valid = false;
        let rawx = UbxRxmRawx {
            rcv_tow: 100000.0,
            week: 2100,
            leap_s: 18,
            num_meas: 1,
            rec_stat: 0,
            version: 1,
            measurements: vec![meas],
        };
        let epoch = rawx.into_epoch_obs();
        // No carrier phase if cp_valid is false
        let has_cp = epoch.satellites[0].observations.iter()
            .any(|o| o.code.obs_type == ObsType::CarrierPhase);
        assert!(!has_cp);
    }

    #[test]
    fn test_ubx_into_epoch_obs_sig_id_2_band() {
        let mut meas = make_rawx_measurement(0, 1, 0);
        meas.sig_id = 1; // sig_id != 0 -> freq_band=2
        let rawx = UbxRxmRawx {
            rcv_tow: 100000.0,
            week: 2100,
            leap_s: 18,
            num_meas: 1,
            rec_stat: 0,
            version: 1,
            measurements: vec![meas],
        };
        let epoch = rawx.into_epoch_obs();
        for obs in &epoch.satellites[0].observations {
            assert_eq!(obs.code.signal.freq_band, 2);
        }
    }
