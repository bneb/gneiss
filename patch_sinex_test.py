with open("crates/gneiss-parsers/src/sinex_bia.rs", "r") as f:
    content = f.read()

target = """        let content = r#"%=BIA 1.00
+BIAS/SOLUTION
*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___
 OSB  G002 G02           C1W       2021:123:00000 2021:123:86400 ns            1.0000000000    0.000000
 OSB  G002 G02           L1W       2021:123:00000 2021:123:86400 ns            2.0000000000    0.000000
 OSB  G002 G02           C2W       2021:123:00000 2021:123:86400 ns            3.0000000000    0.000000
 OSB  G002 G02           L2W       2021:123:00000 2021:123:86400 ns            4.0000000000    0.000000
 OSB  E002 E02           C2C       2021:123:00000 2021:123:86400 ns            5.0000000000    0.000000
 OSB  E002 E02           L2C       2021:123:00000 2021:123:86400 ns            6.0000000000    0.000000
-BIAS/SOLUTION
"#;
        let cursor = Cursor::new(content);
        let bias = SinexBias::parse(cursor).unwrap();
        let t = GpsTime::new(2156, 129600.0); // 2021:123 at 12:00:00
        
        let g02 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let e02 = SatelliteId { constellation: Constellation::Galileo, prn: 2 };

        // Test GPS Fallbacks
        assert_eq!(bias.get_bias(g02, ObsCode::from_str("C1C").unwrap(), t), Some(1.0));
        assert_eq!(bias.get_bias(g02, ObsCode::from_str("L1C").unwrap(), t), Some(2.0));
        assert_eq!(bias.get_bias(g02, ObsCode::from_str("C2L").unwrap(), t), Some(3.0));
        assert_eq!(bias.get_bias(g02, ObsCode::from_str("L2X").unwrap(), t), Some(4.0));

        // Test non-GPS Fallbacks
        assert_eq!(bias.get_bias(e02, ObsCode::from_str("C2X").unwrap(), t), Some(5.0));
        assert_eq!(bias.get_bias(e02, ObsCode::from_str("L2S").unwrap(), t), Some(6.0));
        
        // Ensure no fallback if direct is present
        assert_eq!(bias.get_exact_bias(g02, ObsCode::from_str("C1C").unwrap(), t), None);"""

replacement = """        let content = r#"%=BIA 1.00\n+BIAS/SOLUTION
*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___
 OSB  G002 G02           C1W       2021:123:00000 2021:123:86400 ns            1.0000000000    0.000000
 OSB  G002 G02           L1W       2021:123:00000 2021:123:86400 ns            2.0000000000    0.000000
 OSB  G002 G02           C2W       2021:123:00000 2021:123:86400 ns            3.0000000000    0.000000
 OSB  G002 G02           L2W       2021:123:00000 2021:123:86400 ns            4.0000000000    0.000000
 OSB  E002 E02           C2C       2021:123:00000 2021:123:86400 ns            5.0000000000    0.000000
 OSB  E002 E02           L2C       2021:123:00000 2021:123:86400 ns            6.0000000000    0.000000
-BIAS/SOLUTION\n"#;
        let bias = SinexBias::parse(Cursor::new(content)).unwrap();
        let t = GpsTime::new(2156, 129600.0);
        let g02 = SatelliteId { constellation: Constellation::Gps, prn: 2 };
        let e02 = SatelliteId { constellation: Constellation::Galileo, prn: 2 };
        assert_eq!(bias.get_bias(g02, ObsCode::from_str("C1C").unwrap(), t), Some(1.0));
        assert_eq!(bias.get_bias(g02, ObsCode::from_str("L1C").unwrap(), t), Some(2.0));
        assert_eq!(bias.get_bias(g02, ObsCode::from_str("C2L").unwrap(), t), Some(3.0));
        assert_eq!(bias.get_bias(g02, ObsCode::from_str("L2X").unwrap(), t), Some(4.0));
        assert_eq!(bias.get_bias(e02, ObsCode::from_str("C2X").unwrap(), t), Some(5.0));
        assert_eq!(bias.get_bias(e02, ObsCode::from_str("L2S").unwrap(), t), Some(6.0));
        assert_eq!(bias.get_exact_bias(g02, ObsCode::from_str("C1C").unwrap(), t), None);"""

if target in content:
    content = content.replace(target, replacement)
    with open("crates/gneiss-parsers/src/sinex_bia.rs", "w") as f:
        f.write(content)
