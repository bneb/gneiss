use bitvec::prelude::*;
use super::RtcmParseError;

/// RTCM Message 1005/1006: Stationary RTK Reference Station ARP
#[derive(Debug, Clone, PartialEq)]
pub struct StationArp {
    pub message_number: u16,
    pub station_id: u16,
    pub itrf_epoch_year: u8,
    pub ecef_x: f64,
    pub ecef_y: f64,
    pub ecef_z: f64,
    pub antenna_height: Option<f64>,
}

/// Parses an RTCM 1005 or 1006 payload.
pub fn parse_station_arp(payload: &[u8]) -> Result<StationArp, RtcmParseError> {
    let bits = payload.view_bits::<Msb0>();
    if bits.len() < 12 {
        return Err(RtcmParseError::Incomplete);
    }

    let mut cursor = 0;

    let message_number = bits[cursor..cursor + 12].load_be::<u16>();
    cursor += 12;

    if message_number != 1005 && message_number != 1006 {
        // Not a Station ARP message
        return Err(RtcmParseError::UnsupportedMsmType);
    }

    if bits.len() < 152 {
        return Err(RtcmParseError::Incomplete);
    }

    let station_id = bits[cursor..cursor + 12].load_be::<u16>();
    cursor += 12;

    let itrf_epoch_year = bits[cursor..cursor + 6].load_be::<u8>();
    cursor += 6;

    // Skip GPS(1), GLONASS(1), Galileo(1), Ref-Station Indicator(1)
    cursor += 4;

    let ecef_x_bits = bits[cursor..cursor + 38].load_be::<u64>();
    cursor += 38;
    
    // sign extend 38 bits to 64
    let shift = 64 - 38;
    let ecef_x_int = (ecef_x_bits << shift) as i64 >> shift;
    let ecef_x = ecef_x_int as f64 * 0.0001;

    // Skip Oscillator Indicator(1), Reserved(1)
    cursor += 2;

    let ecef_y_bits = bits[cursor..cursor + 38].load_be::<u64>();
    cursor += 38;
    
    let ecef_y_int = (ecef_y_bits << shift) as i64 >> shift;
    let ecef_y = ecef_y_int as f64 * 0.0001;

    // Skip Quarter Cycle Indicator(2)
    cursor += 2;

    let ecef_z_bits = bits[cursor..cursor + 38].load_be::<u64>();
    cursor += 38;
    
    let ecef_z_int = (ecef_z_bits << shift) as i64 >> shift;
    let ecef_z = ecef_z_int as f64 * 0.0001;

    let mut antenna_height = None;

    if message_number == 1006 {
        if bits.len() < 168 {
            return Err(RtcmParseError::Incomplete);
        }
        let height_bits = bits[cursor..cursor + 16].load_be::<u16>();
        antenna_height = Some(height_bits as f64 * 0.0001);
    }

    Ok(StationArp {
        message_number,
        station_id,
        itrf_epoch_year,
        ecef_x,
        ecef_y,
        ecef_z,
        antenna_height,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_1005() {
        // Hex representation of a 1005 message payload (without preamble/length/crc)
        // Message 1005 is 19 bytes long (152 bits)
        // Let's create a dummy payload.
        let mut bits = bitvec![u8, Msb0; 0; 152];
        
        // message_number = 1005 (0x3ED)
        bits[0..12].store_be(1005_u16);
        // station_id = 1234
        bits[12..24].store_be(1234_u16);
        // itrf_epoch_year = 20
        bits[24..30].store_be(20_u8);
        
        // ecef_x = -2689639.506 m -> -26,896,395,060 (0.1 mm units)
        let x_int = -26896395060_i64;
        bits[34..72].store_be(x_int as u64); // 38 bits
        
        let y_int = -42904386360_i64;
        bits[74..112].store_be(y_int as u64);
        
        let z_int = 38650509560_i64;
        bits[114..152].store_be(z_int as u64);

        let payload = bits.into_vec();
        
        let arp = parse_station_arp(&payload).unwrap();
        assert_eq!(arp.message_number, 1005);
        assert_eq!(arp.station_id, 1234);
        assert_eq!(arp.itrf_epoch_year, 20);
        assert!((arp.ecef_x - (-2689639.506)).abs() < 1e-4);
        assert!((arp.ecef_y - (-4290438.636)).abs() < 1e-4);
        assert!((arp.ecef_z - 3865050.956).abs() < 1e-4);
        assert_eq!(arp.antenna_height, None);
    }
    #[test]
    fn test_parse_1006() {
        let mut bits = bitvec![u8, Msb0; 0; 168];
        
        // message_number = 1006 (0x3EE)
        bits[0..12].store_be(1006_u16);
        // station_id = 1234
        bits[12..24].store_be(1234_u16);
        // itrf_epoch_year = 21
        bits[24..30].store_be(21_u8);
        
        // ecef_x = 1000.0 m -> 10,000,000 (0.1 mm units)
        let x_int = 10000000_i64;
        bits[34..72].store_be(x_int as u64); // 38 bits
        
        let y_int = 20000000_i64;
        bits[74..112].store_be(y_int as u64);
        
        let z_int = 30000000_i64;
        bits[114..152].store_be(z_int as u64);

        // Antenna height = 1.524 m -> 15240 (0.1 mm units)
        bits[152..168].store_be(15240_u16);

        let payload = bits.into_vec();
        
        let arp = parse_station_arp(&payload).unwrap();
        assert_eq!(arp.message_number, 1006);
        assert_eq!(arp.station_id, 1234);
        assert_eq!(arp.itrf_epoch_year, 21);
        assert!((arp.ecef_x - 1000.0).abs() < 1e-4);
        assert!((arp.ecef_y - 2000.0).abs() < 1e-4);
        assert!((arp.ecef_z - 3000.0).abs() < 1e-4);
        assert_eq!(arp.antenna_height, Some(1.524));
    }
}
