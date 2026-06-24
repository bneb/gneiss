use super::{sign_extend_i32, RtcmParseError};
use bitvec::prelude::*;

const RES_ORBIT_RADIAL: f64 = 0.0001;
const RES_ORBIT_TRACK: f64 = 0.0004;
const RES_ORBIT_DOT_RADIAL: f64 = 0.000001;
const RES_ORBIT_DOT_TRACK: f64 = 0.000004;

const RES_CLOCK_C0: f64 = 0.0001;
const RES_CLOCK_C1: f64 = 0.000001;
const RES_CLOCK_C2: f64 = 0.00000002;

const RES_CODE_BIAS: f64 = 0.01;

#[derive(Debug, Clone, PartialEq)]
pub struct SsrHeader {
    pub message_number: u16,
    pub epoch_time: u32,
    pub update_interval: u8,
    pub multiple_message_indicator: bool,
    pub satellite_reference_datum: bool,
    pub iod_ssr: u8,
    pub provider_id: u16,
    pub solution_id: u8,
    pub num_satellites: u8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SsrOrbitSat {
    pub sat_id: u8,
    pub iode: u8,
    pub delta_radial: f64,
    pub delta_along_track: f64,
    pub delta_cross_track: f64,
    pub dot_delta_radial: f64,
    pub dot_delta_along_track: f64,
    pub dot_delta_cross_track: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SsrClockSat {
    pub sat_id: u8,
    pub delta_clock_c0: f64,
    pub delta_clock_c1: f64,
    pub delta_clock_c2: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SsrCodeBiasSat {
    pub sat_id: u8,
    pub num_biases: u8,
    pub biases: Vec<SsrSignalBias>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SsrSignalBias {
    pub signal_and_tracking_mode: u8,
    pub bias: f64,
}

pub fn parse_ssr_header(
    bits: &BitSlice<u8, Msb0>,
) -> Result<(&BitSlice<u8, Msb0>, SsrHeader), RtcmParseError> {
    if bits.len() < 68 {
        return Err(RtcmParseError::Incomplete);
    }

    let header = SsrHeader {
        message_number: bits[0..12].load_be::<u16>(),
        epoch_time: bits[12..32].load_be::<u32>(),
        update_interval: bits[32..36].load_be::<u8>(),
        multiple_message_indicator: bits[36],
        satellite_reference_datum: bits[37],
        iod_ssr: bits[38..42].load_be::<u8>(),
        provider_id: bits[42..58].load_be::<u16>(),
        solution_id: bits[58..62].load_be::<u8>(),
        num_satellites: bits[62..68].load_be::<u8>(),
    };

    Ok((&bits[68..], header))
}

pub fn parse_ssr_orbit(payload: &[u8]) -> Result<(SsrHeader, Vec<SsrOrbitSat>), RtcmParseError> {
    let bits = payload.view_bits::<Msb0>();
    let (mut bits, header) = parse_ssr_header(bits)?;

    let mut sats = Vec::with_capacity(header.num_satellites as usize);
    for _ in 0..header.num_satellites {
        let (next_bits, sat) = parse_orbit_sat(bits)?;
        sats.push(sat);
        bits = next_bits;
    }

    Ok((header, sats))
}

fn parse_orbit_sat(
    bits: &BitSlice<u8, Msb0>,
) -> Result<(&BitSlice<u8, Msb0>, SsrOrbitSat), RtcmParseError> {
    if bits.len() < 135 {
        return Err(RtcmParseError::Incomplete);
    }

    let sat = SsrOrbitSat {
        sat_id: bits[0..6].load_be::<u8>(),
        iode: bits[6..14].load_be::<u8>(),
        delta_radial: sign_extend_i32(bits[14..36].load_be::<u32>(), 22) as f64 * RES_ORBIT_RADIAL,
        delta_along_track: sign_extend_i32(bits[36..56].load_be::<u32>(), 20) as f64
            * RES_ORBIT_TRACK,
        delta_cross_track: sign_extend_i32(bits[56..76].load_be::<u32>(), 20) as f64
            * RES_ORBIT_TRACK,
        dot_delta_radial: sign_extend_i32(bits[76..97].load_be::<u32>(), 21) as f64
            * RES_ORBIT_DOT_RADIAL,
        dot_delta_along_track: sign_extend_i32(bits[97..116].load_be::<u32>(), 19) as f64
            * RES_ORBIT_DOT_TRACK,
        dot_delta_cross_track: sign_extend_i32(bits[116..135].load_be::<u32>(), 19) as f64
            * RES_ORBIT_DOT_TRACK,
    };

    Ok((&bits[135..], sat))
}

pub fn parse_ssr_clock(payload: &[u8]) -> Result<(SsrHeader, Vec<SsrClockSat>), RtcmParseError> {
    let bits = payload.view_bits::<Msb0>();
    let (mut bits, header) = parse_ssr_header(bits)?;

    let mut sats = Vec::with_capacity(header.num_satellites as usize);
    for _ in 0..header.num_satellites {
        let (next_bits, sat) = parse_clock_sat(bits)?;
        sats.push(sat);
        bits = next_bits;
    }

    Ok((header, sats))
}

fn parse_clock_sat(
    bits: &BitSlice<u8, Msb0>,
) -> Result<(&BitSlice<u8, Msb0>, SsrClockSat), RtcmParseError> {
    if bits.len() < 76 {
        return Err(RtcmParseError::Incomplete);
    }

    let sat = SsrClockSat {
        sat_id: bits[0..6].load_be::<u8>(),
        delta_clock_c0: sign_extend_i32(bits[6..28].load_be::<u32>(), 22) as f64 * RES_CLOCK_C0,
        delta_clock_c1: sign_extend_i32(bits[28..49].load_be::<u32>(), 21) as f64 * RES_CLOCK_C1,
        delta_clock_c2: sign_extend_i32(bits[49..76].load_be::<u32>(), 27) as f64 * RES_CLOCK_C2,
    };

    Ok((&bits[76..], sat))
}

pub fn parse_ssr_code_bias(
    payload: &[u8],
) -> Result<(SsrHeader, Vec<SsrCodeBiasSat>), RtcmParseError> {
    let bits = payload.view_bits::<Msb0>();
    let (mut bits, header) = parse_ssr_header(bits)?;

    let mut sats = Vec::with_capacity(header.num_satellites as usize);
    for _ in 0..header.num_satellites {
        let (next_bits, sat) = parse_code_bias_sat(bits)?;
        sats.push(sat);
        bits = next_bits;
    }

    Ok((header, sats))
}

fn parse_code_bias_sat(
    bits: &BitSlice<u8, Msb0>,
) -> Result<(&BitSlice<u8, Msb0>, SsrCodeBiasSat), RtcmParseError> {
    if bits.len() < 11 {
        return Err(RtcmParseError::Incomplete);
    }

    let sat_id = bits[0..6].load_be::<u8>();
    let num_biases = bits[6..11].load_be::<u8>();
    let mut current_bits = &bits[11..];

    let mut biases = Vec::with_capacity(num_biases as usize);
    for _ in 0..num_biases {
        if current_bits.len() < 19 {
            return Err(RtcmParseError::Incomplete);
        }
        biases.push(SsrSignalBias {
            signal_and_tracking_mode: current_bits[0..5].load_be::<u8>(),
            bias: sign_extend_i32(current_bits[5..19].load_be::<u32>(), 14) as f64 * RES_CODE_BIAS,
        });
        current_bits = &current_bits[19..];
    }

    Ok((
        current_bits,
        SsrCodeBiasSat {
            sat_id,
            num_biases,
            biases,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pack (num_bits, value) pairs into a byte buffer in Msb0 order.
    fn pack_bits(pairs: &[(usize, u64)]) -> Vec<u8> {
        let total_bits: usize = pairs.iter().map(|(b, _)| *b).sum();
        let mut bytes = vec![0u8; (total_bits + 7) / 8];
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
}
