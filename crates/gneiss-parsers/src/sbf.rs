//! Septentrio Binary Format (SBF) streaming parser.
//!
//! SBF is the native binary protocol for Septentrio GNSS receivers (AsteRx, mosaic).
//! Each block begins with a sync header (`$@`), a 16-bit CRC-CCITT, a 16-bit block ID/revision,
//! and a 16-bit length (multiple of 4 bytes).

/// Common SBF Block IDs.
pub mod block_ids {
    pub const PVT_GEODETIC: u16 = 4007;
    pub const POS_COV_GEODETIC: u16 = 4008;
    pub const ATT_EULER: u16 = 5938;
    pub const ATT_COV_EULER: u16 = 5939;
    pub const MEAS_EPOCH: u16 = 4027;
    pub const RECEIVER_STATUS: u16 = 4014;
    pub const BASE_VECTOR_GEOD: u16 = 4028;
}

/// A parsed raw SBF block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SbfBlock<'a> {
    /// 13-bit block number (e.g., 4007 for PVTGeodetic).
    pub block_number: u16,
    /// 3-bit block revision number.
    pub block_revision: u8,
    /// Payload slice (excluding the 8-byte SBF header).
    pub payload: &'a [u8],
}

/// Parsing errors for SBF streams.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SbfParseError {
    /// Not enough bytes in input buffer to complete the frame.
    Incomplete,
    /// First two bytes do not match the SBF sync pattern (`$@` / `0x24 0x40`).
    InvalidSync,
    /// Frame length is less than 8 or not a multiple of 4 bytes.
    InvalidLength,
    /// CRC-CCITT checksum validation failed.
    CrcMismatch,
}

/// Precomputed CRC-CCITT lookup table (polynomial 0x1021).
const CRC16_TABLE: [u16; 256] = [
    0x0000, 0x1021, 0x2042, 0x3063, 0x4084, 0x50a5, 0x60c6, 0x70e7,
    0x8108, 0x9129, 0xa14a, 0xb16b, 0xc18c, 0xd1ad, 0xe1ce, 0xf1ef,
    0x1231, 0x0210, 0x3273, 0x2252, 0x52b5, 0x4294, 0x72f7, 0x62d6,
    0x9339, 0x8318, 0xb37b, 0xa35a, 0xd3bd, 0xc39c, 0xf3ff, 0xe3de,
    0x2462, 0x3443, 0x0420, 0x1401, 0x64e6, 0x74c7, 0x44a4, 0x5485,
    0xa56a, 0xb54b, 0x8528, 0x9509, 0xe5ee, 0xf5cf, 0xc5ac, 0xd58d,
    0x3653, 0x2672, 0x1611, 0x0630, 0x76d7, 0x66f6, 0x5695, 0x46b4,
    0xb75b, 0xa77a, 0x9719, 0x8738, 0xf7df, 0xe7fe, 0xd79d, 0xc7bc,
    0x48c4, 0x58e5, 0x6886, 0x78a7, 0x0840, 0x1861, 0x2802, 0x3823,
    0xc9cc, 0xd9ed, 0xe98e, 0xf9af, 0x8948, 0x9969, 0xa90a, 0xb92b,
    0x5af5, 0x4ad4, 0x7ab7, 0x6a96, 0x1a71, 0x0a50, 0x3a33, 0x2a12,
    0xdbfd, 0xcbdc, 0xfbbf, 0xeb9e, 0x9b79, 0x8b58, 0xbb3b, 0xab1a,
    0x6ca6, 0x7c87, 0x4ce4, 0x5cc5, 0x2c22, 0x3c03, 0x0c60, 0x1c41,
    0xedae, 0xfd8f, 0xcdec, 0xddcd, 0xad2a, 0xbd0b, 0x8d68, 0x9d49,
    0x7e97, 0x6eb6, 0x5ed5, 0x4ef4, 0x3e13, 0x2e32, 0x1e51, 0x0e70,
    0xff9f, 0xefbe, 0xdfdd, 0xcffc, 0xbf1b, 0xaf3a, 0x9f59, 0x8f78,
    0x9188, 0x81a9, 0xb1ca, 0xa1eb, 0xd10c, 0xc12d, 0xf14e, 0xe16f,
    0x1080, 0x00a1, 0x30c2, 0x20e3, 0x5004, 0x4025, 0x7046, 0x6067,
    0x83b9, 0x9398, 0xa3fb, 0xb3da, 0xc33d, 0xd31c, 0xe37f, 0xf35e,
    0x02b1, 0x1290, 0x22f3, 0x32d2, 0x4235, 0x5214, 0x6277, 0x7256,
    0xb5ea, 0xa5cb, 0x95a8, 0x8589, 0xf56e, 0xe54f, 0xd52c, 0xc50d,
    0x34e2, 0x24c3, 0x14a0, 0x0481, 0x7466, 0x6447, 0x5424, 0x4405,
    0xa7db, 0xb7fa, 0x8799, 0x97b8, 0xe75f, 0xf77e, 0xc71d, 0xd73c,
    0x26d3, 0x36f2, 0x0691, 0x16b0, 0x6657, 0x7676, 0x4615, 0x5634,
    0xd94c, 0xc96d, 0xf90e, 0xe92f, 0x99c8, 0x89e9, 0xb98a, 0xa9ab,
    0x5844, 0x4865, 0x7806, 0x6827, 0x18c0, 0x08e1, 0x3882, 0x28a3,
    0xcb7d, 0xdb5c, 0xeb3f, 0xfb1e, 0x8bf9, 0x9bd8, 0xabbb, 0xbb9a,
    0x4a75, 0x5a54, 0x6a37, 0x7a16, 0x0af1, 0x1ad0, 0x2ab3, 0x3a92,
    0xfd2e, 0xed0f, 0xdd6c, 0xcd4d, 0xbdaa, 0xad8b, 0x9de8, 0x8dc9,
    0x7c26, 0x6c07, 0x5c64, 0x4c45, 0x3ca2, 0x2c83, 0x1ce0, 0x0cc1,
    0xef1f, 0xff3e, 0xcf5d, 0xdf7c, 0xaf9b, 0xbfba, 0x8fd9, 0x9ff8,
    0x6e17, 0x7e36, 0x4e55, 0x5e74, 0x2e93, 0x3eb2, 0x0ed1, 0x1ef0,
];

/// Computes the Septentrio SBF 16-bit CRC-CCITT over a byte buffer.
#[inline]
pub fn sbf_crc16(data: &[u8]) -> u16 {
    let mut crc: u16 = 0;
    for &byte in data {
        crc = (crc << 8) ^ CRC16_TABLE[((crc >> 8) ^ (byte as u16)) as usize];
    }
    crc
}

/// Scan a buffer to locate the next `$@` sync pattern index.
#[inline]
pub fn find_sbf_sync(input: &[u8]) -> Option<usize> {
    if input.len() < 2 {
        return None;
    }
    input.windows(2).position(|w| w[0] == 0x24 && w[1] == 0x40)
}

/// Parses a single SBF block from an input byte slice.
///
/// Returns the unconsumed remaining bytes and the decoded `SbfBlock`.
pub fn parse_sbf_block(input: &[u8]) -> Result<(&[u8], SbfBlock<'_>), SbfParseError> {
    if input.len() < 8 {
        return Err(SbfParseError::Incomplete);
    }
    if input[0] != 0x24 || input[1] != 0x40 {
        return Err(SbfParseError::InvalidSync);
    }

    let expected_crc = u16::from_le_bytes([input[2], input[3]]);
    let block_id_raw = u16::from_le_bytes([input[4], input[5]]);
    let block_number = block_id_raw & 0x1FFF;
    let block_revision = ((block_id_raw >> 13) & 0x07) as u8;
    let length = u16::from_le_bytes([input[6], input[7]]) as usize;

    if length < 8 || !length.is_multiple_of(4) {
        return Err(SbfParseError::InvalidLength);
    }
    if input.len() < length {
        return Err(SbfParseError::Incomplete);
    }

    let calculated_crc = sbf_crc16(&input[4..length]);
    if calculated_crc != expected_crc {
        return Err(SbfParseError::CrcMismatch);
    }

    let payload = &input[8..length];
    let remaining = &input[length..];
    Ok((remaining, SbfBlock { block_number, block_revision, payload }))
}

/// Stream scanner that parses all valid SBF blocks from a potentially noisy buffer,
/// recovering automatically from sync slips.
pub fn parse_sbf_stream_resilient(mut input: &[u8]) -> (Vec<SbfBlock<'_>>, usize) {
    let mut blocks = Vec::new();
    let mut valid_bytes = 0;

    while !input.is_empty() {
        match parse_sbf_block(input) {
            Ok((rem, block)) => {
                let consumed = input.len() - rem.len();
                valid_bytes += consumed;
                blocks.push(block);
                input = rem;
            }
            Err(SbfParseError::InvalidSync) | Err(SbfParseError::CrcMismatch) | Err(SbfParseError::InvalidLength) => {
                // Advance to next sync
                if let Some(pos) = find_sbf_sync(&input[1..]) {
                    input = &input[1 + pos..];
                } else {
                    break;
                }
            }
            Err(SbfParseError::Incomplete) => {
                break;
            }
        }
    }

    (blocks, valid_bytes)
}

/// Helper function to construct a valid SBF byte buffer for testing or synthetic generation.
pub fn build_sbf_frame(block_number: u16, block_revision: u8, payload: &[u8]) -> Vec<u8> {
    let length = 8 + payload.len();
    assert!(length.is_multiple_of(4), "SBF block total length must be a multiple of 4 bytes");

    let block_id = (block_number & 0x1FFF) | (((block_revision & 0x07) as u16) << 13);
    let mut header_and_payload = Vec::with_capacity(length - 4);
    header_and_payload.extend_from_slice(&block_id.to_le_bytes());
    header_and_payload.extend_from_slice(&(length as u16).to_le_bytes());
    header_and_payload.extend_from_slice(payload);

    let crc = sbf_crc16(&header_and_payload);

    let mut frame = Vec::with_capacity(length);
    frame.push(0x24); // '$'
    frame.push(0x40); // '@'
    frame.extend_from_slice(&crc.to_le_bytes());
    frame.extend_from_slice(&header_and_payload);
    frame
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sbf_crc_known_vector() {
        let test_data = [0x07, 0x0F, 0x10, 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08];
        let crc = sbf_crc16(&test_data);
        assert_ne!(crc, 0);
    }

    #[test]
    fn test_sbf_parse_valid_block() {
        let payload = [0x11, 0x22, 0x33, 0x44];
        let frame = build_sbf_frame(block_ids::PVT_GEODETIC, 1, &payload);

        let (rem, block) = parse_sbf_block(&frame).expect("parse should succeed");
        assert!(rem.is_empty());
        assert_eq!(block.block_number, block_ids::PVT_GEODETIC);
        assert_eq!(block.block_revision, 1);
        assert_eq!(block.payload, &payload);
    }

    #[test]
    fn test_sbf_parse_incomplete() {
        let payload = [0x11, 0x22, 0x33, 0x44];
        let frame = build_sbf_frame(block_ids::MEAS_EPOCH, 0, &payload);

        let res = parse_sbf_block(&frame[..6]);
        assert_eq!(res, Err(SbfParseError::Incomplete));

        let res = parse_sbf_block(&frame[..10]);
        assert_eq!(res, Err(SbfParseError::Incomplete));
    }

    #[test]
    fn test_sbf_parse_invalid_sync() {
        let mut frame = build_sbf_frame(block_ids::PVT_GEODETIC, 0, &[0; 4]);
        frame[0] = 0xFF;

        let res = parse_sbf_block(&frame);
        assert_eq!(res, Err(SbfParseError::InvalidSync));
    }

    #[test]
    fn test_sbf_parse_crc_mismatch() {
        let mut frame = build_sbf_frame(block_ids::ATT_EULER, 2, &[1, 2, 3, 4]);
        frame[9] ^= 0xFF;

        let res = parse_sbf_block(&frame);
        assert_eq!(res, Err(SbfParseError::CrcMismatch));
    }

    #[test]
    fn test_sbf_parse_invalid_length() {
        let mut frame = build_sbf_frame(block_ids::PVT_GEODETIC, 0, &[0; 4]);
        frame[6] = 9;
        frame[7] = 0;

        let res = parse_sbf_block(&frame);
        assert_eq!(res, Err(SbfParseError::InvalidLength));
    }

    #[test]
    fn test_sbf_stream_sequential_blocks() {
        let p1 = [0xAA; 4];
        let p2 = [0xBB; 8];
        let f1 = build_sbf_frame(block_ids::PVT_GEODETIC, 1, &p1);
        let f2 = build_sbf_frame(block_ids::ATT_EULER, 2, &p2);

        let mut stream = Vec::new();
        stream.extend_from_slice(&f1);
        stream.extend_from_slice(&f2);

        let (rem1, b1) = parse_sbf_block(&stream).expect("b1");
        assert_eq!(b1.block_number, block_ids::PVT_GEODETIC);
        assert_eq!(b1.payload, &p1);

        let (rem2, b2) = parse_sbf_block(rem1).expect("b2");
        assert_eq!(b2.block_number, block_ids::ATT_EULER);
        assert_eq!(b2.payload, &p2);
        assert!(rem2.is_empty());
    }

    #[test]
    fn test_sbf_stream_resilient_sync_recovery() {
        let p1 = [0xAA; 4];
        let f1 = build_sbf_frame(block_ids::PVT_GEODETIC, 1, &p1);

        // Prepend garbage bytes before valid frame
        let mut noisy_stream = vec![0xFF, 0xDE, 0xAD, 0xBE, 0xEF];
        noisy_stream.extend_from_slice(&f1);
        noisy_stream.extend_from_slice(&[0x12, 0x34]); // trailing incomplete

        let (blocks, valid_bytes) = parse_sbf_stream_resilient(&noisy_stream);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].block_number, block_ids::PVT_GEODETIC);
        assert_eq!(valid_bytes, f1.len());
    }
}
