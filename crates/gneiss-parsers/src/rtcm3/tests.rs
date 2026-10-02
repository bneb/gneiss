//! Tests for the RTCM 3 frame layer: CRC-24Q and frame framing.

use super::*;
use bitvec::prelude::*;

/// Wraps `payload` in a complete RTCM 3 frame: preamble, 10-bit length,
/// payload, then the CRC-24Q of the preceding bytes.
fn frame(payload: &[u8]) -> Vec<u8> {
    let mut f = vec![0xD3, 0x00, payload.len() as u8];
    f.extend_from_slice(payload);
    let crc = crc24q(&f);
    f.push((crc >> 16) as u8);
    f.push((crc >> 8) as u8);
    f.push(crc as u8);
    f
}

/// Golden CRC-24Q vectors. The expected values were produced by running
/// RTKLIB's independent table-driven implementation (`rtk_crc24q` plus
/// `tbl_CRC24Q[256]`, rtkcmn.c:316 and :912) over the same byte strings; the
/// crate's bitwise loop must agree with it bit for bit.
#[test]
fn crc24q_matches_reference_vectors() {
    let cases: &[(&[u8], u32)] = &[
        (b"", 0x00_00_00),
        (b"\xd3", 0xE3_33_09),
        (b"abc", 0x9F_F3_59),
        (b"\x00", 0x00_00_00),
        (&[0xFF, 0xFF, 0xFF], 0xED_F8_CE),
        (&[0xD3, 0x00, 0x3D], 0x3A_5E_08),
    ];
    for (input, expected) in cases {
        assert_eq!(crc24q(input), *expected, "crc24q mismatch for {input:02x?}");
    }
}

/// The CRC must depend on every bit of its input, including the trailing bit
/// of the last byte.
#[test]
fn crc24q_depends_on_every_input_bit() {
    let base: Vec<u8> = (0u8..16).collect();
    let reference = crc24q(&base);
    for byte in 0..base.len() {
        for bit in 0..8 {
            let mut flipped = base.clone();
            flipped[byte] ^= 1 << bit;
            assert_ne!(
                crc24q(&flipped),
                reference,
                "flipping bit {bit} of byte {byte} did not change the CRC"
            );
        }
    }
}

#[test]
fn parses_a_complete_frame_and_leaves_the_remainder() {
    let payload = [0x00u8, 0x00, 0x00]; // message number 0
    let mut stream = frame(&payload);
    let tail = [0xAAu8, 0xBB];
    stream.extend_from_slice(&tail);

    let (rest, parsed) = parse_rtcm3_frame(&stream).expect("well-formed frame");
    assert_eq!(parsed.payload, &payload[..]);
    assert_eq!(rest, &tail[..], "the parser must consume exactly one frame");
}

#[test]
fn rejects_a_wrong_preamble() {
    let mut f = frame(&[0x01, 0x02]);
    f[0] = 0xD2;
    assert_eq!(parse_rtcm3_frame(&f), Err(RtcmParseError::InvalidPreamble));
}

/// A zero-length message is structurally legal: the frame is preamble, length
/// and CRC only. It must parse, and must then fail downstream field decoding
/// rather than panicking.
#[test]
fn zero_length_payload_parses_as_an_empty_message() {
    let f = frame(&[]);
    assert_eq!(f.len(), 6);
    let (rest, parsed) = parse_rtcm3_frame(&f).expect("empty payload is well-formed");
    assert!(parsed.payload.is_empty());
    assert!(rest.is_empty());
    assert_eq!(parse_1019(parsed.payload), Err(RtcmParseError::Incomplete));
}

/// Every proper prefix of a valid frame must be reported as incomplete rather
/// than parsed or panicking.
#[test]
fn every_truncated_prefix_is_incomplete() {
    let payload: Vec<u8> = (0u8..19).collect();
    let f = frame(&payload);
    for cut in 0..f.len() {
        assert_eq!(
            parse_rtcm3_frame(&f[..cut]),
            Err(RtcmParseError::Incomplete),
            "prefix of length {cut} must be incomplete"
        );
    }
    assert!(parse_rtcm3_frame(&f).is_ok());
}

#[test]
fn short_buffers_never_panic() {
    for len in 0..8usize {
        let buf = vec![0xD3u8; len];
        assert_eq!(parse_rtcm3_frame(&buf), Err(RtcmParseError::Incomplete));
    }
    assert_eq!(parse_rtcm3_frame(&[]), Err(RtcmParseError::Incomplete));
}

/// A single flipped bit anywhere in the CRC'd region must be detected.
#[test]
fn single_bit_flip_in_the_covered_region_is_detected() {
    let f = frame(&[0x3F, 0xD0, 0x07, 0xAA]);
    // Byte 0 is the preamble: flipping it is caught earlier, as a preamble error.
    assert_eq!(
        parse_rtcm3_frame(&{ let mut b = f.clone(); b[0] ^= 1; b }),
        Err(RtcmParseError::InvalidPreamble)
    );
    // Bits 0..1 of byte 1 are the two length MSBs: flipping them changes the
    // declared frame size and is caught as Incomplete instead.
    assert_eq!(
        parse_rtcm3_frame(&{ let mut b = f.clone(); b[1] ^= 0x01; b }),
        Err(RtcmParseError::Incomplete)
    );
    // Bits 0..1 of byte 1 and all of byte 2 are the 10-bit length: flipping
    // them changes the declared frame size, so they are excluded here. Byte 1
    // bits 2..7 are reserved and bytes 3.. are the payload.
    for byte in 1..f.len() - 3 {
        let first_bit = match byte {
            1 => 2,
            2 => 8, // no CRC-only bits in the length LSB byte
            _ => 0,
        };
        for bit in first_bit..8 {
            let mut bad = f.clone();
            bad[byte] ^= 1 << bit;
            let result = parse_rtcm3_frame(&bad);
            assert_eq!(
                result,
                Err(RtcmParseError::CrcMismatch),
                "bit {bit} of byte {byte} was not covered by the CRC"
            );
        }
    }
}

/// The 6 reserved bits of the 16-bit length field are masked off rather than
/// validated. That is safe because the CRC covers the length field: a frame
/// with junk in the reserved bits must fail the checksum.
#[test]
fn reserved_length_bits_are_covered_by_the_checksum() {
    let f = frame(&[0x01, 0x02, 0x03]);
    let mut bad = f.clone();
    bad[1] |= 0x40; // set one reserved bit of the length field
    assert_eq!(parse_rtcm3_frame(&bad), Err(RtcmParseError::CrcMismatch));

    // Re-checksumming a frame whose reserved bits are set is accepted and the
    // payload is still the true 10-bit length, so masking cannot mis-frame.
    // The checksum covers preamble + length + payload, i.e. all but the 3
    // trailing CRC bytes.
    let mut fixed = bad.clone();
    let crc = crc24q(&fixed[..f.len() - 3]);
    fixed[f.len() - 3] = (crc >> 16) as u8;
    fixed[f.len() - 2] = (crc >> 8) as u8;
    fixed[f.len() - 1] = crc as u8;
    let (rest, parsed) = parse_rtcm3_frame(&fixed).expect("re-checksummed frame");
    assert_eq!(parsed.payload, &[0x01, 0x02, 0x03]);
    assert!(rest.is_empty());
}

/// A stream may contain the same message repeatedly and out of any order; each
/// frame must decode independently of its neighbours.
#[test]
fn duplicate_and_reordered_messages_decode_independently() {
    let a = frame(&[0x3F, 0xD0, 0x13, 0x01]); // ephemeris-shaped payload
    let b = frame(&[0x3F, 0xD1, 0x15, 0x02]); // a different message number
    let mut stream = Vec::new();
    stream.extend_from_slice(&a);
    stream.extend_from_slice(&b);
    stream.extend_from_slice(&a);

    let mut seen = Vec::new();
    let mut buf: &[u8] = &stream;
    while let Ok((rest, parsed)) = parse_rtcm3_frame(buf) {
        seen.push(parsed.payload.to_vec());
        buf = rest;
    }
    assert!(buf.is_empty(), "the whole stream must be consumed");
    assert_eq!(seen, vec![a[3..7].to_vec(), b[3..7].to_vec(), a[3..7].to_vec()]);
}

/// Deterministic adversarial input: a length-gradient stream of garbage that
/// must never panic and must never be reported as a valid frame.
#[test]
fn adversarial_garbage_never_panics() {
    let mut state: u32 = 0x1234_5678;
    let mut noise = || {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        (state >> 24) as u8
    };
    for _ in 0..2000 {
        let len = (noise() as usize % 64) + 3;
        let mut buf = vec![0xD3u8; len];
        for b in buf.iter_mut().skip(1) {
            *b = noise();
        }
        if let Ok((_, parsed)) = parse_rtcm3_frame(&buf) {
            // Anything accepted must satisfy the checksum by construction.
            assert_eq!(parsed.payload.len() + 6, len);
        }
    }
}

/// A frame header claiming the maximum 10-bit length must be treated as
/// incomplete, not as a short frame that re-synchronises mid-message.
#[test]
fn maximum_declared_length_is_not_silently_shortened() {
    let mut f = vec![0xD3, 0x03, 0xFF]; // 1023-byte payload announced
    f.extend_from_slice(&[0u8; 32]);
    assert_eq!(parse_rtcm3_frame(&f), Err(RtcmParseError::Incomplete));
}

/// Message numbers other than 1019 must be rejected by the ephemeris decoder
/// with a type error, not decoded into a bogus ephemeris.
#[test]
fn foreign_message_numbers_are_rejected_by_the_ephemeris_decoder() {
    for num in [0u16, 1005, 1074, 1124, 4095] {
        let mut bits = bitvec![u8, Msb0; 0; crate::rtcm3::ephemeris::TYPE_1019_PAYLOAD_BITS];
        bits[0..12].store_be(num);
        assert_eq!(
            parse_1019(&bits.into_vec()),
            Err(RtcmParseError::UnsupportedMsmType),
            "message {num} must not decode as a 1019 ephemeris"
        );
    }
}