//! GPS LNAV (Legacy Navigation) subframe decoder.
//!
//! Decodes GPS L1 C/A broadcast navigation messages from raw subframe words
//! (as delivered by UBX-RXM-SFRBX) into `gneiss_core::ephemeris::Ephemeris`.
//! Reference: IS-GPS-200N, Section 20.3.3.

use gneiss_core::ephemeris::{Ephemeris, GpsEphemeris};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;

use std::f64::consts::PI as GPS_PI;

/// Scale factor for signed 2's complement values from GPS LNAV.
fn sign_scale(raw: u32, bits: u32, scale: f64) -> f64 {
    let mask = (1u64 << bits) - 1;
    let val = (raw as u64) & mask;
    let half = 1u64 << (bits - 1);
    let signed = if val >= half {
        // Sign-extend: val is between half and mask, represents negative number
        (val as i64).wrapping_sub(1i64 << bits)
    } else {
        val as i64
    };
    signed as f64 * scale
}

/// Extract a bit field from the 10 subframe words.
/// Words are 30-bit values stored in the lower 30 bits of each u32.
fn extract_bits(words: &[u32], start_bit: usize, num_bits: usize) -> u32 {
    let mut result: u32 = 0;
    for i in 0..num_bits {
        let bit_idx = start_bit + i;
        let word_idx = bit_idx / 30;
        let bit_in_word = bit_idx % 30;
        if word_idx < words.len() && (words[word_idx] >> bit_in_word) & 1 == 1 {
            result |= 1 << i;
        }
    }
    result
}

/// Parse HOW (Handover Word) from word 2 of any subframe.
/// Returns (tow, subframe_id).
fn parse_how(word2: u32) -> (u32, u8) {
    // HOW format (IS-GPS-200N, Section 20.3.3.1):
    // Bits 0-16: truncated TOW count (17 bits)
    // Bits 17-18: flag bits
    // Bits 19-21: subframe ID (3 bits)
    // Bits 22-27: parity (6 bits)
    let tow = word2 & 0x1FFFF; // 17 bits
    let sf_id = ((word2 >> 19) & 0x7) as u8;
    (tow, sf_id)
}

/// Decode GPS LNAV subframe 1 (clock + health).
fn decode_sf1(words: &[u32], eph: &mut GpsEphemeris, week: u32) {
    // IS-GPS-200N Table 20-I
    eph.iodc = (extract_bits(words, 210, 2) << 8) | extract_bits(words, 60, 8);
    // WN is 10 bits at word 3 bit 6 (bit offset 60+6=66?)
    // Actually: WN is bits 60-69 (word 3 bits 6-15)
    // L2 codes, etc.
    let toc_raw = extract_bits(words, 210 + 2 + 2 + 10, 16);
    eph.toc = GpsTime::new(week, toc_raw as f64 * 16.0);
    eph.af2 = sign_scale(extract_bits(words, 240, 8), 8, 2.0_f64.powi(-55));
    eph.af1 = sign_scale(extract_bits(words, 248, 16), 16, 2.0_f64.powi(-43));
    eph.af0 = sign_scale(extract_bits(words, 264, 22), 22, 2.0_f64.powi(-31));
    eph.tgd = sign_scale(extract_bits(words, 192, 8), 8, 2.0_f64.powi(-31));
}

/// Decode GPS LNAV subframe 2 (ephemeris part 1).
fn decode_sf2(words: &[u32], eph: &mut GpsEphemeris) {
    // IS-GPS-200N Table 20-II
    eph.iode = extract_bits(words, 60, 8);
    eph.crs = sign_scale(extract_bits(words, 68, 16), 16, 2.0_f64.powi(-5));
    eph.delta_n = sign_scale(extract_bits(words, 84, 16), 16, 2.0_f64.powi(-43)) * GPS_PI;
    eph.m0 = sign_scale(extract_bits(words, 100, 32), 32, 2.0_f64.powi(-31)) * GPS_PI;
    eph.cuc = sign_scale(extract_bits(words, 150, 16), 16, 2.0_f64.powi(-29));
    eph.e = (extract_bits(words, 166, 32) as f64) * 2.0_f64.powi(-33);
    eph.cus = sign_scale(extract_bits(words, 210, 16), 16, 2.0_f64.powi(-29));
    eph.sqrt_a = (extract_bits(words, 226, 32) as f64) * 2.0_f64.powi(-19);
    // toe is bits 270-285 (16 bits) in word 10
    eph.toe = GpsTime::new(0, extract_bits(words, 270, 16) as f64 * 16.0);
}

/// Decode GPS LNAV subframe 3 (ephemeris part 2).
fn decode_sf3(words: &[u32], eph: &mut GpsEphemeris) {
    // IS-GPS-200N Table 20-III
    eph.cic = sign_scale(extract_bits(words, 60, 16), 16, 2.0_f64.powi(-29));
    eph.omega0 = sign_scale(extract_bits(words, 76, 32), 32, 2.0_f64.powi(-31)) * GPS_PI;
    eph.cis = sign_scale(extract_bits(words, 108, 16), 16, 2.0_f64.powi(-29));
    eph.i0 = sign_scale(extract_bits(words, 124, 32), 32, 2.0_f64.powi(-31)) * GPS_PI;
    eph.crc = sign_scale(extract_bits(words, 168, 16), 16, 2.0_f64.powi(-5));
    eph.omega = sign_scale(extract_bits(words, 196, 32), 32, 2.0_f64.powi(-31)) * GPS_PI;
    eph.omega_dot = sign_scale(extract_bits(words, 240, 24), 24, 2.0_f64.powi(-43)) * GPS_PI;
    eph.idot = sign_scale(extract_bits(words, 276, 14), 14, 2.0_f64.powi(-43)) * GPS_PI;
}

/// Build GPS ephemerides from a set of UBX-RXM-SFRBX messages.
/// `sfrbx_messages` should be all SFRBX messages for GPS (gnss_id=0).
/// Returns a vector of Ephemeris.
pub fn build_ephemeris_from_sfrbx(
    sfrbx_messages: &[crate::ubx::UbxRxmSfrbx],
) -> Vec<Ephemeris> {
    use std::collections::HashMap;

    // Group words by (sv_id, subframe_id), keeping the latest
    // For GPS LNAV, we need subframes 1, 2, 3 for each satellite
    let mut sf_data: HashMap<(u8, u8), Vec<u32>> = HashMap::new();

    for msg in sfrbx_messages {
        if msg.gnss_id != 0 { continue; } // GPS only
        if msg.words.len() < 10 { continue; } // Need complete subframe

        let (_tow, sf_id) = parse_how(msg.words[1]);
        let key = (msg.sv_id, sf_id);

        // Only keep the latest version of each subframe
        sf_data.insert(key, msg.words.clone());
    }

    // Build ephemeris for each satellite that has all 3 subframes
    let mut ephemerides = Vec::new();
    for sv_id in 1..=32u8 {
        let w1 = match sf_data.get(&(sv_id, 1)) { Some(w) => w.clone(), None => continue };
        let w2 = match sf_data.get(&(sv_id, 2)) { Some(w) => w.clone(), None => continue };
        let w3 = match sf_data.get(&(sv_id, 3)) { Some(w) => w.clone(), None => continue };

        let mut eph = GpsEphemeris {
            sat: SatelliteId { constellation: Constellation::Gps, prn: sv_id },
            toe: GpsTime::new(0, 0.0),
            toc: GpsTime::new(0, 0.0),
            af0: 0.0, af1: 0.0, af2: 0.0,
            crs: 0.0, crc: 0.0, cuc: 0.0, cus: 0.0, cic: 0.0, cis: 0.0,
            m0: 0.0, e: 0.0, sqrt_a: 0.0,
            delta_n: 0.0,
            omega0: 0.0, omega_dot: 0.0, i0: 0.0, idot: 0.0, omega: 0.0,
            tgd: 0.0,
            iode: 0, iodc: 0,
        };

        // Extract week from SF1 (bits 60-69 = word 3 bits 6-15)
        let week = extract_bits(&w1, 60, 10);
        eph.toe = GpsTime::new(week, 0.0);
        eph.toc = GpsTime::new(week, 0.0);

        decode_sf1(&w1, &mut eph, week);
        decode_sf2(&w2, &mut eph);
        decode_sf3(&w3, &mut eph);

        // Fix toe/toc GPS week (SF1 provides it)
        let toe_week = week;
        eph.toe = GpsTime::new(toe_week, eph.toe.tow);
        eph.toc = GpsTime::new(week, eph.toc.tow);

        ephemerides.push(Ephemeris::Gps(eph));
    }

    ephemerides
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_how() {
        // Example HOW: tow=100, sf_id=1
        let how = 100u32 | (1u32 << 19);
        let (tow, sf_id) = parse_how(how);
        assert_eq!(tow, 100);
        assert_eq!(sf_id, 1);
    }

    #[test]
    fn test_extract_bits_simple() {
        let words = vec![0b1011u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32];
        assert_eq!(extract_bits(&words, 0, 4), 0b1011);
        assert_eq!(extract_bits(&words, 1, 2), 0b01);
    }

    #[test]
    fn test_sign_scale_positive() {
        // 3-bit, value=3, scale=1.0 -> 3.0
        assert_eq!(sign_scale(3, 3, 1.0), 3.0);
    }

    #[test]
    fn test_sign_scale_negative() {
        // 3-bit, value=4 (100b) should be -4 in 2's complement
        assert_eq!(sign_scale(4, 3, 1.0), -4.0);
    }
}
