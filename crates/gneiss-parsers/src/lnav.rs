//! GPS LNAV (Legacy Navigation) subframe decoder.
//!
//! Decodes GPS L1 C/A broadcast navigation messages from raw subframe words
//! (as delivered by UBX-RXM-SFRBX) into `gneiss_core::ephemeris::Ephemeris`.
//! Reference: IS-GPS-200N, Section 20.3.3.

use gneiss_core::ephemeris::{Ephemeris, GpsEphemeris};
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;

use std::f64::consts::PI as GPS_PI;

/// Reference week used to resolve the 10-bit subframe week field against.
///
/// GPS week 2048 began 2019-04-06, the first LNAV rollover after the 2019
/// handover. Callers that know the observation time should resolve the week
/// themselves; this constant only removes the 1024-week ambiguity that would
/// otherwise silently place a 2020 ephemeris in 1981.
pub const WEEK_REFERENCE: u32 = 2048;

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

/// Number of information bits carried by one 30-bit LNAV word.
///
/// The remaining 6 bits are non-information (parity) bits. Proven by the
/// reference decoder shipped with this repository, `rtkcmn.c:950` `decode_word`:
/// it masks the information field with `0x3FFFFFC0` (word bits 6..29), shifts the
/// Hamming masks right by 6 before folding, and compares the folded result
/// against `word & 0x3F` for the parity.
pub const INFO_BITS_PER_WORD: usize = 24;

/// Extract a field from the concatenated stream of *information* bits.
///
/// `words` holds one 30-bit LNAV word per entry (two parity bits plus 24 data
/// bits plus four reserved bits is the IS-GPS-200 shape; UBX-RXM-SFRBX delivers
/// the raw 30-bit word). `start` counts information bits from the first
/// transmitted bit of the first word, so a field can never alias or straddle a
/// non-information bit.
///
/// Information bit `j` of a word sits at word bit `29 - (j % 24)`. The first
/// transmitted bit of a field is that field's *most significant* bit, so bit
/// `num_bits - 1 - k` of the result is bit `start + k` of the stream.
fn extract_info_bits(words: &[u32], start: usize, num_bits: usize) -> u32 {
    let mut result: u32 = 0;
    for i in 0..num_bits {
        let idx = start + i;
        let (word_idx, within) = (idx / INFO_BITS_PER_WORD, idx % INFO_BITS_PER_WORD);
        if word_idx < words.len() && (words[word_idx] >> (29 - within)) & 1 == 1 {
            // The first transmitted bit of a field is its most significant bit.
            result |= 1 << (num_bits - 1 - i);
        }
    }
    result
}

/// Extract a bit field from the 10 subframe words, counting raw word bits.
///
/// Bits 0..5 of every word are the non-information parity bits (see
/// [`INFO_BITS_PER_WORD`]), so this accessor is only safe for fields that a
/// caller has already proven to avoid them. New code should use
/// [`extract_info_bits`], which cannot read a parity bit by construction.
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
fn parse_how(words: &[u32]) -> (u32, u8) {
    // HOW format (IS-GPS-200N, Section 20.3.3.1), in information-bit indices
    // within word 2 (word 1 occupies information bits 24..47 of the subframe):
    //   0-16: truncated TOW count (17 bits)
    //   17-18: flag bits
    //   19-21: subframe ID (3 bits)
    // Word 2's information bits start at index 24.
    let tow = extract_info_bits(words, 24, 17);
    let sf_id = extract_info_bits(words, 24 + 19, 3) as u8;
    (tow, sf_id)
}

/// The three subframe decoders below still index fields with [`extract_bits`],
/// whose offsets are expressed in raw word bits and therefore include the six
/// non-information parity bits of every word.
///
/// # Known defect, not yet repaired here
///
/// `decode_sf1`, `decode_sf2` and `decode_sf3` currently return wrong values:
/// of the 23 fields they read, 18 either start on a parity bit or run across
/// one (every 32-bit field - `m0`, `e`, `sqrt_a`, `omega0`, `i0`, `omega`,
/// `af0` - includes six parity bits). Only `crs`, `tgd`, `toc`, `af1` and
/// `idot` happen to avoid them.
///
/// Repairing the offsets needs the per-word field layout of IS-GPS-200
/// Tables 20-I, 20-II and 20-III, which is not available in this repository
/// and could not be verified here. Guessing the offsets and hard-coding the
/// guess would launder a bug into a passing test, so the offsets were left
/// alone and the hazard documented instead. [`build_ephemeris_from_sfrbx`] has
/// no production caller, so nothing downstream consumes these values yet.
///
/// `parse_how` (fixed) and `extract_info_bits` are the verified-correct
/// accessors; converting the three decoders is mechanical once the table is
/// available: express each field as an information-bit index into the 192-bit
/// stream of the eight data words and read it with [`extract_info_bits`].
///
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

        let (_tow, sf_id) = parse_how(&msg.words);
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

        // The SF1 week field is 10 bits wide, so it counts modulo 1024 weeks and
        // must be resolved before it can be used as a `GpsTime` week. Without
        // this, a 2020 broadcast week of 2110 is stored as 62 -> 1981.
        let week = crate::rtcm3::ephemeris::adjust_gps_week(extract_bits(&w1, 60, 10), WEEK_REFERENCE);
        eph.toe = GpsTime::new(week, 0.0);
        eph.toc = GpsTime::new(week, 0.0);

        decode_sf1(&w1, &mut eph, week);
        decode_sf2(&w2, &mut eph);
        decode_sf3(&w3, &mut eph);

        // Re-attach the resolved week: the subframe decoders only produce the
        // time of week.
        eph.toe = GpsTime::new(week, eph.toe.tow);
        eph.toc = GpsTime::new(week, eph.toc.tow);

        ephemerides.push(Ephemeris::Gps(eph));
    }

    ephemerides
}

#[cfg(test)]
mod tests;
