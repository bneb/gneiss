//! Tests for the GPS LNAV subframe decoder.
//!
//! The authority for word layout is the reference decoder shipped with this
//! repository, `rtkcmn.c:950` `decode_word`, which proves that of each 30-bit
//! LNAV word only bits 6..29 carry information (bits 0..5 are the parity bits
//! compared against `word & 0x3F`), and that the first transmitted information
//! bit is the *top* of the word (`data[0] = word >> 22`).

use super::*;

/// A subframe is 10 words; only the information bits are used in these vectors.
fn blank_subframe() -> Vec<u32> {
    vec![0u32; 10]
}

/// Writes `value` into the `num_bits`-wide information field starting at
/// information index `start` of word 1 (the handover word).
fn set_how_field(words: &mut [u32], start: usize, num_bits: usize, value: u32) {
    for k in 0..num_bits {
        let idx = INFO_BITS_PER_WORD + start + k;
        let (w, within) = (idx / INFO_BITS_PER_WORD, idx % INFO_BITS_PER_WORD);
        // The first transmitted bit of a field is its most significant bit.
        if (value >> (num_bits - 1 - k)) & 1 == 1 {
            words[w] |= 1 << (29 - within);
        }
    }
}

/// Regression: the handover decoder used to read the TOW-COUNT and subframe ID
/// out of raw word bits 0..21, which overlaps the six parity bits at bits 0..5.
/// Changing only the parity of a word must not change the decoded values.
#[test]
fn parse_how_is_invariant_to_the_parity_bits() {
    let mut words = blank_subframe();
    set_how_field(&mut words, 0, 17, 0x1_5A3C); // TOW-COUNT
    set_how_field(&mut words, 19, 3, 3); // subframe ID 3
    let (tow, sf_id) = parse_how(&words);
    assert_eq!((tow, sf_id), (0x1_5A3C, 3));

    for parity in 0..64u32 {
        let mut w = words.clone();
        w[1] = (w[1] & !0x3F) | parity;
        assert_eq!(
            parse_how(&w),
            (tow, sf_id),
            "parity pattern {parity:#06b} leaked into the decoded handover word"
        );
    }
}

/// The first transmitted information bit is the most significant bit of a field
/// that starts there, so a lone bit there must decode to the field's MSB.
/// 8-bit field, MSB set -> 0b1000_0000 = 128.
#[test]
fn info_bit_accessor_places_the_first_transmitted_bit_at_the_field_msb() {
    let mut words = blank_subframe();
    words[0] |= 1 << 29; // first information bit of word 1
    assert_eq!(extract_info_bits(&words, 0, 8), 128);
    assert_eq!(extract_info_bits(&words, 0, 1), 1);

    // The last information bit of a word is word bit 6 and is the field LSB.
    words[0] = 1 << 6;
    assert_eq!(extract_info_bits(&words, 0, 24), 1);
    assert_eq!(extract_info_bits(&words, 23, 1), 1);
}

/// The information accessor must be blind to the six parity bits, while the
/// legacy raw-word accessor is not. This pins the defect the fix removed.
#[test]
fn parity_bits_are_invisible_to_the_info_accessor() {
    let words = vec![0x3Fu32; 10]; // every bit set, i.e. all parity
    for start in (0..200).step_by(7) {
        assert_eq!(
            extract_info_bits(&words, start, 4),
            0,
            "info accessor read parity at start {start}"
        );
    }
    assert_ne!(
        extract_bits(&words, 0, 4),
        0,
        "the legacy raw-bit accessor is expected to read parity here"
    );
}

/// Information indices map onto word bits 6..29 and nowhere else.
#[test]
fn info_bit_accessor_covers_exactly_the_information_window() {
    let mut words = blank_subframe();
    for word_bit in 0..30usize {
        words[0] = 1 << word_bit;
        let field = extract_info_bits(&words, 0, INFO_BITS_PER_WORD);
        let touches_parity = word_bit < 6;
        if touches_parity {
            assert_eq!(field, 0, "word bit {word_bit} is a parity bit and must be ignored");
        } else {
            // Word bit 29 is information index 0 (the field MSB) and word bit 6
            // is index 23 (the field LSB).
            assert_eq!(field, 1 << (word_bit - 6), "word bit {word_bit} is information");
        }
    }
}

/// A field wider than one word must continue at the *top* of the next word.
#[test]
fn info_bit_fields_span_word_boundaries() {
    let mut words = blank_subframe();
    // 32-bit field starting at information index 20: information indices
    // 20..52, i.e. 4 bits of word 1 (20..23), 24 bits of word 2 (24..47) and
    // 4 bits of word 3 (48..51). The field LSB is at word 1 bit 29-20 = 9;
    // the field MSB (bit 31) is at word 3 bit 29-3 = 26.
    words[0] |= 1 << 9;
    words[2] |= 1 << 26;
    assert_eq!(extract_info_bits(&words, 20, 32), 0x8000_0001);
    assert_eq!(extract_info_bits(&words, 20, 4), 0b1000, "MSB of the sub-field");
    assert_eq!(extract_info_bits(&words, 24, 24), 0);
    assert_eq!(extract_info_bits(&words, 48, 4), 1);
}

#[test]
fn sign_scale_handles_the_extremes_of_a_field() {
    // 3-bit field, scale 1: 3 -> 3.0, 4 (0b100) -> -4.0, 7 -> -1.0
    assert_eq!(sign_scale(3, 3, 1.0), 3.0);
    assert_eq!(sign_scale(4, 3, 1.0), -4.0);
    assert_eq!(sign_scale(7, 3, 1.0), -1.0);
    // 8-bit field: 0x80 -> -128, 0x7F -> 127
    assert_eq!(sign_scale(0x80, 8, 1.0), -128.0);
    assert_eq!(sign_scale(0x7F, 8, 1.0), 127.0);
    // 32-bit field at the LSB carries the full signed range of the ephemeris
    // angles: 0x8000_0000 -> -2147483648, 0x7FFF_FFFF -> 2147483647
    assert_eq!(sign_scale(0x8000_0000, 32, 1.0), -2147483648.0);
    assert_eq!(sign_scale(0x7FFF_FFFF, 32, 1.0), 2147483647.0);
}

#[test]
fn sign_scale_ignores_bits_above_the_field_width() {
    // Callers pass values already sliced to the field width; the mask makes the
    // helper robust if a wider value is ever handed in.
    assert_eq!(sign_scale(0xFFFF_FFFF, 8, 1.0), -1.0);
    assert_eq!(sign_scale(0x8000_0001, 8, 1.0), 1.0);
}

/// The 10-bit subframe week counts modulo 1024; a 2020 broadcast must not be
/// stored as week 62 (which is 1981).
#[test]
fn ten_bit_week_is_resolved_before_use() {
    use crate::rtcm3::ephemeris::adjust_gps_week;
    assert_eq!(adjust_gps_week(62, 2110), 2110);
    assert_eq!(adjust_gps_week(2110, 2110), 2110);
    // WEEK_REFERENCE is 2048, so a raw 62 maps to 1086.
    assert_eq!(adjust_gps_week(62, WEEK_REFERENCE), 2110); // nearest to 2048
}
