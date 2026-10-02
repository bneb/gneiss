//! Tests for the RTCM 3 type 1019 GPS ephemeris decoder.
//!
//! The payload layout asserted here was derived by hand from RTKLIB
//! `rtcm3.c` `decode_type1019`, which reads the same fields starting at frame
//! bit 24 (= payload bit 0). Expected values are hand arithmetic written next
//! to each assertion, never a re-evaluation of the implementation expression.

use super::*;

/// (payload bit offset, width, decoded-field slot) for every ephemeris field.
const LAYOUT: &[(usize, usize, usize)] = &[
    (18, 10, S_WEEK),
    (34, 14, S_IDOT),
    (48, 8, S_IODE),
    (56, 16, S_TOC_S),
    (72, 8, S_AF2),
    (80, 16, S_AF1),
    (96, 22, S_AF0),
    (118, 10, S_IODC),
    (128, 16, S_CRS),
    (144, 16, S_DELTA_N),
    (160, 32, S_M0),
    (192, 16, S_CUC),
    (208, 32, S_E),
    (240, 16, S_CUS),
    (256, 32, S_SQRT_A),
    (288, 16, S_TOE_S),
    (304, 16, S_CIC),
    (320, 32, S_OMEGA0),
    (352, 16, S_CIS),
    (368, 32, S_I0),
    (400, 16, S_CRC),
    (416, 32, S_OMEGA),
    (448, 24, S_OMEGA_DOT),
    (472, 8, S_TGD),
];

const S_WEEK: usize = 1;
const S_IDOT: usize = 2;
const S_IODE: usize = 3;
const S_TOC_S: usize = 4;
const S_AF2: usize = 5;
const S_AF1: usize = 6;
const S_AF0: usize = 7;
const S_IODC: usize = 8;
const S_CRS: usize = 9;
const S_DELTA_N: usize = 10;
const S_M0: usize = 11;
const S_CUC: usize = 12;
const S_E: usize = 13;
const S_CUS: usize = 14;
const S_SQRT_A: usize = 15;
const S_TOE_S: usize = 16;
const S_CIC: usize = 17;
const S_OMEGA0: usize = 18;
const S_CIS: usize = 19;
const S_I0: usize = 20;
const S_CRC: usize = 21;
const S_OMEGA: usize = 22;
const S_OMEGA_DOT: usize = 23;
const S_TGD: usize = 24;
const N_SLOTS: usize = 25;

/// GPS week 2110 starts 2020-12-20 00:00:00 GPST; 2110 = 2*1024 + 62.
fn ref_week() -> GpsTime {
    GpsTime::new(2110, 0.0)
}

fn bits_with(set_bit: Option<usize>) -> BitVec<u8, Msb0> {
    let mut bits = bitvec![u8, Msb0; 0; TYPE_1019_PAYLOAD_BITS];
    bits[0..12].store_be(1019_u16);
    if let Some(off) = set_bit {
        bits.set(off, true);
    }
    bits
}

/// Builds a 488-bit payload, writing `raw` into every field listed as
/// `(field slot, raw value)`. Raw values are masked to the field width, so
/// negative values may be supplied as `i64 as u64`.
fn payload_with(fields: &[(usize, u64)]) -> Vec<u8> {
    let mut bits = bits_with(None);
    for &(slot, raw) in fields {
        let &(off, width, _) = LAYOUT.iter().find(|e| e.2 == slot).expect("known field slot");
        bits[off..off + width].store_be(raw & ((1u64 << width) - 1));
    }
    bits.into_vec()
}

/// Flattens an ephemeris into a slot vector so tests can assert which field a
/// bit landed in. Integer fields are widened to `f64`; every field is exactly
/// zero for an all-zero payload except the resolved week.
fn decoded_slots(e: &GpsEphemeris) -> [f64; N_SLOTS] {
    [
        e.sat.prn as f64,
        e.toe.week as f64,
        e.idot,
        e.iode as f64,
        e.toc.tow,
        e.af2,
        e.af1,
        e.af0,
        e.iodc as f64,
        e.crs,
        e.delta_n,
        e.m0,
        e.cuc,
        e.e,
        e.cus,
        e.sqrt_a,
        e.toe.tow,
        e.cic,
        e.omega0,
        e.cis,
        e.i0,
        e.crc,
        e.omega,
        e.omega_dot,
        e.tgd,
    ]
}

/// Golden vector: every raw field value below is written explicitly and each
/// expected float is the closed form shown in the comment.
#[test]
fn decodes_hand_encoded_golden_vector() {
    let payload = payload_with(&[
        (S_WEEK, 62), // 2110 mod 1024, resolved against reference week 2110
        (S_IDOT, (-3i64) as u64),
        (S_IODE, 42),
        (S_TOC_S, 32400), // 32400 * 16 s = 518400.0 s
        (S_AF2, (-4i64) as u64),
        (S_AF1, 120),
        (S_AF0, 1_000_000),
        (S_IODC, 900),
        (S_CRS, 32), // 32 * 2^-5 = 1.0 m
        (S_DELTA_N, (-100i64) as u64),
        (S_M0, 1_000_000_000),
        (S_CUC, 20),
        (S_E, 100), // 100 * 2^-33 = 1.1641532182693481e-8
        (S_CUS, (-30i64) as u64),
        (S_SQRT_A, 2_701_846_400), // * 2^-19 = 5153.363037109375 m^0.5
        (S_TOE_S, 32400),          // 32400 * 16 s = 518400.0 s
        (S_CIC, 5),
        (S_OMEGA0, (-2_000_000_000i64) as u64),
        (S_CIS, (-7i64) as u64),
        (S_I0, 950_000_000),
        (S_CRC, 100), // 100 * 2^-5 = 3.125 m
        (S_OMEGA, 123_456_789),
        (S_OMEGA_DOT, (-1_000_000i64) as u64),
        (S_TGD, (-16i64) as u64),
    ]);

    let e = parse_1019_at(&payload, ref_week()).expect("golden 1019 must decode");
    assert_eq!(e.sat.constellation, Constellation::Gps);
    assert_eq!(e.sat.prn, 0, "PRN is left at zero by this vector");
    assert_eq!(e.toe.week, 2110);
    assert_eq!(e.toc.week, 2110);
    assert_eq!(e.toe.tow, 518400.0);
    assert_eq!(e.toc.tow, 518400.0);
    assert_eq!(e.iode, 42);
    assert_eq!(e.iodc, 900);
    // IDOT = -3 * 2^-43 * pi
    assert!((e.idot - (-1.0714732025882517e-12)).abs() < 1e-24);
    // af2 = -4 * 2^-55
    assert!((e.af2 - (-1.1102230246251565e-16)).abs() < 1e-30);
    // af1 = 120 * 2^-43
    assert!((e.af1 - 1.3642420526593924e-11).abs() < 1e-23);
    // af0 = 1000000 * 2^-31
    assert!((e.af0 - 0.00046566128730773926).abs() < 1e-18);
    // C_rs = 32 * 2^-5
    assert!((e.crs - 1.0).abs() < 1e-15);
    // delta_n = -100 * 2^-43 * pi
    assert!((e.delta_n - (-3.571577341960839e-11)).abs() < 1e-23);
    // M0 = 1e9 * 2^-31 * pi
    assert!((e.m0 - 1.4629180792671597).abs() < 1e-15);
    // C_uc = 20 * 2^-29
    assert!((e.cuc - 3.725290298461914e-08).abs() < 1e-22);
    // e = 100 * 2^-33
    assert!((e.e - 1.1641532182693481e-08).abs() < 1e-22);
    // C_us = -30 * 2^-29
    assert!((e.cus - (-5.587935447692871e-08)).abs() < 1e-22);
    // sqrt(A) = 2701846400 * 2^-19 = 5153.363037109375 (NOT the semi-major axis,
    // which is sqrt_a^2 = 26557150.59 m)
    assert!((e.sqrt_a - 5153.363037109375).abs() < 1e-9);
    // C_ic = 5 * 2^-29
    assert!((e.cic - 9.313225746154785e-09).abs() < 1e-21);
    // Omega_0 = -2e9 * 2^-31 * pi
    assert!((e.omega0 - (-2.9258361585343193)).abs() < 1e-15);
    // C_is = -7 * 2^-29
    assert!((e.cis - (-1.30385160446167e-08)).abs() < 1e-22);
    // i0 = 9.5e8 * 2^-31 * pi
    assert!((e.i0 - 1.3897721753038017).abs() < 1e-15);
    // C_rc = 100 * 2^-5
    assert!((e.crc - 3.125).abs() < 1e-15);
    // omega = 123456789 * 2^-31 * pi
    assert!((e.omega - 0.180607168636371).abs() < 1e-15);
    // Omega_dot = -1e6 * 2^-43 * pi
    assert!((e.omega_dot - (-3.571577341960839e-07)).abs() < 1e-19);
    // T_gd = -16 * 2^-31
    assert!((e.tgd - (-7.450580596923828e-09)).abs() < 1e-22);
}

/// The 10-bit week used to be read and discarded, so every toe and toc came
/// back as GPS week 0 (1980-01-06) no matter what the message said.
#[test]
fn week_is_carried_into_toe_and_toc() {
    let payload = payload_with(&[(S_WEEK, 62)]);
    let e = parse_1019_at(&payload, ref_week()).expect("must decode");
    assert_eq!(e.toe.week, 2110, "toe week 0 would be 1980, 40 years early");
    assert_eq!(e.toc.week, 2110, "toc week 0 would be 1980, 40 years early");
}

/// The identical truncated week resolves to a different full week per epoch —
/// exactly the rollover the 10-bit field hides.
#[test]
fn week_rollover_tracks_the_reference_epoch() {
    let payload = payload_with(&[(S_WEEK, 62)]);
    for full in [1086u32, 2110, 3134] {
        let e = parse_1019_at(&payload, GpsTime::new(full, 0.0)).expect("must decode");
        assert_eq!(e.toe.week, full, "raw week 62 must resolve to {full}");
    }
}

#[test]
fn adjust_gps_week_properties() {
    for reference in 1024..4200u32 {
        for raw in 0..1024u32 {
            let w = adjust_gps_week(raw, reference);
            assert_eq!(w % 1024, raw,  "raw {raw} must survive the mapping");
            let dist = (w as i64 - reference as i64).abs();
            assert!(dist <= 512, "w={w} ref={reference} dist={dist} exceeds half a cycle");
        }
    }
    // Closed form: the unique week within half a cycle of the reference.
    // Below week 1024 the correct week can precede the GPS epoch, so the
    // helper clamps at 0 instead of wrapping into a negative week.
    assert_eq!(adjust_gps_week(513, 0), 0);
    assert_eq!(adjust_gps_week(0, 0), 0);
    assert_eq!(adjust_gps_week(1023, 2100), 2047); // 2047 % 1024 == 1023
    assert_eq!(adjust_gps_week(100, 5000), 5220); // 5220 % 1024 == 100
    assert_eq!(adjust_gps_week(62, 2110), 2110); // 2110 = 2*1024 + 62
}

/// Any single set bit must change exactly one decoded field. A wrong shift or
/// mask width merges two fields or leaves one unreachable, and this catches it.
#[test]
fn single_bit_changes_exactly_one_field() {
    let base = decoded_slots(&parse_1019_at(&bits_with(None).into_vec(), ref_week()).unwrap());
    for &(off, width, slot) in LAYOUT {
        let e = parse_1019_at(&bits_with(Some(off)).into_vec(), ref_week()).expect("must decode");
        let got = decoded_slots(&e);
        let changed: Vec<usize> = (0..N_SLOTS).filter(|&i| got[i] != base[i]).collect();
        assert_eq!(
            changed,
            vec![slot],
            "bit at payload offset {off} (w={width}) must land only in field {slot}"
        );
    }
}

#[test]
fn prn_field_is_parsed() {
    let mut bits = bits_with(None);
    bits[12..18].store_be(21_u8);
    let e = parse_1019_at(&bits.into_vec(), ref_week()).expect("must decode");
    assert_eq!(e.sat.prn, 21);
}

/// Every truncation shorter than the full 488-bit payload must be rejected
/// rather than silently decoded from whatever happens to be in the buffer.
#[test]
fn rejects_every_truncation_below_the_full_payload() {
    let payload = payload_with(&[(S_CRS, 32)]);
    for len in 0..TYPE_1019_PAYLOAD_BITS {
        assert_eq!(
            parse_1019_at(&payload[..len / 8], ref_week()),
            Err(RtcmParseError::Incomplete),
            "{len} bits must be rejected, not zero-filled"
        );
    }
    assert!(parse_1019_at(&payload, ref_week()).is_ok());
}

#[test]
fn rejects_foreign_message_number() {
    let mut bits = bits_with(None);
    bits[0..12].store_be(1005_u16);
    assert_eq!(
        parse_1019_at(&bits.into_vec(), ref_week()),
        Err(RtcmParseError::UnsupportedMsmType)
    );
}

/// Each raw below is the most negative value of its field, so a decoder that
/// forgets to sign-extend (or sign-extends at the wrong width) fails here.
#[test]
fn negative_scalars_sign_extend_across_field_widths() {
    let payload = payload_with(&[
        (S_AF2, 0x80),           // 8-bit  -> -128    * 2^-55
        (S_AF1, 0x8000),         // 16-bit -> -32768  * 2^-43
        (S_AF0, 0x200000),       // 22-bit -> -2097152 * 2^-31
        (S_OMEGA_DOT, 0x800000), // 24-bit -> -8388608 * 2^-43 * pi
        (S_IDOT, 0x2000),        // 14-bit -> -8192    * 2^-43 * pi
        (S_TGD, 0x80),           // 8-bit  -> -128     * 2^-31
    ]);
    let e = parse_1019_at(&payload, ref_week()).expect("must decode");
    assert!((e.af2 - (-128.0 * 2.7755575615628914e-17)).abs() < 1e-30);
    assert!((e.af1 - (-32768.0 * 1.1368683772161603e-13)).abs() < 1e-23);
    assert!((e.af0 - (-2097152.0 * 4.656612873077393e-10)).abs() < 1e-18);
    assert!((e.omega_dot - (-8388608.0 * 1.1368683772161603e-13 * core::f64::consts::PI)).abs() < 1e-19);
    assert!((e.idot - (-8192.0 * 1.1368683772161603e-13 * core::f64::consts::PI)).abs() < 1e-22);
    assert!((e.tgd - (-128.0 * 4.656612873077393e-10)).abs() < 1e-22);
}

/// The no-argument wrapper must not silently reintroduce the discarded-week bug:
/// for a payload whose raw week is 62 it must return the week nearest its
/// documented reference, never week 0.
#[test]
fn convenience_wrapper_resolves_the_week() {
    let e = parse_1019(&payload_with(&[(S_WEEK, 62)])).expect("must decode");
    assert_ne!(e.toe.week, 0, "week 0 is 1980-01-06 and can never be a live ephemeris");
    assert_eq!(e.toe.week, e.toc.week);
}
