//! Parser for the CODE/AIUB Bernese monthly and daily GNSS DCB text files
//! (`P1C1YYMM.DCB`, `P1P2YYMM_ALL.DCB`, `comDDDYY.dcb`).
//!
//! Layout, verified against the real files in `datasets/profile_d_f9p/`
//! (`P1C12011.DCB`, `P1P22011_ALL.DCB`, `com21374.dcb`):
//!
//! ```text
//! DIFFERENTIAL (P1-C1) CODE BIASES FOR SATELLITES AND RECEIVERS:
//! PRN / STATION NAME  VALUE (NS)  RMS (NS)  yyyy mm dd hh mm ss  yyyy mm dd hh mm ss
//! G01                     1.614     0.005  2020 12 24 00 00 00  2020 12 25 00 00 00
//! G     ZECK 12351M001    7.962     0.035
//! ```
//!
//! Satellite rows carry a 3-character PRN and a blank station field;
//! receiver rows start with the constellation letter plus a station name.
//! Only satellite biases are emitted.
//!
//! Sign: `DCB(X-Y) = OSB(X) - OSB(Y)`, and the downstream consumer
//! *subtracts* the returned bias from the raw observation, so the emitted
//! `Osb` is `-DCB` with the first signal taken as the zero reference.
//! Cross-checked against `datasets/profile_d_f9p/COD0MGXFIN_..._OSB.BIA`
//! (same centre, same day 359/2020): G01 DCB `+1.614` ns, G01 `C1C` OSB
//! `-1.6145` ns; G02 `-1.585` / `+1.5851`; G05 `+1.134` / `-1.1342`.

use std::io::BufRead;
use std::str::FromStr;
use gneiss_core::sat::{Constellation, SatelliteId};
use gneiss_core::time::GpsTime;
use gneiss_core::obs::ObsCode;
use crate::sinex_bia::{BiasRecord, BiasType};

/// Byte windows of a data line, 0-based half-open, as laid out above.
const PRN_COL: (usize, usize) = (0, 3);
const STATION_COL: (usize, usize) = (3, 20);
const VALUE_COL: (usize, usize) = (20, 35);
const RMS_COL: (usize, usize) = (35, 48);
/// The first `yyyy mm dd hh mm ss` stamp starts at byte 51 and the second
/// at byte 72 in a 91-byte line, so these windows hold exactly one stamp
/// each (19 characters) and no fragment of its neighbour.
const START_COL: (usize, usize) = (48, 70);
const END_COL: (usize, usize) = (71, 91);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DcbDiffType {
    P1C1,
    P2C2,
    P1P2,
    P1C2,
}

/// Classify a `DIFFERENTIAL (...)` block header.
///
/// The real CODE headers name exactly one signal pair, e.g.
/// `DIFFERENTIAL (P1-C1) CODE BIASES FOR SATELLITES AND RECEIVERS:`.
/// Headers the caller cannot act on — the daily files also carry a
/// `DIFFERENTIAL (INTER-FREQ)` receiver block — return `None` so the
/// block is skipped instead of mislabelled.
pub fn parse_dcb_diff_type(line: &str) -> Option<DcbDiffType> {
    if line.contains("(P1-C1)") {
        Some(DcbDiffType::P1C1)
    } else if line.contains("(P2-C2)") {
        Some(DcbDiffType::P2C2)
    } else if line.contains("(P1-P2)") {
        Some(DcbDiffType::P1P2)
    } else if line.contains("(P1-C2)") {
        Some(DcbDiffType::P1C2)
    } else {
        None
    }
}

/// Parse a 3-character PRN field (`G01`, `R24`) into a satellite id.
///
/// Operates on bytes, so a non-ASCII field can never slice a `str` at a
/// non-boundary. The Bernese/CODE files always write exactly one letter
/// plus two digits, so anything shorter or longer is rejected.
pub fn parse_sat_id(s: &str) -> Option<SatelliteId> {
    let b = s.trim().as_bytes();
    if b.len() != 3 || !b[1].is_ascii_digit() || !b[2].is_ascii_digit() {
        return None;
    }
    let constellation = match b[0] {
        b'G' => Constellation::Gps,
        b'R' => Constellation::Glonass,
        b'E' => Constellation::Galileo,
        b'C' => Constellation::Beidou,
        b'J' => Constellation::Qzss,
        _ => return None,
    };
    Some(SatelliteId { constellation, prn: (b[1] - b'0') * 10 + (b[2] - b'0') })
}

/// Decode one `yyyy mm dd hh mm ss` validity stamp.
fn parse_ymd_hms(s: &str) -> Option<GpsTime> {
    let parts: Vec<&str> = s.split_whitespace().collect();
    if parts.len() < 6 {
        return None;
    }
    let y = parts[0].parse::<i32>().ok()?;
    let m = parts[1].parse::<i32>().ok()?;
    let d = parts[2].parse::<i32>().ok()?;
    let hh = parts[3].parse::<i32>().ok()?;
    let mm = parts[4].parse::<i32>().ok()?;
    let ss = parts[5].parse::<f64>().ok()?;
    Some(GpsTime::from_calendar(y, m, d, hh, mm, ss))
}

/// Validity interval declared by a monthly solution header, e.g.
/// `CODE'S MONTHLY GNSS P1-P2 DCB SOLUTION, YEAR 2020, MONTH 11`.
///
/// Monthly files carry no per-row stamps, so without this every record
/// would inherit the parser's 50-year fallback window.
fn header_span(line: &str) -> Option<(GpsTime, GpsTime)> {
    let after = line.split_once("YEAR ")?.1;
    // The real header writes `YEAR 2020, MONTH 11`, so the year token
    // carries the separating comma.
    let year = after.split_whitespace().next()?.trim_end_matches(',').parse::<i32>().ok()?;
    let month = after.split_once("MONTH")?.1.split_whitespace().next()?.parse::<i32>().ok()?;
    if !(1..=12).contains(&month) || !(1980..=2100).contains(&year) {
        return None;
    }
    let start = GpsTime::from_calendar(year, month, 1, 0, 0, 0.0);
    let (next_year, next_month) = if month == 12 { (year + 1, 1) } else { (year, month + 1) };
    Some((start, GpsTime::from_calendar(next_year, next_month, 1, 0, 0, 0.0)))
}

/// Validity interval of one data row, falling back to the file's own
/// interval when the row carries no decodable stamps.
fn line_span(line: &str, file_span: (GpsTime, GpsTime)) -> (GpsTime, GpsTime) {
    let stamps = line
        .get(START_COL.0..START_COL.1)
        .and_then(parse_ymd_hms)
        .zip(line.get(END_COL.0..END_COL.1).and_then(parse_ymd_hms));
    stamps.unwrap_or(file_span)
}

fn map_diff_type_to_bias(
    diff_type: DcbDiffType,
    val_ns: f64,
) -> Option<(BiasType, ObsCode, Option<ObsCode>, f64)> {
    match diff_type {
        DcbDiffType::P1C1 => {
            let obs1 = ObsCode::from_str("C1C").ok()?;
            Some((BiasType::Osb, obs1, None, -val_ns))
        }
        DcbDiffType::P2C2 => {
            let obs1 = ObsCode::from_str("C2C").ok()?;
            Some((BiasType::Osb, obs1, None, -val_ns))
        }
        DcbDiffType::P1P2 => {
            let obs1 = ObsCode::from_str("C1W").ok()?;
            let obs2 = ObsCode::from_str("C2W").ok()?;
            Some((BiasType::Dcb, obs1, Some(obs2), val_ns))
        }
        DcbDiffType::P1C2 => {
            let obs1 = ObsCode::from_str("C2W").ok()?;
            Some((BiasType::Osb, obs1, None, -val_ns))
        }
    }
}

/// Decode one satellite row, or `None` for header, rule and receiver rows.
pub fn parse_dcb_sat_line(
    line: &str,
    diff_type: DcbDiffType,
    file_span: (GpsTime, GpsTime),
) -> Option<BiasRecord> {
    if line.get(STATION_COL.0..STATION_COL.1)?.trim() != "" {
        return None;
    }
    let sat = parse_sat_id(line.get(PRN_COL.0..PRN_COL.1)?)?;
    let rms_end = RMS_COL.1.min(line.len());
    let val_ns = line.get(VALUE_COL.0..VALUE_COL.1)?.trim().parse::<f64>().ok()?;
    let std_dev = line
        .get(RMS_COL.0..rms_end)
        .and_then(|s| s.trim().parse::<f64>().ok())
        .unwrap_or(0.0);
    let (bias_type, obs1, obs2, value) = map_diff_type_to_bias(diff_type, val_ns)?;
    let (start_time, end_time) = line_span(line, file_span);
    Some(BiasRecord {
        bias_type,
        sat,
        station: None,
        obs1,
        obs2,
        start_time,
        end_time,
        unit: "ns".to_string(),
        value,
        std_dev,
    })
}

pub fn parse_bernese_dcb<R: BufRead>(mut reader: R) -> Result<Vec<BiasRecord>, String> {
    let mut records = Vec::new();
    let mut current_diff_type: Option<DcbDiffType> = None;
    let mut file_span: Option<(GpsTime, GpsTime)> = None;
    let mut line_buf = String::new();
    while reader.read_line(&mut line_buf).map_err(|e| e.to_string())? > 0 {
        let trimmed = line_buf.trim_end();
        let span = *file_span.get_or_insert_with(|| header_span(trimmed).unwrap_or(fallback_span()));
        if trimmed.starts_with("DIFFERENTIAL") {
            current_diff_type = parse_dcb_diff_type(trimmed);
        } else if let Some(diff_type) = current_diff_type {
            if let Some(rec) = parse_dcb_sat_line(trimmed, diff_type, span) {
                records.push(rec);
            }
        }
        line_buf.clear();
    }
    Ok(records)
}

/// Window used when a file declares no epoch of its own.
fn fallback_span() -> (GpsTime, GpsTime) {
    (
        GpsTime::from_calendar(2000, 1, 1, 0, 0, 0.0),
        GpsTime::from_calendar(2050, 1, 1, 0, 0, 0.0),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sinex_bia::SinexBias;
    use std::io::Cursor;
    use std::path::PathBuf;

    const TOL: f64 = 1e-9;

    /// Real 91-byte daily row, copied verbatim from
    /// `datasets/profile_d_f9p/com21374.dcb` line 7.
    const REAL_DAILY_ROW: &str = "G01                           1.614       0.005    2020 12 24 00 00 00  2020 12 25 00 00 00";
    /// Real 47-byte monthly row, `P1C12011.DCB` line 7.
    const REAL_MONTHLY_ROW: &str = "G01                           1.496       0.005";
    /// Real header of `P1C12011.DCB` line 0.
    const REAL_MONTHLY_HEADER: &str =
        "CODE'S MONTHLY GPS P1-C1 DCB SOLUTION, YEAR 2020, MONTH 11       04-DEC-20 06:46";


    fn real_file(name: &str) -> Option<std::fs::File> {
        let rel = format!("../../datasets/profile_d_f9p/{name}");
        let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel);
        p.exists().then(|| std::fs::File::open(p).expect("real DCB dataset must open"))
    }

    fn parse_real(name: &str) -> Option<Vec<BiasRecord>> {
        Some(parse_bernese_dcb(std::io::BufReader::new(real_file(name)?)).expect("parse"))
    }

    fn assert_time_eq(actual: GpsTime, y: i32, m: i32, d: i32, what: &str) {
        let want = GpsTime::from_calendar(y, m, d, 0, 0, 0.0);
        assert!(
            (actual - want).abs() < 1.0,
            "{what}: got week {} tow {:.3}, want week {} ({y}-{m}-{d})",
            actual.week,
            actual.tow,
            want.week
        );
    }

    // --- Golden vectors decoded from the real CODE files. ---

    #[test]
    fn real_daily_row_decodes_stamps_from_their_true_columns() {
        let span = fallback_span();
        let rec = parse_dcb_sat_line(REAL_DAILY_ROW, DcbDiffType::P1C1, span)
            .expect("91-byte real row must decode");
        // "2020 12 24 00 00 00" occupies bytes 51..70 of the row, and the
        // closing stamp "2020 12 25 00 00 00" bytes 72..91.
        assert_time_eq(rec.start_time, 2020, 12, 24, "start stamp");
        assert_time_eq(rec.end_time, 2020, 12, 25, "end stamp");
        assert_eq!(rec.sat, SatelliteId { constellation: Constellation::Gps, prn: 1 });
        // DCB(P1-C1) = +1.614 ns; C1W is the zero reference, so the C1C OSB
        // is -1.614 (CODE's own OSB file lists -1.6145 ns for G01/C1C).
        assert!((rec.value - -1.614).abs() < TOL, "value {}", rec.value);
        assert!((rec.std_dev - 0.005).abs() < TOL, "rms {}", rec.std_dev);
        assert_eq!(rec.obs1, ObsCode::from_str("C1C").unwrap());
        assert_eq!(rec.unit, "ns");
    }

    #[test]
    fn real_daily_file_bias_is_valid_only_inside_its_own_day() {
        let Some(recs) = parse_real("com21374.dcb") else { return };
        let bia = SinexBias::new(recs);
        let g01 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let c1c = ObsCode::from_str("C1C").unwrap();
        let at = |d: i32| GpsTime::from_calendar(2020, 12, d, 12, 0, 0.0);
        let (before, inside, after) = (at(20), at(24), at(26));
        // The file declares 2020-12-24 00:00 .. 2020-12-25 00:00 for G01.
        assert_eq!(bia.get_bias(g01, c1c, inside), Some(-1.614));
        assert_eq!(bia.get_bias(g01, c1c, before), None, "row must not apply 4 days early");
        assert_eq!(bia.get_bias(g01, c1c, after), None, "row must not apply 1 day late");
    }

    #[test]
    fn real_monthly_file_takes_its_validity_from_the_header() {
        let Some(recs) = parse_real("P1C12011.DCB") else { return };
        let g01 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let row = recs.iter().find(|r| r.sat == g01).expect("G01 row");
        // Header says YEAR 2020, MONTH 11 -> [2020-11-01, 2020-12-01).
        assert_time_eq(row.start_time, 2020, 11, 1, "monthly start");
        assert_time_eq(row.end_time, 2020, 12, 1, "monthly end");
        // Header says +1.496 ns for G01 (real file, line 7).
        assert!((row.value - -1.496).abs() < TOL, "value {}", row.value);
        assert!((row.std_dev - 0.005).abs() < TOL, "rms {}", row.std_dev);
    }

    #[test]
    fn header_span_handles_year_boundary_and_rejects_junk() {
        let (s, e) = header_span(REAL_MONTHLY_HEADER).expect("real header");
        assert_time_eq(s, 2020, 11, 1, "nov start");
        assert_time_eq(e, 2020, 12, 1, "nov end");
        let dec = "CODE'S MONTHLY GPS P1-C1 DCB SOLUTION, YEAR 2020, MONTH 12";
        assert_time_eq(header_span(dec).expect("december header").0, 2020, 12, 1, "dec start");
        assert_time_eq(header_span(dec).expect("december header").1, 2021, 1, 1, "year rolls");
        // The daily header carries a day-of-year, not a month; the legacy
        // 2008 "YEAR-MONTH 08-01" spelling is not claimed either.
        assert!(header_span("CODE'S MGEX (OSB) RESULTS FOR DAY 359, 2020").is_none());
        assert!(header_span("DIFFERENTIAL (P1-C1) CODE BIASES").is_none());
        assert!(header_span("YEAR 2020, MONTH 13").is_none());
        assert!(header_span("YEAR 0000, MONTH 01").is_none());
    }

    #[test]
    fn monthly_row_without_own_stamps_uses_the_file_span() {
        let (s, e) = header_span(REAL_MONTHLY_HEADER).unwrap();
        let rec = parse_dcb_sat_line(REAL_MONTHLY_ROW, DcbDiffType::P1C1, (s, e)).unwrap();
        assert_time_eq(rec.start_time, 2020, 11, 1, "inherited start");
        assert_time_eq(rec.end_time, 2020, 12, 1, "inherited end");
    }

    // --- Block header classification. ---

    /// Real block headers, verbatim from the four products in datasets/.
    const HEADERS: [(&str, Option<DcbDiffType>); 4] = [
        ("DIFFERENTIAL (P1-C1) CODE BIASES FOR SATELLITES AND RECEIVERS:", Some(DcbDiffType::P1C1)),
        ("DIFFERENTIAL (P2-C2) CODE BIASES FOR SATELLITES AND RECEIVERS:", Some(DcbDiffType::P2C2)),
        ("DIFFERENTIAL (P1-P2) CODE BIASES FOR SATELLITES AND RECEIVERS:", Some(DcbDiffType::P1P2)),
        // The daily file's second block names no supported signal pair.
        ("DIFFERENTIAL (INTER-FREQ) CODE BIASES FOR SATELLITES AND RECEIVERS:", None),
    ];

    #[test]
    fn real_block_headers_classify_as_written() {
        for (line, want) in HEADERS {
            assert_eq!(parse_dcb_diff_type(line), want, "{line}");
        }
        assert_eq!(parse_dcb_diff_type("PRN / STATION NAME        VALUE (NS)  RMS (NS)"), None);
    }

    #[test]
    fn unsupported_block_emits_nothing_from_a_real_file() {
        let Some(recs) = parse_real("com21374.dcb") else { return };
        // com21374.dcb holds 3059 data rows. The supported P1-C1 block is
        // G01..G32 (lines 7..38); the remaining 3027 rows sit in the
        // trailing `DIFFERENTIAL (INTER-FREQ)` receiver block, which names
        // no supported signal pair and must contribute nothing at all.
        assert_eq!(recs.len(), 32, "only the P1-C1 satellite block may emit");
        assert!(recs.iter().all(|r| r.obs1 == ObsCode::from_str("C1C").unwrap()));
        assert!(recs.iter().all(|r| r.sat.constellation == Constellation::Gps));
        assert!((recs[0].value - -1.614).abs() < TOL, "G01 {}", recs[0].value);
    }

    // --- PRN field decoding. ---

    #[test]
    fn sat_id_accepts_only_letter_plus_two_digits() {
        let g = parse_sat_id("G01").expect("G01");
        assert_eq!((g.constellation, g.prn), (Constellation::Gps, 1));
        let r = parse_sat_id(" R24 ").expect("R24 with padding");
        assert_eq!((r.constellation, r.prn), (Constellation::Glonass, 24));
        // Real CODE files never write a one-digit or four-digit PRN.
        assert!(parse_sat_id("G1").is_none());
        assert!(parse_sat_id("G001").is_none());
        assert!(parse_sat_id("G").is_none());
        assert!(parse_sat_id("g01").is_none(), "lowercase is not a Bernese PRN");
        assert!(parse_sat_id("S01").is_none(), "SBAS is not in the supported set");
        assert!(parse_sat_id("***").is_none());
        // Multi-byte UTF-8 must be rejected, never sliced mid-character.
        assert!(parse_sat_id("\u{03a9}01").is_none());
        assert!(parse_sat_id("G\u{00b0}1").is_none());
        assert!(parse_sat_id("\u{4f60}\u{597d}").is_none());
    }

    // --- Fixed-column slicing. ---

    #[test]
    fn real_receiver_rows_are_never_emitted_as_satellites() {
        let span = fallback_span();
        // P1P22011_ALL.DCB receiver row: constellation letter + station.
        assert!(parse_dcb_sat_line("G     ZECK 12351M001          7.962       0.035", DcbDiffType::P1P2, span).is_none());
        // com21374.dcb receiver row whose columns 1..3 read "R20", a
        // syntactically valid GLONASS PRN: only the non-blank station
        // field (bytes 3..20) keeps it out of the record list.
        let r20 = "R20   ZIM3 14001M008        -17.637       0.075    2020 12 24 00 00 00  2020 12 25 00 00 00";
        assert_eq!(r20.get(PRN_COL.0..PRN_COL.1), Some("R20"));
        assert!(parse_dcb_sat_line(r20, DcbDiffType::P1C1, span).is_none());
    }

    #[test]
    fn short_and_boundary_rows_behave_predictably() {
        let span = fallback_span();
        // One byte short of VALUE_COL.1: the value field is incomplete, so
        // the row is rejected rather than read as a truncated number.
        assert_eq!(REAL_MONTHLY_ROW[..VALUE_COL.1 - 1].len(), 34);
        assert!(parse_dcb_sat_line(&REAL_MONTHLY_ROW[..34], DcbDiffType::P1C1, span).is_none());
        // Exactly at VALUE_COL.1: the row parses and keeps every digit.
        let cut = &REAL_MONTHLY_ROW[..VALUE_COL.1];
        let rec = parse_dcb_sat_line(cut, DcbDiffType::P1C1, span).expect("35-byte row parses");
        assert!((rec.value - -1.496).abs() < TOL, "truncated value {}", rec.value);
        // A missing RMS is reported as 0.0 rather than borrowed from the
        // value column.
        let no_rms = "G01                           1.496";
        let rec = parse_dcb_sat_line(no_rms, DcbDiffType::P1C1, span).expect("value-only row");
        assert!((rec.value - -1.496).abs() < TOL);
        assert_eq!(rec.std_dev, 0.0, "absent RMS must not inherit the value");
        // A value shorter than the real field still decodes: the field is
        // right-justified, so "  1.49" is the number 1.49, not 1.496.
        let short = "G01                            1.49       0.005";
        let rec = parse_dcb_sat_line(short, DcbDiffType::P1C1, span).expect("short value field");
        assert!((rec.value - -1.49).abs() < TOL, "short value {}", rec.value);
        assert!((rec.std_dev - 0.005).abs() < TOL);
        // A '+' sign makes the token six bytes; right-justified so it
        // ends at byte 35 it still fits the 15-byte value window.
        let plus = format!("G01{:>12}{:>20}{:>13}", "", "+1.496", "0.005");
        assert_eq!(plus.get(20..35).map(str::trim), Some("+1.496"));
        assert_eq!(plus.get(35..48).map(str::trim), Some("0.005"));
        assert!((parse_dcb_sat_line(&plus, DcbDiffType::P1C1, span).unwrap().value - -1.496).abs() < TOL);
        assert!(parse_dcb_sat_line("\u{00b0}01                           1.496       0.005", DcbDiffType::P1C1, span).is_none());
    }

    #[test]
    fn negative_values_keep_their_sign_through_the_whole_file() {
        // P1P22011_ALL.DCB line 7: G01 DCB(P1-P2) = -6.858 ns. A DCB row
        // keeps its sign (it is already a difference of two signals).
        let text = format!(
            "{}\nDIFFERENTIAL (P1-P2) CODE BIASES FOR SATELLITES AND RECEIVERS:\n\
             PRN / STATION NAME        VALUE (NS)  RMS (NS)\n{}\n",
            REAL_MONTHLY_HEADER.replace("GPS P1-C1", "GNSS P1-P2"),
            "G01                          -6.858       0.008"
        );
        let recs = parse_bernese_dcb(Cursor::new(text)).unwrap();
        assert_eq!(recs.len(), 1);
        let r = &recs[0];
        assert_eq!(r.bias_type, BiasType::Dcb);
        assert_eq!(r.obs1, ObsCode::from_str("C1W").unwrap());
        assert_eq!(r.obs2, Some(ObsCode::from_str("C2W").unwrap()));
        assert!((r.value - -6.858).abs() < TOL, "value {}", r.value);
        assert!((r.std_dev - 0.008).abs() < TOL, "rms {}", r.std_dev);
    }

    #[test]
    fn crlf_and_blank_input_are_handled() {
        let block = format!(
            "DIFFERENTIAL (P1-C1) CODE BIASES FOR SATELLITES AND RECEIVERS:\n{row}\n",
            row = REAL_MONTHLY_ROW
        );
        let lf = parse_bernese_dcb(Cursor::new(&block)).unwrap();
        let crlf = parse_bernese_dcb(Cursor::new(block.replace('\n', "\r\n"))).unwrap();
        assert_eq!(lf.len(), 1, "one satellite row in the block");
        assert_eq!(crlf.len(), 1, "CRLF must parse identically");
        assert!((lf[0].value - -1.496).abs() < TOL);
        assert!((crlf[0].value - lf[0].value).abs() < TOL);
        // Header and rule lines alone carry no records.
        let header_only = format!("{h}\n{h}\n", h = REAL_MONTHLY_HEADER);
        assert_eq!(parse_bernese_dcb(Cursor::new(&header_only)).unwrap().len(), 0);
        assert_eq!(parse_bernese_dcb(Cursor::new("")).unwrap().len(), 0);
        assert_eq!(parse_bernese_dcb(Cursor::new("   \n\t\n")).unwrap().len(), 0);
    }

    #[test]
    fn real_p1p2_file_yields_one_dcb_row_per_listed_satellite() {
        let Some(recs) = parse_real("P1P22011_ALL.DCB") else { return };
        // 32 GPS PRNs plus 21 GLONASS slots; the 345 station rows and the
        // 80-column rule lines must not contribute.
        assert_eq!(recs.len(), 53, "satellite row count");
        assert!(recs.iter().all(|r| r.bias_type == BiasType::Dcb));
        assert!(recs.iter().all(|r| r.station.is_none()));
        let g01 = recs.iter().find(|r| r.sat.prn == 1 && r.sat.constellation == Constellation::Gps);
        assert!((g01.expect("G01").value - -6.858).abs() < TOL);
        let r24 = recs.iter().find(|r| r.sat.prn == 24 && r.sat.constellation == Constellation::Glonass);
        assert!((r24.expect("R24").value - 0.794).abs() < TOL, "R24 value");
    }

    /// DEFECT (outside this module's scope): `SinexBias::get_exact_bias`
    /// only returns `BiasType::Osb`, so every record this parser emits
    /// for a P1-P2 product is indexed and then unreachable. Re-enable when
    /// `sinex_bia::get_exact_bias` serves DCB rows too.
    #[ignore = "known defect: P1-P2 DCB rows are unreachable via get_exact_bias"]
    #[test]
    fn real_p1p2_rows_are_reachable_through_get_exact_bias() {
        let Some(recs) = parse_real("P1P22011_ALL.DCB") else { return };
        let bia = SinexBias::new(recs);
        let g01 = SatelliteId { constellation: Constellation::Gps, prn: 1 };
        let t = GpsTime::from_calendar(2020, 11, 15, 0, 0, 0.0);
        assert_eq!(bia.get_exact_bias(g01, ObsCode::from_str("C1W").unwrap(), t), Some(-6.858));
    }
}
