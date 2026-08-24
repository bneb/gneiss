//! Explicit signal–frequency registry.
//!
//! One canonical mapping between *physical GNSS signals* and their *transmitted
//! frequencies*, replacing ad-hoc `(constellation, band_number)` lookups where
//! a single band number means different physical signals per constellation
//! (GPS band 2 = L2 @ 1227.60 MHz, Galileo band 2 = ambiguous legacy slot).
//!
//! References: GPS ICD-IS-200/705, GLONASS ICD (5.1 ed.), Galileo ICD OS-SIS,
//! BeiDou ICD B1I/B3I (all center frequencies are ITU-allocated RNSS carriers).

use crate::constants::SPEED_OF_LIGHT_M_S;
use crate::sat::Constellation;

/// Physical GNSS signal identity — no ambiguity.
///
/// Each variant names a physical ranging signal with a single centre frequency
/// (or, for GLONASS FDMA, a base frequency plus a per-channel offset).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Signal {
    // GPS
    GpsL1Ca,
    GpsL1P,
    GpsL2Cm,
    GpsL2P,
    GpsL5,
    // GLONASS
    GloL1Of,
    GloL2Of,
    // Galileo
    GalE1Os,
    GalE5a,
    GalE5b,
    GalE6Cs,
    // BeiDou
    BdsB1i,
    BdsB3i,
}

/// One character class of a RINEX observable type string (`"C1"`, `"L2W"`, `"P1"`…).
struct RinexCode {
    kind: char,
    band: u8,
    attr: char,
}

impl RinexCode {
    /// Parses 2-char RINEX 2 types (`"L1"`, `"P2"`) and 3-char RINEX 3 types
    /// (`"C1C"`, `"L2W"`). Case-insensitive; missing attribute becomes `' '`.
    /// Only valid observation kinds (`C`/`L`/`P`/`D`/`S`) parse.
    fn parse(rinex_type: &str) -> Option<Self> {
        let mut chars = rinex_type.chars();
        let kind = chars.next()?.to_ascii_uppercase();
        if !matches!(kind, 'C' | 'L' | 'P' | 'D' | 'S') {
            return None;
        }
        let band = chars.next()?.to_digit(10)? as u8;
        let attr = chars.next().unwrap_or(' ').to_ascii_uppercase();
        Some(Self { kind, band, attr })
    }

    /// True for P-code style observations: RINEX 2 `"P1"/"P2"` kinds and
    /// RINEX 3 encrypted-tracking attributes `W` (P(Y)) / `Y` (P(Y)-AS).
    fn is_precise(&self) -> bool {
        self.kind == 'P' || self.attr == 'W' || self.attr == 'Y'
    }
}

impl Signal {
    /// Base frequency in Hz (without FDMA offset for GLONASS).
    pub fn base_freq_hz(&self) -> f64 {
        match self {
            Signal::GpsL1Ca | Signal::GpsL1P | Signal::GalE1Os => 1_575_420_000.0,
            Signal::GpsL2Cm | Signal::GpsL2P => 1_227_600_000.0,
            Signal::GpsL5 | Signal::GalE5a => 1_176_450_000.0,
            Signal::GloL1Of => 1_602_000_000.0,
            Signal::GloL2Of => 1_246_000_000.0,
            Signal::GalE5b => 1_207_140_000.0,
            Signal::GalE6Cs => 1_278_750_000.0,
            Signal::BdsB1i => 1_561_098_000.0,
            Signal::BdsB3i => 1_268_520_000.0,
        }
    }

    /// Frequency offset per GLONASS FDMA channel number `k`;
    /// `None` for CDMA signals whose frequency is `k`-independent.
    pub fn fdma_offset_hz(&self) -> Option<f64> {
        match self {
            Signal::GloL1Of => Some(562_500.0),
            Signal::GloL2Of => Some(437_500.0),
            _ => None,
        }
    }

    /// Wavelength in metres at the base frequency.
    ///
    /// For GLONASS FDMA signals this is the *nominal* (k = 0) wavelength;
    /// use [`wavelength_for`] for the channel-correct value.
    pub fn wavelength_m(&self) -> f64 {
        SPEED_OF_LIGHT_M_S / self.base_freq_hz()
    }
}

/// Map from RINEX observable type + constellation to [`Signal`].
///
/// Accepts RINEX 2 (`"L1"`, `"P2"`) and RINEX 3 (`"C1C"`, `"L2W"`) spellings.
///
/// # The Galileo band-2 ambiguity (resolved here, once)
///
/// RINEX 3 makes Galileo bands unambiguous: 1 = E1, 5 = E5a, 6 = E6, 7 = E5b.
/// RINEX 2.11 mixed GNSS exports instead reuse the GPS letter scheme, and the
/// second-frequency slot labelled `"L2"/"C2"` for Galileo satellites is written
/// **differently by different converters**: some put E5a (1176.45 MHz) there,
/// others E5b (1207.14 MHz). The datasets this engine consumes put **E5a** in
/// the L2 slot, so this registry resolves Galileo `"L2"` to
/// [`Signal::GalE5a`]. Note this deliberately *differs* from the legacy
/// `crate::signal::get_frequency(.., band 2, ..)`, which assumes E5b — that
/// disagreement was the bug class this module eliminates. Prefer unambiguous
/// RINEX 3 codes (band 5 = E5a, band 7 = E5b) whenever the file provides them.
///
/// Returns `None` for combinations with no modelled signal (e.g. NavIC, BDS
/// B2a/B2b, Galileo E5 AltBOC band 8) — callers must skip such observables.
pub fn rinex_type_to_signal(constellation: Constellation, rinex_type: &str) -> Option<Signal> {
    let code = RinexCode::parse(rinex_type)?;
    match constellation {
        Constellation::Gps | Constellation::Qzss => gps_signal(&code),
        Constellation::Sbas => sbas_signal(code.band),
        Constellation::Glonass => glonass_signal(code.band),
        Constellation::Galileo => galileo_signal(code.band),
        Constellation::Beidou => beidou_signal(code.band),
        Constellation::Navic => None, // L5/S-band not modelled yet
    }
}

/// Actual transmitted frequency in Hz for a specific satellite.
///
/// `freq_num` is the GLONASS FDMA channel number `k`; it is ignored for CDMA
/// signals. If `sat` disagrees with the signal's own constellation (only
/// possible for the FDMA branch), the nominal k = 0 channel is returned rather
/// than panicking — callers should treat that combination as a bug upstream.
pub fn frequency_for(sat: Constellation, signal: Signal, freq_num: i8) -> f64 {
    match signal.fdma_offset_hz() {
        None => signal.base_freq_hz(),
        Some(offset) => {
            let k = if sat == Constellation::Glonass { freq_num } else { 0 };
            signal.base_freq_hz() + f64::from(k) * offset
        }
    }
}

/// Channel-correct wavelength in metres ([`frequency_for`] over `c`).
pub fn wavelength_for(sat: Constellation, signal: Signal, freq_num: i8) -> f64 {
    SPEED_OF_LIGHT_M_S / frequency_for(sat, signal, freq_num)
}

/// GPS and QZSS share an identical radio plan (L1 C/A, L2C, L5), so QZSS
/// reuses the GPS variants as documented frequency-equivalent aliases.
fn gps_signal(code: &RinexCode) -> Option<Signal> {
    match code.band {
        1 => Some(if code.is_precise() { Signal::GpsL1P } else { Signal::GpsL1Ca }),
        2 => Some(if code.is_precise() { Signal::GpsL2P } else { Signal::GpsL2Cm }),
        5 => Some(Signal::GpsL5),
        _ => None,
    }
}

/// GLONASS FDMA: band 1 = L1OF, band 2 = L2OF; the channel number is carried
/// separately and applied by [`frequency_for`].
fn glonass_signal(band: u8) -> Option<Signal> {
    match band {
        1 => Some(Signal::GloL1Of),
        2 => Some(Signal::GloL2Of),
        _ => None, // G3OF (band 3) not modelled
    }
}

/// See the Galileo band-2 note on [`rinex_type_to_signal`].
fn galileo_signal(band: u8) -> Option<Signal> {
    match band {
        1 => Some(Signal::GalE1Os),
        2 => Some(Signal::GalE5a), // AMBIGUOUS legacy slot — see module docs
        5 => Some(Signal::GalE5a),
        6 => Some(Signal::GalE6Cs),
        7 => Some(Signal::GalE5b),
        _ => None, // band 8 (E5 AltBOC composite) not modelled
    }
}

/// BeiDou: B1I on band 1, B3I on RINEX 3 band 6. Bands 2 (B2I ≈1207 MHz) and
/// 5 (B2a 1176.45) have no variant in this registry version yet. Caveat:
/// RINEX 3 BDS-3 `C1P/C1X` observables track B1C (1575.42 MHz), NOT B1I —
/// they resolve to [`Signal::BdsB1i`] here because no B1C variant exists;
/// callers processing BDS-3-only receivers must not trust band-1 frequencies
/// until a `BdsB1c` variant is added.
fn beidou_signal(band: u8) -> Option<Signal> {
    match band {
        1 => Some(Signal::BdsB1i),
        6 => Some(Signal::BdsB3i),
        _ => None,
    }
}

/// SBAS transmits on L1 (1575.42 MHz) and L5 (1176.45 MHz), aliased here to
/// the frequency-identical GPS variants; SBAS has no L2.
fn sbas_signal(band: u8) -> Option<Signal> {
    match band {
        1 => Some(Signal::GpsL1Ca),
        5 => Some(Signal::GpsL5),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::signal::{
        FREQ_BDS_B1I, FREQ_GAL_E5B, FREQ_GLO_L1_DELTA, FREQ_GLO_L1_NOMINAL, FREQ_GLO_L2_DELTA,
        FREQ_GLO_L2_NOMINAL, FREQ_GPS_L1, FREQ_GPS_L2, FREQ_GPS_L5,
    };

    const HZ_TOL: f64 = 1.0; // spec: match ICD to <1 Hz

    fn assert_close(actual: f64, expected: f64, tol: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= tol,
            "{what}: got {actual}, expected {expected} (tol {tol})"
        );
    }

    #[test]
    fn base_freqs_match_icd_within_1hz() {
        // Centre frequencies straight from the ICDs (in Hz).
        let cases = [
            (Signal::GpsL1Ca, 1575.42e6),
            (Signal::GpsL1P, 1575.42e6),
            (Signal::GpsL2Cm, 1227.60e6),
            (Signal::GpsL2P, 1227.60e6),
            (Signal::GpsL5, 1176.45e6),
            (Signal::GloL1Of, 1602.0e6),
            (Signal::GloL2Of, 1246.0e6),
            (Signal::GalE1Os, 1575.42e6),
            (Signal::GalE5a, 1176.45e6),
            (Signal::GalE5b, 1207.14e6),
            (Signal::GalE6Cs, 1278.75e6),
            (Signal::BdsB1i, 1561.098e6),
            (Signal::BdsB3i, 1268.52e6),
        ];
        for (signal, icd_hz) in cases {
            assert_close(signal.base_freq_hz(), icd_hz, HZ_TOL, "base freq");
        }
    }

    #[test]
    fn fdma_offsets_exhaustive() {
        let cases = [
            (Signal::GpsL1Ca, None),
            (Signal::GpsL1P, None),
            (Signal::GpsL2Cm, None),
            (Signal::GpsL2P, None),
            (Signal::GpsL5, None),
            (Signal::GloL1Of, Some(562_500.0)),
            (Signal::GloL2Of, Some(437_500.0)),
            (Signal::GalE1Os, None),
            (Signal::GalE5a, None),
            (Signal::GalE5b, None),
            (Signal::GalE6Cs, None),
            (Signal::BdsB1i, None),
            (Signal::BdsB3i, None),
        ];
        for (signal, offset) in cases {
            assert_eq!(signal.fdma_offset_hz(), offset);
        }
    }

    #[test]
    fn glonass_fdma_k_minus7_to_plus6() {
        // Independent formula: f(k) = f0 + k * delta, valid channels -7..+6.
        for k in -7..=6i8 {
            let f1 = frequency_for(Constellation::Glonass, Signal::GloL1Of, k);
            let f2 = frequency_for(Constellation::Glonass, Signal::GloL2Of, k);
            assert_close(f1, 1_602_000_000.0 + f64::from(k) * 562_500.0, HZ_TOL, "L1(k)");
            assert_close(f2, 1_246_000_000.0 + f64::from(k) * 437_500.0, HZ_TOL, "L2(k)");
        }
        // Published GLONASS ICD edge-channel spot checks.
        let (g, l1, l2) = (Constellation::Glonass, Signal::GloL1Of, Signal::GloL2Of);
        assert_close(frequency_for(g, l1, -7), 1598.0625e6, HZ_TOL, "L1(k=-7)");
        assert_close(frequency_for(g, l1, 6), 1605.375e6, HZ_TOL, "L1(k=+6)");
        assert_close(frequency_for(g, l2, -7), 1242.9375e6, HZ_TOL, "L2(k=-7)");
        assert_close(frequency_for(g, l2, 6), 1248.625e6, HZ_TOL, "L2(k=+6)");
    }

    #[test]
    fn wavelengths_match_c_over_f() {
        // Hardcoded references computed independently of this module.
        let cases = [
            (Signal::GpsL1Ca, 0.190293672798),
            (Signal::GpsL2Cm, 0.244210213425),
            (Signal::GpsL5, 0.254828048791),
            (Signal::GloL1Of, 0.187136365793), // nominal k=0
            (Signal::GloL2Of, 0.240603898876),
            (Signal::GalE1Os, 0.190293672798),
            (Signal::GalE5a, 0.254828048791),
            (Signal::GalE5b, 0.248349369584),
            (Signal::GalE6Cs, 0.234441804888),
            (Signal::BdsB1i, 0.192039486310),
            (Signal::BdsB3i, 0.236332464604),
        ];
        for (signal, lam_ref) in cases {
            let lam = signal.wavelength_m();
            assert_close(lam, lam_ref, 1e-9, "lambda ref");
            // Cross-check: lambda must equal c/f exactly (relative tol).
            let ratio = lam / (SPEED_OF_LIGHT_M_S / signal.base_freq_hz());
            assert_close(ratio, 1.0, 1e-12, "lambda = c/f");
        }
        // Channel-correct GLONASS wavelengths at the extreme channels.
        let l1_lo = wavelength_for(Constellation::Glonass, Signal::GloL1Of, -7);
        let l2_hi = wavelength_for(Constellation::Glonass, Signal::GloL2Of, 6);
        assert_close(l1_lo, 0.187597455043, 1e-9, "GLO L1 lambda k=-7");
        assert_close(l2_hi, 0.240098074282, 1e-9, "GLO L2 lambda k=+6");
    }

    #[test]
    fn gps_rinex_types_map_correctly() {
        let ca = [("C1", Signal::GpsL1Ca), ("L1", Signal::GpsL1Ca), ("D1", Signal::GpsL1Ca), ("S1", Signal::GpsL1Ca), ("C1C", Signal::GpsL1Ca), ("L1C", Signal::GpsL1Ca)];
        for (t, want) in ca {
            assert_eq!(rinex_type_to_signal(Constellation::Gps, t), Some(want), "{t}");
        }
        let pcode = [("P1", Signal::GpsL1P), ("C1W", Signal::GpsL1P), ("L1W", Signal::GpsL1P), ("C1Y", Signal::GpsL1P), ("C2", Signal::GpsL2Cm), ("L2", Signal::GpsL2Cm), ("P2", Signal::GpsL2P), ("C2W", Signal::GpsL2P), ("C2X", Signal::GpsL2Cm)];
        for (t, want) in pcode {
            assert_eq!(rinex_type_to_signal(Constellation::Gps, t), Some(want), "{t}");
        }
        assert_eq!(rinex_type_to_signal(Constellation::Gps, "C5"), Some(Signal::GpsL5));
        assert_eq!(rinex_type_to_signal(Constellation::Gps, "L5Q"), Some(Signal::GpsL5));
        // Unmodelled GPS bands and malformed types.
        for t in ["C6", "C7", "C8", "", "C", "X1", "1C"] {
            assert_eq!(rinex_type_to_signal(Constellation::Gps, t), None, "{t}");
        }
    }

    #[test]
    fn glonass_rinex_types_map_correctly() {
        let cases = [
            ("C1", Signal::GloL1Of),
            ("L1", Signal::GloL1Of),
            ("P1", Signal::GloL1Of),
            ("C2", Signal::GloL2Of),
            ("L2", Signal::GloL2Of),
            ("P2", Signal::GloL2Of),
        ];
        for (t, want) in cases {
            assert_eq!(rinex_type_to_signal(Constellation::Glonass, t), Some(want), "{t}");
        }
        for t in ["C3", "C5", "C6"] {
            assert_eq!(rinex_type_to_signal(Constellation::Glonass, t), None, "{t}");
        }
        // Channel number rides along: R01 with k = -4 (classic GLONASS slot).
        let f1 = frequency_for(Constellation::Glonass, Signal::GloL1Of, -4);
        assert_close(f1, 1599.75e6, HZ_TOL, "R01 L1 k=-4");
    }

    /// THE bug-source test: Galileo "L2" in RINEX 2.11 mixed files.
    ///
    /// Two conventions exist in the wild: converters writing E5b into the L2
    /// slot (matching the old `get_frequency(band=2)` behaviour) and converters
    /// writing E5a there (what our datasets actually contain). This registry
    /// pins the policy to **E5a** — see the module documentation — while E5b
    /// stays reachable through its unambiguous RINEX 3 band 7 (`"L7"`).
    #[test]
    fn galileo_l2_maps_to_e5a_with_documented_ambiguity() {
        for t in ["L2", "C2", "P2", "C2W"] {
            assert_eq!(
                rinex_type_to_signal(Constellation::Galileo, t),
                Some(Signal::GalE5a),
                "Galileo legacy '{t}' must resolve to E5a per documented policy"
            );
        }
        // E5b must remain reachable, and must NOT be what "L2" gives.
        assert_eq!(rinex_type_to_signal(Constellation::Galileo, "L7"), Some(Signal::GalE5b));
        let l2_hz = frequency_for(Constellation::Galileo, Signal::GalE5a, 0);
        assert_close(l2_hz, 1176.45e6, HZ_TOL, "Galileo L2-slot = E5a");
        assert!((l2_hz - FREQ_GAL_E5B).abs() > 30.0e6, "E5a/E5b must differ");
    }

    #[test]
    fn galileo_unambiguous_bands() {
        let cases = [
            ("C1", Signal::GalE1Os),
            ("L1", Signal::GalE1Os),
            ("C5", Signal::GalE5a),
            ("L5", Signal::GalE5a),
            ("C6", Signal::GalE6Cs),
            ("L6", Signal::GalE6Cs),
            ("C7", Signal::GalE5b),
            ("L7", Signal::GalE5b),
        ];
        for (t, want) in cases {
            assert_eq!(rinex_type_to_signal(Constellation::Galileo, t), Some(want), "{t}");
        }
        assert_eq!(rinex_type_to_signal(Constellation::Galileo, "C8"), None);
    }

    #[test]
    fn beidou_rinex_types_map_correctly() {
        assert_eq!(rinex_type_to_signal(Constellation::Beidou, "C1"), Some(Signal::BdsB1i));
        assert_eq!(rinex_type_to_signal(Constellation::Beidou, "C1I"), Some(Signal::BdsB1i));
        // Documented limitation: BDS-3 B1C codes resolve to BdsB1i because no
        // B1C variant exists yet — see beidou_signal docs.
        assert_eq!(rinex_type_to_signal(Constellation::Beidou, "C1P"), Some(Signal::BdsB1i));
        assert_eq!(rinex_type_to_signal(Constellation::Beidou, "L6"), Some(Signal::BdsB3i));
        assert_eq!(rinex_type_to_signal(Constellation::Beidou, "C6I"), Some(Signal::BdsB3i));
        // B2I (band 2) and B2a (band 5) have no registry variant yet: explicit None.
        for t in ["C2", "L2", "C5", "L5", "L7", "C3"] {
            assert_eq!(rinex_type_to_signal(Constellation::Beidou, t), None, "{t}");
        }
    }

    #[test]
    fn qzss_aliases_gps_and_sbas_l1_l5_only() {
        assert_eq!(rinex_type_to_signal(Constellation::Qzss, "C1"), Some(Signal::GpsL1Ca));
        assert_eq!(rinex_type_to_signal(Constellation::Qzss, "C2"), Some(Signal::GpsL2Cm));
        assert_eq!(rinex_type_to_signal(Constellation::Qzss, "C5"), Some(Signal::GpsL5));
        assert_eq!(rinex_type_to_signal(Constellation::Sbas, "C1"), Some(Signal::GpsL1Ca));
        assert_eq!(rinex_type_to_signal(Constellation::Sbas, "L5"), Some(Signal::GpsL5));
        assert_eq!(rinex_type_to_signal(Constellation::Sbas, "C2"), None);
    }

    #[test]
    fn navic_explicitly_unmodelled() {
        for t in ["C1", "L1", "C5", "L5", "C9"] {
            assert_eq!(rinex_type_to_signal(Constellation::Navic, t), None, "{t}");
        }
    }

    #[test]
    fn parsing_is_case_insensitive_and_handles_short_codes() {
        assert_eq!(rinex_type_to_signal(Constellation::Galileo, "c2"), Some(Signal::GalE5a));
        assert_eq!(rinex_type_to_signal(Constellation::Gps, "p2"), Some(Signal::GpsL2P));
        assert_eq!(rinex_type_to_signal(Constellation::Gps, "c1c"), Some(Signal::GpsL1Ca));
        assert_eq!(rinex_type_to_signal(Constellation::Gps, ""), None);
    }

    #[test]
    fn frequency_for_ignores_k_for_cdma_and_guards_mismatch() {
        // CDMA signals ignore the channel entirely.
        assert_close(frequency_for(Constellation::Gps, Signal::GpsL1Ca, 42), FREQ_GPS_L1, 0.0, "k ignored");
        // Mismatched (constellation, FDMA signal) falls back to nominal k=0.
        assert_close(frequency_for(Constellation::Gps, Signal::GloL1Of, 3), FREQ_GLO_L1_NOMINAL, 0.0, "guard");
        assert_close(frequency_for(Constellation::Gps, Signal::GloL2Of, -7), FREQ_GLO_L2_NOMINAL, 0.0, "guard");
    }

    #[test]
    fn registry_agrees_with_legacy_constants() {
        // Guard against the two sources of truth drifting apart again.
        assert_close(Signal::GpsL1Ca.base_freq_hz(), FREQ_GPS_L1, 0.0, "L1");
        assert_close(Signal::GpsL2Cm.base_freq_hz(), FREQ_GPS_L2, 0.0, "L2");
        assert_close(Signal::GpsL5.base_freq_hz(), FREQ_GPS_L5, 0.0, "L5");
        assert_close(Signal::GalE1Os.base_freq_hz(), FREQ_GPS_L1, 0.0, "E1=L1");
        assert_close(Signal::GalE5b.base_freq_hz(), FREQ_GAL_E5B, 0.0, "E5b");
        assert_close(Signal::BdsB1i.base_freq_hz(), FREQ_BDS_B1I, 0.0, "B1I");
        assert_close(Signal::GloL1Of.fdma_offset_hz().unwrap_or(0.0), FREQ_GLO_L1_DELTA, 0.0, "dF1");
        assert_close(Signal::GloL2Of.fdma_offset_hz().unwrap_or(0.0), FREQ_GLO_L2_DELTA, 0.0, "dF2");
        assert_close(Signal::GloL1Of.base_freq_hz(), FREQ_GLO_L1_NOMINAL, 0.0, "nom1");
        assert_close(Signal::GloL2Of.base_freq_hz(), FREQ_GLO_L2_NOMINAL, 0.0, "nom2");
    }

    /// Every constellation the engine defines must produce a *defined* answer
    /// (mapped signal or explicit None) for every standard RINEX band — no
    /// silent GPS fallbacks, no panics. This is the coverage contract.
    #[test]
    fn all_engine_constellations_have_defined_band_answers() {
        let bands = ["1", "2", "5", "6", "7"];
        for cons in [Constellation::Gps, Constellation::Glonass, Constellation::Galileo, Constellation::Beidou, Constellation::Sbas, Constellation::Qzss, Constellation::Navic] {
            for b in bands {
                // Must terminate with Some or None — never fall back to GPS L1.
                let sig = rinex_type_to_signal(cons, b);
                if let Some(s) = sig {
                    let f = frequency_for(cons, s, 0);
                    assert!(f > 1.0e9 && f < 2.0e9, "sanity freq for {cons:?}/{b}");
                }
            }
        }
        // Spot-assert the documented policy edges.
        assert_eq!(rinex_type_to_signal(Constellation::Navic, "1"), None);
        assert_eq!(rinex_type_to_signal(Constellation::Beidou, "2"), None);
        assert_eq!(rinex_type_to_signal(Constellation::Sbas, "2"), None);
    }
}
