//! Explicit signal–frequency registry.
//!
//! One canonical mapping between *physical GNSS signals* and their *transmitted
//! frequencies*, replacing ad-hoc `(constellation, band_number)` lookups where
//! a single band number means different physical signals per constellation
//! (GPS band 2 = L2 @ 1227.60 MHz, Galileo band 2 = ambiguous legacy slot).
//!
//! References: GPS ICD-IS-200/705, GLONASS ICD (5.1 ed.), Galileo ICD OS-SIS,
//! BeiDou ICD B1I/B3I (all centre frequencies are ITU-allocated RNSS carriers).
//!
//! # Integration status
//!
//! The legacy [`crate::signal`] lookups remain wired into the estimators; this
//! registry is the reference they must converge onto. Do not migrate only one
//! call site while others keep legacy band-2 semantics — that recreates the
//! very cross-path inconsistency this module exists to remove.

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
    BdsB2i,
    BdsB3i,
}

impl Signal {
    /// Base frequency in Hz (without FDMA offset for GLONASS).
    pub fn base_freq_hz(&self) -> f64 {
        match self {
            Signal::GpsL1Ca | Signal::GpsL1P | Signal::GalE1Os => 1575.42e6,
            Signal::GpsL2Cm | Signal::GpsL2P => 1227.60e6,
            Signal::GpsL5 | Signal::GalE5a => 1176.45e6,
            Signal::GloL1Of => 1602.0e6,
            Signal::GloL2Of => 1246.0e6,
            Signal::GalE5b => 1207.14e6,
            Signal::GalE6Cs => 1278.75e6,
            Signal::BdsB1i => 1561.098e6,
            Signal::BdsB2i => 1207.14e6,
            Signal::BdsB3i => 1268.52e6,
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

/// One character class of a RINEX observable type string (`"C1"`, `"L2W"`, `"P1"`…).
struct RinexCode {
    kind: char,
    band: u8,
    attr: char,
}

impl RinexCode {
    /// Parses 2-char RINEX 2 types (`"L1"`, `"P2"`) and 3-char RINEX 3 types
    /// (`"C1C"`, `"L2W"`). Case-insensitive; missing attribute becomes `' '`;
    /// characters beyond the third are ignored. Only valid observation kinds
    /// (`C`/`L`/`P`/`D`/`S`) parse.
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

/// Map from RINEX observable type + constellation to [`Signal`].
///
/// Accepts RINEX 2 (`"L1"`, `"P2"`) and RINEX 3 (`"C1C"`, `"L2W"`) spellings.
///
/// # The Galileo band-2 ambiguity (pinned here, once)
///
/// RINEX 3 makes Galileo bands unambiguous: 1 = E1, 5 = E5a, 6 = E6, 7 = E5b.
/// RINEX 2.11 mixed GNSS exports instead reuse the GPS letter scheme, and the
/// second-frequency slot labelled `"L2"/"C2"` for Galileo satellites is written
/// **differently by different converters**: some put E5a (1176.45 MHz) there,
/// others E5b (1207.14 MHz). Per the project's dataset directive this registry
/// resolves Galileo `"L2"` to [`Signal::GalE5a`], which deliberately *differs*
/// from the legacy `crate::signal::get_frequency(.., band 2, ..)` (= E5b);
/// before migrating a call site, confirm which convention your files use.
/// Unambiguous RINEX 3 codes (band 5 = E5a, band 7 = E5b) are always preferred.
///
/// # BeiDou band numbers are version-dependent
///
/// In RINEX 3, B1I rides band 2 (`C2I`, confirmed by MGEX headers) while band
/// 1 carries BDS-3 B1C @1575.42 MHz (`C1P/C1X`). B1C has no variant here yet,
/// so band-1-with-B1C-attribute resolves to `None` rather than risking a
/// 14.3 MHz error. Legacy RINEX 2 files put B1I in the plain band-1 slot.
///
/// Returns `None` for combinations with no modelled signal (NavIC, BDS
/// B1C/B2a/B2b, Galileo E5 AltBOC band 8) — callers must skip such observables.
pub fn rinex_type_to_signal(constellation: Constellation, rinex_type: &str) -> Option<Signal> {
    let code = RinexCode::parse(rinex_type)?;
    match constellation {
        Constellation::Gps | Constellation::Qzss => gps_signal(&code),
        Constellation::Sbas => sbas_signal(code.band),
        Constellation::Glonass => glonass_signal(code.band),
        Constellation::Galileo => galileo_signal(code.band),
        Constellation::Beidou => beidou_signal(&code),
        Constellation::Navic => None, // L5/S-band not modelled yet
    }
}

/// Actual transmitted frequency in Hz for a specific satellite.
///
/// `freq_num` is the GLONASS FDMA channel number `k` (ICD operational range
/// −7..+6; out-of-range values pass through arithmetically — callers must
/// validate); it is ignored for CDMA signals. If `sat` disagrees with the
/// signal's own constellation (only possible for the FDMA branch), the nominal
/// k = 0 channel is returned rather than panicking; treat that combination as
/// an upstream bug.
pub fn frequency_for(sat: Constellation, signal: Signal, freq_num: i8) -> f64 {
    match signal.fdma_offset_hz() {
        None => signal.base_freq_hz(),
        Some(offset) => {
            let k = if sat == Constellation::Glonass {
                freq_num
            } else {
                0
            };
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
        1 => Some(if code.is_precise() {
            Signal::GpsL1P
        } else {
            Signal::GpsL1Ca
        }),
        2 => Some(if code.is_precise() {
            Signal::GpsL2P
        } else {
            Signal::GpsL2Cm
        }),
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

/// BeiDou — band numbers are RINEX-version-dependent, resolved via attributes:
/// RINEX 3 band 2 (`C2I`) carries B1I @1561.098; RINEX 3 band 1 carries BDS-3
/// B1C @1575.42 (`C1P/C1X`) which has NO variant here, so it maps to None
/// rather than risking a 14.3 MHz error. Legacy RINEX 2 files put B1I in the
/// plain band-1 slot (`"C1"`), kept as [`Signal::BdsB1i`]. Band 6 = B3I.
/// Bands 5 (B2a) and 7 (B2b) have no variant in this registry version yet.
fn beidou_signal(code: &RinexCode) -> Option<Signal> {
    let b1c_attr = matches!(code.attr, 'P' | 'X' | 'D');
    match code.band {
        1 if !b1c_attr => Some(Signal::BdsB1i), // RINEX 2 legacy slot only
        2 => Some(Signal::BdsB1i),              // RINEX 3 B1I (e.g. "C2I")
        6 => Some(Signal::BdsB3i),
        7 => Some(Signal::BdsB2i),              // RINEX 3 B2I (e.g. "C7I")
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


/// Engine policy for the secondary (dual-frequency) band per constellation.
///
/// Resolves the historical E5a/E5b ambiguity at the type level: callers
/// receive an explicit `(rinex_band, Signal)` pair instead of guessing.
/// GPS/GLONASS use RINEX band 2; Galileo exports carry E5a on band 5.
pub fn secondary_signal(c: crate::sat::Constellation) -> Option<(u8, Signal)> {
    use crate::sat::Constellation;
    match c {
        Constellation::Gps | Constellation::Qzss => Some((2, Signal::GpsL2Cm)),
        Constellation::Glonass => Some((2, Signal::GloL2Of)),
        Constellation::Galileo => Some((5, Signal::GalE5a)),
        Constellation::Beidou => Some((7, Signal::BdsB2i)),
        _ => None,
    }
}

/// Primary-band (band 1) signal for a constellation.
pub fn primary_signal(c: crate::sat::Constellation) -> Option<Signal> {
    use crate::sat::Constellation;
    match c {
        Constellation::Gps | Constellation::Qzss => Some(Signal::GpsL1Ca),
        Constellation::Glonass => Some(Signal::GloL1Of),
        Constellation::Galileo => Some(Signal::GalE1Os),
        Constellation::Beidou => Some(Signal::BdsB1i),
        _ => None,
    }
}


/// Resolve a RINEX band number to its Signal for a constellation.
///
/// This is the authoritative band→signal mapping; callers reading phase
/// from band N MUST pair it with this signal's wavelength, independent
/// of any preferred-secondary policy.
pub fn signal_for_band(c: crate::sat::Constellation, band: u8) -> Option<Signal> {
    use crate::sat::Constellation;
    match (c, band) {
        (Constellation::Gps | Constellation::Qzss, 1) => Some(Signal::GpsL1Ca),
        (Constellation::Gps | Constellation::Qzss, 2) => Some(Signal::GpsL2Cm),
        (Constellation::Gps | Constellation::Qzss, 5) => Some(Signal::GpsL5),
        (Constellation::Glonass, 1) => Some(Signal::GloL1Of),
        (Constellation::Glonass, 2) => Some(Signal::GloL2Of),
        (Constellation::Galileo, 1) => Some(Signal::GalE1Os),
        (Constellation::Galileo, 5) => Some(Signal::GalE5a),
        // RINEX 3 Galileo band numbering: 5=E5a, 6=E6, 7=E5b.
        // These two were swapped once already (ledger row 11); pinned by
        // frequency_parity tests across ALL observed bands.
        (Constellation::Galileo, 6) => Some(Signal::GalE6Cs),
        (Constellation::Galileo, 7) => Some(Signal::GalE5b),
        (Constellation::Beidou, 1) => Some(Signal::BdsB1i),
        (Constellation::Beidou, 2) => Some(Signal::BdsB1i),
        (Constellation::Beidou, 5) => Some(Signal::BdsB3i),
        (Constellation::Beidou, 7) => Some(Signal::BdsB2i),
        _ => None,
    }
}

/// Frequency for a constellation's band through the Track C registry
/// (`signal_for_band` + `frequency_for`), falling back to the legacy
/// per-band table (`crate::signal::get_frequency`) for the
/// constellation/band combinations `signal_for_band` doesn't cover yet.
/// See this module's doc comment: both paths must agree wherever they
/// overlap, so the fallback is deliberate, not a workaround to remove.
pub fn track_c_frequency(c: Constellation, band: u8, glo_k: i8) -> f64 {
    match signal_for_band(c, band) {
        Some(sig) => frequency_for(c, sig, glo_k),
        None => crate::signal::get_frequency(
            crate::sat::SatelliteId { constellation: c, prn: 0 },
            band,
            glo_k,
        ),
    }
}


#[cfg(test)]
mod tests;
