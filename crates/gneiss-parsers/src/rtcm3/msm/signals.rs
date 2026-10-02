//! MSM signal index to RINEX observation code mappings and conversion scale constants.

use gneiss_core::constants::SPEED_OF_LIGHT_M_S;
use gneiss_core::sat::Constellation;

pub(crate) const RANGE_MS: f64 = SPEED_OF_LIGHT_M_S * 0.001;
pub(crate) const P2_10: f64 = 1.0 / 1024.0;
pub(crate) const P2_24: f64 = 1.0 / 16_777_216.0;
pub(crate) const P2_29: f64 = 1.0 / 536_870_912.0;
pub(crate) const P2_31: f64 = 1.0 / 2_147_483_648.0;

pub(crate) fn msm_signal_rinex_code(c: Constellation, sig_id: u8) -> Option<&'static str> {
    const GPS: [&str; 32] = [
        "", "1C", "1P", "1W", "1Y", "1M", "", "2C", "2P", "2W", "2Y", "2M",
        "", "", "2S", "2L", "2X", "", "", "", "", "5I", "5Q", "5X",
        "", "", "", "", "", "1S", "1L", "1X",
    ];
    const GAL: [&str; 32] = [
        "", "1C", "1A", "1B", "1X", "1Z", "", "6C", "6A", "6B", "6X", "6Z",
        "", "7I", "7Q", "7X", "", "8I", "8Q", "8X", "", "5I", "5Q", "5X",
        "", "", "", "", "", "", "", "",
    ];
    const BDS: [&str; 32] = [
        "", "1I", "1Q", "1X", "", "", "", "6I", "6Q", "6X", "", "",
        "", "7I", "7Q", "7X", "", "", "", "", "", "", "", "",
        "", "", "", "", "", "", "", "",
    ];
    const GLO: [&str; 32] = [
        "", "1C", "1P", "", "", "", "", "2C", "2P", "", "3I", "3Q",
        "3X", "", "", "", "", "", "", "", "", "", "", "",
        "", "", "", "", "", "", "", "",
    ];
    const QZS: [&str; 32] = [
        "", "1C", "", "", "", "", "", "", "6S", "6L", "6X", "",
        "", "", "2S", "2L", "2X", "", "", "", "", "5I", "5Q", "5X",
        "", "", "", "", "", "1S", "1L", "1X",
    ];
    let table = match c {
        Constellation::Gps => &GPS,
        Constellation::Galileo => &GAL,
        Constellation::Beidou => &BDS,
        Constellation::Glonass => &GLO,
        Constellation::Qzss => &QZS,
        _ => return None,
    };
    let idx = sig_id as usize;
    if idx == 0 || idx > 31 {
        return None;
    }
    let code = table[idx];
    if code.is_empty() {
        None
    } else {
        Some(code)
    }
}

pub(crate) fn split_rinex_code(code: &str) -> Option<(u8, char)> {
    let mut chars = code.chars();
    let band = chars.next()?.to_digit(10)? as u8;
    let attribute = chars.next()?;
    Some((band, attribute))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    /// Transcribes the MSM signal-and-tracking-mode tables from the reference
    /// decoder RTKLIB `rtcm3.c` (`msm_sig_gps`, `msm_sig_glo`, `msm_sig_gal`,
    /// `msm_sig_qzs`, lines 64-87), whose entries are annotated `/* 1-12 */` ..
    /// `/* 25-32 */` so that array index equals the RTCM signal and tracking
    /// mode number and index 0 is unused. Signal and tracking mode 1 is
    /// "GPS L1 C/A", so index 1 must be `1C`.
    const REF_GPS: [&str; 32] = [
        "", "1C", "1P", "1W", "1Y", "1M", "", "2C", "2P", "2W", "2Y", "2M",
        "", "", "2S", "2L", "2X", "", "", "", "", "5I", "5Q", "5X",
        "", "", "", "", "", "1S", "1L", "1X",
    ];
    const REF_GLO: [&str; 32] = [
        "", "1C", "1P", "", "", "", "", "2C", "2P", "", "3I", "3Q",
        "3X", "", "", "", "", "", "", "", "", "", "", "",
        "", "", "", "", "", "", "", "",
    ];
    const REF_GAL: [&str; 32] = [
        "", "1C", "1A", "1B", "1X", "1Z", "", "6C", "6A", "6B", "6X", "6Z",
        "", "7I", "7Q", "7X", "", "8I", "8Q", "8X", "", "5I", "5Q", "5X",
        "", "", "", "", "", "", "", "",
    ];
    const REF_QZS: [&str; 32] = [
        "", "1C", "", "", "", "", "", "", "6S", "6L", "6X", "",
        "", "", "2S", "2L", "2X", "", "", "", "", "5I", "5Q", "5X",
        "", "", "", "", "", "1S", "1L", "1X",
    ];

    fn check(sys: Constellation, reference: &[&str; 32]) {
        for (id, expected) in reference.iter().enumerate() {
            assert_eq!(
                msm_signal_rinex_code(sys, id as u8),
                if expected.is_empty() { None } else { Some(*expected) },
                "{sys:?} MSM signal and tracking mode {id}"
            );
        }
    }

    #[test]
    fn gps_tables_match_the_reference() {
        check(Constellation::Gps, &REF_GPS);
    }

    #[test]
    fn glonass_tables_match_the_reference() {
        check(Constellation::Glonass, &REF_GLO);
    }

    #[test]
    fn galileo_tables_match_the_reference() {
        check(Constellation::Galileo, &REF_GAL);
    }

    #[test]
    fn qzss_tables_match_the_reference() {
        check(Constellation::Qzss, &REF_QZS);
    }

    /// Signal id 0 and 32+ are out of range for every constellation.
    #[test]
    fn out_of_range_signal_ids_are_rejected() {
        for id in [0u8, 32, 33, 255] {
            assert_eq!(msm_signal_rinex_code(Constellation::Gps, id), None);
            assert_eq!(msm_signal_rinex_code(Constellation::Galileo, id), None);
        }
    }

    /// Every constellation outside the supported set must return None rather
    /// than fall through to another system's table.
    #[test]
    fn unsupported_constellations_have_no_mapping() {
        for sys in [Constellation::Sbas, Constellation::Navic] {
            assert_eq!(msm_signal_rinex_code(sys, 1), None);
        }
    }

    #[test]
    fn split_rinex_code_separates_band_and_attribute() {
        assert_eq!(split_rinex_code("1C"), Some((1, 'C')));
        assert_eq!(split_rinex_code("5X"), Some((5, 'X')));
        // Zero-length and one-character codes cannot yield both parts.
        assert_eq!(split_rinex_code(""), None);
        assert_eq!(split_rinex_code("1"), None);
        // A non-numeric band is not a valid RINEX observation code.
        assert_eq!(split_rinex_code("AC"), None);
    }

    /// The band handed to `push_carrier_phase_obs` -> `track_c_frequency` is
    /// parsed out of the MSM code string, not the raw MSM signal id
    /// (decoder.rs:155-156). The BeiDou table yields only bands `{1, 6, 7}` —
    /// it has no band-2 entry at all.
    ///
    /// This is what makes the BeiDou band-2 frequency change inert on the RTCM
    /// MSM path: no BeiDou MSM message can ask for band 2, so the frequency
    /// `track_c_frequency` returns there cannot change. Exhaustive over all 32
    /// signal ids, not a hand-picked subset.
    #[test]
    fn beidou_msm_table_can_never_select_band_2() {
        let bands: BTreeSet<u8> = (0u8..32)
            .filter_map(|id| msm_signal_rinex_code(Constellation::Beidou, id))
            .filter_map(|code| split_rinex_code(code).map(|(band, _)| band))
            .collect();
        assert!(
            !bands.contains(&2),
            "MSM must never reach BeiDou band 2; a new code entry would do it"
        );
        assert_eq!(bands, BTreeSet::from([1u8, 6, 7]));
    }
}
