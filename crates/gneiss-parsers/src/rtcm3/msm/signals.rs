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
