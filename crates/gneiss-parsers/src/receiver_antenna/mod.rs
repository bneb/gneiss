//! Receiver antenna PCV models from an ANTEX database.
//!
//! Satellite blocks in ANTEX carry the PRN in the serial field; receiver
//! blocks carry the antenna family in field 1 and the radome in field 2 of
//! `TYPE / SERIAL NO` (e.g. `TRM59800.80     SCIT`). Both parse into
//! [`AntennaPcv`] identically; this module adds the receiver-side lookup,
//! public PCO/PCV fields, zenith-angle interpolation of the NOAZI grid,
//! and the double-difference PCV correction consumed by the DD pipeline.
//!
//! Grid convention: PCV node `i` sits at zenith angle
//! `zen1_deg + i * dzen_deg`; values are millimetres. Interpolation is
//! linear between neighbouring nodes and clamps outside the tabulated
//! range (the elevation mask keeps normal operation well inside it).
//!
//! Sign convention (mirrors phase windup): the observed carrier phase
//! *contains* the antenna signature, `phi_obs = rho/lam + N + phi_PCV`, so
//! the correction returned by [`compute_dd_pcv_correction`] is
//! **subtracted** from the DD phase observation.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use gneiss_core::sat::Constellation;

use crate::antex::{AntennaPcv, AntexDatabase};

/// Read a station's antenna family + radome from a RINEX (2 or 3)
/// observation header's `ANT # / TYPE` line, for ANTEX lookup.
///
/// Returns `None` if the file can't be opened or has no such line in its
/// first 80 lines (the header always precedes any observation record).
pub fn rinex_ant_type(rinex_path: &Path) -> Option<(String, String)> {
    let f = std::fs::File::open(rinex_path).ok()?;
    for line in std::io::BufReader::new(f).lines().take(80).map_while(Result::ok) {
        if line.len() >= 60 && line[60..].trim() == "ANT # / TYPE" {
            let fields: Vec<&str> = line[20..40].split_whitespace().collect();
            let fam = fields.first()?.to_string();
            let rad = fields.get(1).copied().unwrap_or("NONE").to_string();
            return Some((fam, rad));
        }
    }
    None
}

/// One receiver-antenna calibration extracted from an [`AntexDatabase`].
#[derive(Debug, Clone)]
pub struct ReceiverAntenna {
    /// Antenna family, e.g. `"TRM59800.80"`.
    pub ant_type: String,
    /// Radome type (`"SCIT"`), or `"NONE"` for unradomed calibrations.
    pub radome: String,
    /// L1 phase centre offset in millimetres, `[east, north, up]`.
    /// ANTEX tabulates north/east/up; the components are reordered here.
    pub pco_enu_mm: [f64; 3],
    /// L1 NOAZI PCV values (mm) over the zenith grid.
    pub pcv_l1_grid: Vec<f64>,
    /// L2 NOAZI PCV values (mm) over the zenith grid.
    pub pcv_l2_grid: Vec<f64>,
    /// Zenith angle (deg) of the first grid node.
    pub zen1_deg: f64,
    /// Zenith grid step (deg).
    pub dzen_deg: f64,
    /// All NOAZI grids keyed by ANTEX frequency code (`"G01"`, `"R02"`,
    /// ...) so constellations outside GPS L1/L2 stay reachable through
    /// [`ReceiverAntenna::pcv_mm`].
    grids: HashMap<String, Vec<f64>>,
}

impl ReceiverAntenna {
    /// Load a receiver calibration by antenna family and radome.
    ///
    /// Matches the whitespace-split TYPE / SERIAL NO field, so both
    /// `"TRM59800.80     SCIT"` (family + radome) and single-token entries
    /// like `"AOAD/M_T"` (radome `NONE`) resolve. Returns `None` when the
    /// database has no such receiver entry.
    ///
    /// The L1/L2 band grids prefer the GPS codes `G01`/`G02` and fall back
    /// to the GLONASS k=0 codes `R01`/`R02` (the convention IGS ANTEX
    /// files use for receiver GLONASS tables).
    pub fn from_antex(db: &AntexDatabase, ant_type: &str, radome: &str) -> Option<Self> {
        let model = find_block(db, ant_type, radome)?;
        Some(from_model(model))
    }

    /// Alias of [`ReceiverAntenna::from_antex`] retained for engine call
    /// sites that predate the public-field API.
    pub fn lookup(db: &AntexDatabase, family: &str, radome: &str) -> Option<Self> {
        Self::from_antex(db, family, radome)
    }

    /// Display name, e.g. `"TRM59800.80 SCIT"`.
    pub fn antenna_type(&self) -> String {
        format!("{} {}", self.ant_type, self.radome)
    }

    /// NOAZI PCV in millimetres at a zenith angle in degrees for a
    /// frequency band (`1` = L1 grid, `2` = L2 grid, anything else is
    /// `None`). Clamps outside the tabulated range.
    pub fn pcv_mm_at_zenith(&self, freq_band: u8, zenith_deg: f64) -> Option<f64> {
        let grid = match freq_band {
            1 => &self.pcv_l1_grid,
            2 => &self.pcv_l2_grid,
            _ => return None,
        };
        interp_grid(grid, self.zen1_deg, self.dzen_deg, zenith_deg)
    }


    /// NOAZI PCV in millimetres at a zenith angle for an explicit ANTEX
    /// frequency code (`"G01"`, `"R02"`, ...). `None` when the calibration
    /// has no table for that code.
    pub fn pcv_mm(&self, freq_code: &str, zenith_deg: f64) -> Option<f64> {
        interp_grid(self.grids.get(freq_code)?, self.zen1_deg, self.dzen_deg, zenith_deg)
    }
}

/// Find a receiver calibration by family and radome.
fn find_block<'a>(db: &'a AntexDatabase, family: &str, radome: &str) -> Option<&'a AntennaPcv> {
    db.antennas.iter().find(|a| {
        let tokens: Vec<&str> = a.antenna_type.split_whitespace().collect();
        match tokens.as_slice() {
            [t, r] => *t == family && *r == radome,
            [t] => *t == family && radome == "NONE",
            _ => false,
        }
    })
}

/// Extract the public calibration fields from one ANTEX block.
///
/// Band grids prefer GPS `G01`/`G02` and fall back to the GLONASS k=0
/// codes `R01`/`R02`; the PCO comes from the L1 frequency. ANTEX stores
/// the PCO as north/east/up, reordered here to east/north/up.
fn from_model(model: &AntennaPcv) -> ReceiverAntenna {
    let noazi = |codes: &[&str]| {
        codes
            .iter()
            .find_map(|c| model.frequencies.get(*c))
            .map(|f| f.noazi.clone())
            .unwrap_or_default()
    };
    let l1 = model.frequencies.get("G01").or_else(|| model.frequencies.get("R01"));
    let pco_enu_mm = l1.map_or([0.0; 3], |f| [f.pco.y, f.pco.x, f.pco.z]);
    let grids = model
        .frequencies
        .iter()
        .map(|(code, freq)| (code.clone(), freq.noazi.clone()))
        .collect();
    let tokens: Vec<&str> = model.antenna_type.split_whitespace().collect();
    let (ant_type, radome) = match tokens.as_slice() {
        [t] => ((*t).to_string(), "NONE".to_string()),
        [t, r] => ((*t).to_string(), (*r).to_string()),
        _ => (model.antenna_type.clone(), String::new()),
    };
    ReceiverAntenna {
        ant_type,
        radome,
        pco_enu_mm,
        pcv_l1_grid: noazi(&["G01", "R01"]),
        pcv_l2_grid: noazi(&["G02", "R02"]),
        zen1_deg: model.zen1,
        dzen_deg: model.dzen,
        grids,
    }
}

/// Map a pipeline (constellation, frequency band, GLONASS channel) to its
/// ANTEX frequency code. Receiver GLONASS tables in IGS14 are labelled with
/// the k=0 channel (`R01`/`R02`) regardless of the satellite's actual k;
/// BeiDou calibrations are not consumed by this engine yet.
pub fn frequency_code(constellation: Constellation, band: u8) -> Option<String> {
    let prefix = match constellation {
        Constellation::Gps | Constellation::Sbas => 'G',
        Constellation::Qzss => 'J',
        Constellation::Glonass => 'R',
        Constellation::Galileo => 'E',
        Constellation::Beidou | Constellation::Navic => return None,
    };
    if !(1..=9).contains(&band) {
        return None;
    }
    Some(format!("{prefix}{band:0>2}"))
}

/// Linearly interpolate a NOAZI PCV grid (mm) at `zenith_deg`.
///
/// Node `i` sits at `zen1_deg + i*dzen_deg`; clamps outside the range.
fn interp_grid(grid: &[f64], zen1_deg: f64, dzen_deg: f64, zenith_deg: f64) -> Option<f64> {
    if grid.is_empty() || dzen_deg <= 0.0 {
        return None;
    }
    let last = grid.len() - 1;
    let t = (zenith_deg - zen1_deg) / dzen_deg;
    if t <= 0.0 {
        return Some(grid[0]);
    }
    if t >= last as f64 {
        return Some(grid[last]);
    }
    let i = t.floor() as usize;
    Some(grid[i] + (t - i as f64) * (grid[i + 1] - grid[i]))
}

/// Differential receiver PCV embedded in one DD phase observation, metres.
///
/// For a DD formed between satellites at elevations `el_sat`/`el_ref`
/// (radians, rover frame) with different rover/base antennas:
/// `dd_phi contains ([PCV_rov(z_s) - PCV_rov(z_r)] - [PCV_base(z_s) -
/// PCV_base(z_r)])` where `z = 90 deg - el`. Same-family pairs give ~0.
/// The caller subtracts this from the DD phase cycles (`cp - corr/lambda`),
/// matching the windup convention. Returns 0.0 when either calibration
/// lacks the requested frequency: a one-sided correction would be worse
/// than none.
pub fn compute_dd_pcv_correction(
    rover_ant: &ReceiverAntenna,
    base_ant: &ReceiverAntenna,
    freq_code: &str,
    el_sat: f64,
    el_ref: f64,
) -> f64 {
    let zen_sat = 90.0 - el_sat.to_degrees();
    let zen_ref = 90.0 - el_ref.to_degrees();
    let (Some(rov_s), Some(rov_r)) = (
        rover_ant.pcv_mm(freq_code, zen_sat),
        rover_ant.pcv_mm(freq_code, zen_ref),
    ) else {
        return 0.0;
    };
    let (Some(bas_s), Some(bas_r)) = (
        base_ant.pcv_mm(freq_code, zen_sat),
        base_ant.pcv_mm(freq_code, zen_ref),
    ) else {
        return 0.0;
    };
    ((rov_s - rov_r) - (bas_s - bas_r)) / 1000.0
}


#[cfg(test)]
mod tests;
