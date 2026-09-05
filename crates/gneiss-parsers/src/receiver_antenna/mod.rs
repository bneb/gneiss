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
    /// Azimuth grid step (deg), 0.0 if NOAZI only.
    pub dazi_deg: f64,
    /// All NOAZI grids keyed by ANTEX frequency code (`"G01"`, `"R02"`,
    /// ...) so constellations outside GPS L1/L2 stay reachable through
    /// [`ReceiverAntenna::pcv_mm`].
    grids: HashMap<String, Vec<f64>>,
    /// Azimuth grids keyed by ANTEX frequency code.
    azi_grids: HashMap<String, Vec<Vec<f64>>>,
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

    /// PCV in millimetres at azimuth and zenith angles (degrees) for an explicit ANTEX
    /// frequency code (`"G01"`, `"R02"`, ...).
    /// Uses 2D bilinear interpolation if an azimuth grid is present;
    /// otherwise falls back to 1D NOAZI interpolation.
    pub fn pcv_mm_az_zen(&self, freq_code: &str, az_deg: f64, zen_deg: f64) -> Option<f64> {
        if let Some(azi) = self.azi_grids.get(freq_code) {
            if self.dazi_deg > 0.0 && !azi.is_empty() && self.dzen_deg > 0.0 {
                return Some(interp_2d_grid(
                    azi,
                    self.dazi_deg,
                    self.zen1_deg,
                    self.dzen_deg,
                    az_deg,
                    zen_deg,
                ));
            }
        }
        self.pcv_mm(freq_code, zen_deg)
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
    let mut grids = HashMap::new();
    let mut azi_grids = HashMap::new();
    for (code, freq) in &model.frequencies {
        grids.insert(code.clone(), freq.noazi.clone());
        if let Some(ref azi) = freq.azi {
            azi_grids.insert(code.clone(), azi.clone());
        }
    }
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
        dazi_deg: model.dazi,
        grids,
        azi_grids,
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

/// 2D bilinear interpolation of an azimuth-zenith grid with 360-deg azimuth wrapping.
fn interp_2d_grid(
    azi: &[Vec<f64>],
    dazi_deg: f64,
    zen_start_deg: f64,
    zen_step_deg: f64,
    az_deg: f64,
    zen_deg: f64,
) -> f64 {
    let az_norm = az_deg.rem_euclid(360.0);
    let n_azi = azi.len();
    let n_zen = azi[0].len();
    if n_zen == 0 || n_azi == 0 {
        return 0.0;
    }

    let az_pos = az_norm / dazi_deg;
    let az_i = (az_pos.floor() as usize) % n_azi;
    let az_next = (az_i + 1) % n_azi;
    let az_frac = az_pos - az_pos.floor();

    let zen_last = (n_zen - 1) as f64;
    let zen_pos = ((zen_deg - zen_start_deg) / zen_step_deg).clamp(0.0, zen_last);
    let zen_i = zen_pos.floor() as usize;
    let zen_next = (zen_i + 1).min(n_zen - 1);
    let zen_frac = zen_pos - zen_i as f64;

    let v00 = azi[az_i][zen_i];
    let v01 = azi[az_i][zen_next];
    let v10 = azi[az_next][zen_i];
    let v11 = azi[az_next][zen_next];

    let v0 = v00 * (1.0 - zen_frac) + v01 * zen_frac;
    let v1 = v10 * (1.0 - zen_frac) + v11 * zen_frac;

    v0 * (1.0 - az_frac) + v1 * az_frac
}

/// Differential receiver PCV with 2D azimuth and zenith interpolation, metres.
///
/// Sky azimuth and elevation angles are given in radians (local topocentric frame).
/// Rover antenna azimuth is rotated by `rover_heading_rad` relative to North;
/// base antenna is assumed North-aligned (heading 0).
#[allow(clippy::too_many_arguments)]
pub fn compute_dd_pcv_correction_2d(
    rover_ant: &ReceiverAntenna,
    base_ant: &ReceiverAntenna,
    freq_code: &str,
    az_sat_rad: f64,
    el_sat_rad: f64,
    az_ref_rad: f64,
    el_ref_rad: f64,
    rover_heading_rad: f64,
) -> f64 {
    let zen_sat = 90.0 - el_sat_rad.to_degrees();
    let zen_ref = 90.0 - el_ref_rad.to_degrees();
    let az_sat_rov = (az_sat_rad - rover_heading_rad).to_degrees();
    let az_ref_rov = (az_ref_rad - rover_heading_rad).to_degrees();
    let az_sat_bas = az_sat_rad.to_degrees();
    let az_ref_bas = az_ref_rad.to_degrees();
    let (Some(rov_s), Some(rov_r)) = (
        rover_ant.pcv_mm_az_zen(freq_code, az_sat_rov, zen_sat),
        rover_ant.pcv_mm_az_zen(freq_code, az_ref_rov, zen_ref),
    ) else {
        return 0.0;
    };
    let (Some(bas_s), Some(bas_r)) = (
        base_ant.pcv_mm_az_zen(freq_code, az_sat_bas, zen_sat),
        base_ant.pcv_mm_az_zen(freq_code, az_ref_bas, zen_ref),
    ) else {
        return 0.0;
    };
    ((rov_s - rov_r) - (bas_s - bas_r)) / 1000.0
}

/// Differential receiver PCV embedded in one DD phase observation, metres.
///
/// Retained for 1D elevation-only callers using NOAZI.
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
