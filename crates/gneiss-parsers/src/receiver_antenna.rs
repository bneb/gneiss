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
mod tests {
    use super::*;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU32, Ordering};

    // ------------------------------------------------------------------
    // Fixtures
    // ------------------------------------------------------------------

    fn ant_line(data: &str, label: &str) -> String {
        let mut line = String::from(data);
        while line.len() < 60 {
            line.push(' ');
        }
        line.push_str(label);
        line.push('\n');
        line
    }

    static DB_COUNTER: AtomicU32 = AtomicU32::new(0);

    fn db_from(content: &str) -> AntexDatabase {
        let n = DB_COUNTER.fetch_add(1, Ordering::SeqCst);
        let dir = std::env::temp_dir();
        let path = dir.join(format!("test_recv_pcv_{}_{n}.atx", std::process::id()));
        std::fs::write(&path, content).unwrap();
        let db = AntexDatabase::parse(&path).unwrap();
        std::fs::remove_file(&path).ok();
        db
    }

    /// Receiver-style block whose NOAZI row is a linear ramp
    /// `slope * zen` over `0..=18 deg` step 1, with optional extra
    /// `(code, slope)` frequency blocks after the G01 one.
    fn ramp_block(family: &str, radome: &str, g01_slope: f64, extra: &[(&str, f64)]) -> String {
        let type_field = if radome == "NONE" {
            family.to_string()
        } else {
            format!("{family} {radome}")
        };
        let mut b = String::new();
        b.push_str(&ant_line("", "START OF ANTENNA"));
        b.push_str(&ant_line(&format!("{type_field:<40}"), "TYPE / SERIAL NO"));
        b.push_str(&ant_line("     0.0", "DAZI"));
        b.push_str(&ant_line("     0.0  18.0   1.0", "ZEN1 / ZEN2 / DZEN"));
        for (code, slope) in [("G01", g01_slope)].into_iter().chain(extra.iter().copied()) {
            b.push_str(&ant_line(&format!("   {code}"), "START OF FREQUENCY"));
            b.push_str(&ant_line("      0.00      0.00      5.00", "NORTH / EAST / UP"));
            let row: Vec<String> = (0..=18).map(|i| format!("{:>6.2}", slope * i as f64)).collect();
            b.push_str(&format!("   NOAZI{}\n", row.join(" ")));
            b.push_str(&ant_line("", "END OF FREQUENCY"));
        }
        b.push_str(&ant_line("", "END OF ANTENNA"));
        b
    }

    /// Block with distinguishable PCO columns (N=1, E=2, U=3 mm).
    fn pco_block(family: &str) -> String {
        let mut b = String::new();
        b.push_str(&ant_line("", "START OF ANTENNA"));
        b.push_str(&ant_line(&format!("{family:<40}"), "TYPE / SERIAL NO"));
        b.push_str(&ant_line("     0.0", "DAZI"));
        b.push_str(&ant_line("     0.0  18.0   1.0", "ZEN1 / ZEN2 / DZEN"));
        b.push_str(&ant_line("   G01", "START OF FREQUENCY"));
        b.push_str(&ant_line("      1.00      2.00      3.00", "NORTH / EAST / UP"));
        b.push_str(&ant_line("   NOAZI   0.00  -1.00", ""));
        b.push_str(&ant_line("", "END OF FREQUENCY"));
        b.push_str(&ant_line("", "END OF ANTENNA"));
        b
    }

    /// Block whose only frequency carries a PCO but no NOAZI table.
    fn pco_only_block(family: &str) -> String {
        let mut b = String::new();
        b.push_str(&ant_line("", "START OF ANTENNA"));
        b.push_str(&ant_line(&format!("{family:<40}"), "TYPE / SERIAL NO"));
        b.push_str(&ant_line("     0.0", "DAZI"));
        b.push_str(&ant_line("     0.0  18.0   1.0", "ZEN1 / ZEN2 / DZEN"));
        b.push_str(&ant_line("   G01", "START OF FREQUENCY"));
        b.push_str(&ant_line("      1.00      2.00      3.00", "NORTH / EAST / UP"));
        b.push_str(&ant_line("", "END OF FREQUENCY"));
        b.push_str(&ant_line("", "END OF ANTENNA"));
        b
    }

    /// Real IGS14 database; `None` skips dependent tests when the dataset
    /// is not checked out.
    fn real_db() -> Option<AntexDatabase> {
        let path = PathBuf::from("../../datasets/igs14.atx");
        if !path.exists() {
            return None;
        }
        AntexDatabase::parse(&path).ok()
    }

    // ------------------------------------------------------------------
    // Real-data lookups (datasets/igs14.atx)
    // ------------------------------------------------------------------

    #[test]
    fn from_antex_ash701945bm_scit_matches_igs14_values() {
        let Some(db) = real_db() else { return };
        let ant = ReceiverAntenna::from_antex(&db, "ASH701945B_M", "SCIT")
            .expect("igs14 must contain ASH701945B_M SCIT");
        assert_eq!(ant.ant_type, "ASH701945B_M");
        assert_eq!(ant.radome, "SCIT");
        // ANTEX NORTH/EAST/UP = 0.76 / -0.43 / 89.24 mm.
        assert!((ant.pco_enu_mm[0] - (-0.43)).abs() < 1e-9, "east {:?}", ant.pco_enu_mm);
        assert!((ant.pco_enu_mm[1] - 0.76).abs() < 1e-9, "north {:?}", ant.pco_enu_mm);
        assert!((ant.pco_enu_mm[2] - 89.24).abs() < 1e-9, "up {:?}", ant.pco_enu_mm);
        assert!((ant.zen1_deg - 0.0).abs() < 1e-12);
        assert!((ant.dzen_deg - 5.0).abs() < 1e-12);
        assert_eq!(ant.pcv_l1_grid.len(), 17); // 0..=80 deg step 5
        assert_eq!(ant.pcv_l2_grid.len(), 17);
    }

    #[test]
    fn from_antex_trm59800_80_scit_exists() {
        let Some(db) = real_db() else { return };
        let ant = ReceiverAntenna::from_antex(&db, "TRM59800.80", "SCIT")
            .expect("igs14 must contain TRM59800.80 SCIT");
        assert_eq!(ant.radome, "SCIT");
        // ANTEX N/E/U = 1.32 / 0.88 / 84.85 mm.
        assert_eq!(ant.pco_enu_mm, [0.88, 1.32, 84.85]);
        assert_eq!(ant.pcv_l1_grid.len(), 19); // 0..=90 deg step 5
    }

    #[test]
    fn from_antex_leiar20_leim_exists() {
        let Some(db) = real_db() else { return };
        let ant = ReceiverAntenna::from_antex(&db, "LEIAR20", "LEIM")
            .expect("igs14 must contain LEIAR20 LEIM");
        assert_eq!(ant.ant_type, "LEIAR20");
        assert_eq!(ant.radome, "LEIM");
        assert!(ant.dzen_deg > 0.0);
        assert!(!ant.pcv_l1_grid.is_empty());
    }

    /// Mission-specified interpolation spot-check at zenith 10 deg: an
    /// exact grid node for both bands of TRM59800.80 SCIT.
    #[test]
    fn pcv_at_zenith_10_matches_tabulated_node() {
        let Some(db) = real_db() else { return };
        let ant = ReceiverAntenna::from_antex(&db, "TRM59800.80", "SCIT").unwrap();
        assert!((ant.pcv_mm_at_zenith(1, 10.0).unwrap() - (-1.51)).abs() < 1e-9);
        assert!((ant.pcv_mm_at_zenith(2, 10.0).unwrap() - (-0.52)).abs() < 1e-9);
        // Off-node midpoint between the 10 deg (-1.51) and 15 deg (-3.05)
        // L1 nodes.
        assert!((ant.pcv_mm_at_zenith(1, 12.5).unwrap() - (-2.28)).abs() < 1e-9);
    }

    #[test]
    fn interpolation_clamps_outside_grid_range() {
        let Some(db) = real_db() else { return };
        let trm = ReceiverAntenna::from_antex(&db, "TRM59800.80", "SCIT").unwrap();
        assert_eq!(trm.pcv_mm_at_zenith(1, -5.0).unwrap(), 0.00); // first node
        assert_eq!(trm.pcv_mm_at_zenith(1, 95.0).unwrap(), 15.24); // last node
        // ASH701945B_M SCIT stops at 80 deg; beyond that it holds +3.03.
        let ash = ReceiverAntenna::from_antex(&db, "ASH701945B_M", "SCIT").unwrap();
        assert_eq!(ash.pcv_mm_at_zenith(1, 82.5).unwrap(), 3.03);
    }

    /// Different radomes are different calibrations for the same family.
    #[test]
    fn different_radomes_produce_different_pcvs() {
        let Some(db) = real_db() else { return };
        let scit = ReceiverAntenna::from_antex(&db, "ASH701945B_M", "SCIT").unwrap();
        let none = ReceiverAntenna::from_antex(&db, "ASH701945B_M", "NONE").unwrap();
        assert_ne!(scit.pcv_l1_grid, none.pcv_l1_grid);
        assert_ne!(scit.pco_enu_mm, none.pco_enu_mm);
        // A wrong radome does not silently fall back to another entry.
        assert!(ReceiverAntenna::from_antex(&db, "ASH701945B_M", "LEIM").is_none());
    }

    #[test]
    fn from_antex_unknown_type_or_bad_db_entry_is_none() {
        let Some(db) = real_db() else { return };
        assert!(ReceiverAntenna::from_antex(&db, "NOT_AN_ANTENNA", "SCIT").is_none());
        let empty = AntexDatabase::new(Vec::new());
        assert!(ReceiverAntenna::from_antex(&empty, "TRM59800.80", "SCIT").is_none());
    }

    #[test]
    fn lookup_alias_agrees_with_from_antex() {
        let Some(db) = real_db() else { return };
        let via_alias = ReceiverAntenna::lookup(&db, "TRM59800.80", "SCIT").unwrap();
        let direct = ReceiverAntenna::from_antex(&db, "TRM59800.80", "SCIT").unwrap();
        assert_eq!(via_alias.ant_type, direct.ant_type);
        assert_eq!(via_alias.pcv_l1_grid, direct.pcv_l1_grid);
        assert_eq!(via_alias.pco_enu_mm, direct.pco_enu_mm);
        assert_eq!(via_alias.antenna_type(), "TRM59800.80 SCIT");
    }

    // ------------------------------------------------------------------
    // Band selection and grid extraction (synthetic fixtures)
    // ------------------------------------------------------------------

    #[test]
    fn band_outside_1_and_2_is_none() {
        let db = db_from(&ramp_block("TEST", "NONE", -1.0, &[("G02", -2.0)]));
        let ant = ReceiverAntenna::from_antex(&db, "TEST", "NONE").unwrap();
        assert!(ant.pcv_mm_at_zenith(3, 10.0).is_none());
        assert!(ant.pcv_mm_at_zenith(0, 10.0).is_none());
        assert!(ant.pcv_mm_at_zenith(5, 10.0).is_none());
    }

    #[test]
    fn l2_grid_differs_from_l1_for_dual_band_calibration() {
        let db = db_from(&ramp_block("DUAL", "NONE", -1.0, &[("G02", -2.0)]));
        let ant = ReceiverAntenna::from_antex(&db, "DUAL", "NONE").unwrap();
        assert!((ant.pcv_mm_at_zenith(1, 10.0).unwrap() - (-10.0)).abs() < 1e-9);
        assert!((ant.pcv_mm_at_zenith(2, 10.0).unwrap() - (-20.0)).abs() < 1e-9);
    }

    #[test]
    fn glonass_only_calibration_fills_band_grids_from_r01_r02() {
        let mut b = String::new();
        b.push_str(&ant_line("", "START OF ANTENNA"));
        b.push_str(&ant_line(&format!("{:<40}", "GLO_ONLY"), "TYPE / SERIAL NO"));
        b.push_str(&ant_line("     0.0", "DAZI"));
        b.push_str(&ant_line("     0.0  18.0   1.0", "ZEN1 / ZEN2 / DZEN"));
        for (code, slope) in [("R01", -1.0), ("R02", -2.0)] {
            b.push_str(&ant_line(&format!("   {code}"), "START OF FREQUENCY"));
            b.push_str(&ant_line("      0.00      0.00      5.00", "NORTH / EAST / UP"));
            let row: Vec<String> = (0..=18).map(|i| format!("{:>6.2}", slope * i as f64)).collect();
            b.push_str(&format!("   NOAZI{}\n", row.join(" ")));
            b.push_str(&ant_line("", "END OF FREQUENCY"));
        }
        b.push_str(&ant_line("", "END OF ANTENNA"));
        let db = db_from(&b);
        let ant = ReceiverAntenna::from_antex(&db, "GLO_ONLY", "NONE").unwrap();
        assert!((ant.pcv_mm_at_zenith(1, 5.0).unwrap() - (-5.0)).abs() < 1e-9);
        assert!((ant.pcv_mm_at_zenith(2, 5.0).unwrap() - (-10.0)).abs() < 1e-9);
    }

    /// ANTEX tabulates NORTH/EAST/UP; `pco_enu_mm` must be reordered to
    /// east/north/up.
    #[test]
    fn pco_enu_reorders_antex_north_east_up_columns() {
        let db = db_from(&pco_block("PCOORD"));
        let ant = ReceiverAntenna::from_antex(&db, "PCOORD", "NONE").unwrap();
        assert_eq!(ant.pco_enu_mm, [2.0, 1.0, 3.0]);
    }

    /// A calibration with only a PCO (no NOAZI row) loads fine and its
    /// zenith interpolation is `None` rather than a wrong value.
    #[test]
    fn calibration_without_pcv_table_yields_none_interpolation() {
        let db = db_from(&pco_only_block("PCOONLY"));
        let ant = ReceiverAntenna::from_antex(&db, "PCOONLY", "NONE").unwrap();
        assert_eq!(ant.pco_enu_mm, [2.0, 1.0, 3.0]);
        assert!(ant.pcv_l1_grid.is_empty());
        assert!(ant.pcv_mm_at_zenith(1, 10.0).is_none());
        assert!(interp_grid(&[], 0.0, 5.0, 10.0).is_none());
    }

    // ------------------------------------------------------------------
    // Frequency-code mapping and DD correction math
    // ------------------------------------------------------------------

    #[test]
    fn rinex_ant_type_parses_synthetic_header() {
        let content = "\
     2.11           OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE
0220366860          TRM59800.80     SCIT                    ANT # / TYPE
                                                            END OF HEADER
";
        let dir = std::env::temp_dir();
        let path = dir.join(format!("test_ant_type_{}.20o", std::process::id()));
        std::fs::write(&path, content).unwrap();
        let result = rinex_ant_type(&path);
        std::fs::remove_file(&path).ok();
        assert_eq!(result, Some(("TRM59800.80".to_string(), "SCIT".to_string())));
    }

    #[test]
    fn rinex_ant_type_defaults_radome_to_none_for_single_token() {
        let content = "\
                                                            RINEX VERSION / TYPE
0220366860          AOAD/M_T                                ANT # / TYPE
                                                            END OF HEADER
";
        let dir = std::env::temp_dir();
        let path = dir.join(format!("test_ant_type_none_{}.20o", std::process::id()));
        std::fs::write(&path, content).unwrap();
        let result = rinex_ant_type(&path);
        std::fs::remove_file(&path).ok();
        assert_eq!(result, Some(("AOAD/M_T".to_string(), "NONE".to_string())));
    }

    #[test]
    fn rinex_ant_type_missing_file_or_header_is_none() {
        assert!(rinex_ant_type(Path::new("/nonexistent/path.20o")).is_none());

        let dir = std::env::temp_dir();
        let path = dir.join(format!("test_ant_type_missing_{}.20o", std::process::id()));
        std::fs::write(&path, "no antenna header here\nEND OF HEADER\n").unwrap();
        let result = rinex_ant_type(&path);
        std::fs::remove_file(&path).ok();
        assert!(result.is_none());
    }

    /// Real CORS RINEX2 header, gated on the dataset being checked out.
    #[test]
    fn rinex_ant_type_reads_real_cors_header() {
        let path = PathBuf::from("../../datasets/cors_short_baseline/p1811350.20o");
        if !path.exists() {
            return;
        }
        let (fam, rad) = rinex_ant_type(&path).expect("real CORS file must have ANT # / TYPE");
        assert_eq!(fam, "TRM59800.80");
        assert_eq!(rad, "SCIT");
    }

    #[test]
    fn frequency_code_mapping() {
        assert_eq!(frequency_code(Constellation::Gps, 1).as_deref(), Some("G01"));
        assert_eq!(frequency_code(Constellation::Gps, 2).as_deref(), Some("G02"));
        assert_eq!(frequency_code(Constellation::Gps, 5).as_deref(), Some("G05"));
        assert_eq!(frequency_code(Constellation::Glonass, 1).as_deref(), Some("R01"));
        assert_eq!(frequency_code(Constellation::Glonass, 2).as_deref(), Some("R02"));
        assert_eq!(frequency_code(Constellation::Galileo, 5).as_deref(), Some("E05"));
        assert_eq!(frequency_code(Constellation::Qzss, 1).as_deref(), Some("J01"));
        assert_eq!(frequency_code(Constellation::Beidou, 1), None);
        assert_eq!(frequency_code(Constellation::Gps, 42), None);
    }

    #[test]
    fn single_token_type_requires_none_radome() {
        let db = db_from(&ramp_block("AOAD/M_T", "NONE", -1.0, &[]));
        assert!(ReceiverAntenna::from_antex(&db, "AOAD/M_T", "NONE").is_some());
        assert!(ReceiverAntenna::from_antex(&db, "AOAD/M_T", "SCIT").is_none());
    }

    #[test]
    fn dd_correction_zero_for_same_antenna() {
        let db = db_from(&ramp_block("SAME", "NONE", -1.0, &[]));
        let a = ReceiverAntenna::from_antex(&db, "SAME", "NONE").unwrap();
        let corr = compute_dd_pcv_correction(&a, &a, "G01", 0.2618, std::f64::consts::FRAC_PI_4);
        assert!(corr.abs() < 1e-12);
    }

    #[test]
    fn dd_correction_is_antisymmetric_in_station_order() {
        let content = ramp_block("ANT_A", "NONE", -1.0, &[])
            + &ramp_block("ANT_B", "NONE", -2.0, &[]);
        let db = db_from(&content);
        let a = ReceiverAntenna::from_antex(&db, "ANT_A", "NONE").unwrap();
        let b = ReceiverAntenna::from_antex(&db, "ANT_B", "NONE").unwrap();
        // 80 deg elevation -> zenith 10 (grid node), 50 deg -> zenith 40
        // clamps to the last node at 18 deg. Slopes -1 vs -2 mm/deg:
        // [-10-(-18)] - [-20-(-36)] = 8 - 16 = -8 mm.
        let el_s = 80f64.to_radians();
        let el_r = 50f64.to_radians();
        let ab = compute_dd_pcv_correction(&a, &b, "G01", el_s, el_r);
        let ba = compute_dd_pcv_correction(&b, &a, "G01", el_s, el_r);
        assert!((ab - (-0.008)).abs() < 1e-9, "ab={:.6}", ab);
        assert!((ab + ba).abs() < 1e-12);
    }

    #[test]
    fn dd_correction_vanishes_when_elevations_equal() {
        let content = ramp_block("ANT_A", "NONE", -1.0, &[])
            + &ramp_block("ANT_B", "NONE", -2.0, &[]);
        let db = db_from(&content);
        let a = ReceiverAntenna::from_antex(&db, "ANT_A", "NONE").unwrap();
        let b = ReceiverAntenna::from_antex(&db, "ANT_B", "NONE").unwrap();
        assert_eq!(compute_dd_pcv_correction(&a, &b, "G01", 0.5, 0.5), 0.0);
    }

    #[test]
    fn dd_correction_zero_when_frequency_missing_on_one_side() {
        // ANT_B has no G02 table; requesting G02 must yield 0 rather than a
        // one-sided rover-only correction.
        let content = ramp_block("ANT_A", "NONE", -1.0, &[("G02", -1.0)])
            + &ramp_block("ANT_B", "NONE", -2.0, &[]);
        let db = db_from(&content);
        let a = ReceiverAntenna::from_antex(&db, "ANT_A", "NONE").unwrap();
        let b = ReceiverAntenna::from_antex(&db, "ANT_B", "NONE").unwrap();
        assert_eq!(compute_dd_pcv_correction(&a, &b, "G02", 0.2618, std::f64::consts::FRAC_PI_4), 0.0);
    }

    // ------------------------------------------------------------------
    // Real-data cross-family behaviour
    // ------------------------------------------------------------------

    /// Cross-family TRM59800.00 SCIT vs ASH701945B_M SCIT on GPS L1:
    /// non-zero but millimetre-scale, matching hand-computed grid values.
    #[test]
    fn real_cross_family_corrections_match_hand_computed_values() {
        let Some(db) = real_db() else { return };
        let trm = ReceiverAntenna::from_antex(&db, "TRM59800.00", "SCIT").unwrap();
        let ash = ReceiverAntenna::from_antex(&db, "ASH701945B_M", "SCIT").unwrap();
        let el_ref = 40f64.to_radians();
        let d10 = compute_dd_pcv_correction(&trm, &ash, "G01", 10f64.to_radians(), el_ref);
        let d30 = compute_dd_pcv_correction(&trm, &ash, "G01", 30f64.to_radians(), el_ref);
        let d60 = compute_dd_pcv_correction(&trm, &ash, "G01", 60f64.to_radians(), el_ref);
        // Hand-computed from the tabulated NOAZI grids (5 deg nodes):
        // SD(zen) = PCV_trm(zen) - PCV_ash(zen); ref sat at 40 deg el ->
        // zen50, SD(zen50) = -8.82 - (-9.18) = +0.36 mm.
        // SD(zen80)=+1.79, SD(zen60)=+0.71, SD(zen30)=-0.20 mm.
        assert!((d10 - 0.00143).abs() < 1e-4, "d10={:.6}", d10);
        assert!((d30 - 0.00035).abs() < 1e-4, "d30={:.6}", d30);
        assert!((d60 - (-0.00056)).abs() < 1e-4, "d60={:.6}", d60);
        for d in [d10, d30, d60] {
            assert!(d.abs() > 1e-4 && d.abs() < 0.02, "magnitude {d:.6}");
        }
    }

    #[test]
    fn real_leiar_vs_trm_single_difference() {
        let Some(db) = real_db() else { return };
        let trm = ReceiverAntenna::from_antex(&db, "TRM59800.00", "SCIT").unwrap();
        let lei = ReceiverAntenna::from_antex(&db, "LEIAR20", "LEIM").unwrap();
        // Single difference at 60 deg elevation (zen 30):
        // -2.26 - (-7.38) = +5.12 mm.
        let sd = lei.pcv_mm("G01", 30.0).unwrap() - trm.pcv_mm("G01", 30.0).unwrap();
        assert!((sd - 5.12).abs() < 1e-3, "sd={:.3} mm", sd);
    }

    /// Frequency dependence: the L2 correction differs from L1 for the
    /// same geometry (both calibrations tabulate G02).
    #[test]
    fn real_frequency_dependence_l1_vs_l2() {
        let Some(db) = real_db() else { return };
        let trm = ReceiverAntenna::from_antex(&db, "TRM59800.00", "SCIT").unwrap();
        let ash = ReceiverAntenna::from_antex(&db, "ASH701945B_M", "SCIT").unwrap();
        let els = (10f64.to_radians(), 40f64.to_radians());
        let l1 = compute_dd_pcv_correction(&trm, &ash, "G01", els.0, els.1);
        let l2 = compute_dd_pcv_correction(&trm, &ash, "G02", els.0, els.1);
        assert!((l1 - l2).abs() > 1e-4, "L1={:.6} L2={:.6}", l1, l2);
    }
}
