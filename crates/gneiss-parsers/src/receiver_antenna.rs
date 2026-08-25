//! Receiver antenna PCV models from an ANTEX database.
//!
//! Satellite blocks in ANTEX carry the PRN in the serial field; receiver
//! blocks carry the antenna family in field 1 and the radome in field 2 of
//! `TYPE / SERIAL NO` (e.g. `TRM59800.00     SCIT`). Both parse into
//! [`AntennaPcv`] identically; this module adds the receiver-side lookup,
//! zenith-angle interpolation of the NOAZI PCV grid, and the
//! double-difference PCV correction consumed by the DD pipeline.
//!
//! Sign convention (mirrors phase windup): the observed carrier phase
//! *contains* the antenna signature, `phi_obs = rho/lam + N + phi_PCV`, so
//! the correction returned by [`compute_dd_pcv_correction`] is **subtracted**
//! from the DD phase observation.

use gneiss_core::sat::Constellation;

use crate::antex::{AntennaPcv, AntexDatabase};

/// A receiver antenna calibration selected from an [`AntexDatabase`].
#[derive(Debug, Clone)]
pub struct ReceiverAntenna {
    model: AntennaPcv,
}

impl ReceiverAntenna {
    /// Find a receiver calibration by antenna family and radome.
    ///
    /// Matches against the whitespace-split TYPE / SERIAL NO field, so both
    /// `"TRM59800.80     SCIT"` (family + radome) and single-token entries
    /// like `"AOAD/M_T"` (radome `NONE`) resolve. The first structural match
    /// wins; IGS14 lists each family/radome combination exactly once.
    pub fn lookup(db: &AntexDatabase, family: &str, radome: &str) -> Option<Self> {
        let model = db.antennas.iter().find(|a| {
            let parts: Vec<&str> = a.antenna_type.split_whitespace().collect();
            match parts.as_slice() {
                [t, r] => *t == family && *r == radome,
                [t] => *t == family && radome == "NONE",
                _ => false,
            }
        })?;
        Some(Self { model: model.clone() })
    }

    /// The underlying calibration record (PCO vectors, validity, grids).
    pub fn model(&self) -> &AntennaPcv {
        &self.model
    }

    pub fn antenna_type(&self) -> &str {
        &self.model.antenna_type
    }

    /// NOAZI PCV in millimetres at a zenith angle in degrees.
    ///
    /// Returns `None` when the antenna has no tabulation for `freq_code`.
    /// Zenith angles outside the tabulated range clamp to the end values
    /// (the elevation mask keeps normal operation well inside the grid).
    pub fn pcv_mm(&self, freq_code: &str, zenith_deg: f64) -> Option<f64> {
        interp_noazi_mm(&self.model, freq_code, zenith_deg)
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
    // GLONASS receiver tables only expose the k=0 pair R01/R02; other
    // constellations share the band number directly (G01/G02/G05, E01/E05/E07).
    Some(format!("{prefix}{band:0>2}"))
}

/// Linearly interpolate the NOAZI PCV grid (mm) at `zenith_deg`.
///
/// Grid value `i` sits at zenith angle `zen1 + i*dzen`; interpolation is
/// linear between neighbouring nodes and clamps outside `[zen1, zen2]`.
fn interp_noazi_mm(model: &AntennaPcv, freq_code: &str, zenith_deg: f64) -> Option<f64> {
    let freq = model.frequencies.get(freq_code)?;
    let v = freq.noazi.as_slice();
    let step = model.dzen;
    if v.is_empty() || step <= 0.0 {
        return None;
    }
    let last = v.len() - 1;
    let t = (zenith_deg - model.zen1) / step;
    if t <= 0.0 {
        return Some(v[0]);
    }
    if t >= last as f64 {
        return Some(v[last]);
    }
    let i = t.floor() as usize;
    let frac = t - i as f64;
    Some(v[i] + frac * (v[i + 1] - v[i]))
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

    /// Atomically-increasing suffix so parallel `db_from` callers never
    /// collide on one temp path.
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

    /// Receiver-style block with a linear NOAZI ramp: PCV(zen) = -1 mm/deg.
    /// Radome `NONE` yields the single-token TYPE form used by unradomed
    /// calibrations.
    fn linear_block(family: &str, radome: &str) -> String {
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
        b.push_str(&ant_line("   G01", "START OF FREQUENCY"));
        b.push_str(&ant_line("      0.00      0.00      5.00", "NORTH / EAST / UP"));
        let row: Vec<String> = (0..=18).map(|i| format!("{:>6.2}", -(i as f64))).collect();
        b.push_str(&format!("   NOAZI{}\n", row.join("")));
        b.push_str(&ant_line("", "END OF FREQUENCY"));
        b.push_str(&ant_line("", "END OF ANTENNA"));
        b
    }

    fn ant_line(data: &str, label: &str) -> String {
        let mut line = String::from(data);
        while line.len() < 60 {
            line.push(' ');
        }
        line.push_str(label);
        line.push('\n');
        line
    }

    #[test]
    fn test_lookup_family_radome() {
        let db = db_from(&linear_block("TRM59800.00", "SCIT"));
        let ant = ReceiverAntenna::lookup(&db, "TRM59800.00", "SCIT").unwrap();
        assert_eq!(ant.antenna_type(), "TRM59800.00 SCIT");
    }

    #[test]
    fn test_lookup_single_token_requires_none_radome() {
        let db = db_from(&linear_block("AOAD/M_T", "NONE"));
        assert!(ReceiverAntenna::lookup(&db, "AOAD/M_T", "NONE").is_some());
        assert!(ReceiverAntenna::lookup(&db, "AOAD/M_T", "SCIT").is_none());
    }

    #[test]
    fn test_lookup_missing_returns_none() {
        let db = db_from(&linear_block("TRM59800.00", "SCIT"));
        assert!(ReceiverAntenna::lookup(&db, "LEIAR20", "LEIM").is_none());
    }

    #[test]
    fn test_interp_exact_grid_nodes() {
        let db = db_from(&linear_block("TEST", "NONE"));
        let ant = ReceiverAntenna::lookup(&db, "TEST", "NONE").unwrap();
        assert!((ant.pcv_mm("G01", 0.0).unwrap() - 0.0).abs() < 1e-12);
        assert!((ant.pcv_mm("G01", 5.0).unwrap() - (-5.0)).abs() < 1e-12);
        assert!((ant.pcv_mm("G01", 18.0).unwrap() - (-18.0)).abs() < 1e-12);
    }

    #[test]
    fn test_interp_midpoint() {
        let db = db_from(&linear_block("TEST", "NONE"));
        let ant = ReceiverAntenna::lookup(&db, "TEST", "NONE").unwrap();
        assert!((ant.pcv_mm("G01", 12.5).unwrap() - (-12.5)).abs() < 1e-12);
        assert!((ant.pcv_mm("G01", 2.3).unwrap() - (-2.3)).abs() < 1e-12);
    }

    #[test]
    fn test_interp_clamps_outside_range() {
        let db = db_from(&linear_block("TEST", "NONE"));
        let ant = ReceiverAntenna::lookup(&db, "TEST", "NONE").unwrap();
        assert_eq!(ant.pcv_mm("G01", -10.0).unwrap(), 0.0);
        assert_eq!(ant.pcv_mm("G01", 500.0).unwrap(), -18.0);
    }

    #[test]
    fn test_pcv_unknown_frequency_is_none() {
        let db = db_from(&linear_block("TEST", "NONE"));
        let ant = ReceiverAntenna::lookup(&db, "TEST", "NONE").unwrap();
        assert!(ant.pcv_mm("G02", 30.0).is_none());
    }

    #[test]
    fn test_frequency_code_mapping() {
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
    fn test_dd_correction_zero_for_same_antenna() {
        let db = db_from(&linear_block("SAME", "NONE"));
        let a = ReceiverAntenna::lookup(&db, "SAME", "NONE").unwrap();
        let corr = compute_dd_pcv_correction(&a, &a, "G01", 0.2618, 0.7854);
        assert!(corr.abs() < 1e-12);
    }

    #[test]
    fn test_dd_correction_antisymmetric() {
        let db = db_from(&format!("{}{}", linear_block("ANT_A", "NONE"), linear_block("ANT_B", "NONE")));
        let a = ReceiverAntenna::lookup(&db, "ANT_A", "NONE").unwrap();
        let b = ReceiverAntenna::lookup(&db, "ANT_B", "NONE").unwrap();
        let ab = compute_dd_pcv_correction(&a, &b, "G01", 0.2618, 0.7854);
        let ba = compute_dd_pcv_correction(&b, &a, "G01", 0.2618, 0.7854);
        assert!(ab.abs() > 1e-6);
        assert!((ab + ba).abs() < 1e-12);
    }

    #[test]
    fn test_dd_correction_zero_when_elevations_equal() {
        let db = db_from(&format!("{}{}", linear_block("ANT_A", "NONE"), linear_block("ANT_B", "NONE")));
        let a = ReceiverAntenna::lookup(&db, "ANT_A", "NONE").unwrap();
        let b = ReceiverAntenna::lookup(&db, "ANT_B", "NONE").unwrap();
        // Identical zenith angles cancel in every bracket of the DD.
        assert_eq!(compute_dd_pcv_correction(&a, &b, "G01", 0.5, 0.5), 0.0);
    }

    #[test]
    fn test_dd_correction_zero_when_frequency_missing_on_one_side() {
        // ANT_B block has only G01; requesting G02 must yield 0, not a
        // one-sided rover-only correction.
        let content = format!(
            "{}{}",
            linear_block("ANT_A", "NONE"),
            linear_block("ANT_B", "NONE")
        );
        let db = db_from(&content);
        let a = ReceiverAntenna::lookup(&db, "ANT_A", "NONE").unwrap();
        let b = ReceiverAntenna::lookup(&db, "ANT_B", "NONE").unwrap();
        assert_eq!(compute_dd_pcv_correction(&a, &b, "G02", 0.2618, 0.7854), 0.0);
    }

    /// Real-calibration smoke test: cross-family TRM59800.00 SCIT vs
    /// ASH701945B_M SCIT on GPS L1 must be non-zero but millimetre-scale.
    #[test]
    fn test_real_igs14_cross_family_magnitudes() {
        let path = PathBuf::from("../../datasets/igs14.atx");
        if !path.exists() {
            return; // dataset not available in this checkout
        }
        let db = AntexDatabase::parse(&path).unwrap();
        let trm = ReceiverAntenna::lookup(&db, "TRM59800.00", "SCIT").expect("TRM59800.00 SCIT");
        let ash = ReceiverAntenna::lookup(&db, "ASH701945B_M", "SCIT").expect("ASH701945B_M SCIT");

        // Reference satellite held at 40 deg elevation.
        let d10 = compute_dd_pcv_correction(&trm, &ash, "G01", 10f64.to_radians(), 40f64.to_radians());
        let d30 = compute_dd_pcv_correction(&trm, &ash, "G01", 30f64.to_radians(), 40f64.to_radians());
        let d60 = compute_dd_pcv_correction(&trm, &ash, "G01", 60f64.to_radians(), 40f64.to_radians());

        // Reference satellite held at 40 deg elevation. Hand-computed from
        // the tabulated NOAZI grids: single difference SD(zen) =
        // PCV_trm(zen) - PCV_ash(zen), e.g. SD(zen50) = -8.82 - (-9.18)
        // = +0.36 mm; DD(el) = SD(90deg-el) - SD(zen50).
        assert!(d10.abs() > 1e-4 && d10.abs() < 0.02, "d10={:.6} m", d10);
        assert!(d30.abs() > 1e-4 && d30.abs() < 0.02, "d30={:.6} m", d30);
        assert!(d60.abs() > 1e-4 && d60.abs() < 0.02, "d60={:.6} m", d60);

        // Exact grid-node values (metres).
        assert!((d10 - 0.00143).abs() < 1e-4, "d10={:.6}", d10);
        assert!((d30 - 0.00035).abs() < 1e-4, "d30={:.6}", d30);
        assert!((d60 - 0.00057).abs() < 1e-4, "d60={:.6}", d60);
    }

    /// LEIAR20 LEIM vs TRM59800.00 SCIT diverges much more at high
    /// elevation than the Ashtech pair (~4 mm single difference).
    #[test]
    fn test_real_igs14_leiar_vs_trm() {
        let path = PathBuf::from("../../datasets/igs14.atx");
        if !path.exists() {
            return;
        }
        let db = AntexDatabase::parse(&path).unwrap();
        let trm = ReceiverAntenna::lookup(&db, "TRM59800.00", "SCIT").unwrap();
        let lei = ReceiverAntenna::lookup(&db, "LEIAR20", "LEIM").unwrap();
        // Single difference at 60 deg elevation (zen 30):
        // -2.26 - (-7.38) = +5.12 mm.
        let sd = lei.pcv_mm("G01", 30.0).unwrap() - trm.pcv_mm("G01", 30.0).unwrap();
        assert!((sd - 5.12).abs() < 1e-3, "sd={:.3} mm", sd);
    }

    /// Frequency dependence: L2 correction differs from L1 for the same
    /// geometry (both calibrations tabulate G02).
    #[test]
    fn test_real_igs14_frequency_dependence() {
        let path = PathBuf::from("../../datasets/igs14.atx");
        if !path.exists() {
            return;
        }
        let db = AntexDatabase::parse(&path).unwrap();
        let trm = ReceiverAntenna::lookup(&db, "TRM59800.00", "SCIT").unwrap();
        let ash = ReceiverAntenna::lookup(&db, "ASH701945B_M", "SCIT").unwrap();
        let l1 = compute_dd_pcv_correction(&trm, &ash, "G01", 10f64.to_radians(), 40f64.to_radians());
        let l2 = compute_dd_pcv_correction(&trm, &ash, "G02", 10f64.to_radians(), 40f64.to_radians());
        assert!((l1 - l2).abs() > 1e-4, "L1={:.6} L2={:.6}", l1, l2);
    }
}
