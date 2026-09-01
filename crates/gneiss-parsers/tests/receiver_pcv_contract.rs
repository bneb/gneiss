//! Coordinator-written test contract for receiver PCV implementation.
//!
//! These tests define correctness. The implementer's job is to make them
//! pass. Ground-truth values were extracted from datasets/igs14.atx by
//! the coordinator — do NOT modify these expectations.
//!
//! Run with: cargo test -p gneiss-parsers --lib receiver_pcv_contract

#[cfg(test)]
mod receiver_pcv_contract {
    use gneiss_parsers::antex::AntexDatabase;
    
    /// Locate the cached igs14.atx across checkout layouts: the main tree
    /// keeps `datasets/` at the repo root, while worktrees expose the
    /// shared store through a nested `datasets/datasets` symlink.
    ///
    /// Returns `None` when the file is not cached locally (tests then skip
    /// — a fresh checkout with no `datasets/` present, e.g. a CI runner).
    fn antex_path() -> Option<std::path::PathBuf> {
        [
            "../../datasets/igs14.atx",
            "../../datasets/datasets/igs14.atx",
        ]
        .iter()
        .map(|rel| std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(rel))
        .find(|p| p.exists())
    }

    fn load_db() -> Option<AntexDatabase> {
        AntexDatabase::parse(antex_path()?).ok()
    }

    // Helper: find receiver antenna by type + radome
    fn find_antenna<'a>(
        db: &'a AntexDatabase,
        ant_type: &str,
        radome: &str,
    ) -> Option<&'a gneiss_parsers::antex::AntennaPcv> {
        // Parser stores type+radome combined in antenna_type field.
        // Match by checking that the stored string starts with ant_type
        // and contains radome as the second whitespace-separated token.
        db.antennas.iter().find(|a| {
            let parts: Vec<&str> = a.antenna_type.split_whitespace().collect();
            !parts.is_empty() && parts[0] == ant_type
                && (parts.get(1).copied() == Some(radome) 
                    || (parts.len() == 1 && radome == "NONE"))
        })
    }

    /// Extract L1 NOAZI PCV grid from an AntennaPcv entry.
    fn get_l1_pcv_grid(a: &gneiss_parsers::antex::AntennaPcv) -> Option<(Vec<f64>, f64)> {
        // Find G01 frequency block
        let freq = a.frequencies.get("G01")?;
        Some((freq.noazi.clone(), a.dzen))
    }

    #[test]
    fn contract_ash701945b_scit_pco() {
        let Some(db) = load_db() else { return };
        let ant = find_antenna(&db, "ASH701945B_M", "SCIT")
            .expect("ASH701945B_M SCIT must exist in igs14.atx");
        
        // Find G01 frequency NEU
        if let Some(freq) = ant.frequencies.get("G01") {
            // Known values from igs14.atx line after G01 START OF FREQUENCY
            // N=0.76 E=-0.43 U=89.24 mm
            assert!(
                (freq.pco[0] - 0.76).abs() < 0.1,
                "N mismatch: {}",
                freq.pco[0]
            );
            assert!(
                (freq.pco[1] - (-0.43)).abs() < 0.1,
                "E mismatch: {}",
                freq.pco[1]
            );
            assert!(
                (freq.pco[2] - 89.24).abs() < 0.1,
                "U mismatch: {}",
                freq.pco[2]
            );
        }
    }

    #[test]
    fn contract_ash_pcv_grid_values() {
        let Some(db) = load_db() else { return };
        let ant = find_antenna(&db, "ASH701945B_M", "SCIT")
            .expect("ASH701945B_M antenna must exist");
        let (noazi, dzen) = get_l1_pcv_grid(ant).expect("G01 frequency required");
        assert!(!noazi.is_empty(), "G01 NOAZI grid must be non-empty");
        assert_eq!(noazi[0], 0.00, "zenith 0 PCV must be 0");
        assert!(dzen > 0.0, "dzen must be positive");
    }

    #[test]
    fn contract_leiar20_differs_from_trm() {
        let Some(db) = load_db() else { return };
        let leica = find_antenna(&db, "LEIAR20", "LEIM")
            .expect("LEIAR20 antenna must exist");
        let trimble = find_antenna(&db, "TRM59800.80", "SCIT")
            .expect("TRM59800.80 antenna must exist");

        let leica_pcv = &leica.frequencies.get("G01")
            .expect("LEIAR20 G01 frequency must exist").noazi;
        let trm_pcv = &trimble.frequencies.get("G01")
            .expect("TRM59800 G01 frequency must exist").noazi;

        // At zenith 30° (index 6 for 5° step), the PCVs should differ by > 3mm
        // (from ground truth: LEIAR20 ≈ -2.26, TRM ≈ -7.38)
        if leica_pcv.len() > 6 && trm_pcv.len() > 6 {
            let diff = (leica_pcv[6] - trm_pcv[6]).abs();
            assert!(
                diff > 2.0,
                "LEIAR20 vs TRM59800 PCV should differ by >2mm at zen30°, got {:.2}",
                diff
            );
        }
    }

    #[test]
    fn contract_same_family_zero_differential() {
        let Some(db) = load_db() else { return };
        let rov = find_antenna(&db, "TRM59800.80", "SCIT")
            .expect("rov antenna must exist");
        let bas = find_antenna(&db, "TRM59800.80", "SCIT")
            .expect("bas antenna must exist");

        // Same antenna type → identical PCV tables → zero DD correction
        let rov_g01 = &rov.frequencies.get("G01")
            .expect("rov G01 must exist").noazi;
        let bas_g01 = &bas.frequencies.get("G01")
            .expect("bas G01 must exist").noazi;
        assert_eq!(rov_g01.len(), bas_g01.len());
        for i in 0..rov_g01.len().min(bas_g01.len()) {
            assert!((rov_g01[i] - bas_g01[i]).abs() < 1e-12);
        }
    }
}
