//! Print receiver-antenna differential PCV corrections from a real ANTEX
//! file — the numbers quoted in the DD-PCV design notes.
//!
//! Usage: cargo run -p gneiss-parsers --example pcv_check [path/to/igs14.atx]
//!
//! For each antenna pair and GPS L1, prints the single-difference PCV at
//! several satellite elevations and the double-difference correction
//! against a 40 deg reference satellite (the quantity embedded in one DD
//! phase observation).

use gneiss_parsers::antex::AntexDatabase;
use gneiss_parsers::receiver_antenna::{compute_dd_pcv_correction, ReceiverAntenna};

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "datasets/igs14.atx".to_string());
    let db = match AntexDatabase::parse(&path) {
        Ok(db) => db,
        Err(e) => {
            eprintln!("failed to parse {path}: {e:?}");
            std::process::exit(1);
        }
    };
    println!("loaded {} antenna blocks from {path}", db.antennas.len());

    let pairs = [
        ("TRM59800.00 SCIT", "TRM59800.00", "SCIT"),
        ("ASH701945B_M SCIT", "ASH701945B_M", "SCIT"),
        ("LEIAR20 LEIM", "LEIAR20", "LEIM"),
        ("TRM59800.80 NONE", "TRM59800.80", "NONE"),
        ("TRM29659.00 SCIT", "TRM29659.00", "SCIT"),
    ];
    let mut models = Vec::new();
    for (label, fam, rad) in pairs {
        match ReceiverAntenna::lookup(&db, fam, rad) {
            Some(m) => models.push((label, m)),
            None => {
                eprintln!("missing calibration: {label}");
                std::process::exit(1);
            }
        }
    }

    let els_deg = [10.0f64, 30.0, 60.0];
    let el_ref = 40.0f64.to_radians();

    for (name_a, _, _) in pairs {
        // index lookup is trivial at this size
        let (a_label, a) = models.iter().find(|(l, _)| *l == name_a).unwrap();
        for (b_label, b) in &models {
            if b_label == a_label {
                continue;
            }
            let code = "G01";
            print!("{a_label} vs {b_label} [{code}]  DD(el; ref=40deg):");
            for el in els_deg {
                let dd = compute_dd_pcv_correction(a, b, code, el.to_radians(), el_ref);
                print!("  {:>2.0}deg: {:+7.2} mm", el, dd * 1000.0);
            }
            println!();
        }
    }

    // Single-difference profile of the key cross-family pair.
    let (_, trm) = &models[0];
    let (_, ash) = &models[1];
    println!("\nSingle difference PCV_rov - PCV_base [mm], TRM59800.00 SCIT - ASH701945B_M SCIT (G01):");
    for zen in [0.0f64, 15.0, 30.0, 45.0, 60.0, 75.0, 85.0] {
        let s = trm.pcv_mm("G01", zen).unwrap() - ash.pcv_mm("G01", zen).unwrap();
        println!("  zen {:>4.0} deg (el {:>4.0}): {:+6.2}", zen, 90.0 - zen, s);
    }
}
