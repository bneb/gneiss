//! Receiver antenna PCO/PCV loading shared by every offline entry point
//! (benchmark binaries, `gneiss-cli`): RINEX header -> ANTEX lookup ->
//! the types [`super::PostProcessOptions`] and its base-position callers
//! consume. Relocated from `eval_network_ppk.rs` so the CLI isn't stuck
//! re-implementing an already-measured correction (docs/NETWORK_RTK_NEXT_STEPS.md).

use std::path::Path;
use std::sync::Arc;

use nalgebra::Vector3;

use gneiss_parsers::antex::AntexDatabase;
use gneiss_parsers::receiver_antenna::{rinex_ant_type, ReceiverAntenna};

use super::ReceiverPcvPair;

/// Receiver L1 PCO (ECEF, m) for a station: RINEX header antenna
/// type/radome looked up in an ANTEX file, translated ENU -> ECEF at the
/// given antenna reference point.
pub fn station_recv_pco_ecef(rinex_path: &Path, antex_path: &str, arp: Vector3<f64>) -> Option<Vector3<f64>> {
    let (fam, rad) = rinex_ant_type(rinex_path)?;
    let db = match AntexDatabase::parse(antex_path) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("DEBUG: ANTEX parse failed: {:?}", e);
            return None;
        }
    };
    let ant = ReceiverAntenna::lookup(&db, &fam, &rad)?;
    // Parser reorders the ANTEX north/east/up columns to east/north/up.
    let [east_mm, north_mm, up_mm] = ant.pco_enu_mm;

    let llh = gneiss_core::coords::ecef_to_llh(arp);
    let (lat, lon) = (llh.x, llh.y);
    let (slat, clat) = (lat.sin(), lat.cos());
    let (slon, clon) = (lon.sin(), lon.cos());
    let east = Vector3::new(-slon, clon, 0.0);
    let north = Vector3::new(-slat * clon, -slat * slon, clat);
    let up = Vector3::new(clat * clon, clat * slon, slat);
    let m = 1e-3;
    Some(north * (north_mm * m) + east * (east_mm * m) + up * (up_mm * m))
}

/// Receiver antenna PCV models (rover, base) from the ANTEX database.
/// Both headers must resolve to calibrations; otherwise `None` keeps the
/// legacy uncorrected path (a one-sided correction would be worse than
/// none).
pub fn load_receiver_pcv(rover_path: &Path, base_path: &Path, antex_path: &str) -> Option<Arc<ReceiverPcvPair>> {
    let db = match AntexDatabase::parse(antex_path) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("DEBUG: ANTEX parse failed: {:?}", e);
            return None;
        }
    };
    let (rfam, rrad) = match rinex_ant_type(rover_path) {
        Some(x) => x,
        None => {
            eprintln!("DEBUG: ant type not found in {}", rover_path.display());
            return None;
        }
    };
    let rover = ReceiverAntenna::lookup(&db, &rfam, &rrad)?;
    let (bfam, brad) = rinex_ant_type(base_path)?;
    let base = ReceiverAntenna::lookup(&db, &bfam, &brad)?;
    println!(
        "RECV-PCV enabled: rover [{}] base [{}] ({})",
        rover.antenna_type(),
        base.antenna_type(),
        antex_path
    );
    Some(Arc::new(ReceiverPcvPair {
        rover: Arc::new(rover),
        base: Arc::new(base),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn cors_dir() -> PathBuf {
        PathBuf::from("../../datasets/cors_short_baseline")
    }

    fn igs14() -> PathBuf {
        PathBuf::from("../../datasets/igs14.atx")
    }

    /// Real CORS + real ANTEX, gated on the dataset being checked out.
    #[test]
    fn load_receiver_pcv_resolves_real_cross_family_pair() {
        let (rover, base, antex) = (cors_dir().join("p2241350.20o"), cors_dir().join("capo1350.20o"), igs14());
        if !rover.exists() || !base.exists() || !antex.exists() {
            return;
        }
        let pair = load_receiver_pcv(&rover, &base, antex.to_str().unwrap())
            .expect("both CAPO and the P224 rover resolve against igs14.atx");
        assert_eq!(pair.base.ant_type, "LEIAR20");
        assert_eq!(pair.rover.ant_type, "TRM59800.00");
    }

    #[test]
    fn load_receiver_pcv_none_when_antex_missing() {
        let (rover, base) = (cors_dir().join("p2241350.20o"), cors_dir().join("capo1350.20o"));
        if !rover.exists() || !base.exists() {
            return;
        }
        assert!(load_receiver_pcv(&rover, &base, "/nonexistent/igs14.atx").is_none());
    }

    #[test]
    fn station_recv_pco_ecef_none_when_antex_missing() {
        let rover = cors_dir().join("p2241350.20o");
        if !rover.exists() {
            return;
        }
        assert!(station_recv_pco_ecef(&rover, "/nonexistent/igs14.atx", Vector3::zeros()).is_none());
    }

    /// Same station, real calibration: PCO magnitude must be sane (a few
    /// centimetres, never kilometres) and non-zero for a real antenna.
    #[test]
    fn station_recv_pco_ecef_real_capo_is_millimetre_scale() {
        let base = cors_dir().join("capo1350.20o");
        let antex = igs14();
        if !base.exists() || !antex.exists() {
            return;
        }
        let arp = Vector3::new(-2693675.7831, -4273829.9413, 3880383.2888);
        let pco = station_recv_pco_ecef(&base, antex.to_str().unwrap(), arp)
            .expect("CAPO's LEIAR20 LEIM resolves against igs14.atx");
        assert!(pco.norm() < 1.0, "PCO magnitude {} m is not millimetre-scale", pco.norm());
        assert!(pco.norm() > 1e-4, "PCO should be non-zero for a real antenna");
    }
}
