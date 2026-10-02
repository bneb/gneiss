//! Receiver antenna phase-centre variation (PCV) for double-difference RTK.
//!
//! Wraps one receiver calibration from an [`AntexDatabase`] (matched by
//! antenna family + radome) and exposes linear interpolation of its L1
//! NOAZI PCV table plus the differential double-difference correction.
//!
//! Grid convention: PCV node `i` sits at zenith angle
//! `zen_start_deg + i * zen_step_deg`; values are millimetres. Requests
//! outside the tabulated range clamp to the nearest end node.
//!
//! Sign convention: the observed carrier phase *contains* the antenna
//! signature, so [`ReceiverPcv::dd_correction_m`] is **subtracted** from
//! the phase observation, mirroring the phase-windup convention.

use crate::antex::{AntennaPcv, AntexDatabase, FrequencyPcv};
use nalgebra::Vector3;

type Vec3 = Vector3<f64>;

/// Preferred ANTEX frequency code for GPS L1 receiver calibrations.
const L1_PRIMARY: &str = "G01";
/// Legacy two-character L1 code, tried when `G01` is absent.
const L1_FALLBACK: &str = "G1";
/// Millimetres per metre.
const MM_PER_M: f64 = 1000.0;

/// Receiver antenna phase centre variation model.
///
/// `from_antex` fills the fields below; `ant_type`/`radome` echo the
/// requested identifiers so callers can label corrections in logs.
#[derive(Debug, Clone)]
pub struct ReceiverPcv {
    /// Antenna family as requested, e.g. `"ASH701945B_M"`.
    pub ant_type: String,
    /// Radome code as requested (`"SCIT"`), or `"NONE"` when unradomed.
    pub radome: String,
    /// L1 phase centre offset, north/east/up, millimetres (ANTEX order).
    pub pco_neu_mm: nalgebra::Vector3<f64>,
    /// Elevation-dependent L1 PCV values (mm), indexed by zenith grid.
    pub pcv_grid_mm: Vec<f64>,
    /// Azimuth-dependent PCV grid (mm), if present in ANTEX.
    pub azi_grid_mm: Option<Vec<Vec<f64>>>,
    /// Azimuth grid spacing (deg).
    pub dazi_deg: f64,
    /// Zenith angle (deg) of the first grid node.
    pub zen_start_deg: f64,
    /// Zenith grid spacing (deg).
    pub zen_step_deg: f64,
}

impl ReceiverPcv {
    /// Load from ANTEX database by matching type+radome.
    ///
    /// Matches the whitespace-split `TYPE / SERIAL NO` field, so both
    /// `"TRM59800.80     SCIT"` (family + radome) and single-token entries
    /// like `"AOAD/M_T"` (radome `"NONE"` or `""`) resolve. The L1 table
    /// prefers code `G01` and falls back to legacy `G1`. Returns `None`
    /// when no such receiver entry exists or its L1 grid is unusable
    /// (missing, empty, or non-positive zenith step).
    pub fn from_antex(db: &AntexDatabase, ant_type: &str, radome: &str) -> Option<Self> {
        let model = db
            .antennas
            .iter()
            .find(|a| type_radome_matches(&a.antenna_type, ant_type, radome))?;
        let freq = l1_frequency(model)?;
        if model.dzen <= 0.0 || freq.noazi.is_empty() {
            return None;
        }
        Some(Self {
            ant_type: ant_type.to_string(),
            radome: radome.to_string(),
            pco_neu_mm: freq.pco,
            pcv_grid_mm: freq.noazi.clone(),
            azi_grid_mm: freq.azi.clone(),
            dazi_deg: model.dazi,
            zen_start_deg: model.zen1,
            zen_step_deg: model.dzen,
        })
    }

    /// Interpolate PCV (mm) at a zenith angle in degrees (elevation-only NOAZI).
    ///
    /// Uses linear interpolation between grid points, clamps outside range.
    /// Returns 0.0 for an (unusable but constructible) empty grid.
    pub fn interpolate(&self, zenith_deg: f64) -> f64 {
        let grid = &self.pcv_grid_mm;
        if grid.is_empty() || self.zen_step_deg <= 0.0 {
            return 0.0;
        }
        let last = grid.len() - 1;
        let pos = (zenith_deg - self.zen_start_deg) / self.zen_step_deg;
        if pos <= 0.0 {
            return grid[0];
        }
        if pos >= last as f64 {
            return grid[last];
        }
        let i = pos.floor() as usize; // < last, so i + 1 stays in bounds
        let frac = pos - i as f64;
        grid[i] * (1.0 - frac) + grid[i + 1] * frac
    }

    /// Interpolate PCV (mm) at azimuth and zenith angles (degrees).
    ///
    /// Uses 2D bilinear interpolation if a well-formed azimuth grid is
    /// present; otherwise falls back to 1D elevation-dependent
    /// interpolation, which is also the fallback for a truncated grid.
    pub fn interpolate_az_zen(&self, az_deg: f64, zen_deg: f64) -> f64 {
        if let Some(ref azi) = self.azi_grid_mm {
            if self.dazi_deg > 0.0 && self.zen_step_deg > 0.0 {
                if let Some(v) = self.az_zen_grid_mm(azi, az_deg, zen_deg) {
                    return v;
                }
            }
        }
        self.interpolate(zen_deg)
    }

    /// Bilinear sample of the RO-AZI grid, or `None` when the grid cannot
    /// be indexed safely.
    ///
    /// `None` covers an empty grid and rows of unequal length: a short row
    /// would make `azi[az_i][zen_i]` read past its end, and the RO-AZI
    /// table is then less trustworthy than the NOAZI table, so the caller
    /// degrades to the elevation-only model instead of guessing.
    fn az_zen_grid_mm(&self, azi: &[Vec<f64>], az_deg: f64, zen_deg: f64) -> Option<f64> {
        let n_azi = azi.len();
        let n_zen = azi.first()?.len();
        if n_zen == 0 || azi.iter().any(|row| row.len() != n_zen) {
            return None;
        }

        let az_norm = az_deg.rem_euclid(360.0);
        let az_pos = az_norm / self.dazi_deg;
        let az_i = (az_pos.floor() as usize) % n_azi;
        let az_next = (az_i + 1) % n_azi;
        let az_frac = az_pos - az_pos.floor();

        let zen_last = (n_zen - 1) as f64;
        let zen_pos = ((zen_deg - self.zen_start_deg) / self.zen_step_deg).clamp(0.0, zen_last);
        let zen_i = zen_pos.floor() as usize;
        let zen_next = (zen_i + 1).min(n_zen - 1);
        let zen_frac = zen_pos - zen_i as f64;

        let v0 = azi[az_i][zen_i] * (1.0 - zen_frac) + azi[az_i][zen_next] * zen_frac;
        let v1 = azi[az_next][zen_i] * (1.0 - zen_frac) + azi[az_next][zen_next] * zen_frac;
        Some(v0 * (1.0 - az_frac) + v1 * az_frac)
    }

    /// Compute DIFFERENTIAL DD correction in metres:
    /// `[PCV_rov(zen_s) - PCV_rov(zen_r)] - [PCV_base(zen_s) - PCV_base(zen_r)]`.
    ///
    /// Returns the correction to subtract from carrier-phase observations.
    /// Same-type antennas cancel to ~0; cross-type pairs differ by the
    /// millimetre-scale antenna signature.
    ///
    /// SAFETY: the four angles are bare `f64`s, so a caller can pass a
    /// geometrically impossible pair (a reference satellite below the
    /// horizon at one station and above it at the other, or sat/ref
    /// swapped between rover and base). Prefer
    /// [`Self::dd_correction_from_geometry`], which derives all four from
    /// shared satellite positions and cannot represent that state.
    pub fn dd_correction_m(
        rov: &ReceiverPcv,
        bas: &ReceiverPcv,
        zen_rov_sat_deg: f64,
        zen_rov_ref_deg: f64,
        zen_bas_sat_deg: f64,
        zen_bas_ref_deg: f64,
    ) -> f64 {
        let rov_diff_mm = rov.interpolate(zen_rov_sat_deg) - rov.interpolate(zen_rov_ref_deg);
        let bas_diff_mm = bas.interpolate(zen_bas_sat_deg) - bas.interpolate(zen_bas_ref_deg);
        (rov_diff_mm - bas_diff_mm) / MM_PER_M
    }

    /// DD correction derived from shared satellite geometry.
    ///
    /// Zenith and azimuth angles are computed internally from station and
    /// satellite positions, making inconsistent geometry unrepresentable.
    /// `*_llh` are geodetic radians + ellipsoidal height; `sat_pos` and
    /// `ref_sat_pos` are ECEF metres.
    pub fn dd_correction_from_geometry(
        rov: &ReceiverPcv,
        bas: &ReceiverPcv,
        rov_llh: Vec3,
        bas_llh: Vec3,
        sat_pos: Vec3,
        ref_sat_pos: Vec3,
    ) -> f64 {
        let rov_diff_mm = receiver_dd_mm(rov, rov_llh, sat_pos, ref_sat_pos);
        let bas_diff_mm = receiver_dd_mm(bas, bas_llh, sat_pos, ref_sat_pos);
        (rov_diff_mm - bas_diff_mm) / MM_PER_M
    }

    /// Compute standalone receiver PCO correction in metres projected along LOS towards satellite.
    pub fn pco_correction_m(&self, az_deg: f64, el_deg: f64) -> f64 {
        let az_rad = az_deg.to_radians();
        let el_rad = el_deg.to_radians();
        let north_m = self.pco_neu_mm.x * 1e-3;
        let east_m = self.pco_neu_mm.y * 1e-3;
        let up_m = self.pco_neu_mm.z * 1e-3;
        let u_east = az_rad.sin() * el_rad.cos();
        let u_north = az_rad.cos() * el_rad.cos();
        let u_up = el_rad.sin();
        east_m * u_east + north_m * u_north + up_m * u_up
    }

    /// Compute standalone receiver PCV correction in metres along LOS.
    pub fn pcv_correction_m(&self, az_deg: f64, el_deg: f64) -> f64 {
        let zen_deg = (90.0 - el_deg).max(0.0);
        self.interpolate_az_zen(az_deg, zen_deg) / MM_PER_M
    }

    /// Total standalone receiver antenna correction in metres (PCO + PCV) along LOS.
    pub fn total_correction_m(&self, az_deg: f64, el_deg: f64) -> f64 {
        self.pco_correction_m(az_deg, el_deg) + self.pcv_correction_m(az_deg, el_deg)
    }
}

/// PCV difference (mm) between the observation satellite and the
/// reference satellite, seen from one station.
///
/// `llh` supplies the local ENU basis and `sat`/`ref_sat` are ECEF. The
/// two frames are different objects: `gneiss_core::coords::az_el` needs
/// the geodetic position *and* the ECEF position, so the ECEF origin is
/// derived here from the single `llh` argument rather than being passed
/// as a parallel vector that could disagree with it.
fn receiver_dd_mm(pcv: &ReceiverPcv, llh: Vec3, sat: Vec3, ref_sat: Vec3) -> f64 {
    let ecef = gneiss_core::coords::llh_to_ecef(llh);
    let (az_s, zen_s) = az_zen_deg(llh, ecef, sat);
    let (az_r, zen_r) = az_zen_deg(llh, ecef, ref_sat);
    pcv.interpolate_az_zen(az_s, zen_s) - pcv.interpolate_az_zen(az_r, zen_r)
}

/// Azimuth (deg, [0, 360)) and zenith angle (deg) of `sat` as seen from a
/// station given in both geodetic and ECEF coordinates.
fn az_zen_deg(llh: Vec3, ecef: Vec3, sat: Vec3) -> (f64, f64) {
    let (az, el) = gneiss_core::coords::az_el(llh, ecef, sat);
    (az.to_degrees().rem_euclid(360.0), (90.0 - el.to_degrees()).max(0.0))
}

/// Match a stored `TYPE / SERIAL NO` string against family + radome.
/// Receiver blocks store `"FAMILY  RADOME"`; single-token entries
/// (e.g. `"AOAD/M_T"`) denote unradomed calibrations and resolve for the
/// radome codes `"NONE"` or `""`.
fn type_radome_matches(stored: &str, family: &str, radome: &str) -> bool {
    let mut tokens = stored.split_whitespace();
    if tokens.next() != Some(family) {
        return false;
    }
    match tokens.next() {
        Some(stored_radome) => stored_radome == radome,
        None => radome == "NONE" || radome.is_empty(),
    }
}

/// Locate the L1 NOAZI table of one calibration block.
fn l1_frequency(model: &AntennaPcv) -> Option<&FrequencyPcv> {
    model
        .frequencies
        .get(L1_PRIMARY)
        .or_else(|| model.frequencies.get(L1_FALLBACK))
}


#[cfg(test)]
mod tests;
