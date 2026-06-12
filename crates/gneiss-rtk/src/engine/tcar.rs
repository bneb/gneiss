use nalgebra::Vector3;
use gneiss_core::obs::SatObs;
use gneiss_core::sat::SatelliteId;

const C: f64 = gneiss_core::constants::SPEED_OF_LIGHT_M_S;

/// Frequencies for GPS L1, L2, L5
const F1: f64 = 1575.42e6;
const F2: f64 = 1227.60e6;
const F5: f64 = 1176.45e6;

const LAMBDA_EWL: f64 = C / (F2 - F5);
const LAMBDA_WL: f64 = C / (F1 - F2);
const LAMBDA_L1: f64 = C / F1;

pub struct TcarResult {
    pub sat: SatelliteId,
    pub n_ewl: Option<i64>,
    pub n_wl: Option<i64>,
    pub n_nl: Option<i64>,
    pub n1: Option<i64>,
    pub n2: Option<i64>,
}

/// Computes the Extra-Wide-Lane (L2-L5) ambiguity using the Narrow-Lane pseudorange.
pub fn resolve_ewl(sat_obs: &SatObs) -> Option<i64> {
    let cp2 = sat_obs.get_observable_phase(2)?;
    let cp5 = sat_obs.get_observable_phase(5)?;
    let pr2 = sat_obs.get_observable(2)?;
    let pr5 = sat_obs.get_observable(5)?;

    let cp_ewl = cp2 - cp5; // in cycles
    // Narrow-lane pseudorange
    let pr_nl = (F2 * pr2 + F5 * pr5) / (F2 + F5); // in meters

    // N_EWL = \phi_EWL - P_NL / \lambda_EWL
    let float_amb = cp_ewl - (pr_nl / LAMBDA_EWL);
    Some(float_amb.round() as i64)
}

/// Computes the Wide-Lane (L1-L2) ambiguity using the fixed EWL ambiguity.
/// Assumes short baseline (ionosphere negligible).
pub fn resolve_wl(sat_obs: &SatObs, n_ewl_fixed: i64) -> Option<i64> {
    let cp1 = sat_obs.get_observable_phase(1)?;
    let cp2 = sat_obs.get_observable_phase(2)?;
    let cp5 = sat_obs.get_observable_phase(5)?;

    let cp_wl = cp1 - cp2; // in cycles
    let cp_ewl = cp2 - cp5;

    // Range derived from fixed EWL
    let range_ewl = (cp_ewl - n_ewl_fixed as f64) * LAMBDA_EWL;

    let float_amb = cp_wl - (range_ewl / LAMBDA_WL);
    Some(float_amb.round() as i64)
}

/// Fallback Wide-Lane resolution using Melbourne-Wubbena if L5 is not present.
pub fn resolve_wl_mw(sat_obs: &SatObs) -> Option<i64> {
    let cp1 = sat_obs.get_observable_phase(1)?;
    let cp2 = sat_obs.get_observable_phase(2)?;
    let pr1 = sat_obs.get_observable(1)?;
    let pr2 = sat_obs.get_observable(2)?;

    let cp_wl = cp1 - cp2; // in cycles
    let pr_nl = (F1 * pr1 + F2 * pr2) / (F1 + F2); // in meters

    let float_amb = cp_wl - (pr_nl / LAMBDA_WL);
    Some(float_amb.round() as i64)
}

/// Computes the Narrow-Lane (L1) ambiguity using the fixed WL ambiguity.
/// This typically requires geometry (a prior position) or very low ionosphere.
pub fn resolve_nl(sat_obs: &SatObs, _n_wl_fixed: i64, rover_pos: &Vector3<f64>, sat_pos: &Vector3<f64>, rcv_clk: f64, sat_clk: f64) -> Option<i64> {
    let cp1 = sat_obs.get_observable_phase(1)?;
    let _pr1 = sat_obs.get_observable(1)?;
    
    // Geometric range
    let geo_range = (sat_pos - rover_pos).norm() + rcv_clk - sat_clk;
    
    // N1 = cp1 - geo_range / LAMBDA_L1
    let float_amb = cp1 - (geo_range / LAMBDA_L1);
    Some(float_amb.round() as i64)
}

/// Full TCAR pipeline for a single epoch.
pub fn process_tcar_epoch(
    rover_obs: &gneiss_core::obs::EpochObs,
    _base_obs: Option<&gneiss_core::obs::EpochObs>,
    rover_pos: &Vector3<f64>,
    rcv_clk: f64,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris]
) -> Vec<TcarResult> {
    let mut results = Vec::new();
    
    for r_sat in &rover_obs.satellites {
        // Attempt EWL
        let n_ewl = resolve_ewl(r_sat);
        let n_wl;
        let mut n_nl = None;

        if let Some(ewl) = n_ewl {
            n_wl = resolve_wl(r_sat, ewl);
        } else {
            // Fallback to Geometry-Free Melbourne-Wubbena WL
            n_wl = resolve_wl_mw(r_sat);
        }

        // To resolve NL, we need sat pos/clk
        if let Some(wl) = n_wl {
            if let Some(eph) = ephemerides.iter().find(|e| e.sat() == r_sat.sat) {
                // Calculate transmit time (approx)
                let pr1 = r_sat.get_observable(1).unwrap_or(0.0);
                if pr1 > 0.0 {
                    let tx_time = rover_obs.time - (pr1 / C);
                    let (sat_pos, _sat_vel, sat_clk, _sat_drift) = eph.position(tx_time);
                    n_nl = resolve_nl(
                        r_sat, wl, rover_pos, 
                        &sat_pos,
                        rcv_clk, sat_clk * C
                    );
                }
            }
        }

        // Recover original N1 and N2 from combinations
        let mut n1 = None;
        let mut n2 = None;
        if let (Some(nl), Some(wl)) = (n_nl, n_wl) {
            n1 = Some(nl);
            n2 = Some(nl - wl); // Since N_WL = N1 - N2 => N2 = N1 - N_WL
        }

        results.push(TcarResult {
            sat: r_sat.sat,
            n_ewl,
            n_wl,
            n_nl,
            n1,
            n2
        });
    }

    results
}
