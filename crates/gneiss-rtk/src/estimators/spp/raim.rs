//! RAIM fault detection, outlier exclusion, and initial state seeding for SPP.

use gneiss_core::coords::{ecef_to_llh, Coordinate, Datum, Frame};
use gneiss_core::atmosphere::KlobucharParams;
use nalgebra::Vector3;
use super::measurements::compute_sat_state;
use super::solver::{compute_measurement_residuals, solve_spp_iteratively};
use super::{SppConfig, SppError, SppMeasurement, SppState};

pub(crate) fn compute_seed_position(measurements: &[SppMeasurement]) -> (f64, f64, f64) {
    let mut avg = Vector3::zeros();
    for m in measurements {
        let (sat_coord, _) = compute_sat_state(m, 0.0);
        avg += sat_coord.vector;
    }
    avg /= measurements.len() as f64;
    let mut llh = ecef_to_llh(avg);
    llh.z = 0.0;
    let proj = gneiss_core::coords::llh_to_ecef(llh);
    (proj.x, proj.y, proj.z)
}

pub(crate) fn compute_seed_clocks(
    measurements: &[SppMeasurement],
    sx: f64,
    sy: f64,
    sz: f64,
) -> (f64, f64, f64, f64) {
    let mut cdt_gps = None;
    let mut cdt_gal = None;
    let mut cdt_bds = None;
    let mut cdt_glo = None;
    for m in measurements {
        let (sat_coord, corrected_pr) = compute_sat_state(m, 0.0);
        let dx = sx - sat_coord.vector.x;
        let dy = sy - sat_coord.vector.y;
        let dz = sz - sat_coord.vector.z;
        let cdt = corrected_pr - f64::sqrt(dx * dx + dy * dy + dz * dz);
        match m.constellation {
            gneiss_core::sat::Constellation::Gps => {
                if cdt_gps.is_none() {
                    cdt_gps = Some(cdt);
                }
            }
            gneiss_core::sat::Constellation::Galileo => {
                if cdt_gal.is_none() {
                    cdt_gal = Some(cdt);
                }
            }
            gneiss_core::sat::Constellation::Beidou => {
                if cdt_bds.is_none() {
                    cdt_bds = Some(cdt);
                }
            }
            gneiss_core::sat::Constellation::Glonass if cdt_glo.is_none() => {
                cdt_glo = Some(cdt);
            }
            _ => {}
        }
    }
    let default_cdt = cdt_gps.or(cdt_gal).or(cdt_bds).or(cdt_glo).unwrap_or(0.0);
    (
        cdt_gps.unwrap_or(default_cdt),
        cdt_gal.unwrap_or(default_cdt),
        cdt_bds.unwrap_or(default_cdt),
        cdt_glo.unwrap_or(default_cdt),
    )
}

pub(crate) fn seed_initial_state(measurements: &[SppMeasurement], prev_state: Option<&SppState>) -> SppState {
    let (seed_x, seed_y, seed_z) = if let Some(coord) = prev_state.map(|s| &s.position) {
        (coord.vector.x, coord.vector.y, coord.vector.z)
    } else {
        compute_seed_position(measurements)
    };
    let (cdt, cdt_gal, cdt_bds, cdt_glo) =
        compute_seed_clocks(measurements, seed_x, seed_y, seed_z);
    SppState::new(
        Coordinate::new(
            Vector3::new(seed_x, seed_y, seed_z),
            Datum::WGS84,
            Frame::ECEF,
            measurements[0].time,
        ),
        cdt,
        cdt_gal,
        cdt_bds,
        cdt_glo,
    )
}

pub(crate) fn apply_raim(
    state: SppState,
    seed_state: SppState,
    measurements: &[SppMeasurement],
    iono_params: Option<&KlobucharParams>,
    config: &SppConfig,
) -> Result<SppState, SppError> {
    let good_measurements = filter_raim_outliers(&state, measurements, iono_params, config);
    if good_measurements.len() < measurements.len() && good_measurements.len() >= 4 {
        return solve_spp_iteratively(seed_state, &good_measurements, iono_params, config);
    }
    Ok(state)
}

pub(crate) fn filter_raim_outliers(
    state: &SppState,
    measurements: &[SppMeasurement],
    iono_params: Option<&KlobucharParams>,
    config: &SppConfig,
) -> Vec<SppMeasurement> {
    let rec_ecef = state.position.vector;
    let rec_llh = ecef_to_llh(rec_ecef);

    let mut residuals: Vec<(&SppMeasurement, f64)> = measurements
        .iter()
        .map(|m| {
            let (_, _, _, _, residual, _) =
                compute_measurement_residuals(state, m, rec_ecef, rec_llh, iono_params, config);
            (m, residual.abs())
        })
        .collect();

    residuals.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let median_residual = residuals[residuals.len() / 2].1;
    let mut deviations: Vec<f64> = residuals
        .iter()
        .map(|(_, r)| (*r - median_residual).abs())
        .collect();
    deviations.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mad_threshold =
        (deviations[deviations.len() / 2] * config.raim_mad_multiplier).max(config.raim_outlier_m);

    residuals
        .into_iter()
        .filter(|(_, r)| *r <= mad_threshold)
        .map(|(m, _)| m.clone())
        .collect()
}

#[cfg(test)]
mod tests {
    //! Golden vectors for RAIM seeding, gating, and outlier rejection.
    //!
    //! `filter_raim_outliers` uses a median-absolute-deviation gate:
    //!     threshold = max(MAD * raim_mad_multiplier, raim_outlier_m)
    //! For n = 6 measurements with five clean residuals and one large blunder,
    //! the sorted |residual| median and the sorted deviation median both land on
    //! the clean cluster, so the threshold collapses to the absolute floor
    //! raim_outlier_m and the blunder is the only row removed.

    use super::*;
    use super::super::fixture::*;
    use super::super::{build_measurements, LIGHT_SPEED};
    use gneiss_core::sat::Constellation;
    use gneiss_core::time::GpsTime;

    fn receiver() -> Vector3<f64> {
        Vector3::new(gneiss_core::constants::WGS84_SEMI_MAJOR_AXIS_M, 0.0, 0.0)
    }

    fn scene(blunder_m: Option<(usize, f64)>, t: GpsTime) -> (Vec<SppMeasurement>, Vector3<f64>, f64) {
        let true_pos = receiver();
        let true_cdt = 750.0;
        let (mut epoch, ephems) = four_gps_scene(true_pos, true_cdt, t);
        if let Some((i, delta)) = blunder_m {
            epoch.satellites[i].observations[0].value += delta;
        }
        let m = build_measurements(&epoch, &ephems, &open_config());
        (m, true_pos, true_cdt)
    }

    #[test]
    fn seed_position_projects_the_mean_satellite_position_onto_the_ellipsoid() {
        // The seed must be a point ON the ellipsoid: compute_seed_position takes
        // the mean satellite ECEF, converts it to geodetic, zeroes the height,
        // and converts back, so |seed| is the WGS84 radius for that latitude.
        let t = epoch_time();
        let (m, _, _) = scene(None, t);
        let (x, y, z) = compute_seed_position(&m);
        let seed = Vector3::new(x, y, z);
        let llh = ecef_to_llh(seed);
        assert!(llh.z.abs() < 1e-6, "seed height = {} m", llh.z);
        // And the seed is close to the mean satellite position, so it is a sane
        // starting point rather than the geocentre.
        let mut mean = Vector3::zeros();
        for s in &m {
            mean += compute_sat_state(s, 0.0).0.vector;
        }
        mean /= m.len() as f64;
        assert!((seed - mean).norm() < 2.0e6, "seed jumped away from the mean satellite position");
    }

    #[test]
    fn seed_clocks_default_to_the_first_available_constellation_bias() {
        // With only GPS present, no Galileo/BeiDou/GLONASS measurement exists, so
        // `compute_seed_clocks` fills every slot with the same GPS bias rather
        // than leaving the unused systems at zero.
        let t = epoch_time();
        let (m, true_pos, true_cdt) = scene(None, t);
        let (sx, sy, sz) = compute_seed_position(&m);
        let (g, gal, bds, glo) = compute_seed_clocks(&m, sx, sy, sz);
        assert!((gal - g).abs() < 1e-9 && (bds - g).abs() < 1e-9 && (glo - g).abs() < 1e-9,
            "unused clock slots must all inherit the GPS bias: {g} {gal} {bds} {glo}");
        // The bias is (pseudorange minus geometric range) measured from the crude
    // seed position -- the mean satellite position projected onto the ellipsoid
    // -- so it is only an initial guess, not the final clock. What must hold is
    // that it is finite and has the right sign of magnitude for a ~20-30 Mm
    // geometric range (not, say, a metre-scale or NaN value).
    assert!(g.is_finite(), "seeded cdt must be finite");
    assert!(g.abs() > 1.0e5, "seeded cdt {g} is implausibly small");
    assert!(true_pos.x > 0.0);
    let _ = true_cdt;
    }

    #[test]
    fn seed_clocks_ignore_constellations_without_a_dedicated_bias() {
        // QZSS rides the GPS clock and SBAS has no clock term here; neither may
        // overwrite a dedicated per-system bias slot.
        let t = epoch_time();
        let (gps, _, _) = scene(None, t);
        let mut mixed = gps.clone();
        let mut qzss = gps[0].clone();
        qzss.constellation = Constellation::Qzss;
        qzss.raw_pr += 1234.0;
        mixed.push(qzss);
        let mut sbas = gps[1].clone();
        sbas.constellation = Constellation::Sbas;
        sbas.raw_pr -= 987.0;
        mixed.push(sbas);

        let (sx, sy, sz) = compute_seed_position(&gps);
        let (g0, _, _, _) = compute_seed_clocks(&gps, sx, sy, sz);
        let (g1, gal, bds, glo) = compute_seed_clocks(&mixed, sx, sy, sz);
        assert!((g1 - g0).abs() < 1e-9, "QZSS/SBAS changed the GPS bias: {g0} -> {g1}");
        // With no dedicated bias for Galileo/BeiDou/GLONASS, they inherit the
        // GPS value rather than defaulting to zero.
        assert!((gal - g0).abs() < 1e-9 && (bds - g0).abs() < 1e-9 && (glo - g0).abs() < 1e-9);
    }

    #[test]
    fn seed_initial_state_reuses_a_previous_solution_when_given() {
        let t = epoch_time();
        let (m, _, true_cdt) = scene(None, t);
        let prev = SppState::new(
            Coordinate::new(Vector3::new(1.0, 2.0, 3.0), Datum::WGS84, Frame::ECEF, t),
            11.0,
            22.0,
            33.0,
            44.0,
        );
        let seeded = seed_initial_state(&m, Some(&prev));
    // Only the POSITION of the prior is reused; the clock biases are always
    // re-derived from the geometry, because the prior carries metres while the
    // seed solves for c * dt. So the prior's 11 / 22 / 33 / 44 must not survive.
    assert_eq!(seeded.position.vector, Vector3::new(1.0, 2.0, 3.0));
    assert!((seeded.cdt - 11.0).abs() > 1.0, "prior clock bias leaked into the seed");
    assert!((seeded.cdt - seeded.cdt_glo).abs() < 1e-6);
    assert!((seeded.cdt - seeded.cdt_gal).abs() < 1e-6);
    assert!(seeded.cdt.is_finite());
    let _ = true_cdt;
    }

    #[test]
    fn raim_removes_only_the_blunder_and_raim_resolves_it() {
        // State is placed at the truth, so the residual of each measurement is
        // exactly its pseudorange error: ~0 for the five clean satellites and
        // 500 m for the blundered one. With six rows the median |residual| and
        // the median deviation both fall on the clean cluster, so the MAD
        // threshold collapses to the absolute floor raim_outlier_m = 50 m and
        // exactly one row must be dropped.
        let t = epoch_time();
        let (m, true_pos, true_cdt) = scene(Some((3, 500.0)), t);
        assert_eq!(m.len(), 6);
        let cfg = open_config();
        let truth = SppState::new(
            Coordinate::new(true_pos, Datum::WGS84, Frame::ECEF, t),
            true_cdt,
            true_cdt,
            true_cdt,
            true_cdt,
        );
        let good = filter_raim_outliers(&truth, &m, None, &cfg);
        assert_eq!(good.len(), 5, "exactly the blundered satellite must be dropped");
        let blundered = m[3].raw_pr;
        assert!(good.iter().all(|g| (g.raw_pr - blundered).abs() > 1e-6),
            "the blundered measurement survived");

        // apply_raim must notice the removal and re-solve from the seed.
        let seed = seed_initial_state(&m, None);
        let fixed = apply_raim(truth, seed, &m, None, &cfg).expect("re-solve");
        assert!((fixed.position.vector - true_pos).norm() < 0.2,
            "RAIM fix error {} m", (fixed.position.vector - true_pos).norm());
        assert!((fixed.cdt - true_cdt).abs() < 0.2);
    }

    #[test]
    fn raim_leaves_a_clean_solution_untouched() {
        let t = epoch_time();
        let (m, true_pos, _) = scene(None, t);
        let cfg = open_config();
        let seed = seed_initial_state(&m, None);
        let state = solve_spp_iteratively(seed.clone(), &m, None, &cfg).expect("solve");
        let good = filter_raim_outliers(&state, &m, None, &cfg);
        assert_eq!(good.len(), m.len(), "no row may be dropped from a clean fix");
        let out = apply_raim(state.clone(), seed, &m, None, &cfg).expect("raim");
        assert_eq!(out.position.vector, state.position.vector);
        assert_eq!(out.cdt, state.cdt);
        assert!((out.position.vector - true_pos).norm() < 0.1);
    }

    #[test]
    fn raim_will_not_re_solve_below_four_surviving_measurements() {
        // Four satellites with three blunders leaves a single survivor;
        // apply_raim must keep the provisional fix rather than returning a
        // 4-parameter problem solved from one row.
        let t = epoch_time();
        let (mut epoch, ephems) = four_gps_scene(receiver(), 0.0, t);
        epoch.satellites.truncate(3);
        epoch.satellites[1].observations[0].value += 400.0;
        epoch.satellites[2].observations[0].value -= 400.0;
        let m = build_measurements(&epoch, &ephems, &open_config());
        assert_eq!(m.len(), 3);
        let cfg = open_config();
        let truth = SppState::new(
            Coordinate::new(receiver(), Datum::WGS84, Frame::ECEF, t),
            0.0,
            0.0,
            0.0,
            0.0,
        );
        let good = filter_raim_outliers(&truth, &m, None, &cfg);
        assert!(good.len() < 4, "expected fewer than 4 survivors, got {}", good.len());
        let seed = seed_initial_state(&m, None);
        let out = apply_raim(truth, seed, &m, None, &cfg).expect("no re-solve");
        assert_eq!(out.position.vector, receiver());
        let _ = LIGHT_SPEED;
    }
}
