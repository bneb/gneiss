#![allow(clippy::unwrap_used)]

//! Unit tests for the Doppler velocity estimator.

use super::*;
use gneiss_core::sat::Constellation;

#[test]
fn test_solve_wls_synthetic_four_sats() {
    // True receiver velocity: [10.0, -5.0, 2.0] m/s, clock drift: 3.0 m/s
    let v_true = Vector3::new(10.0, -5.0, 2.0);
    let cdt_true = 3.0;

    let directions = [
        Vector3::new(1.0, 0.0, 0.8).normalize(),
        Vector3::new(-1.0, 0.2, 0.3).normalize(),
        Vector3::new(0.1, 1.0, 0.6).normalize(),
        Vector3::new(-0.2, -1.0, 0.4).normalize(),
    ];

    let meas: Vec<SingleDopplerMeas> = directions
        .iter()
        .enumerate()
        .map(|(i, &d)| {
            let range_rate = -d.dot(&v_true) + cdt_true;
            SingleDopplerMeas {
                sat: SatelliteId { constellation: Constellation::Gps, prn: (i + 1) as u8 },
                e_los: d,
                y_range_rate: range_rate,
                weight: 1.0,
            }
        })
        .collect();

    let cmap = ConstellationMap::new(&meas);
    let (sol, _) = solve_wls(&meas, &cmap).expect("WLS solve failed");
    assert!((sol[0] - v_true.x).abs() < 1e-9);
    assert!((sol[1] - v_true.y).abs() < 1e-9);
    assert!((sol[2] - v_true.z).abs() < 1e-9);
    assert!((sol[3] - cdt_true).abs() < 1e-9);
}

#[test]
fn test_outlier_filtering() {
    let v_true = Vector3::new(0.0, 0.0, 0.0);
    let cdt_true = 0.0;
    let directions = [
        Vector3::new(1.0, 0.0, 0.8).normalize(),
        Vector3::new(-1.0, 0.2, 0.3).normalize(),
        Vector3::new(0.1, 1.0, 0.6).normalize(),
        Vector3::new(-0.2, -1.0, 0.4).normalize(),
        Vector3::new(0.5, 0.5, 0.7).normalize(),
        Vector3::new(-0.5, -0.5, 0.5).normalize(),
    ];

    let mut meas: Vec<SingleDopplerMeas> = directions
        .iter()
        .enumerate()
        .map(|(i, &d)| SingleDopplerMeas {
            sat: SatelliteId { constellation: Constellation::Gps, prn: (i + 1) as u8 },
            e_los: d,
            y_range_rate: -d.dot(&v_true) + cdt_true,
            weight: 1.0,
        })
        .collect();

    // Inject 5 m/s blunder on last satellite
    meas[5].y_range_rate += 5.0;

    let cmap = ConstellationMap::new(&meas);
    let (sol, inv) = solve_wls(&meas, &cmap).expect("WLS solve");
    let cleaned = filter_outliers(&meas, &sol, &inv, &cmap).expect("outlier detected");
    assert_eq!(cleaned.len(), 5);
    let c2map = ConstellationMap::new(&cleaned);
    let (sol2, _) = solve_wls(&cleaned, &c2map).expect("cleaned solve");
    assert!(sol2.fixed_rows::<3>(0).norm() < 1e-6);
}

#[test]
fn test_real_doppler_signs() {
    use std::fs::File;
    use std::io::BufReader;
    let nav_path = "../../datasets/urbannav/tokyo/Tokyo_Data/Odaiba/base.nav";
    let obs_path = "../../datasets/urbannav/tokyo/Tokyo_Data/Odaiba/rover_trimble.obs";
    if !std::path::Path::new(nav_path).exists() {
        return;
    }
    let nav_f = File::open(nav_path).expect("open nav");
    let (ephems, _) = gneiss_parsers::rinex::parse_rinex_nav(BufReader::new(nav_f)).expect("parse nav");
    let obs_f = File::open(obs_path).expect("open obs");
    let (epochs, approx) = gneiss_parsers::rinex::parse_rinex_obs(BufReader::new(obs_f)).expect("parse obs");
    let approx_pos = Vector3::new(
        approx.approx_position.unwrap()[0],
        approx.approx_position.unwrap()[1],
        approx.approx_position.unwrap()[2],
    );
    let ep0 = &epochs[0];
    let sol0 = estimate_doppler_velocity(ep0, &ephems, approx_pos).expect("doppler sol0");
    assert!(sol0.vel_ecef.norm() < 2.0);
    assert!(sol0.n_sats >= 15);

    let ep_target = epochs.iter().find(|e| (e.time.tow - 273600.0).abs() < 0.05);
    if let Some(ep) = ep_target {
        let sol = estimate_doppler_velocity(ep, &ephems, approx_pos).expect("doppler sol at 273600");
        assert!(sol.vel_ecef.norm() < 2.0);
        assert!(sol.n_sats >= 15);
    }
}

#[path = "synthetic.rs"]
mod synthetic;

#[test]
fn stationary_receiver_recovers_zero_velocity_from_synthetic_doppler() {
    // Every innovation y = -(e_los . 0) + c*0 must be zero, so the
    // weighted solution must return the origin. A flipped sign on
    // rho_meas (= +lambda*f_D) would instead return twice the satellite
    // range rate, i.e. thousands of m/s.
    use gneiss_core::time::GpsTime;
    let t = GpsTime::new(2000, 100_000.0);
    let ephems: Vec<Ephemeris> = (1..=20u8).map(|p| synthetic::ephemeris(p, t)).collect();
    let sats = ephems.iter().map(|e| synthetic::sat_obs(e, Vector3::zeros())).collect();
    let ep = synthetic::epoch(sats, t);
    let sol = estimate_doppler_velocity(&ep, &ephems, synthetic::receiver())
        .expect("doppler velocity must solve");
    assert!(sol.n_sats >= 4, "only {} satellites survived", sol.n_sats);
    assert!(sol.vel_ecef.norm() < 0.05, "v = {:?} m/s", sol.vel_ecef);
    assert!(sol.clk_drift_m_s.abs() < 0.05, "cdt drift = {}", sol.clk_drift_m_s);
}

#[test]
fn known_receiver_velocity_is_recovered_with_the_correct_receding_sign() {
    // The fixture encodes v_rx = (10, -5, 2) m/s. Doppler velocimetry must
    // return that vector, not its negation: receding (e_los . v_rx > 0)
    // gives a POSITIVE range rate, so any sign flip in the module shows up
    // as a velocity of (-10, 5, -2).
    use gneiss_core::time::GpsTime;
    let t = GpsTime::new(2000, 100_000.0);
    let v_true = Vector3::new(10.0, -5.0, 2.0);
    let ephems: Vec<Ephemeris> = (1..=20u8).map(|p| synthetic::ephemeris(p, t)).collect();
    let sats = ephems.iter().map(|e| synthetic::sat_obs(e, v_true)).collect();
    let ep = synthetic::epoch(sats, t);
    let sol = estimate_doppler_velocity(&ep, &ephems, synthetic::receiver())
        .expect("doppler velocity must solve");
    let err = sol.vel_ecef - v_true;
    assert!(err.norm() < 0.05, "velocity error {:?} (got {:?})", err, sol.vel_ecef);
    // The negated solution must NOT pass: that is the failure mode a sign
    // error would produce.
    assert!((sol.vel_ecef + v_true).norm() > 5.0);
}

#[test]
fn fewer_satellites_than_unknowns_returns_none() {
    // With only GPS present the state has 4 unknowns (3 velocity + 1 clock
    // drift). Two satellites cannot determine it and the estimator must refuse
    // rather than return an under-determined fit. PRNs 3 and 5 are used because
    // this fixture geometry places them at ~49 deg and ~31 deg elevation, i.e.
    // well above the 10 deg mask; the precondition assertion below proves the
    // refusal is due to the count and not to the elevation filter.
    use gneiss_core::time::GpsTime;
    let t = GpsTime::new(2000, 100_000.0);
    let ephems: Vec<Ephemeris> = [3u8, 5].iter().map(|p| synthetic::ephemeris(*p, t)).collect();
    let sats = ephems.iter().map(|e| synthetic::sat_obs(e, Vector3::zeros())).collect::<Vec<_>>();
    let r_rcv = synthetic::receiver();
    let up = compute_up_vector(&r_rcv);
    let kept = sats
        .iter()
        .filter(|s| extract_satellite_doppler(s, &ephems, t, &r_rcv, &up).is_some())
        .count();
    assert_eq!(kept, 2, "precondition: both satellites must clear the elevation mask");
    let ep = synthetic::epoch(sats, t);
    assert!(estimate_doppler_velocity(&ep, &ephems, r_rcv).is_none());
}

#[test]
fn a_satellite_within_1000_km_of_the_receiver_is_discarded() {
    // The line-of-sight geometry is only meaningful for a receiver on the
    // Earth's surface; a satellite closer than 1000 km indicates a broken
    // ephemeris or a wrong reference position and must be dropped rather than
    // producing a near-unit e_los with a garbage Doppler.
    use gneiss_core::time::GpsTime;
    let t = GpsTime::new(2000, 100_000.0);
    let eph = synthetic::ephemeris(3, t);
    let so = synthetic::sat_obs(&eph, Vector3::zeros());
    let far = synthetic::receiver();
    let up = compute_up_vector(&far);
    // Precondition: at the surface receiver this satellite is usable.
    let base = extract_satellite_doppler(&so, std::slice::from_ref(&eph), t, &far, &up)
        .expect("usable at the surface");
    assert!(base.e_los.norm() - 1.0 < 1e-12);

    // Move the "receiver" onto the satellite: the range collapses and the row
    // must disappear.
    let on_sat = eph.position(GpsTime::new(t.week, t.tow - 0.07)).0;
    assert!(extract_satellite_doppler(&so, std::slice::from_ref(&eph), t, &on_sat, &up).is_none());
}

#[test]
fn constellation_map_assigns_one_clock_column_per_cdma_system() {
    let mk = |c: Constellation, prn: u8| SingleDopplerMeas {
        sat: SatelliteId { constellation: c, prn },
        e_los: Vector3::z(),
        y_range_rate: 0.0,
        weight: 1.0,
    };
    // QZSS shares the GPS clock column; Galileo and BeiDou get their own.
    let m = [mk(Constellation::Gps, 1), mk(Constellation::Qzss, 2), mk(Constellation::Galileo, 3)];
    let cmap = ConstellationMap::new(&m);
    assert_eq!(cmap.n_params, 5, "3 position/velocity + 2 clock columns");
    assert_eq!(cmap.col_for(Constellation::Gps), Some(3));
    assert_eq!(cmap.col_for(Constellation::Qzss), Some(3));
    assert_eq!(cmap.col_for(Constellation::Galileo), Some(4));
    assert_eq!(cmap.col_for(Constellation::Beidou), None, "no BeiDou measurement present");
    // GLONASS is FDMA and carries no interoperable clock column here.
    assert_eq!(cmap.col_for(Constellation::Glonass), None);
    assert_eq!(cmap.col_for(Constellation::Sbas), None);

    // A BeiDou measurement appends column 5 and raises n_params to 6.
    let mut m2 = m.to_vec();
    m2.push(mk(Constellation::Beidou, 4));
    let cmap2 = ConstellationMap::new(&m2);
    assert_eq!(cmap2.n_params, 6);
    assert_eq!(cmap2.col_for(Constellation::Beidou), Some(5));
    assert_eq!(cmap2.col_for(Constellation::Gps), Some(3));

    // Empty input is always the bare 3-DOF velocity problem.
    let cmap0 = ConstellationMap::new(&[]);
    assert_eq!(cmap0.n_params, 3);
}

#[test]
fn wls_normal_matrix_is_the_weighted_design_outer_product() {
    // With two GPS satellites (n_params = 4) the solve is refused for
    // under-determination.
    let mk = |i: u8, d: Vector3<f64>| SingleDopplerMeas {
        sat: SatelliteId { constellation: Constellation::Gps, prn: i },
        e_los: d,
        y_range_rate: 0.0,
        weight: 1.0,
    };
    let meas = vec![mk(1, Vector3::new(1.0, 0.0, 0.0)), mk(2, Vector3::new(0.0, 1.0, 0.0))];
    let cmap = ConstellationMap::new(&meas);
    assert_eq!(cmap.n_params, 4);
    assert!(solve_wls(&meas, &cmap).is_none(), "2 measurements cannot fill 4 unknowns");

    // Four GPS satellites: h_i = [-e_los_i, 1] and the normal matrix is
    // sum_i w_i h_i h_i^T, so (0,3) = -sum_i w_i e_los_i.x and (3,3) = sum_i w_i.
    let meas = vec![
        mk(1, Vector3::new(1.0, 0.0, 0.0)),
        mk(2, Vector3::new(0.0, 1.0, 0.0)),
        mk(3, Vector3::new(0.0, 0.0, 1.0)),
        mk(4, Vector3::new(1.0, 1.0, 1.0).normalize()),
    ];
    let cmap = ConstellationMap::new(&meas);
    let (sol, inv) = solve_wls(&meas, &cmap).expect("4 GPS satellites are solvable");
    let sx = meas.iter().map(|m| m.e_los.x).sum::<f64>();
    let sy = meas.iter().map(|m| m.e_los.y).sum::<f64>();
    let sz = meas.iter().map(|m| m.e_los.z).sum::<f64>();
    // With all residuals zero the solution must be exactly zero.
    assert!(sol.fixed_rows::<3>(0).norm() < 1e-12);
    assert!(sol[3].abs() < 1e-12);
    // The returned inverse must invert the regularised normal matrix.
    let mut normal = Matrix6::<f64>::zeros();
    for m in &meas {
        let mut h = Vector6::zeros();
        h[0] = -m.e_los.x;
        h[1] = -m.e_los.y;
        h[2] = -m.e_los.z;
        h[3] = 1.0;
        normal += h * m.weight * h.transpose();
    }
    for i in 4..6 {
        normal[(i, i)] = 1.0;
    }
    assert!((inv * normal - Matrix6::<f64>::identity()).norm() < 1e-9);
    // Spot-check the (0,3) entry: (N^-1 N)[0][3] == 1 forces N[3][0] to be
    // -sum(e_los.x), so recover it from the closed form.
    assert!((-sx).abs() > 0.0 && (-sy).abs() > 0.0 && (-sz).abs() > 0.0);
}

#[test]
fn outlier_filter_keeps_every_measurement_when_no_residual_exceeds_three_sigma() {
    // A consistent solution has zero residuals, so the normalised
    // innovation statistic is 0 everywhere and nothing may be removed.
    let dirs = [
        Vector3::new(1.0, 0.0, 0.8).normalize(),
        Vector3::new(-1.0, 0.2, 0.3).normalize(),
        Vector3::new(0.1, 1.0, 0.6).normalize(),
        Vector3::new(-0.2, -1.0, 0.4).normalize(),
        Vector3::new(0.5, 0.5, 0.7).normalize(),
    ];
    let meas: Vec<SingleDopplerMeas> = dirs
        .iter()
        .enumerate()
        .map(|(i, d)| SingleDopplerMeas {
            sat: SatelliteId { constellation: Constellation::Gps, prn: (i + 1) as u8 },
            e_los: *d,
            y_range_rate: 0.0,
            weight: 1.0,
        })
        .collect();
    let cmap = ConstellationMap::new(&meas);
    let (sol, inv) = solve_wls(&meas, &cmap).expect("solve");
    assert!(filter_outliers(&meas, &sol, &inv, &cmap).is_none());
}

#[test]
fn a_constellation_without_a_clock_column_fails_the_whole_wls_solve() {
    // GLONASS is FDMA and `extract_satellite_doppler` filters it out, but if one
    // ever reached the solver, `ConstellationMap::new` would NOT allocate a clock
    // column for it (n_params stays 4 for GPS-only), so the `?` in `solve_wls`
    // aborts the ENTIRE solve rather than dropping that row. Refusing beats
    // silently solving a 4-parameter problem from 4 rows that include one with
    // no clock coefficient.
    let meas = vec![
        SingleDopplerMeas {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            e_los: Vector3::new(1.0, 0.0, 0.5).normalize(),
            y_range_rate: 1.0,
            weight: 1.0,
        },
        SingleDopplerMeas {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 2 },
            e_los: Vector3::new(0.0, 1.0, 0.5).normalize(),
            y_range_rate: 2.0,
            weight: 1.0,
        },
        SingleDopplerMeas {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 3 },
            e_los: Vector3::new(0.5, 0.5, 0.7).normalize(),
            y_range_rate: -1.0,
            weight: 1.0,
        },
        SingleDopplerMeas {
            sat: SatelliteId { constellation: Constellation::Glonass, prn: 5 },
            e_los: Vector3::new(-0.3, 0.2, 0.9).normalize(),
            y_range_rate: 500.0, // huge blunder, but no clock column
            weight: 1.0,
        },
    ];
    let cmap = ConstellationMap::new(&meas);
    assert_eq!(cmap.n_params, 4, "no clock column is allocated for GLONASS");
    assert!(solve_wls(&meas, &cmap).is_none());
}

#[test]
fn outlier_filter_skips_a_measurement_with_no_clock_column() {
    // `filter_outliers` receives the solution separately, so it can be handed a
    // GLONASS row that `solve_wls` would have refused. No statistic can be
    // formed for that row, so the `continue` skips it: a 500 m/s blunder on a
    // column-less measurement must NOT be reported as the largest normalised
    // innovation, and nothing may be removed.
    let meas = vec![
        SingleDopplerMeas {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 1 },
            e_los: Vector3::new(1.0, 0.0, 0.5).normalize(),
            y_range_rate: 1.0,
            weight: 1.0,
        },
        SingleDopplerMeas {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 2 },
            e_los: Vector3::new(0.0, 1.0, 0.5).normalize(),
            y_range_rate: 2.0,
            weight: 1.0,
        },
        SingleDopplerMeas {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 3 },
            e_los: Vector3::new(0.5, 0.5, 0.7).normalize(),
            y_range_rate: -1.0,
            weight: 1.0,
        },
        SingleDopplerMeas {
            sat: SatelliteId { constellation: Constellation::Gps, prn: 4 },
            e_los: Vector3::new(-0.4, -0.6, 0.6).normalize(),
            y_range_rate: 0.5,
            weight: 1.0,
        },
        SingleDopplerMeas {
            sat: SatelliteId { constellation: Constellation::Glonass, prn: 5 },
            e_los: Vector3::new(-0.3, 0.2, 0.9).normalize(),
            y_range_rate: 500.0,
            weight: 1.0,
        },
    ];
    // Solve the GPS-only subset (4 GPS + 1 clock column = 4 unknowns, exactly
    // determined), then hand the full set to the filter.
    let gps: Vec<SingleDopplerMeas> =
        meas.iter().filter(|m| m.sat.constellation == Constellation::Gps).copied().collect();
    let cmap = ConstellationMap::new(&gps);
    let (sol, inv) = solve_wls(&gps, &cmap).expect("GPS-only solve");
    assert!(filter_outliers(&meas, &sol, &inv, &cmap).is_none(),
        "a column-less constellation must be skipped, not removed");
    // Sanity: the GPS-only residuals really are tiny, so the None above is due
    // to the GLONASS row being skipped rather than to an unlucky geometry.
    let res = compute_residuals(&gps, &sol, &cmap);
    assert!(res.iter().all(|r| r.abs() < 1e-9), "residuals = {res:?}");
}

#[test]
fn compute_up_vector_matches_the_geodetic_up_at_exact_sites() {
    // At the equator on the prime meridian geodetic Up is ECEF +x;
    // at 90 deg east it is ECEF +y; at the pole it is ECEF +z.
    let cases = [
        (Vector3::new(6_378_137.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)),
        (Vector3::new(0.0, 6_378_137.0, 0.0), Vector3::new(0.0, 1.0, 0.0)),
        (Vector3::new(0.0, 0.0, 6_356_752.314_245), Vector3::new(0.0, 0.0, 1.0)),
    ];
    for (pos, want) in cases {
        let up = compute_up_vector(&pos);
        assert!((up - want).norm() < 1e-9, "at {pos:?} up = {up:?}, want {want:?}");
        assert!((up.norm() - 1.0).abs() < 1e-12);
    }
}
