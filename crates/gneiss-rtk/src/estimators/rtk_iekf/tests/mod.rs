#![allow(clippy::unwrap_used)]

use std::sync::Arc;
use super::*;
use super::ar_gate::SEED_VARIANCE_SAFETY_MARGIN;
use crate::sim::generator::{generate_simulation_dataset, SimulationConfig};
use gneiss_core::constants::SPEED_OF_LIGHT_M_S;

mod test_ar_gate;
mod test_clk;
mod test_quad_const;

fn test_engine(start: GpsTime) -> GnssRtkIekf {
    GnssRtkIekf::new(Vector3::new(1.0, 2.0, 3.0), start, 1.0)
}

fn dd_key(sat: u16, band: u8) -> DoubleDiffKey {
    DoubleDiffKey {
        constellation_id: 0,
        sat,
        ref_sat: 1,
        freq_band: band,
    }
}

#[test]
fn seed_variance_matches_margin_squared_times_code_sigma_in_cycles() {
    let pr_var_m2: f64 = 0.64;
    let lambda = 0.1903;
    let sigma_cycles = pr_var_m2.sqrt() / lambda;
    let expected = (SEED_VARIANCE_SAFETY_MARGIN * sigma_cycles).powi(2);
    assert!((seed_ambiguity_variance_cycles2(pr_var_m2, lambda) - expected).abs() < 1e-9);
}

#[test]
fn seed_variance_is_looser_for_noisier_low_elevation_code() {
    let lambda = 0.1903;
    let high_el_pr_var = 0.05;
    let low_el_pr_var = 2.4;
    let high = seed_ambiguity_variance_cycles2(high_el_pr_var, lambda);
    let low = seed_ambiguity_variance_cycles2(low_el_pr_var, lambda);
    assert!(low > high, "low-elevation seed variance ({low}) should exceed high-elevation ({high})");
}

fn run_sim(
    engine: &mut GnssRtkIekf,
    sim: &crate::sim::generator::SimulationDataset,
    base: Vector3<f64>,
) -> Vec<(bool, Vector3<f64>)> {
    let mut out = Vec::new();
    for i in 0..sim.rover_epochs.len() {
        let s = engine
            .process_epoch(&sim.rover_epochs[i], &sim.base_epochs[i], base, &sim.ephemerides)
            .expect("epoch must process");
        out.push((s.is_fixed, s.position_ecef));
    }
    out
}

#[test]
fn test_gnss_rtk_iekf_runs_on_simulated_dataset() {
    let cfg = SimulationConfig {
        duration_s: 10.0,
        ..Default::default()
    };
    let sim = generate_simulation_dataset(&cfg);
    let mut engine = GnssRtkIekf::new(cfg.base_ecef, sim.rover_epochs[0].time, 1.0);

    let mut fixed_count = 0;
    for i in 0..sim.rover_epochs.len() {
        let sol = engine.process_epoch(
            &sim.rover_epochs[i],
            &sim.base_epochs[i],
            cfg.base_ecef,
            &sim.ephemerides,
        );
        assert!(sol.is_ok());
        let s = sol.unwrap();
        let err = (s.position_ecef - sim.truth_positions[i].1).norm();
        if s.is_fixed {
            fixed_count += 1;
            assert!(err < 0.05, "Fixed epoch error should be < 5cm, got {:.4}m", err);
        }
    }

    assert!(fixed_count >= 5, "RTK engine should fix at least 5 epochs");
    let smoothed = engine.smooth();
    assert_eq!(smoothed.len(), 10);
}

#[test]
fn test_select_constellations_drops_glonass_by_default() {
    use gneiss_core::sat::{Constellation, SatelliteId};
    let mk = |c: Constellation, prn: u8| {
        (SatelliteId { constellation: c, prn }, Vector3::zeros())
    };
    let sat_info = vec![
        mk(Constellation::Gps, 6u8),
        mk(Constellation::Glonass, 8u8),
        mk(Constellation::Galileo, 3u8),
        mk(Constellation::Gps, 12u8),
    ];
    let got = GnssRtkIekf::select_constellations(&sat_info, false);
    assert_eq!(got, vec![
        Constellation::Gps as u8,
        Constellation::Galileo as u8,
    ]);
}

#[test]
fn receiver_dd_pcv_m_requires_loaded_pair() {
    let Ok(db) = gneiss_parsers::antex::AntexDatabase::parse("../../datasets/igs14.atx") else {
        return;
    };
    use gneiss_parsers::receiver_antenna::{compute_dd_pcv_correction_2d, ReceiverAntenna};
    let trm = Arc::new(ReceiverAntenna::lookup(&db, "TRM59800.00", "SCIT").expect("igs14 TRM"));
    let ash = Arc::new(ReceiverAntenna::lookup(&db, "ASH701945B_M", "SCIT").expect("igs14 ASH"));
    let mut eng = GnssRtkIekf::new(
        Vector3::new(-3961904.43, 3348994.27, 3698211.71),
        GpsTime::new(2000, 100.0),
        1.0,
    );
    let sat_pos = eng.state.pos_ecef + Vector3::new(1.0e7, 5.0e6, 2.0e7);
    let ref_pos = eng.state.pos_ecef + Vector3::new(0.0, 0.0, 2.4e7);
    let sid = gneiss_core::sat::SatelliteId {
        constellation: gneiss_core::sat::Constellation::Gps,
        prn: 3,
    };

    eng.receiver_pcv = None;
    assert_eq!(eng.receiver_dd_pcv_m(sid, 1, sat_pos, ref_pos), 0.0);

    eng.receiver_pcv = Some((trm.clone(), ash.clone()));
    let dd = eng.receiver_dd_pcv_m(sid, 1, sat_pos, ref_pos);
    let llh = gneiss_core::coords::ecef_to_llh(eng.state.pos_ecef);
    let (az_s, el_s) = gneiss_core::coords::az_el(llh, eng.state.pos_ecef, sat_pos);
    let (az_r, el_r) = gneiss_core::coords::az_el(llh, eng.state.pos_ecef, ref_pos);
    let expected = compute_dd_pcv_correction_2d(
        &trm,
        &ash,
        "G01",
        az_s,
        el_s,
        az_r,
        el_r,
        eng.rover_heading_rad,
    );
    assert!(dd.abs() > 1e-6, "cross-family correction must be non-zero: {dd}");
    assert!((dd - expected).abs() < 1e-12, "dd={dd} expected={expected}");
}
