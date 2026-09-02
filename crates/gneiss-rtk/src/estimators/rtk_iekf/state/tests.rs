#![allow(clippy::unwrap_used)]

use super::*;
use gneiss_core::time::GpsTime;
use nalgebra::Vector3;

#[test]
fn test_gradient_layout_and_roundtrip() {
    let mut state = RtkState::new(Vector3::zeros(), GpsTime::new(2000, 0.0));
    assert_eq!(state.dim(), 6);
    state.enable_zwd(0.0225);
    assert_eq!(state.dim(), 7);
    assert_eq!(state.amb_offset(), 7);
    state.enable_gradients(4e-6);
    // Layout: pos(3) vel(3) zwd@6 gN@7 gE@8 ambs@9..
    assert_eq!(state.dim(), 9);
    assert_eq!(state.amb_offset(), 9);
    assert_eq!(state.grad_idx(), Some((7, 8)));

    let k = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
    state.ensure_ambiguity(k, 5.5, 100.0);
    assert_eq!(state.get_amb_idx(&k), Some(9));

    state.grad_n_m = 0.002;
    state.grad_e_m = -0.001;
    let v = state.to_dvector();
    assert!((v[6] - state.zwd_m).abs() < 1e-15);
    assert!((v[7] - 0.002).abs() < 1e-15 && (v[8] + 0.001).abs() < 1e-15);
    assert!((v[9] - 5.5).abs() < 1e-15);

    state.grad_n_m = 0.0;
    state.update_from_dvector(&v);
    assert!((state.grad_n_m - 0.002).abs() < 1e-15);
    assert!((state.grad_e_m + 0.001).abs() < 1e-15);
}

#[test]
fn test_gradients_disabled_keeps_legacy_layout() {
    let mut state = RtkState::new(Vector3::zeros(), GpsTime::new(2000, 0.0));
    state.enable_gradients(4e-6); // without ZWD: columns at 6,7
    assert_eq!(state.dim(), 8);
    assert_eq!(state.amb_offset(), 8);
    assert_eq!(state.grad_idx(), Some((6, 7)));
}

#[test]
fn test_rtk_state_lifecycle() {
    let mut state = RtkState::new(Vector3::new(100.0, 200.0, 300.0), GpsTime::new(2000, 100.0));
    assert_eq!(state.dim(), 6);

    let k1 = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
    let k2 = DoubleDiffKey { constellation_id: 0, sat: 3, ref_sat: 1, freq_band: 1 };

    state.ensure_ambiguity(k1, 5.0, 100.0);
    state.ensure_ambiguity(k2, 10.0, 100.0);
    assert_eq!(state.dim(), 8);
    assert_eq!(state.get_amb_idx(&k1), Some(6));
    assert_eq!(state.get_amb_idx(&k2), Some(7));

    let vec = state.to_dvector();
    assert_eq!(vec[6], 5.0);
    assert_eq!(vec[7], 10.0);

    state.retain_active_ambiguities(&[k2]);
    assert_eq!(state.dim(), 7);
    assert_eq!(state.get_amb_idx(&k2), Some(6));
}

fn dd_key(sv: u16, band: u8) -> DoubleDiffKey {
    DoubleDiffKey { constellation_id: 0, sat: sv, ref_sat: 1, freq_band: band }
}

#[test]
fn test_iono_disabled_zero_impact() {
    let t = GpsTime::new(2100, 0.0);
    let st = RtkState::new(Vector3::zeros(), t);
    assert!(!st.iono_enabled);
    assert_eq!(st.ionos.len(), 0);
    assert_eq!(st.dim(), 6);
}

#[test]
fn test_iono_dim_and_offset() {
    let t = GpsTime::new(2100, 0.0);
    let mut st = RtkState::new(Vector3::zeros(), t);
    st.iono_enabled = true;
    st.ensure_ambiguity(dd_key(5, 1), 100.0, 100.0);
    st.ensure_iono(dd_key(5, 1), 4.0);
    assert_eq!(st.dim(), 8);
    assert_eq!(st.iono_offset(), st.amb_offset() + 1);
    assert_eq!(st.get_iono_idx(&dd_key(5, 1)), Some(7));
}

#[test]
fn test_iono_roundtrip_preserves_values() {
    let t = GpsTime::new(2100, 0.0);
    let mut st = RtkState::new(Vector3::zeros(), t);
    st.zwd_enabled = true;
    st.iono_enabled = true;
    st.ensure_ambiguity(dd_key(5, 1), 50.0, 100.0);
    st.ensure_iono(dd_key(5, 1), 4.0);
    st.ionos[0].1 = -1.234;
    let v = st.to_dvector();
    let mut st2 = st.clone();
    st2.update_from_dvector(&v);
    assert!((st2.ionos[0].1 - (-1.234)).abs() < 1e-12);
}

#[test]
fn test_ensure_iono_idempotent() {
    let t = GpsTime::new(2100, 0.0);
    let mut st = RtkState::new(Vector3::zeros(), t);
    st.iono_enabled = true;
    st.ensure_iono(dd_key(5, 1), 4.0);
    st.ensure_iono(dd_key(5, 1), 4.0);
    assert_eq!(st.ionos.len(), 1);
    assert_eq!(st.dim(), 7);
}

fn key(sat: u16) -> DoubleDiffKey {
    DoubleDiffKey { constellation_id: 0, sat, ref_sat: 1, freq_band: 1 }
}

#[test]
fn retain_compacts_iono_states_and_keeps_cov_consistent() {
    let t = GpsTime::new(2100, 0.0);
    let mut st = RtkState::new(Vector3::zeros(), t);
    st.iono_enabled = true;
    st.ensure_ambiguity(key(5), 10.0, 100.0);
    st.ensure_ambiguity(key(7), 20.0, 100.0);
    st.ensure_iono(key(5), 4.0);
    st.ensure_iono(key(7), 4.0);
    assert_eq!(st.dim(), st.cov.nrows(), "pre-condition");

    st.retain_active_ambiguities(&[key(5)]);

    assert_eq!(
        st.dim(),
        st.cov.nrows(),
        "dim/cov desync after retain — the storm-day crash"
    );
    assert_eq!(st.ambiguities.len(), 1);
    assert_eq!(st.ionos.len(), 1, "ionos must compact with ambs");
    assert_eq!(st.get_iono_idx(&key(5)), Some(st.iono_offset()));
}

#[test]
fn sat_iono_states_keyed_by_satellite_and_shared_across_pairs() {
    let t = GpsTime::new(2100, 0.0);
    let mut st = RtkState::new(Vector3::zeros(), t);
    st.sat_iono_enabled = true;
    st.ensure_sat_iono(5);
    st.ensure_sat_iono(7);
    st.ensure_sat_iono(9);
    assert_eq!(st.dim(), 6 + 3);
    assert_eq!(st.get_sat_iono_idx(9), Some(8));
    assert_ne!(st.get_sat_iono_idx(5), st.get_sat_iono_idx(7));
}

#[test]
fn sat_iono_multi_constellation_keys() {
    let t = GpsTime::new(2100, 0.0);
    let mut st = RtkState::new(Vector3::zeros(), t);
    st.sat_iono_enabled = true;
    st.ensure_sat_iono_key(0, 10);
    st.ensure_sat_iono_key(2, 10);
    assert_eq!(st.dim(), 6 + 2);
    let idx_gps = st.get_sat_iono_key_idx(0, 10);
    let idx_gal = st.get_sat_iono_key_idx(2, 10);
    assert!(idx_gps.is_some());
    assert!(idx_gal.is_some());
    assert_ne!(idx_gps, idx_gal);
}

#[test]
fn test_transfer_reference_satellite_math() {
    let t = GpsTime::new(2100, 0.0);
    let mut st = RtkState::new(Vector3::zeros(), t);
    let k21 = DoubleDiffKey { constellation_id: 0, sat: 2, ref_sat: 1, freq_band: 1 };
    let k31 = DoubleDiffKey { constellation_id: 0, sat: 3, ref_sat: 1, freq_band: 1 };
    st.ensure_ambiguity(k21, 5.0, 0.04);
    st.ensure_ambiguity(k31, 8.0, 0.04);

    // Transfer reference satellite from 1 to 2
    let transferred = st.transfer_reference_satellite(0, 1, 1, 2);
    assert!(transferred);

    let k12 = DoubleDiffKey { constellation_id: 0, sat: 1, ref_sat: 2, freq_band: 1 };
    let k32 = DoubleDiffKey { constellation_id: 0, sat: 3, ref_sat: 2, freq_band: 1 };

    let idx12 = st.get_amb_idx(&k12).unwrap();
    let idx32 = st.get_amb_idx(&k32).unwrap();

    let x = st.to_dvector();
    assert_eq!(x[idx12], -5.0, "N_12 = -N_21 = -5.0");
    assert_eq!(x[idx32], 3.0, "N_32 = N_31 - N_21 = 8.0 - 5.0 = 3.0");
    assert!(st.cov[(idx12, idx12)] > 0.0);
    assert!(st.cov[(idx32, idx32)] > 0.0);
}
