//! Tests for the iono-free fixed-position update.
//! Child module of `iono_free`, matching the existing `mw` and `formation`
//! splits, so the parent stays under the 500-line budget.

use super::*;
use gneiss_core::obs::{ObsCode, ObsType, Observation, SignalCode};
use gneiss_core::time::GpsTime;

    const F1: f64 = 1575.42e6;
    const F2: f64 = 1227.60e6;

    fn make_key(sat: u16) -> DoubleDiffKey {
        DoubleDiffKey { constellation_id: 0, sat, ref_sat: 1, freq_band: 1 }
    }

    fn phase_obs(sat: SatelliteId, band: u8, value: f64) -> SatObs {
        let code = ObsCode {
            obs_type: ObsType::CarrierPhase,
            signal: SignalCode { freq_band: band, attribute: 'X' },
        };
        let o = Observation { code, value, lock_time: None, lli: Some(0) };
        SatObs { sat, observations: vec![o] }
    }

    /// Build one station's triple: band 1, band 6, band 5.
    fn station(sat: SatelliteId, phase: f64) -> SatObs {
        let mut s = phase_obs(sat, 1, phase);
        for (band, scale) in [(6u8, 1.02), (5, 1.30)] {
            let extra = phase_obs(sat, band, phase * scale).observations.remove(0);
            s.observations.push(extra);
        }
        s
    }

    /// Fixed geometry shared by the band-selection tests.
    fn combine(rs: &SatObs, bs: &SatObs, rr: &SatObs, br: &SatObs) -> Option<IonoFreeMeasurement> {
        form_iono_free_dd(
            rs.sat, rs, bs, rr, br,
            Vector3::new(10_000.0, 20_000.0, 20_000.0),
            Vector3::new(30_000.0, 5_000.0, 10_000.0),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(100.0, 200.0, 300.0),
            make_key(3),
            0,
        )
    }

    /// A degenerate first candidate must not abort the whole measurement.
    ///
    /// GPS `C6X`/L6C (RINEX 3 band 6) is present on all four stations here, so
    /// band 6 wins the preference order — but the legacy table resolves GPS
    /// band 6 to L1, i.e. exactly `f1`, so `|f1 - f2| = 0` and the iono-free
    /// combination would divide by zero. Band 5 (L5 @1176.45 MHz) is observed
    /// and perfectly usable.
    ///
    /// Correct behaviour: skip the degenerate candidate and use band 5. The
    /// code returned `None` instead, discarding a pair it could have combined.
    #[test]
    fn degenerate_secondary_band_falls_through_to_the_next_candidate() {
        let g = |prn: u8| SatelliteId { constellation: Constellation::Gps, prn };
        let s = |prn, phase| station(g(prn), phase);
        let m = combine(&s(3, 2.05e8), &s(9, 2.05e8 - 7.0e3), &s(4, 2.05e8 + 1.0e4), &s(8, 2.05e8 + 2.0e4))
            .expect("a usable secondary band exists; the pair must not be discarded");
        assert_eq!(m.b2, 5, "must fall through band 6 to band 5");
        assert_eq!(m.f1_hz, F1);
        // L5 = 1150 x 1.023 MHz; exact integer identity, no tolerance.
        assert_eq!(m.f2_hz, 1_176_450_000.0);
    }

    /// The same selection must still prefer the earliest candidate when that
    /// candidate is usable — the fall-through must not reorder the preference.
    #[test]
    fn usable_first_candidate_is_still_preferred() {
        let g = |prn: u8| SatelliteId { constellation: Constellation::Gps, prn };
        let s = |prn| {
            let mut s = phase_obs(g(prn), 1, 2.05e8);
            for b in [2u8, 6, 5] {
                let extra = phase_obs(g(prn), b, 2.05e8 + f64::from(b) * 1.0e5).observations.remove(0);
                s.observations.push(extra);
            }
            s
        };
        let m = combine(&s(3), &s(9), &s(4), &s(8)).expect("band 2 is observed and usable");
        assert_eq!(m.b2, 2);
        assert_eq!(m.f2_hz, F2);
    }

    #[test]
    fn test_iono_free_cancels_ionosphere_exactly() {
        let range = 23_000_000.0;
        let n1 = 5.0;
        let n2 = -3.0;
        let iono_m = 0.25; // L1 path iono delay (m)
        // Phase iono is an advance: -iono/λ on each band, L2 scaled by (f1/f2)^2.
        let phi1 = range * F1 / SPEED_OF_LIGHT_M_S + n1 - iono_m * F1 / SPEED_OF_LIGHT_M_S;
        let phi2 = range * F2 / SPEED_OF_LIGHT_M_S + n2 - iono_m * (F1 / F2).powi(2) * F2 / SPEED_OF_LIGHT_M_S;

        let phi_if = combine_iono_free(F1, F2, phi1, phi2);
        let n_if = (F1 * n1 - F2 * n2) / (F1 - F2);
        let expect = range / lambda_iono_free(F1, F2) + n_if;
        assert!((phi_if - expect).abs() < 1e-9, "IF combination must cancel iono, got {:.12e}", phi_if - expect);
    }

    #[test]
    fn test_fixed_iono_free_position_recovers_truth_with_iono() {
        let true_pos = Vector3::new(100.0, 200.0, 300.0);
        let base_pos = Vector3::new(0.0, 0.0, 0.0);
        let state = RtkState::new(true_pos + Vector3::new(0.8, -0.4, 0.2), GpsTime::new(2000, 100.0));

        let sats: [(Vector3<f64>, f64, f64, f64); 6] = [
            (Vector3::new(10_000.0, 20_000.0, 20_000.0), 1.0, 2.0, 0.30),
            (Vector3::new(25_000.0, 5_000.0, 18_000.0), 3.0, -1.0, 0.22),
            (Vector3::new(8_000.0, 30_000.0, 12_000.0), -2.0, 4.0, 0.18),
            (Vector3::new(20_000.0, 12_000.0, 25_000.0), 2.0, 0.0, 0.26),
            (Vector3::new(15_000.0, 22_000.0, 15_000.0), 0.0, -2.0, 0.20),
            (Vector3::new(22_000.0, 18_000.0, 22_000.0), 4.0, 1.0, 0.24),
        ];
        let ref_pos = Vector3::new(30_000.0, 5_000.0, 10_000.0);

        let mut meas = Vec::new();
        for (i, (sat_pos, n1, n2, iono_m)) in sats.iter().enumerate() {
            let base_dd = (sat_pos - base_pos).norm() - (ref_pos - base_pos).norm();
            let geom = (sat_pos - true_pos).norm() - (ref_pos - true_pos).norm() - base_dd;
            let lambda1 = SPEED_OF_LIGHT_M_S / F1;
            let lambda2 = SPEED_OF_LIGHT_M_S / F2;
            let iono_l2 = iono_m * (F1 / F2).powi(2);
            let dd1 = geom / lambda1 + n1 - iono_m / lambda1;
            let dd2 = geom / lambda2 + n2 - iono_l2 / lambda2;
            meas.push(IonoFreeMeasurement {
                key: make_key(2 + i as u16),
                sat_pos: *sat_pos,
                ref_pos,
                base_pos,
                lambda_if: lambda_iono_free(F1, F2),
                b2: 2,
                f1_hz: F1,
                f2_hz: F2,
                dd_phase_if_cycles: combine_iono_free(F1, F2, dd1, dd2),
                variance_cycles2: 1e-4,
                baseline_m: (base_pos - true_pos).norm() + 1000.0,
                dgrad_n_rov: 0.0,
                dgrad_e_rov: 0.0,
            });
        }
        let fixed: Vec<(DoubleDiffKey, f64)> = sats.iter().enumerate()
            .flat_map(|(i, (_, n1, n2, _))| vec![(make_key(2 + i as u16), *n1), (DoubleDiffKey { freq_band: 2, ..make_key(2 + i as u16) }, *n2)])
            .collect();
        let ar = ArResult {
            position_ecef: true_pos,
            cov_position: Matrix3::identity(),
            ratio: 3.0,
            is_fixed: true,
            num_ambiguities: 10,
            fixed_ambiguities: fixed,
        };

        let (pos, _) = match apply_fixed_iono_free(&state, &meas, &ar) {
            IonoFreeOutcome::Solution(p, c) => (p, c),
            _ => panic!("IF stage should produce a solution"),
        };
        let err = (pos - true_pos).norm();
        assert!(err < 0.01, "IF fixed position should recover truth under iono, got {:.4}m", err);
    }

    #[test]
    fn test_rejected_when_integers_shifted_by_common_mode_slip() {
        // Same fixture as above but every N1 and N2 shifted +1: the common
        // mode cancels in N1-N2 (widelane-blind) yet biases the iono-free
        // combination by a full cycle, which the residual gate must reject.
        let true_pos = Vector3::new(100.0, 200.0, 300.0);
        let base_pos = Vector3::new(0.0, 0.0, 0.0);
        let state = RtkState::new(true_pos, GpsTime::new(2000, 100.0));

        let sats: [(Vector3<f64>, f64, f64); 6] = [
            (Vector3::new(10_000.0, 20_000.0, 20_000.0), 1.0, 2.0),
            (Vector3::new(25_000.0, 5_000.0, 18_000.0), 3.0, -1.0),
            (Vector3::new(8_000.0, 30_000.0, 12_000.0), -2.0, 4.0),
            (Vector3::new(20_000.0, 12_000.0, 25_000.0), 2.0, 0.0),
            (Vector3::new(15_000.0, 22_000.0, 15_000.0), 0.0, -2.0),
            (Vector3::new(22_000.0, 18_000.0, 22_000.0), 4.0, 1.0),
        ];
        let ref_pos = Vector3::new(30_000.0, 5_000.0, 10_000.0);
        let lambda1 = SPEED_OF_LIGHT_M_S / F1;
        let lambda2 = SPEED_OF_LIGHT_M_S / F2;

        let mut meas = Vec::new();
        for (i, (sat_pos, n1, n2)) in sats.iter().enumerate() {
            let geom = (*sat_pos - true_pos).norm() - (ref_pos - true_pos).norm();
            let dd1 = geom / lambda1; // integer-free phases: truth-consistent
            let dd2 = geom / lambda2;
            meas.push((make_key(2 + i as u16), *sat_pos, combine_iono_free(F1, F2, dd1, dd2), *n1 + 1.0, *n2 + 1.0));
        }
        let if_meas: Vec<IonoFreeMeasurement> = meas.iter().map(|(k, sat_pos, phase_if, _, _)| {
            IonoFreeMeasurement {
                key: *k,
                sat_pos: *sat_pos,
                ref_pos,
                base_pos,
                lambda_if: lambda_iono_free(F1, F2),
                b2: 2,
                f1_hz: F1,
                f2_hz: F2,
                dd_phase_if_cycles: *phase_if,
                variance_cycles2: 1e-4,
                dgrad_n_rov: 0.0,
                dgrad_e_rov: 0.0,
                baseline_m: 1000.0,
            }
        }).collect();
        let fixed: Vec<(DoubleDiffKey, f64)> = meas.iter()
            .flat_map(|(k, _, _, n1, n2)| vec![(*k, *n1), (DoubleDiffKey { freq_band: 2, ..*k }, *n2)])
            .collect();
        let ar = ArResult {
            position_ecef: true_pos,
            cov_position: Matrix3::identity(),
            ratio: 3.0,
            is_fixed: true,
            num_ambiguities: 12,
            fixed_ambiguities: fixed,
        };
        assert!(matches!(
            apply_fixed_iono_free(&state, &if_meas, &ar),
            IonoFreeOutcome::Rejected
        ));
    }
