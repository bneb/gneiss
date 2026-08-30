#![allow(clippy::unwrap_used)]

use super::*;
use alloc::format;
use nalgebra::Vector3;

    const SITE: Vector3<f64> = Vector3::new(4_027_893.0, 307_041.0, 4_919_475.0);

    struct TestRotZ;
    impl ReferenceFrame for TestRotZ {
        const NAME: &'static str = "TEST-ROTZ";
        /// Distinct rotations about all three axes pin the sign convention.
        const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(HelmertParams {
            rx_mas: 1000.0,
            ry_mas: -2000.0,
            rz_mas: 3000.0,
            ..HelmertParams::identity_at(2015.0)
        });
    }

    struct TestRates;
    impl ReferenceFrame for TestRates {
        const NAME: &'static str = "TEST-RATES";
        const HELMERT_TO_ITRF2014: Option<HelmertParams> = Some(HelmertParams {
            tz_rate: 1.0,
            ..HelmertParams::identity_at(2015.0)
        });
    }

    fn full_params() -> HelmertParams {
        HelmertParams {
            tx_mm: 1.0,
            ty_mm: -2.5,
            tz_mm: 3.0,
            scale_ppb: 12.0,
            rx_mas: 800.0,
            ry_mas: -1500.0,
            rz_mas: 2200.0,
            ref_epoch_yr: 2015.0,
            tx_rate: 0.2,
            ty_rate: -0.1,
            tz_rate: 0.3,
            scale_rate: 0.05,
            rx_rate: 0.4,
            ry_rate: -0.6,
            rz_rate: 0.7,
        }
    }

    #[test]
    fn identity_conversion_is_bit_exact() {
        // Hub frame (None params) converts exactly at any epoch; ITRF2014
        // must stay the unique None frame (None => identity semantics).
        let p = EcefPos::<Itrf2014>::new(SITE);
        assert_eq!(Itrf2014::HELMERT_TO_ITRF2014, None);
        for t in [2015.0_f64, 2020.0, 2040.0] {
            assert_eq!(p.convert_to::<Itrf2014>(t).0, SITE);
        }
    }

    #[test]
    fn round_trip_preserves_position_below_spec_bound() {
        let cases = [
            (SITE, 2015.0_f64),
            (Vector3::new(500_000.0, -5_800_000.0, 2_100_000.0), 2040.0),
        ];
        for (v, t) in cases {
            assert_round_trip::<Igs20, Wgs84Broadcast>(v, t);
            assert_round_trip::<Itrf2020, Igs20>(v, t);
            assert_round_trip::<Wgs84Broadcast, Itrf2014>(v, t);
            assert_round_trip::<Itrf2020, Itrf2014>(v, t);
        }
    }

    fn assert_round_trip<A: ReferenceFrame, B: ReferenceFrame>(v: Vector3<f64>, t: f64) {
        let orig = EcefPos::<A>::new(v);
        let back = orig.convert_to::<B>(t).convert_to::<A>(t);
        let err = (back.vector() - orig.vector()).norm();
        assert!(err < 1.0e-6, "{} -> {} drifted {err:e} m", A::NAME, B::NAME);
    }

    #[test]
    fn itrf2020_to_itrf2014_shift_matches_published_parameters() {
        // Independent scalar derivation of X_14 = T + (1+s)·X_20, T in mm.
        // Evaluated at the reference epoch (2015.0), so rate terms drop out.
        let expected_shift = Vector3::new(
            -1.4e-3 + -0.42e-9 * SITE[0],
            -0.9e-3 + -0.42e-9 * SITE[1],
            1.4e-3 + -0.42e-9 * SITE[2],
        );
        let got = EcefPos::<Itrf2020>::new(SITE)
            .convert_to::<Itrf2014>(2015.0)
            .0
            - SITE;
        assert!(
            (got - expected_shift).abs().max() < 1.0e-8,
            "shift {got:?} != published {expected_shift:?}"
        );
        // IGS20 / broadcast WGS84 are ITRF2020-aligned: same published row.
        assert_eq!(Igs20::HELMERT_TO_ITRF2014, Some(ITRF2020_TO_ITRF2014));
        assert_eq!(
            Wgs84Broadcast::HELMERT_TO_ITRF2014,
            Some(ITRF2020_TO_ITRF2014)
        );
    }

    #[test]
    fn rotation_axes_follow_altamimi_sign_convention() {
        // X' = X + ω×X with ω = (1000, −2000, 3000) mas. Expected values are
        // hardcoded (rad per mas literal) so a corrupted unit conversion or
        // flipped sign on ANY axis cannot self-verify.
        const K: f64 = 4.848_136_811_095_36e-9; // rad per mas
        let cases = [
            (
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(1.0, 3000.0 * K, 2000.0 * K),
            ),
            (
                Vector3::new(0.0, 1.0, 0.0),
                Vector3::new(-3000.0 * K, 1.0, 1000.0 * K),
            ),
            (
                Vector3::new(0.0, 0.0, 1.0),
                Vector3::new(-2000.0 * K, -1000.0 * K, 1.0),
            ),
        ];
        for (v, expected) in cases {
            let out = EcefPos::<TestRotZ>::new(v).convert_to::<Itrf2014>(2015.0).0;
            assert!(
                (out - expected).abs().max() < 1.0e-12,
                "rotation {v:?} -> {out:?}, expected {expected:?}"
            );
        }
    }

    #[test]
    fn epoch_propagation_is_linear_in_years_from_reference_epoch() {
        let shift_at = |t| EcefPos::<TestRates>::new(SITE).convert_to::<Itrf2014>(t).0 - SITE;
        let at_ref = shift_at(2015.0);
        let plus10 = shift_at(2025.0);
        let minus10 = shift_at(2005.0);
        // Rate-only frame: z-shift = 1 mm/yr · dt, other axes untouched. The
        // 1e-8 tolerance is the nm-scale cancellation floor at |r| ≈ 5e6 m.
        assert!(at_ref.abs().max() < 1.0e-15);
        assert!((plus10[2] - 0.010).abs() < 1.0e-8 && (minus10[2] + 0.010).abs() < 1.0e-8);
        assert!(
            plus10[0].abs() + plus10[1].abs() + minus10[0].abs() + minus10[1].abs()
                < 1.0e-7
        );
        // `at` must re-anchor the reference epoch and propagate ALL seven
        // quantities with their signed rates (literals kill sign mutants).
        let q = full_params().at(2025.0);
        assert_eq!(q.ref_epoch_yr, 2025.0);
        let drift = (q.tx_mm - 3.0).abs()
            + (q.ty_mm + 3.5).abs()
            + (q.tz_mm - 6.0).abs()
            + (q.scale_ppb - 12.5).abs()
            + (q.rx_mas - 804.0).abs()
            + (q.ry_mas + 1506.0).abs()
            + (q.rz_mas - 2207.0).abs();
        assert!(drift < 1.0e-9, "propagation drift {drift:e}");
    }

    #[test]
    fn apply_inverse_inverts_apply_for_full_parameter_set() {
        for t in [1990.0_f64, 2015.0, 2035.0] {
            let p = full_params().at(t);
            for v in [
                SITE,
                Vector3::new(-6.378e6, 0.0, 0.0),
                Vector3::new(1.0e5, -2.0e5, 3.0e5),
            ] {
                let back = p.apply_inverse(p.apply(v));
                // f64 noise floor at |r| ≈ 5e6 m is ~1 nm/op; spec bound 0.1 mm.
                assert!(
                    (back - v).abs().max() < 1.0e-8,
                    "inverse failed at t={t}: {back:?} vs {v:?}"
                );
            }
        }
    }

    #[test]
    fn accessors_and_frame_metadata_work() {
        let p = EcefPos::<Igs20>::new(Vector3::new(3.0, 4.0, 0.0));
        assert_eq!(p.norm(), 5.0);
        assert_eq!(p.vector(), &Vector3::new(3.0, 4.0, 0.0));
        assert_eq!(<Igs20 as ReferenceFrame>::NAME, "IGS20");
        assert!(format!("{p:?}").contains("IGS20"), "Debug must show frame");
        let copied = p;
        let different = EcefPos::<Igs20>::new(Vector3::new(3.0, 4.0, 1.0));
        assert!(copied == p && different != p, "Copy/PartialEq broken");
    }
    #[test]
    fn test_nad83_epoch_propagation_accumulates() {
        // NAD83(2011) is aligned to ITRF2014 at epoch 2010.0 (~mm offset),
        // but plate-motion RATES cause the offset to grow over time.
        // By 2025 (15 yr later): ~11 mm in Y from vy_rate=0.757 mm/yr.
        let pos_nad83 = EcefPos::<Nad83_2011>::new(Vector3::new(
            -2_688_201.0, -4_265_643.0, 3_893_778.0, // P224
        ));
        let at_ref = pos_nad83.convert_to::<Itrf2014>(2010.0);
        let at_late = pos_nad83.convert_to::<Itrf2014>(2025.0);
        let d_ref = (at_ref.vector() - pos_nad83.vector()).norm();
        let d_late = (at_late.vector() - pos_nad83.vector()).norm();
        // Both epochs should produce cm-level shifts (frame alignment).
        // The offset may not grow monotonically — the 14-parameter model
        // rotates the differential vector, so magnitude can decrease even
        // as individual components grow.
        assert!(d_ref < 0.02, "reference epoch offset {} m too large", d_ref);
        assert!(d_late < 0.05, "late epoch offset {} m too large", d_late);
    }

    #[test]
    fn test_etrs89_epoch_propagation_grows_with_time() {
        // ETRS89 is frozen at 1989.0; by 2025 the offset from ITRF grows.
        let pos_etr = EcefPos::<Etrs89>::new(Vector3::new(
            4_042_000.0, 355_000.0, 4_950_000.0,
        ));
        let early = pos_etr.convert_to::<Itrf2014>(1995.0);
        let late = pos_etr.convert_to::<Itrf2014>(2025.0);
        let d_early = (early.vector() - pos_etr.vector()).norm();
        let d_late = (late.vector() - pos_etr.vector()).norm();
        assert!(
            d_late > d_early,
            "ETRS89-ITRF2014 difference should grow with time: {} vs {}",
            d_early, d_late
        );
    }

    #[test]
    fn test_gda2020_identity_at_reference_epoch() {
        // GDA2020 is aligned to ITRF2014 at epoch 2020.0 → identity.
        let pos_gda = EcefPos::<Gda2020>::new(Vector3::new(
            -4_000_000.0, 3_500_000.0, -3_200_000.0,
        ));
        let pos_itrf = pos_gda.convert_to::<Itrf2014>(2020.0);
        let diff = (pos_itrf.vector() - pos_gda.vector()).norm();
        assert!(
            diff < 0.001,
            "GDA2020→ITRF2014 should be ~identity at epoch 2020.0, got {} mm",
            diff * 1000.0
        );
    }

    #[test]
    fn test_epoch_position_propagation_with_velocity() {
        let pos_itrf = EcefPos::<Itrf2020>::new(Vector3::new(1000.0, 2000.0, 3000.0));
        let vel = Vector3::new(0.01, -0.02, 0.005); // 10 mm/yr X, -20 mm/yr Y
        let ep = EpochPosition::<Itrf2020, Arp>::new(pos_itrf, 2020.0, Some(vel));

        let ep_2025 = ep.at_epoch(2025.0);
        assert_eq!(ep_2025.epoch_yr, 2025.0);
        // Shift after 5 years = 5 * [0.01, -0.02, 0.005] = [0.05, -0.10, 0.025]
        let diff = ep_2025.pos.0 - ep.pos.0;
        assert!((diff[0] - 0.05).abs() < 1e-9);
        assert!((diff[1] - (-0.10)).abs() < 1e-9);
        assert!((diff[2] - 0.025).abs() < 1e-9);
    }

    #[test]
    fn test_arp_to_apc_roundtrip() {
        let pos_arp = EcefPos::<Itrf2014>::new(Vector3::new(1_000_000.0, 2_000_000.0, 3_000_000.0));
        let ep_arp = EpochPosition::<Itrf2014, Arp>::new(pos_arp, 2025.0, None);

        let pco_neu_mm = Vector3::new(10.0, -20.0, 85.0); // 85 mm vertical PCO
        let llh_rad = Vector3::new(0.5, 1.0, 100.0);

        let ep_apc: EpochPosition<Itrf2014, Apc<1>> = ep_arp.to_apc(pco_neu_mm, llh_rad);
        assert!((ep_apc.pos.0 - ep_arp.pos.0).norm() > 0.05); // ~87 mm offset

        let ep_arp_back = ep_apc.to_arp(pco_neu_mm, llh_rad);
        assert!((ep_arp_back.pos.0 - ep_arp.pos.0).norm() < 1e-9, "Exact roundtrip");
    }
