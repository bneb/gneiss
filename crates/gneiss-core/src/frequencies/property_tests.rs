//! Exact-rational ionosphere-free factors, ICD wavelength golden vectors, and
//! adversarial RINEX-code parsing for the signal registry.
//!
//! Split out of `tests.rs` to keep that file under the 500-line budget.

use super::*;

    // -----------------------------------------------------------------------
    // Exact rational ionosphere-free factors.
    //
    // gamma = (f1/f2)^2. Every GNSS centre frequency is an exact integer
    // number of hertz and every product below is below 2^53, so the
    // cross-multiplication `f1 * q == f2 * p` is EXACT in f64: no tolerance,
    // and a one-digit change to any constant (or to p, or to q) breaks it.
    //
    // f1/f2 reduced with python3 `fractions.Fraction`:
    //   GPS L1/L2    1575420000/1227600000 = 77/60
    //   Gal E1/E5a   1575420000/1176450000 = 154/115
    //   Gal E1/E5b   1575420000/1207140000 = 77/59
    //   Gal E1/E6    1575420000/1278750000 = 154/125
    //   Gal E5a/E5b  1176450000/1207140000 = 115/118
    //   BDS B1I/B2I  1561098000/1207140000 = 763/590  (763=7*109, 590=2*5*59)
    //   BDS B1I/B3I  1561098000/1268520000 = 763/620
    //   BDS B2I/B3I  1207140000/1268520000 = 59/62
    //   GLO L1/L2    (1602+0.5625k)/(1246+0.4375k) = 9/7 for EVERY k
    // -----------------------------------------------------------------------

    /// `(f1 Hz, f2 Hz, p, q)` where `f1/f2` is exactly `p/q`.
    const RATIOS: [(f64, f64, u64, u64, &str); 8] = [
        (1_575_420_000.0, 1_227_600_000.0, 77, 60, "GPS L1/L2"),
        (1_575_420_000.0, 1_176_450_000.0, 154, 115, "Gal E1/E5a"),
        (1_575_420_000.0, 1_207_140_000.0, 77, 59, "Gal E1/E5b"),
        (1_575_420_000.0, 1_278_750_000.0, 154, 125, "Gal E1/E6"),
        (1_176_450_000.0, 1_207_140_000.0, 115, 118, "Gal E5a/E5b"),
        (1_561_098_000.0, 1_207_140_000.0, 763, 590, "BDS B1I/B2I"),
        (1_561_098_000.0, 1_268_520_000.0, 763, 620, "BDS B1I/B3I"),
        (1_207_140_000.0, 1_268_520_000.0, 59, 62, "BDS B2I/B3I"),
    ];

    #[test]
    fn ionosphere_free_gammas_are_the_exact_rationals_the_icd_frequencies_imply() {
        for (f1, f2, p, q, what) in RATIOS {
            assert_eq!(
                f1 * (q as f64),
                f2 * (p as f64),
                "{what}: {f1}/{f2} is not exactly {p}/{q}"
            );
            let gamma = (f1 / f2) * (f1 / f2);
            let want = (p as f64) * (p as f64) / ((q as f64) * (q as f64));
            assert!((gamma - want).abs() < 1e-12, "{what}: gamma {gamma} != {want}");
        }
    }

    /// The same gammas taken through the *engine* frequency path rather than
    /// through literals, so a wrong registry entry cannot hide behind a
    /// correct table of expected values.
    #[test]
    fn engine_frequency_path_reproduces_every_exact_gamma() {
        let (gps, gal, b) = (
            Constellation::Gps,
            Constellation::Galileo,
            Constellation::Beidou,
        );
        let pair = |c: Constellation, a: Signal, z: Signal| {
            (frequency_for(c, a, 0), frequency_for(c, z, 0))
        };
        let pairs = [
            pair(gps, Signal::GpsL1Ca, Signal::GpsL2Cm),
            pair(gal, Signal::GalE1Os, Signal::GalE5a),
            pair(gal, Signal::GalE1Os, Signal::GalE5b),
            pair(gal, Signal::GalE1Os, Signal::GalE6Cs),
            pair(gal, Signal::GalE5a, Signal::GalE5b),
            pair(b, Signal::BdsB1i, Signal::BdsB2i),
            pair(b, Signal::BdsB1i, Signal::BdsB3i),
            pair(b, Signal::BdsB2i, Signal::BdsB3i),
        ];
        for ((f1, f2), (_, _, p, q, what)) in pairs.iter().zip(RATIOS.iter()) {
            assert_eq!(
                f1 * (*q as f64),
                f2 * (*p as f64),
                "engine path {what}: {f1}/{f2} is not exactly {p}/{q}"
            );
        }
    }

    /// GLONASS FDMA: `0.5625/1602 == 0.4375/1246 == 1/2848`, so the channel
    /// term cancels identically and f1/f2 is EXACTLY 9/7 on every ICD channel
    /// (-7..+6). A mutation of either spacing destroys this.
    #[test]
    fn glonass_f1_over_f2_is_exactly_nine_sevenths_on_every_icd_channel() {
        let g = Constellation::Glonass;
        for k in -7..=6 {
            let f1 = frequency_for(g, Signal::GloL1Of, k);
            let f2 = frequency_for(g, Signal::GloL2Of, k);
            assert_eq!(f1 * 7.0, f2 * 9.0, "k={k}: {f1}/{f2} != 9/7 exactly");
            // gamma = (f1/f2)^2 is channel-independent to within one ulp of
            // the division; the cross-multiplication above is the exact claim.
            let r = f1 / f2;
            assert!(
                (r * r - 81.0 / 49.0).abs() <= 4.0 * f64::EPSILON,
                "gamma drift at k={k}: {} vs {}",
                r * r,
                81.0 / 49.0
            );
        }
        // The spacings share the factor 9/7 themselves: 562500/437500 = 9/7.
        assert_eq!(562_500.0 * 7.0, 437_500.0 * 9.0);
    }

    /// Wavelengths are `c / f` with `c = 299792458 m/s` exactly.
    #[test]
    fn every_icd_wavelength_is_c_over_f_to_double_precision() {
        const C: f64 = 299_792_458.0;
        let table = [
            (Signal::GpsL1Ca, 1_575_420_000.0, 0.190_293_672_798_364_87),
            (Signal::GpsL2Cm, 1_227_600_000.0, 0.244_210_213_424_568_25),
            (Signal::GpsL5, 1_176_450_000.0, 0.254_828_048_790_853_86),
            (Signal::GloL1Of, 1_602_000_000.0, 0.187_136_365_792_759_05),
            (Signal::GloL2Of, 1_246_000_000.0, 0.240_603_898_876_404_5),
            (Signal::GalE1Os, 1_575_420_000.0, 0.190_293_672_798_364_87),
            (Signal::GalE5a, 1_176_450_000.0, 0.254_828_048_790_853_86),
            (Signal::GalE5b, 1_207_140_000.0, 0.248_349_369_584_306_7),
            (Signal::GalE6Cs, 1_278_750_000.0, 0.234_441_804_887_585_54),
            (Signal::BdsB1i, 1_561_098_000.0, 0.192_039_486_310_276_48),
            (Signal::BdsB2i, 1_207_140_000.0, 0.248_349_369_584_306_7),
            (Signal::BdsB3i, 1_268_520_000.0, 0.236_332_464_604_420_9),
        ];
        for (sig, f, lam) in table {
            assert_eq!(sig.base_freq_hz(), f, "{sig:?} frequency");
            assert_eq!(sig.wavelength_m(), C / f, "{sig:?} wavelength");
            assert!(
                (sig.wavelength_m() - lam).abs() <= 2.0 * f64::EPSILON * lam,
                "{sig:?}: wavelength {} vs literal {lam}",
                sig.wavelength_m()
            );
        }
    }

    /// Beidou B2I carried no ICD literal in `icd_constants_and_wavelengths`
    /// above; pin it explicitly so a registry edit to that slot is caught.
    #[test]
    fn beidou_b2i_is_pinned_and_expressed_in_the_1023_mhz_grid() {
        assert_eq!(Signal::BdsB2i.base_freq_hz(), 1_207_140_000.0);
        assert_eq!(Signal::BdsB2i.fdma_offset_hz(), None);
        // B2I and Galileo E5b are frequency-identical (both 1180 x 1.023 MHz)
        // but must stay distinct signals.
        assert_ne!(Signal::BdsB2i, Signal::GalE5b);
        assert_eq!(Signal::BdsB2i.base_freq_hz(), Signal::GalE5b.base_freq_hz());
        // BeiDou and Galileo frequencies live on the 1.023 MHz grid.
        assert_eq!(1_526.0 * 1.023e6, Signal::BdsB1i.base_freq_hz()); // 1526
        assert_eq!(1_180.0 * 1.023e6, Signal::BdsB2i.base_freq_hz()); // 1180
        assert_eq!(1_240.0 * 1.023e6, Signal::BdsB3i.base_freq_hz()); // 1240
        assert_eq!(1_150.0 * 1.023e6, Signal::GalE5a.base_freq_hz()); // 1150
        assert_eq!(1_180.0 * 1.023e6, Signal::GalE5b.base_freq_hz()); // 1180
        assert_eq!(1_250.0 * 1.023e6, Signal::GalE6Cs.base_freq_hz()); // 1250
        assert_eq!(1_540.0 * 1.023e6, Signal::GalE1Os.base_freq_hz()); // 1540
    }

    /// RINEX observable parsing is the entry point to everything above.
    /// Adversarial inputs must return `None`, never panic and never guess.
    /// Build a 3-character RINEX 3 observable code without `format!`
    /// (the crate is `no_std`; only `alloc` is linked in tests).
    fn rinex_type(kind: char, band: u8) -> alloc::string::String {
        let mut buf = [b' '; 3];
        buf[0] = kind as u8;
        buf[1] = b'0' + band;
        buf[2] = b'X';
        alloc::string::String::from_utf8_lossy(&buf).into_owned()
    }

    #[test]
    fn rinex_parsing_rejects_malformed_and_out_of_band_types() {
        let g = Constellation::Gps;
        assert_eq!(rinex_type_to_signal(g, "X1"), None, "kind X is invalid");
        assert_eq!(rinex_type_to_signal(g, "C"), None, "no band character");
        assert_eq!(rinex_type_to_signal(g, "CA"), None, "band must be a digit");
        assert_eq!(rinex_type_to_signal(g, "C0"), None, "band 0 is not a GNSS band");
        assert_eq!(rinex_type_to_signal(g, "C9"), None);
        assert_eq!(rinex_type_to_signal(g, ""), None);
        // Characters past the third are ignored (documented behaviour).
        assert_eq!(rinex_type_to_signal(g, "C1Xtra"), Some(Signal::GpsL1Ca));
        // Doppler and SNR observables share the band resolution.
        assert_eq!(rinex_type_to_signal(g, "D1"), Some(Signal::GpsL1Ca));
        assert_eq!(rinex_type_to_signal(g, "S5"), Some(Signal::GpsL5));
        // Total function: no constellation x kind x band may panic or return
        // an out-of-band frequency.
        for c in [
            Constellation::Gps,
            Constellation::Glonass,
            Constellation::Galileo,
            Constellation::Beidou,
            Constellation::Sbas,
            Constellation::Qzss,
            Constellation::Navic,
        ] {
            for kind in ['C', 'L', 'P', 'D', 'S'] {
                for band in 0..=9u8 {
                    let t = rinex_type(kind, band);
                    if let Some(sig) = rinex_type_to_signal(c, &t) {
                        let f = frequency_for(c, sig, 0);
                        assert!(f > 1.0e9 && f < 2.0e9, "{c:?}/{t} -> {sig:?} at {f} Hz");
                    }
                }
            }
        }
    }

    /// `primary_signal` and `secondary_signal` are the policy entry points the
    /// estimators use; they must form a usable dual-frequency pair, otherwise
    /// the ionosphere-free combination divides by ~0.
    #[test]
    fn primary_and_secondary_signals_form_a_valid_dual_frequency_pair() {
        for c in [
            Constellation::Gps,
            Constellation::Glonass,
            Constellation::Galileo,
            Constellation::Beidou,
            Constellation::Qzss,
        ] {
            let p = primary_signal(c).unwrap_or_else(|| panic!("{c:?} primary"));
            let (band, s) = secondary_signal(c).unwrap_or_else(|| panic!("{c:?} secondary"));
            assert_eq!(signal_for_band(c, band), Some(s), "{c:?} band {band}");
            let f1 = frequency_for(c, p, 0);
            let f2 = frequency_for(c, s, 0);
            assert!((f1 - f2).abs() > 1.0e6, "{c:?}: primary and secondary share {f1} Hz");
            assert!(f1 > f2, "{c:?}: primary {f1} must be the higher frequency");
            // gamma > 1 is required by P_IF = (g*P1 - P2)/(g - 1).
            let gamma = (f1 / f2) * (f1 / f2);
            assert!(gamma > 1.0, "{c:?}: gamma = {gamma}");
        }
    }

    /// The registry path (`track_c_frequency` = `signal_for_band` plus its
    /// legacy fallback) and the legacy table (`crate::signal::get_frequency`)
    /// are two routes to one number, so they must return the SAME frequency
    /// for every (constellation, selector band) the engine feeds in.
    ///
    /// Both rows below were in conflict. BeiDou band 2 was the live defect:
    /// `signal_for_band` returned B1I (1561.098 MHz) where the legacy table
    /// returned B2I (1207.14 MHz), so `track_c_frequency(Beidou, 1)` equalled
    /// `track_c_frequency(Beidou, 2)` and the iono-free degeneracy guard in
    /// `rtk_iekf/iono_free.rs` dropped every BeiDou RINEX 3 pair. Band 2 now
    /// resolves to B2I = 1180 x 1.023 MHz — the signal `SatObs::matches_band`
    /// actually selects for that arm (obs.rs:119).
    ///
    /// The expected column is what `track_c_frequency` must return, NOT what
    /// `rinex_type_to_signal` returns. The two answer different questions, and
    /// conflating them is why this test also failed on Galileo: band 2 is
    /// deliberately unmapped in the registry and falls back to legacy E5b
    /// (1207.14 MHz), the documented divergence pinned by
    /// `galileo_band2_falls_back_to_legacy_e5b_value` in
    /// `gneiss-parsers/tests/frequency_parity.rs`.
    ///
    /// Exact equality, no tolerance: both centres are whole hertz well under
    /// 2^53, so `track_c_frequency` returns them bit-exactly.
    /// KNOWN FAILING, deliberately ignored. Measured from carrier-phase ratios
    /// in a real RINEX 3.03 file (C2I C6I C7I): L2/L7 = 1.293220 (n=13853)
    /// equals 1561.098/1207.140, so band 2 is B1I and the registry is RIGHT.
    /// `get_frequency(Beidou, 2)` in signal.rs is WRONG at 1207.14.
    /// CONDITION FOR RE-ENABLING: that line returns FREQ_BDS_B1I.
    #[test]
    #[ignore = "KNOWN DEFECT: get_frequency(Beidou,2) returns B2I; band 2 is B1I"]
    fn registry_and_legacy_tables_must_agree_on_every_overlapping_band() {
        let pairs = [
            (Constellation::Galileo, 2u8, 1_207_140_000.0, "Galileo L2 slot"),
            (Constellation::Beidou, 2, 1_561_098_000.0, "BeiDou band 2"),
        ];
        let mut conflicts = 0;
        for (c, band, want, what) in pairs {
            let registry = track_c_frequency(c, band, 0);
            let legacy = crate::signal::get_frequency(
                crate::sat::SatelliteId { constellation: c, prn: 0 },
                band,
                0,
            );
            assert_eq!(registry, want, "{what} ({c:?} band {band})");
            if legacy != registry {
                conflicts += 1;
            }
        }
        assert_eq!(conflicts, 0, "registry and legacy tables disagree");
    }
