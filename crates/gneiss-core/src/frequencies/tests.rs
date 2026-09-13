#![allow(clippy::unwrap_used)]

    use super::*;
    use crate::signal::{
        FREQ_BDS_B1I, FREQ_GAL_E5B, FREQ_GLO_L1_DELTA, FREQ_GLO_L1_NOMINAL, FREQ_GLO_L2_DELTA,
        FREQ_GLO_L2_NOMINAL, FREQ_GPS_L1, FREQ_GPS_L2, FREQ_GPS_L5,
    };

    const HZ_TOL: f64 = 1.0; // spec: match ICD to <1 Hz

    fn assert_close(actual: f64, expected: f64, tol: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= tol,
            "{what}: got {actual}, expected {expected} (tol {tol})"
        );
    }

    /// Asserts every whitespace-separated RINEX type maps to `want`.
    fn assert_maps(cons: Constellation, types: &str, want: Option<Signal>) {
        for t in types.split_whitespace() {
            assert_eq!(rinex_type_to_signal(cons, t), want, "type {t:?}");
        }
    }

    #[test]
    fn icd_constants_and_wavelengths() {
        // (signal, ICD centre frequency Hz, lambda m); refs computed independently.
        let cases = [
            (Signal::GpsL1Ca, 1575.42e6, 0.190293672798),
            (Signal::GpsL1P, 1575.42e6, 0.190293672798),
            (Signal::GpsL2Cm, 1227.60e6, 0.244210213425),
            (Signal::GpsL2P, 1227.60e6, 0.244210213425),
            (Signal::GpsL5, 1176.45e6, 0.254828048791),
            (Signal::GloL1Of, 1602.0e6, 0.187136365793), // nominal k=0
            (Signal::GloL2Of, 1246.0e6, 0.240603898876),
            (Signal::GalE1Os, 1575.42e6, 0.190293672798),
            (Signal::GalE5a, 1176.45e6, 0.254828048791),
            (Signal::GalE5b, 1207.14e6, 0.248349369584),
            (Signal::GalE6Cs, 1278.75e6, 0.234441804888),
            (Signal::BdsB1i, 1561.098e6, 0.192039486310),
            (Signal::BdsB3i, 1268.52e6, 0.236332464604),
        ];
        for (signal, icd_hz, icd_lam) in cases {
            assert_close(signal.base_freq_hz(), icd_hz, HZ_TOL, "base freq");
            let lam = signal.wavelength_m();
            assert_close(lam, icd_lam, 1e-9, "lambda ref");
            // Cross-path check through wavelength_for: FDMA signals at k=0,
            // CDMA signals with an arbitrary channel, must both agree.
            let via_for = if signal.fdma_offset_hz().is_some() {
                wavelength_for(Constellation::Glonass, signal, 0)
            } else {
                wavelength_for(Constellation::Gps, signal, 3)
            };
            assert_close(via_for, lam, 1e-12, "wavelength_for cross-path");
        }
    }

    #[test]
    fn fdma_offsets_exhaustive() {
        assert_eq!(Signal::GloL1Of.fdma_offset_hz(), Some(562_500.0));
        assert_eq!(Signal::GloL2Of.fdma_offset_hz(), Some(437_500.0));
        // One CDMA representative per constellation: all share the same
        // `_ => None` arm, so a per-variant sweep adds no mutant coverage.
        let cdma = [Signal::GpsL1Ca, Signal::GalE5b, Signal::BdsB3i];
        for signal in cdma {
            assert_eq!(signal.fdma_offset_hz(), None, "{signal:?}");
        }
    }

    #[test]
    fn glonass_fdma_channel_sweep() {
        let g = Constellation::Glonass;
        // Independent formula: f(k) = f0 + k * delta, valid channels −7..+6.
        for k in -7..=6i8 {
            let f1 = frequency_for(g, Signal::GloL1Of, k);
            let f2 = frequency_for(g, Signal::GloL2Of, k);
            assert_close(f1, 1602.0e6 + f64::from(k) * 562_500.0, HZ_TOL, "L1(k)");
            assert_close(f2, 1246.0e6 + f64::from(k) * 437_500.0, HZ_TOL, "L2(k)");
        }
        // Published GLONASS ICD edge-channel spot checks.
        let (l1, l2) = (Signal::GloL1Of, Signal::GloL2Of);
        assert_close(frequency_for(g, l1, -7), 1598.0625e6, HZ_TOL, "L1(k=-7)");
        assert_close(frequency_for(g, l1, 6), 1605.375e6, HZ_TOL, "L1(k=+6)");
        assert_close(frequency_for(g, l2, -7), 1242.9375e6, HZ_TOL, "L2(k=-7)");
        assert_close(frequency_for(g, l2, 6), 1248.625e6, HZ_TOL, "L2(k=+6)");
        // Channel-correct GLONASS wavelengths at the extreme channels.
        let w1 = wavelength_for(g, l1, -7);
        let w2 = wavelength_for(g, l2, 6);
        assert_close(w1, 0.187597455043, 1e-9, "L1 lambda k=-7");
        assert_close(w2, 0.240098074282, 1e-9, "L2 lambda k=+6");
    }

    #[test]
    fn gps_rinex_types_map_correctly() {
        let g = Constellation::Gps;
        assert_maps(g, "C1 L1 D1 S1 C1C L1C c1", Some(Signal::GpsL1Ca));
        assert_maps(g, "P1 p1 C1W L1W C1Y", Some(Signal::GpsL1P));
        assert_maps(g, "C2 L2 C2X c2x l2", Some(Signal::GpsL2Cm));
        assert_maps(g, "P2 C2W L2W", Some(Signal::GpsL2P));
        assert_maps(g, "C5 L5 C5Q l5q", Some(Signal::GpsL5));
        // Unmodelled bands plus malformed / invalid-kind strings.
        assert_maps(g, "C6 C7 C8 X1 1C Q1", None);
        assert_eq!(rinex_type_to_signal(g, ""), None);
        assert_eq!(rinex_type_to_signal(g, "C"), None);
    }

    #[test]
    fn glonass_rinex_types_map_correctly() {
        let r = Constellation::Glonass;
        assert_maps(r, "C1 L1 P1 D1 d1", Some(Signal::GloL1Of));
        assert_maps(r, "C2 L2 P2 C2P c2p", Some(Signal::GloL2Of)); // P shares OF carrier
        assert_maps(r, "C3 C5 C6", None);
        // Channel number rides along: R01 with k=-4 (classic GLONASS slot).
        let f1 = frequency_for(r, Signal::GloL1Of, -4);
        assert_close(f1, 1599.75e6, HZ_TOL, "R01 L1 k=-4");
    }

    /// THE bug-source test: Galileo "L2" in RINEX 2.11 mixed files.
    ///
    /// Two conventions exist in the wild: converters writing E5b into the L2
    /// slot (matching the old `get_frequency(band=2)` behaviour) and converters
    /// writing E5a there (the project dataset directive). This registry pins
    /// the policy to **E5a** — see module documentation — while E5b stays
    /// reachable through its unambiguous RINEX 3 band 7 (`"L7"`).
    #[test]
    fn galileo_l2_maps_to_e5a_with_documented_ambiguity() {
        assert_maps(Constellation::Galileo, "L2 C2 P2 C2W", Some(Signal::GalE5a));
        // E5b must remain reachable, and must NOT be what "L2" gives.
        assert_eq!(
            rinex_type_to_signal(Constellation::Galileo, "L7"),
            Some(Signal::GalE5b)
        );
        let l2_hz = frequency_for(Constellation::Galileo, Signal::GalE5a, 0);
        assert_close(l2_hz, 1176.45e6, HZ_TOL, "Galileo L2-slot = E5a");
        assert!((l2_hz - FREQ_GAL_E5B).abs() > 30.0e6, "E5a/E5b must differ");
    }

    #[test]
    fn galileo_unambiguous_bands() {
        let e = Constellation::Galileo;
        assert_maps(e, "C1 L1 c1a", Some(Signal::GalE1Os));
        assert_maps(e, "C5 L5 C5Q", Some(Signal::GalE5a));
        assert_maps(e, "C6 L6 C6A", Some(Signal::GalE6Cs));
        assert_maps(e, "C7 L7 C7X", Some(Signal::GalE5b));
        assert_maps(e, "C8 L8", None); // E5 AltBOC composite unmodelled
    }

    #[test]
    fn beidou_rinex_types_map_correctly() {
        let c = Constellation::Beidou;
        // RINEX 3: B1I lives on band 2 ("C2I" per real MGEX headers); B3I on 6.
        assert_maps(c, "C2 C2I L2I D2I c2i", Some(Signal::BdsB1i));
        let b1i = frequency_for(c, Signal::BdsB1i, 0);
        assert_close(b1i, 1561.098e6, HZ_TOL, "B1I");
        assert_maps(c, "C6 L6 C6I", Some(Signal::BdsB3i));
        // RINEX 2 legacy: plain band-1 slot carried B1I.
        assert_maps(c, "C1 L1", Some(Signal::BdsB1i));
        // RINEX 3 band 1 is BDS-3 B1C @1575.42 (attrs P/X/D) — NO variant here;
        // must resolve None, never B1I @1561.098 (14.3 MHz error otherwise).
        assert_maps(c, "C1P C1X L1D", None);
        // Band 7 carries BeiDou B2I @1207.14 MHz.
        assert_maps(c, "C7 C7D L7", Some(Signal::BdsB2i));
        // B2a (band 5), B1C (band 1) unmodelled.
        assert_maps(c, "C5 L5 C3 C8", None);
    }

    #[test]
    fn other_constellations_follow_documented_policy() {
        let (q, s) = (Constellation::Qzss, Constellation::Sbas);
        // QZSS aliases the GPS radio plan (identical frequencies); no L6.
        assert_maps(q, "C1 L1", Some(Signal::GpsL1Ca));
        assert_maps(q, "C2 L2", Some(Signal::GpsL2Cm));
        assert_maps(q, "C5 L5", Some(Signal::GpsL5));
        assert_eq!(rinex_type_to_signal(q, "C6"), None);
        // SBAS: L1 + L5 only.
        assert_maps(s, "C1 L1", Some(Signal::GpsL1Ca));
        assert_maps(s, "C5 L5", Some(Signal::GpsL5));
        assert_maps(s, "C2 L2 C6", None);
        // NavIC (L5/S-band) explicitly unmodelled.
        assert_maps(Constellation::Navic, "C1 L1 C5 L5 S5 C9", None);
        // Case-insensitive parsing.
        assert_maps(Constellation::Galileo, "c2", Some(Signal::GalE5a));
        assert_maps(Constellation::Gps, "c1c", Some(Signal::GpsL1Ca));
    }

    #[test]
    fn frequency_for_ignores_k_for_cdma_and_guards_mismatch() {
        let (gps, glo) = (Constellation::Gps, Constellation::Glonass);
        // CDMA signals ignore the channel entirely — including a GLONASS sat
        // argument, whose k must not leak into a CDMA signal's frequency.
        let f = frequency_for(gps, Signal::GpsL1Ca, 42);
        assert_close(f, FREQ_GPS_L1, 0.0, "k ignored");
        let f = frequency_for(glo, Signal::GpsL1Ca, -4);
        assert_close(f, FREQ_GPS_L1, 0.0, "CDMA on GLO");
        // Mismatched (constellation, FDMA signal) falls back to nominal k=0.
        let f = frequency_for(gps, Signal::GloL1Of, 3);
        assert_close(f, FREQ_GLO_L1_NOMINAL, 0.0, "guard");
        let f = frequency_for(gps, Signal::GloL2Of, -7);
        assert_close(f, FREQ_GLO_L2_NOMINAL, 0.0, "guard");
        // Out-of-ICD-range channels pass through arithmetically (no clamping);
        // callers validate the operational −7..+6 window.
        let f = frequency_for(glo, Signal::GloL1Of, 12);
        let want = 1602.0e6 + f64::from(12) * 562_500.0;
        assert_close(f, want, HZ_TOL, "k passthrough");
    }

    #[test]
    fn registry_agrees_with_legacy_constants() {
        // Guard against the two sources of truth drifting apart again.
        assert_close(Signal::GpsL1Ca.base_freq_hz(), FREQ_GPS_L1, 0.0, "L1");
        assert_close(Signal::GpsL2Cm.base_freq_hz(), FREQ_GPS_L2, 0.0, "L2");
        assert_close(Signal::GpsL5.base_freq_hz(), FREQ_GPS_L5, 0.0, "L5");
        assert_close(Signal::GalE1Os.base_freq_hz(), FREQ_GPS_L1, 0.0, "E1=L1");
        assert_close(Signal::GalE5b.base_freq_hz(), FREQ_GAL_E5B, 0.0, "E5b");
        assert_close(Signal::BdsB1i.base_freq_hz(), FREQ_BDS_B1I, 0.0, "B1I");
        let n1 = Signal::GloL1Of.base_freq_hz();
        let n2 = Signal::GloL2Of.base_freq_hz();
        assert_close(n1, FREQ_GLO_L1_NOMINAL, 0.0, "nom1");
        assert_close(n2, FREQ_GLO_L2_NOMINAL, 0.0, "nom2");
        let d1 = Signal::GloL1Of.fdma_offset_hz().unwrap_or_default();
        let d2 = Signal::GloL2Of.fdma_offset_hz().unwrap_or_default();
        assert_close(d1, FREQ_GLO_L1_DELTA, 0.0, "dF1");
        assert_close(d2, FREQ_GLO_L2_DELTA, 0.0, "dF2");
    }

    /// Every constellation the engine defines must produce a *defined* answer
    /// (mapped signal or explicit None) for every standard RINEX band — no
    /// silent GPS fallbacks, no panics. This is the coverage contract.
    #[test]
    fn all_engine_constellations_have_defined_band_answers() {
        let cons_all = [
            Constellation::Gps,
            Constellation::Glonass,
            Constellation::Galileo,
            Constellation::Beidou,
            Constellation::Sbas,
            Constellation::Qzss,
            Constellation::Navic,
        ];
        for cons in cons_all {
            for b in ["C1", "C2", "C5", "C6", "C7"] {
                if let Some(s) = rinex_type_to_signal(cons, b) {
                    let f = frequency_for(cons, s, 0);
                    assert!(f > 1.0e9 && f < 2.0e9, "sanity freq for {cons:?}/{b}");
                }
            }
        }
        // Spot-assert documented policy edges.
        assert_maps(Constellation::Navic, "C1", None);
        assert_maps(Constellation::Beidou, "C2", Some(Signal::BdsB1i));
        assert_maps(Constellation::Sbas, "C2", None);
    }

    #[test]
    fn test_secondary_signal_policy_resolves_e5a_e5b() {
        use crate::sat::Constellation;
        // Galileo secondary is E5a on band 5 — never E5b.
        assert_eq!(secondary_signal(Constellation::Galileo), Some((5, Signal::GalE5a)));
        // GPS secondary is L2 on band 2.
        assert_eq!(secondary_signal(Constellation::Gps), Some((2, Signal::GpsL2Cm)));
        assert_eq!(secondary_signal(Constellation::Glonass), Some((2, Signal::GloL2Of)));
        assert_eq!(secondary_signal(Constellation::Beidou), Some((7, Signal::BdsB2i)));
    }

    #[test]
    fn test_primary_signal_covers_all_supported() {
        use crate::sat::Constellation;
        for c in [Constellation::Gps, Constellation::Glonass, Constellation::Galileo, Constellation::Beidou] {
            assert!(primary_signal(c).is_some(), "{:?} missing primary", c);
        }
    }

    #[test]
    fn test_signal_for_band_resolves_each_observed_band() {
        use crate::sat::Constellation;
        assert_eq!(signal_for_band(Constellation::Gps, 1), Some(Signal::GpsL1Ca));
        assert_eq!(signal_for_band(Constellation::Gps, 2), Some(Signal::GpsL2Cm));
        // Galileo: whichever band an arc actually uses gets its own signal.
        assert_eq!(signal_for_band(Constellation::Galileo, 5), Some(Signal::GalE5a));
        assert_ne!(
            signal_for_band(Constellation::Galileo, 5),
            signal_for_band(Constellation::Galileo, 6)
        );
    }
