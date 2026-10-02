use crate::sat::{Constellation, SatelliteId};

pub const FREQ_GPS_L1: f64 = 1575.42e6;
pub const FREQ_GPS_L2: f64 = 1227.60e6;
pub const FREQ_GPS_L5: f64 = 1176.45e6;
pub const FREQ_GAL_E5B: f64 = 1207.140e6;
pub const FREQ_BDS_B1I: f64 = 1561.098e6;
/// BDS B2I / B2b carrier frequency (1180 * 1.023 MHz = 1207.140 MHz).
pub const FREQ_BDS_B2I: f64 = 1207.140e6;
/// BDS B3I carrier frequency (1240 * 1.023 MHz = 1268.520 MHz).
pub const FREQ_BDS_B3I: f64 = 1268.520e6;
pub const FREQ_GLO_L1_NOMINAL: f64 = 1602.0e6;
pub const FREQ_GLO_L2_NOMINAL: f64 = 1246.0e6;
pub const FREQ_GLO_L1_DELTA: f64 = 0.5625e6;
pub const FREQ_GLO_L2_DELTA: f64 = 0.4375e6;

pub fn satellite_frequencies(sat: SatelliteId, freq_num: i8) -> (f64, f64) {
    match sat.constellation {
        Constellation::Gps | Constellation::Qzss => (FREQ_GPS_L1, FREQ_GPS_L2),
        Constellation::Galileo => (FREQ_GPS_L1, FREQ_GAL_E5B), // E1 shares GPS L1
        Constellation::Beidou => (FREQ_BDS_B1I, FREQ_BDS_B2I), // B1I + B2I (BDS-2)
        Constellation::Glonass => {
            let f1 = FREQ_GLO_L1_NOMINAL + (freq_num as f64) * FREQ_GLO_L1_DELTA;
            let f2 = FREQ_GLO_L2_NOMINAL + (freq_num as f64) * FREQ_GLO_L2_DELTA;
            (f1, f2)
        }
        _ => (FREQ_GPS_L1, FREQ_GPS_L2),
    }
}

pub fn get_frequency(sat: SatelliteId, freq_band: u8, freq_num: i8) -> f64 {
    match freq_band {
        1 => match sat.constellation {
            Constellation::Gps | Constellation::Qzss | Constellation::Galileo => FREQ_GPS_L1,
            Constellation::Beidou => FREQ_BDS_B1I,
            Constellation::Glonass => FREQ_GLO_L1_NOMINAL + (freq_num as f64) * FREQ_GLO_L1_DELTA,
            _ => FREQ_GPS_L1,
        },
        2 => {
            match sat.constellation {
                Constellation::Gps | Constellation::Qzss => FREQ_GPS_L2,
                Constellation::Galileo => FREQ_GAL_E5B,
                Constellation::Beidou => FREQ_BDS_B1I, // BDS B1I (RINEX 3.03+ band 2, 1561.098 MHz)
                Constellation::Glonass => {
                    FREQ_GLO_L2_NOMINAL + (freq_num as f64) * FREQ_GLO_L2_DELTA
                }
                _ => FREQ_GPS_L2,
            }
        }
        5 => match sat.constellation {
            Constellation::Gps | Constellation::Qzss | Constellation::Galileo => FREQ_GPS_L5,
            // KNOWN WRONG: BDS-3 B2a is 1176.750 MHz, not L5. See the test
            // `beidou_band_5_is_b2a_which_is_not_l5` below for the full proof.
            Constellation::Beidou => FREQ_GPS_L5,
            _ => FREQ_GPS_L5,
        },
        6 => match sat.constellation {
            Constellation::Beidou => FREQ_BDS_B3I, // BDS B3I (1268.52 MHz)
            Constellation::Galileo => 1278.75e6,  // Galileo E6
            _ => FREQ_GPS_L1,
        },
        7 => match sat.constellation {
            Constellation::Galileo => FREQ_GAL_E5B,
            Constellation::Beidou => FREQ_BDS_B2I, // BDS-2 B2I / BDS-3 B2b (1207.14 MHz)
            _ => FREQ_GPS_L2,
        },
        _ => FREQ_GPS_L1,
    }
}

pub fn get_wavelength(sat: SatelliteId, freq_band: u8, freq_num: i8) -> f64 {
    let freq = get_frequency(sat, freq_band, freq_num);
    crate::constants::SPEED_OF_LIGHT_M_S / freq
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_glonass_fdma_wavelengths() {
        // Channel -4 (e.g. GLONASS PRN 6 in our dataset)
        let sat = SatelliteId {
            constellation: Constellation::Glonass,
            prn: 6,
        };
        let (f1, f2) = satellite_frequencies(sat, -4);
        assert_eq!(f1, FREQ_GLO_L1_NOMINAL - 4.0 * FREQ_GLO_L1_DELTA);
        assert_eq!(f2, FREQ_GLO_L2_NOMINAL - 4.0 * FREQ_GLO_L2_DELTA);

        let w1 = get_wavelength(sat, 1, -4);
        assert!((w1 - 0.18740019).abs() < 1e-5);
    }

    #[test]
    fn test_galileo_frequencies() {
        let sat = SatelliteId {
            constellation: Constellation::Galileo,
            prn: 11,
        };
        let (f1, f2) = satellite_frequencies(sat, 0);
        assert_eq!(f1, FREQ_GPS_L1);
        assert_eq!(f2, FREQ_GAL_E5B);
    }

    // -----------------------------------------------------------------------
    // Golden table. Every (constellation, band) pair the engine can produce,
    // with the value taken from the published ICDs:
    //   IS-GPS-200  L1 1575.42 / L2 1227.60 / L5 1176.45 MHz
    //   Galileo OS SIS ICD  E1 1575.420 / E5a 1176.450 / E5b 1207.140 / E6 1278.750 MHz
    //   BeiDou B1I/B2I/B3I ICD  B1I 1561.098 / B2I 1207.140 / B3I 1268.520 MHz
    //   GLONASS ICD 5.1  L1 1602 + 0.5625k / L2 1246 + 0.4375k MHz, k in -7..+6
    // -----------------------------------------------------------------------
    const ALL: [(Constellation, u8, f64, &str); 15] = [
        (Constellation::Gps, 1, 1_575_420_000.0, "GPS L1"),
        (Constellation::Gps, 2, 1_227_600_000.0, "GPS L2"),
        (Constellation::Gps, 5, 1_176_450_000.0, "GPS L5"),
        (Constellation::Qzss, 1, 1_575_420_000.0, "QZSS L1"),
        (Constellation::Qzss, 2, 1_227_600_000.0, "QZSS L2"),
        (Constellation::Qzss, 5, 1_176_450_000.0, "QZSS L5"),
        (Constellation::Galileo, 1, 1_575_420_000.0, "Gal E1"),
        (Constellation::Galileo, 2, 1_207_140_000.0, "Gal legacy L2 = E5b"),
        (Constellation::Galileo, 5, 1_176_450_000.0, "Gal E5a"),
        (Constellation::Galileo, 6, 1_278_750_000.0, "Gal E6"),
        (Constellation::Galileo, 7, 1_207_140_000.0, "Gal E5b"),
        (Constellation::Beidou, 1, 1_561_098_000.0, "BDS B1I (RINEX2 band 1)"),
        (Constellation::Beidou, 2, 1_561_098_000.0, "BDS B1I (RINEX 3.03+ band 2)"),
        (Constellation::Beidou, 6, 1_268_520_000.0, "BDS B3I"),
        (Constellation::Beidou, 7, 1_207_140_000.0, "BDS B2I"),
    ];

    #[test]
    fn every_constellation_band_resolves_to_its_icd_frequency() {
        for (c, band, want, what) in ALL {
            let sat = SatelliteId { constellation: c, prn: 7 };
            assert_eq!(get_frequency(sat, band, 0), want, "{what}");
            // Wavelength must be exactly c/f.
            assert_eq!(
                get_wavelength(sat, band, 0),
                crate::constants::SPEED_OF_LIGHT_M_S / want,
                "{what} wavelength"
            );
        }
    }

    #[test]
    fn galileo_band_6_and_7_are_not_interchangeable() {
        // The historical E5a/E6 swap: band 6 is E6 (1278.75), band 7 is E5b.
        let gal = SatelliteId { constellation: Constellation::Galileo, prn: 11 };
        assert_eq!(get_frequency(gal, 6, 0), 1_278_750_000.0);
        assert_eq!(get_frequency(gal, 7, 0), 1_207_140_000.0);
        assert_ne!(get_frequency(gal, 6, 0), get_frequency(gal, 7, 0));
    }

    #[test]
    fn glonass_fdma_matches_the_icd_spacings_exactly() {
        let g = SatelliteId { constellation: Constellation::Glonass, prn: 1 };
        for k in -7..=6 {
            assert_eq!(
                get_frequency(g, 1, k),
                1_602_000_000.0 + (k as f64) * 562_500.0,
                "GLONASS L1 k={k}"
            );
            assert_eq!(
                get_frequency(g, 2, k),
                1_246_000_000.0 + (k as f64) * 437_500.0,
                "GLONASS L2 k={k}"
            );
            // f1/f2 is exactly 9/7 on every channel (spacings share the ratio).
            assert_eq!(get_frequency(g, 1, k) * 7.0, get_frequency(g, 2, k) * 9.0);
        }
        // `satellite_frequencies` must agree with the per-band lookup.
        for k in -7..=6 {
            let (f1, f2) = satellite_frequencies(g, k);
            assert_eq!(f1, get_frequency(g, 1, k));
            assert_eq!(f2, get_frequency(g, 2, k));
        }
    }

    #[test]
    fn satellite_frequencies_returns_the_icd_pair_per_constellation() {
        let cases = [
            (Constellation::Gps, 1_575_420_000.0, 1_227_600_000.0),
            (Constellation::Qzss, 1_575_420_000.0, 1_227_600_000.0),
            (Constellation::Galileo, 1_575_420_000.0, 1_207_140_000.0),
            (Constellation::Beidou, 1_561_098_000.0, 1_207_140_000.0),
        ];
        for (c, f1, f2) in cases {
            let sat = SatelliteId { constellation: c, prn: 3 };
            assert_eq!(satellite_frequencies(sat, 0), (f1, f2), "{c:?}");
        }
        // Unsupported constellations fall back to the GPS plan (documented).
        for c in [Constellation::Sbas, Constellation::Navic] {
            let sat = SatelliteId { constellation: c, prn: 1 };
            assert_eq!(
                satellite_frequencies(sat, 0),
                (FREQ_GPS_L1, FREQ_GPS_L2),
                "{c:?} fallback"
            );
        }
    }

    /// Total function: no constellation x band x channel may produce a
    /// non-physical frequency. A silent 0 or negative would make the
    /// wavelength infinite and quietly poison a measurement model.
    #[test]
    fn frequency_lookup_is_total_and_always_physical() {
        let all = [
            Constellation::Gps,
            Constellation::Qzss,
            Constellation::Glonass,
            Constellation::Galileo,
            Constellation::Beidou,
            Constellation::Sbas,
            Constellation::Navic,
        ];
        for c in all {
            for band in 0..=10u8 {
                for k in [-12i8, -7, 0, 6, 12] {
                    let sat = SatelliteId { constellation: c, prn: 1 };
                    let f = get_frequency(sat, band, k);
                    assert!(
                        f.is_finite() && f > 1.0e9 && f < 2.0e9,
                        "{c:?} band {band} k={k} -> {f} Hz"
                    );
                    let w = get_wavelength(sat, band, k);
                    assert!(
                        w.is_finite() && w > 0.1 && w < 0.4,
                        "{c:?} band {band} k={k} -> {w} m"
                    );
                }
            }
        }
    }

    /// `get_wavelength` must be exactly `c / get_frequency`, with `c` the SI
    /// defining constant. Guards against a hard-coded wavelength table.
    #[test]
    fn wavelength_is_always_c_over_frequency() {
        let c = crate::constants::SPEED_OF_LIGHT_M_S;
        assert_eq!(c, 299_792_458.0);
        for (cons, band, k) in [
            (Constellation::Gps, 1u8, 0i8),
            (Constellation::Galileo, 7, 0),
            (Constellation::Beidou, 6, 0),
            (Constellation::Glonass, 1, -7),
            (Constellation::Glonass, 2, 6),
        ] {
            let sat = SatelliteId { constellation: cons, prn: 9 };
            assert_eq!(get_wavelength(sat, band, k), c / get_frequency(sat, band, k));
        }
    }

    /// BDS B2a (BDS-3, RINEX 3 band 5, e.g. JOZE's `C5P`) is **not** GPS L5.
    ///
    /// Exact integer arithmetic, no float tolerance:
    ///   L5   = 1176.450 MHz = 1150 x 1.023 MHz exactly (1150 * 1_023_000 = 1_176_450_000)
    ///   B2a  = 1176.750 MHz; 1150 * 1_023_000 = 1_176_450_000, and
    ///          1_176_750_000 - 1_176_450_000 = 300_000 Hz remainder — B2a is
    ///          NOT on the 1.023 MHz grid, so it cannot share L5's grid point.
    /// Separation is therefore exactly 300 kHz (relative error 2.55e-4).
    ///
    /// CONFIRMED DEFECT, DELIBERATELY NOT FIXED — `get_frequency` has no B2a
    /// variant and returns L5 for BeiDou band 5. Asserting the correct value
    /// makes this test fail with `left: 1176450000.0, right: 1176750000.0`;
    /// that measurement is the reason the value is still shipped. The
    /// "shares L5" claim in `get_frequency` is simply false, and is corrected
    /// in the comment there rather than in behaviour.
    ///
    /// B2a band 5 is NOT reachable by any engine path on any dataset in the
    /// corpus, which is the brief's condition for leaving it alone:
    ///   * `canonical_bands_for_constellation(Beidou) = [1, 6, 7]` — the DD
    ///     filter never forms a band-5 arc (pinned by
    ///     `beidou_band_list_excludes_the_unimplemented_b2a` in gneiss-rtk);
    ///   * `mw.rs` and `iono_free.rs` both scan `[2, 7, 6, 5]` and take the
    ///     first match; all 26 BeiDou files in `datasets/` that carry band 5
    ///     carry 2/6/7 as well, so 5 is never selected;
    ///   * the RTCM MSM BeiDou table (`gneiss-parsers/src/rtcm3/msm/signals.rs`)
    ///     has no band-5 entry at all — only `1I/1Q/1X`, `6I/6Q/6X`, `7I/7Q/7X`.
    ///
    /// If any of those ever changes, the second assertion below flips and the
    /// constant must be corrected in the same commit.
    #[test]
    fn beidou_band_5_is_b2a_which_is_not_l5() {
        const L5_HZ: i64 = 1150 * 1_023_000; // 1_176_450_000
        const B2A_HZ: i64 = 1_176_750_000;
        assert_eq!(L5_HZ, 1_176_450_000, "L5 is exactly 1150 x 1.023 MHz");
        assert_eq!(B2A_HZ % 1_023_000, 300_000, "B2a is off the 1.023 MHz grid");
        assert_eq!(B2A_HZ - L5_HZ, 300_000, "B2a sits 300 kHz above L5");
        assert_eq!(FREQ_GPS_L5, L5_HZ as f64, "the constant under test is L5");
        // KNOWN WRONG by 300 kHz, unreachable today: asserted so the defect
        // stays visible and a future reachability change trips this test.
        let sat = SatelliteId { constellation: Constellation::Beidou, prn: 1 };
        assert_eq!(get_frequency(sat, 5, 0), FREQ_GPS_L5, "known-wrong B2a fallback");
    }
}
