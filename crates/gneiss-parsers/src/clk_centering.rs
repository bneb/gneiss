//! Constellation-median centering of RINEX clock products, with spread
//! gating, for precise-clock DD corrections.
//!
//! Precise clock products carry an arbitrary per-epoch timescale datum:
//! every satellite's absolute bias can sit hundreds of microseconds off the
//! true system time, and different products (or reference-switch events in
//! downstream pairing) expose that datum differently. Centering each bias by
//! the MEDIAN across its constellation mates removes the common mode before
//! any differencing.
//!
//! After centering we compute the inter-satellite spread (max - min of the
//! centered biases). Healthy products collapse to sub-microsecond spreads;
//! a spread above [`MAX_CENTERED_SPREAD_S`] marks the product epoch as
//! pathological and callers must suppress the differential correction.
//!
//! Note an invariant that guides the design: subtracting a per-epoch
//! constant from every satellite never changes a pairwise difference. What
//! centering buys is (a) a meaningful spread statistic for gating and
//! (b) degenerate-set handling (< [`MIN_CENTERING_MATES`] mates yields no
//! correction instead of an ungated raw differential).

use crate::rinex_clk::RinexClock;
use gneiss_core::sat::SatelliteId;
use gneiss_core::time::GpsTime;

/// Minimum number of valid constellation biases (including the queried
/// satellite's own) required to trust centering. Below this the "median"
/// degenerates to one or two satellites and neither centering nor spread
/// gating is statistically meaningful; [`RinexClock::centered_clock`]
/// returns `None` so callers fall back to no correction.
pub const MIN_CENTERING_MATES: usize = 3;

/// Maximum tolerated inter-satellite spread (max - min) of centered biases,
/// in seconds. Healthy rapid products stay far below this after centering;
/// anything wider means the product's per-satellite biases disagree by more
/// than ~30 km of range and the DD correction must be suppressed for the
/// epoch (`c · 100 µs ≈ 30 km`).
pub const MAX_CENTERED_SPREAD_S: f64 = 100e-6;

/// Centered clock state of one satellite at one instant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CenteredClock {
    /// Satellite bias minus the constellation median bias (seconds).
    pub bias_s: f64,
    /// max - min of the centered biases across all valid constellation
    /// mates (seconds). Datum-invariant: identical to the raw
    /// inter-satellite spread because the median cancels.
    pub spread_s: f64,
}

/// Standard median of a slice: middle element for odd counts, MEAN OF THE
/// TWO MIDDLE elements for even counts (documented policy — the textbook
/// median minimising `Σ|x − m|`, not the lower median). Returns `None` for
/// empty input. The choice only shifts the removed common mode by at most
/// half the central gap; pairwise differences are identical either way.
pub fn sorted_median(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("clock biases are never NaN"));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 1 {
        Some(sorted[mid])
    } else {
        Some((sorted[mid - 1] + sorted[mid]) / 2.0)
    }
}

impl RinexClock {
    /// Bias of `sat` at time `t`, centered by the median over all
    /// constellation mates with valid records near `t`.
    ///
    /// Returns `None` when the satellite itself has no valid record or
    /// fewer than [`MIN_CENTERING_MATES`] constellation biases are valid.
    /// Validity follows [`RinexClock::get_clock_bias`] semantics (records
    /// more than 900 s away do not count).
    pub fn centered_clock(&self, sat: SatelliteId, t: GpsTime) -> Option<CenteredClock> {
        let own = self.get_clock_bias(sat, t)?;
        let biases: Vec<f64> = self
            .satellites
            .keys()
            .filter(|mate| mate.constellation == sat.constellation)
            .filter_map(|mate| self.get_clock_bias(*mate, t))
            .collect();
        if biases.len() < MIN_CENTERING_MATES {
            return None;
        }
        let lo = biases.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = biases.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let median = sorted_median(&biases)?;
        Some(CenteredClock { bias_s: own - median, spread_s: hi - lo })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rinex_clk::ClockRecord;
    use gneiss_core::sat::Constellation;
    use std::collections::HashMap;

    const WEEK: u32 = 2382;

    fn gps(prn: u8) -> SatelliteId {
        SatelliteId { constellation: Constellation::Gps, prn }
    }

    fn gal(prn: u8) -> SatelliteId {
        SatelliteId { constellation: Constellation::Galileo, prn }
    }

    fn record(tow: f64, bias_us: f64) -> ClockRecord {
        ClockRecord { time: GpsTime::new(WEEK, tow), bias: bias_us * 1e-6 }
    }

    /// Product with one record per satellite, all stamped at tow = 0,
    /// biases given in microseconds.
    fn clk_with(biases_us: &[(SatelliteId, f64)]) -> RinexClock {
        RinexClock {
            satellites: biases_us
                .iter()
                .map(|(sat, us)| (*sat, vec![record(0.0, *us)]))
                .collect::<HashMap<_, _>>(),
        }
    }

    fn t0() -> GpsTime {
        GpsTime::new(WEEK, 0.0)
    }

    // ---- median policy ---------------------------------------------------

    #[test]
    fn median_odd_count_picks_middle_element() {
        assert_eq!(sorted_median(&[3.0, 1.0, 2.0]), Some(2.0));
        assert_eq!(sorted_median(&[5.0]), Some(5.0));
        assert_eq!(sorted_median(&[-7.0, 7.0, 0.0]), Some(0.0));
    }

    #[test]
    fn median_even_count_averages_two_middle_elements() {
        // Documented even-count policy: average of the two central order
        // statistics, NOT the lower median ({1,2,3,4} -> 2.5, not 2.0).
        assert_eq!(sorted_median(&[4.0, 1.0, 3.0, 2.0]), Some(2.5));
        assert_eq!(sorted_median(&[-600.0, -590.0, 590.0, 600.0]), Some(0.0));
        assert_eq!(sorted_median(&[]), None);
    }

    // ---- mission case 1: healthy cluster under a large common mode -------

    /// Four satellites sharing a ~+500 µs product datum: after centering
    /// every bias collapses within ±25 µs, the pairwise DD delta between
    /// two sats is tiny, and the spread gate passes.
    #[test]
    fn centering_collapses_common_mode_and_gate_passes() {
        let rc = clk_with(&[
            (gps(28), 500.0),
            (gps(27), 520.0),
            (gps(5), 480.0),
            (gps(10), 510.0),
        ]);
        for prn in [28u8, 27, 5, 10] {
            let c = rc.centered_clock(gps(prn), t0()).expect("4 healthy mates");
            // 1 ps slack absorbs f64 rounding at the exact boundary value.
            assert!(
                c.bias_s.abs() <= 25e-6 + 1e-12,
                "PRN {prn} centered to {:.15} us, expected |x| <= 25 us",
                c.bias_s * 1e6
            );
            assert!(
                c.spread_s <= MAX_CENTERED_SPREAD_S,
                "spread {} us must pass the gate",
                c.spread_s * 1e6
            );
        }
        let a = rc
            .centered_clock(gps(27), t0())
            .expect("healthy cluster")
            .bias_s;
        let b = rc
            .centered_clock(gps(28), t0())
            .expect("healthy cluster")
            .bias_s;
        assert!((a - b).abs() <= 25e-6, "pair delta {} us must be small", (a - b) * 1e6);
        // Datum-invariance: centering cannot alter a pairwise delta.
        let raw_delta = (520.0 - 500.0) * 1e-6;
        assert!((raw_delta - (a - b)).abs() < 1e-15);
    }

    /// The literal mission set keeps one negative outlier {-480 µs}: NO
    /// median can pull it within ±25 µs of the cluster (the example numbers
    /// are internally inconsistent under any median definition). The honest
    /// outcome for that set is exactly what the gate exists for: the
    /// outlier pushes the spread past 100 µs.
    #[test]
    fn literal_mission_set_with_outlier_blows_past_gate() {
        let rc = clk_with(&[
            (gps(28), 500.0),
            (gps(27), 520.0),
            (gps(5), -480.0),
            (gps(10), 510.0),
        ]);
        let outlier = rc.centered_clock(gps(5), t0()).expect("mates present");
        assert!(
            outlier.bias_s.abs() > 25e-6,
            "no median can center -480 into the 480..520 cluster"
        );
        assert!(outlier.spread_s > MAX_CENTERED_SPREAD_S);
    }

    // ---- mission case 2: pathological symmetric set trips the gate -------

    #[test]
    fn pathological_symmetric_biases_trip_spread_gate() {
        let rc = clk_with(&[
            (gps(1), 600.0),
            (gps(2), -600.0),
            (gps(3), 590.0),
            (gps(4), -590.0),
        ]);
        let c = rc.centered_clock(gps(1), t0()).expect("4 mates");
        // Even-count median of {-600,-590,590,600} is (-590+590)/2 = 0.
        assert!((c.bias_s - 600e-6).abs() < 1e-15);
        assert!((c.spread_s - 1200e-6).abs() < 1e-12, "max-min = 1200 us");
        assert!(c.spread_s > MAX_CENTERED_SPREAD_S, "gate must trip");
    }

    // ---- degenerate sets (red-team answers codified) ---------------------

    #[test]
    fn fewer_than_three_valid_mates_returns_none() {
        let pair = clk_with(&[(gps(1), 100.0), (gps(2), 200.0)]);
        assert_eq!(pair.centered_clock(gps(1), t0()), None);
        let solo = clk_with(&[(gps(7), 100.0)]);
        assert_eq!(solo.centered_clock(gps(7), t0()), None);
    }

    #[test]
    fn missing_own_record_returns_none() {
        let rc = clk_with(&[(gps(1), 100.0), (gps(2), 200.0), (gps(3), 300.0)]);
        assert_eq!(rc.centered_clock(gps(9), t0()), None);
    }

    #[test]
    fn other_constellations_are_not_mates() {
        let rc = clk_with(&[
            (gps(1), 100.0),
            (gps(2), 140.0),
            (gal(1), 0.0),
            (gal(2), 10.0),
            (gal(3), 20.0),
            (gal(4), 30.0),
        ]);
        // GPS has only 2 valid mates -> below minimum -> None ...
        assert_eq!(rc.centered_clock(gps(1), t0()), None);
        // ... while Galileo centers normally on its own 4 members.
        let g = rc.centered_clock(gal(2), t0()).expect("4 galileo mates");
        assert!(g.spread_s <= MAX_CENTERED_SPREAD_S);
    }

    #[test]
    fn stale_records_outside_900s_window_are_not_mates() {
        let mut rc = RinexClock::default();
        rc.satellites.insert(gps(1), vec![record(0.0, 50.0)]);
        // These two records are 3600 s from the query time: get_clock_bias
        // rejects them, dropping the valid-mate count below the minimum.
        rc.satellites.insert(gps(2), vec![record(3600.0, 50.0)]);
        rc.satellites.insert(gps(3), vec![record(3600.0, 50.0)]);
        assert_eq!(rc.centered_clock(gps(1), GpsTime::new(WEEK, 0.0)), None);
    }

    // ---- datum invariance -------------------------------------------------

    #[test]
    fn spread_is_invariant_to_common_mode_datum_shift() {
        let base = [600.0f64, -600.0, 590.0, -590.0];
        let mk = |shift: f64| {
            clk_with(
                &base
                    .iter()
                    .enumerate()
                    .map(|(i, b)| (gps(i as u8 + 1), b + shift))
                    .collect::<Vec<_>>(),
            )
        };
        let a = mk(0.0)
            .centered_clock(gps(1), t0())
            .expect("4 mates")
            .spread_s;
        let b = mk(1234.0)
            .centered_clock(gps(1), t0())
            .expect("4 mates")
            .spread_s;
        assert!((a - b).abs() < 1e-15, "shift must not change spread");
        assert!((a - 1200e-6).abs() < 1e-12);
    }
}

#[cfg(test)]
mod gfz_real_file_tests {
    use super::*;
    use crate::rinex_clk::gfz_clk_cached_path;
    use crate::rinex_clk::RinexClock;
    use gneiss_core::sat::Constellation;

    /// On the real DOY160 GFZ rapid product: for a sampled GPS epoch where
    /// centering succeeds, the centered pairwise delta between two
    /// satellites equals their raw delta exactly (the wiring keeps raw
    /// differences), and the reported spread is non-negative and equals the
    /// raw inter-satellite spread. Informational stats printed once.
    #[test]
    fn gfz_centering_preserves_pairwise_deltas_and_measures_spread() {
        let Some(path) = gfz_clk_cached_path() else {
            eprintln!("skipping: GFZ CLK product not cached locally");
            return;
        };
        let content = std::fs::read_to_string(path).expect("GFZ CLK readable");
        let rc = RinexClock::parse(&content);
        let mut gps_sats: Vec<SatelliteId> = rc
            .satellites
            .keys()
            .filter(|s| s.constellation == Constellation::Gps)
            .copied()
            .collect();
        gps_sats.sort_by_key(|s| s.prn);
        assert!(gps_sats.len() >= 20, "expected the full GPS constellation");

        let first_rec = gps_sats
            .iter()
            .filter_map(|s| rc.satellites[s].first())
            .map(|r| r.time.tow)
            .fold(f64::NAN, f64::max);
        let t = GpsTime::new(
            rc.satellites[&gps_sats[0]][0].time.week,
            first_rec,
        );

        let a = gps_sats[0];
        let b = gps_sats[1];
        let (ra, rb) = (
            rc.get_clock_bias(a, t).expect("record exists"),
            rc.get_clock_bias(b, t).expect("record exists"),
        );
        let ca = rc.centered_clock(a, t).expect("full constellation");
        let cb = rc.centered_clock(b, t).expect("full constellation");
        let raw_pair_delta = ra - rb;
        let centered_pair_delta = ca.bias_s - cb.bias_s;
        assert!(
            (raw_pair_delta - centered_pair_delta).abs() < 1e-18,
            "centering must preserve pairwise deltas"
        );
        let raw_spread = {
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for s in &gps_sats {
                if let Some(v) = rc.get_clock_bias(*s, t) {
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
            }
            hi - lo
        };
        assert!((raw_spread - ca.spread_s).abs() < 1e-15);
        assert!(ca.spread_s >= 0.0);
        println!(
            "GFZ DOY160 @ tow {:.0}: GPS spread {:.1} us (gate {:.0} us), \
             pair {:?}{:02}-{:02} delta {:.1} us",
            t.tow,
            ca.spread_s * 1e6,
            MAX_CENTERED_SPREAD_S * 1e6,
            a.constellation,
            a.prn,
            b.prn,
            centered_pair_delta * 1e6
        );
    }
}
