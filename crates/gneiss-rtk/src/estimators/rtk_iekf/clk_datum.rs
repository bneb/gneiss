//! Per-arc precise-clock datum bookkeeping for DD ambiguities.
//!
//! The engine subtracts `c·(dt_sat − dt_ref)` model-side every epoch
//! ([`GnssRtkIekf::precise_clock_dd_m`]). A DD arc's float ambiguity is
//! only consistent with that correction while the arc's `(sat, ref)`
//! pairing — and therefore the correction's datum — stays fixed. This
//! module records, per [`DoubleDiffKey`], the clock datum
//! `c·(dt_sat − dt_ref)` (metres) in force when the arc's ambiguity was
//! last (re)seeded, and keeps the float continuous across events that
//! change the datum:
//!
//! - **Reference switch** (`ref_sat` changes while the satellite keeps
//!   tracking, no cycle slip): the converged float under the old key is
//!   transferred to the new key and stepped by
//!   `(new_datum − old_datum) / λ` cycles via
//!   [`RtkState::adjust_ambiguity`] *before* any innovation is formed,
//!   so the arc does not restart from a noisy code-minus-phase seed.
//! - **Cycle slip** (LLI or GF arc break): the ambiguity resets raw, so
//!   the stored datum resets with it.
//!
//! Without a precise-clock product the whole map stays empty: every
//! correction is 0.0 and this bookkeeping is inert by construction.

use super::state::DoubleDiffKey;
use super::GnssRtkIekf;
use gneiss_core::constants::SPEED_OF_LIGHT_M_S;

/// Pure pair decision for the centered precise-clock DD correction.
pub(crate) fn centered_pair_correction(
    sat: Option<gneiss_parsers::clk_centering::CenteredClock>,
    reference: Option<gneiss_parsers::clk_centering::CenteredClock>,
) -> (f64, bool) {
    use gneiss_parsers::clk_centering::MAX_CENTERED_SPREAD_S;
    const C: f64 = SPEED_OF_LIGHT_M_S;
    match (sat, reference) {
        (Some(a), Some(b)) => {
            let tripped = a.spread_s > MAX_CENTERED_SPREAD_S || b.spread_s > MAX_CENTERED_SPREAD_S;
            if tripped {
                (0.0, true)
            } else {
                (C * (a.bias_s - b.bias_s), false)
            }
        }
        _ => (0.0, false),
    }
}

/// `GNEISS_CLK_TRACE` diagnostic: first 10 evaluations show the centered
/// biases, spreads, and resulting pair correction (metres).
pub(crate) fn clk_centering_trace(
    tow: f64,
    sat_id: gneiss_core::sat::SatelliteId,
    ref_sv: gneiss_core::sat::SatelliteId,
    cs: Option<gneiss_parsers::clk_centering::CenteredClock>,
    cr: Option<gneiss_parsers::clk_centering::CenteredClock>,
    corr_m: f64,
) {
    static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    if n >= 10 {
        return;
    }
    let fmt = |c: Option<gneiss_parsers::clk_centering::CenteredClock>| match c {
        Some(x) => format!("{:+.1}us/spr{:.1}us", x.bias_s * 1e6, x.spread_s * 1e6),
        None => "none".to_string(),
    };
    eprintln!(
        "CLKDD[{n}] tow={:.0} {:?}{:02}-{:02} sat={} ref={} -> {:+.3} m",
        tow,
        sat_id.constellation,
        sat_id.prn,
        ref_sv.prn,
        fmt(cs),
        fmt(cr),
        corr_m
    );
}

/// Seed variance (m^2) for the iono state of a transferred arc; matches
/// the fresh-seed value used in `update_dd_ambiguity`.
const TRANSPLANT_IONO_VAR_M2: f64 = 4.0;

impl GnssRtkIekf {
    /// Clock-datum bookkeeping for one DD pair. Call during measurement
    /// formation, BEFORE the ambiguity update touches the state.
    ///
    /// `datum_m` is the epoch's `c·(dt_sat − dt_ref)` correction in
    /// metres (0.0 without a clock product), `lambda` the pair's
    /// wavelength, `lli_slip` whether the arc just broke.
    pub(crate) fn apply_clock_datum(
        &mut self,
        key: DoubleDiffKey,
        datum_m: f64,
        lambda: f64,
        lli_slip: bool,
    ) {
        if self.precise_clocks.is_none() {
            return; // Inert without the product: no datum exists to track.
        }
        if lli_slip {
            // Genuine cycle slip: the arc restarts from a raw
            // code-minus-phase seed, so its datum restarts with it.
            self.clk_datum_m.insert(key, datum_m);
            return;
        }
        if self.has_live_arc(key) {
            // Established arc: the model-side correction tracks the clocks
            // epoch by epoch; nothing for the float to absorb.
            return;
        }
        self.transplant_across_ref_switch(key, datum_m, lambda);
        self.clk_datum_m.insert(key, datum_m);
    }

    /// True when `key` has both a recorded datum and a live float, i.e.
    /// the arc is established and must not be disturbed.
    fn has_live_arc(&self, key: DoubleDiffKey) -> bool {
        self.clk_datum_m.contains_key(&key) && self.state.get_amb_idx(&key).is_some()
    }

    /// Predecessor arc for the same satellite/band under a DIFFERENT
    /// reference satellite, provided it still carries a live float we can
    /// transfer. At most one such key can exist at a time.
    fn ref_switch_predecessor(&self, key: DoubleDiffKey) -> Option<DoubleDiffKey> {
        self.clk_datum_m
            .keys()
            .find(|k| {
                k.constellation_id == key.constellation_id
                    && k.sat == key.sat
                    && k.freq_band == key.freq_band
                    && k.ref_sat != key.ref_sat
                    && self.state.get_amb_idx(k).is_some()
            })
            .copied()
    }

    /// Move a converged float across a reference-satellite switch.
    ///
    /// The new key's ambiguity is seeded from the old arc's float and
    /// stepped by the datum change so the predicted phase stays
    /// continuous across the switch; the old arc's bookkeeping entry is
    /// consumed. No-op when no live predecessor exists (plain fresh seed).
    fn transplant_across_ref_switch(&mut self, key: DoubleDiffKey, datum_m: f64, lambda: f64) {
        let Some(old_key) = self.ref_switch_predecessor(key) else {
            return;
        };
        let Some(old_datum) = self.clk_datum_m.get(&old_key).copied() else {
            return;
        };
        let Some(old_idx) = self.state.get_amb_idx(&old_key) else {
            return;
        };
        let rel = old_idx - self.state.amb_offset();
        let n_prev = self.state.ambiguities[rel].1;
        let var_prev = self.state.cov[(old_idx, old_idx)];
        let step_cycles = (datum_m - old_datum) / lambda;
        self.state.ensure_ambiguity(key, n_prev, var_prev);
        self.state.adjust_ambiguity(&key, step_cycles);
        self.state.ensure_iono(key, TRANSPLANT_IONO_VAR_M2);
        self.clk_datum_m.remove(&old_key);
        self.clk_ref_switches += 1;
        if std::env::var("GNEISS_CLK_TRACE").is_ok() {
            eprintln!(
                "CLKSWITCH[{}] tow={:.0} cons={} sat={:02} b{} ref {}->{} datum {:+.3}->{:+.3} m step {:+.4} cyc",
                self.clk_ref_switches, self.state.time.tow, key.constellation_id,
                key.sat, key.freq_band, old_key.ref_sat, key.ref_sat,
                old_datum, datum_m, step_cycles
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::estimators::rtk_iekf::state::RtkState;
    use gneiss_core::time::GpsTime;
    use nalgebra::Vector3;

    const LAMBDA: f64 = 0.19029367279836515;

    fn engine_with_clocks() -> GnssRtkIekf {
        let mut eng = GnssRtkIekf::new(
            Vector3::new(1.0, 2.0, 3.0),
            GpsTime::new(2200, 100.0),
            1.0,
        );
        eng.precise_clocks =
            Some(std::sync::Arc::new(gneiss_parsers::rinex_clk::RinexClock::parse("")));
        eng
    }

    fn key(sat: u16, ref_sat: u16, band: u8) -> DoubleDiffKey {
        DoubleDiffKey { constellation_id: 0, sat, ref_sat, freq_band: band }
    }

    #[test]
    fn ref_switch_transfers_float_stepped_by_exactly_the_datum_change() {
        let mut eng = engine_with_clocks();
        let (old, new) = (key(5, 1, 1), key(5, 2, 1));
        let (n_prev, d_a, d_b) = (50.25_f64, 30.0_f64, 45.0_f64);
        eng.state.ensure_ambiguity(old, n_prev, 100.0);
        eng.clk_datum_m.insert(old, d_a);
        eng.apply_clock_datum(new, d_b, LAMBDA, false);

        // Stored ambiguity increased by exactly (dB − dA)/λ ...
        let idx = eng.state.get_amb_idx(&new).expect("transferred float must exist");
        let got = eng.state.ambiguities[idx - eng.state.amb_offset()].1;
        let expected = n_prev + (d_b - d_a) / LAMBDA;
        assert!((got - expected).abs() < 1e-12, "got {got}, want {expected}");
        // ... variance inherited (continuity, not re-randomisation) ...
        assert!((eng.state.cov[(idx, idx)] - 100.0).abs() < 1e-12);
        // ... datum updated and predecessor consumed.
        assert_eq!(eng.clk_datum_m.get(&new), Some(&d_b));
        assert!(!eng.clk_datum_m.contains_key(&old));
        // ... and this engine counted exactly one switch.
        assert_eq!(eng.clk_ref_switches, 1, "one transfer must be counted");
    }

    #[test]
    fn ref_switch_transfer_keeps_phase_innovation_continuous_by_construction() {
        // Same geometry, same underlying integer N_int; observations flip
        // with the pairing's datum exactly as the model does. The
        // pre-switch and post-switch innovations must be identical.
        let mut eng = engine_with_clocks();
        let (old, new) = (key(5, 1, 1), key(5, 2, 1));
        let (n_int, d_a, d_b) = (50.25_f64, 30.0_f64, 45.0_f64);
        // Arc A converged under its own datum-consistent obs/model pair.
        eng.state.ensure_ambiguity(old, n_int + d_a / LAMBDA, 100.0);
        eng.clk_datum_m.insert(old, d_a);

        let innov = |datum: f64, amb: f64| n_int + datum / LAMBDA - amb;
        let innov_before = innov(d_a, n_int + d_a / LAMBDA);
        eng.apply_clock_datum(new, d_b, LAMBDA, false);
        let idx = eng.state.get_amb_idx(&new).unwrap();
        let innov_after = innov(d_b, eng.state.ambiguities[idx - eng.state.amb_offset()].1);
        assert!(
            (innov_after - innov_before).abs() < 1e-12,
            "innovation jumped {innov_before} -> {innov_after}"
        );
        assert!(innov_after.abs() < 1e-9, "no spike possible: {innov_after}");
    }

    #[test]
    fn fresh_seed_without_predecessor_only_records_datum() {
        let mut eng = engine_with_clocks();
        let k = key(7, 1, 1);
        eng.clk_datum_m.insert(k, 3.0);
        // Datum recorded but the float was dropped (sat re-rise): the arc
        // restarts, so only the datum refreshes — no phantom transplant.
        eng.apply_clock_datum(k, 4.0, LAMBDA, false);
        assert_eq!(eng.clk_datum_m.get(&k), Some(&4.0));
        assert!(eng.state.get_amb_idx(&k).is_none(), "bookkeeping never seeds floats");
        eng.apply_clock_datum(k, 5.0, LAMBDA, false);
        assert_eq!(eng.clk_datum_m.get(&k), Some(&5.0));
        assert_eq!(eng.clk_ref_switches, 0);
    }

    #[test]
    fn cycle_slip_resets_both_ambiguity_and_datum() {
        let mut eng = engine_with_clocks();
        let (old, new) = (key(5, 1, 1), key(5, 2, 1));
        eng.state.ensure_ambiguity(old, 50.25, 100.0);
        eng.clk_datum_m.insert(old, 30.0);

        // Slip ON THE NEW PAIRING (ref already switched): raw seed path,
        // datum recorded, never a transplant.
        eng.state.ensure_ambiguity(new, 60.0, 100.0);
        eng.apply_clock_datum(new, 45.0, LAMBDA, true);
        assert_eq!(eng.clk_datum_m.get(&new), Some(&45.0));

        // Slip on an ESTABLISHED arc: datum refreshed alongside the raw
        // re-seed performed by update_dd_ambiguity (verified there).
        eng.apply_clock_datum(old, 31.0, LAMBDA, true);
        assert_eq!(eng.clk_datum_m.get(&old), Some(&31.0));
        assert_eq!(
            eng.clk_datum_m.get(&new),
            Some(&45.0),
            "datum-only entry survives; no float moved"
        );
        eng.apply_clock_datum(old, 32.0, LAMBDA, false);
        assert_eq!(eng.clk_ref_switches, 0, "slips are not switches");
    }

    #[test]
    fn no_clock_product_leaves_bookkeeping_inert() {
        let mut eng = GnssRtkIekf::new(
            Vector3::new(1.0, 2.0, 3.0),
            GpsTime::new(2200, 100.0),
            1.0,
        );
        assert!(eng.precise_clocks.is_none());
        let (old, new) = (key(5, 1, 1), key(5, 2, 1));
        eng.state.ensure_ambiguity(old, 50.25, 100.0);

        eng.apply_clock_datum(new, 45.0, LAMBDA, false);
        eng.apply_clock_datum(old, 31.0, LAMBDA, true);

        assert!(eng.clk_datum_m.is_empty(), "map untouched without clocks");
        assert!(!eng.state.get_amb_idx(&new).is_some());
        let idx = eng.state.get_amb_idx(&old).unwrap();
        let v = eng.state.ambiguities[idx - eng.state.amb_offset()].1;
        assert_eq!(v, 50.25, "float untouched");
        eng.apply_clock_datum(new, 46.0, LAMBDA, false);
        assert_eq!(eng.clk_ref_switches, 0);
    }

    #[test]
    fn adjust_ambiguity_adds_to_stored_float_without_touching_covariance() {
        let mut st = RtkState::new(Vector3::zeros(), GpsTime::new(2100, 0.0));
        let k = key(5, 1, 1);
        st.ensure_ambiguity(k, 10.0, 25.0);
        assert!(st.adjust_ambiguity(&k, 1.5));
        assert_eq!(st.ambiguities[0].1, 11.5);
        assert!((st.cov[(st.get_amb_idx(&k).unwrap(), st.get_amb_idx(&k).unwrap())] - 25.0).abs() < 1e-12);
        assert!(!st.adjust_ambiguity(&key(9, 1, 1), 1.0), "missing key reports no-op");
    }
}
