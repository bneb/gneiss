//! Carrier-phase, ambiguity-arc and slant-ionosphere tests for the UDUC
//! factor builder.  Split from `uduc_builder_tests.rs` to stay under the
//! 500-line file limit; the shared fixtures live there.
#![allow(clippy::unwrap_used)]

use super::tests::*;
use super::*;
use std::collections::HashMap;

// ---- carrier phase ----------------------------------------------------------

#[test]
fn float_ambiguity_is_the_code_minus_the_windup_in_cycles() {
    // With pr = RANGE + IONO and cp = (pr + 5*lambda1)/lambda1 the true integer
    // is 5, and
    //   float_amb = (cp*lambda - windup_m - pr)/lambda
    //             = 5 - windup_rad/(2 pi)      [windup_m/lambda = rad/(2 pi)]
    // The satellite sits on the +Y axis here because the +X geometry makes the
    // receiver dipole projection identically zero (windup == 0), which would
    // leave the windup term untested.  The seeding algebra itself is
    // geometry independent.
    let sat_pos = Vector3::new(0.0, SAT_X, 0.0);
    let mut g = Graph::new(&[0]);
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());
    let mut o = obs_at(4, 1.0, sat_pos);
    o.cp_l1 = Some(cp_cycles(o.pr_l1, LAMBDA1, 5.0));
    run(&mut g, &[o], 0, None, 1.0, &mut slips, &mut trackers);

    let w = reference_windup_rad(&sat_pos);
    assert!(w.abs() > 0.01, "windup is ~0 ({w}); this assertion would be vacuous");
    let ids = g.ambiguities();
    assert_eq!(ids.len(), 1);
    let amb = g.solver.graph.variables[&ids[0].0].value[0];
    let want = 5.0 - w / std::f64::consts::TAU;
    // Tolerance 1e-6 cycles (1.9e-7 m).  `cp` is ~1.06e8 cycles, so re-forming
    // cp*lambda in the implementation carries ulp(2.0e7 m) = 4e-9 m of round
    // off, i.e. ~2e-8 cycles once divided by lambda.  Any real defect (wrong
    // sign, wrong wavelength, missing /lambda) moves this by 0.2 cycles or
    // more, so the bound is ~4 orders of magnitude tighter than needed.
    assert!((amb - want).abs() < 1e-6, "float ambiguity = {amb}, want {want}");
}

#[test]
fn l1_carrier_residual_equals_twice_the_slant_ionosphere_delay() {
    // Substituting the builder's own float_amb into the phase model cancels
    // every term:
    //   r = (cp*lambda - w) - (RANGE - IONO + lambda*float_amb)
    //     = (cp*lambda - w) - RANGE + IONO + cp*lambda - w - pr
    //     = 2*IONO - (pr - RANGE) = 2*IONO          [pr = RANGE + IONO]
    // so the residual is exactly 2 x 10 m = 20 m regardless of geometry,
    // wavelength or windup.  The phase model advances the carrier free of the
    // ionosphere (-I) while the float ambiguity is seeded from the code
    // (delayed by +I), so the two differ by 2I.
    let mut g = Graph::new(&[0]);
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());
    let mut o = obs(5, 1.0);
    o.cp_l1 = Some(cp_cycles(o.pr_l1, LAMBDA1, 5.0));
    run(&mut g, &[o], 0, None, 1.0, &mut slips, &mut trackers);

    let (_, cp) = factor_indices(&g);
    assert_eq!(cp.len(), 1, "one carrier factor expected");
    let v = g.values();
    let r = g.solver.graph.factors[cp[0]].residual(&v)[0];
    assert!((r - 2.0 * IONO).abs() < 1e-6, "carrier residual = {r}, want {}", 2.0 * IONO);
}

#[test]
fn l2_carrier_residual_uses_the_squared_frequency_ratio() {
    // The same cancellation holds for L2, but the phase is seeded from pr_l2
    // and the model scales the slant delay by gamma, so the residual is
    //   2 * gamma * IONO = 2 * 5929/3600 * 10 = 32.938 888 8... m
    // while L1 stays at 2 * IONO = 20 m.  This needs the +X geometry, where
    // the geometric range is exactly RANGE, so that pr - range vanishes.
    let mut g = Graph::new(&[0]);
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());
    let mut o = obs(6, 1.0);
    let pr_l2 = RANGE + GAMMA * IONO;
    o.pr_l2 = Some(pr_l2);
    o.cp_l1 = Some(cp_cycles(o.pr_l1, LAMBDA1, 5.0));
    o.cp_l2 = Some(cp_cycles(pr_l2, LAMBDA2, 7.0));
    run(&mut g, &[o], 0, None, 1.0, &mut slips, &mut trackers);

    let (_, cp) = factor_indices(&g);
    assert_eq!(cp.len(), 2, "L1 and L2 carrier factors expected");
    let v = g.values();
    let mut residuals: Vec<f64> = cp.iter().map(|&i| g.solver.graph.factors[i].residual(&v)[0]).collect();
    residuals.sort_by(f64::total_cmp);
    assert!(
        (residuals[0] - 2.0 * IONO).abs() < 1e-6,
        "L1 carrier residual = {}, want {}",
        residuals[0],
        2.0 * IONO
    );
    assert!(
        (residuals[1] - 2.0 * GAMMA * IONO).abs() < 1e-6,
        "L2 carrier residual = {}, want {}",
        residuals[1],
        2.0 * GAMMA * IONO
    );
}

#[test]
fn l2_ambiguity_is_seeded_with_its_own_wavelength() {
    // Each frequency gets its own ambiguity state, seeded from its own code
    // range through its own wavelength:
    //   float_amb1 = 5 - w/(2 pi)      (lambda1, integer 5)
    //   float_amb2 = 7 - w/(2 pi)      (lambda2, integer 7)
    // Seeding L2 with lambda1 would give 7 * lambda2/lambda1 = 8.988 instead.
    let sat_pos = Vector3::new(0.0, SAT_X, 0.0);
    let mut g = Graph::new(&[0]);
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());
    let mut o = obs_at(6, 1.0, sat_pos);
    let pr_l2 = RANGE + GAMMA * IONO;
    o.pr_l2 = Some(pr_l2);
    o.cp_l1 = Some(cp_cycles(o.pr_l1, LAMBDA1, 5.0));
    o.cp_l2 = Some(cp_cycles(pr_l2, LAMBDA2, 7.0));
    run(&mut g, &[o], 0, None, 1.0, &mut slips, &mut trackers);

    let w = reference_windup_rad(&sat_pos);
    let ids = g.ambiguities();
    assert_eq!(ids.len(), 2, "L1 and L2 must get separate ambiguity states");
    assert_ne!(ids[0].0, ids[1].0);
    assert_eq!((ids[0].1, ids[1].1), (1000, 2000), "frequencies 1 and 2, arc 0");
    let values: Vec<f64> = ids.iter().map(|(id, _)| g.solver.graph.variables[id].value[0]).collect();
    assert!((values[0] - (5.0 - w / std::f64::consts::TAU)).abs() < 1e-6, "L1 amb = {}", values[0]);
    assert!((values[1] - (7.0 - w / std::f64::consts::TAU)).abs() < 1e-6, "L2 amb = {}", values[1]);
}

#[test]
fn phase_information_is_the_elevation_weighted_noise_model() {
    // cp_var = (0.003 m / sin(elevation))^2, inverted by the factor:
    //   30 deg: sin = 0.5 -> 0.006^2 = 3.6e-5 m^2 -> W = 27 777.78
    //    0 deg: sin = 0 clamped to 0.1 -> 0.03^2 = 9e-4 m^2 -> W = 1 111.11
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());
    let mut g = Graph::new(&[0]);
    let mut o = obs(7, 30.0_f64.to_radians());
    o.cp_l1 = Some(cp_cycles(o.pr_l1, LAMBDA1, 5.0));
    run(&mut g, &[o], 0, None, 1.0, &mut slips, &mut trackers);
    let (_, cp) = factor_indices(&g);
    let w30 = g.solver.graph.factors[cp[0]].information()[(0, 0)];
    assert!((w30 - 1.0 / 3.6e-5).abs() < 1e-6, "W(30 deg) = {w30}, want {}", 1.0 / 3.6e-5);

    let mut g = Graph::new(&[0]);
    let mut o = obs(7, 0.0);
    o.cp_l1 = Some(cp_cycles(o.pr_l1, LAMBDA1, 5.0));
    run(&mut g, &[o], 0, None, 1.0, &mut slips, &mut trackers);
    let (_, cp) = factor_indices(&g);
    let w0 = g.solver.graph.factors[cp[0]].information()[(0, 0)];
    assert!((w0 - 1.0 / 9e-4).abs() < 1e-6, "W(0 deg) = {w0}, want {}", 1.0 / 9e-4);
}

// ---- cycle slips, arcs, and the ionosphere random walk ----------------------

#[test]
fn cycle_slip_starts_a_new_arc_with_a_new_ambiguity() {
    let mut g = Graph::new(&[0, 1]);
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());

    // Epoch 0: LLI bit 0 clear -> no slip, ambiguity in arc 0.
    let mut o0 = obs(8, 1.0);
    o0.cp_l1 = Some(cp_cycles(o0.pr_l1, LAMBDA1, 5.0));
    o0.cp_l1_lli = Some(0);
    run(&mut g, &[o0], 0, None, 1.0, &mut slips, &mut trackers);
    assert_eq!(slips.get(&(0, 8)).copied().unwrap_or(0), 0, "LLI = 0 is not a slip");

    // Epoch 1: LLI bit 0 set (half-cycle ambiguity flag) -> arc 1.
    let mut o1 = obs(8, 1.0);
    o1.cp_l1 = Some(cp_cycles(o1.pr_l1, LAMBDA1, 5.0));
    o1.cp_l1_lli = Some(1);
    run(&mut g, &[o1], 1, None, 1.0, &mut slips, &mut trackers);
    assert_eq!(slips[&(0, 8)], 1, "one cycle slip must be counted");

    let ids = g.ambiguities();
    assert_eq!(ids.len(), 2, "a slip must open a new ambiguity arc");
    assert_eq!(ids[0].1 % 1000, 0, "first ambiguity is arc 0");
    assert_eq!(ids[1].1 % 1000, 1, "post-slip ambiguity is arc 1");
    assert_eq!(ids[0].1 / 1000, 1, "both ambiguities are L1");
}

#[test]
fn ambiguity_is_never_re_initialised_once_seeded() {
    let mut g = Graph::new(&[0, 1]);
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());
    let mut o0 = obs(9, 1.0);
    o0.cp_l1 = Some(cp_cycles(o0.pr_l1, LAMBDA1, 5.0));
    run(&mut g, &[o0], 0, None, 1.0, &mut slips, &mut trackers);
    let seeded = g.ambiguities()[0].0;
    let after_epoch_0 = g.solver.graph.variables[&seeded].value[0];

    let mut o1 = obs(9, 1.0);
    o1.cp_l1 = Some(cp_cycles(o1.pr_l1, LAMBDA1, 5.0));
    run(&mut g, &[o1], 1, None, 1.0, &mut slips, &mut trackers);

    assert_eq!(g.ambiguities().len(), 1, "the same arc must reuse its ambiguity state");
    let after_epoch_1 = g.solver.graph.variables[&seeded].value[0];
    assert!(
        (after_epoch_1 - after_epoch_0).abs() < 1e-15,
        "ambiguity was re-seeded: {after_epoch_0} -> {after_epoch_1}"
    );
    // One scalar prior per arc, not one per epoch.
    let priors = g.solver.graph.factors.iter().filter(|f| f.variables() == [seeded]).count();
    assert_eq!(priors, 1, "the float prior must not be re-added every epoch");
}

#[test]
fn iono_random_walk_links_the_previous_epoch_slant_state() {
    // variance = max(4e-4 * dt, 1e-5): with dt = 1.0 that is 4e-4 m^2, so the
    // information is 2 500.  With dt = 0.0 the floor inside the factor bites:
    // 4e-4 * 0.1 = 4e-5 m^2, information 25 000.
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());
    let mut g = Graph::new(&[10]);
    let prev = g.solver.graph.add_variable(VariableKind::IonosphereSlant {
        epoch: 9,
        constellation_id: 0,
        satellite: 10,
    });
    run(&mut g, &[obs(10, 1.0)], 10, Some(9), 1.0, &mut slips, &mut trackers);

    let cur = g.iono_var(10).unwrap();
    let rw: Vec<_> = g.solver.graph.factors.iter().filter(|f| f.variables() == [prev, cur]).collect();
    assert_eq!(rw.len(), 1, "one random walk factor expected");
    assert!((rw[0].information()[(0, 0)] - 2500.0).abs() < 1e-9, "W = {}", rw[0].information()[(0, 0)]);

    let mut g2 = Graph::new(&[10]);
    let prev2 = g2.solver.graph.add_variable(VariableKind::IonosphereSlant {
        epoch: 9,
        constellation_id: 0,
        satellite: 10,
    });
    run(&mut g2, &[obs(10, 1.0)], 10, Some(9), 0.0, &mut slips, &mut trackers);
    let cur2 = g2.iono_var(10).unwrap();
    let rw2: Vec<_> = g2.solver.graph.factors.iter().filter(|f| f.variables() == [prev2, cur2]).collect();
    assert_eq!(rw2.len(), 1);
    assert!((rw2[0].information()[(0, 0)] - 25_000.0).abs() < 1e-9, "W = {}", rw2[0].information()[(0, 0)]);
}

#[test]
fn no_random_walk_without_a_previous_epoch_slant_state() {
    let mut g = Graph::new(&[10]);
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());
    run(&mut g, &[obs(11, 1.0)], 10, Some(9), 1.0, &mut slips, &mut trackers);
    assert_eq!(g.solver.graph.n_factors(), 1, "only the pseudorange factor");
    assert_eq!(g.n_slant(), 1, "the previous-epoch state must not be conjured up");
}

#[test]
fn every_built_factor_references_a_variable_that_exists() {
    // Structural invariant: the builder may not leave a dangling edge, which
    // would panic inside a factor's residual().
    let mut g = Graph::new(&[10]);
    g.solver.graph.add_variable(VariableKind::IonosphereSlant {
        epoch: 9,
        constellation_id: 0,
        satellite: 12,
    });
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());
    let mut o = obs(12, 1.0);
    let pr_l2 = RANGE + GAMMA * IONO;
    o.pr_l2 = Some(pr_l2);
    o.cp_l1 = Some(cp_cycles(o.pr_l1, LAMBDA1, 5.0));
    o.cp_l2 = Some(cp_cycles(pr_l2, LAMBDA2, 7.0));
    run(&mut g, &[o], 10, Some(9), 1.0, &mut slips, &mut trackers);

    let values = g.values();
    assert_eq!(
        g.solver.graph.n_factors(),
        7,
        "2 pseudorange + 2 carrier + 2 ambiguity priors + 1 iono random walk"
    );
    for f in &g.solver.graph.factors {
        for v in f.variables() {
            assert!(values.index_of(*v).is_some(), "factor references missing variable {v:?}");
        }
        assert!(f.residual(&values).iter().all(|r| r.is_finite()), "non-finite residual in {f:?}");
    }
}

#[test]
fn l2_carrier_phase_without_l2_code_falls_back_to_the_l1_range() {
    // A satellite tracked on the L2 carrier but not the L2 code is a normal
    // situation (L2 is weaker and often lost first).  The ambiguity is then
    // seeded from the L1 code range instead, i.e.
    //   float_amb2 = (cp2*lambda2 - windup_m2 - pr_l1) / lambda2.
    // With cp2 chosen for integer 7 against pr_l2 = RANGE + GAMMA*IONO the
    // fallback seeds from the smaller L1 range instead, which RAISES the
    // seed: 7 + (gamma - 1) * IONO / lambda2 - windup_m2 / lambda2.
    let sat_pos = Vector3::new(0.0, SAT_X, 0.0);
    let mut g = Graph::new(&[0]);
    let (mut slips, mut trackers) = (HashMap::new(), HashMap::new());
    let mut o = obs_at(13, 1.0, sat_pos);
    o.cp_l1 = Some(cp_cycles(o.pr_l1, LAMBDA1, 5.0));
    o.cp_l2 = Some(cp_cycles(RANGE + GAMMA * IONO, LAMBDA2, 7.0));
    run(&mut g, &[o], 0, None, 1.0, &mut slips, &mut trackers);

    let (pr, cp) = factor_indices(&g);
    assert_eq!(pr.len(), 1, "there is no L2 code, so only the L1 pseudorange factor");
    assert_eq!(cp.len(), 2, "both carrier phases are still usable");
    let w = reference_windup_rad(&sat_pos);
    let ids = g.ambiguities();
    assert_eq!(ids.len(), 2);
    let l1 = g.solver.graph.variables[&ids[0].0].value[0];
    let l2 = g.solver.graph.variables[&ids[1].0].value[0];
    assert!((l1 - (5.0 - w / std::f64::consts::TAU)).abs() < 1e-6, "L1 amb = {l1}");
    // (GAMMA - 1) * IONO / lambda2 = 6.4694444 / 0.24420532 = 26.491 cycles.
    let shift = (GAMMA - 1.0) * IONO / LAMBDA2;
    let want = 7.0 + shift - w / std::f64::consts::TAU;
    assert!((l2 - want).abs() < 1e-4, "L2 amb = {l2}, want {want} (fallback shift {shift} cycles)");
}
