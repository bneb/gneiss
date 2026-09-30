# TEST_READY — Frame Safety & Epoch Alignment Refactoring

**Date**: 2026-09-25  
**Author**: Test Writer E2E (`teamwork_preview_test_writer_e2e`)  
**Status**: COMPLETE (118/118 tests passing, 0 clippy warnings)  
**Target Package**: `gneiss-tests` (Integration test suite: `test_frame_safety_e2e`)  
**Specification References**:  
- `ORIGINAL_REQUEST.md` (2026-09-25T21:01:47Z)  
- `PROJECT.md` (Features 1–21, Milestones 1–5)  
- `AGENTS.md` (Code Quality Standards & Relational Coupling Guidelines)  

---

## 1. Executive Summary

A comprehensive, requirement-driven, opaque-box 4-tier End-to-End (E2E) integration test suite has been designed, implemented, and verified for the Gneiss Frame Safety & Epoch Alignment Refactoring.

All 118 test cases compile with zero warnings and pass with a 100% success rate:
- **Total Tests**: 118
- **Passed**: 118
- **Failed**: 0
- **Ignored**: 0
- **Clippy Warnings**: 0 (under `cargo clippy -p gneiss-tests --all-targets -- -D warnings`)
- **Execution Time**: ~0.01 seconds

---

## 2. Test Architecture & File Manifest

The test suite is structured modularly under `tests/tests/test_frame_safety_e2e/` with an explicit integration root test harness `tests/tests/test_frame_safety_e2e.rs`.

| File Path | Description | LOC | Test Count | Max Fn Length | Nesting |
|---|---|:---:|:---:|:---:|:---:|
| `tests/tests/test_frame_safety_e2e.rs` | Test root harness module aggregator | 29 | 0 | N/A | 0 |
| `tests/tests/test_frame_safety_e2e/types.rs` | Opaque-box typestates, scale markers, contracts | 326 | 0 | 18 LOC | < 3 |
| `tests/tests/test_frame_safety_e2e/common.rs` | Mathematical oracles, geodetic helpers, geometry | 163 | 0 | 22 LOC | < 3 |
| `tests/tests/test_frame_safety_e2e/tier1_spatial.rs` | Tier 1: Spatial & datum typestates (F1–F5) | 269 | 25 | 18 LOC | < 3 |
| `tests/tests/test_frame_safety_e2e/tier1_temporal.rs`| Tier 1: Temporal typestates & scales (F6–F10) | 232 | 25 | 16 LOC | < 3 |
| `tests/tests/test_frame_safety_e2e/tier1_geometry.rs`| Tier 1: Relational DD geometry & PCV (F11–F13) | 161 | 15 | 18 LOC | < 3 |
| `tests/tests/test_frame_safety_e2e/tier1_estimator.rs`| Tier 1: Estimator state & lever arms (F14–F17)| 204 | 20 | 20 LOC | < 3 |
| `tests/tests/test_frame_safety_e2e/tier2_boundaries.rs`| Tier 2: Boundary, rollover & singularity tests | 218 | 20 | 20 LOC | < 3 |
| `tests/tests/test_frame_safety_e2e/tier3_pairwise.rs` | Tier 3: Cross-feature pairwise interactions | 114 | 8 | 24 LOC | < 3 |
| `tests/tests/test_frame_safety_e2e/tier4_scenarios.rs`| Tier 4: Real-world operational mission scenarios| 76 | 5 | 22 LOC | < 3 |

*All files strictly comply with AGENTS.md standards: strictly < 500 lines per file, strictly <= 32 lines per function, nesting depth < 3 levels.*

---

## 3. Tier Coverage Breakdown

### Tier 1: Feature Functionality (85 Tests)
- **Features 1–5: Spatial Primitives & Reference Frames (25 Tests)**:
  - `test_f01_point3_affine_translation_consistency`: Point translation $p_2 = p_1 + v \implies p_2 - p_1 = v$.
  - `test_f01_spatial_vector_addition_commutativity`: Vector commutativity in local NED frame.
  - `test_f01_spatial_velocity_norm_and_components`: Vector velocity Euclidean norm invariance.
  - `test_f01_ned_covariance_std_extraction`: Extraction of North, East, Down standard deviations.
  - `test_f01_affine_subtraction_roundtrip`: Reversible affine point-vector displacements.
  - `test_f02_reference_frame_names_and_distinctness`: Reference frame tags (ITRF2014, ITRF2020, WGS84, NAD83, JGD2011, PZ-90).
  - `test_f02_pz90_helmert_translation_parameters`: Translation shifts $T_x = 3.0$ mm, $T_y = -1.0$ mm, $T_z = 0.0$ mm.
  - `test_f02_pz90_helmert_rotation_parameters`: Small-angle rotation verification ($R_x = 0.019$, $R_y = -0.042$, $R_z = 0.002$ mas).
  - `test_f02_pz90_coordinate_shift_magnitude`: Bounded PZ-90 to ITRF coordinate shift ($2\text{--}5$ mm).
  - `test_f02_itrf2020_to_itrf2014_epoch_continuity`: Metric shift continuity across ITRF realizations.
  - `test_f03_antenna_lever_arm_initialization`: Body FRD lever arm component tracking.
  - `test_f03_attitude_rotation_body_to_ecef`: Rotation via $C_b^e$ direction cosine matrix.
  - `test_f03_attitude_rotation_90deg_yaw`: +90° yaw maps vehicle forward (+X) to East (+Y in NED).
  - `test_f03_attitude_inverse_roundtrip`: Exact inverse rotation $(C_b^e)^{-1} = C_e^b$.
  - `test_f03_attitude_covariance_rotation_trace_invariance`: Matrix trace preservation under orthogonal rotation.
  - `test_f04_explicit_coords_access`: Explicit `.coords()` accessor without leaky deref.
  - `test_f04_owned_coords_extraction`: Clean `.into_coords()` value extraction.
  - `test_f04_spatial_vector_norm_invariant_under_pure_rotation`: Vector $L_2$ norm invariance under $SO(3)$.
  - `test_f04_null_translation_preserves_position`: Identity translation $p + 0 = p$.
  - `test_f04_vector_addition_associativity`: Vector addition associativity $(v_1 + v_2) + v_3 = v_1 + (v_2 + v_3)$.
  - `test_f05_local_tangent_plane_origin_projection`: Origin projects to $(0, 0, 0)$ in ENU/NED.
  - `test_f05_local_tangent_plane_roundtrip`: Reversible ECEF $\leftrightarrow$ ENU projection.
  - `test_f05_local_tangent_plane_metric_distance_preservation`: Metric Euclidean distance invariance for regional points.
  - `test_f05_local_tangent_plane_covariance_projection_positive_definite`: Positive definiteness preserved in ENU covariance.
  - `test_f05_local_tangent_plane_vertical_up_matches_zenith`: Equator local Up aligns with ECEF $+X$.

- **Features 6–10: Temporal Primitives & Epoch Alignment (25 Tests)**:
  - `test_f06_timescale_name_identifiers`: GPST, BDT, GST, GLONASST, UTC identification.
  - `test_f06_epoch_creation_from_week_and_nanos`: Microsecond/nanosecond exact representation.
  - `test_f06_epoch_normalization_overflow_carries_week`: `tow_nanos >= 604_800s` rolls into next week.
  - `test_f06_timedelta_unit_constructors`: Exact scaling across nanos, millis, and seconds.
  - `test_f06_epoch_monotonic_ordering`: Strict lexicographic ordering by `(week, tow_nanos)`.
  - `test_f07_subtraction_within_same_week`: Integer nanosecond difference within week.
  - `test_f07_subtraction_across_week_boundary`: Exact difference across Saturday/Sunday boundary.
  - `test_f07_addition_forward_time_delta`: Forward epoch advancement.
  - `test_f07_addition_negative_time_delta_across_week`: Backward epoch rollback into preceding week.
  - `test_f07_exact_nanosecond_precision_preservation`: 0 float precision loss across repeated additions.
  - `test_f08_bdt_to_gpst_adds_exact_14_seconds`: BeiDou BDT $\to$ GPST adds exact $14.0\text{ s}$ ($14 \times 10^9$ ns).
  - `test_f08_gpst_to_bdt_subtracts_exact_14_seconds`: GPST $\to$ BDT subtracts exact $14.0\text{ s}$.
  - `test_f08_bdt_gpst_roundtrip_nanosecond_invariance`: Exact roundtrip preservation.
  - `test_f08_bdt_to_gpst_across_week_boundary`: Rollover-safe scale conversion.
  - `test_f08_beidou_clock_correction_alignment`: Eliminates 14-second clock bias artifact.
  - `test_f09_utc_to_gpst_explicit_18_seconds`: Explicit 18-second leap parameter requirement.
  - `test_f09_gpst_to_utc_explicit_18_seconds`: GPST $\to$ UTC deduction.
  - `test_f09_utc_gpst_roundtrip`: Reversible UTC/GPST alignment.
  - `test_f09_future_leap_second_19s_parameter`: Scalable future leap second support.
  - `test_f09_glonass_epoch_relationship`: GLONASS UTC(SU)+3h offset and leap second coupling.
  - `test_f10_continuous_ms_monotonicity`: Linear continuous millisecond counter across week rollovers.
  - `test_f10_epoch_key_equality_and_hashing`: Collision-free hashing and equivalence.
  - `test_f10_is_within_tolerance_positive`: Tolerance-based epoch matching within threshold.
  - `test_f10_is_within_tolerance_negative`: Rejection outside tolerance threshold.
  - `test_f10_sub_millisecond_clock_jitter_handling`: Clock jitter tolerance between base and rover streams.

- **Features 11–13: Relational DD Geometry & Antenna Coupling (15 Tests)**:
  - `test_f11_dd_geometry_construction_and_unit_los_vectors`: Unit line-of-sight extraction.
  - `test_f11_dd_geometry_independent_elevations_and_azimuths`: Independent angles for rover and base.
  - `test_f11_dd_geometry_zero_baseline_vanishes_dd_range`: Geometric DD vanishes on zero baseline.
  - `test_f11_dd_geometry_satellite_range_magnitude`: Valid MEO orbital range ($20{,}000\text{--}30{,}000$ km).
  - `test_f11_dd_geometry_delta_u_rov_for_eskf_jacobian`: Single-difference delta unit vector for ESKF.
  - `test_f12_llh_vs_ecef_units_and_magnitude_distinction`: Catches LLH-passed-as-ECEF bug.
  - `test_f12_proper_ecef_origin_line_of_sight_norm`: True line-of-sight range from geocentric origin.
  - `test_f12_elevation_angle_valid_bounds`: Elevation angle bounded within $[-\pi/2, \pi/2]$.
  - `test_f12_zenith_angle_complement`: Zenith and elevation sum to exactly 90°.
  - `test_f12_receiver_pcv_correction_boundedness`: Millimetric PCV phase corrections ($< 50$ mm).
  - `test_f13_regional_baseline_elevation_divergence`: Elevation angle discrepancy on 50 km baseline ($> 0.05^\circ$).
  - `test_f13_base_elevation_independent_pcv_lookup`: PCV lookup uses true base elevation.
  - `test_f13_zenith_satellite_minimal_dd_pcv_difference`: High-elevation satellite pair stability.
  - `test_f13_double_difference_symmetry_inversion`: Satellite permutation inverts DD range sign.
  - `test_f13_base_and_rover_range_difference_triangle_inequality`: Spatial range bound.

- **Features 14–17: Estimator State & Pipelines (20 Tests)**:
  - `test_f14_eskf_15state_partitioning_and_dimension`: 15-state ESKF vector structure.
  - `test_f14_attitude_quaternion_normalization`: Error quaternion unit norm invariance.
  - `test_f14_body_biases_typed_in_body_frd`: Sensor biases strictly in Body FRD.
  - `test_f14_initial_covariance_symmetry_and_positive_diagonal`: Error covariance symmetry and definiteness.
  - `test_f14_state_propagation_position_step`: Kinematic state transition.
  - `test_f15_gnss_pos_innovation_with_rotated_lever_arm`: Innovation with rotated body lever arm.
  - `test_f15_attitude_error_jacobian_skew_symmetric`: Attitude Jacobian $H_\theta = -[\mathbf{l}^e \times]$.
  - `test_f15_zero_lever_arm_attitude_jacobian_vanishes`: Vanishing attitude coupling when $\mathbf{l}^b = 0$.
  - `test_f15_doppler_velocity_lever_arm_angular_rate_coupling`: Angular rate cross-product $\boldsymbol{\omega}_b \times \mathbf{l}^b$.
  - `test_f15_doppler_gyro_bias_jacobian`: Gyro bias Jacobian $\frac{\partial v}{\partial b_g} = -C_b^e [\mathbf{l}^b \times]$.
  - `test_f16_postprocess_options_typed_base_position`: Typed base position option.
  - `test_f16_postprocess_options_typed_lever_arm`: Typed lever arm option.
  - `test_f16_swfg_pose_datum_consistency`: Factor graph initial pose datum consistency.
  - `test_f16_swfg_dd_factor_lever_arm_projection`: Factor graph lever arm projection.
  - `test_f16_swfg_epoch_delta_nanosecond_precision`: Nanosecond exactness for graph node spacing.
  - `test_f17_odaiba_base_arp_reference_position`: Baseline Tsukuba ARP reference position.
  - `test_f17_odaiba_antenna_height_offset_projection`: Antenna height NED to ECEF projection.
  - `test_f17_odaiba_ned_covariance_to_ecef_projection`: Local variance to ECEF covariance projection.
  - `test_f17_odaiba_zero_lever_arm_invariant`: Odaiba zero lever arm baseline.
  - `test_f17_odaiba_innovation_chi_square_consistency`: Chi-square innovation gating bound.

### Tier 2: Boundary & Corner Cases (20 Tests)
- Antimeridian crossing ($\pm 180^\circ$ longitude continuity).
- North and South polar singularities ($\pm 90^\circ$ latitude height extraction).
- Zero-length lever arm null projection.
- Extreme 10,000 km intercontinental chord distance.
- Midnight week rollover exact nanosecond boundary (`WEEK_NANOS - 1` to `0`).
- Sub-millisecond tolerance exact window edge (5.000 ms inside, 5.001 ms outside).
- Leap second transition instant (+1.0 s step).
- Negative subtraction order anti-symmetry ($t_1 - t_2 = -(t_2 - t_1)$).
- Zero duration delta identity ($t - t = 0$).
- Zenith satellite elevation limit ($\theta_{el} \to \pi/2$).
- Horizon satellite elevation limit ($\theta_{el} \to 0$).
- Co-located base/rover zero double difference.
- Collinear satellite line-of-sight unit vector cancellation.
- Azimuth $[0, 2\pi)$ normalization.
- Zero error quaternion identity.
- $\pi$ ($180^\circ$) yaw rotation boundary.
- High angular rate centrifugal coupling ($10\text{ rad/s}$ turn).
- Near-singular covariance positive eigenvalue stability.
- Infinitesimal time step propagation ($\Delta t = 10^{-6}\text{ s}$).

### Tier 3: Cross-Feature Pairwise Interactions (8 Tests)
- **Spatial + Temporal**: Satellite orbit propagation with typed epoch & datum.
- **Spatial + Geometry**: Double-difference baseline displacement with typed coordinates.
- **Temporal + Estimator**: ESKF state propagation across Saturday/Sunday week rollover.
- **Lever Arm + Attitude + Velocity**: Vehicle turn with lever arm Doppler innovation.
- **Tangent Plane + Covariance**: ECEF to NED/ENU covariance projection.
- **Multi-Constellation**: GPS + BeiDou + UTC joint epoch alignment and scale conversion.
- **Lever Arm + Double-Difference**: DD carrier phase attitude coupling Jacobian $H_\theta = \Delta\mathbf{u}^T [\mathbf{l}^e \times]$.
- **Relational PCV + Double-Difference**: Zero-baseline elevation and azimuth equality.

### Tier 4: Real-World Application Scenarios (5 Scenarios)
- **Scenario 1: Tokyo Odaiba Kinematic INS with Offset Antenna Lever Arm**: High-speed highway vehicle run with lever arm $[0.25, 0.10, -0.85]$ m under roll/pitch vehicle dynamics.
- **Scenario 2: Multi-Station CORS Baseline Across Distinct Datums**: Tsukuba (JGD2011) to Tokyo (ITRF2014) regional CORS baseline (~54 km).
- **Scenario 3: Saturday Midnight Week-Rollover Continuous Kinematic**: 50 Hz continuous kinematic trajectory traversing week rollover with strictly monotonic `EpochKey` generation.
- **Scenario 4: BeiDou B1I/B3I Clock-Bias-Free Baseline Solution**: Joint GPS + BDS processing with 14-second BDT scale conversion eliminating clock bias errors.
- **Scenario 5: High-Dynamic UAV Pitch/Roll with Doppler & DD Updates**: 45° bank turn at 30 m/s with 1.5 rad/s yaw rate and antenna lever arm perturbation.

---

## 4. Verification Commands & Outputs

```bash
# 1. Run the frame safety integration test suite
cargo test -p gneiss-tests --test test_frame_safety_e2e
# Output: test result: ok. 118 passed; 0 failed; 0 ignored; finished in 0.01s

# 2. Run clippy verification across all targets
cargo clippy -p gneiss-tests --all-targets -- -D warnings
# Output: Finished dev profile, 0 warnings

# 3. Check line counts of all test files (strictly < 500 LOC)
wc -l tests/tests/test_frame_safety_e2e.rs tests/tests/test_frame_safety_e2e/*.rs TEST_INFRA.md
# Output: All files <= 326 LOC
```

---

## 5. Conclusion & Readiness

The test infrastructure `TEST_INFRA.md` and complete test suite `test_frame_safety_e2e` are fully operational, verified, and ready. All opaque-box contracts and structural invariants specified in `ORIGINAL_REQUEST.md` and `PROJECT.md` are protected against regressions.
