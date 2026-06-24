# Exploration Report: Bug 22 and Bug 3 Investigation

This report documents findings on Bug 22 (Saastamoinen Zenith Delay Pressure Scaling) and Bug 3 (L2C Phase Shift Bias) in the `gneiss` GNSS engine.

---

## Bug 22: Saastamoinen Zenith Delay and Inconsistent Scaling

### 1. Issue Description
- **Lack of Altitude Scaling in `tropo_nmf`**:
  `AtmosphereModel::tropo_nmf` is used in SPP estimators (e.g. `spp.rs` and `spp_tight.rs`). It is called with `TropoParams::default()` (yielding standard sea level pressure of 1013.25 hPa and temperature of 288.15 K). However, inside `tropo_nmf`, pressure (`params.press_hpa`) and temperature (`params.temp_k`) are **not scaled by the receiver altitude/height** (`pos_llh.z`). This causes the SPP algorithm to compute the tropospheric dry and wet delays as if the receiver is always at sea level, leading to severe delay overestimation at higher geodetic heights.
- **Inconsistent Barometric Constants**:
  `compute_tropo_dry` in `crates/gneiss-rtk/src/engine/ppp_math.rs` implements height scaling using barometric constants ($0.0000226$ and $5.225$) that are inconsistent with the more precise physically-derived constants ($0.000022557$ and $5.2568$) used in `AtmosphereModel::tropo_rtklib_saastamoinen` in `crates/gneiss-core/src/atmosphere.rs`.

### 2. File Paths and Line Numbers
- **File**: `crates/gneiss-core/src/atmosphere.rs`
- **Line Numbers**: 675–683
- **Code Snippet**:
  ```rust
  // Zenith dry and wet delays (simplified Saastamoinen)
  let z_dry = 0.0022768 * params.press_hpa
      / (1.0 - 0.00266 * libm::cos(2.0 * pos_llh.x) - 0.00028 * hgt / 1000.0);

  let e = 6.108
      * libm::exp((17.15 * params.temp_k - 4684.0) / (params.temp_k - 38.45))
      * params.hum_rel;
  let z_wet = 0.002277 * (1255.0 / params.temp_k + 0.05) * e;
  ```

- **File**: `crates/gneiss-rtk/src/engine/ppp_math.rs`
- **Line Numbers**: 271–272
- **Code Snippet**:
  ```rust
  let press_scale = libm::pow(1.0 - 0.0000226 * alt_m, 5.225);
  let p_z = tp.press_hpa * press_scale;
  ```

---

## Bug 3: L2C Phase Shift Bias and Unit Test Fallbacks

### 1. Issue Description
- **Missing L2C Phase Shift Fallback**:
  For GPS, civilian L2C signals have a $+0.25$ cycle phase offset relative to the military L2P(Y) (or codeless `L2W`) signals. When the bias file lacks explicit L2C biases, fallback logic maps observations to the `L2W` bias record. However, since the L2C phase observations are shifted by $+0.25$ cycles, we must apply a correction of $-0.25$ cycles (subtracting $0.25$) when using the `L2W` bias. This shift is completely missing in `apply_osb_corrections` in `ppp_math.rs`.
- **Inaccurate Test Assertions**:
  The unit tests `test_apply_osb_corrections_cp2_fb` (in `ppp_math.rs`) and `test_apply_osb_shift` (in `ppp.rs`) pass only because they assert raw un-shifted values.

### 2. File Paths and Line Numbers
- **File**: `crates/gneiss-rtk/src/engine/ppp_math.rs`
- **Line Numbers**: 69–74
- **Code Snippet**:
  ```rust
  if let Some(v) = out.cp1.as_mut() {
      *v -= out.osb_cp1 / (LIGHT_SPEED / f1);
  }
  if let Some(v) = out.cp2.as_mut() {
      *v -= out.osb_cp2 / (LIGHT_SPEED / f2);
  }
  ```

- **File**: `crates/gneiss-rtk/src/engine/ppp_math.rs` (Unit Test `test_apply_osb_corrections_cp2_fb`)
- **Line Numbers**: 640–644 (specifically line 641)
- **Code Snippet**:
  ```rust
  if let Some(cp2) = out.cp2 {
      assert_eq!(cp2, 300.0, "failed for code L2C");
  }
  ```
  *(Also see lines 590-617 check loop, where `L2L`, `L2S`, `L2X` are asserted to `297.6` rather than `297.6 - 0.25 = 297.35`)*

- **File**: `crates/gneiss-rtk/src/engine/ppp.rs` (Unit Test `test_apply_osb_shift`)
- **Line Numbers**: 811
- **Code Snippet**:
  ```rust
  assert!((res.cp2.unwrap() - (40.0 - (4.0 * 1e-9 * LIGHT_SPEED) / wl2)).abs() < 1e-6);
  ```
  *(This should be `assert!((res.cp2.unwrap() - (40.0 - (4.0 * 1e-9 * LIGHT_SPEED) / wl2 - 0.25)).abs() < 1e-6)`)*
