//! Shared deterministic fixtures for sidereal submodule tests.
//! Test-only: never compiled into production builds.

use std::f64::consts::PI;

use nalgebra::Vector3;

use gneiss_core::time::GpsTime;

use super::phase::{sidereal_phase, Sample};
use super::SmoothedEpoch;

/// Deterministic LCG + Box-Muller so statistical tests never flake.
pub(crate) struct Lcg(pub(crate) u64);

impl Lcg {
    pub(crate) fn uniform(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 11) as f64) / (1u64 << 53) as f64
    }

    pub(crate) fn normal(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-12);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }
}

/// 30 s epochs with a sinusoid of amplitude `amp` locked to sidereal phase
/// plus Gaussian noise of standard deviation `sigma`.
pub(crate) fn synthetic_samples(
    n: usize,
    tow0: f64,
    amp: f64,
    sigma: f64,
    seed: u64,
) -> Vec<Sample> {
    let mut rng = Lcg(seed);
    (0..n)
        .map(|i| {
            let tow = tow0 + i as f64 * 30.0;
            let v = amp * (2.0 * PI * sidereal_phase(tow)).sin() + sigma * rng.normal();
            Sample { tow_s: tow, value: v }
        })
        .collect()
}

/// Truth at lat=0/lon=0 so +ECEF X = up, +Y = east, +Z = north there:
/// axis semantics are readable directly in the assertions.
pub(crate) const TRUTH: Vector3<f64> = Vector3::new(6_378_137.0, 0.0, 0.0);

pub(crate) fn truth_at(_tow: u32) -> Option<Vector3<f64>> {
    Some(TRUTH)
}

pub(crate) fn mk_epoch(tow: f64, ecef: Vector3<f64>) -> SmoothedEpoch {
    SmoothedEpoch {
        time: GpsTime::new(2000, tow),
        position_ecef: ecef,
        velocity_ecef: None,
        attitude: None,
        cov_position: nalgebra::Matrix3::zeros(),
        std_east: 0.0,
        std_north: 0.0,
        std_up: 0.0,
        separation_3d: 0.0,
        quality: 1,
        n_satellites: 8,
    }
}
