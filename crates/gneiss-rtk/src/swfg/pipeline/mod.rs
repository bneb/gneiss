//! Unified measurement pipeline — stacked corrections that transform raw
//! observations into factor graph factors.
//!
//! SPP, PPP, and RTK differ only in which correction passes are stacked.

pub mod aux_factors;
pub mod bias_factor;
pub mod dd_factors;
pub mod factors;
pub mod passes;

pub use aux_factors::{DualAntennaHeadingFactor, OdometerVelocityFactor};
pub use bias_factor::ImuBiasPriorFactor;
pub use dd_factors::{
    elevation_cp_variance, elevation_pr_variance, DdCarrierPhaseFactor, DdDopplerFactor,
    DdPseudorangeFactor, WidelaneConstraintFactor,
};
pub use factors::{
    build_carrier_phase_factor, build_pseudorange_factor, CarrierPhaseFactor,
    DopplerVelocityFactor, NhcFactor, PseudorangeFactor,
};
pub use passes::{
    BroadcastClockCorrection, CorrectionPass, CorrectedObservation, KlobucharIono,
    PreciseClockCorrection, RawObservation, ReceiverState, SaastamoinenTropo, SnrVarianceModel,
};

/// The measurement pipeline transforms raw observations into factors
/// by applying a stack of correction passes.
#[derive(Debug)]
pub struct MeasurementPipeline {
    passes: Vec<Box<dyn CorrectionPass>>,
}

impl MeasurementPipeline {
    pub fn new() -> Self {
        Self { passes: Vec::new() }
    }

    pub fn add_pass(&mut self, pass: Box<dyn CorrectionPass>) {
        self.passes.push(pass);
    }

    /// Set Klobuchar ionosphere parameters. Applied in the next `process` call.
    pub fn set_klobuchar(&mut self, alpha: [f64; 4], beta: [f64; 4]) {
        // Remove any existing KlobucharIono pass and replace with configured one.
        self.passes.retain(|p| p.name() != "klobuchar_iono");
        self.add_pass(Box::new(KlobucharIono { alpha, beta }));
    }

    /// Process raw observations through all correction passes.
    pub fn process(
        &self,
        raw: &[RawObservation],
        rx_state: &ReceiverState,
    ) -> Vec<CorrectedObservation> {
        raw.iter()
            .map(|r| {
                let mut work = r.clone();
                for pass in &self.passes {
                    pass.apply(&mut work, rx_state);
                }
                CorrectedObservation {
                    satellite: work.satellite,
                    constellation_id: work.constellation_id,
                    pr_l1: work.pr_l1,
                    pr_l2: work.pr_l2,
                    cp_l1: work.cp_l1,
                    cp_l1_lli: work.cp_l1_lli,
                    cp_l2: work.cp_l2,
                    doppler: work.doppler,
                    snr_dbhz: work.snr_dbhz,
                    sat_pos_ecef: work.sat_pos_ecef,
                    sat_clock_m: work.sat_clock_m,
                    f1: work.f1,
                    f2: work.f2,
                    freq_num: work.freq_num,
                    elevation_rad: work.elevation_rad,
                    tropo_dry_m: work.tropo_dry_m,
                    tropo_map_wet: work.tropo_map_wet,
                    iono_l1_m: work.iono_l1_m,
                    variance_m2: work.variance_m2,
                    cp_variance_m2: work.cp_variance_m2,
                }
            })
            .collect()
    }

    /// Build a PPP pipeline with precise products, tropo, iono, and variance.
    pub fn ppp_mode() -> Self {
        let mut pipeline = Self::new();
        pipeline.add_pass(Box::new(PreciseClockCorrection {
            clock_biases_m: Vec::new(),
            constellation_id: 0,
        }));
        pipeline.add_pass(Box::new(SaastamoinenTropo));
        pipeline.add_pass(Box::new(KlobucharIono::default()));
        pipeline.add_pass(Box::new(SnrVarianceModel::default()));
        pipeline
    }

    /// Build an RTK pipeline — for short baselines, iono cancels in DD
    /// so we skip the iono correction.
    pub fn rtk_mode() -> Self {
        let mut pipeline = Self::new();
        pipeline.add_pass(Box::new(SnrVarianceModel::default()));
        pipeline
    }

    /// Build an SPP pipeline — minimal corrections.
    pub fn spp_mode() -> Self {
        let mut pipeline = Self::new();
        pipeline.add_pass(Box::new(BroadcastClockCorrection));
        pipeline.add_pass(Box::new(SaastamoinenTropo));
        pipeline.add_pass(Box::new(KlobucharIono::default()));
        pipeline.add_pass(Box::new(SnrVarianceModel::default()));
        pipeline
    }
}

impl Default for MeasurementPipeline {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::swfg::variables::{VariableId, VariableKind, VariableNode, VariableValues};
    use nalgebra::Vector3;
    use std::collections::BTreeMap;

    fn make_test_vars(epoch: u32) -> (BTreeMap<VariableId, VariableNode>, VariableId, VariableId) {
        let mut vars = BTreeMap::new();
        let pose_id = VariableId::new(0);
        let zwd_id = VariableId::new(1);
        let mut node_pose = VariableNode::new(pose_id, VariableKind::Pose { epoch });
        node_pose.set_value(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let mut node_zwd = VariableNode::new(zwd_id, VariableKind::TropoZwd { epoch });
        node_zwd.set_value(&[0.1]);
        vars.insert(pose_id, node_pose);
        vars.insert(zwd_id, node_zwd);
        (vars, pose_id, zwd_id)
    }

    #[test]
    fn pseudorange_factor_residual_zero_at_truth() {
        let (vars, pose_id, zwd_id) = make_test_vars(0);
        let sat_pos = Vector3::new(100.0, 0.0, 0.0);
        let obs = CorrectedObservation {
            satellite: 1,
            constellation_id: 0,
            pr_l1: 100.0,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            cp_l1_lli: None,
            doppler: 0.0,
            snr_dbhz: 40.0,
            sat_pos_ecef: sat_pos,
            sat_clock_m: 0.0,
            f1: 1575.42e6,
            f2: 1227.60e6,
            freq_num: 0,
            elevation_rad: std::f64::consts::FRAC_PI_2,
            tropo_dry_m: 0.0,
            tropo_map_wet: 1.0,
            iono_l1_m: 0.0,
            variance_m2: 1.0,
            cp_variance_m2: 9e-6,
        };

        let factor = build_pseudorange_factor(&obs, 0, pose_id, None, Some(zwd_id), None);
        let values = VariableValues::build(&vars);
        let r = factor.residual(&values);
        assert!((r[0] + 0.1).abs() < 1e-6);
    }

    #[test]
    fn carrier_phase_factor_subtracts_windup_correctly() {
        let (vars, pose_id, zwd_id) = make_test_vars(0);
        let amb_id = VariableId::new(2);
        let mut vars_mut = vars.clone();
        let mut amb_node = VariableNode::new(
            amb_id,
            VariableKind::Ambiguity {
                satellite: 1,
                frequency: 1,
            },
        );
        amb_node.set_value(&[0.0]);
        vars_mut.insert(amb_id, amb_node);

        let lambda = gneiss_core::constants::SPEED_OF_LIGHT_M_S / 1575.42e6;
        let windup_m = 0.05;
        let cp_cycles = 50.0 / lambda;

        let obs = CorrectedObservation {
            satellite: 1,
            constellation_id: 0,
            pr_l1: 50.0,
            pr_l2: None,
            cp_l1: Some(cp_cycles),
            cp_l2: None,
            cp_l1_lli: None,
            doppler: 0.0,
            snr_dbhz: 45.0,
            sat_pos_ecef: Vector3::new(50.0, 0.0, 0.0),
            sat_clock_m: 0.0,
            f1: 1575.42e6,
            f2: 1227.60e6,
            freq_num: 0,
            elevation_rad: std::f64::consts::FRAC_PI_2,
            tropo_dry_m: 0.0,
            tropo_map_wet: 1.0,
            iono_l1_m: 0.0,
            variance_m2: 1.0,
            cp_variance_m2: 9e-6,
        };

        let factor = build_carrier_phase_factor(&obs, 0, pose_id, None, Some(zwd_id), amb_id, None, windup_m);
        let values = VariableValues::build(&vars_mut);
        let r = factor.residual(&values);
        assert!((r[0] - (-windup_m - 0.1)).abs() < 1e-6);
    }

    #[test]
    fn factor_directional_derivatives_strictly_obey_physics() {
        let (vars, pose_id, zwd_id) = make_test_vars(0);
        let clock_id = VariableId::new(2);
        let mut vars_mut = vars.clone();
        vars_mut
            .get_mut(&pose_id)
            .unwrap()
            .set_value(&[1_000_000.0, 2_000_000.0, 6_378_000.0, 0.0, 0.0, 0.0]);
        vars_mut.insert(
            clock_id,
            VariableNode::new(
                clock_id,
                VariableKind::ClockBias {
                    epoch: 0,
                    constellation_id: 0,
                },
            ),
        );

        let sat_pos = Vector3::new(15_000_000.0, 20_000_000.0, 25_000_000.0);
        let geom_range = (sat_pos - Vector3::new(1_000_000.0, 2_000_000.0, 6_378_000.0)).norm();

        let obs1 = CorrectedObservation {
            satellite: 1,
            constellation_id: 0,
            pr_l1: geom_range,
            pr_l2: None,
            cp_l1: None,
            cp_l2: None,
            cp_l1_lli: None,
            doppler: 0.0,
            snr_dbhz: 45.0,
            sat_pos_ecef: sat_pos,
            sat_clock_m: 0.0,
            f1: 1575.42e6,
            f2: 1227.60e6,
            freq_num: 0,
            elevation_rad: 0.5,
            tropo_dry_m: 0.0,
            tropo_map_wet: 0.0,
            iono_l1_m: 0.0,
            variance_m2: 0.25,
            cp_variance_m2: 9e-6,
        };
        let mut obs2 = obs1.clone();
        obs2.sat_clock_m = 10.0;

        let f1 = build_pseudorange_factor(&obs1, 0, pose_id, Some(clock_id), Some(zwd_id), None);
        let f2 = build_pseudorange_factor(&obs2, 0, pose_id, Some(clock_id), Some(zwd_id), None);
        let values = VariableValues::build(&vars_mut);

        let r1 = f1.residual(&values)[0];
        let r2 = f2.residual(&values)[0];

        assert!((r2 - r1 - 10.0).abs() < 1e-9);

        vars_mut
            .get_mut(&clock_id)
            .unwrap()
            .set_value(&[10.0, 0.0, 0.0]);
        let values_clk10 = VariableValues::build(&vars_mut);
        let r_clk10 = f1.residual(&values_clk10)[0];

        assert!((r_clk10 - r1 + 10.0).abs() < 1e-9);
    }
}
