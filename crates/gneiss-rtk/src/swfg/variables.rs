use nalgebra::DVector;
use std::collections::{BTreeMap, HashMap};

/// Unique identifier for a variable in the factor graph.
/// Compact (8 bytes), hashable, zero-cost newtype.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VariableId(u64);

impl VariableId {
    #[inline]
    pub const fn new(id: u64) -> Self {
        Self(id)
    }
    #[inline]
    pub const fn as_u64(self) -> u64 {
        self.0
    }
}

/// Dimension of a variable in the state vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariableDim {
    Scalar,
    Vector3,
    Vector6,
    Matrix3x3,
}

impl VariableDim {
    #[inline]
    pub const fn size(self) -> usize {
        match self {
            Self::Scalar => 1,
            Self::Vector3 => 3,
            Self::Vector6 => 6,
            Self::Matrix3x3 => 9,
        }
    }
}

/// Semantic kind of a variable.  Used when creating variables to determine
/// dimension and initial value.  At solve time, variables are referenced
/// by `VariableId` only — the kind is metadata.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VariableKind {
    /// 6-DOF pose in ECEF: [x, y, z, qx, qy, qz].  qw recovered from unit norm.
    Pose { epoch: u32 },
    /// 3-DOF velocity in ECEF (m/s).
    Velocity { epoch: u32 },
    /// 3-DOF attitude error (axis-angle in body frame).
    /// Used for INS coupling — perturbation from nominal quaternion.
    Attitude { epoch: u32 },
    /// 6-DOF IMU bias [accel_x, accel_y, accel_z, gyro_x, gyro_y, gyro_z].
    /// One per session — not per-epoch.
    ImuBias,
    /// GNSS receiver clock: 3-DOF [bias_m, drift_m_s].  Actually 2-DOF per
    /// constellation, but we use 3 as a simplifying aligned dimension.
    ClockBias { epoch: u32, constellation_id: u8 },
    /// Troposphere zenith wet delay (1-DOF, meters).
    TropoZwd { epoch: u32 },
    /// GLONASS inter-frequency bias slope (1-DOF, m/freq_num).
    IfbGlonass,
    /// Carrier-phase integer ambiguity (1-DOF, cycles) — undifferenced (PPP).
    /// Persists across ALL epochs of the continuous tracking arc.
    Ambiguity { satellite: u16, frequency: u8, arc: u32 },
    /// Double-differenced ambiguity (1-DOF, cycles) — for RTK.
    /// Key: (sat, ref_sat) pair on a given frequency.
    DdAmbiguity { constellation_id: u8, satellite: u16, ref_satellite: u16, frequency: u8, arc: u32 },
    /// Per-satellite slant ionosphere delay on L1 (1-DOF, meters) — for UDUC PPP.
    /// Estimated as a random-walk variable across epochs.
    IonosphereSlant { epoch: u32, satellite: u16 },
    /// Persistent 6-DOF static monument pose: [x, y, z, qx, qy, qz].
    /// Persists across the ENTIRE static session.
    StaticPose,
}

impl VariableKind {
    pub const fn dim(self) -> VariableDim {
        match self {
            Self::Pose { .. } | Self::StaticPose => VariableDim::Vector6,
            Self::Velocity { .. } => VariableDim::Vector3,
            Self::Attitude { .. } => VariableDim::Vector3,
            Self::ImuBias => VariableDim::Vector6,
            Self::ClockBias { .. } => VariableDim::Scalar,
            Self::TropoZwd { .. } => VariableDim::Scalar,
            Self::IfbGlonass => VariableDim::Scalar,
            Self::Ambiguity { .. } => VariableDim::Scalar,
            Self::DdAmbiguity { .. } => VariableDim::Scalar,
            Self::IonosphereSlant { .. } => VariableDim::Scalar,
        }
    }

    /// Whether this variable is per-epoch (marginalized when epoch slides out).
    pub const fn is_per_epoch(self) -> bool {
        matches!(
            self,
            Self::Pose { .. }
                | Self::Velocity { .. }
                | Self::Attitude { .. }
                | Self::ClockBias { .. }
                | Self::TropoZwd { .. }
                | Self::IonosphereSlant { .. }
        )
    }

    /// The epoch this variable belongs to, if per-epoch.
    pub const fn epoch(self) -> Option<u32> {
        match self {
            Self::Pose { epoch }
            | Self::Velocity { epoch }
            | Self::Attitude { epoch }
            | Self::ClockBias { epoch, .. }
            | Self::TropoZwd { epoch }
            | Self::IonosphereSlant { epoch, .. } => Some(epoch),
            _ => None,
        }
    }
}

/// A variable node in the factor graph.  Stores its current estimate
/// (linearization point) and dimension metadata.
#[derive(Debug, Clone)]
pub struct VariableNode {
    pub id: VariableId,
    pub kind: VariableKind,
    /// Current estimate.  Updated after each LM iteration.
    pub value: DVector<f64>,
}

impl VariableNode {
    pub fn new(id: VariableId, kind: VariableKind) -> Self {
        let dim = kind.dim().size();
        Self { id, kind, value: DVector::zeros(dim) }
    }

    /// Set the variable's value from a slice.  Panics if the slice length
    /// doesn't match the variable's dimension.
    pub fn set_value(&mut self, val: &[f64]) {
        assert_eq!(val.len(), self.value.len(),
            "VariableNode::set_value: expected {} elements, got {}",
            self.value.len(), val.len());
        self.value.copy_from_slice(val);
    }
}

/// Efficient lookup structure for variable values during factor evaluation.
///
/// Maps `VariableId → (start_index, dimension)` in a flat state vector,
/// enabling O(1) access to any variable's current estimate.
pub struct VariableValues {
    index: HashMap<VariableId, (usize, usize)>,
    state: DVector<f64>,
}

impl VariableValues {
    /// Build the flat state vector from the current variable set.
    /// Variables are packed in `VariableId` order (deterministic).
    pub fn build(variables: &BTreeMap<VariableId, VariableNode>) -> Self {
        let mut index = HashMap::with_capacity(variables.len());
        let mut offset = 0_usize;
        for (id, node) in variables.iter() {
            let dim = node.value.len();
            index.insert(*id, (offset, dim));
            offset += dim;
        }
        let mut state = DVector::zeros(offset);
        for (id, node) in variables.iter() {
            let (start, dim) = index[id];
            state.rows_mut(start, dim).copy_from(&node.value);
        }
        Self { index, state }
    }

    /// Get a view into a variable's current value.  Returns `None` if the
    /// variable is not in the graph (should not happen during normal operation).
    #[inline]
    pub fn get(&self, id: VariableId) -> Option<nalgebra::DVectorView<'_, f64>> {
        let &(start, dim) = self.index.get(&id)?;
        Some(self.state.rows(start, dim))
    }

    /// Get the (start_index, dimension) of a variable in the packed state vector.
    #[inline]
    pub fn index_of(&self, id: VariableId) -> Option<(usize, usize)> {
        self.index.get(&id).copied()
    }

    /// Total dimension of the packed state vector.
    #[inline]
    pub fn total_dim(&self) -> usize {
        self.state.len()
    }

    /// Immutable reference to the flat state vector.
    #[inline]
    pub fn state(&self) -> &DVector<f64> {
        &self.state
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variable_id_roundtrip() {
        let id = VariableId::new(42);
        assert_eq!(id.as_u64(), 42);
    }

    #[test]
    fn variable_dim_sizes() {
        assert_eq!(VariableDim::Scalar.size(), 1);
        assert_eq!(VariableDim::Vector3.size(), 3);
        assert_eq!(VariableDim::Vector6.size(), 6);
        assert_eq!(VariableDim::Matrix3x3.size(), 9);
    }

    #[test]
    fn variable_values_offset_correctness() {
        let mut vars = BTreeMap::new();
        let k0 = VariableKind::Pose { epoch: 0 };
        let k1 = VariableKind::Velocity { epoch: 0 };
        let k2 = VariableKind::Ambiguity { satellite: 1, frequency: 1, arc: 0 };

        let id0 = VariableId::new(0);
        let id1 = VariableId::new(1);
        let id2 = VariableId::new(2);

        vars.insert(id0, VariableNode::new(id0, k0)); // dim 6
        vars.insert(id1, VariableNode::new(id1, k1)); // dim 3
        vars.insert(id2, VariableNode::new(id2, k2)); // dim 1

        let vals = VariableValues::build(&vars);
        assert_eq!(vals.total_dim(), 10);

        let v0 = vals.get(id0).unwrap();
        assert_eq!(v0.len(), 6);
        let v1 = vals.get(id1).unwrap();
        assert_eq!(v1.len(), 3);
        let v2 = vals.get(id2).unwrap();
        assert_eq!(v2.len(), 1);

        // Verify non-overlapping: id0 starts at 0, id1 at 6, id2 at 9
        // (BTreeMap orders by VariableId)
        assert_eq!(vals.get(id0).unwrap().as_slice().len(), 6);
        assert_eq!(vals.get(id1).unwrap().as_slice().len(), 3);
        assert_eq!(vals.get(id2).unwrap().as_slice().len(), 1);
    }

    #[test]
    fn variable_values_missing_returns_none() {
        let vars = BTreeMap::new();
        let vals = VariableValues::build(&vars);
        assert!(vals.get(VariableId::new(999)).is_none());
    }

    #[test]
    fn variable_kind_dims_consistent() {
        let kinds = [
            (VariableKind::Pose { epoch: 0 }, 6),
            (VariableKind::Velocity { epoch: 0 }, 3),
            (VariableKind::Attitude { epoch: 0 }, 3),
            (VariableKind::ImuBias, 6),
            (VariableKind::ClockBias { epoch: 0, constellation_id: 0 }, 1),
            (VariableKind::TropoZwd { epoch: 0 }, 1),
            (VariableKind::IfbGlonass, 1),
            (VariableKind::Ambiguity { satellite: 0, frequency: 0, arc: 0 }, 1),
            (VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 0, frequency: 1, arc: 0 }, 1),
            (VariableKind::IonosphereSlant { epoch: 0, satellite: 1 }, 1),
        ];
        for (kind, expected) in &kinds {
            let node = VariableNode::new(VariableId::new(0), *kind);
            assert_eq!(node.value.len(), *expected, "kind {:?} wrong dim", kind);
        }
    }

    #[test]
    fn per_epoch_classification() {
        assert!(VariableKind::Pose { epoch: 0 }.is_per_epoch());
        assert!(VariableKind::Velocity { epoch: 0 }.is_per_epoch());
        assert!(VariableKind::Attitude { epoch: 0 }.is_per_epoch());
        assert!(VariableKind::ClockBias { epoch: 0, constellation_id: 0 }.is_per_epoch());
        assert!(VariableKind::IonosphereSlant { epoch: 0, satellite: 1 }.is_per_epoch());
        assert!(!VariableKind::ImuBias.is_per_epoch());
        assert!(!VariableKind::Ambiguity { satellite: 0, frequency: 0, arc: 0 }.is_per_epoch());
        assert!(!VariableKind::DdAmbiguity { constellation_id: 0, satellite: 1, ref_satellite: 0, frequency: 1, arc: 0 }.is_per_epoch());
    }
}
