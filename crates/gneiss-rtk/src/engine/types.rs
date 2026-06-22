use serde::{Deserialize, Serialize};

/// Defines the core positioning mode
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
pub enum EngineMode {
    Spp,
    SppIns,
    SppInsLooselyCoupled,
    #[default]
    Rtk,
    RtkIns,
    RtkInsLooselyCoupled,
    Ppp,
    PppIns,
    PppInsLooselyCoupled,
    PppIekf,
    PppInsIekf,
    RtkInsIekf,
}

impl EngineMode {
    pub fn is_tightly_coupled(&self) -> bool {
        matches!(
            self,
            Self::SppIns | Self::RtkIns | Self::PppIns | Self::PppInsIekf | Self::RtkInsIekf
        )
    }

    pub fn is_ppp(&self) -> bool {
        matches!(
            self,
            Self::Ppp
                | Self::PppIns
                | Self::PppInsLooselyCoupled
                | Self::PppIekf
                | Self::PppInsIekf
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum DynamicsModel {
    Static,
    Pedestrian,
    Automotive,
    Marine,
    Airborne,
}

#[derive(Debug, Clone, PartialEq)]
pub enum EngineError {
    NoObservations,
    InitialSppFailed,
    StateDisappeared,
    InsufficientSatellites,
    MissingBasePosition,
    GeodeticMismatch(&'static str),
}

impl std::fmt::Display for EngineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EngineError::NoObservations => write!(f, "No observations available"),
            EngineError::InitialSppFailed => write!(f, "Initial SPP failed"),
            EngineError::StateDisappeared => write!(f, "EKF state disappeared mid-execution"),
            EngineError::InsufficientSatellites => {
                write!(f, "Insufficient satellites for EKF update")
            }
            EngineError::MissingBasePosition => write!(f, "Base station position must be provided"),
            EngineError::GeodeticMismatch(msg) => write!(f, "Geodetic Gatekeeper Failed: {}", msg),
        }
    }
}
impl std::error::Error for EngineError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_engine_mode_default_is_rtk() {
        assert_eq!(EngineMode::default(), EngineMode::Rtk);
    }

    #[test]
    fn test_is_tightly_coupled() {
        assert!(EngineMode::SppIns.is_tightly_coupled());
        assert!(EngineMode::RtkIns.is_tightly_coupled());
        assert!(EngineMode::PppIns.is_tightly_coupled());
        assert!(EngineMode::PppInsIekf.is_tightly_coupled());
        assert!(EngineMode::RtkInsIekf.is_tightly_coupled());

        assert!(!EngineMode::Spp.is_tightly_coupled());
        assert!(!EngineMode::SppInsLooselyCoupled.is_tightly_coupled());
        assert!(!EngineMode::Rtk.is_tightly_coupled());
        assert!(!EngineMode::RtkInsLooselyCoupled.is_tightly_coupled());
        assert!(!EngineMode::Ppp.is_tightly_coupled());
        assert!(!EngineMode::PppInsLooselyCoupled.is_tightly_coupled());
        assert!(!EngineMode::PppIekf.is_tightly_coupled());
    }

    #[test]
    fn test_is_ppp() {
        assert!(EngineMode::Ppp.is_ppp());
        assert!(EngineMode::PppIns.is_ppp());
        assert!(EngineMode::PppInsLooselyCoupled.is_ppp());
        assert!(EngineMode::PppIekf.is_ppp());
        assert!(EngineMode::PppInsIekf.is_ppp());

        assert!(!EngineMode::Spp.is_ppp());
        assert!(!EngineMode::SppIns.is_ppp());
        assert!(!EngineMode::SppInsLooselyCoupled.is_ppp());
        assert!(!EngineMode::Rtk.is_ppp());
        assert!(!EngineMode::RtkIns.is_ppp());
        assert!(!EngineMode::RtkInsLooselyCoupled.is_ppp());
        assert!(!EngineMode::RtkInsIekf.is_ppp());
    }

    #[test]
    fn test_engine_error_display() {
        assert_eq!(
            format!("{}", EngineError::NoObservations),
            "No observations available"
        );
        assert_eq!(
            format!("{}", EngineError::InitialSppFailed),
            "Initial SPP failed"
        );
        assert_eq!(
            format!("{}", EngineError::StateDisappeared),
            "EKF state disappeared mid-execution"
        );
        assert_eq!(
            format!("{}", EngineError::InsufficientSatellites),
            "Insufficient satellites for EKF update"
        );
        assert_eq!(
            format!("{}", EngineError::MissingBasePosition),
            "Base station position must be provided"
        );
    }

    #[test]
    fn test_engine_error_geodetic_mismatch_display() {
        let err = EngineError::GeodeticMismatch("datum mismatch");
        assert_eq!(
            format!("{}", err),
            "Geodetic Gatekeeper Failed: datum mismatch"
        );
    }

    #[test]
    fn test_engine_error_implements_std_error() {
        fn assert_error<E: std::error::Error>() {}
        assert_error::<EngineError>();
    }

    #[test]
    fn test_engine_mode_partial_eq() {
        assert_eq!(EngineMode::Spp, EngineMode::Spp);
        assert_ne!(EngineMode::Spp, EngineMode::Rtk);
    }

    #[test]
    fn test_dynamics_model_variants() {
        assert_eq!(DynamicsModel::Static as u8, 0);
        assert_eq!(DynamicsModel::Pedestrian as u8, 1);
        assert_eq!(DynamicsModel::Automotive as u8, 2);
        assert_eq!(DynamicsModel::Marine as u8, 3);
        assert_eq!(DynamicsModel::Airborne as u8, 4);
    }

    #[test]
    fn test_dynamics_model_debug_and_clone() {
        let model = DynamicsModel::Automotive;
        let cloned = model;
        assert_eq!(format!("{:?}", cloned), "Automotive");
    }
}
