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
    PppRtklib,
    RtkInsIekf,
    /// Multi-epoch sliding-window factor graph PPP.
    /// Breaks the single-epoch IEKF accuracy floor (~5m) by jointly
    /// optimizing position, clock, tropo, and ambiguities across N epochs
    /// with between-epoch dynamics constraints.
    PppMultiEpoch,
}

impl EngineMode {
    pub fn is_tightly_coupled(&self) -> bool {
        matches!(
            self,
            Self::SppIns | Self::RtkIns | Self::PppIns | Self::PppInsIekf | Self::RtkInsIekf | Self::PppRtklib
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
                | Self::PppRtklib
                | Self::PppMultiEpoch
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
#[derive(Default)]
pub enum DynamicsModel {
    Static,
    Pedestrian,
    #[default]
    Automotive,
    Marine,
    Airborne,
}


#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
#[derive(Default)]
pub enum IonosphereModel {
    /// Klobuchar model — free, broadcast parameters, 1–3 m accuracy.
    #[default]
    Klobuchar,
    /// IONEX grid maps — requires downloaded file, 1–5 cm accuracy.
    Ionex,
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
        assert!(EngineMode::PppMultiEpoch.is_ppp());

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
    fn test_engine_mode_debug_clone_copy() {
        // Verify Debug, Clone, and Copy work
        let mode = EngineMode::PppIekf;
        let _debug = format!("{:?}", mode);
        let cloned = mode;
        assert_eq!(mode, cloned);
    }

    #[test]
    fn test_engine_error_debug() {
        let err = EngineError::NoObservations;
        let debug = format!("{:?}", err);
        assert!(debug.contains("NoObservations"));
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
    fn test_dynamics_model_default() {
        assert_eq!(DynamicsModel::default(), DynamicsModel::Automotive);
    }

    #[test]
    fn test_dynamics_model_debug_and_clone() {
        let model = DynamicsModel::Automotive;
        let cloned = model;
        assert_eq!(format!("{:?}", cloned), "Automotive");
    }

    #[test]
    fn test_ionosphere_model_variants() {
        assert_eq!(IonosphereModel::Klobuchar, IonosphereModel::Klobuchar);
        assert_eq!(IonosphereModel::Ionex, IonosphereModel::Ionex);
        assert_ne!(IonosphereModel::Klobuchar, IonosphereModel::Ionex);
    }

    #[test]
    fn test_ionosphere_model_default() {
        assert_eq!(IonosphereModel::default(), IonosphereModel::Klobuchar);
    }

    #[test]
    fn test_ionosphere_model_debug_and_clone() {
        let model = IonosphereModel::Klobuchar;
        let cloned = model;
        assert_eq!(format!("{:?}", cloned), "Klobuchar");
    }

    #[test]
    fn test_ionosphere_model_partial_eq() {
        assert_eq!(IonosphereModel::Ionex, IonosphereModel::Ionex);
        assert_ne!(IonosphereModel::Klobuchar, IonosphereModel::Ionex);
    }

    #[test]
    fn test_engine_mode_serde_traits() {
        // Compile-time check: EngineMode implements Serialize + Deserialize
        fn assert_serde<T: serde::Serialize + serde::de::DeserializeOwned>() {}
        assert_serde::<EngineMode>();
    }

    #[test]
    fn test_dynamics_model_serde_traits() {
        fn assert_serde<T: serde::Serialize + serde::de::DeserializeOwned>() {}
        assert_serde::<DynamicsModel>();
    }

    #[test]
    fn test_ionosphere_model_serde_traits() {
        fn assert_serde<T: serde::Serialize + serde::de::DeserializeOwned>() {}
        assert_serde::<IonosphereModel>();
    }

    #[test]
    fn test_engine_mode_is_ppp_multiepoch() {
        assert!(EngineMode::PppMultiEpoch.is_ppp());
    }

    #[test]
    fn test_engine_mode_is_tightly_coupled_all_variants() {
        let tight: &[EngineMode] = &[
            EngineMode::SppIns,
            EngineMode::RtkIns,
            EngineMode::PppIns,
            EngineMode::PppInsIekf,
            EngineMode::RtkInsIekf,
        ];
        for mode in tight {
            assert!(mode.is_tightly_coupled(), "{:?} should be tightly coupled", mode);
        }
        let loose: &[EngineMode] = &[
            EngineMode::Spp,
            EngineMode::SppInsLooselyCoupled,
            EngineMode::Rtk,
            EngineMode::RtkInsLooselyCoupled,
            EngineMode::Ppp,
            EngineMode::PppInsLooselyCoupled,
            EngineMode::PppIekf,
            EngineMode::PppMultiEpoch,
        ];
        for mode in loose {
            assert!(!mode.is_tightly_coupled(), "{:?} should not be tightly coupled", mode);
        }
    }
}
