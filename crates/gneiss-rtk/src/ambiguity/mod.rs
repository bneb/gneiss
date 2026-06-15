pub mod ffrt;
pub mod lambda;
pub mod par;

pub struct AmbiguityResolutionResult {
    pub fixed_state: crate::estimators::ekf::filter::RtkState,
    pub z_dd: nalgebra::DVector<f64>,
    pub d_full: nalgebra::DMatrix<f64>,
    pub ratio_test: f64,
    pub subset_size: usize,
}

pub struct AmbiguityFixResult {
    pub fixed_state: crate::estimators::ekf::filter::RtkState,
    pub z_dd: nalgebra::DVector<f64>,
    pub d_full: nalgebra::DMatrix<f64>,
}
