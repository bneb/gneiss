pub mod covariance;
pub mod inversion;
pub mod thresholding;

pub type CovMatrix = nalgebra::DMatrix<f64>;
pub type JacMatrix = nalgebra::DMatrix<f64>;
