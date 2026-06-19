//! Factor Graph Optimization (FGO) Engine
//!
//! This module provides a high-performance, sliding-window Factor Graph Optimization
//! backend for deeply coupled GNSS/INS integration. It replaces the traditional
//! Extended Kalman Filter (EKF) by optimizing over a trajectory of states simultaneously,
//! allowing for better handling of non-linearities and re-linearization of past states.
//!
//! # Architecture
//!
//! ```mermaid
//! graph TD
//!     subgraph Variables
//!         P[Pose x_i] --> V[Vel v_i]
//!         V --> B[Biases b_i]
//!     end
//!     
//!     subgraph Factors
//!         IMU(IMU Pre-integration)
//!         PR(GNSS Pseudorange)
//!         CP(GNSS Carrier Phase)
//!         DOP(GNSS Doppler)
//!         PRIOR(Marginalization Prior)
//!     end
//!     
//!     PRIOR -.->|constrains| P
//!     IMU -.->|links| P
//!     IMU -.->|links| P_next[Pose x_{i+1}]
//!     PR -.->|constrains| P
//!     CP -.->|constrains| P
//!     DOP -.->|constrains| V
//! ```
//!
//! # Components
//!
//! - **Graph (`graph.rs`)**: Manages the collection of Variables and Factors.
//! - **Solver (`solver.rs`)**: Levenberg-Marquardt non-linear optimizer.
//! - **Schur (`schur.rs`)**: Schur complement marginalization algebra.
//! - **Loss (`loss.rs`)**: Robust loss functions (Huber, Cauchy) for outlier rejection.
//! - **Multipath (`multipath.rs`)**: Retroactive GNSS multipath mitigation.

pub mod factors;
pub mod graph;
pub mod loss;
pub mod multipath;
pub mod schur;
pub mod solver;
pub mod variable;
