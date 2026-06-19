//! Machine Learning (ML) Engine for GNSS RAIM
//!
//! This module implements a Graph Neural Network (GNN) for Receiver Autonomous Integrity
//! Monitoring (RAIM). It replaces traditional residual-based fault detection (like
//! Chi-Squared testing) with a self-attention model capable of learning complex, non-linear
//! multipath and spoofing signatures from raw measurement features.
//!
//! # Architecture
//!
//! ```mermaid
//! graph TD
//!     subgraph Inputs
//!         SNR[SNR/CNO]
//!         EL[Elevation]
//!         AZ[Azimuth]
//!         DOP[Doppler]
//!         LOCK[Lock Time]
//!     end
//!
//!     subgraph Graph Neural Network
//!         MLP1[Embedding MLP]
//!         ATTN[Multi-Head Self-Attention]
//!         RES[Residual + ReLU]
//!         MLP2[Output Head MLP]
//!     end
//!     
//!     subgraph Output
//!         SIGMA[Predicted Log Variance]
//!     end
//!
//!     SNR & EL & AZ & DOP & LOCK -->|Features| MLP1
//!     MLP1 -->|Node Embeddings| ATTN
//!     ATTN -->|Message Passing| RES
//!     RES --> MLP2
//!     MLP2 --> SIGMA
//! ```
//!
//! # Components
//!
//! - **GNN RAIM (`gnn_raim.rs`)**: Core model architecture using `candle_core`.
//! - **Dataset (`dataset.rs`)**: Log-variance mapping and Heteroscedastic NLL loss.
//! - **Dataset Loader (`dataset_loader.rs`)**: Parsing and batching of CSV training data.

pub mod dataset;
pub mod dataset_loader;
pub mod gnn_raim;
