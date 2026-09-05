//! Interactive Visual GUI and Web Diagnostic Workspace.
//!
//! Provides a lightweight, embedded local web server delivering:
//! - Interactive Canvas and Leaflet trajectory mapping with color-coded fix status
//! - Multi-channel carrier phase residual, DOP, and satellite tracking inspector
//! - Polar azimuth-elevation skyplot viewer
//! - Automated CORS discovery and one-click mission processing

pub mod assets;
pub mod handlers;
pub mod server;

pub use server::{run_gui_server, GuiArgs};
