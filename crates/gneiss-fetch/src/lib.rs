#![cfg_attr(test, allow(clippy::unwrap_used))]

pub mod hatanaka;
pub mod provider;
pub mod sources;
pub mod uri;

pub use uri::ResourceResolver;
