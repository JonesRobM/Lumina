//! OpenAI-compatible chat client for LuminaCDA.

pub mod config;
pub mod context;
pub mod client;
pub mod sse;
pub mod types;

pub use config::ChatConfig;
pub use client::ChatClient;
pub use types::{ChatEvent, Message, Role};
pub use context::ContextSnapshot;
