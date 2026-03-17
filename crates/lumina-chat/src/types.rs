//! Core message and event types for the chat client.
//!
//! Note: `ChatTask` (the internal worker task) lives in `client.rs` because it
//! depends on `ChatConfig` — keeping it there avoids a circular module dependency.

/// Role of a message participant. Serialised to lowercase for the OpenAI wire format.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Role {
    System,
    User,
    Assistant,
}

/// A single message in the conversation history.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Message {
    pub role: Role,
    pub content: String,
}

/// Events emitted by the background streaming task to the GUI thread.
#[derive(Debug)]
pub enum ChatEvent {
    /// A token chunk from the streaming response.
    Token(String),
    /// The stream has completed successfully.
    Done,
    /// An error occurred; the string is a human-readable description.
    Error(String),
}