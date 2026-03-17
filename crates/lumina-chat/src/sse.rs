//! Server-Sent Events parser for OpenAI-compatible streaming responses.

use crate::types::ChatEvent;

/// Parse a single SSE `data:` line payload.
///
/// Return type: `Option<ChatEvent>`
/// - `"[DONE]"` → `Some(ChatEvent::Done)`
/// - Valid JSON with non-empty `choices[0].delta.content` → `Some(ChatEvent::Token(text))`
/// - Valid JSON with null/missing/empty content → `None` (no-op; don't emit empty tokens)
/// - Invalid JSON → `Some(ChatEvent::Error(description))`
pub fn parse_data_line(payload: &str) -> Option<ChatEvent> {
    if payload == "[DONE]" {
        return Some(ChatEvent::Done);
    }

    match serde_json::from_str::<serde_json::Value>(payload) {
        Err(_) => Some(ChatEvent::Error(format!("Malformed SSE payload: {payload}"))),
        Ok(v) => {
            // Extract choices[0].delta.content — null or missing → None
            let content = v
                .get("choices")
                .and_then(|c| c.get(0))
                .and_then(|ch| ch.get("delta"))
                .and_then(|d| d.get("content"))
                .and_then(|c| c.as_str());

            content
                .filter(|s| !s.is_empty())
                .map(|s| ChatEvent::Token(s.to_string()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn done_signal() {
        let ev = parse_data_line("[DONE]").unwrap();
        assert!(matches!(ev, ChatEvent::Done));
    }

    #[test]
    fn valid_token() {
        let line = r#"{"choices":[{"delta":{"content":"Hello"}}]}"#;
        let ev = parse_data_line(line).unwrap();
        assert!(matches!(ev, ChatEvent::Token(t) if t == "Hello"));
    }

    #[test]
    fn empty_content_returns_none() {
        // Role-only delta at start of stream — no content key
        let line = r#"{"choices":[{"delta":{"role":"assistant"}}]}"#;
        assert!(parse_data_line(line).is_none());
    }

    #[test]
    fn null_content_returns_none() {
        let line = r#"{"choices":[{"delta":{"content":null}}]}"#;
        assert!(parse_data_line(line).is_none());
    }

    #[test]
    fn malformed_json_returns_error() {
        let ev = parse_data_line("not json at all").unwrap();
        assert!(matches!(ev, ChatEvent::Error(_)));
    }
}
