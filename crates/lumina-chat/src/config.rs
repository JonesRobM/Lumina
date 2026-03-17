//! Provider configuration and persistence.

use std::path::PathBuf;

/// Configuration for the chat provider. Serialisable to TOML.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ChatConfig {
    /// Provider identifier: "groq", "ollama", "openrouter", or "custom".
    pub provider: String,
    /// Base URL without trailing slash. The client appends `/chat/completions`.
    pub base_url: String,
    /// Model identifier.
    pub model: String,
    /// API key. None for providers that don't require one (e.g. local Ollama).
    pub api_key: Option<String>,
}

impl Default for ChatConfig {
    fn default() -> Self {
        Self::from_preset("groq")
    }
}

impl ChatConfig {
    /// Construct config from a named preset. Unknown names fall back to Groq.
    pub fn from_preset(name: &str) -> Self {
        match name {
            "ollama" => Self {
                provider: "ollama".into(),
                base_url: "http://localhost:11434/v1".into(),
                model: "llama3.2".into(),
                api_key: None,
            },
            "openrouter" => Self {
                provider: "openrouter".into(),
                base_url: "https://openrouter.ai/api/v1".into(),
                model: "meta-llama/llama-3.3-70b-instruct:free".into(),
                api_key: None,
            },
            _ => Self {
                provider: "groq".into(),
                base_url: "https://api.groq.com/openai/v1".into(),
                model: "llama-3.3-70b-versatile".into(),
                api_key: None,
            },
        }
    }

    /// Platform config file path. Returns None if the OS config dir is unavailable.
    pub fn config_path() -> Option<PathBuf> {
        dirs::config_dir().map(|d| d.join("lumina").join("chat.toml"))
    }

    /// Load from the platform config file. Returns default config on any failure.
    pub fn load() -> Self {
        let Some(path) = Self::config_path() else {
            return Self::default();
        };
        let Ok(text) = std::fs::read_to_string(&path) else {
            return Self::default();
        };
        toml::from_str(&text).unwrap_or_else(|e| {
            log::warn!("Failed to parse chat config ({path:?}): {e}. Using defaults.");
            Self::default()
        })
    }

    /// Save to the platform config file. Logs and ignores errors.
    pub fn save(&self) {
        let Some(path) = Self::config_path() else {
            log::warn!("Could not determine config dir; chat config not saved.");
            return;
        };
        if let Some(parent) = path.parent() {
            if let Err(e) = std::fs::create_dir_all(parent) {
                log::warn!("Could not create config dir {parent:?}: {e}");
                return;
            }
        }
        match toml::to_string_pretty(self) {
            Ok(text) => {
                let commented = format!(
                    "# LuminaCDA chat configuration.\n\
                     # WARNING: API key stored in plaintext. Do not commit this file.\n\n\
                     {text}"
                );
                if let Err(e) = std::fs::write(&path, commented) {
                    log::warn!("Could not write chat config to {path:?}: {e}");
                }
            }
            Err(e) => log::warn!("Could not serialise chat config: {e}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_toml() {
        // Tests the TOML serialisation/deserialisation path in memory.
        // Does NOT call load() or save() to avoid touching the real config dir.
        let cfg = ChatConfig {
            provider: "groq".to_string(),
            base_url: "https://api.groq.com/openai/v1".to_string(),
            model: "llama-3.3-70b-versatile".to_string(),
            api_key: Some("test_key".to_string()),
        };
        let toml_str = toml::to_string(&cfg).unwrap();
        let restored: ChatConfig = toml::from_str(&toml_str).unwrap();
        assert_eq!(cfg.provider, restored.provider);
        assert_eq!(cfg.base_url, restored.base_url);
        assert_eq!(cfg.api_key, restored.api_key);
    }

    #[test]
    fn preset_groq_url() {
        let cfg = ChatConfig::from_preset("groq");
        assert!(cfg.base_url.ends_with("/openai/v1"));
        assert!(!cfg.base_url.ends_with('/'));
    }

    #[test]
    fn preset_ollama_no_key() {
        let cfg = ChatConfig::from_preset("ollama");
        assert!(cfg.api_key.is_none());
        assert!(cfg.base_url.contains("11434"));
    }

    #[test]
    fn malformed_toml_falls_back_to_default() {
        // Confirm toml rejects malformed input (load() would return Self::default() here).
        let result = toml::from_str::<ChatConfig>("not valid toml !!!");
        assert!(result.is_err());
    }
}
