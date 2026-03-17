//! Background streaming HTTP client for OpenAI-compatible providers.

use std::sync::mpsc;

use crate::config::ChatConfig;
use crate::sse::parse_data_line;
use crate::types::{ChatEvent, Message, Role};

/// Static system prompt prepended to every request.
pub const SYSTEM_PROMPT: &str = "\
You are an AI assistant embedded in LuminaCDA, a Coupled Dipole Approximation \
tool for computing optical responses of nanostructures. Help the user understand \
simulation results, interpret spectra, and configure their geometry. \
Be concise and scientifically accurate.";

/// Internal task sent from the main thread to the background worker.
pub(crate) struct ChatTask {
    pub messages: Vec<Message>,
    pub config: ChatConfig,
    #[cfg(feature = "gui")]
    pub egui_ctx: Option<egui::Context>,
}

/// Streaming chat client. Sends HTTP requests on a dedicated background thread.
pub struct ChatClient {
    config: ChatConfig,
    task_tx: mpsc::SyncSender<ChatTask>,
    event_rx: mpsc::Receiver<ChatEvent>,
    #[cfg(feature = "gui")]
    egui_ctx: Option<egui::Context>,
    _worker: std::thread::JoinHandle<()>,
}

impl ChatClient {
    pub fn new(config: ChatConfig) -> Self {
        // Bounded channel: capacity 1 — one request in flight at a time.
        let (task_tx, task_rx) = mpsc::sync_channel::<ChatTask>(1);
        let (event_tx, event_rx) = mpsc::channel::<ChatEvent>();

        let worker = std::thread::Builder::new()
            .name("lumina-chat-worker".into())
            .spawn(move || {
                let rt = tokio::runtime::Builder::new_current_thread()
                    .enable_all()
                    .build()
                    .expect("failed to build tokio runtime");
                while let Ok(task) = task_rx.recv() {
                    rt.block_on(stream_request(task, &event_tx));
                }
                // task_tx dropped on ChatClient drop → recv() returns Err → thread exits cleanly
            })
            .expect("failed to spawn chat worker thread");

        Self {
            config,
            task_tx,
            event_rx,
            #[cfg(feature = "gui")]
            egui_ctx: None,
            _worker: worker,
        }
    }

    /// Store the egui context for repaint requests.
    /// Call this on the first `update()` frame from the GUI.
    #[cfg(feature = "gui")]
    pub fn set_egui_ctx(&mut self, ctx: egui::Context) {
        self.egui_ctx = Some(ctx);
    }

    /// Send a conversation to the provider.
    /// Prepends a system message and dispatches to the background thread.
    /// Non-blocking from the caller's perspective.
    pub fn send(&mut self, history: Vec<Message>, system_prompt: &str) {
        let mut messages = vec![Message {
            role: Role::System,
            content: system_prompt.to_string(),
        }];
        messages.extend(history);

        let task = ChatTask {
            messages,
            config: self.config.clone(),
            #[cfg(feature = "gui")]
            egui_ctx: self.egui_ctx.clone(),
        };

        // try_send: if the channel is full (a request is in flight), drop silently.
        // The Send button is disabled while streaming = true, so this is rare.
        if let Err(e) = self.task_tx.try_send(task) {
            log::warn!("Chat task dropped (channel full or disconnected): {e}");
        }
    }

    /// Poll for the next available event. Non-blocking.
    pub fn poll(&self) -> Option<ChatEvent> {
        self.event_rx.try_recv().ok()
    }

    /// Replace the active config. Takes effect on the next `send()` call.
    pub fn update_config(&mut self, config: ChatConfig) {
        self.config = config;
    }
}

/// Async task: POST to the provider, parse SSE, emit ChatEvents.
async fn stream_request(task: ChatTask, event_tx: &mpsc::Sender<ChatEvent>) {
    let url = format!(
        "{}/chat/completions",
        task.config.base_url.trim_end_matches('/')
    );

    let client = match reqwest::Client::builder()
        .timeout(std::time::Duration::from_secs(60))
        .build()
    {
        Ok(c) => c,
        Err(e) => {
            let _ = event_tx.send(ChatEvent::Error(format!("HTTP client error: {e}")));
            return;
        }
    };

    let body = serde_json::json!({
        "model": task.config.model,
        "messages": task.messages,
        "stream": true,
    });

    let mut req = client.post(&url).json(&body);
    if let Some(key) = &task.config.api_key {
        req = req.bearer_auth(key);
    }

    let response = match req.send().await {
        Ok(r) => r,
        Err(e) => {
            let _ = event_tx.send(ChatEvent::Error(format!("Request failed: {e}")));
            return;
        }
    };

    if !response.status().is_success() {
        let status = response.status();
        let body = response.text().await.unwrap_or_default();
        let _ = event_tx.send(ChatEvent::Error(format!("HTTP {status}: {body}")));
        return;
    }

    // Read SSE stream line by line.
    use futures_util::StreamExt;
    let mut stream = response.bytes_stream();
    let mut buffer = String::new();

    while let Some(chunk) = stream.next().await {
        let chunk = match chunk {
            Ok(b) => b,
            Err(e) => {
                let _ = event_tx.send(ChatEvent::Error(format!("Stream read error: {e}")));
                return;
            }
        };

        buffer.push_str(&String::from_utf8_lossy(&chunk));

        // Process all complete lines in the buffer.
        while let Some(newline_pos) = buffer.find('\n') {
            let line = buffer[..newline_pos].trim_end_matches('\r').to_string();
            buffer.drain(..=newline_pos);

            if let Some(payload) = line.strip_prefix("data: ") {
                if let Some(event) = parse_data_line(payload) {
                    let is_done = matches!(event, ChatEvent::Done);
                    #[cfg(feature = "gui")]
                    if let Some(ctx) = &task.egui_ctx {
                        ctx.request_repaint();
                    }
                    let _ = event_tx.send(event);
                    if is_done {
                        return;
                    }
                }
            }
        }
    }

    // Stream ended without explicit [DONE] — treat as completion.
    #[cfg(feature = "gui")]
    if let Some(ctx) = &task.egui_ctx {
        ctx.request_repaint();
    }
    let _ = event_tx.send(ChatEvent::Done);
}
