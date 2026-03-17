//! Chat panel: provider settings and streaming conversation UI.

use eframe::egui::{self, Color32, RichText, ScrollArea, TextEdit};
use lumina_chat::client::SYSTEM_PROMPT;
use lumina_chat::context::ContextSnapshot;
use lumina_chat::{ChatClient, ChatConfig, ChatEvent, Message, Role};

/// Which sub-tab is active within the chat panel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChatTab {
    Settings,
    Chat,
}

impl Default for ChatTab {
    fn default() -> Self {
        Self::Settings
    }
}

/// All state for the chat panel.
pub struct ChatPanel {
    pub tab: ChatTab,
    pub client: ChatClient,
    /// Editable copy of the config shown in the Settings tab.
    pub config_draft: ChatConfig,
    /// Conversation history.
    pub history: Vec<Message>,
    /// Current text input buffer.
    pub input: String,
    /// True while a streaming response is in progress.
    pub streaming: bool,
    /// Set after Clear — events from any in-flight request are discarded until Done/Error.
    pub cleared: bool,
    /// Set to true on the first update() call so we only call set_egui_ctx once.
    pub egui_ctx_set: bool,
    /// Set by the "Attach context" button; polled by LuminaApp to call attach_context().
    pub attach_requested: bool,
}

impl Default for ChatPanel {
    fn default() -> Self {
        let config = ChatConfig::load();
        Self {
            tab: ChatTab::Settings,
            client: ChatClient::new(config.clone()),
            config_draft: config,
            history: Vec::new(),
            input: String::new(),
            streaming: false,
            cleared: false,
            egui_ctx_set: false,
            attach_requested: false,
        }
    }
}

impl ChatPanel {
    /// Main egui render. Call each frame from `LuminaApp::update`.
    pub fn ui(&mut self, ui: &mut egui::Ui) {
        // Lazily inject egui context on the first frame.
        if !self.egui_ctx_set {
            self.client.set_egui_ctx(ui.ctx().clone());
            self.egui_ctx_set = true;
        }

        // Drain streaming events from the background thread.
        self.drain_events();

        // Top bar: tab selector + Clear button.
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.tab, ChatTab::Settings, "Settings");
            ui.selectable_value(&mut self.tab, ChatTab::Chat, "Chat");
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.button("Clear").clicked() {
                    self.history.clear();
                    self.input.clear();
                    self.streaming = false;
                    self.cleared = true;
                }
            });
        });
        ui.separator();

        match self.tab {
            ChatTab::Settings => self.settings_ui(ui),
            ChatTab::Chat => self.chat_ui(ui),
        }
    }

    /// Inject a pre-built context snapshot as a user message and trigger a send.
    /// Called from `LuminaApp::attach_context`.
    pub fn inject_context_snapshot(&mut self, snapshot: ContextSnapshot) {
        let content = snapshot.to_message_content();
        self.history.push(Message { role: Role::User, content });
        self.history.push(Message { role: Role::Assistant, content: String::new() });
        self.streaming = true;
        self.cleared = false;
        let send_history = self.history[..self.history.len() - 1].to_vec();
        self.client.send(send_history, SYSTEM_PROMPT);
    }

    /// Returns true and resets the flag if "Attach context" was clicked this frame.
    pub fn take_attach_requested(&mut self) -> bool {
        std::mem::take(&mut self.attach_requested)
    }

    // ── Private helpers ───────────────────────────────────────────────────────

    fn drain_events(&mut self) {
        while let Some(event) = self.client.poll() {
            if self.cleared {
                // Discard events from any request that was in flight when Clear was pressed.
                if matches!(event, ChatEvent::Done | ChatEvent::Error(_)) {
                    self.cleared = false;
                    self.streaming = false;
                }
                continue;
            }
            match event {
                ChatEvent::Token(token) => {
                    if let Some(last) = self.history.last_mut() {
                        if last.role == Role::Assistant {
                            last.content.push_str(&token);
                        }
                    }
                }
                ChatEvent::Done => {
                    self.streaming = false;
                }
                ChatEvent::Error(msg) => {
                    self.streaming = false;
                    self.history.push(Message {
                        role: Role::Assistant,
                        content: format!("⚠ Error: {msg}"),
                    });
                }
            }
        }
    }

    fn settings_ui(&mut self, ui: &mut egui::Ui) {
        egui::Grid::new("chat_settings_grid")
            .num_columns(2)
            .spacing([8.0, 6.0])
            .show(ui, |ui| {
                ui.label("Provider:");
                egui::ComboBox::from_id_salt("provider_combo")
                    .selected_text(&self.config_draft.provider)
                    .show_ui(ui, |ui| {
                        for preset in &["groq", "ollama", "openrouter", "custom"] {
                            if ui
                                .selectable_label(
                                    self.config_draft.provider == *preset,
                                    *preset,
                                )
                                .clicked()
                            {
                                if *preset != "custom" {
                                    self.config_draft = ChatConfig::from_preset(preset);
                                } else {
                                    self.config_draft.provider = "custom".into();
                                }
                            }
                        }
                    });
                ui.end_row();

                ui.label("Base URL:");
                ui.add(
                    TextEdit::singleline(&mut self.config_draft.base_url)
                        .desired_width(300.0),
                );
                ui.end_row();

                ui.label("Model:");
                ui.add(
                    TextEdit::singleline(&mut self.config_draft.model)
                        .desired_width(300.0),
                );
                ui.end_row();

                ui.label("API Key:");
                let key_str = self.config_draft.api_key.get_or_insert_with(String::new);
                ui.add(
                    TextEdit::singleline(key_str)
                        .password(true)
                        .desired_width(300.0),
                );
                // Normalise empty string to None.
                if self.config_draft.api_key.as_deref() == Some("") {
                    self.config_draft.api_key = None;
                }
                ui.end_row();
            });

        ui.add_space(8.0);
        if ui.button("Save").clicked() {
            self.config_draft.save();
            self.client.update_config(self.config_draft.clone());
        }

        ui.add_space(8.0);
        ui.separator();
        ui.label(RichText::new("Free providers:").weak());
        ui.label("• Ollama — local, no key needed (ollama.com)");
        ui.label("• Groq — free tier key at console.groq.com");
        ui.label("• OpenRouter — free models at openrouter.ai");
    }

    fn chat_ui(&mut self, ui: &mut egui::Ui) {
        let available = ui.available_height();
        let history_height = (available - 60.0).max(100.0);

        ScrollArea::vertical()
            .max_height(history_height)
            .auto_shrink([false; 2])
            .stick_to_bottom(true)
            .show(ui, |ui| {
                for msg in &self.history {
                    match msg.role {
                        Role::User => {
                            ui.horizontal_wrapped(|ui| {
                                ui.label(RichText::new("You:").strong());
                                ui.label(&msg.content);
                            });
                        }
                        Role::Assistant => {
                            ui.horizontal_wrapped(|ui| {
                                ui.label(
                                    RichText::new("AI:")
                                        .strong()
                                        .color(Color32::from_rgb(100, 180, 255)),
                                );
                                if msg.content.is_empty() && self.streaming {
                                    ui.spinner();
                                } else {
                                    ui.label(&msg.content);
                                }
                            });
                        }
                        Role::System => {}
                    }
                    ui.add_space(4.0);
                }
            });

        ui.separator();

        // Input row.
        let send_triggered = std::cell::Cell::new(false);
        ui.add_enabled_ui(!self.streaming, |ui| {
            ui.horizontal(|ui| {
                let response = ui.add(
                    TextEdit::singleline(&mut self.input)
                        .desired_width(ui.available_width() - 60.0)
                        .hint_text("Ask anything…"),
                );
                if response.lost_focus()
                    && ui.input(|i| i.key_pressed(egui::Key::Enter))
                    && !self.input.trim().is_empty()
                {
                    send_triggered.set(true);
                }
                if ui.button("→").clicked() && !self.input.trim().is_empty() {
                    send_triggered.set(true);
                }
            });
        });

        if send_triggered.get() {
            self.do_send();
        }

        ui.add_space(4.0);
        ui.add_enabled_ui(!self.streaming, |ui| {
            if ui.button("Attach context").clicked() {
                self.attach_requested = true;
            }
        });
    }

    fn do_send(&mut self) {
        let content = std::mem::take(&mut self.input);
        self.history.push(Message { role: Role::User, content });
        self.history.push(Message { role: Role::Assistant, content: String::new() });
        self.streaming = true;
        self.cleared = false;
        let send_history = self.history[..self.history.len() - 1].to_vec();
        self.client.send(send_history, SYSTEM_PROMPT);
    }
}
