//! Lumina GUI application entry point.

mod app;
mod panels;

fn main() -> eframe::Result {
    env_logger::init();

    let options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 800.0])
            .with_min_inner_size([800.0, 600.0]),
        ..Default::default()
    };

    eframe::run_native(
        "Lumina",
        options,
        Box::new(|_cc| Ok(Box::new(app::LuminaApp::default()))),
    )
}
