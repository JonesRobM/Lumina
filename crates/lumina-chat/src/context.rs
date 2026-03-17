//! Serialisation of current simulation state for injection into the chat.

use lumina_core::types::CrossSections;
use lumina_geometry::scene::SceneSpec;

/// An owned snapshot of the simulation state at the moment the user clicks
/// "Attach context". Passed to the AI as a formatted user message.
pub struct ContextSnapshot {
    /// SceneSpec serialised as JSON (owned; avoids borrow conflicts in egui update).
    pub scene_json: String,
    /// Human-readable summary of the latest results, or None if no run has completed.
    pub spectra_summary: Option<String>,
}

impl ContextSnapshot {
    /// Construct a snapshot. All data is cloned/serialised immediately.
    ///
    /// - `spectra`: completed sweep results from `ResultsPanel.spectra`
    /// - `shg_spectra`: `(fundamental_nm, shg_intensity_nm6)` pairs from
    ///   `ResultsPanel.shg_spectra`; pass `None` if SHG was not computed
    pub fn from_scene(
        scene: &SceneSpec,
        spectra: Option<&[CrossSections]>,
        shg_spectra: Option<&[(f64, f64)]>,
    ) -> Self {
        let scene_json = serde_json::to_string_pretty(scene)
            .unwrap_or_else(|_| "{}".to_string());

        let spectra_summary = spectra.map(|s| build_spectra_summary(s, shg_spectra));

        Self { scene_json, spectra_summary }
    }

    /// Format the snapshot as a markdown user message.
    pub fn to_message_content(&self) -> String {
        let results_section = match &self.spectra_summary {
            Some(summary) => summary.clone(),
            None => "No simulation results available yet.".to_string(),
        };

        format!(
            "## Current simulation context\n\n\
             {results_section}\n\n\
             <details>\n<summary>Full scene JSON</summary>\n\n\
             ```json\n{}\n```\n</details>",
            self.scene_json,
        )
    }
}

fn build_spectra_summary(
    spectra: &[CrossSections],
    shg_spectra: Option<&[(f64, f64)]>,
) -> String {
    if spectra.is_empty() {
        return "Results array is empty.".to_string();
    }

    // Find peak extinction.
    let peak = spectra
        .iter()
        .filter(|cs| cs.extinction.is_finite())
        .max_by(|a, b| a.extinction.partial_cmp(&b.extinction).unwrap());

    let mut lines = Vec::new();

    if let Some(p) = peak {
        let abs_ratio = if p.extinction > 0.0 {
            p.absorption / p.extinction
        } else {
            0.0
        };
        lines.push(format!(
            "**Latest results:** peak C_ext = {:.3e} nm² at {:.1} nm, \
             C_abs/C_ext = {:.2}",
            p.extinction, p.wavelength_nm, abs_ratio,
        ));

        if let Some(cd) = p.circular_dichroism {
            lines.push(format!("CD at peak = {cd:.3e} nm²"));
        }
    }

    if let Some(shg) = shg_spectra {
        if !shg.is_empty() {
            let peak_shg = shg
                .iter()
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
            if let Some((wl, i)) = peak_shg {
                lines.push(format!(
                    "**SHG:** peak intensity = {i:.3e} nm⁶ at λ_fund = {wl:.1} nm"
                ));
            }
        }
    }

    lines.join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_results_includes_note() {
        let scene = lumina_geometry::scene::SceneSpec::default();
        let snap = ContextSnapshot::from_scene(&scene, None, None);
        let content = snap.to_message_content();
        assert!(content.contains("No simulation results available yet"));
        assert!(content.contains("simulation context"));
    }

    #[test]
    fn scene_json_is_valid() {
        let scene = lumina_geometry::scene::SceneSpec::default();
        let snap = ContextSnapshot::from_scene(&scene, None, None);
        // The embedded JSON should parse cleanly.
        serde_json::from_str::<serde_json::Value>(&snap.scene_json).unwrap();
    }
}
