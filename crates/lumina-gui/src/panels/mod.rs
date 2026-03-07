//! GUI panels for the Lumina dashboard.

pub mod results;
pub mod scene;
pub mod simulation;

/// Canonical ordered list of material identifiers available in the GUI.
/// Keeps all material ComboBoxes consistent — update here when adding materials.
pub const AVAILABLE_MATERIALS: &[&str] =
    &["Au_JC", "Ag_JC", "Cu_JC", "TiO2_Palik", "SiO2_Palik"];
