//! TOML configuration deserialisation for simulation jobs.

use lumina_geometry::scene::SceneSpec;
use serde::Deserialize;

/// Top-level job configuration.
#[derive(Debug, Deserialize)]
pub struct JobConfig {
    pub simulation: SimulationConfig,
    /// Multi-object scene: deserialises `[[geometry.object]]` arrays.
    pub geometry: SceneSpec,
    #[serde(default)]
    pub output: OutputConfig,
}

/// Simulation parameters from TOML.
#[derive(Debug, Deserialize)]
pub struct SimulationConfig {
    pub wavelengths: WavelengthSpec,
    #[serde(default = "default_env_n")]
    pub environment_n: f64,
    #[serde(default = "default_solver_tolerance")]
    pub solver_tolerance: f64,
    #[serde(default = "default_max_iterations")]
    pub max_iterations: usize,
    /// Compute backend: "auto", "cpu", or "gpu". Default: "auto".
    #[serde(default = "default_backend")]
    pub backend: String,
}

fn default_backend() -> String {
    "auto".into()
}

fn default_env_n() -> f64 {
    1.0
}
fn default_solver_tolerance() -> f64 {
    1e-6
}
fn default_max_iterations() -> usize {
    1000
}

/// Wavelength specification: either a range or explicit list.
#[derive(Debug, Deserialize)]
#[serde(untagged)]
pub enum WavelengthSpec {
    Range {
        range: [f64; 2],
        points: usize,
    },
    List {
        values: Vec<f64>,
    },
}


/// Output configuration.
#[derive(Debug, Deserialize)]
pub struct OutputConfig {
    /// Output directory (default: "./output").
    #[serde(default = "default_output_dir")]
    pub directory: String,
    /// Whether to save spectra as CSV (default: true).
    #[serde(default = "default_true")]
    pub save_spectra: bool,
    /// Whether to also save spectra as JSON (default: false).
    #[serde(default)]
    pub save_json: bool,
    /// Whether to compute and save the near-field map at peak extinction (default: false).
    #[serde(default)]
    pub save_near_field: bool,
}

impl Default for OutputConfig {
    fn default() -> Self {
        Self {
            directory: default_output_dir(),
            save_spectra: true,
            save_json: false,
            save_near_field: false,
        }
    }
}

fn default_output_dir() -> String {
    "./output".into()
}
fn default_true() -> bool {
    true
}

/// Load and parse a TOML job configuration file.
pub fn load_config(path: &std::path::Path) -> anyhow::Result<JobConfig> {
    let content = std::fs::read_to_string(path)?;
    let config: JobConfig = toml::from_str(&content)?;
    Ok(config)
}
