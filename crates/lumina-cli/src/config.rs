//! TOML configuration deserialisation for simulation jobs.

use lumina_core::solver::substrate::SubstrateSpec;
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
    /// Optional nonlinear optics configuration. Absent means linear-only.
    #[serde(default)]
    pub nonlinear: Option<NonlinearConfig>,
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
    /// Optional planar substrate below the nanostructure.
    ///
    /// Example: `[simulation.substrate]` with `z_interface = 0.0` and
    /// `material = "SiO2_Palik"`.
    #[serde(default)]
    pub substrate: Option<SubstrateSpec>,
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

/// Configuration for nonlinear optical response calculations.
///
/// Example `job.toml` section:
/// ```toml
/// [nonlinear]
/// enable_shg = true
/// symmetry = "isotropic_surface"
/// chi_zzz = [1.0, 0.0]   # [real, imag] in nm³
/// chi_zxx = [0.3, 0.0]
/// far_field = false
///
/// enable_thg = true
/// chi3_symmetry = "isotropic_bulk"
/// chi3_xxxx = [1.0, 0.0]   # [real, imag] in nm⁶
/// chi3_xxyy = [0.25, 0.0]
/// ```
#[derive(Debug, Deserialize, Default)]
pub struct NonlinearConfig {
    /// Enable second-harmonic generation (SHG) calculation. Default: false.
    #[serde(default)]
    pub enable_shg: bool,
    /// χ^(2) symmetry class applied to every dipole.
    /// Valid values: `"isotropic_surface"` (C_∞v, surface normal along **z**),
    /// `"zero"` (centrosymmetric, default).
    #[serde(default = "default_shg_symmetry")]
    pub symmetry: String,
    /// χ_zzz component `[real, imag]` in nm³. Required for `isotropic_surface`.
    pub chi_zzz: Option<[f64; 2]>,
    /// χ_zxx component `[real, imag]` in nm³. Required for `isotropic_surface`.
    pub chi_zxx: Option<[f64; 2]>,
    /// Compute far-field radiation pattern at $2\omega$. Default: false.
    #[serde(default)]
    pub far_field: bool,
    /// Enable third-harmonic generation (THG) calculation. Default: false.
    #[serde(default)]
    pub enable_thg: bool,
    /// χ^(3) symmetry class applied to every dipole.
    /// Valid values: `"isotropic_bulk"` (Kleinman-symmetric isotropic bulk),
    /// `"zero"` (default).
    #[serde(default = "default_chi3_symmetry")]
    pub chi3_symmetry: String,
    /// χ_xxxx component `[real, imag]` in nm⁶. Required for `isotropic_bulk`.
    pub chi3_xxxx: Option<[f64; 2]>,
    /// χ_xxyy component `[real, imag]` in nm⁶. Required for `isotropic_bulk`.
    pub chi3_xxyy: Option<[f64; 2]>,
}

fn default_shg_symmetry() -> String {
    "zero".into()
}

fn default_chi3_symmetry() -> String {
    "zero".into()
}

/// Load and parse a TOML job configuration file.
pub fn load_config(path: &std::path::Path) -> anyhow::Result<JobConfig> {
    let content = std::fs::read_to_string(path)?;
    let config: JobConfig = toml::from_str(&content)?;
    Ok(config)
}
