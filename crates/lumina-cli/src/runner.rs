//! Simulation runner: ties together geometry, materials, and solver.

use std::path::Path;

use anyhow::{Context, Result};
use num_complex::Complex64;

use lumina_core::solver::cda::CdaSolver;
use lumina_core::solver::OpticalSolver;
use lumina_core::types::{
    clausius_mossotti, radiative_correction, CrossSections, Dipole, SimulationParams,
};
use lumina_geometry::discretise::discretise_primitive;
use lumina_geometry::primitives::{Primitive, Sphere};
use lumina_materials::johnson_christy::JohnsonChristyMaterial;
use lumina_materials::provider::MaterialProvider;

use crate::config::{JobConfig, ShapeConfig, WavelengthSpec};

/// Run a full simulation from a parsed job configuration.
pub fn run_simulation(job: &JobConfig) -> Result<Vec<CrossSections>> {
    // Build wavelength grid
    let wavelengths = match &job.simulation.wavelengths {
        WavelengthSpec::Range { range, points } => {
            let start = range[0];
            let end = range[1];
            (0..*points)
                .map(|i| start + (end - start) * i as f64 / (*points - 1).max(1) as f64)
                .collect::<Vec<_>>()
        }
        WavelengthSpec::List { values } => values.clone(),
    };

    let params = SimulationParams {
        wavelength_range: [
            *wavelengths.first().unwrap_or(&400.0),
            *wavelengths.last().unwrap_or(&900.0),
        ],
        num_wavelengths: wavelengths.len(),
        environment_n: job.simulation.environment_n,
        solver_tolerance: job.simulation.solver_tolerance,
        max_iterations: job.simulation.max_iterations,
    };

    // Build dipoles from all geometry objects
    let mut all_positions: Vec<[f64; 3]> = Vec::new();
    let mut all_materials: Vec<String> = Vec::new();
    let mut all_spacings: Vec<f64> = Vec::new();

    for obj in &job.geometry.object {
        let primitive = build_primitive(&obj.shape, &obj.name)?;
        let lattice = discretise_primitive(&primitive, obj.dipole_spacing);

        println!(
            "  Object '{}': {} dipoles (spacing={} nm, material={})",
            obj.name,
            lattice.len(),
            obj.dipole_spacing,
            obj.material
        );

        for point in &lattice {
            all_positions.push(point.position);
            all_materials.push(obj.material.clone());
            all_spacings.push(obj.dipole_spacing);
        }
    }

    let total_dipoles = all_positions.len();
    println!("Total dipoles: {}", total_dipoles);
    if total_dipoles == 0 {
        anyhow::bail!("No dipoles generated — check geometry configuration");
    }

    // Load material providers
    let gold = JohnsonChristyMaterial::gold();

    // Run simulation across wavelengths
    let solver = CdaSolver::default();
    let mut all_spectra = Vec::with_capacity(wavelengths.len());

    for (wi, &wl) in wavelengths.iter().enumerate() {
        // Build dipoles with wavelength-dependent polarisability
        let k = 2.0 * std::f64::consts::PI * params.environment_n / wl;
        let epsilon_m = params.environment_n * params.environment_n;

        let mut dipoles = Vec::with_capacity(total_dipoles);
        for i in 0..total_dipoles {
            let epsilon = resolve_material_epsilon(&all_materials[i], wl, &gold)?;
            let volume = all_spacings[i].powi(3);
            let alpha_cm = clausius_mossotti(volume, epsilon, epsilon_m);
            let alpha = radiative_correction(alpha_cm, k);
            dipoles.push(Dipole::isotropic(all_positions[i], alpha));
        }

        let cs = solver
            .compute_cross_sections(&dipoles, wl, &params)
            .map_err(|e| anyhow::anyhow!("Solver error at λ={:.1} nm: {}", wl, e))?;

        if (wi + 1) % 10 == 0 || wi == 0 || wi == wavelengths.len() - 1 {
            println!(
                "  [{}/{}] λ={:.1} nm: C_ext={:.2e}, C_abs={:.2e}, C_sca={:.2e}",
                wi + 1,
                wavelengths.len(),
                wl,
                cs.extinction,
                cs.absorption,
                cs.scattering
            );
        }

        all_spectra.push(cs);
    }

    Ok(all_spectra)
}

/// Build a Primitive from the TOML shape configuration.
fn build_primitive(shape: &ShapeConfig, name: &str) -> Result<Primitive> {
    match shape {
        ShapeConfig::Primitive { shape_type, params } => match shape_type.as_str() {
            "sphere" => {
                let centre = extract_f64_array(params, "centre")
                    .unwrap_or([0.0, 0.0, 0.0]);
                let radius = params
                    .get("radius")
                    .and_then(|v| v.as_float())
                    .context(format!("Object '{}': sphere requires 'radius'", name))?;
                Ok(Primitive::Sphere(Sphere { centre, radius }))
            }
            other => anyhow::bail!("Unsupported shape type '{}' for object '{}'", other, name),
        },
        ShapeConfig::File { geometry_file } => {
            anyhow::bail!(
                "File-based geometry ('{}') not yet implemented for object '{}'",
                geometry_file,
                name
            )
        }
    }
}

fn extract_f64_array(params: &toml::Value, key: &str) -> Option<[f64; 3]> {
    let arr = params.get(key)?.as_array()?;
    if arr.len() != 3 {
        return None;
    }
    Some([
        arr[0].as_float()?,
        arr[1].as_float()?,
        arr[2].as_float()?,
    ])
}

/// Resolve the dielectric function for a material identifier at a given wavelength.
fn resolve_material_epsilon(
    material_id: &str,
    wavelength_nm: f64,
    gold: &JohnsonChristyMaterial,
) -> Result<Complex64> {
    match material_id {
        "Au_JC" => gold
            .dielectric_function(wavelength_nm)
            .map_err(|e| anyhow::anyhow!("Material error: {}", e)),
        _ => anyhow::bail!("Unknown material: '{}'", material_id),
    }
}

/// Write cross-section spectra to a CSV file.
pub fn write_spectra_csv(spectra: &[CrossSections], path: &Path) -> Result<()> {
    use std::io::Write;

    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut file = std::fs::File::create(path)?;
    writeln!(file, "wavelength_nm,extinction_nm2,absorption_nm2,scattering_nm2")?;

    for cs in spectra {
        writeln!(
            file,
            "{:.2},{:.6e},{:.6e},{:.6e}",
            cs.wavelength_nm, cs.extinction, cs.absorption, cs.scattering
        )?;
    }

    println!("Spectra written to: {}", path.display());
    Ok(())
}
