//! Simulation runner: ties together geometry, materials, and solver.

use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use num_complex::Complex64;

use lumina_compute::ComputeBackend;
use lumina_core::solver::cda::CdaSolver;
use lumina_core::solver::{NearFieldPlane, OpticalSolver};
use lumina_core::types::{
    clausius_mossotti, radiative_correction, CrossSections, Dipole, NearFieldMap, SimulationParams,
};
use lumina_geometry::discretise::discretise_primitive;
use lumina_geometry::primitives::{
    Cuboid, Cylinder, Ellipsoid, Helix, Primitive, Sphere,
};
use lumina_materials::johnson_christy::JohnsonChristyMaterial;
use lumina_materials::palik::PalikMaterial;
use lumina_materials::provider::MaterialProvider;

use crate::config::{JobConfig, ShapeConfig, WavelengthSpec};

/// Results from a simulation run.
pub struct SimulationOutput {
    pub spectra: Vec<CrossSections>,
    pub near_field: Option<NearFieldMap>,
}

/// Run a full simulation from a parsed job configuration.
pub fn run_simulation(job: &JobConfig) -> Result<SimulationOutput> {
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

    // Load all material providers
    let gold   = JohnsonChristyMaterial::gold();
    let silver = JohnsonChristyMaterial::silver();
    let copper = JohnsonChristyMaterial::copper();
    let tio2   = PalikMaterial::tio2();
    let sio2   = PalikMaterial::sio2();

    // Select compute backend based on config.
    let backend: Arc<dyn ComputeBackend> = create_backend(&job.simulation.backend);

    // Run simulation across wavelengths
    let solver = CdaSolver {
        backend,
        ..Default::default()
    };
    let mut all_spectra = Vec::with_capacity(wavelengths.len());

    // Track the peak extinction wavelength for near-field computation
    let mut peak_ext = 0.0_f64;
    let mut peak_wl = wavelengths[0];
    let mut peak_dipoles: Vec<Dipole> = Vec::new();

    for (wi, &wl) in wavelengths.iter().enumerate() {
        let k = 2.0 * std::f64::consts::PI * params.environment_n / wl;
        let epsilon_m = params.environment_n * params.environment_n;

        let mut dipoles = Vec::with_capacity(total_dipoles);
        for i in 0..total_dipoles {
            let epsilon = resolve_material_epsilon(
                &all_materials[i], wl, &gold, &silver, &copper, &tio2, &sio2,
            )?;
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

        if cs.extinction > peak_ext {
            peak_ext = cs.extinction;
            peak_wl = wl;
            peak_dipoles = dipoles;
        }

        all_spectra.push(cs);
    }

    // Optionally compute near-field map at peak wavelength
    let near_field = if job.output.save_near_field && !peak_dipoles.is_empty() {
        println!("Computing near-field map at λ={:.1} nm (peak extinction)...", peak_wl);
        match solver.solve_dipoles(&peak_dipoles, peak_wl, &params) {
            Ok(response) => {
                // Estimate bounding radius of the structure for the observation plane
                let half_extent = all_positions
                    .iter()
                    .map(|p| {
                        (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
                    })
                    .fold(0.0_f64, f64::max);
                let plane_half = (half_extent * 2.0).max(20.0);

                let plane = NearFieldPlane {
                    centre: [0.0, 0.0, 0.0],
                    normal: [0.0, 0.0, 1.0],
                    half_width: plane_half,
                    half_height: plane_half,
                    nx: 80,
                    ny: 80,
                };
                match solver.compute_near_field(&peak_dipoles, &response, &plane) {
                    Ok(nf) => Some(nf),
                    Err(e) => {
                        eprintln!("Warning: near-field computation failed: {}", e);
                        None
                    }
                }
            }
            Err(e) => {
                eprintln!("Warning: dipole solve for near-field failed: {}", e);
                None
            }
        }
    } else {
        None
    };

    Ok(SimulationOutput { spectra: all_spectra, near_field })
}

/// Build a Primitive from the TOML shape configuration.
fn build_primitive(shape: &ShapeConfig, name: &str) -> Result<Primitive> {
    match shape {
        ShapeConfig::Primitive { shape_type, params } => match shape_type.as_str() {
            "sphere" => {
                let centre = extract_f64_array(params, "centre").unwrap_or([0.0, 0.0, 0.0]);
                let radius = params
                    .get("radius")
                    .and_then(|v| v.as_float())
                    .context(format!("Object '{}': sphere requires 'radius'", name))?;
                Ok(Primitive::Sphere(Sphere { centre, radius }))
            }
            "cylinder" => {
                let base_centre = extract_f64_array(params, "base_centre")
                    .unwrap_or([0.0, 0.0, 0.0]);
                let axis = extract_f64_array(params, "axis").unwrap_or([0.0, 0.0, 1.0]);
                let radius = params
                    .get("radius")
                    .and_then(|v| v.as_float())
                    .context(format!("Object '{}': cylinder requires 'radius'", name))?;
                let length = params
                    .get("length")
                    .and_then(|v| v.as_float())
                    .context(format!("Object '{}': cylinder requires 'length'", name))?;
                Ok(Primitive::Cylinder(Cylinder { base_centre, axis, radius, length }))
            }
            "cuboid" => {
                let centre = extract_f64_array(params, "centre").unwrap_or([0.0, 0.0, 0.0]);
                let half_extents = extract_f64_array(params, "half_extents")
                    .context(format!("Object '{}': cuboid requires 'half_extents = [hx, hy, hz]'", name))?;
                Ok(Primitive::Cuboid(Cuboid { centre, half_extents }))
            }
            "ellipsoid" => {
                let centre = extract_f64_array(params, "centre").unwrap_or([0.0, 0.0, 0.0]);
                let semi_axes = extract_f64_array(params, "semi_axes")
                    .context(format!("Object '{}': ellipsoid requires 'semi_axes = [a, b, c]'", name))?;
                Ok(Primitive::Ellipsoid(Ellipsoid { centre, semi_axes }))
            }
            "helix" => {
                let base_centre = extract_f64_array(params, "base_centre")
                    .unwrap_or([0.0, 0.0, 0.0]);
                let axis = extract_f64_array(params, "axis").unwrap_or([0.0, 0.0, 1.0]);
                let radius = params
                    .get("radius")
                    .and_then(|v| v.as_float())
                    .context(format!("Object '{}': helix requires 'radius'", name))?;
                let pitch = params
                    .get("pitch")
                    .and_then(|v| v.as_float())
                    .context(format!("Object '{}': helix requires 'pitch'", name))?;
                let turns = params
                    .get("turns")
                    .and_then(|v| v.as_float())
                    .context(format!("Object '{}': helix requires 'turns'", name))?;
                let wire_radius = params
                    .get("wire_radius")
                    .and_then(|v| v.as_float())
                    .context(format!("Object '{}': helix requires 'wire_radius'", name))?;
                Ok(Primitive::Helix(Helix {
                    base_centre,
                    axis,
                    radius,
                    pitch,
                    turns,
                    wire_radius,
                }))
            }
            other => anyhow::bail!(
                "Unsupported shape type '{}' for object '{}'. Valid types: sphere, cylinder, cuboid, ellipsoid, helix",
                other, name
            ),
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
    silver: &JohnsonChristyMaterial,
    copper: &JohnsonChristyMaterial,
    tio2: &PalikMaterial,
    sio2: &PalikMaterial,
) -> Result<Complex64> {
    let provider: &dyn MaterialProvider = match material_id {
        "Au_JC"       => gold,
        "Ag_JC"       => silver,
        "Cu_JC"       => copper,
        "TiO2_Palik"  => tio2,
        "SiO2_Palik"  => sio2,
        _ => anyhow::bail!(
            "Unknown material '{}'. Valid identifiers: Au_JC, Ag_JC, Cu_JC, TiO2_Palik, SiO2_Palik",
            material_id
        ),
    };
    provider
        .dielectric_function(wavelength_nm)
        .map_err(|e| anyhow::anyhow!("Material '{}' at {:.1} nm: {}", material_id, wavelength_nm, e))
}

/// Write cross-section spectra to a CSV file with a metadata header.
pub fn write_spectra_csv(
    spectra: &[CrossSections],
    path: &Path,
    job: &JobConfig,
) -> Result<()> {
    use std::io::Write;

    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut file = std::fs::File::create(path)?;

    // Metadata header
    writeln!(file, "# Lumina CDA Solver — Cross-Section Spectra")?;
    writeln!(file, "# Version: {}", env!("CARGO_PKG_VERSION"))?;
    writeln!(file, "# environment_n: {}", job.simulation.environment_n)?;

    // List object materials and spacings
    for obj in &job.geometry.object {
        writeln!(
            file,
            "# object '{}': material={}, dipole_spacing={} nm",
            obj.name, obj.material, obj.dipole_spacing
        )?;
    }
    writeln!(file, "#")?;

    // Add CD column only if present
    let has_cd = spectra.iter().any(|cs| cs.circular_dichroism.is_some());
    if has_cd {
        writeln!(file, "wavelength_nm,extinction_nm2,absorption_nm2,scattering_nm2,circular_dichroism_nm2")?;
    } else {
        writeln!(file, "wavelength_nm,extinction_nm2,absorption_nm2,scattering_nm2")?;
    }

    for cs in spectra {
        if has_cd {
            writeln!(
                file,
                "{:.2},{:.6e},{:.6e},{:.6e},{:.6e}",
                cs.wavelength_nm,
                cs.extinction,
                cs.absorption,
                cs.scattering,
                cs.circular_dichroism.unwrap_or(0.0),
            )?;
        } else {
            writeln!(
                file,
                "{:.2},{:.6e},{:.6e},{:.6e}",
                cs.wavelength_nm, cs.extinction, cs.absorption, cs.scattering
            )?;
        }
    }

    println!("Spectra written to: {}", path.display());
    Ok(())
}

/// Write cross-section spectra to a JSON file.
pub fn write_spectra_json(
    spectra: &[CrossSections],
    path: &Path,
) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let json = serde_json::to_string_pretty(spectra)
        .map_err(|e| anyhow::anyhow!("JSON serialisation error: {}", e))?;
    std::fs::write(path, json)?;

    println!("Spectra (JSON) written to: {}", path.display());
    Ok(())
}

/// Write a near-field map to a CSV file.
pub fn write_near_field_csv(
    near_field: &NearFieldMap,
    path: &Path,
) -> Result<()> {
    use std::io::Write;

    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut file = std::fs::File::create(path)?;
    writeln!(file, "# Lumina CDA Solver — Near-Field Intensity Map")?;
    writeln!(file, "# Grid: {}x{}", near_field.nx, near_field.ny)?;
    writeln!(
        file,
        "# Extent: x=[{:.2}, {:.2}] y=[{:.2}, {:.2}] nm",
        near_field.extent[0], near_field.extent[1],
        near_field.extent[2], near_field.extent[3],
    )?;
    writeln!(file, "#")?;
    writeln!(file, "x_nm,y_nm,z_nm,|E|_sq")?;

    for (pos, intensity) in near_field.positions.iter().zip(near_field.field_intensity.iter()) {
        writeln!(
            file,
            "{:.4},{:.4},{:.4},{:.6e}",
            pos[0], pos[1], pos[2], intensity
        )?;
    }

    println!("Near-field map written to: {}", path.display());
    Ok(())
}

/// Create a compute backend based on the user's preference string.
///
/// - `"gpu"` — attempt GPU, fail if unavailable.
/// - `"cpu"` — always use CPU.
/// - `"auto"` (default) — try GPU, fall back to CPU.
fn create_backend(preference: &str) -> Arc<dyn ComputeBackend> {
    match preference {
        "cpu" => {
            println!("Backend: CPU");
            Arc::new(lumina_compute::CpuBackend::new())
        }
        "gpu" => {
            #[cfg(feature = "gpu")]
            {
                match lumina_compute::GpuBackend::new_blocking() {
                    Ok(gpu) => {
                        println!("Backend: {}", gpu.device_info().name);
                        Arc::new(gpu)
                    }
                    Err(e) => {
                        eprintln!("GPU requested but unavailable: {}. Aborting.", e);
                        std::process::exit(1);
                    }
                }
            }
            #[cfg(not(feature = "gpu"))]
            {
                eprintln!("GPU requested but binary was built without --features gpu. Aborting.");
                std::process::exit(1);
            }
        }
        _ => {
            // "auto" or any unrecognised value
            #[cfg(feature = "gpu")]
            {
                match lumina_compute::GpuBackend::new_blocking() {
                    Ok(gpu) => {
                        println!("Backend: {} (auto-detected)", gpu.device_info().name);
                        return Arc::new(gpu);
                    }
                    Err(e) => {
                        println!("GPU not available ({}), using CPU", e);
                    }
                }
            }
            Arc::new(lumina_compute::CpuBackend::new())
        }
    }
}
