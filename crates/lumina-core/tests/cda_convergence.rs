//! Test CDA convergence vs discretisation for metallic particles.

use lumina_core::mie::mie_cross_sections;
use lumina_core::solver::cda::CdaSolver;
use lumina_core::solver::OpticalSolver;
use lumina_core::types::{clausius_mossotti, radiative_correction, Dipole, SimulationParams};
use lumina_geometry::discretise::discretise_primitive;
use lumina_geometry::primitives::{Primitive, Sphere};
use num_complex::Complex64;

#[test]
fn test_convergence_metallic() {
    let radius = 8.0;
    let wavelength = 600.0;
    let n_medium = 1.0;
    let epsilon = Complex64::new(-10.0, 1.5); // highly metallic
    let epsilon_m = n_medium * n_medium;
    let k = 2.0 * std::f64::consts::PI * n_medium / wavelength;

    let (mie_ext, _, _) = mie_cross_sections(radius, wavelength, epsilon, n_medium);

    let solver = CdaSolver::default();
    let params = SimulationParams {
        environment_n: n_medium,
        ..Default::default()
    };

    eprintln!("=== Convergence test: ε=(-10, 1.5), R=8nm, λ=600nm ===");
    eprintln!("Mie extinction = {:.4e}", mie_ext);
    eprintln!("{:>8} {:>6} {:>12} {:>12} {:>8}", "spacing", "N", "CDA_ext", "Mie_ext", "err%");

    for &spacing in &[4.0, 3.0, 2.5, 2.0, 1.5] {
        let sphere = Primitive::Sphere(Sphere {
            centre: [0.0, 0.0, 0.0],
            radius,
        });
        let lattice = discretise_primitive(&sphere, spacing);
        let positions: Vec<[f64; 3]> = lattice.iter().map(|p| p.position).collect();

        if positions.is_empty() {
            continue;
        }

        let dipoles: Vec<Dipole> = positions
            .iter()
            .map(|&pos| {
                let volume = spacing.powi(3);
                let alpha_cm = clausius_mossotti(volume, epsilon, epsilon_m);
                let alpha = radiative_correction(alpha_cm, k);
                Dipole::isotropic(pos, alpha)
            })
            .collect();

        let cda = solver
            .compute_cross_sections(&dipoles, wavelength, &params)
            .unwrap();
        let err = (cda.extinction - mie_ext).abs() / mie_ext;

        eprintln!(
            "{:8.1} {:6} {:12.4e} {:12.4e} {:7.1}",
            spacing,
            dipoles.len(),
            cda.extinction,
            mie_ext,
            err * 100.0,
        );
    }
}

#[test]
fn test_convergence_dielectric() {
    let radius = 8.0;
    let wavelength = 600.0;
    let n_medium = 1.0;
    let epsilon = Complex64::new(4.0, 0.5); // lossy dielectric
    let epsilon_m = n_medium * n_medium;
    let k = 2.0 * std::f64::consts::PI * n_medium / wavelength;

    let (mie_ext, _, _) = mie_cross_sections(radius, wavelength, epsilon, n_medium);

    let solver = CdaSolver::default();
    let params = SimulationParams {
        environment_n: n_medium,
        ..Default::default()
    };

    eprintln!("\n=== Convergence test: ε=(4, 0.5), R=8nm, λ=600nm ===");
    eprintln!("Mie extinction = {:.4e}", mie_ext);
    eprintln!("{:>8} {:>6} {:>12} {:>12} {:>8}", "spacing", "N", "CDA_ext", "Mie_ext", "err%");

    for &spacing in &[4.0, 3.0, 2.5, 2.0, 1.5] {
        let sphere = Primitive::Sphere(Sphere {
            centre: [0.0, 0.0, 0.0],
            radius,
        });
        let lattice = discretise_primitive(&sphere, spacing);
        let positions: Vec<[f64; 3]> = lattice.iter().map(|p| p.position).collect();

        if positions.is_empty() {
            continue;
        }

        let dipoles: Vec<Dipole> = positions
            .iter()
            .map(|&pos| {
                let volume = spacing.powi(3);
                let alpha_cm = clausius_mossotti(volume, epsilon, epsilon_m);
                let alpha = radiative_correction(alpha_cm, k);
                Dipole::isotropic(pos, alpha)
            })
            .collect();

        let cda = solver
            .compute_cross_sections(&dipoles, wavelength, &params)
            .unwrap();
        let err = (cda.extinction - mie_ext).abs() / mie_ext;

        eprintln!(
            "{:8.1} {:6} {:12.4e} {:12.4e} {:7.1}",
            spacing,
            dipoles.len(),
            cda.extinction,
            mie_ext,
            err * 100.0,
        );
    }
}
