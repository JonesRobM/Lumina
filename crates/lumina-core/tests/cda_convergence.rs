//! Test CDA convergence: compare RRCM vs LDR for metallic and dielectric particles.

use lumina_core::mie::mie_cross_sections;
use lumina_core::solver::cda::CdaSolver;
use lumina_core::solver::OpticalSolver;
use lumina_core::types::{clausius_mossotti, ldr_correction, Dipole, SimulationParams};
use lumina_geometry::discretise::discretise_primitive;
use lumina_geometry::primitives::{Primitive, Sphere};
use num_complex::Complex64;

mod common;

fn run_convergence(
    label: &str,
    radius: f64,
    wavelength: f64,
    epsilon: Complex64,
    n_medium: f64,
    use_ldr: bool,
) {
    let epsilon_m = n_medium * n_medium;
    let k = 2.0 * std::f64::consts::PI * n_medium / wavelength;
    let propagation = [0.0, 0.0, 1.0];
    let polarisation = [1.0, 0.0, 0.0];

    let (mie_ext, _, _) = mie_cross_sections(radius, wavelength, epsilon, n_medium);

    // Each dipole carries its own FCD cell size, so a spacing sweep is measured
    // against the matching cell automatically (debt item C5).
    let params = SimulationParams {
        environment_n: n_medium,
        ..Default::default()
    };

    let method = if use_ldr { "LDR" } else { "RRCM" };
    eprintln!("\n=== {} [{}]: eps=({:.1},{:.1}), R={}nm, lambda={}nm ===",
        label, method, epsilon.re, epsilon.im, radius, wavelength);
    eprintln!("Mie ext = {:.4e}", mie_ext);
    eprintln!("{:>6} {:>5} {:>11} {:>11} {:>7} {:>11} {:>11}",
        "d", "N", "CDA_ext", "CDA_abs", "err%", "CDA_sca", "sca_frac%");

    for &spacing in &[4.0, 3.0, 2.5, 2.0, 1.5] {
        let sphere = Primitive::Sphere(Sphere {
            centre: [0.0, 0.0, 0.0],
            radius,
        });
        let lattice = discretise_primitive(&sphere, spacing);
        let positions: Vec<[f64; 3]> = lattice.iter().map(|p| p.position).collect();
        let fill_fractions: Vec<f64> = lattice.iter().map(|p| p.fill_fraction).collect();

        if positions.is_empty() {
            continue;
        }

        let solver = CdaSolver::default();

        let dipoles: Vec<Dipole> = if use_ldr {
            // The LDR correction is not offered by `lumina_core::build_dipoles`, which
            // is RRCM-only, so this branch necessarily builds its own dipoles. It is the
            // one deliberate exception to routing tests through the production path --
            // keep it in step if the polarisability pipeline changes.
            positions
                .iter()
                .zip(&fill_fractions)
                .map(|(&pos, &f)| {
                    let alpha_cm = clausius_mossotti(spacing.powi(3), epsilon, epsilon_m);
                    let alpha =
                        ldr_correction(alpha_cm, k, spacing, epsilon, &propagation, &polarisation);
                    Dipole::isotropic(pos, alpha * f, spacing)
                })
                .collect()
        } else {
            common::dipoles_with_constant_eps(&lattice, spacing, wavelength, n_medium, epsilon)
        };

        let cda = solver
            .compute_cross_sections(&dipoles, wavelength, &params)
            .unwrap();
        let err = (cda.extinction - mie_ext).abs() / mie_ext;
        let sca_frac = if cda.extinction > 0.0 { cda.scattering / cda.extinction * 100.0 } else { 0.0 };

        eprintln!(
            "{:6.1} {:5} {:11.4e} {:11.4e} {:6.1} {:11.4e} {:10.1}",
            spacing, dipoles.len(), cda.extinction, cda.absorption,
            err * 100.0, cda.scattering, sca_frac,
        );
    }
}

#[test]
fn test_convergence_comparison() {
    let eps_dielectric = Complex64::new(4.0, 0.5);
    let eps_metallic = Complex64::new(-10.0, 1.5);

    // Dielectric: both methods should work
    run_convergence("Dielectric", 8.0, 600.0, eps_dielectric, 1.0, false);
    run_convergence("Dielectric", 8.0, 600.0, eps_dielectric, 1.0, true);

    // Metallic: compare RRCM vs LDR
    run_convergence("Metallic", 8.0, 600.0, eps_metallic, 1.0, false);
    run_convergence("Metallic", 8.0, 600.0, eps_metallic, 1.0, true);
}
