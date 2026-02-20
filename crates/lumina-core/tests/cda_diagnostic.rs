//! Diagnostic tests to isolate CDA vs Mie discrepancies.

use lumina_core::mie::mie_cross_sections;
use lumina_core::solver::cda::CdaSolver;
use lumina_core::solver::OpticalSolver;
use lumina_core::types::{clausius_mossotti, radiative_correction, Dipole, SimulationParams};
use lumina_geometry::discretise::discretise_primitive;
use lumina_geometry::primitives::{Primitive, Sphere};
use num_complex::Complex64;

/// Helper: run CDA for a sphere with given epsilon and compare to Mie.
fn cda_vs_mie(
    radius: f64,
    spacing: f64,
    wavelength: f64,
    epsilon: Complex64,
    n_medium: f64,
) -> (f64, f64, f64) {
    let k = 2.0 * std::f64::consts::PI * n_medium / wavelength;
    let epsilon_m = n_medium * n_medium;

    let sphere = Primitive::Sphere(Sphere {
        centre: [0.0, 0.0, 0.0],
        radius,
    });
    let lattice = discretise_primitive(&sphere, spacing);
    let positions: Vec<[f64; 3]> = lattice.iter().map(|p| p.position).collect();

    let dipoles: Vec<Dipole> = positions
        .iter()
        .map(|&pos| {
            let volume = spacing.powi(3);
            let alpha_cm = clausius_mossotti(volume, epsilon, epsilon_m);
            let alpha = radiative_correction(alpha_cm, k);
            Dipole::isotropic(pos, alpha)
        })
        .collect();

    let solver = CdaSolver::default();
    let params = SimulationParams {
        environment_n: n_medium,
        ..Default::default()
    };

    let cda = solver
        .compute_cross_sections(&dipoles, wavelength, &params)
        .unwrap();
    let (mie_ext, _, _) = mie_cross_sections(radius, wavelength, epsilon, n_medium);

    let rel_err = (cda.extinction - mie_ext).abs() / mie_ext;

    eprintln!(
        "  ε=({:.2},{:.2}), N={}, CDA={:.4e}, Mie={:.4e}, err={:.1}%",
        epsilon.re,
        epsilon.im,
        dipoles.len(),
        cda.extinction,
        mie_ext,
        rel_err * 100.0,
    );

    (cda.extinction, mie_ext, rel_err)
}

/// Test CDA for a DIELECTRIC sphere (should converge well).
#[test]
fn test_cda_dielectric_sphere() {
    let radius = 10.0;
    let spacing = 3.0;
    let wavelength = 600.0;
    let n_medium = 1.0;

    eprintln!("=== Dielectric sphere (lossless) ===");
    let epsilon = Complex64::new(4.0, 0.0); // n=2
    let (_, _, err) = cda_vs_mie(radius, spacing, wavelength, epsilon, n_medium);
    assert!(
        err < 0.5,
        "Dielectric sphere error ({:.1}%) should be < 50%",
        err * 100.0
    );

    eprintln!("\n=== Dielectric sphere (lossy) ===");
    let epsilon = Complex64::new(4.0, 0.5);
    let (_, _, err) = cda_vs_mie(radius, spacing, wavelength, epsilon, n_medium);
    assert!(
        err < 0.5,
        "Lossy dielectric error ({:.1}%) should be < 50%",
        err * 100.0
    );
}

/// Test CDA for increasingly metallic particles to find the breakdown.
#[test]
fn test_cda_metallic_progression() {
    let radius = 10.0;
    let spacing = 3.0;
    let wavelength = 600.0;
    let n_medium = 1.0;

    let test_cases: Vec<(Complex64, &str)> = vec![
        (Complex64::new(4.0, 0.0), "dielectric n=2"),
        (Complex64::new(2.0, 1.0), "weakly absorbing"),
        (Complex64::new(-1.0, 3.0), "near resonance"),
        (Complex64::new(-2.0, 5.0), "Au-like 450nm"),
        (Complex64::new(-5.0, 2.0), "Au-like 530nm"),
        (Complex64::new(-10.0, 1.5), "Au-like 600nm"),
        (Complex64::new(-17.0, 1.0), "Au-like 700nm"),
    ];

    eprintln!("=== Metallic progression: R={}, d={}, λ={} ===", radius, spacing, wavelength);
    for (eps, label) in &test_cases {
        eprint!("{:20}: ", label);
        cda_vs_mie(radius, spacing, wavelength, *eps, n_medium);
    }
}

/// Test CDA for a single dipole (no interactions) — should match Rayleigh exactly.
#[test]
fn test_cda_single_dipole() {
    let wavelength = 600.0;
    let n_medium = 1.0;
    let k = 2.0 * std::f64::consts::PI * n_medium / wavelength;
    let epsilon = Complex64::new(-10.0, 1.5);
    let epsilon_m = n_medium * n_medium;

    // Single dipole representing the whole sphere
    let radius: f64 = 10.0;
    let volume = 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);
    let alpha_cm = clausius_mossotti(volume, epsilon, epsilon_m);
    let alpha = radiative_correction(alpha_cm, k);

    let dipole = Dipole::isotropic([0.0, 0.0, 0.0], alpha);
    let solver = CdaSolver::default();
    let params = SimulationParams {
        environment_n: n_medium,
        ..Default::default()
    };

    let cda = solver
        .compute_cross_sections(&[dipole], wavelength, &params)
        .unwrap();

    // Rayleigh result
    let c_ext_rayleigh = k * alpha.im;

    // Mie result
    let (mie_ext, _, _) = mie_cross_sections(radius, wavelength, epsilon, n_medium);

    eprintln!("=== Single dipole test (ε={:.1}+{:.1}i) ===", epsilon.re, epsilon.im);
    eprintln!("  CDA (1 dipole): {:.4e}", cda.extinction);
    eprintln!("  Rayleigh:       {:.4e}", c_ext_rayleigh);
    eprintln!("  Mie:            {:.4e}", mie_ext);
    eprintln!(
        "  CDA vs Rayleigh err: {:.2}%",
        (cda.extinction - c_ext_rayleigh).abs() / c_ext_rayleigh * 100.0
    );
    eprintln!(
        "  CDA vs Mie err:      {:.2}%",
        (cda.extinction - mie_ext).abs() / mie_ext * 100.0
    );

    // Single dipole CDA should match Rayleigh exactly
    assert!(
        (cda.extinction - c_ext_rayleigh).abs() / c_ext_rayleigh < 0.01,
        "Single dipole CDA must match Rayleigh to within 1%"
    );
}
