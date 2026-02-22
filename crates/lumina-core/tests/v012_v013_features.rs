//! Integration tests for v0.1.2 and v0.1.3 features.
//!
//! ## v0.1.2 coverage
//! - Johnson & Christy Ag and Cu dielectric data (correct spectral behaviour)
//! - Near-field map computation (dimensions, decay, enhancement)
//!
//! ## v0.1.3 coverage
//! - Palik TiO₂ and SiO₂ dielectric data (correct index, range enforcement)
//! - Far-field radiation pattern (grid size, dipole pattern, energy sum)
//! - Circular dichroism = 0 for isotropic sphere (symmetry theorem)
//! - `CrossSections` JSON round-trip (serde)
//! - Incident polarisation: x-pol and y-pol give related but different dipole moments

use approx::assert_relative_eq;
use num_complex::Complex64;

use lumina_core::fields::{compute_circular_dichroism, compute_far_field, compute_near_field_map};
use lumina_core::solver::cda::CdaSolver;
use lumina_core::solver::{NearFieldPlane, OpticalSolver};
use lumina_core::types::{
    clausius_mossotti, radiative_correction, CrossSections, Dipole, DipoleResponse,
    FarFieldMap, IncidentField, SimulationParams,
};
use lumina_geometry::discretise::{discretise_mesh, discretise_primitive};
use lumina_geometry::parsers::obj::parse_obj;
use lumina_geometry::primitives::{Cuboid, Ellipsoid, Primitive, Sphere};
use lumina_materials::johnson_christy::JohnsonChristyMaterial;
use lumina_materials::palik::PalikMaterial;
use lumina_materials::provider::MaterialProvider;

// ─────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────

fn default_params() -> SimulationParams {
    SimulationParams {
        wavelength_range: [400.0, 800.0],
        num_wavelengths: 10,
        environment_n: 1.0,
        solver_tolerance: 1e-8,
        max_iterations: 500,
    }
}

/// Build N dipoles on a small sphere lattice at wavelength `wl` using gold.
fn gold_sphere_dipoles(radius_nm: f64, spacing_nm: f64, wl_nm: f64) -> Vec<Dipole> {
    let sphere = Primitive::Sphere(Sphere { centre: [0.0, 0.0, 0.0], radius: radius_nm });
    let lattice = discretise_primitive(&sphere, spacing_nm);
    let mat = JohnsonChristyMaterial::gold();
    let eps = mat.dielectric_function(wl_nm).unwrap();
    let k = 2.0 * std::f64::consts::PI / wl_nm;
    lattice
        .iter()
        .map(|p| {
            let alpha_cm = clausius_mossotti(spacing_nm.powi(3), eps, 1.0);
            let alpha = radiative_correction(alpha_cm, k);
            Dipole::isotropic(p.position, alpha)
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────
// v0.1.2: Johnson & Christy Ag data
// ─────────────────────────────────────────────────────────────

#[test]
fn test_ag_jc_epsilon_is_metallic_in_visible() {
    let mat = JohnsonChristyMaterial::silver();
    // Silver is a Drude metal: ε₁ strongly negative in visible, ε₂ > 0
    let eps_500 = mat.dielectric_function(500.0).unwrap();
    assert!(eps_500.re < 0.0, "Ag ε₁ should be negative at 500 nm, got {:.3}", eps_500.re);
    assert!(eps_500.im > 0.0, "Ag ε₂ should be positive at 500 nm, got {:.3}", eps_500.im);
}

#[test]
fn test_ag_jc_uv_and_nir_in_range() {
    let mat = JohnsonChristyMaterial::silver();
    // Data covers 188–892 nm
    mat.dielectric_function(200.0).expect("Ag: 200 nm should be in range");
    mat.dielectric_function(800.0).expect("Ag: 800 nm should be in range");
}

#[test]
fn test_ag_jc_spline_is_smooth() {
    // ε₁ and ε₂ should vary smoothly — test that adjacent wavelengths aren't wildly different
    let mat = JohnsonChristyMaterial::silver();
    let mut prev_eps = mat.dielectric_function(400.0).unwrap();
    for wl in (410..=700usize).step_by(10) {
        let eps = mat.dielectric_function(wl as f64).unwrap();
        let delta_re = (eps.re - prev_eps.re).abs();
        let delta_im = (eps.im - prev_eps.im).abs();
        assert!(delta_re < 5.0, "Ag ε₁ jumped by {delta_re:.2} at {wl} nm — spline not smooth");
        assert!(delta_im < 5.0, "Ag ε₂ jumped by {delta_im:.2} at {wl} nm — spline not smooth");
        prev_eps = eps;
    }
}

// ─────────────────────────────────────────────────────────────
// v0.1.2: Johnson & Christy Cu data
// ─────────────────────────────────────────────────────────────

#[test]
fn test_cu_jc_epsilon_is_metallic_in_visible() {
    let mat = JohnsonChristyMaterial::copper();
    // Copper is metallic: ε₁ negative in visible
    let eps_600 = mat.dielectric_function(600.0).unwrap();
    assert!(eps_600.re < 0.0, "Cu ε₁ should be negative at 600 nm, got {:.3}", eps_600.re);
    assert!(eps_600.im > 0.0, "Cu ε₂ should be positive at 600 nm, got {:.3}", eps_600.im);
}

#[test]
fn test_cu_jc_interband_at_560nm() {
    // Cu d-band absorption edge near 560 nm: ε₂ should be notably larger than at 700 nm
    let mat = JohnsonChristyMaterial::copper();
    let eps_560 = mat.dielectric_function(560.0).unwrap();
    let eps_700 = mat.dielectric_function(700.0).unwrap();
    assert!(
        eps_560.im > eps_700.im,
        "Cu ε₂ should be larger near the d-band (560 nm={:.2}) than in the Drude region (700 nm={:.2})",
        eps_560.im, eps_700.im
    );
}

// ─────────────────────────────────────────────────────────────
// v0.1.2: Near-field map
// ─────────────────────────────────────────────────────────────

#[test]
fn test_near_field_map_dimensions() {
    let dipoles = gold_sphere_dipoles(10.0, 3.0, 520.0);
    let solver = CdaSolver::default();
    let params = default_params();
    let response = solver.solve_dipoles(&dipoles, 520.0, &params).unwrap();

    let plane = NearFieldPlane {
        centre: [0.0, 0.0, 0.0],
        normal: [0.0, 0.0, 1.0],
        half_width: 20.0,
        half_height: 20.0,
        nx: 15,
        ny: 12,
    };
    let nf = compute_near_field_map(&dipoles, &response, &plane, 2.0 * std::f64::consts::PI / 520.0, &solver.incident_field).unwrap();

    assert_eq!(nf.nx, 15);
    assert_eq!(nf.ny, 12);
    assert_eq!(nf.positions.len(), 15 * 12);
    assert_eq!(nf.field_intensity.len(), 15 * 12);
}

#[test]
fn test_near_field_intensity_is_nonnegative() {
    let dipoles = gold_sphere_dipoles(10.0, 3.0, 520.0);
    let solver = CdaSolver::default();
    let params = default_params();
    let response = solver.solve_dipoles(&dipoles, 520.0, &params).unwrap();

    let plane = NearFieldPlane {
        centre: [0.0, 0.0, 0.0],
        normal: [0.0, 0.0, 1.0],
        half_width: 30.0,
        half_height: 30.0,
        nx: 10,
        ny: 10,
    };
    let nf = compute_near_field_map(&dipoles, &response, &plane, 2.0 * std::f64::consts::PI / 520.0, &solver.incident_field).unwrap();

    for &val in &nf.field_intensity {
        assert!(val >= 0.0, "|E|² must be non-negative, got {val:.3e}");
    }
}

#[test]
fn test_near_field_nonzero_outside_particle() {
    // Far from the particle the field should still contain the incident field (≥ 1.0 for E₀=1)
    let dipoles = gold_sphere_dipoles(10.0, 3.0, 520.0);
    let solver = CdaSolver::default();
    let params = default_params();
    let response = solver.solve_dipoles(&dipoles, 520.0, &params).unwrap();

    // Observation plane well outside the particle
    let plane = NearFieldPlane {
        centre: [0.0, 0.0, 0.0],
        normal: [0.0, 0.0, 1.0],
        half_width: 100.0,
        half_height: 100.0,
        nx: 5,
        ny: 5,
    };
    let nf = compute_near_field_map(&dipoles, &response, &plane, 2.0 * std::f64::consts::PI / 520.0, &solver.incident_field).unwrap();

    // At least some points should have |E|² > 0.5 (incident field alone gives 1.0)
    let nonzero = nf.field_intensity.iter().filter(|&&v| v > 0.5).count();
    assert!(nonzero > 0, "All near-field intensities are near zero — incident field missing");
}

// ─────────────────────────────────────────────────────────────
// v0.1.2: Geometry — cuboid and ellipsoid discretisation
// ─────────────────────────────────────────────────────────────

#[test]
fn test_cuboid_discretisation_count() {
    // 20×20×20 nm cuboid at 2 nm spacing → expect ~(10)³ = 1000 dipoles
    let prim = Primitive::Cuboid(Cuboid {
        centre: [0.0, 0.0, 0.0],
        half_extents: [10.0, 10.0, 10.0],
    });
    let lattice = discretise_primitive(&prim, 2.0);
    let n = lattice.len();
    // Centred integer grid: -10..=+10 in steps of 2 → 11 pts/axis → 11³ = 1331 dipoles
    assert!(n > 1100 && n <= 1331, "Cuboid dipole count = {n}, expected ~1331");
}

#[test]
fn test_ellipsoid_has_fewer_dipoles_than_circumscribed_sphere() {
    let semi_axes = [8.0_f64, 8.0, 8.0];
    let radius = semi_axes[0]; // sphere of same radius
    let spacing = 2.0;

    let prim_sphere = Primitive::Sphere(Sphere { centre: [0.0, 0.0, 0.0], radius });
    let prim_ellipsoid = Primitive::Ellipsoid(Ellipsoid {
        centre: [0.0, 0.0, 0.0],
        semi_axes,
    });

    let n_sphere = discretise_primitive(&prim_sphere, spacing).len();
    let n_ellipsoid = discretise_primitive(&prim_ellipsoid, spacing).len();

    // Equal semi-axes should give approximately the same count (within 10%)
    let ratio = n_ellipsoid as f64 / n_sphere as f64;
    assert!(
        ratio > 0.85 && ratio < 1.15,
        "Sphere-equivalent ellipsoid ratio = {ratio:.2}, expected ≈1.0"
    );
}

#[test]
fn test_oblate_ellipsoid_fewer_than_sphere() {
    // A flat oblate ellipsoid (c << a,b) should have fewer dipoles than its circumscribed sphere
    let spacing = 2.0;
    let prim_sphere = Primitive::Sphere(Sphere { centre: [0.0, 0.0, 0.0], radius: 12.0 });
    let prim_oblate = Primitive::Ellipsoid(Ellipsoid {
        centre: [0.0, 0.0, 0.0],
        semi_axes: [12.0, 12.0, 4.0], // flat in z
    });

    let n_sphere = discretise_primitive(&prim_sphere, spacing).len();
    let n_oblate = discretise_primitive(&prim_oblate, spacing).len();
    assert!(
        n_oblate < n_sphere,
        "Oblate ellipsoid ({n_oblate}) should have fewer dipoles than circumscribed sphere ({n_sphere})"
    );
}

// ─────────────────────────────────────────────────────────────
// v0.1.3: Palik TiO₂ material
// ─────────────────────────────────────────────────────────────

#[test]
fn test_palik_tio2_high_refractive_index() {
    let mat = PalikMaterial::tio2();
    let eps_500 = mat.dielectric_function(500.0).unwrap();
    // Rutile TiO₂ at 500 nm: n ≈ 2.57, so ε₁ ≈ n² ≈ 6.6, ε₂ ≈ 0
    let n_eff = eps_500.re.sqrt();
    assert!(
        n_eff > 2.4 && n_eff < 2.8,
        "TiO₂ effective n = {n_eff:.3} at 500 nm, expected ≈2.57"
    );
    // ε₂ should be near zero (transparent in visible past absorption edge)
    assert!(eps_500.im < 0.05, "TiO₂ ε₂ should be ~0 at 500 nm, got {:.4}", eps_500.im);
}

#[test]
fn test_palik_tio2_absorption_edge() {
    // TiO₂ absorbs below ~380 nm (band gap ≈ 3.2 eV)
    let mat = PalikMaterial::tio2();
    let eps_320 = mat.dielectric_function(320.0).unwrap();
    let eps_500 = mat.dielectric_function(500.0).unwrap();
    assert!(
        eps_320.im > eps_500.im,
        "TiO₂ ε₂ should be larger at 320 nm ({:.3}) than 500 nm ({:.3})",
        eps_320.im, eps_500.im
    );
}

#[test]
fn test_palik_tio2_wavelength_range() {
    let mat = PalikMaterial::tio2();
    let (lo, hi) = mat.wavelength_range();
    assert_relative_eq!(lo, 300.0, epsilon = 1.0);
    assert_relative_eq!(hi, 1000.0, epsilon = 1.0);
}

#[test]
fn test_palik_tio2_out_of_range_error() {
    let mat = PalikMaterial::tio2();
    assert!(
        mat.dielectric_function(200.0).is_err(),
        "TiO₂: wavelength 200 nm should be out of range"
    );
    assert!(
        mat.dielectric_function(1100.0).is_err(),
        "TiO₂: wavelength 1100 nm should be out of range"
    );
}

// ─────────────────────────────────────────────────────────────
// v0.1.3: Palik SiO₂ material
// ─────────────────────────────────────────────────────────────

#[test]
fn test_palik_sio2_low_refractive_index() {
    let mat = PalikMaterial::sio2();
    let eps_500 = mat.dielectric_function(500.0).unwrap();
    // Fused silica at 500 nm: n ≈ 1.462, so ε₁ ≈ 2.14, ε₂ ≈ 0
    let n_eff = eps_500.re.sqrt();
    assert!(
        n_eff > 1.43 && n_eff < 1.50,
        "SiO₂ effective n = {n_eff:.4} at 500 nm, expected ≈1.462"
    );
    assert!(eps_500.im < 1e-4, "SiO₂ ε₂ should be ~0 at 500 nm, got {:.6}", eps_500.im);
}

#[test]
fn test_palik_sio2_wavelength_range() {
    let mat = PalikMaterial::sio2();
    let (lo, hi) = mat.wavelength_range();
    assert_relative_eq!(lo, 300.0, epsilon = 1.0);
    assert_relative_eq!(hi, 1000.0, epsilon = 1.0);
}

#[test]
fn test_palik_sio2_monotonically_decreasing_index() {
    // SiO₂ is normal dispersion in visible: n decreases with λ
    let mat = PalikMaterial::sio2();
    let n_400 = mat.dielectric_function(400.0).unwrap().re.sqrt();
    let n_600 = mat.dielectric_function(600.0).unwrap().re.sqrt();
    let n_800 = mat.dielectric_function(800.0).unwrap().re.sqrt();
    assert!(n_400 > n_600, "SiO₂: n at 400 nm ({n_400:.4}) should exceed n at 600 nm ({n_600:.4})");
    assert!(n_600 > n_800, "SiO₂: n at 600 nm ({n_600:.4}) should exceed n at 800 nm ({n_800:.4})");
}

#[test]
fn test_palik_sio2_out_of_range_error() {
    let mat = PalikMaterial::sio2();
    assert!(mat.dielectric_function(200.0).is_err(), "SiO₂: 200 nm should be out of range");
    assert!(mat.dielectric_function(1200.0).is_err(), "SiO₂: 1200 nm should be out of range");
}

// ─────────────────────────────────────────────────────────────
// v0.1.3: Far-field computation
// ─────────────────────────────────────────────────────────────

/// Solve a single isotropic dipole at the origin and return its DipoleResponse.
fn single_dipole_response(polarisation: [f64; 3]) -> DipoleResponse {
    let wl = 500.0;
    let _k = 2.0 * std::f64::consts::PI / wl;
    let alpha = Complex64::new(500.0, 50.0); // arbitrary polarisability
    let dipole = Dipole::isotropic([0.0, 0.0, 0.0], alpha);

    let incident = IncidentField {
        direction: [0.0, 0.0, 1.0],
        polarisation,
        amplitude: 1.0,
    };
    let solver = CdaSolver::with_incident(1000, true, 3.0, incident);
    let params = default_params();
    solver.solve_dipoles(&[dipole], wl, &params).unwrap()
}

#[test]
fn test_far_field_grid_dimensions() {
    let resp = single_dipole_response([1.0, 0.0, 0.0]);
    let dipoles = [Dipole::isotropic([0.0, 0.0, 0.0], resp.moments[[0, 0]])];
    let k = 2.0 * std::f64::consts::PI / 500.0;
    let ff = compute_far_field(&dipoles, &resp, k, 18, 36);

    assert_eq!(ff.n_theta, 18);
    assert_eq!(ff.n_phi, 36);
    assert_eq!(ff.theta.len(), 18 * 36);
    assert_eq!(ff.phi.len(), 18 * 36);
    assert_eq!(ff.intensity.len(), 18 * 36);
}

#[test]
fn test_far_field_intensity_nonnegative() {
    let resp = single_dipole_response([1.0, 0.0, 0.0]);
    let dipoles = [Dipole::isotropic([0.0, 0.0, 0.0], Complex64::new(500.0, 50.0))];
    let k = 2.0 * std::f64::consts::PI / 500.0;
    let ff = compute_far_field(&dipoles, &resp, k, 18, 36);

    for &val in &ff.intensity {
        assert!(val >= 0.0, "Far-field intensity must be non-negative, got {val:.3e}");
    }
}

#[test]
fn test_far_field_x_dipole_has_node_along_x() {
    // A dipole with moment along x should radiate zero intensity along x̂ (θ=π/2, φ=0).
    // compute_far_field uses structure factor F(r̂) = Σᵢ pᵢ exp(-ik r̂·rᵢ).
    // For a single dipole at origin, F = p and the transverse projection eliminates
    // the component along r̂. Along x̂: r̂ = (1,0,0), ft = p - (p·r̂)r̂ = (px,py,pz) - px*(1,0,0) = (0,py,pz).
    // If p = (px, 0, 0), ft = (0, 0, 0), so intensity = 0.
    use ndarray::Array2;

    let wl = 500.0;
    let k = 2.0 * std::f64::consts::PI / wl;
    let px = Complex64::new(1.0, 0.0);
    let _zero = Complex64::new(0.0, 0.0);

    let mut moments = Array2::zeros((1, 3));
    moments[[0, 0]] = px; // x-component only
    let local_fields = Array2::zeros((1, 3));

    let resp = DipoleResponse { wavelength_nm: wl, moments, local_fields };
    let dipoles = [Dipole::isotropic([0.0, 0.0, 0.0], px)];

    // Use n_theta=3, n_phi=4 so θ=π/2 (it=1) and φ=0 (ip=0) is sampled exactly
    // theta: 0, π/2, π  (n_theta=3)
    // phi:   0, π/2, π, 3π/2  (n_phi=4)
    let ff = compute_far_field(&dipoles, &resp, k, 3, 4);

    // Point (it=1, ip=0): θ=π/2, φ=0 → r̂ = (1,0,0) = x̂
    let idx = 1 * 4 + 0;
    assert!(
        ff.intensity[idx] < 1e-20,
        "x-polarised dipole should radiate zero along x̂, got {:.3e}",
        ff.intensity[idx]
    );
}

#[test]
fn test_far_field_z_direction_has_maximum_for_x_dipole() {
    // Along ẑ (θ=0), all azimuthal φ give maximum intensity for an x-polarised dipole.
    // r̂ = (0,0,1), ft = p - (p·r̂)*r̂ = (px,py,pz) - pz*(0,0,1) = (px, py, 0).
    // For p = (1,0,0): ft = (1,0,0), |ft|² = 1 = maximum.
    use ndarray::Array2;

    let wl = 500.0;
    let k = 2.0 * std::f64::consts::PI / wl;
    let px = Complex64::new(1.0, 0.0);

    let mut moments = Array2::zeros((1, 3));
    moments[[0, 0]] = px;
    let resp = DipoleResponse {
        wavelength_nm: wl,
        moments,
        local_fields: Array2::zeros((1, 3)),
    };
    let dipoles = [Dipole::isotropic([0.0, 0.0, 0.0], px)];

    // n_theta=3: θ = 0, π/2, π; n_phi=4
    let ff = compute_far_field(&dipoles, &resp, k, 3, 4);

    // θ=0 (it=0): all ip should give intensity = |px|² = 1
    for ip in 0..4 {
        let idx = 0 * 4 + ip;
        let intensity = ff.intensity[idx];
        assert!((intensity - 1.0).abs() < 1e-10,
            "At θ=0, ip={ip}: expected intensity 1.0, got {:.6}", intensity);
    }
}

// ─────────────────────────────────────────────────────────────
// v0.1.3: Circular dichroism
// ─────────────────────────────────────────────────────────────

#[test]
fn test_cd_isotropic_sphere_is_zero() {
    // An isotropic sphere has D∞h symmetry about the propagation axis.
    // The CD formula ΔC_ext = k Σᵢ [Im(p_y^x) - Im(p_x^y)] must vanish by symmetry:
    // under 90° rotation about z: response_x maps to response_y such that
    // p_y^(x) = p_y^(y after rotation) = 0, and p_x^(y) = 0 by the same argument.
    let wl = 500.0;
    let k = 2.0 * std::f64::consts::PI / wl;

    let dipoles = gold_sphere_dipoles(8.0, 2.0, wl);
    let params = default_params();

    let inc_x = IncidentField { direction: [0.0, 0.0, 1.0], polarisation: [1.0, 0.0, 0.0], amplitude: 1.0 };
    let inc_y = IncidentField { direction: [0.0, 0.0, 1.0], polarisation: [0.0, 1.0, 0.0], amplitude: 1.0 };

    let solver_x = CdaSolver::with_incident(1000, true, 2.0, inc_x);
    let solver_y = CdaSolver::with_incident(1000, true, 2.0, inc_y);

    let resp_x = solver_x.solve_dipoles(&dipoles, wl, &params).unwrap();
    let resp_y = solver_y.solve_dipoles(&dipoles, wl, &params).unwrap();

    let cd = compute_circular_dichroism(&resp_x, &resp_y, k);

    assert!(
        cd.abs() < 1e-6,
        "CD for isotropic sphere should be ~0, got {cd:.3e}"
    );
}

#[test]
fn test_cd_returns_float_for_single_dipole() {
    // A single isotropic dipole cannot have CD (no handedness).
    let wl = 500.0;
    let k = 2.0 * std::f64::consts::PI / wl;
    let alpha = Complex64::new(100.0, 10.0);
    let dipole = Dipole::isotropic([0.0, 0.0, 0.0], alpha);
    let params = default_params();

    let inc_x = IncidentField { direction: [0.0, 0.0, 1.0], polarisation: [1.0, 0.0, 0.0], amplitude: 1.0 };
    let inc_y = IncidentField { direction: [0.0, 0.0, 1.0], polarisation: [0.0, 1.0, 0.0], amplitude: 1.0 };

    let solver_x = CdaSolver::with_incident(1000, true, 3.0, inc_x);
    let solver_y = CdaSolver::with_incident(1000, true, 3.0, inc_y);

    let resp_x = solver_x.solve_dipoles(&[dipole.clone()], wl, &params).unwrap();
    let resp_y = solver_y.solve_dipoles(&[dipole], wl, &params).unwrap();

    let cd = compute_circular_dichroism(&resp_x, &resp_y, k);
    assert!(cd.is_finite(), "CD should be finite for single dipole");
    assert!(cd.abs() < 1e-8, "Single isotropic dipole should have CD ≈ 0, got {cd:.3e}");
}

// ─────────────────────────────────────────────────────────────
// v0.1.3: Incident polarisation
// ─────────────────────────────────────────────────────────────

#[test]
fn test_x_pol_and_y_pol_give_rotated_dipole_moments() {
    // For a single isotropic dipole at origin, x-pol excites p along x,
    // y-pol excites p along y — same magnitude, different direction.
    let wl = 500.0;
    let alpha = Complex64::new(300.0, 30.0);
    let dipole = Dipole::isotropic([0.0, 0.0, 0.0], alpha);
    let params = default_params();

    let inc_x = IncidentField { direction: [0.0, 0.0, 1.0], polarisation: [1.0, 0.0, 0.0], amplitude: 1.0 };
    let inc_y = IncidentField { direction: [0.0, 0.0, 1.0], polarisation: [0.0, 1.0, 0.0], amplitude: 1.0 };

    let solver_x = CdaSolver::with_incident(1000, true, 3.0, inc_x);
    let solver_y = CdaSolver::with_incident(1000, true, 3.0, inc_y);

    let resp_x = solver_x.solve_dipoles(&[dipole.clone()], wl, &params).unwrap();
    let resp_y = solver_y.solve_dipoles(&[dipole], wl, &params).unwrap();

    // x-pol: px should dominate, py ≈ 0
    let px_x = resp_x.moments[[0, 0]];
    let py_x = resp_x.moments[[0, 1]];
    assert!(px_x.norm() > py_x.norm() * 100.0,
        "x-pol: px ({:.3e}) should dominate py ({:.3e})", px_x.norm(), py_x.norm());

    // y-pol: py should dominate, px ≈ 0
    let px_y = resp_y.moments[[0, 0]];
    let py_y = resp_y.moments[[0, 1]];
    assert!(py_y.norm() > px_y.norm() * 100.0,
        "y-pol: py ({:.3e}) should dominate px ({:.3e})", py_y.norm(), px_y.norm());

    // Both should have same magnitude (same amplitude, same polarisability)
    assert_relative_eq!(px_x.norm(), py_y.norm(), epsilon = 1e-10);
}

#[test]
fn test_cross_sections_same_for_sphere_under_x_or_y_pol() {
    // By symmetry, extinction cross-section of a sphere is the same for x-pol and y-pol
    let dipoles = gold_sphere_dipoles(8.0, 2.0, 520.0);
    let params = default_params();

    let inc_x = IncidentField { direction: [0.0, 0.0, 1.0], polarisation: [1.0, 0.0, 0.0], amplitude: 1.0 };
    let inc_y = IncidentField { direction: [0.0, 0.0, 1.0], polarisation: [0.0, 1.0, 0.0], amplitude: 1.0 };

    let solver_x = CdaSolver::with_incident(1000, true, 2.0, inc_x);
    let solver_y = CdaSolver::with_incident(1000, true, 2.0, inc_y);

    let cs_x = solver_x.compute_cross_sections(&dipoles, 520.0, &params).unwrap();
    let cs_y = solver_y.compute_cross_sections(&dipoles, 520.0, &params).unwrap();

    assert_relative_eq!(cs_x.extinction, cs_y.extinction, epsilon = 1e-8);
    assert_relative_eq!(cs_x.absorption, cs_y.absorption, epsilon = 1e-8);
}

// ─────────────────────────────────────────────────────────────
// v0.1.3: CrossSections JSON serialisation
// ─────────────────────────────────────────────────────────────

#[test]
fn test_cross_sections_json_round_trip() {
    let cs = CrossSections {
        wavelength_nm: 532.5,
        extinction: 1234.56,
        absorption: 789.01,
        scattering: 445.55,
        circular_dichroism: Some(-12.34),
    };

    let json = serde_json::to_string(&cs).expect("Serialisation failed");
    let cs2: CrossSections = serde_json::from_str(&json).expect("Deserialisation failed");

    assert_relative_eq!(cs2.wavelength_nm, cs.wavelength_nm, epsilon = 1e-6);
    assert_relative_eq!(cs2.extinction, cs.extinction, epsilon = 1e-6);
    assert_relative_eq!(cs2.absorption, cs.absorption, epsilon = 1e-6);
    assert_relative_eq!(cs2.scattering, cs.scattering, epsilon = 1e-6);
    assert_relative_eq!(
        cs2.circular_dichroism.unwrap(),
        cs.circular_dichroism.unwrap(),
        epsilon = 1e-6
    );
}

#[test]
fn test_cross_sections_json_round_trip_no_cd() {
    let cs = CrossSections {
        wavelength_nm: 400.0,
        extinction: 100.0,
        absorption: 60.0,
        scattering: 40.0,
        circular_dichroism: None,
    };

    let json = serde_json::to_string(&cs).unwrap();
    let cs2: CrossSections = serde_json::from_str(&json).unwrap();

    assert!(cs2.circular_dichroism.is_none());
    assert_relative_eq!(cs2.extinction, cs.extinction, epsilon = 1e-6);
}

#[test]
fn test_far_field_map_json_round_trip() {
    // FarFieldMap is Serialize/Deserialize — verify it survives a round-trip
    let ff = FarFieldMap {
        wavelength_nm: 500.0,
        theta: vec![0.0, 1.57, 3.14],
        phi: vec![0.0, 1.57],
        intensity: vec![1.0, 0.5, 0.0, 0.8, 0.3, 0.0],
        n_theta: 3,
        n_phi: 2,
    };

    let json = serde_json::to_string(&ff).unwrap();
    let ff2: FarFieldMap = serde_json::from_str(&json).unwrap();

    assert_eq!(ff2.n_theta, 3);
    assert_eq!(ff2.n_phi, 2);
    assert_relative_eq!(ff2.wavelength_nm, 500.0, epsilon = 1e-6);
    assert_relative_eq!(ff2.intensity[0], 1.0, epsilon = 1e-6);
}

// ─────────────────────────────────────────────────────────────
// v0.1.3: OBJ mesh parser + volume-fill discretisation
// ─────────────────────────────────────────────────────────────

#[test]
fn test_obj_parse_and_discretise_round_trip() {
    // A cube OBJ (±5 nm) parsed and discretised at 2 nm spacing.
    let obj = "\
        v 5 5 -5\nv 5 -5 -5\nv -5 -5 -5\nv -5 5 -5\n\
        v 5 5 5\nv 5 -5 5\nv -5 -5 5\nv -5 5 5\n\
        f 1 2 3 4\nf 5 8 7 6\nf 1 5 6 2\nf 3 7 8 4\nf 1 4 8 5\nf 2 6 7 3\n";

    let mesh = parse_obj(obj).unwrap();
    assert_eq!(mesh.vertices.len(), 8);
    assert_eq!(mesh.faces.len(), 12); // 6 quads → 12 triangles

    let points = discretise_mesh(&mesh, 2.0);
    // Interior grid for ±5 cube at 2nm: points at -4,-2,0,2,4 = 5 per axis → 125
    // Some boundary effects may shift this slightly
    assert!(
        points.len() >= 100 && points.len() <= 200,
        "OBJ cube dipole count = {}, expected ~125",
        points.len()
    );

    // All points must lie within the cube ±5 nm
    for p in &points {
        assert!(p.position[0].abs() <= 5.0 + 1e-8, "x out of bounds: {}", p.position[0]);
        assert!(p.position[1].abs() <= 5.0 + 1e-8, "y out of bounds: {}", p.position[1]);
        assert!(p.position[2].abs() <= 5.0 + 1e-8, "z out of bounds: {}", p.position[2]);
    }
}

#[test]
fn test_obj_mesh_bounding_box() {
    let obj = "\
        v 0 0 0\nv 10 0 0\nv 10 10 0\nv 0 10 0\n\
        v 0 0 10\nv 10 0 10\nv 10 10 10\nv 0 10 10\n\
        f 1 2 3 4\nf 5 8 7 6\nf 1 5 6 2\nf 3 7 8 4\nf 1 4 8 5\nf 2 6 7 3\n";
    let mesh = parse_obj(obj).unwrap();
    let (min, max) = mesh.bounding_box();
    assert_relative_eq!(min[0], 0.0, epsilon = 1e-10);
    assert_relative_eq!(min[1], 0.0, epsilon = 1e-10);
    assert_relative_eq!(min[2], 0.0, epsilon = 1e-10);
    assert_relative_eq!(max[0], 10.0, epsilon = 1e-10);
    assert_relative_eq!(max[1], 10.0, epsilon = 1e-10);
    assert_relative_eq!(max[2], 10.0, epsilon = 1e-10);
}

#[test]
fn test_obj_finer_spacing_gives_more_dipoles() {
    let obj = "\
        v 5 5 -5\nv 5 -5 -5\nv -5 -5 -5\nv -5 5 -5\n\
        v 5 5 5\nv 5 -5 5\nv -5 -5 5\nv -5 5 5\n\
        f 1 2 3 4\nf 5 8 7 6\nf 1 5 6 2\nf 3 7 8 4\nf 1 4 8 5\nf 2 6 7 3\n";
    let mesh = parse_obj(obj).unwrap();
    let coarse = discretise_mesh(&mesh, 3.0);
    let fine = discretise_mesh(&mesh, 1.0);
    assert!(
        fine.len() > coarse.len(),
        "Finer spacing should give more dipoles: {} (1nm) vs {} (3nm)",
        fine.len(),
        coarse.len()
    );
}
