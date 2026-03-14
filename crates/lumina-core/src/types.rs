//! Core types shared across the Lumina framework.
//!
//! This module defines the fundamental data structures used throughout the
//! simulation pipeline: dipoles, simulation parameters, and result containers.

use ndarray::Array2;
use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// A single point dipole in the simulation domain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Dipole {
    /// Position in 3D space (nanometres).
    pub position: [f64; 3],
    /// Complex polarisability tensor (3x3) at the current frequency.
    /// Stored as a flat [9] array in row-major order.
    pub polarisability: [Complex64; 9],
}

impl Dipole {
    /// Create a dipole with isotropic (scalar) polarisability.
    pub fn isotropic(position: [f64; 3], alpha: Complex64) -> Self {
        let zero = Complex64::new(0.0, 0.0);
        Self {
            position,
            polarisability: [
                alpha, zero, zero,
                zero, alpha, zero,
                zero, zero, alpha,
            ],
        }
    }
}

/// Incident plane wave specification.
#[derive(Debug, Clone)]
pub struct IncidentField {
    /// Propagation direction (unit vector).
    pub direction: [f64; 3],
    /// Polarisation vector (unit vector, perpendicular to direction).
    pub polarisation: [f64; 3],
    /// Amplitude (V/m). Typically set to 1.0 for cross-section calculations.
    pub amplitude: f64,
}

impl Default for IncidentField {
    fn default() -> Self {
        Self {
            direction: [0.0, 0.0, 1.0],    // propagating along z
            polarisation: [1.0, 0.0, 0.0],  // x-polarised
            amplitude: 1.0,
        }
    }
}

impl IncidentField {
    /// Evaluate the incident field at a given position.
    ///
    /// $\mathbf{E}_{\text{inc}}(\mathbf{r}) = E_0 \hat{\mathbf{e}} \exp(i \mathbf{k} \cdot \mathbf{r})$
    pub fn at_position(&self, position: &[f64; 3], k: f64) -> [Complex64; 3] {
        let kdotr = k * (
            self.direction[0] * position[0]
            + self.direction[1] * position[1]
            + self.direction[2] * position[2]
        );
        let phase = Complex64::new(0.0, kdotr).exp();
        let e0 = Complex64::from(self.amplitude);
        [
            e0 * self.polarisation[0] * phase,
            e0 * self.polarisation[1] * phase,
            e0 * self.polarisation[2] * phase,
        ]
    }
}

/// Parameters defining a simulation run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationParams {
    /// Wavelength range [start, end] in nanometres.
    pub wavelength_range: [f64; 2],
    /// Number of wavelength points to evaluate.
    pub num_wavelengths: usize,
    /// Refractive index of the surrounding medium.
    pub environment_n: f64,
    /// Solver tolerance for iterative methods (e.g. GMRES).
    pub solver_tolerance: f64,
    /// Maximum iterations for iterative solvers.
    pub max_iterations: usize,
    /// Bloch wavevector **k**∥ in nm⁻¹ (in-plane, for periodic structures).
    /// Use `[0.0, 0.0, 0.0]` for the Γ point (normal incidence).
    #[serde(default)]
    pub k_bloch: [f64; 3],
    /// Runtime substrate correction at the current wavelength.
    /// `None` = free-space (no substrate correction).
    ///
    /// Not parsed from TOML; set per-wavelength by the runner / GUI after
    /// resolving the substrate material's dielectric function. Contains either
    /// a scalar Fresnel factor or a full Sommerfeld integrator.
    #[serde(skip, default)]
    pub substrate_runtime: Option<crate::solver::substrate::SubstrateRuntime>,
}

impl Default for SimulationParams {
    fn default() -> Self {
        Self {
            wavelength_range: [400.0, 900.0],
            num_wavelengths: 100,
            environment_n: 1.0,
            solver_tolerance: 1e-6,
            max_iterations: 1000,
            k_bloch: [0.0, 0.0, 0.0],
            substrate_runtime: None,
        }
    }
}

/// Optical cross-sections at a single wavelength.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CrossSections {
    /// Wavelength (nm).
    pub wavelength_nm: f64,
    /// Extinction cross-section (nm²).
    pub extinction: f64,
    /// Absorption cross-section (nm²).
    pub absorption: f64,
    /// Scattering cross-section (nm²).
    pub scattering: f64,
    /// Circular dichroism ΔC_ext = C_ext(LCP) − C_ext(RCP) (nm²).
    /// `None` unless CD mode is enabled.
    pub circular_dichroism: Option<f64>,
}

/// Optical cross-sections at a single (wavelength, **k**∥) point.
///
/// Extends [`CrossSections`] with the in-plane Bloch wavevector for
/// periodic-array dispersion calculations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlochCrossSections {
    /// Wavelength (nm).
    pub wavelength_nm: f64,
    /// In-plane Bloch wavevector **k**∥ (nm⁻¹).
    pub k_parallel: [f64; 3],
    /// Extinction cross-section (nm²).
    pub extinction: f64,
    /// Absorption cross-section (nm²).
    pub absorption: f64,
    /// Scattering cross-section (nm²).
    pub scattering: f64,
    /// Circular dichroism ΔC_ext (nm²). `None` unless CD mode is active.
    pub circular_dichroism: Option<f64>,
}

/// Full dispersion map: cross-sections sampled on a 2D (ω, **k**∥) grid.
///
/// # Layout
/// `data[i * k_points.len() + j]` is the response at
/// `wavelengths_nm[i]` and `k_points[j]`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DispersionMap {
    /// Wavelength grid (nm).
    pub wavelengths_nm: Vec<f64>,
    /// In-plane **k**∥ points (nm⁻¹).
    pub k_points: Vec<[f64; 3]>,
    /// Flattened array of cross-sections; length = wavelengths × k_points.
    pub data: Vec<BlochCrossSections>,
}

/// Far-field radiation pattern sampled on the unit sphere.
///
/// The differential scattering intensity at each direction r̂ is proportional to
/// $|(\mathbf{I} - \hat{\mathbf{r}}\hat{\mathbf{r}}) \cdot \mathbf{F}(\hat{\mathbf{r}})|^2$
/// where $\mathbf{F}(\hat{\mathbf{r}}) = \sum_i \mathbf{p}_i \exp(-ik\hat{\mathbf{r}} \cdot \mathbf{r}_i)$.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FarFieldMap {
    /// Wavelength at which the pattern was computed (nm).
    pub wavelength_nm: f64,
    /// Polar angles θ ∈ [0, π] (radians), one per sample.
    pub theta: Vec<f64>,
    /// Azimuthal angles φ ∈ [0, 2π) (radians), one per sample.
    pub phi: Vec<f64>,
    /// Differential scattering intensity (nm⁴, unnormalised) at each (θ, φ).
    pub intensity: Vec<f64>,
    /// Number of θ samples.
    pub n_theta: usize,
    /// Number of φ samples.
    pub n_phi: usize,
}

/// The solved dipole moments for a given excitation.
#[derive(Debug, Clone)]
pub struct DipoleResponse {
    /// Wavelength (nm) at which this was computed.
    pub wavelength_nm: f64,
    /// Complex dipole moments, shape (N, 3).
    pub moments: Array2<Complex64>,
    /// Local electric field at each dipole, shape (N, 3).
    pub local_fields: Array2<Complex64>,
}

/// Near-field map on a 2D plane.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NearFieldMap {
    /// Grid positions, shape (Nx * Ny, 3).
    pub positions: Vec<[f64; 3]>,
    /// |E|^2 values at each grid point.
    pub field_intensity: Vec<f64>,
    /// Number of points along x.
    pub nx: usize,
    /// Number of points along y.
    pub ny: usize,
    /// Spatial extent: [x_min, x_max, y_min, y_max] in nm.
    pub extent: [f64; 4],
}

/// Complete results from a simulation run.
#[derive(Debug, Clone)]
pub struct SimulationResult {
    /// Cross-section spectra across all wavelengths.
    pub spectra: Vec<CrossSections>,
    /// Dipole responses at each wavelength (optional, can be large).
    pub dipole_responses: Option<Vec<DipoleResponse>>,
    /// Near-field maps (optional).
    pub near_field_maps: Option<Vec<NearFieldMap>>,
}

/// Compute the Clausius-Mossotti polarisability for a small volume element.
///
/// $\alpha_{\text{CM}} = 3 V \epsilon_0 \frac{\epsilon - \epsilon_m}{\epsilon + 2\epsilon_m}$
///
/// Note: we work in Gaussian-like units where $\epsilon_0$ is absorbed, so the
/// returned polarisability has units of volume (nm^3).
///
/// # Arguments
/// * `volume_nm3` - Volume of the dipole element in nm^3 (typically $d^3$).
/// * `epsilon` - Complex dielectric function of the material.
/// * `epsilon_m` - Dielectric constant of the surrounding medium ($n_m^2$).
pub fn clausius_mossotti(volume_nm3: f64, epsilon: Complex64, epsilon_m: f64) -> Complex64 {
    let eps_m = Complex64::from(epsilon_m);
    3.0 * volume_nm3 * (epsilon - eps_m) / (epsilon + 2.0 * eps_m)
}

/// Apply the Draine radiative correction to the Clausius-Mossotti polarisability.
///
/// This ensures optical theorem consistency (energy conservation):
/// $\alpha_{\text{rad}} = \frac{\alpha_{\text{CM}}}{1 - \frac{i k^3}{6\pi} \alpha_{\text{CM}}}$
///
/// # Arguments
/// * `alpha_cm` - Clausius-Mossotti polarisability.
/// * `k` - Wavenumber in the medium (nm^{-1}).
pub fn radiative_correction(alpha_cm: Complex64, k: f64) -> Complex64 {
    let correction = Complex64::new(0.0, k.powi(3) / (6.0 * std::f64::consts::PI));
    alpha_cm / (Complex64::from(1.0) - correction * alpha_cm)
}

/// Compute the full 3×3 polarisability tensor for an ellipsoidal dipole cell.
///
/// Uses the generalised Clausius-Mossotti formula with per-axis depolarisation
/// factors, then applies the Draine radiative reaction correction per component,
/// and finally rotates the diagonal tensor to the lab (world) frame:
///
/// $$\boldsymbol{\alpha} = R\,\text{diag}(\alpha_x, \alpha_y, \alpha_z)\,R^T$$
///
/// The per-axis CM polarisability is:
/// $$\alpha_i^{\text{CM}} = V \frac{\varepsilon - \varepsilon_m}
///   {\varepsilon_m + L_i(\varepsilon - \varepsilon_m)}$$
///
/// Reduces to the standard isotropic formula when $L_i = 1/3$ and $R = I$.
///
/// # Arguments
/// * `volume_nm3`  — Dipole cell volume in nm³ (typically $d^3$ for spacing $d$).
/// * `depol`       — Depolarisation factors $[L_x, L_y, L_z]$ in the object-local frame.
/// * `epsilon`     — Complex dielectric function of the material at this wavelength.
/// * `epsilon_m`   — Surrounding medium dielectric constant $\varepsilon_m = n_m^2$.
/// * `k`           — Wavenumber in the medium in nm⁻¹.
/// * `rot`         — 3×3 rotation matrix (row-major `[f64; 9]`) from the
///                   object-local principal-axis frame to the lab frame.
///
/// # Returns
/// Full 3×3 complex polarisability tensor as `[Complex64; 9]` in row-major order.
pub fn ellipsoid_polarisability_tensor(
    volume_nm3: f64,
    depol: [f64; 3],
    epsilon: Complex64,
    epsilon_m: f64,
    k: f64,
    rot: [f64; 9],
) -> [Complex64; 9] {
    let eps_m = Complex64::from(epsilon_m);
    let delta_eps = epsilon - eps_m;

    // Compute the three principal-axis polarisabilities.
    let alphas: [Complex64; 3] = std::array::from_fn(|i| {
        let denom = eps_m + Complex64::from(depol[i]) * delta_eps;
        let alpha_cm = volume_nm3 * delta_eps / denom;
        radiative_correction(alpha_cm, k)
    });

    // Rotate to lab frame: α = R · diag(αx, αy, αz) · Rᵀ
    // Element (i,j) = Σ_k rot[i*3+k] * alphas[k] * rot[j*3+k]
    let mut tensor = [Complex64::from(0.0); 9];
    for i in 0..3 {
        for j in 0..3 {
            let mut s = Complex64::from(0.0);
            for k in 0..3 {
                s += alphas[k] * rot[i * 3 + k] * rot[j * 3 + k];
            }
            tensor[i * 3 + j] = s;
        }
    }
    tensor
}

/// Second-order nonlinear susceptibility tensor χ^(2).
///
/// A rank-3 tensor with components $\chi^{(2)}_{abc}$, where index $a$ is the
/// output polarisation direction and $b$, $c$ are the input polarisation
/// directions of the two fundamental photons.
///
/// The nonlinear dipole moment at $2\omega$ is:
/// $$p_a(2\omega) = \sum_{b,c} \chi^{(2)}_{abc}\,
///   E_{\text{loc},b}(\omega)\, E_{\text{loc},c}(\omega)$$
///
/// In our unit convention (ε₀ absorbed, lengths in nm), χ^(2) has units of nm³.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Chi2Tensor {
    /// Tensor components in row-major order.
    /// Component $(a,b,c)$ is at index `a*9 + b*3 + c`.
    pub components: [Complex64; 27],
}

impl Chi2Tensor {
    /// Zero tensor — no second-order response (centrosymmetric material).
    pub fn zero() -> Self {
        Self { components: [Complex64::new(0.0, 0.0); 27] }
    }

    /// Construct from explicit 27 components (row-major).
    pub fn from_components(components: [Complex64; 27]) -> Self {
        Self { components }
    }

    /// Contract with a local field vector: $p_a = \sum_{b,c} \chi_{abc} E_b E_c$.
    ///
    /// # Arguments
    /// * `e` — Local electric field $[E_x, E_y, E_z]$ at frequency $\omega$.
    ///
    /// # Returns
    /// Nonlinear dipole moment $[p_x, p_y, p_z]$ at frequency $2\omega$.
    pub fn contract(&self, e: &[Complex64; 3]) -> [Complex64; 3] {
        let mut p = [Complex64::new(0.0, 0.0); 3];
        for a in 0..3 {
            for b in 0..3 {
                for c in 0..3 {
                    p[a] += self.components[a * 9 + b * 3 + c] * e[b] * e[c];
                }
            }
        }
        p
    }

    /// Isotropic surface tensor with surface normal along **z** (symmetry group C_∞v).
    ///
    /// The independent non-zero elements are:
    /// - $\chi_{zzz}$ (normal–normal–normal)
    /// - $\chi_{zxx} = \chi_{zyy}$ (normal–tangential–tangential)
    /// - $\chi_{xxz} = \chi_{xzx} = \chi_{yyz} = \chi_{yzy}$ (Kleinman symmetry)
    ///
    /// # Arguments
    /// * `chi_zzz` — Normal–normal–normal component.
    /// * `chi_zxx` — Normal–tangential–tangential component.
    ///   Kleinman symmetry gives $\chi_{xxz} = \chi_{zxx}$.
    pub fn isotropic_surface(chi_zzz: Complex64, chi_zxx: Complex64) -> Self {
        let mut t = Self::zero();
        // χ_zzz: (2,2,2)
        t.components[2 * 9 + 2 * 3 + 2] = chi_zzz;
        // χ_zxx = χ_zyy: (2,0,0) and (2,1,1)
        t.components[2 * 9 + 0 * 3 + 0] = chi_zxx;
        t.components[2 * 9 + 1 * 3 + 1] = chi_zxx;
        // χ_xxz = χ_xzx: (0,0,2) and (0,2,0)  [Kleinman: χ_xxz = χ_zxx]
        t.components[0 * 9 + 0 * 3 + 2] = chi_zxx;
        t.components[0 * 9 + 2 * 3 + 0] = chi_zxx;
        // χ_yyz = χ_yzy: (1,1,2) and (1,2,1)
        t.components[1 * 9 + 1 * 3 + 2] = chi_zxx;
        t.components[1 * 9 + 2 * 3 + 1] = chi_zxx;
        t
    }
}

/// Result from a second-harmonic generation (SHG) calculation.
///
/// Produced by [`crate::nonlinear::compute_shg_response`].
#[derive(Debug, Clone)]
pub struct ShgResult {
    /// Fundamental frequency wavelength $\omega$ (nm).
    pub fundamental_nm: f64,
    /// Second-harmonic wavelength $\lambda/2$ (nm).
    pub harmonic_nm: f64,
    /// Nonlinear source moments at $2\omega$ before inter-dipole coupling,
    /// shape $(N, 3)$.
    ///
    /// $p^{\text{NL}}_{i,a} = \sum_{b,c} \chi^{(2)}_{abc}\,
    ///   E_{\text{loc},b}(\omega)\, E_{\text{loc},c}(\omega)$
    pub source_moments: Array2<Complex64>,
    /// Self-consistently driven dipole moments at $2\omega$, shape $(N, 3)$.
    ///
    /// Solution of $\mathbf{A}(2\omega)\,\mathbf{p}(2\omega) = \mathbf{p}^{\text{NL}}$.
    pub driven_moments: Array2<Complex64>,
    /// Total SHG intensity $\propto \sum_i |\mathbf{p}_i(2\omega)|^2$ (nm⁶, unnormalised).
    ///
    /// Scales as $|E_{\text{inc}}|^4$. Use this for spectral sweeps and
    /// relative comparisons between structures or wavelengths.
    pub shg_intensity: f64,
    /// Far-field radiation pattern at $2\omega$. `None` unless `calc_far_field`
    /// was `true` when calling [`crate::nonlinear::compute_shg_response`].
    pub far_field: Option<FarFieldMap>,
}

/// Third-order nonlinear susceptibility tensor χ^(3).
///
/// A rank-4 tensor with components $\chi^{(3)}_{abcd}$, where index $a$ is the
/// output polarisation direction and $b$, $c$, $d$ are the input polarisation
/// directions of the three fundamental photons.
///
/// The nonlinear dipole moment at $3\omega$ is:
/// $$p_a(3\omega) = \sum_{b,c,d} \chi^{(3)}_{abcd}\,
///   E_{\text{loc},b}(\omega)\, E_{\text{loc},c}(\omega)\, E_{\text{loc},d}(\omega)$$
///
/// In our unit convention (ε₀ absorbed, lengths in nm), χ^(3) has units of nm⁶.
#[derive(Debug, Clone)]
pub struct Chi3Tensor {
    /// Tensor components in row-major order.
    /// Component $(a,b,c,d)$ is at index `a*27 + b*9 + c*3 + d`.
    pub components: [Complex64; 81],
}

impl Chi3Tensor {
    /// Zero tensor — no third-order response (e.g. vanishing χ^(3)).
    pub fn zero() -> Self {
        Self { components: [Complex64::new(0.0, 0.0); 81] }
    }

    /// Construct from explicit 81 components (row-major).
    pub fn from_components(components: [Complex64; 81]) -> Self {
        Self { components }
    }

    /// Contract with a local field vector:
    /// $p_a = \sum_{b,c,d} \chi_{abcd} E_b E_c E_d$.
    ///
    /// # Arguments
    /// * `e` — Local electric field $[E_x, E_y, E_z]$ at frequency $\omega$.
    ///
    /// # Returns
    /// Nonlinear dipole moment $[p_x, p_y, p_z]$ at frequency $3\omega$.
    pub fn contract(&self, e: &[Complex64; 3]) -> [Complex64; 3] {
        let mut out = [Complex64::new(0.0, 0.0); 3];
        for a in 0..3 {
            for b in 0..3 {
                for c in 0..3 {
                    for d in 0..3 {
                        out[a] += self.components[a * 27 + b * 9 + c * 3 + d]
                            * e[b]
                            * e[c]
                            * e[d];
                    }
                }
            }
        }
        out
    }

    /// Isotropic bulk (centrosymmetric, Kleinman-symmetric) tensor.
    ///
    /// For an isotropic medium the independent components are $\chi_{xxxx}$ and
    /// $\chi_{xxyy}$.  Under Kleinman symmetry
    /// $\chi_{xyxy} = \chi_{xyyx} = \chi_{xxyy}$.
    ///
    /// The full tensor:
    /// $$\chi_{abcd} = \chi_{xxyy}(\delta_{ab}\delta_{cd}
    ///   + \delta_{ac}\delta_{bd} + \delta_{ad}\delta_{bc})
    ///   + (\chi_{xxxx} - 3\chi_{xxyy})\,\delta_{ab}\delta_{ac}\delta_{ad}$$
    ///
    /// # Arguments
    /// * `chi_xxxx` — Diagonal component.
    /// * `chi_xxyy` — Off-diagonal component (Kleinman: $\chi_{xyxy} = \chi_{xxyy}$).
    pub fn isotropic_bulk(chi_xxxx: Complex64, chi_xxyy: Complex64) -> Self {
        let mut t = Self::zero();
        let extra_diag = chi_xxxx - Complex64::from(3.0) * chi_xxyy;
        for a in 0..3usize {
            for b in 0..3usize {
                for c in 0..3usize {
                    for d in 0..3usize {
                        let mut val = Complex64::new(0.0, 0.0);
                        // δ_ab δ_cd term
                        if a == b && c == d {
                            val += chi_xxyy;
                        }
                        // δ_ac δ_bd term
                        if a == c && b == d {
                            val += chi_xxyy;
                        }
                        // δ_ad δ_bc term
                        if a == d && b == c {
                            val += chi_xxyy;
                        }
                        // Diagonal correction: only when a=b=c=d
                        if a == b && b == c && c == d {
                            val += extra_diag;
                        }
                        t.components[a * 27 + b * 9 + c * 3 + d] = val;
                    }
                }
            }
        }
        t
    }
}

/// Result from a third-harmonic generation (THG) calculation.
///
/// Produced by [`crate::nonlinear::compute_thg_response`].
#[derive(Debug, Clone)]
pub struct ThgResult {
    /// Fundamental frequency wavelength $\omega$ (nm).
    pub fundamental_nm: f64,
    /// Third-harmonic wavelength $\lambda/3$ (nm).
    pub harmonic_nm: f64,
    /// Nonlinear source moments at $3\omega$ before inter-dipole coupling,
    /// shape $(N, 3)$.
    ///
    /// $p^{\text{NL}}_{i,a} = \sum_{b,c,d} \chi^{(3)}_{abcd}\,
    ///   E_{\text{loc},b}(\omega)\, E_{\text{loc},c}(\omega)\, E_{\text{loc},d}(\omega)$
    pub source_moments: Array2<Complex64>,
    /// Self-consistently driven dipole moments at $3\omega$, shape $(N, 3)$.
    ///
    /// Solution of $\mathbf{A}(3\omega)\,\mathbf{p}(3\omega) = \mathbf{p}^{\text{NL}}$.
    pub driven_moments: Array2<Complex64>,
    /// Total THG intensity $\propto \sum_i |\mathbf{p}_i(3\omega)|^2$ (nm⁹, unnormalised).
    ///
    /// Scales as $|E_{\text{inc}}|^6$. Use this for spectral sweeps and
    /// relative comparisons between structures or wavelengths.
    pub thg_intensity: f64,
    /// Far-field radiation pattern at $3\omega$. `None` unless `calc_far_field`
    /// was `true` when calling [`crate::nonlinear::compute_thg_response`].
    pub far_field: Option<FarFieldMap>,
}

/// Apply the Lattice Dispersion Relation (LDR) correction to the polarisability.
///
/// The LDR correction (Draine & Goodman, *Astrophys. J.* **405**, 685, 1993)
/// modifies the Clausius-Mossotti polarisability so that a cubic lattice of
/// point dipoles correctly reproduces the bulk dielectric dispersion relation.
/// This is critical for metallic particles where the standard CM prescription
/// gives poor convergence.
///
/// The corrected inverse polarisability is:
/// $$\alpha_{\text{LDR}}^{-1} = \alpha_{\text{CM}}^{-1}
///   - \frac{ik^3}{6\pi}
///   - \frac{k^2}{d}\bigl(b_1 + b_2 \epsilon + b_3 \epsilon S\bigr)$$
///
/// where $S$ depends on the propagation and polarisation directions.
///
/// # Arguments
/// * `alpha_cm` - Clausius-Mossotti polarisability (nm^3).
/// * `k` - Wavenumber in the medium (nm^{-1}).
/// * `d` - Lattice spacing (nm).
/// * `epsilon` - Complex dielectric function of the particle material.
/// * `propagation` - Propagation direction (unit vector).
/// * `polarisation` - Polarisation direction (unit vector).
pub fn ldr_correction(
    alpha_cm: Complex64,
    k: f64,
    d: f64,
    epsilon: Complex64,
    propagation: &[f64; 3],
    polarisation: &[f64; 3],
) -> Complex64 {
    // Draine & Goodman (1993) LDR coefficients
    let b1 = Complex64::new(-1.891_531_6, 0.0);
    let b2 = Complex64::new(0.164_846_9, 0.0);
    let b3 = Complex64::new(-1.770_000_4, 0.0);

    // S = sum_i (a_i * k_hat_i)^2 where a = polarisation, k_hat = propagation
    // For x-polarised, z-propagating: S = 0 (since a_x*k_x = 1*0 = 0, etc.)
    let s: f64 = (0..3).map(|i| {
        let ak = polarisation[i] * propagation[i];
        ak * ak
    }).sum();
    let s = Complex64::from(s);

    // Radiative reaction: ik^3/(6π) in our convention (already includes 1/(4π))
    let rad_reaction = Complex64::new(0.0, k.powi(3) / (6.0 * std::f64::consts::PI));

    // LDR lattice correction (Draine & Goodman 1993).
    // The D&G coefficients are for the α_DG convention where α_DG = α_our/(4π),
    // so the correction must be divided by 4π to convert to our convention.
    let lattice_term = (k * k / d) * (b1 + b2 * epsilon + b3 * epsilon * s)
        / (4.0 * std::f64::consts::PI);

    // α_LDR⁻¹ = α_CM⁻¹ − ik³/(6π) − (k²/d)(b₁ + b₂ε + b₃εS)/(4π)
    let inv_alpha_cm = Complex64::from(1.0) / alpha_cm;
    let inv_alpha_ldr = inv_alpha_cm - rad_reaction - lattice_term;

    Complex64::from(1.0) / inv_alpha_ldr
}
