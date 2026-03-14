//! Substrate image-dipole correction for the CDA.
//!
//! Models the reflection of dipole fields from a planar, non-magnetic substrate
//! at z = z_interface using an image-dipole approximation.
//!
//! # Physical model
//!
//! For a nanostructure above a substrate with dielectric function ε_sub(ω),
//! each real dipole **p**_j at position **r**_j induces an image dipole at the
//! mirror position **r**_j′ = (x_j, y_j, 2z_interface − z_j):
//!
//! **p**_j′ = r · T · **p**_j,  T = diag(−1, −1, +1)
//!
//! where the reflection factor `r` is either:
//!
//! - **Quasi-static** (`use_retarded = false`):
//!   Δε = (ε_sub − ε_env)/(ε_sub + ε_env)
//! - **Retarded** (`use_retarded = true`, default):
//!   r_amp = (n_sub − n_env)/(n_sub + n_env)
//!   where n = √ε is the complex refractive index, evaluated at k∥ = 0
//!   (normal-incidence Fresnel amplitude reflection coefficient).
//!
//! The retarded form accounts for radiation retardation via the full
//! frequency-dependent refractive index rather than the dielectric constant.
//! At low frequencies the two forms are equivalent; they differ significantly
//! for metallic substrates at optical frequencies where Im(ε) is large.
//!
//! The image contribution to the interaction matrix is:
//!
//! A_ij += −r · T · **G**(r_i, r_j′)
//!
//! (including i = j, the self-image term).
//!
//! # Limitations
//!
//! Both forms are single-image approximations (k∥ = 0). Full Sommerfeld
//! integral treatment with k∥-dependent Fresnel coefficients is deferred
//! to v0.5.

use std::sync::Arc;
use num_complex::Complex64;
use serde::{Deserialize, Serialize};

use super::sommerfeld::SommerfeldGreens;

fn default_use_retarded() -> bool { true }

/// Specifies a planar substrate below the nanostructure.
///
/// All dipoles must satisfy `z > z_interface` for the image-dipole formula to
/// be physically valid.
///
/// # TOML example
/// ```toml
/// [simulation.substrate]
/// z_interface = 0.0
/// material    = "SiO2_Palik"
/// use_retarded  = true    # optional, default true
/// use_sommerfeld = false  # optional, default false; overrides use_retarded
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubstrateSpec {
    /// z-coordinate of the substrate–medium interface in nm.
    /// Defaults to 0.0 (interface at the origin).
    #[serde(default)]
    pub z_interface: f64,
    /// Material identifier for the substrate. Resolved via the same material
    /// providers used for particle materials (e.g. `"SiO2_Palik"`).
    pub material: String,
    /// Use the retarded (refractive-index based) Fresnel amplitude coefficient
    /// instead of the quasi-static Δε.  Defaults to `true`.
    ///
    /// Ignored when `use_sommerfeld = true`.
    #[serde(default = "default_use_retarded")]
    pub use_retarded: bool,
    /// Use the full Sommerfeld integral treatment (k∥-dependent Fresnel
    /// coefficients). More accurate than the scalar image-dipole approximation,
    /// especially for particles close to the substrate (z < λ/4).
    ///
    /// Defaults to `false`. When `true`, overrides `use_retarded`.
    #[serde(default)]
    pub use_sommerfeld: bool,
}

/// Runtime substrate correction, constructed once per wavelength.
///
/// Implements `Debug` manually because `SommerfeldGreens` only derives `Debug`
/// and `Arc<SommerfeldGreens>` does too, so the derive would work — but we keep
/// the manual impl explicit for clarity.
///
/// `SubstrateRuntime` is what the solver receives — it holds the
/// wavelength-resolved reflection information rather than material names.
///
/// # Construction
///
/// Call [`SubstrateRuntime::from_spec`] in the runner / GUI per-wavelength loop.
#[derive(Clone, Debug)]
pub enum SubstrateRuntime {
    /// Scalar image-dipole: r · T · G(r_i, r_j′).
    ///
    /// Corresponds to `use_retarded = false` (quasi-static Δε) or
    /// `use_retarded = true` (retarded Fresnel amplitude).
    Fresnel {
        z_interface: f64,
        /// Scalar reflection factor (quasi-static Δε or retarded r_amp).
        delta_eps: Complex64,
    },
    /// Full Sommerfeld integral: **G**_refl(**r**_i, **r**_j) from k∥ integration.
    Sommerfeld {
        z_interface: f64,
        /// Pre-initialised integrator (cheap to clone — contains only a few f64/Complex64 fields).
        greens: Arc<SommerfeldGreens>,
    },
}

impl SubstrateRuntime {
    /// Build the runtime substrate from a spec and the wavelength-resolved
    /// substrate dielectric function.
    ///
    /// # Arguments
    /// * `spec`      — The substrate configuration.
    /// * `eps_sub`   — Substrate ε at the current wavelength.
    /// * `eps_env`   — Environment ε = n_env² at the current wavelength.
    /// * `wavelength_nm` — Current wavelength in nm (needed for Sommerfeld).
    pub fn from_spec(spec: &SubstrateSpec, eps_sub: Complex64, eps_env: f64, wavelength_nm: f64) -> Self {
        if spec.use_sommerfeld {
            SubstrateRuntime::Sommerfeld {
                z_interface: spec.z_interface,
                greens: Arc::new(SommerfeldGreens::new(wavelength_nm, eps_sub, eps_env)),
            }
        } else {
            SubstrateRuntime::Fresnel {
                z_interface: spec.z_interface,
                delta_eps: substrate_reflection_factor(eps_sub, eps_env, spec.use_retarded),
            }
        }
    }

    /// The z-coordinate of the interface, regardless of variant.
    pub fn z_interface(&self) -> f64 {
        match self {
            SubstrateRuntime::Fresnel { z_interface, .. } => *z_interface,
            SubstrateRuntime::Sommerfeld { z_interface, .. } => *z_interface,
        }
    }
}

/// Compute the quasi-static Fresnel reflection factor.
///
/// Δε = (ε_sub − ε_env) / (ε_sub + ε_env)
///
/// This is the quasi-static (long-wavelength) limit valid when the
/// nanostructure is much smaller than the wavelength. For metallic substrates
/// at optical frequencies, prefer [`fresnel_reflection_amplitude`] which
/// accounts for retardation.
///
/// # Arguments
/// * `eps_sub`  — Substrate dielectric function at the current wavelength.
/// * `eps_env`  — Environment dielectric constant ε_env = n_env².
pub fn fresnel_delta_eps(eps_sub: Complex64, eps_env: f64) -> Complex64 {
    let env = Complex64::from(eps_env);
    (eps_sub - env) / (eps_sub + env)
}

/// Compute the retarded normal-incidence Fresnel amplitude reflection coefficient.
///
/// r_amp = (n_sub − n_env) / (n_sub + n_env)
///
/// where n_sub = √ε_sub (principal complex square root) and
/// n_env = √ε_env (real, positive).
///
/// This is the retarded version of the image-dipole reflection factor. It
/// reduces to `fresnel_delta_eps` in the quasi-static limit but correctly
/// accounts for retardation at optical frequencies via the complex refractive
/// index. For a lossless dielectric (ε_sub real, positive) this equals the
/// standard Fresnel formula for normal incidence.
///
/// # Arguments
/// * `eps_sub`  — Substrate dielectric function at the current wavelength.
/// * `eps_env`  — Environment dielectric constant ε_env = n_env².
pub fn fresnel_reflection_amplitude(eps_sub: Complex64, eps_env: f64) -> Complex64 {
    let n_sub = eps_sub.sqrt(); // principal branch: Re(sqrt) >= 0
    let n_env = Complex64::from(eps_env.sqrt());
    (n_sub - n_env) / (n_sub + n_env)
}

/// Select the appropriate reflection factor based on `use_retarded`.
///
/// Convenience wrapper used by the CLI runner and GUI app so that call sites
/// only need to change the `SubstrateSpec` field rather than their own logic.
pub fn substrate_reflection_factor(eps_sub: Complex64, eps_env: f64, use_retarded: bool) -> Complex64 {
    if use_retarded {
        fresnel_reflection_amplitude(eps_sub, eps_env)
    } else {
        fresnel_delta_eps(eps_sub, eps_env)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fresnel_delta_eps_glass_in_air() {
        // SiO2 in air: ε_sub ≈ 2.25, ε_env = 1.0 → Δε = 1.25/3.25 ≈ 0.385
        let eps_sub = Complex64::new(2.25, 0.0);
        let delta = fresnel_delta_eps(eps_sub, 1.0);
        let expected = 1.25 / 3.25;
        assert!(
            (delta.re - expected).abs() < 1e-10,
            "Δε.re = {:.4}, expected {:.4}",
            delta.re, expected
        );
        assert!(delta.im.abs() < 1e-12, "Δε for lossless dielectric should be real");
    }

    #[test]
    fn test_fresnel_delta_eps_same_medium() {
        // Substrate and environment identical → Δε = 0
        let eps = Complex64::new(2.0, 0.0);
        let delta = fresnel_delta_eps(eps, 2.0);
        assert!(delta.norm() < 1e-12, "Δε should be 0 when substrate = environment");
    }

    #[test]
    fn test_fresnel_delta_eps_metallic_substrate() {
        // Metal substrate: large |ε_sub| → Δε → 1
        let eps_sub = Complex64::new(-100.0, 5.0); // metal-like
        let delta = fresnel_delta_eps(eps_sub, 1.0);
        // As ε_sub → ∞, Δε → 1
        assert!(delta.norm() > 0.95, "metallic substrate should give |Δε| close to 1");
    }

    #[test]
    fn test_fresnel_reflection_amplitude_glass_in_air() {
        // SiO2 in air: n_sub = 1.5, n_env = 1.0 → r = 0.5/2.5 = 0.2
        let eps_sub = Complex64::new(2.25, 0.0); // n = 1.5
        let r = fresnel_reflection_amplitude(eps_sub, 1.0);
        let expected = (1.5 - 1.0) / (1.5 + 1.0); // 0.2
        assert!(
            (r.re - expected).abs() < 1e-10,
            "r_amp.re = {:.6}, expected {:.6}",
            r.re, expected
        );
        assert!(r.im.abs() < 1e-12, "r_amp for lossless dielectric should be real");
    }

    #[test]
    fn test_fresnel_reflection_amplitude_same_medium() {
        // n_sub == n_env → r = 0
        let eps = Complex64::new(2.25, 0.0);
        let r = fresnel_reflection_amplitude(eps, 2.25);
        assert!(r.norm() < 1e-12, "r_amp should be 0 when substrate = environment");
    }

    #[test]
    fn test_fresnel_reflection_amplitude_differs_from_delta_eps() {
        // For a dielectric the two factors differ: r_amp uses n, Δε uses ε.
        // n=1.5, ε=2.25: r_amp=0.2, Δε=0.385
        let eps_sub = Complex64::new(2.25, 0.0);
        let r_amp = fresnel_reflection_amplitude(eps_sub, 1.0);
        let delta = fresnel_delta_eps(eps_sub, 1.0);
        assert!(
            (r_amp.re - 0.2).abs() < 1e-10,
            "r_amp = {:.6}, expected 0.2", r_amp.re
        );
        assert!(
            (delta.re - 1.25 / 3.25).abs() < 1e-10,
            "delta_eps = {:.6}, expected {:.6}", delta.re, 1.25_f64 / 3.25_f64
        );
        // They should be different
        assert!(
            (r_amp.re - delta.re).abs() > 0.1,
            "r_amp and Δε should differ for dielectrics"
        );
    }

    #[test]
    fn test_substrate_reflection_factor_dispatch() {
        let eps_sub = Complex64::new(2.25, 0.0);
        let eps_env = 1.0_f64;
        let r_retarded = substrate_reflection_factor(eps_sub, eps_env, true);
        let r_quasi = substrate_reflection_factor(eps_sub, eps_env, false);
        assert_eq!(r_retarded, fresnel_reflection_amplitude(eps_sub, eps_env));
        assert_eq!(r_quasi, fresnel_delta_eps(eps_sub, eps_env));
    }

    #[test]
    fn test_substrate_spec_defaults() {
        // Check that the serde defaults are correct
        let spec: SubstrateSpec = serde_json::from_str(
            r#"{"z_interface": 0.0, "material": "SiO2_Palik"}"#
        ).unwrap();
        assert!(spec.use_retarded, "use_retarded should default to true");
        assert!(!spec.use_sommerfeld, "use_sommerfeld should default to false");
    }

    #[test]
    fn test_substrate_runtime_fresnel() {
        let spec = SubstrateSpec {
            z_interface: 0.0,
            material: "SiO2_Palik".into(),
            use_retarded: true,
            use_sommerfeld: false,
        };
        let eps_sub = Complex64::new(2.25, 0.0);
        let rt = SubstrateRuntime::from_spec(&spec, eps_sub, 1.0, 600.0);
        match rt {
            SubstrateRuntime::Fresnel { delta_eps, .. } => {
                // n_sub = 1.5, n_env = 1.0 → r = 0.5/2.5 = 0.2
                assert!((delta_eps.re - 0.2).abs() < 1e-10,
                    "Fresnel runtime delta_eps.re = {:.6}", delta_eps.re);
            }
            SubstrateRuntime::Sommerfeld { .. } => panic!("Expected Fresnel variant"),
        }
    }

    #[test]
    fn test_substrate_runtime_sommerfeld() {
        let spec = SubstrateSpec {
            z_interface: 5.0,
            material: "SiO2_Palik".into(),
            use_retarded: true,
            use_sommerfeld: true,
        };
        let eps_sub = Complex64::new(2.25, 0.0);
        let rt = SubstrateRuntime::from_spec(&spec, eps_sub, 1.0, 600.0);
        match rt {
            SubstrateRuntime::Sommerfeld { z_interface, .. } => {
                assert_eq!(z_interface, 5.0);
            }
            SubstrateRuntime::Fresnel { .. } => panic!("Expected Sommerfeld variant"),
        }
    }
}
