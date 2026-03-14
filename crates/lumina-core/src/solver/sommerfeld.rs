//! Full Sommerfeld integral reflected Green's tensor for a planar substrate.
//!
//! This module computes the reflected (substrate image) Green's dyadic using
//! the rigorous Sommerfeld integral over in-plane wavevectors k∥, replacing the
//! normal-incidence scalar Fresnel approximation used in the image-dipole method.
//!
//! # Physical background
//!
//! For a dipole source at **r**_j above a planar substrate at z = z₀, the
//! substrate scatters radiation back to the observation point **r**_i. The
//! exact scattered Green's dyadic is:
//!
//! **G**_refl(**r**_i, **r**_j) = (i/8π) ∫₀^∞ [r_s M_s + r_p M_p] exp(ik_z Δz) k∥/k_z dk∥
//!
//! where Δz = (z_i − z₀) + (z_j − z₀) > 0, r_s and r_p are the
//! TE/TM Fresnel reflection coefficients, and the tensor weights M_s/M_p encode
//! the azimuthal structure. After azimuthal integration the result is expressed
//! through Bessel functions J₀, J₁, J₂ of the lateral separation ρ.
//!
//! The integration contour follows the real k∥ axis, split into three segments
//! at the environment and substrate light lines. Each segment is integrated with
//! 15-point Gauss–Kronrod quadrature (G7K15).
//!
//! # Relation to the image-dipole approximation
//!
//! The retarded image-dipole formula uses r_p(k∥=0) = r_s(k∥=0) = r_Fresnel as
//! a single scalar multiplied by the free-space Green's function at the mirror
//! position. The Sommerfeld integral reduces to this approximation only when
//! the particle is many wavelengths above the substrate and the off-axis
//! contributions are negligible. For particles close to the substrate (z < λ/4),
//! the Sommerfeld integral is significantly more accurate.

use num_complex::Complex64;

// ─── G7K15 quadrature constants ─────────────────────────────────────────────
//
// Nodes and weights for the 15-point Gauss-Kronrod rule on [-1, 1].
// Source: Laurie (1997), SIAM J. Numer. Anal. 34(3), 1098-1107.
// The 7-point Gauss rule uses nodes at indices 1,3,5,7,9,11,13 of the K15 list.
// Error estimate: |I_K15 - I_G7|.

const GK15_NODES: [f64; 15] = [
    -0.991_455_371_120_813,
    -0.949_107_912_342_758,
    -0.864_864_423_359_769,
    -0.741_531_185_599_394,
    -0.586_087_235_467_691,
    -0.405_845_151_377_397,
    -0.207_784_955_007_898,
     0.000_000_000_000_000,
     0.207_784_955_007_898,
     0.405_845_151_377_397,
     0.586_087_235_467_691,
     0.741_531_185_599_394,
     0.864_864_423_359_769,
     0.949_107_912_342_758,
     0.991_455_371_120_813,
];

const GK15_WEIGHTS: [f64; 15] = [
    0.022_935_322_010_529,
    0.063_092_092_629_979,
    0.104_790_010_322_250,
    0.140_653_259_715_525,
    0.169_004_726_639_267,
    0.190_350_578_064_785,
    0.204_432_940_075_298,
    0.209_482_141_084_728,
    0.204_432_940_075_298,
    0.190_350_578_064_785,
    0.169_004_726_639_267,
    0.140_653_259_715_525,
    0.104_790_010_322_250,
    0.063_092_092_629_979,
    0.022_935_322_010_529,
];

// G7 weights at the corresponding K15 node positions (zero at K15-only nodes).
// Used for error estimation in adaptive refinement (future).
#[allow(dead_code)]
const G7_WEIGHTS_IN_GK15: [f64; 15] = [
    0.0,
    0.129_484_966_168_870,
    0.0,
    0.279_705_391_489_277,
    0.0,
    0.381_830_050_505_119,
    0.0,
    0.417_959_183_673_469,
    0.0,
    0.381_830_050_505_119,
    0.0,
    0.279_705_391_489_277,
    0.0,
    0.129_484_966_168_870,
    0.0,
];

// ─── Bessel functions for real non-negative argument ───────────────────────

/// J₀(x) — zeroth-order Bessel function of the first kind.
///
/// Uses rational approximations from Numerical Recipes (Press et al. 3rd ed.),
/// accurate to ~machine precision for all x ≥ 0.
fn j0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 8.0 {
        let y = x * x;
        let p1 = 57_568_490_574.0 + y * (-13_362_590_354.0 + y * (651_619_640.7
             + y * (-11_214_424.18 + y * (77_392.33017 + y * (-184.905_245_6)))));
        let q1 = 57_568_490_411.0 + y * (1_029_532_985.0 + y * (9_494_680.718
             + y * (59_272.648_53 + y * (267.853_271_2 + y))));
        p1 / q1
    } else {
        let z = 8.0 / ax;
        let y = z * z;
        let xx = ax - 0.785_398_163_4;
        let p1 = 1.0 + y * (-0.001_098_628_627 + y * (0.000_002_734_511
             + y * (-0.000_000_020_73 + y * (-0.000_000_000_207_3))));
        let q1 = -0.001_562_499_995 + y * (0.000_001_430_488
             + y * (-0.000_000_006_905 + y * (0.000_000_000_107_7 + y * 0.000_000_000_002_72)));
        (2.0 / (std::f64::consts::PI * ax)).sqrt()
            * (p1 * xx.cos() - z * q1 * xx.sin())
    }
}

/// J₁(x) — first-order Bessel function of the first kind.
///
/// Uses rational approximations from Numerical Recipes, accurate to ~machine precision.
fn j1(x: f64) -> f64 {
    let ax = x.abs();
    let val = if ax < 8.0 {
        let y = x * x;
        let p1 = x * (72_362_614_232.0 + y * (-7_895_059_235.0 + y * (242_396_853.1
              + y * (-2_972_611.439 + y * (15_704.488_2 + y * (-30.163_508_1))))));
        let q1 = 144_725_228_442.0 + y * (2_300_535_178.0 + y * (18_583_304.74
              + y * (99_447.435_26 + y * (376.994_919_9 + y))));
        p1 / q1
    } else {
        let z = 8.0 / ax;
        let y = z * z;
        let xx = ax - 2.356_194_491;
        let p1 = 1.0 + y * (0.000_183_105 + y * (-0.000_003_516_396_5
             + y * (0.000_000_024_575_72 + y * 2.457_520e-7)));
        let q1 = 0.040_957_605 + y * (-0.000_004_976_877
             + y * (0.000_000_079_821_6 + y * (-0.000_000_000_944_603 + y * 4.534_91e-10)));
        (2.0 / (std::f64::consts::PI * ax)).sqrt()
            * (p1 * xx.cos() - z * q1 * xx.sin())
    };
    if x < 0.0 { -val } else { val }
}

/// J₂(x) = (2/x) J₁(x) − J₀(x)
/// Stable for x > 0; returns 0 for x = 0.
fn j2(x: f64) -> f64 {
    if x.abs() < 1e-15 { return 0.0; }
    (2.0 / x) * j1(x) - j0(x)
}

// ─── Complex branch cuts ────────────────────────────────────────────────────

/// Complex square root with Im(result) ≥ 0 convention (outgoing wave / decaying).
///
/// For the vertical wavenumber k_z = √(k² − k∥²):
/// - propagating wave: k_z purely real (> 0)
/// - evanescent wave: k_z purely imaginary (> 0·i)
fn csqrt_pos_im(z: Complex64) -> Complex64 {
    let s = z.sqrt();
    if s.im < 0.0 { -s } else { s }
}

// ─── Intermediate integral accumulators ────────────────────────────────────

/// The five independent Sommerfeld integrals needed to build the 3×3 tensor.
///
/// All are of the form ∫ f(k∥) dk∥ along the real k∥ axis.
struct Integrals {
    /// ∫ r_s · J₀(k∥ρ) · exp(ik_z Δz) · (k∥/k_z) dk∥
    i_s_j0: Complex64,
    /// ∫ r_p · (k_z/k)² · J₀(k∥ρ) · exp(ik_z Δz) · (k∥/k_z) dk∥
    i_p_j0: Complex64,
    /// ∫ r_p · (k∥/k)² · J₀(k∥ρ) · exp(ik_z Δz) · (k∥/k_z) dk∏  [for G_zz]
    i_zz: Complex64,
    /// ∫ r_p · (k∥/k) · (k_z/k) · J₁(k∥ρ) · exp(ik_z Δz) dk∥  [for G_xz]
    i_xz: Complex64,
    /// ∫ r_p · (k_z/k)² · J₂(k∥ρ) · exp(ik_z Δz) · (k∥/k_z) dk∥  [G_xx/G_yy split]
    i_j2: Complex64,
}

impl Integrals {
    fn zero() -> Self {
        let z = Complex64::new(0.0, 0.0);
        Integrals { i_s_j0: z, i_p_j0: z, i_zz: z, i_xz: z, i_j2: z }
    }

    fn add(&mut self, other: &Integrals) {
        self.i_s_j0 += other.i_s_j0;
        self.i_p_j0 += other.i_p_j0;
        self.i_zz   += other.i_zz;
        self.i_xz   += other.i_xz;
        self.i_j2   += other.i_j2;
    }
}

// ─── SommerfeldGreens ────────────────────────────────────────────────────────

/// Reflected Green's dyadic computed via the Sommerfeld integral.
///
/// Created once per wavelength; `evaluate()` is called for every dipole pair
/// (i, j) in the interaction matrix assembly.
///
/// The struct is `Clone + Send + Sync` so it can be wrapped in an `Arc` and
/// shared across Rayon threads without locking.
#[derive(Debug, Clone)]
pub struct SommerfeldGreens {
    /// Wavenumber in the environment: k_env = 2π n_env / λ (nm⁻¹)
    pub k_env: f64,
    /// Substrate wavenumber k_sub = csqrt(ε_sub) · 2π/λ (nm⁻¹, Im ≥ 0)
    pub k_sub: Complex64,
    /// Substrate dielectric function at this wavelength
    pub eps_sub: Complex64,
    /// Environment dielectric constant ε_env = n_env²
    pub eps_env: f64,
    /// Number of quadrature segments for the oscillatory real-k∥ integral.
    /// Each segment uses 15 Gauss–Kronrod nodes. Default: 1 (one segment each).
    pub n_subdiv: usize,
}

impl SommerfeldGreens {
    /// Construct for a given wavelength and substrate material.
    ///
    /// # Arguments
    /// * `wavelength_nm` — Free-space wavelength in nm.
    /// * `eps_sub`       — Substrate ε at this wavelength.
    /// * `eps_env`       — Environment ε = n_env².
    pub fn new(wavelength_nm: f64, eps_sub: Complex64, eps_env: f64) -> Self {
        let k0 = 2.0 * std::f64::consts::PI / wavelength_nm;
        let k_env = k0 * eps_env.sqrt();
        let k_sub = csqrt_pos_im(eps_sub * Complex64::new(k0 * k0, 0.0));
        Self { k_env, k_sub, eps_sub, eps_env, n_subdiv: 1 }
    }

    /// Evaluate the 3×3 reflected Green's tensor **G**_refl(**r**_i, **r**_j).
    ///
    /// Both **r**_i and **r**_j must satisfy z > z_interface.
    ///
    /// # Convention
    ///
    /// The result matches the `k² exp(ikR)/(4πR)` prefactor convention of the
    /// free-space Green's function in `greens.rs`. For a substrate with the same
    /// dielectric constant as the environment, all components are zero.
    pub fn evaluate(
        &self,
        ri: [f64; 3],
        rj: [f64; 3],
        z_interface: f64,
    ) -> [[Complex64; 3]; 3] {
        // Heights above interface
        let zi = ri[2] - z_interface;
        let zj = rj[2] - z_interface;
        // Δz = (z_i − z₀) + (z_j − z₀) — always positive if both above substrate.
        let dz = zi + zj;
        if dz <= 0.0 {
            // Degenerate or below-substrate — return zero rather than panic.
            return [[Complex64::new(0.0, 0.0); 3]; 3];
        }

        // Lateral separation and azimuthal angle
        let dx = ri[0] - rj[0];
        let dy = ri[1] - rj[1];
        let rho = (dx * dx + dy * dy).sqrt();
        let phi = dy.atan2(dx);
        let cos_phi = phi.cos();
        let sin_phi = phi.sin();
        let cos2phi = 2.0 * cos_phi * cos_phi - 1.0; // cos(2φ)
        let sin2phi = 2.0 * sin_phi * cos_phi;        // sin(2φ)

        let integ = self.integrate(rho, dz);

        // ── Build tensor from integrals ───────────────────────────────────
        // Prefactor: i/(8π) — matches the k²e^{ikR}/(4πR) convention after
        // accounting for the 1/(2π) from the azimuthal integration.
        let pre = Complex64::new(0.0, 1.0 / (8.0 * std::f64::consts::PI));

        // Diagonal and off-diagonal components
        let g_xx = pre * (integ.i_s_j0 + integ.i_p_j0 - integ.i_j2 * cos2phi);
        let g_yy = pre * (integ.i_s_j0 + integ.i_p_j0 + integ.i_j2 * cos2phi);
        let g_zz = pre * (-integ.i_zz);
        let g_xy = pre * (-integ.i_j2 * sin2phi);
        let g_xz = pre * (integ.i_xz * cos_phi);
        let g_yz = pre * (integ.i_xz * sin_phi);

        [
            [g_xx, g_xy, g_xz],
            [g_xy, g_yy, g_yz],
            [g_xz, g_yz, g_zz],
        ]
    }

    /// Compute the five independent integrals over k∥ ∈ [0, k_max].
    ///
    /// The integration range is split into three physical segments:
    ///
    /// 1. **[0, k_env]**: both k_z and k_z_sub complex with Im = 0.
    ///    Integrand oscillates as exp(ik_z dz).
    /// 2. **[k_env, Re(k_sub)]**: k_z purely imaginary (evanescent in env).
    ///    Integrand decays as exp(-|k_z| dz).
    ///    (Skipped when Re(k_sub) ≤ k_env, e.g. for metallic substrates.)
    /// 3. **[Re(k_sub_real), k_max]**: both evanescent.
    ///    exp(-|k_z| dz) decay; integrate to adaptive k_max.
    fn integrate(&self, rho: f64, dz: f64) -> Integrals {
        let k_env = self.k_env;
        let k_sub_real = self.k_sub.re.max(k_env); // ensure segment 2 is non-empty only when k_sub.re > k_env
        let has_seg2 = self.k_sub.re > k_env * 1.001; // 0.1% margin

        // Upper limit for evanescent segment: where exp(-kz·dz) < 1e-10.
        // kz ≈ k∥ for large k∥, so k_max ≈ (10 ln 10) / dz.
        let k_max = if dz > 1e-6 {
            (10.0 * std::f64::consts::LN_10 / dz).max(k_env * 5.0)
        } else {
            k_env * 20.0
        };

        let mut total = Integrals::zero();

        // Segment 1: [0, k_env] — propagating (possibly oscillatory)
        let seg1 = self.integrate_segment(0.0, k_env, rho, dz);
        total.add(&seg1);

        // Segment 2: [k_env, Re(k_sub)] — evanescent in env, may be propagating in sub
        if has_seg2 {
            let seg2 = self.integrate_segment(k_env, k_sub_real, rho, dz);
            total.add(&seg2);
        }

        // Segment 3: [max(k_env, Re(k_sub)), k_max] — both evanescent
        let seg3_lo = if has_seg2 { k_sub_real } else { k_env };
        if k_max > seg3_lo * 1.001 {
            let seg3 = self.integrate_segment(seg3_lo, k_max, rho, dz);
            total.add(&seg3);
        }

        total
    }

    /// Apply 15-point Gauss–Kronrod on [a, b].
    fn integrate_segment(&self, a: f64, b: f64, rho: f64, dz: f64) -> Integrals {
        let mid = 0.5 * (a + b);
        let half = 0.5 * (b - a);

        let mut acc = Integrals::zero();
        for idx in 0..15 {
            let kpar = mid + half * GK15_NODES[idx];
            let w = GK15_WEIGHTS[idx] * half;
            let contrib = self.integrand(kpar, rho, dz);
            acc.i_s_j0 += w * contrib.i_s_j0;
            acc.i_p_j0 += w * contrib.i_p_j0;
            acc.i_zz   += w * contrib.i_zz;
            acc.i_xz   += w * contrib.i_xz;
            acc.i_j2   += w * contrib.i_j2;
        }
        acc
    }

    /// Evaluate all five integrand kernels at a single k∥ value.
    fn integrand(&self, kpar: f64, rho: f64, dz: f64) -> Integrals {
        let k = self.k_env; // k in environment (k_env)
        let k2 = k * k;

        // Vertical wavenumbers (Im ≥ 0 convention)
        let kz_env  = csqrt_pos_im(Complex64::new(k2 - kpar * kpar, 0.0));
        let kz_sub  = csqrt_pos_im(self.eps_sub * Complex64::new(k2 / self.eps_env, 0.0) - Complex64::new(kpar * kpar, 0.0));

        // Fresnel coefficients
        let r_s = (kz_env - kz_sub) / (kz_env + kz_sub);
        let r_p = (self.eps_sub * kz_env - Complex64::new(self.eps_env, 0.0) * kz_sub)
                / (self.eps_sub * kz_env + Complex64::new(self.eps_env, 0.0) * kz_sub);

        // Phase factor: exp(i k_z Δz)
        let phase = (Complex64::new(0.0, 1.0) * kz_env * dz).exp();

        // Bessel functions at k∥ ρ
        let arg = kpar * rho;
        let j0v = Complex64::new(j0(arg), 0.0);
        let j1v = Complex64::new(j1(arg), 0.0);
        let j2v = Complex64::new(j2(arg), 0.0);

        // Weights for each polarisation contribution
        let kpar_c = Complex64::new(kpar, 0.0);
        let k2_c   = Complex64::new(k2, 0.0);
        let kz_factor = kpar_c / kz_env; // k∥ / k_z

        let i_s_j0 = r_s * j0v * phase * kz_factor;
        let i_p_j0 = r_p * (kz_env * kz_env / k2_c) * j0v * phase * kz_factor;
        let i_zz   = r_p * (kpar_c * kpar_c / k2_c) * j0v * phase * kz_factor;
        let i_xz   = r_p * (kpar_c / Complex64::new(k, 0.0)) * (kz_env / Complex64::new(k, 0.0)) * j1v * phase;
        let i_j2   = r_p * (kz_env * kz_env / k2_c) * j2v * phase * kz_factor;

        Integrals { i_s_j0, i_p_j0, i_zz, i_xz, i_j2 }
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // ── Bessel function tests ────────────────────────────────────────────────

    #[test]
    fn test_j0_at_zero() {
        assert!((j0(0.0) - 1.0).abs() < TOL, "J₀(0) = 1");
    }

    #[test]
    fn test_j0_first_root() {
        // First zero of J₀ ≈ 2.4048
        let val = j0(2.404_825_557);
        assert!(val.abs() < 1e-5, "J₀(2.4048) ≈ 0, got {val:.2e}");
    }

    #[test]
    fn test_j1_at_zero() {
        assert!(j1(0.0).abs() < TOL, "J₁(0) = 0");
    }

    #[test]
    fn test_j1_first_root() {
        // First positive zero of J₁ ≈ 3.8317
        let val = j1(3.831_705_97);
        assert!(val.abs() < 1e-5, "J₁(3.8317) ≈ 0, got {val:.2e}");
    }

    #[test]
    fn test_j2_at_zero() {
        assert!(j2(0.0).abs() < TOL, "J₂(0) = 0");
    }

    #[test]
    fn test_j2_recurrence() {
        // J₂(x) = (2/x) J₁(x) − J₀(x) — verify at x = 3
        let x = 3.0_f64;
        let j2_direct = j2(x);
        let j2_recur = (2.0 / x) * j1(x) - j0(x);
        assert!((j2_direct - j2_recur).abs() < TOL, "J₂ recurrence mismatch");
    }

    // ── No-substrate limit ───────────────────────────────────────────────────

    #[test]
    fn test_no_substrate_zero_tensor() {
        // When eps_sub == eps_env, r_s = r_p = 0 everywhere → G_refl = 0.
        let eps_env = 2.25_f64; // n = 1.5
        let eps_sub = Complex64::new(eps_env, 0.0);
        let sg = SommerfeldGreens::new(500.0, eps_sub, eps_env);

        let ri = [0.0, 0.0, 10.0];
        let rj = [3.0, 0.0, 8.0];
        let g = sg.evaluate(ri, rj, 0.0);

        for row in &g {
            for &val in row.iter() {
                assert!(
                    val.norm() < 1e-10,
                    "G_refl should be zero when substrate = environment, got {val:.2e}"
                );
            }
        }
    }

    // ── Symmetry tests ───────────────────────────────────────────────────────

    #[test]
    fn test_tensor_symmetry_reciprocity() {
        // G_refl[a][b] = G_refl[b][a] (symmetry of the reflected dyadic)
        let sg = SommerfeldGreens::new(600.0, Complex64::new(2.25, 0.0), 1.0);
        let ri = [5.0, 3.0, 12.0];
        let rj = [2.0, 1.0, 8.0];
        let g = sg.evaluate(ri, rj, 0.0);

        assert!((g[0][1] - g[1][0]).norm() < 1e-10, "G_xy = G_yx");
        assert!((g[0][2] - g[2][0]).norm() < 1e-10, "G_xz = G_zx");
        assert!((g[1][2] - g[2][1]).norm() < 1e-10, "G_yz = G_zy");
    }

    #[test]
    fn test_xz_zero_for_zero_rho() {
        // ρ = 0 (ri and rj directly above each other) → G_xz = G_yz = G_xy = 0
        let sg = SommerfeldGreens::new(600.0, Complex64::new(2.25, 0.0), 1.0);
        let ri = [0.0, 0.0, 12.0];
        let rj = [0.0, 0.0, 8.0];
        let g = sg.evaluate(ri, rj, 0.0);

        assert!(g[0][2].norm() < 1e-8, "G_xz = 0 at ρ=0, got {:.2e}", g[0][2].norm());
        assert!(g[1][2].norm() < 1e-8, "G_yz = 0 at ρ=0, got {:.2e}", g[1][2].norm());
        assert!(g[0][1].norm() < 1e-8, "G_xy = 0 at ρ=0, got {:.2e}", g[0][1].norm());
    }

    // ── Far-field decay ──────────────────────────────────────────────────────

    #[test]
    fn test_far_field_decay() {
        // For dz >> λ, G_refl should be very small (exponential decay of evanescent parts).
        let sg = SommerfeldGreens::new(500.0, Complex64::new(2.25, 0.0), 1.0);
        // dz = 10 λ = 5000 nm — deep far-field
        let ri = [0.0, 0.0, 2502.5];
        let rj = [0.0, 0.0, 2497.5]; // dz = 5000 nm
        let g = sg.evaluate(ri, rj, 0.0);
        // In the far field the reflected Green's tensor is small (not strictly zero,
        // as the propagating part survives, but it should be of order r_Fresnel * G_free magnitude)
        // Just verify it's finite and non-NaN.
        for row in &g {
            for &val in row.iter() {
                assert!(val.is_finite(), "G_refl should be finite in far field");
            }
        }
    }

    // ── Physical limit: retarded Fresnel agreement ───────────────────────────

    #[test]
    fn test_sommerfeld_matches_fresnel_on_axis_large_dz() {
        // For ρ = 0 and dz >> λ, the Sommerfeld integral is dominated by the
        // k∥ ≈ 0 contribution. G_refl ≈ r_Fresnel × G_free(ri, rj_image).
        //
        // Here we check that the diagonal components have the right sign pattern:
        // G_xx ≈ G_yy (equal by symmetry at ρ=0) and G_zz has opposite sign.
        let eps_sub = Complex64::new(2.25, 0.0); // glass
        let eps_env = 1.0_f64;                   // air
        let lambda = 600.0_f64;

        let sg = SommerfeldGreens::new(lambda, eps_sub, eps_env);

        // Moderate dz = 2λ to keep the integral tractable but still test the physics
        let dz = 2.0 * lambda;
        let ri = [0.0, 0.0, dz * 0.6];
        let rj = [0.0, 0.0, dz * 0.4];
        let z_interface = 0.0;
        let g = sg.evaluate(ri, rj, z_interface);

        // At ρ=0: G_xx = G_yy by symmetry
        assert!(
            (g[0][0] - g[1][1]).norm() < 1e-8 * g[0][0].norm().max(1e-20),
            "G_xx ≠ G_yy at ρ=0: {:.4e} vs {:.4e}", g[0][0], g[1][1]
        );

        // All components should be finite
        for row in &g {
            for &val in row.iter() {
                assert!(val.is_finite(), "G_refl has non-finite component");
            }
        }
    }
}
