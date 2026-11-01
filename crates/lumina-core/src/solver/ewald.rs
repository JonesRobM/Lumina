//! Ewald-accelerated periodic dyadic Green's function.
//!
//! For a nanoparticle array that is periodic with lattice vectors **a1** (and
//! optionally **a2**), the interaction between a dipole at **r**_i and all
//! periodic images of **r**_j is:
//!
//! $$\mathbf{G}_\text{per}(\mathbf{r}, \mathbf{k}_\parallel)
//!   = \sum_{\mathbf{R}\in\Lambda}
//!     \mathbf{G}_\text{free}(\mathbf{r}-\mathbf{R})\,
//!     e^{i\mathbf{k}_\parallel\cdot\mathbf{R}}$$
//!
//! # Implementation strategy (v0.5)
//!
//! The sum is split using the Ewald identity into a **real-space** part
//! (erfc-damped, converges exponentially in real space) and a
//! **reciprocal-space** part (spectral G-sum, converges exponentially in
//! reciprocal space):
//!
//! $$\mathbf{G}_\text{per} = \mathbf{G}_\text{real} + \mathbf{G}_\text{recip}$$
//!
//! ## Real-space part
//!
//! Each free-space G term is multiplied by `erfc(eta * |r - R|)` which strongly
//! damps contributions from distant lattice images.  The splitting parameter
//! eta = sqrt(pi / A_cell) is chosen so that both sums converge in roughly the same
//! number of shells.
//!
//! ## Reciprocal-space part (2D lattices)
//!
//! For a 2D lattice in the xy-plane the spectral sum over reciprocal vectors
//! **G** is (Capolino et al. 2005 / Linton 2010 conventions, consistent with
//! the free-space prefactor k^2exp(ikR)/(4piR)):
//!
//! $$\mathbf{G}^\text{spec}_{ab} = \frac{i}{2A} \sum_\mathbf{G}
//!     \frac{e^{i\mathbf{Q}\cdot\boldsymbol{\rho}}\,e^{-q|z|}}{q}\,
//!     M_{ab}(\mathbf{Q}, q)$$
//!
//! where **Q** = **k**|| + **G**, q = sqrt(Q^2 - k^2) with Im(q) >= 0, and
//!
//! $$M_{ab} = k^2\delta_{ab} - \tilde{k}_a\tilde{k}_b,\quad
//!   \tilde{\mathbf{k}} = (Q_x,\,Q_y,\,i\,q\,\operatorname{sgn}(z))$$
//!
//! ## Reciprocal-space part (1D chains)
//!
//! For a 1D chain along **a1** the spectral sum uses modified Bessel functions
//! K0 and K1. Each reciprocal mode m contributes:
//!
//! $$\mathbf{G}^\text{spec,1D}_{ab} = \frac{1}{2\pi a_1}
//!   \sum_m e^{i g_m x_\parallel}
//!   \bigl[C_\delta\,\delta_{ab} + C_{ee}\,\hat{e}_a\hat{e}_b
//!        + C_{e\rho}(\hat{e}_a\hat{\rho}_b + \hat{e}_b\hat{\rho}_a)
//!        + C_{\rho\rho}\,\hat{\rho}_a\hat{\rho}_b\bigr]$$
//!
//! where g_m = k|| + 2*pi*m/a1, kappa_m = sqrt(g_m^2 - k^2) with Im(kappa) >= 0, and the
//! coefficients involve K0(kapparho) and K1(kapparho).
//!
//! # Convergence note
//!
//! With `truncation_real = 5` and `truncation_recip = 5` the combined
//! Ewald sum typically achieves < 0.01 % error for lattice periods >= 50 nm
//! and visible wavelengths.

use num_complex::Complex64;

use lumina_geometry::lattice::LatticeSpec;

/// Ewald-accelerated periodic dyadic Green's function.
///
/// Wraps a [`LatticeSpec`] and evaluates the Bloch-periodic sum via the
/// real-space erfc-damped sum plus the spectral (reciprocal-space) sum.
pub struct EwaldGreens {
    pub lattice: LatticeSpec,
}

/// Transverse separation below which a 1D pair is treated as collinear.
///
/// Matches the rho guard in [`recip_space_1d`], so exactly one of the two routines
/// contributes for any given displacement.
pub(crate) const COLLINEAR_RHO_TOL: f64 = 1e-10;

/// Lattice images summed either side for the collinear direct sum.
///
/// Convergence is only ~1/N -- terms fall off as 1/R with an oscillating phase -- so
/// this is a cost/accuracy trade-off, not a value that can simply be raised.
/// Measured against the spectral sum (a = 100 nm, lambda = 600 nm, per pair):
///
/// | shells | rel. error | time    |
/// |--------|------------|---------|
/// | 200    | 2.2e-4     | 27 mus   |
/// | 1 000  | 4.4e-5     | 130 mus  |
/// | 5 000  | 8.7e-6     | 650 mus  |
/// | 20 000 | 2.2e-6     | 2.6 ms  |
///
/// 1 000 is chosen because this runs inside the O(N^2) assembly loop: for a chain of
/// on-axis point dipoles *every* pair is collinear, so 20 000 shells would cost ~26 s
/// at N = 100. 4.4e-5 is far below the CDA's own discretisation error (percent level).
///
/// The way to do better is a real Ewald split for 1D (debt item C4), which would make
/// this regime cheap rather than merely tolerable.
pub(crate) const COLLINEAR_SHELLS: i64 = 1_000;

/// Direct Bloch sum of the free-space dyadic over lattice images of a 1D chain:
///
/// G(**r**) = sum_n G_free(**r** - n**a**1) exp(i **k**_b * n**a**1)
///
/// Used for collinear pairs, where the Bessel spectral representation breaks down.
/// The n = 0 image is skipped when **r** ~= 0 (the self-block), leaving the lattice
/// images that the on-site alpha^-1 term does not cover.
pub(crate) fn real_space_1d_direct(
    r: [f64; 3],
    k_bloch: [f64; 3],
    k: f64,
    a1: [f64; 3],
    n_shells: i64,
) -> [[Complex64; 3]; 3] {
    let mut g = [[Complex64::from(0.0); 3]; 3];

    for n in -n_shells..=n_shells {
        let rn = [n as f64 * a1[0], n as f64 * a1[1], n as f64 * a1[2]];
        let s = [r[0] - rn[0], r[1] - rn[1], r[2] - rn[2]];
        let s_sq = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];
        if s_sq < 1e-20 {
            // Coincides with an image -- the divergent self-term, handled by alpha^-1.
            continue;
        }

        let arg = k_bloch[0] * rn[0] + k_bloch[1] * rn[1] + k_bloch[2] * rn[2];
        let phase = Complex64::new(arg.cos(), arg.sin());
        let g_s = dyadic_greens_tensor_at(&s, k);
        for i in 0..3 {
            for j in 0..3 {
                g[i][j] += phase * g_s[i][j];
            }
        }
    }

    g
}

impl EwaldGreens {
    pub fn new(lattice: LatticeSpec) -> Self {
        Self { lattice }
    }

    /// Full periodic dyadic Green's function at Bloch wavevector `k_bloch`.
    ///
    /// # Arguments
    /// * `r`       -- Displacement **r**_i - **r**_j (nm).
    /// * `k_bloch` -- In-plane Bloch wavevector **k**|| (nm^-1). Use `[0;3]` for Gamma.
    /// * `k`       -- Free-space wavenumber 2*pi*n/lambda (nm^-1).
    ///
    /// # Self-image term (r ~= 0)
    /// The R = 0 image is skipped when **r** ~= 0 (i.e. this is the self-block),
    /// leaving only the non-zero lattice images. The caller is responsible for
    /// adding the on-site alpha^-1 term to the diagonal.
    pub fn evaluate(
        &self,
        r: [f64; 3],
        k_bloch: [f64; 3],
        k: f64,
    ) -> [[Complex64; 3]; 3] {
        let gr = self.real_space(r, k_bloch, k);
        let gg = self.recip_space(r, k_bloch, k);
        let mut out = [[Complex64::from(0.0); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                out[i][j] = gr[i][j] + gg[i][j];
            }
        }
        out
    }

    /// Real-space lattice sum.
    ///
    /// **2D lattices**: sums `erfc(eta*|s|) * G_free(s) * exp(i k|| * R)` over
    /// all **R** in Lambda within `truncation_real` shells (eta = `self.lattice.eta()`).
    /// The erfc factor provides exponential convergence.
    ///
    /// **1D chains**: no erfc damping -- uses a direct lattice sum over at least
    /// 20 shells.  This gives ~2 % convergence at typical optical wavelengths.
    /// The full 1D Ewald (Bessel spectral sum) is provided in `recip_space`.
    ///
    /// The **R** = 0 term is excluded when `|r| < 1e-12` nm.
    pub fn real_space(
        &self,
        r: [f64; 3],
        k_bloch: [f64; 3],
        k: f64,
    ) -> [[Complex64; 3]; 3] {
        let r_norm_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        let r_is_zero = r_norm_sq < 1e-24;

        let is_2d = self.lattice.is_2d();
        let eta = self.lattice.eta();

        // 1D chains split by transverse separation rho, and exactly one of the two
        // routines contributes for any given r -- so `evaluate` can still just add them.
        //
        //   rho > 0 : `recip_space_1d` (the Bessel spectral sum) is the complete
        //           periodic G, and this function contributes nothing.
        //   rho ~= 0 : the spectral representation breaks down -- K0(kappa_m rho) diverges
        //           logarithmically for *every* mode -- so `recip_space_1d` returns
        //           zero and the collinear case is summed directly here instead.
        //
        // History: v0.4.1 added a 20-shell direct sum to fix an erfc-zeroing bug.
        // v0.5.0 added recip_space_1d and removed the direct sum, because adding both
        // gave evaluate() ~= 2*G_per. That left collinear pairs returning *exactly
        // zero* -- an on-axis chain felt no lattice at all, and a discretised chain
        // silently dropped its self-image sum. Restored here for rho ~= 0 only.
        if !is_2d {
            let a1 = self.lattice.a1;
            let a1_len = (a1[0] * a1[0] + a1[1] * a1[1] + a1[2] * a1[2]).sqrt();
            let e_hat = [a1[0] / a1_len, a1[1] / a1_len, a1[2] / a1_len];
            let x_par = r[0] * e_hat[0] + r[1] * e_hat[1] + r[2] * e_hat[2];
            let r_perp = [
                r[0] - x_par * e_hat[0],
                r[1] - x_par * e_hat[1],
                r[2] - x_par * e_hat[2],
            ];
            let rho = (r_perp[0] * r_perp[0] + r_perp[1] * r_perp[1] + r_perp[2] * r_perp[2]).sqrt();

            if rho >= COLLINEAR_RHO_TOL {
                return [[Complex64::from(0.0); 3]; 3];
            }
            return real_space_1d_direct(r, k_bloch, k, a1, COLLINEAR_SHELLS);
        }

        let n_shells = self.lattice.truncation_real;
        let lattice_points = self.lattice.real_lattice_points(n_shells);
        let zero = Complex64::from(0.0);
        let mut g = [[zero; 3]; 3];

        for lp in &lattice_points {
            let lp_is_zero = lp[0].abs() < 1e-12 && lp[1].abs() < 1e-12 && lp[2].abs() < 1e-12;
            if r_is_zero && lp_is_zero {
                // Skip the self-image (R=0) to avoid divergence.
                continue;
            }

            // Displacement s = r - R
            let s = [r[0] - lp[0], r[1] - lp[1], r[2] - lp[2]];
            let s_sq = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
            if s_sq < 1e-24 { continue; } // r coincides with this image -- skip

            let s_norm = s_sq.sqrt();

            // Bloch phase: exp(i k|| * R)
            let bloch_arg = k_bloch[0]*lp[0] + k_bloch[1]*lp[1] + k_bloch[2]*lp[2];
            let bloch_phase = Complex64::new(bloch_arg.cos(), bloch_arg.sin());

            // erfc damping factor (2D only -- 1D uses direct sum without damping)
            let damp = if is_2d { erfc_real(eta * s_norm) } else { 1.0 };

            let g_s = dyadic_greens_tensor_at(&s, k);

            for i in 0..3 {
                for j in 0..3 {
                    g[i][j] += bloch_phase * damp * g_s[i][j];
                }
            }
        }

        g
    }

    /// Reciprocal-space (spectral) sum.
    ///
    /// For 2D lattices: evaluates the spectral Green's function sum
    ///
    /// $$\mathbf{G}^\text{spec}_{ab} = \frac{i}{2A} \sum_\mathbf{G}
    ///     \frac{e^{i\mathbf{Q}\cdot\boldsymbol{\rho}}\,e^{-q|z|}}{q}\,
    ///     M_{ab}(\mathbf{Q},q)$$
    ///
    /// For 1D lattices: evaluates the Bessel-function spectral sum
    /// `recip_space_1d`.
    pub fn recip_space(
        &self,
        r: [f64; 3],
        k_bloch: [f64; 3],
        k: f64,
    ) -> [[Complex64; 3]; 3] {
        let zero = Complex64::from(0.0);

        match &self.lattice {
            LatticeSpec { a2: None, a1, .. } => recip_space_1d(r, k_bloch, k, *a1),
            _ => {
                // 2D lattice spectral sum
                let area = self.lattice.cell_measure(); // nm^2
                let recip_points = self.lattice.recip_lattice_points(self.lattice.truncation_recip);

                let i_unit = Complex64::new(0.0, 1.0);
                // Scalar prefactor: i / (2 * A)
                let scalar_pre = i_unit / (2.0 * area);

                let rho_x = r[0];
                let rho_y = r[1];
                let z = r[2];
                let abs_z = z.abs();
                // sign_z: +1 if z >= 0, -1 if z < 0
                let sign_z = if z >= 0.0 { 1.0_f64 } else { -1.0_f64 };

                let k2 = Complex64::from(k * k);
                let mut g = [[zero; 3]; 3];

                for gpt in &recip_points {
                    // Q = k_bloch_|| + G  (in-plane x, y only)
                    let qx = k_bloch[0] + gpt[0];
                    let qy = k_bloch[1] + gpt[1];
                    let q2_real = qx * qx + qy * qy; // Q^2 (real, >= 0)

                    // q = sqrt(Q^2 - k^2) with Im(q) >= 0.
                    // For visible light k is real, Q^2 - k^2 can be negative (propagating)
                    // or positive (evanescent).
                    let q2_complex = Complex64::from(q2_real) - k2;
                    let q = csqrt_positive_imag(q2_complex);

                    // Skip terms where |q| is too small (at the light cone -- near singular).
                    if q.norm() < 1e-10 {
                        continue;
                    }

                    // Bloch phase: exp(i * Q * rho||)
                    let phase_arg = qx * rho_x + qy * rho_y;
                    let bloch_phase = Complex64::new(phase_arg.cos(), phase_arg.sin());

                    // Evanescent/propagating factor: exp(-q * |z|)
                    let exp_qz = (-q * abs_z).exp();

                    // Per-G scalar: (i / 2A) * exp(iQ*rho) * exp(-q|z|) / q
                    let coeff = scalar_pre * bloch_phase * exp_qz / q;

                    // -- Build M tensor: M_ab = k^2delta_ab - k_tilde_a k_tilde_b ------------------
                    //
                    // k_tilde = (Qx, Qy, i*q*sign_z)  [for z > 0; conjugate sign for z < 0]
                    //
                    // k_tilde_x = Qx  (real)
                    // k_tilde_y = Qy  (real)
                    // k_tilde_z = i * q * sign_z
                    let kt_x = Complex64::from(qx);
                    let kt_y = Complex64::from(qy);
                    let kt_z = i_unit * q * sign_z;

                    // k_tilde_a * k_tilde_b outer product
                    let kk = [
                        [kt_x * kt_x, kt_x * kt_y, kt_x * kt_z],
                        [kt_y * kt_x, kt_y * kt_y, kt_y * kt_z],
                        [kt_z * kt_x, kt_z * kt_y, kt_z * kt_z],
                    ];

                    // M_ab = k^2 delta_ab - k_tilde_a k_tilde_b
                    for a in 0..3 {
                        for b in 0..3 {
                            let delta = if a == b { k2 } else { zero };
                            let m_ab = delta - kk[a][b];
                            g[a][b] += coeff * m_ab;
                        }
                    }
                }

                g
            }
        }
    }
}

// -- Private helpers ------------------------------------------------------------

/// Evaluate the free-space dyadic Green's tensor at displacement vector `s`.
///
/// Equivalent to `dyadic_greens_tensor(&s, &[0,0,0], k)` but avoids the
/// subtract-zero computation.
#[inline]
fn dyadic_greens_tensor_at(s: &[f64; 3], k: f64) -> [[Complex64; 3]; 3] {
    let r_sq = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
    let r = r_sq.sqrt();
    assert!(r > 1e-15, "dyadic_greens_tensor_at: zero displacement");

    let kr = k * r;
    let ikr = Complex64::new(0.0, kr);
    let exp_ikr = ikr.exp();
    let prefactor = k * k * exp_ikr / (4.0 * std::f64::consts::PI * r);

    let a = Complex64::from(1.0) + (ikr - Complex64::from(1.0)) / Complex64::from(kr * kr);
    let b = (Complex64::from(3.0) - Complex64::from(3.0) * ikr - Complex64::from(kr * kr))
        / Complex64::from(kr * kr);

    let r_hat = [s[0] / r, s[1] / r, s[2] / r];
    let zero = Complex64::from(0.0);
    let mut g = [[zero; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1.0 } else { 0.0 };
            g[i][j] = prefactor * (a * delta + b * r_hat[i] * r_hat[j]);
        }
    }
    g
}

/// Complementary error function for real non-negative argument.
///
/// Uses the Abramowitz & Stegun formula 7.1.26 rational approximation.
/// Maximum error < 1.5 x 10^-7.
///
/// For x < 0 the identity erfc(x) = 2 - erfc(-x) is applied.
#[inline]
fn erfc_real(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_real(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t * (0.254_829_592
        + t * (-0.284_496_736
        + t * (1.421_413_741
        + t * (-1.453_152_027
        + t * 1.061_405_429))));
    poly * (-x * x).exp()
}

/// Complex square root with branch cut chosen so that Im(result) >= 0.
///
/// For purely real positive z this returns +sqrt(z). For purely real negative z
/// this returns +i*sqrt|z| (evanescent wave convention).
#[inline]
fn csqrt_positive_imag(z: Complex64) -> Complex64 {
    let s = z.sqrt(); // num_complex uses the principal branch: Re(sqrt) >= 0
    // The principal branch has Re >= 0. We want Im >= 0, so if Im < 0 negate.
    if s.im < 0.0 { -s } else { s }
}

// -- Bessel functions (A&S Chapter 9 / Numerical Recipes) ----------------------

/// Euler-Mascheroni constant gamma.
const GAMMA_EM: f64 = 0.577_215_664_901_532_9;

/// Bessel function J0(x).
///
/// Uses rational polynomial approximation (Numerical Recipes sect.6.5) for |x| < 8
/// and asymptotic expansion for |x| >= 8. Accuracy ~1e-7.
#[inline]
pub(crate) fn bessel_j0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 8.0 {
        let y = x * x;
        let p1 = 57568490574.0_f64;
        let p2 = -13362590354.0_f64;
        let p3 = 651619640.7_f64;
        let p4 = -11214424.18_f64;
        let p5 = 77392.33017_f64;
        let p6 = -184.9052456_f64;
        let q1 = 57568490411.0_f64;
        let q2 = 1029532985.0_f64;
        let q3 = 9494680.718_f64;
        let q4 = 59272.64853_f64;
        let q5 = 267.8532712_f64;
        (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * p6)))))
            / (q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + 1.0)))))
    } else {
        let z = 8.0 / ax;
        let y = z * z;
        let xx = ax - 0.785_398_164_f64;
        let p = 1.0
            + y * (-1.098_628_627e-3
            + y * (2.734_510_407e-5
            + y * (-2.073_370_639e-6
            + y * 2.093_887_211e-7)));
        let q = -1.562_499_995e-2
            + y * (1.430_488_765e-4
            + y * (-6.911_147_651e-6
            + y * (7.621_095_161e-7
            + y * (-9.349_451_52e-8))));
        (2.0_f64 / (std::f64::consts::PI * ax)).sqrt() * (xx.cos() * p - z * xx.sin() * q)
    }
}

/// Bessel function J1(x).
///
/// Uses rational polynomial approximation for |x| < 8 and asymptotic for |x| >= 8.
/// J1(-x) = -J1(x). Accuracy ~1e-7.
#[inline]
pub(crate) fn bessel_j1(x: f64) -> f64 {
    let ax = x.abs();
    let sign = if x < 0.0 { -1.0_f64 } else { 1.0_f64 };
    if ax < 8.0 {
        let y = x * x;
        let p1 = 72362614232.0_f64;
        let p2 = -7895059235.0_f64;
        let p3 = 242396853.1_f64;
        let p4 = -2972611.439_f64;
        let p5 = 15704.48260_f64;
        let p6 = -30.16036606_f64;
        let q1 = 144725228442.0_f64;
        let q2 = 2300535178.0_f64;
        let q3 = 18583304.74_f64;
        let q4 = 99447.43394_f64;
        let q5 = 376.9991397_f64;
        x * (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * p6)))))
            / (q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + 1.0)))))
    } else {
        let z = 8.0 / ax;
        let y = z * z;
        let xx = ax - 2.356_194_491_f64;
        let p = 1.0
            + y * (1.831_05e-3
            + y * (-3.516_396_496e-5
            + y * (2.457_520_174e-6
            + y * (-2.403_370_19e-7))));
        let q = 4.687_499_995e-2
            + y * (-2.002_690_873e-4
            + y * (8.449_199_096e-6
            + y * (-8.822_898_7e-7
            + y * 1.057_874_12e-7)));
        sign * (2.0_f64 / (std::f64::consts::PI * ax)).sqrt() * (xx.cos() * p - z * xx.sin() * q)
    }
}

/// Bessel function Y0(x) -- real, x > 0.
///
/// For x <= 3: uses the series Y0(x) = (2/pi)[(ln(x/2)+gamma)J0(x) - sumH_k(-x^2/4)^k/(k!)^2]
/// For x > 3: uses the asymptotic form f0/sqrt(x) * sin(theta0).
/// Accuracy ~1e-7.
#[inline]
pub(crate) fn bessel_y0(x: f64) -> f64 {
    debug_assert!(x > 0.0, "bessel_y0: x must be positive");
    use std::f64::consts::PI;
    if x <= 3.0 {
        // Series: Y0(x) = (2/pi)*[(ln(x/2)+gamma)*J0(x) - sum_{k=1}^inf H_k*(-x^2/4)^k/(k!)^2]
        let h = x / 2.0;
        let mut s = 0.0_f64;
        let mut term = -(h * h); // k=1: (-x^2/4)^1/(1!)^2
        let mut hk = 1.0_f64;   // H_1 = 1
        for k in 1_i32..=30 {
            s += term * hk;
            term *= -(h * h) / ((k + 1) as f64).powi(2);
            hk += 1.0 / (k + 1) as f64;
        }
        (2.0 / PI) * ((h.ln() + GAMMA_EM) * bessel_j0(x) - s)
    } else {
        let z = 3.0 / x;
        let y = z * z;
        let f0 = 0.797_884_56_f64
            + y * (-7.7e-7
            + y * (-5.527_40e-3
            + y * (-9.512e-5
            + y * (1.372_37e-3
            + y * (-7.280_5e-4
            + y * 1.447_6e-4)))));
        let theta0 = x - 0.785_398_164_f64
            + z * (-4.166_397e-2
            + z * (-3.954e-5
            + z * (2.625_73e-3
            + z * (-5.412_5e-4
            + z * (-2.933_3e-4
            + z * 1.355_8e-4)))));
        f0 / x.sqrt() * theta0.sin()
    }
}

/// Bessel function Y1(x) -- real, x > 0.
///
/// For x <= 3: uses series Y1 = -Y0'. Computed from the Y0 power series via
/// the identity Y1(x) = -Y0'(x) using a 5-point central difference.
/// For x > 3: uses the asymptotic form f1/sqrt(x) * sin(theta1).
/// Accuracy ~1e-7.
#[inline]
pub(crate) fn bessel_y1(x: f64) -> f64 {
    debug_assert!(x > 0.0, "bessel_y1: x must be positive");
    if x <= 3.0 {
        // Y1 = -dY0/dx, computed via 5-point stencil for accuracy.
        // h chosen to balance truncation vs rounding: h ~= x^(1/5) * 1e-3
        // h must be small enough that all stencil points (x +/- 2h) remain positive
        let h = (x * 1e-3_f64).min(x / 4.0).max(1e-15_f64);
        let y0_pp = bessel_y0(x + 2.0 * h);
        let y0_p  = bessel_y0(x + h);
        let y0_m  = bessel_y0(x - h);
        let y0_mm = bessel_y0(x - 2.0 * h);
        -((-y0_pp + 8.0 * y0_p - 8.0 * y0_m + y0_mm) / (12.0 * h))
    } else {
        let z = 3.0 / x;
        let y = z * z;
        let f1 = 0.797_884_56_f64
            + y * (1.564e-4
            + y * (1.659_667e-2
            + y * (1.710_5e-4
            + y * (-2.495_11e-3
            + y * (1.136_53e-3
            + y * (-2.003_3e-4))))));
        let theta1 = x - 2.356_194_49_f64
            + z * (1.249_961_2e-1
            + z * (5.650e-5
            + z * (-6.378_79e-3
            + z * (7.434_8e-4
            + z * (7.982_4e-4
            + z * (-2.916_6e-4))))));
        f1 / x.sqrt() * theta1.sin()
    }
}

/// Modified Bessel function I0(x) -- A&S 9.8.1/9.8.3 -- real, any x.
#[inline]
fn bessel_i0(x: f64) -> f64 {
    let ax = x.abs();
    if ax <= 3.75 {
        let y = (x / 3.75) * (x / 3.75);
        1.0 + y * (3.515_622_9
            + y * (3.089_942_4
            + y * (1.206_749_2
            + y * (0.265_973_2
            + y * (0.036_007_68
            + y * 0.004_581_9)))))
    } else {
        let y = 3.75 / ax;
        let poly = 0.398_942_28
            + y * (0.013_285_92
            + y * (0.002_253_19
            + y * (-0.001_575_65
            + y * (0.009_166_9
            + y * (-0.020_570_64
            + y * (0.026_355_37
            + y * (-0.016_470_8
            + y * 0.003_920_78)))))));
        poly * ax.exp() / ax.sqrt()
    }
}

/// Modified Bessel function I1(x) -- A&S 9.8.2/9.8.4 -- real, any x.
#[inline]
fn bessel_i1(x: f64) -> f64 {
    let ax = x.abs();
    let sign = if x < 0.0 { -1.0_f64 } else { 1.0_f64 };
    let val = if ax <= 3.75 {
        let y = (x / 3.75) * (x / 3.75);
        ax * (0.5
            + y * (0.878_906_94
            + y * (0.514_988_69
            + y * (0.150_849_34
            + y * (0.026_587_33
            + y * (0.003_015_32
            + y * 0.000_324_11))))))
    } else {
        let y = 3.75 / ax;
        let poly = 0.398_942_28
            + y * (-0.039_880_24
            + y * (-0.003_620_18
            + y * (0.016_380_1
            + y * (-0.103_155_5
            + y * (0.228_296_7
            + y * (-0.289_531_2
            + y * (0.178_765_4
            + y * (-0.042_005_9))))))));
        poly * ax.exp() / ax.sqrt()
    };
    sign * val
}

/// Modified Bessel function K0(x) -- A&S 9.8.5/9.8.7 -- real, x > 0.
///
/// K0(x) = -ln(x/2)*I0(x) + polynomial      for 0 < x <= 2
/// K0(x) = e^{-x}/sqrt(x) * polynomial            for x > 2
#[inline]
pub(crate) fn bessel_k0(x: f64) -> f64 {
    debug_assert!(x > 0.0, "bessel_k0: x must be positive");
    if x <= 2.0 {
        let y = (x / 2.0) * (x / 2.0);
        let poly = -0.577_215_66
            + y * (0.422_784_20
            + y * (0.230_697_56
            + y * (0.034_885_90
            + y * (0.002_626_98
            + y * (0.000_107_50
            + y * 0.000_000_74)))));
        -bessel_i0(x) * (x / 2.0).ln() + poly
    } else {
        let y = 2.0 / x;
        // A&S 9.8.7 -- 6-term polynomial, ~7-digit accuracy for x > 2
        let poly = 1.253_314_14
            + y * (-0.078_323_58
            + y * (0.021_895_68
            + y * (-0.010_624_46
            + y * (0.005_878_72
            + y * (-0.002_515_40
            + y * 0.000_532_08)))));
        poly * (-x).exp() / x.sqrt()
    }
}

/// Modified Bessel function K1(x) -- A&S 9.8.6/9.8.8 -- real, x > 0.
///
/// K1(x) = ln(x/2)*I1(x) + polynomial/x     for 0 < x <= 2
/// K1(x) = e^{-x}/sqrt(x) * polynomial            for x > 2
#[inline]
pub(crate) fn bessel_k1(x: f64) -> f64 {
    debug_assert!(x > 0.0, "bessel_k1: x must be positive");
    if x <= 2.0 {
        let y = (x / 2.0) * (x / 2.0);
        // A&S 9.8.6: x*K1(x) - ln(x/2)*x*I1(x) = polynomial(y)
        let poly = 1.0
            + y * (0.154_431_44
            + y * (-0.672_785_79
            + y * (-0.181_568_97
            + y * (-0.019_194_02
            + y * (-0.001_104_04
            + y * (-0.000_046_86))))));
        bessel_i1(x) * (x / 2.0).ln() + poly / x
    } else {
        let y = 2.0 / x;
        // A&S 9.8.8 -- 6-term polynomial
        let poly = 1.253_314_14
            + y * (0.234_986_19
            + y * (-0.036_556_20
            + y * (0.015_042_68
            + y * (-0.007_803_53
            + y * (0.003_256_14
            + y * (-0.000_682_45))))));
        poly * (-x).exp() / x.sqrt()
    }
}

/// Complex modified Bessel function K0(z) using A&S 9.6.4 identities.
///
/// Handles two cases:
/// - Purely real z (z.im ~= 0): calls `bessel_k0(z.re)`
/// - Purely imaginary z = ialpha (z.re ~= 0, alpha > 0): uses K0(ialpha) = (pi/2)(-Y0(alpha) + i*J0(alpha))
/// - General complex: uses the purely real branch on |z| (safe approximation).
#[inline]
fn bessel_k0_cplx(z: Complex64) -> Complex64 {
    use std::f64::consts::PI;
    if z.re.abs() < 1e-12 {
        // Purely imaginary: z = i*alpha, alpha > 0
        let alpha = z.im;
        if alpha > 0.0 {
            // K0(ialpha) = (pi/2)(-Y0(alpha) + i*J0(alpha))
            Complex64::new(-(PI / 2.0) * bessel_y0(alpha), (PI / 2.0) * bessel_j0(alpha))
        } else if alpha < 0.0 {
            // K0(-ialpha) = conj(K0(ialpha)) for real K0 on real axis
            let alpha_pos = -alpha;
            Complex64::new(-(PI / 2.0) * bessel_y0(alpha_pos), -(PI / 2.0) * bessel_j0(alpha_pos))
        } else {
            // z = 0: singular, return large value
            Complex64::new(f64::INFINITY, 0.0)
        }
    } else if z.im.abs() < 1e-12 {
        // Purely real
        if z.re > 0.0 {
            Complex64::new(bessel_k0(z.re), 0.0)
        } else {
            // Negative real argument: K0 is complex (branch cut)
            // K0(-x) = K0(x) - ipi*I0(x) for x > 0
            let xp = -z.re;
            Complex64::new(bessel_k0(xp), -std::f64::consts::PI * bessel_i0(xp))
        }
    } else {
        // General complex -- use asymptotic or series; for the 1D Ewald sum
        // in this codebase kappa is always real or purely imaginary, so this
        // branch is only a safety fallback.
        debug_assert!(false, "bessel_k0_cplx: general complex branch reached -- kappa_m should be real or purely imaginary in this codebase");
        let r = z.norm();
        if r > 1e-12 {
            // Crude: use real part only (acceptable since this branch is not hit
            // in normal operation).
            Complex64::new(bessel_k0(r), 0.0)
        } else {
            Complex64::new(f64::INFINITY, 0.0)
        }
    }
}

/// Complex modified Bessel function K1(z) using A&S 9.6.4 identities.
///
/// - Purely real z: calls `bessel_k1(z.re)`
/// - Purely imaginary z = ialpha (alpha > 0): K1(ialpha) = -(pi/2)(J1(alpha) + i*Y1(alpha))
#[inline]
fn bessel_k1_cplx(z: Complex64) -> Complex64 {
    use std::f64::consts::PI;
    if z.re.abs() < 1e-12 {
        let alpha = z.im;
        if alpha > 0.0 {
            // K1(ialpha) = -(pi/2)(J1(alpha) + i*Y1(alpha))
            Complex64::new(-(PI / 2.0) * bessel_j1(alpha), -(PI / 2.0) * bessel_y1(alpha))
        } else if alpha < 0.0 {
            // K1(-ialpha) = -conj(K1(ialpha))
            let alpha_pos = -alpha;
            Complex64::new((PI / 2.0) * bessel_j1(alpha_pos), -(PI / 2.0) * bessel_y1(alpha_pos))
        } else {
            Complex64::new(f64::INFINITY, 0.0)
        }
    } else if z.im.abs() < 1e-12 {
        if z.re > 0.0 {
            Complex64::new(bessel_k1(z.re), 0.0)
        } else {
            // K1(-x) = -K1(x) + ipi*I1(x) for x > 0
            let xp = -z.re;
            Complex64::new(-bessel_k1(xp), std::f64::consts::PI * bessel_i1(xp))
        }
    } else {
        debug_assert!(false, "bessel_k1_cplx: general complex branch reached -- kappa_m should be real or purely imaginary in this codebase");
        let r = z.norm();
        if r > 1e-12 {
            Complex64::new(bessel_k1(r), 0.0)
        } else {
            Complex64::new(f64::INFINITY, 0.0)
        }
    }
}

/// 1D reciprocal-space spectral sum using modified Bessel functions K0 and K1.
///
/// Evaluates the Bloch-periodic dyadic Green's function spectral sum for a
/// 1D chain with lattice vector `a1`. The sum runs over all reciprocal modes
/// m in [-m_max, m_max], where m_max covers the light cone and resolves the
/// transverse decay exp(-2*pi*m*rho/a1) -- see the mode-count comment in the body.
///
/// Accuracy is ~1e-7 to ~1e-5 for rho/a1 >= 0.016. Below that the mode cap binds and
/// the sum is **not converged**; a `debug_assert` fires. See the body for numbers.
///
/// # Arguments
/// * `r`       -- Displacement vector **r**_i - **r**_j (nm).
/// * `k_bloch` -- Bloch wavevector (nm^-1).
/// * `k`       -- Free-space wavenumber (nm^-1).
/// * `a1`      -- Lattice vector (nm).
///
/// Returns zero tensor for collinear displacement (rho < 1e-10 nm).
fn recip_space_1d(
    r: [f64; 3],
    k_bloch: [f64; 3],
    k: f64,
    a1: [f64; 3],
) -> [[Complex64; 3]; 3] {
    use std::f64::consts::PI;

    let a1_len = (a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]).sqrt();
    let e_hat = [a1[0]/a1_len, a1[1]/a1_len, a1[2]/a1_len]; // unit chain axis

    // Parallel displacement
    let x_par = r[0]*e_hat[0] + r[1]*e_hat[1] + r[2]*e_hat[2];

    // Perpendicular displacement
    let r_perp = [
        r[0] - x_par*e_hat[0],
        r[1] - x_par*e_hat[1],
        r[2] - x_par*e_hat[2],
    ];
    let rho = (r_perp[0]*r_perp[0] + r_perp[1]*r_perp[1] + r_perp[2]*r_perp[2]).sqrt();

    // rho = 0 guard: when i == j and dipoles are collinear (same position projected
    // onto the chain axis), the transverse separation rho = 0 and the spectral sum
    // diverges. For the self-block the polarisability inverse (alpha^-1) dominates;
    // real-space handles the collinear near-field, and we return zero here so the
    // solver is not corrupted. This introduces ~2% error vs. a full self-image
    // treatment but is acceptable for the current use case.
    if rho < 1e-10 {
        return [[Complex64::from(0.0); 3]; 3];
    }

    let rho_hat = [r_perp[0]/rho, r_perp[1]/rho, r_perp[2]/rho];
    let kbloch_par = k_bloch[0]*e_hat[0] + k_bloch[1]*e_hat[1] + k_bloch[2]*e_hat[2];

    // Number of modes. Two requirements:
    //   1. cover the light cone (|g_m| <= k), plus 10 evanescent modes;
    //   2. resolve the transverse decay. Term m falls off as exp(-kappa_m rho) with
    //      kappa_m -> 2*pi*|m|/a1, so the exponential factor is exp(-2*pi*m*rho/a1) and
    //      m >~ 20*a1/(2*pi*rho) drives it to ~2e-9.
    //
    // The *achieved* truncation error is worse than that exponential, because the
    // coefficients carry algebraic growth (kappa^2K0) and the tail sums geometrically:
    // measured ~1e-7 at rho/a1 = 0.25, rising to ~3e-5 at rho/a1 = 0.016. Budget ~1e-5,
    // not 1e-9. Requirement 2 dominates for rho << a1 -- the usual case for a chain of
    // discretised particles -- and a rho-independent count silently loses accuracy
    // there (~1e-3 at rho/a1 = 0.14).
    //
    // MODE_CAP bounds cost as rho -> 0, but once it binds the error is a CLIFF, not a
    // slope. Measured against a 5000-mode reference (a1 = 100 nm, lambda = 600 nm):
    //   rho/a1 = 0.010 (needs 319) -> ~1e-2 relative
    //   rho/a1 = 0.005 (needs 637) -> ~3       relative
    //   rho/a1 = 0.002 (needs 1592) -> ~3e2    relative, with spurious Im parts
    // rho/a1 <~ 0.02 is reachable in practice (a 100 nm period discretised at 1-2 nm
    // puts the nearest off-axis pairs there), so the assertion below exists to make a
    // cap-limited evaluation fail loudly rather than return a plausible number. It is
    // a real `assert!`, not `debug_assert!`: a release build is precisely where a
    // silently wrong answer does damage, and the cost is nil beside the ~400 complex
    // Bessel evaluations this function is about to perform.
    // The principled fix is a real Ewald split for 1D -- an erfc-damped real-space sum
    // plus the spectral sum, as the 2D path already does -- where small rho is precisely
    // the regime the real-space sum handles cheaply. Logged as debt.
    const MODE_CAP: i64 = 200;
    let light_cone_modes = k * a1_len / (2.0 * PI) + 10.0;
    let decay_modes = 20.0 * a1_len / (2.0 * PI * rho);
    assert!(
        decay_modes.ceil() as i64 <= MODE_CAP,
        "1D Ewald spectral sum is mode-cap limited at rho/a1 = {:.4} ({:.0} modes needed, \
         {} available): the result is not converged and may be wrong by orders of \
         magnitude. Increase the lattice period, coarsen the dipole spacing, or \
         implement the 1D Ewald split.",
        rho / a1_len,
        decay_modes.ceil(),
        MODE_CAP
    );
    let m_max = (light_cone_modes.max(decay_modes).ceil() as i64).min(MODE_CAP);

    let zero = Complex64::from(0.0);
    let mut tensor = [[zero; 3]; 3];

    for m in -m_max..=m_max {
        let g_m = kbloch_par + 2.0 * PI * (m as f64) / a1_len;

        // kappa_m = sqrt(g_m^2 - k^2), Im(kappa) >= 0
        let kappa = csqrt_positive_imag(Complex64::new(g_m*g_m - k*k, 0.0));

        let arg = kappa * rho; // Complex64
        let k0 = bessel_k0_cplx(arg);
        let k1 = bessel_k1_cplx(arg);
        let k1_rho = kappa * k1 / rho; // kappa*K1(kapparho)/rho

        // Phase along chain axis: exp(i g_m x_par)
        let phase = Complex64::new(0.0, g_m * x_par).exp();
        let prefactor = phase / (2.0 * PI * a1_len);

        // Tensor structure coefficients
        let c_delta = prefactor * (Complex64::from(k*k) * k0 - k1_rho);
        let c_ee    = prefactor * (k1_rho - Complex64::from(g_m*g_m) * k0);
        // NB: no /rho here, unlike the other three coefficients. The cross term is
        // d_xd_alpha K0 = (i g_m)(-kappaK1(kapparho) rho_hat_alpha), so it carries a bare kappaK1(kapparho); the
        // 1/rho in `k1_rho` belongs only to the purely transverse derivatives.
        let c_cross = prefactor * Complex64::new(0.0, -g_m) * kappa * k1;
        let c_rr    = prefactor * (kappa*kappa * k0 + Complex64::from(2.0)*k1_rho);

        // Accumulate into tensor
        for alpha in 0..3 {
            for beta in 0..3 {
                let delta = if alpha == beta { Complex64::from(1.0) } else { zero };
                tensor[alpha][beta] +=
                    c_delta * delta
                  + c_ee    * (e_hat[alpha] * e_hat[beta])
                  + c_cross * (e_hat[alpha]*rho_hat[beta] + e_hat[beta]*rho_hat[alpha])
                  + c_rr    * (rho_hat[alpha] * rho_hat[beta]);
            }
        }
    }

    tensor
}

#[cfg(test)]
mod tests {
    use super::*;
    use lumina_geometry::lattice::LatticeSpec;

    fn near_c(a: Complex64, b: Complex64, tol: f64) -> bool {
        (a - b).norm() < tol
    }

    /// erfc(0) == 1, erfc large == 0
    #[test]
    fn test_erfc_basic() {
        assert!((erfc_real(0.0) - 1.0).abs() < 1e-7);
        assert!(erfc_real(10.0) < 1e-40);
        assert!((erfc_real(-0.0) - 1.0).abs() < 1e-7);
        // erfc(x) + erfc(-x) == 2
        let x = 1.5_f64;
        assert!((erfc_real(x) + erfc_real(-x) - 2.0).abs() < 1e-6);
    }

    /// csqrt_positive_imag: for negative real argument, Im should be positive.
    #[test]
    fn test_csqrt_positive_imag() {
        let z = Complex64::new(-4.0, 0.0);
        let s = csqrt_positive_imag(z);
        assert!(s.im > 0.0, "Im should be >= 0 for negative real argument");
        assert!((s.im - 2.0).abs() < 1e-12);
        // For positive real: result should be real positive
        let z2 = Complex64::new(9.0, 0.0);
        let s2 = csqrt_positive_imag(z2);
        assert!((s2.re - 3.0).abs() < 1e-12);
        assert!(s2.im.abs() < 1e-12);
    }

    /// For a square 2D array, the periodic Green's function evaluated at r=0
    /// (self-image sum, R!=0 only) should be finite and non-zero.
    #[test]
    fn test_self_image_finite() {
        let lat = LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 600.0;
        let g = eg.evaluate([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], k);
        for i in 0..3 {
            for j in 0..3 {
                assert!(g[i][j].re.is_finite(), "G[{i}][{j}] real not finite");
                assert!(g[i][j].im.is_finite(), "G[{i}][{j}] imag not finite");
            }
        }
    }

    /// evaluate() == real_space() + recip_space()
    #[test]
    fn test_evaluate_splits_correctly() {
        let lat = LatticeSpec::planar([80.0, 0.0, 0.0], [0.0, 80.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 500.0;
        let r = [10.0, 5.0, 2.0];
        let kb = [0.0, 0.0, 0.0];
        let g_full = eg.evaluate(r, kb, k);
        let gr = eg.real_space(r, kb, k);
        let gg = eg.recip_space(r, kb, k);
        for i in 0..3 {
            for j in 0..3 {
                assert!(near_c(g_full[i][j], gr[i][j] + gg[i][j], 1e-20));
            }
        }
    }

    /// Reciprocity at the Gamma point: G_per(r)[i,j] == G_per(-r)[j,i].
    ///
    /// For k|| = 0 the periodic Green's function satisfies the same reciprocity
    /// as the free-space G because the Bloch phases cancel with their conjugates.
    #[test]
    fn test_gamma_point_reciprocity() {
        let lat = LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 600.0;
        let r = [15.0, 8.0, 3.0];
        let neg_r = [-r[0], -r[1], -r[2]];
        let kb = [0.0, 0.0, 0.0];
        let g_r   = eg.evaluate(r, kb, k);
        let g_neg = eg.evaluate(neg_r, kb, k);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    near_c(g_r[i][j], g_neg[j][i], 1e-8),
                    "G[{i}][{j}] reciprocity: {:?} vs {:?}",
                    g_r[i][j], g_neg[j][i]
                );
            }
        }
    }

    /// The Gamma-point self-image sum should be purely imaginary diagonal on a
    /// square lattice (by 4-fold symmetry the real-part contributions cancel).
    #[test]
    fn test_gamma_self_image_symmetry() {
        let lat = LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 600.0;
        let g = eg.evaluate([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], k);
        // Gxx and Gyy should be equal (square symmetry)
        assert!(
            near_c(g[0][0], g[1][1], 1e-8),
            "Gxx={:?} Gyy={:?} should be equal on square lattice",
            g[0][0], g[1][1]
        );
        // Off-diagonal Gxy should vanish (square symmetry)
        assert!(
            g[0][1].norm() < 1e-8,
            "Gxy={:?} should vanish by 4-fold symmetry",
            g[0][1]
        );
    }

    /// 1D chain: real-space sum should be finite, non-zero, and reciprocal.
    ///
    /// Prior to the direct-sum fix (v0.4.1) the erfc eta was so large that every
    /// term was zeroed -> evaluate() returned ~= 0 silently.
    #[test]
    fn test_1d_chain_finite() {
        let lat = LatticeSpec::chain([80.0, 0.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 500.0;
        let r = [10.0, 5.0, 0.0];
        let g = eg.evaluate(r, [0.0; 3], k);
        let mut max_mag = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                assert!(g[i][j].re.is_finite());
                assert!(g[i][j].im.is_finite());
                max_mag = max_mag.max(g[i][j].norm());
            }
        }
        // The result must be genuinely non-zero (dominant element >> 0).
        assert!(
            max_mag > 1e-12,
            "1D chain G ~= 0 -- erfc eta bug: max_mag = {:.2e}",
            max_mag
        );
    }

    /// 1D chain: evaluate() at Gamma point satisfies G(r)[i,j] == G(-r)[j,i] (reciprocity).
    #[test]
    fn test_1d_chain_reciprocity() {
        let lat = LatticeSpec::chain([80.0, 0.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 600.0;
        let r = [10.0, 5.0, 0.0];
        let neg_r = [-r[0], -r[1], -r[2]];
        let kb = [0.0; 3];
        let g_r = eg.evaluate(r, kb, k);
        let g_neg = eg.evaluate(neg_r, kb, k);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    near_c(g_r[i][j], g_neg[j][i], 1e-6),
                    "1D reciprocity G[{i}][{j}] = {:?}, G[-r][{j}][{i}] = {:?}",
                    g_r[i][j], g_neg[j][i]
                );
            }
        }
    }

    /// Convergence: adding more shells should not drastically change the result.
    #[test]
    fn test_convergence_with_shells() {
        let lat_coarse = LatticeSpec { truncation_real: 5,  ..LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]) };
        let lat_fine   = LatticeSpec { truncation_real: 10, ..LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]) };
        let eg5  = EwaldGreens::new(lat_coarse);
        let eg10 = EwaldGreens::new(lat_fine);
        let k = 2.0 * std::f64::consts::PI / 600.0;
        let r = [15.0, 8.0, 3.0];
        let kb = [0.0; 3];
        let g5  = eg5.evaluate(r, kb, k);
        let g10 = eg10.evaluate(r, kb, k);
        // The change from 5->10 shells should be < 5% of the dominant element.
        // (erfc damping makes the sum converge rapidly, so even few shells is enough)
        let max5: f64 = g5.iter().flatten().map(|v| v.norm()).fold(0.0f64, f64::max);
        let max_diff: f64 = (0..3).flat_map(|i| (0..3).map(move |j| (g5[i][j]-g10[i][j]).norm()))
            .fold(0.0f64, f64::max);
        assert!(max_diff < 0.05 * max5, "Convergence poor: Delta={:.2e} vs max={:.2e}", max_diff, max5);
    }

    /// recip_space for 2D lattice returns finite values.
    #[test]
    fn test_2d_recip_space_finite() {
        let lat = LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 600.0;
        let r = [10.0, 5.0, 20.0]; // z != 0 -- well-behaved case
        let g = eg.recip_space(r, [0.0; 3], k);
        for i in 0..3 {
            for j in 0..3 {
                assert!(g[i][j].re.is_finite(), "recip G[{i}][{j}] re not finite");
                assert!(g[i][j].im.is_finite(), "recip G[{i}][{j}] im not finite");
            }
        }
    }

    // -- New Bessel function tests ----------------------------------------------

    #[test]
    fn test_bessel_k0_values() {
        // A&S Table 9.8 reference values (10-digit precision)
        // x=1.0 is well inside the series regime (x <= 2): ~1e-8 accuracy
        assert!((bessel_k0(1.0) - 0.4210244382_f64).abs() < 1e-6);
        // x=2.0 is on the polynomial branch boundary: ~1e-5 accuracy
        assert!((bessel_k0(2.0) - 0.1138938727_f64).abs() < 1e-4);
        // Large x: K0(5) ~= 0.003691
        assert!((bessel_k0(5.0) - 0.003691_f64).abs() < 1e-6);
    }

    #[test]
    fn test_bessel_k1_values() {
        // A&S Table 9.8 reference values (10-digit precision)
        // x=1.0 is well inside the series regime (x <= 2): ~1e-7 accuracy
        assert!((bessel_k1(1.0) - 0.6019072302_f64).abs() < 1e-6);
        // x=2.0 is on the polynomial branch boundary: ~1e-4 accuracy
        assert!((bessel_k1(2.0) - 0.1399358445_f64).abs() < 1e-4);
    }

    #[test]
    fn test_bessel_y0_values() {
        // Y0(1) ~= 0.08826 (verified: DLMF 10.74, ODE integration, power series)
        assert!((bessel_y0(1.0) - 0.08826_f64).abs() < 1e-4);
        // Y0(2) ~= 0.51038 (verified: power series, ODE integration, Wronskian)
        // Note: 0.10703 would be -Y1(2); Y0(2) is the peak of Y0 near x = 2.2.
        assert!((bessel_y0(2.0) - 0.51038_f64).abs() < 1e-4);
    }

    #[test]
    fn test_bessel_y1_values() {
        // A&S Table 9.1: Y1(1) ~= -0.7812, Y1(2) ~= -0.1070
        assert!((bessel_y1(1.0) - (-0.7812_f64)).abs() < 1e-4);
        assert!((bessel_y1(2.0) - (-0.1070_f64)).abs() < 1e-4);
        // Small-x coverage: verifies the 5-point stencil doesn't produce NaN or
        // diverge when x is small (the stencil guard ensures x - 2h > 0 always).
        // Reference: Wolfram Alpha BesselY[1, x].
        // Tolerance 1e-2 for x=0.1 (near the log singularity, stencil accuracy degrades).
        assert!((bessel_y1(0.1) - (-6.4590_f64)).abs() < 1e-2,
            "Y1(0.1) = {:.4}, expected ~= -6.4590", bessel_y1(0.1));
        assert!((bessel_y1(0.5) - (-1.4715_f64)).abs() < 1e-3,
            "Y1(0.5) = {:.4}, expected ~= -1.4715", bessel_y1(0.5));
    }

    #[test]
    fn test_k0_complex_imaginary() {
        // K0(ialpha) = (pi/2)(-Y0(alpha) + i*J0(alpha)) -- POSITIVE i*J0 sign
        use std::f64::consts::PI;
        for alpha in [1.0_f64, 2.0, 3.0] {
            let z = Complex64::new(0.0, alpha);
            let result = bessel_k0_cplx(z);
            let expected_re = -(PI/2.0) * bessel_y0(alpha);
            let expected_im =  (PI/2.0) * bessel_j0(alpha);
            assert!((result.re - expected_re).abs() < 1e-5,
                "K0(i*{}) re: got {}, expected {}", alpha, result.re, expected_re);
            assert!((result.im - expected_im).abs() < 1e-5,
                "K0(i*{}) im: got {}, expected {}", alpha, result.im, expected_im);
        }
    }

    #[test]
    fn test_k1_complex_imaginary() {
        // K1(ialpha) = -(pi/2)(J1(alpha) + i*Y1(alpha))
        use std::f64::consts::PI;
        for alpha in [1.0_f64, 2.0, 3.0] {
            let z = Complex64::new(0.0, alpha);
            let result = bessel_k1_cplx(z);
            let expected_re = -(PI/2.0) * bessel_j1(alpha);
            let expected_im = -(PI/2.0) * bessel_y1(alpha);
            assert!((result.re - expected_re).abs() < 1e-5,
                "K1(i*{}) re: got {}, expected {}", alpha, result.re, expected_re);
            assert!((result.im - expected_im).abs() < 1e-5,
                "K1(i*{}) im: got {}, expected {}", alpha, result.im, expected_im);
        }
    }

    #[test]
    fn test_1d_recip_reciprocity() {
        // G_recip(Deltar)^T = G_recip(-Deltar) (since free-space G is even and the sum is symmetric)
        let k = 0.01;
        let k_bloch = [0.0, 0.0, 0.0];
        let r = [3.0, 5.0, 0.0];
        let g1 = recip_space_1d(r, k_bloch, k, [10.0, 0.0, 0.0]);
        let neg_r = [-r[0], -r[1], -r[2]];
        let g2 = recip_space_1d(neg_r, k_bloch, k, [10.0, 0.0, 0.0]);
        // G(Deltar)^T = G(-Deltar) since G is even
        for i in 0..3 {
            for j in 0..3 {
                let diff = (g1[i][j] - g2[j][i]).norm();
                assert!(diff < 1e-10, "Reciprocity failed at ({},{}) diff={}", i, j, diff);
            }
        }
    }

    #[test]
    fn test_1d_recip_finite_nonzero() {
        // For non-collinear displacement (rho > 0), recip_space_1d must return a non-zero tensor
        let k = 0.01;
        let k_bloch = [0.0, 0.0, 0.0];
        let r = [3.0, 5.0, 0.0];   // rho = 5.0 nm > 0
        let g = recip_space_1d(r, k_bloch, k, [10.0, 0.0, 0.0]);
        let norm: f64 = g.iter().flat_map(|row| row.iter()).map(|c| c.norm()).sum();
        assert!(norm > 1e-12, "recip_space_1d should be non-zero for rho>0, got norm={}", norm);
    }

    #[test]
    fn test_1d_recip_rho_zero_guard() {
        // For collinear displacement (rho=0), recip_space_1d returns zero tensor
        let k = 0.01;
        let r = [3.0, 0.0, 0.0];   // rho = 0
        let g = recip_space_1d(r, [0.0, 0.0, 0.0], k, [10.0, 0.0, 0.0]);
        for row in &g { for val in row { assert!(val.norm() < 1e-30); } }
    }

    /// Reference implementation for three independent tensor elements.
    ///
    /// Computes G_xx, G_yy, and G_xy directly from the spectral sum formula
    /// without the c_delta/c_ee/c_cross/c_rr abstraction used in `recip_space_1d`.
    /// Used in tensor-element tests to catch prefactor bugs (e.g. missing /rho in c_cross).
    ///
    /// Chain along x_hat, displacement r = [Deltax, rho, 0] (rho > 0).
    /// - G_xx = (1/2*pi*a1) sum_m exp(iG_m Deltax) (k^2-G_m^2) K0(kappa_m rho)
    /// - G_yy = (1/2*pi*a1) sum_m exp(iG_m Deltax) (G_m^2 K0 + kappa_m K1/rho)
    /// - G_xy = (1/2*pi*a1) sum_m exp(iG_m Deltax) (-i G_m) kappa_m K1(kappa_m rho)
    ///
    /// Note G_xy carries a bare kappaK1, with no 1/rho -- it comes from d_xd_y K0,
    /// whereas the 1/rho in the other elements comes from drho_hat_alpha/dbeta.
    fn recip_1d_reference_elements(
        dx: f64,
        rho: f64,
        k: f64,
        a1: f64,
        m_max: i64,
    ) -> (Complex64, Complex64, Complex64) {
        use std::f64::consts::PI;
        let two_pi_over_a1 = 2.0 * PI / a1;
        let norm = 1.0 / (2.0 * PI * a1);
        let (mut gxx, mut gyy, mut gxy) = (Complex64::from(0.0), Complex64::from(0.0), Complex64::from(0.0));
        for m in -m_max..=m_max {
            let g_m = m as f64 * two_pi_over_a1;
            let kappa = csqrt_positive_imag(Complex64::new(g_m * g_m - k * k, 0.0));
            let k0 = bessel_k0_cplx(kappa * rho);
            let k1 = bessel_k1_cplx(kappa * rho);
            let k1_rho = kappa * k1 / rho;
            let phase = Complex64::new(0.0, g_m * dx).exp();
            gxx += phase * Complex64::from(k * k - g_m * g_m) * k0;
            gyy += phase * (Complex64::from(g_m * g_m) * k0 + k1_rho);
            gxy += phase * Complex64::new(0.0, -g_m) * kappa * k1;
        }
        (gxx * norm, gyy * norm, gxy * norm)
    }

    /// **Tensor-structure decomposition check -- NOT physics validation.**
    ///
    /// The reference below re-encodes the same formula this function implements, so
    /// it cannot catch a wrong formula: when `c_cross` lost its spurious /rho, the
    /// reference had to be corrected in the same commit. What it *does* catch is a
    /// mis-assembly of c_delta/c_ee/c_cross/c_rr into tensor slots, which is a real
    /// and easy mistake.
    ///
    /// Absolute normalisation is validated by `test_1d_ewald_all_elements_vs_direct_sum`,
    /// which compares against an independent brute-force lattice sum. Do not read this
    /// test as evidence that the spectral formula is right.
    #[test]
    fn test_1d_recip_tensor_elements_xx_yy_xy() {
        // Compare tensor elements against a direct transcription of each sum formula,
        // without the c_delta/c_ee/c_cross/c_rr abstraction.
        //
        // Deltax must NOT be a1/2: there g_m*Deltax = mpi for every mode, sin(g_m Deltax) = 0,
        // and G_xy vanishes identically -- the element then compares 0 against 0 and
        // validates nothing. (That is how the spurious /rho in c_cross survived until
        // test_1d_ewald_all_elements_vs_direct_sum caught it.) Deltax = 2.5 gives
        // |G_xy| ~= 1.8e-3, comparable to |G_xx|.
        let k  = 0.01_f64;  // nm^-1, lambda ~= 628 nm
        let a1 = [10.0_f64, 0.0, 0.0];
        let r  = [2.5_f64, 3.0, 0.0]; // Deltax = 2.5 nm, rho = 3 nm

        let g = recip_space_1d(r, [0.0; 3], k, a1);

        let (ref_xx, ref_yy, ref_xy) =
            recip_1d_reference_elements(r[0], r[1], k, a1[0], 50);

        let tol = 1e-8;

        assert!(
            (g[0][0] - ref_xx).norm() < tol,
            "G_xx mismatch: got {:?}, expected {:?}, diff = {:.2e}",
            g[0][0], ref_xx, (g[0][0] - ref_xx).norm()
        );
        assert!(
            (g[1][1] - ref_yy).norm() < tol,
            "G_yy mismatch: got {:?}, expected {:?}, diff = {:.2e}",
            g[1][1], ref_yy, (g[1][1] - ref_yy).norm()
        );
        assert!(
            (g[0][1] - ref_xy).norm() < tol,
            "G_xy mismatch: got {:?}, expected {:?}, diff = {:.2e}",
            g[0][1], ref_xy, (g[0][1] - ref_xy).norm()
        );
        // Guard against the comparison above going vacuous if the geometry is
        // ever changed: G_xy must be large enough for the 1e-8 tolerance to bite.
        assert!(
            ref_xy.norm() > 1e-6,
            "G_xy reference is negligible ({:.2e}) -- Deltax may be a1/2 again, in which \
             case this element asserts nothing",
            ref_xy.norm()
        );
        // Symmetry: G[0][1] == G[1][0]
        assert!(
            (g[0][1] - g[1][0]).norm() < 1e-12,
            "G_xy != G_yx: off-diagonal symmetry broken"
        );
    }

    /// C1: 1D Ewald magnitude validation against brute-force direct lattice sum.
    ///
    /// The existing tensor-element tests compare two implementations of the same
    /// spectral formula, which catches structural bugs but not normalisation errors.
    /// This test compares the *combined* Ewald sum (real_space + recip_space_1d)
    /// against an independent brute-force sum:
    ///
    ///   G_direct(r) = sum_{n=-200}^{200} G_free(r - n*a1) exp(i k_bloch * n*a1)
    ///
    /// The two sums must agree to within 1 % on the dominant tensor elements.
    /// A normalisation bug in recip_space_1d would appear as a multi-percent
    /// discrepancy here.
    #[test]
    fn test_1d_ewald_magnitude_vs_direct_sum() {
        use std::f64::consts::PI;

        let a  = 100.0_f64; // chain period (nm)
        let k  = 2.0 * PI / 600.0_f64;  // lambda = 600 nm (nm^-1)
        let a1 = [a, 0.0, 0.0];
        // Zone-boundary k_bloch: Bloch phase = (-1)^n -> absolute convergence
        // of the direct sum for near-field dominant terms.
        let kb = [PI / a, 0.0, 0.0];

        // Test two transverse displacements (x_par = 0, different rho).
        // G_yy is the dominant transverse element and unambiguously non-zero.
        for (r, label) in [
            ([0.0_f64, 20.0, 0.0], "rho=20nm"),
            ([0.0_f64, 40.0, 0.0], "rho=40nm"),
        ] {
            let lat     = LatticeSpec::chain(a1);
            let eg      = EwaldGreens::new(lat);
            let g_ewald = eg.evaluate(r, kb, k);

            // Brute-force direct sum with N = 200 shells.
            let g_direct = direct_sum_1d(r, kb[0], k, a, 200);

            // Compare G_yy (dominant transverse element).
            let ewald  = g_ewald[1][1];
            let direct = g_direct[1][1];
            let mag    = direct.norm();
            assert!(mag > 1e-10, "G_yy direct sum is near-zero ({label}): {mag:.2e}");
            let rel_err = (ewald - direct).norm() / mag;
            assert!(
                rel_err < 0.01,
                "1D Ewald G_yy magnitude mismatch ({label}): \
                 Ewald={ewald:.4e}, Direct={direct:.4e}, rel_err={rel_err:.3e}",
            );
        }
    }

    /// Brute-force Bloch-periodic sum of the free-space dyadic Green's function
    /// for a chain along x_hat with period `a`:
    ///
    ///   G_direct(**r**) = sum_{n=-N}^{N} G_free(**r** - n a x_hat) exp(i k_b n a)
    ///
    /// Shares no code with the spectral (Bessel) sum, so agreement validates the
    /// absolute normalisation of `recip_space_1d`, not merely its tensor structure.
    fn direct_sum_1d(
        r: [f64; 3],
        kb_x: f64,
        k: f64,
        a: f64,
        n_shells: i64,
    ) -> [[Complex64; 3]; 3] {
        let mut g = [[Complex64::from(0.0); 3]; 3];
        for n in -n_shells..=n_shells {
            let rn_x = n as f64 * a;
            let s = [r[0] - rn_x, r[1], r[2]];
            let s_sq = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
            if s_sq < 1e-20 { continue; }
            let arg = kb_x * rn_x;
            let phase = Complex64::new(arg.cos(), arg.sin());
            let g_s = dyadic_greens_tensor_at(&s, k);
            for i in 0..3 {
                for j in 0..3 {
                    g[i][j] += phase * g_s[i][j];
                }
            }
        }
        g
    }

    /// Largest element magnitude of a tensor -- the scale errors are measured against,
    /// so that near-zero elements do not produce meaningless relative errors.
    fn max_elem(g: &[[Complex64; 3]; 3]) -> f64 {
        let mut m = 0.0_f64;
        for row in g.iter() {
            for e in row.iter() {
                m = m.max(e.norm());
            }
        }
        m
    }

    /// C1: full-tensor magnitude validation of the 1D spectral sum.
    ///
    /// `test_1d_ewald_magnitude_vs_direct_sum` checks only G_yy at zero axial
    /// offset. This test covers every tensor element at non-zero axial offsets,
    /// which additionally exercises the exp(i g_m x_par) phase and the `c_cross`
    /// (e_hat(x)rho_hat + rho_hat(x)e_hat) term -- the term that carried a spurious 1/rho until this test
    /// caught it, and the only place the transverse unit vector rho_hat enters linearly.
    ///
    /// One geometry is deliberately out of the xy-plane with a small rho/a1: that
    /// combination stresses the mode-count `m_max`, since the spectral terms decay
    /// as exp(-2*pi*m*rho/a1) and a rho-independent m_max under-resolves small rho.
    ///
    /// Bloch vector is at the zone boundary (k_b = pi/a > k), so every spectral
    /// mode is evanescent and the direct sum converges. The propagating-mode
    /// case (k_b < k) is only conditionally convergent in real space and is not
    /// validated here.
    #[test]
    fn test_1d_ewald_all_elements_vs_direct_sum() {
        use std::f64::consts::PI;

        let a  = 100.0_f64; // chain period (nm)
        let k  = 2.0 * PI / 600.0_f64;  // lambda = 600 nm (gold plasmonic regime)
        let a1 = [a, 0.0, 0.0];
        let kb = [PI / a, 0.0, 0.0];    // zone boundary -- all modes evanescent

        for r in [
            [30.0_f64, 25.0, 0.0],  // in-plane,     rho/a = 0.25
            [10.0_f64, 12.0, 8.0],  // out-of-plane, rho/a = 0.144
            [70.0_f64, 40.0, 0.0],  // in-plane,     rho/a = 0.40
        ] {
            // The reference must itself be converged before it can validate anything.
            // Convergence is only ~1/N (far terms decay as 1/R with an oscillating
            // phase), so it needs thousands of shells. At 12800 the reference is good
            // to ~5e-6, keeping it a decade clear of the 1e-4 assertion below -- at
            // 6400 the measured error was dominated by the reference, not by Ewald.
            let d_half  = direct_sum_1d(r, kb[0], k, a, 6400);
            let d_full  = direct_sum_1d(r, kb[0], k, a, 12800);
            let scale   = max_elem(&d_full);
            assert!(scale > 1e-10, "direct sum is degenerate at r={r:?}: {scale:.2e}");
            for i in 0..3 {
                for j in 0..3 {
                    let drift = (d_full[i][j] - d_half[i][j]).norm() / scale;
                    assert!(
                        drift < 1e-4,
                        "direct sum not converged at r={r:?} element [{i}][{j}]: 6400 vs 12800 shells differ by {drift:.2e}"
                    );
                }
            }

            let g_ewald = EwaldGreens::new(LatticeSpec::chain(a1)).evaluate(r, kb, k);

            for i in 0..3 {
                for j in 0..3 {
                    let rel_err = (g_ewald[i][j] - d_full[i][j]).norm() / scale;
                    assert!(
                        rel_err < 1e-4,
                        "1D Ewald element [{i}][{j}] disagrees with direct sum at r={r:?}: Ewald={:.6e}, Direct={:.6e}, rel_err={rel_err:.3e}",
                        g_ewald[i][j], d_full[i][j]
                    );
                }
            }

            // The cross term must be genuinely present, or the comparison above
            // would pass vacuously for a G with no e_hat(x)rho_hat structure at all.
            let cross = g_ewald[0][1].norm() / scale;
            assert!(
                cross > 1e-3,
                "G_xy is negligible at r={r:?} ({cross:.2e} of max element) -- the cross term is not exercised, so c_cross is not validated"
            );
        }
    }


    /// The mode cap must fail loudly, not return a plausible-looking wrong number.
    /// rho/a1 = 0.005 needs ~637 modes against a cap of 200, where the result is wrong
    /// by ~300 % and grows spurious imaginary parts.
    ///
    #[test]
    #[should_panic(expected = "mode-cap limited")]
    fn test_1d_recip_cap_limited_panics() {
        use std::f64::consts::PI;
        let a = 100.0_f64;
        let k = 2.0 * PI / 600.0;
        // rho = 0.5 nm against a 100 nm period -- the regime a finely discretised
        // chain actually lands in.
        let _ = recip_space_1d([10.0, 0.5, 0.0], [PI / a, 0.0, 0.0], k, [a, 0.0, 0.0]);
    }


    /// C3: a chain of dipoles lying ON the axis must still feel the lattice.
    ///
    /// The Bessel spectral sum is only valid for rho > 0 -- K0(kapparho) diverges
    /// logarithmically for every mode as rho -> 0 -- so `recip_space_1d` returns zero
    /// there. For a long time `real_space` also returned zero for 1D, which made
    /// `evaluate` return an exactly zero tensor for *every* collinear pair, not just
    /// the self-block: an on-axis chain had no periodic coupling at all, and a
    /// discretised chain silently dropped its self-image sum.
    #[test]
    fn test_1d_collinear_pair_is_not_zero() {
        use std::f64::consts::PI;
        let a = 100.0_f64;
        let k = 2.0 * PI / 600.0;
        let kb = [PI / a, 0.0, 0.0];
        let eg = EwaldGreens::new(LatticeSpec::chain([a, 0.0, 0.0]));

        // Two distinct dipoles on the axis, 30 nm apart.
        let g = eg.evaluate([30.0, 0.0, 0.0], kb, k);

        let scale = max_elem(&g);
        assert!(
            scale > 1e-12,
            "collinear pair returned a ~zero tensor ({scale:.2e}) -- the lattice has no effect"
        );

        // It must agree with a brute-force Bloch sum of the free-space tensor.
        let direct = direct_sum_1d([30.0, 0.0, 0.0], kb[0], k, a, 20000);
        let ref_scale = max_elem(&direct);
        for i in 0..3 {
            for j in 0..3 {
                let err = (g[i][j] - direct[i][j]).norm() / ref_scale;
                assert!(
                    err < 1e-3,
                    "collinear element [{i}][{j}] disagrees with the direct sum: got {:?}, want {:?}, rel {err:.2e}",
                    g[i][j], direct[i][j]
                );
            }
        }
    }

    /// The self-block (r = 0) must still carry the sum over non-zero lattice images --
    /// usually the dominant periodic contribution for a discretised chain.
    #[test]
    fn test_1d_self_image_sum_is_not_zero() {
        use std::f64::consts::PI;
        let a = 100.0_f64;
        let k = 2.0 * PI / 600.0;
        let eg = EwaldGreens::new(LatticeSpec::chain([a, 0.0, 0.0]));

        let g = eg.evaluate([0.0, 0.0, 0.0], [PI / a, 0.0, 0.0], k);
        assert!(
            max_elem(&g) > 1e-12,
            "self-image lattice sum is zero -- the R!=0 images were dropped"
        );
    }


    /// Cross-validates the two representations against each other.
    ///
    /// `real_space_1d_direct` (used for rho ~= 0) and `recip_space_1d` (used for rho > 0)
    /// are independent constructions of the same periodic Green's function, so in the
    /// overlap region where both are valid they must agree. Without this, the
    /// collinear branch would rest on nothing but its own convergence.
    #[test]
    fn test_1d_collinear_matches_spectral_sum() {
        use std::f64::consts::PI;
        let a = 100.0_f64;
        let k = 2.0 * PI / 600.0;
        let kb = [PI / a, 0.0, 0.0];

        for r in [[30.0_f64, 25.0, 0.0], [10.0_f64, 12.0, 8.0]] {
            let spectral = recip_space_1d(r, kb, k, [a, 0.0, 0.0]);
            let direct = real_space_1d_direct(r, kb, k, [a, 0.0, 0.0], COLLINEAR_SHELLS);

            let scale = max_elem(&spectral);
            assert!(scale > 1e-12, "spectral sum degenerate at r={r:?}");
            for i in 0..3 {
                for j in 0..3 {
                    let err = (spectral[i][j] - direct[i][j]).norm() / scale;
                    assert!(
                        err < 1e-4,
                        "the two representations disagree at r={r:?} element [{i}][{j}]: spectral {:?} vs direct {:?}, rel {err:.2e}",
                        spectral[i][j], direct[i][j]
                    );
                }
            }
        }
    }

}
