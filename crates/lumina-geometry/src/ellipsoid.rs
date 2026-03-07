//! Ellipsoid depolarisation factors via Gauss–Legendre quadrature.
//!
//! The depolarisation factors $L_x, L_y, L_z$ for an ellipsoid with semi-axes
//! $[a, b, c]$ are defined by:
//!
//! $$L_i = \frac{abc}{2} \int_0^\infty
//!   \frac{ds}{(s_i^2 + s)\sqrt{(a^2+s)(b^2+s)(c^2+s)}}$$
//!
//! where $s_x = a^2$, $s_y = b^2$, $s_z = c^2$ (the squared semi-axes).
//!
//! They satisfy $L_x + L_y + L_z = 1$ and reduce to $L_i = 1/3$ for a sphere.
//!
//! # Numerical method
//!
//! The substitution $s = s_i^2(1 - u^2)/u^2$ (equivalently $u = s_i/\sqrt{s_i^2 + s}$)
//! maps $u \in [0, 1] \to s \in [\infty, 0)$ and transforms the integrand to
//! the remarkably clean form:
//!
//! $$h_i(u) = \frac{2u^2}{\sqrt{q_a(u)\,q_b(u)\,q_c(u)}}$$
//!
//! where $q_j(u) = s_i^2 + (a_j^2 - s_i^2)\,u^2$.
//!
//! For a sphere this reduces to $h(u) = 2u^2/a^3$, a polynomial, so 16-point
//! Gauss–Legendre quadrature gives the exact answer. For general ellipsoids
//! the integrand is smooth and convergence is rapid.

/// Eight positive Gauss–Legendre nodes on $[-1, 1]$ for the symmetric 16-point rule.
/// By symmetry the negative nodes are $-x_k$, both with weight $w_k$.
/// Source: Abramowitz & Stegun, Table 25.4.
const GL16_NODES_POS: [f64; 8] = [
    0.095_012_509_836_062_3,
    0.281_603_550_779_258_9,
    0.458_016_777_657_227_4,
    0.617_876_244_402_643_8,
    0.755_404_408_355_003_0,
    0.865_631_202_383_565_0,
    0.944_575_023_073_232_6,
    0.989_400_934_991_649_9,
];

/// Weights corresponding to the positive nodes above (same for $\pm x_k$).
const GL16_WEIGHTS: [f64; 8] = [
    0.189_450_610_455_068_5,
    0.182_603_415_044_923_6,
    0.169_156_519_395_002_5,
    0.149_595_988_816_576_7,
    0.124_628_971_255_533_9,
    0.095_158_511_682_492_8,
    0.062_253_523_938_647_9,
    0.027_152_459_411_754_1,
];

/// Compute the three depolarisation factors $[L_x, L_y, L_z]$ for an ellipsoid
/// with semi-axes `[a, b, c]`.
///
/// All three satisfy $0 < L_i < 1$ and $L_x + L_y + L_z = 1$.
/// For a sphere ($a = b = c$): all equal $1/3$.
/// For a prolate spheroid ($a > b = c$): $L_x < 1/3 < L_y = L_z$.
/// For an oblate spheroid ($a = b > c$): $L_x = L_y < 1/3 < L_z$.
///
/// # Arguments
/// * `semi_axes` — `[a, b, c]` in nm (or any consistent unit).
///
/// # Panics
/// Panics in debug mode if any semi-axis is ≤ 0.
pub fn ellipsoid_depol_factors(semi_axes: [f64; 3]) -> [f64; 3] {
    let [a, b, c] = semi_axes;
    debug_assert!(a > 0.0 && b > 0.0 && c > 0.0, "semi-axes must be positive");

    let a2 = a * a;
    let b2 = b * b;
    let c2 = c * c;
    let prefactor = a * b * c * 0.5;

    let lx = prefactor * gl_integrate(a2, b2, c2, a2);
    let ly = prefactor * gl_integrate(a2, b2, c2, b2);
    let lz = prefactor * gl_integrate(a2, b2, c2, c2);
    [lx, ly, lz]
}

/// Evaluate $\int_0^\infty \frac{ds}{(s_i^2 + s)\sqrt{(a^2+s)(b^2+s)(c^2+s)}}$
/// using the substitution $s = s_i^2(1-u^2)/u^2$ followed by 16-point GL.
///
/// The transformed integrand $h(u) = 2u^2/\sqrt{q_a(u)\,q_b(u)\,q_c(u)}$,
/// $q_j(u) = s_i^2 + (a_j^2 - s_i^2)\,u^2$, is smooth and polynomial-like
/// near $u=0$, yielding rapid convergence.
fn gl_integrate(a2: f64, b2: f64, c2: f64, si2: f64) -> f64 {
    let mut sum = 0.0;
    for k in 0..8 {
        let x = GL16_NODES_POS[k];
        // Factor of 0.5: maps the [-1,1] GL rule to [0,1].
        let w = GL16_WEIGHTS[k] * 0.5;
        sum += w * gl_integrand(a2, b2, c2, si2, (1.0 + x) * 0.5);
        sum += w * gl_integrand(a2, b2, c2, si2, (1.0 - x) * 0.5);
    }
    sum
}

/// Evaluate $h_i(u) = 2u^2 / \sqrt{q_a \cdot q_b \cdot q_c}$
/// where $q_j(u) = s_i^2 + (a_j^2 - s_i^2)\,u^2$.
#[inline]
fn gl_integrand(a2: f64, b2: f64, c2: f64, si2: f64, u: f64) -> f64 {
    let u2 = u * u;
    let qa = si2 + (a2 - si2) * u2;
    let qb = si2 + (b2 - si2) * u2;
    let qc = si2 + (c2 - si2) * u2;
    2.0 * u2 / (qa * qb * qc).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(a: f64, b: f64, eps: f64, msg: &str) {
        assert!(
            (a - b).abs() < eps,
            "{msg}: got {a:.12}, expected {b:.12}, diff {:.2e}",
            (a - b).abs()
        );
    }

    #[test]
    fn test_depol_sphere() {
        // For a sphere, GL integrates u² exactly → result is exact.
        let [lx, ly, lz] = ellipsoid_depol_factors([10.0, 10.0, 10.0]);
        assert_close(lx, 1.0 / 3.0, 1e-12, "sphere Lx");
        assert_close(ly, 1.0 / 3.0, 1e-12, "sphere Ly");
        assert_close(lz, 1.0 / 3.0, 1e-12, "sphere Lz");
    }

    #[test]
    fn test_depol_sum_unity() {
        // Note: 16-point GL is accurate to ~1 ppm for aspect ratios up to ~20:1.
        // Beyond that (e.g. 100:1) a composite rule would be needed, but such
        // extreme ellipsoids are not physically meaningful for CDA dipoles.
        for axes in [
            [10.0_f64, 5.0, 3.0],
            [20.0, 20.0, 5.0],
            [4.0, 8.0, 16.0],
            [8.0, 8.0, 8.0],
        ] {
            let [lx, ly, lz] = ellipsoid_depol_factors(axes);
            // Each L_i uses a different substitution; errors accumulate rather
            // than cancel. ~1 ppm accuracy on the sum is fine for CDA physics.
            assert_close(lx + ly + lz, 1.0, 1e-6, &format!("sum for {axes:?}"));
        }
    }

    #[test]
    fn test_depol_prolate_ordering() {
        // Prolate: a (long axis) has smaller L than the transverse axes.
        let [lx, ly, lz] = ellipsoid_depol_factors([30.0, 10.0, 10.0]);
        assert!(lx < 1.0 / 3.0, "prolate Lx={lx:.8} should be < 1/3");
        assert!(ly > 1.0 / 3.0, "prolate Ly={ly:.8} should be > 1/3");
        assert_close(ly, lz, 1e-8, "prolate Ly == Lz by symmetry");
    }

    #[test]
    fn test_depol_oblate_ordering() {
        // Oblate: c (short axis) has larger L than the in-plane axes.
        let [lx, ly, lz] = ellipsoid_depol_factors([20.0, 20.0, 5.0]);
        assert!(lz > 1.0 / 3.0, "oblate Lz={lz:.8} should be > 1/3");
        assert_close(lx, ly, 1e-8, "oblate Lx == Ly by symmetry");
    }

    #[test]
    fn test_depol_prolate_analytical() {
        // Prolate spheroid a=20, b=c=10.
        // Closed form (Bohren & Huffman, Appendix B):
        //   L_a = (1 - e²) / e² * (artanh(e) / e − 1)
        //   where e = sqrt(1 − (b/a)²)
        let a = 20.0_f64;
        let b = 10.0_f64;
        let e2 = 1.0 - (b / a).powi(2);
        let e = e2.sqrt();
        let artanh_e = 0.5 * ((1.0 + e) / (1.0 - e)).ln();
        let l_a_analytical = (1.0 - e2) / e2 * (artanh_e / e - 1.0);
        let [lx, _, _] = ellipsoid_depol_factors([a, b, b]);
        assert_close(lx, l_a_analytical, 1e-8, "prolate L_a analytical");
    }
}
