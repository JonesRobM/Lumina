# Mie Theory

Mie theory provides the exact analytical solution for electromagnetic scattering by a homogeneous sphere. It serves as the primary validation benchmark for the CDA solver.

## Formulation

For a sphere of radius \\(a\\) with complex refractive index \\(m = n_{\text{sphere}} / n_{\text{medium}}\\) and size parameter \\(x = k a\\), the scattering coefficients are:

\\[
a_n = \frac{m \psi_n(mx) \psi_n'(x) - \psi_n(x) \psi_n'(mx)}{m \psi_n(mx) \xi_n'(x) - \xi_n(x) \psi_n'(mx)}
\\]

\\[
b_n = \frac{\psi_n(mx) \psi_n'(x) - m \psi_n(x) \psi_n'(mx)}{\psi_n(mx) \xi_n'(x) - m \xi_n(x) \psi_n'(mx)}
\\]

where \\(\psi_n\\) and \\(\xi_n\\) are Riccati-Bessel functions.

## Cross-Sections

\\[
C_{\text{ext}} = \frac{2\pi}{k^2} \sum_{n=1}^{\infty} (2n + 1) \, \text{Re}(a_n + b_n)
\\]

\\[
C_{\text{sca}} = \frac{2\pi}{k^2} \sum_{n=1}^{\infty} (2n + 1) (|a_n|^2 + |b_n|^2)
\\]

The series is truncated when \\(|a_n|\\) and \\(|b_n|\\) fall below a tolerance (typically \\(10^{-12}\\)).

## Validation Protocol

The CDA is validated against Mie theory using a gold nanosphere (\\(R = 10\\) nm). The expected accuracy depends on the spectral region and discretisation:

| Regime | Wavelength | Spacing | Expected Error |
|--------|-----------|---------|---------------|
| Interband (Au) | 400–520 nm | \\(d = 3\\) nm | \\(< 30\%\\) |
| Drude (Au) | \\(> 550\\) nm | \\(d = 3\\) nm | \\(> 50\%\\) (known limitation) |
| Dielectric | any | \\(d = 3\\) nm | \\(< 15\%\\) |
| Dielectric | any | \\(d \leq 2\\) nm | \\(< 5\%\\) |

The 30% tolerance at coarse spacing (\\(d = 3\\) nm) is appropriate for a sphere with only \\(\sim 170\\) dipoles. Finer discretisation significantly improves accuracy for dielectrics. For highly metallic particles, the staircase surface artefact limits convergence regardless of spacing; see [v0.1 Validation](./validation.md) for details.

## Reference

Bohren, C. F. & Huffman, D. R., *Absorption and Scattering of Light by Small Particles*, Wiley (1983).
