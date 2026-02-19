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

The CDA result for a 20 nm gold nanosphere should reproduce the Mie theory extinction spectrum to within 5% at all wavelengths, provided the dipole spacing is sufficiently fine (\\(d \leq 2\\) nm).

## Reference

Bohren, C. F. & Huffman, D. R., *Absorption and Scattering of Light by Small Particles*, Wiley (1983).
