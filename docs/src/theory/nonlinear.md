# Nonlinear Optics (SHG/THG)

> **Status**: Planned for post-v0.1.

For materials with broken inversion symmetry, the second-order nonlinear polarisation generates light at twice the fundamental frequency (second-harmonic generation, SHG).

## Second-Harmonic Generation

The nonlinear dipole moment at frequency \\(2\omega\\) is:

\\[
\mathbf{p}_i(2\omega) = \epsilon_0 \chi^{(2)}_{ijk} : \mathbf{E}_{\text{loc},j}(\omega) \mathbf{E}_{\text{loc},k}(\omega)
\\]

where \\(\chi^{(2)}\\) is the second-order susceptibility tensor, which depends on the material's symmetry group.

## Implementation Approach

1. Solve the linear CDA system at frequency \\(\omega\\) to obtain local fields.
2. Compute the nonlinear source terms \\(\mathbf{p}_i(2\omega)\\) from the local fields.
3. Re-solve the CDA system at frequency \\(2\omega\\) with these sources.
4. Compute far-field radiation from the \\(2\omega\\) dipoles.
