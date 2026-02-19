# Coupled Dipole Approximation

The Coupled Dipole Approximation (CDA), also known as the Discrete Dipole Approximation (DDA), models an arbitrary nanostructure as an ensemble of \\(N\\) polarisable point dipoles. Each dipole interacts with all others through the electromagnetic field, and the self-consistent response is obtained by solving a \\(3N \times 3N\\) complex linear system.

## Formulation

### Self-Consistent Equation

The dipole moment \\(\mathbf{p}_i\\) at site \\(i\\) satisfies:

\\[
\mathbf{p}_i = \alpha_i \left( \mathbf{E}_{\text{inc}}(\mathbf{r}_i) + \sum_{j \neq i} \mathbf{G}(\mathbf{r}_i, \mathbf{r}_j) \cdot \mathbf{p}_j \right)
\\]

where:
- \\(\alpha_i\\) is the polarisability of dipole \\(i\\),
- \\(\mathbf{E}_{\text{inc}}\\) is the incident electric field,
- \\(\mathbf{G}\\) is the dyadic Green's function (see [Dyadic Green's Function](./greens-function.md)).

### Matrix Form

Rearranging into matrix form:

\\[
\mathbf{A} \mathbf{p} = \mathbf{E}_{\text{inc}}
\\]

where the \\(3 \times 3\\) blocks of \\(\mathbf{A}\\) are:
- Diagonal: \\(\mathbf{A}_{ii} = \alpha_i^{-1}\\)
- Off-diagonal: \\(\mathbf{A}_{ij} = -\mathbf{G}(\mathbf{r}_i, \mathbf{r}_j)\\)

### Polarisability

For a dipole representing a small volume \\(V = d^3\\) of material with dielectric function \\(\epsilon(\omega)\\) in a medium with \\(\epsilon_m\\), the Clausius-Mossotti polarisability is:

\\[
\alpha_{\text{CM}} = 3V\epsilon_0 \frac{\epsilon - \epsilon_m}{\epsilon + 2\epsilon_m}
\\]

Radiative corrections (Draine & Goodman, 1993) are applied to ensure energy conservation.

## Validity Criteria

The CDA is accurate provided:
1. The dipole spacing \\(d\\) satisfies \\(|m| k d < 1\\), where \\(m\\) is the complex refractive index.
2. A sufficient number of dipoles is used to resolve the geometry (typically \\(N > 10^3\\) for quantitative accuracy).

## Cross-Section Computation

From the solved dipole moments, the optical cross-sections are:

- **Extinction**: \\(C_{\text{ext}} = \frac{4\pi k}{|E_0|^2} \sum_i \text{Im}(\mathbf{E}_{\text{inc},i}^* \cdot \mathbf{p}_i)\\)
- **Absorption**: \\(C_{\text{abs}} = \frac{4\pi k}{|E_0|^2} \sum_i \left[ \text{Im}(\mathbf{p}_i \cdot (\alpha_i^{-1})^* \cdot \mathbf{p}_i^*) - \frac{2}{3} k^3 |\mathbf{p}_i|^2 \right]\\)
- **Scattering**: \\(C_{\text{sca}} = C_{\text{ext}} - C_{\text{abs}}\\)
