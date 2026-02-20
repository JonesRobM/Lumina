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

## Unit Convention

Lumina uses a Gaussian-like convention where \\(\epsilon_0\\) is absorbed into the definitions. The key quantities and their relationships to the DDSCAT convention are:

| Quantity | Lumina convention | DDSCAT convention |
|----------|------------------|-------------------|
| Green's function | \\(G = \frac{k^2 e^{ikR}}{4\pi R} [\ldots]\\) | \\(G_{\text{DD}} = \frac{k^2 e^{ikR}}{R} [\ldots]\\) |
| Polarisability | \\(\alpha = 3V\frac{\epsilon - \epsilon_m}{\epsilon + 2\epsilon_m}\\) | \\(\alpha_{\text{DD}} = \frac{3V}{4\pi}\frac{\epsilon - \epsilon_m}{\epsilon + 2\epsilon_m}\\) |
| Cross-section | \\(C_{\text{ext}} = \frac{k}{\|E_0\|^2} \sum_i \operatorname{Im}(\mathbf{E}^*_i \cdot \mathbf{p}_i)\\) | \\(C_{\text{ext}} = \frac{4\pi k}{\|E_0\|^2} \sum_i \operatorname{Im}(\mathbf{E}^*_i \cdot \mathbf{p}_i)\\) |

The \\(4\pi\\) factors in \\(\alpha\\) and \\(G\\) cancel in the system \\(\mathbf{A}\mathbf{p} = \mathbf{E}\\), so the cross-section formula uses \\(k / \|E_0\|^2\\) rather than \\(4\pi k / \|E_0\|^2\\). Both conventions yield identical physical cross-sections.

## Polarisability Prescriptions

### Clausius-Mossotti (CM)

For a dipole representing a small volume \\(V = d^3\\) of material with dielectric function \\(\epsilon(\omega)\\) in a medium with \\(\epsilon_m\\):

\\[
\alpha_{\text{CM}} = 3V \frac{\epsilon - \epsilon_m}{\epsilon + 2\epsilon_m}
\\]

This is the bare electrostatic polarisability. It does not satisfy the optical theorem and must be corrected for radiative effects.

### Radiative Reaction Corrected (RRCM)

The radiative reaction correction (Draine, 1988) ensures energy conservation:

\\[
\alpha_{\text{RRCM}}^{-1} = \alpha_{\text{CM}}^{-1} - \frac{ik^3}{6\pi}
\\]

or equivalently:

\\[
\alpha_{\text{RRCM}} = \frac{\alpha_{\text{CM}}}{1 - \frac{ik^3}{6\pi} \alpha_{\text{CM}}}
\\]

This is the default prescription used in Lumina v0.1. It ensures that the extinction of a single dipole satisfies the optical theorem exactly, and produces accurate results for dielectric and moderately metallic particles.

### Lattice Dispersion Relation (LDR)

The LDR correction (Draine & Goodman, 1993) modifies the Clausius-Mossotti polarisability so that a cubic lattice of point dipoles correctly reproduces the bulk dielectric dispersion relation:

\\[
\alpha_{\text{LDR}}^{-1} = \alpha_{\text{CM}}^{-1} - \frac{ik^3}{6\pi} - \frac{k^2}{d} (b_1 + b_2 \epsilon + b_3 \epsilon S)
\\]

where \\(b_1 = -1.8915316\\), \\(b_2 = 0.1648469\\), \\(b_3 = -1.7700004\\), and \\(S = \sum_i (\hat{a}_i \hat{k}_i)^2\\) depends on the propagation and polarisation directions.

> **Note:** In Lumina's convention, the LDR lattice correction term includes a \\(1/(4\pi)\\) factor relative to the original Draine & Goodman coefficients, because those coefficients were derived for the DDSCAT convention where \\(\alpha_{\text{DD}} = \alpha / (4\pi)\\). In practice, the LDR correction is numerically negligible for the particle sizes tested in v0.1 (the correction terms are \\(\sim 10^{-6}\\) of \\(\alpha^{-1}\\)). Both RRCM and LDR give equivalent results.

## Cross-Section Computation

From the solved dipole moments, the optical cross-sections are:

- **Extinction**: \\(C_{\text{ext}} = \frac{k}{|E_0|^2} \sum_i \operatorname{Im}(\mathbf{E}_{\text{inc},i}^* \cdot \mathbf{p}_i)\\)
- **Absorption**: \\(C_{\text{abs}} = \frac{k}{|E_0|^2} \sum_i \left[ \operatorname{Im}(\mathbf{p}_i \cdot (\alpha_i^{-1})^* \cdot \mathbf{p}_i^*) - \frac{2}{3} k^3 |\mathbf{p}_i|^2 \right]\\)
- **Scattering**: \\(C_{\text{sca}} = C_{\text{ext}} - C_{\text{abs}}\\)

See [Validation](./validation.md) for accuracy benchmarks against Mie theory.

## Validity Criteria

The CDA is accurate provided:
1. The dipole spacing \\(d\\) satisfies \\(|m| k d < 1\\), where \\(m\\) is the complex refractive index.
2. A sufficient number of dipoles is used to resolve the geometry (typically \\(N > 10^3\\) for quantitative accuracy).
3. The lattice is centred on the shape centre to avoid asymmetric staircase artefacts (see [Geometry & Shapes](../user-guide/geometry.md)).

## References

- Purcell, E. M. & Pennypacker, C. R., *Astrophys. J.* **186**, 705 (1973).
- Draine, B. T. & Flatau, P. J., *J. Opt. Soc. Am. A* **11**, 1491 (1994).
- Draine, B. T. & Goodman, J., *Astrophys. J.* **405**, 685 (1993).
