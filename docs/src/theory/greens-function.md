# Dyadic Green's Function

The dyadic Green's function \\(\mathbf{G}(\mathbf{r}, \mathbf{r}')\\) describes the electric field at position \\(\mathbf{r}\\) due to an oscillating point dipole at \\(\mathbf{r}'\\).

## Free-Space Form

\\[
\mathbf{G}(\mathbf{r}, \mathbf{r}') = \frac{k^2 e^{ikR}}{4\pi \epsilon_0 R}
\left[ \left(1 + \frac{ikR - 1}{k^2 R^2}\right)\mathbf{I}
+ \frac{3 - 3ikR - k^2 R^2}{k^2 R^2} \frac{\mathbf{R}\mathbf{R}^T}{R^2} \right]
\\]

where \\(\mathbf{R} = \mathbf{r} - \mathbf{r}'\\) and \\(R = |\mathbf{R}|\\).

## Physical Regimes

The Green's function contains contributions from three regimes:

### Near field (\\(kR \ll 1\\))
The \\(1/R^3\\) terms dominate. The field has the character of a static dipole field and is responsible for near-field enhancement effects (e.g. SERS, tip-enhanced spectroscopy).

### Intermediate field (\\(kR \sim 1\\))
The \\(1/R^2\\) terms contribute. This is the transition region between near and far field.

### Far field (\\(kR \gg 1\\))
The \\(1/R\\) terms dominate. The field is a propagating transverse wave. This regime governs scattering and extinction cross-sections.

## Symmetry Properties

The free-space Green's function obeys reciprocity:

\\[
\mathbf{G}(\mathbf{r}, \mathbf{r}') = \mathbf{G}^T(\mathbf{r}', \mathbf{r})
\\]

This is verified in the test suite as a fundamental consistency check.

## Implementation Notes

In Lumina, the Green's function is implemented in `lumina-core/src/solver/cda/greens.rs`. The factor \\(\epsilon_0\\) is absorbed into the polarisability definition for cleaner numerics.
