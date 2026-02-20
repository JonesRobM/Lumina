# v0.1 Validation

This chapter documents the validation of the Lumina CDA solver against Mie theory for homogeneous spheres. All results use the v0.1 implementation: direct LU factorisation (`faer`), Clausius-Mossotti polarisability with radiative reaction correction (RRCM), and centred cubic lattice discretisation.

## Methodology

### Reference: Mie Theory

Mie theory provides the exact analytical cross-sections for scattering by a homogeneous sphere (see [Mie Theory](./mie-theory.md)). It serves as the ground-truth benchmark. The Lumina implementation uses Riccati-Bessel functions with the series truncated when \\(|a_n|, |b_n| < 10^{-12}\\).

### Test Geometry

All validation tests use a sphere of radius \\(R\\) discretised into a centred cubic lattice with spacing \\(d\\). The number of dipoles scales as:

\\[
N \approx \frac{4}{3}\pi \left(\frac{R}{d}\right)^3
\\]

For the standard test case (\\(R = 10\\) nm, \\(d = 3\\) nm), \\(N \approx 170\\) dipoles.

### Material Data

Gold optical constants are taken from Johnson & Christy (1972), covering 188–892 nm (43 data points). The tabulated \\((n, k)\\) values are converted to the complex dielectric function:

\\[
\epsilon_1 = n^2 - k^2, \qquad \epsilon_2 = 2nk
\\]

and interpolated with natural cubic splines.

## Test Results

### Single-Dipole Rayleigh Limit

A single dipole with polarisability \\(\alpha_{\text{RRCM}}\\) for a sphere of volume \\(V = \frac{4}{3}\pi R^3\\) reproduces the Mie Rayleigh limit to machine precision:

| Material | \\(\epsilon\\) | CDA \\(C_{\text{ext}}\\) | Mie \\(C_{\text{ext}}\\) | Error |
|----------|---------------|--------------------------|--------------------------|-------|
| Dielectric | \\(4.0 + 0.5i\\) | \\(3.14 \times 10^1\\) | \\(3.14 \times 10^1\\) | 0.00% |

This confirms that the cross-section formula, polarisability prescription, and unit conventions are mutually consistent.

### Dielectric Sphere

A dielectric sphere (\\(R = 10\\) nm, \\(d = 3\\) nm, \\(N \approx 170\\)) is tested with several dielectric constants across three wavelengths. All errors are well below 30%.

| Material | \\(\epsilon\\) | \\(\lambda\\) (nm) | CDA \\(C_{\text{ext}}\\) | Mie \\(C_{\text{ext}}\\) | Error |
|----------|---------------|---------------------|--------------------------|--------------------------|-------|
| TiO\\(_2\\)-like | \\(4.0 + 0.0i\\) | 500, 600, 700 | — | — | \\(\lesssim 15\%\\) |
| Lossy dielectric | \\(4.0 + 0.5i\\) | 500, 600, 700 | — | — | \\(\lesssim 15\%\\) |
| High-index lossy | \\(6.0 + 1.0i\\) | 500, 600, 700 | — | — | \\(\lesssim 15\%\\) |
| Near resonance | \\(-1.0 + 3.0i\\) | 500, 600, 700 | — | — | \\(\lesssim 13\%\\) |

Finer discretisation (\\(d \leq 2\\) nm) reduces errors to \\(<5\%\\) for all dielectric cases.

### Gold Sphere — Interband Region

Gold in the interband region (400–520 nm) has large \\(\epsilon_2\\), giving moderate \\(|\epsilon_1|/\epsilon_2\\) ratios. The CDA is well-behaved here.

| \\(\lambda\\) (nm) | \\(\epsilon_1\\) | \\(\epsilon_2\\) | \\(|\epsilon_1|/\epsilon_2\\) | CDA Error |
|---------------------|-----------------|-----------------|-------------------------------|-----------|
| 420 | \\(-1.48\\) | \\(5.28\\) | 0.28 | \\(\sim 22\%\\) |
| 450 | \\(-1.66\\) | \\(5.29\\) | 0.31 | \\(\sim 18\%\\) |
| 480 | \\(-2.82\\) | \\(4.10\\) | 0.69 | \\(\sim 25\%\\) |
| 510 | \\(-4.83\\) | \\(2.28\\) | 2.12 | \\(\sim 28\%\\) |

All pass the 30% tolerance threshold at \\(d = 3\\) nm.

### Gold Sphere — Drude Region (Known Limitation)

At wavelengths above \\(\sim 550\\) nm, gold enters the Drude regime where \\(|\epsilon_1| \gg \epsilon_2\\). The CDA with RRCM fails to converge for these cases:

| \\(\lambda\\) (nm) | \\(\epsilon_1\\) | \\(\epsilon_2\\) | \\(|\epsilon_1|/\epsilon_2\\) | CDA Error |
|---------------------|-----------------|-----------------|-------------------------------|-----------|
| 550 | \\(-7.16\\) | \\(1.84\\) | 3.9 | \\(\sim 50\%\\) |
| 600 | \\(-9.64\\) | \\(1.38\\) | 7.0 | \\(> 100\%\\) |
| 700 | \\(-16.7\\) | \\(1.24\\) | 13.5 | \\(> 200\%\\) |

This is a well-known limitation of the DDA/CDA at coarse discretisation; see [Known Limitations](#known-limitations) below.

### Convergence with Spacing

For a dielectric sphere (\\(\epsilon = 4.0 + 0.5i\\), \\(R = 8\\) nm, \\(\lambda = 550\\) nm), the CDA converges monotonically as the lattice is refined:

| Spacing \\(d\\) (nm) | \\(N\\) | CDA \\(C_{\text{ext}}\\) | Mie \\(C_{\text{ext}}\\) | Error |
|-----------------------|--------|--------------------------|--------------------------|-------|
| 4.0 | \\(\sim 30\\) | — | — | \\(\sim 10\%\\) |
| 3.0 | \\(\sim 80\\) | — | — | \\(\sim 5\%\\) |
| 2.5 | \\(\sim 120\\) | — | — | \\(\sim 3\%\\) |
| 2.0 | \\(\sim 260\\) | — | — | \\(\sim 1.5\%\\) |

For metallic particles (\\(\epsilon = -10 + 1.5i\\)), convergence is non-monotonic and the errors remain large (\\(> 50\%\\)) even at \\(d = 1.5\\) nm with \\(N \approx 600\\) dipoles. This is characteristic of the staircase surface artefact described below.

## v0.1.1 Results: Filtered Coupled Dipole (FCD)

### FCD Method

The v0.1.1 release introduces the Filtered Coupled Dipole (FCD) method, also known as the Integrated Green's Tensor (IGT) approach (Yurkin & Hoekstra, JQSRT 2007). For near-field interactions (\\(|r_i - r_j| \leq 2d\\)), the point-dipole Green's function is replaced by its volume average over the source cell using 3D Gauss-Legendre quadrature (\\(3^3 = 27\\) evaluation points per cell):

\\[
G_{\text{FCD}}(\mathbf{r}_i, \mathbf{r}_j) = \frac{1}{V_j} \int_{V_j} G(\mathbf{r}_i, \mathbf{r}') \, d\mathbf{r}'
\\]

This smooths the \\(1/R^3\\) near-field singularity that amplifies surface staircase artefacts for metallic particles.

### GMRES Iterative Solver

v0.1.1 also adds a restarted GMRES(m) solver for systems with \\(N > 1000\\) dipoles, enabling finer discretisation. GMRES agrees with the direct LU solver to machine precision (\\(\sim 10^{-13}\\) relative error) and is activated automatically when \\(N\\) exceeds the iterative threshold (default 1000).

### FCD Results — Gold Sphere

The FCD method with \\(d = 2\\) nm (\\(N \approx 515\\) dipoles) shows improved accuracy in the interband region compared to the standard point-dipole CDA:

| \\(\lambda\\) (nm) | \\(\epsilon_1\\) | \\(\epsilon_2\\) | \\(|\epsilon_1|/\epsilon_2\\) | Point CDA Error | FCD Error |
|---------------------|-----------------|-----------------|-------------------------------|-----------------|-----------|
| 420 | \\(-1.70\\) | \\(5.71\\) | 0.3 | \\(\sim 25\%\\) | \\(\sim 16\%\\) |
| 480 | \\(-1.78\\) | \\(4.55\\) | 0.4 | \\(\sim 22\%\\) | \\(\sim 13\%\\) |
| 510 | \\(-3.15\\) | \\(3.06\\) | 1.0 | \\(\sim 28\%\\) | \\(\sim 24\%\\) |
| 550 | \\(-5.93\\) | \\(2.10\\) | 2.8 | \\(\sim 50\%\\) | \\(\sim 81\%\\) |
| 600 | \\(-9.42\\) | \\(1.50\\) | 6.3 | \\(> 100\%\\) | \\(\sim 112\%\\) |

The FCD provides a meaningful improvement (\\(5\text{–}10\\) percentage points) in the interband region (\\(|\epsilon_1|/\epsilon_2 < 3\\)). However, in the Drude region (\\(|\epsilon_1|/\epsilon_2 > 5\\)), the errors remain large. This confirms that the staircase artefact is not fully resolved by volume-averaging alone — future releases will investigate surface averaging and smooth boundary representations.

### Geometry Completions

v0.1.1 adds containment and bounding-box implementations for:

- **Cylinder**: Defined by base centre, axis direction, length, and radius. Containment tests via axis projection and perpendicular distance.
- **Helix**: Defined by base centre, axis, helix radius, pitch, turns, and wire radius. Containment tests via nearest-point search on the helical centreline.

Both shapes produce correct dipole lattices via the existing centred cubic discretisation.

## Known Limitations

### Surface Staircase Artefacts

When a smooth surface (e.g. a sphere) is discretised onto a cubic lattice, dipoles near the surface occupy cubic cells that straddle the boundary. This creates a "staircase" surface with artificially sharp edges at the lattice scale. For highly metallic particles (\\(|\epsilon_1| \gg \epsilon_2\\)):

1. The strong near-field coupling (\\(\propto 1/R^3\\)) amplifies errors at adjacent surface dipoles.
2. The staircase geometry creates artificial polarisation charges at surface steps.
3. These artefacts do not diminish quickly with finer discretisation because each refinement creates a new staircase pattern.

### Accuracy Boundaries

The CDA with RRCM is reliable when the metallicity ratio \\(|\epsilon_1|/\epsilon_2 \lesssim 5\\). A practical guide:

| \\(|\epsilon_1|/\epsilon_2\\) | Accuracy | Example |
|-------------------------------|----------|---------|
| \\(< 1\\) | Excellent (\\(< 15\%\\)) | Dielectrics, Au at 420 nm |
| \\(1 \text{–} 3\\) | Good (\\(15\text{–}30\%\\)) | Au at 510 nm |
| \\(3 \text{–} 5\\) | Fair (\\(30\text{–}50\%\\)) | Au at 550 nm |
| \\(> 5\\) | Poor (\\(> 100\%\\)) | Au at 600+ nm, Ag in the UV |

### Mitigation Strategies (v0.2+)

The following improvements are planned for future releases:

- **Filtered Coupled Dipole (FCD)**: Replaces the free-space Green's function with a lattice-filtered version, eliminating aliasing artefacts from the staircase surface.
- **Surface averaging**: Assigns reduced polarisabilities to dipoles that partially overlap the surface boundary.
- **Smooth surface methods**: Interpolation-based boundary representations that avoid the staircase entirely.
- **GMRES iterative solver**: Enables much finer discretisation (\\(N > 10^4\\)) that partially compensates for the staircase.

## Test Suite

The v0.1 validation tests are located in `crates/lumina-core/tests/` and can be run with:

```bash
cargo test --workspace
```

### `mie_validation.rs`

Three integration tests:

1. **`test_cda_vs_mie_gold_sphere`** — Validates CDA against Mie theory for a gold sphere (\\(R = 10\\) nm, \\(d = 3\\) nm) at four wavelengths in the interband region (420–510 nm). Tolerance: 30%.

2. **`test_cda_vs_mie_dielectric_sphere`** — Validates CDA against Mie theory for a dielectric sphere with four different \\(\epsilon\\) values at three wavelengths (500–700 nm). Tolerance: 30%.

3. **`test_finer_discretisation_improves_accuracy`** — Confirms that refining the lattice spacing from 4 nm to 2 nm monotonically reduces the CDA error for a dielectric sphere (\\(\epsilon = 4.0 + 0.5i\\), \\(R = 8\\) nm).

### `cda_diagnostic.rs`

Diagnostic tests for debugging and characterisation:

1. **Single-dipole Rayleigh limit** — Verifies exact agreement with Mie theory for a single dipole (zero CDA error).

2. **Dielectric sphere** — Tests a well-conditioned dielectric case (\\(\epsilon = 4.0 + 0.5i\\)) where the CDA should converge.

3. **Metallic progression** — Sweeps \\(|\epsilon_1|/\epsilon_2\\) from 2 to 10, demonstrating the systematic increase in CDA error for increasingly metallic particles.

### `cda_convergence.rs`

Convergence comparison between RRCM and LDR polarisability prescriptions for dielectric and metallic particles at spacings from 4 nm to 1.5 nm. Confirms that both prescriptions give equivalent results in Lumina's unit convention.

## References

- Johnson, P. B. & Christy, R. W., *Phys. Rev. B* **6**, 4370 (1972).
- Draine, B. T. & Flatau, P. J., *J. Opt. Soc. Am. A* **11**, 1491 (1994).
- Draine, B. T. & Goodman, J., *Astrophys. J.* **405**, 685 (1993).
- Yurkin, M. A. & Hoekstra, A. G., *J. Quant. Spectrosc. Radiat. Transf.* **106**, 558 (2007).
