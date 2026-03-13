# Periodic Structures & Ewald Summation

> **Status (v0.4.1):** 2D lattices — full Ewald acceleration (erfc real-space + spectral reciprocal sum). 1D chains — direct lattice sum with ≥ 20 shells (no erfc, ~2 % convergence at optical wavelengths). Bessel-function 1D spectral sum deferred.

## Motivation

Many technologically important structures — semiconductor slabs, photonic crystals, metasurfaces, and gratings — are periodic in one or two spatial dimensions. The CDA must account for interactions between a dipole and all its periodic images.

## The Convergence Problem

The naive lattice sum of the dyadic Green's function over periodic images:

\\[
\mathbf{G}_{\text{lat}}(\mathbf{r}, \mathbf{r}') = \sum_{\mathbf{R}_n} \mathbf{G}(\mathbf{r}, \mathbf{r}' + \mathbf{R}_n)
\\]

is **conditionally convergent** — it converges extremely slowly (if at all) in its original form. Direct summation over many shells is feasible only for small unit cells.

## Ewald Summation

Ewald's method resolves this by splitting the sum into two rapidly convergent parts:

1. **Real-space sum**: Short-range interactions computed in direct space, damped by a Gaussian screening function \\(e^{-\eta^2 R^2}\\).
2. **Reciprocal-space sum**: Long-range interactions computed in Fourier space, which converge exponentially.

The splitting parameter \\(\eta\\) controls the balance between the two sums. Both sums converge exponentially, making the total computation \\(O(N)\\) per unit cell.

## Physical Nuance

Ewald summation does not approximate the physics — it is an exact mathematical rearrangement of the lattice sum. The result is identical to the infinite direct sum; only the computational strategy differs.

The method naturally handles:
- **1D periodicity** (nanowires, gratings)
- **2D periodicity** (slabs, metasurfaces)
- **3D periodicity** (bulk photonic crystals)

Oblique incidence is accommodated via the Bloch phase factor \\(e^{i\mathbf{k}_\parallel \cdot \mathbf{R}_n}\\).

## Current Implementation (v0.3.0)

### LatticeSpec

Defined in `lumina-geometry/src/lattice.rs`. Specifies the Bravais lattice vectors:

```toml
# 1D chain along x, period 50 nm
[lattice]
type = "chain"
a1 = [50.0, 0.0, 0.0]

# 2D square lattice, 50 nm × 50 nm
[lattice]
type = "planar"
a1 = [50.0, 0.0, 0.0]
a2 = [0.0, 50.0, 0.0]
```

`LatticeSpec` exposes:
- `real_lattice_points(n_shells)` — all lattice vectors within `n_shells` shells
- `reciprocal_vectors()` — primitive reciprocal lattice vectors
- `eta()` — auto-selected splitting parameter (\\(\eta = \sqrt{\pi} / L_{\text{cell}}\\))

### EwaldGreens

`lumina-core::solver::ewald::EwaldGreens` implements the full Ewald-accelerated periodic interaction tensor. The evaluate method returns:

\\[
\mathbf{G}_{\text{per}} = \mathbf{G}_{\text{real}} + \mathbf{G}_{\text{recip}}
\\]

**Real-space part (2D)**: Each free-space Green's tensor term is multiplied by `erfc(η·|r−R|)`, which exponentially damps distant images. With η = √(π/A), 5 shells give < 0.01 % error.

**Real-space part (1D chains)**: No erfc damping is applied. A direct lattice sum over at least 20 shells is used instead:

\\[
\mathbf{G}_{\text{real}}^{1D}(\mathbf{r}, k_\parallel) = \sum_{n = -N_\text{shells}}^{N_\text{shells}} \mathbf{G}_\text{free}(\mathbf{r} - n\mathbf{a}_1)\, e^{i n k_\parallel a_1}
\\]

With 20 shells this gives ~2 % accuracy at visible wavelengths for chains with period ≥ 50 nm. Use `truncation_real = 50` in the TOML for < 0.5 % error at the cost of more computation.

> **Why no erfc for 1D?** The optimal Ewald splitting parameter η = √(π/|**a**₁|) was derived for 2D. Applied to 1D chains it gives η·|**a**₁| ≈ 1.77, making `erfc(η·|**a**₁|) ≈ 0.004` — 99.6 % of the nearest-image interaction would be pushed into the reciprocal sum. The 1D reciprocal sum requires modified Bessel functions K₀/K₁ with complex arguments and is deferred to a future release. The direct sum is a correct, if slower-converging, alternative.

**Reciprocal-space part (2D lattices)**: Evaluates the spectral sum

\\[
\mathbf{G}^{\text{spec}}_{ab} = \frac{i}{2A} \sum_{\mathbf{G}}
    \frac{e^{i\mathbf{Q}\cdot\boldsymbol{\rho}}\,e^{-q|z|}}{q}\,
    M_{ab}(\mathbf{Q}, q)
\\]

where **Q** = **k**∥ + **G**, \\(q = \sqrt{Q^2 - k^2}\\) with Im(q) ≥ 0, and \\(M_{ab} = k^2\delta_{ab} - \tilde{k}_a\tilde{k}_b\\).

**Reciprocal-space part (1D chains)**: Returns zero. The direct real-space sum carries the full interaction.

### Convergence for 1D chains

With the default `truncation_real = 5` (auto-clamped to 20 for 1D), the direct sum error is dominated by the truncated tail. For a chain with period 80 nm at λ = 500 nm:

| `truncation_real` | Approx. error | Wall time (relative) |
|-------------------|---------------|----------------------|
| 20 | ~2 % | 1× |
| 50 | ~0.5 % | 2.5× |
| 100 | ~0.2 % | 5× |

### Bloch Boundary Conditions

Set the Bloch wavevector \\(\mathbf{k}_\parallel\\) in the simulation section:

```toml
[simulation]
k_bloch = [0.01, 0.0, 0.0]   # in nm⁻¹
```

This multiplies each lattice-image Green's tensor by the phase factor \\(e^{i \mathbf{k}_\parallel \cdot \mathbf{R}_n}\\), implementing exact Bloch boundary conditions for oblique incidence.

### Dispersion Maps

When sweeping over multiple \\(\mathbf{k}_\parallel\\) values, the output type is `DispersionMap`:

| Field | Description |
|-------|-------------|
| `wavelengths` | λ grid (nm) |
| `k_points` | \\(k_\parallel\\) grid (nm⁻¹) |
| `extinction` | C\_ext(λ, k) matrix |

## Known Limitations (v0.4.1)

| Scenario | Status |
|----------|--------|
| 2D periodic array, any unit cell size | ✓ Full Ewald (erfc real-space + spectral recip) |
| 1D chain | ✓ Direct lattice sum, ≥ 20 shells, ~2 % accuracy |
| Oblique incidence (arbitrary \\(\mathbf{k}_\parallel\\)) | ✓ Bloch phases on all sums |
| 3D periodicity | Not supported |
| 1D Bessel spectral sum (K₀/K₁) | Deferred — needed for < 1 % on-axis chains |
