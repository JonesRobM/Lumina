# Periodic Structures & Ewald Summation

> **Status (v0.4.0):** Full Ewald acceleration implemented — erfc-damped real-space sum + 2D spectral reciprocal-space sum.

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

**Real-space part**: Each free-space Green's tensor term is multiplied by `erfc(η·|r−R|)`, which exponentially damps distant images. With the auto-selected η = √(π/A), 5 shells are sufficient for < 0.01 % error.

**Reciprocal-space part** (2D lattices): Evaluates the spectral sum

\\[
\mathbf{G}^{\text{spec}}_{ab} = \frac{i}{2A} \sum_{\mathbf{G}}
    \frac{e^{i\mathbf{Q}\cdot\boldsymbol{\rho}}\,e^{-q|z|}}{q}\,
    M_{ab}(\mathbf{Q}, q)
\\]

where **Q** = **k**∥ + **G**, \\(q = \sqrt{Q^2 - k^2}\\) with Im(q) ≥ 0, and \\(M_{ab} = k^2\delta_{ab} - \tilde{k}_a\tilde{k}_b\\).

For 1D chains the reciprocal sum returns zero (the erfc-damped direct sum converges adequately).

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

## Known Limitations (v0.4.0)

| Scenario | Status |
|----------|--------|
| 2D periodic array, any unit cell size | ✓ Full Ewald (v0.4) |
| 1D chain | ✓ erfc-damped real-space sum; reciprocal sum not needed |
| Oblique incidence (arbitrary \\(\mathbf{k}_\parallel\\)) | ✓ Correct with Bloch phases |
| 3D periodicity | Not supported |
| Reciprocal sum for 1D chain | Deferred (direct sum adequate) |
