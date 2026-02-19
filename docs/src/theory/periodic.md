# Periodic Structures & Ewald Summation

> **Status**: Planned for v0.3+.

## Motivation

Many technologically important structures — semiconductor slabs, photonic crystals, metasurfaces, and gratings — are periodic in one or two spatial dimensions. The CDA must account for interactions between a dipole and all its periodic images.

## The Convergence Problem

The naive lattice sum of the dyadic Green's function over periodic images:

\\[
\mathbf{G}_{\text{lat}}(\mathbf{r}, \mathbf{r}') = \sum_{\mathbf{R}_n} \mathbf{G}(\mathbf{r}, \mathbf{r}' + \mathbf{R}_n)
\\]

is **conditionally convergent** — it converges extremely slowly (if at all) in its original form. Direct summation is computationally intractable.

## Ewald Summation

Ewald's method resolves this by splitting the sum into two rapidly convergent parts:

1. **Real-space sum**: Short-range interactions computed in direct space, damped by a Gaussian screening function \\(e^{-\eta^2 R^2}\\).
2. **Reciprocal-space sum**: Long-range interactions computed in Fourier space, which converge exponentially.

The splitting parameter \\(\eta\\) controls the balance between the two sums. Both sums converge exponentially, making the total computation \\(O(N)\\) per unit cell.

## Physical Nuance

It is essential to understand that Ewald summation does not approximate the physics — it is an exact mathematical rearrangement of the lattice sum. The result is identical to the infinite direct sum; only the computational strategy differs.

The method naturally handles:
- **1D periodicity** (nanowires, gratings)
- **2D periodicity** (slabs, metasurfaces)
- **3D periodicity** (bulk photonic crystals)

Oblique incidence is accommodated via the Bloch phase factor \\(e^{i\mathbf{k}_\parallel \cdot \mathbf{R}_n}\\).
