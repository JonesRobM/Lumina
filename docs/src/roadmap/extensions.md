# Future Extensions

This page logs planned extensions beyond the core CDA functionality.

## Alternative Solver Methods

The `OpticalSolver` trait is designed to accommodate fundamentally different electromagnetic simulation approaches:

### Boundary Element Method (BEM)
Solves Maxwell's equations on the surface of a structure rather than its volume. More efficient than CDA for large, smooth objects but requires surface meshing.

### T-Matrix Method
Expands fields in vector spherical harmonics. Highly efficient for clusters of spherical or spheroidal particles. Exact for spheres (reduces to Mie theory for a single sphere).

### Finite-Difference Time-Domain (FDTD)
Propagates E and B fields on a spatial grid in the time domain. Handles arbitrary geometries and nonlinear materials. May be implemented as a standalone project that shares the `OpticalSolver` interface.

## Material Classes Beyond Plasmonics

- **Ceramics** — High-refractive-index dielectrics (TiO2, Si3N4) for Mie resonances.
- **Composites** — Effective medium theories (Maxwell-Garnett, Bruggeman) for mixed materials.
- **2D Materials** — Surface conductivity models for graphene, MoS2, and other layered materials.
- **Active Media** — Gain materials for lasing and loss compensation studies.

## Advanced Physics

- **Magnetic dipoles** — Extend beyond electric dipoles to include magnetic dipole and electric quadrupole contributions.
- **Substrate effects** — Image dipole method or layered-medium Green's function for particles on or near surfaces.
- **Thermal effects** — Temperature-dependent material properties and photothermal heating.
- **Quantum corrections** — Nonlocal dielectric response and quantum tunnelling for sub-nanometre gaps.

## Workflow & Infrastructure

- **Python bindings** via PyO3 for integration with the scientific Python ecosystem.
- **Jupyter widget** for interactive notebook-based simulations.
- **Cloud deployment** — containerised solver for cloud HPC (AWS, Azure, GCP).
- **Database** — Materials and results database with search and comparison tools.
