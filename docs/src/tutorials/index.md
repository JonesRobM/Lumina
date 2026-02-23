# Tutorials

These step-by-step tutorials take you from a blank configuration file to a fully interpreted result. Each tutorial is self-contained; you can work through them in order or jump to whichever topic interests you.

## Prerequisites

- Lumina built in release mode (`cargo build --release`)
- Basic familiarity with the command line
- (Optional) GPU build for Tutorial 6: `cargo build --release --features gpu`
- (Optional) Python + matplotlib for post-processing plots

## Tutorial List

| # | Tutorial | What you learn |
|---|----------|---------------|
| 1 | [Gold Nanosphere Spectrum](./01-gold-nanosphere.md) | Basic CLI workflow, extinction/absorption/scattering |
| 2 | [Nanorod & Plasmon Red-Shift](./02-nanorod.md) | Cylinder geometry, aspect ratio dependence |
| 3 | [Near-Field Intensity Map](./03-near-field.md) | Near-field output, \|E\|² interpretation |
| 4 | [Chiral Helix](./04-chiral-helix.md) | Helix geometry, multi-object configurations |
| 5 | [Custom Dielectric Material](./05-custom-material.md) | Constant-ε material, convergence study |
| 6 | [GPU-Accelerated Large System](./06-gpu-acceleration.md) | GPU backend, large dipole counts, benchmarking |

## Conventions Used Throughout

- All lengths are in **nanometres (nm)** unless stated otherwise.
- Cross-sections are in **nm²**.
- Wavelength always refers to the free-space wavelength λ₀ = 2π/k.
- Example TOML files are also available in the `examples/` directory at the repository root.
