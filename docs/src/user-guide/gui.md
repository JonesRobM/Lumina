# GUI Dashboard

The Lumina GUI provides an interactive environment for configuring, running, and visualising simulations.

## Launching

```bash
cargo run --release -p lumina-gui
```

## Panels

### Geometry
Configure shapes, import files, adjust dipole spacing, and preview the dipole lattice in the 3D viewport.

### Materials
Select from built-in materials or define custom dielectric functions. View the interpolated \\(\epsilon(\lambda)\\) plot.

### Simulation
Set the wavelength range, medium refractive index, and solver parameters. Launch and monitor simulation progress.

### Results
View extinction/absorption/scattering spectra, near-field intensity maps, and export raw data.
