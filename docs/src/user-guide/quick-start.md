# Quick Start

This guide walks through computing the extinction spectrum of a gold nanosphere.

## 1. Create a Configuration File

Save the following as `my_sphere.toml`:

```toml
[simulation]
wavelengths = { range = [400.0, 900.0], points = 100 }
environment_n = 1.0

[[geometry.object]]
name = "Gold_Sphere"
type = "sphere"
centre = [0.0, 0.0, 0.0]
radius = 20.0
material = "Au_JC"
dipole_spacing = 2.0

[output]
directory = "./output"
save_spectra = true
```

## 2. Run via CLI

```bash
lumina-cli run my_sphere.toml
```

## 3. Or Use the GUI

Launch the GUI and configure the same parameters interactively:

```bash
cargo run -p lumina-gui
```

## 4. Examine the Results

The output directory will contain:
- `spectra.csv` — wavelength, extinction, absorption, and scattering cross-sections.
- `near_field/` — field intensity maps (if enabled).
