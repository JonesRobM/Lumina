# Configuration (TOML)

Lumina uses TOML files for simulation configuration. This page documents all available options.

## Sections

### `[simulation]`

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `wavelengths` | Range or list | Required | Wavelength specification |
| `environment_n` | float | `1.0` | Refractive index of surrounding medium |
| `solver_tolerance` | float | `1e-6` | Convergence tolerance for iterative solver |
| `max_iterations` | int | `1000` | Maximum GMRES iterations |

### `[[geometry.object]]`

Each object in the simulation is defined as a TOML array entry.

| Key | Type | Description |
|-----|------|-------------|
| `name` | string | Human-readable label |
| `material` | string | Material identifier (e.g. `"Au_JC"`) |
| `dipole_spacing` | float | Lattice spacing in nm |
| `type` | string | Shape type: `sphere`, `cylinder`, `cuboid`, `ellipsoid` |
| `geometry_file` | string | Path to `.xyz`, `.cif`, or `.obj` file (alternative to `type`) |

Shape-specific parameters are documented in the [Geometry](./geometry.md) section.

### `[output]`

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `directory` | string | `"./output"` | Output directory path |
| `save_spectra` | bool | `true` | Save cross-section spectra as CSV |
| `save_near_field` | bool | `false` | Compute and save near-field maps |
| `save_dipoles` | bool | `false` | Save raw dipole moments |
