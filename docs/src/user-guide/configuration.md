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
| `backend` | string | `"auto"` | Compute backend: `"auto"`, `"cpu"`, or `"gpu"` |
| `k_bloch` | `[float; 3]` | `[0,0,0]` | Bloch wavevector in nm⁻¹ (for periodic structures) |

The `backend` field selects the compute backend for GMRES matrix-vector products:

- `"auto"` (default) — uses GPU if available, falls back to CPU
- `"cpu"` — always use CPU (Rayon)
- `"gpu"` — require GPU (wgpu); errors if no GPU is available

> **Note:** GPU acceleration requires building with `--features gpu`. When built without this feature, the `backend` field is ignored and CPU is always used.

### `[[geometry.object]]`

Each object in the simulation is defined as a TOML array entry. Objects are merged into a single dipole list in definition order.

| Key | Type | Description |
|-----|------|-------------|
| `name` | string | Human-readable label |
| `material` | string | Material identifier (e.g. `"Au_JC"`) — not used for `coreshell` |
| `dipole_spacing` | float | Lattice spacing in nm |
| `type` | string | Shape type: `sphere`, `cylinder`, `cuboid`, `ellipsoid`, `helix`, `coreshell`, `file` |

**Transform** (optional, all fields default to identity):

```toml
[geometry.object.transform]
position    = [0.0, 0.0, 0.0]   # nm
rotation_deg = [0.0, 0.0, 0.0]  # Euler XYZ in degrees
scale       = 1.0
```

**CoreShell** fields (when `type = "coreshell"`):

| Key | Type | Description |
|-----|------|-------------|
| `core_material` | string | Material of the innermost solid region |
| `[core_shape]` | table | Shape of the core (same keys as a primitive) |
| `[[shells]]` | array | Ordered list of `{ thickness, material }` entries |

**File import** (when `type = "file"`):

| Key | Type | Description |
|-----|------|-------------|
| `path` | string | Path to `.xyz` or `.obj` file |
| `species_map` | table | Maps element symbols to material strings (`.xyz` only) |

Shape-specific parameters are documented in the [Geometry](./geometry.md) section.

### `[lattice]` (v0.3.0+)

Enables periodic boundary conditions. Omit this section for isolated (non-periodic) structures.

```toml
# 1D chain, period 50 nm along x
[lattice]
type = "chain"
a1 = [50.0, 0.0, 0.0]

# 2D square lattice
[lattice]
type = "planar"
a1 = [50.0, 0.0, 0.0]
a2 = [0.0,  50.0, 0.0]
```

See [Periodic Structures](../theory/periodic.md) for the underlying theory.

### `[simulation.substrate]` (v0.3.0+)

Models a flat dielectric substrate below the structure using an image-dipole approximation. Omit this section for a free-standing structure.

```toml
[simulation.substrate]
material     = "SiO2_Palik"   # substrate material identifier
z_interface  = -25.0          # z-coordinate of the interface in nm (structure is above)
use_retarded = true           # optional — default true (v0.4.0+)
```

The substrate is also available in the GUI via the **Substrate** checkbox in the Simulation panel.

Each real dipole **p**_j induces an image dipole at its mirror position below the interface. The image contribution is scaled by the reflection factor:

- `use_retarded = true` (default, v0.4.0+): retarded normal-incidence Fresnel amplitude r = (n_sub − n_env)/(n_sub + n_env), where n = √ε is the complex refractive index.
- `use_retarded = false` (quasi-static, v0.3.0 behaviour): Δε = (ε_sub − ε_env)/(ε_sub + ε_env).

For a glass substrate (n = 1.5) in air, these give r = 0.20 vs Δε = 0.385 — a factor of ~2. The retarded form is physically correct for optical-frequency calculations.

All dipoles must satisfy z > z_interface.

### `[nonlinear]` (v0.4.0+)

Enables second-harmonic generation (SHG) computation. Omit this section for linear-only simulations.

```toml
[nonlinear]
enable_shg = true
symmetry   = "isotropic_surface"   # "isotropic_surface" | "zero"
chi_zzz    = [1.0, 0.0]           # χ_zzz [real, imag] in nm³
chi_zxx    = [0.3, 0.0]           # χ_zxx [real, imag] in nm³
far_field  = false                 # compute far-field pattern at 2ω
```

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `enable_shg` | bool | `false` | Enable SHG calculation |
| `symmetry` | string | `"zero"` | χ^(2) symmetry class: `"isotropic_surface"` or `"zero"` |
| `chi_zzz` | `[float; 2]` | — | χ_zzz = [real, imag] in nm³ (required for isotropic_surface) |
| `chi_zxx` | `[float; 2]` | — | χ_zxx = [real, imag] in nm³ (required for isotropic_surface) |
| `far_field` | bool | `false` | Compute far-field radiation pattern at 2ω |
| `enable_thg` | bool | `false` | Enable THG calculation |
| `chi3_symmetry` | string | `"zero"` | χ^(3) symmetry class: `"isotropic_bulk"` or `"zero"` |
| `chi3_xxxx` | `[float; 2]` | — | χ_xxxx = [real, imag] in nm⁶ (required for isotropic_bulk) |
| `chi3_xxyy` | `[float; 2]` | — | χ_xxyy = [real, imag] in nm⁶ (required for isotropic_bulk) |

Outputs: `output/shg_spectrum.csv` (`fundamental_nm`, `harmonic_nm`, `shg_intensity_nm6`) and/or `output/thg_spectrum.csv` (`fundamental_nm`, `harmonic_nm`, `thg_intensity_nm9`).

See [Nonlinear Optics](../theory/nonlinear.md) for the underlying theory.

### `[output]`

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `directory` | string | `"./output"` | Output directory path |
| `save_spectra` | bool | `true` | Save cross-section spectra as CSV |
| `save_json` | bool | `false` | Also save spectra as JSON |
| `save_near_field` | bool | `false` | Compute and save near-field maps at peak extinction |
