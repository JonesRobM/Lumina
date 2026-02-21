# Tutorial 5: Custom Dielectric Material

**Goal:** Define a material with a user-specified complex dielectric constant ε = ε₁ + iε₂ and use it in a simulation. Then run a convergence study showing how the CDA result improves as the dipole spacing decreases.

**Time:** ~15 minutes

---

## Background

Lumina ships with Johnson & Christy gold, silver, and copper optical data. For any other material — a semiconductor, a custom glass, an experimentally measured thin film — you can supply either:

1. **A constant complex dielectric:** `"Custom(er, ei)"` where `er = ε₁` and `ei = ε₂`
2. **A tabulated CSV file:** `"CSV:path/to/epsilon.csv"` (v0.1.3+)

This tutorial uses the constant-ε option, which is ideal for benchmarking and sensitivity analysis.

---

## Step 1: Dielectric Sphere (reference case)

A dielectric sphere with ε = 4.0 + 0.0i (a TiO₂-like insulator) has an exact Mie theory solution and small CDA error (< 15%). This makes it an ideal convergence benchmark.

Create `examples/dielectric_sphere.toml`:

```toml
# Tutorial 5 — dielectric sphere, ε = 4 + 0i, convergence study
[simulation]
wavelengths   = { range = [400.0, 800.0], points = 50 }
environment_n = 1.0

[[geometry.object]]
name          = "Dielectric_Sphere"
type          = "sphere"
centre        = [0.0, 0.0, 0.0]
radius        = 20.0           # nm
material      = "Custom(4.0, 0.0)"   # ε₁=4.0, ε₂=0.0 (lossless)
dipole_spacing = 4.0           # nm — start coarse

[output]
directory     = "./output/tutorial_05/d4"
save_spectra  = true
```

Run:

```bash
cargo run --release --bin lumina-cli -- examples/dielectric_sphere.toml
```

---

## Step 2: Convergence Study

The CDA result should converge as the dipole spacing d → 0. Run the same geometry at four spacings:

```bash
for D in 4 3 2 1; do
  sed "s/dipole_spacing = 4.0/dipole_spacing = ${D}.0/" \
      examples/dielectric_sphere.toml \
    | sed "s|d4|d${D}|" \
    > /tmp/sphere_d${D}.toml
  cargo run --release --bin lumina-cli -- /tmp/sphere_d${D}.toml
done
```

Expected dipole counts (R=20 nm sphere):

| Spacing (nm) | N (approx) | Solver |
|-------------|-----------|--------|
| 4.0 | ~520 | Direct LU |
| 3.0 | ~1 250 | GMRES |
| 2.0 | ~4 100 | GMRES |
| 1.0 | ~33 500 | GMRES |

> **Warning:** d=1 nm at R=20 nm generates ~33 000 dipoles. GMRES with N=33k takes ~5–10 minutes on a modern CPU.

---

## Step 3: Compare to Mie Theory

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Mie extinction for a dielectric sphere (ε=4+0i, R=20nm, vacuum)
# Generated via the lumina-core mie_validation test suite
# You can also use miepython: pip install miepython
try:
    import miepython
    wl = np.linspace(400, 800, 50)
    m  = np.sqrt(4.0 + 0j)   # refractive index for ε=4
    x  = 2 * np.pi * 20.0 / wl  # size parameter
    qext, qsca, qback, g = miepython.mie(m, x)
    cext_mie = qext * np.pi * 20.0**2  # nm²
except ImportError:
    print("Install miepython for Mie reference: pip install miepython")
    cext_mie = None

spacings = [4, 3, 2, 1]
colours  = ["#F44336", "#FF9800", "#4CAF50", "#2196F3"]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: spectra for each dipole spacing
for d, colour in zip(spacings, colours):
    df = pd.read_csv(f"output/tutorial_05/d{d}/spectra.csv")
    axes[0].plot(df["wavelength_nm"], df["extinction_nm2"],
                 label=f"d = {d} nm", color=colour)

if cext_mie is not None:
    axes[0].plot(wl, cext_mie, "k--", linewidth=2, label="Mie (exact)")

axes[0].set_xlabel("Wavelength (nm)")
axes[0].set_ylabel("Extinction Cross-Section (nm²)")
axes[0].set_title("Convergence: dielectric sphere ε=4")
axes[0].legend()

# Right: relative error at λ=600 nm vs dipole spacing
if cext_mie is not None:
    ref_idx = np.argmin(np.abs(wl - 600))
    ref_val = cext_mie[ref_idx]
    errors = []
    for d in spacings:
        df = pd.read_csv(f"output/tutorial_05/d{d}/spectra.csv")
        cda_val = np.interp(600, df["wavelength_nm"], df["extinction_nm2"])
        errors.append(abs(cda_val - ref_val) / ref_val * 100)

    axes[1].loglog(spacings, errors, "ko-")
    axes[1].set_xlabel("Dipole spacing d (nm)")
    axes[1].set_ylabel("Relative error at λ=600 nm (%)")
    axes[1].set_title("CDA convergence rate")
    # Expected: roughly O(d²) for dielectrics
    d_range = np.array([0.8, 5.0])
    axes[1].loglog(d_range, 5 * d_range**2, "r--", label="O(d²)")
    axes[1].legend()
    axes[1].invert_xaxis()

plt.tight_layout()
plt.savefig("dielectric_convergence.png", dpi=150)
```

---

## Step 4: Lossy Dielectric

To add absorption, set ε₂ > 0:

```toml
material = "Custom(4.0, 0.5)"   # ε = 4.0 + 0.5i  (lossy dielectric)
```

With non-zero ε₂, absorption becomes non-trivial and `C_abs` is no longer zero. The Mie series still converges, so you can benchmark:

```python
m_lossy = np.sqrt(4.0 + 0.5j)
qext_l, qsca_l, qback_l, g_l = miepython.mie(m_lossy, x)
```

Expected result: CDA error < 16% across 500–700 nm (see [Validation table](../theory/validation.md)).

---

## Step 5: Using the Rust API Directly

For programmatic access to the custom material:

```rust
use lumina_materials::provider::MaterialProvider;
use lumina_materials::custom::ConstantMaterial;
use num_complex::Complex64;

// Lossless dielectric ε = 4 + 0i
let material = ConstantMaterial::new(Complex64::new(4.0, 0.0));

let epsilon = material.dielectric_function(600.0).unwrap();
// epsilon = 4.0 + 0.0i
```

Or via the convenience alias from `lumina_core::types`:

```rust
use num_complex::Complex64;

let epsilon = Complex64::new(4.0, 0.5);
let alpha_cm = clausius_mossotti(spacing.powi(3), epsilon, 1.0);
```

---

## Choosing a Dipole Spacing

As a rule of thumb:

| Material type | Recommended d/R | Notes |
|---------------|----------------|-------|
| Lossless dielectric | ≤ 0.15 | Error < 5% |
| Lossy dielectric | ≤ 0.12 | Error < 10% |
| Au interband (420–510 nm) | ≤ 0.10 | Error < 25% with FCD |
| Au Drude (> 550 nm) | ≤ 0.05 | Still 30–80% error; use BEM |

For a 20 nm sphere, d = 0.10 R = 2 nm is the practical minimum on a desktop CPU without GPU acceleration.

---

## Next Steps

- Read the [Validation](../theory/validation.md) chapter for a systematic accuracy analysis
- See the [Roadmap](../roadmap/versions.md) for planned v0.2.0 GPU and surface-averaging improvements
- Return to the [Tutorials index](./index.md) for a summary of all tutorials
