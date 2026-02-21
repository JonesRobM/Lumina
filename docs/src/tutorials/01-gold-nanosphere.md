# Tutorial 1: Gold Nanosphere Spectrum

**Goal:** Compute the extinction, absorption, and scattering cross-sections of a 20 nm gold nanosphere in vacuum across 400–800 nm, then compare the result to analytical Mie theory.

**Time:** ~5 minutes

---

## Background

The gold nanosphere is the canonical test case for any CDA implementation. Gold nanoparticles exhibit a localised surface plasmon resonance (LSPR) in the visible region. For a 20 nm sphere the LSPR peak sits near 520 nm. Mie theory provides an analytical solution — Lumina's validation suite checks against this.

The key physics:

- Johnson & Christy (1972) dielectric data, interpolated via cubic spline
- Clausius-Mossotti polarisability with radiative reaction correction (RRCM)
- Filtered Coupled Dipole (FCD) for improved accuracy near-surface

---

## Step 1: Write the Configuration File

Create `examples/gold_sphere_tutorial.toml`:

```toml
# Tutorial 1 — 20 nm gold sphere in vacuum
[simulation]
wavelengths    = { range = [400.0, 800.0], points = 80 }
environment_n  = 1.0          # vacuum
solver_tolerance  = 1e-6
max_iterations    = 500

[[geometry.object]]
name           = "Au_Sphere"
type           = "sphere"
centre         = [0.0, 0.0, 0.0]
radius         = 20.0         # nm
material       = "Au_JC"      # Johnson & Christy gold
dipole_spacing = 2.0          # nm  → ~4189 dipoles

[output]
directory      = "./output/tutorial_01"
save_spectra   = true
save_near_field = false
save_dipoles   = false
```

> **Dipole count:** With `radius = 20.0` and `dipole_spacing = 2.0`, Lumina discretises the sphere into approximately 4 000–4 500 dipoles. This triggers the GMRES solver automatically (threshold: N > 1 000).

---

## Step 2: Run the Simulation

```bash
cargo run --release --bin lumina-cli -- examples/gold_sphere_tutorial.toml
```

Expected terminal output:

```
[Lumina] Geometry: Au_Sphere — 4189 dipoles (sphere R=20nm, d=2nm)
[Lumina] Solver: GMRES(m=30), tolerance=1e-6
[Lumina] Computing 80 wavelengths…
[Lumina] λ=400 nm … C_ext=1.42e2 nm²   C_abs=1.24e2 nm²   C_sca=1.78e1 nm²
…
[Lumina] λ=520 nm … C_ext=9.83e2 nm²   C_abs=7.61e2 nm²   C_sca=2.22e2 nm²   ← resonance
…
[Lumina] Done. Results written to ./output/tutorial_01/spectra.csv
```

---

## Step 3: Inspect the Output

The CSV file has the following structure:

```
wavelength_nm,extinction_nm2,absorption_nm2,scattering_nm2
400.0,142.3,124.1,18.2
405.1,158.9,137.4,21.5
…
```

Quick check with Python:

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("output/tutorial_01/spectra.csv")

plt.figure(figsize=(8, 5))
plt.plot(df["wavelength_nm"], df["extinction_nm2"],  label="Extinction",  color="black")
plt.plot(df["wavelength_nm"], df["absorption_nm2"],  label="Absorption",  color="red",   linestyle="--")
plt.plot(df["wavelength_nm"], df["scattering_nm2"],  label="Scattering",  color="blue",  linestyle=":")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Cross-Section (nm²)")
plt.title("Au sphere, R=20 nm, vacuum")
plt.legend()
plt.tight_layout()
plt.savefig("gold_sphere.png", dpi=150)
```

---

## Step 4: Understand the Result

You should observe:

| Wavelength | Feature |
|-----------|---------|
| 420–510 nm | Interband absorption (d→sp transitions in Au) |
| ~520 nm | LSPR peak in extinction and scattering |
| >550 nm | Drude tail — scattering falls off sharply |

**Conservation check:** `C_ext = C_abs + C_sca` should hold to numerical precision (< 0.1% deviation).

```python
df["check"] = (df["absorption_nm2"] + df["scattering_nm2"]) / df["extinction_nm2"]
print(df["check"].describe())   # mean ≈ 1.000, std < 0.001
```

---

## Step 5: Compare to Mie Theory (optional)

Lumina includes a Mie theory calculator in its test suite. Run it with `--nocapture` to print the comparison:

```bash
cargo test --release -p lumina-core --test mie_validation -- --nocapture
```

The output table shows CDA error vs Mie at each wavelength. For R=10 nm at d=2 nm expect:

| Region | Typical error |
|--------|--------------|
| 420–510 nm (interband) | < 25% |
| 520–550 nm (near-LSPR) | < 30% |
| > 550 nm (Drude) | 80–300% |

The Drude region error is a known limitation of the CDA on coarse grids — see [Validation](../theory/validation.md) and the v0.2.0 roadmap for surface-averaging improvements.

---

## GUI Alternative

All of the above can be done interactively:

```bash
cargo run --release --bin lumina-gui
```

1. **Geometry panel** → Shape: Sphere, Radius: 20, Spacing: 2
2. **Material panel** → Au (Johnson & Christy)
3. **Simulation panel** → 400–800 nm, 80 points, environment n=1.0
4. **Run** → spectra appear in the Results panel with zoom/pan

---

## Next Steps

- [Tutorial 2: Nanorod & Plasmon Red-Shift](./02-nanorod.md) — change the geometry to a cylinder and observe the longitudinal plasmon
- [Tutorial 3: Near-Field Map](./03-near-field.md) — visualise |E|² around the resonance
