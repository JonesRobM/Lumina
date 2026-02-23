# Materials

Lumina provides frequency-dependent material properties through the `MaterialProvider` trait.

## Built-in Materials

### Metals (Johnson & Christy)

| Identifier | Description | Wavelength Range | Data Points | Source |
|-----------|-------------|-----------------|-------------|--------|
| `Au_JC` | Gold | 188–892 nm | 43 | Johnson & Christy (1972) |
| `Ag_JC` | Silver | 188–892 nm | 43 | Johnson & Christy (1972) |
| `Cu_JC` | Copper | 188–892 nm | 43 | Johnson & Christy (1972) |

### Dielectrics (Palik)

| Identifier | Description | Wavelength Range | Source |
|-----------|-------------|-----------------|--------|
| `TiO2_Palik` | Titanium dioxide | 300–1000 nm | Palik (1998) |
| `SiO2_Palik` | Silicon dioxide | 300–1000 nm | Palik (1998) |

### Johnson & Christy Metal Data

The metal datasets contain the full 43-point tabulations from Johnson & Christy (1972), as reproduced in the widely-used [refractiveindex.info](https://refractiveindex.info) database. The original \\((n, k)\\) data is converted to the complex dielectric function on import:

\\[
\epsilon_1 = n^2 - k^2, \qquad \epsilon_2 = 2nk
\\]

Key spectral regions for gold:

| Region | \\(\lambda\\) (nm) | Character | \\(|\epsilon_1|/\epsilon_2\\) |
|--------|---------------------|-----------|-------------------------------|
| UV | 188–300 | Interband transitions | \\(< 1\\) |
| Visible (interband) | 400–520 | Large \\(\epsilon_2\\) | \\(0.3\text{–}2\\) |
| Visible (Drude) | 550–892 | \\(\|\epsilon_1\| \gg \epsilon_2\\) | \\(4\text{–}30\\) |

The CDA is most accurate in the interband region where \\(|\epsilon_1|/\epsilon_2\\) is moderate. See [Validation](../theory/validation.md) for accuracy benchmarks.

### Palik Dielectric Data

TiO\\(_2\\) and SiO\\(_2\\) data are taken from the Palik (1998) handbook. These are high-refractive-index and low-refractive-index dielectrics commonly used as substrates and coatings in nanophotonics.

- **TiO\\(_2\\):** \\(n \approx 2.4\text{–}2.7\\) across 300–1000 nm, with an absorption edge near 350 nm.
- **SiO\\(_2\\):** \\(n \approx 1.45\text{–}1.47\\) across 300–1000 nm, essentially lossless in the visible.

## Custom Materials

### Constant dielectric

Specify a constant complex dielectric function in the TOML configuration:

```toml
material = "Custom(4.0, 0.5)"   # ε = 4.0 + 0.5i
```

This is useful for benchmarking, convergence studies, and materials where the wavelength dependence is negligible over the simulated range.

### Tabulated CSV (v0.1.3+)

Provide a CSV file with wavelength-dependent data:

```toml
material = "CSV:path/to/epsilon.csv"
```

The CSV should have columns: `wavelength_nm`, `epsilon_real`, `epsilon_imag`.

## Interpolation

Tabulated data is interpolated using natural cubic splines to provide smooth, continuous dielectric functions at arbitrary wavelengths within the data range. Extrapolation beyond the data range is not supported and will return an error.

## Viewing \\(\epsilon(\lambda)\\) Curves

The GUI materials panel displays the real and imaginary parts of the dielectric function as interactive plots. Select a material from the dropdown to see its \\(\epsilon_1(\lambda)\\) and \\(\epsilon_2(\lambda)\\) curves across the full data range.

## References

- Johnson, P. B. & Christy, R. W., "Optical Constants of the Noble Metals", *Phys. Rev. B* **6**, 4370 (1972).
- Palik, E. D., *Handbook of Optical Constants of Solids*, Academic Press (1998).
