# Materials

Lumina provides frequency-dependent material properties through the `MaterialProvider` trait.

## Built-in Materials

| Identifier | Description | Wavelength Range | Data Points | Source |
|-----------|-------------|-----------------|-------------|--------|
| `Au_JC` | Gold | 188–892 nm | 43 | Johnson & Christy (1972) |
| `Ag_JC` | Silver | 400–800 nm | — | Johnson & Christy (1972) |
| `Cu_JC` | Copper | 400–800 nm | — | Johnson & Christy (1972) |

### Johnson & Christy Gold Data

The gold dataset contains the full 43-point tabulation from Johnson & Christy (1972), as reproduced in the widely-used [refractiveindex.info](https://refractiveindex.info) database. The original \\((n, k)\\) data is converted to the complex dielectric function on import:

\\[
\epsilon_1 = n^2 - k^2, \qquad \epsilon_2 = 2nk
\\]

Key spectral regions for gold:

| Region | \\(\lambda\\) (nm) | Character | \\(|\epsilon_1|/\epsilon_2\\) |
|--------|---------------------|-----------|-------------------------------|
| UV | 188–300 | Interband transitions | \\(< 1\\) |
| Visible (interband) | 400–520 | Large \\(\epsilon_2\\) | \\(0.3\text{–}2\\) |
| Visible (Drude) | 550–892 | \\(\|\epsilon_1\| \gg \epsilon_2\\) | \\(4\text{–}30\\) |

The CDA is most accurate in the interband region where \\(|\epsilon_1|/\epsilon_2\\) is moderate. See [v0.1 Validation](../theory/validation.md) for accuracy benchmarks.

## Custom Materials

Specify a constant refractive index in the TOML configuration or provide a CSV file with wavelength-dependent data.

## Interpolation

Tabulated data is interpolated using natural cubic splines to provide smooth, continuous dielectric functions at arbitrary wavelengths within the data range. Extrapolation beyond the data range is not supported and will return an error.

## Reference

Johnson, P. B. & Christy, R. W., "Optical Constants of the Noble Metals", *Phys. Rev. B* **6**, 4370 (1972).
