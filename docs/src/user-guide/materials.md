# Materials

Lumina provides frequency-dependent material properties through the `MaterialProvider` trait.

## Built-in Materials

| Identifier | Description | Wavelength Range | Source |
|-----------|-------------|-----------------|--------|
| `Au_JC` | Gold | 400–800 nm | Johnson & Christy (1972) |
| `Ag_JC` | Silver | 400–800 nm | Johnson & Christy (1972) |
| `Cu_JC` | Copper | 400–800 nm | Johnson & Christy (1972) |

## Custom Materials

Specify a constant refractive index in the TOML configuration or provide a CSV file with wavelength-dependent data.

## Interpolation

Tabulated data is interpolated using natural cubic splines to provide smooth, continuous dielectric functions at arbitrary wavelengths within the data range. Extrapolation beyond the data range is not supported and will return an error.
