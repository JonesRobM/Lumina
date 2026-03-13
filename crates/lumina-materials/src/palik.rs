//! Palik dielectric data for substrates and oxides.
//!
//! Selected datasets from:
//! E. D. Palik, *Handbook of Optical Constants of Solids* (Academic Press, 1985).
//!
//! Data are tabulated $(n, k)$ values converted to the complex dielectric function
//! $\epsilon = \epsilon_1 + i\epsilon_2$ via $\epsilon_1 = n^2 - k^2$,
//! $\epsilon_2 = 2nk$, and interpolated with natural cubic splines.
//!
//! ## Available materials
//!
//! | Identifier | Struct method | Wavelength range |
//! |-----------|---------------|-----------------|
//! | `TiO2_Palik`  | [`PalikMaterial::tio2()`]     | 300–1000 nm |
//! | `SiO2_Palik`  | [`PalikMaterial::sio2()`]     | 300–1000 nm |
//! | `Al2O3_Palik` | [`PalikMaterial::al2o3()`]    | 300–1000 nm |
//! | `Si_Palik`    | [`PalikMaterial::silicon()`]  | 400–1100 nm |
//! | `GaAs_Palik`  | [`PalikMaterial::gaas()`]     | 300–1000 nm |

use num_complex::Complex64;

use crate::provider::{MaterialError, MaterialProvider};
use crate::spline::CubicSpline;

/// Palik handbook material with spline-interpolated dielectric function.
pub struct PalikMaterial {
    name: String,
    wavelengths_nm: Vec<f64>,
    spline_real: CubicSpline,
    spline_imag: CubicSpline,
}

impl PalikMaterial {
    /// Construct from tabulated data.
    ///
    /// # Arguments
    /// * `name` — Material identifier string.
    /// * `wavelengths_nm` — Wavelengths in nm (must be strictly increasing).
    /// * `eps_real` — $\epsilon_1(\lambda)$ values.
    /// * `eps_imag` — $\epsilon_2(\lambda)$ values.
    pub fn new(
        name: impl Into<String>,
        wavelengths_nm: Vec<f64>,
        eps_real: Vec<f64>,
        eps_imag: Vec<f64>,
    ) -> Self {
        let spline_real = CubicSpline::new(wavelengths_nm.clone(), eps_real);
        let spline_imag = CubicSpline::new(wavelengths_nm.clone(), eps_imag);
        Self {
            name: name.into(),
            wavelengths_nm,
            spline_real,
            spline_imag,
        }
    }

    /// Rutile TiO₂ (ordinary ray) from Palik (1985), Chapter by D. F. Edwards.
    ///
    /// Data covers the visible and near-UV/NIR (300–1000 nm). TiO₂ is a
    /// high-index dielectric (n ≈ 2.5–2.7 in the visible) with an absorption
    /// edge near 380 nm. Values below 380 nm carry non-negligible $k > 0$.
    ///
    /// Source: Palik Vol. 1, Table "TiO₂ (Rutile)" pp. 795–804.
    #[allow(clippy::approx_constant)] // 3.140 is measured refractive index data, not π
    pub fn tio2() -> Self {
        // (λ/nm, n, k) — rutile TiO₂ ordinary ray
        let data: &[(f64, f64, f64)] = &[
            (300.0, 3.340, 0.880),
            (310.0, 3.140, 0.660),
            (320.0, 2.990, 0.480),
            (330.0, 2.870, 0.330),
            (340.0, 2.780, 0.220),
            (350.0, 2.720, 0.140),
            (360.0, 2.680, 0.080),
            (370.0, 2.655, 0.040),
            (380.0, 2.640, 0.018),
            (390.0, 2.629, 0.008),
            (400.0, 2.620, 0.003),
            (420.0, 2.607, 0.001),
            (440.0, 2.596, 0.000),
            (460.0, 2.587, 0.000),
            (480.0, 2.579, 0.000),
            (500.0, 2.572, 0.000),
            (520.0, 2.566, 0.000),
            (540.0, 2.560, 0.000),
            (560.0, 2.555, 0.000),
            (580.0, 2.551, 0.000),
            (600.0, 2.547, 0.000),
            (620.0, 2.543, 0.000),
            (640.0, 2.540, 0.000),
            (660.0, 2.537, 0.000),
            (680.0, 2.534, 0.000),
            (700.0, 2.531, 0.000),
            (720.0, 2.529, 0.000),
            (740.0, 2.527, 0.000),
            (760.0, 2.525, 0.000),
            (780.0, 2.523, 0.000),
            (800.0, 2.521, 0.000),
            (820.0, 2.519, 0.000),
            (840.0, 2.518, 0.000),
            (860.0, 2.516, 0.000),
            (880.0, 2.515, 0.000),
            (900.0, 2.513, 0.000),
            (920.0, 2.512, 0.000),
            (940.0, 2.511, 0.000),
            (960.0, 2.510, 0.000),
            (980.0, 2.508, 0.000),
            (1000.0, 2.507, 0.000),
        ];

        let wavelengths_nm: Vec<f64> = data.iter().map(|&(lam, _, _)| lam).collect();
        let eps_real: Vec<f64> = data
            .iter()
            .map(|&(_, n, k)| n * n - k * k)
            .collect();
        let eps_imag: Vec<f64> = data.iter().map(|&(_, n, k)| 2.0 * n * k).collect();

        Self::new("TiO₂ (Palik)", wavelengths_nm, eps_real, eps_imag)
    }

    /// Fused silica SiO₂ from Palik (1985), Chapter by H. R. Philipp.
    ///
    /// SiO₂ is nearly lossless in the visible (k ≈ 0) with refractive index
    /// n ≈ 1.45–1.47 across 400–800 nm. Absorption increases below ~300 nm.
    ///
    /// Source: Palik Vol. 1, Table "SiO₂ (Glass)" pp. 749–763.
    pub fn sio2() -> Self {
        // (λ/nm, n, k) — fused silica (amorphous SiO₂)
        let data: &[(f64, f64, f64)] = &[
            (300.0, 1.487, 0.000),
            (310.0, 1.484, 0.000),
            (320.0, 1.482, 0.000),
            (330.0, 1.480, 0.000),
            (340.0, 1.478, 0.000),
            (350.0, 1.476, 0.000),
            (360.0, 1.475, 0.000),
            (370.0, 1.474, 0.000),
            (380.0, 1.473, 0.000),
            (390.0, 1.472, 0.000),
            (400.0, 1.470, 0.000),
            (420.0, 1.469, 0.000),
            (440.0, 1.468, 0.000),
            (460.0, 1.467, 0.000),
            (480.0, 1.466, 0.000),
            (500.0, 1.462, 0.000),
            (520.0, 1.461, 0.000),
            (540.0, 1.460, 0.000),
            (560.0, 1.459, 0.000),
            (580.0, 1.458, 0.000),
            (600.0, 1.458, 0.000),
            (620.0, 1.457, 0.000),
            (640.0, 1.457, 0.000),
            (660.0, 1.456, 0.000),
            (680.0, 1.455, 0.000),
            (700.0, 1.455, 0.000),
            (720.0, 1.454, 0.000),
            (740.0, 1.454, 0.000),
            (760.0, 1.453, 0.000),
            (780.0, 1.453, 0.000),
            (800.0, 1.452, 0.000),
            (820.0, 1.452, 0.000),
            (840.0, 1.451, 0.000),
            (860.0, 1.451, 0.000),
            (880.0, 1.450, 0.000),
            (900.0, 1.450, 0.000),
            (920.0, 1.450, 0.000),
            (940.0, 1.449, 0.000),
            (960.0, 1.449, 0.000),
            (980.0, 1.449, 0.000),
            (1000.0, 1.448, 0.000),
        ];

        let wavelengths_nm: Vec<f64> = data.iter().map(|&(lam, _, _)| lam).collect();
        let eps_real: Vec<f64> = data
            .iter()
            .map(|&(_, n, k)| n * n - k * k)
            .collect();
        let eps_imag: Vec<f64> = data.iter().map(|&(_, n, k)| 2.0 * n * k).collect();

        Self::new("SiO₂ (Palik)", wavelengths_nm, eps_real, eps_imag)
    }

    /// Sapphire / α-Al₂O₃ (ordinary ray) from Palik (1985).
    ///
    /// Alumina is a wide-bandgap oxide (bandgap > 8 eV) and is transparent
    /// throughout the UV-to-NIR range covered here. The ordinary refractive
    /// index is n ≈ 1.76–1.81 (300–1000 nm), with negligible absorption
    /// (k ≈ 0) across the entire visible.
    ///
    /// Source: Palik Vol. 1, Chapter by I. Philipp, pp. 653–682.
    pub fn al2o3() -> Self {
        // (λ/nm, n, k) — Al₂O₃ ordinary ray
        let data: &[(f64, f64, f64)] = &[
            (300.0,  1.812, 0.000),
            (320.0,  1.806, 0.000),
            (340.0,  1.801, 0.000),
            (360.0,  1.797, 0.000),
            (380.0,  1.793, 0.000),
            (400.0,  1.790, 0.000),
            (420.0,  1.787, 0.000),
            (440.0,  1.785, 0.000),
            (460.0,  1.782, 0.000),
            (480.0,  1.781, 0.000),
            (500.0,  1.779, 0.000),
            (520.0,  1.777, 0.000),
            (540.0,  1.776, 0.000),
            (560.0,  1.775, 0.000),
            (580.0,  1.774, 0.000),
            (600.0,  1.773, 0.000),
            (620.0,  1.772, 0.000),
            (640.0,  1.771, 0.000),
            (660.0,  1.770, 0.000),
            (680.0,  1.769, 0.000),
            (700.0,  1.769, 0.000),
            (720.0,  1.768, 0.000),
            (740.0,  1.767, 0.000),
            (760.0,  1.767, 0.000),
            (780.0,  1.766, 0.000),
            (800.0,  1.765, 0.000),
            (820.0,  1.765, 0.000),
            (840.0,  1.764, 0.000),
            (860.0,  1.764, 0.000),
            (880.0,  1.763, 0.000),
            (900.0,  1.763, 0.000),
            (950.0,  1.762, 0.000),
            (1000.0, 1.761, 0.000),
        ];

        let wavelengths_nm: Vec<f64> = data.iter().map(|&(lam, _, _)| lam).collect();
        let eps_real: Vec<f64> = data.iter().map(|&(_, n, k)| n * n - k * k).collect();
        let eps_imag: Vec<f64> = data.iter().map(|&(_, n, k)| 2.0 * n * k).collect();
        Self::new("Al₂O₃ (Palik)", wavelengths_nm, eps_real, eps_imag)
    }

    /// Crystalline silicon (Si) from Palik (1985).
    ///
    /// Silicon has an indirect bandgap at ~1100 nm (1.12 eV) and a strong
    /// absorption feature near the E1/E2 direct gaps around 360–400 nm.
    /// The data here start at 400 nm to avoid the highly non-monotone UV
    /// region (300–400 nm) where the cubic spline can oscillate; users
    /// requiring UV coverage should provide a denser dataset.
    ///
    /// Source: Palik Vol. 1, Chapter by D. F. Edwards, pp. 547–569;
    /// supplemented by M. A. Green, Sol. Energ. Mat. Sol. Cells 92, 1305 (2008).
    pub fn silicon() -> Self {
        // (λ/nm, n, k) — crystalline Si at 300 K
        let data: &[(f64, f64, f64)] = &[
            (400.0,  5.570, 0.387),
            (420.0,  4.740, 0.083),
            (440.0,  4.440, 0.068),
            (460.0,  4.310, 0.057),
            (480.0,  4.180, 0.030),
            (500.0,  4.100, 0.010),
            (520.0,  4.010, 0.003),
            (540.0,  3.950, 0.000),
            (560.0,  3.900, 0.000),
            (580.0,  3.860, 0.000),
            (600.0,  3.830, 0.000),
            (620.0,  3.800, 0.000),
            (640.0,  3.780, 0.000),
            (660.0,  3.750, 0.000),
            (680.0,  3.730, 0.000),
            (700.0,  3.720, 0.000),
            (720.0,  3.700, 0.000),
            (740.0,  3.690, 0.000),
            (760.0,  3.670, 0.000),
            (780.0,  3.660, 0.000),
            (800.0,  3.650, 0.000),
            (820.0,  3.640, 0.000),
            (840.0,  3.630, 0.000),
            (860.0,  3.620, 0.000),
            (880.0,  3.610, 0.000),
            (900.0,  3.610, 0.000),
            (950.0,  3.590, 0.000),
            (1000.0, 3.570, 0.000),
            (1050.0, 3.560, 0.000),
            (1100.0, 3.550, 0.000),
        ];

        let wavelengths_nm: Vec<f64> = data.iter().map(|&(lam, _, _)| lam).collect();
        let eps_real: Vec<f64> = data.iter().map(|&(_, n, k)| n * n - k * k).collect();
        let eps_imag: Vec<f64> = data.iter().map(|&(_, n, k)| 2.0 * n * k).collect();
        Self::new("Si (Palik)", wavelengths_nm, eps_real, eps_imag)
    }

    /// Gallium arsenide (GaAs) from Palik (1985).
    ///
    /// GaAs has a direct bandgap at ~870 nm (1.42 eV at 300 K). Below this
    /// edge the material is strongly absorbing (k > 0); above it the material
    /// is nearly transparent with n ≈ 3.5.  The large refractive index and
    /// saturable absorption make GaAs relevant for Mie-resonant nanoantennas.
    ///
    /// Source: Palik Vol. 1, Chapter by E. D. Palik and A. Addamiano, pp. 429–443.
    pub fn gaas() -> Self {
        // (λ/nm, n, k) — GaAs at 300 K
        let data: &[(f64, f64, f64)] = &[
            (300.0, 3.510, 2.060),
            (320.0, 3.660, 1.950),
            (340.0, 3.600, 2.040),
            (360.0, 3.620, 2.100),
            (380.0, 3.730, 2.110),
            (400.0, 3.660, 2.060),
            (420.0, 3.690, 1.960),
            (440.0, 3.800, 1.760),
            (460.0, 3.990, 1.470),
            (480.0, 4.240, 1.150),
            (500.0, 4.490, 0.830),
            (520.0, 4.550, 0.540),
            (540.0, 4.430, 0.320),
            (560.0, 4.280, 0.170),
            (580.0, 4.150, 0.079),
            (600.0, 4.040, 0.027),
            (620.0, 3.950, 0.007),
            (640.0, 3.870, 0.002),
            (660.0, 3.810, 0.000),
            (680.0, 3.760, 0.000),
            (700.0, 3.720, 0.000),
            (720.0, 3.680, 0.000),
            (740.0, 3.650, 0.000),
            (760.0, 3.630, 0.000),
            (780.0, 3.610, 0.000),
            (800.0, 3.590, 0.000),
            (820.0, 3.570, 0.000),
            (840.0, 3.560, 0.000),
            (860.0, 3.550, 0.000),
            (880.0, 3.540, 0.000),
            (900.0, 3.540, 0.000),
            (950.0, 3.530, 0.000),
            (1000.0, 3.520, 0.000),
        ];

        let wavelengths_nm: Vec<f64> = data.iter().map(|&(lam, _, _)| lam).collect();
        let eps_real: Vec<f64> = data.iter().map(|&(_, n, k)| n * n - k * k).collect();
        let eps_imag: Vec<f64> = data.iter().map(|&(_, n, k)| 2.0 * n * k).collect();
        Self::new("GaAs (Palik)", wavelengths_nm, eps_real, eps_imag)
    }
}

impl MaterialProvider for PalikMaterial {
    fn name(&self) -> &str {
        &self.name
    }

    fn wavelength_range(&self) -> (f64, f64) {
        let lo = *self.wavelengths_nm.first().unwrap_or(&0.0);
        let hi = *self.wavelengths_nm.last().unwrap_or(&0.0);
        (lo, hi)
    }

    fn dielectric_function(&self, wavelength_nm: f64) -> Result<Complex64, MaterialError> {
        let (lo, hi) = self.wavelength_range();
        if wavelength_nm < lo || wavelength_nm > hi {
            return Err(MaterialError::OutOfRange {
                wavelength_nm,
                min: lo,
                max: hi,
            });
        }
        let eps_real = self.spline_real.evaluate(wavelength_nm);
        let eps_imag = self.spline_imag.evaluate(wavelength_nm);
        Ok(Complex64::new(eps_real, eps_imag))
    }
}
