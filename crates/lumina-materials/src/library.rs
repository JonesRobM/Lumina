//! Lookup of material providers by the identifiers used in job files.
//!
//! Job files and `species_map` entries name materials with short identifiers
//! such as `"Au_JC"` or `"Pt_LD"`. [`MaterialLibrary`] builds every provider
//! once and resolves those identifiers.

use num_complex::Complex64;

use crate::johnson_christy::JohnsonChristyMaterial;
use crate::lorentz_drude::LorentzDrudeMaterial;
use crate::palik::PalikMaterial;
use crate::provider::{MaterialError, MaterialProvider};

/// All built-in materials, keyed by job-file identifier.
pub struct MaterialLibrary {
    entries: Vec<(&'static str, Box<dyn MaterialProvider>)>,
}

impl MaterialLibrary {
    /// Build every built-in material.
    pub fn standard() -> Self {
        let entries: Vec<(&'static str, Box<dyn MaterialProvider>)> = vec![
            ("Au_JC", Box::new(JohnsonChristyMaterial::gold())),
            ("Ag_JC", Box::new(JohnsonChristyMaterial::silver())),
            ("Cu_JC", Box::new(JohnsonChristyMaterial::copper())),
            ("TiO2_Palik", Box::new(PalikMaterial::tio2())),
            ("SiO2_Palik", Box::new(PalikMaterial::sio2())),
            ("Al2O3_Palik", Box::new(PalikMaterial::al2o3())),
            ("Si_Palik", Box::new(PalikMaterial::silicon())),
            ("GaAs_Palik", Box::new(PalikMaterial::gaas())),
            ("Ag_LD", Box::new(LorentzDrudeMaterial::silver())),
            ("Al_LD", Box::new(LorentzDrudeMaterial::aluminium())),
            ("Au_LD", Box::new(LorentzDrudeMaterial::gold())),
            ("Be_LD", Box::new(LorentzDrudeMaterial::beryllium())),
            ("Cr_LD", Box::new(LorentzDrudeMaterial::chromium())),
            ("Cu_LD", Box::new(LorentzDrudeMaterial::copper())),
            ("Ni_LD", Box::new(LorentzDrudeMaterial::nickel())),
            ("Pd_LD", Box::new(LorentzDrudeMaterial::palladium())),
            ("Pt_LD", Box::new(LorentzDrudeMaterial::platinum())),
            ("Ti_LD", Box::new(LorentzDrudeMaterial::titanium())),
            ("W_LD", Box::new(LorentzDrudeMaterial::tungsten())),
        ];
        Self { entries }
    }

    /// Iterate over `(identifier, provider)` pairs in display order.
    pub fn iter(&self) -> impl Iterator<Item = (&'static str, &dyn MaterialProvider)> + '_ {
        self.entries.iter().map(|(id, p)| (*id, p.as_ref()))
    }

    /// Look up a provider by identifier.
    pub fn get(&self, id: &str) -> Result<&dyn MaterialProvider, MaterialError> {
        self.iter()
            .find(|(name, _)| *name == id)
            .map(|(_, p)| p)
            .ok_or_else(|| {
                let valid: Vec<&str> = self.iter().map(|(name, _)| name).collect();
                MaterialError::NotFound(format!("'{id}'. Valid: {}", valid.join(", ")))
            })
    }

    /// Dielectric function of material `id` at `wavelength_nm`.
    pub fn dielectric_function(
        &self,
        id: &str,
        wavelength_nm: f64,
    ) -> Result<Complex64, MaterialError> {
        self.get(id)?.dielectric_function(wavelength_nm)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolves_every_identifier() {
        let lib = MaterialLibrary::standard();
        for (id, provider) in lib.iter() {
            let (lo, hi) = provider.wavelength_range();
            let mid = 0.5 * (lo + hi);
            assert!(lib.dielectric_function(id, mid).is_ok(), "{id} at {mid} nm");
        }
    }

    #[test]
    fn lorentz_drude_names_resolve() {
        let lib = MaterialLibrary::standard();
        for id in [
            "Ag_LD", "Al_LD", "Au_LD", "Be_LD", "Cr_LD", "Cu_LD", "Ni_LD", "Pd_LD", "Pt_LD",
            "Ti_LD", "W_LD",
        ] {
            assert!(lib.dielectric_function(id, 600.0).is_ok(), "{id}");
        }
    }

    #[test]
    fn unknown_identifier_lists_valid_names() {
        let lib = MaterialLibrary::standard();
        let err = lib
            .get("Pt")
            .err()
            .expect("unknown id must fail")
            .to_string();
        assert!(
            err.contains("'Pt'") && err.contains("Pt_LD") && err.contains("Au_JC"),
            "{err}"
        );
    }
}
