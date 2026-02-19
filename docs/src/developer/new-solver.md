# Adding a New Solver

Lumina is designed to support multiple electromagnetic simulation methods beyond the CDA. This guide explains how to add a new solver.

## Step 1: Implement `OpticalSolver`

Create a new module under `lumina-core/src/solver/` and implement the trait:

```rust
pub trait OpticalSolver {
    fn compute_cross_sections(&self, dipoles: &[Dipole], wavelength_nm: f64, params: &SimulationParams) -> Result<CrossSections, SolverError>;
    fn solve_dipoles(&self, dipoles: &[Dipole], wavelength_nm: f64, params: &SimulationParams) -> Result<DipoleResponse, SolverError>;
    fn compute_near_field(&self, dipoles: &[Dipole], response: &DipoleResponse, plane: &NearFieldPlane) -> Result<NearFieldMap, SolverError>;
    fn method_name(&self) -> &str;
}
```

## Step 2: Register the Solver

Add the new solver to the solver selection logic in the CLI and GUI so that users can choose it via configuration.

## Step 3: Validate

Every new solver should be validated against known analytical solutions or established numerical codes before use.
