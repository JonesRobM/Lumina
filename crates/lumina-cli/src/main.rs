//! Lumina command-line interface.
//!
//! Run simulations from TOML configuration files:
//! ```sh
//! lumina-cli run job.toml
//! lumina-cli validate job.toml
//! lumina-cli materials
//! ```

mod config;
mod runner;

use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "lumina-cli")]
#[command(about = "Lumina: Coupled Dipole Approximation Framework")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run a simulation from a TOML configuration file.
    Run {
        /// Path to the job configuration file.
        config: PathBuf,
        /// Output directory (overrides config file setting).
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Validate a configuration file without running the simulation.
    Validate {
        /// Path to the job configuration file.
        config: PathBuf,
    },
    /// Display information about available materials.
    Materials,
}

fn main() -> anyhow::Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Run { config, output } => {
            println!("Lumina CDA Solver");
            println!("=================");
            let job = config::load_config(&config)?;
            println!("Configuration: {}", config.display());

            let result = runner::run_simulation(&job)?;

            // Determine output directory
            let out_dir = output.unwrap_or_else(|| PathBuf::from(&job.output.directory));

            // CSV spectra (default on)
            if job.output.save_spectra {
                let csv_path = out_dir.join("spectra.csv");
                runner::write_spectra_csv(&result.spectra, &csv_path, &job)?;
            }

            // JSON spectra (optional)
            if job.output.save_json {
                let json_path = out_dir.join("spectra.json");
                runner::write_spectra_json(&result.spectra, &json_path)?;
            }

            // Near-field map (optional)
            if let Some(nf) = &result.near_field {
                let nf_path = out_dir.join("near_field.csv");
                runner::write_near_field_csv(nf, &nf_path)?;
            }

            println!("Simulation complete.");
            Ok(())
        }
        Commands::Validate { config } => {
            let _job = config::load_config(&config)?;
            println!("Configuration is valid: {}", config.display());
            Ok(())
        }
        Commands::Materials => {
            println!("Available materials:");
            println!();
            println!("  Johnson & Christy (1972) metals:");
            println!("    Au_JC      — Gold,   188–892 nm");
            println!("    Ag_JC      — Silver, 188–892 nm");
            println!("    Cu_JC      — Copper, 188–892 nm");
            println!();
            println!("  Palik Handbook dielectrics:");
            println!("    TiO2_Palik — Rutile TiO₂, 300–1000 nm");
            println!("    SiO2_Palik — Fused silica SiO₂, 300–1000 nm");
            Ok(())
        }
    }
}
