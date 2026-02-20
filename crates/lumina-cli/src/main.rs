//! Lumina command-line interface.
//!
//! Run simulations from TOML configuration files:
//! ```sh
//! lumina-cli run job.toml
//! lumina-cli validate job.toml
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

            let spectra = runner::run_simulation(&job)?;

            // Determine output path
            let out_dir = output
                .unwrap_or_else(|| PathBuf::from(&job.output.directory));
            let csv_path = out_dir.join("spectra.csv");

            runner::write_spectra_csv(&spectra, &csv_path)?;
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
            println!("  Au_JC  - Gold (Johnson & Christy, 400-800 nm)");
            println!("  Ag_JC  - Silver (Johnson & Christy) [placeholder]");
            println!("  Cu_JC  - Copper (Johnson & Christy) [placeholder]");
            Ok(())
        }
    }
}
