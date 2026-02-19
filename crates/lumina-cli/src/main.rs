//! Lumina command-line interface.
//!
//! Run simulations from TOML configuration files:
//! ```sh
//! lumina-cli run job.toml
//! lumina-cli validate job.toml
//! ```

mod config;

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
        Commands::Run { config } => {
            let job = config::load_config(&config)?;
            println!("Loaded configuration: {:?}", config);
            println!("Simulation: {:?}", job.simulation);
            println!("Objects: {}", job.geometry.object.len());
            // TODO: Build geometry, create solver, run simulation, write output.
            println!("Simulation engine not yet implemented.");
            Ok(())
        }
        Commands::Validate { config } => {
            let _job = config::load_config(&config)?;
            println!("Configuration is valid: {:?}", config);
            Ok(())
        }
        Commands::Materials => {
            println!("Available materials:");
            println!("  Au_JC  - Gold (Johnson & Christy)");
            println!("  Ag_JC  - Silver (Johnson & Christy)");
            println!("  Cu_JC  - Copper (Johnson & Christy)");
            // TODO: List from material registry.
            Ok(())
        }
    }
}
