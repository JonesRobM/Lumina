# GUI Dashboard

The Lumina GUI provides an interactive environment for configuring, running, and visualising simulations. As of v0.2.1, wavelengths are solved in parallel across all CPU cores, making the GUI competitive with the CLI for large systems.

## Launching

```bash
# CPU only
cargo run --release -p lumina-gui

# With GPU acceleration toggle
cargo run --release -p lumina-gui --features gpu
```

> **Important:** Always use `--release` for production simulations. Debug mode is 10-50x slower for numerical workloads.

## Panels

### Geometry

Configure shapes (sphere, cylinder, cuboid, ellipsoid, helix), import files (.xyz, .obj), adjust dipole spacing, and preview the dipole lattice as a 2D scatter plot (XY/XZ/YZ projections).

#### Directory-Scoped File Browser (v0.2.1)

For `.xyz` and `.obj` imports, the geometry panel offers a folder-based workflow:

1. Click **Set folder...** to select a directory containing structure files
2. A dropdown populates with all `.xyz` and `.obj` files in that directory
3. Select a file from the dropdown to load it immediately
4. Click **Refresh** to re-scan the directory after adding or removing files

The traditional **Open .xyz file...** and **Open .obj file...** buttons remain available below the dropdown for one-off file picks.

### Materials

Select from built-in materials (Au, Ag, Cu via Johnson & Christy; TiO\(_2\), SiO\(_2\) via Palik) or define custom dielectric functions. View the interpolated \(\epsilon(\lambda)\) plot showing \(\epsilon_1\) and \(\epsilon_2\) curves.

### Simulation

Set the wavelength range, medium refractive index, and solver parameters. Additional controls:

- **Incident polarisation** — select x-pol, y-pol, or circular
- **Circular dichroism** — checkbox to compute \(\Delta C_\text{ext}\) (requires two solver calls per wavelength)
- **GPU acceleration** — checkbox to enable GPU matvec (only visible when built with `--features gpu`)
- **Debug output** — checkbox to show per-wavelength solver diagnostics (v0.2.1)

Launch the simulation and monitor progress via the progress bar. Wavelengths are solved in parallel using all available CPU cores.

#### Debug Output (v0.2.1)

When **Debug output** is enabled, the simulation panel displays a scrollable monospace log with:

- Backend type (CPU or GPU device name)
- Wavelength range, environment refractive index, and material
- For XYZ imports: auto-detected nearest-neighbour distance and FCD threshold
- Dipole count, matrix dimensions, and solver method (Direct LU or GMRES)
- Per-wavelength: \(\epsilon(\lambda)\), \(C_\text{ext}\), \(C_\text{abs}\), and solve time in milliseconds
- Total simulation time

The log updates in real time as wavelengths complete (note: completion order may differ from wavelength order due to parallel execution).

### Results

View interactive plots powered by `egui_plot`:

- **Spectra** — extinction, absorption, and scattering cross-sections vs wavelength (with optional CD curve)
- **Near-field heatmap** — \(|E|^2 / |E_0|^2\) on the XY plane at peak extinction, log-scale colourmap
- **Far-field polar plot** — E-plane and H-plane radiation pattern cuts

Export results to CSV or JSON via the file dialog (native file picker via `rfd`).

### Chat

An AI assistant panel powered by any OpenAI-compatible provider. Configure it under the **Settings** tab, then switch to the **Chat** tab to ask questions about your simulation.

#### Settings tab

| Field | Description |
|-------|-------------|
| Provider | Preset selector: Ollama (local), Groq, OpenRouter, or Custom |
| Base URL | Full versioned URL prefix (e.g. `https://api.groq.com/openai/v1`) |
| Model | Model identifier passed verbatim to the provider |
| API Key | Optional; stored in plaintext in the platform config directory |

Click **Save** to persist settings to `<config_dir>/lumina/chat.toml`.

#### Chat tab

- Type a message and press **Enter** to send (Shift+Enter inserts a newline).
- Tokens stream in real time; the input is disabled while a response is arriving.
- **Attach context** — injects the current scene geometry and simulation results as a structured message, then sends immediately. Use this before asking questions about your simulation.
- **Clear** — wipes the conversation history. Any in-flight stream completes silently in the background.

#### Free providers

| Provider | Free tier | Setup |
|----------|-----------|-------|
| [Ollama](https://ollama.com) | Fully free, runs locally | `ollama serve` + pull a model |
| [Groq](https://console.groq.com) | Free API key, generous rate limits | Sign up, copy key |
| [OpenRouter](https://openrouter.ai) | Free-tier models available | Sign up, copy key |

> **Note:** API keys are stored in plaintext. Avoid using paid keys with high blast radius; free-tier keys are the intended use case.
