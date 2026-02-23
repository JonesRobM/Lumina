# CLI Usage

The Lumina CLI (`lumina-cli`) enables headless simulation runs from TOML configuration files.

## Commands

### `run`

Execute a simulation:

```bash
lumina-cli run job.toml
```

### `validate`

Check a configuration file for errors without running:

```bash
lumina-cli validate job.toml
```

### `materials`

List available built-in materials:

```bash
lumina-cli materials
```

## GPU Acceleration

When built with `--features gpu`, the CLI supports GPU-accelerated GMRES:

```bash
cargo run --release -p lumina-cli --features gpu -- run job.toml
```

Control the backend via the TOML configuration:

```toml
[simulation]
backend = "auto"   # "cpu", "gpu", or "auto" (default)
```

- `"auto"` — tries GPU, falls back to CPU if unavailable
- `"cpu"` — always use CPU
- `"gpu"` — require GPU (errors if unavailable)

If the binary was built without `--features gpu`, the `backend` field is ignored and CPU is always used.

## Output Formats

The CLI writes results to the directory specified in `[output]`:

| File | Description | Controlled by |
|------|-------------|---------------|
| `spectra.csv` | Wavelength, C\_ext, C\_abs, C\_sca (and CD if computed) | `save_spectra = true` (default) |
| `spectra.json` | Same data in JSON format | `save_json = true` |
| `near_field/near_field_xy_peak.csv` | Near-field \\(\|E\|^2 / \|E_0\|^2\\) at peak extinction | `save_near_field = true` |

## Environment Variables

| Variable | Description |
|----------|-------------|
| `RUST_LOG` | Set logging level (`error`, `warn`, `info`, `debug`, `trace`) |
| `LUMINA_THREADS` | Override the number of CPU threads (default: all available) |
