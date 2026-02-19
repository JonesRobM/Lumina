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

## Environment Variables

| Variable | Description |
|----------|-------------|
| `RUST_LOG` | Set logging level (`error`, `warn`, `info`, `debug`, `trace`) |
| `LUMINA_THREADS` | Override the number of CPU threads (default: all available) |
