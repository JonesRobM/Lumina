# Installation

## Prerequisites

- **Rust** 1.75 or later ([rustup.rs](https://rustup.rs))
- A C/C++ compiler (MSVC on Windows, gcc/clang on Linux/macOS)

### GPU acceleration (optional)

To use GPU-accelerated GMRES, you also need:

- A GPU with Vulkan, DirectX 12, or Metal support
- Up-to-date GPU drivers (NVIDIA 535+, AMD Adrenalin 23.7+, or macOS 13+)

Most discrete GPUs from 2018 onwards and integrated GPUs from 2020 onwards are supported via wgpu.

## Building from Source

### CPU only (default)

```bash
git clone https://github.com/YOUR_USERNAME/Lumina.git
cd Lumina
cargo build --release
```

### With GPU acceleration

```bash
cargo build --release --features gpu
```

This enables the `gpu` feature across the workspace (`lumina-compute`, `lumina-core`, `lumina-gui`, and `lumina-cli`).

## Running the GUI

```bash
# CPU only
cargo run --release -p lumina-gui

# With GPU toggle enabled
cargo run --release -p lumina-gui --features gpu
```

When built with `--features gpu`, the simulation panel shows a **GPU acceleration** checkbox.

## Running the CLI

```bash
# CPU only
cargo run --release -p lumina-cli -- run examples/gold_sphere.toml

# With GPU support
cargo run --release -p lumina-cli --features gpu -- run examples/gold_sphere.toml
```

The CLI reads the `backend` field from the TOML `[simulation]` section:

```toml
[simulation]
backend = "auto"   # "cpu", "gpu", or "auto" (default)
```

- `"auto"` — uses GPU if available, falls back to CPU
- `"cpu"` — always CPU
- `"gpu"` — requires GPU (errors if unavailable)

## Running Tests

```bash
# All tests (CPU)
cargo test --workspace --exclude lumina-gui

# Including GPU tests
cargo test --workspace --exclude lumina-gui --features gpu

# GPU benchmark (release mode, with output)
cargo test -p lumina-core --features gpu --release -- gpu_benchmark --nocapture
```

## Troubleshooting

### GPU not detected

If GPU initialisation fails, Lumina falls back to CPU automatically (when `backend = "auto"`). Check:

1. **Drivers:** Run `vulkaninfo` (Linux/Windows) or check System Preferences > GPU (macOS) to verify driver support.
2. **wgpu compatibility:** Lumina uses wgpu v24, which requires Vulkan 1.1+ (Linux/Windows) or Metal (macOS).
3. **Headless servers:** GPU compute works without a display — wgpu can use headless Vulkan.
