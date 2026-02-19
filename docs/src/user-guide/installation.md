# Installation

## Prerequisites

- **Rust** 1.75 or later ([rustup.rs](https://rustup.rs))
- A C/C++ compiler (MSVC on Windows, gcc/clang on Linux/macOS)

## Building from Source

```bash
git clone https://github.com/YOUR_USERNAME/Lumina.git
cd Lumina
cargo build --release
```

## Running the GUI

```bash
cargo run --release -p lumina-gui
```

## Running the CLI

```bash
cargo run --release -p lumina-cli -- run examples/gold_sphere.toml
```

## Optional Features

GPU compute support (requires a compatible GPU and drivers):

```bash
cargo build --release --features lumina-compute/gpu
```
