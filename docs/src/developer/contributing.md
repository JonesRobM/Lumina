# Contributing

## Getting Started

1. Fork the repository and clone your fork.
2. Install Rust 1.75+ via [rustup](https://rustup.rs).
3. Run `cargo test --workspace --exclude lumina-gui` to ensure everything passes.
4. (Optional) Run `cargo test --workspace --exclude lumina-gui --features gpu` to include GPU tests.

## Code Style

- Follow standard `rustfmt` formatting (`cargo fmt`).
- Run `cargo clippy --workspace` (and `cargo clippy --workspace --features gpu` if working on GPU code) and address all warnings.
- All public functions must include documentation comments.
- Use British English spelling in documentation and comments.

## Testing

- Unit tests live alongside the code in `#[cfg(test)]` modules.
- Integration tests live in the workspace-level `tests/` directory.
- Property-based tests (using `proptest`) are encouraged for mathematical functions.

## Pull Requests

- Keep PRs focused on a single change.
- Include tests for new functionality.
- Update documentation if the public API changes.
