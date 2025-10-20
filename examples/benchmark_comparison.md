# Benchmark Comparison Example

This document shows how to run and interpret benchmark comparisons.

## Quick Start

### 1. Run basic benchmarks (rug only)
```bash
cargo bench
```

This will benchmark the default `rug` backend.

### 2. Compare with pure Rust alternatives
```bash
cargo bench --features bench-num-bigint,bench-ibig
```

This compares:
- rug (GMP-based)
- num-bigint (pure Rust)
- ibig (pure Rust)

### 3. Full comparison including all GMP backends
```bash
cargo bench --features bench-num-bigint,bench-ibig,bench-rust-gmp
```

This includes all available backends.

## Understanding Results

### Mulmod (Modular Multiplication)
Lower is better. Example output:
```
mulmod/rug              time:   [241.20 ns ...]
mulmod/num-bigint       time:   [418.34 ns ...]
mulmod/ibig             time:   [361.86 ns ...]
mulmod/rust-gmp         time:   [246.56 ns ...]
```

**Interpretation**: rug and rust-gmp are fastest (GMP-based), ibig is ~1.5x slower, num-bigint is ~1.7x slower.

### Powmod (Modular Exponentiation)
Lower is better. Example output:
```
powmod/rug              time:   [2.6124 µs ...]
powmod/num-bigint       time:   [21.021 µs ...]
powmod/ibig             time:   [6.6168 µs ...]
powmod/rust-gmp         time:   [2.6221 µs ...]
```

**Interpretation**: GMP backends are ~8x faster than num-bigint, ~2.5x faster than ibig for this operation.

## Discrete Log Benchmarks

To focus on discrete log algorithm performance:
```bash
cargo bench --bench discrete_log
```

This tests:
- trial_mul (exhaustive search)
- shanks_steps (baby-step giant-step)
- pollard_rho (Pollard's rho)
- pohlig_hellman (Pohlig-Hellman)
- Auto algorithm selection

## Advanced: HTML Reports

Criterion generates HTML reports in `target/criterion/`:
```bash
cargo bench
# Then open target/criterion/report/index.html in a browser
```

## Performance Tips

1. **For maximum accuracy**: Close other applications and run benchmarks multiple times
2. **For quick comparisons**: Use `--quick` flag: `cargo bench -- --quick`
3. **For specific benchmarks**: Use filtering: `cargo bench mulmod`

## Platform-Specific Notes

### Linux
- Install GMP: `sudo apt-get install libgmp-dev`
- Best performance on x86_64

### macOS
- Install GMP: `brew install gmp`
- M1/M2 Macs show excellent performance

### Windows
- GMP installation can be complex
- Consider pure Rust alternatives (ibig, num-bigint)
- See rug documentation for Windows setup

### WebAssembly
- GMP backends (rug, rust-gmp) won't work
- Use `ibig` or `num-bigint`
- ibig generally performs better

## Example: Choosing a Backend for Your Project

### High-Performance Server
```toml
[dependencies]
discrete-logarithm = "1.0"
# Uses rug by default - fastest option
```

### WebAssembly Target
```toml
[dependencies]
discrete-logarithm = { version = "1.0", default-features = false }
# You would need to modify the library to support ibig/num-bigint as backend
# Currently only benchmarks support multiple backends
```

## Note on Library vs Benchmarks

**Important**: The discrete-logarithm library itself uses `rug` as its bignum backend. The optional features (`bench-num-bigint`, `bench-ibig`, `bench-rust-gmp`) are **only for benchmarking** - they don't change what the library uses internally.

If you need a discrete logarithm library using a different bignum backend, you would need to:
1. Fork the library
2. Replace `rug::Integer` with your preferred backend
3. Update all algorithms to use the new API

The benchmarks demonstrate the performance characteristics of different bignum libraries to help inform such decisions.
