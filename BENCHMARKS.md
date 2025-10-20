# Bignum Library Benchmarks

This document provides information about benchmarking different bignum libraries with the discrete-logarithm crate.

## Supported Libraries

The following bignum libraries are supported for benchmarking:

### GMP-based Libraries
- **rug** (default): Rust bindings to GMP via MPFR. Provides comprehensive arbitrary-precision arithmetic with excellent performance.
- **rust-gmp**: Direct Rust bindings to GMP. Slightly simpler API than rug, similar performance.

### Pure Rust Libraries
- **num-bigint**: Popular pure Rust implementation. Portable but slower than GMP-based libraries.
- **ibig**: Modern pure Rust big integer library with focus on performance. Faster than num-bigint.

### Not Included
- **ramp**: Requires nightly Rust due to use of unstable features. Cannot be used with stable Rust.

## Running Benchmarks

### Default Benchmarks (rug only)
```bash
cargo bench
```

### With Specific Backend
```bash
# Benchmark with num-bigint
cargo bench --features bench-num-bigint

# Benchmark with ibig
cargo bench --features bench-ibig

# Benchmark with rust-gmp
cargo bench --features bench-rust-gmp
```

### Compare All Backends
```bash
cargo bench --features bench-num-bigint,bench-ibig,bench-rust-gmp
```

## Benchmark Suites

### 1. Basic Operations (`bignum_ops`)
Tests fundamental bignum operations:
- **mulmod**: Modular multiplication `(g * h) mod p`
- **powmod**: Modular exponentiation `g^exp mod p`

### 2. Discrete Logarithm Algorithms (`discrete_log`)
Tests the actual discrete logarithm algorithms:
- **trial_mul**: Exhaustive search (for small orders)
- **shanks_steps**: Baby-step giant-step algorithm
- **pollard_rho**: Pollard's rho algorithm
- **pohlig_hellman**: Pohlig-Hellman algorithm

## Performance Notes

### Expected Performance Characteristics

1. **GMP-based libraries (rug, rust-gmp)**:
   - Fastest overall performance
   - Highly optimized C library underneath
   - Best for production use with large numbers
   - Requires system GMP installation

2. **Pure Rust libraries (ibig, num-bigint)**:
   - Portable, no external dependencies
   - Slower than GMP-based libraries
   - ibig is typically 2-3x faster than num-bigint
   - Good for WebAssembly and embedded targets

### Sample Results

Based on typical benchmark runs on modern hardware:

**Modular Multiplication (mulmod)**:
- rug: ~240 ns
- rust-gmp: ~245 ns
- ibig: ~360 ns
- num-bigint: ~410 ns

**Modular Exponentiation (powmod, exp=65537)**:
- rug: ~2.6 µs
- rust-gmp: ~2.6 µs
- ibig: ~6.6 µs
- num-bigint: ~21 µs

## Choosing a Backend

### Use rug (default) if:
- You need maximum performance
- You're building for desktop/server environments
- GMP dependencies are acceptable
- You want the most comprehensive API

### Use rust-gmp if:
- You need GMP performance with a simpler API
- You prefer direct GMP bindings
- Performance is critical

### Use ibig if:
- You need pure Rust implementation
- You're targeting WebAssembly
- You need portable code without external dependencies
- Performance is still important but not critical

### Use num-bigint if:
- You need pure Rust implementation
- You want the most mature pure Rust option
- Performance is less critical than ecosystem compatibility
- You're already using other num-* crates

## Contributing Benchmark Results

If you run benchmarks on your system, consider sharing the results! This helps the community understand performance characteristics across different platforms and architectures.

### Platform Information to Include:
- CPU model and architecture
- Operating system
- Rust version
- GMP version (for rug/rust-gmp)
- Benchmark command used
- Full benchmark output

## Implementation Details

### Core Library
The discrete-logarithm library itself uses **rug** as its core bignum implementation. The benchmarks allow comparison of different backends but the library code is not abstracted to support multiple backends at runtime.

### Why Not Multi-Backend?
Supporting multiple bignum backends at runtime would require:
1. Significant abstraction layer overhead
2. Loss of type safety and zero-cost abstractions
3. Increased maintenance burden
4. Reduced performance due to dynamic dispatch

The current approach provides benchmarks for comparison while keeping the library focused on a single, high-performance backend (rug).

### Future Considerations
If there's strong demand for supporting alternative backends in the library itself (not just benchmarks), this could be reconsidered. However, it would likely be implemented as separate feature flags that select the backend at compile time, not runtime.

## See Also
- [Criterion.rs](https://github.com/bheisler/criterion.rs) - The benchmarking framework used
- [rug documentation](https://docs.rs/rug/)
- [num-bigint documentation](https://docs.rs/num-bigint/)
- [ibig documentation](https://docs.rs/ibig/)
- [rust-gmp documentation](https://docs.rs/rust-gmp/)
