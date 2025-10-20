# Discrete Logarithm Solver

[![Build](https://github.com/skyf0l/discrete-logarithm/actions/workflows/ci.yml/badge.svg)](https://github.com/skyf0l/discrete-logarithm/actions/workflows/ci.yml)
[![Crate.io](https://img.shields.io/crates/v/discrete-logarithm.svg)](https://crates.io/crates/discrete-logarithm)
[![codecov](https://codecov.io/gh/skyf0l/discrete-logarithm/branch/main/graph/badge.svg)](https://codecov.io/gh/skyf0l/discrete-logarithm)

Fast discrete logarithm solver in Rust.

The code is based on the [sympy](https://github.com/sympy/sympy) implementation and translated to Rust.

Based on [rug](https://crates.io/crates/rug), it can use [arbitrary-precision numbers (aka BigNum)](https://en.wikipedia.org/wiki/Arbitrary-precision_arithmetic).

## Algorithm

This library solves the discrete logarithm problem: given `b`, `a`, and `n`, find the smallest non-negative integer `x` such that `b^x ≡ a (mod n)`.

The main `discrete_log` function intelligently selects the most efficient algorithm based on the characteristics of the input, specifically the order of the group. The following algorithms are implemented:

| Algorithm | Complexity | Memory | Use Case | Description |
|-----------|------------|--------|----------|-------------|
| **Trial Multiplication** | O(order) | O(1) | Very small orders (< 1,000) | Exhaustive search testing each exponent sequentially |
| **Baby-Step Giant-Step** | O(√order) | O(√order) | Prime orders when memory usage is acceptable | Time-memory tradeoff algorithm that precomputes a table of values |
| **Pollard's Rho** | O(√order) | O(1) | Large prime orders where memory is constrained | Randomized algorithm with minimal memory requirements, same expected time as Shanks |
| **Pohlig-Hellman** | O(∑ e_i(log(n) + √p_i)) | O(log order) | Composite orders (non-prime) | Reduces the problem to smaller subproblems using the factorization of the group order |
| **Index Calculus** | O(exp(2√(log(n)log(log(n))))) | O(B) | Very large prime orders where exp(2√(log(n)log(log(n)))) < √order | Most efficient for very large primes, uses smooth numbers and linear algebra |

### Algorithm Selection Logic

The library automatically selects the optimal algorithm:

1. If order < 1,000: use **Trial Multiplication**
2. If order is prime (or probably prime):
   - If 4√(log(n)log(log(n))) < log(order) - 10: use **Index Calculus**
   - Else if order < 10^12: use **Baby-Step Giant-Step**
   - Else: use **Pollard's Rho**
3. If order is composite: use **Pohlig-Hellman**

This automatic selection ensures optimal performance across different problem sizes and characteristics.

## License

Licensed under either of

- Apache License, Version 2.0
  ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
- MIT license
  ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
