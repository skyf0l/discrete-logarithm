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

| Algorithm | Complexity | Memory | Use Case |
|-----------|------------|--------|----------|
| **Trial Multiplication**<br>Exhaustive search testing each exponent sequentially | O(order) | O(1) | Very small orders (< 1,000) |
| **Baby-Step Giant-Step**<br>Time-memory tradeoff algorithm that precomputes a table of values | O(√order) | O(√order) | Prime orders when memory usage is acceptable |
| **Pollard's Rho**<br>Randomized algorithm with minimal memory requirements, same expected time as Shanks | O(√order) | O(1) | Large prime orders where memory is constrained |
| **Pohlig-Hellman**<br>Reduces the problem to smaller subproblems using the factorization of the group order | O(∑ e_i(log(n) + √p_i)) | O(log(order)) | Composite orders (non-prime) |
| **Index Calculus**<br>Most efficient for very large primes, uses smooth numbers and linear algebra | O(exp(2√(log(n)log(log(n))))) | O(B) | Very large prime orders where exp(2√(log(n)log(log(n)))) < √order |

### Algorithm Selection Logic

The library automatically selects the optimal algorithm:

1. If order < 1,000: use **Trial Multiplication**
2. If order is prime (or probably prime):
   - If 4√(log(n)log(log(n))) < log(order) - 10: use **Index Calculus**
   - Else if order < 10^12: use **Baby-Step Giant-Step**
   - Else: use **Pollard's Rho**
3. If order is composite: use **Pohlig-Hellman**

This automatic selection ensures optimal performance across different problem sizes and characteristics.

### Group Order

The order of `b` and its prime factorization are computed in one pass: the order of the group of
units (the totient of `n`) is built from the factorization of `n` and of every `p - 1`, then the
primes the order of `b` does not need are divided out. Pohlig-Hellman, which works on that
factorization, does not have to compute it again, and the order needs no primality test.

Factoring is done by trial division by the primes below 65536, then by Brent's variant of
Pollard's rho on the cofactor; the root of a perfect power is factored instead of the power
itself. Moduli with two prime factors above ~10^12 are out of reach, as factoring them is the
general factoring problem: pass what you know with `discrete_log_with_factors`,
`discrete_log_with_order` or `discrete_log_with_prime_order` to skip this step.

`n` must be positive, and modulo 1 the logarithm is always 0. Bases that share a factor with the
modulus are not rejected: their powers are not a subgroup of the units, so the search is bounded
by the order of the group instead.

The randomized algorithms (Pollard's rho, index calculus) have seeded variants
(`discrete_log_pollard_rho_with_seed`, `discrete_log_index_calculus_with_seed`): the same seed
always makes the same choices, and a retry with another seed makes other ones.

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
