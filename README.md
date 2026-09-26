# discrete-logarithm

[![crates.io](https://img.shields.io/crates/v/discrete-logarithm.svg)](https://crates.io/crates/discrete-logarithm)
[![docs.rs](https://img.shields.io/docsrs/discrete-logarithm)](https://docs.rs/discrete-logarithm)
[![MSRV](https://img.shields.io/crates/msrv/discrete-logarithm)](https://crates.io/crates/discrete-logarithm)
[![License](https://img.shields.io/crates/l/discrete-logarithm)](#license)

Fast discrete logarithms modulo any integer, in Rust, following [sympy](https://github.com/sympy/sympy)'s `ntheory` and built on [rug](https://crates.io/crates/rug) (GMP) for arbitrary-precision integers.

The algorithms, their selection and the exact results are sympy's; the implementations no longer are, see [Credits and references](#credits-and-references).

## Algorithm

This library solves the discrete logarithm problem: given `b`, `a`, and `n`, find the smallest non-negative integer `x` such that `b^x ≡ a (mod n)`.

The main `discrete_log` function intelligently selects the most efficient algorithm based on the characteristics of the input, specifically the order of the group. The following algorithms are implemented:

| Algorithm | Complexity | Memory | Use Case |
|-----------|------------|--------|----------|
| **Trial Multiplication**<br>Exhaustive search testing each exponent sequentially | O(order) | O(1) | Very small orders (< 1,000) |
| **Baby-Step Giant-Step**<br>Time-memory tradeoff, interleaving baby and giant steps so a small logarithm is found in O(√x) | O(√order) | O(√order) | Prime orders below 2^31, where its table stays around a megabyte |
| **Pollard's Rho**<br>Randomized r-adding walk (Teske) with Brent cycle detection, no memory at all | O(√order) | O(1) | Prime orders above 2^31, and any order when memory is constrained |
| **Pohlig-Hellman**<br>Reduces the problem to smaller subproblems using the factorization of the group order | O(∑ e_i(log(n) + √p_i)) | O(log(order)) | Composite orders (non-prime) |
| **Index Calculus**<br>Smooth relations and linear algebra, with the single large prime variation | O(exp(2√(log(n)log(log(n))))) | O((B/log B)²), capped at 32 MiB | Large prime orders, from about a 45 bit modulus on a safe prime |

### Algorithm Selection Logic

The library automatically selects the optimal algorithm:

1. If order < 1,000: use **Trial Multiplication**
2. If order is prime (or probably prime):
   - If order < 2^31: use **Baby-Step Giant-Step**
   - Else if 4√(log(n)log(log(n))) < log(order) + 11: use **Index Calculus**
   - Else: use **Pollard's Rho**
3. If order is composite: use **Pohlig-Hellman**

The two boundaries were measured on this implementation, not inherited: baby-step giant-step stops
at 2^31 because above it the cache misses of its table cost more than the extra instructions of
Pollard's rho (at 2^31 the two are within 5%, and rho needs no memory at all), and index calculus
takes over where it is measured about twice as cheap as rho, which on a safe prime is a modulus of
about 45 bits (three times as cheap at 48 bits, and the gap keeps growing). The index calculus
boundary is deliberately conservative: unlike the other algorithms it can give up,
so it is only chosen where it wins by a wide margin. 10^12 remains the memory cap of
`discrete_log_shanks_steps` itself, which refuses orders at or above it.

### Group Order

The order of `b` and its prime factorization are computed in one pass: Carmichael's lambda, the
exponent of the group of units, is built from the factorization of `n` and of every `p - 1`, then
the primes the order of `b` does not need are divided out. A base that shares a factor with the
modulus gets the Euler totient instead, since its powers are not a subgroup of the units and the
value is then only a search bound. Pohlig-Hellman, which works on that
factorization, does not have to compute it again, and the order needs no primality test.

Factoring is done by trial division by the primes below 65536, then by Brent's variant of
Pollard's rho on the cofactor; the root of a perfect power is factored instead of the power
itself. Two prime factors of about 2^50 each take a couple of seconds and 2^56 each about half a
minute, and past that see [What is out of reach](#what-is-out-of-reach): pass
what you know with `discrete_log_with_factors`,
`discrete_log_with_order` or `discrete_log_with_prime_order` to skip this step.

### Arithmetic

Moduli that fit in a machine word are computed with Montgomery multiplication on `u64`, with
`u128` intermediates and no allocation; larger ones use GMP through
[rug](https://crates.io/crates/rug), with the operations done in place so a loop allocates
nothing. Every algorithm above is written once against that layer and runs on either. This is
where most of the speed comes from: on word-size inputs the bignum library and the allocator used
to cost more than the algorithms themselves.

`n` must be positive, and modulo 1 the logarithm is always 0. Bases that share a factor with the
modulus are not rejected: their powers are not a subgroup of the units, so the search is bounded
by the order of the group instead. One exponent can escape that bound, since such a base can run
through `totient(n) + 1` distinct powers, and the last of them is then reported as having no
logarithm. sympy answers the same way, and the cases are rare (3,503 of the 180,000 base and
modulus pairs below 600).

The randomized algorithms (Pollard's rho, index calculus) have seeded variants
(`discrete_log_pollard_rho_with_seed`, `discrete_log_index_calculus_with_seed`): the same seed
always makes the same choices, and a retry with another seed makes other ones.

The optional `parallel` feature adds `discrete_log_pollard_rho_parallel`, the parallel collision
search of van Oorschot and Wiener: several threads walk independently and meet at distinguished
points. Measured at a 2^47 order, it is about 3x faster on 2 threads and 5x to 8x on 8 threads;
the length of a single walk is long-tailed, so those are averages over many instances, and part of
the gain on few threads is that distinguished points stop a walk earlier than cycle detection
does. It is opt-in because it is the only entry point that spawns threads, and it needs no
dependency beyond `std`.

### What is out of reach

Three limits are worth knowing before reaching for a cryptographic size:

- A large prime order modulo a modulus above about 88 bits is refused immediately with
  `Error::OutOfReach`: index calculus is the only algorithm that could take it, and its relation
  matrix would pass the 32 MiB the crate is willing to spend (2.1 GB at 128 bits, 125 TB at 256).
  A smaller prime order modulo the same modulus goes to Pollard's rho instead, which has no such
  limit but needs about √order group operations.
- Index calculus can also give up on an instance it does accept, returning
  `Error::LogDoesNotExist` for a logarithm that exists. It is the only algorithm here that can,
  which is why it is chosen conservatively; no give-up was observed in about 3,000 instances from
  30 to 56 bits.
- A modulus that is the product of two prime factors of about 2^56 or more **does not return**.
  Factoring it is the general factoring problem, and the factoriser has no time budget: it is not
  an error, the call runs until the process is stopped. Pass a known factorization or a known
  order to skip it.

## Credits and references

The crate began as a translation of sympy's `ntheory.residue_ntheory` (BSD-3-Clause, SymPy
Development Team), and still follows it where it is visible to callers: the same algorithms, the
same selection rules, the same logarithms and the same corner cases, including bases that share a
factor with the modulus. Correctness is checked against sympy directly, on 75,610 differential
cases (every base and target for every modulus up to 60, then random moduli up to 10^11).

The implementations are no longer translations. What they follow instead:

- P. L. Montgomery, *Modular multiplication without trial division* (1985), for the word-size
  arithmetic every algorithm runs on.
- E. Teske, *Speeding up Pollard's rho method for computing discrete logarithms* (1998), for the
  r-adding walk, which replaces the three-branch walk with a squaring step.
- R. P. Brent, *An improved Monte Carlo factorization algorithm* (1980), for the cycle detection
  and for the rho used by the factoriser.
- P. C. van Oorschot and M. J. Wiener, *Parallel collision search with cryptanalytic
  applications* (1999), for the optional parallel search.
- The interleaved form of baby-step giant-step, which grows both sequences together so that a
  small logarithm costs O(sqrt(x)) instead of O(sqrt(order)).
- The single large prime variation, standard in the factoring and index calculus literature, for
  the relation search of index calculus.
- A. J. Menezes, P. C. van Oorschot and S. A. Vanstone, *Handbook of Applied Cryptography*
  (1997), the reference sympy itself cites for these algorithms.

The selection thresholds are measured on this implementation rather than inherited, with the
benchmarks in [benches](benches/README.md).

## License

Licensed under either of

- Apache License, Version 2.0
  ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
- MIT license
  ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.

The parts listed in [LICENSE-BSD-SYMPY](LICENSE-BSD-SYMPY) began as a translation of sympy and
remain subject to its 3-clause BSD licence, whose notice that file keeps.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
