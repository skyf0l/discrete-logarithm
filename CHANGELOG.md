# Changelog

## 2.0.0

### Breaking

- `Error` is `#[non_exhaustive]`: a later version can add a variant without that being a breaking
  change, and a match on it now needs a catch-all arm.
- `Error` has three new variants.
  - `InvalidModulus`. A modulus below 1 returns it instead of `NotRelativelyPrime` (`n_order`) or
    `LogDoesNotExist` (`discrete_log_with_order`). `discrete_log_index_calculus`, which used to
    validate nothing, returns it too.
  - `InvalidOrder`. An order that is not positive, or a factorization that cannot be the prime
    factorization of the number it is given for, returns it instead of panicking or running without
    end: see the fixes below.
  - `OutOfReach`, for an instance beyond what the chosen algorithm can do, where `LogDoesNotExist`
    used to claim that no logarithm existed. `discrete_log_shanks_steps` returns it for an order at
    or above `MAX_ORDER`, and `discrete_log_index_calculus` for a modulus whose relation matrix
    would be above its memory budget.
- `LogDoesNotExist` from `discrete_log_pollard_rho` and `discrete_log_index_calculus` is documented
  for what it is: no logarithm found within the budget of the search, not a proof that none exists.
- The algorithm selection tests the bounds on the order before the index calculus comparison. A
  prime order below `2**31` now goes to baby-step giant-step, where the comparison alone also chose
  index calculus for the small moduli (`n` from about 2000 to 1.16 million), 3.6 to 6.6 times slower
  there than the baby-step giant-step it displaced.
- Modulo 1 every residue is 0, so every algorithm returns `Ok(0)` for `n == 1`.
  `discrete_log_pollard_rho(1, 0, 0, Some(1000))` returned an arbitrary exponent and
  `discrete_log_index_calculus(1, ...)` an error.
- `discrete_log_pollard_rho_parallel` clamps `threads` to the new `MAX_THREADS`, and runs one of the
  walks in the calling thread, so `threads` walks cost `threads - 1` spawned threads.
- `discrete_log` and `discrete_log_with_factors` no longer reject a base that shares a factor with
  the modulus. They compute the order of the group and search, as sympy does, so
  `discrete_log(9, 0, 3)` returns 2 where it used to return an error.
- `discrete_log_pollard_rho(227, 3**7, 5)` returns 132 instead of failing: the logarithm exists,
  and the walk used until now could not find it.
- The minimum supported Rust version is declared: 1.86, and the crate is on edition 2024.

### Added

- `MAX_THREADS`, the cap `discrete_log_pollard_rho_parallel` clamps its thread count to (`parallel`
  feature).
- `discrete_log_with_prime_order`, for an order already known to be prime.
- `discrete_log_pohlig_hellman_with_factors`, for an order whose factorization is known.
- `n_order_with_factors` is now exported.
- `discrete_log_pollard_rho_with_seed` and `discrete_log_index_calculus_with_seed`: the randomized
  algorithms with a seed, so a run can be reproduced or retried differently.
- Optional `parallel` feature, adding `discrete_log_pollard_rho_parallel`, the parallel collision
  search of van Oorschot and Wiener. Off by default, `std` threads only.
- A benchmark suite (`benches/`) measuring instruction counts, and the differential test against
  sympy that gates the results.

### Fixed

- Factoring no longer assumes that what trial division leaves behind is prime. Moduli with two
  prime factors above 15.5 million gave wrong group orders, failed logarithms, and could send
  Pohlig-Hellman into unbounded recursion.
- Index calculus no longer fails on small instances whose relation budget was spent on relations
  that carried no new information, and no longer loops forever when the target is 0.
- Pollard's rho no longer panics for an order below 4 and draws its exponents from the whole range.
- Trial multiplication no longer overflows its counter above 2^31 steps. It also validates the
  modulus like every other entry point (`discrete_log_trial_mul(0, ...)` returned a panicking
  division by zero, a negative modulus a meaningless answer), answers 0 modulo 1, and treats a
  negative `order` as an empty search instead of walking the whole of a `usize`.
- Baby-step giant-step returns `NotRelativelyPrime` instead of panicking when the base cannot be
  inverted.
- `n_order` of an order that is a perfect power of a large prime now works, and `n_order(a, 1)`
  returning 1 is documented.
- The entry points taking a factorization now check it, and return `InvalidOrder` instead of
  panicking, hanging or allocating without bound. `discrete_log_pohlig_hellman_with_factors` panicked
  with a division by zero on a prime of 0 and with "a non negative exponent" on a negative order or
  a negative prime; `n_order_with_factors` and `discrete_log_with_factors` panicked on an exponent of
  0 (a subtraction that overflowed, or `3**4294967295` in release), ran for minutes on an exponent of
  `usize::MAX`, and looped for ever while growing the order on a factorization belonging to another
  number (`n_order_with_factors(3, 4, {3: 1})`), the base not being a unit of a prime power that does
  not divide the modulus. Exponents are also refused before anything casts them to a `u32`, which
  used to truncate silently.
- `discrete_log_pollard_rho_parallel` no longer panics on a large thread count: 100 000 threads was
  "failed to spawn thread: WouldBlock" and `usize::MAX` a capacity overflow. A spawn the operating
  system refuses now leaves the search to the walks already started.
- `discrete_log_index_calculus` stops the walk that looks for the first relation as soon as it
  reaches zero, which is where a base that is not invertible sends it: `(4, 3, 0, Some(10**12))`
  returned an error after hours of multiplying zero by zero, and now returns at once.
- `discrete_log_index_calculus` refuses a modulus whose dense relation matrix would be above 32 MiB
  with `OutOfReach`, which is from about 88 bits on. `discrete_log` sends a prime modulus of 128 bits
  or more with a prime order to index calculus, whose factor base (17 859 primes, 2.4 GiB of rows at
  128 bits, 49 GiB at 160) is far beyond memory: it used to spin over candidates and never return,
  and now fails in under a millisecond. A modulus of 72 bits, 836 primes and 5.3 MiB of rows, is
  still solved.
- Brent's cycle length in the factorization no longer wraps a `u64` after 63 doublings, which turned
  the walk into a loop of no steps that never ended (an overflow panic in debug).
- `ModRing::pow` refuses a negative exponent the same way in both rings, and in release builds as in
  debug ones: the word-size ring used to raise to the absolute value of it where GMP panicked.

### Performance

Against 1.1.0, on the same machine: small and medium moduli 500x to 8500x faster (a real
factoriser with an early exit instead of dividing by the first million primes every time), a
2^47 prime order 31x (word-size Montgomery arithmetic, an r-adding walk with Brent cycle
detection), the 108 digit smooth case 130x (Carmichael lambda, one pass for the order and its
factorization, a shared table of baby steps per prime), index calculus 6.8x (ring arithmetic, the
single large prime variation, walk-based candidate generation). The algorithm selection thresholds
are now measured on this implementation instead of inherited.

## 1.1.0 and earlier

No changelog was kept.
